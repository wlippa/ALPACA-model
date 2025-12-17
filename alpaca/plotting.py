from pathlib import Path
from importlib import resources
import warnings
import argparse
import math

import pandas as pd
import numpy as np
import seaborn as sns
import plotly.graph_objects as go
import networkx as nx
from plotly.subplots import make_subplots

from alpaca.utils import read_tree_json, flat_list
from alpaca.plotting_helpers import (
    get_chr_table,
    remove_duplicates_preserve_order,
    get_unique_lists,
    get_tree_edges,
    clean_output,
    find_parent,
)


_DEFAULT_CHR_TABLE = "hg19_chromosome_lengths.csv"


def _format_rgb(rgb_tuple):
    return f'rgb{tuple(int(channel) for channel in rgb_tuple)}'


def build_copy_number_palette(max_state):
    """Generate a discrete palette with 0=blue, 1=white, >=2 shades of red."""
    if max_state is None or math.isnan(max_state):
        max_state = 0
    max_state = int(math.ceil(max(0, max_state)))

    palette = {0: (0, 0, 255)}
    if max_state >= 1:
        palette[1] = (255, 255, 255)

    red_states = max(0, max_state - 1)
    if red_states:
        reds = sns.color_palette("Reds", red_states)
        for idx, color in enumerate(reds, start=2):
            palette[idx] = tuple(int(round(channel * 255)) for channel in color)
    return palette


def load_chr_table(custom_path=None):
    """Load the chromosome table from either a user supplied path or the packaged resource."""
    if custom_path:
        return get_chr_table(custom_path)

    package_candidate = None
    if hasattr(resources, "files"):
        try:
            package_candidate = resources.files("alpaca.resources").joinpath(_DEFAULT_CHR_TABLE)
        except ModuleNotFoundError:
            package_candidate = None
    if package_candidate and package_candidate.is_file():
        with package_candidate.open("rb") as handle:
            return get_chr_table(handle)
    if not hasattr(resources, "files"):
        try:
            with resources.path("alpaca.resources", _DEFAULT_CHR_TABLE) as resource_path:
                return get_chr_table(resource_path)
        except (FileNotFoundError, ModuleNotFoundError):
            pass

    fallback = Path(__file__).resolve().parents[1] / "resources" / _DEFAULT_CHR_TABLE
    if fallback.exists():
        return get_chr_table(fallback)

    raise FileNotFoundError(
        "Unable to locate chromosome length table. Provide a path via --chr-table."
    )


def _read_table(table_path):
    sep = "\t" if str(table_path).lower().endswith((".tsv", ".txt")) else ","
    return pd.read_csv(table_path, sep=sep)


def load_mutation_table(input_dir, explicit_path=None):
    """Read the driver mutation table if it is present."""
    if explicit_path:
        path = Path(explicit_path)
        if not path.exists():
            raise FileNotFoundError(f"Mutation table not found: {path}")
        return _read_table(path)

    input_dir = Path(input_dir)
    candidates = (
        "driver_mutations.csv",
        "driver_mutations.tsv",
        "mutations.csv",
        "mutations.tsv",
    )
    for candidate in candidates:
        candidate_path = input_dir / candidate
        if candidate_path.exists():
            return _read_table(candidate_path)
    return None


def prepare_driver_mutations(mutation_df, tumour_id, chr_table):
    """Return driver mutations annotated with absolute coordinates for the selected tumour."""
    if mutation_df is None or mutation_df.empty:
        return None

    driver_df = mutation_df.copy()
    if "tumour_id" in driver_df.columns:
        driver_df = driver_df[driver_df["tumour_id"] == tumour_id]
        if driver_df.empty:
            return None

    if "clone" not in driver_df.columns:
        warnings.warn("Driver mutation table is missing 'clone' column; skipping mutation overlay.")
        return None

    if "abs_position" not in driver_df.columns:
        chr_column = next((c for c in ("chr", "chromosome", "chrom") if c in driver_df.columns), None)
        pos_column = next((c for c in ("position", "pos", "start", "bp") if c in driver_df.columns), None)
        if chr_column is None or pos_column is None:
            warnings.warn(
                "Driver mutation table requires chromosome and position columns to compute absolute coordinates; skipping mutation overlay."
            )
            return None
        chr_lookup = chr_table.set_index("chr")["shift"]
        driver_df["chr"] = driver_df[chr_column].astype(str).apply(
            lambda val: val if val.startswith("chr") else f"chr{val}"
        )
        driver_df["abs_position"] = driver_df["chr"].map(chr_lookup) + driver_df[pos_column].astype(float)
        driver_df = driver_df.dropna(subset=["abs_position"])
        if driver_df.empty:
            return None

    return driver_df.reset_index(drop=True)


def find_alpaca_output_file(output_dir):
    matches = sorted(Path(output_dir).glob("ALPACA_output*.csv"))
    if not matches:
        raise FileNotFoundError(
            f"Could not find ALPACA output CSV in {output_dir}. Expected a file named 'ALPACA_output*.csv'."
        )
    return matches[0]


def plot_heatmap_with_tree(tree, alpaca_output, cp_table, chr_table, driver_mutations=None, allele='A', max_cpn_cap=8):
    tumour_id = alpaca_output.tumour_id.iloc[0]
    mrca = tree[0][0]
    
    # drop diploid:
    alpaca_output = alpaca_output[alpaca_output.clone != 'diploid']
    # add coords:
    alpaca_output['chr'] = 'chr' + alpaca_output.segment.str.split('_', expand=True)[0]
    alpaca_output['Start'] = alpaca_output.segment.str.split('_', expand=True)[1].astype(int)
    alpaca_output['End'] = alpaca_output.segment.str.split('_', expand=True)[2].astype(int)
    # modify segment positions to absolute:
    alpaca_output = alpaca_output.merge(chr_table, on='chr', how='left')
    alpaca_output['abs_start'] = alpaca_output['Start'] + alpaca_output['shift']
    alpaca_output['abs_end'] = alpaca_output['End'] + alpaca_output['shift']
    alpaca_output = alpaca_output.sort_values(['abs_start'], ascending=False)
    alpaca_output.loc[alpaca_output['pred_CN_A'] > max_cpn_cap, 'pred_CN_A'] = max_cpn_cap
    alpaca_output.loc[alpaca_output['pred_CN_B'] > max_cpn_cap, 'pred_CN_B'] = max_cpn_cap
    max_cp_state = int(math.ceil(max(alpaca_output['pred_CN_A'].max(), alpaca_output['pred_CN_B'].max())))
    copy_number_palette = build_copy_number_palette(max_cp_state)
    
    number_of_clones = len(alpaca_output.clone.unique())
    max_levels = max([len(b) for b in tree])
    tree_with_levels = [dict(zip(b, range(0, len(b)))) for b in tree]
    tree_with_levels =   pd.concat([pd.DataFrame(tree_with_levels[x], index=[0]).transpose() for x, _ in enumerate(tree_with_levels)]).reset_index().drop_duplicates().rename(
        columns={'index': 'clone', 0: 'level'})
    
    # make empty y_loc df
    clone_y_location = dict(zip(alpaca_output.clone.unique(), range(0, number_of_clones)))
    clone_y_location = pd.DataFrame(clone_y_location, index=[0]).transpose().reset_index().rename(columns={'index': 'clone', 0: 'y_loc'})
    clone_y_location['y_loc'] = 100
    # find sections:
    # section is a part of a path that requires its own horizontal space on the final graph
    ori_tree = tree.copy()
    sections = [[tree[0][0]]]
    while len(flat_list(ori_tree)) > 1:
        for i, branch in enumerate(ori_tree):
            if (branch != []) and (remove_duplicates_preserve_order(flat_list(sections)) != remove_duplicates_preserve_order(flat_list(tree))):
                branching_clones = list(pd.Series(flat_list(ori_tree)).value_counts()[pd.Series(flat_list(ori_tree)).value_counts() > 1].index)
                if branching_clones == []:
                    branching_clones = [mrca]
                section_start = max([branch.index(x) for x in branching_clones if x in branch]) + 1
                section = branch[section_start:]
                ori_tree[ori_tree.index(branch)] = branch[:section_start]
                ori_tree = get_unique_lists(ori_tree)
                if section != []:
                    sections.append(section)
    # order sections according to proximity
    # start with on of the longest paths
    section_termini = [x[-1] for x in sections]
    initial_node = [s for s in sections if len(s) == max([len(x) for x in sections])][0][-1]
    # simplify tree graph:
    simple_tree = [[clone for clone in branch if clone in section_termini] for branch in tree]
    
    edges = get_tree_edges(simple_tree)
    
    G = nx.Graph()
    for edge in edges:
        G.add_edge(edge[0], edge[1], weight=0)
    
    nodes = [initial_node] + [x for x in list(G.nodes) if x is not initial_node]
    distance_to_nodes = dict(nx.all_pairs_shortest_path_length(G))
    processed_nodes = [[s for s in sections if s[-1] == initial_node][0][-1]]
    
    while len(nodes) > 0:
        node = processed_nodes[-1]
        neighbours = pd.DataFrame(distance_to_nodes[node], index=['val']).transpose()
        neighbours.drop(inplace=True, index=processed_nodes)
        try:
            neighbours.drop(inplace=True, index=node)
        except KeyError:
            pass
        neighbours = neighbours.reset_index().rename(columns={'index': 'clone'}).merge(tree_with_levels, how='left', on='clone')
        neighbours['level'] = neighbours['level'].astype(int)
        if len(neighbours) > 0:
            # if there are to neighbours equaly close, choose the one which has higher level (i.e. is deeper in the tree)
            # to do that, multiply proximity by negaive level:
            closest_candidates = neighbours[neighbours.val == min(neighbours.val)]
            closest_neighbours = closest_candidates[closest_candidates.level == max(closest_candidates.level)]
            closest_neighbour = closest_neighbours.clone.values[0]
            if closest_neighbour not in processed_nodes:
                processed_nodes.append(closest_neighbour)
        nodes = [n for n in nodes if n != node]
    
    sorted_sections_raw = []
    for node in processed_nodes:
        sorted_sections_raw.append([s for s in sections if s[-1] == node][0])
    
    # join single-clone sections below MRCA with their descendants:
    sorted_sections = []
    below_MRCA = True
    skip_element = False
    for ss in sorted_sections_raw:
        if skip_element:
            skip_element = False
            continue
        if mrca in ss:
            below_MRCA = False
        if below_MRCA:
            if len(ss) == 1:
                if sorted_sections_raw.index(ss) < len(sorted_sections_raw) - 1:
                    joined = ss + sorted_sections_raw[sorted_sections_raw.index(ss) + 1]
                    sorted_sections.append(joined)
                else:
                    sorted_sections[-1] = sorted_sections[-1] + ss
                skip_element = True
            else:
                sorted_sections.append(ss)
        else:
            sorted_sections.append(ss)
    
    # assign y location on the plot to sorted sections:
    available_y_locs = range(0, number_of_clones)
    below_MRCA = True
    for section in sorted_sections:
        if mrca in section:
            below_MRCA = False
        locs_for_this_section = list(available_y_locs[:len(section)])
        if below_MRCA:
            locs_for_this_section = list(reversed(locs_for_this_section))
        locs_for_this_section_dict = dict(zip(section, locs_for_this_section))
        available_y_locs = available_y_locs[len(section):]
        for n in locs_for_this_section_dict.keys():
            clone_y_location.loc[clone_y_location.clone == n, 'y_loc'] = locs_for_this_section_dict[n]
    
    tree_graph_df = pd.merge(tree_with_levels, clone_y_location)
    
    total_plot_height = max(1000, 75 * number_of_clones)
    clone_prop_title = 'Clone proportions in regions'

    s = [[{"type": "xy", "rowspan": number_of_clones}, {"type": "xy", "rowspan": number_of_clones}, {"type": "xy"}]]
    for c in range(number_of_clones - 1):
        s.append([None, None, {"type": "xy"}])
    
    fig = make_subplots(
        rows=number_of_clones, cols=3,
        column_widths=[0.1, 0.8, 0.1],
        specs=s, horizontal_spacing=0.02, vertical_spacing=0.01,
        subplot_titles=('', '', clone_prop_title))
    
    for clone_pos in tree_graph_df.y_loc:
        hline = go.Scatter(showlegend=False,
                           x=[-0.3, max_levels],
                           y=[clone_pos + 0.5, clone_pos + 0.5],
                           mode='lines',
                           line=dict(color='Green', dash='dot'))
        
        fig.append_trace(hline, row=1, col=1)
    
    # *** plot tree ***
    for branch in tree:
        branch_df = tree_graph_df[tree_graph_df.clone.isin(branch)]
        fig.append_trace(go.Scatter(
            showlegend=False,
            name='tree',
            x=branch_df['level'],
            y=branch_df['y_loc'],
            mode='lines+markers',
            marker=dict(
                symbol='circle',
                color='purple',
                size=10,
                line=dict(
                    color='purple',
                    width=2)),
            text=branch_df['clone']),
            row=1, col=1
        )
    fig.update_yaxes(
        showgrid=True,
        tickmode='array',
        tickvals=list(tree_graph_df.sort_values('y_loc').y_loc),
        ticktext=list(tree_graph_df.sort_values('y_loc').clone),
        range=[-0.5, number_of_clones - 0.5],
        showticklabels=True, zeroline=False, row=1, col=1)
    fig.update_xaxes(
        showgrid=True, zeroline=True, row=1, col=1
    )
    fig.update_yaxes(
        showgrid=True,
        tickmode='array',
        tickvals=list(tree_graph_df.sort_values('y_loc').y_loc),
        ticktext=list(tree_graph_df.sort_values('y_loc').clone),
        range=[-0.5, number_of_clones - 0.5],
        showticklabels=True, zeroline=False, row=1, col=2)
    fig.update_xaxes(
        showgrid=True, zeroline=True, row=1, col=2
    )
    
    # *** plot clones ***
    df = alpaca_output
    y_limit = df['pred_CN_A'].max()
    shapes = []
    chr_len = chr_table.copy()
    fig = plot_heat_map(
        alpaca_output.copy(),
        allele,
        fig,
        tree_graph_df,
        copy_number_palette,
        chr_table,
        driver_mutations,
    )
    for clone in df.clone.unique():
        i = int(tree_graph_df[tree_graph_df.clone == clone]['y_loc'].iloc[0])
        i = list(reversed(range(number_of_clones)))[i]
        
        fig.update_xaxes(
            showgrid=False, row=i + 1, col=2
        )
        # add proportions in regions:
        clone_cp = cp_table.loc[[clone]]
        
        showscale = i == len(df.clone.unique()) - 1
        clone_proportion_heatmap = go.Heatmap(
            z=clone_cp.values,
            x=clone_cp.columns,
            y=clone_cp.index,
            text=np.round(clone_cp.values,2),
            texttemplate="%{text}",
            textfont={"size":12},
            colorscale='Blues',
            showscale=False,
            colorbar=dict(
                tickfont=dict(size=12),
                orientation='h',
                x=0.89,
                y=-0.1,
                len=0.1,
                thickness=20),
            hoverinfo='z',
            zauto=False,
            zmin=0,
            zmax=1,
        )
        fig.add_trace(clone_proportion_heatmap, row=i + 1, col=3)
        
        
        # cleanup axes:
        
        fig.update_yaxes(showticklabels=False, row=i + 1, col=3)
        fig.update_yaxes(
            showticklabels=False, row=i + 1, col=2)
        sample_names = clone_cp.columns
        if i != len(df.clone.unique()) - 1:
            fig.update_xaxes(showticklabels=False, row=i + 1, col=1)
            
            fig.update_xaxes(showticklabels=False, row=i + 1, col=3)
        else:
            sample_names = clone_cp.columns
            # if sample names are in the long format, with tumour_id, split them:
            if tumour_id in sample_names[0]:
                sample_names = [x.split(f'{tumour_id}_')[1] for x in sample_names]
            fig.update_xaxes(
                tickmode='array',
                ticktext=sample_names
                ,showticklabels=True, row=i + 1, col=3)
        
        if i == 0:
            fig.update_xaxes(
                tickmode='array',
                tickvals=chr_len['cumsum'] - (chr_len['len'] / 2),
                ticktext=[str(x) for x in list(range(1, 23))]
                , showticklabels=True, row=i + 1, col=2)
    # subtitle font size:
    fig.update_annotations(font_size=12)
    fig.update_layout(
        title=f'{tumour_id}<br>Allele: {allele}',
        plot_bgcolor='rgba(255,255,255,0)',
        autosize=False,
        width=1600,
        height=total_plot_height,
        legend_tracegroupgap=10,
        legend=dict(
            orientation='h',
            yanchor="top",
            y=1.4,
            xanchor="left",
            x=0.2
        ))
    
    # build legend:
    
    max_palette_state = max(copy_number_palette.keys())
    legend_items = [
        (str(state), _format_rgb(copy_number_palette.get(state, copy_number_palette[max_palette_state])))
        for state in range(0, max_palette_state + 1)
    ]
    for label, colour in legend_items:
        fig.add_trace(go.Scatter(
            legendgroup='copy-number',
            x=[None],
            y=[None],
            mode='markers',
            name=label,
            marker=dict(
                color=colour,
                size=10,
                line=dict(color='black', width=1)
            ),
            showlegend=True
        ), row=1, col=1)
    
    return fig


def plot_cpn_per_clone(tree, alpaca_output, cp_table, chr_table, driver_mutations=None, max_cpn_cap=8):
    tumour_id = alpaca_output.tumour_id.iloc[0]
    mrca = tree[0][0]
    
    # drop diploid:
    alpaca_output = alpaca_output[alpaca_output.clone != 'diploid']
    # add coords:
    alpaca_output['chr'] = 'chr' + alpaca_output.segment.str.split('_', expand=True)[0]
    alpaca_output['Start'] = alpaca_output.segment.str.split('_', expand=True)[1].astype(int)
    alpaca_output['End'] = alpaca_output.segment.str.split('_', expand=True)[2].astype(int)
    # modify segment positions to absolute:
    alpaca_output = alpaca_output.merge(chr_table, on='chr', how='left')
    alpaca_output['abs_start'] = alpaca_output['Start'] + alpaca_output['shift']
    alpaca_output['abs_end'] = alpaca_output['End'] + alpaca_output['shift']
    alpaca_output = alpaca_output.sort_values(['abs_start'], ascending=False)
    alpaca_output.loc[alpaca_output['pred_CN_A'] > max_cpn_cap, 'pred_CN_A'] = max_cpn_cap
    alpaca_output.loc[alpaca_output['pred_CN_B'] > max_cpn_cap, 'pred_CN_B'] = max_cpn_cap
    
    number_of_clones = len(alpaca_output.clone.unique())
    max_levels = max([len(b) for b in tree])
    tree_with_levels = [dict(zip(b, range(0, len(b)))) for b in tree]
    tree_with_levels = pd.concat([pd.DataFrame(tree_with_levels[x], index=[0]).transpose() for x, _ in enumerate(tree_with_levels)]).reset_index().drop_duplicates().rename(
        columns={'index': 'clone', 0: 'level'})
    
    # make empty y_loc df
    clone_y_location = dict(zip(alpaca_output.clone.unique(), range(0, number_of_clones)))
    clone_y_location = pd.DataFrame(clone_y_location, index=[0]).transpose().reset_index().rename(columns={'index': 'clone', 0: 'y_loc'})
    clone_y_location['y_loc'] = 100
    # find sections:
    # section is a part of a path that requires its own horizontal space on the final graph
    ori_tree = tree.copy()
    sections = [[tree[0][0]]]
    while len(flat_list(ori_tree)) > 1:
        for i, branch in enumerate(ori_tree):
            if (branch != []) and (set(flat_list(sections)) != set(flat_list(tree))):
                branching_clones = list(pd.Series(flat_list(ori_tree)).value_counts()[pd.Series(flat_list(ori_tree)).value_counts() > 1].index)
                if branching_clones == []:
                    branching_clones = [mrca]
                section_start = max([branch.index(x) for x in branching_clones if x in branch]) + 1
                section = branch[section_start:]
                ori_tree[ori_tree.index(branch)] = branch[:section_start]
                ori_tree = get_unique_lists(ori_tree)
                if section != []:
                    sections.append(section)
    # order sections according to proximity
    # start with on of the longest paths
    section_termini = [x[-1] for x in sections]
    initial_node = [s for s in sections if len(s) == max([len(x) for x in sections])][0][-1]
    # simplify tree graph:
    simple_tree = [[clone for clone in branch if clone in section_termini] for branch in tree]
    
    edges = get_tree_edges(simple_tree)
    
    G = nx.Graph()
    for edge in edges:
        G.add_edge(edge[0], edge[1], weight=0)
    
    nodes = [initial_node] + [x for x in list(G.nodes) if x is not initial_node]
    distance_to_nodes = dict(nx.all_pairs_shortest_path_length(G))
    processed_nodes = [[s for s in sections if s[-1] == initial_node][0][-1]]
    
    while len(nodes) > 0:
        node = processed_nodes[-1]
        neighbours = pd.DataFrame(distance_to_nodes[node], index=['val']).transpose()
        neighbours.drop(inplace=True, index=processed_nodes)
        try:
            neighbours.drop(inplace=True, index=node)
        except KeyError:
            pass
        neighbours = neighbours.reset_index().rename(columns={'index': 'clone'}).merge(tree_with_levels, how='left', on='clone')
        neighbours['level'] = neighbours['level'].astype(int)
        if len(neighbours) > 0:
            # if there are to neighbours equaly close, choose the one which has higher level (i.e. is deeper in the tree)
            # to do that, multiply proximity by negaive level:
            closest_candidates = neighbours[neighbours.val == min(neighbours.val)]
            closest_neighbours = closest_candidates[closest_candidates.level == max(closest_candidates.level)]
            closest_neighbour = closest_neighbours.clone.values[0]
            if closest_neighbour not in processed_nodes:
                processed_nodes.append(closest_neighbour)
        nodes = [n for n in nodes if n != node]
    
    sorted_sections_raw = []
    for node in processed_nodes:
        sorted_sections_raw.append([s for s in sections if s[-1] == node][0])
    
    # join single-clone sections below MRCA with their descendants:
    sorted_sections = []
    below_MRCA = True
    skip_element = False
    for ss in sorted_sections_raw:
        if skip_element:
            skip_element = False
            continue
        if mrca in ss:
            below_MRCA = False
        if below_MRCA:
            if len(ss) == 1:
                joined = ss + sorted_sections_raw[sorted_sections_raw.index(ss) + 1]
                sorted_sections.append(joined)
                skip_element = True
            else:
                sorted_sections.append(ss)
        else:
            sorted_sections.append(ss)
    
    # assign y location on the plot to sorted sections:
    available_y_locs = range(0, number_of_clones)
    below_MRCA = True
    for section in sorted_sections:
        if mrca in section:
            below_MRCA = False
        locs_for_this_section = list(available_y_locs[:len(section)])
        if below_MRCA:
            locs_for_this_section = list(reversed(locs_for_this_section))
        locs_for_this_section_dict = dict(zip(section, locs_for_this_section))
        available_y_locs = available_y_locs[len(section):]
        for n in locs_for_this_section_dict.keys():
            clone_y_location.loc[clone_y_location.clone == n, 'y_loc'] = locs_for_this_section_dict[n]
    
    tree_graph_df = pd.merge(tree_with_levels, clone_y_location)
    
    clone_prop_title = 'Clone proportions in regions'

    s = [[{"type": "xy", "rowspan": number_of_clones}, {"type": "xy"}, {"type": "xy"}]]
    for c in range(number_of_clones - 1):
        s.append([None, {"type": "xy"}, {"type": "xy"}])
    
    fig = make_subplots(
        rows=number_of_clones, cols=3,
        column_widths=[0.1, 0.8, 0.1],
        specs=s, horizontal_spacing=0.02, vertical_spacing=0.01,
        subplot_titles=('', '', clone_prop_title))
    
    for clone_pos in tree_graph_df.y_loc:
        hline = go.Scatter(
            showlegend=False,
            x=[-0.3, max_levels],
            y=[clone_pos + 0.5, clone_pos + 0.5],
            mode='lines',
            line=dict(color='Green', dash='dot'))
        
        fig.append_trace(hline, row=1, col=1)
    
    # *** plot tree ***
    for branch in tree:
        branch_df = tree_graph_df[tree_graph_df.clone.isin(branch)]
        fig.append_trace(go.Scatter(
            showlegend=False,
            name='tree',
            x=branch_df['level'],
            y=branch_df['y_loc'],
            mode='lines+markers',
            marker=dict(
                symbol='circle',
                color='purple',
                size=10,
                line=dict(
                    color='purple',
                    width=2)),
            text=branch_df['clone']),
            row=1, col=1
        )
    fig.update_yaxes(
        showgrid=False,
        tickmode='array',
        tickvals=list(tree_graph_df.sort_values('y_loc').y_loc),
        ticktext=list(tree_graph_df.sort_values('y_loc').clone),
        range=[-0.3, number_of_clones - 0.7],
        showticklabels=True, zeroline=False, row=1, col=1)
    fig.update_xaxes(
        showgrid=False, zeroline=False, row=1, col=1
    )
    total_plot_height = max(1000, 75 * number_of_clones)
    # *** plot clones ***
    df = alpaca_output
    shapes = []
    chr_len = chr_table.copy()
    for clone in df.clone.unique():
        i = int(tree_graph_df[tree_graph_df.clone == clone]['y_loc'].iloc[0])
        i = list(reversed(range(number_of_clones)))[i]
        y_limit = df['pred_CN_A'].max()
        clone_df = df[df.clone == clone]
        clone_df = clone_df[['abs_start', 'abs_end', 'pred_CN_A', 'pred_CN_B']].drop_duplicates()
        clone_df = clone_df.sort_values('abs_start', ascending=True)
        clone_df['space'] = None
        x = flat_list([[v[1]['abs_start'], v[1]['abs_end'], v[1]['space']] for v in clone_df.iterrows()])
        ya = flat_list([[v[1]['pred_CN_A'], v[1]['pred_CN_A'], v[1]['space']] for v in clone_df.iterrows()])
        yb = flat_list([[v[1]['pred_CN_B'], v[1]['pred_CN_B'], v[1]['space']] for v in clone_df.iterrows()])
        clone_df_sameAB = clone_df[clone_df['pred_CN_A'] == clone_df['pred_CN_B']]
        if clone_df_sameAB.empty:
            yab = []
        else:
            yab = flat_list([[v[1]['pred_CN_B'], v[1]['pred_CN_B'], v[1]['space']] for v in clone_df_sameAB.iterrows()])
        clone_mutations = driver_mutations[driver_mutations.clone == clone] if driver_mutations is not None else pd.DataFrame([])
        
        fig.add_trace(go.Scatter(
            showlegend=False,
            x=x,
            y=ya,
            mode='lines',
            line=dict(
                color='rgb(255, 164, 0)',
                width=1)),
            row=i + 1, col=2
        )
        
        fig.add_trace(go.Scatter(
            showlegend=False,
            x=x,
            y=yb,
            mode='lines',
            line=dict(
                color='rgb(0, 128, 128)',
                width=1)),
            row=i + 1, col=2
        )
        
        fig.add_trace(go.Scatter(
            showlegend=False,
            x=x,
            y=yab,
            mode='lines',
            line=dict(
                color='rgb(187, 206, 186)',
                width=1)),
            row=i + 1, col=2
        )
        
        # add loss/gain areas:
        ancestor_clone = find_parent(clone, tree)
        if ancestor_clone == 'diploid':
            ancestor_df = clone_df.copy()
            ancestor_df.loc[:, 'pred_CN_A'] = 1
            ancestor_df.loc[:, 'pred_CN_B'] = 1
        else:
            ancestor_df = df[df.clone == ancestor_clone]
            ancestor_df = ancestor_df[['abs_start', 'abs_end', 'pred_CN_A', 'pred_CN_B']].drop_duplicates()
            ancestor_df = ancestor_df.sort_values('abs_start', ascending=True)
            ancestor_df['space'] = None
        gain_loss_df = pd.merge(clone_df, ancestor_df, on=['abs_start', 'abs_end', 'space'])
        gain_loss_df['gain'] = (gain_loss_df['pred_CN_A_x'] > gain_loss_df['pred_CN_A_y']) | (gain_loss_df['pred_CN_B_x'] > gain_loss_df['pred_CN_B_y'])
        gain_loss_df['loss'] = (gain_loss_df['pred_CN_A_x'] < gain_loss_df['pred_CN_A_y']) | (gain_loss_df['pred_CN_B_x'] < gain_loss_df['pred_CN_B_y'])
        
        gain_and_loss_df = gain_loss_df[(gain_loss_df['loss'] == True) & (gain_loss_df['gain'] == True)]
        gain_df = gain_loss_df[(gain_loss_df['loss'] == False) & (gain_loss_df['gain'] == True)]
        loss_df = gain_loss_df[(gain_loss_df['loss'] == True) & (gain_loss_df['gain'] == False)]
        D = {'gain_and_loss_df': gain_and_loss_df, 'gain_df': gain_df, 'loss_df': loss_df}
        
        xxyy = {}
        
        for d in D.items():
            if len(d[1] > 0):
                xxyy[f'x_{d[0]}'] = flat_list([[v[1]['abs_start'], v[1]['abs_start'], v[1]['abs_end'], v[1]['abs_end'], v[1]['abs_start'], v[1]['space']] for v in d[1].iterrows()])
                xxyy[f'y_{d[0]}'] = flat_list([[0, y_limit, y_limit, 0, 0, v[1]['space']] for v in d[1].iterrows()])
        
        try:
            fig.add_trace(go.Scatter(
                showlegend=False,
                x=xxyy['x_gain_df'],
                y=xxyy['y_gain_df'],
                mode='lines',
                fill="toself",
                line=dict(
                    color='rgb(255, 230, 229,0.5)',
                    width=1)),
                row=i + 1, col=2
            )
        except KeyError:
            pass
        try:
            fig.add_trace(go.Scatter(
                showlegend=False,
                x=xxyy['x_loss_df'],
                y=xxyy['y_loss_df'],
                mode='lines',
                fill="toself",
                line=dict(
                    color='rgb(230, 229, 255,0.5)',
                    width=1)),
                row=i + 1, col=2
            )
        except KeyError:
            pass
        try:
            fig.add_trace(go.Scatter(
                showlegend=False,
                x=xxyy['x_gain_and_loss_df'],
                y=xxyy['y_gain_and_loss_df'],
                mode='lines',
                fill="toself",
                line=dict(
                    color='rgb(230, 207, 232,0.5)',
                    width=1)),
                row=i + 1, col=2
            )
        except KeyError:
            pass
        if not clone_mutations.empty:
            fig.add_trace(go.Scatter(
                x=clone_mutations.abs_position,
                y=[0] * len(clone_mutations),
                hovertext=clone_mutations.get('gene'),
                mode='markers',
                showlegend=False,
                marker=dict(color='teal')
            ), row=i + 1, col=2)
        #add static annotations:
        '''
        for index, row in clone_mutations.iterrows():
            fig.add_annotation(
                x=row['abs_position'],
                y=[0]*len(clone_mutations),
                text=row['gene'],
                font=dict(size=20),
                showarrow=False,
                yshift=0,
                xref=f"x{i+1}",
                yref=f"y2"
            )
         '''   
        fig.update_xaxes(
            showgrid=False, row=i + 1, col=2
        )
        
        # add proportions in regions:
        clone_cp = cp_table.loc[[clone]]
        
        showscale = i == len(df.clone.unique()) - 1
        clone_proportion_heatmap = go.Heatmap(
            z=clone_cp.values,
            x=clone_cp.columns,
            y=clone_cp.index,
            text=np.round(clone_cp.values,2),
            texttemplate="%{text}",
            textfont={"size":12},
            colorscale='Blues',
            showscale=False,
            colorbar=dict(
                tickfont=dict(size=12),
                orientation='h',
                x=0.89,
                y=-0.1,
                len=0.1,
                thickness=20),
            hoverinfo='z',
            zauto=False,
            zmin=0,
            zmax=1,
        )
        fig.add_trace(clone_proportion_heatmap, row=i + 1, col=3)
        
        # cleanup axes:
        for chromosome_line in chr_len['cumsum']:
            fig.add_trace(go.Scatter(
                x=[chromosome_line, chromosome_line],
                y=[0, y_limit],
                mode='lines', line=dict(color='black', width=1, dash='dot'), showlegend=False), row=i + 1, col=2)
        fig.update_yaxes(showticklabels=False, row=i + 1, col=3)
        #fig.update_yaxes(showticklabels=False, row=i + 1, col=2)
        sample_names = clone_cp.columns
        if i != len(df.clone.unique()) - 1:
            fig.update_xaxes(showticklabels=False, row=i + 1, col=1)
            fig.update_xaxes(showticklabels=False, row=i + 1, col=2)
            fig.update_xaxes(showticklabels=False, row=i + 1, col=3)
        else:
            fig.update_xaxes(
                tickmode='array',
                tickvals=chr_len['cumsum'] - (chr_len['len'] / 2),
                ticktext=[str(x) for x in list(range(1, 23))]
                , showticklabels=True, row=i + 1, col=2)
            sample_names = clone_cp.columns
            # if sample names are in the long format, with tumour_id, split them:
            if tumour_id in sample_names[0]:
                sample_names = [x.split(f'{tumour_id}_')[1] for x in sample_names]
            fig.update_xaxes(
                tickmode='array',
                ticktext=sample_names
                ,showticklabels=True, row=i + 1, col=3)
        

            
    
    # subtitle font size:
    fig.update_annotations(font_size=12)
    fig.update_layout(
        title=f'{tumour_id}',
        plot_bgcolor='rgb(255,255,255)',
        autosize=False,
        width=1600,
        height=total_plot_height,
        legend_tracegroupgap=10,
        legend=dict(
            orientation='h',
            yanchor="top",
            y=1.4,
            xanchor="left",
            x=0.2
        ))
    
    # build legend:
    # gain loss areas:
    fig.add_trace(go.Scatter(
        legendgroup='1',
        x=[None],
        y=[None],
        mode="markers",
        name="gain relative to parent",
        marker=dict(size=10, color='rgb(255, 230, 229)', symbol='square'),
    ), row=1, col=1)
    fig.add_trace(go.Scatter(
        legendgroup='1',
        x=[None],
        y=[None],
        mode="markers",
        name="loss relative to parent",
        marker=dict(size=10, color='rgb(230, 229, 255)', symbol='square'),
    ), row=1, col=1)
    fig.add_trace(go.Scatter(
        legendgroup='1',
        x=[None],
        y=[None],
        mode="markers",
        name="gain and loss relative to parent",
        marker=dict(size=10, color='rgb(230, 207, 232)', symbol='square'),
    ), row=1, col=1)
    # cpn lines:
    fig.add_trace(go.Scatter(
        legendgroup='2',
        x=[None],
        y=[None],
        mode="markers+lines",
        name="A allele",
        marker=dict(size=10, color='rgb(255, 164, 0)', symbol='line-ew'),
    ), row=1, col=1)
    fig.add_trace(go.Scatter(
        legendgroup='2',
        x=[None],
        y=[None],
        mode="markers+lines",
        name="B allele",
        marker=dict(size=10, color='rgb(0, 128, 128)', symbol='line-ew'),
    ), row=1, col=1)
    fig.add_trace(go.Scatter(
        legendgroup='2',
        x=[None],
        y=[None],
        mode="markers+lines",
        name="A and B allele",
        marker=dict(size=10, color='rgb(187, 206, 186)', symbol='line-ew'),
    ), row=1, col=1)
    
    return fig


def plot_heat_map(patient_output, allele, fig, tree_graph_df, cpn_palette, chr_table, driver_mutations=None):
    clones = tree_graph_df.sort_values('y_loc', ascending=True).clone
    clone_number = len(clones)
    chromosome_table = chr_table.copy()
    tumour_id = patient_output.tumour_id.iloc[0]
    patient_output = clean_output(patient_output)
    patient_output['predicted_cpn'] = patient_output[f'pred_CN_{allele}']
    # patient_output['fractional_cpn'] = patient_output[f'mphase{allele}_allele']
    patient_output = patient_output.merge(tree_graph_df)
    patient_output = patient_output.sort_values('y_loc')
    max_palette_state = max(cpn_palette.keys())
    
    def getColour(cp_state):
        try:
            cp_key = int(round(cp_state))
        except (TypeError, ValueError):
            cp_key = 0
        cp_key = max(0, min(cp_key, max_palette_state))
        colour = cpn_palette.get(cp_key, cpn_palette[max_palette_state])
        return _format_rgb(colour)
    
    # fig = go.Figure(layout_xaxis_range=[0, 3.2 * 10 ** 9], layout_yaxis_range=[0, len(clones) + 1])
    for clone_index, clone_name in enumerate(clones):
        clone_df = patient_output[patient_output['clone'] == clone_name]
        for cp_state in clone_df.predicted_cpn.unique():
            clone_df_cp_state = clone_df[clone_df.predicted_cpn == cp_state]
            segments_predicted = [[tuple([row[1]['abs_start'], clone_index - 0.5]), tuple([row[1]['abs_end'], clone_index + 0.5])] for row in clone_df_cp_state.iterrows()]
            segments_predicted_unique = []
            for x in segments_predicted:
                if x not in segments_predicted_unique:
                    segments_predicted_unique.append(x)
            cpn_color = getColour(cp_state)
            for rectangle in segments_predicted_unique:
                segment = clone_df_cp_state[clone_df_cp_state['abs_start'] == rectangle[0][0]].segment.unique()[0]
                fig.add_trace(go.Scatter(
                    showlegend=False,
                    x=[rectangle[0][0], rectangle[0][0], rectangle[1][0], rectangle[1][0]],
                    y=[rectangle[0][1], rectangle[1][1], rectangle[1][1], rectangle[0][1]],
                    # y=[rectangle[0][1], rectangle[0][1] + 1, rectangle[0][1] + 1, rectangle[0][1]],
                    fill='toself',
                    mode='lines',
                    fillcolor=cpn_color,
                    line_color=cpn_color,
                    name=f'clone: {rectangle[0][1]}, seg: {segment}'), row=1, col=2)
        if driver_mutations is not None:
            clone_mutations = driver_mutations[driver_mutations.clone == clone_name]
            if not clone_mutations.empty:
                fig.add_trace(go.Scatter(
                    x=clone_mutations.abs_position,
                    y=[clone_index] * len(clone_mutations),
                    mode='markers',
                    showlegend=False,
                    marker=dict(color='teal')
                ), row=1, col=2)
        
    chromosomes = [[tuple([row[1]['cumsum'], 0]), tuple([row[1]['cumsum'], clone_number])] for row in chromosome_table.iterrows()]
    for chromosome_line in chromosomes:
        fig.add_trace(go.Scatter(
            x=[chromosome_line[0][0], chromosome_line[0][0]],
            y=[-0.5, clone_number + 0.5],
            mode='lines', line=dict(color='black', width=1, dash='dash'), showlegend=False), row=1, col=2)

    return fig


def main():
    parser = argparse.ArgumentParser(description="Plot ALPACA heatmap with tree")
    parser.add_argument('--output_directory', type=str, required=True, help='Output directory')
    parser.add_argument('--input_directory', type=str, required=True, help='Input directory')
    parser.add_argument('--chr-table', dest='chr_table', type=str, default=None,
                        help='Optional path to a chromosome length table (defaults to the packaged hg19 table).')
    parser.add_argument('--mutation-table', dest='mutation_table', type=str, default=None,
                        help='Optional path to a driver mutation table (CSV/TSV).')
    args = parser.parse_args()
    output_directory = Path(args.output_directory).expanduser().resolve()
    input_directory = Path(args.input_directory).expanduser().resolve()

    chr_table = load_chr_table(args.chr_table)

    tree_path = input_directory / 'tree_paths.json'
    if not tree_path.exists():
        raise FileNotFoundError(f"Missing tree_paths.json in {input_directory}")
    tree = read_tree_json(str(tree_path))

    alpaca_output_file = find_alpaca_output_file(output_directory)
    alpaca_output = pd.read_csv(alpaca_output_file)

    cp_table_path = input_directory / 'cp_table.csv'
    if not cp_table_path.exists():
        raise FileNotFoundError(f"Missing cp_table.csv in {input_directory}")
    cp_table = pd.read_csv(cp_table_path, index_col='clone')

    mutation_table = load_mutation_table(input_directory, args.mutation_table)
    tumour_id = alpaca_output.tumour_id.iloc[0]
    driver_mutations = prepare_driver_mutations(mutation_table, tumour_id, chr_table)

    for allele in ('A', 'B'):
        plot = plot_heatmap_with_tree(
            tree,
            alpaca_output,
            cp_table,
            chr_table,
            driver_mutations,
            allele=allele,
        )
        plot.write_image(str(output_directory / f"{tumour_id}_{allele}_heatmap.pdf"))

    cn_changes = plot_cpn_per_clone(
        tree,
        alpaca_output,
        cp_table,
        chr_table,
        driver_mutations,
        max_cpn_cap=8,
    )
    cn_changes.write_image(str(output_directory / f"{tumour_id}_cn_changes_per_clone.pdf"))

if __name__ == "__main__":
    main()

