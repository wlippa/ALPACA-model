import json

from alpaca.utils import read_tree_json, newick_to_paths


def test_newick_to_paths_parses_branch_lengths():
    newick = "((clone2:0.1,clone3:0.2)clone1:0.3,clone4:0.5)root;"
    expected = [
        ["root", "clone1", "clone2"],
        ["root", "clone1", "clone3"],
        ["root", "clone4"],
    ]
    assert newick_to_paths(newick) == expected


def test_read_tree_json_prefers_json_over_newick(tmp_path):
    json_tree = [["root", "a"], ["root", "b"]]
    json_path = tmp_path / "tree_paths.json"
    json_path.write_text(json.dumps(json_tree))
    nwk_path = tmp_path / "tree_paths.nwk"
    nwk_path.write_text("(a,b)root;")

    result = read_tree_json(str(json_path))
    assert result == json_tree


def test_read_tree_json_falls_back_to_newick(tmp_path):
    nwk_path = tmp_path / "tree_paths.nwk"
    nwk_path.write_text("(clone2,clone3)clone1;")

    result = read_tree_json(str(tmp_path / "tree_paths.json"))
    assert result == [["clone1", "clone2"], ["clone1", "clone3"]]


def test_newick_handles_trailing_wrapped_parents():
    newick = "((clone14)clone1,((clone7,clone8,(((clone20)clone4)clone16)clone17)clone6,((clone15)clone11)clone9)clone2)clone3;"
    expected = [
        ['clone3', 'clone1', 'clone14'], 
        ['clone3', 'clone2', 'clone6', 'clone7'], 
        ['clone3', 'clone2', 'clone6', 'clone8'], 
        ['clone3', 'clone2', 'clone6', 'clone17', 'clone16', 'clone4', 'clone20'], 
        ['clone3', 'clone2', 'clone9', 'clone11', 'clone15']
        ]
    assert newick_to_paths(newick) == expected
