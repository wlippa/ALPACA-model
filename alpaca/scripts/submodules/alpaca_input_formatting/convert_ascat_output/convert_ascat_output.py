#!/usr/bin/env python3
"""

This script reads TSVs exported by `unpack_RDS_ascat_helper.R` (combined files
containing all samples) and computes bootstrapped CIs for `nAraw`/`nBraw`,
recentering intervals on ASCAT's original fractional copy numbers.

Next, it writes output in ALPACA input format.

Usage:
  python extract_and_bootstrap.py \
    --tumour_id LTX0000_Tumour1 \
    --output_dir '.' \
    --segments out_segments_all.tsv \
    --snps out_snp_data_all.tsv \
    --purity_ploidy out_purity_ploidy_all.tsv \
    --n_boot 1000 \
    --gamma 1


Output:
  Writes ALPACA inputs files to `<output_dir>`

Requirements:
  pip install pandas numpy
  Obtain gamma paramter from your ASCAT run (typically 1 or 0.55) and set `--gamma` accordingly.
"""

import argparse
import sys
import numpy as np
import pandas as pd


def compute_nAB_raw(logr, baf, chr_arr, rho, psi, gamma, haploid_chrs, null_chrs):
    logr = np.asarray(logr, dtype=float)
    baf = np.asarray(baf, dtype=float)
    chr_arr = np.asarray(chr_arr, dtype=object)

    diploid = ~np.isin(chr_arr, haploid_chrs)
    is_null = np.isin(chr_arr, null_chrs)
    scaling = (1 - rho) * 2 + rho * psi
    powered = np.power(2.0, logr / float(gamma))

    nAraw = np.where(
        diploid,
        (rho - 1.0 - (baf - 1.0) * powered * scaling) / rho,
        np.where(is_null, 0.0, (rho - 1.0 + scaling * powered) / rho),
    )
    nBraw = np.where(
        diploid,
        (rho - 1.0 + baf * powered * scaling) / rho,
        0.0,
    )

    total = nAraw + nBraw
    neg_both = total < 0
    neg_A = (~neg_both) & (nAraw < 0)
    neg_B = (~neg_both) & (nBraw < 0)

    nAraw[neg_both] = 0.0
    nBraw[neg_both] = 0.0
    nBraw[neg_A] = total[neg_A]
    nAraw[neg_A] = 0.0
    nAraw[neg_B] = total[neg_B]
    nBraw[neg_B] = 0.0

    return nAraw, nBraw


def main():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--tumour_id", required=True, help="Tumour ID for the sample"
    )
    p.add_argument(
        "--segments", required=True, help="TSV with combined segments (all samples)"
    )
    p.add_argument(
        "--snps", required=True, help="TSV with combined per-SNP data (all samples)"
    )
    p.add_argument(
        "--purity_ploidy",
        required=True,
        help="TSV with combined purity/ploidy/meta (all samples)",
    )
    p.add_argument("--gamma", type=float, default=1)
    p.add_argument("--n_boot", type=int, default=1000)
    p.add_argument("--output_dir", default=".", help="Directory to write outputs")
    p.add_argument(
        "--ci_level",
        type=float,
        default=0.95,
        help="Confidence interval level between 0 and 1 (e.g., 0.95)",
    )
    p.add_argument(
        "--sample",
        default=None,
        help="Optional: process only this sample name",
    )
    args = p.parse_args()

    # Validate CI level and compute percentiles (as percentages for np.nanpercentile)
    if args.ci_level <= 0.0 or args.ci_level >= 1.0:
        print("Invalid --ci_level; must be between 0 and 1 (exclusive).")
        sys.exit(1)
    ci_level = float(args.ci_level)
    lower_pct = (1.0 - ci_level) / 2.0 * 100.0 
    upper_pct = (1.0 + ci_level) / 2.0 * 100.0

    seg_df = pd.read_csv(args.segments, sep="\t")
    snp_df = pd.read_csv(args.snps, sep="\t")
    pur_df = pd.read_csv(args.purity_ploidy, sep="\t")

    # Build sample list to process using only the `sample` column
    samples_info = []
    if args.sample is not None:
        sarg = str(args.sample)
        # prefer purity/ploidy file if present
        if "sample" in pur_df.columns:
            prow = pur_df[pur_df["sample"].astype(str) == sarg]
            if prow.shape[0] == 0:
                # fallback to segments or snps
                if "sample" in seg_df.columns:
                    uniq = seg_df[seg_df["sample"].astype(str) == sarg][["sample"]].drop_duplicates()
                    if uniq.shape[0] == 0:
                        print("Sample name not found in purity or segments TSVs:", sarg)
                        sys.exit(1)
                    samples_info = [ {"sample": r["sample"]} for i, r in uniq.reset_index(drop=True).iterrows() ]
                elif "sample" in snp_df.columns:
                    uniq = snp_df[snp_df["sample"].astype(str) == sarg][["sample"]].drop_duplicates()
                    if uniq.shape[0] == 0:
                        print("Sample name not found in purity or snps TSVs:", sarg)
                        sys.exit(1)
                    samples_info = [ {"sample": r["sample"]} for i, r in uniq.reset_index(drop=True).iterrows() ]
                else:
                    print("Sample name not found in purity TSV:", sarg)
                    sys.exit(1)
            else:
                samples_info = prow[["sample"]].drop_duplicates().to_dict("records")
        elif "sample" in seg_df.columns:
            uniq = seg_df[seg_df["sample"].astype(str) == sarg][["sample"]].drop_duplicates()
            if uniq.shape[0] == 0:
                print("Sample name not found in segments TSV:", sarg)
                sys.exit(1)
            samples_info = [ {"sample": r["sample"]} for i, r in uniq.reset_index(drop=True).iterrows() ]
        elif "sample" in snp_df.columns:
            uniq = snp_df[snp_df["sample"].astype(str) == sarg][["sample"]].drop_duplicates()
            if uniq.shape[0] == 0:
                print("Sample name not found in snps TSV:", sarg)
                sys.exit(1)
            samples_info = [ {"sample": r["sample"]} for i, r in uniq.reset_index(drop=True).iterrows() ]
        else:
            print("Cannot resolve sample name:", sarg)
            sys.exit(1)
    else:
        # default: gather unique samples from purity, then segments, then snps
        if "sample" in pur_df.columns:
            uniq = pur_df["sample"].astype(str).unique()
            samples_info = [{"sample": s} for s in uniq]
        elif "sample" in seg_df.columns:
            uniq = seg_df["sample"].astype(str).unique()
            samples_info = [{"sample": s} for s in uniq]
        elif "sample" in snp_df.columns:
            uniq = snp_df["sample"].astype(str).unique()
            samples_info = [{"sample": s} for s in uniq]
        else:
            print("Could not determine samples from purity/segments/snp TSVs")
            sys.exit(1)

    ci_all = []

    for info in samples_info:
        sample_name = str(info.get("sample"))
        print("Processing sample:", sample_name)

        if "sample" in seg_df.columns:
            seg_df_sample = seg_df[seg_df["sample"].astype(str) == sample_name].reset_index(
                drop=True
            )
        else:
            print("Segments TSV missing 'sample' column; cannot select sample:", sample_name)
            sys.exit(1)

        if "sample" in snp_df.columns:
            snp_df_sample = snp_df[snp_df["sample"].astype(str) == sample_name].reset_index(
                drop=True
            )
        else:
            print("SNP TSV missing 'sample' column; cannot select sample:", sample_name)
            sys.exit(1)

        if seg_df_sample.shape[0] == 0:
            print("No segments for sample", sample_name, "- skipping")
            continue

        # Ensure snp_df_sample has segment_id
        if "segment_id" not in snp_df_sample.columns:
            snp_df_sample["segment_id"] = np.nan
            for i, r in seg_df_sample.iterrows():
                mask = (
                    (snp_df_sample["chr"] == r["chr"])
                    & (snp_df_sample["pos"] >= r["startpos"])
                    & (snp_df_sample["pos"] <= r["endpos"])
                )
                snp_df_sample.loc[mask, "segment_id"] = i + 1

        # purity / psi / ploidy / gender (use `sample` column only)
        if "sample" in pur_df.columns:
            prow = pur_df[pur_df["sample"].astype(str) == sample_name]
            if prow.shape[0] == 0:
                prow = pur_df.iloc[0:1]
        else:
            prow = pur_df.iloc[0:1]

        rho = float(prow["purity"].iloc[0])
        try:
            psi_val = float(prow["psi"].iloc[0])
        except Exception:
            psi_val = float(prow.get("psi", pd.Series([1.0])).iloc[0])
        ploidy = float(prow["ploidy"].iloc[0]) if "ploidy" in prow.columns else None
        gender = str(prow["gender"].iloc[0]) if "gender" in prow.columns else "XX"
        sexchromosomes = (
            str(prow["sexchromosomes"].iloc[0]).split(",")
            if "sexchromosomes" in prow.columns
            else ["X", "Y"]
        )

        # best-effort haploid/null chromosomes
        haploid_chrs = list(set([gender[0], gender[1]])) if len(gender) >= 2 else []
        if len(gender) >= 2 and gender[0] == gender[1]:
            haploid_chrs = [x for x in haploid_chrs if x != gender[0]]
        null_chrs = [
            c
            for c in sexchromosomes
            if c not in (gender[0:2] if len(gender) >= 2 else [])
        ]

        # bootstrap
        n_seg = len(seg_df_sample)
        boot_nA = np.full((args.n_boot, n_seg), np.nan)
        boot_nB = np.full((args.n_boot, n_seg), np.nan)
        rng = np.random.default_rng(42)

        for b in range(args.n_boot):
            if (b + 1) % 100 == 0:
                print("iteration", b + 1, "/", args.n_boot, "for sample", sample_name)
            for s in range(n_seg):
                seg_snps = snp_df_sample[snp_df_sample["segment_id"] == (s + 1)]
                if seg_snps.shape[0] == 0:
                    continue
                idx = rng.choice(
                    seg_snps.index.values, size=seg_snps.shape[0], replace=True
                )
                boot_logr = seg_snps.loc[idx, "LogR"].values
                boot_baf = seg_snps.loc[idx, "BAF"].values
                boot_chr = seg_snps.loc[idx, "chr"].values
                valid = (~np.isnan(boot_logr)) & (~np.isnan(boot_baf))
                if valid.sum() == 0:
                    continue
                mean_logr = float(np.mean(boot_logr[valid]))
                mean_baf = float(np.mean(boot_baf[valid]))
                seg_chr_val = boot_chr[valid][0]
                nA_val, nB_val = compute_nAB_raw(
                    [mean_logr],
                    [mean_baf],
                    [seg_chr_val],
                    rho,
                    psi_val,
                    args.gamma,
                    haploid_chrs,
                    null_chrs,
                )
                boot_nA[b, s] = nA_val[0]
                boot_nB[b, s] = nB_val[0]

        # center and compute CI
        if "nAraw" in seg_df_sample.columns:
            orig_nAraw = seg_df_sample["nAraw"].astype(float).values
        else:
            orig_nAraw = np.zeros(n_seg, dtype=float)
        if "nBraw" in seg_df_sample.columns:
            orig_nBraw = seg_df_sample["nBraw"].astype(float).values
        else:
            orig_nBraw = np.zeros(n_seg, dtype=float)

        boot_nA_mean = np.nanmean(boot_nA, axis=0)
        boot_nB_mean = np.nanmean(boot_nB, axis=0)
        boot_nA_lo = np.nanpercentile(boot_nA, lower_pct, axis=0)
        boot_nA_hi = np.nanpercentile(boot_nA, upper_pct, axis=0)
        boot_nB_lo = np.nanpercentile(boot_nB, lower_pct, axis=0)
        boot_nB_hi = np.nanpercentile(boot_nB, upper_pct, axis=0)

        shifted_nAraw_ci_lo = np.maximum(0.0, orig_nAraw + (boot_nA_lo - boot_nA_mean))
        shifted_nAraw_ci_hi = np.maximum(0.01, orig_nAraw + (boot_nA_hi - boot_nA_mean))
        shifted_nBraw_ci_lo = np.maximum(0.0, orig_nBraw + (boot_nB_lo - boot_nB_mean))
        shifted_nBraw_ci_hi = np.maximum(0.01, orig_nBraw + (boot_nB_hi - boot_nB_mean))

        ci_df = pd.DataFrame(
            {
                "sample": sample_name,
                "segment_id": np.arange(1, n_seg + 1),
                "chr": seg_df_sample["chr"].values,
                "startpos": seg_df_sample["startpos"].values,
                "endpos": seg_df_sample["endpos"].values,
                "nAraw": orig_nAraw,
                "nBraw": orig_nBraw,
                "nAraw_mean_boot": boot_nA_mean,
                "nBraw_mean_boot": boot_nB_mean,
                "nAraw_ci_lo": shifted_nAraw_ci_lo,
                "nAraw_ci_hi": shifted_nAraw_ci_hi,
                "nBraw_ci_lo": shifted_nBraw_ci_lo,
                "nBraw_ci_hi": shifted_nBraw_ci_hi,
                "n_snps": [
                    int(np.sum(snp_df_sample["segment_id"] == (s + 1)))
                    for s in range(n_seg)
                ],
            }
        )

        ci_all.append(ci_df)

    if len(ci_all) == 0:
        print("No CI results generated (no samples processed).")
        sys.exit(1)

    result_df = pd.concat(ci_all, ignore_index=True)
    
    # set expected column names:
    result_df.rename(columns={
        "nAraw_ci_lo": "lower_CI_A",
        "nAraw_ci_hi": "upper_CI_A",
        "nBraw_ci_lo": "lower_CI_B",
        "nBraw_ci_hi": "upper_CI_B",
        "nAraw": "cpnA",
        "nBraw": "cpnB",
    }, inplace=True)
    
    result_df['tumour_id'] = args.tumour_id
    result_df['chr'] = result_df['chr'].astype(str).replace("chr", "", regex=False).replace("Chr", "", regex=False).replace("CHR", "", regex=False).replace("X", "23").replace("Y", "24")
    result_df['segment'] = result_df['chr'] + "_" + result_df['startpos'].astype(str) + "_" + result_df['endpos'].astype(str)
    result_df['ci_value'] = ci_level
    ci_table = result_df.copy()[["segment", "sample", "lower_CI_A", "upper_CI_A", "lower_CI_B", "upper_CI_B", "tumour_id", "ci_value"]]
    input_table = result_df.copy()[["tumour_id", "sample", "segment", "cpnA", "cpnB"]]
    ci_table_file = f"{args.output_dir}/ci_table.csv"
    input_table_file = f"{args.output_dir}/ALPACA_input_table.csv"
    ci_table.to_csv(ci_table_file, index=False)
    input_table.to_csv(input_table_file, index=False)
    print("Wrote", ci_table_file, "and", input_table_file)


if __name__ == "__main__":
    main()
