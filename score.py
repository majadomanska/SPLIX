#!/usr/bin/env python3

import argparse
import sys
import csv
import subprocess
import re

import pandas as pd


csv.field_size_limit(sys.maxsize)


def norm_seq(s):
    """Normalize sequence."""
    if s is None:
        return ""
    if pd.isna(s):
        return ""
    return str(s).strip().upper()


def to_rna(seq):
    """Convert T to U for pairing logic."""
    return norm_seq(seq).replace("T", "U")


def cag_ctg(intron):
    """Find CAG motif near the 5' end 
    and CTG motif near the 3' end of the intron.
    """
    intron = norm_seq(intron)

    if not intron:
        return 0

    has_cag = "CAG" in intron[:10]
    has_ctg = "CTG" in intron[-10:] if len(intron) >= 3 else False
    return int(has_cag and has_ctg)


PURINES = {"A", "G"}
PYRIMIDINES = {"C", "T", "U"}


def is_purine(base):
    """Return 1 if the base is a purine, otherwise return 0."""
    base = norm_seq(base)

    if not base:
        return 0
    
    return int(base in PURINES)


def is_pyrimidine(base):
    """Return 1 if the base is a pyrimidine, otherwise return 0."""
    base = norm_seq(base)

    if not base:
        return 0
    
    return int(base in PYRIMIDINES)


def last_base(seq):
    """ Return last nucleotide of the sequence."""
    seq = norm_seq(seq)

    if not seq:
        return ""
    
    return seq[-1]


def first_base(seq):
    """ Return first nucleotide of the sequence."""
    seq = norm_seq(seq)

    if not seq:
        return ""
    
    return seq[0]


# Allowed base pairs:
PAIR_OK = [
    {"A", "U"},
    {"C", "G"},
    {"G", "U"},
]


def base_in_pos(seq, pos):
    """Return base in position (from start)."""
    seq = norm_seq(seq)

    if pos <= 0 or len(seq) < pos:
        return ""
    
    return seq[pos - 1]


def base_from_end_pos(seq, k_from_end):
    """Return base in k-position from end."""
    seq = norm_seq(seq)

    if k_from_end <= 0 or len(seq) < k_from_end:
        return ""
    
    return seq[-k_from_end]

def pairs(b1, b2):
    """Return 1 if two bases can be paired."""
    b1 = to_rna(b1)
    b2 = to_rna(b2)

    if not b1 or not b2:
        return 0
    
    return int({b1, b2} in PAIR_OK)


def is_cg_pair(b1, b2):
    """Return 1 if bases are C-G pair."""
    b1 = norm_seq(b1)
    b2 = norm_seq(b2)
    return int(b1 == "C" and b2 == "G")


def is_at_pair(b1, b2):
    """Return 1 if bases are A-T pair."""
    b1 = norm_seq(b1).replace("U", "T")
    b2 = norm_seq(b2).replace("U", "T")
    return int(b1 == "A" and b2 == "T")


def is_gc_pair(b1, b2):
    """Return 1 if bases are G-C pair."""
    b1 = norm_seq(b1)
    b2 = norm_seq(b2)
    return int(b1 == "G" and b2 == "C")


def compute_pair_cols(intron_seq):
    """Compute structural pairing features."""
    b4 = base_in_pos(intron_seq, 4)
    bm6 = base_from_end_pos(intron_seq, 6)

    b5 = base_in_pos(intron_seq, 5)
    bm7 = base_from_end_pos(intron_seq, 7)

    b6 = base_in_pos(intron_seq, 6)
    bm8 = base_from_end_pos(intron_seq, 8)

    return {
        "4_-6": pairs(b4, bm6),
        "5_-7": pairs(b5, bm7),
        "6_-8": pairs(b6, bm8),
        "CG_4_-6": is_cg_pair(b4, bm6),
        "AT_5_-7": is_at_pair(b5, bm7),
        "GC_6_-8": is_gc_pair(b6, bm8),
    }


def compute_boundary_cols(prev_exon_last5, intron_full, next_exon_first5):
    """Compute nucleotide class features at exon-intron boundaries."""
    return {
        "prev_exon_Y": is_pyrimidine(last_base(prev_exon_last5)),
        "5'_intron_R": is_purine(first_base(intron_full)),
        "3'_intron_Y": is_pyrimidine(last_base(intron_full)),
        "next_exon_R": is_purine(first_base(next_exon_first5)),
    }


def compute_conventional_cols(intron_seq):
    """Compute features typical for conventional introns."""
    intron_seq = norm_seq(intron_seq)

    b1 = base_in_pos(intron_seq, 1)
    b2 = base_in_pos(intron_seq, 2)
    bm1 = base_from_end_pos(intron_seq, 1)
    bm2 = base_from_end_pos(intron_seq, 2)
    bm3 = base_from_end_pos(intron_seq, 3)

    pyrimidine_positions = [
        base_from_end_pos(intron_seq, k)
        for k in range(3, 13)
    ]

    n_pyrimidines = sum(
        1 for base in pyrimidine_positions
        if base in PYRIMIDINES
    )

    return {
        "G_1": int(b1 == "G"),
        "T_2": int(b2 == "T"),
        "C_2": 0.5 if b2 == "C" else 0,
        "G_-1": int(bm1 == "G"),
        "A_-2": int(bm2 == "A"),
        "C_-3": int(bm3 == "C"),
        "pyrimidines_-3_to_-12": round(n_pyrimidines / 10, 2),
    }


def make_secondary_seq(intron_first20, intron_last20, spacer_len=20): #spacer =20 n DOL
    """Build artificial sequence for RNAfold."""
    left_intron = to_rna(intron_first20)
    right_intron = to_rna(intron_last20)
    spacer = "N" * spacer_len

    return left_intron + spacer + right_intron


def run_rnafold(seq):
    """Run RNAfold and return Vienna structure."""
    if not seq:
        return ""

    try:
        result = subprocess.run(
            ["RNAfold", "--noPS"],
            input=seq + "\n",
            text=True,
            capture_output=True,
            check=True,
        )

        lines = result.stdout.strip().splitlines()

        if len(lines) < 2:
            return ""

        second_line = lines[1].strip()
        m = re.match(r"^([().]+)\s+\(", second_line)

        if not m:
            return ""

        return m.group(1)

    except Exception:
        return ""


def parse_dot_bracket(structure):
    """Return pairing map for paired position."""
    stack = []
    pairs_map = {}

    for i, ch in enumerate(structure):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            if stack:
                j = stack.pop()
                pairs_map[i] = j
                pairs_map[j] = i

    return pairs_map


def get_cross_pairs(structure, left_len=20, spacer_len=20, right_len=20):  #TODO
    """
    Return base pairs connecting the left and right intron regions.
    Output pairs are tuples: (left_pos, right_pos), 0-based.
    """
    if not structure:
        return []

    pairs_map = parse_dot_bracket(structure)

    left_start = 0
    left_end = left_len - 1
    right_start = left_len + spacer_len
    right_end = right_start + right_len - 1

    cross_pairs = []
    seen = set()

    for i, j in pairs_map.items():
        pair = tuple(sorted((i, j)))

        if pair in seen:
            continue

        seen.add(pair)

        if left_start <= i <= left_end and right_start <= j <= right_end:
            cross_pairs.append((i, j))
        elif left_start <= j <= left_end and right_start <= i <= right_end:
            cross_pairs.append((j, i))

    cross_pairs.sort(key=lambda x: x[0])

    return cross_pairs


def longest_consecutive_pairs(cross_pairs):
    """
    Find longest run of consecutive pairs:
    left +1 and right -1.
    """
    if not cross_pairs:
        return []

    best = [cross_pairs[0]]
    current = [cross_pairs[0]]

    for k in range(1, len(cross_pairs)):
        prev_left, prev_right = cross_pairs[k - 1]
        curr_left, curr_right = cross_pairs[k]

        if curr_left == prev_left + 1 and curr_right == prev_right - 1:
            current.append(cross_pairs[k])
        else:
            if len(current) > len(best):
                best = current
            current = [cross_pairs[k]]

    if len(current) > len(best):
        best = current

    return best


def pair_score(b1, b2):
    """Score RNA base pair: CG/GC = 3, AU/UA = 2, GU/UG = 1."""
    b1 = to_rna(b1)
    b2 = to_rna(b2)
    pair = {b1, b2}

    if pair == {"C", "G"}:
        return 3
    if pair == {"A", "U"}:
        return 2
    if pair == {"G", "U"}:
        return 1
    
    return 0


def score_cross_pairs(seq, cross_pairs):
    """Score selected cross-pairs."""
    if not seq or not cross_pairs:
        return 0

    total = 0

    for left_pos, right_pos in cross_pairs:
        total += pair_score(seq[left_pos], seq[right_pos])

    return total


def normalize_score(value, max_value):
    """Normalize score to range 0-1."""
    try:
        value = float(value)
    except Exception:
        return 0.0

    if max_value <= 0:
        return 0.0

    normalized = value / max_value

    if normalized < 0:
        return 0.0
    if normalized > 1:
        return 1.0

    return round(normalized, 2)


def compute_secondary_structure_cols(
    intron_first20,
    intron_last20,
    max_pairs_score_max=48,
    total_pairs_score_max=50,
):
    """Compute RNAfold-based secondary structure columns."""
    intron_first20 = norm_seq(intron_first20)
    intron_last20 = norm_seq(intron_last20)

    if len(intron_first20) < 20 or len(intron_last20) < 20:
        return {
            "max_pairs": 0,
            "max_pairs_score": 0,
            "total_pairs": 0,
            "total_pairs_score": 0,
            "max_pairs_normalized": 0.0,
            "max_pairs_score_normalized": 0.0,
            "total_pairs_normalized": 0.0,
            "total_pairs_score_normalized": 0.0,
            "secondary_structure": "",
        }

    seq = make_secondary_seq(intron_first20, intron_last20, spacer_len=20)
    structure = run_rnafold(seq)
    cross_pairs = get_cross_pairs(structure, left_len=20, spacer_len=20, right_len=20)
    max_pair_run = longest_consecutive_pairs(cross_pairs)

    max_pairs = len(max_pair_run)
    max_pairs_score = score_cross_pairs(seq, max_pair_run)
    total_pairs = len(cross_pairs)
    total_pairs_score = score_cross_pairs(seq, cross_pairs)

    return {
        "max_pairs": max_pairs,
        "max_pairs_score": max_pairs_score,
        "total_pairs": total_pairs,
        "total_pairs_score": total_pairs_score,
        "max_pairs_normalized": normalize_score(max_pairs, 20),
        "max_pairs_score_normalized": normalize_score(max_pairs_score, max_pairs_score_max),
        "total_pairs_normalized": normalize_score(total_pairs, 20),
        "total_pairs_score_normalized": normalize_score(total_pairs_score, total_pairs_score_max),
        "secondary_structure": structure,
    }


def run_scoring(
    inp_csv, 
    out_csv, 
    include_full=False, 
    filter_c_lub_n=None, 
    filter_konstyt=None, 
    short_out_csv=None, 
    short_intron_length=30,
    max_pairs_score_max=48,
    total_pairs_score_max=50,    
    ):
    """Run intron scoring and save output to CSV."""
    df = pd.read_csv(inp_csv, sep=None, engine="python", dtype=str)
    df.columns = [str(c).strip() for c in df.columns]

    for col in [
        "prev_exon_last5",
        "next_exon_first5",
        "intron_first20",
        "intron_last20",
        "intron_full",
    ]:
        df[col] = df[col].apply(norm_seq)

    df["length"] = df["intron_full"].apply(len)

    df["conv_or_nonconv"] = df["conv_or_nonconv"].fillna("").astype(str).str.strip()
    df["const_or_alt"] = df["const_or_alt"].fillna("").astype(str).str.strip()

    if filter_c_lub_n is not None:
        df = df[df["conv_or_nonconv"] == str(filter_c_lub_n)]

    if filter_konstyt is not None:
        df = df[df["const_or_alt"] == str(filter_konstyt)]

    df = df.reset_index(drop=True)

    base_cols = [
        "scaffold",
        "id",
        "name",
        "prev_exon_Y",
        "5'_intron_R",
        "3'_intron_Y",
        "next_exon_R",
        "CAG_CTG",
        "4_-6",
        "5_-7",
        "6_-8",
        "CG_4_-6",
        "AT_5_-7",
        "GC_6_-8",
        "SUM_sequence",
        "max_pairs",
        "max_pairs_score",
        "total_pairs",
        "total_pairs_score",
        "max_pairs_normalized",
        "max_pairs_score_normalized",
        "total_pairs_normalized",
        "total_pairs_score_normalized",
        "SUM_structure",
        "G_1",
        "T_2",
        "C_2",
        "G_-1",
        "A_-2",
        "C_-3",
        "pyrimidines_-3_to_-12",
        "SUM_conventional",
        "conv_or_nonconv",
        "const_or_alt",
        "intron_full",
    ]

    full_cols = [
        "length",
        "prev_exon_last5",
        "next_exon_first5",
        "intron_first20",
        "intron_last20",
        "secondary_structure",
    ]

    if include_full:
        cols = base_cols[:3] + full_cols + base_cols[3:]
    else:
        cols = base_cols

    if df.empty:
        out = pd.DataFrame(columns=cols)
        out.to_csv(out_csv, index=False)

        return {
            "saved_rows": 0,
            "short_rows": 0,
            "out_csv": str(out_csv),
            "short_out_csv": None,
        }

    out = pd.DataFrame()
    out["scaffold"] = df["scaffold"].astype(str)
    out["id"] = df["id"].astype(str)
    out["name"] = df["name"].astype(str)

    boundary_df = df.apply(
        lambda row: compute_boundary_cols(
            row["prev_exon_last5"],
            row["intron_full"],
            row["next_exon_first5"],
        ),
        axis=1,
    ).apply(pd.Series)
    out = pd.concat([out, boundary_df], axis=1)

    out["CAG_CTG"] = df["intron_full"].apply(cag_ctg)

    pair_df = df["intron_full"].apply(compute_pair_cols).apply(pd.Series)
    out = pd.concat([out, pair_df], axis=1)

    secondary_df = df.apply(
        lambda row: compute_secondary_structure_cols(
            row["intron_first20"],
            row["intron_last20"],
            max_pairs_score_max=max_pairs_score_max,
            total_pairs_score_max=total_pairs_score_max,
        ),
        axis=1,
    ).apply(pd.Series)

    observed_max_pairs_score = secondary_df["max_pairs_score"].max()
    observed_total_pairs_score = secondary_df["total_pairs_score"].max()

    if observed_max_pairs_score > max_pairs_score_max:
        print(
            "WARNING: observed max_pairs_score "
            f"({observed_max_pairs_score}) exceeds "
            f"configured max_pairs_score_max "
            f"({max_pairs_score_max}). "
            "Normalized values above the limit "
            "were clipped to 1."
        )

    if observed_total_pairs_score > total_pairs_score_max:
        print(
            "WARNING: observed total_pairs_score "
            f"({observed_total_pairs_score}) exceeds "
            f"configured total_pairs_score_max "
            f"({total_pairs_score_max}). "
            "Normalized values above the limit "
            "were clipped to 1."
        )

    out = pd.concat([out, secondary_df], axis=1)

    sequence_score_cols = [
        "prev_exon_Y",
        "5'_intron_R",
        "3'_intron_Y",
        "next_exon_R",
        "CAG_CTG",
        "4_-6",
        "5_-7",
        "6_-8",
        "CG_4_-6",
        "AT_5_-7",
        "GC_6_-8",
    ]

    out["SUM_sequence"] = out[sequence_score_cols].sum(axis=1)

    structure_score_cols = [
        "max_pairs_normalized",
        "max_pairs_score_normalized",
        "total_pairs_normalized",
        "total_pairs_score_normalized",
    ]

    out["SUM_structure"] = out[structure_score_cols].sum(axis=1).round(2)

    conventional_df = df["intron_full"].apply(compute_conventional_cols).apply(pd.Series)
    out = pd.concat([out, conventional_df], axis=1)

    conventional_score_cols = [
        "G_1",
        "T_2",
        "C_2",
        "G_-1",
        "A_-2",
        "C_-3",
        "pyrimidines_-3_to_-12",
    ]

    out["SUM_conventional"] = out[conventional_score_cols].sum(axis=1).round(2)

    out["conv_or_nonconv"] = df["conv_or_nonconv"].astype(str)
    out["const_or_alt"] = df["const_or_alt"].astype(str)
    out["intron_full"] = df["intron_full"].astype(str)
    out["length"] = df["length"].astype(int)

    if include_full:
        out["prev_exon_last5"] = df["prev_exon_last5"].astype(str)
        out["next_exon_first5"] = df["next_exon_first5"].astype(str)
        out["intron_first20"] = df["intron_first20"].astype(str)
        out["intron_last20"] = df["intron_last20"].astype(str)

    short_df = out[out["length"] < short_intron_length].copy()
    out = out[out["length"] >= short_intron_length].copy()

    out = out[cols]

    if not short_df.empty:
        short_df = short_df[cols]

    out.to_csv(out_csv, index=False)

    if short_out_csv is not None and not short_df.empty:
        short_df.to_csv(short_out_csv, index=False)

    return {
        "saved_rows": len(out),
        "short_rows": len(short_df),
        "out_csv": str(out_csv),
        "short_out_csv": str(short_out_csv)
        if short_out_csv is not None and not short_df.empty
        else None,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inp")
    ap.add_argument("-o", "--out", default="scored.csv")
    args = ap.parse_args()

    run_scoring(args.inp, args.out)


if __name__ == "__main__":
    main()
