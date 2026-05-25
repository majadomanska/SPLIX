import pandas as pd


FEATURE_DETAIL_COLUMNS = [
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
    "G_1",
    "T_2",
    "C_2",
    "G_-1",
    "A_-2",
    "C_-3",
]

FRACTION_DETAIL_COLUMNS = [
    "pyrimidines_-3_to_-12",
    "max_pairs_normalized",
    "max_pairs_score_normalized",
    "total_pairs_normalized",
    "total_pairs_score_normalized",
]

RAW_DETAIL_COLUMNS = [
    "max_pairs",
    "max_pairs_score",
    "total_pairs",
    "total_pairs_score",
]

DETAIL_COLUMNS = (
    FEATURE_DETAIL_COLUMNS
    + FRACTION_DETAIL_COLUMNS
    + RAW_DETAIL_COLUMNS
)


def _normalize(x):
    """Normalize value to stripped string."""
    if pd.isna(x):
        return ""
    return str(x).strip()


def _prepare_df(path):
    """Load and prepare scoring dataframe."""
    df = pd.read_csv(path, sep=None, engine="python", dtype=str)
    df.columns = [str(c).strip() for c in df.columns]

    required = [
        "SUM_sequence_normalized",
        "SUM_structure_normalized",
        "SUM_conventional_normalized",
        "conv_or_nonconv",
        "const_or_alt",
    ]

    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns: {missing}")

    df["SUM_sequence_normalized"] = pd.to_numeric(df["SUM_sequence_normalized"], errors="coerce")
    df["SUM_structure_normalized"] = pd.to_numeric(df["SUM_structure_normalized"], errors="coerce")
    df["SUM_conventional_normalized"] = pd.to_numeric(df["SUM_conventional_normalized"], errors="coerce")

    df["conv_or_nonconv"] = df["conv_or_nonconv"].map(_normalize).str.upper()
    df["const_or_alt"] = df["const_or_alt"].map(_normalize).str.lower()

    return df


def _is_alt(s):
    """Select alternative introns."""
    return s.isin(["alt"])


def _is_const(s):
    """Select constitutive introns."""
    return s.isin(["const"])


def _normalized_percent(subdf, score_col):
    """Compute mean normalized score percentage."""
    if len(subdf) == 0:
        return None

    values = pd.to_numeric(subdf[score_col], errors="coerce").dropna()

    if len(values) == 0:
        return None

    return values.mean() * 100.0


def _share_percent(part, whole):
    """Compute percentage share."""
    if whole == 0:
        return None
    return (part / whole) * 100.0


def _fmt(x):
    """Format percentage value."""
    if x is None:
        return "NA"
    return f"{x:.2f}%"


def _collect_groups(df):
    """Split dataframe into intron groups."""
    alt = df[_is_alt(df["const_or_alt"])]
    const = df[_is_const(df["const_or_alt"])]

    conv = df[df["conv_or_nonconv"] == "C"]
    nonconv = df[df["conv_or_nonconv"] == "N"]

    conv_alt = conv[_is_alt(conv["const_or_alt"])]
    conv_const = conv[_is_const(conv["const_or_alt"])]

    nonconv_alt = nonconv[_is_alt(nonconv["const_or_alt"])]
    nonconv_const = nonconv[_is_const(nonconv["const_or_alt"])]

    return {
        "alt": alt,
        "const": const,
        "conv": conv,
        "nonconv": nonconv,
        "conv_alt": conv_alt,
        "conv_const": conv_const,
        "nonconv_alt": nonconv_alt,
        "nonconv_const": nonconv_const,
    }


def _collect_scoring(groups, score_col):
    """Collect scoring statistics for all intron groups."""
    return {
        "conv_vs_nonconv": (
            _normalized_percent(groups["conv"], score_col),
            _normalized_percent(groups["nonconv"], score_col),
        ),
        "within_conv": (
            _normalized_percent(groups["conv_alt"], score_col),
            _normalized_percent(groups["conv_const"], score_col),
        ),
        "within_nonconv": (
            _normalized_percent(groups["nonconv_alt"], score_col),
            _normalized_percent(groups["nonconv_const"], score_col),
        ),
    }


def _collect_composition(groups):
    """Collect composition statistics for all intron groups."""
    n_conv = len(groups["conv"])
    n_nonconv = len(groups["nonconv"])
    n_conv_alt = len(groups["conv_alt"])
    n_conv_const = len(groups["conv_const"])
    n_nonconv_alt = len(groups["nonconv_alt"])
    n_nonconv_const = len(groups["nonconv_const"])

    return {
        "conv_vs_nonconv": (
            _share_percent(n_conv, n_conv + n_nonconv),
            _share_percent(n_nonconv, n_conv + n_nonconv),
        ),
        "within_conv": (
            _share_percent(n_conv_alt, n_conv_alt + n_conv_const),
            _share_percent(n_conv_const, n_conv_alt + n_conv_const),
        ),
        "within_nonconv": (
            _share_percent(n_nonconv_alt, n_nonconv_alt + n_nonconv_const),
            _share_percent(n_nonconv_const, n_nonconv_alt + n_nonconv_const),
        ),
    }


def _build_txt(sequence_scoring, structure_scoring, conventional_scoring, composition):
    """Build final TXT report.""" 

    return f"""SEQUENCE SCORING

Conventional vs Non-conventional
  Conventional:      {_fmt(sequence_scoring["conv_vs_nonconv"][0])}
  Non-conventional:  {_fmt(sequence_scoring["conv_vs_nonconv"][1])}

Conventional introns
  Alternative:   {_fmt(sequence_scoring["within_conv"][0])}
  Constitutive:  {_fmt(sequence_scoring["within_conv"][1])}

Non-conventional introns
  Alternative:   {_fmt(sequence_scoring["within_nonconv"][0])}
  Constitutive:  {_fmt(sequence_scoring["within_nonconv"][1])}


STRUCTURE SCORING

Conventional vs Non-conventional
  Conventional:      {_fmt(structure_scoring["conv_vs_nonconv"][0])}
  Non-conventional:  {_fmt(structure_scoring["conv_vs_nonconv"][1])}

Conventional introns
  Alternative:   {_fmt(structure_scoring["within_conv"][0])}
  Constitutive:  {_fmt(structure_scoring["within_conv"][1])}

Non-conventional introns
  Alternative:   {_fmt(structure_scoring["within_nonconv"][0])}
  Constitutive:  {_fmt(structure_scoring["within_nonconv"][1])}

  
CONVENTIONAL SCORING

Conventional vs Non-conventional
  Conventional:      {_fmt(conventional_scoring["conv_vs_nonconv"][0])}
  Non-conventional:  {_fmt(conventional_scoring["conv_vs_nonconv"][1])}

Conventional introns
  Alternative:   {_fmt(conventional_scoring["within_conv"][0])}
  Constitutive:  {_fmt(conventional_scoring["within_conv"][1])}

Non-conventional introns
  Alternative:   {_fmt(conventional_scoring["within_nonconv"][0])}
  Constitutive:  {_fmt(conventional_scoring["within_nonconv"][1])}

  
COMPOSITION

Conventional vs Non-conventional
  Conventional:      {_fmt(composition["conv_vs_nonconv"][0])}
  Non-conventional:  {_fmt(composition["conv_vs_nonconv"][1])}

Conventional introns
  Alternative:   {_fmt(composition["within_conv"][0])}
  Constitutive:  {_fmt(composition["within_conv"][1])}

Non-conventional introns
  Alternative:   {_fmt(composition["within_nonconv"][0])}
  Constitutive:  {_fmt(composition["within_nonconv"][1])}
"""


def _detail_value(subdf, col):
    """Compute detailed statistic for one column."""
    if len(subdf) == 0 or col not in subdf.columns:
        return None

    values = pd.to_numeric(subdf[col], errors="coerce").dropna()

    if len(values) == 0:
        return None

    if col in FEATURE_DETAIL_COLUMNS:
        return (values > 0).mean() * 100.0

    if col in FRACTION_DETAIL_COLUMNS:
        return values.mean() * 100.0

    if col in RAW_DETAIL_COLUMNS:
        return values.mean()

    return None


def _fmt_detail(x, col):
    """Format detailed statistic."""
    if x is None:
        return "NA"

    if col in RAW_DETAIL_COLUMNS:
        return f"{x:.2f}"

    return f"{x:.2f}%"


def _build_detailed_table(title, groups_dict):
    """Build detailed statistics table for selected groups."""
    lines = [
        "",
        title,
        "",
    ]

    group_names = list(groups_dict.keys())
    header = "Feature".ljust(32)

    for name in group_names:
        header += name.rjust(18)

    lines.append(header)
    lines.append("-" * len(header))

    detail_sections = [
        ("SEQUENCE", [
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
        ]),
        ("STRUCTURE", [
            "max_pairs_normalized",
            "max_pairs_score_normalized",
            "total_pairs_normalized",
            "total_pairs_score_normalized",
            "max_pairs",
            "max_pairs_score",
            "total_pairs",
            "total_pairs_score",
        ]),
        ("CONVENTIONAL", [
            "G_1",
            "T_2",
            "C_2",
            "G_-1",
            "A_-2",
            "C_-3",
            "pyrimidines_-3_to_-12",
        ])
    ]

    for section_name, columns in detail_sections:
        if not columns:
            continue

        lines.append("")
        lines.append(section_name)

        for col in columns:
            row = col.ljust(32)

            for name in group_names:
                value = _detail_value(groups_dict[name], col)
                row += _fmt_detail(value, col).rjust(18)

            lines.append(row)

    return "\n".join(lines)


def _build_detailed_txt(groups):
    """Build detailed statistics TXT tables."""
    sections = []

    sections.append(
        _build_detailed_table(
            "DETAILED STATS: Conventional vs Non-conventional",
            {
                "Conventional": groups["conv"],
                "Non-conventional": groups["nonconv"],
            },
        )
    )

    sections.append(
        _build_detailed_table(
            "DETAILED STATS: Conventional introns",
            {
                "Alternative": groups["conv_alt"],
                "Constitutive": groups["conv_const"],
            },
        )
    )

    sections.append(
        _build_detailed_table(
            "DETAILED STATS: Non-conventional introns",
            {
                "Alternative": groups["nonconv_alt"],
                "Constitutive": groups["nonconv_const"],
            },
        )
    )

    return "\n\n".join(sections)


def run_stats(inp_csv, out_txt, detailed=True):
    """Run statistics analysis and save TXT report."""
    df = _prepare_df(inp_csv)
    groups = _collect_groups(df)

    sequence_scoring = _collect_scoring(groups, "SUM_sequence")
    structure_scoring = _collect_scoring(groups, "SUM_structure")
    conventional_scoring = _collect_scoring(groups, "SUM_conventional")

    composition = _collect_composition(groups)

    txt = _build_txt(sequence_scoring, structure_scoring, conventional_scoring, composition)

    if detailed:
        txt += "\n\n"
        txt += _build_detailed_txt(groups)

    with open(out_txt, "w", encoding="utf-8") as f:
        f.write(txt)

    return {
        "out_txt": str(out_txt),
        "sequence_scoring": sequence_scoring,
        "structure_scoring": structure_scoring,
        "conventional_scoring": conventional_scoring,
        "composition": composition,
        "detailed": detailed,
    }