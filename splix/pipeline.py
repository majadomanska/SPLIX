#!/usr/bin/env python3

import argparse
import time
from pathlib import Path

import yaml

from splix.extract import run_extract
from splix.score import run_scoring
from splix.stats import run_stats


def load_config(config_path):
    defaults = {
        "outputs": {
            "save_scored": True,
            "save_extracted": False,
        },
        "scored": {
            "include_full": True,
            "short_intron_length": 40,
        },
        "filters": {
            "conv_or_nonconv": None,
            "const_or_alt": None,
        },
        "stats": {
            "detailed": True,
        },
        "structure_scoring": {
            "max_pairs_score_max": None,
            "total_pairs_score_max": None,
        },
    }

    with open(config_path, "r", encoding="utf-8") as f:
        user_cfg = yaml.safe_load(f) or {}

    for section, values in user_cfg.items():
        if section in defaults and isinstance(values, dict):
            defaults[section].update(values)
        else:
            defaults[section] = values

    return defaults


def main():
    start_time = time.time()
    default_config = Path(__file__).resolve().parent / "config.yaml"

    ap = argparse.ArgumentParser(description="Score introns")
    ap.add_argument("--gff", required=True, help="GFF file")
    ap.add_argument("--fasta", required=True, help="FASTA file")
    ap.add_argument("--config", default=str(default_config), help="Configuration file")
    ap.add_argument("--name", default=None, help="Custom base name for output files")
    args = ap.parse_args()

    gff = Path(args.gff)
    fasta = Path(args.fasta)
    cfg = load_config(args.config)

    base = args.name if args.name else gff.stem
    result_dir = Path(f"{base}_results")
    result_dir.mkdir(parents=True, exist_ok=True)

    extracted_csv = result_dir / f"{base}_extracted.csv"
    scored_csv = result_dir / f"{base}_scored.csv"
    short_scored_csv = result_dir / f"{base}_short_introns.csv"
    log_txt = result_dir / f"{base}_log.txt"

    # Step 1: Exctract
    print("Step 1: EXTRACT")

    extract_stats = run_extract(
        gff_path=str(gff),
        fasta_path=str(fasta),
        out_csv=str(extracted_csv),
    )

    print("Step 1: DONE")

    # Step 2: Score
    print("Step 2: SCORE")

    # Remove old short introns file.
    if short_scored_csv.exists():
        short_scored_csv.unlink()

    # Structure scoring configuration.
    structure_cfg = cfg.get("structure_scoring", {})

    max_pairs_score_max = structure_cfg.get("max_pairs_score_max", None,)

    total_pairs_score_max = structure_cfg.get("total_pairs_score_max", None,)

    # Validate configuration values.
    if (max_pairs_score_max is not None and not 1 <= max_pairs_score_max <= 60):
        raise ValueError("max_pairs_score_max must be between 1 and 60")

    if (total_pairs_score_max is not None and not 1 <= total_pairs_score_max <= 60):
        raise ValueError("total_pairs_score_max must be between 1 and 60")

    score_stats = run_scoring(
        inp_csv=str(extracted_csv),
        out_csv=str(scored_csv),
        include_full=cfg["scored"]["include_full"],
        filter_c_lub_n=cfg["filters"]["conv_or_nonconv"],
        filter_konstyt=cfg["filters"]["const_or_alt"],
        short_out_csv=str(short_scored_csv),
        short_intron_length=cfg["scored"]["short_intron_length"],
        max_pairs_score_max=max_pairs_score_max,
        total_pairs_score_max=total_pairs_score_max,
    )

    print("Step 2: DONE")

    # Step 3: Stats
    print("Step 3: STATS")

    stats_txt = result_dir / f"{base}_stats.txt"

    stats_stats = run_stats(
        inp_csv=str(scored_csv),
        out_txt=str(stats_txt),
        detailed=cfg["stats"]["detailed"],
    )

    print("Step 3: DONE")
    
    # Handle output files based on confgiguration.
    if not cfg["outputs"]["save_extracted"] and extracted_csv.exists():
        extracted_csv.unlink()

    if not cfg["outputs"]["save_scored"]:
        if scored_csv.exists():
            scored_csv.unlink()
        if short_scored_csv.exists():
            short_scored_csv.unlink()

    # LOG
    elapsed = time.time() - start_time
    mode = "full" if cfg["scored"]["include_full"] else "base"

    log_lines = [
        "Step 1: EXTRACT",
        f"FASTA file: {fasta.name}",
        f"GFF file: {gff.name}",
        f"Config file: {Path(args.config).name}",
        f"Mode: {mode}",
        f"Short intron threshold: {cfg['scored']['short_intron_length']} nt",
        f"Filter: Conventional or Nonconventional: {cfg['filters']['conv_or_nonconv']}",
        f"Filter: Constitutive or Alternative: {cfg['filters']['const_or_alt']}",
        f"Number of scaffolds in FASTA: {extract_stats['num_scaffolds']}",
        f"Number of scaffolds in GFF: {extract_stats['num_scaffolds_gff']}",
        f"Total GFF records: {extract_stats['total_gff_records']}",
        f"Identified introns: {extract_stats['identified_introns']}",
        f"Unique introns after merging: {extract_stats['unique_introns']}",
        f"Saved extracted rows: {extract_stats['saved_rows']}",
        f"Extracted CSV: {extract_stats['out_csv']}",
        f"Scaffolds missing in FASTA: {len(extract_stats['missing_scaffolds_in_fasta'])}",
        "",
        "Step 2: SCORE",
        f"Input CSV: {extracted_csv}",
        f"Saved scored rows: {score_stats['saved_rows']}",
        f"Saved short introns: {score_stats['short_rows']}",
        f"Scored CSV: {score_stats['out_csv']}",
        f"max_pairs_score_max: {score_stats['max_pairs_score_max']}",
        f"total_pairs_score_max: {score_stats['total_pairs_score_max']}",
    ]

    if score_stats["short_out_csv"] is not None:
        log_lines.append(f"Short introns CSV: {score_stats['short_out_csv']}")
    else:
        log_lines.append("Short introns CSV: None")

    all_warnings = []

    all_warnings.extend(extract_stats.get("warnings", []))
    all_warnings.extend(score_stats.get("warnings", []))

    if all_warnings:
        log_lines.append("")
        log_lines.append("WARNINGS")

        for warning in all_warnings:
            log_lines.append(warning)

    log_lines.extend([
        "",
        "Step 3: STATS",
        f"Stats TXT: {stats_stats['out_txt']}",
        f"Detailed stats: {cfg['stats']['detailed']}",
        "",
        "Done!",
        f"Results in: {result_dir}",
        f"Elapsed time: {elapsed:.2f} s",
    ])

    with open(log_txt, "w", encoding="utf-8") as f:
        f.write("\n".join(str(x) for x in log_lines) + "\n")

    print(f"\nSuccess! Results in: {result_dir}")
    print(f"Log file: {log_txt}")


if __name__ == "__main__":
    main()
