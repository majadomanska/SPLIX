# NAZWA #TODO

A bioinformatics pipeline for extraction, classification, and scoring of introns from GFF and FASTA files.

This tool is designed for analysis of intron features associated with conventional and non-conventional introns, including sequence composition, splice-site motifs, and RNA secondary structure. 

Developed and tested on Linux.

The pipeline consists of three main steps:
1. EXTRACT
2. SCORE
3. STATS

This pipeline:
1. extracts introns from genomic annotations,
2. classifies introns,
3. computes sequence-based scores,
4. predicts RNA secondary structure using RNAfold,
5. generates summary statistics.

---

# Features

- intron extraction from GFF annotations
- strand-aware sequence extraction
- support for alternative and constitutive introns
- support for conventional and non-conventional introns
- splice-site motif scoring
- RNA secondary structure analysis using RNAfold
- configurable filtering system
- automatic statistics generation
- detailed pipeline logs

---

# Requirements

- Python >= 3.10
- RNAfold from the ViennaRNA package

Python dependencies are listed in `requirements.txt`.

---

# Installation

Clone repository:

```bash
git clone https://github.com/majadomanska/introns.git
cd introns
```

(Optional but recommended) create virtual environment:

```bash
python3 -m venv venv
source venv/bin/activate
```

Install Python dependencies:

```bash
pip install -r requirements.txt
```

---

# RNAfold Dependency

This project requires RNAfold from the ViennaRNA package.

Ubuntu / Debian:

```bash
sudo apt install vienna-rna     #NIE WIEM CZY DZIALA TODO
```

Check installation:

```bash
RNAfold --version
```

---

# Input Files

This pipeline requires:

- FASTA file containing genomic sequences
- GFF file containing genome annotations

---

# Input Requirements

Scaffold names in FASTA and GFF files must match exactly.

Example:

FASTA:
```text
>scaffold_1
```

GFF:
```text
scaffold_1
```

Input FASTA and GFF files are expected to be located in the current working directory unless full or relative paths are provided.

Differences in scaffold naming may cause introns to be skipped during extraction.

The FASTA file may contain additional scaffolds not present in the GFF file.
Only scaffolds present in the GFF annotation are analyzed.

---

# Usage

Basic usage:

```bash
python pipeline.py \ # z pyproject.toml bedzie ladniej #TODO
    --fasta genome.fasta \
    --gff annotations.gff
```

Using custom config:

```bash
python pipeline.py \
    --fasta genome.fasta \
    --gff annotations.gff \
    --config config.yaml
```

Using custom output name:

```bash
python pipeline.py \
    --fasta genome.fasta \
    --gff annotations.gff \
    --name example_run
```

---

# Configuration

Example configuration file:

```yaml
outputs:
  save_scored: true
  save_extracted: false

scored:
  include_full: true
  short_intron_length: 30

filters:
  conv_or_nonconv: null       # null / "C" / "N"
  const_or_alt: null          # null / "const" / "alt"

stats:
  detailed: true
  
structure_scoring:
  max_pairs_score_max: 48 # range: 0-60
  total_pairs_score_max: 50 # range: 0-60 
```

---

# Configuration Parameters

## outputs

| Parameter | Description |
|---|---|
| save_scored | Save scored intron table |
| save_extracted | Save extracted introns before scoring |

---

## scored

| Parameter | Description |
|---|---|
| include_full | Include additional sequence and structure columns in output |
| short_intron_length | Minimum intron length required for scoring |

---

## structure_scoring

| Parameter | Description |
|---|---|
| max_pairs_score_max | Maximum value used to normalize `max_pairs_score` |
| total_pairs_score_max | Maximum value used to normalize `total_pairs_score` |

Default values:

```yaml
structure_scoring:
  max_pairs_score_max: 48
  total_pairs_score_max: 50
```

These defaults were selected based on observed score distributions in the *E. longa* genome dataset.

If observed scores exceed configured normalization limits, the pipeline reports a warning and values above the limit are clipped to `1` during normalization.

Example warning:

```text
WARNING: observed total_pairs_score (40) exceeds configured total_pairs_score_max (10). Normalized values above the limit were clipped to 1.
```

---

## filters

| Parameter | Description |
|---|---|
| conv_or_nonconv | Filter introns by type (`C` or `N`) |
| const_or_alt | Filter introns by class (`const` or `alt`) |

---

## stats

| Parameter | Description |
|---|---|
| detailed | Generate additional detailed statistics |

---

# Output Files

This pipeline creates a results directory containing:

| File | Description |
|---|---|
| *_extracted.csv | Extracted introns before scoring |
| *_scored.csv | Final scored introns |
| *_short_introns.csv | Introns filtered out due to short length |
| *_stats.txt | Summary statistics |
| *_log.txt | Pipeline execution log |

If `save_extracted: false`, the extracted intron table is not saved.

---

# Output Columns

The following columns are additionally included when:

```yaml
scored:
  include_full: true
```

Additional columns:

| Column | Description |
|---|---|
| length | Intron length |
| prev_exon_last5 | Last 5 nt of upstream exon |
| next_exon_first5 | First 5 nt of downstream exon |
| intron_first20 | First 20 nt of intron |
| intron_last20 | Last 20 nt of intron |
| secondary_structure | RNAfold secondary structure prediction |

---

## General Columns

| Column | Description |
|---|---|
| scaffold | Scaffold name |
| id | Unique intron identifier based on scaffold, coordinates, and strand |
| name | Original intron or transcript name |
| length | Intron length |
| conv_or_nonconv | Conventional (`C`) or non-conventional (`N`) intron |
| const_or_alt | Constitutive (`const`) or alternative (`alt`) intron |
| intron_full | Full intron sequence |

---

## Sequence Context Columns

| Column | Description |
|---|---|
| prev_exon_last5 | Last 5 nt of upstream exon |
| next_exon_first5 | First 5 nt of downstream exon |
| intron_first20 | First 20 nt of intron |
| intron_last20 | Last 20 nt of intron |

---

## Sequence Scoring Columns

These columns describe sequence-based splice-site and motif features.

Most sequence features are binary:
- `1` = feature present,
- `0` = feature absent.

| Column | Description |
|---|---|
| prev_exon_Y | Pyrimidine at upstream exon boundary (0 or 1) |
| `5'_intron_R` | Purine at intron 5' boundary (0 or 1) |
| `3'_intron_Y` | Pyrimidine at intron 3' boundary (0 or 1) |
| next_exon_R | Purine at downstream exon boundary (0 or 1) |
| CAG_CTG | Presence of CAG/CTG motifs (0 or 1) |
| 4_-6 | Pairing between positions +4 and -6 (0 or 1) |
| 5_-7 | Pairing between positions +5 and -7 (0 or 1) |
| 6_-8 | Pairing between positions +6 and -8 (0 or 1) |
| CG_4_-6 | C-G pairing at +4/-6 (0 or 1) |
| AT_5_-7 | A-T pairing at +5/-7 (0 or 1) |
| GC_6_-8 | G-C pairing at +6/-8 (0 or 1) |
| SUM_sequence | Combined sequence-based score |

---

## Structure Scoring Columns

These columns describe RNA secondary structure features predicted using RNAfold.

RNA base pairs are weighted according to pairing strength:
- `C-G` or `G-C` = 3 points
- `A-U` or `U-A` = 2 points
- `G-U` or `U-G` = 1 point

| Column | Description |
|---|---|
| secondary_structure | RNAfold secondary structure prediction |
| max_pairs | Maximum consecutive base-pair stretch (0–20) |
| max_pairs_score | Score derived from maximum pairing stretch based on RNA pair types (0–48) |
| total_pairs | Total number of detected pairings (0–20) |
| total_pairs_score | Total score of all detected pairings based on RNA pair types (0–50) |
| max_pairs_normalized | Normalized maximum pairing value scaled to 0–1 |
| max_pairs_score_normalized | Normalized maximum pairing score scaled to 0–1 |
| total_pairs_normalized | Normalized total pairing value scaled to 0–1 |
| total_pairs_score_normalized | Normalized total pairing score scaled to 0–1 |
| SUM_structure | Combined structure-based score |

---

## Conventional Intron Scoring Columns

These columns describe sequence features associated with conventional introns.

Most features are binary:
- `1` = feature present,
- `0` = feature absent.

Exceptions:
- `C_2` contributes `0.5` when present,
- `pyrimidines_-3_to_-12` is normalized to a value between `0` and `1` based on pyrimidine content.

| Column | Description |
|---|---|
| G_1 | Presence of G at position +1 (0 or 1) |
| T_2 | Presence of T at position +2 (0 or 1) |
| C_2 | Presence of C at position +2 (0 or 0.5) |
| G_-1 | Presence of G at position -1 (0 or 1) |
| A_-2 | Presence of A at position -2 (0 or 1) |
| C_-3 | Presence of C at position -3 (0 or 1) |
| pyrimidines_-3_to_-12 | Pyrimidine-rich region near the 3' splice site (0–1) |
| SUM_conventional | Combined conventional intron score |

---

# Scoring

The pipeline computes three major feature groups:

1. sequence-based scoring,
2. RNA secondary structure scoring,
3. scoring of sequence features associated with conventional introns.

These feature groups are summarized using:
- `SUM_sequence`
- `SUM_structure`
- `SUM_conventional`

---

# Statistics

The pipeline generates summary statistics for:
- alternative vs constitutive introns,
- conventional vs non-conventional introns,
- alternative vs constitutive introns within conventional introns,
- alternative vs constitutive introns within non-conventional introns,
- sequence scoring,
- structure scoring,
- scoring of sequence features associated with conventional introns.

When:

```yaml
stats:
  detailed: true
```

statistics for individual feature columns are additionally reported.

---

# Project Structure

```text
project/
├── pipeline.py
├── extract.py
├── score.py
├── stats.py
├── config.yaml
├── requirements.txt
└── README.md
```

---

# Example Workflow

```bash
python pipeline.py \
    --fasta example.fasta \
    --gff example.gff \
    --config config.yaml
```

Example output:

```text
example_results/
├── example_extracted.csv
├── example_scored.csv
├── example_short_introns.csv
├── example_stats.txt
└── example_log.txt
```

---

# Pipeline Log

The pipeline automatically generates a detailed log file containing execution parameters and summary statistics.

Example log file:

```text
Step 1: EXTRACT
FASTA file: 50scaffolds.fasta
GFF file: the_eight.gff
Config file: config.yaml
Mode: full
Short intron threshold: 150 nt
Filter: Conventional or Nonconventional: None
Filter: Constitutive or Alternative: None
Number of scaffolds in FASTA: 50
Number of scaffolds in GFF: 8
Total GFF records: 1826
Identified introns: 874
Unique introns after merging: 175
Saved extracted rows: 175
Extracted CSV: the_eight_results/the_eight_extracted.csv

Step 2: SCORE
Input CSV: the_eight_results/the_eight_extracted.csv
Saved scored rows: 161
Saved short introns: 14
Scored CSV: the_eight_results/the_eight_scored.csv

Step 3: STATS
Stats TXT: the_eight_results/the_eight_stats.txt
Detailed stats: True
Short introns CSV: the_eight_results/the_eight_short_introns.csv

Done!
Results in: the_eight_results
Elapsed time: 0.97 s
```

---

# Notes

- Introns shorter than the configured minimum length are excluded from scoring.
- RNAfold must be available in the system PATH.
- All sequence operations are strand-aware.

---

# License

#TODO

---

# Citation

#TODO
