# SPLIX

A bioinformatics pipeline for extraction, classification, and scoring of introns from GFF and FASTA files.

SPLIX is a tool designed for analysis of intron features associated with conventional and non-conventional introns, including sequence composition, splice-site motifs and RNA pairing potential between intron termini estimated using RNAfold.

Developed and tested on Linux.  #TODO test na innych - dodac 

This pipeline:
1. extracts introns from genomic annotations,
2. classifies introns,
3. computes sequence-based scores,
4. estimates splice-site pairing potential using RNAfold,
5. computes RNAfold-based pairing scores,
6. generates summary statistics.

---

# Features

- intron extraction from GFF annotations
- strand-aware sequence extraction
- support for alternative and constitutive introns
- support for conventional and non-conventional introns
- splice-site motif scoring
- RNAfold-based splice-site pairing analysis
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
git clone https://github.com/majadomanska/SPLIX.git
cd SPLIX
```

(Optional but recommended) create virtual environment:

```bash
python3 -m venv splix_env
source splix_env/bin/activate
```

Install SPLIX:

```bash
pip install -e . 
```


Alternatively, SPLIX may also be used directly from the source code without installation:

```bash
python splix/pipeline.py \
    --fasta genome.fasta \
    --gff annotations.gff
```

In this case, Python dependencies may be installed using:

```bash
pip install -r requirements.txt
```
---

# RNAfold Dependency

This project requires RNAfold from the ViennaRNA package.

Ubuntu / Debian:

```bash
sudo apt install vienna-rna
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

Input FASTA and GFF files may be specified using relative or absolute paths.

Differences in scaffold naming will cause introns to be skipped during extraction.

The program is strand-aware and processes both + and - strand introns correctly.
However, the strand orientation in the GFF annotation must match the orientation of sequences in the FASTA file.
If reversed annotations are used, a correspondingly reversed FASTA file must also be provided.

The FASTA file may contain additional scaffolds not present in the GFF file.
Only scaffolds present in the GFF annotation are analyzed.

---

# Usage

Basic usage:

```bash
splix \
    --fasta genome.fasta \
    --gff annotations.gff
```

Using custom config:

```bash
splix \
    --fasta genome.fasta \
    --gff annotations.gff \
    --config config.yaml
```

Using custom output name:

```bash
splix \
    --fasta genome.fasta \
    --gff annotations.gff \
    --name example_run
```

---

# Configuration


By default, SPLIX uses the bundled `config.yaml` file located inside the package.
Custom configuration files may be provided using `--config`.

Example configuration file:

```yaml
outputs:
  save_scored: true
  save_extracted: false

scored:
  include_full: true
  short_intron_length: 40

filters:
  conv_or_nonconv: null       # null / "C" / "N"
  const_or_alt: null          # null / "const" / "alt"

stats:
  detailed: true
  
structure_scoring:
  max_pairs_score_max: null # range: 1-60
  total_pairs_score_max: null # range: 1-60 
```

---

# Configuration Parameters

## outputs

| Parameter | Description |
|---|---|
| save_scored | Save scored intron table - main output, `default: true`|
| save_extracted | Save extracted introns before scoring, `default: false`|

---

## scored

| Parameter | Description |
|---|---|
| include_full | Include additional sequence and structure columns in output, `default: true` |
| short_intron_length | Minimum intron length required for scoring, `default: 40` |

---

## structure_scoring

| Parameter | Description |
|---|---|
| max_pairs_score_max | Maximum value used to normalize `max_pairs_score`|
| total_pairs_score_max | Maximum value used to normalize `total_pairs_score`|

Default values:

```yaml
structure_scoring:
  max_pairs_score_max: null
  total_pairs_score_max: null
```

When set to `null`, the pipeline automatically uses the highest observed score value from the dataset for normalization.
This allows normalization to adapt automatically to datasets with different score ranges without requiring manual parameter tuning.

Users may also define fixed normalization limits manually in the config file (range: 1–60).

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
| intron_first20 | First 20 nt of intron |
| intron_last20 | Last 20 nt of intron |
| next_exon_first5 | First 5 nt of downstream exon |
| secondary_structure | RNAfold secondary structure prediction |
| secondary_structure_16 | Rescue RNAfold pairing structure calculated using 16 nt terminal windows when no pairing was detected with the default 20 nt windows |

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
| intron_first20 | First 20 nt of intron |
| intron_last20 | Last 20 nt of intron |
| next_exon_first5 | First 5 nt of downstream exon |

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
| SUM_sequence | Combined sequence-based score (0-11) |
| SUM_sequence_normalized | Normalized combined sequence-based score scaled to 0–1|

---

## RNAfold Pairing Strategy

The pipeline does not predict the full intron secondary structure.

Instead, it estimates pairing potential between intron termini using constrained local folding with RNAfold.

For each intron:
- the first 20 nt of the intron are extracted,
- the last 20 nt of the intron are extracted,
- a spacer of 20 `N` nucleotides is inserted between them.

The analyzed sequence therefore has the following structure:

```text
[first 20 nt] + NNNNNNNNNNNNNNNNNNNN + [last 20 nt]
```
If no pairing is detected using the default 20 nt terminal windows (`max_pairs = 0` and `total_pairs = 0`), SPLIX performs an additional rescue folding step using 16 nt terminal windows.

---

## Structure Scoring Columns

These columns describe splice-site pairing features estimated using RNAfold-based local folding analysis.

RNA base pairs are weighted according to pairing strength:
- `C-G` or `G-C` = 3 points
- `A-U` or `U-A` = 2 points
- `G-U` or `U-G` = 1 point

| Column | Description |
|---|---|
| secondary_structure | RNAfold dot-bracket pairing prediction |
| max_pairs | Maximum consecutive base-pair stretch (0–20) |
| max_pairs_score | Score derived from maximum pairing stretch based on RNA pair types (0–60) |
| total_pairs | Total number of detected pairings (0–20) |
| total_pairs_score | Total score of all detected pairings based on RNA pair types (0–60) |
| max_pairs_normalized | Normalized maximum pairing value scaled to 0–1 |
| max_pairs_score_normalized | Normalized maximum pairing score scaled to 0–1 |
| total_pairs_normalized | Normalized total pairing value scaled to 0–1 |
| total_pairs_score_normalized | Normalized total pairing score scaled to 0–1 |
| SUM_structure | Combined structure-based score (0-4)|
| SUM_structure_normalized | Combined structure-based score scaled to 0–1|


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
| pyrimidines_-3_to_-12 | Pyrimidine-rich region near the 3' splice site - number of pyrimidines in positions from -3 to -12 normalized to a value between 0 and 1 (0–1) |
| SUM_conventional | Combined conventional intron score (0-6)|
| SUM_conventional_normalized | Combined conventional intron score scaled to 0–1|

---

# Scoring

The pipeline computes three major feature groups:

1. sequence-based scoring,
2. RNAfold-based splice-site pairing scoring,
3. scoring of sequence features associated with conventional introns.

These feature groups are summarized using:
- `SUM_sequence` and `SUM_sequence_normalized`
- `SUM_structure` and `SUM_structure_normalized`
- `SUM_conventional` and `SUM_conventional_normalized`

---

# Statistics

The pipeline generates summary statistics for four main result groups:

- sequence scoring,
- RNAfold-based structure scoring,
- conventional intron feature scoring,
- intron composition.

Summary statistics are reported for:
- conventional vs non-conventional introns,
- alternative vs constitutive introns within conventional introns,
- alternative vs constitutive introns within non-conventional introns.

In summary scoring sections, percentages represent the mean normalized score for a given intron group, multiplied by 100.

For example:

```text
SEQUENCE SCORING

Conventional vs Non-conventional
  Conventional:      28.62%
  Non-conventional:  61.02%
```

means that conventional introns have an average `SUM_sequence_normalized` value of `0.2862`, while non-conventional introns have an average value of `0.6102`.

The `COMPOSITION` section is different. Its percentages describe the proportion of introns in each group, not scoring values.

When:

```yaml
stats:
  detailed: true
```

statistics for individual feature columns are additionally reported.

In detailed statistics:
- binary features are reported as percentages, showing how often a feature is present in a given group,
- normalized score columns are also reported as percentages,
- raw numeric columns are reported as mean values, not percentages.

For example:
- `G_1 = 100.00%` means that all introns in the group have `G` at position +1,
- `max_pairs_normalized = 16.54%` means that the mean normalized value is `0.1654`,
- `max_pairs = 3.31` means that the average maximum consecutive pairing stretch is `3.31` base pairs.

Detailed statistics are grouped into:
- `SEQUENCE` features,
- `STRUCTURE` features,
- `CONVENTIONAL` intron features.

---

# Project Structure

```text
SPLIX/
├── pyproject.toml
├── README.md
├── requirements.txt
└── splix/
    ├── __init__.py
    ├── config.yaml
    ├── pipeline.py
    ├── extract.py
    ├── score.py
    └── stats.py
```

---

# Example Workflow

```bash
splix \
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
Scaffolds missing in FASTA: 0

Step 2: SCORE
Input CSV: the_eight_results/the_eight_extracted.csv
Saved scored rows: 161
Saved short introns: 14
Scored CSV: the_eight_results/the_eight_scored.csv

Step 3: STATS
Stats TXT: the_eight_results/the_eight_stats.txt
Detailed stats: True
max_pairs_score_max: 48
total_pairs_score_max: 50
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

