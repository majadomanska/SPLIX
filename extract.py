#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# Data class representing one parsed GFF record.
@dataclass
class GFFRecord:
    scaffold: str
    source: str
    feature: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    attributes_raw: str
    attrs: Dict[str, str]
    rec_id: str = ""
    parent: str = ""
    name: str = ""


# Data class representing one processed intron candidate.
@dataclass
class IntronCandidate:
    scaffold: str
    strand: str
    start: int
    end: int
    original_id: str
    original_name: str
    feature: str
    parent: str

    prev_exon_last5: str = ""
    next_exon_first5: str = ""
    intron_first20: str = ""
    intron_last20: str = ""
    intron_full: str = ""

    conv_or_nonconv: str = ""
    const_or_alt: str = ""

    warnings: List[str] = field(default_factory=list)
    support_count: int = 1

    intron_group_no: Optional[int] = None
    intron_class: str = ""  # single / alt
    alt_label: str = ""  # alt1 / alt2 / ...


    @property
    def unique_key(self):
        """Return a key used to identify duplicate introns."""
        return (self.scaffold, self.strand, self.start, self.end)


    @property
    def stable_id(self):
        """Return a stable intron identifier based on genomic coordinates."""
        return f"{self.scaffold}:{self.start}-{self.end}:{self.strand}"


def parse_gff_attributes(attr_text):
    """Parse GFF attribute column into a dictionary."""
    attrs = {}
    text = attr_text.strip()

    if not text:
        return attrs

    for item in text.split(";"):
        item = item.strip()

        if not item:
            continue

        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key.strip()] = value.strip().strip('"')
        else:
            m = re.match(r'^(\S+)\s+"?(.*?)"?$', item)
            if m:
                key = m.group(1).strip()
                value = m.group(2).strip().strip('"')
                attrs[key] = value

    return attrs


def load_fasta(fasta_path):
    """Load FASTA sequences into a dictionary indexed by scaffold name."""
    genome = {}
    current_name = None
    seq_chunks = []

    with open(fasta_path) as fh:
        for line in fh:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if current_name is not None:
                    genome[current_name] = "".join(seq_chunks).upper()

                current_name = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)

    if current_name is not None:
        genome[current_name] = "".join(seq_chunks).upper()

    return genome


def read_gff(gff_path):
    """Read GFF file and convert each valid row into a GFFRecord object."""
    records = []

    with open(gff_path) as fh:
        for line in fh:
            line = line.rstrip("\n")

            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")

            if len(parts) != 9:
                continue

            scaffold, source, feature, start, end, score, strand, phase, attributes = parts
            attrs = parse_gff_attributes(attributes)

            records.append(
                GFFRecord(
                    scaffold=scaffold,
                    source=source,
                    feature=feature,
                    start=int(start),
                    end=int(end),
                    score=score,
                    strand=strand,
                    phase=phase,
                    attributes_raw=attributes,
                    attrs=attrs,
                    rec_id=attrs.get("ID", ""),
                    parent=attrs.get("Parent", ""),
                    name=attrs.get("Name", ""),
                )
            )

    return records


def is_intron_feature(feature):
    """Check whether a GFF feature represents an intron."""
    return "intron" in feature.strip().lower()


def extract_introns(records):
    """Extract only intron records."""
    return [rec for rec in records if is_intron_feature(rec.feature)]


def is_exon_feature(feature):
    """Check whether a GFF feature represents an exon."""
    return "exon" in feature.strip().lower()


def extract_exons(records):
    """Extract only exon records."""
    return [rec for rec in records if is_exon_feature(rec.feature)]


def extract_c_or_n(feature):
    """Extract type C or N from the feature name based on its suffix."""
    f = feature.strip()
    m = re.search(r'intron[_-]?([cn])$', f, flags=re.IGNORECASE)
    if m:
        return m.group(1).upper()
    return ""


def reverse_complement(seq):
    """Return the reverse-complemented DNA sequence."""
    trans = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(trans)[::-1]


def slice_genome_sequence(genome, scaffold, start, end):
    """Extract an exact genomic sequence fragment."""
    if scaffold not in genome:
        raise KeyError(f"Scaffold {scaffold!r} not found in FASTA")

    seq = genome[scaffold]

    if start < 1 or end > len(seq) or start > end:
        raise ValueError(
            f"Wrong coordinates for scaffold={scaffold}: "
            f"start={start}, end={end}, len={len(seq)}" 
        )

    return seq[start - 1:end].upper()


def safe_slice(genome, scaffold, start, end):
    """Extract a genomic fragment, clipping coordinates to sequence bounds."""
    if scaffold not in genome:
        return ""

    seq = genome[scaffold]
    start = max(1, start)
    end = min(len(seq), end)

    if start > end:
        return ""

    return seq[start - 1:end].upper()


def build_intron_candidate(intron, genome):
    """Build a processed intron candidate with oriented sequence fragments."""
    warnings = []

    intron_seq = slice_genome_sequence(genome, intron.scaffold, intron.start, intron.end)
    c_or_n = extract_c_or_n(intron.feature)

    if not c_or_n:
        warnings.append("missing_C_or_N")

    if intron.strand == "+":
        prev_exon_seq = safe_slice(genome, intron.scaffold, intron.start - 5, intron.start - 1)
        next_exon_seq = safe_slice(genome, intron.scaffold, intron.end + 1, intron.end + 5)
        intron_oriented = intron_seq

    elif intron.strand == "-":
        left_genomic = safe_slice(genome, intron.scaffold, intron.start - 5, intron.start - 1)
        right_genomic = safe_slice(genome, intron.scaffold, intron.end + 1, intron.end + 5)

        prev_exon_seq = reverse_complement(right_genomic)
        next_exon_seq = reverse_complement(left_genomic)
        intron_oriented = reverse_complement(intron_seq)
    else:
        warnings.append("unknown_strand")
        return None

    prev_exon_last5 = prev_exon_seq[-5:] if len(prev_exon_seq) >= 5 else prev_exon_seq
    next_exon_first5 = next_exon_seq[:5] if len(next_exon_seq) >= 5 else next_exon_seq
    intron_first20 = intron_oriented[:20] if len(intron_oriented) >= 20 else intron_oriented
    intron_last20 = intron_oriented[-20:] if len(intron_oriented) >= 20 else intron_oriented

    if len(prev_exon_seq) < 5:
        warnings.append("short_left_flank")
    if len(next_exon_seq) < 5:
        warnings.append("short_right_flank")

    return IntronCandidate(
        scaffold=intron.scaffold,
        strand=intron.strand,
        start=intron.start,
        end=intron.end,
        original_id=intron.rec_id,
        original_name=intron.name,
        feature=intron.feature,
        parent=intron.parent,
        prev_exon_last5=prev_exon_last5,
        next_exon_first5=next_exon_first5,
        intron_first20=intron_first20,
        intron_last20=intron_last20,
        intron_full=intron_oriented,
        conv_or_nonconv=c_or_n,
        warnings=warnings,
    )


def deduplicate_introns(candidates):
    """Merge introns with identical coordinates into single entries."""
    merged = {}

    for cand in candidates:
        key = cand.unique_key
        
        if key not in merged:
            merged[key] = cand
        else:
            existing = merged[key]
            existing.support_count += 1
            existing.warnings = list(dict.fromkeys(existing.warnings + cand.warnings))

            if not existing.original_name and cand.original_name:
                existing.original_name = cand.original_name
            if not existing.original_id and cand.original_id:
                existing.original_id = cand.original_id
            if not existing.conv_or_nonconv and cand.conv_or_nonconv:
                existing.conv_or_nonconv = cand.conv_or_nonconv

    return list(merged.values())


def overlaps(a_start, a_end, b_start, b_end):
    """Return True if two intervals overlap."""
    return not (a_end < b_start or b_end < a_start)


def build_exon_index(exons):
    """Group exons by scaffold and strand."""
    exon_index = defaultdict(list)

    for ex in exons:
        exon_index[(ex.scaffold, ex.strand)].append(ex)

    return exon_index


def has_alternative_exon_context(cand, exon_index):
    """Check whether a single intron is in alternative exon context."""
    same_region_exons = exon_index.get((cand.scaffold, cand.strand), [])

    for ex in same_region_exons:
        # Exon exactly ending before intron = normal left exon
        if ex.end == cand.start - 1:
            continue

        # Exon exactly starting after intron = normal right exon
        if ex.start == cand.end + 1:
            continue

        # Exon overlapping intron body suggests alternative exon structure
        if overlaps(ex.start, ex.end, cand.start, cand.end):
            return True

    return False


def assign_overlap_groups(candidates, exon_index):
    """Group overlapping introns on the same scaffold and strand."""
    by_scaffold_strand = defaultdict(list)

    for cand in candidates:
        by_scaffold_strand[(cand.scaffold, cand.strand)].append(cand)

    group_no = 0

    for _, group_candidates in by_scaffold_strand.items():
        ordered = sorted(group_candidates, key=lambda x: (x.start, x.end))

        current_cluster = []
        current_cluster_max_end = None

        for cand in ordered:
            if not current_cluster:
                current_cluster = [cand]
                current_cluster_max_end = cand.end
                continue

            if cand.start <= current_cluster_max_end:
                current_cluster.append(cand)
                current_cluster_max_end = max(current_cluster_max_end, cand.end)
            else:
                group_no += 1
                finalize_overlap_cluster(current_cluster, group_no, exon_index)
                current_cluster = [cand]
                current_cluster_max_end = cand.end

        if current_cluster:
            group_no += 1
            finalize_overlap_cluster(current_cluster, group_no, exon_index)

    return candidates


def finalize_overlap_cluster(cluster, group_no, exon_index):
    """Assign group-based labels and classify introns as constitutive or alternative."""
    ordered = sorted(cluster, key=lambda x: (x.start, x.end))

    for cand in ordered:
        cand.intron_group_no = group_no

    if len(ordered) == 1:
        only = ordered[0]
        only.intron_class = "single"
        only.alt_label = ""

        if only.support_count == 1 and has_alternative_exon_context(only, exon_index):
            only.const_or_alt = "alt"
        else:
            only.const_or_alt = "const"

    else:
        for i, cand in enumerate(ordered, start=1):
            cand.intron_class = "alt"
            cand.alt_label = f"alt{i}"
            cand.const_or_alt = "alt"
            

def write_ready_csv(candidates, out_csv):
    """Write the final intron table to CSV file."""
    fieldnames = [
        "scaffold",
        "id",
        "name",
        "conv_or_nonconv",
        "const_or_alt",
        "prev_exon_last5",
        "next_exon_first5",
        "intron_first20",
        "intron_last20",
        "intron_full",
    ]

    with open(out_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()

        for cand in candidates:
            writer.writerow(
                {
                    "scaffold": cand.scaffold,
                    "id": cand.stable_id,
                    "name": cand.original_name or cand.original_id or "",
                    "conv_or_nonconv": cand.conv_or_nonconv,
                    "const_or_alt": cand.const_or_alt,
                    "prev_exon_last5": cand.prev_exon_last5,
                    "next_exon_first5": cand.next_exon_first5,
                    "intron_first20": cand.intron_first20,
                    "intron_last20": cand.intron_last20,
                    "intron_full": cand.intron_full,
                }
            )



def run_extract(gff_path, fasta_path, out_csv):
    """Run full intron extraction workflow."""
    genome = load_fasta(fasta_path) 
    records = read_gff(gff_path)
    introns = extract_introns(records) 
    exons = extract_exons(records)
    exon_index = build_exon_index(exons)

    candidates = []

    for intron in introns:
        try:
            cand = build_intron_candidate(intron, genome)

            if cand is not None:
                candidates.append(cand)
        except Exception as e:
            print(f"Skipping intron {intron.rec_id or intron.name or intron.feature}: {e}")
            continue

    deduped = deduplicate_introns(candidates)
    annotated = assign_overlap_groups(deduped, exon_index)
    annotated = sorted(annotated, key=lambda x: (x.scaffold, x.strand, x.start, x.end))

    Path(out_csv).parent.mkdir(parents=True, exist_ok=True) 
    write_ready_csv(annotated, out_csv)

    return {
        "num_scaffolds": len(genome),
        "total_gff_records": len(records),
        "num_scaffolds_gff": len({r.scaffold for r in records}),
        "identified_introns": len(introns),
        "unique_introns": len(deduped),
        "saved_rows": len(annotated),
        "out_csv": str(out_csv),
    }


def main():
    ap = argparse.ArgumentParser(description="Extract introns")
    ap.add_argument("--gff", required=True, help="GFF file")
    ap.add_argument("--fasta", required=True, help="FASTA file")
    ap.add_argument("--out", required=True, help="output CSV file")
    args = ap.parse_args() 
    
    run_extract(
        gff_path=args.gff,
        fasta_path=args.fasta,
        out_csv=args.out,
    )


if __name__ == "__main__":
    main()
