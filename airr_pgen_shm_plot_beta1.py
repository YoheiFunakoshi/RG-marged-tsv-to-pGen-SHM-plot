#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""AIRR outfmt19 TSV -> pGen, SHM, and pGen-SHM plot.

This is the beta_1 implementation for RG merged-read AIRR TSV files.
Input is an IgBLAST AIRR outfmt 19 TSV, not FASTA.
"""

from __future__ import annotations

import argparse
import csv
import io
import math
import os
import re
import sys
import zipfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

try:
    from openpyxl import Workbook
except Exception:  # pragma: no cover - handled by check_setup / runtime error
    Workbook = None


PGEN_EDGES = [1e-5, 1e-10, 1e-15, 1e-20, 1e-25, 1e-30, 1e-35, 1e-40]
PGEN_LABELS = [
    "1e-5~1e-10",
    "1e-10~1e-15",
    "1e-15~1e-20",
    "1e-20~1e-25",
    "1e-25~1e-30",
    "1e-30~1e-35",
    "1e-35~1e-40",
    "1e-40~",
]
SHM_LABELS = [
    "0~2%",
    "2~4%",
    "4~6%",
    "6~8%",
    "8~10%",
    "10~12%",
    "12~14%",
    "14~16%",
    "16~18%",
    "18~20%",
    "20~%",
]
ROW_LEVEL_COLUMNS = [
    "sequence_id",
    "junction_aa",
    "shm",
    "pgen",
    "log10_pgen",
    "junction",
    "v_identity",
    "locus",
    "productive",
    "v_call",
    "j_call",
    "same_xy_count",
    "plot_weight",
]
AA_ALPHABET = set("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
TRUTHY = {"T", "TRUE", "1", "Y", "YES"}


LogFn = Callable[[str], None]


@dataclass
class AnalysisConfig:
    input_path: Path
    output_dir: Path
    sample: str
    cache_path: Path
    use_duplicate_count: bool = False
    min_v_align_len: int = 0
    locus: str = "IGH"
    xlim: tuple[float, float] = (-30.0, -5.0)
    ylim: tuple[float, float] = (0.0, 15.0)
    bw_factor: float = 0.8
    prefix: str | None = None


def log_default(message: str) -> None:
    print(message)


def safe_sample_name(path: Path) -> str:
    name = path.name
    for suffix in (".igblast.airr.tsv", ".airr.tsv", ".tsv", ".txt", ".zip"):
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
            break
    safe = re.sub(r'[<>:"/\\|?*\x00-\x1f]+', "_", name).strip(" ._")
    return safe[:90] or "sample"


def is_truthy(value: str | None) -> bool:
    return str(value or "").strip().upper() in TRUTHY


def is_acgt(seq: str | None) -> bool:
    seq = str(seq or "").strip()
    return bool(seq) and re.fullmatch(r"[ACGTacgt]+", seq) is not None


def aa_ok(aa: str | None) -> bool:
    aa = str(aa or "").strip().upper()
    if not aa:
        return False
    if "*" in aa or "X" in aa:
        return False
    return all(ch in AA_ALPHABET for ch in aa)


def to_int(value: str | None, default: int = 0) -> int:
    try:
        text = str(value or "").strip()
        if not text:
            return default
        return int(float(text))
    except Exception:
        return default


def to_float(value: str | None) -> float | None:
    try:
        text = str(value or "").strip()
        if not text:
            return None
        return float(text)
    except Exception:
        return None


def pgen_bin_label(pgen: float) -> str:
    if pgen <= 0:
        return "1e-40~"
    for i in range(len(PGEN_EDGES) - 1):
        hi = PGEN_EDGES[i]
        lo = PGEN_EDGES[i + 1]
        if lo <= pgen < hi:
            return PGEN_LABELS[i]
    if pgen < PGEN_EDGES[-1]:
        return "1e-40~"
    return "1e-5~1e-10"


def shm_bin_label(shm: float) -> str:
    if shm < 0:
        shm = 0.0
    if shm >= 20.0:
        return "20~%"
    start = int(math.floor(shm / 2.0)) * 2
    start = max(0, min(18, start))
    return f"{start}~{start + 2}%"


def shm_from_identity(v_identity: float) -> float:
    if v_identity <= 1.2:
        shm = (1.0 - v_identity) * 100.0
    else:
        shm = 100.0 - v_identity
    return max(0.0, shm)


def ungapped_len(seq: str | None) -> int:
    return len(re.sub(r"[-.\s]", "", str(seq or "")))


def open_airr_tsv(path: Path) -> Iterable[dict[str, str]]:
    """Yield AIRR rows from a plain TSV or a ZIP containing one TSV."""
    if path.suffix.lower() == ".zip":
        with zipfile.ZipFile(path, "r") as zf:
            tsv_names = [n for n in zf.namelist() if n.lower().endswith((".tsv", ".txt"))]
            if not tsv_names:
                raise ValueError("ZIP does not contain a TSV/TXT file.")
            if len(tsv_names) > 1:
                raise ValueError("ZIP contains multiple TSV/TXT files. Please unzip and choose one.")
            with zf.open(tsv_names[0], "r") as fb:
                text = io.TextIOWrapper(fb, encoding="utf-8-sig", errors="replace", newline="")
                reader = csv.DictReader(text, delimiter="\t")
                if reader.fieldnames is None:
                    raise ValueError("Input TSV has no header.")
                for row in reader:
                    yield row
        return

    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Input TSV has no header.")
        for row in reader:
            yield row


def get_fieldnames(path: Path) -> list[str]:
    if path.suffix.lower() == ".zip":
        with zipfile.ZipFile(path, "r") as zf:
            tsv_names = [n for n in zf.namelist() if n.lower().endswith((".tsv", ".txt"))]
            if len(tsv_names) != 1:
                raise ValueError("ZIP must contain exactly one TSV/TXT file.")
            with zf.open(tsv_names[0], "r") as fb:
                header = fb.readline().decode("utf-8-sig", errors="replace").rstrip("\n\r")
        return header.split("\t")
    with path.open("r", encoding="utf-8-sig", errors="replace") as handle:
        return handle.readline().rstrip("\n\r").split("\t")


def v_alignment_len(row: dict[str, str], fields: set[str]) -> int:
    if "v_sequence_alignment" in fields:
        return ungapped_len(row.get("v_sequence_alignment"))
    if "v_alignment_start" in fields and "v_alignment_end" in fields:
        start = to_int(row.get("v_alignment_start"), default=0)
        end = to_int(row.get("v_alignment_end"), default=0)
        if start > 0 and end >= start:
            return end - start + 1
    return 0


def load_pgen_cache(cache_path: Path) -> dict[str, float]:
    cache: dict[str, float] = {}
    if not cache_path.exists():
        return cache
    with cache_path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            aa = (row.get("junction_aa") or row.get("AA_JUNCTION") or "").strip().upper()
            pgen = to_float(row.get("pgen"))
            if aa and pgen is not None:
                cache[aa] = pgen
    return cache


def compute_pgen_for_aas(aas: Iterable[str], cache_path: Path, log: LogFn) -> dict[str, float]:
    cache = load_pgen_cache(cache_path)
    todo = [aa for aa in sorted(set(aas)) if aa not in cache]
    if not todo:
        log(f"pGen cache hit: {len(cache):,} entries; no new OLGA calculation.")
        return cache

    try:
        import olga
        import olga.generation_probability as generation_probability
        import olga.load_model as load_model
    except Exception as exc:
        raise RuntimeError("OLGA is not installed. Install with: pip install olga") from exc

    olga_dir = Path(olga.__file__).resolve().parent
    model_dir = olga_dir / "default_models" / "human_B_heavy"
    if not model_dir.exists():
        raise RuntimeError(f"OLGA human_B_heavy model directory not found: {model_dir}")

    genomic_data = load_model.GenomicDataVDJ()
    genomic_data.load_igor_genomic_data(
        str(model_dir / "model_params.txt"),
        str(model_dir / "V_gene_CDR3_anchors.csv"),
        str(model_dir / "J_gene_CDR3_anchors.csv"),
    )
    gen_model = load_model.GenerativeModelVDJ()
    gen_model.load_and_process_igor_model(str(model_dir / "model_marginals.txt"))
    pgen_model = generation_probability.GenerationProbabilityVDJ(gen_model, genomic_data)

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    need_header = not cache_path.exists()
    with cache_path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        if need_header:
            writer.writerow(["junction_aa", "pgen"])
        for i, aa in enumerate(todo, start=1):
            try:
                pgen = float(pgen_model.compute_aa_CDR3_pgen(aa, print_warnings=False))
            except TypeError:
                try:
                    pgen = float(pgen_model.compute_aa_CDR3_pgen(aa))
                except Exception:
                    pgen = 0.0
            except Exception:
                pgen = 0.0
            cache[aa] = pgen
            writer.writerow([aa, f"{pgen:.17g}"])
            if i % 100 == 0 or i == len(todo):
                log(f"pGen computed {i:,}/{len(todo):,} new AA sequences.")
    return cache


def read_and_aggregate(config: AnalysisConfig, log: LogFn):
    fields = set(get_fieldnames(config.input_path))
    required = {"productive", "junction", "junction_aa", "v_identity"}
    missing = sorted(required - fields)
    if missing:
        raise ValueError(f"Missing required AIRR column(s): {', '.join(missing)}")

    has_vj = "vj_in_frame" in fields
    has_stop = "stop_codon" in fields
    has_locus = "locus" in fields
    has_dup = "duplicate_count" in fields
    has_vlen = bool({"v_sequence_alignment", "v_alignment_start", "v_alignment_end"} & fields)

    stats = Counter()
    aa_counts = Counter()
    j_to_shm: dict[str, list[float]] = defaultdict(list)
    j_to_aa: dict[str, Counter] = defaultdict(Counter)
    j_to_read_count = Counter()
    j_to_weighted_read_count = Counter()
    j_to_vlen: dict[str, list[int]] = defaultdict(list)
    row_records: list[dict[str, object]] = []

    for row in open_airr_tsv(config.input_path):
        stats["rows_total"] += 1

        locus = (row.get("locus") or "").strip()
        if has_locus and config.locus:
            if locus and locus != config.locus:
                stats["drop_non_IGH_locus"] += 1
                continue
            if not locus:
                stats["kept_empty_locus"] += 1

        if not is_truthy(row.get("productive")):
            stats["drop_nonproductive"] += 1
            continue
        if has_vj and not is_truthy(row.get("vj_in_frame")):
            stats["drop_vj_outframe"] += 1
            continue
        if has_stop and is_truthy(row.get("stop_codon")):
            stats["drop_stopcodon"] += 1
            continue

        junction = (row.get("junction") or "").strip().upper()
        if not is_acgt(junction):
            stats["drop_bad_junction_nt"] += 1
            continue

        aa = (row.get("junction_aa") or "").strip().upper()
        if not aa_ok(aa):
            stats["drop_bad_junction_aa"] += 1
            continue

        v_identity = to_float(row.get("v_identity"))
        if v_identity is None:
            stats["drop_missing_v_identity"] += 1
            continue

        vlen = v_alignment_len(row, fields) if has_vlen else 0
        if config.min_v_align_len > 0:
            if not has_vlen:
                raise ValueError("min V alignment length was set, but no V alignment length column is available.")
            if vlen < config.min_v_align_len:
                stats["drop_short_v_alignment"] += 1
                continue

        shm = shm_from_identity(v_identity)
        weight = 1
        if config.use_duplicate_count and has_dup:
            weight = max(1, to_int(row.get("duplicate_count"), default=1))

        sequence_id = (row.get("sequence_id") or "").strip()
        if not sequence_id:
            sequence_id = f"row_{stats['rows_total']}"
            stats["missing_sequence_id_replaced"] += 1

        aa_counts[aa] += weight
        j_to_shm[junction].append(shm)
        j_to_aa[junction][aa] += 1
        j_to_read_count[junction] += 1
        j_to_weighted_read_count[junction] += weight
        if vlen:
            j_to_vlen[junction].append(vlen)
        stats["kept_weighted_reads"] += weight
        row_records.append(
            {
                "sequence_id": sequence_id,
                "junction_aa": aa,
                "shm": shm,
                "pgen": 0.0,
                "log10_pgen": math.nan,
                "junction": junction,
                "v_identity": v_identity,
                "locus": locus,
                "productive": (row.get("productive") or "").strip(),
                "v_call": (row.get("v_call") or "").strip(),
                "j_call": (row.get("j_call") or "").strip(),
                "same_xy_count": 0,
                "plot_weight": weight,
            }
        )

    stats["kept_reads"] = sum(j_to_read_count.values())
    stats["kept_unique_junction_nt"] = len(j_to_read_count)
    stats["kept_unique_junction_aa"] = len(aa_counts)
    stats["has_vj_in_frame_column"] = int(has_vj)
    stats["has_stop_codon_column"] = int(has_stop)
    stats["has_locus_column"] = int(has_locus)
    stats["has_duplicate_count_column"] = int(has_dup)
    stats["has_v_alignment_length_data"] = int(has_vlen)

    if stats["kept_reads"] == 0:
        raise ValueError("No rows left after beta_1 filtering. Check AIRR columns and filters.")
    log(f"Rows kept after beta_1 filtering: {stats['kept_reads']:,}")
    log(f"Unique junction(nt) points: {stats['kept_unique_junction_nt']:,}")
    log(f"Unique junction_aa for pGen: {stats['kept_unique_junction_aa']:,}")
    return aa_counts, j_to_shm, j_to_aa, j_to_read_count, j_to_weighted_read_count, j_to_vlen, row_records, stats


def write_qc_summary(stats: Counter, out_tsv: Path) -> None:
    with out_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        for key in sorted(stats):
            writer.writerow([key, stats[key]])


def write_pgen_bins(aa_counts: Counter, aa_to_pgen: dict[str, float], out_tsv: Path) -> tuple[list[float], list[float]]:
    unique_counts = Counter()
    weighted_counts = Counter()

    for aa, weight in aa_counts.items():
        label = pgen_bin_label(float(aa_to_pgen.get(aa, 0.0)))
        unique_counts[label] += 1
        weighted_counts[label] += int(weight)

    total_unique = sum(unique_counts.values()) or 1
    total_weighted = sum(weighted_counts.values()) or 1
    frac_unique = [unique_counts.get(label, 0) / total_unique for label in PGEN_LABELS]
    frac_weighted = [weighted_counts.get(label, 0) / total_weighted for label in PGEN_LABELS]

    with out_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["bin", "count_uniqueAA", "fraction_uniqueAA", "count_weightedReads", "fraction_weightedReads"])
        for label, fu, fw in zip(PGEN_LABELS, frac_unique, frac_weighted):
            writer.writerow([label, unique_counts.get(label, 0), f"{fu:.12g}", weighted_counts.get(label, 0), f"{fw:.12g}"])
    return frac_unique, frac_weighted


def plot_barh(labels: list[str], fractions: list[float], ylabel: str, xlabel: str, title: str, out_png: Path) -> None:
    y = np.arange(len(labels))[::-1]
    plt.figure(figsize=(8.4, 5.0))
    plt.barh(y, fractions[::-1], color="#2b7fb8")
    plt.yticks(y, labels[::-1])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    xmax = max(fractions) * 1.12 if fractions and max(fractions) > 0 else 1.0
    plt.xlim(0, xmax)
    plt.grid(True, axis="x", alpha=0.25)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def write_points(
    junctions: Iterable[str],
    j_to_shm: dict[str, list[float]],
    j_to_aa: dict[str, Counter],
    j_to_read_count: Counter,
    j_to_weighted_read_count: Counter,
    j_to_vlen: dict[str, list[int]],
    aa_to_pgen: dict[str, float],
    out_tsv: Path,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    total_reads = sum(j_to_read_count.values()) or 1
    total_weighted_reads = sum(j_to_weighted_read_count.values()) or 1
    with out_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "junction_nt",
            "junction_aa",
            "read_count",
            "shm_median",
            "pgen",
            "log10_pgen",
            "v_seq_len_median",
            "aa_candidate_count",
            "read_fraction",
            "weighted_read_count",
            "weighted_read_fraction",
        ])
        for junction in junctions:
            shm_values = j_to_shm.get(junction, [])
            if not shm_values:
                continue
            aa_counter = j_to_aa.get(junction, Counter())
            representative_aa = aa_counter.most_common(1)[0][0] if aa_counter else ""
            pgen = float(aa_to_pgen.get(representative_aa, 0.0))
            log10_pgen = math.log10(pgen) if pgen > 0 else math.nan
            vlen_values = j_to_vlen.get(junction, [])
            vlen_median = int(np.median(vlen_values)) if vlen_values else 0
            read_count = int(j_to_read_count.get(junction, 0))
            weighted_read_count = int(j_to_weighted_read_count.get(junction, read_count))
            row = {
                "junction_nt": junction,
                "junction_aa": representative_aa,
                "read_count": read_count,
                "read_fraction": read_count / total_reads,
                "weighted_read_count": weighted_read_count,
                "weighted_read_fraction": weighted_read_count / total_weighted_reads,
                "shm_median": float(np.median(shm_values)),
                "pgen": pgen,
                "log10_pgen": log10_pgen,
                "v_seq_len_median": vlen_median,
                "aa_candidate_count": len(aa_counter),
            }
            rows.append(row)
            writer.writerow([
                row["junction_nt"],
                row["junction_aa"],
                row["read_count"],
                f"{row['shm_median']:.12g}",
                f"{pgen:.17g}",
                "" if math.isnan(log10_pgen) else f"{log10_pgen:.12g}",
                row["v_seq_len_median"],
                row["aa_candidate_count"],
                f"{row['read_fraction']:.12g}",
                row["weighted_read_count"],
                f"{row['weighted_read_fraction']:.12g}",
            ])
    return rows


def add_pgen_to_row_records(row_records: list[dict[str, object]], aa_to_pgen: dict[str, float]) -> None:
    xy_counts = Counter()
    for row in row_records:
        aa = str(row["junction_aa"])
        pgen = float(aa_to_pgen.get(aa, 0.0))
        log10_pgen = math.log10(pgen) if pgen > 0 else math.nan
        row["pgen"] = pgen
        row["log10_pgen"] = log10_pgen
        if pgen > 0 and not math.isnan(log10_pgen):
            xy_counts[(log10_pgen, float(row["shm"]))] += 1

    for row in row_records:
        log10_pgen = float(row["log10_pgen"])
        if float(row["pgen"]) > 0 and not math.isnan(log10_pgen):
            row["same_xy_count"] = int(xy_counts[(log10_pgen, float(row["shm"]))])
        else:
            row["same_xy_count"] = 0


def clean_excel_value(value: object) -> object:
    if isinstance(value, float) and math.isnan(value):
        return None
    return value


def write_row_level_excel(row_records: list[dict[str, object]], out_xlsx: Path) -> None:
    if Workbook is None:
        raise RuntimeError("openpyxl is required to write Excel output. Install it with: pip install openpyxl")
    wb = Workbook(write_only=True)
    ws = wb.create_sheet("accepted_rows")
    ws.append(ROW_LEVEL_COLUMNS)
    for row in row_records:
        ws.append([clean_excel_value(row.get(column, "")) for column in ROW_LEVEL_COLUMNS])
    wb.save(out_xlsx)


def write_row_level_tsv(row_records: list[dict[str, object]], out_tsv: Path) -> None:
    with out_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(ROW_LEVEL_COLUMNS)
        for row in row_records:
            values = []
            for column in ROW_LEVEL_COLUMNS:
                value = clean_excel_value(row.get(column, ""))
                if value is None:
                    value = ""
                values.append(value)
            writer.writerow(values)


def compact_xy_points(
    rows: list[dict[str, object]],
    config: AnalysisConfig,
    weight_key: str | None = None,
) -> list[dict[str, float]]:
    compacted: dict[tuple[float, float], dict[str, float]] = {}
    for row in rows:
        pgen = float(row["pgen"])
        log10_pgen = float(row["log10_pgen"])
        shm = float(row["shm"])
        if (
            pgen <= 0
            or math.isnan(log10_pgen)
            or not (config.xlim[0] <= log10_pgen <= config.xlim[1])
            or not (config.ylim[0] <= shm <= config.ylim[1])
        ):
            continue
        key = (log10_pgen, shm)
        if key not in compacted:
            compacted[key] = {
                "log10_pgen": log10_pgen,
                "shm": shm,
                "same_xy_count": 0.0,
                "plot_weight": 0.0,
            }
        compacted[key]["same_xy_count"] += 1.0
        if weight_key:
            compacted[key]["plot_weight"] += max(0.0, float(row.get(weight_key, 1.0)))
        else:
            compacted[key]["plot_weight"] += 1.0
    return list(compacted.values())


def write_shm_hist(
    points: list[dict[str, object]],
    out_tsv: Path,
    weight_key: str | None = None,
    shm_key: str = "shm_median",
) -> list[float]:
    counts = Counter()
    for row in points:
        weight = float(row.get(weight_key, 1.0)) if weight_key else 1.0
        counts[shm_bin_label(float(row[shm_key]))] += weight
    total = sum(counts.values()) or 1
    fractions = [counts.get(label, 0) / total for label in SHM_LABELS]
    with out_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        if weight_key:
            writer.writerow(["bin", "count_weightedReads", "fraction_weightedReads"])
        else:
            writer.writerow(["bin", "count", "fraction"])
        for label, fraction in zip(SHM_LABELS, fractions):
            count_value = counts.get(label, 0)
            if float(count_value).is_integer():
                count_value = int(count_value)
            writer.writerow([label, count_value, f"{fraction:.12g}"])
    return fractions


def plot_kde(
    points: list[dict[str, object]],
    config: AnalysisConfig,
    out_png: Path,
    log: LogFn,
    weight_key: str | None = None,
) -> int:
    filtered = [
        row
        for row in points
        if float(row["pgen"]) > 0
        and not math.isnan(float(row["log10_pgen"]))
        and config.xlim[0] <= float(row["log10_pgen"]) <= config.xlim[1]
        and config.ylim[0] <= float(row["shm_median"]) <= config.ylim[1]
    ]
    if len(filtered) < 5:
        log("Too few valid points for KDE; writing scatter fallback instead.")
        plot_scatter(filtered, config, out_png)
        return len(filtered)

    xs = np.array([float(row["log10_pgen"]) for row in filtered], dtype=float)
    ys = np.array([float(row["shm_median"]) for row in filtered], dtype=float)
    weights = None
    if weight_key:
        weights = np.array([max(0.0, float(row.get(weight_key, 0.0))) for row in filtered], dtype=float)
        if weights.sum() <= 0:
            weights = None
    try:
        kde_args = {
            "bw_method": lambda obj: obj.scotts_factor() * config.bw_factor,
        }
        if weights is not None:
            kde_args["weights"] = weights
        kde = gaussian_kde(np.vstack([xs, ys]), **kde_args)
        xi = np.linspace(config.xlim[0], config.xlim[1], 250)
        yi = np.linspace(config.ylim[0], config.ylim[1], 250)
        x_grid, y_grid = np.meshgrid(xi, yi)
        density = kde(np.vstack([x_grid.ravel(), y_grid.ravel()])).reshape(x_grid.shape)
    except TypeError as exc:
        if weights is not None:
            log(f"Weighted KDE is not supported by this SciPy ({exc}); writing unweighted KDE instead.")
            return plot_kde(points, config, out_png, log, weight_key=None)
        log(f"KDE failed ({exc}); writing scatter fallback instead.")
        plot_scatter(filtered, config, out_png)
        return len(filtered)
    except Exception as exc:
        log(f"KDE failed ({exc}); writing scatter fallback instead.")
        plot_scatter(filtered, config, out_png)
        return len(filtered)

    plt.figure(figsize=(6.2, 6.0))
    plt.contourf(x_grid, y_grid, density, levels=12, cmap="YlOrRd")
    plt.xlabel("pGen (log10)")
    plt.ylabel("%Mutation")
    suffix = "weighted KDE" if weight_key else "KDE"
    plt.title(f"{config.sample} pGen-SHM {suffix}")
    plt.xlim(config.xlim)
    plt.ylim(config.ylim)
    plt.grid(True, alpha=0.35)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()
    return len(filtered)


def plot_weighted_scatter(points: list[dict[str, object]], config: AnalysisConfig, out_png: Path) -> int:
    filtered = [
        row
        for row in points
        if float(row["pgen"]) > 0
        and not math.isnan(float(row["log10_pgen"]))
        and config.xlim[0] <= float(row["log10_pgen"]) <= config.xlim[1]
        and config.ylim[0] <= float(row["shm_median"]) <= config.ylim[1]
    ]
    xs = np.array([float(row["log10_pgen"]) for row in filtered], dtype=float)
    ys = np.array([float(row["shm_median"]) for row in filtered], dtype=float)
    weights = np.array([max(1.0, float(row.get("weighted_read_count", row.get("read_count", 1)))) for row in filtered], dtype=float)
    if len(filtered) == 0:
        sizes = np.array([], dtype=float)
    else:
        log_weights = np.log10(weights + 1.0)
        if float(log_weights.max()) > float(log_weights.min()):
            sizes = 12.0 + 80.0 * (log_weights - log_weights.min()) / (log_weights.max() - log_weights.min())
        else:
            sizes = np.full_like(log_weights, 28.0)

    plt.figure(figsize=(6.2, 6.0))
    plt.scatter(xs, ys, s=sizes, alpha=0.38, color="#2b7fb8", edgecolors="none")
    plt.xlabel("pGen (log10)")
    plt.ylabel("%Mutation")
    plt.title(f"{config.sample} pGen-SHM scatter weighted by reads")
    plt.xlim(config.xlim)
    plt.ylim(config.ylim)
    plt.grid(True, alpha=0.35)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()
    return len(filtered)


def plot_xy_kde(
    xy_points: list[dict[str, float]],
    config: AnalysisConfig,
    out_png: Path,
    log: LogFn,
    title_suffix: str,
) -> int:
    if len(xy_points) < 5:
        log("Too few valid row-level points for KDE; writing scatter fallback instead.")
        plot_xy_scatter(xy_points, config, out_png, title_suffix)
        return int(sum(point.get("plot_weight", 1.0) for point in xy_points))

    xs = np.array([float(point["log10_pgen"]) for point in xy_points], dtype=float)
    ys = np.array([float(point["shm"]) for point in xy_points], dtype=float)
    weights = np.array([max(0.0, float(point.get("plot_weight", 1.0))) for point in xy_points], dtype=float)
    if weights.sum() <= 0:
        weights = None
    try:
        kde_args = {
            "bw_method": lambda obj: obj.scotts_factor() * config.bw_factor,
        }
        if weights is not None:
            kde_args["weights"] = weights
        kde = gaussian_kde(np.vstack([xs, ys]), **kde_args)
        xi = np.linspace(config.xlim[0], config.xlim[1], 250)
        yi = np.linspace(config.ylim[0], config.ylim[1], 250)
        x_grid, y_grid = np.meshgrid(xi, yi)
        density = kde(np.vstack([x_grid.ravel(), y_grid.ravel()])).reshape(x_grid.shape)
    except TypeError as exc:
        log(f"Weighted KDE is not supported by this SciPy ({exc}); writing scatter fallback instead.")
        plot_xy_scatter(xy_points, config, out_png, title_suffix)
        return int(sum(point.get("plot_weight", 1.0) for point in xy_points))
    except Exception as exc:
        log(f"KDE failed ({exc}); writing scatter fallback instead.")
        plot_xy_scatter(xy_points, config, out_png, title_suffix)
        return int(sum(point.get("plot_weight", 1.0) for point in xy_points))

    plt.figure(figsize=(6.2, 6.0))
    plt.contourf(x_grid, y_grid, density, levels=12, cmap="YlOrRd")
    plt.xlabel("pGen (log10)")
    plt.ylabel("%Mutation")
    plt.title(f"{config.sample} pGen-SHM {title_suffix}")
    plt.xlim(config.xlim)
    plt.ylim(config.ylim)
    plt.grid(True, alpha=0.35)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()
    return int(sum(point.get("plot_weight", 1.0) for point in xy_points))


def plot_xy_scatter(
    xy_points: list[dict[str, float]],
    config: AnalysisConfig,
    out_png: Path,
    title_suffix: str,
) -> int:
    xs = np.array([float(point["log10_pgen"]) for point in xy_points], dtype=float)
    ys = np.array([float(point["shm"]) for point in xy_points], dtype=float)
    weights = np.array([max(1.0, float(point.get("plot_weight", 1.0))) for point in xy_points], dtype=float)
    if len(xy_points) == 0:
        sizes = np.array([], dtype=float)
    else:
        log_weights = np.log10(weights + 1.0)
        if float(log_weights.max()) > float(log_weights.min()):
            sizes = 10.0 + 75.0 * (log_weights - log_weights.min()) / (log_weights.max() - log_weights.min())
        else:
            sizes = np.full_like(log_weights, 24.0)

    plt.figure(figsize=(6.2, 6.0))
    plt.scatter(xs, ys, s=sizes, alpha=0.35, color="#2b7fb8", edgecolors="none")
    plt.xlabel("pGen (log10)")
    plt.ylabel("%Mutation")
    plt.title(f"{config.sample} pGen-SHM {title_suffix}")
    plt.xlim(config.xlim)
    plt.ylim(config.ylim)
    plt.grid(True, alpha=0.35)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()
    return int(sum(point.get("plot_weight", 1.0) for point in xy_points))


def plot_scatter(points: list[dict[str, object]], config: AnalysisConfig, out_png: Path) -> None:
    xs = [float(row["log10_pgen"]) for row in points if float(row["pgen"]) > 0]
    ys = [float(row["shm_median"]) for row in points if float(row["pgen"]) > 0]
    plt.figure(figsize=(6.2, 6.0))
    plt.scatter(xs, ys, s=10, alpha=0.4, color="#2b7fb8", edgecolors="none")
    plt.xlabel("pGen (log10)")
    plt.ylabel("%Mutation")
    plt.title(f"{config.sample} pGen-SHM")
    plt.xlim(config.xlim)
    plt.ylim(config.ylim)
    plt.grid(True, alpha=0.35)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def write_run_log(lines: list[str], out_path: Path) -> None:
    with out_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(lines))
        handle.write("\n")


def run_analysis(config: AnalysisConfig, log: LogFn = log_default) -> dict[str, Path]:
    config.output_dir.mkdir(parents=True, exist_ok=True)
    prefix = config.prefix or safe_sample_name(Path(config.sample))
    outputs = {
        "qc_summary": config.output_dir / f"{prefix}_qc_summary.tsv",
        "pgen_bins": config.output_dir / f"{prefix}_pgen_bins.tsv",
        "pgen_bins_unique_png": config.output_dir / f"{prefix}_pgen_bins_unique.png",
        "pgen_bins_weighted_png": config.output_dir / f"{prefix}_pgen_bins_weighted.png",
        "shm_hist": config.output_dir / f"{prefix}_shm_hist.tsv",
        "shm_hist_png": config.output_dir / f"{prefix}_shm_hist.png",
        "shm_hist_weighted": config.output_dir / f"{prefix}_shm_hist_weighted.tsv",
        "shm_hist_weighted_png": config.output_dir / f"{prefix}_shm_hist_weighted.png",
        "rows_xlsx": config.output_dir / f"{prefix}_pgen_shm_rows.xlsx",
        "points": config.output_dir / f"{prefix}_pgen_shm_points.tsv",
        "kde_png": config.output_dir / f"{prefix}_pgen_shm_kde_unweighted.png",
        "kde_weighted_png": config.output_dir / f"{prefix}_pgen_shm_kde_weighted.png",
        "scatter_reads_png": config.output_dir / f"{prefix}_pgen_shm_scatter_reads.png",
        "scatter_weighted_png": config.output_dir / f"{prefix}_pgen_shm_scatter_weighted.png",
        "run_log": config.output_dir / f"{prefix}_run_log.txt",
    }
    run_log_lines: list[str] = []

    def log_both(message: str) -> None:
        run_log_lines.append(message)
        log(message)

    log_both(f"Input AIRR TSV: {config.input_path}")
    log_both(f"Output folder: {config.output_dir}")
    log_both(f"Sample: {config.sample}")
    log_both(f"pGen cache: {config.cache_path}")
    log_both(f"Locus policy: keep {config.locus}; empty locus is kept with a QC count.")

    aa_counts, j_to_shm, j_to_aa, j_to_read_count, j_to_weighted_read_count, j_to_vlen, row_records, stats = read_and_aggregate(
        config,
        log_both,
    )
    write_qc_summary(stats, outputs["qc_summary"])
    log_both(f"Saved QC summary: {outputs['qc_summary']}")

    aa_all = sorted(aa_counts.keys())
    aa_to_pgen = compute_pgen_for_aas(aa_all, config.cache_path, log_both)
    add_pgen_to_row_records(row_records, aa_to_pgen)

    frac_unique, frac_weighted = write_pgen_bins(aa_counts, aa_to_pgen, outputs["pgen_bins"])
    plot_barh(PGEN_LABELS, frac_unique, "pGen", "Frequency", f"{config.sample} pGen bins (unique AA)", outputs["pgen_bins_unique_png"])
    plot_barh(PGEN_LABELS, frac_weighted, "pGen", "Frequency", f"{config.sample} pGen bins (weighted by reads)", outputs["pgen_bins_weighted_png"])
    log_both(f"Saved pGen bins: {outputs['pgen_bins']}")
    log_both(f"Saved pGen plots: {outputs['pgen_bins_unique_png']} ; {outputs['pgen_bins_weighted_png']}")

    write_row_level_excel(row_records, outputs["rows_xlsx"])
    write_row_level_tsv(row_records, outputs["points"])
    zero_points = sum(1 for row in row_records if float(row["pgen"]) <= 0)
    duplicate_xy_rows = sum(1 for row in row_records if int(row.get("same_xy_count", 0)) > 1)
    unique_xy_points = len({(float(row["log10_pgen"]), float(row["shm"])) for row in row_records if float(row["pgen"]) > 0})
    log_both(f"Saved row-level pGen-SHM Excel: {outputs['rows_xlsx']}")
    log_both(f"Saved row-level pGen-SHM points TSV: {outputs['points']}")
    log_both(f"pGen=0 rows retained in Excel/points TSV and excluded from plots: {zero_points:,}")
    log_both(f"Row-level plot rows: {len(row_records):,}; unique plotted x-y coordinates: {unique_xy_points:,}; rows sharing an x-y coordinate: {duplicate_xy_rows:,}")

    shm_frac = write_shm_hist(row_records, outputs["shm_hist"], shm_key="shm")
    plot_barh(SHM_LABELS, shm_frac, "%Mutation", "Frequency", f"{config.sample} SHM histogram", outputs["shm_hist_png"])
    log_both(f"Saved SHM histogram: {outputs['shm_hist']} ; {outputs['shm_hist_png']}")

    shm_weighted_frac = write_shm_hist(row_records, outputs["shm_hist_weighted"], weight_key="plot_weight", shm_key="shm")
    plot_barh(
        SHM_LABELS,
        shm_weighted_frac,
        "%Mutation",
        "Frequency",
        f"{config.sample} SHM histogram (weighted by reads)",
        outputs["shm_hist_weighted_png"],
    )
    log_both(f"Saved weighted SHM histogram: {outputs['shm_hist_weighted']} ; {outputs['shm_hist_weighted_png']}")

    xy_points = compact_xy_points(row_records, config)
    weighted_xy_points = compact_xy_points(row_records, config, weight_key="plot_weight")
    scatter_points = plot_xy_scatter(xy_points, config, outputs["scatter_reads_png"], "row-level scatter")
    log_both(f"Saved row-level pGen-SHM scatter plot: {outputs['scatter_reads_png']}")
    log_both(f"Row-level scatter plotted rows: {scatter_points:,}; plotted x-y coordinates: {len(xy_points):,}")

    kde_points = plot_xy_kde(xy_points, config, outputs["kde_png"], log_both, "row-level KDE")
    log_both(f"Saved row-level pGen-SHM KDE plot: {outputs['kde_png']}")
    log_both(f"Row-level KDE plotted rows: {kde_points:,}")

    weighted_kde_points = plot_xy_kde(weighted_xy_points, config, outputs["kde_weighted_png"], log_both, "row-level weighted KDE")
    log_both(f"Saved row-level weighted pGen-SHM KDE plot: {outputs['kde_weighted_png']}")
    log_both(f"Row-level weighted KDE plotted weight: {weighted_kde_points:,}")

    scatter_points = plot_xy_scatter(weighted_xy_points, config, outputs["scatter_weighted_png"], "row-level weighted scatter")
    log_both(f"Saved weighted pGen-SHM scatter plot: {outputs['scatter_weighted_png']}")
    log_both(f"Weighted scatter plotted weight: {scatter_points:,}; plotted x-y coordinates: {len(weighted_xy_points):,}")

    write_run_log(run_log_lines, outputs["run_log"])
    log(f"Saved run log: {outputs['run_log']}")
    return outputs


def parse_range(value: str, label: str) -> tuple[float, float]:
    try:
        left, right = [float(x.strip()) for x in value.split(",", 1)]
    except Exception as exc:
        raise argparse.ArgumentTypeError(f"{label} must be two comma-separated numbers, e.g. -30,-5") from exc
    if right <= left:
        raise argparse.ArgumentTypeError(f"{label}: max must be greater than min.")
    return left, right


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="AIRR outfmt19 TSV -> pGen / SHM / pGen-SHM plot (beta_1)")
    parser.add_argument("--input", required=True, help="Input AIRR outfmt19 TSV, TXT, or ZIP containing one TSV.")
    parser.add_argument("--outdir", default="", help="Output folder. Default: same folder as input.")
    parser.add_argument("--sample", default="", help="Sample name for plot titles and output prefix.")
    parser.add_argument("--pgen-cache", default="", help="pGen cache TSV. Default: <outdir>/pgen_cache.tsv")
    parser.add_argument("--use-duplicate-count", action="store_true", help="Use duplicate_count as pGen weighted count if present.")
    parser.add_argument("--min-v-align-len", type=int, default=0, help="Optional V alignment length filter. Default: 0 disabled.")
    parser.add_argument("--locus", default="IGH", help="Expected locus. Default: IGH. Empty disables locus filtering.")
    parser.add_argument("--xlim", default="-30,-5", help="KDE x-axis log10 pGen range, e.g. -30,-5.")
    parser.add_argument("--ylim", default="0,15", help="KDE y-axis SHM range, e.g. 0,15.")
    parser.add_argument("--bw-factor", type=float, default=0.8, help="KDE bandwidth multiplier. Default: 0.8.")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    input_path = Path(args.input).resolve()
    if not input_path.exists():
        parser.error(f"input not found: {input_path}")
    output_dir = Path(args.outdir).resolve() if args.outdir else input_path.parent
    sample = args.sample.strip() or safe_sample_name(input_path)
    cache_path = Path(args.pgen_cache).resolve() if args.pgen_cache else output_dir / "pgen_cache.tsv"
    config = AnalysisConfig(
        input_path=input_path,
        output_dir=output_dir,
        sample=sample,
        cache_path=cache_path,
        use_duplicate_count=bool(args.use_duplicate_count),
        min_v_align_len=max(0, int(args.min_v_align_len)),
        locus=args.locus.strip(),
        xlim=parse_range(args.xlim, "xlim"),
        ylim=parse_range(args.ylim, "ylim"),
        bw_factor=float(args.bw_factor),
        prefix=safe_sample_name(Path(sample)),
    )
    run_analysis(config)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
