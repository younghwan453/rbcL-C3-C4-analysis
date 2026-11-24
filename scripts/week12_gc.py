"""Week 12: 검증된 rbcL CDS에서 염기 조성을 계산하고 표와 그림을 저장한다."""

from collections import Counter
from pathlib import Path

from Bio import SeqIO
from Bio.Data import CodonTable
import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
PHOTOSYNTHESIS = ["C3", "C4"]
TAXONOMY = ["Monocot", "Eudicot"]
METRICS = ["at_pct", "gc_pct", "gc1_pct", "gc2_pct", "gc3_pct"]
GC_POSITIONS = ["gc1_pct", "gc2_pct", "gc3_pct"]
COLORS = {"C3": "#0072B2", "C4": "#D55E00"}
EXPECTED_COUNTS = {("C3", "Monocot"): 8, ("C3", "Eudicot"): 7,
                   ("C4", "Monocot"): 8, ("C4", "Eudicot"): 7}


def nucleotide_metrics(sequence, translation_table):
    """완전 CDS와 번역표에서 원시 개수와 종결 코돈 제외 백분율을 반환한다."""
    sequence = str(sequence).upper()
    if len(sequence) < 6 or len(sequence) % 3 or set(sequence) - set("ACGT"):
        raise ValueError("Expected an unambiguous complete CDS with at least two codons.")
    table = CodonTable.unambiguous_dna_by_id[translation_table]
    if sequence[-3:] not in table.stop_codons:
        raise ValueError("CDS lacks a terminal stop codon in its translation table.")
    sense = sequence[:-3]
    if any(sense[i:i + 3] in table.stop_codons for i in range(0, len(sense), 3)):
        raise ValueError("CDS contains an internal stop codon.")
    # 원시 개수에는 종결 코돈을 포함하고 모든 백분율에서 공통으로 제외한다.
    full_counts = Counter(sequence)
    sense_counts = Counter(sense)
    result = {"cds_length": len(sequence)}
    result.update({f"{base}_count": full_counts[base] for base in "ATGC"})
    result.update({
        "composition_nt_count": len(sense),
        "sense_codon_count": len(sense) // 3,
        "at_pct": 100 * (sense_counts["A"] + sense_counts["T"]) / len(sense),
        "gc_pct": 100 * (sense_counts["G"] + sense_counts["C"]) / len(sense),
    })
    for index, key in enumerate(GC_POSITIONS):
        position_bases = sense[index::3]
        result[key] = 100 * sum(base in "GC" for base in position_bases) / len(position_bases)
    return result


def load_species_metrics(processed_dir):
    """FASTA·메타데이터·검증표를 종 ID로 결합하고 종별 GC 결과를 반환한다."""
    metadata = pd.read_csv(processed_dir / "species_metadata.csv")
    validation = pd.read_csv(processed_dir / "sequence_validation.csv")
    with (processed_dir / "rbcL_nucleotide.fasta").open(encoding="utf-8") as handle:
        fasta = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    if metadata["species_id"].duplicated().any() or validation["species_id"].duplicated().any():
        raise ValueError("Species IDs must be unique in each CSV.")
    if not metadata["validation_status"].eq("validated").all() or not validation["validation_status"].eq("validated").all():
        raise ValueError("The final GC analysis requires validated species records.")
    if set(metadata["species_id"]) != set(validation["species_id"]) or set(metadata["species_id"]) != set(fasta):
        raise ValueError("FASTA, metadata and validation must contain the same species IDs.")
    checks = validation.set_index("species_id")
    rows = []
    for row in metadata.itertuples(index=False):
        check = checks.loc[row.species_id]
        metrics = nucleotide_metrics(fasta[row.species_id].seq, int(check["translation_table"]))
        if (metrics["cds_length"] != row.cds_length or metrics["cds_length"] != check["cds_length"]
                or str(fasta[row.species_id].seq[-3:]).upper() != check["stop_codon"]
                or row.nucleotide_accession != check["nucleotide_accession"]):
            raise ValueError(f"Sequence and record metadata disagree for {row.species_id}.")
        rows.append({
            "species_id": row.species_id, "scientific_name": row.scientific_name,
            "taxonomic_group": row.taxonomic_group, "photosynthesis_type": row.photosynthesis_type,
            "nucleotide_accession": row.nucleotide_accession, **metrics,
        })
    return pd.DataFrame(rows)


def group_summary(species):
    """종별 결과를 받아 세 분류 기준에서 n·평균·표본 표준편차를 반환한다."""
    long = species.melt(
        id_vars=["species_id", "photosynthesis_type", "taxonomic_group"],
        value_vars=METRICS, var_name="metric", value_name="percent",
    )
    tables = []
    for label, groups in (
        ("Photosynthesis type", ["photosynthesis_type"]),
        ("Taxonomic group", ["taxonomic_group"]),
        ("Photosynthesis type and taxonomic group", ["photosynthesis_type", "taxonomic_group"]),
    ):
        summary = (
            long.groupby([*groups, "metric"], sort=False)["percent"]
            .agg(n="count", mean_pct="mean", sd_pct="std").reset_index()
        )
        summary.insert(0, "grouping", label)
        tables.append(summary)
    return pd.concat(tables, ignore_index=True).reindex(
        columns=["grouping", "photosynthesis_type", "taxonomic_group", "metric", "n", "mean_pct", "sd_pct"]
    )


def save_figure(fig, path):
    """Matplotlib 그림을 300 dpi PNG로 저장하고 메모리에서 닫는다."""
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def make_figures(species, figures_dir):
    """종별 지표 DataFrame에서 요청한 네 가지 PNG를 생성한다."""
    sns.set_theme(style="whitegrid", palette=list(COLORS.values()))
    ordered = species.sort_values("gc_pct", ascending=False, kind="stable")
    fig, ax = plt.subplots(figsize=(8.2, max(4.0, len(ordered) * 0.30)), layout="constrained")
    sns.barplot(
        data=ordered, y="scientific_name", x="gc_pct", hue="photosynthesis_type",
        hue_order=PHOTOSYNTHESIS, palette=COLORS, dodge=False, errorbar=None, ax=ax,
    )
    fig.suptitle("rbcL GC Content by Species")
    ax.set(xlabel="GC content (%)", ylabel="Species")
    ax.set_xlim(0, ordered["gc_pct"].max() * 1.13)
    for bars in ax.containers:
        ax.bar_label(bars, fmt="%.2f", padding=2, fontsize=7)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncols=2, frameon=False)
    save_figure(fig, figures_dir / "gc_by_species.png")

    fig, ax = plt.subplots(figsize=(6.2, 4.5), layout="constrained")
    sns.boxplot(
        data=species, x="photosynthesis_type", y="gc_pct", hue="photosynthesis_type",
        order=PHOTOSYNTHESIS, hue_order=PHOTOSYNTHESIS, palette=COLORS,
        legend=False, showfliers=False, ax=ax,
    )
    sns.swarmplot(data=species, x="photosynthesis_type", y="gc_pct",
                  order=PHOTOSYNTHESIS, color="#212121", size=4, ax=ax)
    ax.set(title="rbcL GC Content in C3 and C4 Plants", xlabel="Photosynthesis type", ylabel="GC content (%)")
    save_figure(fig, figures_dir / "gc_c3_vs_c4.png")

    positions = species.melt(id_vars="photosynthesis_type", value_vars=GC_POSITIONS,
                             var_name="position", value_name="position_gc_pct")
    positions["position"] = positions["position"].map(dict(zip(GC_POSITIONS, ["GC1", "GC2", "GC3"])))
    fig, ax = plt.subplots(figsize=(7.2, 4.8), layout="constrained")
    sns.boxplot(
        data=positions, x="position", y="position_gc_pct", hue="photosynthesis_type",
        order=["GC1", "GC2", "GC3"], hue_order=PHOTOSYNTHESIS, palette=COLORS, ax=ax,
    )
    ax.set(title="rbcL GC Content by Codon Position", xlabel="Codon position", ylabel="GC content (%)")
    ax.legend(title="Photosynthesis")
    save_figure(fig, figures_dir / "gc_by_codon_position.png")

    fig, ax = plt.subplots(figsize=(7.4, 5.0), layout="constrained")
    labeled = species.rename(columns={
        "photosynthesis_type": "Photosynthesis type", "taxonomic_group": "Taxonomic group",
    })
    sns.scatterplot(
        data=labeled, x="gc3_pct", y="gc_pct", hue="Photosynthesis type",
        style="Taxonomic group", hue_order=PHOTOSYNTHESIS, style_order=TAXONOMY,
        palette=COLORS, s=66, alpha=0.83, ax=ax,
    )
    ax.set(title="rbcL GC3 vs. Overall GC Content", xlabel="GC3 (%)", ylabel="GC content (%)")
    save_figure(fig, figures_dir / "gc3_vs_gc.png")


def main():
    """검증된 30종을 확인하고 Week 12 결과 CSV와 PNG를 저장한다."""
    species = load_species_metrics(ROOT / "data/processed")
    actual = species.groupby(["photosynthesis_type", "taxonomic_group"]).size().to_dict()
    if actual != EXPECTED_COUNTS:
        raise ValueError(f"Unexpected photosynthesis/taxonomy sample counts: {actual}")
    results_dir = ROOT / "week12/results"
    figures_dir = ROOT / "week12/figures"
    results_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    species.to_csv(results_dir / "species_gc.csv", index=False)
    group_summary(species).to_csv(results_dir / "group_gc_summary.csv", index=False)
    make_figures(species, figures_dir)
    print(f"GC analysis: {len(species)} species, 2 CSV tables, 4 PNG figures")


if __name__ == "__main__":
    main()
