"""Week 13: 검증된 rbcL CDS와 단백질의 종별 조성 및 그룹 요약을 계산한다."""

from collections import Counter
from pathlib import Path

from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd

from scripts.week12_gc import COLORS, EXPECTED_COUNTS, PHOTOSYNTHESIS, plt, sns


ROOT = Path(__file__).resolve().parents[1]
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
GROUPS = ["photosynthesis_type", "taxonomic_group"]
PROPERTIES = {
    "protein_length": "Protein length (aa)",
    "molecular_weight_da": "Molecular weight (Da)",
    "isoelectric_point": "Isoelectric point (pI)",
    "gravy": "GRAVY",
    "aromaticity": "Aromaticity (fraction)",
}


def sequence_composition(cds, protein, translation_table):
    """완전 CDS와 단백질을 대조한 뒤 종 한 개의 조성과 특성을 계산한다."""
    cds = cds.upper()
    protein = protein.upper()
    if set(cds) - set("ACGT"):
        raise ValueError("CDS contains ambiguous bases.")
    if not protein or set(protein) - set(AMINO_ACIDS):
        raise ValueError("Protein must contain only standard amino acids.")
    if str(cds.translate(table=translation_table, cds=True)) != protein:
        raise ValueError("CDS translation differs from protein FASTA.")

    analysis = ProteinAnalysis(protein)
    amino_acids = [
        {"amino_acid": aa, "count": count, "percent": 100 * count / len(protein)}
        for aa, count in sorted(analysis.count_amino_acids().items())
    ]
    properties = {
        "protein_length": len(protein),
        "molecular_weight_da": analysis.molecular_weight(),
        "isoelectric_point": analysis.isoelectric_point(),
        "gravy": analysis.gravy(),
        "aromaticity": analysis.aromaticity(),
    }

    # 종결 코돈을 제외하고, 시작 위치를 포함한 실제 DNA 코돈을 센다.
    forward = CodonTable.unambiguous_dna_by_id[translation_table].forward_table
    counts = Counter(str(cds[i:i + 3]) for i in range(0, len(cds) - 3, 3))
    family_size = Counter(forward.values())
    family_count = Counter()
    for codon, count in counts.items():
        family_count[forward[codon]] += count

    codons = []
    for codon, aa in sorted(forward.items()):
        codons.append({
            "codon": codon,
            "amino_acid": aa,
            "count": counts[codon],
            "frequency_pct": 100 * counts[codon] / len(protein),
            "synonymous_family_size": family_size[aa],
            "synonymous_family_count": family_count[aa],
            "rscu": counts[codon] * family_size[aa] / family_count[aa] if family_count[aa] else None,
        })
    return amino_acids, properties, codons


def load_composition(processed):
    """종 ID로 메타데이터, 검증표 및 두 FASTA를 일대일 연결한다."""
    metadata, validation = [
        pd.read_csv(processed / name, keep_default_na=False).set_index("species_id")
        for name in ("species_metadata.csv", "sequence_validation.csv")
    ]
    with (processed / "rbcL_nucleotide.fasta").open(encoding="utf-8") as handle:
        cds_records = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    with (processed / "rbcL_protein.fasta").open(encoding="utf-8") as handle:
        protein_records = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    if not metadata.index.is_unique or not validation.index.is_unique:
        raise ValueError("Species IDs must be unique in each CSV.")
    ids = set(metadata.index)
    if "" in ids or ids != set(validation.index) or ids != set(cds_records) or ids != set(protein_records):
        raise ValueError("The four inputs must contain the same nonempty species IDs.")
    if metadata.groupby(GROUPS).size().to_dict() != EXPECTED_COUNTS:
        raise ValueError("Unexpected species counts in the four study groups.")

    amino_acids, properties, codons = [], [], []
    for species_id, record in metadata.iterrows():
        check = validation.loc[species_id]
        cds = cds_records[species_id].seq.upper()
        protein = str(protein_records[species_id].seq).upper()
        shared = ["scientific_name", "validation_status", "nucleotide_accession", "protein_accession",
                  "cds_length", "protein_length"]
        if (not record[shared].equals(check[shared]) or record["validation_status"] != "validated"
                or str(check["annotation_translation_match"]).lower() != "true"
                or len(cds) != record["cds_length"] or len(protein) != record["protein_length"]
                or str(cds[-3:]) != check["stop_codon"]):
            raise ValueError(f"Sequence or validation metadata disagree for {species_id}.")

        table = int(check["translation_table"])
        try:
            species_aa, species_properties, species_codons = sequence_composition(cds, protein, table)
        except (ValueError, KeyError, CodonTable.TranslationError) as error:
            raise ValueError(f"{species_id}: {error}") from error
        label = {"species_id": species_id, **record[["scientific_name", *GROUPS]].to_dict()}
        amino_acids.extend({**label, **row} for row in species_aa)
        properties.append({**label, **species_properties})
        codons.extend({**label, "translation_table": table, **row} for row in species_codons)
    return pd.DataFrame(amino_acids), pd.DataFrame(properties), pd.DataFrame(codons)


def group_summary(rows):
    """C3/C4 전체 및 분류군별로 종 단위 평균과 표본 표준편차를 계산한다."""
    summaries = []
    for group_keys in ([GROUPS[0]], GROUPS):
        grouped = rows.groupby([*group_keys, "measure", "feature"], sort=True)["value"]
        summary = grouped.agg(n_species="size", n_available="count", mean="mean", sd="std").reset_index()
        if len(group_keys) == 1:
            summary["taxonomic_group"] = "All"
        summaries.append(summary)
    columns = [*GROUPS, "measure", "feature", "n_species", "n_available", "mean", "sd"]
    return pd.concat(summaries, ignore_index=True)[columns]


def save_figure(fig, path):
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_results(amino_acids, properties, codons, summary, output):
    """종별 히트맵과 광합성 유형·분류군별 요약 그림을 저장한다."""
    species_order = properties.sort_values([*GROUPS, "scientific_name"])
    species_ids = species_order["species_id"]
    species_labels = [
        f"{row.photosynthesis_type} {row.taxonomic_group[0]} | {row.scientific_name}"
        for row in species_order.itertuples()
    ]
    specifications = [
        (amino_acids, "amino_acid", "percent", "Amino acid composition", "Amino acids (%)",
         "amino_acid_heatmap", "viridis", None),
        (codons, "codon", "frequency_pct", "Sense codon usage", "Sense codons (%)",
         "codon_frequency_heatmap", "viridis", None),
        (codons, "codon", "rscu", "Relative synonymous codon usage", "RSCU",
         "rscu_heatmap", "vlag", 1),
    ]
    for frame, feature, value, title, scale, filename, palette, center in specifications:
        # 코돈 열은 동의 코돈을 같은 아미노산 아래에 모아서 표시한다.
        order = list(AMINO_ACIDS)
        if feature == "codon":
            order = frame.sort_values(["amino_acid", "codon"])["codon"].drop_duplicates().tolist()
        matrix = frame.pivot(index="species_id", columns=feature, values=value).loc[species_ids, order]
        matrix.index = species_labels
        fig, ax = plt.subplots(figsize=(18 if feature == "codon" else 10, 10))
        ax.set_facecolor("#d9d9d9")
        sns.heatmap(matrix, ax=ax, cmap=palette, center=center, cbar_kws={"label": scale},
                    xticklabels=True, yticklabels=True)
        for boundary in species_order.groupby(GROUPS, sort=False).size().cumsum().iloc[:-1]:
            ax.axhline(boundary, color="white", linewidth=1.3)
        ax.set(title=title, xlabel="DNA codon" if feature == "codon" else "Amino acid",
               ylabel="Species (M: monocot, E: eudicot)")
        save_figure(fig, output / f"{filename}.png")

    aa_groups = summary[summary["measure"].eq("amino_acid_pct")]
    fig, axes = plt.subplots(3, 1, figsize=(12, 9), sharex=True, layout="constrained")
    for ax, (taxon, subset) in zip(axes, aa_groups.groupby("taxonomic_group", sort=True)):
        sns.barplot(data=subset, x="feature", y="mean", hue="photosynthesis_type", ax=ax,
                    order=list(AMINO_ACIDS), hue_order=PHOTOSYNTHESIS, palette=COLORS, errorbar=None)
        counts = subset.groupby("photosynthesis_type")["n_species"].first()
        ax.set(title=f"{taxon} (C3 n={counts['C3']}, C4 n={counts['C4']})",
               xlabel="", ylabel="Mean amino acid (%)")
        ax.legend(title="Photosynthesis")
    axes[-1].set_xlabel("Amino acid")
    save_figure(fig, output / "amino_acid_group_bars.png")

    fig, axes = plt.subplots(2, 3, figsize=(13, 8), layout="constrained")
    box_options = {"data": properties, "x": "taxonomic_group", "hue": "photosynthesis_type",
                   "order": ["Monocot", "Eudicot"], "hue_order": PHOTOSYNTHESIS}
    for ax, (metric, label) in zip(axes.flat, PROPERTIES.items()):
        sns.boxplot(y=metric, ax=ax, palette=COLORS, showfliers=False, **box_options)
        sns.stripplot(y=metric, ax=ax, dodge=True, jitter=False, palette={"C3": "black", "C4": "black"},
                      size=2.7, alpha=0.6, legend=False, **box_options)
        ax.set(xlabel="", ylabel=label)
    axes.flat[-1].set_visible(False)
    save_figure(fig, output / "protein_properties_boxplots.png")


def run(processed, output):
    amino_acids, properties, codons = load_composition(processed)
    labels = ["species_id", *GROUPS]
    aa_long = amino_acids.rename(columns={"amino_acid": "feature", "percent": "value"}).assign(
        measure="amino_acid_pct")
    property_long = properties.melt(
        id_vars=labels, value_vars=list(PROPERTIES), var_name="feature").assign(measure="protein_property")
    codon_long = codons.melt(
        id_vars=[*labels, "codon"], value_vars=["frequency_pct", "rscu"], var_name="measure"
    ).rename(columns={"codon": "feature"})
    summary = group_summary(pd.concat([aa_long, property_long, codon_long], ignore_index=True))

    results = output / "results"
    figures = output / "figures"
    results.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)
    tables = {
        "species_amino_acid_composition": amino_acids,
        "species_protein_properties": properties,
        "species_codon_frequency": codons.drop(columns=["rscu"]),
        "species_rscu": codons.drop(columns=["frequency_pct"]),
        "group_composition_summary": summary,
    }
    for name, table in tables.items():
        table.to_csv(results / f"{name}.csv", index=False, float_format="%.10f")
    plot_results(amino_acids, properties, codons, summary, figures)
    print(f"Analyzed {len(properties)} species. Wrote {len(tables)} CSV files and 5 PNG figures.")


if __name__ == "__main__":
    run(ROOT / "data" / "processed", ROOT / "week13")
