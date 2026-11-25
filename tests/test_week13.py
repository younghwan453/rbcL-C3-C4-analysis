"""13주차의 손 계산, 경계 조건, 입력 연결 및 저장 결과를 확인한다."""

import math
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory
import unittest

from Bio.Data.CodonTable import TranslationError
from Bio.Seq import Seq
from PIL import Image
import pandas as pd

from scripts.week13_composition import PROPERTIES, group_summary, load_composition, run, sequence_composition


ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"


class CompositionTests(unittest.TestCase):
    def test_hand_calculation_and_missing_family(self):
        amino, properties, codons = sequence_composition(Seq("ATGGCTGCTGCCTAA"), "MAAA", 11)
        aa = {row["amino_acid"]: row for row in amino}
        usage = {row["codon"]: row for row in codons}
        self.assertEqual(len(aa), 20)
        self.assertEqual((aa["A"]["count"], aa["A"]["percent"]), (3, 75))
        self.assertEqual((aa["Y"]["count"], aa["Y"]["percent"]), (0, 0))
        self.assertEqual(properties["protein_length"], 4)
        self.assertAlmostEqual(properties["gravy"], (1.9 + 3 * 1.8) / 4)
        self.assertEqual(properties["aromaticity"], 0)
        self.assertEqual((usage["GCT"]["count"], usage["GCT"]["frequency_pct"]), (2, 50))
        self.assertAlmostEqual(usage["GCT"]["rscu"], 8 / 3)
        self.assertAlmostEqual(usage["GCC"]["rscu"], 4 / 3)
        self.assertEqual(usage["GCA"]["rscu"], 0)
        self.assertIsNone(usage["TGG"]["rscu"])
        self.assertNotIn("TAA", usage)

    def test_alternative_start_counts_original_codon(self):
        amino, _, codons = sequence_composition(Seq("GTGTAG"), "M", 11)
        self.assertEqual(sum(row["count"] for row in codons), 1)
        self.assertEqual(next(row["count"] for row in amino if row["amino_acid"] == "M"), 1)
        gtg = next(row for row in codons if row["codon"] == "GTG")
        self.assertEqual((gtg["amino_acid"], gtg["rscu"]), ("V", 4))

    def test_invalid_inputs(self):
        cases = [
            ("ATGGCNTAA", "MA", ValueError), ("ATGGCCTAA", "MX", ValueError),
            ("ATGGCCTAA", "MK", ValueError), ("ATGTAAGCCTAG", "MAA", TranslationError),
            ("ATGGCC", "M", TranslationError), ("ATGGTAA", "MV", TranslationError),
        ]
        for cds, protein, error in cases:
            with self.subTest(cds=cds, protein=protein), self.assertRaises(error):
                sequence_composition(Seq(cds), protein, 11)

    def test_group_summary_is_species_weighted_and_skips_missing(self):
        rows = pd.DataFrame({
            "photosynthesis_type": ["C3"] * 3, "taxonomic_group": ["Monocot", "Monocot", "Eudicot"],
            "measure": ["rscu"] * 3, "feature": ["GCA"] * 3, "value": [0, math.nan, 4],
        })
        row = group_summary(rows).query("taxonomic_group == 'All'").iloc[0]
        self.assertEqual((row.n_species, row.n_available, row["mean"]), (3, 2, 2))
        self.assertAlmostEqual(row.sd, math.sqrt(8))

    def test_rejects_unlinked_protein_id(self):
        with TemporaryDirectory() as directory:
            processed = Path(directory)
            for name in ("species_metadata.csv", "sequence_validation.csv",
                         "rbcL_nucleotide.fasta", "rbcL_protein.fasta"):
                shutil.copyfile(PROCESSED / name, processed / name)
            fasta = processed / "rbcL_protein.fasta"
            fasta.write_text(fasta.read_text().replace(">oryza_sativa", ">unlinked"), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "same nonempty species IDs"):
                load_composition(processed)

    def test_real_dataset_and_output(self):
        with TemporaryDirectory() as directory:
            output = Path(directory)
            run(PROCESSED, output)
            aa = pd.read_csv(output / "results/species_amino_acid_composition.csv")
            codons = pd.read_csv(output / "results/species_codon_frequency.csv")
            rscu = pd.read_csv(output / "results/species_rscu.csv")
            properties = pd.read_csv(output / "results/species_protein_properties.csv").set_index("species_id")
            summary = pd.read_csv(output / "results/group_composition_summary.csv")
            self.assertEqual(len(properties), 30)
            self.assertEqual(set(properties.columns) & set(PROPERTIES), set(PROPERTIES))
            self.assertTrue((aa.groupby("species_id").size() == 20).all())
            self.assertTrue((codons.groupby("species_id").size() == 61).all())
            for table, percentage in [(aa, "percent"), (codons, "frequency_pct")]:
                for species_id, group in table.groupby("species_id"):
                    self.assertEqual(group["count"].sum(), properties.loc[species_id, "protein_length"])
                    self.assertAlmostEqual(group[percentage].sum(), 100)
            for _, group in rscu.groupby(["species_id", "amino_acid"]):
                if group["synonymous_family_count"].iloc[0] > 0:
                    self.assertAlmostEqual(group["rscu"].sum(), len(group))
                else:
                    self.assertTrue(group["rscu"].isna().all())
            self.assertEqual(len(summary), 6 * (20 + len(PROPERTIES) + 61 * 2))
            c4_leucine = aa.query("photosynthesis_type == 'C4' and amino_acid == 'L'")
            reported = summary.query(
                "photosynthesis_type == 'C4' and taxonomic_group == 'All' "
                "and measure == 'amino_acid_pct' and feature == 'L'"
            ).iloc[0]
            self.assertEqual(reported.n_species, 15)
            self.assertAlmostEqual(reported["mean"], c4_leucine["percent"].mean())
            figures = sorted((output / "figures").glob("*.png"))
            self.assertEqual(len(figures), 5)
            for figure in figures:
                with Image.open(figure) as image:
                    image.verify()
                with Image.open(figure) as image:
                    self.assertAlmostEqual(image.info["dpi"][0], 300, delta=0.01)
                    self.assertAlmostEqual(image.info["dpi"][1], 300, delta=0.01)


if __name__ == "__main__":
    unittest.main()
