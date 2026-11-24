"""11주차 후보 CSV의 구조와 모집 설계를 확인한다."""

import csv
from collections import Counter
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
FIELDS = (
    "species_id scientific_name common_name family taxonomic_group "
    "photosynthesis_type photosynthesis_reference ncbi_taxonomy_id "
    "nucleotide_accession protein_accession cds_length protein_length "
    "download_date source_url validation_status"
).split()


class CandidateMetadataTests(unittest.TestCase):
    """입력은 주차 문서가 참조하는 CSV이며 출력은 unittest 결과이다."""

    @classmethod
    def setUpClass(cls):
        """후보 CSV를 읽고 열 이름과 종별 행을 테스트에 제공한다."""
        with (ROOT / "data/processed/species_metadata.csv").open(
            encoding="utf-8", newline=""
        ) as handle:
            reader = csv.DictReader(handle)
            cls.fields = reader.fieldnames
            cls.rows = list(reader)

    def test_schema_and_row_width(self):
        self.assertEqual(self.fields, FIELDS)
        self.assertTrue(all(list(row) == FIELDS for row in self.rows))
        self.assertTrue(all(value is not None for row in self.rows for value in row.values()))

    def test_unique_species_and_identifiers(self):
        self.assertEqual(len(self.rows), 30)
        for field in ("scientific_name", "species_id"):
            self.assertEqual(len({row[field] for row in self.rows}), 30)
        for row in self.rows:
            self.assertEqual(row["species_id"], row["scientific_name"].lower().replace(" ", "_"))

    def test_four_group_counts(self):
        counts = Counter(
            (row["photosynthesis_type"], row["taxonomic_group"]) for row in self.rows
        )
        self.assertEqual(counts, {
            ("C3", "Monocot"): 8,
            ("C3", "Eudicot"): 7,
            ("C4", "Monocot"): 8,
            ("C4", "Eudicot"): 7,
        })

    def test_pending_rows_leave_sequence_fields_empty(self):
        for row in self.rows:
            self.assertTrue(row["validation_status"])
            if row["validation_status"] == "pending_sequence_validation":
                self.assertTrue(all(row[field] == "" for field in FIELDS[7:-1]))

    def test_every_species_has_a_source(self):
        for row in self.rows:
            self.assertTrue(row["photosynthesis_reference"].startswith("https://"), row["species_id"])

    def test_flaveria_pairing_and_document_link(self):
        c3 = {"Flaveria pringlei", "Flaveria robusta"}
        c4 = {"Flaveria bidentis", "Flaveria trinervia"}
        actual = {
            row["scientific_name"]: row["photosynthesis_type"] for row in self.rows
            if row["scientific_name"].startswith("Flaveria ")
        }
        self.assertEqual(actual, dict.fromkeys(c3, "C3") | dict.fromkeys(c4, "C4"))
        readme = (ROOT / "week11/README.md").read_text(encoding="utf-8")
        self.assertIn("(../data/processed/species_metadata.csv)", readme)

    def test_four_pilot_species_cover_all_groups(self):
        expected = {
            "Oryza sativa": ("C3", "Monocot"),
            "Flaveria pringlei": ("C3", "Eudicot"),
            "Zea mays": ("C4", "Monocot"),
            "Flaveria bidentis": ("C4", "Eudicot"),
        }
        actual = {
            row["scientific_name"]: (row["photosynthesis_type"], row["taxonomic_group"])
            for row in self.rows if row["scientific_name"] in expected
        }
        self.assertEqual(actual, expected)
        plan = (ROOT / "week11/research_plan.md").read_text(encoding="utf-8")
        self.assertTrue(all(species in plan for species in expected))
