"""Week 12 단위·파일 연결 테스트. NCBI 호출 없이 합성 서열을 사용한다."""

from io import StringIO
import logging
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import BeforePosition, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
from PIL import Image

from scripts.week12_collect import (
    RecordExcluded, check_source, collect_species, validate_cds, write_outputs,
)
from scripts.week12_gc import group_summary, load_species_metrics, make_figures, nucleotide_metrics


TOY_CDS = "ATGGCTGCATAA"
TOY_PROTEIN = "MAA"
TOY_TAXON = {"ScientificName": "Example species", "TaxId": "12345"}
TOY_ROW = {"species_id": "example_species", "scientific_name": "Example species"}


def make_record(sequence=TOY_CDS, protein=TOY_PROTEIN, strand=1, table=None, partial=False):
    """검증 사례에 필요한 실제 GenBank feature 구조의 합성 레코드를 만든다."""
    bases = Seq(sequence)
    record = SeqRecord(
        bases if strand == 1 else bases.reverse_complement(), id="TEST000001.1",
        description="Example species chloroplast, complete genome",
    )
    record.annotations["molecule_type"] = "DNA"
    record.annotations["organism"] = TOY_TAXON["ScientificName"]
    source = {"organism": [TOY_TAXON["ScientificName"]], "db_xref": ["taxon:12345"]}
    qualifiers = {"gene": ["rbcL"], "protein_id": ["TEST000002.1"], "translation": [protein]}
    if table is not None:
        qualifiers["transl_table"] = [str(table)]
    start = BeforePosition(0) if partial else 0
    record.features = [
        SeqFeature(SimpleLocation(0, len(bases), strand=1), type="source", qualifiers=source),
        SeqFeature(SimpleLocation(start, len(bases), strand=strand), type="CDS", qualifiers=qualifiers),
    ]
    return record


class CdsTests(unittest.TestCase):
    def test_reverse_strand_and_default_table(self):
        record = make_record(strand=-1)
        cds, protein, data = validate_cds(record)
        self.assertEqual((cds, protein), (TOY_CDS, TOY_PROTEIN))
        self.assertEqual((data["strand"], data["translation_table"]), (-1, 11))

    def test_explicit_translation_table(self):
        record = make_record("ATGTGATAA", "MW", table=4)
        _, protein, data = validate_cds(record)
        self.assertEqual((protein, data["translation_table_source"]), ("MW", "GenBank"))

    def test_standalone_legacy_annotation_requires_both_title_and_product(self):
        record = make_record()
        record.description = "Example species chloroplast rbcL gene"
        feature = record.features[1]
        feature.qualifiers.pop("gene")
        feature.qualifiers["note"] = ["unnamed protein product; rubisco large subunit"]
        with self.assertRaises(RecordExcluded):
            validate_cds(record)
        _, _, data = validate_cds(record, allow_standalone_annotation=True)
        self.assertEqual(data["cds_identity_evidence"], "standalone_rbcL_title_and_large_subunit_annotation")
        feature.qualifiers["note"] = ["unrelated large subunit"]
        with self.assertRaises(RecordExcluded):
            validate_cds(record, allow_standalone_annotation=True)

    def test_partial_and_ambiguous_cds_are_excluded(self):
        for case in (make_record(partial=True), make_record("ATGGCNGCATAA"), make_record("ATGGCATA")):
            with self.subTest(record=case.seq), self.assertRaises(RecordExcluded):
                validate_cds(case)

    def test_internal_stop_and_annotation_mismatch_are_excluded(self):
        for case in (make_record("ATGTAAGCATAA"), make_record(protein="MAV"), make_record("ATGGCTGCAAAA")):
            with self.subTest(record=case.seq), self.assertRaises(RecordExcluded):
                validate_cds(case)

    def test_source_taxonomy_must_match(self):
        record = make_record()
        self.assertEqual(check_source(record, TOY_TAXON), ("Example species", "12345"))
        record.features[0].qualifiers["organism"] = ["Example other"]
        with self.assertRaisesRegex(RecordExcluded, "source name"):
            check_source(record, TOY_TAXON)

    @patch("scripts.week12_collect.request")
    @patch("scripts.week12_collect.identify_species", return_value=TOY_TAXON)
    def test_refseq_selection_raw_preservation_and_duplicate_accession(self, identify, request):
        record = make_record()
        buffer = StringIO()
        SeqIO.write(record, buffer, "genbank")
        raw = buffer.getvalue()
        request.side_effect = [
            {"Count": "1", "IdList": [record.id]}, raw,
            {"Count": "1", "IdList": [record.id]}, raw,
            *[{"Count": "0", "IdList": []}] * 3,
        ]
        with TemporaryDirectory() as folder:
            logger = logging.getLogger("test_week12")
            used = set()
            path = Path(folder)
            with self.assertLogs(logger, level="INFO") as logs:
                pair, result = collect_species(TOY_ROW, path, used, logger)
                other, excluded = collect_species(TOY_ROW, path, used, logger)
            self.assertEqual(pair, (TOY_CDS, TOY_PROTEIN))
            self.assertTrue(identify.called)
            self.assertEqual(result["source_tier"], "refseq_plastome")
            self.assertEqual((path / f"{record.id}.gb").read_text(encoding="utf-8"), raw)
            self.assertEqual(used, {"TEST000001"})
            self.assertIsNone(other)
            self.assertEqual(excluded["validation_status"], "excluded")
            self.assertTrue(any("already selected" in message for message in logs.output))


class GcTests(unittest.TestCase):
    def test_hand_calculation(self):
        result = nucleotide_metrics(TOY_CDS, 11)
        self.assertEqual([result[f"{base}_count"] for base in "ATGC"], [4, 3, 3, 2])
        self.assertEqual((result["cds_length"], result["composition_nt_count"]), (12, 9))
        for key, expected in {"at_pct": 400/9, "gc_pct": 500/9, "gc1_pct": 200/3,
                              "gc2_pct": 200/3, "gc3_pct": 100/3}.items():
            with self.subTest(metric=key):
                self.assertAlmostEqual(result[key], expected)

    def test_one_sense_codon_is_a_valid_boundary(self):
        result = nucleotide_metrics("ATGTAA", 11)
        self.assertEqual(result["sense_codon_count"], 1)
        self.assertEqual((result["gc1_pct"], result["gc2_pct"], result["gc3_pct"]), (0, 0, 100))

    def test_invalid_sequences_raise_errors(self):
        for sequence in ("", "ATG", "ATGTAAN", "ATGCCC", "ATGTAATAA"):
            with self.subTest(sequence=sequence), self.assertRaises(ValueError):
                nucleotide_metrics(sequence, 11)

    def make_linked_files(self, path):
        """서로 연결된 합성 메타데이터·검증 CSV 및 FASTA 파일을 만든다."""
        validation, metadata, sequences = [], [], {}
        for index, (photo, taxon) in enumerate(
            (("C3", "Monocot"), ("C3", "Eudicot"), ("C4", "Monocot"), ("C4", "Eudicot")), 1
        ):
            species_id = f"test_{index}"
            accession = f"TEST00000{index}.1"
            row = {"species_id": species_id, "scientific_name": f"Test species {index}",
                   "nucleotide_accession": accession, "cds_length": len(TOY_CDS),
                   "validation_status": "validated"}
            metadata.append({**row, "photosynthesis_type": photo, "taxonomic_group": taxon})
            validation.append({**row, "translation_table": 11, "stop_codon": "TAA"})
            sequences[species_id] = (TOY_CDS, TOY_PROTEIN)
        write_outputs(path, list(metadata[0]), metadata, validation, sequences)

    def test_file_integration_group_tables_and_300_dpi_figures(self):
        with TemporaryDirectory() as folder:
            path = Path(folder)
            self.make_linked_files(path)
            species = load_species_metrics(path)
            self.assertEqual(len(species), 4)
            self.assertTrue(species["gc_pct"].eq(500 / 9).all())
            figures = path / "figures"
            figures.mkdir()
            make_figures(species, figures)
            tables = group_summary(species)
            by_photo = tables.loc[tables["grouping"].eq("Photosynthesis type")]
            self.assertTrue(by_photo["n"].eq(2).all())
            self.assertTrue(by_photo["sd_pct"].eq(0).all())
            self.assertEqual(len(list(figures.glob("*.png"))), 4)
            for file in figures.glob("*.png"):
                with Image.open(file) as image:
                    self.assertTrue(all(abs(dpi - 300) < 0.1 for dpi in image.info["dpi"]))

    def test_missing_fasta_id_is_rejected(self):
        with TemporaryDirectory() as folder:
            path = Path(folder)
            self.make_linked_files(path)
            (path / "rbcL_nucleotide.fasta").write_text(">extra\nATGTAA\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "same species IDs"):
                load_species_metrics(path)


if __name__ == "__main__":
    unittest.main()
