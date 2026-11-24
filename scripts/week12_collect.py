"""Week 12: NCBI GenBank에서 종당 대표 rbcL CDS를 수집·검증한다."""

import argparse
import csv
from datetime import date
from io import StringIO
import logging
import os
from pathlib import Path
import re
import socket
from tempfile import TemporaryDirectory
import time
from urllib.error import HTTPError

from Bio import Entrez, SeqIO
from Bio.Data import CodonTable
from Bio.SeqFeature import ExactPosition


ROOT = Path(__file__).resolve().parents[1]
PILOT = ("oryza_sativa", "flaveria_pringlei", "zea_mays", "flaveria_bidentis")
VALIDATION_COLUMNS = (
    "species_id", "scientific_name", "ncbi_taxonomy_id", "nucleotide_accession",
    "protein_accession", "record_organism", "record_taxonomy_id", "source_tier",
    "cds_identity_evidence", "cds_location", "strand", "translation_table", "translation_table_source",
    "cds_length", "protein_length", "stop_codon", "annotation_translation_match",
    "validation_status", "reason",
)
ACCESSION = re.compile(r"[A-Z][A-Z0-9_]+\.\d+")
MAX_CANDIDATES = 30


class RecordExcluded(ValueError):
    """후보 종 또는 레코드가 명시된 검증 기준을 통과하지 못했다."""


def request(operation, parser, **parameters):
    """Entrez 요청과 파서를 받아 결과를 반환한다. 한 번에 한 요청을 보낸다."""
    for attempt in range(3):
        try:
            # 재시도를 포함해 요청 시작을 최소 0.37초 간격으로 유지한다.
            time.sleep(0.37)
            with operation(**parameters) as response:
                return parser(response)
        except HTTPError as error:
            failure = f"HTTP {error.code}"
            if error.code != 429 and error.code < 500:
                break
        except (OSError, ValueError) as error:
            failure = type(error).__name__
        if attempt < 2:
            time.sleep(2 ** attempt)
    # 오류 객체에는 이메일이나 키가 들어간 URL이 포함될 수 있다.
    raise RuntimeError(f"NCBI request failed ({failure})") from None


def taxon_records(*ids):
    """NCBI Taxonomy ID를 받아 레코드 목록을 반환한다."""
    return request(Entrez.efetch, Entrez.read, db="taxonomy", id=",".join(ids), retmode="xml")


def taxon_names(record):
    """Taxonomy 레코드의 현재 학명과 학명 동의어를 소문자 집합으로 반환한다."""
    names = {str(record["ScientificName"]).casefold()}
    other = record.get("OtherNames", {})
    for key in ("Synonym", "EquivalentName", "GenbankSynonym"):
        values = other.get(key, [])
        if isinstance(values, str):
            values = [values]
        names.update(str(value).casefold() for value in values)
    return names


def identify_species(name):
    """입력 학명을 Taxonomy의 유일한 종 단위 레코드와 대조한다."""
    hits = request(
        Entrez.esearch, Entrez.read, db="taxonomy", term=f'"{name}"[All Names]',
        retmax=20,
    )
    if not hits["IdList"]:
        raise RecordExcluded("no NCBI Taxonomy search hit")
    matches = [
        record for record in taxon_records(*hits["IdList"])
        if record["Rank"] == "species" and name.casefold() in taxon_names(record)
    ]
    if len(matches) != 1:
        raise RecordExcluded("no unique NCBI species match")
    return matches[0]


def candidate_queries(tax_id):
    """종 ID에서 우선순위가 있는 NCBI Nucleotide 검색식 네 개를 반환한다."""
    organism = f"txid{tax_id}[ORGN]"
    plastome = (
        f"{organism} AND (chloroplast[Title] OR plastid[Title] OR organelle[Title]) "
        "AND 50000:300000[SLEN]"
    )
    cds = (
        f'{organism} AND (rbcL[Title] OR "ribulose-1,5-bisphosphate carboxylase"[Title]) '
        'AND 1300:2300[SLEN] NOT "partial cds"[Title]'
    )
    return (
        ("refseq_plastome", plastome + " AND srcdb_refseq[PROP]"),
        ("genbank_plastome", plastome + " NOT srcdb_refseq[PROP]"),
        ("refseq_rbcL_cds", cds + " AND srcdb_refseq[PROP]"),
        ("genbank_rbcL_cds", cds + " NOT srcdb_refseq[PROP]"),
    )


def check_source(record, species):
    """GenBank source의 종명을 확인하고 실제 레코드의 종명·taxid를 반환한다."""
    sources = [feature for feature in record.features if feature.type == "source"]
    if len(sources) != 1:
        raise RecordExcluded("expected one GenBank source feature")
    source = sources[0].qualifiers
    organism = source.get("organism", [""])[0]
    tax_ids = [value[6:] for value in source.get("db_xref", []) if value.startswith("taxon:")]
    if len(tax_ids) != 1 or not organism:
        raise RecordExcluded("missing or ambiguous source organism/taxid")
    source_id = tax_ids[0]
    if source_id == str(species["TaxId"]):
        if organism.casefold() not in taxon_names(species):
            raise RecordExcluded("source name differs from NCBI Taxonomy")
    else:
        # 재배군, 품종, 아종의 레코드는 동일한 종 조상과 실제 이름을 확인한다.
        if not any(organism.casefold().startswith(name + " ") for name in taxon_names(species)):
            raise RecordExcluded("record organism is outside the species")
        taxa = taxon_records(source_id)
        if len(taxa) != 1 or organism.casefold() not in taxon_names(taxa[0]):
            raise RecordExcluded("record organism could not be confirmed in Taxonomy")
        lineage = [str(node["TaxId"]) for node in taxa[0].get("LineageEx", [])]
        if str(species["TaxId"]) not in lineage or re.search(r"\s[x×]\s", organism):
            raise RecordExcluded("source is not an infraspecific record of the target species")
    return organism, source_id


def validate_cds(record, *, allow_standalone_annotation=False):
    """GenBank 레코드에서 완전 rbcL CDS, 단백질과 검증 메타데이터를 반환한다."""
    cds_features = [feature for feature in record.features if feature.type == "CDS"]
    features = [
        feature for feature in cds_features
        if any(
            gene.casefold() == "rbcl" for gene in feature.qualifiers.get("gene", [])
        )
    ]
    identity_evidence = "gene=rbcL"
    # 오래된 독립 레코드에서는 rbcL이 표제에만 있고 CDS에 산물 주석이 있다.
    if (not features and allow_standalone_annotation and len(cds_features) == 1
            and re.search(r"\brbcl\b", record.description, re.IGNORECASE)
            and not cds_features[0].qualifiers.get("gene")):
        annotations = cds_features[0].qualifiers
        labels = " ".join(annotations.get("product", []) + annotations.get("note", []))
        if re.search(r"(rubisco|ribulose[^;]*)\s+large subunit", labels, re.IGNORECASE):
            features = cds_features
            identity_evidence = "standalone_rbcL_title_and_large_subunit_annotation"
    if len(features) != 1:
        raise RecordExcluded("expected one CDS with verified rbcL identity")
    feature = features[0]
    location = feature.location
    qualifiers = feature.qualifiers
    if (
        not location or location.strand not in (-1, 1) or "partial" in qualifiers
        or any(type(part.start) is not ExactPosition or type(part.end) is not ExactPosition
               or part.ref or part.ref_db or part.strand != location.strand
               or not 0 <= int(part.start) < int(part.end) <= len(record)
               for part in location.parts)
    ):
        raise RecordExcluded("CDS has a partial or invalid location")
    if qualifiers.get("codon_start", ["1"]) != ["1"]:
        raise RecordExcluded("complete CDS must have codon_start=1")
    # feature.extract는 역방향과 복합 위치를 GenBank 주석대로 처리한다.
    cds = str(feature.extract(record.seq)).upper()
    if set(cds) - set("ACGT"):
        raise RecordExcluded("CDS contains ambiguous nucleotides")
    if len(cds) < 6 or len(cds) % 3:
        raise RecordExcluded("CDS length is invalid or not a multiple of three")
    try:
        table = int(qualifiers.get("transl_table", ["11"])[0])
        protein = str(feature.extract(record.seq).translate(table=table, cds=True))
    except (CodonTable.TranslationError, ValueError, KeyError) as error:
        raise RecordExcluded(f"invalid complete CDS translation ({type(error).__name__})") from None
    annotated = "".join("".join(qualifiers.get("translation", [])).split())
    if not annotated or protein != annotated:
        raise RecordExcluded("missing or mismatching GenBank translation")
    protein_id = qualifiers.get("protein_id", [""])[0]
    if not ACCESSION.fullmatch(protein_id):
        raise RecordExcluded("missing versioned protein accession")
    details = {
        "protein_accession": protein_id,
        "cds_identity_evidence": identity_evidence,
        "cds_location": str(location),
        "strand": location.strand,
        "translation_table": table,
        "translation_table_source": "GenBank" if "transl_table" in qualifiers else "plant_plastid_default",
        "cds_length": len(cds),
        "protein_length": len(protein),
        "stop_codon": cds[-3:],
        "annotation_translation_match": "true",
    }
    return cds, protein, details


def collect_species(row, raw_dir, used_accessions, logger):
    """후보 한 행에서 선택한 서열 쌍과 검증 정보를 반환하고 원본을 보존한다."""
    details = {key: row[key] for key in ("species_id", "scientific_name")}
    species_id = row["species_id"]
    try:
        species = identify_species(row["scientific_name"])
    except RecordExcluded as error:
        logger.warning("%s: %s", species_id, error)
        return None, {**details, "validation_status": "excluded", "reason": str(error)}

    tax_id = str(species["TaxId"])
    details["ncbi_taxonomy_id"] = tax_id
    logger.info("%s: confirmed NCBI species %s (taxid %s)", species_id, species["ScientificName"], tax_id)
    seen = set()
    for tier, query in candidate_queries(tax_id):
        hits = request(
            Entrez.esearch, Entrez.read, db="nuccore", term=query,
            retmax=MAX_CANDIDATES, idtype="acc",
        )
        logger.info("%s: %s %s hits, first %d returned", species_id, tier, hits["Count"], len(hits["IdList"]))
        for accession in hits["IdList"]:
            if accession in seen:
                continue
            seen.add(accession)
            try:
                if not ACCESSION.fullmatch(accession):
                    raise RecordExcluded("non-versioned NCBI search identifier")
                raw = request(
                    Entrez.efetch, lambda response: response.read(), db="nuccore",
                    id=accession, rettype="gbwithparts", retmode="text",
                )
                record = SeqIO.read(StringIO(raw), "genbank")
                if not ACCESSION.fullmatch(record.id) or record.id != accession:
                    raise RecordExcluded("GenBank accession/version differs from search hit")
                if tier.endswith("_plastome"):
                    title = record.description.casefold()
                    if "chloroplast" not in title and "plastid" not in title:
                        raise RecordExcluded("record title does not identify a plastome")
                    if "complete genome" not in title and record.annotations.get("topology") != "circular":
                        raise RecordExcluded("plastome assembly has no confirmed circular topology")
                if record.id.split(".")[0] in used_accessions:
                    raise RecordExcluded("nucleotide accession already selected for another species")
                organism, source_id = check_source(record, species)
                cds, protein, features = validate_cds(
                    record, allow_standalone_annotation=tier.endswith("_cds")
                )
            except (ValueError, LookupError) as error:
                reason = str(error) if isinstance(error, RecordExcluded) else type(error).__name__
                logger.warning("%s: excluded %s (%s)", species_id, accession, reason)
                continue

            raw_path = raw_dir / f"{record.id}.gb"
            if raw_path.exists():
                if raw_path.read_text(encoding="utf-8") != raw.replace("\r\n", "\n"):
                    raise RuntimeError(f"Archived GenBank differs from fetched record {record.id}")
            else:
                raw_path.write_text(raw, encoding="utf-8", newline="\n")
            used_accessions.add(record.id.split(".")[0])
            logger.info("%s: selected %s (%s)", species_id, record.id, tier)
            return (cds, protein), {
                **details, **features, "nucleotide_accession": record.id,
                "record_organism": organism, "record_taxonomy_id": source_id,
                "source_tier": tier, "validation_status": "validated", "reason": "",
            }

    reason = "no acceptable full-length rbcL CDS in returned candidates"
    logger.warning("%s: %s", species_id, reason)
    return None, {**details, "validation_status": "excluded", "reason": reason}


def write_csv(path, fields, rows):
    """지정한 열과 행을 UTF-8 CSV로 기록한다."""
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_outputs(output_dir, metadata_fields, rows, validation, sequences):
    """수집 단계의 메타데이터, 감사 CSV 및 FASTA를 기록한다."""
    write_csv(output_dir / "species_metadata.csv", metadata_fields, rows)
    write_csv(output_dir / "sequence_validation.csv", VALIDATION_COLUMNS, validation)
    for index, name in enumerate(("rbcL_nucleotide.fasta", "rbcL_protein.fasta")):
        with (output_dir / name).open("w", encoding="utf-8", newline="\n") as handle:
            for species_id, pair in sequences.items():
                sequence = pair[index]
                handle.write(f">{species_id}\n")
                for start in range(0, len(sequence), 70):
                    handle.write(sequence[start:start + 70] + "\n")


def main():
    """환경 연락처와 수집 단계를 읽어 파일을 만들고 성공 여부를 반환한다."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase", required=True, choices=("pilot", "all"))
    args = parser.parse_args()
    email = os.environ.get("NCBI_EMAIL", "").strip()
    if not re.fullmatch(r"[^\s@]+@[^\s@]+\.[^\s@]+", email):
        parser.error("Set NCBI_EMAIL to a valid contact email in the environment.")
    Entrez.email = email
    Entrez.api_key = os.environ.get("NCBI_API_KEY") or None
    Entrez.tool = "rbcl_c3_c4_analysis"
    Entrez.max_tries = 1
    socket.setdefaulttimeout(30)
    processed_dir = ROOT / "data/processed"
    pilot_dir = ROOT / "week12/results/pilot"
    if args.phase == "all":
        pilot_report = pilot_dir / "sequence_validation.csv"
        if not pilot_report.exists():
            parser.error("Complete --phase pilot before --phase all.")
        with pilot_report.open(encoding="utf-8", newline="") as handle:
            pilot_checks = {row["species_id"]: row["validation_status"] for row in csv.DictReader(handle)}
        if pilot_checks != dict.fromkeys(PILOT, "validated"):
            parser.error("All four pilot species must pass before full collection.")
    with (processed_dir / "species_metadata.csv").open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames
        rows = list(reader)
    if len(rows) != 30 or len({row["species_id"] for row in rows}) != 30:
        parser.error("Expected 30 unique species IDs in species_metadata.csv.")
    if args.phase == "pilot":
        by_id = {row["species_id"]: row for row in rows}
        rows = [by_id[species_id] for species_id in PILOT]
    raw_dir = ROOT / "data/raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    log_buffer = StringIO()
    logging.basicConfig(
        level=logging.INFO, format="%(levelname)s %(message)s",
        handlers=[logging.StreamHandler(), logging.StreamHandler(log_buffer)],
    )
    logger = logging.getLogger("week12")
    today = date.today().isoformat()
    logger.info("NCBI collection phase=%s date=%s", args.phase, today)
    sequences = {}
    used_accessions = set()
    validation = []
    try:
        for row in rows:
            pair, details = collect_species(row, raw_dir, used_accessions, logger)
            validation.append(details)
            for key in ("ncbi_taxonomy_id", "nucleotide_accession", "protein_accession",
                        "cds_length", "protein_length", "validation_status"):
                row[key] = details.get(key, "")
            row["download_date"] = today if pair else ""
            row["source_url"] = (
                f"https://www.ncbi.nlm.nih.gov/nuccore/{details['nucleotide_accession']}" if pair else ""
            )
            if pair:
                sequences[row["species_id"]] = pair
    except RuntimeError as error:
        logger.error("Collection interrupted: %s", error)
    logger.info("Validated %d of %d species", len(sequences), len(rows))
    complete = len(validation) == len(rows) == len(sequences)
    output_dir = pilot_dir if args.phase == "pilot" else processed_dir
    if args.phase == "all" and not complete:
        output_dir = ROOT / "week12/results/failed_collection"
    output_dir.mkdir(parents=True, exist_ok=True)
    # 수집이 끝난 파일만 공개 경로에 연결하고 전체 수집 실패는 별도 보존한다.
    with TemporaryDirectory(prefix=".collection-", dir=output_dir) as folder:
        stage = Path(folder)
        (stage / "download_log.txt").write_text(log_buffer.getvalue(), encoding="utf-8", newline="\n")
        names = ["sequence_validation.csv", "download_log.txt"]
        if args.phase == "all" and not complete:
            write_csv(stage / "sequence_validation.csv", VALIDATION_COLUMNS, validation)
        else:
            write_outputs(stage, fields, rows, validation, sequences)
            names = ["rbcL_nucleotide.fasta", "rbcL_protein.fasta", *names, "species_metadata.csv"]
        for name in names:
            (stage / name).replace(output_dir / name)
    return 0 if complete else 1


if __name__ == "__main__":
    raise SystemExit(main())
