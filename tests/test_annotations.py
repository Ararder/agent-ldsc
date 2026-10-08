import numpy as np
import pandas as pd
import pytest

from agent_ldsc_worker.annotations import normalize, parse_bed
from agent_ldsc_worker.intervals import merge_intervals, snp_membership
from agent_ldsc_worker.util import WorkerError, normalize_chrom


def write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return p


def bed_spec(**kw):
    return {"id": "a", "type": "bed", "file": "x.bed", "species": "human", "build": "GRCh37", **kw}


def gene_spec(**kw):
    return {"id": "g", "type": "genes", "file": "g.txt", "species": "human", "id_type": "symbol",
            "gene_model": "gene_span", "flank_bp": 0, **kw}


# ---------------------------------------------------------------- interval primitives

def test_membership_off_by_one_matches_make_annot():
    merged = pd.DataFrame({"chrom": [1], "start": [100], "end": [200]})
    bp = np.array([100, 101, 200, 201])  # 1-based SNP positions
    # [100,200) covers 1-based 101..200
    assert snp_membership(merged, 1, bp).tolist() == [0, 1, 1, 0]
    assert snp_membership(merged, 2, bp).tolist() == [0, 0, 0, 0]


def test_merge_overlapping_and_bookended_only_within_chromosome():
    df = pd.DataFrame({"chrom": [2, 1, 1, 1, 10], "start": [5, 10, 15, 30, 0], "end": [9, 20, 30, 40, 5]})
    m = merge_intervals(df)
    assert m.values.tolist() == [[1, 10, 40], [2, 5, 9], [10, 0, 5]]  # numeric chrom order: 2 before 10


def test_chrom_aliases():
    assert [normalize_chrom(c) for c in ["chr1", "1", "chr22", "chrX", "X", "chr01", "chr23", "chr6_ssto_hap7"]] == \
        [1, 1, 22, None, None, None, None, None]


# ---------------------------------------------------------------- BED

def test_bed_parse_rejects_malformed(tmp_path):
    p = write(tmp_path, "x.bed", "track name=x\nchr1\t10\t5\nchr1\t-1\t5\nchr1\ta\t5\nchr1\t5\n")
    with pytest.raises(WorkerError) as e:
        parse_bed(p)
    assert e.value.code == "INPUT_INVALID"
    assert [d["line"] for d in e.value.details] == [2, 3, 4, 5]


def test_bed_grch37_direct_merges_overlapping_peaks(tmp_path, bundle, strict_profile):
    p = write(tmp_path, "x.bed", "chr1\t100000\t100500\tp1\nchr1\t100400\t101000\tp2\n1\t200000\t200100\n")
    ann = normalize(bed_spec(), p, bundle, strict_profile)
    assert ann.intervals.values.tolist() == [[1, 100000, 101000], [1, 200000, 200100]]
    assert ann.summary["kept_bp_before_merge"] == 500 + 600 + 100
    assert ann.summary["merged_bp"] == 1000 + 100


def test_bed_non_autosomal_requires_allow_drop(tmp_path, bundle, strict_profile):
    p = write(tmp_path, "x.bed", "chr1\t100000\t100500\nchrX\t1\t100\n")
    with pytest.raises(WorkerError) as e:
        normalize(bed_spec(), p, bundle, strict_profile)
    assert e.value.code == "MAPPING_LOSS"
    ann = normalize(bed_spec(allow_drop=["non_autosomal"]), p, bundle, strict_profile)
    assert ann.summary["records_dropped"] == {"non_autosomal": 1}
    assert len(ann.audit) == 2


def test_bed_undeclared_build_rejected(tmp_path, bundle, strict_profile):
    p = write(tmp_path, "x.bed", "chr1\t100000\t100500\n")
    spec = bed_spec()
    del spec["build"]
    with pytest.raises(WorkerError) as e:
        normalize(spec, p, bundle, strict_profile)
    assert e.value.code == "BUILD_UNSUPPORTED"


def test_bed_past_chromosome_end_rejected(tmp_path, bundle, strict_profile):
    p = write(tmp_path, "x.bed", "chr1\t2999990\t3000010\n")
    with pytest.raises(WorkerError) as e:
        normalize(bed_spec(), p, bundle, strict_profile)
    assert e.value.code == "INPUT_INVALID"


@pytest.mark.liftover
def test_liftover_maps_and_records_deleted(tmp_path, bundle, strict_profile):
    # chain: hg38 chr1 [0,500k) -> +1000 ; [500k,510k) deleted ; chr2 identity
    p = write(tmp_path, "x.bed", "chr1\t100000\t100500\tok\nchr1\t502000\t503000\tgone\nchr2\t5000\t6000\tid\n")
    with pytest.raises(WorkerError) as e:
        normalize(bed_spec(build="GRCh38"), p, bundle, strict_profile)
    assert e.value.code == "MAPPING_LOSS"
    ann = normalize(bed_spec(build="hg38", allow_drop=["liftover_unmapped"]), p, bundle, strict_profile)
    assert ann.intervals.values.tolist() == [[1, 101000, 101500], [2, 5000, 6000]]
    assert ann.summary["records_dropped"] == {"liftover_deleted": 1}
    assert ann.summary["liftover"]["min_match"] == 0.95
    row = ann.audit.set_index("input_name").loc["gone"]
    assert (row.status, row.reason) == ("dropped", "liftover_deleted")


@pytest.mark.liftover
def test_liftover_partial_overlap_of_deletion_is_dropped(tmp_path, bundle, strict_profile):
    # 50% of this interval lies in the deleted segment -> below minMatch 0.95
    p = write(tmp_path, "x.bed", "chr1\t499000\t501000\n")
    ann_err = None
    try:
        normalize(bed_spec(build="GRCh38", allow_drop=["liftover_unmapped"]), p, bundle, strict_profile)
    except WorkerError as e:
        ann_err = e
    assert ann_err is not None and ann_err.code == "INPUT_INVALID"  # nothing left


# ---------------------------------------------------------------- genes

def test_gene_symbol_alias_prev_and_flank(tmp_path, bundle, strict_profile):
    p = write(tmp_path, "g.txt", "# my genes\nGENEA\nOLDB\nGENEA\n")
    ann = normalize(gene_spec(flank_bp=1000), p, bundle, strict_profile)
    a = ann.audit.set_index("input_line")
    assert a.loc[2, "gene_id"] == "ENSG00000000001"
    assert a.loc[3, "gene_id"] == "ENSG00000000002"  # unique previous symbol
    assert a.loc[4, "status"] == "duplicate"
    # GENEA 1-based 100001..140000 -> BED [100000,140000) +- 1000
    assert ann.intervals.values.tolist() == [[1, 99000, 141000], [2, 899000, 921000]]


def test_gene_flank_clipped_to_chromosome(tmp_path, bundle, strict_profile):
    p = write(tmp_path, "g.txt", "GENEE\n")
    ann = normalize(gene_spec(flank_bp=100_000), p, bundle, strict_profile)
    assert ann.intervals.values.tolist() == [[22, 2_750_000, 3_000_000]]


def test_gene_ambiguous_alias_fails_with_candidates(tmp_path, bundle, strict_profile):
    p = write(tmp_path, "g.txt", "AMBIG\n")
    with pytest.raises(WorkerError) as e:
        normalize(gene_spec(allow_drop=["gene_unmapped"]), p, bundle, strict_profile)
    assert e.value.code == "GENE_AMBIGUOUS"
    assert e.value.details[0]["candidates"] == "HGNC:1:GENEA|HGNC:2:GENEB"


def test_gene_losses_need_allow_drop(tmp_path, bundle, strict_profile):
    p = write(tmp_path, "g.txt", "GENEA\nGENEC\nGENED\nGENEF\nNOPE\ngenea\n")
    with pytest.raises(WorkerError) as e:
        normalize(gene_spec(), p, bundle, strict_profile)
    assert e.value.code == "MAPPING_LOSS"
    ann = normalize(gene_spec(allow_drop=["gene_unmapped", "non_autosomal"]), p, bundle, strict_profile)
    assert ann.summary["records_dropped"] == {
        "gene_no_ensembl_id": 1, "gene_remap_partial_1": 1, "gene_symbol_not_found": 2, "non_autosomal": 1}
    assert ann.summary["genes_mapped"] == 1


def test_gene_ensembl_ids_strip_version(tmp_path, bundle, strict_profile):
    p = write(tmp_path, "g.txt", "gene_id\nENSG00000000001.5\nENSG00000000002\nENSG123\n")
    ann = normalize(gene_spec(id_type="ensembl_gene", header=True, allow_drop=["gene_unmapped"]), p, bundle,
                    strict_profile)
    assert ann.audit.gene_id.tolist()[:2] == ["ENSG00000000001", "ENSG00000000002"]
    assert ann.summary["records_dropped"] == {"gene_invalid_id": 1}


def test_genes_and_equivalent_bed_give_identical_canonical_annotation(tmp_path, bundle, strict_profile):
    g = normalize(gene_spec(flank_bp=500), write(tmp_path, "g.txt", "GENEA\n"), bundle, strict_profile)
    b = normalize(bed_spec(), write(tmp_path, "x.bed", "chr1\t99500\t140500\n"), bundle, strict_profile)
    assert g.content_hash == b.content_hash
