import csv
import gzip
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/prepare_fenland_coloc.py"
spec = importlib.util.spec_from_file_location("prepare_fenland_coloc", SCRIPT)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def tsv(path, header, records):
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=header, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(records)


class PrepareFenlandColocTests(unittest.TestCase):
    def test_selection_filters_conflicts_and_mvp_stream_match(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            genes = d / "genes.tsv"
            tsv(genes, ["gene_id", "symbol", "available_followup", "author_reference"], [
                dict(gene_id="G1", symbol="ONE", available_followup="False", author_reference="True")])
            targets = d / "targets.tsv"
            target_header = ["pheno", "chr", "region_start", "region_end", "Assay", "Uniprot_ID", "ensembl_gene_id", "id", "assay_multi_protein"]
            target_rows = [
                dict(pheno="x", chr="1", region_start="1", region_end="9", Assay="A", Uniprot_ID="U1", ensembl_gene_id="G1", id="a.1", assay_multi_protein="False"),
                dict(pheno="y", chr="1", region_start="10", region_end="19", Assay="B", Uniprot_ID="U1", ensembl_gene_id="G1", id="b.1", assay_multi_protein="False"),
            ]
            tsv(targets, target_header, target_rows)
            indices = d / "indices.tsv"
            index_header = ["id", "rsid", "log10p", "sdY_estimated"]
            tsv(indices, index_header, [dict(id="a.1", rsid="rs1", log10p="2", sdY_estimated="2"),
                                        dict(id="b.1", rsid="rs2", log10p="3", sdY_estimated="4")])
            regions = d / "regions.tsv.gz"
            region_header = ["chr", "rsid", "pos", "Allele2", "Allele1", "Freq1", "info", "Effect", "StdErr", "log10p", "id"]
            data = [
                dict(chr="1", rsid="rs10", pos="10", Allele2="A", Allele1="C", Freq1=".2", info=".8", Effect="8", StdErr="2", log10p="4", id="a.1"),
                dict(chr="1", rsid="rs20", pos="20", Allele2="A", Allele1="C", Freq1=".2", info=".8", Effect="12", StdErr="3", log10p="7", id="b.1"),
                # Conflict group: same rsID, incompatible position/alleles, both must be removed.
                dict(chr="1", rsid="rs30", pos="30", Allele2="A", Allele1="C", Freq1=".3", info=".8", Effect="1", StdErr="1", log10p="1", id="b.1"),
                dict(chr="1", rsid="rs30", pos="31", Allele2="G", Allele1="T", Freq1=".3", info=".8", Effect="1", StdErr="1", log10p="1", id="b.1"),
                dict(chr="1", rsid="rs40", pos="40", Allele2="A", Allele1="C", Freq1=".0009", info=".8", Effect="1", StdErr="1", log10p="1", id="b.1"),
            ]
            with gzip.open(regions, "wt", newline="", encoding="utf-8") as f:
                w = csv.DictWriter(f, fieldnames=region_header, delimiter="\t", lineterminator="\n")
                w.writeheader(); w.writerows(data)
            mvp = d / "mvp.ma"
            mvp.write_text("SNP\tA1\tA2\tfreq\tb\tse\tp\tn\nrs20\tC\tA\t.2\t.1\t.1\t.5\t100\nrs30\tC\tA\t.3\t.1\t.1\t.5\t100\nrs99\tA\tG\t.4\t.1\t.1\t.5\t100\n", encoding="utf-8")
            out = d / "out"
            record = mod.build({"genes": genes, "targets": targets, "indices": indices, "regions": regions, "mvp": mvp}, out)
            with gzip.open(out / "full_cis.tsv.gz", "rt", encoding="utf-8") as f:
                cis = list(csv.DictReader(f, delimiter="\t"))
            self.assertEqual([r["SNP"] for r in cis], ["rs20"])
            self.assertEqual(cis[0]["beta"], "3.0")
            self.assertEqual(record["conflict_rsids_excluded"], 1)
            self.assertEqual(record["requested_ids_missing_from_mvp"], 0)
            with gzip.open(out / "mvp_matched.ma.gz", "rt", encoding="utf-8") as f:
                self.assertEqual(sum(1 for _ in f), 2)
            self.assertEqual(json.loads((out / "run.json").read_text())["selected_assays"], 1)


if __name__ == "__main__":
    unittest.main()
