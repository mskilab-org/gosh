import csv
import tempfile
import unittest
from pathlib import Path

from gosh_cli.core.outputs import Outputs


class StructuralVariantOutputsTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.results = self.root / "results"
        self.samplesheet = self.root / "samplesheet.csv"
        self.samplesheet.write_text("patient,sample,status\nP,T,1\nP,N,0\n")

    def publish(self, relative):
        path = self.results / "P" / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()
        return str(path)

    def test_esvee_resume_uses_unfiltered_not_assembly_or_germline(self):
        for name in ("raw", "ref_depth", "germline"):
            self.publish(f"esvee/T.esvee.{name}.vcf.gz")
        somatic = self.publish("esvee/T.esvee.somatic.vcf.gz")
        unfiltered = self.publish("esvee/T.esvee.unfiltered.vcf.gz")
        somatic_ffpe = self.publish(
            "sv_chimera_filter/tumor/T.esvee.somatic.ffpe_filtered.vcf.gz"
        )
        raw_ffpe = self.publish(
            "sv_chimera_filter/tumor/T.esvee.unfiltered.ffpe_filtered.vcf.gz"
        )
        self.publish("esvee/T.esvee.somatic.vcf.gz.tbi")
        outputs = Outputs(str(self.results), str(self.samplesheet))
        expected = {
            "structural_variants": somatic,
            "structural_variants_unfiltered": unfiltered,
            "structural_variants_raw": unfiltered,
            "structural_variants_chimera_filtered": somatic_ffpe,
            "structural_variants_raw_chimera_filtered": raw_ffpe,
        }
        record = outputs.outputs[0]
        self.assertEqual({key: record[key] for key in expected}, expected)
        sheet = self.root / "resume.csv"
        outputs.emit_samplesheet_csv(str(sheet))
        with sheet.open() as stream:
            rows = list(csv.DictReader(stream))
        tumor = next(row for row in rows if row["status"] == "1")
        self.assertEqual(tumor["vcf"], somatic)
        self.assertEqual(tumor["vcf_raw"], unfiltered)
        self.assertEqual(tumor["structural_variants_chimera_filtered"], somatic_ffpe)
        self.assertEqual(tumor["structural_variants_raw_chimera_filtered"], raw_ffpe)

    def test_somatic_ffpe_cannot_substitute_for_unfiltered_ffpe(self):
        self.publish("sv_chimera_filter/tumor/T.esvee.somatic.ffpe_filtered.vcf.gz")
        self.publish("esvee/T.esvee.raw.vcf.gz")
        self.publish("esvee/T.esvee.germline.vcf.gz")
        record = Outputs(str(self.results), str(self.samplesheet)).outputs[0]
        self.assertEqual(record["structural_variants_raw_chimera_filtered"], "")
        self.assertEqual(record["structural_variants_raw"], "")
        self.assertEqual(record["structural_variants"], "")

    def test_gridss_precedes_esvee_and_samplesheet_precedes_discovery(self):
        self.publish("esvee/T.esvee.somatic.vcf.gz")
        self.publish("esvee/T.esvee.unfiltered.vcf.gz")
        self.publish("sv_chimera_filter/tumor/T.esvee.somatic.ffpe_filtered.vcf.gz")
        self.publish("sv_chimera_filter/tumor/T.esvee.unfiltered.ffpe_filtered.vcf.gz")
        expected = {
            "structural_variants": self.publish("gridss/T.gripss.filtered.vcf.gz"),
            "structural_variants_unfiltered": self.publish("gridss/T.gripss.vcf.gz"),
            "structural_variants_raw": self.publish("gridss/T.gridss.vcf.gz"),
            "structural_variants_chimera_filtered": self.publish(
                "sv_chimera_filter/tumor/T.gridss.filtered.ffpe_filtered.vcf.gz"
            ),
            "structural_variants_raw_chimera_filtered": self.publish(
                "sv_chimera_filter/tumor/T.gridss.ffpe_filtered.vcf.gz"
            ),
        }
        record = Outputs(str(self.results), str(self.samplesheet)).outputs[0]
        self.assertEqual({key: record[key] for key in expected}, expected)
        self.samplesheet.write_text(
            "patient,sample,status,vcf\nP,T,1,/provided/somatic.vcf.gz\nP,N,0,\n"
        )
        record = Outputs(str(self.results), str(self.samplesheet)).outputs[0]
        self.assertEqual(record["structural_variants"], "/provided/somatic.vcf.gz")


if __name__ == "__main__":
    unittest.main()
