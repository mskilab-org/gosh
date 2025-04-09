import os
import re
import glob
import csv

# Define the output keys (modify in one place if needed)
OUTPUT_KEYS = [
    "patient_id",
    "sample_ids",
    "tumor_type",
    "disease",
    "primary_site",
    "sex",
    "tumor_bam",
    "normal_bam",
    "msisensorpro",
    "structural_variants",
    "structural_variants_unfiltered",
    "tumor_coverage",
    "snvs_somatic",
    "snvs_somatic_unfiltered",
    "snvs_germline",
    "het_pileups",
    "purple_pp_range",
    "purple_pp_best_fit",
    "multiplicity",
    "multiplicity_germline",
    "multiplicity_hetsnps",
    "variant_annotations_somatic",
    "variant_annotations_germline",
    "oncokb_snv",
    "oncokb_cna",
    "oncokb_fusions",
    "jabba_gg",
    "karyograph",
    "jabba_gg_balanced",
    "events",
    "fusions",
    "jabba_gg_allelic",
    "signatures_activities_sbs",
    "signatures_matrix_sbs",
    "signatures_decomposed_sbs",
    "signatures_activities_indel",
    "signatures_matrix_indel",
    "signatures_decomposed_indel",
    "hrdetect",
    "onenesstwoness",
    "qc_dup_rate",
    "qc_alignment_summary",
    "qc_insert_size",
    "qc_coverage_metrics",
    "qc_coverage_metrics_tumor",
    "qc_coverage_metrics_normal",
]

# Map each output key to its file regex pattern(s)
OUTPUT_FILES_MAPPING = {
    "msisensorpro": r"msisensorpro/.*_somatic$",
    "structural_variants": [
        r"gridss.*/.*/.*high_confidence_somatic\.vcf\.bgz$",
        r"tumor_only_junction_filter/.*/.*somatic\.filtered\.sv\.rds$"
    ],
    "structural_variants_unfiltered": r"gridss.*/.*\.gridss\.filtered\.vcf\.gz$",
    "tumor_coverage": r"dryclean/tumor/.*/drycleaned\.cov\.rds$",
    "snvs_somatic": r"sage/somatic/tumor_only_filter/.*/.*\.sage\.pass_filtered\.tumoronly\.vcf\.gz$",
    "snvs_somatic_unfiltered": r"sage/somatic/.*/.*sage\.somatic\.vcf\.gz$",
    "snvs_germline": r"sage/germline/.*/.*sage\.germline\.vcf\.gz$",
    "het_pileups": r"amber/sites\.txt$",
    "purple_pp_range": r"purple/.*purple\.purity\.range\.tsv$",
    "purple_pp_best_fit": r"purple/.*purple\.purity\.tsv$",
    "multiplicity": r"snv_multiplicity/.*/.*est_snv_cn_somatic\.rds$",
    "multiplicity_germline": r"snv_multiplicity/.*/.*est_snv_cn_germline\.rds$",
    "multiplicity_hetsnps": r"snv_multiplicity/.*/.*est_snv_cn_hetsnps\.rds$",
    "variant_annotations_somatic": r"snpeff/somatic/.*/.*ann\.bcf$",
    "variant_annotations_germline": r"snpeff/germline/.*/.*ann\.bcf$",
    "oncokb_snv": r"oncokb/.*/oncokb_snv\.rds$",
    "oncokb_cna": r"oncokb/.*/oncokb_cna\.rds$",
    "oncokb_fusions": r"oncokb/.*/oncokb_fusions\.rds$",
    "jabba_gg": r"jabba/.*/jabba\.simple\.gg\.rds$",
    "karyograph": r"jabba/.*/karyograph\.rds$",
    "jabba_gg_balanced": r"non_integer_balance/.*/non_integer\.balanced\.gg\.rds$",
    "events": r"events/.*/complex\.rds$",
    "fusions": r"fusions/.*/fusions\.rds$",
    "jabba_gg_allelic": r"lp_phased_balance/.*/lp_phased\.balanced\.gg\.rds$",
    "signatures_activities_sbs": r"sigprofilerassignment/sbs_results/Assignment_Solution/Activities/sbs_Assignment_Solution_Activities\.txt",
    "signatures_matrix_sbs": r"sigprofilerassignment/sig_inputs/output/SBS/sigmat_results\.SBS96\.all",
    "signatures_decomposed_sbs": r"sigprofilerassignment/sbs_results/Assignment_Solution/Activities/Decomposed_MutationType_Probabilities\.txt",
    "signatures_activities_indel": r"sigprofilerassignment/indel_results/Assignment_Solution/Activities/indel_Assignment_Solution_Activities\.txt",
    "signatures_matrix_indel": r"sigprofilerassignment/sig_inputs/output/ID/sigmat_results\.ID83\.all",
    "signatures_decomposed_indel": r"sigprofilerassignment/indel_results/.*/Decomposed_MutationType_Probabilities\.txt",
    "hrdetect": r"hrdetect/.*/hrdetect_results\.rds",
    "onenesstwoness": r"oneness_twoness/.*/oneness_twoness\.rds$",
    "qc_dup_rate": r"gatk_qc/.*/.*metrics",
    "qc_alignment_summary": r"picard_qc/.*/.*alignment_summary_metrics",
    "qc_insert_size": r"picard_qc/.*/.*insert_size_metrics",
    "qc_coverage_metrics": r"picard_qc/.*/.*coverage_metrics",
    "qc_coverage_metrics_tumor": r"picard_qc/tumor/.*/.*coverage_metrics",
    "qc_coverage_metrics_normal": r"picard_qc/normal/.*/.*coverage_metrics",
}

class Outputs:
    def __init__(self, outputs_dir: str, samplesheet: str):
        self.outputs_dir = outputs_dir
        self.samplesheet = samplesheet
        self.samples_data = self._read_samplesheet()
        self.outputs = self._collect_outputs()

    def _read_samplesheet(self) -> dict:
        """
        Read the samplesheet CSV and return a dictionary keyed by patient_id.
        Each value should include:
        - "sample_ids": a list of sample IDs for the patient
        - "metadata": a dict of additional metadata (if provided) that matches output keys.
        """
        patient_data = {}
        with open(self.samplesheet, newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                patient_id = row.get("patient", "").strip()
                if not patient_id:
                    continue
                if patient_id not in patient_data:
                    patient_data[patient_id] = {
                        "patient_id": patient_id,
                        "sample_ids": [],
                        "tumor_type": "",
                        "disease": "",
                        "primary_site": "",
                        "sex": "",
                        "tumor_bam": "",
                        "normal_bam": "",
                    }
                # Append the sample_id
                sample_id = row.get("sample_id", "").strip()
                if sample_id:
                    patient_data[patient_id]["sample_ids"].append(sample_id)
                # For every key defined in OUTPUT_KEYS (besides patient_id and sample_ids),
                # if the row contains that column, store it (prefer this over any output file)
                for key in OUTPUT_KEYS:
                    if key in ("patient_id", "sample_ids"):
                        continue
                    if row.get(key):
                        if key in ("tumor_type", "disease", "primary_site", "sex"):
                            # Store these keys directly at the top level
                            patient_data[patient_id][key] = row[key].strip()
        return patient_data

    def _collect_outputs(self) -> list:
        """
        For each patient_id from the samplesheet, scan the outputs directory to find files matching
        the regex patterns defined in OUTPUT_FILES_MAPPING. For each output key, use the value from
        the samplesheet metadata if present; otherwise, set it to the matched file path.
        Use an empty string if nothing is found.
        """
        outputs_list = []
        for patient_id, data in self.samples_data.items():
            # Initialize with empty strings for all keys
            record = {key: "" for key in OUTPUT_KEYS}
            record["patient_id"] = patient_id
            record["sample_ids"] = data.get("sample_ids", [])
            # Overwrite with samplesheet metadata where provided
            for key, value in data.get("metadata", {}).items():
                record[key] = value
            # Overwrite with top-level keys where provided
            for key in ("tumor_type", "disease", "primary_site", "sex"):
                if key in data:
                    record[key] = data[key]

            patient_dir = os.path.join(self.outputs_dir, patient_id)
            # For each output key that is still empty, search in patient's outputs
            for key, pattern in OUTPUT_FILES_MAPPING.items():
                if record.get(key):
                    continue  # prefer samplesheet value if available
                patterns = pattern if isinstance(pattern, list) else [pattern]
                for pat in patterns:
                    # Use glob to recursively search for files under the patient directory
                    search_pattern = os.path.join(patient_dir, "**", "*")
                    for filepath in glob.glob(search_pattern, recursive=True):
                        rel_path = os.path.relpath(filepath, patient_dir)
                        if re.search(pat, rel_path):
                            record[key] = filepath  # assign full path or process as needed
                            break
                    if record.get(key):
                        break
            outputs_list.append(record)
        return outputs_list

    def emit_output_csv(self, csv_path: str):
        """
        Write the collected outputs (self.outputs) to a CSV file at csv_path.
        The CSV should include all keys except 'sample_ids' (to preserve one-to-one mapping).
        For any missing values, output an empty string.
        """
        # Fields for CSV omit 'sample_ids'
        fieldnames = [key for key in OUTPUT_KEYS if key != "sample_ids"]
        with open(csv_path, "w", newline="") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for row in self.outputs:
                # Construct a row excluding 'sample_ids'
                writer.writerow({key: row.get(key, "") for key in fieldnames})
