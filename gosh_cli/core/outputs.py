import os
import re
import glob
import csv
import sys
from typing import Optional, List

# Define the output keys (modify in one place if needed)
OUTPUT_KEYS = [
    "patient_id",
    "sample_ids",
    "tumor_sample",
    "normal_sample",
    "tumor_type",
    "disease",
    "primary_site",
    "sex",
    "bam_tumor",
    "bam_normal",
    "qc_dup_rate",
    "qc_dup_rate_tumor",
    "qc_dup_rate_normal",
    "qc_alignment_summary",
    "qc_alignment_summary_tumor",
    "qc_alignment_summary_normal",
    "qc_insert_size",
    "qc_insert_size_tumor",
    "qc_insert_size_normal",
    "qc_coverage_metrics",
    "qc_coverage_metrics_tumor",
    "qc_coverage_metrics_normal",
    "msisensorpro",
    "msisensorpro_germline",
    "structural_variants",
    "structural_variants_retiered",
    "structural_variants_unfiltered",
    "structural_variants_raw",
    "frag_cov_tumor",
    "frag_cov_normal",
    "coverage_tumor",
    "coverage_normal",
    "snvs_somatic",
    "snvs_somatic_unfiltered",
    "snvs_germline",
    "het_pileups",
    "amber_dir",
    "cobalt_dir",
    "purple_pp_best_fit",
    "purple_pp_range",
    "purity",
    "ploidy",
    "seg",
    "nseg",
    "variant_annotations_somatic",
    "variant_annotations_germline",
    "variant_annotations_somatic_vcf",
    "variant_annotations_germline_vcf",
    "karyograph",
    "jabba_rds",
    "jabba_gg",
    "jabba_gg_balanced",
    "jabba_gg_allelic",
    "events",
    "fusions",
    "fusions_junctions",
    "multiplicity",
    "multiplicity_germline",
    "multiplicity_hetsnps",
    "oncokb_snv",
    "oncokb_cna",
    "oncokb_fusions",
    "signatures_activities_sbs",
    "signatures_matrix_sbs",
    "signatures_decomposed_sbs",
    "signatures_post_prob_sbs",
    "signatures_activities_indel",
    "signatures_matrix_indel",
    "signatures_decomposed_indel",
    "signatures_post_prob_indel",
    "ffpe_impact_vcf",
  	"ffpe_impact_vcf_filtered",
    "hrdetect",
    "onenesstwoness",
    "conpair_concordance",
    "conpair_contamination"
]

# Define the default samplesheet columns
SAMPLESHEET_FIELDNAMES = [
    "patient", "sample", "status", "sex", 
    "bam", 
    "qc_dup_rate", "qc_alignment_summary", "qc_insert_size", "qc_coverage_metrics", 
    "msi", "msi_germline", 
    "hets", "amber_dir", 
    "frag_cov", "dryclean_cov",
    "cobalt_dir", "purple_pp_best_fit", "purple_pp_range", "purity", "ploidy", 
    "seg", "nseg", 
    "vcf", "vcf_raw",
    "structural_variants_retiered",
    "structural_variants_chimera_filtered", "structural_variants_raw_chimera_filtered",
    "jabba_rds", "jabba_gg", "ni_balanced_gg", "lp_balanced_gg", "events", "fusions",
    "snv_somatic_vcf", "snv_germline_vcf", 
    "variant_somatic_ann", "variant_somatic_bcf", 
    "variant_germline_ann", "variant_germline_bcf",
    "snv_multiplicity", 
    "oncokb_maf", "oncokb_fusions", "oncokb_cna",
    "sbs_signatures", "indel_signatures", "signatures_matrix", 
    "signatures_matrix_indel",
    "signatures_decomposed_sbs",
    "signatures_decomposed_indel",
    "hrdetect", "onenesstwoness", "conpair_concordance", "conpair_contamination"
]


OUTPUT_FILES_MAPPING_OLD = {
    "qc_dup_rate": r"qc_metrics/gatk/.*/.*metrics",
    "qc_alignment_summary": r"qc_metrics/picard/.*/.*alignment_summary_metrics",
    "qc_insert_size": r"qc_metrics/picard/.*/.*insert_size_metrics",
    "qc_coverage_metrics": r"qc_metrics/picard/.*/.*coverage_metrics",
    "qc_coverage_metrics_tumor": r"qc_metrics/picard/tumor/.*/.*coverage_metrics",
    "qc_coverage_metrics_normal": r"qc_metrics/picard/normal/.*/.*coverage_metrics",
    "msisensorpro": r"msisensorpro/.*(?<!_dis)(?<!_somatic)(?<!_germline)(?<!\.list)$",
    "structural_variants": [
        r"sv_calling/gridss_somatic/.*/.*high_confidence_somatic\.vcf\.bgz$",
        r"sv_calling/tumor_only_junction_filter/.*/somatic\.filtered\.sv\.rds$",
    ],
    "structural_variants_unfiltered": r"sv_calling/gridss/.*/.*\.gridss\.filtered\.vcf\.gz$",
    "frag_cov_tumor": r"coverage/fragcounter_tumor/.*/cov\.rds$",
    "frag_cov_normal": r"coverage/fragcounter_normal/.*/cov\.rds$",
    "coverage_tumor": r"coverage/dryclean_tumor/.*/drycleaned\.cov\.rds$",
    "coverage_normal": r"coverage/dryclean_normal/.*/drycleaned\.cov\.rds$",
    "snvs_somatic": [
        r"snv_calling/sage/somatic/tumor_only_filter/.*/.*\.sage\.pass_filtered\.tumoronly\.vcf\.gz$",
        r"snv_calling/sage/somatic/.*/.*\.sage\.pass_filtered\.vcf\.gz$"
    ],
    "snvs_somatic_unfiltered": [
        r"snv_calling/sage/somatic/.*/.*sage\.somatic\.vcf\.gz$",
        r"snv_calling/sage/somatic/.*/.*\.sage\.pass_filtered\.vcf\.gz$"
	],
    "snvs_germline": r"snv_calling/sage/germline/.*/.*sage\.germline\.vcf\.gz$",
    "het_pileups": [
        r"amber/.*/sites\.txt$",
        r"hetpileups/.*/sites\.txt$",
    ],
    "amber_dir": r"amber/.*/amber/",
    "cobalt_dir": r"cobalt/.*/cobalt/",
    "purple_pp_range": r"purple/.*/.*purple\.purity\.range\.tsv$",
    "purple_pp_best_fit": r"purple/.*/.*purple\.purity\.tsv$",
    "seg": r"cbs/.*/seg.rds",
    "nseg": r"cbs/.*/nseg.rds",
    "multiplicity": r"snv_multiplicity/.*/.*est_snv_cn_somatic\.rds$",
    "multiplicity_germline": r"snv_multiplicity/.*/.*est_snv_cn_germline\.rds$",
    "multiplicity_hetsnps": r"snv_multiplicity/.*/.*est_snv_cn_hets\.rds$",
    "variant_annotations_somatic": r"variant_annotations/snpeff/somatic/.*/.*ann\.bcf$",
    "variant_annotations_germline": r"variant_annotations/snpeff/germline/.*/.*ann\.bcf$",
    "variant_annotations_somatic_vcf": r"variant_annotations/snpeff/somatic/.*/.*ann\.vcf$",
    "variant_annotations_germline_vcf": r"variant_annotations/snpeff/germline/.*/.*ann\.vcf$",
    "oncokb_snv": r"oncokb/.*/merged_oncokb\.maf$",
    "oncokb_cna": r"oncokb/.*/merged_oncokb_cna\.tsv$",
    "oncokb_fusions": r"oncokb/.*/merged_oncokb_fusions\.tsv$",
    "karyograph": r"jabba/.*/karyograph\.rds$",
    "jabba_rds": r"jabba/.*/jabba\.simple\.rds$",
    "jabba_gg": r"jabba/.*/jabba\.simple\.gg\.rds$",
    "jabba_gg_balanced": r"non_integer_balance/.*/non_integer\.balanced\.gg\.rds$",
    "jabba_gg_allelic": r"lp_phased_balance/.*/lp_phased\.balanced\.gg\.rds$",
    "events": r"events/.*/complex\.rds$",
    "fusions": r"fusions/.*/fusions\.rds$",
    "signatures_activities_sbs": r"signatures/sigprofilerassignment/somatic/.*/sbs_results/Assignment_Solution/Activities/sbs_Assignment_Solution_Activities\.txt",
    "signatures_matrix_sbs": r"signatures/sigprofilerassignment/somatic/.*/sig_inputs/output/SBS/sigmat_results\.SBS96\.all",
    "signatures_decomposed_sbs": r"signatures/sigprofilerassignment/somatic/.*/sbs_results/Assignment_Solution/Activities/Decomposed_MutationType_Probabilities\.txt",
    "signatures_activities_indel": r"signatures/sigprofilerassignment/somatic/.*/indel_results/Assignment_Solution/Activities/indel_Assignment_Solution_Activities\.txt",
    "signatures_matrix_indel": r"signatures/sigprofilerassignment/somatic/.*/sig_inputs/output/ID/sigmat_results\.ID83\.all",
    "signatures_decomposed_indel": r"signatures/sigprofilerassignment/somatic/.*/indel_results/.*/Decomposed_MutationType_Probabilities\.txt",
    "hrdetect": r"hrdetect/.*/hrdetect_results\.rds",
    "onenesstwoness": r"oneness/.*/oneness_twoness\.rds$",
}

# Map each output key to its file regex pattern(s)
OUTPUT_FILES_MAPPING = {
    "bam": [
        r"alignment/.*\.bam$"	
	],
    "bam_tumor": [
        r"alignment/.*\.bam$"	
	],
    "bam_normal": [
        r"alignment/.*\.bam$"	
	],
    "qc_dup_rate": [
        r"gatk_qc/.*/.*metrics",
        r"gatk_qc/.*metrics",
        r"markduplicates/.*metrics",
        r"alignment/.*duplicate-metrics.txt",
        r"parabricks/.*duplicate-metrics.txt"
	],
    "qc_dup_rate_tumor": [
        r"gatk_qc/.*/.*metrics",
        r"gatk_qc/.*metrics",
        r"markduplicates/.*metrics",
        r"alignment/tumor/.*duplicate-metrics.txt",
        r"parabricks/.*duplicate-metrics.txt"
	],
    "qc_dup_rate_normal": [
        r"gatk_qc/.*/.*metrics",
        r"gatk_qc/.*metrics",
        r"markduplicates/.*metrics",
        r"alignment/normal/.*duplicate-metrics.txt",
        r"parabricks/.*duplicate-metrics.txt"
	],
    "qc_alignment_summary": [
        r"picard_qc/.*/.*alignment_summary_metrics",
        r"picard_qc/.*alignment_summary_metrics",
        r"alignment/.*qc_metrics/alignment.txt",
        r"parabricks/.*alignment.txt"
	],
    "qc_alignment_summary_tumor": [
        r"picard_qc/.*/.*alignment_summary_metrics",
        r"picard_qc/.*alignment_summary_metrics",
        r"alignment/tumor/.*qc_metrics/alignment.txt"
        r"parabricks/.*alignment.txt"
	],
    "qc_alignment_summary_normal": [
        r"picard_qc/.*/.*alignment_summary_metrics",
        r"picard_qc/.*alignment_summary_metrics",
        r"alignment/normal/.*qc_metrics/alignment.txt",
        r"parabricks/.*alignment.txt"
	],
    "qc_insert_size": [
        r"picard_qc/.*insert_size_metrics",
        r"picard_qc/.*/.*insert_size_metrics",
        r"alignment/.*qc_metrics/insert_size.txt",
        r"parabricks/.*insert_size.txt"
	],
    "qc_insert_size_tumor": [
        r"picard_qc/.*insert_size_metrics",
        r"alignment/tumor/.*qc_metrics/insert_size.txt",
        r"parabricks/.*insert_size.txt"
	],
    "qc_insert_size_normal": [
        r"picard_qc/.*insert_size_metrics",
        r"alignment/normal/.*qc_metrics/insert_size.txt",
        r"parabricks/.*insert_size.txt"
	],
    "qc_coverage_metrics": [
        r"picard_qc/.*coverage_metrics",
        r"picard_qc/.*/.*coverage_metrics",
        r"parabricks_qc/.*coverage_metrics",
        r"parabricks/.*coverage_metrics"
	],
    "qc_coverage_metrics_tumor": [
        r"picard_qc/.*coverage_metrics",
        r"picard_qc/tumor/.*/.*coverage_metrics",
        r"parabricks_qc/tumor/.*coverage_metrics",
        r"parabricks/.*coverage_metrics"
	],
    "qc_coverage_metrics_normal": [
        r"picard_qc/.*coverage_metrics",
        r"picard_qc/normal/.*/.*coverage_metrics",
        r"parabricks_qc/normal/.*coverage_metrics",
        r"parabricks/.*coverage_metrics"
	],
    "msisensorpro": r"msisensorpro/.*_report$",
    "structural_variants": [
        r"gridss/.*high_confidence_somatic\.vcf\.bgz$",
        r"gridss/.*somatic\.filtered\.sv\.no\.gnomAD\.rds$",
        r"tumor_only_junction_filter/.*/somatic\.filtered\.sv\.rds$",
        r"gridss/.*somatic\.filtered\.sv\.rds$",
    ],
    "structural_variants_retiered": [
        r"jabba/.*somatic\.filtered\.sv.*___tiered\.rds$"
    ],
    "structural_variants_unfiltered": r"gridss.*/.*\.gridss\.filtered\.vcf\.gz$",
    "structural_variants_raw": r"gridss.*/.*\.gridss\.vcf\.gz$",
    "structural_variants_chimera_filtered": r"sv_chimera_filter.*/.*\.filtered.ffpe_filtered\.vcf\.gz$",
    "structural_variants_raw_chimera_filtered": r"sv_chimera_filter.*/(?!.*\.filtered\.).*\.ffpe_filtered\.vcf\.gz$",
    "frag_cov_tumor": r"fragcounter/tumor/cov\.rds$",
    "frag_cov_normal": r"fragcounter/normal/cov\.rds$",
    "coverage_tumor": r"dryclean/tumor/drycleaned\.cov\.rds$",
    "coverage_normal": r"dryclean/normal/drycleaned\.cov\.rds$",
    "snvs_somatic": [
        r"sage/somatic/tumor_only_filter/.*\.sage\.pass_filtered\.tumoronly\.vcf\.gz$",
        r"sage/somatic/.*\.sage\.pass_filtered\.vcf\.gz$"
    ],
    "snvs_somatic_unfiltered": [
        r"sage/somatic/.*sage\.somatic\.vcf\.gz$",
        r"sage/somatic/.*\.sage\.pass_filtered\.vcf\.gz$"
	],
    "snvs_germline": r"sage/germline/.*sage\.germline\.vcf\.gz$",
    "het_pileups": r"amber/sites\.txt$",
    "amber_dir": r"amber/amber/",
    "cobalt_dir": r"cobalt/cobalt/",
    "purple_pp_range": r"purple/.*purple\.purity\.range\.tsv$",
    "purple_pp_best_fit": r"purple/.*purple\.purity\.tsv$",
    "purple_qc": r"purple/.*purple\.qc$",
    "seg": r"cbs/seg.rds",
    "nseg": r"cbs/nseg.rds",
    "multiplicity": r"snv_multiplicity/.*est_snv_cn_somatic\.rds$",
    "multiplicity_germline": r"snv_multiplicity/.*est_snv_cn_germline\.rds$",
    "multiplicity_hetsnps": r"snv_multiplicity/.*est_snv_cn_hets\.rds$",
    "variant_annotations_somatic": r"snpeff/somatic/.*ann\.bcf$",
    "variant_annotations_germline": r"snpeff/germline/.*ann\.bcf$",
    "variant_annotations_somatic_vcf": r"snpeff/somatic/.*ann\.vcf$",
    "variant_annotations_germline_vcf": r"snpeff/germline/.*ann\.vcf$",
    "oncokb_snv": r"oncokb/merged_oncokb\.maf$",
    "oncokb_cna": r"oncokb/merged_oncokb_cna\.tsv$",
    "oncokb_fusions": r"oncokb/merged_oncokb_fusions\.tsv$",
    "karyograph": r"jabba/karyograph\.rds$",
    "jabba_rds": r"jabba/jabba\.simple\.rds$",
    "jabba_gg": r"jabba/jabba\.simple\.gg\.rds$",
    "jabba_gg_balanced": r"non_integer_balance/non_integer\.balanced\.gg\.rds$",
    "jabba_gg_allelic": r"lp_phased_balance/lp_phased\.balanced\.gg\.rds$",
    "events": r"events/complex\.rds$",
    "fusions": r"fusions/fusions\.rds$",
    "fusions_junctions": r"fusions/altedge\.annotations\.tsv$",
    "signatures_activities_sbs": r"sigprofilerassignment/sbs_results/Assignment_Solution/Activities/sbs_Assignment_Solution_Activities\.txt",
    "signatures_post_prob_sbs": r"sigprofilerassignment/sbs_results/Assignment_Solution/Activities/Decomposed_Mutation_Probabilities/Decomposed.*\.txt",
    "signatures_matrix_sbs": r"sigprofilerassignment/sig_inputs/output/SBS/sigmat_results\.SBS96\.all",
    "signatures_decomposed_sbs": r"sigprofilerassignment/sbs_results/Assignment_Solution/Activities/Decomposed_MutationType_Probabilities\.txt",
    "signatures_activities_indel": r"sigprofilerassignment/indel_results/Assignment_Solution/Activities/indel_Assignment_Solution_Activities\.txt",
    "signatures_post_prob_indel": r"sigprofilerassignment/indel_results/Assignment_Solution/Activities/Decomposed_Mutation_Probabilities/Decomposed.*\.txt",
    "signatures_decomposed_indel": r"sigprofilerassignment/indel_results/Assignment_Solution/Activities/Decomposed_Mutation_Probabilities/Decomposed.*\.txt",
    "ffpe_impact_vcf": r"ffpe_impact/.*ffpe_annotated\.vcf\.gz$",
  	"ffpe_impact_vcf_filtered": r"ffpe_impact/.*ffpe_annotated_filtered\.vcf\.gz$",
    "signatures_matrix_indel": r"sigprofilerassignment/sig_inputs/output/ID/sigmat_results\.ID83\.all",
    "signatures_decomposed_indel": r"sigprofilerassignment/indel_results/.*/Decomposed_MutationType_Probabilities\.txt",
    "hrdetect": r"hrdetect/hrdetect_results\.rds",
    "onenesstwoness": r"onenesstwoness/onenesstwoness_results\.rds$",
    "conpair_concordance": r"conpair/concordance\.txt$",
    "conpair_contamination": r"conpair/contamination\.txt$"
}

class Outputs:
    def __init__(self, outputs_dir: str, samplesheet: str, use_old: bool = False, prefer_outputs: bool = False):
        self.outputs_dir = os.path.abspath(outputs_dir)
        self.samplesheet = samplesheet
        self.samples_data = self._read_samplesheet()
        self.outputs = self._collect_outputs(use_old_output_files_mapping=use_old, prefer_outputs = prefer_outputs)

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

            # Maps for columns that need to be split by status (tumor "1" vs normal "0")
            conditional_mapping = {
                "bam": ("bam_tumor", "bam_normal"),
                "qc_dup_rate": ("qc_dup_rate", "qc_dup_rate_normal"),
                "qc_alignment_summary": ("qc_alignment_summary", "qc_alignment_summary_normal"),
                "qc_insert_size": ("qc_insert_size", "qc_insert_size_normal"),
                "qc_coverage_metrics": ("qc_coverage_metrics", "qc_coverage_metrics_normal"),
                "frag_cov": ("frag_cov_tumor", "frag_cov_normal"),
                "dryclean_cov": ("coverage_tumor", "coverage_normal"),
            }

            # Maps for columns that map directly (samplesheet column -> Outputs key)
            direct_mapping = {
                "sex": "sex",
                "msi": "msisensorpro",
                "msi_germline": "msisensorpro_germline",
                "hets": "het_pileups",
                "amber_dir": "amber_dir",
                "cobalt_dir": "cobalt_dir",
                "purity": "purity",
                "ploidy": "ploidy",
                "seg": "seg",
                "nseg": "nseg",
                "vcf": "structural_variants",
                "vcf_raw": "structural_variants_raw",
                "structural_variants_retiered": "structural_variants_retiered",
                "structural_variants_chimera_filtered": "structural_variants_chimera_filtered",
                "structural_variants_raw_chimera_filtered": "structural_variants_raw_chimera_filtered",
                "jabba_rds": "jabba_rds",
                "jabba_gg": "jabba_gg",
                "ni_balanced_gg": "jabba_gg_balanced",
                "lp_balanced_gg": "jabba_gg_allelic",
                "purple_pp_best_fit": "purple_pp_best_fit",
                "purple_pp_range": "purple_pp_range",
                "events": "events",
                "fusions": "fusions",
                "snv_somatic_vcf": "snvs_somatic",
                "snv_germline_vcf": "snvs_germline",
                "variant_somatic_ann": "variant_annotations_somatic_vcf",
                "variant_somatic_bcf": "variant_annotations_somatic",
                "variant_germline_ann": "variant_annotations_germline_vcf",
                "variant_germline_bcf": "variant_annotations_germline",
                "snv_multiplicity": "multiplicity",
                "oncokb_maf": "oncokb_snv",
                "oncokb_fusions": "oncokb_fusions",
                "oncokb_cna": "oncokb_cna",
                "sbs_signatures": "signatures_activities_sbs",
                "indel_signatures": "signatures_activities_indel",
                "signatures_matrix": "signatures_matrix_sbs",
                "signatures_matrix_indel": "signatures_matrix_indel",
                "signatures_decomposed_sbs": "signatures_decomposed_sbs",
                "signatures_decomposed_indel": "signatures_decomposed_indel",
                "hrdetect": "hrdetect",
                "onenesstwoness": "onenesstwoness",
                "conpair_concordance": "conpair_concordance",
                "conpair_contamination": "conpair_contamination"
            }

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
                    }

                status = row.get("status", "").strip()
                sample_id = row.get("sample", "").strip()
                if sample_id and sample_id not in patient_data[patient_id]["sample_ids"]:
                    if status == "1" and patient_data[patient_id]["sample_ids"]:
                        patient_data[patient_id]["sample_ids"].insert(0, sample_id)
                    else:
                        patient_data[patient_id]["sample_ids"].append(sample_id)

                bam = row.get("bam", "").strip()
                if bam and status:
                    if status == "1":  # Tumor sample
                        patient_data[patient_id]["bam_tumor"] = bam
                    elif status == "0":  # Normal sample
                        patient_data[patient_id]["bam_normal"] = bam

                # Process all other columns from the samplesheet row
                for col, value in row.items():
                    value = value.strip() if value else ""
                    if not value or col in ["patient", "sample", "status", "bam"]:
                        continue

                    if col in conditional_mapping:
                        # Assign to tumor or normal key based on status
                        out_key = conditional_mapping[col][0] if status == "1" else (
                                  conditional_mapping[col][1] if status == "0" else None)
                        if out_key and not patient_data[patient_id].get(out_key):
                            patient_data[patient_id][out_key] = value

                    elif col in direct_mapping:
                        out_key = direct_mapping[col]
                        if not patient_data[patient_id].get(out_key):
                            patient_data[patient_id][out_key] = value

                    else:
                        # For any other columns not explicitly mapped, store them using the original column name.
                        if col not in patient_data[patient_id]:
                            patient_data[patient_id][col] = value

        return patient_data

    def _apply_old_mapping(self, record: dict) -> None:
        """Apply the old output files mapping to the record."""
        mapping = OUTPUT_FILES_MAPPING_OLD
        for key, pattern in mapping.items():
            if record.get(key):
                continue  # prefer samplesheet value if available
            patterns = pattern if isinstance(pattern, list) else [pattern]
            for pat in patterns:
                # Derive process directory prefix from the pattern (assume the prefix is the literal part before '/.*/')
                if '/.*/' in pat:
                    process_prefix = pat.split('/.*/')[0]
                else:
                    process_prefix = os.path.dirname(pat)
                # Build the search directory using the process prefix
                search_dir = os.path.join(self.outputs_dir, process_prefix)
                if key == "msisensorpro":
                    # Special case for msisensorpro, which has patient in the filename
                    search_dir = self.outputs_dir
                search_pattern = os.path.join(search_dir, "**", "*")

                for root, dirs, files in os.walk(search_dir):
                    # Remove 'work' from dirs to prevent os.walk from traversing it
                    if 'work' in dirs and os.path.samefile(os.path.join(root, 'work'), 
                                                           os.path.join(search_dir, 'work')):
                        dirs.remove('work')

                    # Process files in the current directory
                    for file in files:
                        filepath = os.path.join(root, file)

                        is_matching_sample = False
                        for sample_id in record.get("sample_ids", []):
                            if sample_id in filepath:
                                is_matching_sample = True
                                break

                        if is_matching_sample:
                            if re.search(pat, filepath):
                                # record absolute path
                                record[key] = filepath
                                if pat.endswith("/"):
                                    record[key] = os.path.dirname(filepath)
                                break
                    if record.get(key):
                        break

    def _collect_outputs(
        self,
        use_old_output_files_mapping = False,
        prefer_outputs = False
    ) -> list:
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
            # use the sample id for old mapping
            record["sample_ids"] = data.get("sample_ids", [])
            ## FIXME: Assuming that sample ids is a length 2 list,
            ## where tumor id always goes first (encoded by _read_samplesheet)
            sample_ids = data.get("sample_ids", [])
            ln = len(sample_ids)
            is_paired = ln == 2
            is_empty = ln == 0
            is_invalid = ln > 2
            
            if is_invalid:
                print(sampleids)
                raise ValueError("More than one sample id found - not supported yet")
            
            if not is_empty:
                record["tumor_sample"] = sample_ids[0]
            if is_paired:
                record["normal_sample"] = sample_ids[1]

            # Overwrite with top-level keys where provided
            for key in data.keys():
                record[key] = data[key]

            if use_old_output_files_mapping:
                self._apply_old_mapping(record)
            else:
                mapping = OUTPUT_FILES_MAPPING
                for key, pattern in mapping.items():
                    is_key_tumor = "_tumor" in key
                    is_key_normal = "_normal" in key
                    is_key_neither_tumor_normal = not is_key_tumor and not is_key_normal
                    # if record.get(key):
                    if not prefer_outputs and record.get(key):
                        continue  # prefer samplesheet value if available
                    patterns = pattern if isinstance(pattern, list) else [pattern]
                    patient_dir = os.path.join(self.outputs_dir, patient_id)
                    is_pattern_filepath_matched = False ## Initializing break conditional
                    ## Pattern finding
                    for pat in patterns:
                        search_pattern = os.path.join(patient_dir, "**", "*")
                        for filepath in glob.glob(search_pattern, recursive=True):
                            rel_path = os.path.relpath(filepath, patient_dir)
                            ## FIXME: hack logic is that if sample_ids is not present, return length 2 list
                            ## with random, unmatchable strings
                            data_sample_ids = data.get(
                                "sample_ids", 
                                [
                                    "VBHWHNhrLQwyX56NDOBoMWBO", 
                                    "dT1Z99GJSU1XT95v1vARKdOt"
                                ]
                            )
                            is_tumor_sample_id_in_path = bool(re.search(data_sample_ids[0], rel_path))
                            is_normal_sample_id_in_path = False
                            if len(data_sample_ids) > 1:
                            	is_normal_sample_id_in_path = bool(re.search(data_sample_ids[1], rel_path))
                            # is_sample_id_in_path = is_tumor_sample_id_in_path or is_normal_sample_id_in_path
                            is_sample_id_in_path = any([bool(re.search(sample_id, rel_path)) for sample_id in data_sample_ids])
                            is_sample_id_absent_in_path = not is_sample_id_in_path
                            is_pattern_present = bool(re.search(pat, rel_path))
                            ## Pattern is in file path, sample_id isn't in file path
                            is_proceed_with_first_file_match = is_pattern_present and is_sample_id_absent_in_path
                            ## Pattern is in file path, sample_id is in file path, and key matches with the sample_id type (tumor vs normal)
                            is_sample_id_file_matched = is_pattern_present and is_sample_id_in_path
                            is_proceed_with_tumor_sample_id_file_match = is_sample_id_file_matched and is_tumor_sample_id_in_path and is_key_tumor
                            is_proceed_with_normal_sample_id_file_match = is_sample_id_file_matched and is_normal_sample_id_in_path and is_key_normal
                            ## Pattern is in file path, sample_id is in file path and is the tumor, but the column is neither tumor/normal specific
                            ## The below will preferentially populate with the tumor (which is what we want in most cases)
                            is_proceed_with_tumor_file_match = is_sample_id_file_matched and is_key_neither_tumor_normal and is_tumor_sample_id_in_path
                            is_filepath_to_be_populated = (
                                is_proceed_with_first_file_match
                                or is_proceed_with_tumor_sample_id_file_match
                                or is_proceed_with_normal_sample_id_file_match
                                or is_proceed_with_tumor_file_match
							)
                            if is_filepath_to_be_populated:
                                record[key] = filepath
                                if pat.endswith("/"):
                                    record[key] = os.path.dirname(filepath)
                            is_key_populated = bool(record.get(key)) 
                            if is_key_populated:
                                is_pattern_filepath_matched = True
                                break
                        if is_pattern_filepath_matched:
                            break

            # New: Populate purity and ploidy from purple.purity.tsv (if available)
            purity_file = record["purple_pp_best_fit"]
            if purity_file:
                with open(purity_file) as pf:
                    lines = pf.read().splitlines()
                    if len(lines) >= 2:
                        headers = lines[0].split("\t")
                        values = lines[1].split("\t")
                        mapping_dict = dict(zip(headers, values))
                        if "purity" in mapping_dict:
                            record["purity"] = mapping_dict["purity"]
                        if "ploidy" in mapping_dict:
                            record["ploidy"] = mapping_dict["ploidy"]

            outputs_list.append(record)
        return outputs_list

    def emit_output_csv(self, csv_path: Optional[str] = None,
                        include_columns: Optional[List[str]] = None,
                        exclude_columns: Optional[List[str]] = None):
        """
        Write the collected outputs (self.outputs) to a CSV file at csv_path.
        The CSV includes keys from OUTPUT_KEYS except 'sample_ids'.
        Columns can be filtered using include_columns or exclude_columns.
        For any missing values, output an empty string.
        """
        base_fieldnames = [key for key in OUTPUT_KEYS if key != "sample_ids"]

        if include_columns:
            # Filter base_fieldnames, maintaining the order of include_columns if they exist in base
            included_set = set(include_columns)
            fieldnames = [col for col in base_fieldnames if col in included_set]
            # Ensure the order matches include_columns for columns present in base_fieldnames
            fieldnames.sort(key=lambda col: include_columns.index(col) if col in include_columns else float('inf'))
        elif exclude_columns:
            excluded_set = set(exclude_columns)
            fieldnames = [col for col in base_fieldnames if col not in excluded_set]
        else:
            fieldnames = base_fieldnames

        if not fieldnames:
            print("Error: No columns selected for output.", file=sys.stderr)
            return

        # remove empty columns
        fieldnames = [col for col in fieldnames if any(row.get(col) for row in self.outputs)]
        if not fieldnames:
            print("Error: No columns selected for output.", file=sys.stderr)
            return

        if csv_path:
            output_stream = open(csv_path, "w", newline="")
        else:
            output_stream = sys.stdout
        writer = csv.DictWriter(output_stream, fieldnames=fieldnames)
        writer.writeheader()
        for row in self.outputs:
            writer.writerow({key: row.get(key, "") for key in fieldnames})
        if csv_path:
            output_stream.close()

    def emit_samplesheet_csv(self, csv_path: Optional[str] = None,
                             include_columns: Optional[List[str]] = None,
                             exclude_columns: Optional[List[str]] = None):
        """
        Write a tall samplesheet CSV where every row corresponds to one sample.
        Columns can be filtered using include_columns or exclude_columns.
        For paired samples, a row is written for the tumor (status "1", bam, frag_cov, dryclean_cov from bam_tumor, frag_cov_tumor, coverage_tumor)
        and another for the normal (status "0", bam, frag_cov, dryclean_cov from bam_normal, frag_cov_normal, coverage_normal).
        Other fields are copied directly from the record (they come from the original samplesheet or outputs).
        The mapping from the output record keys to the samplesheet columns is defined as follows:

            samplesheet_col              -> outputs key (or conditional):
            patient                      -> patient_id
            sample                       -> sample_ids (single value)
            status                       -> "1" for tumor, "0" for normal
            sex                          -> sex
            bam                          -> bam_tumor (if tumor) or bam_normal (if normal)
            qc_dup_rate                  -> GATK EstimateLibraryComplexity
            qc_dup_rate_tumor            -> GATK EstimateLibraryComplexity from tumor
            qc_dup_rate_normal           -> GATK EstimateLibraryComplexity from normal
            qc_insert_size               -> Picard CollectMultipleMetrics 
            qc_insert_size_tumor         -> Picard CollectMultipleMetrics from tumor 
            qc_insert_size_normal        -> Picard CollectMultipleMetrics from normal 
            qc_alignment_summary         -> Picard CollectMultipleMetrics
            qc_alignment_summary_tumor   -> Picard CollectMultipleMetrics from tumor
            qc_alignment_summary_normal  -> Picard CollectMultipleMetrics from normal
            qc_coverage_metrics          -> Picard CollectWGSMetrics
            qc_coverage_metrics_tumor    -> Picard CollectWGSMetrics from tumor
            qc_coverage_metrics_normal   -> Picard CollectWGSMetrics from normal
            msi                          -> msisensorpro
            hets                         -> het_pileups
            amber_dir                    -> amber_dir
            frag_cov                     -> frag_cov_tumor (if tumor) or frag_cov_normal (if normal)
            dryclean_cov                 -> coverage_tumor (if tumor) or coverage_normal (if normal)
            cobalt_dir                   -> cobalt_dir
            purity                       -> purity
            ploidy                       -> ploidy
            seg                          -> seg
            nseg                         -> nseg
            vcf                          -> structural_variants
            jabba_rds                    -> jabba_rds
            jabba_gg                     -> jabba_gg
            ni_balanced_gg               -> jabba_gg_balanced
            lp_balanced_gg               -> jabba_gg_allelic
            events                       -> events
            fusions                      -> fusions
            snv_somatic_vcf              -> snvs_somatic
            snv_germline_vcf             -> snvs_germline
            variant_somatic_ann          -> variant_annotations_somatic_vcf
            variant_somatic_bcf          -> variant_annotations_somatic
            variant_germline_ann         -> variant_annotations_germline_vcf
            variant_germline_bcf         -> variant_annotations_germline
            snv_multiplicity             -> multiplicity
            oncokb_maf                   -> oncokb_snv
            oncokb_fusions               -> oncokb_fusions
            oncokb_cna                   -> oncokb_cna
            sbs_signatures               -> signatures_activities_sbs
            indel_signatures             -> signatures_activities_indel
            signatures_matrix            -> signatures_matrix_sbs
            hrdetect                     -> hrdetect
            onenesstwoness               -> onenesstwoness

        Rows are generated by inspecting each record in self.outputs.
        If both bam_tumor and bam_normal are present, two rows are emitted
        (the first sample in sample_ids corresponds to tumor and the second to normal).
        If only one BAM is available, only a single row is emitted.
        """
        base_fieldnames = SAMPLESHEET_FIELDNAMES

        if include_columns:
            # Filter base_fieldnames, maintaining the order of include_columns if they exist in base
            included_set = set(include_columns)
            fieldnames = [col for col in base_fieldnames if col in included_set]
             # Ensure the order matches include_columns for columns present in base_fieldnames
            fieldnames.sort(key=lambda col: include_columns.index(col) if col in include_columns else float('inf'))
        elif exclude_columns:
            excluded_set = set(exclude_columns)
            fieldnames = [col for col in base_fieldnames if col not in excluded_set]
        else:
            fieldnames = base_fieldnames

        if not fieldnames:
            print("Error: No columns selected for output.", file=sys.stderr)
            return

        if csv_path:
            output_stream = open(csv_path, "w", newline="")
        else:
            output_stream = sys.stdout
        writer = csv.DictWriter(output_stream, fieldnames=fieldnames)
        writer.writeheader()

        for record in self.outputs:
            sample_rows = []
            if len(record.get("sample_ids", [])) >= 2:
                print(f"Found paired sample for patient {record['patient_id']}, {record['sample_ids']}")

                # Tumor row (status "1")
                tumor_sample = record["sample_ids"][0] if record["sample_ids"] else ""
                tumor_row = {
                    "patient": record.get("patient_id", ""),
                    "sample": tumor_sample,
                    "status": "1",
                    "sex": record.get("sex", ""),
                    "bam": record.get("bam_tumor", ""),
                    "qc_dup_rate": record.get("qc_dup_rate", ""),
                    "qc_alignment_summary": record.get("qc_alignment_summary", ""),
                    "qc_insert_size": record.get("qc_insert_size", ""),
                    "qc_coverage_metrics": record.get("qc_coverage_metrics", ""),
                    "msi": record.get("msisensorpro", ""),
                    "msi_germline": record.get("msisensorpro_germline", ""),
                    "hets": record.get("het_pileups", ""),
                    "amber_dir": record.get("amber_dir", ""),
                    "frag_cov": record.get("frag_cov_tumor", ""),
                    "dryclean_cov": record.get("coverage_tumor", ""),
                    "cobalt_dir": record.get("cobalt_dir", ""),
                    "purple_pp_best_fit": record.get("purple_pp_best_fit", ""),
                    "purple_pp_range": record.get("purple_pp_range", ""),
                    "purity": record.get("purity", ""),
                    "ploidy": record.get("ploidy", ""),
                    "seg": record.get("seg", ""),
                    "nseg": record.get("nseg", ""),
                    "vcf": record.get("structural_variants", ""),
                    "vcf_raw": record.get("structural_variants_raw", ""),
                    "structural_variants_retiered": record.get("structural_variants_retiered", ""),
                    "structural_variants_chimera_filtered": record.get("structural_variants_chimera_filtered", ""),
                    "structural_variants_raw_chimera_filtered": record.get("structural_variants_raw_chimera_filtered", ""),
                    "jabba_rds": record.get("jabba_rds", ""),
                    "jabba_gg": record.get("jabba_gg", ""),
                    "ni_balanced_gg": record.get("jabba_gg_balanced", ""),
                    "lp_balanced_gg": record.get("jabba_gg_allelic", ""),
                    "events": record.get("events", ""),
                    "fusions": record.get("fusions", ""),
                    "snv_somatic_vcf": record.get("snvs_somatic", ""),
                    "snv_germline_vcf": record.get("snvs_germline", ""),
                    "variant_somatic_ann": record.get("variant_annotations_somatic_vcf", ""),
                    "variant_somatic_bcf": record.get("variant_annotations_somatic", ""),
                    "variant_germline_ann": record.get("variant_annotations_germline_vcf", ""),
                    "variant_germline_bcf": record.get("variant_annotations_germline", ""),
                    "snv_multiplicity": record.get("multiplicity", ""),
                    "oncokb_maf": record.get("oncokb_snv", ""),
                    "oncokb_fusions": record.get("oncokb_fusions", ""),
                    "oncokb_cna": record.get("oncokb_cna", ""),
                    "sbs_signatures": record.get("signatures_activities_sbs", ""),
                    "indel_signatures": record.get("signatures_activities_indel", ""),
                    "signatures_matrix": record.get("signatures_matrix_sbs", ""),
                    "signatures_matrix_indel": record.get("signatures_matrix_indel", ""),
                    "signatures_decomposed_sbs": record.get("signatures_decomposed_sbs", ""),
                    "signatures_decomposed_indel": record.get("signatures_decomposed_indel", ""),
                    "hrdetect": record.get("hrdetect", ""),
                    "onenesstwoness": record.get("onenesstwoness", ""),
                    "conpair_concordance": record.get("conpair_concordance", ""),
                    "conpair_contamination": record.get("conpair_contamination", "")
                }
                sample_rows.append(tumor_row)

                # Normal row (status "0"): second sample in sample_ids if available
                normal_sample = record["sample_ids"][1] if len(record["sample_ids"]) >= 2 else ""
                normal_row = tumor_row.copy()
                normal_row["sample"] = normal_sample
                normal_row["status"] = "0"
                normal_row["bam"] = record.get("bam_normal", "")
                normal_row["qc_dup_rate"] = record.get("qc_dup_rate_normal", "")
                normal_row["qc_alignment_summary"] = record.get("qc_alignment_summary_normal", "")
                normal_row["qc_insert_size"] = record.get("qc_insert_size_normal", "")
                normal_row["qc_coverage_metrics"] = record.get("qc_coverage_metrics_normal", "")
                normal_row["frag_cov"] = record.get("frag_cov_normal", "")
                normal_row["dryclean_cov"] = record.get("coverage_normal", "")
                sample_rows.append(normal_row)
            else:
                status = "1"
                bam_val = record.get("bam_tumor", "")
                qc_dup_val = record.get("qc_dup_rate_tumor", "")
                qc_insert_val = record.get("qc_insert_size_tumor", "")
                qc_alignment_val = record.get("qc_alignment_summary_tumor", "")
                qc_coverage_val = record.get("qc_coverage_metrics_tumor", "")
                frag_cov = record.get("frag_cov_tumor", "")
                dryclean_cov = record.get("coverage_tumor", "")

                sample_val = record["sample_ids"][0] if record["sample_ids"] else ""
                sample_rows.append({
                    "patient": record.get("patient_id", ""),
                    "sample": sample_val,
                    "status": status,
                    "sex": record.get("sex", ""),
                    "bam": bam_val,
                    "qc_dup_rate": record.get("qc_dup_rate", ""),
                    "qc_alignment_summary": record.get("qc_alignment_summary", ""),
                    "qc_insert_size": record.get("qc_insert_size", ""),
                    "qc_coverage_metrics": record.get("qc_coverage_metrics", ""),
                    "msi": record.get("msisensorpro", ""),
                    "hets": record.get("het_pileups", ""),
                    "amber_dir": record.get("amber_dir", ""),
                    "frag_cov": frag_cov,
                    "dryclean_cov": dryclean_cov,
                    "cobalt_dir": record.get("cobalt_dir", ""),
                    "purple_pp_best_fit": record.get("purple_pp_best_fit", ""),
                    "purple_pp_range": record.get("purple_pp_range", ""),
                    "purity": record.get("purity", ""),
                    "ploidy": record.get("ploidy", ""),
                    "seg": record.get("seg", ""),
                    "nseg": record.get("nseg", ""),
                    "vcf": record.get("structural_variants", ""),
                    "vcf_raw": record.get("structural_variants_raw", ""),
                    "structural_variants_retiered": record.get("structural_variants_retiered", ""),
                    "structural_variants_chimera_filtered": record.get("structural_variants_chimera_filtered", ""),
                    "structural_variants_raw_chimera_filtered": record.get("structural_variants_raw_chimera_filtered", ""),
                    "jabba_rds": record.get("jabba_rds", ""),
                    "jabba_gg": record.get("jabba_gg", ""),
                    "ni_balanced_gg": record.get("jabba_gg_balanced", ""),
                    "lp_balanced_gg": record.get("jabba_gg_allelic", ""),
                    "events": record.get("events", ""),
                    "fusions": record.get("fusions", ""),
                    "snv_somatic_vcf": record.get("snvs_somatic", ""),
                    "snv_germline_vcf": record.get("snvs_germline", ""),
                    "variant_somatic_ann": record.get("variant_annotations_somatic_vcf", ""),
                    "variant_somatic_bcf": record.get("variant_annotations_somatic", ""),
                    "variant_germline_ann": record.get("variant_annotations_germline_vcf", ""),
                    "variant_germline_bcf": record.get("variant_annotations_germline", ""),
                    "snv_multiplicity": record.get("multiplicity", ""),
                    "oncokb_maf": record.get("oncokb_snv", ""),
                    "oncokb_fusions": record.get("oncokb_fusions", ""),
                    "oncokb_cna": record.get("oncokb_cna", ""),
                    "sbs_signatures": record.get("signatures_activities_sbs", ""),
                    "indel_signatures": record.get("signatures_activities_indel", ""),
                    "signatures_matrix": record.get("signatures_matrix_sbs", ""),
                    "signatures_matrix_indel": record.get("signatures_matrix_indel", ""),
                    "signatures_decomposed_sbs": record.get("signatures_decomposed_sbs", ""),
                    "signatures_decomposed_indel": record.get("signatures_decomposed_indel", ""),
                    "hrdetect": record.get("hrdetect", ""),
                    "onenesstwoness": record.get("onenesstwoness", ""),
                    "conpair_concordance": record.get("conpair_concordance", ""),
                    "conpair_contamination": record.get("conpair_contamination", "")
                })

            for row in sample_rows:
                # Filter the row to include only the selected fieldnames
                filtered_row = {key: row.get(key, "") for key in fieldnames}
                writer.writerow(filtered_row)
        if csv_path:
            output_stream.close()
