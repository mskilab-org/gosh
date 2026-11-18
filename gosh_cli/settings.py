GENOME_MAP = {
    "hg19": "GATK.GRCh37",  
    "gatk.grch37": "GATK.GRCh37", 
    "grch37": "GATK.GRCh37",
    "hg38": "GATK.GRCh38", 
    "gatk.grch38": "GATK.GRCh38",
    "grch38": "GATK.GRCh38", 
    "wmg-hg38": "WMG-hg38"
}
# Define the full list of tools available in the pipeline
ALL_TOOLS = [
    "aligner",
    "collect_wgs_metrics",
    "collect_multiple_metrics",
    "estimate_library_complexity",
    "msisensorpro",
    "gridss",
    "amber",
    "fragcounter",
    "dryclean",
    "cbs",
    "sage",
    "purple",
    "jabba",
    "non_integer_balance",
    "lp_phased_balance",
    "events",
    "fusions",
    "snpeff",
    "snv_multiplicity",
    "oncokb",
    "signatures",
    "hrdetect",
    "onenesstwoness",
]
ALL_TOOLS_STR = ",".join(ALL_TOOLS)