import os
import sys
import json
from .samplesheet import check_if_tumor_only

genome_map = {
    'hg19': 'GATK.GRCh37',
    'hg38': 'GATK.GRCh38'
}
def create_params_file(
    preset="default",
    samplesheet="./samplesheet.csv",
    outdir="./results/",
    genome="hg19",
):
    default_input = samplesheet
    default_outdir = outdir
    default_tools = "all"
    default_genome = genome_map[genome]
    default_email = "example_email@gmail.com"

    # Check if default input exists
    input_exists = os.path.exists(default_input)

    # Prompt for input
    if input_exists:
        print(f"Default input samplesheet found at: {default_input}")
        input_path = default_input
    else:
        input_prompt = f"Enter samplesheet CSV file path (Press Enter to use default [./samplesheet.csv]): "
        input_path = input(input_prompt).strip()
        if not input_path:
            if input_exists:
                input_path = default_input
            else:
                print("Error: You must provide a path to the samplesheet.")
                sys.exit(1)
        else:
            if not os.path.exists(input_path):
                print(f"Error: The provided input file '{input_path}' does not exist.")
                sys.exit(1)

    # Determine if tumor-only mode should be enabled
    try:
        is_tumor_only = check_if_tumor_only(input_path)
        mode = "No normals found. Running in tumor-only mode." if is_tumor_only else "Running in paired tumor-normal mode."
        print(mode)
    except (FileNotFoundError, ValueError) as e:
        print(f"Error checking samplesheet: {e}")
        sys.exit(1)

    # Check if default outdir exists
    outdir_exists = os.path.exists(default_outdir)
    outdir_status = "found!" if outdir_exists else "not found!"

    # Prompt for outdir
    if outdir_status == "found!":
        print(f"Default output directory found at: {default_outdir}")
    else:
        outdir_prompt = f"Enter output directory [default: {default_outdir} ({outdir_status})] (Press Enter to use default): "
        outdir = input(outdir_prompt).strip() or default_outdir

    # New preset based configuration.
    presets = {
        "default": "",
        "jabba": "sage,snpeff,snv_multiplicity,signatures,hrdetect",
        "hrd": "non_integer_balance,lp_phased_balance,events,fusions",
        "heme": "msisensorpro,hrdetect,onenesstwoness"
    }

    if preset != "default":
        print(f"Preset automatically set to '{preset}'")
        preset_used = preset
    else:
        preset_prompt = (
            "Available presets:\n"
            " - default (recommended): runs all tools; aligner, bamqc, gridss, amber, fragcounter, dryclean, cbs, sage, purple, jabba, non_integer_balance, lp_phased_balance, events, fusions, snpeff, snv_multiplicity, signatures, hrdetect \n"
            " - jabba: runs all tools necessary for JaBbA outputs (skips tools: sage, snpeff, snv_multiplicity, signatures, hrdetect)\n"
            " - hrd: runs HR deficiency pipeline (skips tools: non_integer_balance, lp_phased_balance, events, fusions)\n"
            " - heme: runs heme pipeline (skips tools: msisensorpro, hrdetect, onenesstwoness)\n"
            "Enter preset option (options: default, jabba, hrd, heme): "
        )
        preset_used = input(preset_prompt).strip().lower() or "default"
        if preset_used not in presets:
            print(f"Warning: Invalid preset '{preset_used}'. Using default.")
            preset_used = "default"

    # Prompt for genome
    if not genome in genome_map.values():
        genome_prompt = (
            f"Enter genome [default: hg19] (options: hg19, hg38) (Press Enter to use default): "
        )
        genome_input = input(genome_prompt).strip() or "hg19"
        genome = genome_map.get(genome_input.lower())

    if not genome == "GATK.GRCh37" and not genome == "GATK.GRCh38":
        print(f"Warning: Invalid genome '{genome_input}'. Using default '{default_genome}'.")
        genome = genome_map[default_genome]

    # Prompt for email
    email_prompt = f"Enter email address [{default_email}] (Press Enter to skip): "
    email = input(email_prompt).strip() or ""

    # Create params dictionary
    params = {
        "input": input_path,
        "outdir": outdir,
        "genome": genome,
        "email": email,
        "tumor_only": is_tumor_only
    }
    if preset_used != "default":
        params["skip_tools"] = presets[preset_used]
    if preset_used == "heme":
        print("Adding heme specific parameters...")
        params.update({
            "is_retier_whitelist_junctions": "true",
            "purple_use_svs": "false",
            "purple_use_smlvs": "false",
            "purple_highly_diploid_percentage": 1.0,
            "purple_min_purity": 0.25,
            "purple_ploidy_penalty_factor": 0.6
        })

    print("Parameters:")
    print(json.dumps(params, indent=4))

    # Write to params.json
    with open("params.json", "w") as f:
        json.dump(params, f, indent=4)
    print("Created 'params.json' with the provided parameters.")
