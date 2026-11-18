import os
import sys
import json
from .presets import PRESET_NAMES, build_params


def create_params_file(
    preset="default",
    samplesheet="./samplesheet.csv",
    outdir="./results/",
    genome="hg19",
):
    default_input = samplesheet
    default_outdir = outdir
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

    # Check if default outdir exists
    outdir_exists = os.path.exists(default_outdir)
    outdir_status = "found!" if outdir_exists else "not found!"

    # Prompt for outdir
    if outdir_status == "found!":
        print(f"Default output directory found at: {default_outdir}")
    else:
        outdir_prompt = f"Enter output directory [default: {default_outdir} ({outdir_status})] (Press Enter to use default): "
        outdir = input(outdir_prompt).strip() or default_outdir

    if preset != "default":
        print(f"Preset automatically set to '{preset}'")
        preset_used = preset
    else:
        preset_prompt = (
            f"Available presets: {', '.join(PRESET_NAMES)}\n"
            "Enter preset option [default]: "
        )
        preset_used = input(preset_prompt).strip().lower() or "default"
        if preset_used not in PRESET_NAMES:
            print(f"Warning: Invalid preset '{preset_used}'. Using default.")
            preset_used = "default"

    # Prompt for genome
    genome_prompt = f"Enter genome [default: {genome}] (options: hg19, hg38, wmg-hg38) (Press Enter to use default): "
    genome = input(genome_prompt).strip() or genome

    # Prompt for email
    email_prompt = f"Enter email address [{default_email}] (Press Enter to skip): "
    email = input(email_prompt).strip() or ""

    params = build_params(preset_used, input_path, outdir, genome, email)
    print(
        "No normals found. Running in tumor-only mode."
        if params["tumor_only"]
        else "Running in paired tumor-normal mode."
    )

    custom_params_prompt = "Enter custom parameters? Y/N (default: N): "
    custom_params = input(custom_params_prompt).strip().lower() or "n"
    if custom_params == "y":

        def parse_value(raw):
            if raw.lower() == "true":
                return True
            if raw.lower() == "false":
                return False
            try:
                return int(raw)
            except ValueError:
                pass
            try:
                return float(raw)
            except ValueError:
                pass
            return raw

        while True:
            custom_params_key = input(
                "Parameter name (Press Enter when done): "
            ).strip()
            if not custom_params_key:
                confirm = (
                    input(
                        "No parameter name entered. Are you finished? Y/N (default: Y): "
                    )
                    .strip()
                    .lower()
                    or "y"
                )
                if confirm == "y":
                    break
                else:
                    continue
            raw_value = input(f"Value for '{custom_params_key}': ").strip()
            params[custom_params_key] = parse_value(raw_value)

    print("Parameters:")
    print(json.dumps(params, indent=4))

    # Write to params.json
    with open("params.json", "w") as f:
        json.dump(params, f, indent=4)
    print("Created 'params.json' with the provided parameters.")
