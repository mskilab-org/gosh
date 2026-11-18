import csv
import os

def check_if_tumor_only(samplesheet_path):
    """Validate all case statuses in one pass and return the uniform run mode."""
    if not samplesheet_path or not os.path.isfile(samplesheet_path):
        raise FileNotFoundError(f"Samplesheet not found at {samplesheet_path}")

    cases = {}
    with open(samplesheet_path, newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        fields = reader.fieldnames or []
        if "status" not in fields:
            raise ValueError("Samplesheet missing required 'status' column.")
        if "pair" not in fields and "patient" not in fields:
            raise ValueError("Samplesheet requires a 'pair' or 'patient' column.")
        for line, row in enumerate(reader, start=2):
            status = (row.get("status") or "").strip()
            if status not in ("0", "1"):
                raise ValueError(f"Samplesheet row {line}: status must be 0 or 1")
            case = (row.get("pair") or "").strip() or (row.get("patient") or "").strip()
            if not case:
                raise ValueError(f"Samplesheet row {line}: pair or patient must be nonempty")
            cases[case] = cases.get(case, 0) | (1 if status == "1" else 2)

    if not cases:
        raise ValueError("Samplesheet must contain at least one sample")
    for case, statuses in cases.items():
        if statuses == 2:
            raise ValueError(f"Samplesheet case {case} contains only normal samples")
    modes = set(cases.values())
    if len(modes) != 1:
        raise ValueError("Samplesheet mixes paired and singleton tumor cases")
    return modes == {1}
