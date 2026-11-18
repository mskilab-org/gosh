"""Operator-managed presets and immutable per-run parameter snapshots."""

import fcntl
import json
import os
import stat
import tempfile
from pathlib import Path

from .samplesheet import check_if_tumor_only
from ..settings import ALL_TOOLS, GENOME_MAP


BUILTIN_PRESETS = {
    "default": {"skip_tools": ""},
    "ffpe": {
        "skip_tools": "",
        "filter_ffpe_impact": True,
        "filter_ffpe_chimera": True,
        "sv_filter_ffpe_chimera": True,
    },
    "tumor-only": {
        "skip_tools": "",
        "purple_use_svs": False,
        "purple_use_smlvs": False,
    },
    "heme": {
        "skip_tools": "msisensorpro,hrdetect,onenesstwoness",
        "is_heme": True,
        "is_retier_whitelist_junctions": True,
        "purple_use_svs": False,
        "purple_use_smlvs": False,
        "purple_highly_diploid_percentage": 1.0,
        "purple_min_purity": 0.25,
        "purple_ploidy_penalty_factor": 0.6,
        "purple_revise_purity_ploidy": True,
    },
    "jabba": {"skip_tools": "sage,snpeff,snv_multiplicity,signatures,hrdetect"},
    "hrd": {"skip_tools": "non_integer_balance,lp_phased_balance,events,fusions"},
    "taps": {"skip_tools": "non_integer_balance,lp_phased_balance"},
}
PRESET_NAMES = tuple(BUILTIN_PRESETS)
_RUN_OWNED = {"input", "outdir", "genome", "email", "tumor_only"}
_SECRET_PARTS = ("secret", "password", "token", "api_key")
_BUILTIN_TYPES = {
    key: type(value)
    for values in BUILTIN_PRESETS.values()
    for key, value in values.items()
}


def _reject_constant(value):
    raise ValueError("Preset values must be valid JSON (no NaN or Infinity)")


def parse_json(value):
    """Decode JSON without Python's non-standard numeric constants."""
    return json.loads(value, parse_constant=_reject_constant)


def _check_name(name):
    if name not in BUILTIN_PRESETS:
        raise ValueError(f"Unknown preset: {name}")


def _check_key(key):
    if not isinstance(key, str) or not key.strip():
        raise ValueError("Preset parameter keys must be nonempty strings")
    normalized = key.lower()
    if any(part in normalized for part in _SECRET_PARTS):
        raise ValueError("Credentials cannot be stored in presets")
    if normalized in _RUN_OWNED:
        raise ValueError(f"Run-owned parameter cannot be set in a preset: {key}")


def _check_nested_secrets(value):
    if isinstance(value, dict):
        for key, child in value.items():
            if any(part in key.lower() for part in _SECRET_PARTS):
                raise ValueError("Credentials cannot be stored in presets")
            _check_nested_secrets(child)
    elif isinstance(value, list):
        for child in value:
            _check_nested_secrets(child)


def _validate_override(key, value):
    _check_key(key)
    expected = _BUILTIN_TYPES.get(key)
    if expected is float:
        valid = type(value) in (int, float)
    else:
        valid = expected is None or type(value) is expected
    if not valid:
        raise ValueError(f"Incorrect JSON type for preset parameter: {key}")
    _check_nested_secrets(value)
    # Also reject non-JSON types and non-finite numbers in programmatic callers.
    json.dumps(value, allow_nan=False)
    if key == "skip_tools":
        tools = {tool.strip() for tool in value.split(",")} if value else set()
        invalid = tools - set(ALL_TOOLS)
        if invalid:
            raise ValueError(f"Invalid skip_tools: {', '.join(sorted(invalid))}")


def _config_path(config_path):
    if config_path is not None:
        return Path(config_path)
    base = os.environ.get("XDG_CONFIG_HOME") or str(Path.home() / ".config")
    return Path(base) / "gosh" / "presets.json"


def _read_overrides(config_path):
    try:
        with config_path.open() as stream:
            data = json.load(stream, parse_constant=_reject_constant)
    except FileNotFoundError:
        return {}
    except (json.JSONDecodeError, ValueError) as exc:
        raise ValueError(f"Invalid preset config: {config_path}") from exc
    if not isinstance(data, dict):
        raise ValueError("Preset config must be a JSON object")
    for name, overrides in data.items():
        _check_name(name)
        if not isinstance(overrides, dict):
            raise ValueError(f"Overrides for {name} must be a JSON object")
        for key, value in overrides.items():
            _validate_override(key, value)
    return data


def resolve_preset(name, config_path=None) -> dict:
    """Return built-ins plus validated persistent overrides for one preset."""
    _check_name(name)
    data = _read_overrides(_config_path(config_path))
    return {**BUILTIN_PRESETS[name], **data.get(name, {})}


def _save_overrides(config_path, data):
    try:
        mode = stat.S_IMODE(config_path.stat().st_mode)
    except FileNotFoundError:
        mode = 0o600
    fd, temporary = tempfile.mkstemp(prefix=".presets-", dir=config_path.parent)
    try:
        with os.fdopen(fd, "w") as stream:
            os.fchmod(stream.fileno(), mode)
            json.dump(data, stream, indent=2, allow_nan=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, config_path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def set_override(name, key, value, config_path=None):
    _check_name(name)
    _validate_override(key, value)
    target = _config_path(config_path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.with_name("presets.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        data = _read_overrides(target)
        data.setdefault(name, {})[key] = value
        _save_overrides(target, data)


def unset_override(name, key, config_path=None):
    _check_name(name)
    _check_key(key)
    target = _config_path(config_path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.with_name("presets.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        data = _read_overrides(target)
        if key not in data.get(name, {}):
            return
        del data[name][key]
        if not data[name]:
            del data[name]
        _save_overrides(target, data)


def build_params(preset, samplesheet, outdir, genome, email="", config_path=None) -> dict:
    """Build a snapshot; sample modes and preset overrides are validated once."""
    effective = resolve_preset(preset, config_path)
    tumor_only = check_if_tumor_only(samplesheet)
    if preset == "tumor-only" and not tumor_only:
        raise ValueError("The tumor-only preset requires every sample to have status 1")
    reference = GENOME_MAP.get(genome.lower())
    if reference is None:
        raise ValueError(f"Unsupported genome reference: {genome}")
    params = {
        "input": os.fspath(samplesheet),
        "outdir": os.fspath(outdir),
        "genome": reference,
        "email": email,
        "tumor_only": tumor_only,
        "purple_use_svs": not tumor_only,
        "purple_use_smlvs": not tumor_only,
    }
    params.update(effective)
    return params
