"""Persistent preset overrides for future parameter snapshots."""

import json

import click

from ..core.presets import (
    PRESET_NAMES,
    parse_json,
    resolve_preset,
    set_override,
    unset_override,
)


@click.group(name="preset")
def preset_cli():
    """Inspect and edit pipeline presets (affects newly prepared runs only)."""


@preset_cli.command(name="list")
def list_presets():
    """List the supported preset names."""
    for name in PRESET_NAMES:
        click.echo(name)


@preset_cli.command(name="show")
@click.argument("name", type=click.Choice(PRESET_NAMES))
def show_preset(name):
    """Show effective preset values as JSON."""
    try:
        effective = resolve_preset(name)
    except (OSError, ValueError) as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(json.dumps(effective, indent=2, allow_nan=False))


@preset_cli.command(name="set")
@click.argument("name", type=click.Choice(PRESET_NAMES))
@click.argument("key")
@click.argument("json_value")
def set_preset(name, key, json_value):
    """Persist a parameter override. JSON_VALUE must be valid JSON."""
    try:
        value = parse_json(json_value)
        set_override(name, key, value)
    except (OSError, ValueError, TypeError) as exc:
        raise click.ClickException(str(exc)) from exc


@preset_cli.command(name="unset")
@click.argument("name", type=click.Choice(PRESET_NAMES))
@click.argument("key")
def unset_preset(name, key):
    """Remove an override, restoring the built-in value when present."""
    try:
        unset_override(name, key)
    except (OSError, ValueError) as exc:
        raise click.ClickException(str(exc)) from exc
