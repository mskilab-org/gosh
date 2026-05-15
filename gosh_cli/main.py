import click
from .cli.run import run_cli
from .cli.debug import debug_cli
from .cli.help import help_cli
from .cli.protect import protect_cli
from .cli.purge import purge_cli
from .cli.cache import cache_cli
from gosh_cli import __version__

@click.group()
@click.version_option(version=__version__, prog_name="gosh")
def cli():
    """gOSh - gOS sHell"""
    pass

# Register command groups
cli.add_command(run_cli, name='run')
cli.add_command(debug_cli, name='debug')
cli.add_command(help_cli, name='help')
cli.add_command(protect_cli, name='protect')
cli.add_command(purge_cli, name='purge')
cli.add_command(cache_cli, name='cache')
