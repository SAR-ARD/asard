import click
import ERS_NRB.processor as process


@click.command()
@click.option('--config-file', '-c', required=True, type=click.Path(),
              help='Full path to an INI-style configuration text file.')
def cli(config_file):
    process.main(config_file=config_file)
