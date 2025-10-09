import click


@click.group(name='asard',
             no_args_is_help=True,
             invoke_without_command=True)
@click.option('--version', is_flag=True,
              help='Print asard version information and exit. Overrides all other arguments.')
def cli(version):
    if version:
        import ERS_NRB
        print(ERS_NRB.__version__)


@cli.group(name='rb',
           no_args_is_help=True,
           invoke_without_command=True)
def rb():
    pass


@rb.command(name='process',
            no_args_is_help=True,
            context_settings=dict(
                ignore_unknown_options=True,
                allow_extra_args=True, )
            )
@click.option('--config-file', '-c', required=False, type=click.Path(),
              help="Full path to an INI-style configuration text file. "
                   "If not defined, the package's default file will be used.")
@click.option('--debug', is_flag=True,
              help='Print debugging information for pyroSAR modules.')
@click.pass_context
def process(ctx, config_file, debug):
    """
    Central asard radar backscatter (rb) processing command.

    A config file can be defined from which all configuration is read.
    If not defined, configuration will be read from the package's default file.
    Additional options can be passed to override individual processing parameters
    in the configuration file. For example, to read all values from the configuration
    file except the acquisition mode and the annotation layers:

    asard rb process -c config.ini --acq_mode IMP --annotation dm,id
    """
    from ERS_NRB import process
    extra = {ctx.args[i][2:]: ctx.args[i + 1] for i in range(0, len(ctx.args), 2)}
    process(config_file=config_file, debug=debug, **extra)
