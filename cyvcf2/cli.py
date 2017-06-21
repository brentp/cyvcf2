import logging
from datetime import datetime

import coloredlogs
import click

from . import VCF

from . import __version__

log = logging.getLogger(__name__)
LOG_LEVELS = ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']


def print_header(vcf_obj, include, exclude):
    """Print the header of a vcf obj
    
    Parameters
    ----------
    vcf_obj: cyvcf2.VCF
    include: tuple
        set of strings with info fields that should be included
    exclude: tuple
        set of strings with info fields that should be excluded
    """
    for header_line in vcf_obj.raw_header.split('\n'):
        if len(header_line) == 0:
            continue
        print_line = True
        if 'INFO' not in header_line:
            pass
        elif include:
            print_line = any('ID='+inc+','for inc in include)
        elif exclude:
            print_line = not any('ID='+exc+','for exc in exclude)
        if print_line:
            click.echo(header_line)

def print_variant(variant, include, exclude):
    """Print a vcf variant
    
    Parameters
    ----------
    variant: cyvcf2.Variant
    include: tuple
        set of strings with info fields that should be included
    include: tuple
        set of strings with info fields that should be included
    """
    if include:
        for info in variant.INFO:
            key = info[0]
            if not key in include:
                del variant.INFO[key]
    if exclude:
        for exc in exclude:
            if variant.INFO.get(exc):
                del variant.INFO[exc]
                
    print_string = str(variant).rstrip()
    click.echo(print_string)
        
    

@click.command()
@click.argument('vcf',
    metavar='<vcf_file> or -',
)
@click.option('-c', '--chrom',
    help="Specify what chromosome to include.",
)
@click.option('-s', '--start',
    type=int,
    help="Specify the start of region.",
)
@click.option('-e', '--end',
    type=int,
    help="Specify the end of the region.",
)
@click.option('--include',
    help="Specify what info field to include.",
    multiple=True,
)
@click.option('--exclude',
    help="Specify what info field to exclude.",
    multiple=True,
)
@click.option('--loglevel',
    default='INFO',
    type=click.Choice(LOG_LEVELS),
    help="Set the level of log output.",
    show_default=True,
)
@click.option('--silent', 
    is_flag=True,
    help='Skip printing of vcf.'
)
@click.pass_context
def cyvcf2(context, vcf, include, exclude, chrom, start, end, loglevel, silent):
    """fast vcf parsing with cython + htslib"""
    coloredlogs.install(log_level=loglevel)
    start_parsing = datetime.now()
    log.info("Running cyvcf2 version %s", __version__)

    if include and exclude:
        log.warning("Can not use include and exclude at the same time")
        context.abort()

    region = ''
    if (chrom or start or end):
        if not (chrom and start and end):
            log.warning("Please specify chromosome, start and end for region")
            context.abort()
        else:
            region = "{0}:{1}-{2}".format(chrom, start, end)

    vcf_obj = VCF(vcf)
    
    for inclusion in include:
        if vcf_obj.contains(inclusion):
            log.info("Including %s in output", inclusion)
        else:
            log.warning("%s does not exist in header", inclusion)
            context.abort()

    for exclusion in exclude:
        if vcf_obj.contains(exclusion):
            log.info("Excluding %s in output", exclusion)
        else:
            log.warning("%s does not exist in header", exclusion)
            context.abort()
    
    if not silent:
        print_header(vcf_obj, include, exclude)
    
    nr_variants = None
    try:
        for nr_variants, variant in enumerate(vcf_obj(region)):
            if not silent:
                print_variant(variant, include, exclude)
    except Exception as err:
        log.warning(err)
        context.abort()
    
    if nr_variants:
        log.info("{0} variants parsed".format(nr_variants+1))
        log.info("Time to parse variants: {0}".format(datetime.now() - start_parsing))
        
    