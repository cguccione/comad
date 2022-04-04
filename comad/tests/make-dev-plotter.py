#!/usr/bin/env python

import click
import os
import pkg_resources
from comad.utils import biom2data_tax
from jinja2 import Environment, FileSystemLoader

ARG_TYPE = click.Path(exists=True, dir_okay=False, file_okay=True)

SUPPORT_FILES = pkg_resources.resource_filename('comad', 'support_files')
TEMPLATES = os.path.join(SUPPORT_FILES, 'templates')

@click.command()
@click.option(
    '--biom',
    required=True,
    type=ARG_TYPE,
    help='TODO')
def main(biom):
     #Convert data from biom to csv files for Neufit
    fnData, fnTaxonomy = biom2data_tax(biom, "test", './data')
    env = Environment(loader=FileSystemLoader(TEMPLATES))
    template = env.get_template('test-template.html')
    with open('development-page.html', 'w') as f:
    	f.write(template.render({
            'base_url': SUPPORT_FILES})

if __name__ == '__main__':
	main()