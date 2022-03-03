import os
import click
from .__init__ import cli
from comad.utils import (biom2data_tax,
                          non_neutral_outliers)
from comad.neufit import comad_pipeline, neufit

@cli.command(name='full_comad')
@click.option(
    '--biom',
    required=True,
    help='TODO')
@click.option(
    '--output_filename',
    required=True,
    help='TODO')
@click.option(
    '--output_folder_path',
    required=True,
    help='TODO')
def standalone_neufit(biom : str,
                      output_filename : str,
                      output_folder_path: str):
    '''Calls all functions needed to create neutral model 
    
    Written by: Caitlin Guccione, 08-25-2021
    
    Parameters
    ----------
    output_filename: str
        Name/nickname of dataset (ex. 'combined'). Will be incorperated into 
        filename of all comad outputs. 
    '''
    # run within wrapper
    comad_pipeline(biom, output_filename,
           output_folder_path)


@cli.command(name='neufit')
@click.option(
    '--fnData',
    required=True,
    help='TODO')
@click.option(
    '--fnTaxonomy',
    required=True,
    help='TODO')
@click.option(
    '--output-filename',
    required=True,
    help='TODO')
@click.option(
    '--output_folder_path',
    required=True,
    help='TODO')
def standalone_neufit(fnData : str,
                      fnTaxonomy : str,
                      output_filename : str,
                      output_folder_path: str):
    '''Calls all functions needed to create neutral model 
    
    Written by: Caitlin Guccione, 08-25-2021
    
    Parameters
    ----------
    output_filename: str
        Name/nickname of dataset (ex. 'combined'). Will be incorperated into 
        filename of all comad outputs. 
    '''
    # run within wrapper
    neufit(fnData, fnTaxonomy, output_filename,
           output_folder_path)