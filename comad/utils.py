import pandas as pd
from biom import load_table 
import os
import numpy as np

def biom2data_tax(biom_filename, output_filename, output_folder_path):
    '''Imports biom file -> pandas dataframe -> data.csv, taxonomy.csv 
    
    Written by: Caitlin Guccione, 08-25-2021
    
    Parameters
    ----------
    biom_filename: str
        The filename of biom file
    output_filename: str
        The filename to be used as the header in _data.csv and
        _taxonomy.csv files
    output_folder_path: str
        Path of all output files

     Returns
     -------
    fnD: str
        The file location of [name]_data.csv
    csv file
        Creates[name]_data.csv, a data.csv file needed to run Neufit 
    fnT: str
        The file location of [name]_taxonomy.csv
    csv file
        Creates [name]_taxonomy.csv, a taxonomy.csv file needed
        to run Neufit
      
     TODO
     ----
     - Insure there isn't a better way to inpur filename
     - Make sure the path direction is clear for all inputs
     - Add more clear description of csv files 
     
      '''

    #Create folder for _data.csv and _tax.csv files for Neufit input
    data_tax_path = output_folder_path + '/data_tax_csv/'
    os.makedirs(data_tax_path, exist_ok=True) 

    
    #Make filename and import data
    featureTable = load_table(biom_filename) 
    #https://biom-format.org/documentation/generated/biom.load_table.html
    
    #Create _data.csv
    pandas_featureTable = pd.DataFrame(featureTable.matrix_data.toarray(),
                                       featureTable.ids('observation'), 
                                       featureTable.ids())
    fnD = data_tax_path +  output_filename + '_data.csv'
    pandas_featureTable.to_csv(fnD, sep='\t')
    
    #Create _taxonomy.csv
    pandas_TaxTable = pd.DataFrame(featureTable.metadata_to_dataframe('observation'))
    pandas_TaxTable.set_axis(['Kingdom', 'Phylum', 'Class', 'Order',
                              'Family', 'Genus', 'Species'], axis=1)
    fnT = data_tax_path + output_filename + '_taxonomy.csv'
    pandas_TaxTable.to_csv(fnT, sep='\t')
    
    return(fnD, fnT)


def non_neutral_outliers(file_header, occurr_freqs, threshold = 0.5):
    ''' Creates the most Non-neutral csv file 
    
    Written by: Caitlin Guccione, 08-25-2021
    
    Parameters
    ----------
    file_header: str, path
        Filepath for all comad outputs. Includes path, data nickname and 
        time stamp.
    occurr_freqs: pandas df
        Df header: otu_id, mean_abundance, occurrence, Kingdom, Phylum, 
        Class, Order, Family, Genus, Species, predicted_occurrence, 
        lower_conf_int, upper_conf_int
    threshold: int, optional
        Autoset to 0.5, but determines which bacteria are considered 
        non-neutral strictly for this csv file
        
     Returns
     -------
     csv file
         [name]_NonNeutralOutliers.csv, holds all the species 
         furthest off the neutral curve
    
    TODO
    ----
    - See if there is a way to eliminate the Domain vs Kingdom 
    differences more cleanly 
    
    '''
    
    #Create standout microbes df with approriate col names
    occurr_freqs_col_names = list(occurr_freqs.columns)
    standoutMicrobes_col_names = ['Difference off Neutral Model']
    start = False
    for i in occurr_freqs_col_names:
        if i == 'occurrence':
            start = True 
        elif start == True:
            if i != 'predicted_occurrence':
                standoutMicrobes_col_names.append(i)
            else:
                break
    standoutMicrobes = pd.DataFrame(columns = standoutMicrobes_col_names)
    
    #Loop and find most non-neutral microbes based upon threshold
    row_count = 0

    for index, row in occurr_freqs.iterrows():
        diff = abs(row['occurrence'] - row['predicted_occurrence'])
        if diff > threshold:
            new_row = [diff]
            for i in standoutMicrobes_col_names[1:]:
                new_row.append(row[i])
            standoutMicrobes.loc[row_count] = new_row
            row_count +=1
    standoutMicrobes = standoutMicrobes.sort_values(by =['Difference off Neutral Model'],
                                                    ascending=False)
    
    #Display and export non-neutral microbes as csv
    print("\nTop NonNeutral Microbes")
    display(standoutMicrobes)
    print('=========================================================\n')
    
    fn = file_header + '_NonNeutral_Outliers.csv'
    
    standoutMicrobes.to_csv(fn, sep = ',', index=False)