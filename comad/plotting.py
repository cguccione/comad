import os
import pandas as pd
import numpy as np
from matplotlib import pyplot
from math import log10
from statsmodels.stats.proportion import proportion_confint
from scipy.stats import beta
from lmfit import fit_report

def beta_cdf(p, N, m):
    # Expected long term distribution under the neutral model (truncated cumulative beta-distribution)
    return beta.cdf(1.0, N*m*p, N*m*(1.0-p)) - beta.cdf(1.0/N, N*m*p, N*m*(1.0-p))

def neufit_plot(occurr_freqs, beta_fit, n_samples, n_reads, r_square, save_plot, HP_color):
    
    #Taken from : https://github.com/cguccione/neufit_gillespie_pipeline/blob/main/Fig4_selectionPlotting.ipynb#enroll-beta
    
    pyplot.cla() #Clears previous plot - to avoid double keys

    #Prepare results plot
    pyplot.xlabel('Mean relative abundance across samples', fontsize=18)
    pyplot.xscale('log')
    x_range = np.logspace(log10(min(occurr_freqs['mean_abundance'])/10), 0, 1000)
    pyplot.xlim(min(x_range), max(x_range))
    pyplot.xticks(fontsize=16)
    pyplot.ylabel('Occurrence frequency in samples', fontsize=18)
    pyplot.ylim(-0.05, 1.05)
    pyplot.yticks(fontsize=16)

    # Plot data points
    pyplot.plot(occurr_freqs['mean_abundance'], occurr_freqs['occurrence'], 'o', markersize=6, fillstyle='full', color='black')

    #>>Calculate lower and upper range
    lower, upper = proportion_confint(beta_cdf(x_range, n_reads, beta_fit.best_values['m'])*n_samples, n_samples, alpha=0.05, method='wilson')

    #>>Calculate/plot the main fit line
    #Orginal plotting colors
    pyplot.plot(x_range, beta_cdf(x_range, n_reads, beta_fit.best_values['m']), '-', lw=5, color='darkred')
    pyplot.plot(x_range, lower, '--', lw=2, color='darkred')
    pyplot.plot(x_range, upper, '--', lw=2, color='darkred')
    pyplot.fill_between(x_range, lower, upper, color='lightgrey')

    #Extract m value
    m = fit_report(beta_fit).split('(fixed)')[1]
    m = 'm = ' + str(round(float(m.strip('\n').split('+')[0].strip(' ').strip('m: ')),3))

    #Plot R^2 and m values
    pyplot.text(0.05, 0.9, m, fontsize=16, transform=pyplot.gca().transAxes)
    pyplot.text(0.05, 0.8, '$R^2 = ' + '{:1.2f}'.format(r_square) + '$', fontsize=16, transform=pyplot.gca().transAxes)
    
    if HP_color != False:
        hp_occurr_freqs = occurr_freqs.loc['k__Bacteria;p__Proteobacteria;c__Epsilonproteobacteria;o__Campylobacterales;f__Helicobacteraceae;g__Helicobacter;s__Helicobacter pylori' , :]
        pyplot.plot(hp_occurr_freqs['mean_abundance'], hp_occurr_freqs['occurrence'], 'o', markersize=6, fillstyle='full', color='magenta')
    
    if save_plot != False:
        #Save plot
        pyplot.tight_layout()
        pyplot.savefig(save_plot + '.pdf')
        
        #Save df
        occurr_freqs.to_csv(save_plot + '.tsv', sep='\t')

    pyplot.show()