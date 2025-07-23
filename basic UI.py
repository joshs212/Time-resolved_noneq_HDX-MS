#import modules 

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from os import path

# import dataframe

dataframe = pd.read_csv('raw/State_data_3_states.csv')

# choose desired out file destination
outpath = 'out'


## here you select which functions you want to run:

#this is only required when generating matplotlib uptake plots, otherwise set to false 
generate_matplotlib_plots = False

# this is only needed if generating k values, turn off otherwise
generate_k_plots = False

#choose desired k value and any peptides to highlight on the grey cluster scatter plot
input_k_value = 7
choose_selected_peptides = ['222-240', '246-252', '244-256', '203-218', '203-219']

#this is only needed if generating cluster scatter plots
generate_cluster_scatter_plots = False

#k value is required for this section to work, generates plots using novel technique
generate_plots_from_new_analysis = True




## This section does basic formatting for the df and generates basic matplotlib uptake plots 
from basic_plotting import (add_columns, gen_uniques, uptake_plots, allstates_uptake_plots)

# these are required for everything
df = add_columns(dataframe)

unique_peptides, unique_states, unique_startend = gen_uniques(df)


if generate_matplotlib_plots == True:
    allstates_uptake_plot = allstates_uptake_plots(df, outpath)
    individual_uptake_plots = uptake_plots(df,outpath)
else:
    print('allstates_uptake_plot is turned off!')
    print('individual_uptake_plots is turned off!')

## This next section is the Sum Uptake clustering technique 

from clusters import (process_cluster_data, choose_k, gen_clusterdata, gen_scatter_plot, grey_scatter_plot)


if generate_k_plots == True:
    delta_dfs = process_cluster_data(df)
    k_selection = choose_k(delta_dfs, outpath, max_clusters = 12)
else:
    print('k_selection is turned off!')



if generate_cluster_scatter_plots == True:
    delta_dfs = process_cluster_data(df)
    #data processing for clustering method
    clusterdata, unique_labels, centroids = gen_clusterdata(delta_dfs, outpath, n_clusters = input_k_value)
    
    #function for generating cluster scatter plot
    sc_plot = gen_scatter_plot(clusterdata, unique_labels, centroids, outpath)

    #you can change the selected peptides by changing their aa sequence
    grey_sc = grey_scatter_plot(clusterdata, outpath, unique_labels, centroids, selected_peptides = choose_selected_peptides)
else:
    print('gen_clusterdata is turned off!')
    print('gen_scatter_plot is turned off!')
    print('grey_scatter_plot is turned off!')


## this is the new code section relating to my new delta at each time point analysis

from integration_of_recursive import (initial_processing, recursive_func, deletion_func, fix_new_issue, grad_generator, gen_gradients, gen_deltas, data_processingcleaning, plot_by_cluster, plot_individual_colored)

if generate_plots_from_new_analysis == True:
    # function that does all data processing for the novel analysis 
    plotting_ready = data_processingcleaning(df, outpath)

    #plotting functions
    plot_indiv_colored = plot_individual_colored(plotting_ready, outpath)
    organise_by_cluster = plot_by_cluster(plotting_ready, outpath)
else:
    print('data_processingcleaning is turned off!')
    print('plot_individual_colored is turned off!')
    print('plot_by_cluster is turned off!')


