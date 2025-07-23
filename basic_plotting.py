#import modules 

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from os import path


#import dataset

df = pd.read_csv('raw/State_data_3_states.csv')
outpath = "out"

# add %uptkae and Start-End columns to whole dataset
def add_columns(df):

    df['%Uptake'] = (df['Uptake']/df['MaxUptake'])*100
    df['startend'] = df['Start'].astype(str) + '-' + df['End'].astype(str)

    
    return df

df = add_columns(df)

# genrate lists of all unique peptides, states and startend 
def gen_uniques(df):

    unique_peptides = pd.unique(df['Sequence'])
    unique_states = pd.unique(df['State'])
    unique_startend = pd.unique(df['startend'])

    
    return unique_peptides, unique_states, unique_startend

unique_peptides, unique_states, unique_startend = gen_uniques(df)


# genetate uptake plots with all states in one plot
def allstates_uptake_plots(df, outpath):

    allstates_uptake_plots_folder = os.path.join(outpath, '3 states uptake plots')

    os.makedirs(allstates_uptake_plots_folder, exist_ok = True)
    
    for p in unique_peptides:
        
        # define figure parameters

        plt.figure(figsize=(15, 15))
        plt.rc('font', size=30)
        plt.rc('lines', linewidth=2)

        for s in unique_states:
            
            # pull relevant data from the df

            y = df.loc[((df['Sequence'] == p) & (df['State'] == s)), '%Uptake']
            x = df.loc[((df['Sequence'] == p) & (df['State'] == s)), 'Exposure']
            startendl = df.loc[((df['Sequence'] == p) & (df['State'] == s)), 'startend']
            startend = startendl.values[0]

            # plot the data

            plt.plot(x,y,alpha=0.9, marker = 'o', label = s)
        
        #change graph visuals
        
        plt.title(f'{startend} - {p}')
        plt.xlabel("Time (ms)", fontsize=28)
        plt.ylabel("Uptake (%)", fontsize=28)
        plt.tick_params(axis='both', which='major', labelsize=12)
        plt.legend(loc='upper left', fontsize=7)
        plt.tight_layout()
        plt.ylim(0,100)
        plt.semilogx()        
        plt.legend()
        plt.savefig(path.join(allstates_uptake_plots_folder, f"{startend} - {p}.png"))
        #plt.show()
        plt.close()

        print('allstates_uptake_plots Complete!')


#generate uptake plots for each peptide in each state individually
def uptake_plots(df, outpath):

    uptake_plots_folder = os.path.join(outpath, 'individual uptake plots')

    os.makedirs(uptake_plots_folder, exist_ok = True)
    
    for p in unique_peptides:

        for s in unique_states:
            
            #pull relevant data

            y = df.loc[((df['Sequence'] == p) & (df['State'] == s)), '%Uptake']
            x = df.loc[((df['Sequence'] == p) & (df['State'] == s)), 'Exposure']
            startendl = df.loc[((df['Sequence'] == p) & (df['State'] == s)), 'startend']
            startend = startendl.values[0]

            #plot and change graph visuals
            plt.plot(x,y,alpha=0.9, marker = 'o', label = s)
        
            plt.title(f'{startend} - {p}({s})')
            plt.xlabel("Time (ms)", fontsize=12)
            plt.ylabel("Uptake (%)", fontsize=12)
            plt.tick_params(axis='both', which='major', labelsize=12)
            plt.tight_layout()
            plt.ylim(0,100)
            plt.semilogx()        
            plt.savefig(path.join(uptake_plots_folder, f"{startend} - {p}({s}).png"))
            #plt.show()
            plt.close()

            print('uptake_plots Complete!')



