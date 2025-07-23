import pandas as pd
import numpy as np
import scipy
from scipy import signal
import matplotlib.pyplot as plt
from basic_plotting import (gen_uniques, add_columns)
from clusters import(process_cluster_data, gen_clusterdata)
from mpl_toolkits import mplot3d
import os 
from os import path
import matplotlib.colors as mcolors
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm

#import df
df = pd.read_csv('raw/State_data_3_states.csv')

df = add_columns(df)
unique_peptides, unique_states, unique_startend = gen_uniques(df)

#set outpath
outpath = 'out'

def initial_processing(df, outpath):

    #obtain clusterdata to enable organisation by cluster
    delta_dfs = process_cluster_data(df)
    clusterdata, unique_labels, centroids = gen_clusterdata(delta_dfs, outpath)

    big_list = []

    for p in unique_peptides:

        df_list = []

        for s in unique_states:
            
            #pull relevant info from df
            y = df.loc[((df['Sequence'] == p) & (df['State'] == s)), 'Uptake']
            x = df.loc[((df['Sequence'] == p) & (df['State'] == s)), 'Exposure']
            startend = df.loc[((df['Sequence'] == p) & (df['State'] == s)), 'startend']

            #convert to series so it can be used
            se = pd.Series(startend)

            #use savgol_filter function to smooth data 
            y2 = signal.savgol_filter(y, 3,1)

            # initiate df to store data, creates a df for each peptide state
            initial_format_df = pd.DataFrame(data = [], columns = ['startend','labels', 'peptide','state','x', 'y'])

            initial_format_df['x'] = x.values
            initial_format_df['y'] = y2
            initial_format_df['peptide'] = p
            initial_format_df['state'] = s
            initial_format_df['startend'] = se.values

            df_list.append(initial_format_df)

        df_df = pd.concat(df_list)

        big_list.append(df_df)

    #now reformatted df using smoothed data and desired columns
    first_format_df = pd.concat(big_list)

    peptide_clusters = clusterdata[['peptide', 'labels']]

    unique_labels = np.sort(unique_labels)

    #map labls from clusterdata to correct peptides
    peptide_to_label = dict(zip(peptide_clusters['peptide'], peptide_clusters['labels']))
    first_format_df['labels'] = first_format_df['peptide'].map(peptide_to_label)

    first_format_df.to_csv(path.join(outpath, 'firstformat_df.csv'))

    return(first_format_df, unique_labels)


#function to delete any uptake values that would result in a negative gradient, things removed using hash can be restored to check if func is working as intended
def recursive_func(df):

    #print("\nStarting recursive function")

    i = 0
    # loop to iterate through each uptake value
    while i < len(df) - 1:
        
        #checks if value will result in negative gradient
        if df['y'].iloc[i + 1] - df['y'].iloc[i] < 0:

            #print(f"❌ Removing row {i+1} (y value: {df['y'].iloc[i + 1]})")

            #deletes value and resets the loop to 0
            df = df.drop(df.index[[i + 1]], axis = 0).reset_index(drop=True)

            i = 0

        else:
            
            #loop continues if value is ok
            i += 1          
    
    #print("\n✅ Final cleaned DataFrame:")
    #print(df)

    return(df)

# function built to specifically check if recursive function is working as intended, only needs to be ran when checking that
def check_rec_func(df):
    
    all_clear = True

    for i in range(len(df) -1):

        if df['y'].iloc[i + 1] - df['y'].iloc[i] < 0:

            print(f"negative number spotted {i + 1}")

            all_clear = False

        if all_clear:

            print("all clear!")         

    return

#func that integrates recursive func
def deletion_func(df, outpath):

    del_list = []

    for p in unique_peptides:
        
        perpep = []

        for s in unique_states:
            
            #pull relevant data
            y = df.loc[((df['peptide'] == p) & (df['state'] == s)), 'y']
            x = df.loc[((df['peptide'] == p) & (df['state'] == s)), 'x']
            startend = df.loc[((df['peptide'] == p) & (df['state'] == s)), 'startend']
            labels = df.loc[((df['peptide'] == p) & (df['state'] == s)), 'labels']

            #start new df with only these columns
            filtered_df = pd.DataFrame(data = [], columns = ['startend','peptide','state','x', 'y', 'labels'])

            filtered_df['x'] = x
            filtered_df['y'] = y
            filtered_df['startend'] = startend
            filtered_df['labels'] = labels
            filtered_df['peptide'] = p
            filtered_df['state'] = s

            #implement recursive func 
            filtered_df = recursive_func(filtered_df)

            #appends df that has had its problematic values removed
            perpep.append(filtered_df)

        perpep_df = pd.concat(perpep)

        del_list.append(perpep_df)
    #returns df that has been fully sorted by recursive func
    deleted_df = pd.concat(del_list)

    deleted_df.to_csv(path.join(outpath, 'post_recursive.csv'))

    return(deleted_df)

# due to the deletion many peptides have different lengths of uptake values between the states, cannot compare the states like this
#this function is used to extend every df to have 31 timepoints and if it has been deleted an NA value is inserted there
def expand_with_nan(unique_timepoints, time, data):

    expanded_time = []
    expanded_data = []

    for i in range(len(unique_timepoints)):

        if unique_timepoints[i] in time:
            
            time = list(time)
            index = time.index(unique_timepoints[i])
            expanded_time.append(unique_timepoints[i])
            expanded_data.append(data[index])

        else:

            expanded_time.append(np.nan)
            expanded_data.append(np.nan)

    return(expanded_time, expanded_data)

    
#function that integrates expanded_with_nan func 
def fix_mismatch_issue(df, outpath):

    unique_timepoints = df['x'].unique()
    unique_timepoints.sort()

    flattened_df_list = []

    for p in unique_peptides:
        
        #create new df with the different states as columns not rows
        newshape_df = pd.DataFrame(data= [], columns = ['startend', 'peptide', 'labels', 'apo_time', 'eq_time', 'noneq_time', 'apo_data', 'eq_data', 'noneq_data'])

        startend = df.loc[((df['peptide'] == p) & (df['state'] == unique_states[0])), 'startend']
        labels = df.loc[((df['peptide'] == p) & (df['state'] == unique_states[0])), 'labels']

        apo_time = df.loc[((df['peptide'] == p) & (df['state'] == unique_states[0])), 'x']
        eq_time = df.loc[((df['peptide'] == p) & (df['state'] == unique_states[1])), 'x']
        noneq_time = df.loc[((df['peptide'] == p) & (df['state'] == unique_states[2])), 'x']

        apo_y = df.loc[((df['peptide'] == p) & (df['state'] == unique_states[0])), 'y']
        eq_y = df.loc[((df['peptide'] == p) & (df['state'] == unique_states[1])), 'y']
        noneq_y = df.loc[((df['peptide'] == p) & (df['state'] == unique_states[2])), 'y']

        #expand time values for all states
        expanded_apo_time, expanded_apo_data = expand_with_nan(unique_timepoints, apo_time, apo_y)
        expanded_eq_time, expanded_eq_data = expand_with_nan(unique_timepoints, eq_time, eq_y)
        expanded_noneq_time, expanded_noneq_data = expand_with_nan(unique_timepoints, noneq_time, noneq_y)

       #these factors kept the same within a peptide
        first_label = labels.iloc[0] if not labels.empty else np.nan
        first_startend = startend.iloc[0] if not startend.empty else np.nan

        expanded_startend = [first_startend] * len(expanded_apo_time)
        expanded_labels = [first_label] * len(expanded_apo_time)

        #new df with the expanded values
        newshape_df['apo_time'] = expanded_apo_time
        newshape_df['eq_time'] = expanded_eq_time
        newshape_df['noneq_time'] = expanded_noneq_time
        newshape_df['apo_data'] = expanded_apo_data
        newshape_df['eq_data'] = expanded_eq_data
        newshape_df['noneq_data'] = expanded_noneq_data
        newshape_df['peptide'] = p
        newshape_df['startend'] = expanded_startend
        newshape_df['labels'] = expanded_labels

        flattened_df_list.append(newshape_df)

    # return df with all peptides having been expanded
    flattened_df = pd.concat(flattened_df_list)

    #removes any row where a value is na so states are comparable
    flattened_df = flattened_df.dropna().reset_index(drop = True)

    #system to remove peptides that are too short
    peptides_to_remove = []

    for p in unique_peptides:

        apo_data = flattened_df.loc[(flattened_df['peptide'] == p), 'apo_data']

        if len(apo_data) < 15:

            peptides_to_remove.append(p)

    
    print(f'removed peptides:{peptides_to_remove}')
    print(f'number of removed peptides: {len(peptides_to_remove)}')

    #remove peptides
    flattened_df = flattened_df[~flattened_df['peptide'].isin(peptides_to_remove)]

    flattened_df.to_csv(path.join(outpath, 'flattened_df.csv'))
    
    return(flattened_df)


#function to generate gradients
def grad_generator(x, y):

    delta_y = np.diff(y)
    delta_x = np.diff(x)

    grads = []

    for dy,dx in zip(delta_y, delta_x):

        grad = dy/dx

        grads.append(grad)

    return(grads, delta_y, delta_x)

#function to implement func to generate gradients
def gen_gradients(df, outpath):

    big_df = []

    for p in unique_peptides:

        se = df.loc[(df['peptide'] == p), 'startend'].values
        labels = df.loc[(df['peptide'] == p), 'labels'].values

        apo_data = df.loc[(df['peptide'] == p), 'apo_data']
        eq_data = df.loc[(df['peptide'] == p), 'eq_data']
        noneq_data = df.loc[(df['peptide'] == p), 'noneq_data']
        
        apo_time = df.loc[(df['peptide'] == p), 'apo_time'].values
        eq_time = df.loc[(df['peptide'] == p), 'eq_time']
        noneq_time = df.loc[(df['peptide'] == p), 'noneq_time']

        apo_grad, apo_dy, apo_dx = grad_generator(apo_time, apo_data)
        eq_grad, eq_dy, eq_dx = grad_generator(eq_time, eq_data)
        noneq_grad, noneq_dy, noneq_dx = grad_generator(noneq_time, noneq_data)

        dydx_df = pd.DataFrame(data = [], columns = ['startend','labels', 'peptide','apo_x', 'dx', 'apo_dy', 'eq_dy', 'noneq_dy', 'apo_grad', 'eq_grad', 'noneq_grad'])

        dydx_df['apo_dy'] = apo_dy
        dydx_df['eq_dy'] = eq_dy
        dydx_df['noneq_dy'] = noneq_dy
        dydx_df['dx'] = apo_dx
        dydx_df['apo_x'] = apo_time[1:] 
        dydx_df['peptide'] = p  
        dydx_df['apo_grad'] = apo_grad
        dydx_df['eq_grad'] = eq_grad
        dydx_df['noneq_grad'] = noneq_grad
        dydx_df['startend'] = se[1:]
        dydx_df['labels'] = labels[1:]

        big_df.append(dydx_df)

    dydx_df = pd.concat(big_df)

    dydx_df.to_csv(path.join(outpath, 'grads_df.csv'))

    return(dydx_df)

# func to generate the difference between the gradients at each time point
def gen_deltas(df, outpath):

    big_df = []

    for p in unique_peptides:

        startend = df.loc[(df['peptide'] == p), 'startend']
        labels = df.loc[(df['peptide'] == p), 'labels']
        apo_x = df.loc[(df['peptide'] == p), 'apo_x']
        dx = df.loc[(df['peptide'] == p), 'dx']

        apo_grad = df.loc[(df['peptide'] == p), 'apo_grad']
        eq_grad = df.loc[(df['peptide'] == p), 'eq_grad']
        noneq_grad = df.loc[(df['peptide'] == p), 'noneq_grad']

        apo_noneq = noneq_grad - apo_grad
        eq_noneq = noneq_grad - eq_grad

        deltas_df = pd.DataFrame(data = [], columns = ['startend', 'labels', 'peptide', 'apo_x', 'dx', 'apo_noneq', 'eq_noneq'])

        deltas_df['apo_noneq'] = apo_noneq
        deltas_df['eq_noneq'] = eq_noneq
        deltas_df['peptide'] = p
        deltas_df['startend'] = startend
        deltas_df['labels'] = labels
        deltas_df['apo_x'] = apo_x
        deltas_df['dx'] = dx

        big_df.append(deltas_df)

    deltas_df = pd.concat(big_df)

    #normalising the data
    max_apo_noneq = deltas_df['apo_noneq'].abs().max()
    max_eq_noneq = deltas_df['eq_noneq'].abs().max()

    deltas_df['norm_apo_noneq'] = deltas_df['apo_noneq']/max_apo_noneq
    deltas_df['norm_eq_noneq'] = deltas_df['eq_noneq']/max_eq_noneq

    deltas_df.to_csv(path.join(outpath, 'deltas_df.csv'))

    return(deltas_df)



# a function to package all data processing func together
def data_processingcleaning(df, outpath):

    first_format_df, unique_labels = initial_processing(df, outpath)

    cleaned_df = deletion_func(first_format_df, outpath)

    flattened_df = fix_mismatch_issue(cleaned_df, outpath)

    grads_df = gen_gradients(flattened_df, outpath)

    plotting_ready = gen_deltas(grads_df, outpath)

    print('data_processingcleaning Complete!')
    return(plotting_ready)

#funcion to plot each peptide indivudually
def plot_by_cluster(df, outpath):

    #sorts them into folders based on the cluster the peptide belongs to
    full_path = os.path.join(outpath, 'plots by cluster')
    os.makedirs(full_path, exist_ok = True)

    unique_labels = df['labels'].unique()

    for l in unique_labels:
        
        selected_peptides = df.loc[(df['labels'] == l), 'peptide']

        filtered_df = df[df['peptide'].isin(selected_peptides)]
        
        filtered_unique_peptides = filtered_df['peptide'].unique()

        for p in filtered_unique_peptides:

            #plot for each peptide
            x = filtered_df.loc[(filtered_df['peptide'] == p), 'norm_apo_noneq']
            y = filtered_df.loc[(filtered_df['peptide'] == p), 'norm_eq_noneq']

            plt.plot(x,y, label = p)
            plt.scatter(x,y)
    
        #plotting them all on the same graph
        plt.title(f'{l}')
        plt.ylim(-1.1,1.1)
        plt.xlim(-1.1,1.1)
        plt.axhline(xmin = -1.1, xmax = 1.1, y = 0, linestyle = '--')
        plt.axvline(ymin = -1.1, ymax = 1.1, x = 0, linestyle = '--')
        plt.xlabel('Difference from apo', fontsize=11)
        plt.ylabel('Difference from eq', fontsize=11)
        plt.legend()

        #plt.show()
        plt.savefig(path.join(full_path, f'{l}.png'))

    print('plot_by_cluster Complete!')



    return()


# same as previous func but adds in colour change based on time 
def plot_individual_colored(df, outpath):

    full_path = os.path.join(outpath, 'colored individual plots by cluster')
    unique_labels = df['labels'].unique()
    os.makedirs(full_path, exist_ok = True)

    for l in unique_labels:

        label_folders = os.path.join(full_path, str(l))

        os.makedirs(label_folders, exist_ok = True)

        selected_peptides = df.loc[(df['labels'] == l), 'peptide']

        filtered_df = df[df['peptide'].isin(selected_peptides)]
        
        filtered_unique_peptides = filtered_df['peptide'].unique()

        for p in filtered_unique_peptides:

            x = df.loc[((df['peptide'] == p) & (df['labels'] == l)), 'norm_apo_noneq']
            y = df.loc[((df['peptide'] == p) & (df['labels'] == l)), 'norm_eq_noneq']
            startendl = df.loc[((df['peptide'] == p) & (df['labels'] == l)), 'startend']
            time = df.loc[((df['peptide'] == p) & (df['labels'] == l)), 'apo_x']
            startend = startendl.values[0]


            # Create segments for the line plot
            points = np.column_stack([x, y])
            segments = np.array([points[:-1], points[1:]]).transpose(1, 0, 2)  # Create line segments

            # Color normalization for time points
            vmin_adjusted = np.percentile(time, 1)
            vmax_adjusted = 10**3
            norm = LogNorm(vmin=vmin_adjusted, vmax=vmax_adjusted)

            cmap = plt.get_cmap('plasma_r')
            lc = LineCollection(segments, cmap=cmap, norm=norm, linewidth=2)
            lc.set_array(time[:-1])

            fig, ax = plt.subplots(figsize=(15, 15))
            plt.rc('font', size=30)
            plt.rc('lines', linewidth=2)
            ax.add_collection(lc)

            # Scatter the individual points (small dots for all, larger for marked times)

            for t in time:

                xt = df.loc[((df['peptide'] == p) & (df['labels'] == l) & (df['apo_x'] == t)), 'norm_apo_noneq']
                yt = df.loc[((df['peptide'] == p) & (df['labels'] == l) & (df['apo_x'] == t)), 'norm_eq_noneq']

                sc = ax.scatter(xt, yt, c=t, cmap=cmap, norm=norm, edgecolors='k', s=50, zorder=11)


            # Add colorbar
            cbar = plt.colorbar(sc, ax=ax)
            cbar.set_label('Time (log scale)')

            # Graph aesthetics
            ax.set_xlim(-1.1, 1.1)
            ax.set_ylim(-1.1, 1.1)
            ax.axhline(0, color='grey', alpha=0.13)
            ax.axvline(0, color='grey', alpha=0.13)
            ax.set_xlabel('Difference in gradient from apo', fontsize=28)
            ax.set_ylabel('Difference in gradient from eq', fontsize=28)
            ax.set_title(f'{startend} - {p}')
            

            plt.savefig(path.join(label_folders, f'{startend} - {p}.png'))
            plt.close()

    print('plot_individual_colored Complete!')



#def plot_3d(df):

    #for l in unique_labels:
        
        #ig = plt.figure()
        #ax = plt.axes(projection='3d')

        #for p in unique_peptides:

            #xs = df.loc[((df['peptide'] == p) & (df['labels'] == l)), 'norm_apo_noneq']
            #zs = df.loc[((df['peptide'] == p) & (df['labels'] == l)), 'norm_eq_noneq']
            #ys = df.loc[((df['peptide'] == p) & (df['labels'] == l)), 'apo_time']

            #scatter_plot = ax.scatter3D(xs,np.log10(ys),zs, color = 'red', cmap = 'autumn', s = 0.1)

            #ax.plot3D(xs,np.log10(ys),zs)

        #ax.set_xlim(-1.1,1.1)
        #ax.set_zlim(-1.1,1.1)

        #ax.set_xlabel('difference from apo')
        #ax.set_zlabel('difference from eq')
        #ax.set_ylabel('time')

        #plt.show()



















