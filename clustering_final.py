import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os
from os import path
from reattemptreal import (gen_uniques, add_columns)
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, silhouette_samples

df = pd.read_csv('raw/State_data_3_states.csv')
df = add_columns(df)
unique_sequences, unique_states, unique_startend = gen_uniques(df)

outpath = 'out'


def process_cluster_data(df):

    delta_dfs_list = []

    peptides_remove = []

    for p in unique_sequences:

        remove_peptide = False

        apo_y = df.loc[((df['Sequence'] == p) & (df['State'] == unique_states[0])), 'Uptake']
        eq_y = df.loc[((df['Sequence'] == p) & (df['State'] == unique_states[1])), 'Uptake']
        noneq_y = df.loc[((df['Sequence'] == p) & (df['State'] == unique_states[2])), 'Uptake']

        apo_x = df.loc[((df['Sequence'] == p) & (df['State'] == unique_states[0])), 'Exposure']
        eq_x = df.loc[((df['Sequence'] == p) & (df['State'] == unique_states[1])), 'Exposure']
        noneq_x = df.loc[((df['Sequence'] == p) & (df['State'] == unique_states[2])), 'Exposure']

        mapo_y = df.loc[((df['Sequence'] == p) & (df['State'] == unique_states[0])), 'MaxUptake']
        meq_y = df.loc[((df['Sequence'] == p) & (df['State'] == unique_states[1])), 'MaxUptake']
        mnoneq_y = df.loc[((df['Sequence'] == p) & (df['State'] == unique_states[2])), 'MaxUptake']

        norm_apo_y = apo_y/mapo_y
        norm_eq_y = eq_y/meq_y
        norm_noneq_y = noneq_y/mnoneq_y

        startend = df.loc[(df['Sequence'] == p), 'StartEnd']
        se = pd.Series(startend)
        se1 = pd.unique(se)

        apo_sumy = sum(norm_apo_y)/len(apo_x)
        eq_sumy = sum(norm_eq_y)/len(eq_x)
        noneq_sumy = sum(norm_noneq_y)/len(noneq_x)

        noneq_apo = noneq_sumy - apo_sumy
        noneq_eq = noneq_sumy - eq_sumy

        delta_df = pd.DataFrame(data = [], columns = ['startend','peptide','noneq-apo', 'noneq-eq'])
        delta_df['startend'] = se1
        delta_df['peptide'] = p
        delta_df['noneq-apo'] = noneq_apo
        delta_df['noneq-eq'] = noneq_eq

        delta_dfs_list.append(delta_df)

        for s in unique_states:

             max_u = df.loc[((df['Sequence'] == p) & (df['State'] == s)), 'MaxUptake']
             e = df.loc[((df['Sequence'] == p) & (df['State'] == s)), 'Exposure']

             if max_u.empty or e.empty:

                remove_peptide = True
                break
             
             if remove_peptide:
                peptides_remove.append(p)

    delta_dfs = pd.concat(delta_dfs_list)

    delta_dfs = delta_dfs[delta_dfs['peptide'] != 'IWGVEPSRQRLPAPDEKIP']

    delta_dfs = delta_dfs[~delta_dfs['peptide'].isin(peptides_remove)]

    max_noneq_apo = delta_dfs['noneq-apo'].abs().max()
    max_noneq_eq = delta_dfs['noneq-eq'].abs().max()

    delta_dfs['norm noneq-apo'] = delta_dfs['noneq-apo']/max_noneq_apo
    delta_dfs['norm noneq-eq'] = delta_dfs['noneq-eq']/max_noneq_eq

    return(delta_dfs)

prints = process_cluster_data(df)

def choose_k(df, outpath, max_clusters = 12):
    

    k_plots_folder = os.path.join(outpath, 'k plots')
    os.makedirs(k_plots_folder, exist_ok = True)

    #select only these columns from df
    Clusterdata = df[['norm noneq-apo', 'norm noneq-eq']]

    silhouette_scores = []
    inertia_scores = []

    #loop through different numbers of clusters
    for i in range(2,max_clusters):

        #create k means clustering for selected no of clusters
        kmeans = KMeans(n_clusters = i, init ='k-means++', max_iter = 1000).fit(Clusterdata)

        #generate silhouette score for tnis number of clusters
        silhouette_avg = silhouette_score(Clusterdata, kmeans.labels_)
        
        inertia_scores.append(kmeans.inertia_)
        silhouette_scores.append(silhouette_avg)

    # set x-axis for both graphs
    x = list(range(2,max_clusters))

    # plot graphs
    plt.plot(x,silhouette_scores, marker = 'o', linewidth=2, color='cadetblue')
    plt.title('The Silhouette method for optimal k')
    plt.xlabel('Number of clusters')
    plt.ylabel('Average Silhouette score')
    plt.xticks(ticks=range(1, max_clusters + 1))
    plt.savefig(path.join(k_plots_folder, "The Silhouette method for optimal k.png"))
    
    plt.show()

    plt.plot(x,inertia_scores, marker = 'o', linewidth=2, color='cadetblue')
    plt.xlabel('Number of clusters')
    plt.ylabel('WCSS')
    plt.title('The elbow method for optimal k')
    plt.savefig(path.join(k_plots_folder, "The elbow method for optimal k.png"))
    
    #plt.show()

    return silhouette_scores

#plot = choose_k(prints)

def gen_clusterdata(df: pd.DataFrame, outpath, n_clusters = 7):
    """
    Makes a scatter plot
    Clusters_r is a value of clusters"""

    Clusterdata = pd.DataFrame(data = [], columns = ['normnoneq_apo', 'normnoneq_eq'])

    Clusterdata['normnoneq_apo'] = df['norm noneq-apo']
    Clusterdata['normnoneq_eq'] = df['norm noneq-eq']

    kmeans = KMeans(n_clusters = n_clusters, init ='k-means++', random_state= 10).fit(Clusterdata)

    Clusterdata['labels'] = kmeans.fit_predict(Clusterdata[['normnoneq_apo', 'normnoneq_eq']])

    Clusterdata['peptide'] = df['peptide']
    
    unique_labels = Clusterdata['labels'].unique()
    unique_labels = np.sort(unique_labels)

    # plot centroids on the graph
    Centroids = kmeans.cluster_centers_

    centroids = kmeans.cluster_centers_
    Centroids = pd.DataFrame(centroids)

    cen_x = [i[0] for i in centroids]
    cen_y = [i[1] for i in centroids]
    Clusterdata['cen_x'] = Clusterdata.labels.map({i: cen_x[i] for i in range(n_clusters)})
    Clusterdata['cen_y'] = Clusterdata.labels.map({i: cen_y[i] for i in range(n_clusters)})

    Clusterdata.to_csv(path.join(outpath, 'clusterdata.csv'))

    return(Clusterdata, unique_labels, centroids)



def gen_scatter_plot(Clusterdata, unique_labels, centroids, outpath):

    
    plt.figure(figsize=(15, 15))
    plt.rc('font', size=30)
    plt.rc('lines', linewidth=2)

    for i, l in enumerate(unique_labels):
        
        x = Clusterdata.loc[(Clusterdata['labels'] == l), 'normnoneq_apo']
        y = Clusterdata.loc[(Clusterdata['labels'] == l), 'normnoneq_eq']

        plt.scatter(x, y, label = l, edgecolors='grey', s=450, alpha=0.3)

        for idx, val in Clusterdata.iterrows():
            x = [val.normnoneq_apo, val.cen_x]
            y = [val.normnoneq_eq, val.cen_y]
            plt.plot(x, y, c='grey', alpha=0.03, linewidth=1)

        plt.annotate(l, xy=(centroids[i]), family='sans-serif',
                     bbox=dict(boxstyle="circle", pad=0.01, facecolor='white', alpha=0.5),
                     horizontalalignment='center', verticalalignment='center', size=35, weight='bold', color='slategrey')

    plt.legend()


    plt.axhline(0, color='grey', alpha=0.13)
    plt.axvline(0, color='grey', alpha=0.13)
    plt.xlabel('Difference from apo', fontsize=28)
    plt.ylabel('Difference from eq', fontsize=28)
    plt.xlim(-1.1,1.1)
    plt.ylim(-1.1,1.1)

    #plt.show()
    plt.savefig(path.join(outpath, 'scatter plot.png'))
    plt.close()

    return

#plotting = gen_scatter_plot(prints)
    




        




    



