# -*- coding: utf-8 -*-

from __future__ import print_function

from sklearn.datasets import make_blobs
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

import pandas as pd 
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler, RobustScaler

import matplotlib.pyplot as plt
from matplotlib import style
style.use("ggplot")
import matplotlib.patches as mpatches
		
from scipy.spatial import distance



def get_distance1(x1):
    boo = centroids[0]
    return distance.euclidean(x1, boo)

def get_distance2(x1):
    boo = centroids[1]
    return distance.euclidean(x1, boo)

def get_distance3(x1):
    boo = centroids[2]
    return distance.euclidean(x1, boo)


data = pd.read_csv("cluster.txt" , sep="\t" , index_col=0 , header=None )


#print (data.head())


data1 = data.dropna()



data2 = data1.astype(int)

#standard_scaler = StandardScaler()

#Norm_data = standard_scaler.fit_transform(data2)

#print (Norm_data[0:,])

#data2.to_csv("test.txt" , sep='\t')


pca = PCA(n_components=2)



new_pca = pca.fit_transform(data2)
##
new_pca = pd.DataFrame(new_pca)
new_pca.index = data2.index
new_pca.columns = ['PC1' , 'PC2']


new_pca.plot.scatter(x= 'PC1' , y = 'PC2')
plt.tight_layout()
plt.savefig('PCA_kmeans.pdf')


## convert to numpy array 

X = new_pca.as_matrix()


kmeans = KMeans(n_clusters=2, init='k-means++').fit(X)

centroids = kmeans.cluster_centers_
labels = kmeans.labels_
distances = kmeans.inertia_

print(centroids)
print(labels)


     
      

    


colors = ["g.","r.","c.","y.", 'k.' , 'm.']

for i in range(len(X)):
    #print("coordinate:",X[i], "label:", labels[i])
    plt.plot(X[i][0], X[i][1],colors[labels[i]], markersize = 10)


plt.scatter(centroids[:, 0],centroids[:, 1], marker = "x", s=150, linewidths = 5, zorder = 10)
plt.tight_layout()
plt.savefig("kmeans_groups.pdf")

results = pd.DataFrame([data2.index,labels, X]).T

results.columns = ['id' , 'group' , 'X']

#print (results)

#esults['dis_c1'] = 

#for i in range(len(centroids)):
#    col =  ("centroid " + str(i))
#    value = centroids[i]
results['Dis_C1'] = results.X.apply(get_distance1)
results['Dis_C2'] = results.X.apply(get_distance2)
#results['Dis_C3'] = results.X.apply(get_distance3)



#    for t in range(len(X)):
#        
#        dst = get_distance(X[t], value)
#        cid = data2.index[i]
#        print (str(cid) + " " + str(dst)) 







results.to_csv("kmeans_results.txt", sep='\t' , header=None , index=None)
#
##print (new_pca.head)
#
##
###print(pca.explained_variance_ratio_)
##
#
#
### cluster some shit
#k_means = cluster.KMeans(n_clusters=3,)
##
#k_means.fit(new_pca)
#
##donkey1 = pd.DataFrame([data2.index,donkey]).T
#
##donkey1.to_csv("transformed_data.txt" , sep='\t' , header = None , index=None)
#
##print (donkey1)
#
##print (new_pca)
#centroids = k_means.cluster_centers_
#
##print (centroids)
#
##numpy.savetxt("cluster_centers.txt" , centroids , delimiter="\t")
#
#
#
##print (centroids[:, 0])
#
#labels = k_means.labels_


## save some shit
#results = pd.DataFrame([data2.index,labels]).T
#results_data = pd.DataFrame([data2.index, Norm_data,]).T 


#results.to_csv("kmeans_results.txt", sep='\t' , header=None , index=None)
#results_data.to_csv("nor_data.txt", sep='\t' , header=None , index=None)
#print (results)

## Plot some shit 
#
#
#colors = ["green" , "red" , "blue"]
#
#
#
#
#
y =  labels                  
#                   
range_n_clusters = [2, 3, 4, 5, 6]

for n_clusters in range_n_clusters:
    # Create a subplot with 1 row and 2 columns
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.set_size_inches(18, 7)

    # The 1st subplot is the silhouette plot
    # The silhouette coefficient can range from -1, 1 but in this example all
    # lie within [-0.1, 1]
    ax1.set_xlim([-0.1, 1])
    # The (n_clusters+1)*10 is for inserting blank space between silhouette
    # plots of individual clusters, to demarcate them clearly.
    ax1.set_ylim([0, len(X) + (n_clusters + 1) * 10])

    # Initialize the clusterer with n_clusters value and a random generator
    # seed of 10 for reproducibility.
    clusterer = KMeans(n_clusters=n_clusters, random_state=10)
    cluster_labels = clusterer.fit_predict(X)

    # The silhouette_score gives the average value for all the samples.
    # This gives a perspective into the density and separation of the formed
    # clusters
    silhouette_avg = silhouette_score(X, cluster_labels)
    print("For n_clusters =", n_clusters,
          "The average silhouette_score is :", silhouette_avg)

    # Compute the silhouette scores for each sample
    sample_silhouette_values = silhouette_samples(X, cluster_labels)

    y_lower = 10
    for i in range(n_clusters):
        # Aggregate the silhouette scores for samples belonging to
        # cluster i, and sort them
        ith_cluster_silhouette_values = \
            sample_silhouette_values[cluster_labels == i]

        ith_cluster_silhouette_values.sort()

        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i

        color = cm.spectral(float(i) / n_clusters)
        ax1.fill_betweenx(np.arange(y_lower, y_upper),
                          0, ith_cluster_silhouette_values,
                          facecolor=color, edgecolor=color, alpha=0.7)

        # Label the silhouette plots with their cluster numbers at the middle
        ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

        # Compute the new y_lower for next plot
        y_lower = y_upper + 10  # 10 for the 0 samples

    ax1.set_title("The silhouette plot for the various clusters.")
    ax1.set_xlabel("The silhouette coefficient values")
    ax1.set_ylabel("Cluster label")

    # The vertical line for average silhouette score of all the values
    ax1.axvline(x=silhouette_avg, color="red", linestyle="--")

    ax1.set_yticks([])  # Clear the yaxis labels / ticks
    ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])

    # 2nd Plot showing the actual clusters formed
    colors = cm.spectral(cluster_labels.astype(float) / n_clusters)
    ax2.scatter(X[:, 0], X[:, 1], marker='.', s=30, lw=0, alpha=0.7,
                c=colors)

    # Labeling the clusters
    centers = clusterer.cluster_centers_
    # Draw white circles at cluster centers
    ax2.scatter(centers[:, 0], centers[:, 1],
                marker='o', c="white", alpha=1, s=200)

    for i, c in enumerate(centers):
        ax2.scatter(c[0], c[1], marker='$%d$' % i, alpha=1, s=50)

    ax2.set_title("The visualization of the clustered data.")
    ax2.set_xlabel("Feature space for the 1st feature")
    ax2.set_ylabel("Feature space for the 2nd feature")

    #plt.suptitle(("Silhouette analysis for KMeans clustering on sample data "
#                  "with n_clusters = %d" % n_clusters),
#                 fontsize=8, fontweight='bold')
    
    name = ("silhouette_" + str(n_clusters) + ".pdf")
    plt.tight_layout()
    plt.savefig(name)