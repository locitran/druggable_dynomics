import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist
from sklearn.metrics import silhouette_score
from sklearn.metrics import silhouette_samples
from scipy.cluster.hierarchy import average, fcluster



'''
The size of the input data(matrix) should be (n_samples, n_features)
'''

###############################################################################
### Do hierarchical clustering.                                             ###
###############################################################################
def hierarchical_clustering(data, n_clust):

    ### calculate pairwise distance(pdist()) and do linkage(average())
    Z = average(pdist(data))
    cluster_ids = fcluster(Z, n_clust, criterion='maxclust')

    return Z, cluster_ids


###############################################################################
### Try 2-15 clusters and calculate the average silhouette score of each n. ###
### Return unique clustering results and the silhouette score of each n.    ###
###############################################################################
def try_n_clusters(data, minclust, maxclust):
        
    maxclust += 1
    ### check min/ max cluster to make sure 2 <= no. of clusters <= n_sample-1
    assert (maxclust <= data.shape[0] and minclust >= 2), \
        "maxclust needs to be smaller then (n_samples-1) and minclust needs to be bigger then 1"
        
    ### store the coresponding cluster ids for each number of clusters in m  
    m = []
    cid_arr = []
    
    for i in range(minclust, maxclust):
        Z, cluster_ids = hierarchical_clustering(data, i)
        cid_arr.append([[i],[cluster_ids]])
        n_clusters = len(set(cluster_ids))
        
    ### Try 2-15 clusters and store the unique clustering results in matrix m
        if i == minclust:
            m.append([n_clusters, cluster_ids])
        elif n_clusters != m[-1][0]:
            m.append([n_clusters, cluster_ids])
        else:
            continue
        
    silhouette_scores = []
    
    for j in range(len(m)):
    
        ### get the number of clusters in the j th clustering result
        n = m[j][0]
        ### get the j th clustering result
        cluster_ids = m[j][1]
        ### get the score for every sample using
        score_arr = silhouette_samples(data, cluster_ids)
        ### get average score using "cluster_ids" as label
        avg_score = silhouette_score(data, cluster_ids)
        # silhouette_scores.append([n, avg_score, score_arr]) #ChihHui
        silhouette_scores.append([n, avg_score]) #Navi
        print(silhouette_scores);

    return silhouette_scores, cid_arr



###############################################################################
### Read the result from try_n_clusters and see which n has the best        ###
### average silhouette score.                                               ###
### Return the optimal number of clusters(n) and the matrix\                ###
### "silhouette_scores", it stores all scores of each sample at different n ###
###############################################################################
def optimal_n_clustering(data, minclust=2, maxclust=15, method='scikit' ):
    

    silhouette_scores, cid_arr = np.array(try_n_clusters(data, minclust, maxclust), dtype=object)
    ### the index of the optimal number of cluster having the highest score
    max_id = np.argmax(silhouette_scores[:, 1])
    ### the optimal number of cluster which has the highest score
    optimal_n, score= silhouette_scores[max_id, 0], silhouette_scores[max_id, 1]        
    
    

    print(f'The optimal number of clusters is {optimal_n}, the average silhouette score is {score:.4f}.')
    
    cluster_ids = cid_arr[max_id][1]

    return optimal_n, cluster_ids, silhouette_scores

def do_clustering_base_on_silhouette_score(universe, pcs, pc_vars, filename, xray, output_selection, new_positions, superimpose_flag):

    def write_pdb_file(universe, frame, cluster_id, path, output_selection, filename, new_positions, superimpose_flag):
        universe.trajectory[frame]
        protein = universe.select_atoms(output_selection)
        if superimpose_flag:
            protein.atoms.positions = new_positions[frame]
        protein.atoms.write(os.path.join(path,f'{filename}_cluster{cluster_id}_{frame+1}.pdb'))

    path = os.path.dirname(filename)
    pc1, pc2 = pcs[0], pcs[1]

    if xray:
        data = np.array([pc1[:-1], pc2[:-1]]).T
        pc_pos = np.array([pc1[:-1], pc2[:-1]])
    else:
        data = np.array([pc1, pc2]).T
        pc_pos = np.array([pc1, pc2])

    ### data size: (n_samples, n_features)
    optimal_n, cluster_ids, silhouette_scores = optimal_n_clustering(data, 2, 10)
    n_clust = [i[0] for i in silhouette_scores[:, :2]]
    n_scores = [round(i[1], 4) for i in silhouette_scores[:, :2]]

    plt.figure(dpi=90)
    plt.plot(n_clust, n_scores, linewidth="2", linestyle="-", marker="o")
    plt.xlabel('number of cluster')
    plt.ylabel('silhouette score')
    for x, y in zip(n_clust, n_scores):
        plt.text(x-0.1, y+0.005, f'{y:.4f}')
    plt.savefig(f'{filename}_silhouette_score.png')

    #find cluster center and plot
    clusters = [i for i in set(cluster_ids[0])]
    
    pc1_var, pc2_var = pc_vars[0], pc_vars[1]

    color_labels = ['lightgreen', 'orange', 'cyan' , 'blue', 'red', 'pink', 'black', 'gray', 'purple', 'green']
    center_color = ['green', 'red', 'lightgreen', 'lightblue', 'pink', 'yellow', 'gray', 'black', 'blue', 'lightgreen']
    legend = []
    represent_ids = []

    plt.figure()
    for i in range(len(clusters)):
    ### plot each cluster
        clust_id = clusters[i]
        clust_pos = pc_pos.T[cluster_ids[0] == clust_id]
        print(f'cluster{clust_id}: {clust_pos.shape[0]} frames')
        plt.scatter(clust_pos.T[0], clust_pos.T[1], c=color_labels[i], marker='.')
        legend.append(f'cluster {clust_id}')

        ### plot centroid
        mean_pos = clust_pos.mean(axis=0)
        plt.scatter(mean_pos[0], mean_pos[1], marker='^', c=center_color[i])
        legend.append(f'centroid {clust_id}')
        represent_id = np.linalg.norm(pc_pos.T - mean_pos, axis=1).argmin()
        represent_ids.append(represent_id)
        print(f'cluster{clust_id}_representative:{represent_id+1}')

    if xray:
        plt.scatter(pc1[-1], pc2[-1], c=color_labels[i+1], marker='.')
        legend.append(f'xray structure')

    plt.legend(legend, loc='upper right', bbox_to_anchor=(1.30, 1.03))
    plt.xlabel(f'pc1 {pc1_var*100:.2f}%')
    plt.ylabel(f'pc2 {pc2_var*100:.2f}%')

    plt.tight_layout()
    
    plt.savefig(f'{filename}_PCA_plot.png')

    for i in clusters:
        write_pdb_file(universe,represent_ids[i-1], i, path, output_selection, filename, new_positions, superimpose_flag)
