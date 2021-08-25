#For vlm need to save linear/non-linear pca,t-SNE,UMAP,  and true pseudotime/gamma, cluster labels

#git clone https://github.com/pachterlab/GFCP_2021.git
#import vis 
#vis.makeEmbeds()....


makeEmbeds(vlm,embeds)
'''
Save embedding objects in embeds lists in vlm

Parameters
----------
vlm : Velocyto loompy object
embeds : List of embeddings e.g. ['PCA','UMAP','t-SNE']

Returns
-------
'''


princCurvePlots(vlm)
'''
Plot principle curve coordinates for linear PCA embedding

Parameters
----------

Returns
-------
'''

gridArrowPlots(vlm,trans,embeds,sim=False)
'''
Plot arrow embeddings for vlm data with defined count transformations

Parameters
----------
sim: boolean to look for true pseudotime or cluster

Returns
-------




'''

getJaccard(embed1,embed2)
'''
Get jaccard distance between embeddings

Parameters
----------

Returns
-------
'''


jaccardPlots(vlm,pairs,n_neigh)
'''
Plot jaccard distances for neighbors between pairs of embeddings

Parameters
----------

Returns
-------

'''

angleDevPlots(vlm,trans,n_neigh,embed)
'''
Plot angle deviations from transformations over varying neighbors for embedding

Parameters
----------

Returns
-------
'''



getImputed(vlm)
'''
Get gamma inference from imputed counts

Parameters
----------

Returns
-------

'''

#Get imputed counts (S_x etc)

#Get gamma inference


plotPhase(vlm)
'''
Plot phase portrait

Parameters
----------

Returns
-------
'''

plotGammaK(vlm,gene_idx)
'''
Plot gamma over k neighbors for gene at gene_idx

Parameters
----------

Returns
-------

'''

phasePlots(vlm,n_neighs,genes)
'''
Plot phase portrais with gamma distributions for various genes across n_neighs

Parameters
----------

Returns
-------
'''

#For genes 1...n

#For k 

#Get inferred parameters

#Plot phase

#Plot inferred gamma + line for upregulate cells





