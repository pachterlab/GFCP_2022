#For vlm need to save linear/non-linear pca,t-SNE,UMAP,  and true pseudotime/gamma, cluster labels

#!git clone https://github.com/pachterlab/GFCP_2021.git
#!cd GFCP_2021
#import vis 
#vis.makeEmbeds()....


makeEmbeds(vlm,embeds):
'''
Save embedding objects in embeds lists in vlm

Parameters
----------
vlm : Velocyto loompy object
embeds : List of embeddings e.g. ['PCA','UMAP','t-SNE']

Returns
-------
'''

getImputed(vlm):
'''
Get gamma inference from imputed counts

Parameters
----------

Returns
-------

'''

#Get imputed counts (S_x etc)

#Get gamma inference

getJaccard(embed1,embed2)
'''
Get jaccard distance between embeddings

Parameters
----------

Returns
-------
'''

# ---------------- Plotting -------------

princCurvePlots(ax,vlm):
'''
Plot principle curve coordinates for linear PCA embedding

Parameters
----------

Returns
-------
'''





plotEmbed(ax,vlm,embed):
'''
Plot given embedding (UMAP, t-SNE, etc)

Parameters
----------

Returns
-------
'''



plotGrid(ax,vlm,embed):
'''
Plot grid with arrows given embedding

Parameters
----------

Returns
-------
'''


gridArrowPlots(vlm,trans,embeds,sim=False):
'''
Plot arrow embeddings for vlm data with defined count transformations

Parameters
----------
sim: boolean to look for true pseudotime or cluster and/or principal curve

Returns
-------

'''



plotJaccard(ax,vlm,pair):
'''
Single jaccard distance plot
'''

jaccardPlots(vlm,pairs,n_neigh):
'''
Plot jaccard distances for neighbors between pairs of embeddings

Parameters
----------

Returns
-------

'''



plotTheta(ax,vlm,embed,baseline):
'''
Single angle deviation plot
'''

angleDevPlots(vlm,trans,n_neigh,embed,baseline):
'''
Plot angle deviations from transformations over varying neighbors for embedding (only compared to baseline)

Parameters
----------

Returns
-------
'''




plotPhase(ax,vlm):
'''
Plot phase portrait

Parameters
----------

Returns
-------
'''

plotGammaK(ax,vlm,gene_idx):
'''
Plot gamma over k neighbors for gene at gene_idx

Parameters
----------

Returns
-------

'''

phasePlots(vlm,n_neighs,genes):
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





