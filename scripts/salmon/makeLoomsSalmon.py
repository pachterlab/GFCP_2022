import numpy as np
import pandas as pd
import scanpy as sc
import anndata
#import scvelo as scv
#import matplotlib
import scipy
import json
import os

mm10_dirs = ["brain_10x_5k","brain_nuc_10x_5k","desai_dmso_v2","desai_idu_v2",
"heart_10k_v3","heart_1k_v3","neuron_10k_v3","neuron_1k_v3"]

hg38_dirs = ["pbmc_10k_v3","pbmc_1k_v3"]

for i in mm10_dirs:

	frydir = i+"_res"
	e2n_path = "geneid_to_name.txt"
	fpath = os.path.sep.join([frydir, "quant.json"])
	#if !os.path.exists(fpath):
	#  fpath = os.path.sep.join([frydir, "meta_info.json"])

	meta_info = json.load(open(fpath))
	ng = meta_info['num_genes']
	usa_mode = meta_info['usa_mode']

	if usa_mode:
		print("processing input in USA mode, will return A+S as the spliced count, and U as the unspliced count")
	else:
		print("please follow previous steps to generate the ount matrix in the USA mode")
		assert(False)

	af_raw = sc.read_mtx(os.path.sep.join([frydir, "alevin", "quants_mat.mtx"]))
	ng = int(ng/3)
	e2n = dict([ l.rstrip().split() for l in open(e2n_path).readlines()])
	var_names = [ l.rstrip() for l in open(os.path.sep.join([frydir, "alevin", "quants_mat_cols.txt"])).readlines()][:ng]
	var_names = [e2n[e] for e in var_names]

	obs_names = [ l.rstrip() for l in open(os.path.sep.join([frydir, "alevin", "quants_mat_rows.txt"])).readlines() ]

	x = af_raw.X
	spliced = x[:,range(0,ng)] + x[:,range(2*ng,3*ng)]
	unspliced = x[:,range(ng, 2*ng)]

	# create AnnData using spliced and unspliced count matrix
	adata = anndata.AnnData(X = spliced, layers = dict(spliced = spliced, unspliced = unspliced),
		obs =  pd.DataFrame({'barcode': np.array(obs_names)},index=np.array(obs_names)),
		var = pd.DataFrame({'gene_name': np.array(var_names)},index=np.array(var_names)))
	#adata.obs_names = pd.DataFrame({'barcode': np.array(barcodes)},index=np.array(obs_names))
	#adata.var_names = pd.DataFrame({'Gene': np.array(geneNames)},index=np.array(geneNames))
	adata.var_names_make_unique()

	fname = i+'.loom'

	adata.write_loom(fname)

for i in hg38_dirs:

	frydir = i+"_res"
	e2n_path = "geneid_to_name_hg.txt"
	fpath = os.path.sep.join([frydir, "quant.json"])
	#if !os.path.exists(fpath):
	#  fpath = os.path.sep.join([frydir, "meta_info.json"])

	meta_info = json.load(open(fpath))
	ng = meta_info['num_genes']
	usa_mode = meta_info['usa_mode']

	if usa_mode:
		print("processing input in USA mode, will return A+S as the spliced count, and U as the unspliced count")
	else:
		print("please follow previous steps to generate the ount matrix in the USA mode")
		assert(False)

	af_raw = sc.read_mtx(os.path.sep.join([frydir, "alevin", "quants_mat.mtx"]))
	ng = int(ng/3)
	e2n = dict([ l.rstrip().split() for l in open(e2n_path).readlines()])
	var_names = [ l.rstrip() for l in open(os.path.sep.join([frydir, "alevin", "quants_mat_cols.txt"])).readlines()][:ng]
	var_names = [e2n[e] for e in var_names]

	obs_names = [ l.rstrip() for l in open(os.path.sep.join([frydir, "alevin", "quants_mat_rows.txt"])).readlines() ]

	x = af_raw.X
	spliced = x[:,range(0,ng)] + x[:,range(2*ng,3*ng)]
	unspliced = x[:,range(ng, 2*ng)]

	# create AnnData using spliced and unspliced count matrix
	adata = anndata.AnnData(X = spliced, layers = dict(spliced = spliced, unspliced = unspliced),
		obs =  pd.DataFrame({'barcode': np.array(obs_names)},index=np.array(obs_names)),
		var = pd.DataFrame({'gene_name': np.array(var_names)},index=np.array(var_names)))
	#adata.obs_names = pd.DataFrame({'barcode': np.array(barcodes)},index=np.array(obs_names))
	#adata.var_names = pd.DataFrame({'Gene': np.array(geneNames)},index=np.array(geneNames))
	adata.var_names_make_unique()

	fname = i+'.loom'

	adata.write_loom(fname)


