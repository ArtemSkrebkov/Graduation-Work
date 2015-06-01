import rpy2.robjects as robjects
import numpy as np
from scipy.spatial.distance import pdist
from scipy.cluster import hierarchy
import matplotlib.pyplot as plt
import math

def init():
	robjects.r.library("Biobase")
	robjects.r.library("GEOquery")	

def getMetaStrInfo(softFile, info):
	result = robjects.r["$"](robjects.r.Meta(softFile), info)

	return result[0]

def getDataTable(softFile):
	sampleCount = int(robjects.r["$"](robjects.r.Meta(softFile), "sample_count")[0])
	featureCount = int(robjects.r["$"](robjects.r.Meta(softFile), "feature_count")[0])
	shift = 2
	
	table = robjects.r.Table(softFile)
	colNames = robjects.r.colnames(table)

	sampleNames = []
	for i in range(shift, shift + sampleCount):
		sampleNames.append(colNames[i])

	featureList = []
	for i in range(0, sampleCount):
		featureList.append([])
		for j in range(0, featureCount):
			try:
				featureList[i].append(float(table[i + shift][j]))
			except ValueError:
				featureList[i].append(0.0)

	data = [map(float, featureList[i]) for i in range(0, sampleCount)]

	return sampleNames, data

def getGeneList(softFile, repeate = False):
	sampleCount = int(robjects.r["$"](robjects.r.Meta(softFile), "sample_count")[0])
	featureCount = int(robjects.r["$"](robjects.r.Meta(softFile), "feature_count")[0])
	shift = 4

	table = robjects.r.Table(softFile)
	geneList = []
	if repeate:
		#print(table[1].levels[table[1][featureCount - 1] - 1])
		genes = list(table[1].levels)
		indexes = list(table[1])
		for i in range(0, featureCount):
			geneList.append(genes[indexes[i] - 1])
	else:
		geneList = list(table[1].levels)
	
	return geneList

def norm(data):
    matrix = np.array(data, 'f')
    len_val = len(matrix[1, :])
    for i in range(len_val):
        local_min = matrix[:, i].min()
        if local_min !=  0.0:
            matrix[:, i] -= local_min
        local_max = matrix[:, i].max()
        if local_max !=  0.0:
            matrix[:, i] /= local_max
    return matrix.tolist()

def hierarchyDraw(Z, labels):
    plt.figure()
    hierarchy.dendrogram(Z, labels=labels, leaf_font_size=8, count_sort=True)
    plt.show()

def getClusters(samples):
	n = len(samples)

	dist = pdist(samples, 'euclidean')

	Z = hierarchy.linkage(dist, method='average')

	clusters = []
	for i in range(0, n):
		clusters.append([i])

	for i in range(0, len(Z)):
		clusters.append(clusters[int(Z[i][0])] + clusters[int(Z[i][1])])

	return clusters

def getListTopGenes(samples, genes, clusters, numCluster1, numCluster2, countTopGenes = 20):
	featureCount = len(samples[0])

	countClusters = len(clusters)
	maxClusterDist = [0] * featureCount
	if numCluster1 != numCluster2:
		ni = len(clusters[numCluster1])
		nj = len(clusters[numCluster2])
		for l in range(0, ni):
			si = samples[clusters[numCluster1][l]]
			for k in range(0, nj):
				sj = samples[clusters[numCluster2][k]]
				if (clusters[numCluster2][k] != samples[clusters[numCluster1][l]]):
					for q in range(0, featureCount):
						m = math.fabs(si[q] - sj[q])
						if m > maxClusterDist[q]:
							maxClusterDist[q] = m
	print(len(genes))
	print(len(maxClusterDist))		
	genesMaxDist = zip(genes, maxClusterDist)
	print(genesMaxDist[0])
	genesMaxDist = sorted(genesMaxDist, key = lambda x: x[1], reverse = True)
	print(genesMaxDist[0])
	return genesMaxDist[0:countTopGenes]

def drawClusters(samples, sampleNames, dataIsNorm):
	if dataIsNorm:
		samples = norm(samples)

	dist = pdist(samples, 'euclidean')
	Z = hierarchy.linkage( dist, method='average' )

	hierarchyDraw(Z, sampleNames)

def hierarchyTest():
	gds4847 = robjects.r.getGEO(filename = "./samples/GDS4847_full.soft.gz")
	sampleNames1, data1 = getDataTable(gds4847)
	
	gds4840 = robjects.r.getGEO(filename = "./samples/GDS4840_full.soft.gz")
	sampleNames2, data2 = getDataTable(gds4840)
	
	data = norm(data1) + norm(data2)
	sampleNames = sampleNames1 + sampleNames2
	
	dist = pdist(data2, 'euclidean')
	
	Z = hierarchy.linkage( dist, method='average' )

	hierarchyDraw(Z, sampleNames2)

init()

gds4840 = robjects.r.getGEO(filename = "./samples/GDS1649_full.soft.gz")
#geneList = getGeneList(gds4840)
sampleNames, samples = getDataTable(gds4840)
#genes = getGeneList(gds4840, True)
#clusters = getClusters(samples)
#topGenes = getListTopGenes(samples, genes, clusters, 6, 7, 100)
#print(clusters)
#print(topGenes)
#f = open('topGenes.txt', 'w')
#for i in range(len(topGenes)):
#	f.write(topGenes[i][0] + "\n")

drawClusters(samples, sampleNames, False)
#f = open('genes_gds4847', 'w')
#for gene in geneList:
#	f.write("%s\n" % gene)
#hierarchyTest()
