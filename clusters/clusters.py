import rpy2.robjects as robjects
import numpy as np
from scipy.spatial.distance import pdist
from scipy.cluster import hierarchy
import matplotlib.pyplot as plt

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
			featureList[i].append(table[i + shift][j])
	data = [map(float, featureList[i]) for i in range(0, sampleCount)]

	return sampleNames, data

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

init()

gds4847 = robjects.r.getGEO(filename = "./samples/GDS4847_full.soft.gz")
sampleNames1, data1 = getDataTable(gds4847)

gds4840 = robjects.r.getGEO(filename = "./samples/GDS4840_full.soft.gz")
sampleNames2, data2 = getDataTable(gds4840)

data = norm(data1) + norm(data2)
sampleNames = sampleNames1 + sampleNames2

dist = pdist(data, 'euclidean')

Z = hierarchy.linkage( dist, method='average' )

hierarchyDraw(Z, sampleNames)
