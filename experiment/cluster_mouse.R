library(GEOquery)
library(Biobase)
library(cluster)
library(openxlsx)
library(dendextend)

geneSymbols <- vector()
genesID <- vector
namesSample <- vector()
typesTissue <- vector()

readGSESamples <- function(gsefilename, gplfilename, isGeneSymbol)
{
	gse <- getGEO(filename = gsefilename)

	probesets <- Table(GPLList(gse)[[1]])$ID

	data.matrix <- do.call('cbind',lapply(GSMList(gse),
							function(x)
							{
								tab <- Table(x)
								mymatch <- match(probesets,tab$ID_REF)
								return(tab$VALUE[mymatch])
							}))
	data.matrix <- apply(data.matrix,2,function(x) {as.numeric(as.character(x))})
	#data.matrix <- log2(data.matrix)
	rownames(data.matrix) <- probesets

	#geneSymbols <- vector()
	#save import information
	if (isGeneSymbol)
	{
		gpl <- getGEO(filename = gplfilename)
		geneSymbols <<- as.vector(Table(gpl)[,c("Gene symbol")])
	}
	else
	{
		geneSymbols <<- as.vector(Table(GPLList(gse)[[1]])$Symbol)
	}

	gsms <- GSMList(gse)
	for (i in 1:length(gsms))
	{
		metaInfo <- Meta(gsms[[i]])
		if ( length(grep("normal", metaInfo$source_name_ch1, ignore.case=TRUE)) != 0 ||
			 length(grep("non-tumor", metaInfo$source_name_ch1, ignore.case=TRUE)) != 0)
		{
			typesTissue <<- c(typesTissue, 1)
		}
		else
		{
			typesTissue <<- c(typesTissue, 2)
		}
	}
	genesID <<- as.vector(probesets)

	
	return (as.table(data.matrix[,]))
}

readGSESamplesWithGenes <- function(gsefilename, gplfilename, genes, isGeneSymbol)
{
	gse <- getGEO(filename = gsefilename)

	probesets <- Table(GPLList(gse)[[1]])$ID

	data.matrix <- do.call('cbind',lapply(GSMList(gse),
							function(x)
							{
								tab <- Table(x)
								mymatch <- match(probesets,tab$ID_REF)
								return(tab$VALUE[mymatch])
							}))
	data.matrix <- apply(data.matrix,2,function(x) {as.numeric(as.character(x))})
	#data.matrix <- log2(data.matrix)
	#data.matrix <- exp(data.matrix * log(2.0))
	rownames(data.matrix) <- probesets

	geneSymbols <- vector()
	if (isGeneSymbol)
	{
		gpl <- getGEO(filename = gplfilename)
		geneSymbols <- as.vector(Table(gpl)[,c("Gene symbol")])
	}
	else
	{
		geneSymbols <- as.vector(Table(GPLList(gse)[[1]])$Symbol)
	}

	gsms <- GSMList(gse)
	for (i in 1:length(gsms))
	{
		metaInfo <- Meta(gsms[[i]])
		if ( length(grep("normal", metaInfo$source_name_ch1, ignore.case=TRUE)) != 0 ||
			 length(grep("non-tumor", metaInfo$source_name_ch1, ignore.case=TRUE)) != 0)
		{
			typesTissue <<- c(typesTissue, 1)
		}
		else
		{
			typesTissue <<- c(typesTissue, 2)
		}
	}
	genesID <<- as.vector(probesets)

	numGenes <- vector()
	n <- length(genes)
	m <- length(geneSymbols)

	for (j in 1:n)
	{
		countСoincidence <- 0
		for (i in 1:m)
		{
			sumExpression <- 0.0
			if (tolower(geneSymbols[i]) == tolower(genes[j]))
			{
				countСoincidence <- countСoincidence + 1 
				if (countСoincidence == 1)
					numGenes <- c(numGenes, i)
				else
				{
					k <- length(numGenes)
					data.matrix[k,] <- data.matrix[k,] + data.matrix[i,]
				}
				#break
			}
		}
		k <- length(numGenes)
		if (countСoincidence != 0)
		{
			data.matrix[k,] <- data.matrix[k,] / countСoincidence
		}
	}

	return (as.table(data.matrix[numGenes,]))
}

GetAccuracy <- function(hier, countClusters, trueAnswers)
{
	dv2 <- as.table(cutree(as.hclust(hier), k = countClusters))

	namesSamples <- rownames(dv2)
	n <- length(namesSamples)
	countError <- 0
	for (i in 1:length(namesSamples))
	{
		if (dv2[namesSamples[i]] == trueAnswers[namesSamples[i]])
		{
			countError <- countError + 1
		}
	}

	return (1.0 - countError / n)
}

getGenesFromID <- function(gplfilename, curIDs)
{
	genes <- vector()

	gpl <- getGEO(filename = gplfilename)
	geneSymbols <- as.vector(Table(gpl)[,c("Gene symbol")])
	IDs <- as.vector(Table(gpl)[,c("ID")])

	n <- length(IDs)
	m <- length(curIDs)

	for (i in 1:m)
	{
		for (j in 1:n)
		{
			if (IDs[j] == curIDs[i])
			{
				genes <- c(genes, geneSymbols[j])
			}
		}
	}

	return (genes)
}

getIntersection <- function(filenameA, filenameB)
{
	fileA <- read.table(filenameA, header = TRUE, sep = "\t")
	fileB <- read.table(filenameB, header = TRUE, sep = "\t")

	genesIDA <- as.vector(fileA[,1])
	genesIDB <- as.vector(fileB[,1])
	genesA <- as.vector(fileA[,2])
	genesB <- as.vector(fileB[,2])
	nA <- length(fileA[,1]); nB <- length(fileB[,1]);

	interGeneName <- vector()
	interGeneIDA <- vector()
	interGeneIDB <- vector()
	for (i in 1:nA)
	{
		for (j in 1:nB)
		{
			if (genesA[i] != "" && genesB != "")
			{
				if (tolower(genesA[i]) == tolower(genesB[j]))
				{
					nInter <- length(interGeneName)
					flag <- TRUE
					if (nInter != 0)
					{
						for (k in 1:nInter)
						{
							if (tolower(genesB[j]) == tolower(interGeneName[k]))
							{
								flag <- FALSE; break;
							}
						}
					}
					if (flag)
					{
						interGeneName <- c(interGeneName, genesB[j])
						interGeneIDA <- c(interGeneIDA, genesIDA[i])
						interGeneIDB <- c(interGeneIDB, genesIDB[j])
						break
					}
				}
			}
		}
	}

	results <- data.frame(Gene.name = interGeneName, Gene.IDA = interGeneIDA, Gene.IDB = interGeneIDB)

	return (results)
}

getIntersectionWithRepeates <- function(filenameA, filenameB)
{
	fileA <- read.table(filenameA, header = TRUE, sep = "\t")
	fileB <- read.table(filenameB, header = TRUE, sep = "\t")

	genesIDA <- as.vector(fileA[,1])
	genesIDB <- as.vector(fileB[,1])
	genesA <- as.vector(fileA[,2])
	genesB <- as.vector(fileB[,2])
	nA <- length(fileA[,1]); nB <- length(fileB[,1]);

	interGeneName <- vector()
	interGeneIDA <- vector()
	interGeneIDB <- vector()
	for (i in 1:nA)
	{
		for (j in 1:nB)
		{
			if (genesA[i] != "" && genesB != "")
			{
				if (tolower(genesA[i]) == tolower(genesB[j]))
				{
					interGeneName <- c(interGeneName, genesB[j])
					interGeneIDA <- c(interGeneIDA, genesIDA[i])
					interGeneIDB <- c(interGeneIDB, genesIDB[j])
				}
			}
		}
	}

	results <- data.frame(Gene.name = interGeneName, Gene.IDA = interGeneIDA, Gene.IDB = interGeneIDB)

	return (results)
}

getBedFile <- function(db, host, organism, platform, genesIDS, filename)
{
	ensembl <- useMart(biomart=db, host = host, path="/biomart/martservice", dataset = organism)

	bed <- getBM(attributes = c('chromosome_name', 'start_position','end_position', platform, 'strand'), 
			     filters = platform, values = genesIDS, mart = ensembl)
	numberTrue <- vector()
	for (i in 1:length(bed$chromosome_name))
	{
		if (nchar(bed$chromosome_name[i]) <= 2)
			numberTrue <- c(numberTrue, i)
	}
	bed <- bed[numberTrue,]
	bed$chromosome_name <- paste("chr",bed$chromosome_name,sep="")
	bed$strand[bed$strand==1] <- "+"
	bed$strand[bed$strand==-1] <- "-"
	bed$score <- 0
	bed <- bed[,c('chromosome_name', 'start_position', 'end_position', platform, 'score', 'strand')]

	if (length(filename) != 0)
	{
		write.table(bed, filename, row.names = FALSE, col.names = FALSE)
	}

	return (bed)
}

HierAnalysis <- function(gsefilename, gplfilename, name = "", listgenes = vector())
{
	expressionGene <- 0
	if (length(listgenes) == 0)
	{
		expressionGene <- readGSESamples(gsefilename, gplfilename, TRUE)
	}
	else
	{
		expressionGene <- readGSESamplesWithGenes(gsefilename, gplfilename, listgenes, TRUE)
	}
	expressionGeneTranspose <- t(expressionGene)
	hier <- diana(expressionGeneTranspose)
	dend <- as.dendrogram(hier)
	
	y <- typesTissue
	par( cex = 0.8)
	colorDend <- color_branches(dend, k = length(y), col = y)

	coef <- coefHier(hier)

	title <- paste(name, "\nDivisive coef = " ,as.character(coef))
	plot(colorDend, main = title)
}

genes <- as.vector(read.table("/home/artem/projects/dimplom/tools/data/geo/experiment/intersec.txt", sep = "\t"))
genes <- as.vector(genes[, 1])
#expressionGene <- readGSESamples("/home/artem/projects/dimplom/tools/data/geo/experiment/mouse/GSE49200_family.soft", 
#								 "/home/artem/projects/dimplom/tools/data/geo/experiment/mouse/GPL81.soft", TRUE)
#expressionGene <- readGSESamples("/home/artem/projects/dimplom/tools/data/geo/experiment/human/GSE32863_family.soft", 
#								 "", FALSE)
#expressionGene <- readGSESamplesWithGenes("/home/artem/projects/dimplom/tools/data/geo/experiment/mouse/GSE49200_family.soft", 
#								 		  "/home/artem/projects/dimplom/tools/data/geo/experiment/mouse/GPL81.soft",
#								           genes, TRUE)
#expressionGene <- readGSESamplesWithGenes("/home/artem/projects/dimplom/tools/data/geo/experiment/human/GSE32863_family.soft", 
#								 		  "/home/artem/projects/dimplom/tools/data/geo/experiment/human/GPL6884.annot",
#								           genes, FALSE)

#expressionGeneTranspose <- t(expressionGene)
#file <- paste("/home/artem/", "human_last_article.xlsx")
#res <- write.xlsx(expressionGene, file)  
#write.csv(expressionGene, "/home/artem/human_last_article.csv")
#pvalues <- WelchTestForGeneExpression(expressionGeneTranspose, 2, 1, 2)
#sortedPvalues <- as.table(pvalues[,order(pvalues)])

#hier <- diana(expressionGeneTranspose)
#plot(hier, type = "fun")

#trueAnswers <- read.table("/home/artem/projects/dimplom/tools/data/geo/experiment/HumanClusters.txt", sep = ",")
#samplesNames <- trueAnswers[,1]
#answers <- trueAnswers[,2]
#names(answers) <- samplesNames
#print(GetAccuracy(hier, 2, answers))

#maxDiffGenes("/home/artem/projects/dimplom/tools/data/geo/experiment/mouse/GSE49200_family.soft", 
#			 "/home/artem/projects/dimplom/tools/data/geo/experiment/mouse/GPL81.soft", 
#			 2, 1, 2)

#HierAnalysis("/home/artem/projects/dimplom/tools/data/geo/experiment/mouse/GSE49200_family.soft", 
#	         "/home/artem/projects/dimplom/tools/data/geo/experiment/mouse/GPL81.soft", "Mouse genes")

HierAnalysis("/home/artem/projects/dimplom/tools/data/geo/experiment/human/GSE32863_family.soft", 
			 "/home/artem/projects/dimplom/tools/data/geo/experiment/human/GPL6884.annot", "Human intersection genes", genes)