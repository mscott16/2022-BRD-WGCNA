library(WGCNA)
library(plyr)
library(dplyr)
library(viridis)
library(edgeR)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()
enableWGCNAThreads(16)

setwd("C:/Users/...")

#metadata
traitData = read.csv("2017_complete.csv")

#genecount data
count<-read.csv("2017_lcpm.csv")

datExpr0 = as.data.frame(t(count[,-c(1)]))
names(datExpr0) =count$X
rownames(datExpr0)=names(count)[-c(1)]

nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)

names(traitData)

expressionSamples = rownames(datExpr0)
traitRows = match(expressionSamples, traitData$ID)
datTraits = traitData[traitRows, -1]
rownames(datTraits) = traitData[traitRows, 1]
collectGarbage()

A=adjacency(t(datExpr0), type = "distance")
k = as.numeric(apply(A, 2, sum)) - 1
z.k = scale(k)

thresholdZ.k = -5
outlierColor = ifelse(z.k < thresholdZ.k, "red", "black")
sampleTree = hclust(as.dist(1 - A), method = "average")
traitColors = data.frame(numbers2colors(datTraits, signed = FALSE, colors = inferno(30)))
dimnames(traitColors)[[2]] = paste(names(datTraits), sep = "")
datColors = data.frame(outlierC = outlierColor, traitColors)

#Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree, groupLabels = names(datColors), colors = datColors, 
                    main = "Sample dendrogram and trait heatmap", addTextGuide = TRUE)


###-----------soft threshold-----------###
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft=pickSoftThreshold(datExpr0,powerVector=powers, networkType = "signed")
# Plot the results:
par(mfrow = c(1, 2))
# SFT index as a function of different powers
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 
     xlab = "Soft Threshold (power)", ylab = "SFT, signed R^2", type = "n", main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 
     labels = powers, col = "darkblue")
# this line corresponds to using an R^2 cut-off of h
abline(h = 0.8, col = "darkorange")
# Mean connectivity as a function of different powers
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], type = "n", xlab = "Soft Threshold (power)", 
     ylab = "Mean Connectivity", main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "darkblue")


###-----------block-wise network construction and module detection-----------###
bwnet2 = blockwiseModules(datExpr0, corType = "pearson", maxBlockSize = 12795, TOMType = "signed",
                          networkType = "signed", power = 8, minModuleSize = 30, mergeCutHeight = .25, 
                          numericLabels = TRUE, saveTOMs = TRUE, pamRespectsDendro = FALSE, saveTOMFileBase = "2017_TOM2")

moduleColorsBlockwise = bwnet2$colors

blockNumber = 1
plotDendroAndColors(bwnet2$dendrograms[[blockNumber]], moduleColorsBlockwise[bwnet2$blockGenes[[blockNumber]]], 
                    "Module colors", main = paste("Dendrogram and module colors in block", blockNumber), 
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)


###-----------module correlation-----------###
#Choose a module assignment
moduleColorsBRD = labels2colors(bwnet2$colors)

#Define numbers of genes and samples
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)

#Recalculate MEs
MEs0 = moduleEigengenes(datExpr0, moduleColorsBRD)$eigengenes
MEsBRD = orderMEs(MEs0)
modTraitCor = cor(MEsBRD, datTraits, use = "p")
modTraitP = corPvalueStudent(modTraitCor, nSamples)

#Correlations and p-values
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", 
                   sep = "")
dim(textMatrix) = dim(modTraitCor)
par(mar = c(6, 8.8, 3, 2.2))

#Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = modTraitCor, xLabels = names(datTraits), yLabels = names(MEsBRD), 
               ySymbols = names(MEsBRD), colorLabels = FALSE, colors = plasma(60), 
               textMatrix = textMatrix, setStdMargins = FALSE, 
               cex.text = 0.5, zlim = c(-1,1), main = paste("Module-trait relationships"))


names(datExpr0)[moduleColorsBRD=="lightgreen"]
names(datExpr0)[moduleColorsBRD=="yellowgreen"]
names(datExpr0)[moduleColorsBRD=="tan"]
names(datExpr0)[moduleColorsBRD=="purple"]


###----------Preservation model-----------###
counts_2019<-read.csv("2019_lcpm.csv")

datExpr2019 = as.data.frame(t(counts_2019[,-c(1)]))
names(datExpr2019)=counts_2019$X
rownames(datExpr2019)=names(counts_2019)[-c(1)]



setLabels = c("X_2017", "X_2019")
multiExpr=list(X_2017=list(data=datExpr0),
               X_2019=list(data=datExpr2019))
moduleColors2017=moduleColorsBRD
multiColor=list(X_2017=moduleColors2017)

nPermutations1=200
set.seed(1)

mp = modulePreservation(multiExpr, multiColor,
                        referenceNetworks = 1, nPermutations = nPermutations1,
                        randomSeed = 1, quickCor = 0, verbose = 3)

ref=1
test = 2

Obs.PreservationStats=mp$preservation$observed[[ref]][[test]]
Z.PreservationStats=mp$preservation$Z[[ref]][[test]]

#Look at the observed preservation statistics
Obs.PreservationStats

#Z statistics from the permutation test analysis
Z.PreservationStats
modColors = rownames(Obs.PreservationStats)
moduleSize = Obs.PreservationStats$moduleSize

#Omit the grey module (background genes)
selectModules = !(modColors %in% c("grey"))

#Text labels for points
point.label = modColors[selectModules]

#Composite preservation statistics
medianRank=Obs.PreservationStats$medianRank.pres
Zsummary=Z.PreservationStats$Zsummary.pres

write.csv(Obs.PreservationStats, "medianRank.csv")
write.csv(Z.PreservationStats, "Zsummary.csv")


par(mfrow=c(1,2),mar = c(4.5,4.5,2.5,1))

#plot medianRank versus module size
plot(moduleSize[selectModules],medianRank[selectModules],col=1,
     bg=modColors[selectModules],pch = 21,main="medianRank Preservation",
     cex = 2, ylab ="medianRank",xlab="Module size", log="x")
labelPoints(moduleSize[selectModules],medianRank[selectModules],point.label,cex=1)
abline(h = 2, col = "blue")

#plot Zsummary versus module size
plot(moduleSize[selectModules],Zsummary[selectModules], col = 1,
     bg=modColors[selectModules],pch = 21,main="Zsummary Preservation",
     cex=2,ylab ="Zsummary", xlab = "Module size", log = "x")
labelPoints(moduleSize[selectModules],Zsummary[selectModules],point.label,cex=1)

#Add threshold lines for Zsummary
abline(h=0); abline(h=2, col = "blue", lty = 2); abline(h=10, col = "red", lty = 2)

plot(medianRank[selectModules],Zsummary[selectModules], col = 1,
     bg=modColors[selectModules],pch = 21,main="Preservation Analysis",
     cex=2,ylab ="Zsummary", xlab = "medianRank")

labelPoints(medianRank[selectModules],Zsummary[selectModules],point.label,cex=1)
abline(h=10, col = "black", lty = 4); abline(v=5, col = "black", lty = 2)

preserved<-read.csv("preservation.csv")

prevcol<-modColors[selectModules]

require("ggrepel")
ggplot(data=preserved, aes(x=medianRank, y = Zsummary, label=module)) + 
        geom_point(color = prevcol, size = 5) +
        geom_text_repel(max.overlaps = Inf, size=3) +
        theme_light() +
        geom_point(shape=1, size =5, color = "black") + 
        geom_hline(yintercept = 10, linetype = 4, size = 1) +
        geom_vline(xintercept = 5, linetype = 4, size = 1)

names(datExpr0)[moduleColorsBRD=="black"]
names(datExpr0)[moduleColorsBRD=="purple"]
names(datExpr0)[moduleColorsBRD=="lightgreen"]
names(datExpr0)[moduleColorsBRD=="tan"]
names(datExpr0)[moduleColorsBRD=="steelblue"]
names(datExpr0)[moduleColorsBRD=="orange"]
names(datExpr0)[moduleColorsBRD=="violet"]
names(datExpr0)[moduleColorsBRD=="royalblue"]

###-----------gene significance and module association-----------###
risk = as.data.frame(datTraits$Treament_Frequency)
names(risk) = "risk"
modNames = substring(names(MEsBRD), 3)

geneModuleMembership = as.data.frame(cor(datExpr0, MEsBRD, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(datExpr0, risk, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(risk), sep="")
names(GSPvalue) = paste("p.GS.", names(risk), sep="")

module = "violet"
column = match(module, modNames)
moduleGenes = moduleColorsBRD==module

sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Treatment Frequency",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

names(datExpr0)
names(datExpr0)[moduleColorsBRD=="violet"]

write.csv(geneModuleMembership, file="ModuleMembership.csv")

y=datTraits$Treament_Frequency
GS1=as.numeric(cor(y,datExpr0, use="p"))
GeneSignificance=abs(GS1)

datKME = signedKME(datExpr0, MEsBRD, outputColumnName = "MM.")
head(datKME)

FilterGenes= abs(GS1)> .3 & abs(datKME$MM.violet)>.7
table(FilterGenes)
dimnames(data.frame(datExpr0))[[2]][FilterGenes]
