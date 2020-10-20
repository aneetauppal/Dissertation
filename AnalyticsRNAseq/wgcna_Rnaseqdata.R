#WGCNA used from dissertation work
#WGCNA written by WGCNA and adapted for RNAseq data 


# Code builds off of normalized count matrix produced by deseq2
library(dplyr)
library(tools)
library(DESeq2)
library("WGCNA")


#calculated by deseq2 from deseq object for use in wgcna 
#Get variance stabilizing transformation for WGCNA using count data
Wgcnaobj2=varianceStabilizingTransformation(countDataMatrix, blind = TRUE, fitType="parametric")

#remove genes that are not "known" Refseq's (XM_ XR_ XP) 
wgcnaob <- Wgcnaobj2[rownames(Wgcnaobj2) %like% "NM_",]
newresults <- t(wgcnaob)

options(stringsAsFactors = FALSE);
no.obs = 30


#check that all genes are good with no missing values 
gsg = goodSamplesGenes(newresults, verbose =3)
gsg$allOK

#not all ok so - 
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(newresults)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(newresults)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  newresults_g = newresults[gsg$goodSamples, gsg$goodGenes]
}
#check again 
#check that all genes are good with no missing values 
gsg2 = goodSamplesGenes(newresults_g, verbose =3)
gsg2$allOK
#allgood

#detect outliers
sampleTree = hclust(dist(newresults), method ="average")
#OR 
sampleTree = hclust(dists, method ="average")

plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

abline(h = 15, col = "red");

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
#did not remove any samples - i don't think there were any tell tale outliers
Tstable2 <- t(stable2)
traitrows = rownames(Wgcnaobj2t)

traitColors= numbers2colors(y, signed = FALSE)
datTraits = stable2[,3]


plotDendroAndColors(sampleTree, traitColors,
                    groupLabels = colnames(stable2)[2],
                    main = "Dendrogram clustered by samples", cex.colorLabels = 1.3, cex.dendroLabels = 1.3, 
                    cex.rowText = 1.3)


save(newresults, datTraits, file = "rnaseqoils-02-dataInput.RData")

allowWGCNAThreads()
lnames = load(file="rnaseqoils-02-dataInput.RData")
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(newresults, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#invoke when using WGCNA
cor <- WGCNA::cor
#to switch back
#cor<-stats::cor
net = blockwiseModules(newresults, power = 6,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "TOM", 
                       verbose = 3)

table(net$colors)

TheLabels = matchLabels(net$colors,moduleLabels)
NewLabels = labels2colors(TheLabels)
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
#since net produced two blockwise modules each will get graphed seperately 
plotDendroAndColors(net$dendrograms[[1]], NewLabels[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors

sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[2]], NewLabels[net$blockGenes[[2]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05)


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "rnaseq-03-networkConstruction-auto.RData")


lnames = load(file = "rnaseqoils-02-dataInput.RData" )
lnames = load(file = "rnaseq-03-networkConstruction-auto.RData")



# Define numbers of genes and samples
nGenes = ncol(newresults)
nSamples = nrow(newresults)
dataexpr = data.frame(newresults)

# here we define the adjacency matrix using soft thresholding with beta=6
ADJ1=abs(cor(dataexpr,use="p"))^6
# When you have relatively few genes (<5000) use the following code
k=as.vector(apply(ADJ1,2,sum, na.rm=T))
# When you have a lot of genes use the following code
k=softConnectivity(datE=dataexpr,power=6)
# Plot a histogram of k and a scale free topology plot
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")


datExpr=dataexpr[, rank(-k,ties.method="first" )<=3600]

dissADJ=1-ADJ1
dissTOM=TOMdist(ADJ1)
collectGarbage()
hierADJ=hclust(as.dist(dissADJ), method="average" )


# Plot the resulting clustering tree together with the true color assignment
sizeGrWindow(10,5);
plotDendroAndColors(hierADJ, colors = data.frame(truemodule), dendroLabels = FALSE, hang = 0.03,
                    main = "Gene hierarchical clustering dendrogram and simulated module colors" )



# Recalculate MEs with color labels
MEs0 = moduleEigengenes(dataexpr, moduleColors)$eigengenes
signif(cor(MEs0, use="p"),2)

dissimME=(1-t(cor(MEs0, method="p")))/2
hclustdatME=hclust(as.dist(dissimME), method="average" )
# Plot the eigengene dendrogram
par(mfrow=c(1,1))
plot(hclustdatME, main="Clustering tree based of the module eigengenes")

sizeGrWindow(8,9)
plotMEpairs(MEs0,y=y)



signif(cor(MEs0, ModuleEigengeneNetwork1[,-1]),2)



MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)


sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));

pdf("module-trait_relationships.pdf", height=10, width=15)  

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors=blueWhiteRed(28296),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#above heatmap giving errors onto next part 3b
# Define variable weight containing the weight column of datTrait
treatments = as.data.frame(datTraits1$datTraits);
names(treatments) = "treatment"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(dataexpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(dataexpr, treatments, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(treatments), sep="");
names(GSPvalue) = paste("p.GS.", names(treatments), sep="");
module = "yellow"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for treatment type",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


names(dataexpr)
names(dataexpr)[moduleColors=="yellow"]


dissTOM = 1-TOMsimilarityFromExpr(dataexpr, power = 6);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")


nSelect = 400
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There’s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, 400 selected genes")

# Recalculate module eigengenes
MEs2 = moduleEigengenes(dataexpr, moduleColors)$eigengenes
# Isolate weight from the clinical traits
weight = as.data.frame(datTraits);
names(weight) = "condition"
# Add the weight to existing module eigengenes
METs = orderMEs(cbind(MEs2, weight))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MEs2, "Eigengene dendogram and eigengene adjancey heatmap", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)

# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MEs2, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE, colorLabels=TRUE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(METs, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)


sizeGrWindow(7,7)
# The following shows the correlations between the top genes
plotNetworkHeatmap(dataexpr, plotGenes = gene.names,
                   networkType="unsigned", useTOM=FALSE,
                   power=1, main="signed correlations")


