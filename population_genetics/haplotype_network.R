# haplotype_network.R, written by Jessica Goodheart
# Based in part on https://arundurvasula.wordpress.com/2016/02/24/haplotype-networks-in-r/
# Last updated 17 July 2019

# Install packages
install.packages("pegas")
install.packages("ape")
install.packages("RColorBrewer")
install.packages("mmod")
install.packages("adegenet")

# Load packages
library(ape, pegas, RColorBrewer, mmod, adegenet)

# Input alignment data
#input <- "VtsujiiCO1.TXT"
input <- "Vtsujii16S.fa"

# Read data into ape - 16S
d <- ape::read.dna(input, format='fasta')
f <- ape::dist.dna(d)

# Populations vector -16S 
rownames(d) <- substr(rownames(d), 0, 2)

# Haplotype network - 16S - Range from 1 to 46 specimens per haplotype
h <- pegas::haplotype(d)
h <- sort(h, what = "label")
(net <- pegas::haploNet(h, getProb=TRUE))
ind.hap<-with(
  stack(setNames(attr(h, "index"), rownames(h))),
  table(hap=ind, pop=rownames(d)[values])
)

# Populations labels
pops.label <- c("CI", "SD", "SP")

# Plot 16S
cols<-c("purple","blue","green")
pdf("haplotypes_16S.pdf",width=6,height=4,paper='special')
plot(net, size=attr(net, "freq"), pie=ind.hap, labels=FALSE, bg=cols, threshold=0, show.mutation=2)
legend(-45, 30, pops.label, col=cols, pch=19)
dev.off()

# DAPC Analysis - first find strata
strata <- data.frame("pops"=substr(labels(d),0,2))
g <- as.genind.DNAbin(d, pops=substr(labels(d),0,2))
g@strata <- strata

# Find clusters - retained 4 PCs and 3 clusters
grp.16S <- find.clusters(g, max.n.clust=3)
grps.16S <- table(pop(g),grp.16S$grp)
table.value(table(pop(g), grp.16S$grp), col.lab=paste("inferred", 1:3),row.lab=paste("original", pops.label))

# Plot clustering results (PDF) - 3 clusters favored
pdf("clusters_16S.pdf",width=6,height=4,paper='special')
table.value(table(pop(g), grp.16S$grp), col.lab=paste("inferred", 1:3),row.lab=paste("original", pops.label))
dev.off()

#################################################
#### COI
#################################################
# Input alignment data
input.COI <- "VtsujiiCO1.fa"

# Read data into ape - COI
d.COI <- ape::read.dna(input.COI, format='fasta')
f.COI <- ape::dist.dna(d.COI)

# Populations vector -COI 
rownames(d.COI) <- substr(rownames(d.COI), 0, 2)

# Haplotype network - COI
h.COI <- pegas::haplotype(d.COI)
h.COI <- sort(h.COI, what = "label")
(net.COI <- pegas::haploNet(h.COI))
ind.hap.COI<-with(
  stack(setNames(attr(h.COI, "index"), rownames(h.COI))),
  table(hap=ind, pop=rownames(d.COI)[values])
)

# Plot COI
cols<-c("purple","blue","green")
pdf("haplotypes_COI.pdf",width=6,height=4,paper='special')
plot(net.COI, size=attr(net.COI, "freq"), pie=ind.hap.COI, labels=FALSE, bg=cols, threshold=0, show.mutation=2)
legend(-25, 7, pops.label, col=cols, pch=19)
dev.off()

#DAPC Analysis - find strata with AMOVA
strata.COI <- data.frame("pops"=substr(labels(d.COI),0,2))
g.COI <- as.genind.DNAbin(d.COI, pops=substr(labels(d.COI),0,2))
g.COI@strata <- strata.COI

# Find clusters - retained 14 PCs and 3 clusters
grp.COI <- find.clusters(g.COI, max.n.clust=3)
grps.COI <- table(pop(g.COI),grp.COI$grp)
table.value(table(pop(g.COI), grp.COI$grp), col.lab=paste("inferred", 1:3),row.lab=paste("original", pops.label))

# Plot clustering results - 3 clusters favored
pdf("clusters_COI.pdf",width=6,height=4,paper='special')
table.value(table(pop(g.COI), grp.COI$grp), col.lab=paste("inferred", 1:3),row.lab=paste("original", pops.label))
dev.off()
