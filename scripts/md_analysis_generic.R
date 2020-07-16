library(bio3d)
library(ggfortify)
library(igraph)

mod.select <- function(x, thres=0.1) {
  remodel <- community.tree(x, rescale = TRUE)
  n.max = length(unique(x$communities$membership))
  ind.max = which(remodel$num.of.comms == n.max)
  v = remodel$modularity[length(remodel$modularity):ind.max]
  v = rev(diff(v))
  fa = which(v>=thres)[1] - 1
  ncomm = ifelse(is.na(fa), min(remodel$num.of.comms), n.max - fa)
  
  ind <- which(remodel$num.of.comms == ncomm)
  network.amendment(x, remodel$tree[ind, ])
}

dcdfile <- "~/ORF2/md_apo/trj_orf2_nomn_3/orf2_CA.dcd" # replace with your dcd file
pdbfile <- "~/ORF2/md_apo/trj_orf2_nomn_3/orf2_CA.pdb" # replace with your pdb file
pdbfile2 <- "~/ORF2/md_apo/trj_orf2_nomn_3/orf2_prot.pdb" # replace with your pdb file (all atoms)

dcd <- read.dcd(dcdfile)
pdb <- read.pdb(pdbfile)
pdb2 <- read.pdb(pdbfile2)

ca.inds <- atom.select(pdb, elety="CA") # makes sure only carbon alphas are being chosen

xyz <- fit.xyz(fixed=pdb$xyz, mobile=dcd, 
               fixed.inds=ca.inds$xyz,
               mobile.inds=ca.inds$xyz)

cij<-dccm(xyz[,ca.inds$xyz])
plot(cij, main="ORF2 En without Mn2+") # plots x-correlation matrix

rd <- rmsd(xyz[1,ca.inds$xyz], xyz[,ca.inds$xyz]) # calculates RMSD

plot(rd, typ="l", ylab="RMSD", xlab="Frame No.")
points(lowess(rd), typ="l", col="red", lty=2, lwd=2) 

hist(rd, breaks=40, freq=FALSE, xlab="RMSD")
lines(density(rd), col="gray", lwd=3) # plots RMSD as histogram

rf <- rmsf(xyz[,ca.inds$xyz]) # calculates RMSF
plot(rf, ylab="RMSF", xlab="Residue Position", typ="l")

pc <- pca.xyz(xyz[,ca.inds$xyz]) # plots PCA
plot(pc, col=bwr.colors(nrow(xyz)) )

cij <- dccm(xyz)
net <- cna(cij) # calculates residue network
plot(net, pdb2)

# whole code below is just to make the network easier to see #

layout3D <- layout.cna(net, pdb2, k=3) 
layout2D <- layout.cna(net, pdb2, k=2)


grp.col <- bwr.colors(20)[net$communities$membership]

grp.col <- list()
for (i in 1:max(net$communities$membership)) {
  grp.tmp <- which(net$communities$membership == i)
  grp.col[[i]] <- grp.tmp
}


colbar.full <- vmd_colors()[net$communities$membership]
colbar.comms <- vmd_colors(max(net$communities$membership), alpha = 0.5)

plot(net, pdb, full = TRUE, mark.groups = grp.col, mark.col = colbar.comms)


