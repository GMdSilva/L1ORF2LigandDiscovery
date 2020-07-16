library(ggfortify)
library(ggplot)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(bio3d)

#

mynamestheme <- theme(plot.title = element_text(family = "Helvetica", face = "bold", size = (15)), 
                      axis.title = element_text(family = "Helvetica", size = (12)),
                      axis.text = element_text(family = "Helvetica", size = (12)))

pca <- read.table("~/ORF2/md_apo/results/bs_pca/orf2_bs_pca1.csv",skip=1, sep="," )

pca$row_num <- seq.int(nrow(pca))
# 1
(d <- ggplot(pca, aes(V1, V2)) +
  geom_point(aes(colour=row_num)))
d <- d + scale_color_distiller(palette="RdYlBu")
d <- d + theme(legend.position = "none")
d <- d + labs(y="PC1",x="PC2", title="ORF2 En with Mn2+ - Binding Site 1")
d <- d+ xlim(-4, 4)
d <- d+ ylim(-4, 4)
d <- d+mynamestheme

pca2 <- read.table("~/ORF2/md_apo/results/bs_pca/orf2_bs_pca2.csv",skip=1, sep=",")

pca2$row_num <- seq.int(nrow(pca2))
# 1
(d2 <- ggplot(pca2, aes(V1, V2)) +
    geom_point(aes(colour=row_num)))
d2 <- d2 + scale_color_distiller(palette="RdYlBu")
d2 <- d2 + theme(legend.position = "none")
d2 <- d2 + labs(y="PC1",x="PC2", title="ORF2 En with Mn2+ - Binding Site 2")
d2 <- d2 + xlim(-4, 4)
d2 <- d2 + ylim(-4, 4)
d2 <- d2+mynamestheme

pca3 <- read.table("~/ORF2/md_apo/results/bs_pca/orf2_bs_pca3.csv",skip=1, sep=",", )

pca3$row_num <- seq.int(nrow(pca3))
# 1
(d3 <- ggplot(pca3, aes(V1, V2)) +
    geom_point(aes(colour=row_num)))
d3 <- d3 + scale_color_distiller(palette="RdYlBu")
d3 <- d3 + theme(legend.position = "none")
d3 <- d3 + labs(y="PC1",x="PC2", title="ORF2 En with Mn2+ - Binding Site 3")
d3 <- d3 + xlim(-4, 4)
d3 <- d3 + ylim(-4, 4)
d3 <- d3+mynamestheme

pca4 <- read.table("~/ORF2/md_apo/results/bs_pca/orf2_nomn_bs_pca1.csv",skip=1, sep=",", )

pca4$row_num <- seq.int(nrow(pca4))
# 1
(d4 <- ggplot(pca4, aes(V1, V2)) +
    geom_point(aes(colour=row_num)))
d4 <- d4 + scale_color_distiller(palette="RdYlBu")
d4 <- d4 + theme(legend.position = "none")
d4 <- d4 + labs(y="PC1",x="PC2", title="ORF2 En without Mn2+ - Binding Site")
d4 <- d4 + xlim(-4, 4)
d4 <- d4 + ylim(-4, 4)
d4 <- d4+mynamestheme

pca5 <- read.table("~/ORF2/md_apo/results/bs_pca/orf2_nomn_bs_pca2.csv",skip=1, sep=",", )

pca5$row_num <- seq.int(nrow(pca5))
# 1
(d5 <- ggplot(pca5, aes(V1, V2)) +
    geom_point(aes(colour=row_num)))
d5 <- d5 + scale_color_distiller(palette="RdYlBu")
d5 <- d5 + theme(legend.position = "none")
d5 <- d5 + labs(y="PC1",x="PC2", title="ORF2 En without Mn2+ - Binding Site 2")
d5 <- d5 + xlim(-4, 4)
d5 <- d5 + ylim(-4, 4)
d5 <- d5+mynamestheme

pca6 <- read.table("~/ORF2/md_apo/results/bs_pca/orf2_nomn_bs_pca3.csv",skip=1, sep=",", )

pca6$row_num <- seq.int(nrow(pca5))
# 1
(d6 <- ggplot(pca6, aes(V1, V2)) +
    geom_point(aes(colour=row_num)))
d6 <- d6 + scale_color_distiller(palette="RdYlBu")
d6 <- d6 + theme(legend.position = "none")
d6 <- d6 + labs(y="PC1",x="PC2", title="ORF2 En without Mn2+ - Binding Site 3")
d6 <- d6 + xlim(-4, 4)
d6 <- d6 + ylim(-4, 4)
d6 <- d6+mynamestheme

grid.arrange(d, d2,d3,d4,d5,d6, nrow = 2)

##

