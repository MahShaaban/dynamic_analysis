library(bio3d)

dcdfile <- "data/HD_Zn.dcd"
pdbfile <- "data/HD_Zn.pdb"

dcd <- read.dcd(dcdfile)
pdb <- read.pdb(pdbfile)

ca.inds <- atom.select(pdb, elety="CA")


xyz <- fit.xyz(fixed=pdb$xyz, mobile=dcd,
               fixed.inds=ca.inds$xyz,
               mobile.inds=ca.inds$xyz)

dim(xyz) == dim(dcd)

rd <- rmsd(xyz[1,ca.inds$xyz], xyz[,ca.inds$xyz])

pdf(file = 'figures/rmsd.pdf')
plot(rd, typ="l", ylab="RMSD", xlab="Frame No.")
points(lowess(rd), typ="l", col="red", lty=2, lwd=2)
dev.off()

pdf(file = 'figures/rmsd_histogram.pdf')
hist(rd, breaks=40, freq=FALSE, main="RMSD Histogram", xlab="RMSD")
lines(density(rd), col="gray", lwd=3)
dev.off()

pdf(file = 'figures/rmsf.pdf')
rf <- rmsf(xyz[,ca.inds$xyz])
plot(rf, ylab="RMSF", xlab="Residue Position", typ="l")
dev.off()

pdf(file = 'figures/pca1.pdf')
pc <- pca.xyz(xyz[,ca.inds$xyz])
plot(pc, col=bwr.colors(nrow(xyz)) )
dev.off()

pdf(file = 'figures/pca2.pdf')
hc <- hclust(dist(pc$z[,1:2]))
grps <- cutree(hc, k=2)
plot(pc, col=grps)
dev.off()

pdf(file = 'figures/residue_position.pdf')
plot.bio3d(pc$au[,1], ylab="PC1 (A)", xlab="Residue Position", typ="l")
points(pc$au[,2], typ="l", col="blue")
dev.off()

p1 <- mktrj.pca(pc, pc=1, b=pc$au[,1], file="output/pc1.pdb")
p2 <- mktrj.pca(pc, pc=2,b=pc$au[,2], file="output/pc2.pdb")

write.ncdf(p1, "output/trj_pc1.nc")

cij <- dccm(xyz[,ca.inds$xyz])

pdf(file = 'figures/residue_cross_correlation.pdf')
plot(cij)
dev.off()

pymol.dccm(cij, pdb, type="launch")
