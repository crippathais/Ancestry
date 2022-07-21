####################
# PCA
# Created by Thais Crippa
# July 19th, 2022
####################

# install packages
pacotes <- c("rgl", "scales", "pca3d", "RColorBrewer", "stats", 
              "epicontacts", "hierfstat", "dplyr", "magick")

if(sum(as.numeric(!pacotes %in% installed.packages())) != 0){
  instalador <- pacotes[!pacotes %in% installed.packages()]
  for(i in 1:length(instalador)) {
    install.packages(instalador, dependencies = T)
    break()}
  sapply(pacotes, require, character = T) 
} else {
  sapply(pacotes, require, character = T) 
}

####################
#Eigenvectors
####################

# read PC eigenvectors from PLINK to R
pca <- read.table("alldatahg381KGP_MAF_LD_IBD.eigenvec",header=T,sep="\t")

# color_all_regions<-recode(levels(pca$FID),"'AFR'='red';'BRA'='black';'AMR'='orange';'EAS'='yellow';'EUR'='blue';'SAS'='purple'")
# c('AYMARAN','MAYAN','NAHUAN','QUECHUAN')='green'

# pca[,'FID']<-if(grepl('^PEL', pca$IID), 'NAT')

####################
# Eigenvalues
####################

# read PC eingenvalues from PLINK to R
pca_eigenvalue <- read.table("alldatahg381KGP_MAF_LD.eigenval",header=F)

# calculate PC eigenvalues proportions
pca_eigenvalue <-data.frame(eingenvalue=pca_eigenvalue$V1,variation=round((pca_eigenvalue$V1/sum(pca_eigenvalue$V1))*100,digits = 1))

####################
# PCA 2D
####################

# create labels for the plot
pcx_eigenvalue <- paste0("PC1 (",pca_eigenvalue$variation[1]," %)")
pcy_eigenvalue <- paste0("PC2 (",pca_eigenvalue$variation[2]," %)")
pcz_eigenvalue <- paste0("PC3 (",pca_eigenvalue$variation[3]," %)")

temp <- as.factor(pca$FID)

myCol <- transp(c("red", # AFR
                  "orange", # AMR
                  "black",# BRA
                  "yellow",# EAS
                  "blue",# EUR
                  "green", # NAT
                  "purple"),# SAS
                .7)[temp]

png(filename = 'pca2D.png',width = 2000, height = 1000,res = 150)
par(mfrow=c(1,3))
plot(pca$PC1,pca$PC2, col=myCol, cex=1, pch=16, main="PC1 x PC2", xlab=pcx_eigenvalue, ylab=pcy_eigenvalue) #PC1 X PC2
plot(pca$PC1,pca$PC3, col=myCol, cex=1, pch=16, main="PC1 x PC3", xlab=pcx_eigenvalue, ylab=pcz_eigenvalu) #PC1 X PC3
plot(pca$PC2,pca$PC3, col=myCol, cex=1, pch=16, main="PC2 x PC3", xlab=pcy_eigenvalue, ylab=pcz_eigenvalue) #PC2 X PC3
dev.off()

# Filter by "eyes"
pca1<-pca[pca$PC3>=0.2 | pca$PC3<=-0.6 | pca$PC1<=-0.2,]
pca1<-pca1[1:2]
write.table(pca1,file="ExcludeSamplesPCA.txt",quote=F, row.names = F)

####################
# PCA 3D
####################
size <- c(10,10, 10, 10, 10, 10, 10) #Number of populations

irisList <- list(pca[pca$Region=="BRA",],
                 pca[pca$Region=="AME",],
                 pca[pca$Region=="AFR",],
                 pca[pca$Region=="EAS",],
                 pca[pca$Region=="EUR",],
                 pca[pca$Region=="SAS",],
                 pca[pca$Region=="NAT",])

color <- c("black", # BRA
        "red",
        "blue",
        "orange",
        "yellow",
        "blue",
        "purple")

with(pca, plot3d(PC1, PC2, PC3,type="s", box=FALSE, size=0))

# Use a separate call to points3d() to plot points of each size
for(i in seq_along(irisList)) {
  with(irisList[[i]], points3d(PC1,PC2, PC3, col=color[[i]], alpha=0.7,size=size[[i]]), shade=0.1)}

snapshotPCA3d(file="PCA3D.png")

# We can indicate the axis and the rotation velocity
play3d( spin3d( axis = c(0, 0, 1), rpm = 20), duration = 10 )

# Save like gif
movie3d(
  movie="PCA3DAnimation", 
  spin3d( axis = c(0, 0, 1), rpm = 7),
  duration = 10, 
  type = "gif", 
  clean = TRUE)



