library(Seurat)
library(ggplot2)
library(Matrix)
library(gplots)
library(scater)
library(dplyr)
#setwd("~/Documents/project_tracking/Goodridge_Helen/PB-6925--05--14--2019/merged_run/ADT_HTO/")
setwd("~/Documents/project_tracking/Engman_David/DO-8312--12--05--2019_Mouse")

#pool_of_blood_and_bone_marrow_HTO_RNA
matrix_dir1 = "./pool_of_blood_and_bone_marrow/HTO_RNA/filtered_feature_bc_matrix/"
barcode.path1 <- paste0(matrix_dir1, "barcodes.tsv.gz")
features.path1 <- paste0(matrix_dir1, "features.tsv.gz")
matrix.path1 <- paste0(matrix_dir1, "matrix.mtx.gz")
mat1 <- readMM(file = matrix.path1)
feature.names1 = read.delim(features.path1, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names1 = read.delim(barcode.path1, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat1) = barcode.names1$V1
rownames(mat1) = feature.names1$V1 
#rownames(mat) = feature.names$V2 
pattern1 <- c("Hashtag_1_WT1_blood","Hashtag_2_WT2_blood","Hashtag_3_SS1_blood","Hashtag_4_SS2_blood","Hashtag_5_WT1_bone_marrow", "Hashtag_6_WT2_bone_marrow", "Hashtag_11_SS1_bone_marrow", "Hashtag_12_SS2_bone_marrow")
mat_tag1 <- mat1[grep(paste(pattern1, collapse="|"), rownames(mat1), value = TRUE), ]
mat_RNA1 <- mat1[1:27998, ]



#pool_of_blood_and_bone_marrow_ADT_RNA
matrix_dir2 = "./pool_of_blood_and_bone_marrow/ADT_RNA/filtered_feature_bc_matrix/"
barcode.path2 <- paste0(matrix_dir2, "barcodes.tsv.gz")
features.path2 <- paste0(matrix_dir2, "features.tsv.gz")
matrix.path2 <- paste0(matrix_dir2, "matrix.mtx.gz")
mat2 <- readMM(file = matrix.path2)
feature.names2 = read.delim(features.path2, 
                            header = FALSE,
                            stringsAsFactors = FALSE)
barcode.names2 = read.delim(barcode.path2, 
                            header = FALSE,
                            stringsAsFactors = FALSE)
colnames(mat2) = barcode.names2$V1
rownames(mat2) = feature.names2$V2
#rownames(mat) = feature.names$V2 

antibody1 <- c("CD184","Ly-6G","CD11b","CD62L")
#mat_RNA2 <- mat[1:27998, ]
mat_AD1 <-  mat2[grep(paste(antibody1, collapse="|"), rownames(mat2), value = TRUE), ]

 
temp_feature <- rbind(feature.names1, feature.names2) #56008
temp_feature1 <- temp_feature[!duplicated(temp_feature[c(1,2)]),]





#pool_of_lung_and_liver_HTO_RNA
matrix_dir3 = "./pool_of_lung_and_liver/HTO_RNA/filtered_feature_bc_matrix/"
barcode.path3 <- paste0(matrix_dir3, "barcodes.tsv.gz")
features.path3 <- paste0(matrix_dir3, "features.tsv.gz")
matrix.path3 <- paste0(matrix_dir3, "matrix.mtx.gz")
mat3 <- readMM(file = matrix.path3)
feature.names3 = read.delim(features.path3, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names3 = read.delim(barcode.path3, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat3) = barcode.names3$V1
rownames(mat3) = feature.names3$V2
#rownames(mat) = feature.names$V2 

pattern2 <- c("Hashtag_1_WT1_Lung","Hashtag_2_WT2_Lung","Hashtag_3_SS1_Lung","Hashtag_4_SS2_Lung","Hashtag_5_WT1_Liver", "Hashtag_6_WT2_Liver", "Hashtag_11_SS1_Liver", "Hashtag_12_SS2_Liver")
mat_tag2 <- mat3[grep(paste(pattern2, collapse="|"), rownames(mat3), value = TRUE), ]
mat_RNA3 <- mat3[1:27998, ]


#pool_of_lung_and_liver_ADT_RNA
matrix_dir4 = "./pool_of_lung_and_liver/ADT_RNA/filtered_feature_bc_matrix/"
barcode.path4 <- paste0(matrix_dir4, "barcodes.tsv.gz")
features.path4 <- paste0(matrix_dir4, "features.tsv.gz")
matrix.path4 <- paste0(matrix_dir4, "matrix.mtx.gz")
mat4 <- readMM(file = matrix.path4)
feature.names4 = read.delim(features.path4, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names4 = read.delim(barcode.path4, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat4) = barcode.names4$V1
rownames(mat4) = feature.names4$V2
#rownames(mat) = feature.names$V2 

antibody2 <- c("CD184","Ly-6G","CD11b","CD62L")
#mat_RNA4 <- mat[1:27998, ]
mat_AD2 <-  mat4[grep(paste(antibody2, collapse="|"), rownames(mat4), value = TRUE), ]


temp_feature <- rbind(feature.names3, feature.names4) #56008
temp_feature2 <- temp_feature[!duplicated(temp_feature[c(1,2)]),] #28010


#RNA, ADT and HTO for pool_of_blood_and_bone_marrow: mat_RNA3, mat_tag2, mat_AD2
#RNA, ADT and HTO for pool_of_lung_and_liver: mat_RNA1, mat_tag1, mat_AD1




# Load in the UMI matrix
#pbmc.umis <- readRDS("pbmc_umi_mtx.rds")
# For generating a hashtag count matrix from FASTQ files, please refer to
# https://github.com/Hoohm/CITE-seq-Count.  Load in the HTO count matrix
#pbmc.htos <- readRDS("pbmc_hto_mtx.rds")
# Select cell barcodes detected by both RNA and HTO In the example datasets we have already
# filtered the cells for you, but perform this step for clarity.
#joint.bcs <- intersect(colnames(pbmc.umis), colnames(pbmc.htos))
# Subset RNA and HTO counts by joint cell barcodes
#pbmc.umis <- pbmc.umis[, joint.bcs]
#pbmc.htos <- as.matrix(pbmc.htos[, joint.bcs])

#RNA, ADT and HTO for pool_of_blood_and_bone_marrow: mat_RNA3, mat_tag2, mat_AD2
#RNA, ADT and HTO for pool_of_lung_and_liver: mat_RNA1, mat_tag1, mat_AD1




# Confirm that the HTO have the correct names
#rownames(pbmc.htos)
rownames(mat_tag1)
rownames(mat_AD1)
# Setup Seurat object
data1.umis <- mat_RNA1

data1 <- CreateSeuratObject(counts = data1.umis)

#library("scater")
#HTO.raw.data <- as.matrix(GetAssayData(data.hashtag@assays$HTO, slot = "counts"))
#example_sce <- SingleCellExperiment(
#  assays = list(counts = HTO.raw.data),
  #colData = sc_example_cell_info)
#example_sce <- calculateQCMetrics(example_sce)
#keep.total_lower <- isOutlier(example_sce$total_counts, nmads=3,type="lower", log=TRUE)
#keep.total_higher <- isOutlier(example_sce$total_counts, nmads=3,type="higher", log=TRUE)

#filtered_lower <- example_sce[,keep.total_lower]
#filtered_higher <- example_sce[,keep.total_higher]

# Normalize RNA data with log normalization
data1 <- NormalizeData(data1)

# Find and scale variable features
data1 <- FindVariableFeatures(data1)
data1 <- ScaleData(data1)

# Add HTO data as a new assay independent from RNA
data1[["HTO"]] <- CreateAssayObject(counts = mat_tag1)
data1[["ADT"]] <- CreateAssayObject(counts = mat_AD1)

#For CITE-seq data, we do not recommend typical
# LogNormalization. Instead, we use a centered log-ratio (CLR) normalization, computed
# independently for each feature.  This is a slightly improved procedure from the original
# publication, and we will release more advanced versions of CITE-seq normalizations soon.
data1 <- NormalizeData(data1, assay = "HTO", normalization.method = "CLR")
data1 <- NormalizeData(data1, assay = "ADT", normalization.method = "CLR")
saveRDS(data1, "pool_of_blood_and_bone_marrow.rds")
data1 <- readRDS("./pool_of_blood_and_bone_marrow.rds")


#RNA, ADT and HTO for pool_of_lung_and_liver: mat_RNA1, mat_tag1, mat_AD1
data2.umis <- mat_RNA3
data2 <- CreateSeuratObject(counts = data2.umis)
data2 <- NormalizeData(data2)
data2 <- FindVariableFeatures(data2)
data2 <- ScaleData(data2)
data2[["HTO"]] <- CreateAssayObject(counts = mat_tag2)
data2[["ADT"]] <- CreateAssayObject(counts = mat_AD2)
data2 <- NormalizeData(data2, assay = "HTO", normalization.method = "CLR")
data2 <- NormalizeData(data2, assay = "ADT", normalization.method = "CLR")
saveRDS(data2, "pool_of_lung_and_liver.rds")
data2 <- readRDS("./pool_of_lung_and_liver.rds")


test <- data1
name <- "pool_of_blood_and_bone_marrow"

test <- data2
name <- "pool_of_lung_and_liver"

#rownames(test@assays$RNA)[grep('^Mt-',rownames(test@assays$RNA))]
#rownames(test@assays$RNA)[grep('^MT-',rownames(test@assays$RNA))]
rownames(test@assays$RNA)[grep('^mt-',rownames(test@assays$RNA))]


test[["percent.mt"]] <- PercentageFeatureSet(test, pattern = "^mt-")
head(test@meta.data, 5)
plot1 <- FeatureScatter(test, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(test, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
test_1 <- subset(test, subset = nFeature_RNA >= 300 & percent.mt <= 15)

plot1_1 <- FeatureScatter(test_1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2_1 <- FeatureScatter(test_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#CombinePlots(plots = list(plot1_1, plot2_1))
pdf(paste(name,"_QC_scRNA.pdf",sep=""),16,12)
test_4 <- subset(test, subset = nFeature_RNA < 300 & percent.mt > 15)
#write.csv(test@meta.data,"pool_of_blood_and_bone_marro_meta_data.csv")
e<-dim(test_4@assays$RNA)[2]
#e=0
test_2 <- subset(test, subset = percent.mt <= 15)
a<-dim(test@assays$RNA)[2]-dim(test_2@assays$RNA)[2]-e
test_3 <- subset(test, subset = nFeature_RNA >= 300)
b<-dim(test@assays$RNA)[2]-dim(test_3@assays$RNA)[2]-e
c<-dim(test_1@assays$RNA)[2]
d<-dim(test@assays$RNA)[2]
text1<-paste("Sample Name:",name,sep=" ")
text2<-paste(a,"cells failed mito% <= 15%",sep=" ")
text3<-paste(b,"cells failed total # expressed genes >= 300.",sep=" ")
text5<-paste(e,"cells failed total # expressed genes >= 300 and mito% <= 15%.",sep=" ")
text4<-paste("There are",c,"out of",d,"cells remained after filtering.",sep=" ")
text<-paste(text1,text2,text3,text5,text4,sep="\n");

textplot(text,halign="center",valign="center",cex=2)
textplot("Before Filtering",halign="center",valign="center",cex=5)
VlnPlot(test, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

CombinePlots(plots = list(plot1, plot2))
textplot("After Filtering",halign="center",valign="center",cex=5)
VlnPlot(test_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
CombinePlots(plots = list(plot1_1, plot2_1))
dev.off()


test_1<-NormalizeData(object = test_1)
expr_raw<-GetAssayData(object = test_1, slot = "counts")
expr_norm<-GetAssayData(object = test_1, slot = "data")

HTO_raw<-GetAssayData(object = test_1, assay = "HTO", slot = "counts")
HTO_norm<-GetAssayData(object = test_1, assay = "HTO", slot = "data")

ADT_raw<-GetAssayData(object = test_1, assay = "ADT", slot = "counts")
ADT_norm<-GetAssayData(object = test_1, assay = "ADT", slot = "data")

expr_ADT_HTO_raw <- rbind(expr_raw, HTO_raw, ADT_raw)

expr_ADT_HTO_norm <- rbind(expr_norm, HTO_norm, ADT_norm)

#temp_feature1
expr_ADT_HTO_raw_final <- cbind(temp_feature1[, c(1,2)], expr_ADT_HTO_raw) 
colnames(expr_ADT_HTO_raw_final)[1] <- "Ensembl_ID"
colnames(expr_ADT_HTO_raw_final)[2] <- "Gene"

expr_ADT_HTO_norm_final <- cbind(temp_feature1[, c(1,2)], expr_ADT_HTO_norm) 
colnames(expr_ADT_HTO_norm_final)[1] <- "Ensembl_ID"
colnames(expr_ADT_HTO_norm_final)[2] <- "Gene"


#temp_feature2
expr_ADT_HTO_raw_final <- cbind(temp_feature2[, c(1,2)], expr_ADT_HTO_raw) 
colnames(expr_ADT_HTO_raw_final)[1] <- "Ensembl_ID"
colnames(expr_ADT_HTO_raw_final)[2] <- "Gene"

expr_ADT_HTO_norm_final <- cbind(temp_feature2[, c(1,2)], expr_ADT_HTO_norm) 
colnames(expr_ADT_HTO_norm_final)[1] <- "Ensembl_ID"
colnames(expr_ADT_HTO_norm_final)[2] <- "Gene"




raw_name<-paste0("./",name,"_Expr_HTO_ADT_raw_QC.csv")
norm_name<-paste0("./",name,"_Expr_HTO_ADT_norm_QC.csv")

write.csv(expr_ADT_HTO_raw_final,raw_name,quote=F,row.names = FALSE)
write.csv(expr_ADT_HTO_norm_final,norm_name,quote=F,row.names = FALSE)


temp<-test_1@assays$RNA
barcode<-colnames(temp)
barcode<-data.frame(Barcode=barcode)
barcode[,1]<-paste(barcode[,1],"-1",sep="")
write.csv(barcode,paste(name,"barcode_filtered.csv",sep="_"),row.names = F,quote = F)
