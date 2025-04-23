#==============================================
#        九十 利用DESeq2分析RNAseq数据        #
#==============================================
#library(BiocManager)
#BiocManager::install("DESeq2")

library("DESeq2")
#==============================================
#            读取readcount矩阵文件              #
#==============================================
countdata <- read.csv("CountMatrix.csv",header = T,row.names = 1)
head(countdata, 10)


#生成分组信息
coldata <- read.csv("sample.csv",header = T)
#关键步骤，生成Deseq对象
dds <- DESeqDataSetFromMatrix(countData = countdata,colData = coldata,design = ~ cell + dex)


#==============================================
#              查看DESeqDataSet对象           #
#==============================================

dim(dds)
assay(dds)
assayNames(dds)
colSums(assay(dds))
rowRanges(dds)
colData(dds)


#过滤没有reads比对上的基因，所有reads数为零
nrow(dds)
dds <- dds[rowSums(counts(dds)) > 1,]
nrow(dds)
colSums(assay(dds))
barplot(colSums(assay(dds)),las=3,col=rep(rainbow(4),each=2))
#==============================================
#               多维数据探索                  #
#==============================================
# 将数据通过rolg方法与vst方法转换，这样可以用于后面计算距离矩阵
## ----rlog方法----------------------------------------------------------------
rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)

## ----vst方法-----------------------------------------------------------------
vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)


#利用转换后的结果计算样品之间距离关系
#方法1--欧氏距离--------------------------
sampleDists <- dist(t(assay(rld)))
sampleDists

library("pheatmap")
library("RColorBrewer")

## ----distheatmap, fig.width = 6.1, fig.height = 4.5----------------------
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$dex, rld$cell, sep = " - " )
colnames(sampleDistMatrix) <- NULL
pheatmap(sampleDistMatrix,clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists)

## 方法2，使用Poisson Distance方法------------------------------------------
library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( rld$dex, rld$cell, sep=" - " )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd)

## 方法3，PCA--------------------------------------------------------------
plotPCA(rld, intgroup = c("dex", "cell"))

## ------------------------------------------------------------------------
pcaData <- plotPCA(rld, intgroup = c( "dex", "cell"), returnData = TRUE)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))
library(ggplot2)
ggplot(pcaData, aes(x = PC1, y = PC2, color = dex, shape = cell)) +
  geom_point(size =3) +  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +  coord_fixed()

## 方法4，MDS plot---------------------------------------------------------
library(dplyr)
mds <- as.data.frame(colData(rld))  %>%   cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = dex, shape = cell)) +  
  geom_point(size = 3) + coord_fixed()
mdsPois <- as.data.frame(colData(dds)) %>%
  cbind(cmdscale(samplePoisDistMatrix))
ggplot(mdsPois, aes(x = `1`, y = `2`, color = dex, shape = cell)) +
  geom_point(size = 3) + coord_fixed()

#==============================================
#                 差异表达计算                 #
#==============================================

dep <- DESeq(dds)
res <- results(dep)
res
write.csv(x = res,file = "des.csv")


#筛选出差异表达基因，log2foldchange >=1或者<=-1,并且q值小于0.05的基因
res <- na.omit(res)
res0.05 <- res[(abs(res$log2FoldChange)>=1 & res$padj <=0.05), ]
nrow(res0.05)

#筛选出差异表达基因，log2foldchange >=2,并且q值小于0.01的基因
res0.01 <- res[(abs(res$log2FoldChange)>=1 & res$padj <=0.01), ]
nrow(res0.01)

write.csv(res0.05,file = "res0.05.csv")
#找出差异表达最显著的基因，按log2差异倍数排序
head(res0.05[order(res0.05$log2FoldChange,decreasing = TRUE), ])
head(res0.05[order(res0.05$log2FoldChange,decreasing = FALSE), ])

#==============================================
#                  结果可视化                 #
#==============================================
#MA图
plotMA(res, ylim = c(-10,10))

#火山图
m <- res
head(m)
m <- na.omit(m)
plot(m$log2FoldChange,m$padj)
plot(m$log2FoldChange,-1*log10(m$padj))
plot(m$log2FoldChange,-1*log10(m$padj),xlim = c(-10,10),ylim=c(0,100))
m <- transform(m,padj=-1*log10(m$padj))
down <- m[m$log2FoldChange<=-1,] 
up <- m[m$log2FoldChange>=1,]
no <- m[m$log2FoldChange>-1 & m$log2FoldChange <1,] 

plot(no$log2FoldChange,no$padj,xlim = c(-10,10),ylim=c(0,100),col="blue",pch=16,
     cex=0.8,main = "Gene Expression",
     xlab = "log2FoldChange",ylab="-log10(Qvalue)")

points(up$log2FoldChange,up$padj,col="red",pch=16,cex=0.8)
points(down$log2FoldChange,down$padj,col="green",pch=16,cex=0.8)
abline(v = c(-1,1),h = 2)


#==============================================
#         九十一 利用edgeR分析RNAseq数据      #
#==============================================
rm(list=ls())
#安装R包
install.packages("BiocM")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("RNAseq123")
library(limma)
library(Glimma)
library(edgeR)

#设置工作目录,xxx为表达矩阵数据目录
#setwd("xxx")

#==============================================
#    方法一：读取featureCounts表达矩阵        #
#==============================================
countdata <- read.table("gene.counts.matrix.S_vs_C.edgeR.count_matrix",header = T,row.names = 1)
head(countdata)
#过滤表达矩阵
nrow(countdata)
countdata <- countdata[rowSums(countdata)>10,]
nrow(countdata)

#读取样品信息
coldata <- read.table("../sample.info.txt",header = T)
head(coldata)
View(coldata)

#生成DGElist对象
dge <- DGEList(counts = countdata, genes = row.names(countdata),group = as.factor(coldata$group))
dge

#添加样品信息
dge$samples$cell <- as.factor(coldata$sample)

#备份一份
dge2 <- dge

#==============================================
#     方法二：读取htseq多个文件表达信息       #
#==============================================
setwd("../htseq/")
files <- dir(pattern = "tsv")
files

#查看文件内容
read.delim(files[1], nrow=5,header = F)
#生成DGElist对象
x <- readDGE(files = files,columns = c(1,2))
x
#查看属性
class(x)
x$samples
x$counts
dim(x)
samplenames <- substring(colnames(x),1,nchar(colnames(x)))
samplenames

colnames(x) <- samplenames

x$samples$group <- as.factor(coldata$dex)
x$samples$cell <- as.factor(coldata$cell)
x$samples

#==============================================
#         DGElist对象探索          #
#==============================================
barplot(colSums(dge$counts),las=2,col = rep(rainbow(4),each=2))
opar <- par(no.readonly = TRUE)
par(mar=c(7.1,4.1,4.1,2.1))
barplot(colSums(dge$counts),las=2,col = rep(c("red","green"),each=3))
par(opar)
#==============================================
#         对数据计算标准化因子            #
#==============================================
# Normalization method: "TMM","TMMwsp","RLE","upperquartile","none"
dge.tmm <- calcNormFactors(dge,method = "TMM")
dge.tmm$samples
plotMDS(dge.tmm)

dge.tmmwsp <- calcNormFactors(dge,method = "TMMwsp")
dge.tmmwsp$samples
plotMDS(dge.tmmwsp)

dge.upperquartile <- calcNormFactors(dge,method = "upperquartile")
dge.upperquartile$samples
plotMDS(dge.upperquartile)

dge.rle <- calcNormFactors(dge,method = "RLE") # DESeq2, cuffdiff
dge.rle$samples
plotMDS(dge.rle)

dge.none <- calcNormFactors(dge,method = "none")
dge.none$samples
plotMDS(dge.none)

dge <- calcNormFactors(dge,method = "TMM")
dge
dge$samples

#检查样本中的异常值，使用plotMDS()作图。
plotMDS(dge)
plotMDS(dge,labels = dge$samples$cell,col=rep(c("red","green"),each=3))
#==============================================
#               差异表达计算               #
#==============================================
#1 生成实验矩阵
model.matrix( ~ group,data = dge$samples)
design.matrix <- model.matrix( ~ group ,data = dge$samples)
rownames(design.matrix) <- rownames(dge$samples)
design.matrix

#2 估计离散系数（dispersion）
dge <- estimateDisp(dge,design = design.matrix)
dge$common.dispersion
dge$tagwise.dispersion

#可以用BCV plot查看离散系数；
plotBCV(dge, cex = 0.8)

# plot var and mean
plotMeanVar(dge, show.raw=TRUE, show.tagwise=TRUE, show.binned=TRUE)

#3 exactTest检验
dge_res <- exactTest(dge)
dge_results <- as.data.frame(topTags(dge_res,n=nrow(dge$counts),sort.by = "logFC"))
View(dge_results)

#3  最大似然法检测
fit <- glmFit(y = dge, design = design.matrix)
lrt <- glmLRT(fit, coef=2)
lrt_results <- as.data.frame(topTags(lrt,n=nrow(dge$counts),sort.by = "logFC"))
View(lrt_results)

#4 筛选差异表达基因
write.csv(x = dge_results,file = "deg_results.csv")


#筛选出差异表达基因，logFC >=1或者<=-1,并且q值小于0.05的基因
res <- na.omit(dge_results)
res0.05 <- res[(abs(res$logFC)>=1 & res$FDR <=0.05), ]
nrow(res0.05)

#筛选出差异表达基因，logFC >=2,并且q值小于0.01的基因
res0.01 <- res[(abs(res$logFC)>=1 & res$FDR <=0.01), ]
nrow(res0.01)

write.csv(res0.05,file = "res0.05.csv")
#找出差异表达最显著的基因，按log2差异倍数排序
head(res0.05[order(res0.05$logFC,decreasing = TRUE), ])
head(res0.05[order(res0.05$FDR,decreasing = TRUE), ])

#5 结果可视化 
#绘制MA plot
select.sign.gene = decideTestsDGE(dge_res, p.value = 0.01) 
select.sign.gene_id = rownames(dge_res)[as.logical(select.sign.gene)]
plotSmear(dge_res, de.tags = select.sign.gene_id, cex = 0.5,ylim=c(-4,4)) 
abline(h = c(-2, 2), col = "blue")

#绘制差异基因散点图（类似火山图）；
plotMD(dge_res,ylim=c(-10,10))
abline(h=c(-1, 1), col="green")





#==============================================
#               九十二 基因功能注释           #
#==============================================
library(dplyr)
library(tidyr)
library(DOSE)
library(GO.db)
library(org.Hs.eg.db)
library(GSEABase)
library(clusterProfiler)

dta <- read.csv("res0.05.csv",header = T,row.names = 1,stringsAsFactors = F)

## ------------------------------------------------------------------------
x <- rownames(dta)
length(x)
pm <- AnnotationDbi::loadDb("../org.Petromyzon_marinus.eg.sqlite")
keytypes(pm.db)
eg <-  bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb=pm.db)
head(eg)
ids <- bitr(x, fromType="SYMBOL", toType="SYMBOL", OrgDb=pm.db)
head(ids)
go <-  bitr(x,fromType = "SYMBOL",toType = c("SYMBOL","GO","ONTOLOGY"),OrgDb = pm.db)
head(go)
#====================================
# GO功能富集分析，输入Gene ID 为GI号 #
#====================================
gene <- eg$ENTREZID
gene.df <- bitr(gene, fromType = "ENTREZID",
                toType = c("SYMBOL"),
                OrgDb = pm.db)
head(gene.df)
ggo <- groupGO(gene     = gene,
               OrgDb    = pm.db,
               ont      = "MF",
               level = 3,
               readable = TRUE)
head(ggo)
View(as.data.frame(ggo))

ego <- enrichGO(gene          = gene,
                OrgDb         = pm.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                readable      = TRUE)
head(ego)
View(as.data.frame(ego))
write.csv(as.data.frame(ego),file = "GO_BP_gene.csv")

## GO功能富集可视化 
barplot(ggo, drop=TRUE, showCategory=12)
barplot(ego, drop=TRUE, showCategory=12)
## ----fig.height=5, fig.width=8-------------------------------------------
barplot(ego, showCategory=15,
        color = "p.adjust")+
  scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))

dotplot(ego,
        x = "GeneRatio",
        font.size=10,
        color = "p.adjust")+
        scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))
goplot(ego)

#====================================
#           KEGG富集分析            #
#====================================
search_kegg_organism("pmrn", by='kegg_code')
hsa <- search_kegg_organism('hsa', by='kegg_code')
search_kegg_organism('sea', by='common_name')
ecoli <- search_kegg_organism('Escherichia coli', by='scientific_name')
dim(hsa)
head(hsa)


kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
head(kk)


## ----eval = FALSE--------------------------------------------------------
mkk <- enrichMKEGG(gene = gene,organism = 'hsa')

browseKEGG(kk, 'hsa04110')
barplot(kk)
dotplot(kk)

#====================================
#      #单独绘制GO条目图             #
#====================================
library(ggplot2)
go <- read.csv("go.csv",header = T)
head(go)
go_sort <- go[order(go$Ontology,-go$Percentage),]
m <- go_sort[go_sort$Ontology=="Molecular function",][1:10,]
c <- go_sort[go_sort$Ontology=="Cellular component",][1:10,]
b <- go_sort[go_sort$Ontology=="Biological process",][1:10,]
slimgo <- rbind(b,c,m)

#首先需要将Trem转换为因子
slimgo$Term=factor(slimgo$Term,levels=slimgo$Term)
head(slimgo)
colnames(slimgo)
p <- ggplot(data = slimgo, mapping = aes(x=Term,y=Percentage,fill=Ontology))
p
p+geom_bar(stat="identity")
p+geom_bar(stat="identity")+coord_flip()
p+geom_bar(stat="identity")+coord_flip()+scale_x_discrete(limits=rev(levels(slimgo$Term)))
p+geom_bar(stat="identity")+coord_flip()+scale_x_discrete(limits=rev(levels(slimgo$Term)))+guides(fill=FALSE)
p+geom_bar(stat="identity")+coord_flip()+scale_x_discrete(limits=rev(levels(slimgo$Term)))+guides(fill=FALSE)+theme_bw()

#====================================
#     单独绘制kegg条目图            #
#====================================
library(ggplot2)

pathway <-  read.csv("kegg.csv",header=T)
head(pathway)
colnames(pathway)

p <-  ggplot(data=pathway,mapping = aes(x=richFactor,y=Pathway))
p
p + geom_point()
p + geom_point(aes(size=AvsB))
p + geom_point(aes(size=AvsB,color=Qvalue))
p + geom_point(aes(size=AvsB,color=Qvalue)) + scale_colour_gradient(low="green",high="red")
p + geom_point(aes(size=AvsB,color=Qvalue)) + scale_colour_gradient(low="green",high="red")+labs(title="Top20 of pathway enrichment",x="Rich factor",y="Pathway Name",color="-log10(Qvalue)",size="Gene Numbers")
p + geom_point(aes(size=AvsB,color=Qvalue)) + scale_colour_gradient(low="green",high="red")+labs(title="Top20 of pathway enrichment",x="Rich factor",y="Pathway Name",color="-log10(Qvalue)",size="Gene Numbers")+theme_bw()

#==============================================
#               GSEA分析                     #
#==============================================
#====================================
#             方法1  GSEAbase       #
#====================================
#BiocManager::install("GSEABase")
library(GSEABase)
library(clusterProfiler)
library(org.Hs.eg.db)

#生成基因排序列表
x <- read.csv("DESeq2.156098.dir/deg_result.csv",row.names = 1)
head(x)

geneList <- x$log2FoldChange
names(geneList) <- toupper(rownames(x))
geneList <- sort(geneList, decreasing = T)
head(geneList)
#names(geneList) <- mapIds(org.Hs.eg.db, keys=names(geneList), column="SYMBOL", keytype="ENSEMBL",multiVals = "first")
head(geneList)
#去除NA
geneList <- geneList[!is.na(names(geneList))]

#读取gmt文件
gmtfile <- "msigdb.v7.5.1.symbols.gmt"
geneset <- read.gmt(gmtfile)

#GSEA分析
egmt <- GSEA(geneList, TERM2GENE = geneset, minGSSize = 1, pvalueCutoff = 0.01, verbose = F,eps = 0)

#查看结果
gsea.out.df <- egmt@result
View(head(gsea.out.df))
gsea.out.df$ID

#绘制GSEA图
library(enrichplot)
options(repr.plot.width=6,repr.plot.height=4)
gseaplot2(egmt, geneSetID = "GOCC_INNER_MITOCHONDRIAL_MEMBRANE_PROTEIN_COMPLEX", pvalue_table = T)

#====================================
#           方法2  fgsea            #
#====================================
BiocManager::install("fgsea")
BiocManager::install("msigdbr")
BiocManager::install("fgsea")
BiocManager::install("fgsea")
BiocManager::install("fgsea")

library(msigdbr)
library(fgsea)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)

## 预定义基因集
mdb_c2 <- msigdbr(species = "Homo sapiens")
mdb_c2
mdb_kegg = mdb_c2 [grep("^KEGG",mdb_c2$gs_name),]
mdb_kegg

fgsea_sets <- mdb_kegg %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_sets

#读取基因差异表达文件
x <- read.csv("res.csv",row.names = 1)
head(x)


#基因按logFC排序
rownames(x)
x$genes <- mapIds(org.Hs.eg.db, keys=rownames(x), column="SYMBOL", keytype="ENSEMBL",multiVals = "first")
x <- na.omit(x)
gene_log2 <- x %>% arrange(desc(log2FoldChange)) %>% dplyr::select(genes,log2FoldChange)

generanks<- deframe(gene_log2)
generanks

#fgsea分析
fgseaRes <- fgsea(pathways = fgsea_sets, 
                  stats = geneList,
                  minSize=15,
                  maxSize=500,
                  nperm=10000)
fgseaRes <- fgseaRes[fgseaRes$padj <= 0.05]

#fgsea绘图
plotEnrichment(fgsea_sets[["GOBP_AMIDE_BIOSYNTHETIC_PROCESS"]],geneList) + 
  labs(title="GOBP_AMIDE_BIOSYNTHETIC_PROCESS")

