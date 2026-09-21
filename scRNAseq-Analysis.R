install.packages("Seurat")
library(Seurat)

library(Seurat)
# 1. حساب نسبة الجينات الميتوكندرية لكل خلية
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# 2. رسم بياني لمعاينة الجودة (QC Plots)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# 1. فلترة البيانات طبقاً لمعايير الجودة
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# 2. Normalization (تطبيع البيانات)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# 3. تحديد الجينات الأكثر تغيراً (Highly Variable Features)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# التأكد من نجاح الفلترة
pbmc
# 1. Scale Data (معايرة البيانات)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# 2. Run PCA (تحليل المكونات الرئيسية)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# 3. Clustering (تجميع الخلايا بناءً على الشبه بينها)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# 4. Run UMAP (إنشاء الخريطة ثنائية الأبعاد)
pbmc <- RunUMAP(pbmc, dims = 1:10)

# 5. رسم خريطة الـ UMAP النهائية!
DimPlot(pbmc, reduction = "umap", label = TRUE)
# 1. إعطاء كل Cluster اسمها البيولوجي الحقيقي
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", 
                     "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

# 2. رسم خريطة الـ UMAP النهائية بالأسماء
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# 3. حفظ الرسمة كصورة عالية الجودة للمشروع
ggsave("scRNAseq_PBMC_UMAP.png", width = 8, height = 6, dpi = 300)
library(ggplot2)
ggsave("scRNAseq_PBMC_UMAP.png", width = 8, height = 6, dpi = 300)

