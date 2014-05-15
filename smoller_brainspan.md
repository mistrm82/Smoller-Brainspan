Smoller Brainspan Analysis
========================================================





Array analysis for Erin Dunn (erindunn@pngu.mgh.harvard.edu) at MGH. Contact Meeta Mistry (mmistry@hsph.harvard.edu) for additional details. Request from client was:

> It doesn't look like BrainCloud is an ideal dataset because of the issues related to the PMI. Given this, I think we should probably focus on BrainSpan and run a parallel set of analyses to see if things look cleaner in that new dataset. 

## Workflow: 
* grab the BrainSpan data set and metadata
* re-run the QC of the metadata and basic expression analysis 
* isolate one brain region only: Orbitofrontal Cortex and Amygdala

## Setup

### Bioconductor and R libraries used


```r
library(ggplot2)
library(gtable)
library(scales)
library(RColorBrewer)
library(GEOquery)
library(affy)
library(arrayQualityMetrics)
library(reshape)
library(xtable)
library(ruv)
library(limma)
library(Biobase)
library(gridExtra)
library(stringr)
library(knitr)
library(png)
library(sva)

source("http://dl.dropboxusercontent.com/u/4253254/Resources/functions.r")
```


### Get variables
- get base directory for analyses
- specify data and results directories
- specify column headers used in metadata file



```r
# Setup directory variables
baseDir <- "."
dataDir <- file.path(baseDir, "data")
metaDir <- file.path(dataDir, "meta")
resultsDir <- file.path(baseDir, "results")
```



### Load the expression data
RMA background corrected data, with quantile normalization and log2-transformation. The median of all probe sets within one gene (transcript cluster) was used as the estimate of gene expression.  In the supplement '...a total of 17,565 mainly protein coding genes were surveyed'


```r

# Load GEO data
gse_gene <- getGEO(filename = file.path(dataDir, "geo/GSE25219-GPL5175_series_matrix.txt.gz"))
# gse_probe <- getGEO(filename=file.path(dataDir,
# 'geo/GSE25219-GPL5188_series_matrix.txt.gz'))
```



### Extract metadata and relevant categories





```r
pheno_gene <- read.delim(file.path(metaDir, "brainspan_samples_metadata.txt"), 
    row.names = 1)
meta_donor <- read.delim(file.path(metaDir, "brainspan_donor_metadata.txt"), 
    row.names = 1)

# Subset data by region and hemisphere
meta_ofc <- pheno_gene[which(pheno_gene$BrainRegion == "OFC" & pheno_gene$Hemisphere == 
    "R"), ]
meta_ofc$PMI <- as.numeric(as.character(meta_ofc$PMI))
meta_ofc$Stage <- factor(meta_ofc$Stage)
meta_amy <- pheno_gene[which(pheno_gene$BrainRegion == "AMY" & pheno_gene$Hemisphere == 
    "R"), ]
```



### Data exploration: focus for now OFC
Phenotype data is loaded in from which we can isolate our two brain regions of interest. In total there are 57 donors, and from each multiple samples were taken from various brain regions. 

Exploring the metadata for consistency, and generating a quick overview:


```r

# Gender distribution
gender <- rbind(table(meta_ofc$Sex), table(meta_amy$Sex))
row.names(gender) <- c("Orbitofrontal cortex", "Amygdala")
kable(gender, format = "html")
```

<table>
 <thead>
  <tr>
   <th>   </th>
   <th> F </th>
   <th> M </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td> Orbitofrontal cortex </td>
   <td> 17 </td>
   <td> 19 </td>
  </tr>
  <tr>
   <td> Amygdala </td>
   <td> 18 </td>
   <td> 19 </td>
  </tr>
</tbody>
</table>

```r

# Age table
age.table <- as.character(unique(meta_donor$Period))
age.table <- strsplit(age.table, ",")
age.table <- do.call("rbind", age.table)
row.names(age.table) <- rep("", nrow(age.table))
colnames(age.table) <- c("Stage", "Description")
kable(data.frame(age.table), format = "html", row.names = FALSE)
```

<table>
 <thead>
  <tr>
   <th> Stage </th>
   <th> Description </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td> 1 </td>
   <td>  Embryonic </td>
  </tr>
  <tr>
   <td> 2 </td>
   <td>  Early fetal </td>
  </tr>
  <tr>
   <td> 3 </td>
   <td>  Early fetal </td>
  </tr>
  <tr>
   <td> 4 </td>
   <td>  Early mid-fetal </td>
  </tr>
  <tr>
   <td> 5 </td>
   <td>  Early mid-fetal </td>
  </tr>
  <tr>
   <td> 6 </td>
   <td>  Late mid-fetal </td>
  </tr>
  <tr>
   <td> 7 </td>
   <td>  Late fetal </td>
  </tr>
  <tr>
   <td> 8 </td>
   <td>  Neonatal and early infancy  </td>
  </tr>
  <tr>
   <td> 9 </td>
   <td>  Late infancy  </td>
  </tr>
  <tr>
   <td> 10 </td>
   <td>  Early childhood  </td>
  </tr>
  <tr>
   <td> 11 </td>
   <td>  Middle and late childhood  </td>
  </tr>
  <tr>
   <td> 12 </td>
   <td>  Adolescence  </td>
  </tr>
  <tr>
   <td> 13 </td>
   <td>  Young adulthood  </td>
  </tr>
  <tr>
   <td> 14 </td>
   <td>  Middle adulthood  </td>
  </tr>
  <tr>
   <td> 15 </td>
   <td>  Late adulthood  </td>
  </tr>
</tbody>
</table>

```r

# Age distribution
ggplot(meta_ofc, aes(Stage)) + geom_bar() + ggtitle("Orbitofrontal Cortex: Age distribution") + 
    xlab("Developmental Stage") + ylab("Number of Samples") + theme(axis.text.x = element_text(colour = "grey20", 
    size = 15, angle = 0, hjust = 0.5, vjust = 0.5, face = "plain"), plot.title = element_text(size = rel(2)), 
    axis.title = element_text(size = rel(1.25)))
```

<img src="figure/testMeta1.png" title="plot of chunk testMeta" alt="plot of chunk testMeta" width="800px" />

```r


# pH measurements
ggplot(na.omit(meta_ofc), aes(x = Stage, y = pH, fill = Stage)) + geom_boxplot() + 
    ggtitle("Orbitofrontal Cortex: pH levels") + xlab("Stage") + guides(fill = FALSE) + 
    theme(plot.title = element_text(size = rel(2)), axis.title = element_text(size = rel(1.5)), 
        axis.text = element_text(size = rel(1.25)))
```

<img src="figure/testMeta2.png" title="plot of chunk testMeta" alt="plot of chunk testMeta" width="800px" />

```r


# Time between death and sample collection; remove non-numeric values
ggplot(na.omit(meta_ofc), aes(x = Stage, y = PMI, fill = Stage)) + geom_boxplot() + 
    ggtitle("Orbitofrontal Cortex: Postmortem Intervals") + xlab("Stage") + 
    ylab("Postmortem Interval") + guides(fill = FALSE) + theme(legend.position = "none", 
    plot.title = element_text(size = rel(2)), axis.title = element_text(size = rel(1.5)), 
    axis.text = element_text(size = rel(1.25)))
```

<img src="figure/testMeta3.png" title="plot of chunk testMeta" alt="plot of chunk testMeta" width="800px" />

```r


# RNA Integrity
ggplot(na.omit(meta_ofc), aes(x = Stage, y = RIN, fill = Stage)) + geom_boxplot() + 
    ggtitle("Orbitofrontal Cortex: RIN") + xlab("Stage") + ylab("Postmortem Interval") + 
    guides(fill = FALSE) + theme(legend.position = "none", plot.title = element_text(size = rel(2)), 
    axis.title = element_text(size = rel(1.5)), axis.text = element_text(size = rel(1.25)))
```

<img src="figure/testMeta4.png" title="plot of chunk testMeta" alt="plot of chunk testMeta" width="800px" />


As we have seen with the Braincloud data, there is a  positive correlation observed with devleopemental stage and PMI and a negative correlations with RIN. Not surprising, fetal samples are more likely to have been obtained faster. The authors also acknolwedge this in the manuscript and account for it by incorporating both as covariates in the model. 

### Quality Control

ArrayQualityMetrics QC report for [GSE25219](./results/report_OFC/index.html)






```r

arrayQualityMetrics(expressionset = eset.ofc, intgroup = c("fetal"), outdir = "./results/report_OFC", 
    force = TRUE, do.logtransform = FALSE)
```



### Clustering of data 

<img src="figure/clusteringDendro.png" title="plot of chunk clusteringDendro" alt="plot of chunk clusteringDendro" width="800px" />



A false color heatmap of the distance between arrays demonstrates a high degree of similarity among fetal samples and likewise with non-fetal samples. Two fetal samples clustering with postnatal, simialr to dendorgram above. Remove these samples.
<img src="figure/image1_.png" title="plot of chunk image1 " alt="plot of chunk image1 " width="800px" style="display: block; margin: auto;" />


Density plots (smoothed histograms) for all arrays follow a similar distribution shape and range.  The shape of distribution is questionable with a fairly heavy right tail.

<img src="figure/image2_.png" title="plot of chunk image2 " alt="plot of chunk image2 " width="800px" style="display: block; margin: auto;" />


### Quick check with raw data
Try checking the quality on a handful of .CEL files and see if we see the same dsitribution. If not, it might be better to work directly from .CEL files.

Reading in : ./data/geo/CEL/GSM703924_HSB92-OFC-R.CEL.gz
Reading in : ./data/geo/CEL/GSM703972_HSB97-OFC-R.CEL.gz
Reading in : ./data/geo/CEL/GSM704048_HSB100-OFC-R.CEL.gz
Reading in : ./data/geo/CEL/GSM704080_HSB102-OFC-R.CEL.gz
Reading in : ./data/geo/CEL/GSM704110_HSB103-OFC-R.CEL.gz
<table>
 <thead>
  <tr>
   <th>   </th>
   <th> SampleName </th>
   <th> BrainCode </th>
   <th> BrainRegion </th>
   <th> Hemisphere </th>
   <th> Sex </th>
   <th> Age </th>
   <th> Stage </th>
   <th> PMI </th>
   <th> pH </th>
   <th> RIN </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td> GSM703924 </td>
   <td> HSB92_OFC_R </td>
   <td> HSB92 </td>
   <td> OFC </td>
   <td> R </td>
   <td> M </td>
   <td> 21 PCW </td>
   <td> 6 </td>
   <td> 4 </td>
   <td> 6.65 </td>
   <td> 9.1 </td>
  </tr>
  <tr>
   <td> GSM703972 </td>
   <td> HSB97_OFC_R </td>
   <td> HSB97 </td>
   <td> OFC </td>
   <td> R </td>
   <td> F </td>
   <td> 17 PCW </td>
   <td> 5 </td>
   <td> 1 </td>
   <td> NA </td>
   <td> 9.9 </td>
  </tr>
  <tr>
   <td> GSM704048 </td>
   <td> HSB100_OFC_R </td>
   <td> HSB100 </td>
   <td> OFC </td>
   <td> R </td>
   <td> F </td>
   <td> 19 PCW </td>
   <td> 6 </td>
   <td> 4 </td>
   <td> 6.56 </td>
   <td> 9.5 </td>
  </tr>
  <tr>
   <td> GSM704080 </td>
   <td> HSB102_OFC_R </td>
   <td> HSB102 </td>
   <td> OFC </td>
   <td> R </td>
   <td> F </td>
   <td> 21 PCW </td>
   <td> 6 </td>
   <td> >13 </td>
   <td> 5.89 </td>
   <td> 9.4 </td>
  </tr>
  <tr>
   <td> GSM704110 </td>
   <td> HSB103_OFC_R </td>
   <td> HSB103 </td>
   <td> OFC </td>
   <td> R </td>
   <td> M </td>
   <td> 13 PCW </td>
   <td> 4 </td>
   <td> 1.5 </td>
   <td> NA </td>
   <td> 9.3 </td>
  </tr>
</tbody>
</table>





The raw data seems better than what we obtained from GEO, with higher signal intensities. Two samples show a slightly wider distribution and another is skewed to the left. Will check how this changes after normalization. 

<img src="figure/image3_.png" title="plot of chunk image3 " alt="plot of chunk image3 " width="800px" style="display: block; margin: auto;" />


The data was normalized for differential gene expression analysis using Robust Multichip Average (RMA) in the oligo BioConductor package. Here, RMA normalizes the intensity values at the probe level, and collapses probes into "core" transcripts based on annotations provided by Affymetrix.


```r
geneSummaries <- rma(affyRaw, target = "core", background = T, normalize = T)
```


### QC after normalization
Repeat the previous QC using the normalized data, and we see that the distributions are similar to what we had found originally pulled from GEO. One option is removing those wide distribution samples.




<img src="figure/image4_.png" title="plot of chunk image4 " alt="plot of chunk image4 " width="800px" style="display: block; margin: auto;" />


### A simple linear model fit including RIN and PMI as covariates
In the Kang et al study, age was evaluated using ANOVA and both RIN and PMI were included as covariates. To stay consistent, we performed an ANCOVA the same on Orbitofrontal cortex data, modeling Age/Developemental Stage as a continuous variable.


```r

# Remove outlier samples
remove <- c(which(label(ddata)$x == 1), which(label(ddata)$x == 2))

# Remove NA values
meta_new <- meta_new[-remove, ]
meta_new <- meta_new[which(!is.na(meta_new$PMI)), ]
data_new <- exprs(eset.ofc)[, rownames(meta_new)]

# Update expression set
exprs(eset.ofc) <- data_new
pData(eset.ofc) <- meta_new

# Model fit
mod <- model.matrix(~Stage + PMI + RIN, pData(eset.ofc))
fit <- lmFit(eset.ofc, mod)
fit <- eBayes(fit)

topStage <- topTable(fit, coef = 2, number = nrow(exprs(eset.ofc)), adjust.method = "BH")
hist(topStage$P.Value, col = "grey", border = F, main = "P-value distribution: Age", 
    xlab = "P-value")
```

<img src="figure/anova1.png" title="plot of chunk anova" alt="plot of chunk anova" width="800px" />

```r

topPMI <- topTable(fit, coef = 3, number = nrow(exprs(eset.ofc)), adjust.method = "BH")
hist(topPMI$P.Value, col = "grey", border = F, main = "P-value distribution: PMI", 
    xlab = "P-value")
```

<img src="figure/anova2.png" title="plot of chunk anova" alt="plot of chunk anova" width="800px" />

```r

topRIN <- topTable(fit, coef = 4, number = nrow(exprs(eset.ofc)), adjust.method = "BH")
hist(topRIN$P.Value, col = "grey", border = F, main = "P-value distribution: RIN", 
    xlab = "P-value")
```

<img src="figure/anova3.png" title="plot of chunk anova" alt="plot of chunk anova" width="800px" />


Alot of differentially expressed genes with Age/Developemntal Stage (3020) even at a quite stringent threshold (padj < 0.001). Suprisingly few significant changes associated with RIN (880) and none associated with postmortem interval.

### Comparing the expression changes of our top hits with age to changes with PMI
Take the top 6 probes that are affected by age. Even though the changes are not identical, we still see a similarity in the trend as we did with the Braincloud data.




<img src="figure/topPlot.png" title="plot of chunk topPlot" alt="plot of chunk topPlot" width="800px" />


### Comparing the expression changes of our top hits with age to changes with RIN

<img src="figure/topPlot2.png" title="plot of chunk topPlot2" alt="plot of chunk topPlot2" width="800px" />


