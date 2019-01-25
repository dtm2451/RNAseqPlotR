# RNAseqPlotR
A set of color blind friendly functions built for analysis and visualization of single cell and bulk RNA-sequencing data

**For a tutorial on how to use these functions, see [RNAseqPlotR_Vignette.html](https://github.com/dtm2451/RNAseqPlotR/blob/master/RNAseqPlotR_Vignette.html)**

It includes various helper and plotting functions for working with RNAseq data. I parsonally use the helper functions (especially `meta()` and `get.metas()`) constantly, but I imagine that the main draw for others will be the plotting functions.

All plotting functions spit out default-themed plots upon minimal coding input for daily analysis needs, but they also allow various manipulations to provide for out-of-the-box submission-quality figures as well.

I built the functions while analyzing single cell RNAseq data with Seurat, but since then, I have started to add functionality for handling bulk RNAseq data as well.  Currently, the bulk capabilities only extend to DESeq-analyzed bulk data, but I plan to add functionality for edgeR in the near future.

NOTE: These functions are currently in a functional though not yet prime-time state.  My goal is that they will eventually be released as part of a package that will have all the typical function documentation.  Currently, documentation is in the form of header comments in the function code.  Hopefully, with those and with this tutorial, all the workings of these functions can be figured out!  If you have questions about how to do something, or would like to suggest new features, please message me!

## To use:

Copy code in RNAseqPlotR.R into a script with that name.  Then, use this code to load the functions into your workspace:

```
source("LOCATION/RNAseqPlotR.R")
```

For a tutorial on how to use these functions, see [RNAseqPlotR_Vignette.html](https://github.com/dtm2451/RNAseqPlotR/blob/master/RNAseqPlotR_Vignette.html)

## Plotting Functions

**`DBDimPlot`** = handles all needs for Seurat TSNEPlot / PCAPlot / DimPlot functions.  Improves on the Seurat functions' capabilities to present continuous (including negative) numerical data, or descrete data (clustering, samples, batches, condition, etc.) in various ways.

**`DBPlot`** = handles needs of Seurat's VlnPlot function. Allows generation of jitter/dot-plot, boxplot, and/or violin-plot representation of numerical data, with order of what's on top easily settable. Data can be expression of particular genes or any numerical metadata like percent.mito, nUMI, and nGene.  Colors and grouping of cells is tunable through discrete inputs.

**`DBBarPlot`** = No analogous function currently in Seurat, which is a bit crazy imho. Most common use: Plotting the cluster breakdown of all cells of each sample. Essentially, it is similar to DBPlot, but for discrete variables. Handles plotting of discrete data on a per-sample or per-condition grouping. 

## Color adjustment functions (note: use these on a color.panel, not on a generated plot)

**`Darken`**: Darkens a color or color.panel by a given amount. Usage = `Darken(color, percent.change = 0.25, relative = T)` where color is a hexadecimal string providing the rgb color values (like the values of MYcolors), percent.change is given as a decimal between [0:1] and represents how much you want to change the colors by, and relative is a logical stating whether the change should be relative to the current color amounts versus absolute (a.k.a. relative to the full color range possible.)

**`Lighten`**: Lightens a color or color.panel by a given amount. Usage = `Lighten(color, percent.change = 0.25, relative = T)` where color is a hexadecimal string providing the rgb color values (like the values of MYcolors), percent.change is given as a decimal between [0:1] and represents how much you want to change the colors by, and relative is a logical stating whether the change should be relative to the current color amounts versus absolute (a.k.a. relative to the full color range possible.)

## Helper functions

These make manipulating Seurat data, and using my plotting functons, easier.

**`get.metas()` and `get.genes()`**: Returns the list of meta slots or the list of genes included in the dataset.  Works exactly like typing `names(object@meta.data)` or `rownames(object@raw.data)`, only easier. Usage = `get.metas(object)`, `get.metas("object")`, or `get.metas()` with `DEFAULT` set.

**`is.meta()` and `is.gene()`**: Returns TRUE or FALSE for whether a "meta-name" or "gene" input is part of the dataset. Usage = `is.meta("meta-name", object)`, `is.meta("meta-name", "object")`, or `is.meta("meta-name")` with `DEFAULT` set.

**`meta()` and `gene()`**: Returns the values a metadata for every cell or the normalized expression data (`@data` slot) for all cells. Usage = `meta("meta-name", object)`, `meta("meta-name", "object")`, or `meta("meta-name")` with `DEFAULT` set.

**`meta.levels()`**: Returns the range of values of metadata. Like running `levels(as.factor(object@meta.data$meta))`. Usage = `meta.levels("meta-name", object)`, `meta.levels("meta-name", "object")`, or `meta.levels("meta-name")` with `DEFAULT` set.  Alternatively, can reurn the counts of each value of the meta if the optional input `table.out` is set to `TRUE`.

**`extDim()`**: extracts the loadings of each cell for a given dimensional reduction space. Usage = extDim(reduction.use, dim, object) where: reduction.use is the all lowercase, quoted, key for the target dr space, "tsne", "pca", "ica", "cca", "cca.aligned"; where dim is the component #, like 1 vs 2 for PC1 vs PC2; where object is the "quoted" name of the seurat object of interest.

**`rankedBarcodes`**: Spits out a CellRanger-websummary-like plot of #UMI vs #cells.  Usage = rankedBarcodes(object) with the actual object being given, not the name in quotes.
