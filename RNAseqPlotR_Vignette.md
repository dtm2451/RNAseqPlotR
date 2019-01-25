# RNAseqPlotR

This is a tutorial for using my RNAseqPlotR package of Seurat/DESeq/EdgeR auxiliary functions.  It is currently in the form of an Rscript, [RNAseqPlotR.R](https://github.com/dtm2451/RNAseqPlotR/RNAseqPlot.R).

RNAseqPlotR.R is a color blind friendly package of functions built for analysis and visualization of rnaseq data.  It includes various helper and plotting functions for working with RNAseq data. I personally use the helper functions (especially meta() and get.metas()) constantly, but I imagine that the main draw for others will be the plotting functions.

All plotting functions spit out default-themed plots upon minimal coding input for daily analysis needs, but they also allow various manipulations to provide for out-of-the-box submission-quality figures as well.

I built the functions while analyzing single cell RNAseq data with Seurat.  Since then, I have started to add functionality for handling bulk RNAseq data as well. Currently, bulk capabilities only extend to DESeq-analyzed bulk data. I plan to add functionality for edgeR in the near future. I will add information about handling bulk data with these functions at that time.

NOTE: These functions are currently in a functional though not yet prime-time state.  My goal is that they will eventually be released as part of a package that will have all the typical function documentation.  Currently, documentation is in the form of header comments in the function code.  Hopefully, with those and with this tutorial, all the workings of these functions can be figured out!  If you have questions about how to do something, or would like to suggest new features, my eamil is <daniel.bunis@ucsf.edu>.

## Colors:

The default colors of this package are color blind friendly.  Source: [Wong B, "Points of view: Color blindness." Nature Methods, 2011.](https://www.nature.com/articles/nmeth.1618)  Currently, when you source my package into your workspace, a variable called MYcolors is created that has the 8 colors refenced in the Points of View paper stored inside.  Lighter and darker versions of these same colors are then appended to make it a 24 color vector.

## Plotting Functions

**DBDimPlot** = handles all needs for Seurat TSNEPlot / PCAPlot / DimPlot functions.  Improves on the Seurat functions' capabilities to present continuous (including negative) numerical data, or descrete data (clustering, samples, batches, condition, etc.) in various ways.

**DBPlot** = handles needs of Seurat's VlnPlot function. Allows generation of jitter/dot-plot, boxplot, and/or violin-plot representation of numerical data, with order of what's on top easily settable. Data can be expression of particular genes or any numerical metadata like percent.mito, nUMI, and nGene.  Colors and grouping of cells is tunable through discrete inputs.

**DBBarPlot** = No analogous function currently in Seurat.  Handles plotting of discrete data on a per-sample or per-condition grouping. Essentially, it is similar to DBPlot, but for discrete variables. Example: cluster makeup of a sample.

Loading in the functions...

```
source("~/Desktop/RNAseqPlotR.R")
```

Now, load in your Seurat dataset, and you'll be ready to get going!  If you are working with bulk RNAseq data, there are extra steps given in a later section.  I will eventually build an entire tutorial for this type of data.

```
#Load in Seurat object used for making all example figures in this vignette:
HSPCs <- readRDS("~/Box Sync/Layering Analysis/10X/Savespot1901/HSPCs_cca-alinged.rds")
```

### Basic use

The basic use of most functions, including all of the plotting functions is `funtion(var, object, other.inputs)` where var is the target variable (often this will be the "name" of a metadata or gene) and object is the Seurat-object target.  One of the most useful other inputs is probably cells.use.

**`var`** - When var is given as a string name of a gene or metadata, default plot titles are generally generated and the functions will automatically grab the relevant values to plot.  A discrete vector can also be provided, but note that even if cells.use is going to be used to subset down to only certain cells, this vector must include data for all the cells.

**`object`** - Object is the Seurat or RNAseq-class the function should call on.  It can be entered in any of 3 ways:

- the actual object
- "quoted" name of the object
- or left out completely after setting `DEFAULT <- "object"`

**`DEFAULT`** - Every function in RNAseqPlotR has `object = DEFAULT` as it's default setup. Thus, if you do not provide an `object` input, it will look for a variable named DEFAULT within your workspace. Try this: Store your quoted object name as a variable named DEFAULT, code = `DEFAULT <- "object"`, to skip the object input altogether! For example, if I start by running `DEFAULT <- "HSPCs"`, then a t-SNE plot with datapoints colored by their age can be generated with just `DBDimPlot("age")` rather than with `DBDimPlot("age","HSPCs")`.

Other Required Inputs, **`group.by` and `color.by`** - DBBarPlot(var, object, group.by) and DBPlot(var, object, group.by, color.by) have 1 and 2 other required inputs.  For both, group.by is required in order to set the x axis groupings.  For DBPlot, color.by is required as well for setting the fill color of the violin-plotting or box-plotting.  The inputs for each of these variables are the "quoted" name of a meta.data of the object.  These should be whatever metadata you wish to have the data grouped by / colored by.

#### Basic use examples

```
#DBDimPlot
DBDimPlot("age", "HSPCs")

#DBPlot
DBPlot("percent.mito", "HSPCs", group.by = "Sample", color.by = "age")

#DBBarPlot
DBBarPlot("new.HSPCcelltype", "HSPCs", group.by = "Sample")

#Set DEFAULT
DEFAULT <- "HSPCs" #After setting this, the object slot can be left out entirely!

#DBDimPlot, but allowing it to find "HSPCs" when it looks for DEFAULT
DBDimPlot("age")

#same for DBPlot
DBPlot("percent.mito", group.by = "Sample", color.by = "age")

#same for DBBarPlot
DBBarPlot("new.HSPCcelltype", group.by = "Sample")
```

### Intuitive inputs

"ident" = If "ident" is provided as the `var` then all functions will look for `object@ident`, which is where clustering gets stored in a Seurat object.

"gene-name" or "meta-name" = If a character string is provided that is not "ident", then helper functions is.meta() and is.gene() will be called to determine how to proceed. If the "quoted" name of a metadata slot is given as `var`, then `object@metadata$var` will be used and the plot will be titled "var" by default.  If the "quoted" name of a gene is given as `var`, then `object@data[gene,]` will be used and the plot will be titled "Expression of var" by default.

```
DBDimPlot("ident") # Default: main plot title will also be set to Ident.  y-axis label too for DBPlot.
DBDimPlot("age") # Default: main plot title will be set to 'var'. y-axis label too for DBPlot.
DBDimPlot("CD34") # Default: main plot will be set to "'var' expression". y-axis label too for DBPlot.
```

### Useful Manipulations Examples

A couple examples before jumping in, to showcase the flexible functionality:

```
DBDimPlot("new.HSPCcelltype",
          main = "Easily plot where cells lie in PC/tsne/cca-space, and set all labels",
          sub = "Easily set titles and all labels too!",
          xlab = "PC1 (50%) (actually tSNE1, but this is an example!)",
          ylab = "PC2 (20%)",
          legend.title = "Celltype",
          do.label = T,
          ellipse = T,
          cells.use = meta("new.HSPCcelltype")!="0", colors = c(2:8))
DBPlot("Score",
       group.by = "Sample",
       color.by = "age",
       plots = c("vlnplot","jitter","boxplot"),
       cells.use = meta("new.HSPCcelltype")=="GMP",
       labels = c("Adult-1", "Adult-2", "Fetal-1", "Fetal-2", "Fetal-3", "Cord-1", "Cord-2"),
       jitter.size = 0.5,
       boxplot.color = "white",
       boxplot.fill = F,
       y.breaks = c(-70, -35, 0, 35, 70),
       hline = c(12, -20),
       jitter.width = 0.35,
       sub = "MUCH MORE FLEXIBLE then Seurat's native DimPlot and VlnPlot"
)
```

(Images to come soon.)

### Key advanced inputs:

**`cells.use`, `show.others`** = in DBPlot (cells.use only) and DBDimPlot (both), for showing or highlighting only certain cells. `cells.use` input is used for either showing/highlighting only certain cells with DBPlot and DBDimPlot. For example, only HSCs below.  In DBDimPlot, non-target cells will still be shown, by default, as light grey dots. `show.others = F` is an optional additional input for DBDimPlot which can be used to turn off the inclusion of other cells. Usage = For cells.use, provide either the list of cell.names subsetted to just the ones you want (the same way it is used in Seurat), OR provide a logical that says whether each cell should be included (a.k.a. `THIS` in: `object@cell.names[THIS]`). For show.others, use either TRUE or FALSE.

```
# Using the list of cell names method:
DBPlot("Score", group.by = "Sample", color.by = "age", cells.use = HSPCs@cell.names[meta("new.HSPCcelltype")=="HSC"])
# Using the logical method:
DBDimPlot("age", cells.use = meta("new.HSPCcelltype")=="HSC")
# Removing the gray 'others' dots.
DBDimPlot("age", cells.use = meta("new.HSPCcelltype")=="HSC", show.others = F)
```

**`do.label`, `label.size`, `highlight.labels`** = in DBDimPlot, for labeling clusters / ages / conditions / any discrete classifier of the cells in the dataset.  Setting `do.label = TRUE` turns on labeling. `label.size = 10` adjusts the label size. `highlight.labels = TRUE` controls whether there will be a white box around the labels.  Default are label.size = 5, and highlight.labels = TRUE

**`ellipse`** = in DBDimPlot, whether an ellipse should be drawn to help orient the clustering of certain groups when plotting discrete variables. Usage: `ellipse = TRUE`

```
# Using the default labeling
DBDimPlot("new.HSPCcelltype", do.label = T)
# Adding in ellipses
DBDimPlot("new.HSPCcelltype", do.label = T, ellipse = T)
# Adjusting the label size and removing the highlight/background
DBDimPlot("new.HSPCcelltype", do.label = T, label.size = 10, highlight.labels = F)
```

**`reduction.use`, `dim.1`, `dim.2`** = in DBDimPlot, for setting which dimensional reduction space and dimensions to use. Usage = For reduction.use, provide the all lowercase, quoted, key for the target dr space of the object (`object@dr$"THIS"`), a.k.a. "tsne", "pca", "ica", "cca", "cca.aligned". For dim.1 and dim.2 provide the component #, like 1 and 2 for PC1 and PC2.  dim.1 sets the x-axis and dim.2 sets the y-axis, just like in Seurat's DimPlot functions.

```
# Switching to cca instead of the default, tsne.
DBDimPlot("new.HSPCcelltype", reduction.use = "cca.aligned", dim.1 = 1, dim.2 = 2)
```

**`main`, `sub`, `xlab`, `ylab`, `legend.title`** = These set the titles for the plot, `main` and `sub`, for the axes, `xlab` and `ylab`, and `legend.title` for the legend.  Except for `legend.title` which is currently only used in DBDimPlot, all others work for all my plotting functions.

```
# Adjusting all titles
DBDimPlot("new.HSPCcelltype", reduction.use = "cca.aligned", dim.1 = 1, dim.2 = 2,
          main = "Where cells lie in CCA-space",
          sub = "After CCA-alignment",
          xlab = "CC1",
          ylab = "CC2",
          legend.title = "Celltype")
```

**`plots`** = in DBPlot, this input is used to set the order of data plotting. Options are "vlnplot" for making violin plots, "boxplot" for making box plots, and "jitter" (dot plotting with random spread added in the x direction for visualization).  Usage = `plots = c("vlnplot", "boxplot", "jitter")`. The order that these options are given sets the bottom to top order that they will be plotted.  So actually giving `c("vlnplot", "boxplot", "jitter")` will put the violin plots on the bottom, with a boxplot on top of that and the jitter plotting on the top layer.  Default: plots = c("jitter","vlnplot") as I have found this to be most useful for 10,000+ cell datasets.

```
# This default plotting code...
DBPlot("Score", group.by = "Sample", color.by = "age", cells.use = meta("new.HSPCcelltype")=="HSC")
# ...will produce the same plot as:  (jitter plot behind a violin plot)
DBPlot("Score", group.by = "Sample", color.by = "age", cells.use = meta("new.HSPCcelltype")=="HSC",
       plots = c("jitter", "vlnplot"))
# Changing the order and adding a boxplot in the middle: (violin plot in the back, boxplot on top, jitter on top of that) 
DBPlot("Score", group.by = "Sample", color.by = "age", cells.use = meta("new.HSPCcelltype")=="HSC",
       plots = c("vlnplot", "boxplot", "jitter"))
```

**`color.panel`, `colors`** for setting the colors to be used. Default = a colorblind friendly colorset with 8 colors! Credit to Wong, B. “Points of view: Color blindness” Nature Methods, 2011.  color.panel sets the color options, and colors can be used to adjust which colors are used for what data.  Usage: (plotting.colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000") which is run when the file is sourced); `color.panel = MYcolors`. This is the default. colors = c(1:8) by default or c(8,2:7) if black is wanted for the first grouping instead 

```
#Default of DBDimPlot("new.HSPCcelltype") is the same as
DBDimPlot("new.HSPCcelltype", color.panel = MYcolors, colors = c(1:8))
#And with black (slot 8) first instead...
DBDimPlot("new.HSPCcelltype", colors = c(8,2:7))
```

**`shape` / `shape.by`, `shapes` / `jitter.shapes`** = in DBDimPlot and DBPlot respectively. *shape.by* in DBPlot acts exactly like `color.by` except that instead of setting the fill colors for `plots = "vlnplot"` or `"boxplot"`, it sets the shaping used for `plot = "jitter"`.  Usage = provide the "quoted" name of a metadata that has the discrete options you are looking to have shapes be set by.  *shape* works very similarly in DBDimPlot except that it changes the shapes of the dots in the dimplot AND in that it can additionally be used to change the shape of all dots by simply setting it to a different number corresponding to a ggplot shape from the defult 16 (= filled circles). Usage = `shape = 16` OR `shape = "meta-name"`. *shapes* and *jitter.shapes* provides the set of shapes that should be used. For a list of possibilities, see here <https://www.rstudio.com/wp-content/uploads/2015/03/ggplot2-cheatsheet.pdf>. Default: `shapes=c(16,15,17,23,25,8)` and same for `jitter.shapes=c(16,15,17,23,25,8)`.

```
DBPlot("Score", group.by = "Sample", color.by = "age", cells.use = (meta("new.HSPCcelltype")=="MEP" | meta("new.HSPCcelltype")=="GMP"), plots = c("vlnplot", "jitter"),
       shape.by = "new.HSPCcelltype")
# To reverse the shapes, I change the order of jitter.shapes (default was c(16,15,17,23,25,8))
DBPlot("Score", group.by = "Sample", color.by = "age", cells.use = (meta("new.HSPCcelltype")=="MEP" | meta("new.HSPCcelltype")=="GMP"), plots = c("vlnplot", "jitter"),
       shape.by = "new.HSPCcelltype",
       jitter.shapes = c(15,16))
```

**`range`, `low`, and `high`** = in DBPlot, `range` sets the bounds of the y axis. In DBDimPlot, `range` sets the cutoffs for lowest and highest values that the `low` and `high` colors should be used for. `low` and `high` set the color scale. Usage: `range = c(low-end, high-end)` with low-end and high-end both being numbers. `low` or `high` = "quoted" color name, e.g. "blue", or "quoted" hex code color representation, e.g. "#F0E442". Defaults are the extents of the data for range and yellow&blue fro low&high.

```
# Adjusting the range of DBPlot or DBBarPlot
DBPlot("Score", group.by = "Sample", color.by = "age", cells.use = meta("new.HSPCcelltype")=="HSC",
       plots = c("vlnplot", "boxplot", "jitter"),
       range = c(-40,50))

# Adjusting the range and colors of a DBDimPloot
DBDimPlot("Score", cells.use = meta("new.HSPCcelltype")=="HSC",
          range = c(-20,10), low = "grey", high = "red")
```

**`theme`** = If you would like to use your own theme, you can provide a ggplot theme here.  Provide the theme object (would be labeled as a 'list' in your workspace).  Note: This particular feature has only been lightly tested.

```
#Make a theme object
prettyplot.theme <- theme(text = element_text(size = 14, color="black"),
                    panel.background = element_rect(fill = "transparent",colour = NA),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    axis.text=element_text(color="black"),
                    legend.title = element_text(colour="black", size=14,face="bold"),
                    legend.text = element_text(colour="black", size = 12, face="plain"),
                    plot.background = element_rect(fill = "transparent",colour = NA))
#Make a plot with this theme
DBDimPlot("new.HSPCcelltype", theme = prettyplot.theme)
#Make a plot by calling theme inside the Plotting function call:
DBDimPlot("new.HSPCcelltype",
          theme = theme(text = element_text(size = 14, color="black"),
                        panel.background = element_rect(fill = "transparent",colour = NA),
                        panel.grid.minor = element_blank(),
                        panel.grid.major = element_blank(),
                        axis.text=element_text(color="black"),
                        legend.title = element_text(colour="black", size=14,face="bold"),
                        legend.text = element_text(colour="black", size = 12, face="plain"),
                        plot.background = element_rect(fill = "transparent",colour = NA))
          )
```

**DBPlot specific adjustments** There are many more.

- `labels`: names to change x labels to.
- `hline`: y value(s) where a dashed horizontal line should go. hline = 0.5 or hline = c(0.5, 0.9, 0.99)
- `hline.linetype`: Any ggplot linetype should work.  Default is "dashed"
- `hline.color`: color(s) of the horizontal line(s). hline.color = "black" or hline.color = c("black", "white", "grey")
- `jitter.size`: the size of the jitter shapes. Default is 1.
- `jitter.width`: the width/spread of the jitter in the x direction.  Default is 0.2.
- `jitter.color`: the color of the jitter shapes. Default is "black".
- `jitter.shapes`: the shapes to use.  Default is c(16,15,17,23,25,8) / the first in there which corresponds to dots.
- `jitter.shape.legend.size`: Changes the size of the shape key in the legend.  Use a number.  OR, set to "none" to remove from the legend completely.
- `boxplot.width`: the width/spread of the boxplot in the x direction.  Default is 0.2.
- `boxplot.color`: the color of the lines of the boxplot. Default is "black".
- `boxplot.show.outliers`: whether outliers should by including in the boxplot. If no jitter is being added, this should be set to TRUE. Default is FALSE in order to not have duplicate dots to what's in the jitter. 
- `box.plot.fill`: whether the boxplot should be filled in or not. Default is TRUE. 
- `y.breaks`: a list of breaks that should be used as major gridlines. Usage: y.breaks = c(break1,break2,break3,etc.). NOTE: The low and highs of this variable will override `range`.
- `title.legend`: whether to leave the title for the plot's legend. Default is FALSE, a.k.a. no legend title.

An example using many of these:
```
#Before
DBPlot("Score", group.by = "Sample", color.by = "age", plots = c("vlnplot","jitter","boxplot"), cells.use = meta("new.HSPCcelltype")=="GMP")
#After
DBPlot("Score", group.by = "Sample", color.by = "age", plots = c("vlnplot","jitter","boxplot"), cells.use = meta("new.HSPCcelltype")=="GMP",
       labels = c("Adult-1", "Adult-2", "Fetal-1", "Fetal-2", "Fetal-3", "Cord-1", "Cord-2"),
       jitter.size = 0.5,
       boxplot.color = "white",
       boxplot.fill = F,
       y.breaks = c(-70, -35, 0, 35, 70),
       hline = c(12, -20),
       jitter.width = 0.35
)
```

(Image to come soon.)

## Helper functions

These make manipulating Seurat data, and using my plotting functons, easier.

**`get.metas()` and `get.genes()`**: Returns the list of meta slots or the list of genes included in the dataset.  Works exactly like typing `names(object@meta.data)` or `rownames(object@raw.data)`, only easier. Usage = `get.metas(object)`, `get.metas("object")`, or `get.metas()` with `DEFAULT` set.

**`is.meta()` and `is.gene()`**: Returns TRUE or FALSE for whether a "meta-name" or "gene" input is part of the dataset. Usage = `is.meta("meta-name", object)`, `is.meta("meta-name", "object")`, or `is.meta("meta-name")` with `DEFAULT` set.

**`meta()` and `gene()`**: Returns the values a metadata for every cell or the normalized expression data (`@data` slot) for all cells. Usage = `meta("meta-name", object)`, `meta("meta-name", "object")`, or `meta("meta-name")` with `DEFAULT` set.

**`meta.levels()`**: Returns the range of values of metadata. Like running `levels(as.factor(object@meta.data$meta))`. Usage = `meta.levels("meta-name", object)`, `meta.levels("meta-name", "object")`, or `meta.levels("meta-name")` with `DEFAULT` set.  Alternatively, can reurn the counts of each value of the meta if the optional input `table.out` is set to `TRUE`.

### Other included helper functions

**`extDim()`**: extracts the loadings of each cell for a given dimensional reduction space. Usage = extDim(reduction.use, dim, object) where: reduction.use is the all lowercase, quoted, key for the target dr space, "tsne", "pca", "ica", "cca", "cca.aligned"; where dim is the component #, like 1 vs 2 for PC1 vs PC2; where object is the "quoted" name of the seurat object of interest.

**`rankedBarcodes`**: Spits out a CellRanger-websummary-like plot of #UMI vs #cells.  Usage = rankedBarcodes(object) with the actual object being given, not the name in quotes.
