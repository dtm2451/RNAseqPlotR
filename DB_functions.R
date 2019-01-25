#These are my (sc)RNAseq functions and some things that they rely on.
# DBDimPlot, DBPlot, and DBBarPlot are plotting functions.
# All others are helper functions that can make Seurat and DESeq analysis easier to code!
# Support for other common bulk RNAseq differential expression packages is planned.

MYcolors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
library(ggplot2, quietly = T)
library(Seurat, quietly = T)
library(colorspace, quietly = T)

################## MAIN FUNCTIONS ###################

################################ DBDimPlot ##################################
DBDimPlot <- function(var="ident", object = DEFAULT, reduction.use = NA, dim.1 = 1, dim.2 = 2, theme = NA,
                      size=1, shape=16, shapes=c(16,15,17,23,25,8),
                      legend.size = 5, shape.legend.size = 5, legend.title = NULL, shape.legend.title = NULL,
                      main = NULL, sub = NULL, xlab="make", ylab ="make", auto.title = T,
                      cells.use = NULL, show.others=TRUE, ellipse = F,
                      do.label = F, label.size = 5, highlight.labels = T, label.by = NULL,
                      rename.groups = NA,
                      low = "#F0E442", high = "#0072B2", range = NULL, color.panel = MYcolors, colors = 1:8){
  #Makes a plot where colored dots are overlayed onto the dim.reduction plot of choice.
  #
  #object                 the Seurat or RNAseq object
  #var                    Target Variable = either values or a metadata (in "quotes"), gene (in "quotes"), or "ident"
  #reduction.use          "pca", "tsne", "ica", etc.  Default = tsne for Seurat objects, and pca for RNAseq objects
  #dim.1                  The component number to use on the x-axis.  Default = 1
  #dim.2                  The component number to use on the y-axis.  Default = 2
  #main                   plot title
  #sub                    plot subtitle
  #cells.use              cells to show: either in the form of a character list of names, or a logical that is the same
                          # length as the number of cells in the object (a.k.a. *THIS*: object@cell.names[*THIS*])
  #show.others            TRUE by default, whether other cells should be shown in the background
  #size                   number for size of all highlighted points
  #shape                  number for setting shape OR name of metadata to use for setting shape
  #shapes                 the shapes to use.  Default is a list of 6.  There are more, but not many of the default
                          # ggplot options are great.  I recommend using colors for any variable with 7+ options.
  #low                    color for lowest values of var/range
  #high                   color for highest values of var/range
  #range                  limits for color scaling if var is a continuous value
  #color.panel            a list of colors to be used for when plotting a discrete var.
  #colors                 indexes / order of colors from color.panel to use
  #theme                  Allows setting of a theme. Uses theme_bw if not provided.
                          # To provide, either say "prettyplot" and the prettyplot.1 declared in
                          # this file will be used.  Or provide a your own theme.
  #xlab & ylab            labels for the x and y axes.
  #do.label               Whether to add text labels at the center (median) of clusters for grouping vars
  #label.size             Size of the the labels text
  #highlight.labels       Whether the labels should have a box behind them
  #legend.size            The Size to increase the plotting of legend shapes to (for discrete variable plotting)
  #legend.label           The title for the legend.  Set to NULL if not wanted.
  #rename.groups          new names for the identities of var.  Change to NULL to remove labeling altogether.
  
  #Establish Defaults
  #If cells.use = NA (was not provided), populate it to be all cells or all samples.
  if (classof(object)=="seurat" & is.null(cells.use)) {cells.use <- eval(expr = parse(text = paste0(object,"@cell.names")))}
  if (classof(object)=="RNAseq" & is.null(cells.use)) {cells.use <- meta("Samples", object = object)}
  #If reduction.use = NA (was not provided), populate it to be tsne or pca.
  if (classof(object)=="seurat" & is.na(reduction.use)) {reduction.use <- "tsne"}
  if (classof(object)=="RNAseq" & is.na(reduction.use)) {reduction.use <- "pca"}

  #Build data for populating dat, the data.frame for plotyting.
  #Determine the identity of the provided 'var' and populate Y, the variable used for coloring.
  if(typeof(var)=="character"){
    #If "ident" pull the @ident object from the seurat object
    if(var == "ident"){Y <- eval(expr = parse(text = paste0(object, "@ident")))}
    #If "is.meta" pull the @meta.data$"var" from the RNAseq or seurat object
    if(is.meta(var, object)){Y <- meta(var, object)}
    #If "is.gene" pull the gene expression data from the RNAseq or seurat object
    if(is.gene(var, object)){Y <- gene(var, object)}
    #Otherwise, var is likely a full set of data already, so just make Y = var
  } else {Y <- var}
  #Determine the identity of the provided 'shape'.
    #If it is a meta.data name, pull the meta.  else (it is a number) just carry it through
  if(typeof(shape)=="character"){Shape <- meta(shape, object)} else{Shape <- shape}
  
  #Populate the data.frame to be used for plotting
  full_dat <- data.frame(Y = Y,
                         dim1 = extDim(reduction.use,dim.1,object),
                         dim2 = extDim(reduction.use,dim.2,object),
                         size = size,
                         shape = Shape,
                         Group = Y)

  #Subset to cells.use
  if (classof(object)=="seurat"){
    if (typeof(cells.use)=="logical"){
      Others_dat <- full_dat[!cells.use,]
      Target_dat <- full_dat[cells.use,]
    } else {
    Others_dat <- full_dat[!(eval(expr = parse(text = paste0(object,"@cell.names"))) %in% cells.use),]
    Target_dat <- full_dat[eval(expr = parse(text = paste0(object,"@cell.names"))) %in% cells.use,]
    }
  } 
  if (classof(object)=="RNAseq"){
    if (typeof(cells.use)=="logical"){
      Others_dat <- full_dat[!cells.use,]
      Target_dat <- full_dat[cells.use,]
    } else {
    Others_dat <- full_dat[!(meta("Samples", object) %in% cells.use),]
    Target_dat <- full_dat[meta("Samples", object) %in% cells.use,]
    }
  } 
  
  ###Start building the plot###
  p <- ggplot()
      
  #If not set to NULL and not changed from "make", add default axis labels (ex. "tsne1" or "pca2")
  if (!(is.null(xlab))) {if (xlab=="make") {xlab <- paste0(reduction.use,dim.1)}}
  if (!(is.null(ylab))) {if (ylab=="make") {ylab <- paste0(reduction.use,dim.2)}}
  p <- p + ylab(ylab) + xlab(xlab)
  
  #Then Add more layers:
  
  ###Add the data###
  #Make gray dots on the bottom layer if show.others = T and cells.use is a subset of all the cells / samples.
  if (show.others & dim(Others_dat)[1]>1) {
    p <- p + geom_point(data=Others_dat, aes(x = dim1, y = dim2), size=0.5, color = "gray90")
  }
  #Overlay the target data on top
  # If 'shape' input was the name of a meta.data, aka type=character, treat shape as an aesthetic for performing grouping.
  # Otherwise it is a number and belongs outside of aes.  
  if (typeof(shape)=="character") {
    p <- p + geom_point(data=Target_dat, aes(x = dim1, y = dim2, colour = Y, shape=shape), size=size) +
      scale_shape_manual(values = shapes[1:length(levels(as.factor(Target_dat$shape)))],
                         labels = levels(as.factor(as.character(Target_dat$shape))))
  }  else {
    p <- p + geom_point(data=Target_dat, aes(x = dim1, y = dim2, colour = Y), shape= Shape, size=size)
  }
  
  ###Add ellipse###
  ### Draw an ellipse if ellipse = T.
  if (ellipse) { p <- p + stat_ellipse(data=Target_dat,
                                       aes(x = dim1, y = dim2, colour = Y),
                                       type = "t",
                                       linetype = 2,
                                       size = 0.5
                                       )}
  
  ###Add titles###
  #If not provided, autogenerate based on the identity of var
  if (is.null(main) & auto.title==T){
    #If var is a meta.data, make the title = var
    if(is.meta(var, object)){main <- var}
    #If var is a gene, make the title "Expression of "var
    if(is.gene(var, object)){main <- paste0("Expression of ", var)}
  }
  #If main was provided, use that as the main title
  if (!is.null(main)) {
    p <- p + ggtitle(main, subtitle = sub)
  }
  #If sub was provided, use that as the subtitle
  # if (!is.null(sub)) {
  #   p <- p + ggtitle(subtitle = sub)
  # }
  
  ### Add Labels ###
  if (do.label) {
    #If something different than the "var" is NOT given to 'label.by'
    if (is.null(label.by)){
      #Make a text plot at the median x and y values for each cluster
      #Determine medians
      cent.1 = sapply(levels(as.factor(Target_dat$Y)), function(X) median(Target_dat$dim1[Target_dat$Y==X]))
      cent.2 = sapply(levels(as.factor(Target_dat$Y)), function(X) median(Target_dat$dim2[Target_dat$Y==X]))
      #Add labels
      if (highlight.labels){
        #Add labels with a white background
        p <- p + 
          geom_label(data = data.frame(x=cent.1, y=cent.2),
                     aes(x = x, y = y, label = if(!(is.na(rename.groups[1]))){rename.groups} else {levels(as.factor(Target_dat$Y))}),
                     size = label.size)
      } else {
        #Add labels without a white background
        p <- p + 
          geom_text(data = data.frame(x=cent.1, y=cent.2),
                    aes(x = x, y = y, label = if(!(is.na(rename.groups[1]))){rename.groups} else {levels(as.factor(Target_dat$Y))}),
                    size = label.size)
      }
    } else { #If a distinct variable to label based on was given, use that instead
        #Error: if (!(is.meta(label.by)))
        #Grab the label.by variable
        label.data <- as.factor(meta(label.by, object)[cells.use])
        #Determine medians
        cent.1 = sapply(levels(as.factor(label.data)), function(X) median(Target_dat$dim1[label.data==X]))
        cent.2 = sapply(levels(as.factor(label.data)), function(X) median(Target_dat$dim2[label.data==X]))
        #Add labels
        if (highlight.labels){
          #Add labels with a white background
          p <- p + 
            geom_label(data = data.frame(x=cent.1, y=cent.2),
                       aes(x = x, y = y, label = if(!(is.na(rename.groups[1]))){rename.groups} else {levels(label.data)}),
                       size = label.size)
        } else {
          #Add labels without a white background
          p <- p + 
            geom_text(data = data.frame(x=cent.1, y=cent.2),
                      aes(x = x, y = y, label = if(!(is.na(rename.groups[1]))){rename.groups} else {levels(label.data)}),
                      size = label.size)
        }
    }
  }
  
  ### Set the colors ###
  ### Also change the size of the dots in the legend if showing groupings ###
  #If var yielded a list of groups for plotting (should be in the form of a list of strings = character, or a factor = integer)
  if (typeof(Y)=="character" | typeof(Y)=="integer"){
    #If the number of levels/groups is less than the length of the color.panel, use the color.panel set.
    if (length(levels(as.factor(as.character(Target_dat$Y))))<=length(color.panel[colors])){
      #If labels.rename was changed from NA, rename the grouping labels to rename.groups values
      if (!(is.na(rename.groups[1]))){
        p <- p+ scale_colour_manual(name = legend.title,
                                    values = color.panel[colors],
                                    labels = rename.groups)
      } else {
        #If not, just set the colors and name for the legend key
        p <- p+ scale_colour_manual(name = legend.title,
                                    values = color.panel[colors])
      }
      #Also change the legend properties.
      if (!is.null(legend.size)){
                        #First, change the size of the dots in the legend, unless legend.size was set to NULL
                        #Also, set the legend title.  Given input if given, or remove if still null.
        p <- p + guides(colour = guide_legend(override.aes = list(size=legend.size), title = legend.title),
                        shape = guide_legend(override.aes = list(size=shape.legend.size), title = shape.legend.title)
                        )
      }
    }
  } else {
    #Otherwise, the data is continous, so set a gradient that goes from 'low' input color to 'high' input color.
    p <- p + scale_color_gradient(low= low, high = high, limits = range,
                                  #Next, set the legend title.  Given input if given, or remove if still null.
                                  name = legend.title)
  }
  
  ### Set the theme ###
  #Use theme_bw if 'theme' = NA (was not provided), or use prettyplot.1 if "prettyplot" was provided, or
  # use provided theme if a full one is provided, aka = a list.
  if (is.na(theme)){
    p <- p + theme_bw()
  } else {
    if (theme=="prettyplot") {p <- p + prettyplot.1}
    if (typeof(theme)=="list") {p <- p + theme}
  }
  return(p)
}

################################# DBPlot ####################################

DBPlot <- function(var, object = DEFAULT, main = NULL, sub = NULL, group.by,
                   color.by, shape.by = "", cells.use = NULL, plots = c("jitter","vlnplot"),
                   color.panel = MYcolors, colors = c(1:8),
                   ylab = "make", y.breaks = NULL, range = NULL,
                   xlab = NULL, labels = NULL, rotate.labels = TRUE,
                   hline=NULL, hline.linetype = "dashed", hline.color = "black",
                   jitter.size=1, jitter.width=0.2, jitter.color = "black", jitter.shapes=c(16,15,17,23,25,8),
                   jitter.shape.legend.size = 3,
                   boxplot.width = 0.2, boxplot.color = "black", boxplot.show.outliers = F, boxplot.fill =T,
                   reorder.x = 1:length(meta.levels(group.by, object)),
                   title.legend = F, auto.title=T){
  #Makes a ggplot plot where color/shapes are overlayed onto the dim.reduction plot of choice.
  #
  #object                 the Seurat or RNAseq Object = name of object in "quotes". REQUIRED.
  #var                    Target Variable = values, OR a gene or metadata in "quotes". REQUIRED.
  #group.by               "metadata" to use for separating values.  Default is by sample. REQUIRED.
  #color.by               "metadata" to use for coloring. Affects boxplot and vlnplot fills. REQUIRED for both.
  #shape.by               "metadata" to use for setting the shape of jitter.  Default = just dots. Ignored if not a metadata
  #cells.use              Cells to include: either in the form of a character list of names, or a logical that is the same
                          # length as the number of cells in the object (a.k.a. *USE* in object@cell.names[*USE*])
  #main                   plot title
  #sub                    plot subtitle
  #color.panel            the set of colors to draw from
  #colors                 indexes / or order, of colors from color.panel to actual use
  #plots                  types of plots to include: possibilities = "jitter", "boxplot", "vlnplot"
  #labels                 names to change x labels to.
  #rotate.labels          Logical, whether the labels should be rotated
  #ylab                   y axis label, default is 
  #hline                  y value(s) where a dashed horizontal line should go
  #hline.linetype         Type of line.  Any ggplot linetype should work.  Defaults to "dashed"
  #hline.color            color(s) of the horizontal line(s)
  #jitter.size            the size of the jitter shapes.
  #jitter.width           the width/spread of the jitter in the x direction
  #jitter.color           the color of the jitter shapes
  #jitter.shapes          the shapes to use.  Default is the first in the list, which corresponds to dots.
  #jitter.shape.legend.size     Changes the size of the shape key in the legend.  Use a number.  OR, set to "none"
                                # to remove from the legend completely
  #boxplot.width          the width/spread of the boxplot in the x direction
  #boxplot.color          the color of the lines of the boxplot
  #boxplot.show.outliers  whether outliers should by including in the boxplot. If no jitter is being added, this should
                          # be set to TRUE.  Default is FALSE to not have duplicate dots to what's in the jitter. 
  #box.plot.fill          whether the boxplot should be filled in or not. 
  #reorder.x              sequence of numbers from 1:length(meta.levels(group.by)) for providing a new order
  # for the samples.  Default = alphabetical then numerical.
  #y.breaks               a list of breaks that should be used as major gridlines. c(break1,break2,break3,etc.)
                          # NOTE: The low and highs of this variable will override `range`.
  #title.legend           whether to leave the title for the plot's legend
  
  #If cells.use = NA (was not provided), populate it to be all cells or all samples.
  if (classof(object)=="seurat" & is.null(cells.use)) {cells.use <- eval(expr = parse(text = paste0(object,"@cell.names")))}
  if (classof(object)=="RNAseq" & is.null(cells.use)) {cells.use <- meta("Samples", object)}
  
  ###Determine what the y-axis should be.
  #non-direct input options are: the name of a metadata or a gene (in "quotes").
  # Both these options also get a default axis title.
  if (is.meta(var, object)){Y <- meta(var, object)}
  if (is.gene(var, object)){Y <- gene(var, object)}
  if (!(is.meta(var, object)|is.gene(var, object))){Y <- var}
  #Unless ylab has been changed from default "make" to NULL, set default y-labels if var is "metadata" or "gene".
  if (!(is.null(ylab))){
    if (ylab=="make" & is.meta(var, object)) {ylab <- var}
    if (ylab=="make" & is.gene(var,object)) {ylab <- paste0(var," expression")}
    if (ylab=="make") {ylab <- NULL}
  }
  #Update a holder 'Shape' variable if 'shape.by' is a meta
  if (is.meta(shape.by, object)){Shape <- meta(shape.by, object)} else {Shape = NA}
  
  ###Make dataframe for storing the plotting data:
  #The data =Y, how to group the data = sample, how to color the groupings = color, and what the shape should be if there is a "jitter" made.
  full_dat <- data.frame(Y = Y,
                         sample = meta(group.by, object),
                         color = meta(color.by, object),
                         shape = Shape)
  #Subset the data.frame to only the cells in cell.use.
  if (classof(object)=="seurat"){
    if (typeof(cells.use)=="logical"){
      Target_dat <- full_dat[cells.use,]
    } else {
      Target_dat <- full_dat[eval(expr = parse(text = paste0(object,"@cell.names"))) %in% cells.use,]
    }
  }
  if (classof(object)=="RNAseq"){
    if (typeof(cells.use)=="logical"){
      Target_dat <- full_dat[cells.use,]
    } else {
      Target_dat <- full_dat[meta("Samples", object) %in% cells.use,]
    }
  }
  
  #Reorder x groupings (steps 1 and 2)
  #1-Rename the grouping (=Target_dat$sample) labels in order to set their order.
  #2-Store originals in orig.names.
  #3-Names will be set labels or back to orig.names further down in the code.
  if (typeof(reorder.x)=="integer" | typeof(reorder.x)=="double"){
    #This step is necessary becuase ggplot orders by "character" which would put 10 after 1 instead of 9.
    reorder.x <- as.character(unlist(sapply(reorder.x, function(X) ifelse(X<10,paste0("0",X),X))))
  }
  orig.names <- levels(as.factor(as.character(Target_dat$sample)))
  Target_dat$sample <- as.factor(as.character(Target_dat$sample))
  levels(Target_dat$sample) <- reorder.x
  Target_dat$sample <- as.factor(as.character(Target_dat$sample))
  
  #####Start making the plot
  p <- ggplot(Target_dat, aes(x=sample, y=Y, fill=color))
  #Set the legend to not have a title, unless requested.
  if (!title.legend) {p <- p + theme(legend.title=element_blank())}
  #Add the y label
  p<- p + ylab(ylab)
  #Set the y-axis limits if a range is given.
  if (!is.null(range)){p <- p + ylim(range)}
  
  ###Add data based on what is requested in plots, *ordered by their order*
  for (i in 1:length(plots)){
    #If next request is "boxplot", make a boxplot.
    if (plots[i] == "boxplot") {
      if (boxplot.show.outliers) {
        p <- p + geom_boxplot(width=boxplot.width, color = boxplot.color,
                              alpha = ifelse(boxplot.fill, 1, 0))
      } else {
        p <- p + geom_boxplot(width=boxplot.width, color = boxplot.color,
                              alpha = ifelse(boxplot.fill, 1, 0),
                              outlier.shape = NA)
      }
    }
    #If next request is "jitter", make a jitter.  a.k.a. add dots with y=Y, and randomized spread in x direction.
    if (plots[i] == "jitter") {
      #If shape.by metadata given, use it. Else, shapes[1] which = dots (16) by default
      if (is.meta(shape.by, object)){
        #Make jitter with shapes
        p <- p + geom_jitter(size=jitter.size,
                             width=jitter.width,
                             height = 0,
                             aes(shape = shape),
                             color = jitter.color)
        #Actually set the shapes to jitter.shapes.
        p <- p + scale_shape_manual(values = jitter.shapes[1:length(levels(as.factor(Target_dat$shape)))],
                                    labels = levels(as.factor(as.character(Target_dat$shape))))
        #Also change the size of the shape key in the legend, unless jitter.shape.legend.size was set to NA or "none".
        if (!is.na(jitter.shape.legend.size) & jitter.shape.legend.size!="none"){
          p <- p + guides(shape = guide_legend(override.aes = list(size=jitter.shape.legend.size)))
        }
        if (jitter.shape.legend.size=="none"){
          p <- p + guides(shape = "none")
        }
      } else {
        p <- p + geom_jitter(size=jitter.size,
                             width=jitter.width,
                             height = 0,
                             shape = jitter.shapes[1],
                             color = jitter.color)
      }
    }
    if (plots[i] == "vlnplot") {
      p <- p + geom_violin()
    }
  }
  
  #Add horizontal lines if given.
  if (!is.null(hline))  {p <- p + geom_hline(yintercept=hline, linetype= hline.linetype, color = hline.color)}
  
  #Set default title if no 'main' was given and if 'autotitle' is set to TRUE
  if (is.null(main) & auto.title==T){
    if(is.meta(var, object)){main <- var}
    if(is.gene(var, object)){main <- paste0("Expression of ", var)}
  }
  #Change the x-axis labels if a labels was given.
  if (!(is.null(labels))){
    p <- p + scale_x_discrete(labels=labels)
  } else {
    # If it was not given, put the orig.names back
    p <- p + scale_x_discrete(labels=orig.names[order(reorder.x)])
  }
  #Change y-axis limits/breaks if requested
  if (!is.null(y.breaks)) {
    p <- p + scale_y_continuous(breaks = y.breaks) + coord_cartesian(ylim=c(min(y.breaks),max(y.breaks)))
  }
  #Rotate Labels if rotate.labels = TRUE
  if (rotate.labels) {p <- p + theme(axis.text.x= element_text(angle=45, hjust = 1.3, vjust = 1.2, size=12))}
  #Set colors to color.panel[colors], set the x-axis formatting, and add titles.
  p <- p + scale_fill_manual(values=color.panel[colors]) +
    ggtitle(main, sub) + xlab(xlab)
  
  #DONE. Return the plot
  return(p)
}

########## DBBarPlot: Builds a stacked bar plot to show the composition of samples / ages / 'group.by' ##########
DBBarPlot <- function(var="ident", object = DEFAULT, group.by = "Sample",
                      cells.use = NULL,
                      color.panel = MYcolors, colors = c(1:length(color.panel)),
                      xlab = NULL, ylab = NA, x.labels = NA, rotate.labels = TRUE,
                      y.labels = c(0,0.5,1), label.by = NA,
                      main = NULL, sub = NULL, rename.groups = NA,
                      legend.title = NULL,
                      reorder.x = 1:length(meta.levels(group.by, object))
                      ){
  #This function will build a barplot colored by counts of discrete identities in 'var', and grouped by sample,
  # age, or any other grouping given to 'group.by'.
  #
  #object                 the Seurat or RNAseq object to draw from = name of object in "quotes". REQUIRED.
  #var                    Target Variable = values, OR a metadata in "quotes". REQUIRED. Length must be same 
                          # as the number of cells in the Seurat object or Samples in the RNAseq object
  #group.by               "metadata" to use for separating values.  Default is by sample. REQUIRED.
  #cells.use              Cells to include: either in the form of a character list of names, or a logical that is the same
                          # length as the number of cells in the object (a.k.a. *USE* in object@cell.names[*USE*])
  #color.panel            the set of colors to draw from
  #colors                 indexes / or order, of colors from color.panel to actual use
  #xlab/ylab              The text titles for the axis.  xlab is blank by default.  Provide ylab = NULL to remove
  #x.labels               Replacement x-axis labels to use instead of the identities of gorup.by
  #y.labels               The numerical labels for the y axis.  This will set the range of the plot!
  #label.by               "metaadata" to use for renaming the groupings.  If there are repeats, "-#" will be added.
  #main/sub               main = plot title, sub = plot subtitle.
  #rename.groups          new names for the identities of var.  Change to NULL to remove labeling altogether.
  #legend.label           Title for the legend
  #reorder.x              sequence of numbers from 1:length(meta.levels(group.by)) for providing a new order
                          # for the samples.  Default = alphabetical then numerical.
  
  #If cells.use = NA (was not provided), populate it to be all cells or all samples.
  if (classof(object)=="seurat" & is.null(cells.use)) {cells.use <- eval(expr = parse(text = paste0(object,"@cell.names")))}
  if (classof(object)=="RNAseq" & is.null(cells.use)) {cells.use <- meta("Samples", object)}
  
  ####Retrieve metas: var, group.by
  #var to y.var
    #If name of meta in "quotes", obtain the meta
  if(length(var)==1 & typeof(var)=="character") {
    if (is.meta(var, object)){
      y.var <- meta(var, object)
    }
  }
  #group.by to x.var
    #If name of meta in "quotes", obtain the meta
  if(length(group.by)==1 & typeof(group.by)=="character") {
    if (is.meta(group.by, object)){
      x.var <- as.factor(meta(group.by, object))
    }
  }
  #Subset the x.var and y.var to only the cells in cell.use.
  if (classof(object)=="seurat"){
    if (typeof(cells.use)=="logical"){
      x.var <- as.factor(as.character(x.var[cells.use]))
      y.var <- as.factor(as.character(y.var[cells.use]))
    } else {
      x.var <- x.var[eval(expr = parse(text = paste0(object,"@cell.names"))) %in% cells.use]
      y.var <- y.var[eval(expr = parse(text = paste0(object,"@cell.names"))) %in% cells.use]
    }
  }
  if (classof(object)=="RNAseq"){
    if (typeof(cells.use)=="logical"){
      x.var <- x.var[cells.use]
      y.var <- y.var[cells.use]
    } else {
      x.var <- x.var[meta("Samples", object) %in% cells.use,]
      y.var <- y.var[meta("Samples", object) %in% cells.use]
    }
  }
  #Reorder x groupings (steps 1 and 2)
  #1-Rename the x.var labels in order to set their order.
  #2-Store originals in orig.names.
  #3-Names will be set back to orig.names or x.labels further down on in the code.
  if (typeof(reorder.x)=="integer"){
    reorder.x <- as.character(unlist(sapply(reorder.x, function(X) ifelse(X<10,paste0("0",X),X))))
  }
  orig.names <- levels(x.var)
  levels(x.var) <- reorder.x
  
  #Build data (Make a dataframe while calculating the percent makeup of x.var groups by y.var identities.)
                    #Generate the x.grouping data (needs to be the identities of x.var each individually repeated
                    # the number of times that there are distinct levels in the var / y.var.)
  dat <- data.frame(grouping = c(sapply(levels(as.factor(x.var)), function(X) rep(X, length(levels(as.factor(y.var)))))),
                    #Label the y.ident that this perentage is for. (Used for coloring)
                    y.ident = rep(levels(as.factor(y.var)), length(levels(as.factor(x.var)))),
                    #Generate the percents
                    y.percents = c(sapply(levels(as.factor(x.var)), function(X)
                      unlist(sapply(levels(as.factor(y.var)), function(Y)
                        #Number of Xs that are Ys, divided by the total number of Xs.
                        sum(y.var==Y & x.var == X)/sum(x.var == X)
                        ))
                      ))
  )
  
  #Build Plot
  p <- ggplot(data=dat, aes(x = grouping)) +
    #Add the bars.
    geom_col(aes(y=y.percents, fill = y.ident))
    #Populate ylab title if not provided / left as NA.  If given as NULL,
    # this will be used and will later remove the label.
    if (is.na(ylab) & is.meta(var,object)){
      ylab <- ifelse(is.null(ylab), NULL, paste0("Percent of ", 
                                                 ifelse(classof(object)=="seurat",
                                                        "cells\ncalled as ",
                                                        "reads\ncalled as"),
                                                 var))
    }
    #Add the y-axis labels and name to the plot
    p <- p + scale_y_continuous(breaks= y.labels,
                                limits = c(0,1),
                                name = ylab)
    ###Set the colors & and rename the groupings
    #If labels.rename was changed from NA, rename the labels to rename.groups
    if (!(is.na(rename.groups[1]))){
      p <- p+ scale_fill_manual(name = legend.title,
                                values = color.panel[colors],
                                labels = rename.groups)
    } else {
    #If not, just set the colors and name for the legend key
      p <- p+ scale_fill_manual(name = legend.title,
                                values = color.panel[colors])
    }
    
    #Add the x axis title, plot title and subtitle, 
    p <- p + xlab(xlab) +
      ggtitle(main, subtitle = sub) +
      theme(axis.text.x= element_text(size=10)) +
      if (rotate.labels) {theme(axis.text.x= element_text(angle=45, hjust = 1.3, vjust = 1.2, size=12))} +
      theme(legend.title=element_text(size=12)) +
      theme(legend.text=element_text(size=14))
    
    #Rename the x-axis labels if an x.labels or label.by was given.
    if (!(is.na(x.labels[1]))){
      p <- p + scale_x_discrete(labels=x.labels)
    } else {
      # If it was not given, put the orig.names back
      p <- p + scale_x_discrete(labels=orig.names[order(reorder.x)])
    }
    
    #DONE return the plot
    return(p)
}

#################### Helper Functions ########################

#### is.meta: Is this the name of a meta.data slot in my dataset? ####
is.meta <- function(test, object=DEFAULT){
  #test = the name, in "quotes", that will be tested for being a metadata slot.
  #object = either a Seurat or RNAseq object OR the name of one in "quotes" 

  #Bypass for ident for Seurat objects
  if(test=="ident" & classof(object)=="seurat") {return(TRUE)}
  
  #For all other cases...
  test %in% get.metas(object)
  }

#### is.gene: Is this the name of a gene in my dataset? ####
is.gene <- function(test, object=DEFAULT){
  #test = the name, in "quotes", that will be tested for being a gene in the dataset.
  #object = either a Seurat or RNAseq object OR the name of one in "quotes"
  if(typeof(object)=="character")
  {
    if (classof(object)=="seurat"){
      return(test %in% rownames(eval(expr = parse(text = paste0(object,"@raw.data")))))
      } 
    if (classof(object)=="RNAseq"){
      return(test %in% rownames(eval(expr = parse(text = paste0(object,"@counts")))))
      } 
  } else {
    if (class(object)=="seurat"){
      return(test %in% rownames(object@raw.data))
    }
    if (class(object)=="RNAseq"){
      return(test %in% rownames(object@counts))
    }
  }
}

#### get.metas: prints the names of all the metadata lists for the object ####
get.metas <- function(object=DEFAULT){
  #object = either a Seurat or RNAseq object OR the name of one in "quotes"
  if(typeof(object)=="character"){
    names(eval(expr = parse(text = paste0(object,"@meta.data"))))
  } else {names(object@meta.data)}
}

#### get.genes: prints the names of all the genes for a Seurat or RNAseq ####
get.genes <- function(object=DEFAULT){
  #object = either a Seurat or RNAseq object OR the name of one in "quotes"
  if (classof(object)=="seurat"){
    if(typeof(object)=="character"){
      return(rownames(eval(expr = parse(text = paste0(object,"@raw.data")))))
    } else {return(rownames(object@raw.data))}
  }
  if (classof(object)=="RNAseq"){
    if(typeof(object)=="character"){
      rownames(eval(expr = parse(text = paste0(object,"@counts"))))
    } else {rownames(object@counts)}
  }
}

#### meta: for extracting the values of a particular metadata for all cells/samples ####
meta <- function(meta = "age", object=DEFAULT){
  #meta = the name of a meta.data slot in "quotes"
  #object = either a Seurat or RNAseq object OR the name of one in "quotes"
  if(typeof(object)=="character"){
    if(meta=="ident"){return(as.character(eval(expr = parse(text = paste0(object,"@ident")))))}
    else{return(eval(expr = parse(text = paste0(object,"@meta.data$",meta))))}    
  } else {
    if(meta=="ident"){return(as.character(object@ident))}
    else{return(eval(expr = parse(text = paste0("object@meta.data$",meta))))}
  }
}

#### gene: for extracting the expression values of a particular gene for all cells/samples ####
gene <- function(gene, object=DEFAULT){
  #gene = the name of a gene, in the dataset, in "quotes"
  #object = either a Seurat or RNAseq object OR the name of one in "quotes"
  if (classof(object)=="seurat"){
    ind <- grep(paste0("^",gene,"$"),rownames(eval(expr = parse(text = paste0(object,"@raw.data")))))
    OUT <- eval(expr = parse(text = paste0(object,"@data[gene,]")))
  }
  if (classof(object)=="RNAseq"){
    ind <- grep(paste0("^",gene,"$"),rownames(eval(expr = parse(text = paste0(object,"@counts")))))
    OUT <- eval(expr = parse(text = paste0(object,"@rlog[ind,]")))
  }
  OUT
}

#### meta.levels: for obtaining the different classifications of a meta.data
meta.levels <- function(meta, object = DEFAULT, table.out = F){
  #meta = the name of a meta.data slot in "quotes"
  #object = either a Seurat or RNAseq object OR the name of one in "quotes"
  if (table.out){
    table(meta(meta, object))
  } else {
    levels(as.factor(meta(meta, object)))
  }
}

#### extDim: for extracting the PC/tsne/IC/dim.reduction loadings of all cells/samples ####
extDim <- function(reduction.use, dim=1, object=DEFAULT){
  #reduction.use = the name of the @dr$ slot if the object is a Seurat object.  This is not used for 'RNAseq'
  # objects because they only have pca.
  #dim = the component number to get.  a.k.a. PC*1* versus PC*2*.
  #object = the name of a Seurat or RNAseq object, in "quotes"
  if (classof(object)=="seurat"){
    OUT <- eval(expr = parse(text = paste0(object,"@dr$",reduction.use,"@cell.embeddings[,",dim,"]")))
  }
  if (classof(object)=="RNAseq"){
    OUT <- eval(expr = parse(text = paste0(object,"@pca[[1]]$x[,",dim,"]")))
  }
  OUT
}

#### classof: for determining if 'object' is a Seurat or RNAseq ####
classof <- function (object = DEFAULT){
  #object = the name of a Seurat or RNAseq object, in "quotes"
  class(eval(expr = parse(text = paste0(object)))) 
}

#### rankedBarcodes: recreates the 10X cells vs UMI plot ###########
# Not as pretty, but it gets the job done if you need it.
rankedBarcodes <- function(object){
  #Takes in a Seurat object and spits out a cellranger-like rankedBarcodes plot for all the cells in the object
  UMI.counts <- object@meta.data$nUMI[order(object@meta.data$nUMI, decreasing=TRUE)]
  Barcodes <- 1:length(UMI.counts)
  plot(log10(Barcodes),log10(UMI.counts), type="s")
}

############################################################################################################

#### COLOR PANEL MODIFICATIONS BLOCK ####

############################################################################################################

#### Darken: For darkening colors ####
Darken <- function(colors, percent.change = 0.25, relative = T){
  # Darkens the given color(s) by a set amount.
  #Utilizes lighten and darken functions of the colorspace package to alter colors
  darken(colors, amount = percent.change, space = "HLS", fixup = TRUE, method = ifelse(relative,"relative","absolute"))
}

#### Lighten: For lightening colors ####
Lighten <- function(colors, percent.change = 0.25, relative = T){
  # Lightens the given color(s) by a set amount.
  #Utilizes lighten and darken functions of the colorspace package to alter colors
  lighten(colors, amount = percent.change, space = "HLS", fixup = TRUE, method = ifelse(relative,"relative","absolute"))
}

################################################################################################################

### Bulk RNA-Seq block ###
# The code in this block creates an object for DESeq2 data with similarish structure to a Seurat object.
# = helpful for analyzing bulk and single cell data together
# = helpful for analyzing bulk data if you are used to Seurat structure
# = adds compatibility to my plotting and helper functions for bulk RNAseq data analyzed with DESeq2

library(DESeq2, quietly = T)

Class <- setClass("RNAseq",
                  representation(
                    counts = "matrix",
                    dds = "DESeqDataSet",
                    rlog = "matrix",
                    meta.data = "data.frame",
                    pca = "list",
                    var.genes = "character",
                    samples = "character",
                    exp.filter = "logical",
                    CVs = "numeric"
                    ),
                  prototype(
                    counts = matrix(),
                    dds = new("DESeqDataSet"),
                    rlog = matrix(),
                    meta.data = data.frame(),
                    pca = list(),
                    var.genes = character(),
                    samples = character(),
                    exp.filter = logical(),
                    CVs = double()
                  )
  )

#### NewRNAseq builds a Seurat-like structure object off of a DESeq input.  Also, can run PCA ####
NewRNAseq <- function(dds, #A DESeq object, *the output of DESeq()*
                      run_PCA = FALSE,#If changed to TRUE, function will:
                                           # auto-populate the rlog, var.genes, and PCA fields.
                      pc.genes = NULL,
                      Ngenes = 2500, #How many genes to use for running PCA, (and how many genes
                                     # will be stored in @var.genes)
                      blind = FALSE, #Whether or not the rlog estimation should be blinded to sample info.
                                     # Run `?rlog` for more info
                      counts = NULL  #Raw Counts data, matrix with columns = genes and rows = samples.
                                     # not required, but can be provided.
                      ){
  #INPUTS
  #dds                The DESeq2 object for your data, *the output of DESeq()*
  #run_PCA            FALSE by default.  Whether @rlog @var.genes and @pca[[1]] population are desired.
  #pc.genes           The genes that should be used for calculating PCA.  If null, a per condition expression filter
                      # will be applied, followed by a selection of Ngenes number of genes that have the highest
                      # coefficient of variation (CV=mean/sd)
  #Ngenes             How many genes to use for the PCA
  #blind              Whether rlog estimation should be blinded to sample info. Run `?rlog` for more info
  #counts

  ########## Create the Object ########################
  #Create the object with whatever inputs were given, a.k.a. creates objects@counts and any other level
  # within str(object).  Will all be NULL unless provided in the function call
  object <- new("RNAseq", dds = dds)

  ########## Run Autopopulations ######################
  # Will run by default because this function requires a dds object to be given.
  ##populate dds
  object@dds <- dds
  ##populate @counts
    #Use the data provided if it was, otherwise, grab from the dds
  if (!(is.null(counts))) {object@counts <- counts
    } else {object@counts <- counts(dds)}
  ##populate @samples
  object@samples <- colnames(object@counts)
  ##populate some of @meta.data
  #   1st add samples, then Nreads.
  object@meta.data <- data.frame(Samples = object@samples,
                                 Nreads = colSums(object@counts))
  rownames(object@meta.data) <- object@samples

  ##Also add colData from dds to @meta.data slot
  #Turn colData into a data.frame, and merge that with current meta.data, BUT do not include any
  # dublicate sets.  For example, Samples will be ignored in colData because it was already grabbed
  # from the counts matrix
  object@meta.data <- cbind(object@meta.data,
                            data.frame(object@dds@colData@listData)[(!duplicated(
                              c(names(data.frame(object@dds@colData@listData)),names(object@meta.data)),
                              fromLast=T
                              ))[1:length(object@dds@colData@listData)]])
  ########## Will run if run_PCA = TRUE ##################
  ##Will populate: rlog, pca, var.genes
  if (run_PCA){
    ##populate rlog
    object@rlog <- assay(rlog(object@dds, blind = blind))
    #Filter rlog to only the genes expressed in at least 75% of samples from test group (ONLY WORKS FOR ONE TEST GROUP)
      test_meta <- strsplit(as.character(object@dds@design), split = "~")[[2]][1]
      #Store this metadata as an easily accessible variable to speed up the next step.
      classes <- meta(test_meta, object)
      #For each gene, return TRUE if... the gene is expressed in >75% of samples from each condition used to build the dds
      ##populate exp.filter
      object@exp.filter <- sapply(1:dim(object@counts)[1], function(X)
        #Ensures that the #of classifications = the number TRUEs in what the nested sapply produces
        length(levels(as.factor(classes)))==sum(
          #For each classification of the test variable, check for >= 75% expression
          #This half sets a variable to each of the saparate classifications,
          # and says to run the next lines on each
          sapply(levels(as.factor(classes)), function(Y)
            #This part of the function determine how many of the samples express the gene
            (sum(object@counts[X, classes==Y]>0))
            #This part compares the above to (the number of samples of the current classification*75%)
             >= (sum(classes==Y)*.75)
          )
        )
      )
    data_for_prcomp <- as.data.frame(object@rlog)[object@exp.filter,]
    #calculate CV by dividing mean by sd
    ## populate CVs
    object@CVs <- apply(X = object@rlog, MARGIN = 1, FUN = sd)/apply(X = object@rlog, MARGIN = 1, FUN = mean)
    #Trim rlog data and RawCV_rlog by expression filter variable = object@exp.filter
    #arrange by CV_rank, higher CVs first
    data_for_prcomp<- data_for_prcomp[order(object@CVs[object@exp.filter], decreasing = T),]
    ##populate var.genes
    object@var.genes <- rownames(data_for_prcomp)[1:(min(Ngenes,dim(data_for_prcomp)[1]))]
    ##populate pca : Run PCA on the top Ngenes CV genes that survive a 75% expression per condition filter
    object@pca <- list(prcomp(t(data_for_prcomp[1:Ngenes,]), center = T, scale = T))
  }
  #OUTPUT: (This is how functions "work" in R.  The final line is what they return.)
  object
}


PCAcalc <- function(object = DEFAULT,
                    genes.use = NULL,
                    Ngenes = 2500, #How many genes to use for running PCA, (and how many genes
                    # will be stored in @var.genes)
                    blind = FALSE #Whether or not the rlog estimation should be blinded to sample info.
                    # Run `?rlog` for more info
                    ){
  #INPUTS
  #object         The RNAseq object to work on
  #genes.use      The genes that should be used for calculating PCA.  If null, a per condition expression filter
                  # will be applied, followed by a selection of Ngenes number of genes that have the highest
                  # coefficient of variation (CV=mean/sd)
  #run_PCA        FALSE by default.  Whether @rlog @var.genes and @pca[[1]] population are desired.
  #Ngenes         How many genes to use for the PCA
  #blind          Whether rlog estimation should be blinded to sample info. Run `?rlog` for more info

  ##populate rlog
  object@rlog <- assay(rlog(object@dds, blind = blind))
  
  ######### IF no genes.use given, use the CVs and ExpFilter to pick genes ###########
  if(is.null(genes.use)){
    #Filter rlog to only the genes expressed in at least 75% of samples from test group (ONLY WORKS FOR ONE TEST GROUP)
    test_meta <- strsplit(as.character(object@dds@design), split = "~")[[2]][1]
    #Store this metadata as an easily accessible variable to speed up the next step.
    classes <- meta(test_meta, object)
    #For each gene, return TRUE if... the gene is expressed in >75% of samples from each condition used to build the dds
    ##populate exp.filter
    object@exp.filter <- sapply(1:dim(object@counts)[1], function(X)
      #Ensures that the #of classifications = the number TRUEs in what the nested sapply produces
      length(levels(as.factor(classes)))==sum(
        #For each classification of the test variable, check for >= 75% expression
        #This half sets a variable to each of the saparate classifications,
        # and says to run the next lines on each
        sapply(levels(as.factor(classes)), function(Y)
          #This part of the function determine how many of the samples express the gene
          (sum(object@counts[X, classes==Y]>0))
          #This part compares the above to (the number of samples of the current classification*75%)
          >= (sum(classes==Y)*.75)
        )
      )
    )
    data_for_prcomp <- as.data.frame(object@rlog)[object@exp.filter,]
    #calculate CV by dividing mean by sd
    ## populate CVs
    object@CVs <- apply(X = object@rlog, MARGIN = 1, FUN = sd)/apply(X = object@rlog, MARGIN = 1, FUN = mean)
    #Trim rlog data and RawCV_rlog by expression filter variable = object@exp.filter
    #arrange by CV_rank, higher CVs first
    data_for_prcomp<- data_for_prcomp[order(object@CVs[object@exp.filter], decreasing = T),]
    ##populate var.genes
    object@var.genes <- rownames(data_for_prcomp)[1:(min(Ngenes,dim(data_for_prcomp)[1]))]
    ##populate pca : Run PCA on the top Ngenes CV genes that survive a 75% expression per condition filter
    object@pca <- list(prcomp(t(data_for_prcomp[1:Ngenes,]), center = T, scale = T))
  }
  ######### IF genes.use given, use them  ###########
  if(!(is.null(genes.use))){
    #Filter rlog to only the genes in genes.use
    data_for_prcomp <- as.data.frame(object@rlog)[genes.use,]
    ##populate pca : Run PCA on the top Ngenes CV genes that survive a 75% expression per condition filter
    object@pca <- list(prcomp(t(data_for_prcomp), center = T, scale = T))
  }
  object
}

MYcolors <- c(MYcolors, Darken(MYcolors, 0.25), Lighten(MYcolors,0.1))
