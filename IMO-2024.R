library(knitr)
library(mixOmics)
knitr::opts_chunk$set(dpi = 100, echo= TRUE, warning=FALSE, message=FALSE, fig.align = 'center',
                      fig.show=TRUE, fig.keep = 'all', fig.height= 6, fig.width=7)

gene <- read.csv("DIABLO_gene.csv",row.names = 1,header = T)
metabolite <- read.csv("DIABLO_metabolite.csv",row.names = 1,header = T)
group <- read.table("group.txt",header = T)

mydata <- list(G= gene, M= metabolite)


lapply(mydata, dim) # check their dimensions

mygroup <- factor(group$group,levels=c("PBAT","PBS","PLA","PE","PS","CK")) # set the response variable as the Y dataframe
summary(mygroup)


## ---- fig.show = "hold", out.width = "33%", fig.cap = "FIGURE 1: Circle Correlation Plots for pairwise PLS models on the breast TCGA data. Only displays the top 25 features for each dimension, subsetting by those with a correlation above 0.5. "----
list.keepX = c(25, 25) # select arbitrary values of features to keep
list.keepY = c(25, 25)

pls <- spls(gene, metabolite, keepX = list.keepX, keepY = list.keepY) # generate three pairwise PLS models

plotVar(pls, cutoff = 0.5, title = "Gene vs Metabolite", legend = c("Gene", "Metabolite"), # plot features of first PLS
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 20), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))

## ---- echo = FALSE--------------------------------------------------------------------------------------------------
pls1 <- spls(gene, metabolite,, ncomp = 1, keepX = 25, keepY = 25)
cor(pls1$variates$X, pls1$variates$Y)

## -------------------------------------------------------------------------------------------------------------------
mydesign  <-  matrix(0.1, ncol = length(mydata), nrow = length(mydata), # for square matrix filled with 0.1s
                     dimnames = list(names(mydata), names(mydata)))
diag(mydesign) = 0 # set diagonal to 0s

mydesign

## -------------------------------------------------------------------------------------------------------------------
basic.diablo.model <-  block.splsda(X = mydata, Y = mygroup, ncomp = 5, design = mydesign) # form basic DIABLO model


## ---- fig.cap = "FIGURE 2: Choosing the number of components in `block.plsda` using `perf()` with 10 × 10-fold CV function in the `breast.TCGA` study. Classification error rates (overall and balanced, see Section 7.3) are represented on the y-axis with respect to the number of components on the x-axis for each prediction distance presented in PLS-DA"----
perf.diablo = perf(basic.diablo.model, validation = 'Mfold', folds = 10, nrepeat = 10) # run component number tuning with repeated CV

plot(perf.diablo) # plot output of tuning


## -------------------------------------------------------------------------------------------------------------------
ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"] # set the optimal ncomp value
perf.diablo$choice.ncomp$WeightedVote # show the optimal choice for ncomp for each dist metric

# set grid of values for each component to test
test.keepX = list (G = c(5,10,15,20,25,50,100,150,161), 
                   M = c(5,10,15,20,25,50,100,150,154))

#######run the feature selection tuning
tune.TCGA = tune.block.splsda(X = mydata, Y = mygroup, ncomp = 3, 
                              test.keepX = test.keepX, design = mydesign,
                              validation = 'Mfold', folds = 9, nrepeat = 1,
                              dist = "centroids.dist")

list.keepX = tune.TCGA$choice.keepX # set the optimal values of features to retain
list.keepX

# set the optimised DIABLO model
final.diablo.model = block.splsda(X = mydata, Y = mygroup, ncomp = 5, 
                                  keepX = list.keepX, design = mydesign)
final.diablo.model$design # design matrix for the final model
# the features selected to form the first component
selectVar(final.diablo.model, block = 'M', comp = 1)$M$name 

###########
plotDiablo(final.diablo.model, ncomp = 1,
           col.per.group =c("#1F77B4FF","#FF7F0EFF", "#2CA02CFF", "#D62728FF", "#9467BDFF", "#8C564BFF"),
           margin = 0.1)

plotIndiv(final.diablo.model, ind.names = FALSE, legend = TRUE, group = mygroup,
          col.per.group =c("#1F77B4FF","#FF7F0EFF", "#2CA02CFF", "#D62728FF", "#9467BDFF", "#8C564BFF"),
          style = "ggplot2",
          ellipse = TRUE,
          ellipse.level = 0.95,
          title = 'DIABLO Sample Plots')

plotArrow(final.diablo.model, ind.names = FALSE, legend = TRUE, 
          col.per.group =c("#1F77B4FF","#FF7F0EFF", "#2CA02CFF", "#D62728FF", "#9467BDFF", "#8C564BFF"),
          title = 'DIABLO')

plotVar(final.diablo.model, var.names = FALSE, 
        style = 'graphics', legend = TRUE,
        pch = c(16, 17), cex = c(2,2), 
        col = c("#1F77B4FF","#FF7F0EFF"))

circosPlot(final.diablo.model, cutoff = 0.8, line = TRUE,
           color.blocks= c("#1F77B4FF","#FF7F0EFF"),
           color.cor = c("chocolate3","grey20"), size.labels = 1.5)

network(final.diablo.model, blocks = c(1,2),
        color.node = c("#1F77B4FF","#FF7F0EFF"), cutoff = 0.8)

library(igraph)
my.network = network(final.diablo.model, blocks = c(1,2),
                     color.node = c('darkorchid', 'brown1'), cutoff = 0.8)
write.graph(my.network$gR, file = "mynetwork.gml", format = "gml")

