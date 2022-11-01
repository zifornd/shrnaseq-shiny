theme_set(theme_classic())

#### Functions ####

##### Quality control #####

indexcounts <- function(data) {
  df <- melt(data.frame(t(colSums(data$x$counts))))
  colnames(df)=c("Sample", "Count")
  ggplotly(
    ggplot(df, aes(x=Sample, y=Count)) +  
      geom_bar(stat = "identity", fill="#b8dbcc", color="black") +
      labs(x="", y = "Counts") +
      theme(plot.title = element_text(hjust = 0.5))
  )
}
guidecounts <- function(data) {
  df <- melt(data.frame(t(rowSums(data$x$counts))))
  colnames(df)=c("Guide RNA", "Count")
  ggplotly(
    ggplot(df, aes(x=Count)) +
      geom_density(fill="#b8dbcc", alpha=0.8) +
      labs(x="Counts", y = "Density") +
      theme(plot.title = element_text(hjust = 0.5)) 
  )
}
bcv <- function(data) {
  df <- data.frame(cbind(data$xglm$AveLogCPM, 
                      data$xglm$trended.dispersion, 
                      data$xglm$common.dispersion, 
                      data$xglm$tagwise.dispersion))
  colnames(df) <- c("AveLogCPM", "Trended dispersion",
                 "Common dispersion", "Tagwise dispersion")
  colors <- c("Trended dispersion" = "#b8dbcc", "Common dispersion" = "red", "Tagwise dispersion" = "#b8dbcc")
  
  ggplotly(
    ggplot(df, aes(x=AveLogCPM)) +
      geom_line(aes(y=`Trended dispersion`, color="Trended dispersion")) + 
      geom_line(aes(y=`Common dispersion`, color="Common dispersion")) + 
      geom_point(aes(y=`Tagwise dispersion`, color="Tagwise dispersion")) + 
      labs(x="Average log CPM",
           y = "Biological coefficient of variation", 
           color = "") +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_color_manual(values = colors)
      )
}
mds <- function(data) {
  data$x$samples$group <- factor(data$x$samples$group, levels = c(unique(data$x$samples$group)))
  mat <- cpm(data$x$counts, log = TRUE, prior.count = 1)
  var <- matrixStats::rowVars(mat)
  num <- min(500, length(var))
  ind <- order(var, decreasing = TRUE)[seq_len(num)]
  dst <- dist(t(mat[ind, ]))
  mds <- cmdscale(as.matrix(dst))
  dat <- data.frame(
    MD1 = mds[, 1], 
    MD2 = mds[, 2], 
    group = data$x$samples$group
  )
  
  ggplotly(
    ggplot(dat, aes(MD1, MD2, colour = group)) + 
      geom_point(size = 3) + 
      labs(x = "MDS 1", y = "MDS 2", colour = "") + 
      labs(title="A") +  
      scale_color_brewer(palette = "Set3")
  )
}
cor_mds <- function(data) {
  data$x$samples$group <- factor(data$x$samples$group, levels = c(unique(data$x$samples$group)))
  mat <- data$corrected
  var <- matrixStats::rowVars(mat)
  num <- min(500, length(var))
  ind <- order(var, decreasing = TRUE)[seq_len(num)]
  dst <- dist(t(mat[ind, ]))
  mds <- cmdscale(as.matrix(dst))
  dat <- data.frame(
    MD1 = mds[, 1], 
    MD2 = mds[, 2], 
    group = data$x$samples$group
  )
  
  ggplotly(
    ggplot(dat, aes(MD1, MD2, colour = group)) + 
      geom_point(size = 3) + 
      labs(x = "MDS 1", y = "MDS 2", colour = "") +
      labs(title="B") + 
      scale_color_brewer(palette = "Set3")
  )
}
pca <- function(data) {
  data$x$samples$group <- factor(data$x$samples$group, levels = c(unique(data$x$samples$group)))
  mat <- cpm(data$x$counts, log=TRUE, prior.count = 1)
  var <- matrixStats::rowVars(mat)
  num <- min(500, length(var))
  ind <- order(var, decreasing = TRUE)[seq_len(num)]
  pca <- prcomp(t(mat[ind, ]))
  pct <- (pca$sdev ^ 2) / sum(pca$sdev ^ 2)
  dat <- data.frame(
    PC1 = pca$x[,1],
    PC2 = pca$x[,2],
    group = data$x$samples$group)
  
  ggplotly(
    ggplot(dat, aes_string(x = "PC1", y = "PC2", color = "group")) + 
      geom_point(size = 3) + 
      xlab(paste0("PC1: ", round(pct[1] * 100), "% variance")) + 
      ylab(paste0("PC2: ", round(pct[2] * 100), "% variance")) + 
      coord_fixed() +
      labs(title="A", color="") + 
      scale_color_brewer(palette="Set3")
  )
}
cor_pca <- function(data) {
  data$x$samples$group <- factor(data$x$samples$group, levels = c(unique(data$x$samples$group)))
  var <- matrixStats::rowVars(data$corrected)
  num <- min(500, length(var))
  ind <- order(var, decreasing = TRUE)[seq_len(num)]
  pca <- prcomp(t(data$corrected[ind, ]))
  pct <- (pca$sdev ^ 2) / sum(pca$sdev ^ 2)
  dat <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], group = data$x$samples$group)
  
  ggplotly(
    ggplot(dat, aes_string(x = "PC1", y = "PC2", color = "group")) + 
      geom_point(size = 3) + 
      xlab(paste0("PC1: ", round(pct[1] * 100), "% variance")) + 
      ylab(paste0("PC2: ", round(pct[2] * 100), "% variance")) + 
      coord_fixed() +
      labs(title="B", color="") +
      scale_color_brewer(palette = "Set3")
  )
}
sampledist <- function(data) {
  mat <- cpm(data$x$counts, log=TRUE, prior.count = 1)
  sampleDists <- dist(t(mat))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(data$x$samples$group, data$x$samples$Replicate, sep = " - " )
  colnames(sampleDistMatrix) <- paste(data$x$samples$group, data$x$samples$Replicate, sep = " - " )
  getPalette <- colorRampPalette(brewer.pal(9, "BuGn"))
  heatmaply(sampleDistMatrix, col=getPalette, main = "A") 
}
cor_sampledist <- function(data) {
  sampleDists <- dist(t(data$corrected))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(data$x$samples$group, data$x$samples$Replicate, sep = " - " )
  colnames(sampleDistMatrix) <- paste(data$x$samples$group, data$x$samples$Replicate, sep = " - " )
  getPalette <- colorRampPalette(brewer.pal(9, "BuGn"))
  heatmaply(sampleDistMatrix, col=getPalette, main = "B")
}

##### Differential expression tab #####

hist <- function(data, inputcontrast) {
  obj <- get('data')[[paste0(inputcontrast, "_lrt")]]
  df <- data.frame(obj$table)
  ggplotly(
    ggplot(df, aes(x=PValue)) + 
      geom_histogram(bins=45,fill="#b8dbcc", color="black") +
      labs(x="Guide RNA P values", y = "Frequency") +
      theme(plot.title = element_text(hjust = 0.5)) 
  )
}
plotsmear <- function(data, inputcontrast, FCthres, FDRthres) {
  obj <- get('data')[[paste0(inputcontrast, "_lrt")]]
  top2 <- topTags(obj, n=Inf)
  top2ids <- top2$table[(top2$table$logFC>FCthres | top2$table$logFC<(-(FCthres))),1]
  df <- data.frame(obj$table)
  df$Guide <- rownames(df)
  colors <- c("FDR sig." = "red")
  
  ggplotly(
    ggplot(df, aes(x=logCPM, y=logFC, text=Guide)) +
      geom_point(color = "#b8dbcc") +
      geom_point(data = df %>% filter(row.names(df) %in% top2ids), color = "#000000") +
      geom_point(data = df %>% filter(row.names(df) %in% row.names(top2)[top2$table$FDR<FDRthres]), aes(color = "FDR sig.")) +
      geom_hline(yintercept=(-(FCthres)), linetype="dashed", color="#b8dbcc") +
      geom_hline(yintercept=0,  color="cornflowerblue") +
      geom_hline(yintercept=(FCthres), linetype="dashed", color="#b8dbcc") +
      labs(color="") +
      scale_color_manual(values = colors)
  )
}
volcano <- function(data, inputcontrast, FDRthres) {
  obj <- get('data')[[paste0(inputcontrast, "_lrt")]]
  top2 <- topTags(obj, n=Inf)
  df <- data.frame(obj$table)
  df$Guide <- rownames(df)
  ggplotly(
    ggplot(df, aes(x=logFC, y=-10*log10(PValue), text=Guide)) + geom_point() +
      geom_point() + 
      geom_point(data = df %>% filter(row.names(df) %in% row.names(top2)[top2$table$FDR<FDRthres]), color = "red") +
      labs( x="M", y = "-10*log(P-value)")
  )
}
de <- function(data, inputcontrast, topguides) {
  obj <- get('data')[[paste0(inputcontrast, "_lrt")]]
  top2 <- data.frame(topTags(obj, n=Inf))
  top2 <- top2[order(top2$logFC),]
  selY <- rownames(top2)[c(1:topguides)]
  
  mat <- cpm(data$x$counts, log=TRUE, prior.count = 1)
  colnames(mat) <- paste(data$x$samples$group, data$x$samples$Replicate, sep = " - " )
  mat <- subset(mat, rownames(mat) %in% selY)
  getPalette <- colorRampPalette(brewer.pal(9, "BuGn"))
  heatmaply(mat, col=getPalette, main="A" )
}
cor_de <- function(data, inputcontrast, topguides) {
  obj <- get('data')[[paste0(inputcontrast, "_lrt")]]
  top2 <- data.frame(topTags(obj, n=Inf))
  top2=top2[order(top2$logFC),]
  selY <- rownames(top2)[c(1:topguides)]
  
  colnames(data$corrected) <- paste(data$x$samples$group, data$x$samples$Replicate, sep = " - " )
  corrected <- subset(data$corrected, rownames(data$corrected) %in% selY)
  getPalette <- colorRampPalette(brewer.pal(9, "BuGn"))
  heatmaply(corrected, col=getPalette, main="B")
}
detable <- function(data, inputcontrast, corrected) {
  if (corrected=="Uncorrected") {
    mat <- cpm(data$x$counts, log=TRUE, prior.count = 1)
    colnames(mat) <- paste(data$x$samples$group, data$x$samples$Replicate, sep = " - " )
    mat
  } else {
    colnames(data$corrected) <- paste(data$x$samples$group, data$x$samples$Replicate, sep = " - " )
    data$corrected
  }
}

##### Gene level #####

camera <- function(data, inputcontrast) {
  obj <- get('data')[[paste0(inputcontrast, "_camera")]]
  colnames(obj) <- c("nGuides", "Direction", "Pvalue", "FDR", "Gene")
  obj[,c(5,1:4)]
}
camerarank <- function(data, inputcontrast, FDRthres, s) {
  obj <- get('data')[[paste0(inputcontrast, "_camera")]]
  colnames(obj) <- c("nGuides", "Direction", "Pvalue", "FDR", "Gene")
  res <- as.data.frame(obj)
  
  if (min(res$FDR) > FDRthres) {
    FDRthres=(min(res$FDR)+0.01)
  }
  
  res$Status <- factor("NS", levels = c("Up", "NS", "Down"))
  res$Status[res$Direction == "Up" & res$FDR < FDRthres] <- "Up"
  res$Status[res$Direction == "Down" & res$FDR < FDRthres] <- "Down"
  res$Pvalue <- -log10(res$Pvalue)
  res$Rank <- rank(res$Pvalue)
  
  col <- c(
    "Up"   = "#FF0000",
    "NS"   = "#B8DBCC",
    "Down" = "#6495ED"
  )

  res.s <- res[s, , drop = FALSE]
  
  plt <- ggplot(res, aes(x = Rank, y = Pvalue, colour = Status, size = nGuides, label = Gene)) + 
    geom_point() + 
    geom_point(data = res.s, shape = 1, color="black") + 
    geom_text(data = res.s, hjust = 0, nudge_x = 10, color="black", size=3) +  
    scale_colour_manual(values = col, breaks = names(col)) + 
    labs(
      x = "Rank",
      y = "-log10(Pvalue)",
      colour = "Status"
    ) 
  
  plt <- ggplotly(plt)

  lab <- c(
    "Up"   = sprintf("Up (%s)", comma(sum(res$Status == "Up"))),
    "NS"   = sprintf("NS (%s)", comma(sum(res$Status == "NS"))),
    "Down" = sprintf("Down (%s)", comma(sum(res$Status == "Down")))
  )
  
  plt$x$layout$legend$title$text <- "Status"
  
  for (i in  1:length(unique(res$Status))) {
    plt$x$data[[i]]$name <- lab[plt$x$data[[i]]$name]
  }
  
  plt <- layout(p = plt, legend = list(itemsizing = "constant"))
  plt

}
genelevel <- function(data, inputcontrast) {
  obj <- get('data')[[paste0(inputcontrast, "_genelevel")]]
  colnames(obj) <- c("Gene", "nGuides", 	"Mean logFC", "IQR logFC",
                  "Direction mean logFC",	"Direction smallest Pvalue",
                  "Stouffer's Pvalue",	"Stouffer's FDR")
  row.names(obj) <- NULL
  obj
}
generank <- function(data, inputcontrast, FCthres, s) {
  obj <- get('data')[[paste0(inputcontrast, "_genelevel")]]
  colnames(obj) <- c("Gene", "nGuides", 	"Mean logFC", "IQR logFC",
                     "Direction mean logFC",	"Direction smallest Pvalue",
                     "Stouffer's Pvalue",	"Stouffer's FDR")
  res <- as.data.frame(obj)
  FC_index=c(paste0("> ", FCthres), paste0("< ", FCthres, " and > -", FCthres), paste0("< -", FCthres))
  
  res$LogFC <- factor(FC_index[2], levels = FC_index)
  
  res$LogFC[res$`Mean logFC`>FCthres] <- FC_index[1]
  
  res$LogFC[res$`Mean logFC`<(-(FCthres))] <- FC_index[3]
  
  res$Rank <- rank(res$`Mean logFC`)
  
  col <- c( "#FF0000","#B8DBCC", "#6495ED")
  names(col)=c(FC_index)
  res.s <- res[s, , drop = FALSE]
  
  plt <- ggplot(res, aes(x = Rank, y = `Mean logFC`, colour = LogFC, size = nGuides, label = Gene)) + 
    geom_point() + 
    geom_point(data = res.s, shape = 1, color="black") + 
    geom_text(data = res.s, hjust = 0, nudge_x = 10, color="black", size=3) +  
    scale_colour_manual(values = col, breaks = names(col)) + 
    labs(
      x = "Rank",
      y = "Mean log fold-change",
      colour = "LogFC"
    )
  
  plt <- ggplotly(plt)
  
  lab <- c(sprintf(paste0(FC_index[1], " (%s)"), comma(sum(res$LogFC == FC_index[1]))),
           sprintf(paste0(FC_index[2], " (%s)"), comma(sum(res$LogFC == FC_index[2]))),
           sprintf(paste0(FC_index[3], " (%s)"), comma(sum(res$LogFC == FC_index[3])))
  )
  
  names(lab)=c(FC_index)
  
  plt$x$layout$legend$title$text <- "Log fold-change"
  
  for (i in  1:length(unique(res$LogFC ))) {
    plt$x$data[[i]]$name <- lab[plt$x$data[[i]]$name]
  }
  
  
  plt <- layout(p = plt, legend = list(itemsizing = "constant"))
  
  plt
}
barcode <- function(data, inputcontrast, gene) {
  obj <- get('data')[[paste0(inputcontrast, "_lrt")]]
  genesymbollist <- list()
  genesymbols <- as.character(data$x$genes$Gene)
  unq <- unique(genesymbols)
  unq <- unq[!is.na(unq)]
  
  for (i in unq) {
    sel <- genesymbols == i & !is.na(genesymbols)
      genesymbollist[[i]]  <- which(sel)
  }
  barcodeplot(obj$table$logFC,index=genesymbollist[[gene]],
              labels=c("Negative logFC", "Positive logFC"),
              quantile=c(-0.5,0.5))
}

##### Enrichment #####

gotable <- function(data, inputcontrast) {
  obj <- get('data')[[paste0(inputcontrast, "_go")]]
  obj <- obj[order(obj$P.Up),]
  obj$P.Up <- signif(obj$P.Up,4)
  obj$P.Down <- signif(obj$P.Down,4)
  colnames(obj) <- c("Term", "Ont", "nGenes", "nGenes.Up", "nGenes.Down", "P.Up", "P.Down")
  obj
}
upgoplot <- function(data, inputcontrast, topgos) {
  obj <- get('data')[[paste0(inputcontrast, "_go")]]
  pos <- obj[order(obj$P.Up),]
  pos$logPvalue <- -log(pos$P.Up)
  subset <- pos[c(1:topgos),]
  subset$Term <- str_wrap(subset$Term, width = 20)
  subset$Term <- factor(subset$Term, levels=rev(subset$Term))
  
  ggplotly(
    ggplot(subset, aes(x=logPvalue, y=Term, group=Up)) +
      geom_point(aes(size=Up),color="red") +
      labs(x="-log P value", y = "")
  , height = 600)
}
downgoplot <- function(data, inputcontrast, topgos) {
  obj <- get('data')[[paste0(inputcontrast, "_go")]]
  neg <- obj[order(obj$P.Down),]
  neg$logPvalue <- -log(neg$P.Down)
  subset <- neg[c(1:topgos),]
  subset$Term <- str_wrap(subset$Term, width = 20)
  subset$Term <- factor(subset$Term, levels=rev(subset$Term))
  
  ggplotly(
    ggplot(subset, aes(x=logPvalue, y=Term)) +
      geom_point(aes(size=Down),color="cornflowerblue") +
      labs(x="-log P value", y = "")  
  , height = 600) 
}
keggtable <- function(data, inputcontrast) {
  obj <- get('data')[[paste0(inputcontrast, "_kegg")]]
  obj <- obj[order(obj$P.Up),]
  obj$P.Up <- signif(obj$P.Up,4)
  obj$P.Down <- signif(obj$P.Down,4)
  colnames(obj) <- c("Pathway", "nGenes", "nGenes.Up", "nGenes.Down", "P.Up", "P.Down")
  obj
}
upkeggplot <- function(data, inputcontrast, topkeggs) {
  obj <- get('data')[[paste0(inputcontrast, "_kegg")]]
  pos <- obj[order(obj$P.Up),]
  pos$logPvalue <- -log(pos$P.Up)
  subset <- pos[c(1:topkeggs),]
  subset$Pathway <- str_wrap(subset$Pathway, width = 20)
  subset$Pathway <- factor(subset$Pathway, levels=rev(subset$Pathway))
  
  ggplotly(
    ggplot(subset, aes(x=logPvalue, y=Pathway, group=Up)) +
      geom_point(aes(size=Up),color="red") +
      labs(x="-log P value", y = "")
  , height = 600)
}
downkeggplot <- function(data, inputcontrast, topkeggs) {
  obj <- get('data')[[paste0(inputcontrast, "_kegg")]]
  neg <- obj[order(obj$P.Down),]
  neg$logPvalue <- -log(neg$P.Down)
  subset <- neg[c(1:topkeggs),]
  subset$Pathway <- str_wrap(subset$Pathway, width = 20)
  subset$Pathway <- factor(subset$Pathway, levels=rev(subset$Pathway))
  
  ggplotly(
    ggplot(subset, aes(x=logPvalue, y=Pathway)) +
      geom_point(aes(size=Down),color="cornflowerblue") +
      labs(x="-log P value", y = "")
  , height = 600)
}

##### Comparing contrasts #####

guidecontrast <- function(data, gene) {
  name <- names(data)[which(str_detect(names(data), "_lrt"))]
  split <- str_split(name, "_")
  contrasts <- as.vector(sapply(split,"[[",1))
  compareFC <- NULL
  comparepval <- NULL
  for (i in contrasts) {
    obj <- get('data')[[paste0(i, "_camera")]]
    colnames(obj) <- c("nGuides", "Direction", "Pvalue", "FDR", "Gene")
    
    obj <- obj[obj$Gene==gene,]
    if (nrow(obj)==0) {
      obj[1,"FDR"]<- "NA"
    }
    comparepval <- data.frame(rbind(comparepval, obj$FDR))
    
    obj <- get('data')[[paste0(i, "_lrt")]]
    compareFC <- data.frame(cbind(compareFC,obj$table$logFC))
  }
  
  colnames(compareFC) <- contrasts
  compareFC$Guides <- rownames(obj$table)
  sel=which(str_detect(compareFC$Guides,  (paste0(gene, "\\."))))
  res <- data.frame(compareFC[sel,])
  d <- melt(res)
  colnames(d) <- c("Guides","Contrast","logFC")
  
  comparepval$Contrast <- unique(d$Contrast)
  comparepval$yloc <- max(d$logFC)+0.5 
  comparepval$fdr <- comparepval[,1]
  comparepval$label=c("")
  for (i in 1:nrow(comparepval)) {
    if (is.na(comparepval[i, "fdr"])) {next}
    if (comparepval[i,"fdr"]<0.05) {
      comparepval[i,"label"]="*"
    } 
    if (comparepval[i,"fdr"]<0.01) {
      comparepval[i,"label"]="**"
    }
    if (comparepval[i,"fdr"]<0.001) {
      comparepval[i,"label"]="***"
    }
    }
 
  ggplotly(
    ggplot(d, aes(x=Contrast, y=logFC)) +
      geom_boxplot(outlier.shape = NA, fill="#b8dbcc", show.legend=FALSE) +
      geom_jitter(aes(colour = Guides), show.legend = TRUE) +
      geom_hline(yintercept=0, linetype="dashed", color = "gray") +
      geom_text(data = comparepval, aes(x=Contrast, y = yloc, label =label), vjust=-0.3)
  )
}
FCcontrast <- function(data, inputcontrast) {
  name <- names(data)[which(str_detect(names(data), "_lrt"))]
  split <- str_split(name, "_")
  contrasts <- as.vector(sapply(split,"[[",1))
  vector <- NULL
  dat <- NULL
  for (gene in unique(data$x$genes$Gene)) {
    for (i in contrasts) {
      obj <- get('data')[[paste0(i, "_lrt")]]
      df <- obj$table[which(str_detect(rownames(obj$table),  gene)),]
      vector <- cbind(vector,mean(df$logFC))
    }
    nguides <- nrow(data$x$genes[data$x$genes$Gene==gene,])
    vector <- cbind(gene, nguides, vector)
    dat <- rbind(dat,vector)
    vector <- NULL
  }
  colnames(dat) <- c("Gene", "nGuides", contrasts)
  dat
}
