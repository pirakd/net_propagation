library(dplyr, warn.conflicts = FALSE)
library(clusterProfiler, warn.conflicts = FALSE, quietly=TRUE)
library(msigdbr, warn.conflicts = FALSE)
library(data.table,warn.conflicts = FALSE)

simplifyEnrichBySimilarUniverseMembership <- function(enrichResultsTable, gmt, groupColumn='group', cutHeight = 0.99, broadest=TRUE)
{
  if (length(unique(enrichResultsTable$ID)) < 2){
    message ("Nothing to simplify")
    return (list (enrichResultsTable, data.frame()))
  }
  setDT(enrichResultsTable); setDT(gmt)

  ##Select Significant Terms
  target_overrep_sig <- enrichResultsTable[p.adjust < 0.05,]#qvalue < 0.01]
  #target_overrep_sig <- enrichResultsTable[pvalue < 0.05,]#qvalue < 0.01]

  ##Prepare Significant GO Term Jaccard Similarity Matrix
  sig_go_terms <- unique(target_overrep_sig$ID)

  message ("Computing universal gene overlap between ", length(sig_go_terms), " significant GO terms from ", length(unique(enrichResultsTable[[groupColumn]])), " ", groupColumn, "(s)")

  gmt.subset <- gmt[term %in% sig_go_terms, .(term=factor(term), gene=factor(gene))]
  termByGeneMat <-  Matrix::sparseMatrix(as.integer(gmt.subset$term), as.integer(gmt.subset$gene),
                                         dimnames=list(levels(gmt.subset$term),
                                                       levels(gmt.subset$gene)))

  go_dist_mat <- dist(termByGeneMat, method="binary")
  hc <- hclust(go_dist_mat) # you can change the method here.

  clusters <- cutree(hc, h=cutHeight)
  clusters <- data.table (cluster = as.numeric(clusters), ID = attributes(clusters)$names )

  message ("GO terms clustered into ", max(clusters$cluster), " clusters")

  ## go gene set lengths
  gmt.setLengths <- gmt[,.(setSize = length(unique(gene))), by = term]
  clusters <- merge (clusters, gmt.setLengths, by.x="ID", by.y="term")

  clusterInfo <- merge (enrichResultsTable, clusters, by = "ID")

  clusterInfo[,maxSet := max (setSize), by = c("cluster", groupColumn)]
  #winners <- clusterInfo[,.SD[which(count == maxSet)], by = .(cluster, Bait)]  #keeps all ties...needs updating

  if (broadest){
    winners <- clusterInfo[p.adjust < 0.05,.SD[which.max(setSize),],by=c("cluster", groupColumn)]  #chooses the first in case of tie breakers
    message (length(unique(winners$ID)), " representative GO terms choosing the BROADEST significant term per GO-cluster per ", groupColumn)
  }else{
    winners <- clusterInfo[p.adjust < 0.05,.SD[which.min(p.adjust),],by=c("cluster", groupColumn)]  #chooses the first in case of tie breakers
    message (length(unique(winners$ID)), " representative GO terms choosing the MOST significant term per GO-cluster per ", groupColumn)
  }
  result <- enrichResultsTable[ID %in% winners$ID,]
  list(simplified = result, clusterInfo = clusterInfo)
}

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=4) {
  stop("Not all parameters were supplied to R script", call.=FALSE)
}
canonical_pathway_path = args[1]
enrichment_scores_path = args[2]
out_scores_path = args[3]
cut_height = args[4]


gmt <- clusterProfiler::read.gmt(canonical_pathway_path)
enrichTable <- as.data.frame(read.delim(enrichment_scores_path, sep = "\t"))
simp <- simplifyEnrichBySimilarUniverseMembership(enrichTable,
                                                  gmt,
                                                  groupColumn="group",
                                                  cutHeight = cut_height,
                                                  broadest=FALSE)


results = simp$clusterInfo
write.table(results, file = file(out_scores_path, open='w'),
           sep = "\t", row.names=FALSE, col.names=TRUE, append=FALSE, quote=FALSE)