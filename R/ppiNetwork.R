#' @importFrom igraph E V E<- V<-
#' @importFrom ggplot2 ggplot guides guide_legend aes rel element_text
#'   element_rect element_blank
#' @import methods
NULL

#' Plot PPI network for gene-set clusters identified using vissE
#'
#' This function plots the protein-protein interaction (PPI) network for a
#' gene-set cluster identified using vissE. The international molecular exchange
#' (IMEx) PPI is used to obtain PPIs for genes present in a gene-set cluster.
#'
#' @param ppidf a data.frame, containing a protein-protein interaction from the
#'   IMEx database. This can be retrieved from the [msigdb::getIMEX()] function.
#' @param threshConfidence a numeric, specifying the confidence threshold to
#'   apply to determine high confidence interactions. This should be a value
#'   between 0 and 1 (default is 0).
#' @param threshFrequency a numeric, specifying the frequency threshold to apply
#'   to determine more frequent genes in the gene-set cluster. The frequecy of a
#'   gene is computed as the proportion of gene-sets to which the gene belongs.
#'   This should be a value between 0 and 1 (default is 0.25).
#' @param threshStatistic a numeric, specifying the threshold to apply to
#'   gene-level statistics (e.g. a log fold-change). This should be a value
#'   between 0 and 1 (default is 0).
#' @param threshUseAbsolute a logical, indicating whether the `threshStatistic`
#'   threshold should be applied to absolute values (default TRUE). This can be
#'   used to threshold on statistics such as the log fold-chage from a
#'   differential expression analysis.
#'
#' @inheritParams plotMsigNetwork
#' @inheritParams plotGeneStats
#'
#' @return a ggplot object with the protein-protein interaction networks plot
#'   for each gene-set cluster.
#' @export
#'
#' @examples
#' data(hgsc)
#' grps = list('early' = 'HALLMARK_ESTROGEN_RESPONSE_EARLY', 'late' = 'HALLMARK_ESTROGEN_RESPONSE_LATE')
#' ppi = msigdb::getIMEX(org = 'hs', inferred = TRUE)
#' plotMsigPPI(ppi, hgsc, grps)
#' 
plotMsigPPI <-
  function(ppidf,
           msigGsc,
           groups,
           geneStat = NULL,
           statName = 'Gene-level statistic',
           threshConfidence = 0,
           threshFrequency = 0.25,
           threshStatistic = 0,
           threshUseAbsolute = TRUE,
           topN = 5,
           nodeSF = 1,
           edgeSF = 1,
           lytFunc = 'graphopt',
           lytParams = list()) {
    checkGroups(groups, names(msigGsc))
    stopifnot(threshConfidence >= 0)
    stopifnot(threshFrequency >= 0)
    stopifnot(nodeSF > 0)
    stopifnot(edgeSF > 0)
    stopifnot(is.null(geneStat) || !is.null(names(geneStat)))
    
    #compute combined graph
    gr = computeMsigGroupPPI(
      ppidf,
      msigGsc,
      groups,
      geneStat,
      threshConfidence,
      threshFrequency,
      threshStatistic,
      threshUseAbsolute,
      topN
    )
    
    #plot base graph
    lytParams = c(list(graph = gr, layout = lytFunc), lytParams)
    p1 = do.call(ggraph::ggraph, lytParams) +
      ggraph::geom_edge_link(
        aes(colour = weight),
        show.legend = FALSE,
        colour = '#666666',
        lwd = 0.5
      ) +
      ggraph::geom_node_point(aes(size = Degree, text = label), colour = '#FFFFFF') +
      ggraph::geom_node_text(aes(label = Label), repel = TRUE) +
      ggraph::geom_node_point(
        aes(fill = geneStat, size = Degree),
        alpha = 0.75,
        shape = 21,
        stroke = 0.2 * nodeSF
      )
    
    #choose colour scale
    if (is.null(geneStat)) {
      p1 = p1 +
        ggplot2::scale_fill_manual(values = c('A' = 1), na.value = 'orchid')
    } else {
      if (threshUseAbsolute & !all(geneStat >= 0)) {
        lims = c(-1, 1) * stats::quantile(abs(geneStat), 0.99)
        palname = 'vik'
      } else {
        lims = stats::quantile(abs(geneStat), c(0.01, 0.99))
        palname = 'tokyo'
      }
      p1 = p1 +
        scico::scale_fill_scico(palette = palname, na.value = 'orchid', limits = lims, oob = scales::squish)
    }
    
    #plot
    p1 = p1 +
      ggplot2::facet_wrap(~ Group, scales = 'free') +
      guides(size = guide_legend(ncol = 2)) +
      ggplot2::scale_size_continuous(range = c(0.1, 6) * nodeSF) +
      bhuvad_theme() +
      ggplot2::theme(
        legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5, size = rel(1.5)),
        panel.background = element_rect(fill = '#FFFFFF', colour = '#000000'),
        strip.background = element_rect(colour = '#000000'),
        axis.text = element_blank(),
        axis.title = element_blank()
      )
    
    return(p1)
  }

computeMsigGroupPPI <- function(ppidf,
                                msigGsc,
                                groups,
                                geneStat = NULL,
                                threshConfidence = 0,
                                threshFrequency = 0.25,
                                threshStatistic = 0,
                                threshUseAbsolute = TRUE,
                                topN = 5) {
  #identify the organism
  idType = msigdb::getMsigIdType(msigGsc)
  org = msigdb::getMsigOrganism(msigGsc, idType)
  
  #check PPI has records for organism of interest
  taxid = c('hs' = '9606', 'mm' = 10090)[org]
  isorg = ppidf$Taxid %in% taxid
  if (any(isorg))
    ppidf = ppidf[isorg, , drop = FALSE]
  else
    stop('No PPIs found for organism, please ensure that the correct PPI database has been loaded.')
  
  #retrieve PPI
  if (is(idType, 'SymbolIdentifier')) {
    colnames(ppidf)[colnames(ppidf) %in% c('SymbolA', 'SymbolB')] = c('from', 'to')
  } else {
    colnames(ppidf)[colnames(ppidf) %in% c('EntrezA', 'EntrezB')] = c('from', 'to')
  }
  
  #compute group-specific PPIs
  ppi_list = lapply(names(groups), function(grpname) {
    gs = groups[[grpname]]
    genes = table(unlist(lapply(msigGsc[gs], GSEABase::geneIds)))
    ig = computeMsigPPI(names(genes), ppidf, threshConfidence)
    
    #remove unconnected nodes
    ig = igraph::induced_subgraph(ig, V(ig)$name[igraph::degree(ig) > 0])
    
    #if graph is not empty
    if (length(V(ig))>  0) {
      V(ig)$label = V(ig)$name
      V(ig)$Freq = genes[V(ig)$name] / length(gs)
      
      #select frequent nodes
      isFreq = V(ig)$Freq >= threshFrequency
      
      #add gene-stats if present
      if (!is.null(geneStat)) {
        V(ig)$geneStat = geneStat[V(ig)$name]
        #filter nodes
        gst = V(ig)$geneStat
        if (threshUseAbsolute)
          gst = abs(gst)
        keepStat = !is.na(gst) & gst >= threshStatistic
      } else {
        V(ig)$geneStat = NA
        keepStat = TRUE #keep all nodes
      }
      
      #create unique node names for each group
      keep = isFreq & keepStat
      V(ig)$name = paste(V(ig)$name, grpname, sep = '_')
      
      if (any(keep)) {
        ig = igraph::induced_subgraph(ig, V(ig)$name[keep])
        #remove unconnected nodes
        ig = igraph::induced_subgraph(ig, igraph::degree(ig) > 0)
      } else {
        ig = igraph::graph_from_data_frame(data.frame(numeric(), numeric()))
      }
    }
    return(ig)
  })
  
  #convert to data.frames
  names(ppi_list) = names(groups)
  nodedf = plyr::ldply(ppi_list, igraph::as_data_frame, what = 'vertices', .id = 'Group')
  edgedf = plyr::ldply(ppi_list, igraph::as_data_frame, what = 'edges', .id = 'Group')
  
  #create NA entry for empty graphs
  nagrps = setdiff(names(groups), nodedf$Group)
  if (length(nagrps) > 0) {
    nodedf = rbind(
      nodedf,
      data.frame(
        'Group' = nagrps,
        'name' = NA,
        'Degree' = NA,
        'label' = NA,
        'Freq' = NA,
        'geneStat' = NA
      )
    )
  }
  
  #add labels
  #identify outliers
  nodedf = plyr::ddply(nodedf, 'Group', function(x) {
    st = rank(x$Degree) * rank(x$Freq)
    if (!is.null(geneStat)) {
      gst = geneStat[x$label]
      st = st * rank(abs(gst))
    }
    islab = !is.na(st) & rank(-st) <= topN
    
    x$Label = ''
    x$Label[islab] = x$label[islab]
    return(x)
  })
  
  #plot base graph
  gr = tidygraph::tbl_graph(nodedf, edgedf, node_key = 'name')
  
  return(gr)
}

computeMsigPPI <- function(genes, ppidf, threshConfidence = 0) {
  colnames(ppidf)[colnames(ppidf) %in% c('Confidence')] = c('weight')
  
  #identify the organism and use the subset PPI
  ppidf = ppidf[ppidf$from %in% genes &
                  ppidf$to %in% genes &
                  ppidf$weight >= threshConfidence, , drop = FALSE]
  ppidf = ppidf[ppidf$from != ppidf$to ,]
  
  #create igraph
  edgecols = c('from', 'to', 'weight')
  ppidf = ppidf[, union(edgecols, colnames(ppidf)), drop = FALSE]
  ppi_ig = igraph::graph_from_data_frame(ppidf)
  igraph::V(ppi_ig)$Degree = igraph::degree(ppi_ig)
  
  return(ppi_ig)
}