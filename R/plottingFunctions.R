#' @importFrom igraph E V E<- V<-
#' @importFrom ggplot2 ggplot guides guide_legend aes rel element_text
#'   element_rect element_blank
#' @import methods
NULL

#' Compute and plot word frequencies for multiple MSigDB collections
#'
#' Given a gene set collection, this function computes the word frequency of
#' gene set names from the Molecular Signatures Database (MSigDB) collection
#' (split by _). Word frequencies are also computed using short descriptions
#' attached with each gene set object.
#'
#' @param groups a named list, of character vectors or numeric indices
#'   specifying node groupings. Each element of the list represent a group and
#'   contains a character vector with node names.
#' @param type a character, specifying the source of text mining. Either gene
#'   set names (`Name`) or descriptions (`Short`) can be used.
#'
#' @inheritParams computeMsigWordFreq
#'
#' @return a ggplot object.
#' @export
#'
#' @examples
#' data("hgsc")
#' groups <- list('g1' = names(hgsc)[1:25], 'g2' = names(hgsc)[26:50])
#' plotMsigWordcloud(hgsc, groups, rmwords = getMsigBlacklist())
#'
plotMsigWordcloud <-
  function(msigGsc,
           groups,
           weight = NULL,
           measure = c('tfidf', 'tf'),
           version = msigdb::getMsigdbVersions(),
           org = c('auto', 'hs', 'mm'),
           rmwords = getMsigBlacklist(),
           type = c('Name', 'Short')) {
    checkGroups(groups, names(msigGsc))
    
    measure = match.arg(measure)
    org = match.arg(org)
    type = match.arg(type)
    
    #add gene set counts for each group
    names(groups) = sapply(names(groups), function(x) {
      paste0(x, ' (n = ', length(groups[[x]]), ')')
    })
    
    #create list of genesets
    msigGsc_list = lapply(groups, function(x) msigGsc[x])
    
    #compute word frequencies
    worddf = plyr::ldply(msigGsc_list, function(x) {
      df = computeMsigWordFreq(x, weight, measure, version, org, rmwords)[[type]]
      df$freq = df$freq / max(df$freq)
      df = df[seq_len(min(30, nrow(df))), ]
      df$angle = sample(c(0, 90), nrow(df), replace = TRUE, prob = c(0.65, 0.35))
      return(df)
    }, .id = 'NodeGroup')
    
    #plot
    p1 = ggplot(worddf, aes(label = word, size = freq, color = freq, angle = angle)) +
      ggwordcloud::geom_text_wordcloud(rm_outside = TRUE, shape = 'circle', eccentricity = 0.65) +
      ggplot2::facet_wrap(~ NodeGroup, scales = 'free') +
      scico::scale_colour_scico(palette = 'acton', direction = -1) +
      ggplot2::scale_size_area(max_size = 6 / log10(1 + length(msigGsc_list))) +
      bhuvad_theme()
    
    return(p1)
  }

#' Plot a gene set overlap network
#'
#' Plots a network of gene set overlap with overlap computed using the
#' [computeMsigOverlap()] and a graph created using [computeMsigNetwork()].
#'
#' @param ig an igraph object, containing a network of gene set overlaps
#'   computed using [computeMsigNetwork()].
#' @param markGroups a named list, of character vectors. Each element of the
#'   list represent a group and contains a character vector with node names. Up
#'   to 12 groups can be visualised in the plot.
#' @param genesetStat a named numeric, statistic to project onto the nodes.
#'   These could be p-values, log fold-changes or gene set score from a
#'   singscore-based analysis.
#' @param nodeSF a numeric, indicating the scaling factor to apply to node
#'   sizes.
#' @param edgeSF a numeric, indicating the scaling factor to apply to edge
#'   widths.
#' @param lytFunc a character, specifying the layout to use (see
#'   `ggraph::create_layout()`).
#' @param lytParams a named list, containing additional parameters needed for
#'   the layout (see `ggraph::create_layout()`).
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' data(hgsc)
#' ovlap <- computeMsigOverlap(hgsc, thresh = 0.15)
#' ig <- computeMsigNetwork(ovlap, hgsc)
#' groups <- list(
#'   'g1' = c("HALLMARK_HYPOXIA", "HALLMARK_GLYCOLYSIS"),
#'   'g2' = c("HALLMARK_INTERFERON_GAMMA_RESPONSE")
#' )
#'
#' plotMsigNetwork(ig, markGroups = groups)
#' 
plotMsigNetwork <-
  function(ig,
           markGroups = NULL,
           genesetStat = NULL,
           nodeSF = 1,
           edgeSF = 1,
           lytFunc = 'graphopt',
           lytParams = list()) {
    stopifnot(nodeSF > 0)
    stopifnot(edgeSF > 0)
    stopifnot(is.null(genesetStat) || !is.null(names(genesetStat)))
    stopifnot(is.null(markGroups) || checkGroups(markGroups, V(ig)$name))
    
    if (length(markGroups) > 12) {
      warning("Only the first 12 components will be plot")
      markGroups = markGroups[seq_len(12)]
    }
    
    #remove unconnected nodes
    ig = igraph::induced_subgraph(ig, V(ig)[igraph::degree(ig) > 0])
    
    #add custom annotation is no category annotated
    if (!all(c('Category') %in% igraph::list.vertex.attributes(ig))) {
      V(ig)$Category = rep('custom', length(V(ig)))
    }
    
    #node colour map
    colmap_nodes = RColorBrewer::brewer.pal(10, 'Set3')
    names(colmap_nodes) = c('archived', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', 'h', 'custom')
    
    #plot base graph
    lytParams = c(list(graph = igraph::as.directed(ig), layout = lytFunc), lytParams)
    p1 = do.call(ggraph::ggraph, lytParams) +
      ggraph::geom_edge_link(
        edge_width = 0.2 * edgeSF,
        alpha = 1 / log10(length(igraph::E(ig))),
        colour = '#66666666'
      ) +
      ggraph::geom_node_point(aes(size = Size), colour = '#FFFFFF')
    
    if (is.null(genesetStat)) {
      p1 = p1 +
        ggraph::geom_node_point(
          aes(fill = Category, size = Size),
          alpha = 0.75,
          shape = 21,
          stroke = 0.2 * nodeSF
        ) +
        ggplot2::scale_fill_manual(values = colmap_nodes) +
        guides(fill = guide_legend(ncol = 4, override.aes = list(size = 5)))
    } else {
      #add stats
      p1$data$genesetStat = genesetStat[p1$data$name]
      
      p1 = p1 +
        ggraph::geom_node_point(
          aes(fill = genesetStat, size = Size),
          alpha = 0.75,
          shape = 21,
          stroke = 0.2 * nodeSF
        )
      if (!all(genesetStat >= 0)) {
        lims = c(-1, 1) * stats::quantile(abs(genesetStat), 0.99)
        palname = 'cork'
        dir = 1
      } else {
        lims = stats::quantile(abs(genesetStat), c(0.01, 0.99))
        palname = 'tokyo'
        dir = -1
      }
      p1 = p1 +
        scico::scale_fill_scico(
          palette = palname,
          na.value = '#AAAAAA',
          limits = lims,
          oob = scales::squish,
          direction = dir
        )
    }
    
    #general theme settings
    p1 = p1 +
      guides(size = guide_legend(ncol = 2)) +
      ggplot2::scale_size_continuous(range = c(0.1, 6) * nodeSF) +
      ggplot2::theme_void() +
      ggplot2::theme(
        legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5, size = rel(1.5))
      )
    
    #mark groups
    if (!is.null(markGroups)) {
      checkGroups(markGroups, p1$data$name)
      
      #add gene set counts for each group
      names(markGroups) = sapply(names(markGroups), function(x) {
        paste0(x, ' (n = ', length(markGroups[[x]]), ')')
      })
      
      #complex hull for groups
      hulldf = plyr::ldply(markGroups, function(x) {
        df = p1$data[p1$data$name %in% x, ]
        df = df[grDevices::chull(df$x, df$y), ]
      }, .id = 'NodeGroup')
      
      p1 = p1 +
        ggforce::geom_shape(
          aes(x, y, colour = NodeGroup),
          fill = NA,
          radius = ggplot2::unit(0.01, 'npc'),
          expand = ggplot2::unit(0.02, 'npc'),
          data = hulldf
        ) +
        ggplot2::scale_colour_brewer(palette = 'Paired') +
        guides(colour = guide_legend(ncol = 4))
    }
    
    return(p1)
  }

#' Plot gene statistics for clusters of gene sets
#'
#' This function plots gene statistics against gene frequencies for any given
#' cluster of gene sets. The plot can be used to identify genes that are
#' over-represented in a cluster of gene-sets (identified based on gene-set
#' overlaps) and have a strong statistic (e.g. log fold-chage or p-value).
#'
#' @param geneStat a named numeric, containing the statistic to be displayed.
#'   The vector must be named with either gene Symbols or Entrez IDs depending
#'   on annotations in msigGsc.
#' @param statName a character, specifying the name of the statistic.
#' @param topN a numeric, specifying the number of genes to label. The top genes
#'   are those with the largest count and statistic.
#'
#' @inheritParams plotMsigWordcloud
#'
#' @return a ggplot object, plotting the gene-level statistic against gene
#'   frequencies in the cluster of gene sets.
#' @export
#'
#' @examples
#' library(GSEABase)
#'
#' data(hgsc)
#' groups <- list('g1' = names(hgsc)[1:25], 'g2' = names(hgsc)[26:50])
#'
#' #create statistics
#' allgenes = unique(unlist(geneIds(hgsc)))
#' gstats = rnorm(length(allgenes))
#' names(gstats) = allgenes
#'
#' #plot
#' plotGeneStats(gstats, hgsc, groups)
#' 
plotGeneStats <- function(geneStat, msigGsc, groups, statName = 'Gene-level statistic', topN = 5) {
  checkGroups(groups, names(msigGsc))
  
  #compute frequencies
  genefreq = plyr::ldply(groups, function (x) {
    gc = table(unlist(lapply(msigGsc[x], GSEABase::geneIds)))
    gc = data.frame('Gene' = names(gc), 'Count' = as.numeric(gc))
    return(gc)
  }, .id = 'Group')
  
  #add stats
  genefreq$GeneStat = geneStat[genefreq$Gene]
  stopifnot(any(!is.na(genefreq$GeneStat)))
  
  #identify outliers
  genefreq = plyr::ddply(genefreq, 'Group', function(x) {
    x = x[!is.na(x$GeneStat) & !is.infinite(x$GeneStat), ]
    st = rank(abs(x$GeneStat)) * rank(x$Count)
    x$rank = rank(-st)
    return(x)
  })
  
  #plot
  p1 = ggplot(genefreq, aes(Count, GeneStat)) +
    ggplot2::geom_jitter(data = genefreq[genefreq$rank > topN, ], shape = '.', colour = 'gray80') +
    ggplot2::geom_jitter(data = genefreq[genefreq$rank <= topN, ]) +
    ggrepel::geom_text_repel(aes(label = Gene), data = genefreq[genefreq$rank <= topN, ]) +
    ggplot2::xlab('Gene-set count') +
    ggplot2::ylab(statName) +
    ggplot2::facet_wrap(~ Group, scales = 'free_x') +
    bhuvad_theme()
  
  return(p1)
}

#' Custom theme
#'
#' @param rl a numeric, scaling factor to apply to text sizes
#'
#' @return a ggplot2 theme
#' @export
#'
#' @examples
#' p1 = ggplot2::ggplot()
#' p1 + bhuvad_theme()
#'
bhuvad_theme = function (rl = 1.1) {
  stopifnot(rl > 0)
  ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.border = element_rect(colour = 'black', fill = NA),
      panel.grid = element_blank(),
      axis.title = element_text(size = rel(rl) * 1.1),
      axis.text = element_text(size = rel(rl)),
      plot.title = element_text(size = rel(rl) * 1.2),
      strip.background = element_rect(fill = NA, colour = 'black'),
      strip.text = element_text(size = rel(rl)),
      legend.text = element_text(size = rel(rl)),
      legend.title = element_text(size = rel(rl), face = 'italic')
    )
}

checkGroups <- function(groups, gscnames) {
  stopifnot(length(groups) > 0)
  stopifnot(all(sapply(groups, length) > 0))
  stopifnot(is.list(groups))
  stopifnot(!is.null(names(groups)))
  lapply(names(groups), function(grpname) {
    if (!all(groups[[grpname]] %in% gscnames))
      stop(sprintf('group "%s" contains unknown members', grpname))
  })
  
  return(TRUE)
}
