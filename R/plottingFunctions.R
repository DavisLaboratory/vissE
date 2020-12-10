#' @importFrom igraph E V E<- V<-
#' @importFrom ggplot2 ggplot guides guide_legend aes rel element_text
#'   element_rect element_blank
#' @import methods
NULL

#' Compute word frequencies for MSigDB collections
#'
#' Given a gene set collection, this function computes the word frequency of
#' gene set names from the Molecular Signatures Database (MSigDB) collection
#' (split by _). Word frequencies are also computed using short descriptions
#' attached with each gene set object.
#'
#' @param msigGsc a GeneSetCollection object, containing gene sets from the
#'   MSigDB. The [GSEABase::getBroadSets()] function can be used to parse XML
#'   files downloaded from MSigDB.
#' @param groups a named list, of character vectors or numeric indices
#'   specifying node groupings. Each element of the list represent a group and
#'   contains a character vector with node names.
#' @param rmwords a character vector, containing a blacklist of words to discard
#'   from the analysis.
#' @param type a character, specifying the source of text mining. Either gene
#'   set names (`Name`) or descriptions (`Short`) can be used.
#'
#' @return a ggplot object.
#' @export
#'
#' @examples
#'
#' data("hgsc")
#' groups <- list('g1' = 1:10, 'g2' = 11:20)
#' plotMsigWordcloud(hgsc, groups, getMsigBlacklist())
#'
plotMsigWordcloud <-
  function(msigGsc,
           groups,
           rmwords = getMsigBlacklist(),
           type = c('Name', 'Short')) {
  stopifnot(is.list(groups))
  type = match.arg(type)

  #add gene set counts for each group
  names(groups) = sapply(names(groups), function(x) {
    paste0(x, ' (n = ', length(groups[[x]]), ')')
  })

  #create list of genesets
  msigGsc_list = lapply(groups, function(x) msigGsc[x])

  #compute word frequencies
  worddf = plyr::ldply(msigGsc_list, function(x) {
    df = computeMsigWordFreq(x, rmwords)[[type]]
    df$col = df$freq / max(df$freq)
    df = df[1:min(30, nrow(df)), ]
    df$angle = sample(c(0, 90), nrow(df), replace = TRUE, prob = c(0.65, 0.35))
    return(df)
  }, .id = 'NodeGroup')

  #plot
  p1 = ggplot(worddf, aes(label = word, size = freq, color = col, angle = angle)) +
    ggwordcloud::geom_text_wordcloud(rm_outside = TRUE, shape = 'circle', eccentricity = 0.65) +
    ggplot2::facet_wrap(~ NodeGroup, scales = 'free') +
    scico::scale_colour_scico(palette = 'acton', direction = -1) +
    ggplot2::scale_size_area(max_size = 12) +
    current_theme()

  return(p1)
}

#' Plot a gene set overlap network
#'
#' Plots a network of gene set overlap with overlap computed using the
#' [computeMsigOverlap()] and a graph created using [computeMsigNetwork()].
#'
#' @param ig an igraph object, containing a network of gene set overlaps
#'   computed using [computeMsigNetwork()].
#' @param markGroups a named list, of character vectors or numeric indices
#'   specifying node groupings. Each element of the list represent a group and
#'   contains a character vector with node names. Up to 12 groups can be
#'   visualised in the plot.
#' @param nodeSF a numeric, indicating the scaling factor to apply to node
#'   sizes.
#' @param edgeSF a numeric, indicating the scaling factor to apply to edge
#'   widths.
#' @param lytFunc a function, that computes layouts and returns a matrix with 2
#'   columns specifying the x and y coordinates of nodes. Layout functions in
#'   the igraph package can be used here.
#' @param lytParams a named list, containing additional parameters to be passed
#'   on to the layout function.
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#'
#' data(hgsc)
#' ovlap <- computeMsigOverlap(hgsc)
#' ig <- computeMsigNetwork(ovlap, hgsc)
#' groups <- list('g1' = c(1, 10, 11, 14, 21), 'g2' = c(3, 5, 13))
#'
#' plotMsigNetwork(ig, markGroups = groups)
#'
plotMsigNetwork <-
  function(ig,
           markGroups = NULL,
           enrichStat = NULL,
           nodeSF = 1,
           edgeSF = 1,
           lytFunc = igraph::layout_with_graphopt,
           lytParams = list()) {
  stopifnot(nodeSF > 0)
  stopifnot(edgeSF > 0)
  stopifnot(is.function(lytFunc))
  stopifnot(is.null(markGroups) | is.list(markGroups))
  stopifnot(is.null(enrichStat) | !is.null(names(enrichStat)))

  if (length(markGroups) > 12) {
    warning("Only the first 12 components will be plot")
    markGroups = markGroups[1:12]
  }

  #remove unconnected nodes
  ig = igraph::induced_subgraph(ig, V(ig)[igraph::degree(ig) > 0])

  #compute layout if none present
  if (!all(c('x', 'y') %in% igraph::list.vertex.attributes(ig))) {
    lytParams = c(list('graph' = ig), lytParams)
    lyt = do.call(lytFunc, lytParams)
    V(ig)$x = lyt[, 1]
    V(ig)$y = lyt[, 2]
  }
  if (!all(c('Category') %in% igraph::list.vertex.attributes(ig))) {
    V(ig)$Category = rep('custom', length(V(ig)))
  }

  #convert to dataframe
  nodedf = igraph::as_data_frame(ig, what = 'vertices')
  edgedf = igraph::as_data_frame(ig, what = 'edges')

  #node colour map
  colmap_nodes = RColorBrewer::brewer.pal(10, 'Set3')
  names(colmap_nodes) = c('archived', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', 'h', 'custom')

  #order nodes
  nodedf = nodedf[order(nodedf$Category), ]

  #add edge coordinates
  edgedf$x_begin = nodedf[edgedf$from, 'x']
  edgedf$x_end = nodedf[edgedf$to, 'x']
  edgedf$y_begin = nodedf[edgedf$from, 'y']
  edgedf$y_end = nodedf[edgedf$to, 'y']

  p1 = ggplot() +
    ggplot2::geom_segment(
      aes(
        x = x_begin,
        y = y_begin,
        xend = x_end,
        yend = y_end
      ),
      lwd = 0.2 * edgeSF,
      alpha = 1 / log10(nrow(edgedf)),
      colour = '#66666666',
      data = edgedf
    ) +
    #to cover up edges
    ggplot2::geom_point(
      aes(x, y, size = Size),
      fill = 'white',
      colour = 'white',
      shape = 21,
      stroke = 0.2 * nodeSF,
      data = nodedf
    ) +
    ggplot2::scale_size_continuous(range = c(0.1, 6) * nodeSF) +
    guides(size = guide_legend(ncol = 2)) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = 'bottom',
      plot.title = element_text(hjust = 0.5, size = rel(1.5))
    )

  if (is.null(enrichStat)) {
    p1 = p1 +
      #plot nodes
      ggplot2::geom_point(
        aes(x, y, fill = Category, size = Size),
        alpha = 0.75,
        shape = 21,
        stroke = 0.2 * nodeSF,
        data = nodedf
      ) +
      ggplot2::scale_fill_manual(values = colmap_nodes) +
      guides(fill = guide_legend(ncol = 4, override.aes = list(size = 5)))
  } else {
    #add stats
    commongs = intersect(nodedf$name, names(enrichStat))
    nodedf$EnrichStat = NA
    nodedf[commongs, 'EnrichStat'] = enrichStat[commongs]

    p1 = p1 +
      #plot nodes
      ggplot2::geom_point(
        aes(x, y, fill = EnrichStat, size = Size),
        alpha = 0.75,
        shape = 21,
        stroke = 0.2 * nodeSF,
        data = nodedf
      ) +
      scico::scale_fill_scico(palette = 'cork', na.value = '#AAAAAA')
  }

  #mark groups
  if (!is.null(markGroups)) {
    #add gene set counts for each group
    names(markGroups) = sapply(names(markGroups), function(x) {
      paste0(x, ' (n = ', length(markGroups[[x]]), ')')
    })

    #complex hull for groups
    hulldf = plyr::ldply(markGroups, function(x) {
      df = nodedf[x, ]
      df = df[grDevices::chull(df$x, df$y), ]
    }, .id = 'NodeGroup')

    p1 = p1 +
      ggnewscale::new_scale_colour() +
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

#' Custom theme
#'
#' @param rl a numeric, scaling factor to apply to text sizes
#'
#' @return a ggplot2 theme
#' @export
#'
#' @examples
#'
current_theme = function (rl = 1.1) {
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
