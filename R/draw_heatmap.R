#' Draw Heatmap for Microbial Sample Taxa Matrix
#'
#' This function creates a heatmap for a given microbial sample taxa matrix. 
#' It is specifically designed for visualizing abundance patterns across 
#' different samples and taxa in microbiome studies. The function applies 
#' an arcsinh transformation to the data for normalization and better visualization 
#' of abundance patterns, especially useful in handling highly skewed microbiome data.
#'
#' The heatmap is generated using `ggplot2` and `reshape2` packages, with taxa 
#' on the x-axis and samples on the y-axis. The color intensity in the heatmap 
#' represents the arcsinh-transformed abundance of each taxa in each sample.
#'
#' @param X A numeric matrix representing microbial sample taxa data, where rows 
#'   represent samples and columns represent taxa.
#' @param title An optional title for the heatmap
#' @return A ggplot object representing the heatmap. This can be further customized 
#'   or directly plotted.
#'
#' @examples
#' # Create a sample matrix representing microbial sample taxa data
#' mat <- matrix(runif(20), nrow = 5)
#' 
#' # Draw heatmap
#' draw_heatmap(mat)
#'
#' @export
#'
#' @importFrom ggplot2 ggplot geom_tile theme_bw labs theme element_text 
#' scale_fill_gradient element_blank
#' @importFrom reshape2 melt
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr filter
#' @importFrom dplyr select
#'
draw_heatmap <- function(X, title="") {
  reshape2::melt(asinh(X)) %>%
    dplyr::rename(sample = Var1, taxa = Var2, asinh.abun = value) %>%
    ggplot2::ggplot(., aes (x = taxa, y = sample, fill = asinh.abun)) +
    ggplot2::geom_tile() + ggplot2::theme_bw() + ggplot2::ggtitle(title) +
    ggplot2::labs(fill = "arcsinh\nabundance") +
    ggplot2::theme(plot.title = element_text(hjust = 0.5),
                   axis.title.x=element_blank(), axis.title.y=element_blank(),
                   axis.text.x = element_text(size=3, angle=90), axis.text.y = element_text(size=4)) +
    ggplot2::scale_fill_gradient(low="yellow", high="red", na.value = "grey") +
    theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),
          axis.title.y = element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
    theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
          panel.background = element_blank(),panel.border = element_blank())+
    theme(legend.position="none") +
    ggplot2::coord_fixed(ratio = 1)  # Fixing the aspect ratio
}
