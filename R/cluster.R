
#' @useDynLib tidygenomics
#' @importFrom Rcpp sourceCpp
NULL

#' Intersect data frames based on chromosome, start and end.
#'
#' @param x A dataframe.
#' @param by A character vector with 3 entries which are the chromosome, start and end column.
#'   For example: \code{by=c("chr", "start", "end")}
#' @return The dataframe with the additional column of the cluster
#' @examples
#' x1 <- data_frame(id = 1:4, bla=letters[1:4],
#'                  chromosome = c("chr1", "chr1", "chr2", "chr1"),
#'                  start = c(100, 120, 300, 260),
#'                  end = c(150, 250, 350, 450))
#' genome_cluster(x1, by=c("chromosome", "start", "end"))
#' genome_cluster(x1, by=c("chromosome", "start", "end"), max_distance=10)
#' @export
genome_cluster <- function(x, by=NULL, max_distance=0, cluster_column_name="cluster_id"){

  x_groups <- dplyr::groups(x)
  x <- dplyr::ungroup(x)
  regroup <- function(d) {
    if (is.null(x_groups)) {
      return(d)
    }
    g <- purrr::map_chr(x_groups, as.character)
    missing <- !(g %in% colnames(d))
    g[missing] <- paste0(g[missing], ".x")
    dplyr::group_by_(d, .dots = g)
  }

  if (is.null(by) | length(by) != 3) {
    stop("genome_cluster must join on exactly three columns")
  }


  cluster_counter <- 0
  ret <- x %>%
    dplyr::group_by_(by[1]) %>%
    # mutate(cluster_column_name=cluster_interval(start, end))
    dplyr::mutate_(.dots=setNames(list(paste0("tidygenomics::cluster_interval(", by[2], ",", by[3], ")")), cluster_column_name)) %>%
    dplyr::ungroup() %>%
    # dplyr::mutate(cluster=(paste0(chromosome,"-", cluster) %>% as.factor %>% as.numeric())-1)
    dplyr::mutate_(.dots=setNames(list(paste0("as.numeric(as.factor(paste0(", by[1],",\"-\",", cluster_column_name, ")))-1")),
                                  cluster_column_name))


  ret <-regroup(ret)
  return(ret)

}
