

#' Join intervals on chromosomes in data frames, to the closest partner
#'
#' @param x A dataframe.
#' @param y A dataframe.
#' @param by A character vector with 3 entries which are used to match the chromosome, start and end column.
#'   For example: \code{by=c("Chromosome"="chr", "Start"="start", "End"="end")}
#' @param mode One of "inner", "full", "left", "right", "semi" or "anti".
#' @param distance_column_name A string that is used as the new column name with the distance.
#' If \code{NULL} no new column is added.
#' @param max_distance The maximum distance that is allowed to join 2 entries.
#' @param select A string that is passed on to \code{IRanges::distanceToNearest}, can either be
#'  all which means that in case that multiple intervals have the same distance all are reported, or
#'  arbitrary which means in that case one would be chosen at random.
#' @param ... Additional arguments parsed on to genome_join_closest.
#' @return The joined dataframe of \code{x} and \code{y}.
#' @examples
#'
#' library(dplyr)
#'
#' x1 <- data.frame(id = 1:4, bla=letters[1:4],
#'                  chromosome = c("chr1", "chr1", "chr2", "chr2"),
#'                  start = c(100, 200, 300, 400),
#'                  end = c(150, 250, 350, 450))
#'
#' x2 <- data.frame(id = 1:4, BLA=LETTERS[1:4],
#'                  chromosome = c("chr1", "chr2", "chr2", "chr1"),
#'                  start = c(140, 210, 400, 300),
#'                  end = c(160, 240, 415, 320))
#' j <- genome_intersect(x1, x2, by=c("chromosome", "start", "end"), mode="both")
#' print(j)
#' @export
genome_join_closest <- function(x, y, by=NULL,  mode = "inner",
                                distance_column_name=NULL, max_distance=Inf, select="all"){

  # Nearly all of this code is copied from https://github.com/dgrtwo/fuzzyjoin

  if (!requireNamespace("IRanges", quietly = TRUE)) {
    stop("genome_join_closest requires the IRanges package: ",
         "https://bioconductor.org/packages/release/bioc/html/IRanges.html")
  }

  select <- match.arg(select, c("all", "arbitrary"))

  by <- dplyr::common_by(by, x, y)
  if (length(by$x) != 3) {
    stop("genome_join_closest must join on exactly three columns")
  }

  f <- function(x, y) {
    # nest around the chromosome column
    x$..index <- seq_len(nrow(x))
    y$..index <- seq_len(nrow(y))

    nested_x <- dplyr::group_by_at(x, 1) %>% tidyr::nest()
    nested_y <- dplyr::group_by_at(y, 1) %>% tidyr::nest()
    by <- c(colnames(nested_y)[1])
    names(by) <- colnames(nested_x)[1]

    joined <- dplyr::inner_join(nested_x, nested_y, by = by)

    # find matching ranges in each
    find_closest <- function(xd, yd) {
      r1 <- IRanges::IRanges(xd[[1]], xd[[2]])
      r2 <- IRanges::IRanges(yd[[1]], yd[[2]])
      o <- as.data.frame(IRanges::distanceToNearest(r1, r2, select=select))

      data.frame(x = xd$..index[o$queryHits], y = yd$..index[o$subjectHits], ..distance=o$distance) %>%
        dplyr::filter(`..distance` < max_distance)
    }

    ret <- purrr::map2_df(joined$data.x, joined$data.y, find_closest)

    if(! is.null(distance_column_name)){
      ret[[distance_column_name]] <- ret$..distance
    }
    ret$..distance <- NULL

    ret
  }

  fuzzyjoin::fuzzy_join(x, y, mode = mode, index_match_fun = f, multi_by = by)

}


#' @rdname genome_join_closest
#' @export
genome_inner_join_closest <- function(x, y, by = NULL, ...) {
  genome_join_closest (x, y,  by, mode = "inner", ...)
}


#' @rdname genome_join_closest
#' @export
genome_left_join_closest <- function(x, y, by = NULL, ...) {
  genome_join_closest (x, y,  by, mode = "left", ...)
}


#' @rdname genome_join_closest
#' @export
genome_right_join_closest <- function(x, y, by = NULL, ...) {
  genome_join_closest (x, y,  by, mode = "right", ...)
}


#' @rdname genome_join_closest
#' @export
genome_full_join_closest <- function(x, y, by = NULL, ...) {
  genome_join_closest (x, y,  by, mode = "full", ...)
}


#' @rdname genome_join_closest
#' @export
genome_semi_join_closest <- function(x, y, by = NULL, ...) {
  genome_join_closest (x, y,  by, mode = "semi", ...)
}


#' @rdname genome_join_closest
#' @export
genome_anti_join_closest <- function(x, y, by = NULL, ...) {
  genome_join_closest (x, y,  by, mode = "anti", ...)
}


