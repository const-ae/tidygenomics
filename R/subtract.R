


#' Subtract one data frame from another based on chromosome, start and end.
#'
#' @param x A dataframe.
#' @param y A dataframe.
#' @param by A character vector with 3 entries which are used to match the chromosome, start and end column.
#'   For example: \code{by=c("Chromosome"="chr", "Start"="start", "End"="end")}
#' @return The subtracted dataframe of \code{x} and \code{y} with the new boundaries.
#' @examples
#'
#' library(dplyr)
#'
#' x1 <- data.frame(id = 1:4, bla=letters[1:4],
#'                  chromosome = c("chr1", "chr1", "chr2", "chr1"),
#'                  start = c(100, 200, 300, 400),
#'                  end = c(150, 250, 350, 450))
#'
#' x2 <- data.frame(id = 1:4, BLA=LETTERS[1:4],
#'                  chromosome = c("chr1", "chr2", "chr1", "chr1"),
#'                  start = c(120, 210, 300, 400),
#'                  end = c(125, 240, 320, 415))
#'
#' j <- genome_subtract(x1, x2, by=c("chromosome", "start", "end"))
#' print(j)
#'
#'
#' @export
genome_subtract <- function(x, y, by=NULL){

  # Much of this code is copied from https://github.com/dgrtwo/fuzzyjoin

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

  by <- dplyr::common_by(by, x, y)

  if (length(by$x) != 3) {
    stop("genome_join must join on exactly three columns")
  }


  f <- function(x,y){
    # nest around the chromosome column
    x$..index <- seq_len(nrow(x))
    y$..index <- seq_len(nrow(y))
    nested_x <- tidyr::nest_(x, "x_data", colnames(x)[-1])
    nested_y <- tidyr::nest_(y, "y_data", colnames(y)[-1])
    by <- c(colnames(nested_y)[1])
    names(by) <- colnames(nested_x)[1]

    joined <- dplyr::inner_join(nested_x, nested_y, by = by)

    # find matching ranges in each
    find_subtractions <- function(xd, yd) {
      r1 <- IRanges::IRanges(xd[[1]], xd[[2]])
      r2 <- IRanges::IRanges(yd[[1]], yd[[2]])

      subtraction <- IRanges::setdiff(r1, r2)

      o <- as.data.frame(IRanges::findOverlaps(subtraction, r1))
      data.frame(x = xd$..index[o$subjectHits],
                 ..start=pmax(IRanges::start(subtraction)[o$queryHits], IRanges::start(r1)[o$subjectHits]),
                 ..end=pmin(IRanges::end(subtraction)[o$queryHits], IRanges::end(r1)[o$subjectHits]))
    }

    ret <- purrr::map2_df(joined$x_data, joined$y_data, find_subtractions)
    ret
  }

  d1 <- x[, by$x, drop = FALSE]
  d2 <- y[, by$y, drop = FALSE]

  matches <- f(d1, d2)
  ret <- x %>%
    dplyr::select(- dplyr::one_of(by$x[-1])) %>%
    dplyr::mutate(..id=seq_len(n())) %>%
    dplyr::inner_join(matches[, c("x", "..start", "..end")], by=c("..id"="x")) %>%
    dplyr::rename(!! by$x[2] := `..start`, !! by$x[3] := `..end`) %>%
    dplyr::select(- `..id`) %>%
    regroup()
  return(ret)

}
