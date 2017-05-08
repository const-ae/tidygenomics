
#' Intersect data frames based on chromosome, start and end.
#'
#' @param x A dataframe.
#' @param y A dataframe.
#' @param by A character vector with 3 entries which are used to match the chromosome, start and end column.
#'   For example: \code{by=c("Chromosome"="chr", "Start"="start", "End"="end")}
#' @param mode One of "both", "left", "right" or "anti".
#' @return The intersected dataframe of \code{x} and \code{y} with the new boundaries.
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
#'
#'
#' @importFrom dplyr "%>%"
#'
#' @export
genome_intersect <- function(x, y, by=NULL, mode= "both"){

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

  mode <- match.arg(mode, c("both", "left", "right", "anti"))

  by <- dplyr::common_by(by, x, y)

  if (length(by$x) != 3) {
    stop("genome_join must join on exactly three columns")
  }


  index_match_fun <- function(x,y){
    # nest around the chromosome column
    x$..index <- seq_len(nrow(x))
    y$..index <- seq_len(nrow(y))
    nested_x <- tidyr::nest_(x, "x_data", colnames(x)[-1])
    nested_y <- tidyr::nest_(y, "y_data", colnames(y)[-1])
    by <- c(colnames(nested_y)[1])
    names(by) <- colnames(nested_x)[1]

    joined <- dplyr::inner_join(nested_x, nested_y, by = by)

    # find matching ranges in each
    find_overlaps <- function(xd, yd) {
      r1 <- IRanges::IRanges(xd[[1]], xd[[2]])
      r2 <- IRanges::IRanges(yd[[1]], yd[[2]])
      o <- as.data.frame(IRanges::findOverlaps(r1, r2))
      intersection <- IRanges::pintersect(r1[o$queryHits], r2[o$subjectHits])
      data.frame(x = xd$..index[o$queryHits], y = yd$..index[o$subjectHits],
                 ..start=IRanges::start(intersection), ..end=IRanges::end(intersection))
    }

    ret <- purrr::map2_df(joined$x_data, joined$y_data, find_overlaps)
    ret
  }

  d1 <- x[, by$x, drop = FALSE]
  d2 <- y[, by$y, drop = FALSE]
  matches <- index_match_fun(d1, d2)

  matches$i <- NULL
  if (mode == "anti") {
    if (nrow(matches) == 0) {
      return(regroup(x))
    }
    return(regroup(x[-sort(unique(matches$x)), ]))
  }
  if (mode == "left") {
    ret <- x %>%
      dplyr::select(- dplyr::one_of(by$x[-1])) %>%
      dplyr::mutate(..id=seq_len(n())) %>%
      dplyr::inner_join(matches[, c("x", "..start", "..end")], by=c("..id"="x")) %>%
      dplyr::rename_(.dots=stats::setNames(c("..start", "..end"), by$x[-1])) %>%
      dplyr::select_(quote(- `..id`)) %>%
      regroup()
    return(ret)
  }
  else if (mode == "right") {
    ret <- y %>%
      dplyr::select(- dplyr::one_of(by$y[-1])) %>%
      dplyr::mutate(..id=seq_len(n())) %>%
      dplyr::inner_join(matches[,c("y", "..start", "..end")], by=c("..id"="y")) %>%
      dplyr::rename_(.dots=stats::setNames(c("..start", "..end"), by$y[-1])) %>%
      dplyr::select_(quote(- `..id`)) %>%
      regroup()
    return(ret)
  }

  matches <- dplyr::arrange(matches, x, y)
  for (n in intersect(colnames(x), colnames(y))) {
    if(! n %in% by$x){
      x <- dplyr::rename_(x, .dots = structure(n, .Names = paste0(n,".x")))
    }
    if(! n %in% by$y){
      y <- dplyr::rename_(y, .dots = structure(n, .Names = paste0(n,".y")))
    }
  }

  ret <- dplyr::bind_cols(x[matches$x, , drop = FALSE] %>% dplyr::select(- dplyr::one_of(by$x[-1])),
                          y[matches$y, , drop = FALSE] %>% dplyr::select(- dplyr::one_of(by$y)))
  if (ncol(matches) > 2) {
    extra_cols <- matches[, -(1:2), drop = FALSE]
    ret <- dplyr::bind_cols(ret, extra_cols) %>%
      dplyr::rename_(.dots=stats::setNames(c("..start", "..end"), by$x[-1]))
  }
  regroup(ret)


}
