



#' Calculates the complement to the intervals covered by the intervals in
#' a data frame. It can optionally take a \code{chromosome_size} data frame
#' that contains 2 or 3 columns, the first the names of chromosome and in case
#' there are 2 columns the size or first the start index and lastly the end index
#' on the chromosome.
#' @param x A data frame for which the complement is calculated
#' @param chromosome_size A dataframe with at least 2 columns that contains
#'  first the chromosome name and then the size of that chromosome. Can be NULL
#'  in which case the largest value per chromosome from \code{x} is used.
#' @param by A character vector with 3 entries which are the chromosome, start and end column.
#'   For example: \code{by=c("chr", "start", "end")}
#' @examples
#'
#' library(dplyr)
#'
#' x1 <- data.frame(id = 1:4, bla=letters[1:4],
#'                  chromosome = c("chr1", "chr1", "chr2", "chr1"),
#'                  start = c(100, 200, 300, 400),
#'                  end = c(150, 250, 350, 450))
#'
#' genome_complement(x1, by=c("chromosome", "start", "end"))
#' @export
genome_complement <- function(x, chromosome_size=NULL, by=NULL){

  if (is.null(by) | length(by) != 3) {
    stop("genome_complement must work on exactly three columns")
  }


  if(is.null(chromosome_size)){
    chromosome_size <- x %>%
      dplyr::group_by(!! sym(by[1])) %>%
      dplyr::summarize(start = 1,
                       end = max(!! sym(by[3])))
  }else if(ncol(chromosome_size) == 2){
    chromosome_size <- cbind(chromosome_size[, 1, drop=FALSE], data.frame(start=1), chromosome_size[, -1, drop=FALSE])
  }

  colnames(chromosome_size)[1:3] <- by

  chromosome_size %>%
    genome_subtract(x, by=by)
}
