
context("genome_complement")

library(dplyr)

x1 <- data_frame(id = 1:4, bla=letters[1:4],
                 chromosome = c("chr1", "chr1", "chr2", "chr1"),
                 start = c(100, 200, 300, 400),
                 end = c(150, 250, 350, 450))

test_that("Calculating the complement of a sequence works", {
  j <- genome_complement(x1, by=c("chromosome", "start", "end"))
  print(j)
  expect_equal(j$chromosome, c("chr1", "chr1", "chr1", "chr2"))
  expect_equal(j$start, c(1,151, 251,1))
  expect_equal(j$end, c(99,199, 399, 299))
})
