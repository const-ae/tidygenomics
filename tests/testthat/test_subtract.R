
context("genome_subtract")


suppressPackageStartupMessages(library(dplyr))

x1 <- tibble(id = 1:4, bla=letters[1:4],
                 chromosome = c("chr1", "chr1", "chr2", "chr1"),
                 start = c(100, 200, 300, 400),
                 end = c(150, 250, 350, 450))

x2 <- tibble(id = 1:4, BLA=LETTERS[1:4],
                 chromosome = c("chr1", "chr2", "chr1", "chr1"),
                 start = c(120, 210, 300, 400),
                 end = c(125, 240, 320, 415))

test_that("Subtraction of 2 data frames works as expected", {
  j <- genome_subtract(x1, x2, by=c("chromosome", "start", "end"))
  # print(j)
  expect_equal(colnames(j), c("id", "bla", "chromosome", "start", "end"))
  expect_equal(j$start, c(100, 126, 200, 300, 416))
  expect_equal(j$end,   c(119, 150, 250, 350, 450))
})



test_that("Edge cases of subtraction of 2 data frames works as expected", {
  x1 <- tibble(id = 1:2, bla=letters[1:2],
                   chromosome = c("chr1", "chr1"),
                   start = c(100, 200),
                   end = c(150, 250))

  x2 <- tibble(id = 1:4, BLA=LETTERS[1:4],
                   chromosome = c("chr1", "chr1", "chr1", "chr1"),
                   start = c(120, 110, 190, 400),
                   end = c(125, 122, 320, 415))

  j <- genome_subtract(x1, x2, by=c("chromosome", "start", "end"))
  print(j)
  expect_equal(colnames(j), c("id", "bla", "chromosome", "start", "end"))
  expect_equal(j$start, c(100, 126))
  expect_equal(j$end,   c(109, 150))
})


test_that("during subtraction the intervals are not unified", {
  x1 <- tibble(id = 1:3, bla=letters[1:3],
                   chromosome = c("chr1", "chr1", "chr1"),
                   start = c(100, 115, 200),
                   end = c(150, 160, 250))

  x2 <- tibble(id = 1, BLA=LETTERS[1],
                   chromosome = c("chr1"),
                   start = c(110),
                   end = c(130))

  j <- genome_subtract(x1, x2, by=c("chromosome", "start", "end"))
  print(j)
  expect_equal(colnames(j), c("id", "bla", "chromosome", "start", "end"))
  expect_equal(j$start, c(100, 131, 131, 200))
  expect_equal(j$end,   c(109, 150, 160, 250))
})
