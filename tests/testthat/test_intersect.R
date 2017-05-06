
context("genome_intersect")

suppressPackageStartupMessages(library(dplyr))

x1 <- data_frame(id = 1:4, bla=letters[1:4],
                 chromosome = c("chr1", "chr1", "chr2", "chr2"),
                 start = c(100, 200, 300, 400),
                 end = c(150, 250, 350, 450))

x2 <- data_frame(id = 1:4, BLA=LETTERS[1:4],
                 chromosome = c("chr1", "chr2", "chr2", "chr1"),
                 start = c(140, 210, 400, 300),
                 end = c(160, 240, 415, 320))

test_that("Intersection (both) of 2 data frames works as expected", {
  j <- genome_intersect(x1, x2, by=c("chromosome", "start", "end"), mode="both")
  # print(j)
  expect_equal(colnames(j), c("id.x", "bla", "chromosome", "id.y", "BLA", "start", "end"))
  expect_equal(j$start, c(140, 400))
  expect_equal(j$end, c(150, 415))
})

test_that("Intersection of 2 data frames works for multi-overlap ranges", {
  x2 <- data_frame(id = 1, BLA=LETTERS[1],
                   chromosome = c("chr1"),
                   start = c(140),
                   end = c(220))
  j <- genome_intersect(x1, x2, by=c("chromosome", "start", "end"), mode="both")
  # print(j)
  expect_equal(colnames(j), c("id.x", "bla", "chromosome", "id.y", "BLA", "start", "end"))
  expect_equal(j$start, c(140, 200))
  expect_equal(j$end, c(150, 220))
  expect_equal(j$id.x, c(1,2))
  expect_equal(j$id.y, c(1,1))

})



test_that("Intersection of 2 data frames works for multi-overlap ranges the other way around", {
  x1 <- data_frame(id = 1, bla=letters[1],
                   chromosome = c("chr1"),
                   start = c(100),
                   end = c(420))
  j <- genome_intersect(x1, x2, by=c("chromosome", "start", "end"), mode="both")
  # print(j)
  expect_equal(colnames(j), c("id.x", "bla", "chromosome", "id.y", "BLA", "start", "end"))
  expect_equal(j$start, c(140, 300))
  expect_equal(j$end, c(160, 320))
  expect_equal(j$id.x, c(1,1))
  expect_equal(j$id.y, c(1,4))

})





