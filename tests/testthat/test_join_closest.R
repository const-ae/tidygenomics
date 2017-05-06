
context("genome_join_closest")

library(dplyr)

x1 <- data_frame(id = 1:4, bla=letters[1:4],
                 chromosome = c("chr1", "chr1", "chr2", "chr3"),
                 start = c(100, 200, 300, 400),
                 end = c(150, 250, 350, 450))

x2 <- data_frame(id = 1:4, BLA=LETTERS[1:4],
                 chromosome = c("chr1", "chr1", "chr1", "chr2"),
                 start = c(220, 210, 300, 400),
                 end = c(225, 240, 320, 415))

test_that("Joining with closest works as expected", {
  j <- genome_join_closest(x1, x2, by=c("chromosome", "start", "end"), distance_column_name="distance", mode="left")
  print(j)
  expect_equal(colnames(j), c("id.x", "bla", "chromosome.x", "start.x", "end.x",
                              "id.y", "BLA", "chromosome.y", "start.y", "end.y", "distance"))
  expect_equal(j$start.y, c(210, 220, 210, 400, NA))
  expect_equal(j$distance, c(59, 0, 0, 49, NA))
})


