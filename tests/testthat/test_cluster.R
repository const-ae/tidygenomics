
context("genome_cluster")

library(dplyr)

x1 <- tibble(id = 1:4, bla=letters[1:4],
                 chromosome = c("chr1", "chr1", "chr2", "chr1"),
                 start = c(100, 120, 300, 260),
                 end = c(150, 250, 350, 450))


test_that("genome_clustering assings that correct clusters", {
  j <- genome_cluster(x1, by=c("chromosome", "start", "end"), max_distance=5)

  print(j)

  expect_equal(j$cluster_id, c(0,0,2,1))
})


test_that("cluster_interval works", {
  starts <- c(50, 100, 120)
  ends <- c(75, 130, 150)
  j <- cluster_interval(starts, ends)
  expect_equal(j, c(0,1,1))
  expect_equal(cluster_interval(starts, ends, max_distance = 24), c(0,1,1))
  expect_equal(cluster_interval(starts, ends, max_distance = 25), c(0,0,0))

  starts <- c(50, 100, 120, 180, 350)
  ends <- c(75, 200, 150, 210, 400)
  expect_equal(cluster_interval(starts, ends), c(0,1,1,1,2))

  starts <- c(500, 300, 150)
  ends <- c(510, 310, 160)
  expect_equal(cluster_interval(starts, ends), c(2,1,0))

  expect_equal(cluster_interval(numeric(0), numeric(0)), numeric(0))
})
