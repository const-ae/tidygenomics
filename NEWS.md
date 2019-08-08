
# v0.1.2

* Remove remaining calls to deprecated tidyr functions to become compatible with tidyr v1.0.0.
  Thanks to @jennybc for the pull request (#6)

# v0.1.1

* Fix issue #5
  - the genome_cluster method assigned all chunks to cluster zero if their end was smaller 
    than the end of the first entry

* Port dplyr calls to new tidyeval API
  - This avoids plenty of deprecation warnings

* Add pkgdown webpage: https://const-ae.github.io/tidygenomics/

# Initial Release (v0.1.0)

First acceptance on CRAN
