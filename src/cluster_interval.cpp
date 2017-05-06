#include <Rcpp.h>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace std;


// Function copied from http://stackoverflow.com/a/12399290/604854
// Modified to handle the RTYPE

template <int RTYPE>
vector<size_t> sort_indexes(const Vector<RTYPE>& v) {

  // initialize original index locations
  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

//' Cluster ranges which are implemented as 2 equal-length numeric vectors.
//' @param starts A numeric vector that defines the starts of each interval
//' @param ends A numeric vector that defines the ends of each interval
//' @examples
//' starts <- c(50, 100, 120)
//' ends <- c(75, 130, 150)
//' j <- cluster_interval(starts, ends)
//' j == c(0,1,1)
//' @export
// [[Rcpp::export]]
IntegerVector cluster_interval(NumericVector starts, NumericVector ends, int max_distance=0) {

  // Require that starts and ends are the same length

  // The implementation is inspired by the bedtools implementation:
  // https://github.com/arq5x/bedtools2/blob/14fbbb8aed5c6a04685da2cee3f11b98d70304a7/src/clusterBed/clusterBed.cpp
  IntegerVector result(starts.size());
  int cluster_id = 0;
  int prev_end = starts[0];
  for (auto i: sort_indexes(clone(starts))) {
    if(starts[i] - prev_end > max_distance){
      cluster_id++;
      prev_end = ends[i];
    }else{
      if(ends[i] > prev_end){
        prev_end = ends[i];
      }
    }
    result[i] = cluster_id;
  }

  return result;
}



