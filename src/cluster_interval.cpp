#include <Rcpp.h>


using namespace Rcpp;
using namespace std;


// The following code was copied from the stan math library
// https://github.com/stan-dev/stan/blob/e118db2b78ed33c40f7b5c774f3ce5b85aa5dfdf/src/stan/math/matrix/sort_indices.hpp

/*
 * Copyright (c) 2011--2015, Stan Developers and their Assignees
 All rights reserved.

 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 * Neither the name of Columbia University nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

template <bool ascending, typename C>
class index_comparator {
  const C& xs_;
public:
  /**
  * Construct an index comparator holding a reference
  * to the specified container.
  *
  * @patam xs Container
  */
  index_comparator(const C& xs) : xs_(xs) { }

  /**
  * Return true if the value at the first index is sorted in
  * front of the value at the second index;  this will depend
  * on the template parameter <code>ascending</code>.
  *
  * @param i Index of first value for comparison
  * @param j Index of second value for comparison
  */
  bool operator()(int i, int j) const {
    if (ascending)
      return xs_[i] < xs_[j];
    else
      return xs_[i] > xs_[j];
  }
};


/**
 * Return an integer array of indices of the specified container
 * sorting the values in ascending or descending order based on
 * the value of the first template prameter.
 *
 * @tparam ascending true if sort is in ascending order
 * @tparam C type of container
 * @param xs Container to sort
 * @return sorted version of container
 */
template <bool ascending, typename C>
std::vector<int> sort_indices(const C& xs) {
  typename C::size_type size = xs.size();
  std::vector<int> idxs;
  idxs.resize(size);
  for (typename C::size_type i = 0; i < size; ++i)
    idxs[i] = i;
  index_comparator<ascending,C> comparator(xs);
  std::sort(idxs.begin(), idxs.end(), comparator);
  return idxs;
}


// [[Rcpp::export]]
IntegerVector sort_indices(NumericVector x){
  return wrap(sort_indices<true>(as<std::vector<double> >(x)));
}


//' Cluster ranges which are implemented as 2 equal-length numeric vectors.
//' @param starts A numeric vector that defines the starts of each interval
//' @param ends A numeric vector that defines the ends of each interval
//' @param max_distance The maximum distance up to which intervals are still considered to be
//'  the same cluster. Default: 0.
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
  int cluster_id = -1;
  int prev_end = std::numeric_limits<int>::min();
  vector<int> indices = sort_indices<true>(as<std::vector<double> >(starts));
  for (int j = 0; j < indices.size(); j++) {
    int i = indices[j];
    Rcpp::checkUserInterrupt();
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



