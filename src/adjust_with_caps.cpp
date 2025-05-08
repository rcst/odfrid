// #include <iostream>
// #include <vector>
#include <random>
// #include <numeric>
#include <algorithm>

#include "Rcpp.h"

using namespace Rcpp;

// [[Rcpp::export]]
std::vector<int> adjustWithCaps(const std::vector<int>& x, const std::vector<int>& z) {
  if (x.size() != z.size()) {
    stop("Vectors x and z must have the same length.");
  }

  int totalSumX = std::accumulate(x.begin(), x.end(), 0);
  int totalCapacity = std::accumulate(z.begin(), z.end(), 0);

  if (totalCapacity < totalSumX) {
    stop("The total capacity defined by z is smaller than the total sum of x. Redistribution is impossible.");
  }

  int n = x.size();
  std::vector<int> result = x;
  int cappedSum = 0;

  // First, cap all values at z_i if they exceed it
  for (int i = 0; i < n; ++i) {
    if (result[i] > z[i]) {
      result[i] = z[i];
    }
    cappedSum += result[i];
  }

  int excess = totalSumX - cappedSum;
  if (excess == 0) return result; // Nothing to redistribute

  // Create a list of indices eligible for receiving excess
  std::vector<int> eligibleIndices;
  for (int i = 0; i < n; ++i) {
    if (result[i] < z[i]) {
      eligibleIndices.push_back(i);
    }
  }

  std::random_device rd;
  std::mt19937 gen(rd());

  // Redistribute the excess randomly within capacity
  while (excess > 0) {
    std::shuffle(eligibleIndices.begin(), eligibleIndices.end(), gen);
    for (int idx : eligibleIndices) {
      if (result[idx] < z[idx]) {
        result[idx]++;
        excess--;
        if (excess == 0) break;
      }
    }
  }

  return result;
}
