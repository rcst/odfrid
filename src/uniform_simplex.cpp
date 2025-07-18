#include "Rcpp.h"
#include <random>
#include <algorithm>
#include "odfrid/uniform_simplex.h"

using namespace Rcpp;

// Generate N values uniformly distributed on the simplex (sum = 1),
// then optionally scale and shift each by C and A respectively.
std::vector<double> uniformSimplexSample(int N, double C, double A) {
  // std::random_device rd;
  // std::mt19937 gen(rd()); // Mersenne Twister RNG
  // std::uniform_real_distribution<> uniform(0.0, 1.0);

  std::vector<double> u(N);
  for (int i = 0; i < N; ++i) {
    double v = (double)runif(1, 0, 1)[0];
    u[i] = -std::log(v); // Exponential(1) sample
  }

  double sumU = std::accumulate(u.begin(), u.end(), 0.0);

  std::vector<double> p(N);
  for (int i = 0; i < N; ++i) {
    p[i] = C * (u[i] / sumU) + A;
  }

  return p;
}

std::vector<unsigned int> adjustWithCaps(const std::vector<unsigned int>& x, const std::vector<unsigned int>& z) {
  if (x.size() != z.size()) {
    stop("Vectors x and z must have the same length.");
  }

  int totalSumX = std::accumulate(x.begin(), x.end(), 0);
  int totalCapacity = std::accumulate(z.begin(), z.end(), 0);

  if (totalCapacity < totalSumX) {
    stop("The total capacity defined by z is smaller than the total sum of x. Redistribution is impossible.");
  }

  int n = x.size();
  std::vector<unsigned int> result = x;
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

// Custom comparator for sorting by difference
bool compareByDifference(const TempElement &a, const TempElement &b) {
    return a.difference > b.difference;
}

// Custom comparator for sorting by original index
bool compareByIndex(const TempElement &a, const TempElement &b) {
    return a.index < b.index;
}

std::vector<unsigned int> roundWithPreservedSum(const std::vector<double>& fn) {
    int n = fn.size();
    std::vector<TempElement> tempArr(n);
    double arraySum = 0;

    // Calculate expected sum
    for (double val : fn)
        arraySum += val;

    int lowerSum = 0;

    // Populate tempArr
    for (int i = 0; i < n; ++i) {
        int floored = std::floor(fn[i]);
        tempArr[i] = {
            floored,
            fn[i] - floored,
            i
        };
        lowerSum += floored;
    }

    // Sort by difference descending
    std::sort(tempArr.begin(), tempArr.end(), compareByDifference);

    int difference = static_cast<int>(std::round(arraySum - lowerSum));

    // Adjust values with largest remainders
    for (int i = 0; i < difference; ++i) {
        tempArr[i].result += 1;
    }

    // Sort back to original index order
    std::sort(tempArr.begin(), tempArr.end(), compareByIndex);

    // Extract result values
    std::vector<unsigned int> roundedResult(n);
    for (int i = 0; i < n; ++i) {
        roundedResult[i] = tempArr[i].result;
    }

    return roundedResult;
}
