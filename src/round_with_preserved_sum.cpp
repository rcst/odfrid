#include "Rcpp.h"
#include "round_with_preserved_sum.h"

using namespace Rcpp;

// Custom comparator for sorting by difference
bool compareByDifference(const TempElement &a, const TempElement &b) {
    return a.difference > b.difference;
}

// Custom comparator for sorting by original index
bool compareByIndex(const TempElement &a, const TempElement &b) {
    return a.index < b.index;
}

// [[Rcpp::export]]
std::vector<int> roundWithPreservedSum(const std::vector<double>& fn) {
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
    std::vector<int> roundedResult(n);
    for (int i = 0; i < n; ++i) {
        roundedResult[i] = tempArr[i].result;
    }

    return roundedResult;
}
