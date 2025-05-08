#ifndef ROUND_WITH_PRESERVED_SUM_H
#define ROUND_WITH_PRESERVED_SUM_H

struct TempElement {
    int result;        // floor value
    double difference; // fractional part
    int index;         // original index
};

bool compareByDifference(const TempElement &a, const TempElement &b);
bool compareByIndex(const TempElement &a, const TempElement &b);
std::vector<int> roundWithPreservedSum(const std::vector<double>& fn);

#endif
