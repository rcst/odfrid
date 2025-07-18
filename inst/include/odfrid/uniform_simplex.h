#ifndef UNIFORM_SIMPLEX_H
#define UNIFORM_SIMPLEX_H

struct TempElement {
    int result;        // floor value
    double difference; // fractional part
    int index;         // original index
};

bool compareByDifference(const TempElement &a, const TempElement &b);
bool compareByIndex(const TempElement &a, const TempElement &b);
std::vector<unsigned int> roundWithPreservedSum(const std::vector<double>& fn);
std::vector<unsigned int> adjustWithCaps(const std::vector<unsigned int>& x, const std::vector<unsigned int>& z);
std::vector<double> uniformSimplexSample(int N, double C = 1.0, double A = 0.0);

#endif
