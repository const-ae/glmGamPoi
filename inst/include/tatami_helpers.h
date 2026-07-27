#ifndef TATAMI_HELP_H
#define TATAMI_HELP_H

#include "tatami/tatami.hpp"

// Oracle sub-class which just always fetches a compile-time constant index
template <tatami::NumericMatrix::index_type IndexValue> class ConstIndexOracle final : public tatami::Oracle<tatami::NumericMatrix::index_type> {
private:
  tatami::PredictionIndex length_;

public:
  ConstIndexOracle(const tatami::NumericMatrix::index_type length) : length_(sanisizer::cast<tatami::PredictionIndex>(length)) {}
  tatami::PredictionIndex total() const { return length_; }
  tatami::NumericMatrix::index_type get(const tatami::PredictionIndex i) const { return IndexValue; }
};

#endif