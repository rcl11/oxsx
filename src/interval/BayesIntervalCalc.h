#ifndef __OXSX_BAYESINTERVALCALC__
#define __OXSX_BAYESINTERVALCALC__
#include <Histogram.h>
#include <PDF.h>

class BayesIntervalCalc{
 public:
    static double UpperBound(Histogram lh_, double cl_); // needs to be a 1D histogram
    static double UpperBound(double expectedCounts_, int observedCounts_, double cl_, double tolerance_, double startValue_ = 0);
    static double UpperBound(PDF& b_prob_, double expectedCounts_, double expectedCountsSigma_, int observedCounts_, double cl_);
    static double OneSidedUpperInterval(double expectedCounts_, int observedCounts_, double upperEdge_);
 private:   
    static double Marginalise(PDF& b_prob, double s, double b, double width, int n);
    static double CumulativePosterior(PDF& b_prob, double signal1, double signal2, double b, double width, int n);
    static double PoissonProb(PDF& b_prob, double s, double b, double width, int n);
};
#endif
