#include <BayesIntervalCalc.h>
#include <PdfExceptions.h>
#include <Histogram.h>
#include <iostream>

double
BayesIntervalCalc::UpperBound(Histogram posterior_, double cl_){
    // Makes no sense for an empty histogram
    if(!posterior_.Integral())
        throw 0;
    
    if(cl_ <= 0)
        throw 0;
    
    // only works for 1D histograms currently
    if(posterior_.GetNDims() != 1)
        throw DimensionError("Unable to produce upper limit on posteriors dim != 1. Marginalise?");

    posterior_.Normalise();

    // Integrate from the minimum to x until the total probability exceeds cl_
    // find the first bin <b> for which the integral  <0>-><b> > cl
    int    critBin  = 0;
    double integral = 0;
    double sum = 0;
    while(true){
        sum = 0;              
        for(size_t i = 0; i < critBin; i++){
            sum += posterior_.GetBinContent(i);
        }
        if (sum < cl_)
            critBin++;
        else
            break;
        
    }

    // now interpolate between the bins    
    double upperEdge = posterior_.GetAxes().GetBinHighEdge(critBin, 0);
    double lowerEdge = posterior_.GetAxes().GetBinLowEdge(critBin, 0);
    double content   = posterior_.GetBinContent(critBin);

    return upperEdge - (upperEdge - lowerEdge) * (sum - cl_)/content;
         
}
