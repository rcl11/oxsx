#include <BayesIntervalCalc.h>
#include <Exceptions.h>
#include <Histogram.h>
#include <TH1D.h>
#include <gsl/gsl_sf_gamma.h> // imcomplete gamma function
#include <iostream>
#include <cmath>

double
BayesIntervalCalc::UpperBound(Histogram posterior_, double cl_){
    // Makes no sense for an empty histogram
    if(!posterior_.Integral())
        throw ValueError("BayesIntervalCalc::Empty histogram passed!");
    
    if(cl_ <= 0)
        throw ValueError(Formatter() << "BayesIntervalCalc:: cl = " << cl_
                         << " , must be >0!");
    
    // only works for 1D histograms currently
    if(posterior_.GetNDims() != 1)
        throw DimensionError("BayesIntervalCal", 1, posterior_.GetNDims(), 
                             "Only implemented for 1D, marginalise?");

    posterior_.Normalise();

    // Integrate from the minimum to x until the total probability exceeds cl_
    // find the first bin <b> for which the integral  <0>-><b> > cl
    int    critBin  = 0;
    double integral = 0;
    double sum = 0;
    while(sum < cl_){
        sum = 0;              
        for(size_t i = 0; i < critBin; i++){
            sum += posterior_.GetBinContent(i);
        }

        critBin++;        
    }

    // now interpolate between the bins    
    double upperEdge = posterior_.GetAxes().GetBinHighEdge(critBin, 0);
    double lowerEdge = posterior_.GetAxes().GetBinLowEdge(critBin, 0);
    double content   = posterior_.GetBinContent(critBin);

    return upperEdge - (upperEdge - lowerEdge) * (sum - cl_)/content;
         
}

double
BayesIntervalCalc::OneSidedUpperInterval(double expectedCounts_, int observedCounts_, double upperEdge_){
    return (gsl_sf_gamma_inc_P(observedCounts_ + 1, expectedCounts_ + upperEdge_ ) - gsl_sf_gamma_inc_P(observedCounts_ +1, expectedCounts_))/ (1 - gsl_sf_gamma_inc_P(observedCounts_ + 1, expectedCounts_));
}

double 
BayesIntervalCalc::UpperBound(double expectedCounts_, int observedCounts_, double cl_, double tolerance_, double startValue_){
    double currentCL = 0;
    double limit = 0;
    while(std::abs(currentCL - cl_) > tolerance_){
        limit += tolerance_;
        currentCL = OneSidedUpperInterval(expectedCounts_, observedCounts_, limit);
    }
    return limit;
}

double 
BayesIntervalCalc::UpperBound(PDF &b_prob_, double expectedCounts_, double expectedCountsSigma_, int observedCounts_, double cl_){
    // Choose some signal points. Try using 25 equally spaced points between 0 and expectedCounts_ (rounded to int)
    int n_sig_points = 25;
    std::vector<double> signals;
    for(int i=0; i<n_sig_points; i++) signals.push_back(double(i*int(expectedCounts_+0.5)/n_sig_points));
    
    AxisCollection axes;
    axes.AddAxis(BinAxis("sig",0,int(expectedCounts_+0.5),n_sig_points,"sig"));
    Histogram mlhs_hist(axes);
    
    std::vector<double> mlhs;
    double sum = 0.0;
    int critBin=0;
    for(size_t i=0; i<signals.size()-1; i++){    
        if(sum > 0.99) sum = 1.0;
        else sum += BayesIntervalCalc::CumulativePosterior(b_prob_, signals[i], signals[i+1], expectedCounts_, expectedCountsSigma_, observedCounts_);
        mlhs.push_back(sum);
        if(sum < cl_) critBin++;
    }
    //Fill histogram object with cumulative posterior values to be interpolated for signal at CL=cl_
    for(int i=0; i<mlhs.size(); i++) {
        mlhs_hist.Fill(i, mlhs[i]);
    }
    
    // now interpolate between the bins    
    double upperEdge = mlhs_hist.GetAxes().GetBinHighEdge(critBin, 0);
    double lowerEdge = mlhs_hist.GetAxes().GetBinLowEdge(critBin, 0);
    double content   = mlhs_hist.GetBinContent(critBin);

    return upperEdge - (upperEdge - lowerEdge) * (sum - cl_)/content;
}

double BayesIntervalCalc::Marginalise(PDF& b_prob_, double s, double b, double width, int n){
   //If width is 0, no need to marginalise
   if(width == 0){
       return BayesIntervalCalc::PoissonProb(b_prob_,s,b,width,n);   
   } else {
       //Integrate over b with fixed s
       //Fix bkg integral to run between 0 and 1000 in 1000 steps for now
       int n_bkg_points = 1000;
       std::vector<double> bkgs;
       for(size_t i=0; i<n_bkg_points; i++) bkgs.push_back(double(i*1000/n_bkg_points));
       TH1D bkg_integral = TH1D("bkg","bkg",n_bkg_points, 0, 1000);
       for(size_t i=0; i<n_bkg_points; i++){
            bkg_integral.SetBinContent(i+1, BayesIntervalCalc::PoissonProb(b_prob_,s,bkgs[i],width,n));
       }
       return bkg_integral.Integral("width");
   }
}

double BayesIntervalCalc::CumulativePosterior(PDF& b_prob_, double signal1, double signal2, double expectedCounts_, double expectedCountsSigma_, int observedCounts_){
    // Fix signal integral to run between signal1 and signal2 in 10000 steps for now
    int n_sig_points = 10000;
    std::vector<double> sigs;
    for(size_t i=0; i<n_sig_points; i++) sigs.push_back(signal1 + double(i*(signal2-signal1)/n_sig_points));
    // Need to compute the total normalisation as well 
    // Compute normalisation by integrating signal from 0 to 1000(inf) in 10000 steps for now
    int n_norm_sig_points = 10000;
    std::vector<double> norm_sigs;
    for(size_t i=0; i<n_norm_sig_points; i++) norm_sigs.push_back(double(i*(1000.0)/n_sig_points));
    
    TH1D sig_integral = TH1D("sig","sig",n_sig_points, signal1, signal2);
    TH1D norm_sig_integral = TH1D("norm_sig","norm_sig",n_norm_sig_points, 0, 1000);

    for(size_t i=0; i<n_norm_sig_points; i++) norm_sig_integral.SetBinContent(i+1, BayesIntervalCalc::Marginalise(b_prob_,norm_sigs[i],expectedCounts_,expectedCountsSigma_,observedCounts_));
    double norm = norm_sig_integral.Integral("width");
    
    //Integrate s from signal1 to signal2
    for(size_t i=0; i<n_sig_points; i++) {
        double result = BayesIntervalCalc::Marginalise(b_prob_,sigs[i],expectedCounts_,expectedCountsSigma_,observedCounts_)/norm;
        sig_integral.SetBinContent(i+1, result );
    }
     
    return sig_integral.Integral("width");    
    
}

double BayesIntervalCalc::PoissonProb(PDF& b_prob_, double s, double b, double width, int n) {
   double computed_value = 1.0;
   //Apply a Poisson probability (s+b)^n / n! * exp(-(s+b)) and whatever background distribution is stored in b_prob
   if(width == 0){
        computed_value *= std::pow(s + b, n);
        computed_value *= 1.0/(gsl_sf_fact(n));
        computed_value *= exp(-1*(s + b));
   } else {
        std::vector<double> values_b_prob;
        values_b_prob.push_back(b);
       
        computed_value *= b_prob_( values_b_prob);
        computed_value *= std::pow(s + b, n);
        computed_value *= 1.0/(gsl_sf_fact(n));
        computed_value *= exp(-1*(s + b));
   }
   return computed_value;
}
