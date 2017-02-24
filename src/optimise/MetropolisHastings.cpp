#include <MetropolisHastings.h>
#include <TestStatistic.h>
#include <Exceptions.h>
#include <HistTools.h>
#include <Formatter.hpp>
#include <Combinations.hpp>
#include <Rand.h>
#include <iostream>
#include <math.h> //sqrt

const std::vector<double>&
MetropolisHastings::GetSigmas() const{
    return fSigmas;
}

void
MetropolisHastings::SetSigmas(const std::vector<double>& sigmas_){
    fSigmas = sigmas_;
}

unsigned
MetropolisHastings::GetBurnIn() const{
    return fBurnIn;
}

void
MetropolisHastings::SetBurnIn(unsigned burnIn_){
    fBurnIn = burnIn_;
}

unsigned
MetropolisHastings::GetThinFactor() const{
    return fThinFactor;
}

void
MetropolisHastings::SetThinFactor(unsigned thinFac_){
    fThinFactor = thinFac_;
}

unsigned
MetropolisHastings::GetMaxIter() const{
    return fMaxIter;
}

void
MetropolisHastings::SetMaxIter(unsigned maxIter_){
    fMaxIter = maxIter_;
}

double
MetropolisHastings::GetRejectionRate() const{
    return fRejectionRate;
}

std::vector<double>
MetropolisHastings::GetMinima() const{
    return fMinima;
}

void
MetropolisHastings::SetMinima(const std::vector<double>& minima_){
    fMinima = minima_;
}


std::vector<double>
MetropolisHastings::GetMaxima() const{
    return fMaxima;
}

void
MetropolisHastings::SetMaxima(const std::vector<double>& maxima_){
    fMaxima = maxima_;
}


void
MetropolisHastings::SetFlipSign(bool b_){
    fFlipSign = b_;
}

bool
MetropolisHastings::GetFlipSign() const{
    return fFlipSign;
}

bool
MetropolisHastings::GetTestStatLogged() const{
    return fTestStatLogged;
}

void
MetropolisHastings::SetTestStatLogged(bool b_){
    fTestStatLogged = b_;
}


bool 
MetropolisHastings::GetSaveFullHistogram() const{
    return fSaveFullHistogram;
}

void
MetropolisHastings::SetSaveFullHistogram(bool b_){
    fSaveFullHistogram = b_;
}

void
MetropolisHastings::SetHistogramAxes(const AxisCollection& template_){
    fHistogramAxes = template_;
}
             
AxisCollection
MetropolisHastings::GetHistogramAxes() const{
    return fHistogramAxes;
}                                                        

void
MetropolisHastings::SetInitialTrial(const std::vector<double>& trial_){
    fInitialTrial = trial_;
}

std::vector<double>
MetropolisHastings::GetInitialTrial() const{
    return fInitialTrial;
}

const FitResult&
MetropolisHastings::Optimise(TestStatistic* testStat_){
    testStat_->RegisterFitComponents();
    pTestStatistic = testStat_;
    fNDims = fMinima.size();
    fRejectionRate = 0;
    fBestFit.resize(fNDims, 0);
    fMaxVal = 0;

    // Check initialisation
    size_t nParams = testStat_ -> GetParameterCount();
    if(   (fInitialTrial.size() != nParams && fInitialTrial.size())
          || fMaxima.size() != nParams
          || fMinima.size() != nParams
          || fSigmas.size() != nParams
          )
        throw LogicError(Formatter()
                         << "Metropolis Hastings initialisation error - Got "
                         << fMinima.size() << " Minima, "
                         << fMaxima.size() << " Maxima and "
                         << fSigmas.size() << " Sigmas"
                         << " - Need one per fit parameter (" << nParams
                         << ")"
                         << "\n If initial trial is specified it should be one per fit parameter"
                         );    

    InitialiseHistograms();

    // 1. Choose a random starting point or the user defined one
    std::vector<double> currentStep;
    if(fInitialTrial.size())
      currentStep = fInitialTrial;
    else
      for(size_t i = 0; i < fMinima.size(); i++){
          currentStep.push_back(fMinima.at(i) + Rand::Uniform() * (fMaxima.at(i) - fMinima.at(i)));
      }

    std::cout << "Metropolis Hastings::Initial Position @ \t\t";
    for(size_t i = 0; i < currentStep.size(); i++)
        std::cout << currentStep.at(i) << ",";
    std::cout << std::endl;

    fSample.empty();
    //    fSample.reserve(fMaxIter - fBurnIn);
   
    // 2. Loop step through the space a fixed number of times and
    for(unsigned i = 0; i < fMaxIter; i++){
        // a. Save the point in question if you are past burn-in phase and according to thinning
        if (i > fBurnIn && !(i%fThinFactor)){
          //            fSample.push_back(currentStep);

            if(fSaveFullHistogram)
                fHist.Fill(currentStep);

            else{
                FillProjections(currentStep);
            }
        }
      
        if(!(i%100000))
            std::cout << i << "  /  " << fMaxIter 
              << "\t" << (i - fRejectionRate) / static_cast<double>(i) << std::endl;

        // b. Propose a new step according to a random jump distribution
        std::vector<double> proposedStep = JumpDraw(currentStep);

        // c. Test for acceptance
        if (StepAccepted(currentStep, proposedStep))
            currentStep = proposedStep;
        else
            fRejectionRate++;
    }
    fRejectionRate/= static_cast<double>(fMaxIter);
    std::cout << "Metropolis Hastings:: acceptance rate = " << 1 - fRejectionRate << std::endl;

    fFitResult.SetBestFit(fBestFit);
    fFitResult.SetParameterNames(pTestStatistic->GetParameterNames());
    fFitResult.SetStatSample(fSample);
    fFitResult.SetStatSpace(fHist);
    fFitResult.SetValid(true);
    if(!fSaveFullHistogram)
      SaveProjections();

    return fFitResult;
}

bool
MetropolisHastings::StepAccepted(const std::vector<double>& thisStep_, 
                                 const std::vector<double>& proposedStep_){

    // dont step outside of the fit region 
    for (size_t i = 0; i < fMinima.size(); i++){
        if (proposedStep_.at(i) < fMinima.at(i) || proposedStep_.at(i) > fMaxima.at(i))
            return false;
    }

    pTestStatistic -> SetParameters(thisStep_);
    double thisVal = pTestStatistic -> Evaluate();
    
    pTestStatistic -> SetParameters(proposedStep_);
    double proposedVal = pTestStatistic -> Evaluate();

    if(fFlipSign){
        thisVal = -thisVal;
        proposedVal = -proposedVal;
    }

    if(thisVal > fMaxVal || fMaxVal == 0.0){
        fMaxVal = thisVal;
        fBestFit = thisStep_;
    }
    
    double acceptanceParam = 0;
    if(fTestStatLogged)
        acceptanceParam = JumpProbRatio(thisStep_, proposedStep_) * exp(proposedVal - thisVal);

    else
        acceptanceParam = JumpProbRatio(thisStep_, proposedStep_) * proposedVal/thisVal;
    
    if (acceptanceParam > 1)
        return true;

    else if (acceptanceParam >= Rand::Uniform())
        return true;

    else{
        return false;
    }
}

std::vector<double> 
MetropolisHastings::JumpDraw(const std::vector<double>& thisStep_) const{
    std::vector<double> newStep(thisStep_);
    for(size_t i = 0; i < newStep.size(); i++){
        if(fMinima.at(i) == fMaxima.at(i))
            newStep[i] = fMinima.at(i);
        else
            newStep[i] = thisStep_.at(i) + Rand::Gaus(0, fSigmas.at(i));
    }
    return newStep;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// std::vector<double>                                                                             //
// MetropolisHastings::JumpDraw(const std::vector<double>& thisStep_) const{                       //
//     double stepSize = Rand::Gaus(0, fSigma);                                                    //
//                                                                                                 //
//     // choose random direction on NParam -1 Sphere according to                                 //
//     // http://mathworld.wolfram.com/HyperspherePointPicking.html (Muller 1959, Marsaglia 1972)  //
//                                                                                                 //
//     std::vector<double> newPoint;                                                               //
//     double norm = 0;                                                                            //
//     for(size_t i = 0; i < fNDims; i++){                                                         //
//         double component = Rand::Gaus();                                                        //
//         newPoint.push_back(component);                                                          //
//         norm += component * component;                                                          //
//     }                                                                                           //
//                                                                                                 //
//     norm = sqrt(norm);                                                                          //
//     for(size_t i = 0; i < fNDims; i++){                                                         //
//         newPoint[i] *= stepSize/norm;                                                           //
//     }                                                                                           //
//                                                                                                 //
//     // Add to orginal vector and return                                                         //
//     for(size_t i = 0; i < fNDims; i++)                                                          //
//         newPoint[i] += thisStep_.at(i);                                                         //
//                                                                                                 //
//     return newPoint;                                                                            //
// }                                                                                               //
/////////////////////////////////////////////////////////////////////////////////////////////////////


void
MetropolisHastings::InitialiseHistograms(){
    // Create histogram(s)
    // If the binning is specified, use that
    AxisCollection histAxes;
    if(fHistogramAxes.GetNDimensions()){
        if(fHistogramAxes.GetNDimensions() != pTestStatistic->GetParameterCount())
            throw DimensionError("MetropolisHastings::InitialiseHistograms",
                                 pTestStatistic->GetParameterCount(),
                                 fHistogramAxes.GetNDimensions(),
                                 "Template axes - one axis for each parameter"
                                 );
        histAxes = fHistogramAxes;
    }

    // otherwise take a guess
    else{
        std::vector<std::string> paramNames = pTestStatistic->GetParameterNames();
        for(size_t i = 0; i < fMinima.size(); i++){
            histAxes.AddAxis(BinAxis(paramNames.at(i), fMinima.at(i), fMaxima.at(i), 
                                 int(pow(fMaxIter, 1./fMinima.size()))));
        }
    }

    // Set up the histogram, either with a big ND histogram
    if(fSaveFullHistogram){
        fHist = Histogram(histAxes);
    }
    
    
    // otherwise just save the visualisable projections    
    else{
        // 2D
        std::vector<size_t> dimensionList = SequentialElements<size_t>(size_t(0), histAxes.GetNDimensions());
        f2DProjectionIndices  = FixedLengthCombinationsNoDuplicates<size_t>(dimensionList, 2);
        f2DProjections        = HistTools::MakeAllHists(histAxes, f2DProjectionIndices);
        
        // 1D
        f1DProjections = HistTools::MakeAllHists(histAxes, 
                                                 FixedLengthCombinationsNoDuplicates<size_t>(dimensionList, 1));
    }
}


void
MetropolisHastings::FillProjections(const std::vector<double>& params_){
    // 1D
    for(size_t i = 0; i < params_.size(); i++){
        f1DProjections[i].Fill(params_.at(i));
    }

    // 2D
    std::vector<double> fillVals(2, 0);
    for(size_t i = 0; i < f2DProjectionIndices.size(); i++){
        const std::vector<size_t>& dimsToFill = f2DProjectionIndices.at(i);
        Histogram& histToFill = f2DProjections[i];
        fillVals[0] = params_.at(dimsToFill.at(0));        
        fillVals[1] = params_.at(dimsToFill.at(1));
        histToFill.Fill(fillVals);
    }
}

void
MetropolisHastings::SaveProjections(){
  std::map<std::string, Histogram> map1d;
  for(size_t i = 0; i < f1DProjections.size(); i++){
    const Histogram& hist = f1DProjections.at(i);
    std::string name = hist.GetAxes().GetAxis(0).GetName();
    map1d[name + "_projection"] = hist;
  }

  std::map<std::string, Histogram> map2d;
  for(size_t i = 0; i < f2DProjections.size(); i++){
    const Histogram& hist = f2DProjections.at(i);
    std::string name0 = hist.GetAxes().GetAxis(0).GetName();
    std::string name1 = hist.GetAxes().GetAxis(1).GetName();

    map2d[name0 + "_" + name1 + "_projection"] = hist;
  }

  fFitResult.Set1DProjections(map1d);
  fFitResult.Set2DProjections(map2d);
  
}
