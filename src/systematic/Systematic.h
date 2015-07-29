/**************************************************************************/
/* Responsible for constructing a response matrix and forwarding the call */
/**************************************************************************/
#ifndef __SYSTEMATIC__
#define __SYSTEMATIC__
#include "ResponseMatrix.h"
#include "../pdf/binned/BinnedPdf.h"

class Systematic{
 public:
    Systematic()   {}
    virtual ~Systematic()  {}

    BinnedPdf 
    operator()(const BinnedPdf& pdf_, const std::vector<size_t>& indicies_) const;
        
    void 
    SetResponse(const std::vector<std::vector<double> >& responseMatrix_);

    const ResponseMatrix& GetResponse() const;
        
 protected:
    ResponseMatrix fResponseMatrix;
};
#endif