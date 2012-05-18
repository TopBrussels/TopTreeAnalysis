#ifndef MyChi2_h
#define MyChi2_h

#include <TMinuitMinimizer.h>
#include "TNtuple.h"
#include "TBranch.h"

#include <string>
#include <vector>
#include <iostream>
#include "TH1F.h"
#include "TMath.h"

const int CosThetaNbins = 15;       

using namespace std; 

namespace Annik
{

class MyChi2 : public ROOT::Math::IMultiGenFunction {
public:
  // Mandatory functions 
  
  virtual MyChi2 * Clone() const {return new MyChi2;};
  virtual unsigned int NDim() const {return ndimen_;}
  
  // Constructors
  MyChi2(){};
  MyChi2(TH1F* datah, TH1F* signalh, TH1F* bkgh,double ff0, double ffl, double ffr) : datah_(datah),signalh_(signalh),bkgh_(bkgh),ff0_(ff0),ffl_(ffl),ffr_(ffr){};  
  void SetNDim(unsigned int nDim) { ndimen_ = nDim; }
  void SetNTuple(TNtuple *genttbarhisto[CosThetaNbins]);

  double DoEval(const double* x) const;

private:
  
  TH1F* datah_;
  TH1F* signalh_;
  TH1F* bkgh_;
  double ff0_;
  double ffl_;
  double ffr_;
  unsigned int ndimen_;
  TNtuple **genttbarhisto_;//[CosThetaNbins];
};

}

#endif
