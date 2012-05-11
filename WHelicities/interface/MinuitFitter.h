#ifndef MinuitFitter_h
#define MinuitFitter_h

#include <TMinuitMinimizer.h>
#include "TNtuple.h"

#include <string>
#include <vector>
#include <iostream>
#include "TH1F.h"

#include "MyChi2.h"

using namespace std; 

class MinuitFitter{

  public:
   MinuitFitter();
   MinuitFitter(TH1F* datah, TH1F* signalh, TH1F* bkgh,double ff0, double ffl, double ffr, TNtuple *genttbarhisto[CosThetaNbins], int ndimen);
   ~MinuitFitter();
   double GetF0Result() { return f0result_; }
   double GetF0Error() { return ef0result_; }
   double GetFLResult() { return flresult_; }
   double GetFLError() { return eflresult_; }
   double GetFRResult() { return frresult_; }
   double GetFRError() { return efrresult_; }

  private:
   Annik::MyChi2 myChi2_;
   TMinuitMinimizer* minimizer_;
   double f0result_;
   double ef0result_;
   double flresult_;
   double eflresult_;
   double frresult_;
   double efrresult_;

};

#endif
