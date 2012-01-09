#ifndef BTagCosThetaCalculation_h
#define BTagCosThetaCalculation_h

#include "TLorentzVector.h"

#include <string>
#include <vector>
#include <iostream>

using namespace std;

class BTagCosThetaCalculation{

  public:
   BTagCosThetaCalculation();
   ~BTagCosThetaCalculation();
   
   float Calculation(TLorentzVector muon, TLorentzVector Neutrino, TLorentzVector leptonicBJet);

   float CalcOrigKins(int BLeptonicIndex, int BHadronicIndex,  TLorentzVector muon, vector<TLorentzVector> selectedJets, float MassW, float MassTop);

};

#endif

