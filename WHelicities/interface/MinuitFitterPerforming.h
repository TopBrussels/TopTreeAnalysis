#ifndef MinuitFitterPerforming_h
#define MinuitFitterPerforming_h

#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include "TH1F.h"
#include <sstream>

//Includes for Minuit
#include <TMinuitMinimizer.h>
#include "MinuitFitter.h"
#include "BTagName.h"


using namespace std;

class MinuitFitterPerforming{

 public:
  MinuitFitterPerforming();
  MinuitFitterPerforming(string Name, int CosThetaBinNumber, int ndimen, int iDataSet, vector<float> CosThetaValuesKinFit[], vector<float> CosThGenKinFit[], vector<float> EventCorrectionWeightKinFit[], string dataSetName, bool semiMuon, bool semiElectron, bool SignalOnly, bool DataResults, bool JESResults, bool JERResults, bool WSystResults, bool TTScalingResults, bool TTMatchingResults, bool PUResults, bool UnclusEnergyResults, bool TopMassResults, bool TriggEvtSelResults, bool bTagResults, int dataSetsSize);

 private:
  std::string CosThetaString;
  std::string CosThetaDataString;
  std::string CosThetaJESPlusString;
  std::string CosThetaJESMinusString;
  std::string CosThetaJERPlusString;
  std::string CosThetaJERMinusString;
  std::string CosThetaWPlusString;
  std::string CosThetaWMinusString;
  std::string CosThetaTTScalingUpString;
  std::string CosThetaTTScalingDownString;
  std::string CosThetaTTMatchingUpString;
  std::string CosThetaTTMatchingDownString;
  std::string CosThetaPUPlusString;
  std::string CosThetaPUMinusString;
  std::string CosThetaUnclusEnergyPlusString;
  std::string CosThetaUnclusEnergyMinusString;
  std::string CosThetaTopMassPlusString;
  std::string CosThetaTopMassMinusString;
  std::string CosThetaTriggEvtSelPlusString;
  std::string CosThetaTriggEvtSelMinusString;
  std::string CosThetaBTagSystPlusString;
  std::string CosThetaBTagSystMinusString;  
  std::string CosThetaSignalString;
  std::string CosThetaBckgString;
  
};

#endif
