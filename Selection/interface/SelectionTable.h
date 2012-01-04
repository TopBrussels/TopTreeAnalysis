#ifndef SelectionTable_h
#define SelectionTable_h
#include "Selection.h"
#include "TopTreeAnalysis/Tools/interface/Efficiency.h"
#include "TopTreeAnalysis/Content/interface/Dataset.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include "TH1F.h"
#include "TFile.h"

using namespace std;

class SelectionTable {

  
 public:
  
  SelectionTable();
  SelectionTable(vector<string> listOfCuts, vector<Dataset*> listOfDatasets);
  SelectionTable(const SelectionTable& psel);
  ~SelectionTable();

  void SetLuminosity(float Lumi) {Lumi_ = Lumi;};

  float Luminosity() const {return Lumi_;};
  vector<string> ListOfCuts() const {return listOfCuts_;};
  vector<Dataset*> ListOfDatasets() const {return listOfDatasets_;};

  float BinomialeError(float p, float n) const {return (sqrt((p*(1-p)/n)*p));};
  float ErrorCalculator( float number, float p, float factor) const { return (sqrt(number*(1-p))*factor); };
  //float ErrorCalculator(float number, float factor) const {return (sqrt(number)*factor);};
  float WilsonScoreIntervalLow(float Non, float Ntot);
  float WilsonScoreIntervalHigh(float Non, float Ntot);
  float FactorCalculator(float nofevents, float xsection) const { return(Lumi_*xsection/nofevents);};


  void Clear(){;};
  void Fill(unsigned int DatasetNumber, vector<float> PassTheCuts);
  void Fill(unsigned int DatasetNumber, unsigned int CutNumber, float value);
  void TableCalculator(bool mergeTT = false, bool mergeQCD = false, bool mergeW = false, bool mergeZ = false, bool mergeST = false);
  void Scale(float Lumi);
  void WriteTable(ofstream& fout, float** listTable_, bool writeMerged = false);
  void WriteTable(ofstream& fout, float** listTable_, float** listTableError_, bool writeMerged = false);
  void WriteTable(ofstream& fout, float** listTable_, float** listTableErrorHigh_, float** listTableErrorLow_, bool writeMerged = false);
  void Write(string filename, bool WithError = true);
  void Write(ofstream& fout, bool WithError = true);

 private:

  float Lumi_; // in 1/pb
  vector<string> listOfCuts_;
  vector<Dataset*> listOfDatasets_;
  //first dimension: cuts
  //second dimension: datasets
  float**  nofEventsRaw_;	//from Fill method
  float**  nofEventsRawError_;
  float**  nofEventsRawErrorHigh_;
  float**  nofEventsRawErrorLow_;
  float**  nofEvents_;	        //rescaled
  float**  nofEventsError_;
  float**  nofEventsExpErrorHigh_;
  float**  nofEventsExpErrorLow_;
  float**  nofEventsMcErrorHigh_;
  float**  nofEventsMcErrorLow_;
  float**  cutEfficiency_;
  float**  cutEfficiencyErrorHigh_;
  float**  cutEfficiencyErrorLow_;
  float**  cutEfficiencyError_;
  float**  totalCutEfficiency_;
  float**  totalCutEfficiencyErrorHigh_;
  float**  totalCutEfficiencyErrorLow_;
  float**  totalCutEfficiencyError_;
  
  //Merged stuff
  vector<Dataset*> listOfDatasetsMerged_;
  float**  nofEventsRawMerged_;
  float**  nofEventsRawErrorMerged_;
  float**  nofEventsRawErrorHighMerged_;
  float**  nofEventsRawErrorLowMerged_;
  float**  nofEventsMerged_;
  float**  nofEventsErrorMerged_;
  float**  nofEventsExpErrorHighMerged_;
  float**  nofEventsExpErrorLowMerged_;
  float**  nofEventsMcErrorHighMerged_;
  float**  nofEventsMcErrorLowMerged_;
  float**  cutEfficiencyMerged_;
  float**  cutEfficiencyErrorMerged_;
  float**  cutEfficiencyErrorHighMerged_;
  float**  cutEfficiencyErrorLowMerged_;
  float**  totalCutEfficiencyMerged_;
  float**  totalCutEfficiencyErrorMerged_;
  float**  totalCutEfficiencyErrorHighMerged_;
  float**  totalCutEfficiencyErrorLowMerged_;

  bool merge_;
};

#endif
