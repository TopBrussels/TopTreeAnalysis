#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

// Root stuff
#include "TROOT.h"
#include "TH1F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TF2.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TBranch.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TDirectory.h"
#include "TNtuple.h"
#include "TProfile.h"

#include "../Tools/interface/PlottingTools.h"
#include "../Tools/interface/MultiSamplePlot.h"
#include "../Content/interface/AnalysisEnvironment.h"
#include "../Content/interface/Dataset.h"
#include "../WHelicities/interface/WTree.h"
#include "../JESMeasurement/interface/FullKinFit.h"
#include "../MCInformation/interface/LumiReWeighting.h"
#include "../MCInformation/interface/Lumi3DReWeighting.h"
#include "../Selection/interface/SelectionTable.h"

#include "Style.C"

//Includes for made classes:
#include "../WHelicities/interface/BTagName.h"
#include "../WHelicities/interface/BTagJetSelection.h"
#include "../WHelicities/interface/BTagCosThetaCalculation.h"
#include "../WHelicities/interface/MinuitFitter.h"
//#include "../WHelicities/interface/KinematicFitClass.h"

//Includes for Minuit
#include <TMinuitMinimizer.h>
//#include "TNtuple.h"  //A simple tree restricted to a list of float variables only. 

const int ndimen = 3;

const int CosThetaBinNumber = 30;       
const int NumberSSVHP=14;
const int NumberSSVHE=14;
const int NumberTCHP=14;
const int NumberTCHE=14;
const int NumberCSV=14;
const int NumberJP=14;
const int NumberJBP=14;
int NumberOfHelicityBins=100; 
TNtuple *genttbarhistoHadr[CosThetaBinNumber];  //This is the vector of ntuples containing the generated values of cos theta* for each cos theta* reconstructed bin
TNtuple *genttbarhistoHadrAndLeptWOnly[CosThetaBinNumber];  //This is the vector of ntuples containing the generated values of cos theta* for each cos theta* reconstructed bin
TNtuple *genttbarhistoHadrAndLept[CosThetaBinNumber];  //This is the vector of ntuples containing the generated values of cos theta* for each cos theta* reconstructed bin

using namespace std;
using namespace reweight;

/// TGraphAsymmErrors
map<string,TGraphAsymmErrors*> graphAsymmErr;
map<string,TGraphErrors*> graphErr;

/// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;

/// MultiSamplePlot
map<string,MultiSamplePlot*> MSPlot;

int main (int argc, char *argv[])
{
  clock_t start = clock();

  cout << "*******************************************************" << endl;
  cout << " Beginning of the program for W Helicities Analysis " << endl;
  cout << "   --> bTag in jet-quark matching ! " << endl;
  cout << "*******************************************************" << endl;
 
  //  setTDRStyle();
  setMyStyle();
  
  string pathPNG = "Plots/bTagUsedInJetMatching/";
  mkdir(pathPNG.c_str(),0777);
  mkdir((pathPNG+"MSPlot/").c_str(),0777);
  mkdir((pathPNG+"CosThetaPlots/").c_str(),0777);
  
  float Luminosity;
  vector<string> inputWTree, nameDataSet;

  //--------------------------------------------------------------
  //     DataSamples needed for 2011 Data for muon channel:     --
  //--------------------------------------------------------------
  
  string decayChannel;
  bool semiMuon = false;
  bool semiElectron = true;
  if(semiMuon == true){decayChannel = "SemiMu";}
  else if(semiElectron == true){decayChannel = "SemiEl";}

  string UsedTrigger;
  bool IsoMu172024Trigger = false;
  bool TriCentralJet30Trigger = true;
  if(IsoMu172024Trigger == true){
    UsedTrigger = "IsoMu172024Trigger";
    Luminosity = 4938.985;
    UsedTrigger = UsedTrigger;
  }
  else if(TriCentralJet30Trigger == true){
    UsedTrigger = "TriCentralJet30Trigger";
    if(semiMuon == true){
      Luminosity = 4656.959;
    }
    else if(semiElectron == true){
      Luminosity = 4966.709;
    }
  }
  else if(TriCentralJet30Trigger == false && IsoMu172024Trigger == false){
    UsedTrigger = "NoTrigger";
  }  
  cout << "Executing the W Helicities analysis for an integrated luminosity of " << Luminosity << " pb^-1" << endl;

  //Booleans to load in different root files
  bool SignalOnly = false;
  bool DataResults = true;
  bool JESResults = true;
  bool JERResults = true;
  bool WSystResults = true;
  bool TTScalingResults = true;
  bool TTMatchingResults = true;
  bool PUResults = true;
  bool UnclusEnergyResults = true;
  bool TopMassResults = true;
  bool TriggEvtSelResults = true;
  cout << " Obtaining results for : Data = " << DataResults << " , JES = " << JESResults << " , JER = " << JERResults << " , WSyst = " << WSystResults << endl;

  //Different KinFit configurations studied:
  bool HadronicOnly = true;
  bool HadronicAndLeptonicWOnly = true;
  bool HadronicAndLeptonic = true;

  //Boolean to put MinuitFitter on and off:
  bool PerformMinuit = true;

  //Boolean to create MSPlots:
  bool CreateMSPlots = false;
 
  //1) Nominal samples:
  if(SignalOnly == false){
    if(DataResults == true){
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_"+"Data_"+decayChannel+".root").c_str());
      nameDataSet.push_back("Data");
    }
    inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ST_SingleTop_tChannel_tbar_"+decayChannel+".root").c_str());
    nameDataSet.push_back("ST_SingleTop_tChannel_tbar");
    //    inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ST_SingleTop_tChannel_t_"+decayChannel+".root").c_str());
    //    nameDataSet.push_back("ST_SingleTop_tChannel_t");
    inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ST_SingleTop_tWChannel_tbar_"+decayChannel+".root").c_str());
    nameDataSet.push_back("ST_SingleTop_tWChannel_tbar");
    inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ST_SingleTop_tWChannel_t_"+decayChannel+".root").c_str());
    nameDataSet.push_back("ST_SingleTop_tWChannel_t");
    inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_TTbarJets_Other_"+decayChannel+".root").c_str());
    nameDataSet.push_back("TTbarJets_Other");
    if(semiMuon == true){
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_TTbarJets_SemiEl_"+decayChannel+".root").c_str());  //In muon channel case SemiEl is considered as background
      nameDataSet.push_back("TTbarJets_SemiEl");
    }
    else if(semiElectron == true){
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_TTbarJets_SemiMuon_"+decayChannel+".root").c_str());  //In electron channel case SemiMu is considered as background
      nameDataSet.push_back("TTbarJets_SemiMuon");
    }
    inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_WJets_"+decayChannel+".root").c_str());
    nameDataSet.push_back("WJets_Nominal");
    inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ZJets_"+decayChannel+".root").c_str());
    nameDataSet.push_back("ZJets");  
    
    //2) JES Plus/Min samples:
    if(JESResults == true ){
      //Consider JESPlus samples:
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_JESPlus_1Sig_ST_SingleTop_tChannel_tbar_"+decayChannel+".root").c_str());
      nameDataSet.push_back("JESPlus_ST_SingleTop_tChannel_tbar");
      //      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_JESPlus_1Sig_ST_SingleTop_tChannel_t_"+decayChannel+".root").c_str());
      //      nameDataSet.push_back("JES_ST_SingleTop_tChannel_t");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_JESPlus_1Sig_ST_SingleTop_tWChannel_tbar_"+decayChannel+".root").c_str());
      nameDataSet.push_back("JESPlus_ST_SingleTop_tWChannel_tbar");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_JESPlus_1Sig_ST_SingleTop_tWChannel_t_"+decayChannel+".root").c_str());
      nameDataSet.push_back("JESPlus_ST_SingleTop_tWChannel_t");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_JESPlus_1Sig_TTbarJets_Other_"+decayChannel+".root").c_str());
      nameDataSet.push_back("JESPlus_TTbarJets_Other");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_JESPlus_1Sig_TTbarJets_SemiEl_"+decayChannel+".root").c_str());
      nameDataSet.push_back("JESPlus_TTbarJets_SemiEl");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_JESPlus_1Sig_WJets_"+decayChannel+".root").c_str());    
      nameDataSet.push_back("JESPlus_WJets");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_JESPlus_1Sig_ZJets_"+decayChannel+".root").c_str());
      nameDataSet.push_back("JESPlus_ZJets");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_JESPlus_1Sig_TTbarJets_SemiMuon_"+decayChannel+".root").c_str());
      nameDataSet.push_back("JESPlus_TTbarJets_SemiMuon");

      //Consider JESMinus samples:
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_JESMinus_1Sig_ST_SingleTop_tChannel_tbar_"+decayChannel+".root").c_str());
      nameDataSet.push_back("JESMinus_ST_SingleTop_tChannel_tbar");
      //      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_JESMinus_1Sig_ST_SingleTop_tChannel_t_"+decayChannel+".root").c_str());
      //      nameDataSet.push_back("JES_ST_SingleTop_tChannel_t");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_JESMinus_1Sig_ST_SingleTop_tWChannel_tbar_"+decayChannel+".root").c_str());
      nameDataSet.push_back("JESMinus_ST_SingleTop_tWChannel_tbar");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_JESMinus_1Sig_ST_SingleTop_tWChannel_t_"+decayChannel+".root").c_str());
      nameDataSet.push_back("JESMinus_ST_SingleTop_tWChannel_t");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_JESMinus_1Sig_TTbarJets_Other_"+decayChannel+".root").c_str());
      nameDataSet.push_back("JESMinus_TTbarJets_Other");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_JESMinus_1Sig_TTbarJets_SemiEl_"+decayChannel+".root").c_str());
      nameDataSet.push_back("JESMinus_TTbarJets_SemiEl");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_JESMinus_1Sig_WJets_"+decayChannel+".root").c_str());    
      nameDataSet.push_back("JESMinus_WJets");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_JESMinus_1Sig_ZJets_"+decayChannel+".root").c_str());
      nameDataSet.push_back("JESMinus_ZJets");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_JESMinus_1Sig_TTbarJets_SemiMuon_"+decayChannel+".root").c_str());
      nameDataSet.push_back("JESMinus_TTbarJets_SemiMuon");
    }
    
    //3) JER Plus/Min samples:
    if(JERResults == true){
      //Consider JERPlus samples:
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_JERPlus_ST_SingleTop_tChannel_tbar_"+decayChannel+".root").c_str());
      nameDataSet.push_back("JERPlus_ST_SingleTop_tChannel_tbar");
      //      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_JERPlus_ST_SingleTop_tChannel_t_"+decayChannel+".root").c_str());
      //      nameDataSet.push_back("JER_ST_SingleTop_tChannel_t");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_JERPlus_ST_SingleTop_tWChannel_tbar_"+decayChannel+".root").c_str());
      nameDataSet.push_back("JERPlus_ST_SingleTop_tWChannel_tbar");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_JERPlus_ST_SingleTop_tWChannel_t_"+decayChannel+".root").c_str());
      nameDataSet.push_back("JERPlus_ST_SingleTop_tWChannel_t");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_JERPlus_TTbarJets_Other_"+decayChannel+".root").c_str());
      nameDataSet.push_back("JERPlus_TTbarJets_Other");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_JERPlus_TTbarJets_SemiEl_"+decayChannel+".root").c_str());
      nameDataSet.push_back("JERPlus_TTbarJets_SemiEl");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_JERPlus_WJets_"+decayChannel+".root").c_str());    
      nameDataSet.push_back("JERPlus_WJets");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_JERPlus_ZJets_"+decayChannel+".root").c_str());
      nameDataSet.push_back("JERPlus_ZJets");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_JERPlus_TTbarJets_SemiMuon_"+decayChannel+".root").c_str());
      nameDataSet.push_back("JERPlus_TTbarJets_SemiMuon");

      //Consider JERMinus samples:
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_JERMinus_ST_SingleTop_tChannel_tbar_"+decayChannel+".root").c_str());
      nameDataSet.push_back("JERMinus_ST_SingleTop_tChannel_tbar");
      //      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_JERMinus_ST_SingleTop_tChannel_t_"+decayChannel+".root").c_str());
      //      nameDataSet.push_back("JER_ST_SingleTop_tChannel_t");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_JERMinus_ST_SingleTop_tWChannel_tbar_"+decayChannel+".root").c_str());
      nameDataSet.push_back("JERMinus_ST_SingleTop_tWChannel_tbar");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_JERMinus_ST_SingleTop_tWChannel_t_"+decayChannel+".root").c_str());
      nameDataSet.push_back("JERMinus_ST_SingleTop_tWChannel_t");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_JERMinus_TTbarJets_Other_"+decayChannel+".root").c_str());
      nameDataSet.push_back("JERMinus_TTbarJets_Other");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_JERMinus_TTbarJets_SemiEl_"+decayChannel+".root").c_str());
      nameDataSet.push_back("JERMinus_TTbarJets_SemiEl");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_JERMinus_WJets_"+decayChannel+".root").c_str());    
      nameDataSet.push_back("JERMinus_WJets");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_JERMinus_ZJets_"+decayChannel+".root").c_str());
      nameDataSet.push_back("JERMinus_ZJets");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_JERMinus_TTbarJets_SemiMuon_"+decayChannel+".root").c_str());
      nameDataSet.push_back("JERMinus_TTbarJets_SemiMuon");
    }
    
    //4) WJets systematics:  
    if(WSystResults == true){
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_WJets_"+decayChannel+".root").c_str());  
      nameDataSet.push_back("WSystPlus_WJets");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_WJets_"+decayChannel+".root").c_str());  
      nameDataSet.push_back("WSystMinus_WJets");
    }

    //5) TT scaling systematics:
    if(TTScalingResults == true){
      if(semiMuon == true){
	inputWTree.push_back("WTree/KinFit_WTree_IsoMu172024Trigger_TTbarJets_ScalingUp_SemiMu.root");
	nameDataSet.push_back("TTScalingUp");
	inputWTree.push_back("WTree/KinFit_WTree_IsoMu172024Trigger_TTbarJets_ScalingDown_SemiMu.root");
	nameDataSet.push_back("TTScalingDown");
      }
      else if(semiElectron == true){
	inputWTree.push_back("WTree/KinFit_WTree_TriCentralJet30Trigger_TTbarJets_ScalingUp_SemiEl.root");
	nameDataSet.push_back("TTScalingUp");
	inputWTree.push_back("WTree/KinFit_WTree_TriCentralJet30Trigger_TTbarJets_ScalingDown_SemiEl.root");
	nameDataSet.push_back("TTScalingDown");
      }
    }
    
    //6) TT matching systematics:
    if(TTMatchingResults == true){
      if(semiMuon == true){
	inputWTree.push_back("WTree/KinFit_WTree_IsoMu172024Trigger_TTbarJets_MatchingUp_SemiMu.root");
	nameDataSet.push_back("TTMatchingUp");
	inputWTree.push_back("WTree/KinFit_WTree_IsoMu172024Trigger_TTbarJets_MatchingDown_SemiMu.root");
	nameDataSet.push_back("TTMatchingDown");
      }
      else if(semiElectron == true){
	inputWTree.push_back("WTree/KinFit_WTree_TriCentralJet30Trigger_TTbarJets_MatchingUp_SemiEl.root");
	nameDataSet.push_back("TTMatchingUp");
	inputWTree.push_back("WTree/KinFit_WTree_TriCentralJet30Trigger_TTbarJets_MatchingDown_SemiEl.root");
	nameDataSet.push_back("TTMatchingDown");
      }
    }
    
    //7) PU systematics:
    if(PUResults == true){
      //Consider PU Plus samples:      
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ST_SingleTop_tChannel_tbar_"+decayChannel+".root").c_str());
      nameDataSet.push_back("PUPlus_ST_SingleTop_tChannel_tbar");
      //      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ST_SingleTop_tChannel_t_"+decayChannel+".root").c_str());
      //      nameDataSet.push_back("PUPlus_ST_SingleTop_tChannel_t");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ST_SingleTop_tWChannel_tbar_"+decayChannel+".root").c_str());
      nameDataSet.push_back("PUPlus_ST_SingleTop_tWChannel_tbar");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ST_SingleTop_tWChannel_t_"+decayChannel+".root").c_str());
      nameDataSet.push_back("PUPlus_ST_SingleTop_tWChannel_t");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_TTbarJets_Other_"+decayChannel+".root").c_str());
      nameDataSet.push_back("PUPlus_TTbarJets_Other");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_TTbarJets_SemiEl_"+decayChannel+".root").c_str());
      nameDataSet.push_back("PUPlus_TTbarJets_SemiEl");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_WJets_"+decayChannel+".root").c_str());    
      nameDataSet.push_back("PUPlus_WJets");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ZJets_"+decayChannel+".root").c_str());
      nameDataSet.push_back("PUPlus_ZJets");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_TTbarJets_SemiMuon_"+decayChannel+".root").c_str());
      nameDataSet.push_back("PUPlus_TTbarJets_SemiMuon");

      //Consider PU Minus samples:
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ST_SingleTop_tChannel_tbar_"+decayChannel+".root").c_str());
      nameDataSet.push_back("PUMinus_ST_SingleTop_tChannel_tbar");
      //      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ST_SingleTop_tChannel_t_"+decayChannel+".root").c_str());
      //      nameDataSet.push_back("PUMinus_ST_SingleTop_tChannel_t");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ST_SingleTop_tWChannel_tbar_"+decayChannel+".root").c_str());
      nameDataSet.push_back("PUMinus_ST_SingleTop_tWChannel_tbar");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ST_SingleTop_tWChannel_t_"+decayChannel+".root").c_str());
      nameDataSet.push_back("PUMinus_ST_SingleTop_tWChannel_t");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_TTbarJets_Other_"+decayChannel+".root").c_str());
      nameDataSet.push_back("PUMinus_TTbarJets_Other");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_TTbarJets_SemiEl_"+decayChannel+".root").c_str());
      nameDataSet.push_back("PUMinus_TTbarJets_SemiEl");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_WJets_"+decayChannel+".root").c_str());    
      nameDataSet.push_back("PUMinus_WJets");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ZJets_"+decayChannel+".root").c_str());
      nameDataSet.push_back("PUMinus_ZJets");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_TTbarJets_SemiMuon_"+decayChannel+".root").c_str());
      nameDataSet.push_back("PUMinus_TTbarJets_SemiMuon");
    }

    //8) Unclustered energy
    if(UnclusEnergyResults == true){
      //Consider UnclusEnergy Plus samples:
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_UnclEnergyPlus10Perc_ST_SingleTop_tChannel_tbar_"+decayChannel+".root").c_str());
      nameDataSet.push_back("UnclusEnergyPlus_ST_SingleTop_tChannel_tbar");
      //      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_UnclEnergyPlus10Perc_ST_SingleTop_tChannel_t_"+decayChannel+".root").c_str());
      //      nameDataSet.push_back("UnclusEnergyPlus_ST_SingleTop_tChannel_t");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_UnclEnergyPlus10Perc_ST_SingleTop_tWChannel_tbar_"+decayChannel+".root").c_str());
      nameDataSet.push_back("UnclusEnergyPlus_ST_SingleTop_tWChannel_tbar");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_UnclEnergyPlus10Perc_ST_SingleTop_tWChannel_t_"+decayChannel+".root").c_str());
      nameDataSet.push_back("UnclusEnergyPlus_ST_SingleTop_tWChannel_t");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_UnclEnergyPlus10Perc_TTbarJets_Other_"+decayChannel+".root").c_str());
      nameDataSet.push_back("UnclusEnergyPlus_TTbarJets_Other");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_UnclEnergyPlus10Perc_TTbarJets_SemiEl_"+decayChannel+".root").c_str());
      nameDataSet.push_back("UnclusEnergyPlus_TTbarJets_SemiEl");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_UnclEnergyPlus10Perc_WJets_"+decayChannel+".root").c_str());    
      nameDataSet.push_back("UnclusEnergyPlus_WJets");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_UnclEnergyPlus10Perc_ZJets_"+decayChannel+".root").c_str());
      nameDataSet.push_back("UnclusEnergyPlus_ZJets");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_UnclEnergyPlus10Perc_TTbarJets_SemiMuon_"+decayChannel+".root").c_str());
      nameDataSet.push_back("UnclusEnergyPlus_TTbarJets_SemiMuon");

      //Consider UnclusEnergy Minus samples:
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_UnclEnergyMinus10Perc_ST_SingleTop_tChannel_tbar_"+decayChannel+".root").c_str());
      nameDataSet.push_back("UnclusEnergyMinus_ST_SingleTop_tChannel_tbar");
      //      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_UnclEnergyMinus10Perc_ST_SingleTop_tChannel_t_"+decayChannel+".root").c_str());
      //      nameDataSet.push_back("UnclusEnergy_ST_SingleTop_tChannel_t");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_UnclEnergyMinus10Perc_ST_SingleTop_tWChannel_tbar_"+decayChannel+".root").c_str());
      nameDataSet.push_back("UnclusEnergyMinus_ST_SingleTop_tWChannel_tbar");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_UnclEnergyMinus10Perc_ST_SingleTop_tWChannel_t_"+decayChannel+".root").c_str());
      nameDataSet.push_back("UnclusEnergyMinus_ST_SingleTop_tWChannel_t");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_UnclEnergyMinus10Perc_TTbarJets_Other_"+decayChannel+".root").c_str());
      nameDataSet.push_back("UnclusEnergyMinus_TTbarJets_Other");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_UnclEnergyMinus10Perc_TTbarJets_SemiEl_"+decayChannel+".root").c_str());
      nameDataSet.push_back("UnclusEnergyMinus_TTbarJets_SemiEl");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_UnclEnergyMinus10Perc_WJets_"+decayChannel+".root").c_str());    
      nameDataSet.push_back("UnclusEnergyMinus_WJets");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_UnclEnergyMinus10Perc_ZJets_"+decayChannel+".root").c_str());
      nameDataSet.push_back("UnclusEnergyMinus_ZJets");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_UnclEnergyMinus10Perc_TTbarJets_SemiMuon_"+decayChannel+".root").c_str());
      nameDataSet.push_back("UnclusEnergyMinus_TTbarJets_SemiMuon");
    }
    //9) Top Mass
    if(TopMassResults == true){
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_TTbarjets_TopMass175_5_"+decayChannel+".root").c_str());
      nameDataSet.push_back("TopMassPlus");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_TTbarjets_TopMass169_5_"+decayChannel+".root").c_str());
      nameDataSet.push_back("TopMassMinus");
    }

    //10) Trigger & Event selection
    if(TriggEvtSelResults == true){
      //Consider TriggEvtSel Plus samples:
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ST_SingleTop_tChannel_tbar_"+decayChannel+".root").c_str());
      nameDataSet.push_back("TriggEvtSelPlus_ST_SingleTop_tChannel_tbar");
      //      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_SingleTop_tChannel_t_"+decayChannel+".root").c_str());
      //      nameDataSet.push_back("TriggEvtSelPlus_ST_SingleTop_tChannel_t");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ST_SingleTop_tWChannel_tbar_"+decayChannel+".root").c_str());
      nameDataSet.push_back("TriggEvtSelPlus_ST_SingleTop_tWChannel_tbar");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ST_SingleTop_tWChannel_t_"+decayChannel+".root").c_str());
      nameDataSet.push_back("TriggEvtSelPlus_ST_SingleTop_tWChannel_t");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_TTbarJets_Other_"+decayChannel+".root").c_str());
      nameDataSet.push_back("TriggEvtSelPlus_TTbarJets_Other");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_TTbarJets_SemiEl_"+decayChannel+".root").c_str());
      nameDataSet.push_back("TriggEvtSelPlus_TTbarJets_SemiEl");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_WJets_"+decayChannel+".root").c_str());    
      nameDataSet.push_back("TriggEvtSelPlus_WJets");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ZJets_"+decayChannel+".root").c_str());
      nameDataSet.push_back("TriggEvtSelPlus_ZJets");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_TTbarJets_SemiMuon_"+decayChannel+".root").c_str());
      nameDataSet.push_back("TriggEvtSelPlus_TTbarJets_SemiMuon");

      //Consider TriggEvtSel Minus samples:
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ST_SingleTop_tChannel_tbar_"+decayChannel+".root").c_str());
      nameDataSet.push_back("TriggEvtSelMinus_ST_SingleTop_tChannel_tbar");
      //      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ST_SingleTop_tChannel_t_"+decayChannel+".root").c_str());
      //      nameDataSet.push_back("TriggEvtSel_ST_SingleTop_tChannel_t");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ST_SingleTop_tWChannel_tbar_"+decayChannel+".root").c_str());
      nameDataSet.push_back("TriggEvtSelMinus_ST_SingleTop_tWChannel_tbar");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ST_SingleTop_tWChannel_t_"+decayChannel+".root").c_str());
      nameDataSet.push_back("TriggEvtSelMinus_ST_SingleTop_tWChannel_t");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_TTbarJets_Other_"+decayChannel+".root").c_str());
      nameDataSet.push_back("TriggEvtSelMinus_TTbarJets_Other");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_TTbarJets_SemiEl_"+decayChannel+".root").c_str());
      nameDataSet.push_back("TriggEvtSelMinus_TTbarJets_SemiEl");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_WJets_"+decayChannel+".root").c_str());    
      nameDataSet.push_back("TriggEvtSelMinus_WJets");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ZJets_"+decayChannel+".root").c_str());
      nameDataSet.push_back("TriggEvtSelMinus_ZJets");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_TTbarJets_SemiMuon_"+decayChannel+".root").c_str());
      nameDataSet.push_back("TriggEvtSelMinus_TTbarJets_SemiMuon");      
    }

  }//End of signalOnly = false loop
  
  //TTbarJets_SemiMuon sample should always be put as latest sample to avoid crash of TMinuitMinimizer !!
  if(semiMuon == true){
    inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_TTbarJets_SemiMuon_"+decayChannel+".root").c_str());
    nameDataSet.push_back("TTbarJets_SemiMuon");
  }
  else if(semiElectron == true){
    inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_TTbarJets_SemiEl_"+decayChannel+".root").c_str());
    nameDataSet.push_back("TTbarJets_SemiEl");
  }
  
  TFile *fout = new TFile (("WHelicities_Analysis_"+UsedTrigger+".root").c_str(), "RECREATE");
  fout->cd();

  //-----------------------------------------------
  //  Output .tex files for presentation/paper
  //-----------------------------------------------
  //mkdir("BTagPerformanceStudy/",0777);

  string bTagPerf = "BTagPerformanceStudy/";
  mkdir(("HadrKinFit/"+bTagPerf).c_str(),0777);
  mkdir(("HadrAndLeptWOnlyKinFit/"+bTagPerf).c_str(),0777);
  mkdir(("HadrAndLeptKinFit/"+bTagPerf).c_str(),0777);
 
  //oooooOOOOOOOOOOOoooooooooooOOOOOOOOOOOOOOOOoooooooooOOOOOOOOOOOoooooooooooooooooOOOOOOOOOOOOOOOoooooooooooooOOOOOOOOOOOoooooooo
  //1) Own results: No btag at event selection, top mass constraint for leptonic top in KinFit and offline muon pt cut of 27
  //oooooOOOOOOOOOOOoooooooooooOOOOOOOOOOOOOOOOoooooooooOOOOOOOOOOOoooooooooooooooooOOOOOOOOOOOOOOOoooooooooooooOOOOOOOOOOOoooooooo

  string PresentationTexTitle;
  PresentationTexTitle =("BTagPerformanceStudy/"+UsedTrigger+"_"+decayChannel+".tex").c_str();
 
  ofstream PresentationTexHadr(("HadrKinFit/"+PresentationTexTitle).c_str());
  ofstream PresentationTexHadrAndLeptWOnly(("HadrAndLeptWOnlyKinFit/"+PresentationTexTitle).c_str());
  ofstream PresentationTexHadrAndLept(("HadrAndLeptKinFit/"+PresentationTexTitle).c_str());

  ofstream FMinusTexHadr("HadrKinFit/FMinus.tex");
  ofstream FPlusTexHadr("HadrKinFit/FPlus.tex");
  ofstream FZeroTexHadr("HadrKinFit/FZero.tex");

  ofstream FMinusTexHadrAndLeptWOnly("HadrAndLeptWOnlyKinFit/FMinus.tex");
  ofstream FPlusTexHadrAndLeptWOnly("HadrAndLeptWOnlyKinFit/FPlus.tex");
  ofstream FZeroTexHadrAndLeptWOnly("HadrAndLeptWOnlyKinFit/FZero.tex");

  ofstream FMinusTexHadrAndLept("HadrAndLeptKinFit/FMinus.tex");
  ofstream FPlusTexHadrAndLept("HadrAndLeptKinFit/FPlus.tex");
  ofstream FZeroTexHadrAndLept("HadrAndLeptKinFit/FZero.tex");

  //-----------------------------------------
  //  Start of filling of .tex files !!
  //-----------------------------------------
  //Hadronic KinFit configuration: 
  PresentationTexHadr << " \\begin{table} " << endl;
  PresentationTexHadr << " \\begin{tiny} " << endl;
  PresentationTexHadr << " \\renewcommand{\\arraystretch}{1.2} " << endl;
  PresentationTexHadr << " \\begin{center} " << endl;
  for(int ii=0; ii<nameDataSet.size(); ii++){
    if(ii==0){PresentationTexHadr << "  \\begin{tabular}{|c|";}
    else if(ii < (nameDataSet.size()-1)){PresentationTexHadr << "c|";}
    else{PresentationTexHadr << "c|c|c|c|c|c|c|c|c|c|c|} " << endl;}
  }
  PresentationTexHadr << " \\hline " << endl;
  for(int ii=0; ii < nameDataSet.size();ii++){
    PresentationTexHadr << " & " << nameDataSet[ii];
  }
  PresentationTexHadr << " & bLept correct & F+ & F+ - SM & $\\delta$ F+ & F- & F- - SM & $\\delta$ F- & F0 & F0 - SM & $\\delta$ F0 \\\\ " << endl;
  PresentationTexHadr << " \\hline " << endl;

  //Hadronic and Leptonic W Only KinFit configuration: 
  PresentationTexHadrAndLeptWOnly << " \\begin{table} " << endl;
  PresentationTexHadrAndLeptWOnly << " \\begin{tiny} " << endl;
  PresentationTexHadrAndLeptWOnly << " \\renewcommand{\\arraystretch}{1.2} " << endl;
  PresentationTexHadrAndLeptWOnly << " \\begin{center} " << endl;
  for(int ii=0; ii<nameDataSet.size(); ii++){
    if(ii==0){PresentationTexHadrAndLeptWOnly << "  \\begin{tabular}{|c|";}
    else if(ii < (nameDataSet.size()-1)){PresentationTexHadrAndLeptWOnly << "c|";}
    else{PresentationTexHadrAndLeptWOnly << "c|c|c|c|c|c|c|c|c|c|c|} " << endl;}
  }
  PresentationTexHadrAndLeptWOnly << " \\hline " << endl;
  for(int ii=0; ii < nameDataSet.size();ii++){
    PresentationTexHadrAndLeptWOnly << " & " << nameDataSet[ii];
  }
  PresentationTexHadrAndLeptWOnly << " & bLept correct & F+ & F+ - SM & $\\delta$ F+ & F- & F- - SM & $\\delta$ F- & F0 & F0 - SM & $\\delta$ F0 \\\\ " << endl;
  PresentationTexHadrAndLeptWOnly << " \\hline " << endl;

  //Hadronic and Leptonic KinFit configuration: 
  PresentationTexHadrAndLept << " \\begin{table} " << endl;
  PresentationTexHadrAndLept << " \\begin{tiny} " << endl;
  PresentationTexHadrAndLept << " \\renewcommand{\\arraystretch}{1.2} " << endl;
  PresentationTexHadrAndLept << " \\begin{center} " << endl;
  for(int ii=0; ii<nameDataSet.size(); ii++){
    if(ii==0){PresentationTexHadrAndLept << "  \\begin{tabular}{|c|";}
    else if(ii < (nameDataSet.size()-1)){PresentationTexHadrAndLept << "c|";}
    else{PresentationTexHadrAndLept << "c|c|c|c|c|c|c|c|c|c|c|} " << endl;}
  }
  PresentationTexHadrAndLept << " \\hline " << endl;
  for(int ii=0; ii < nameDataSet.size();ii++){
    PresentationTexHadrAndLept << " & " << nameDataSet[ii];
  }
  PresentationTexHadrAndLept << " & bLept correct & F+ & F+ - SM & $\\delta$ F+ & F- & F- - SM & $\\delta$ F- & F0 & F0 - SM & $\\delta$ F0 \\\\ " << endl;
  PresentationTexHadrAndLept << " \\hline " << endl;

  // initialize histograms
  cout << "Initializing histograms" << endl;

  histo1D["lumiWeights3D"]= new TH1F("lumiWeights3D","lumiWeights3D",25,0,50);

  Int_t N;
  Double_t *xVec, *yVec, *xErrHi, *xErrLo, *yErrHi, *yErrLo;
  if(semiElectron == true){
    TFile *ElecFile_ = new TFile("JetLegTriggerEfficiencyIsoLepTriJetJetMult4.root","read");
    TGraph* effHists_= (TGraph*) ElecFile_->Get("JetLegTriggerEfficiencyIsoLepTriJetJetMult4")->Clone();
    /// read values from TGraphAsymmerrorSFeles
    N = effHists_->GetN();
    xVec = effHists_->GetX();
    yVec = effHists_->GetY();
    xErrHi = effHists_->GetEXhigh();
    xErrLo = effHists_->GetEXlow();
    yErrHi = effHists_->GetEYhigh();
    yErrLo = effHists_->GetEYlow();
  }
  else if(semiMuon == true){
    TFile *MuonFile_ = new TFile("MuonEffSF2011.root","read");
    TGraph* effHists_  = (TGraph*) MuonFile_->Get("tapAllSFeta")->Clone();
    /// read values from TGraphAsymmErrors
    N      = effHists_->GetN();
    xVec   = effHists_->GetX();
    yVec   = effHists_->GetY();
    xErrHi = effHists_->GetEXhigh();
    xErrLo = effHists_->GetEXlow();
    yErrHi = effHists_->GetEYhigh();
    yErrLo = effHists_->GetEYlow(); 
  } 

  //Standard Model helicity values:
  float SMfrResult = 0.0334506;
  float SMflResult = 0.321241;
  float SMf0Result = 0.64491;
  histo1D["StandardCosThetaFit"]=new TH1F("StCosThetaFit","StCosThetaFit",200,-1,1);   
  histo1D["StandardCosTheta"]=new TH1F("StCosTheta","StCosTheta",200,-1,1);   

  //Zie code Stijn voor alle gebruikte controle plots !

  histo1D["JetCombHadr"] = new TH1F("JetCombHadr","JetCombHadr",14,-0.5,13.5);
  histo1D["JetCombHadrAndLeptWOnly"] = new TH1F("JetCombHadrAndLeptWOnly","JetCombHadrAndLeptWOnly",14,-0.5,13.5);
  histo1D["JetCombHadrAndLept"] = new TH1F("JetCombHaAndLeptdr","JetCombHadrAndLept",14,-0.5,13.5);

  histo1D["LeptWMassHadr"] = new TH1F("LeptWMassHadr","LeptWMassHadr",40,0,150);
  histo1D["LeptTopMassHadr"] = new TH1F("LeptTopMassHadr","LeptTopMassHadr",60,50,350);

  //   //Histograms to obtain ratio of cos theta* for alternative and SM helicities
  //   float HelicityWeight[3];
  //   int SizeArray = 3;
  //   float HelicityFraction[SizeArray][3];  //0:Longitudinal; 1:Righthanded; 2:Lefthanded
  //   float UsedDistributionValue[SizeArray];
  //   std::ostringstream HelicityNumbers[SizeArray][3];
  //   for(int helicityNumbers=0;helicityNumbers<SizeArray;helicityNumbers++){
  //     HelicityWeight[helicityNumbers]=1;
  //     UsedDistributionValue[helicityNumbers]=0;
  //     if(helicityNumbers == 0){
  //       HelicityFraction[helicityNumbers][0]=0.5;
  //       HelicityFraction[helicityNumbers][1]=0.;
  //       HelicityFraction[helicityNumbers][2]=0.5;
  //     }
  //     else if(helicityNumbers==1){
  //       HelicityFraction[helicityNumbers][0]=0.6;
  //       HelicityFraction[helicityNumbers][1]=0.2;
  //       HelicityFraction[helicityNumbers][2]=0.2;
  //     }
  //     else if(helicityNumbers==2){
  //       HelicityFraction[helicityNumbers][0]=0.65;
  //       HelicityFraction[helicityNumbers][1]=0.1;
  //       HelicityFraction[helicityNumbers][2]=0.25;
  //     }
  //     HelicityNumbers[helicityNumbers][0] << HelicityFraction[helicityNumbers][0];
  //     HelicityNumbers[helicityNumbers][1] << HelicityFraction[helicityNumbers][1];
  //     HelicityNumbers[helicityNumbers][2] << HelicityFraction[helicityNumbers][2];
    
  //     std::string HistoName = "CosThetaRight"+HelicityNumbers[helicityNumbers][1].str()+"Long"+HelicityNumbers[helicityNumbers][0].str()+"Left"+HelicityNumbers[helicityNumbers][2].str();
  //     TString THistoName = "CosThetaRight"+HelicityNumbers[helicityNumbers][1].str()+"Long"+HelicityNumbers[helicityNumbers][0].str()+"Left"+HelicityNumbers[helicityNumbers][2].str();
  //     std::string MSSVHEbTagHistoName = "MSSVHEbTagCosThetaRight"+HelicityNumbers[helicityNumbers][1].str()+"Long"+HelicityNumbers[helicityNumbers][0].str()+"Left"+HelicityNumbers[helicityNumbers][2].str();
  //     TString MSSVHEbTagTHistoName = "MSSVHEbTagCosThetaRight"+HelicityNumbers[helicityNumbers][1].str()+"Long"+HelicityNumbers[helicityNumbers][0].str()+"Left"+HelicityNumbers[helicityNumbers][2].str();
  //     std::string TCSVbTagHistoName = "TCSVbTagCosThetaRight"+HelicityNumbers[helicityNumbers][1].str()+"Long"+HelicityNumbers[helicityNumbers][0].str()+"Left"+HelicityNumbers[helicityNumbers][2].str();
  //     TString TCSVbTagTHistoName = "TCSVbTagCosThetaRight"+HelicityNumbers[helicityNumbers][1].str()+"Long"+HelicityNumbers[helicityNumbers][0].str()+"Left"+HelicityNumbers[helicityNumbers][2].str();

  //     histo1D[HistoName] = new TH1F(THistoName,THistoName,200,-1,1);               
  //     histo1D[MSSVHEbTagHistoName] = new TH1F(MSSVHEbTagTHistoName,MSSVHEbTagTHistoName,200,-1,1);               
  //     histo1D[TCSVbTagHistoName] = new TH1F(TCSVbTagTHistoName,TCSVbTagTHistoName,200,-1,1);               
  //   }

  float XSection, EqLumi, NominalNormFactor;
  vector<Dataset*> datasets; // needed for MSPlots
  string dataSetName;
  for(unsigned int iDataSet=0; iDataSet<inputWTree.size(); iDataSet++){
    TFile* inFile = new TFile(inputWTree[iDataSet].c_str(),"READ");
    TTree* inConfigTree = (TTree*) inFile->Get("configTreeWTreeFile");
    TBranch* d_br = (TBranch*) inConfigTree->GetBranch("Dataset");
    TClonesArray* tc_dataset = new TClonesArray("Dataset",0);
    d_br->SetAddress(&tc_dataset);
    inConfigTree->GetEvent(0);
    Dataset* dataSet = (Dataset*) tc_dataset->At(0);
    int color = 0;
    XSection = dataSet->Xsection();
    EqLumi = dataSet->EquivalentLumi();    
    NominalNormFactor = dataSet->NormFactor();
    dataSetName = nameDataSet[iDataSet];
 
    if( dataSetName.find("WSystMinus_WJets") == 0 && WSystResults == true){
      dataSet->SetEquivalentLuminosity( dataSet->EquivalentLumi() / (0.7) );  //WJets Minus 30%
      //dataSet->SetEquivalentLuminosity( dataSet->EquivalentLumi() / (0.0000001) );  //WJets Minus 100%
    }
    if(dataSetName.find("WSystPlus_WJets") == 0 && WSystResults == true){
      dataSet->SetEquivalentLuminosity( dataSet->EquivalentLumi() / (1.3) ); //WJets Plus 30 %
      //dataSet->SetEquivalentLuminosity( dataSet->EquivalentLumi() / (2.) ); //WJets Plus 100 %
    }

    if( dataSetName.find("Data") == 0 && semiElectron == true && TriCentralJet30Trigger == true){
      dataSet->SetEquivalentLuminosity( Luminosity);
      dataSet->SetTitle("Data (4.966 fb^{-1}");
    }
   
    //Nominal samples:
    if( dataSet->Name().find("TTbarJets_SemiMuon") == 0 && dataSet->Name().find("JES") != 0) color = kRed+1;
    if( dataSet->Name().find("TTbarJets_SemiEl") == 0 && dataSet->Name().find("JES") != 0) color = kRed-4;
    if( dataSet->Name().find("TTbarJets_Other") == 0 && dataSet->Name().find("JES") != 0 ) color = kRed-7;
    if( dataSet->Name().find("WJets") == 0 ){
      dataSet->SetTitle("W#rightarrowl#nu");
      color = kGreen-3;
    }
    if( dataSet->Name().find("ZJets") == 0 && dataSet->Name().find("JES") != 0 ){
      dataSet->SetTitle("Z/#gamma*#rightarrowl^{+}l^{-}");
      color = kAzure-2;
    }
    if( dataSet->Name().find("ST") == 0 && dataSet->Name().find("JES") != 0 ) color = kMagenta;
    if( dataSet->Name().find("QCD") == 0 && dataSet->Name().find("JES") !=0 ) color = kBlue;        
      
    Dataset* tmpDS = new Dataset(dataSet->Name(), dataSet->Title(), dataSet->DoIt(), color, dataSet->LineStyle(), dataSet->LineWidth(), dataSet->NormFactor(), dataSet->Xsection());
    tmpDS->SetEquivalentLuminosity( dataSet->EquivalentLumi() );
    datasets.push_back( tmpDS );
  }

  ////////////////////////////////////////////////////////////////////
  //  MSPlots for case with b-tag after jet combination selection:  //
  ////////////////////////////////////////////////////////////////////
  if(CreateMSPlots == true){
    MSPlot["JetCombinationChiSqOnly"]= new MultiSamplePlot(datasets,"JetCombinationChiSqOnly",14,-0.5,13.5,"JetCombinationChiSqOnly");
    //MSPlot["PtLeptonicBOriginalChiSqOnly"]= new MultiSamplePlot(datasets,"PtLeptonicBOriginalChiSqOnly",20,0,100,"PtLeptonicBOriginalChiSqOnly");
    //MSPlot["PtHadronicBOriginalChiSqOnly"]= new MultiSamplePlot(datasets,"PtHadronicBOriginalChiSqOnly",20,0,100,"PtHadronicBOriginalChiSqOnly");
    //MSPlot["PtLight1JetOriginalChiSqOnly"]= new MultiSamplePlot(datasets,"PtLight1JetOriginalChiSqOnly",20,0,100,"PtLight1JetOriginalChiSqOnly");
    //MSPlot["PtLight2JetOriginalChiSqOnly"]= new MultiSamplePlot(datasets,"PtLight2JetOriginalChiSqOnly",20,0,100,"PtLight2JetOriginalChiSqOnly");

    //No b-tag applied
    MSPlot["JetCombCosThetaFoundChiSqOnlyNoBTag"]= new MultiSamplePlot(datasets,"JetCombCosThetaFoundChiSqOnlyNoBTag",13,-0.5,12.5,"JetCombCosThetaFoundChiSqOnlyNoBTag");
    MSPlot["ChiSqKinFitCosThetaFoundChiSqOnlyNoBTag"]= new MultiSamplePlot(datasets,"ChiSqKinFitCosThetaFoundChiSqOnlyNoBTag",40,0,150,"ChiSqKinFitCosThetaFoundChiSqOnlyNoBTag");
    MSPlot["ProbKinFitCosThetaFoundChiSqOnlyNoBTag"]= new MultiSamplePlot(datasets,"ProbKinFitCosThetaFoundChiSqOnlyNoBTag",20,0,1,"ProbKinFitCosThetaFoundChiSqOnlyNoBTag");

    MSPlot["CosThetaChiSqOnlyNoBTag"]= new MultiSamplePlot(datasets,"CosThetaChiSqOnlyNoBTag",CosThetaBinNumber,-1,1,"CosThetaChiSqOnlyNoBTag");
    MSPlot["nPVChiSqOnlyNoBTag"]= new MultiSamplePlot(datasets,"nPVChiSqOnlyNoBTag",25,0,25,"nPVChiSqOnlyNoBTag");
    MSPlot["TransverseMassChiSqOnlyNoBTag"]= new MultiSamplePlot(datasets,"TransverseMassChiSqOnlyNoBTag",30,0,200,"TransverseMassChiSqOnlyNoBTag");
    MSPlot["TransverseMassSelectedParticlesChiSqOnlyNoBTag"]= new MultiSamplePlot(datasets,"TransverseMassSelectedParticlesChiSqOnlyNoBTag",30,0,200,"TransverseMassSelectedParticlesChiSqOnlyNoBTag");
  
    MSPlot["Jet1PtChiSqOnlyNoBTag"]= new MultiSamplePlot(datasets,"Jet1PtChiSqOnlyNoBTag",100,0,500,"Jet1PtChiSqOnlyNoBTag");
    MSPlot["Jet2PtChiSqOnlyNoBTag"]= new MultiSamplePlot(datasets,"Jet2PtChiSqOnlyNoBTag",100,0,500,"Jet2PtChiSqOnlyNoBTag");
    MSPlot["Jet3PtChiSqOnlyNoBTag"]= new MultiSamplePlot(datasets,"Jet3PtChiSqOnlyNoBTag",100,0,500,"Jet3PtChiSqOnlyNoBTag");
    MSPlot["Jet4PtChiSqOnlyNoBTag"]= new MultiSamplePlot(datasets,"Jet4PtChiSqOnlyNoBTag",100,0,500,"Jet4PtChiSqOnlyNoBTag");
    
    MSPlot["PtLeptonicBChiSqOnlyNoBTag"]= new MultiSamplePlot(datasets,"PtLeptonicBChiSqOnlyNoBTag",100,0,500,"PtLeptonicBChiSqOnlyNoBTag");
    MSPlot["PtHadronicBChiSqOnlyNoBTag"]= new MultiSamplePlot(datasets,"PtHadronicBChiSqOnlyNoBTag",100,0,500,"PtHadronicBChiSqOnlyNoBTag");
    MSPlot["PtLightJet1ChiSqOnlyNoBTag"]= new MultiSamplePlot(datasets,"PtLightJet1ChiSqOnlyNoBTag",100,0,500,"PtLightJet1ChiSqOnlyNoBTag");
    MSPlot["PtLightJet2ChiSqOnlyNoBTag"]= new MultiSamplePlot(datasets,"PtLightJet2ChiSqOnlyNoBTag",100,0,500,"PtLightJet2ChiSqOnlyNoBTag");
    
    MSPlot["leptonPhiChiSqOnlyNoBTag"]= new MultiSamplePlot(datasets,"leptonPhiChiSqOnlyNoBTag",30,-4,4,"leptonPhiChiSqOnlyNoBTag");
    MSPlot["leptonPxChiSqOnlyNoBTag"]= new MultiSamplePlot(datasets,"leptonPxChiSqOnlyNoBTag",50,-200,200,"leptonPxChiSqOnlyNoBTag");
    MSPlot["leptonPyChiSqOnlyNoBTag"]= new MultiSamplePlot(datasets,"leptonPyChiSqOnlyNoBTag",50,-200,200,"leptonPyChiSqOnlyNoBTag");
    MSPlot["leptonPzChiSqOnlyNoBTag"]= new MultiSamplePlot(datasets,"leptonPzChiSqOnlyNoBTag",50,-500,500,"leptonPzChiSqOnlyNoBTag");
    MSPlot["leptonPtChiSqOnlyNoBTag"]= new MultiSamplePlot(datasets,"leptonPtChiSqOnlyNoBTag",50,0,500,"leptonPtChiSqOnlyNoBTag");
    MSPlot["leptonMassChiSqOnlyNoBTag"]= new MultiSamplePlot(datasets,"leptonMassChiSqOnlyNoBTag",30,0.1,0.15,"leptonMassChiSqOnlyNoBTag");
    
    MSPlot["METPhiChiSqOnlyNoBTag"]= new MultiSamplePlot(datasets,"METPhiChiSqOnlyNoBTag",30,-4,4,"METPhiChiSqOnlyNoBTag");
    MSPlot["METPxChiSqOnlyNoBTag"]= new MultiSamplePlot(datasets,"METPxChiSqOnlyNoBTag",50,-200,200,"METPxChiSqOnlyNoBTag");
    MSPlot["METPyChiSqOnlyNoBTag"]= new MultiSamplePlot(datasets,"METPyChiSqOnlyNoBTag",50,-200,200,"METPyChiSqOnlyNoBTag");
    MSPlot["METPzChiSqOnlyNoBTag"]= new MultiSamplePlot(datasets,"METPzChiSqOnlyNoBTag",50,-500,500,"METPzChiSqOnlyNoBTag");
    MSPlot["METPtChiSqOnlyNoBTag"]= new MultiSamplePlot(datasets,"METPtChiSqOnlyNoBTag",100,0,500,"METPtChiSqOnlyNoBTag");
    MSPlot["METMassChiSqOnlyNoBTag"]= new MultiSamplePlot(datasets,"METMassChiSqOnlyNoBTag",30,0,0.01,"METMassChiSqOnlyNoBTag");
    
    //Specific b-tag applied
    MSPlot["JetCombCosThetaFoundChiSqOnlySpecificBTag"]= new MultiSamplePlot(datasets,"JetCombCosThetaFoundChiSqOnlySpecificBTag",13,-0.5,12.5,"JetCombCosThetaFoundChiSqOnlySpecificBTag");
    MSPlot["ChiSqKinFitCosThetaFoundChiSqOnlySpecificBTag"]= new MultiSamplePlot(datasets,"ChiSqKinFitCosThetaFoundChiSqOnlySpecificBTag",40,0,150,"ChiSqKinFitCosThetaFoundChiSqOnlySpecificBTag");
    MSPlot["ProbKinFitCosThetaFoundChiSqOnlySpecificBTag"]= new MultiSamplePlot(datasets,"ProbKinFitCosThetaFoundChiSqOnlySpecificBTag",20,0,1,"ProbKinFitCosThetaFoundChiSqOnlySpecificBTag");
    
    MSPlot["CosThetaChiSqOnlySpecificBTag"]= new MultiSamplePlot(datasets,"CosThetaChiSqOnlySpecificBTag",CosThetaBinNumber,-1,1,"CosThetaChiSqOnlySpecificBTag");
    MSPlot["nPVChiSqOnlySpecificBTag"]= new MultiSamplePlot(datasets,"nPVChiSqOnlySpecificBTag",25,0,25,"nPVChiSqOnlySpecificBTag");
    MSPlot["TransverseMassChiSqOnlySpecificBTag"]= new MultiSamplePlot(datasets,"TransverseMassChiSqOnlySpecificBTag",30,0,200,"TransverseMassChiSqOnlySpecificBTag");
    MSPlot["TransverseMassSelectedParticlesChiSqOnlySpecificBTag"]= new MultiSamplePlot(datasets,"TransverseMassSelectedParticlesChiSqOnlySpecificBTag",30,0,200,"TransverseMassSelectedParticlesChiSqOnlySpecificBTag");
    
    MSPlot["Jet1PtChiSqOnlySpecificBTag"]= new MultiSamplePlot(datasets,"Jet1PtChiSqOnlySpecificBTag",100,0,500,"Jet1PtChiSqOnlySpecificBTag");
    MSPlot["Jet2PtChiSqOnlySpecificBTag"]= new MultiSamplePlot(datasets,"Jet2PtChiSqOnlySpecificBTag",100,0,500,"Jet2PtChiSqOnlySpecificBTag");
    MSPlot["Jet3PtChiSqOnlySpecificBTag"]= new MultiSamplePlot(datasets,"Jet3PtChiSqOnlySpecificBTag",100,0,500,"Jet3PtChiSqOnlySpecificBTag");
    MSPlot["Jet4PtChiSqOnlySpecificBTag"]= new MultiSamplePlot(datasets,"Jet4PtChiSqOnlySpecificBTag",100,0,500,"Jet4PtChiSqOnlySpecificBTag");
    
    MSPlot["PtLeptonicBChiSqOnlySpecificBTag"]= new MultiSamplePlot(datasets,"PtLeptonicBChiSqOnlySpecificBTag",100,0,500,"PtLeptonicBChiSqOnlySpecificBTag");
    MSPlot["PtHadronicBChiSqOnlySpecificBTag"]= new MultiSamplePlot(datasets,"PtHadronicBChiSqOnlySpecificBTag",100,0,500,"PtHadronicBChiSqOnlySpecificBTag");
    MSPlot["PtLightJet1ChiSqOnlySpecificBTag"]= new MultiSamplePlot(datasets,"PtLightJet1ChiSqOnlySpecificBTag",100,0,500,"PtLightJet1ChiSqOnlySpecificBTag");
    MSPlot["PtLightJet2ChiSqOnlySpecificBTag"]= new MultiSamplePlot(datasets,"PtLightJet2ChiSqOnlySpecificBTag",100,0,500,"PtLightJet2ChiSqOnlySpecificBTag");
    
    MSPlot["leptonPhiChiSqOnlySpecificBTag"]= new MultiSamplePlot(datasets,"leptonPhiChiSqOnlySpecificBTag",30,-4,4,"leptonPhiChiSqOnlySpecificBTag");
    MSPlot["leptonPxChiSqOnlySpecificBTag"]= new MultiSamplePlot(datasets,"leptonPxChiSqOnlySpecificBTag",50,-200,200,"leptonPxChiSqOnlySpecificBTag");
    MSPlot["leptonPyChiSqOnlySpecificBTag"]= new MultiSamplePlot(datasets,"leptonPyChiSqOnlySpecificBTag",50,-200,200,"leptonPyChiSqOnlySpecificBTag");
    MSPlot["leptonPzChiSqOnlySpecificBTag"]= new MultiSamplePlot(datasets,"leptonPzChiSqOnlySpecificBTag",50,-500,500,"leptonPzChiSqOnlySpecificBTag");
    MSPlot["leptonPtChiSqOnlySpecificBTag"]= new MultiSamplePlot(datasets,"leptonPtChiSqOnlySpecificBTag",100,0,500,"leptonPtChiSqOnlySpecificBTag");
    MSPlot["leptonMassChiSqOnlySpecificBTag"]= new MultiSamplePlot(datasets,"leptonMassChiSqOnlySpecificBTag",30,0.1,0.15,"leptonMassChiSqOnlySpecificBTag");
    
    MSPlot["METPhiChiSqOnlySpecificBTag"]= new MultiSamplePlot(datasets,"METPhiChiSqOnlySpecificBTag",30,-4,4,"METPhiChiSqOnlySpecificBTag");
    MSPlot["METPxChiSqOnlySpecificBTag"]= new MultiSamplePlot(datasets,"METPxChiSqOnlySpecificBTag",50,-200,200,"METPxChiSqOnlySpecificBTag");
    MSPlot["METPyChiSqOnlySpecificBTag"]= new MultiSamplePlot(datasets,"METPyChiSqOnlySpecificBTag",50,-200,200,"METPyChiSqOnlySpecificBTag");
    MSPlot["METPzChiSqOnlySpecificBTag"]= new MultiSamplePlot(datasets,"METPzChiSqOnlySpecificBTag",50,-500,500,"METPzChiSqOnlySpecificBTag");
    MSPlot["METPtChiSqOnlySpecificBTag"]= new MultiSamplePlot(datasets,"METPtChiSqOnlySpecificBTag",100,0,500,"METPtChiSqOnlySpecificBTag");
    MSPlot["METMassChiSqOnlySpecificBTag"]= new MultiSamplePlot(datasets,"METMassChiSqOnlySpecificBTag",30,0,0.01,"METMassChiSqOnlySpecificBTag");
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    ///////////////////////////////////////////////////////////////////////////////
    //  MSPlots for jet combination together with b-tag (Hadr and Lept KinFit):  //
    ///////////////////////////////////////////////////////////////////////////////
    MSPlot["JetCombination"]= new MultiSamplePlot(datasets,"JetCombination",14,-0.5,13.5,"JetCombination");
    
    //No b-tag applied:
    MSPlot["JetCombCosThetaFoundNoBTag"]= new MultiSamplePlot(datasets,"JetCombCosThetaFoundNoBTag",13,-0.5,12.5,"JetCombCosThetaFoundNoBTag");
    MSPlot["ChiSqKinFitCosThetaFoundNoBTag"]= new MultiSamplePlot(datasets,"ChiSqKinFitCosThetaFoundNoBTag",40,0,150,"ChiSqKinFitCosThetaFoundNoBTag");
    MSPlot["ProbKinFitCosThetaFoundNoBTag"]= new MultiSamplePlot(datasets,"ProbKinFitCosThetaFoundNoBTag",20,0,1,"ProbKinFitCosThetaFoundNoBTag");
    
    MSPlot["CosThetaNoBTag"]= new MultiSamplePlot(datasets,"CosThetaNoBTag",CosThetaBinNumber,-1,1,"CosThetaNoBTag");
    MSPlot["nPVNoBTag"]= new MultiSamplePlot(datasets,"nPVNoBTag",25,0,25,"nPVNoBTag");
    MSPlot["TransverseMassNoBTag"]= new MultiSamplePlot(datasets,"TransverseMassNoBTag",30,0,200,"TransverseMassNoBTag");
    MSPlot["TransverseMassSelectedParticlesNoBTag"]= new MultiSamplePlot(datasets,"TransverseMassSelectedParticlesNoBTag",30,0,200,"TransverseMassSelectedParticlesNoBTag");
    
    MSPlot["Jet1PtNoBTag"]= new MultiSamplePlot(datasets,"Jet1PtNoBTag",100,0,500,"Jet1PtNoBTag");
    MSPlot["Jet2PtNoBTag"]= new MultiSamplePlot(datasets,"Jet2PtNoBTag",100,0,500,"Jet2PtNoBTag");
    MSPlot["Jet3PtNoBTag"]= new MultiSamplePlot(datasets,"Jet3PtNoBTag",100,0,500,"Jet3PtNoBTag");
    MSPlot["Jet4PtNoBTag"]= new MultiSamplePlot(datasets,"Jet4PtNoBTag",100,0,500,"Jet4PtNoBTag");
    
    MSPlot["PtLeptonicBNoBTag"]= new MultiSamplePlot(datasets,"PtLeptonicBNoBTag",100,0,500,"PtLeptonicBNoBTag");
    MSPlot["PtHadronicBNoBTag"]= new MultiSamplePlot(datasets,"PtHadronicBNoBTag",100,0,500,"PtHadronicBNoBTag");
    MSPlot["PtLightJet1NoBTag"]= new MultiSamplePlot(datasets,"PtLightJet1NoBTag",100,0,500,"PtLightJet1NoBTag");
    MSPlot["PtLightJet2NoBTag"]= new MultiSamplePlot(datasets,"PtLightJet2NoBTag",100,0,500,"PtLightJet2NoBTag");
    
    MSPlot["leptonPhiNoBTag"]= new MultiSamplePlot(datasets,"leptonPhiNoBTag",30,-4,4,"leptonPhiNoBTag");
    MSPlot["leptonPxNoBTag"]= new MultiSamplePlot(datasets,"leptonPxNoBTag",50,-200,200,"leptonPxNoBTag");
    MSPlot["leptonPyNoBTag"]= new MultiSamplePlot(datasets,"leptonPyNoBTag",50,-200,200,"leptonPyNoBTag");
    MSPlot["leptonPzNoBTag"]= new MultiSamplePlot(datasets,"leptonPzNoBTag",50,-500,500,"leptonPzNoBTag");
    MSPlot["leptonPtNoBTag"]= new MultiSamplePlot(datasets,"leptonPtNoBTag",100,0,500,"leptonPtNoBTag");
    MSPlot["leptonMassNoBTag"]= new MultiSamplePlot(datasets,"leptonMassNoBTag",30,0.1,0.15,"leptonMassNoBTag");
    
    MSPlot["METPhiNoBTag"]= new MultiSamplePlot(datasets,"METPhiNoBTag",30,-4,4,"METPhiNoBTag");
    MSPlot["METPxNoBTag"]= new MultiSamplePlot(datasets,"METPxNoBTag",50,-200,200,"METPxNoBTag");
    MSPlot["METPyNoBTag"]= new MultiSamplePlot(datasets,"METPyNoBTag",50,-200,200,"METPyNoBTag");
    MSPlot["METPzNoBTag"]= new MultiSamplePlot(datasets,"METPzNoBTag",50,-500,500,"METPzNoBTag");
    MSPlot["METPtNoBTag"]= new MultiSamplePlot(datasets,"METPtNoBTag",100,0,500,"METPtNoBTag");
    MSPlot["METMassNoBTag"]= new MultiSamplePlot(datasets,"METMassNoBTag",30,0,0.01,"METMassNoBTag");

    MSPlot["HadronicWMassNoBTag"]= new MultiSamplePlot(datasets,"HadronicWMassNoBTag",50,50,110,"HadronicWMassNoBTag");
    MSPlot["HadronicTopMassNoBTag"]= new MultiSamplePlot(datasets,"HadronicTopMassNoBTag",50,130,210,"HadronicTopMassNoBTag");
    MSPlot["LeptonicWMassNoBTag"]= new MultiSamplePlot(datasets,"LeptonicWMassNoBTag",50,500,110,"LeptonicWMassNoBTag");
    MSPlot["LeptonicTopMassNoBTag"]= new MultiSamplePlot(datasets,"LeptonicTopMassNoBTag",50,130,210,"LeptonicTopMassNoBTag");
    
    //Specific b-tag applied
    MSPlot["JetCombCosThetaFoundSpecificBTag"]= new MultiSamplePlot(datasets,"JetCombCosThetaFoundSpecificBTag",13,-0.5,12.5,"JetCombCosThetaFoundSpecificBTag");
    MSPlot["ChiSqKinFitCosThetaFoundSpecificBTag"]= new MultiSamplePlot(datasets,"ChiSqKinFitCosThetaFoundSpecificBTag",40,0,150,"ChiSqKinFitCosThetaFoundSpecificBTag");
    MSPlot["ProbKinFitCosThetaFoundSpecificBTag"]= new MultiSamplePlot(datasets,"ProbKinFitCosThetaFoundSpecificBTag",20,0,1,"ProbKinFitCosThetaFoundSpecificBTag");

    MSPlot["CosThetaSpecificBTag50Bins"]= new MultiSamplePlot(datasets,"CosThetaSpecificBTag50Bins",50,-1,1,"CosThetaSpecificBTag50Bins");
    MSPlot["CosThetaSpecificBTag30Bins"]= new MultiSamplePlot(datasets,"CosThetaSpecificBTag30Bins",30,-1,1,"CosThetaSpecificBTag30Bins");
    MSPlot["CosThetaSpecificDoubleBTag50Bins"]= new MultiSamplePlot(datasets,"CosThetaSpecificDoubleBTag50Bins",50,-1,1,"CosThetaSpecificDoubleBTag50Bins");
    MSPlot["CosThetaSpecificDoubleBTag30Bins"]= new MultiSamplePlot(datasets,"CosThetaSpecificDoubleBTag30Bins",30,-1,1,"CosThetaSpecificDoubleBTag30Bins");

    MSPlot["nPVSpecificBTag"]= new MultiSamplePlot(datasets,"nPVSpecificBTag",25,0,25,"nPVSpecificBTag");
    MSPlot["TransverseMassSpecificBTag"]= new MultiSamplePlot(datasets,"TransverseMassSpecificBTag",30,0,200,"TransverseMassSpecificBTag");
    MSPlot["TransverseMassSelectedParticlesSpecificBTag"]= new MultiSamplePlot(datasets,"TransverseMassSelectedParticlesSpecificBTag",30,0,200,"TransverseMassSelectedParticlesSpecificBTag");
  
    MSPlot["Jet1PtSpecificBTag"]= new MultiSamplePlot(datasets,"Jet1PtSpecificBTag",100,0,500,"Jet1PtSpecificBTag");
    MSPlot["Jet2PtSpecificBTag"]= new MultiSamplePlot(datasets,"Jet2PtSpecificBTag",100,0,500,"Jet2PtSpecificBTag");
    MSPlot["Jet3PtSpecificBTag"]= new MultiSamplePlot(datasets,"Jet3PtSpecificBTag",100,0,500,"Jet3PtSpecificBTag");
    MSPlot["Jet4PtSpecificBTag"]= new MultiSamplePlot(datasets,"Jet4PtSpecificBTag",100,0,500,"Jet4PtSpecificBTag");

    MSPlot["PtLeptonicBSpecificBTag"]= new MultiSamplePlot(datasets,"PtLeptonicBSpecificBTag",100,0,500,"PtLeptonicBSpecificBTag");
    MSPlot["PtHadronicBSpecificBTag"]= new MultiSamplePlot(datasets,"PtHadronicBSpecificBTag",100,0,500,"PtHadronicBSpecificBTag");
    MSPlot["PtLightJet1SpecificBTag"]= new MultiSamplePlot(datasets,"PtLightJet1SpecificBTag",100,0,500,"PtLightJet1SpecificBTag");
    MSPlot["PtLightJet2SpecificBTag"]= new MultiSamplePlot(datasets,"PtLightJet2SpecificBTag",100,0,500,"PtLightJet2SpecificBTag");

    MSPlot["leptonPhiSpecificBTag"]= new MultiSamplePlot(datasets,"leptonPhiSpecificBTag",30,-4,4,"leptonPhiSpecificBTag");
    MSPlot["leptonPxSpecificBTag"]= new MultiSamplePlot(datasets,"leptonPxSpecificBTag",50,-200,200,"leptonPxSpecificBTag");
    MSPlot["leptonPySpecificBTag"]= new MultiSamplePlot(datasets,"leptonPySpecificBTag",50,-200,200,"leptonPySpecificBTag");
    MSPlot["leptonPzSpecificBTag"]= new MultiSamplePlot(datasets,"leptonPzSpecificBTag",50,-500,500,"leptonPzSpecificBTag");
    MSPlot["leptonPtSpecificBTag"]= new MultiSamplePlot(datasets,"leptonPtSpecificBTag",100,0,500,"leptonPtSpecificBTag");
    MSPlot["leptonMassSpecificBTag"]= new MultiSamplePlot(datasets,"leptonMassSpecificBTag",30,0.1,0.15,"leptonMassSpecificBTag");

    MSPlot["METPhiSpecificBTag"]= new MultiSamplePlot(datasets,"METPhiSpecificBTag",30,-4,4,"METPhiSpecificBTag");
    MSPlot["METPxSpecificBTag"]= new MultiSamplePlot(datasets,"METPxSpecificBTag",50,-200,200,"METPxSpecificBTag");
    MSPlot["METPySpecificBTag"]= new MultiSamplePlot(datasets,"METPySpecificBTag",50,-200,200,"METPySpecificBTag");
    MSPlot["METPzSpecificBTag"]= new MultiSamplePlot(datasets,"METPzSpecificBTag",50,-500,500,"METPzSpecificBTag");
    MSPlot["METPtSpecificBTag"]= new MultiSamplePlot(datasets,"METPtSpecificBTag",100,0,500,"METPtSpecificBTag");
    MSPlot["METMassSpecificBTag"]= new MultiSamplePlot(datasets,"METMassSpecificBTag",30,0,0.01,"METMassSpecificBTag");  

    MSPlot["HadronicWMassSpecificBTag"]= new MultiSamplePlot(datasets,"HadronicWMassSpecificBTag",50,50,110,"HadronicWMassSpecificBTag");
    MSPlot["HadronicTopMassSpecificBTag"]= new MultiSamplePlot(datasets,"HadronicTopMassSpecificBTag",50,130,210,"HadronicTopMassSpecificBTag");
    MSPlot["LeptonicWMassSpecificBTag"]= new MultiSamplePlot(datasets,"LeptonicWMassSpecificBTag",50,500,110,"LeptonicWMassSpecificBTag");
    MSPlot["LeptonicTopMassSpecificBTag"]= new MultiSamplePlot(datasets,"LeptonicTopMassSpecificBTag",50,130,210,"LeptonicTopMassSpecificBTag");

   
    //Check nPrimary vertices for different executed cuts !!
    MSPlot["nPVBeforeCuts"] = new MultiSamplePlot(datasets,"nPVBeforeCuts" , 20, 0, 20, "nPVBeforeCuts");
    MSPlot["nPVAfterTransverseMassCut"] = new MultiSamplePlot(datasets,"nPVAfterTransverseMassCut" , 20, 0, 20, "nPVAfterTransverseMassCut");

    MSPlot["TransverseMassBeforeCut"]=new MultiSamplePlot(datasets,"TransverseMassBeforeCut",50,0,200,"TransverseMassBeforeCut");
    MSPlot["TransverseMassAfterCut"]=new MultiSamplePlot(datasets,"TransverseMassAfterCut",50,0,200,"TransverseMassAfterCut"); 
  
    MSPlot["RelIsoLepton"]=new MultiSamplePlot(datasets,"RelIsoLepton",50,0,0.5,"RelIsoLepton");

    MSPlot["ChiSqHadrComb0"] = new MultiSamplePlot(datasets,"ChiSqHadrComb0",40,0,150,"ChiSqHadrComb0");
    MSPlot["ChiSqHadrComb1"] = new MultiSamplePlot(datasets,"ChiSqHadrComb1",40,0,150,"ChiSqHadrComb1");
    MSPlot["ChiSqHadrComb2"] = new MultiSamplePlot(datasets,"ChiSqHadrComb2",40,0,150,"ChiSqHadrComb2");
    MSPlot["ChiSqHadrComb3"] = new MultiSamplePlot(datasets,"ChiSqHadrComb3",40,0,150,"ChiSqHadrComb3");
    MSPlot["ChiSqHadrComb4"] = new MultiSamplePlot(datasets,"ChiSqHadrComb4",40,0,150,"ChiSqHadrComb4");
    MSPlot["ChiSqHadrComb5"] = new MultiSamplePlot(datasets,"ChiSqHadrComb5",40,0,150,"ChiSqHadrComb5");
    MSPlot["ChiSqHadrComb6"] = new MultiSamplePlot(datasets,"ChiSqHadrComb6",40,0,150,"ChiSqHadrComb6");
    MSPlot["ChiSqHadrComb7"] = new MultiSamplePlot(datasets,"ChiSqHadrComb7",40,0,150,"ChiSqHadrComb7");
    MSPlot["ChiSqHadrComb8"] = new MultiSamplePlot(datasets,"ChiSqHadrComb8",40,0,150,"ChiSqHadrComb8");
    MSPlot["ChiSqHadrComb9"] = new MultiSamplePlot(datasets,"ChiSqHadrComb9",40,0,150,"ChiSqHadrComb9");
    MSPlot["ChiSqHadrComb10"] = new MultiSamplePlot(datasets,"ChiSqHadrComb10",40,0,150,"ChiSqHadrComb10");
    MSPlot["ChiSqHadrComb11"] = new MultiSamplePlot(datasets,"ChiSqHadrComb11",40,0,150,"ChiSqHadrComb11");

    MSPlot["ChiSqHadrAndLeptWOnlyComb0"] = new MultiSamplePlot(datasets,"ChiSqHadrAndLeptWOnlyComb0",40,0,150,"ChiSqHadrAndLeptWOnlyComb0");
    MSPlot["ChiSqHadrAndLeptWOnlyComb1"] = new MultiSamplePlot(datasets,"ChiSqHadrAndLeptWOnlyComb1",40,0,150,"ChiSqHadrAndLeptWOnlyComb1");
    MSPlot["ChiSqHadrAndLeptWOnlyComb2"] = new MultiSamplePlot(datasets,"ChiSqHadrAndLeptWOnlyComb2",40,0,150,"ChiSqHadrAndLeptWOnlyComb2");
    MSPlot["ChiSqHadrAndLeptWOnlyComb3"] = new MultiSamplePlot(datasets,"ChiSqHadrAndLeptWOnlyComb3",40,0,150,"ChiSqHadrAndLeptWOnlyComb3");
    MSPlot["ChiSqHadrAndLeptWOnlyComb4"] = new MultiSamplePlot(datasets,"ChiSqHadrAndLeptWOnlyComb4",40,0,150,"ChiSqHadrAndLeptWOnlyComb4");
    MSPlot["ChiSqHadrAndLeptWOnlyComb5"] = new MultiSamplePlot(datasets,"ChiSqHadrAndLeptWOnlyComb5",40,0,150,"ChiSqHadrAndLeptWOnlyComb5");
    MSPlot["ChiSqHadrAndLeptWOnlyComb6"] = new MultiSamplePlot(datasets,"ChiSqHadrAndLeptWOnlyComb6",40,0,150,"ChiSqHadrAndLeptWOnlyComb6");
    MSPlot["ChiSqHadrAndLeptWOnlyComb7"] = new MultiSamplePlot(datasets,"ChiSqHadrAndLeptWOnlyComb7",40,0,150,"ChiSqHadrAndLeptWOnlyComb7");
    MSPlot["ChiSqHadrAndLeptWOnlyComb8"] = new MultiSamplePlot(datasets,"ChiSqHadrAndLeptWOnlyComb8",40,0,150,"ChiSqHadrAndLeptWOnlyComb8");
    MSPlot["ChiSqHadrAndLeptWOnlyComb9"] = new MultiSamplePlot(datasets,"ChiSqHadrAndLeptWOnlyComb9",40,0,150,"ChiSqHadrAndLeptWOnlyComb9");
    MSPlot["ChiSqHadrAndLeptWOnlyComb10"] = new MultiSamplePlot(datasets,"ChiSqHadrAndLeptWOnlyComb10",40,0,150,"ChiSqHadrAndLeptWOnlyComb10");
    MSPlot["ChiSqHadrAndLeptWOnlyComb11"] = new MultiSamplePlot(datasets,"ChiSqHadrAndLeptWOnlyComb11",40,0,150,"ChiSqHadrAndLeptWOnlyComb11");

    MSPlot["ChiSqHadrAndLeptComb0"] = new MultiSamplePlot(datasets,"ChiSqHadrAndLeptComb0",40,0,150,"ChiSqHadrAndLeptComb0");
    MSPlot["ChiSqHadrAndLeptComb1"] = new MultiSamplePlot(datasets,"ChiSqHadrAndLeptComb1",40,0,150,"ChiSqHadrAndLeptComb1");
    MSPlot["ChiSqHadrAndLeptComb2"] = new MultiSamplePlot(datasets,"ChiSqHadrAndLeptComb2",40,0,150,"ChiSqHadrAndLeptComb2");
    MSPlot["ChiSqHadrAndLeptComb3"] = new MultiSamplePlot(datasets,"ChiSqHadrAndLeptComb3",40,0,150,"ChiSqHadrAndLeptComb3");
    MSPlot["ChiSqHadrAndLeptComb4"] = new MultiSamplePlot(datasets,"ChiSqHadrAndLeptComb4",40,0,150,"ChiSqHadrAndLeptComb4");
    MSPlot["ChiSqHadrAndLeptComb5"] = new MultiSamplePlot(datasets,"ChiSqHadrAndLeptComb5",40,0,150,"ChiSqHadrAndLeptComb5");
    MSPlot["ChiSqHadrAndLeptComb6"] = new MultiSamplePlot(datasets,"ChiSqHadrAndLeptComb6",40,0,150,"ChiSqHadrAndLeptComb6");
    MSPlot["ChiSqHadrAndLeptComb7"] = new MultiSamplePlot(datasets,"ChiSqHadrAndLeptComb7",40,0,150,"ChiSqHadrAndLeptComb7");
    MSPlot["ChiSqHadrAndLeptComb8"] = new MultiSamplePlot(datasets,"ChiSqHadrAndLeptComb8",40,0,150,"ChiSqHadrAndLeptComb8");
    MSPlot["ChiSqHadrAndLeptComb9"] = new MultiSamplePlot(datasets,"ChiSqHadrAndLeptComb9",40,0,150,"ChiSqHadrAndLeptComb9");
    MSPlot["ChiSqHadrAndLeptComb10"] = new MultiSamplePlot(datasets,"ChiSqHadrAndLeptComb10",40,0,150,"ChiSqHadrAndLeptComb10");
    MSPlot["ChiSqHadrAndLeptComb11"] = new MultiSamplePlot(datasets,"ChiSqHadrAndLeptComb11",40,0,150,"ChiSqHadrAndLeptComb11");
  
    MSPlot["KinFitProbabilityHadrComb0"]=new MultiSamplePlot(datasets,"KinFitProbabilityHadrComb0",50,-1,1,"KinFitProbabilityHadrComb0");  
    MSPlot["KinFitProbabilityHadrComb1"]=new MultiSamplePlot(datasets,"KinFitProbabilityHadrComb1",50,-1,1,"KinFitProbabilityHadrComb1");  
    MSPlot["KinFitProbabilityHadrComb2"]=new MultiSamplePlot(datasets,"KinFitProbabilityHadrComb2",50,-1,1,"KinFitProbabilityHadrComb2");  
    MSPlot["KinFitProbabilityHadrComb3"]=new MultiSamplePlot(datasets,"KinFitProbabilityHadrComb3",50,-1,1,"KinFitProbabilityHadrComb3");  
    MSPlot["KinFitProbabilityHadrComb4"]=new MultiSamplePlot(datasets,"KinFitProbabilityHadrComb4",50,-1,1,"KinFitProbabilityHadrComb4");  
    MSPlot["KinFitProbabilityHadrComb5"]=new MultiSamplePlot(datasets,"KinFitProbabilityHadrComb5",50,-1,1,"KinFitProbabilityHadrComb5");  
    MSPlot["KinFitProbabilityHadrComb6"]=new MultiSamplePlot(datasets,"KinFitProbabilityHadrComb6",50,-1,1,"KinFitProbabilityHadrComb6");  
    MSPlot["KinFitProbabilityHadrComb7"]=new MultiSamplePlot(datasets,"KinFitProbabilityHadrComb7",50,-1,1,"KinFitProbabilityHadrComb7");  
    MSPlot["KinFitProbabilityHadrComb8"]=new MultiSamplePlot(datasets,"KinFitProbabilityHadrComb8",50,-1,1,"KinFitProbabilityHadrComb8");  
    MSPlot["KinFitProbabilityHadrComb9"]=new MultiSamplePlot(datasets,"KinFitProbabilityHadrComb9",50,-1,1,"KinFitProbabilityHadrComb9");  
    MSPlot["KinFitProbabilityHadrComb10"]=new MultiSamplePlot(datasets,"KinFitProbabilityHadrComb10",50,-1,1,"KinFitProbabilityHadrComb10");  
    MSPlot["KinFitProbabilityHadrComb11"]=new MultiSamplePlot(datasets,"KinFitProbabilityHadrComb11",50,-1,1,"KinFitProbabilityHadrComb11");  

    MSPlot["KinFitProbabilityHadrAndLeptWOnlyComb0"]=new MultiSamplePlot(datasets,"KinFitProbabilityHadrAndLeptWOnlyComb0",50,-1,1,"KinFitProbabilityHadrAndLeptWOnlyComb0");
    MSPlot["KinFitProbabilityHadrAndLeptWOnlyComb1"]=new MultiSamplePlot(datasets,"KinFitProbabilityHadrAndLeptWOnlyComb1",50,-1,1,"KinFitProbabilityHadrAndLeptWOnlyComb1");
    MSPlot["KinFitProbabilityHadrAndLeptWOnlyComb2"]=new MultiSamplePlot(datasets,"KinFitProbabilityHadrAndLeptWOnlyComb2",50,-1,1,"KinFitProbabilityHadrAndLeptWOnlyComb2");
    MSPlot["KinFitProbabilityHadrAndLeptWOnlyComb3"]=new MultiSamplePlot(datasets,"KinFitProbabilityHadrAndLeptWOnlyComb3",50,-1,1,"KinFitProbabilityHadrAndLeptWOnlyComb3");
    MSPlot["KinFitProbabilityHadrAndLeptWOnlyComb4"]=new MultiSamplePlot(datasets,"KinFitProbabilityHadrAndLeptWOnlyComb4",50,-1,1,"KinFitProbabilityHadrAndLeptWOnlyComb4");
    MSPlot["KinFitProbabilityHadrAndLeptWOnlyComb5"]=new MultiSamplePlot(datasets,"KinFitProbabilityHadrAndLeptWOnlyComb5",50,-1,1,"KinFitProbabilityHadrAndLeptWOnlyComb5");
    MSPlot["KinFitProbabilityHadrAndLeptWOnlyComb6"]=new MultiSamplePlot(datasets,"KinFitProbabilityHadrAndLeptWOnlyComb6",50,-1,1,"KinFitProbabilityHadrAndLeptWOnlyComb6");
    MSPlot["KinFitProbabilityHadrAndLeptWOnlyComb7"]=new MultiSamplePlot(datasets,"KinFitProbabilityHadrAndLeptWOnlyComb7",50,-1,1,"KinFitProbabilityHadrAndLeptWOnlyComb7");
    MSPlot["KinFitProbabilityHadrAndLeptWOnlyComb8"]=new MultiSamplePlot(datasets,"KinFitProbabilityHadrAndLeptWOnlyComb8",50,-1,1,"KinFitProbabilityHadrAndLeptWOnlyComb8");
    MSPlot["KinFitProbabilityHadrAndLeptWOnlyComb9"]=new MultiSamplePlot(datasets,"KinFitProbabilityHadrAndLeptWOnlyComb9",50,-1,1,"KinFitProbabilityHadrAndLeptWOnlyComb9");
    MSPlot["KinFitProbabilityHadrAndLeptWOnlyComb10"]=new MultiSamplePlot(datasets,"KinFitProbabilityHadrAndLeptWOnlyComb10",50,-1,1,"KinFitProbabilityHadrAndLeptWOnlyComb10");
    MSPlot["KinFitProbabilityHadrAndLeptWOnlyComb11"]=new MultiSamplePlot(datasets,"KinFitProbabilityHadrAndLeptWOnlyComb11",50,-1,1,"KinFitProbabilityHadrAndLeptWOnlyComb11");

    MSPlot["KinFitProbabilityHadrAndLeptComb0"]=new MultiSamplePlot(datasets,"KinFitProbabilityHadrAndLeptComb0",50,-1,1,"KinFitProbabilityHadrAndLeptComb0");
    MSPlot["KinFitProbabilityHadrAndLeptComb1"]=new MultiSamplePlot(datasets,"KinFitProbabilityHadrAndLeptComb1",50,-1,1,"KinFitProbabilityHadrAndLeptComb1");
    MSPlot["KinFitProbabilityHadrAndLeptComb2"]=new MultiSamplePlot(datasets,"KinFitProbabilityHadrAndLeptComb2",50,-1,1,"KinFitProbabilityHadrAndLeptComb2");
    MSPlot["KinFitProbabilityHadrAndLeptComb3"]=new MultiSamplePlot(datasets,"KinFitProbabilityHadrAndLeptComb3",50,-1,1,"KinFitProbabilityHadrAndLeptComb3");
    MSPlot["KinFitProbabilityHadrAndLeptComb4"]=new MultiSamplePlot(datasets,"KinFitProbabilityHadrAndLeptComb4",50,-1,1,"KinFitProbabilityHadrAndLeptComb4");
    MSPlot["KinFitProbabilityHadrAndLeptComb5"]=new MultiSamplePlot(datasets,"KinFitProbabilityHadrAndLeptComb5",50,-1,1,"KinFitProbabilityHadrAndLeptComb5");
    MSPlot["KinFitProbabilityHadrAndLeptComb6"]=new MultiSamplePlot(datasets,"KinFitProbabilityHadrAndLeptComb6",50,-1,1,"KinFitProbabilityHadrAndLeptComb6");
    MSPlot["KinFitProbabilityHadrAndLeptComb7"]=new MultiSamplePlot(datasets,"KinFitProbabilityHadrAndLeptComb7",50,-1,1,"KinFitProbabilityHadrAndLeptComb7");
    MSPlot["KinFitProbabilityHadrAndLeptComb8"]=new MultiSamplePlot(datasets,"KinFitProbabilityHadrAndLeptComb8",50,-1,1,"KinFitProbabilityHadrAndLeptComb8");
    MSPlot["KinFitProbabilityHadrAndLeptComb9"]=new MultiSamplePlot(datasets,"KinFitProbabilityHadrAndLeptComb9",50,-1,1,"KinFitProbabilityHadrAndLeptComb9");
    MSPlot["KinFitProbabilityHadrAndLeptComb10"]=new MultiSamplePlot(datasets,"KinFitProbabilityHadrAndLeptComb10",50,-1,1,"KinFitProbabilityHadrAndLeptComb10");
    MSPlot["KinFitProbabilityHadrAndLeptComb11"]=new MultiSamplePlot(datasets,"KinFitProbabilityHadrAndLeptComb11",50,-1,1,"KinFitProbabilityHadrAndLeptComb11");
  
    MSPlot["Jet1Pt"] = new MultiSamplePlot(datasets, "Jet1Pt", 100,0,500,"Jet1Pt");
    MSPlot["Jet2Pt"] = new MultiSamplePlot(datasets, "Jet2Pt", 100,0,500,"Jet2Pt");
    MSPlot["Jet3Pt"] = new MultiSamplePlot(datasets, "Jet3Pt", 100,0,500,"Jet3Pt");
    MSPlot["Jet4Pt"] = new MultiSamplePlot(datasets, "Jet4Pt", 100,0,500,"Jet4Pt");
    
    MSPlot["DeltaRJetLepton"] =new MultiSamplePlot(datasets,"DeltaRJetLepton",50,0,4,"DeltaRJetLepton");
    MSPlot["DeltaRMuonJet"] =new MultiSamplePlot(datasets,"DeltaRMuonJet",50,0,4,"DeltaRMuonJet");
    MSPlot["DeltaRElectronJet"] =new MultiSamplePlot(datasets,"DeltaRElectronJet",50,0,10,"DeltaRElectronJet");
  }

  ////////////////////////////////
  //     Selection Table        // 
  ////////////////////////////////

  vector<string> CutsSelecTableMacro;
  CutsSelecTableMacro.push_back(string("offline cuts"));
  CutsSelecTableMacro.push_back(string("M SSVHE bTag"));
  CutsSelecTableMacro.push_back(string("Muon Pt Cut"));
  CutsSelecTableMacro.push_back(string("TransverseMass Cut"));

  SelectionTable selecTableMacro(CutsSelecTableMacro,datasets);
  selecTableMacro.SetLuminosity(Luminosity);

  // initialize LumiReWeighting stuff  
  cout << " Starting LumiReWeighting stuff " << endl;

  Lumi3DReWeighting Lumi3DWeights;
  Lumi3DWeights = Lumi3DReWeighting("PileUpReweighting/pileup_MC_Fall11.root","PileUpReweighting/pileup_FineBin_2011Data_UpToRun180252.root", "pileup", "pileup");
  Lumi3DWeights.weight3D_init(1.0);

  if(dataSetName.find("PUMinus_") == true && PUResults == true)
    Lumi3DWeights.weight3D_init(0.92);
  if(dataSetName.find("PUPlus_") == true && PUResults == true)
    Lumi3DWeights.weight3D_init(1.08);
  
  PoissonMeanShifter PShiftUp = PoissonMeanShifter(0.6); // PU-systematic
  PoissonMeanShifter PShiftDown = PoissonMeanShifter(-0.6); // PU-systematic

  /////////////////////////////////////////
  // Initializing used variables
  /////////////////////////////////////////
  float MassW=83.6103;
  float MassTop = 172.956;
  float SigmaW=11.1534;  //Obtained from gaussian fit on Top and W distribution with simulated information
  float SigmaTop=18.232;

  std::string bTag[14];
  for(int ii=1;ii<=14;ii++){
    std::stringstream out;
    out << ii;
    bTag[ii-1]  =  out.str();
  }

  int NumberTCHEbTags = 13;
  int NumberTCHPbTags = 13;
  int NumberSSVHEbTags = 13;
  int NumberSSVHPbTags=13;
  int NumberCSVbTags=13;
  int NumberJPbTags=13;
  int NumberJBPbTags=13;
  int TotalNumberbTags = NumberTCHEbTags + 1 + NumberTCHPbTags+1 + NumberSSVHEbTags + 1 + NumberSSVHPbTags + 1 + NumberCSVbTags + 1 + NumberJPbTags + 1 + NumberJBPbTags + 1;
  
  std::string bTagFileOutput[TotalNumberbTags];
  map<int, int[5]> ProcessedbTags;
  std::string PresentationOutput[TotalNumberbTags];
  int NumberRemainingEventsKinFitHadr[TotalNumberbTags][inputWTree.size()];
  int NumberRemainingEventsKinFitHadrAndLeptWOnly[TotalNumberbTags][inputWTree.size()];
  int NumberRemainingEventsKinFitHadrAndLept[TotalNumberbTags][inputWTree.size()];
  int NumberBLeptCorrectEventsKinFitHadr[TotalNumberbTags][inputWTree.size()];
  int NumberBLeptCorrectEventsKinFitHadrAndLeptWOnly[TotalNumberbTags][inputWTree.size()];
  int NumberBLeptCorrectEventsKinFitHadrAndLept[TotalNumberbTags][inputWTree.size()];

  //Initialize:
  for(int ii=0;ii<14;ii++){  //tche
    for(int jj=0;jj<14;jj++){  //tchp
      for(int kk=0;kk<14;kk++){  //ssvhe
	for(int ll=0;ll<14;ll++){  //ssvhp
	  for(int nn=0;nn<14;nn++){   //csv
	    for(int qq=0;qq<14;qq++){   //jp
	      for(int pp=0;pp<14;pp++){   //jbp
		bTagFileOutput[ii+jj+kk+ll+nn+qq+pp]=" Wrong entry chosen";
		PresentationOutput[ii+jj+kk+ll+nn+qq+pp]=" Wrong entry chosen";
		for(int mm=0;mm<inputWTree.size();mm++){
		  NumberRemainingEventsKinFitHadr[ii+jj+kk+ll+nn+qq+pp][mm]=0;
		  NumberRemainingEventsKinFitHadrAndLeptWOnly[ii+jj+kk+ll+nn+qq+pp][mm]=0;
		  NumberRemainingEventsKinFitHadrAndLept[ii+jj+kk+ll+nn+qq+pp][mm]=0;
		  NumberBLeptCorrectEventsKinFitHadr[ii+jj+kk+ll+nn+qq+pp][mm]=0;
		  NumberBLeptCorrectEventsKinFitHadrAndLeptWOnly[ii+jj+kk+ll+nn+qq+pp][mm]=0;
		  NumberBLeptCorrectEventsKinFitHadrAndLept[ii+jj+kk+ll+nn+qq+pp][mm]=0;
		}
	      }
	    }
	  }
	}
      }
    }
  }

  //Count number of events for prob KinFit cut (cut at 0.01) --> Check behavior compared to percentage correct events
  int NumberEventsSignalCosThetaFound = 0;
  int NumberEventsLowProbCosThetaFound = 0;
  int NumberEventsSignalCorrectLeptBCosThetaFound = 0;
  int NumberEventsSignalCorrectHadrBCosThetaFound = 0;
  int NumberEventsSignalCorrectLightCosThetaFound = 0;
  int NumberEventsSignalCorrectHadronicCosThetaFound = 0;
  int NumberEventsLowProbCorrectLeptBCosThetaFound = 0;
  int NumberEventsLowProbCorrectHadrBCosThetaFound = 0;
  int NumberEventsLowProbCorrectLightCosThetaFound = 0;
  int NumberEventsLowProbCorrectHadronicCosThetaFound = 0;
 
  //Count number of events for different b-tag method:
  int NumberEventsChiSqOnly[inputWTree.size()];
  int NumberEventsChiSqOnlyJetCombFound[inputWTree.size()];
  int NumberEventsChiSqOnlyCosThetaFound[inputWTree.size()];
  int NumberEventsChiSqOnlySpecificBTag[inputWTree.size()];

  int NumberEventsChiSqAndBTag[inputWTree.size()];
  int NumberEventsChiSqAndBTagJetCombFound[inputWTree.size()];
  int NumberEventsChiSqAndBTagCosThetaFound[inputWTree.size()];
  int NumberEventsChiSqAndBTagSpecificBTag[inputWTree.size()];

  for(int iCombi =0; iCombi < inputWTree.size(); iCombi++){
    NumberEventsChiSqOnly[iCombi] = 0;
    NumberEventsChiSqOnlyJetCombFound[iCombi] = 0;
    NumberEventsChiSqOnlyCosThetaFound[iCombi] = 0;
    NumberEventsChiSqOnlySpecificBTag[iCombi] = 0;
    
    NumberEventsChiSqAndBTag[iCombi] = 0;
    NumberEventsChiSqAndBTagJetCombFound[iCombi] = 0;
    NumberEventsChiSqAndBTagCosThetaFound[iCombi] = 0;
    NumberEventsChiSqAndBTagSpecificBTag[iCombi] = 0;
  }

  //Call classes made :
  BTagName bTagName = BTagName();  //for bTagFileOutput name giving
  BTagJetSelection bTagJetSelection = BTagJetSelection();
  BTagCosThetaCalculation bTagCosThetaCalculation = BTagCosThetaCalculation();

  int LengthOfPresentationArray=0;
  int NumberSelectedEvents =0;
  int NumberEventsBeforeCuts = 0;
  int NumberSelectedDataEvents = 0;
  int NumberDataEventsBeforeCuts = 0;

  //Order of indices used in Kinematic Fit:
  int BLeptIndex[12]={3,2,3,1,2,1,3,0,2,0,1,0};
  int BHadrIndex[12]={2,3,1,3,1,2,0,3,0,2,0,1};
  int Quark1Index[12]={0,0,0,0,0,0,1,1,1,1,2,2};
  int Quark2Index[12]={1,1,2,2,3,3,2,2,3,3,3,3};
  
  //BTag values:
  float TCHEbTag[3] = {1.7,3.3,10.2};
  float TCHPbTag[3] = {1.19,1.93,3.41};
  float SSVHEbTag[3] = {0.,1.74,3.05};
  float SSVHPbTag[3] = {0.,1.,2.};
  float CSVbTag[3] = {0.244,0.679,0.898};
  float JPbTag[3] = {0.275,0.545,0.790};
  float JBPbTag[3] = {1.33,2.55,3.744};

  /////////////////////////////////////////
  // Loop on datasets
  /////////////////////////////////////////
  
  for(unsigned int iDataSet=0; iDataSet<inputWTree.size(); iDataSet++){

    string dataSetName = nameDataSet[iDataSet];
    std::cout << " dataSetName : " << dataSetName << endl;
    
    TFile* inFile = new TFile(inputWTree[iDataSet].c_str(),"READ");
    std::cout  << " inputWTree[iDataSet].c_str() " << inputWTree[iDataSet].c_str() << endl;

    TTree* inWTreeTree = (TTree*) inFile->Get("WTreeTree");
    TBranch* m_br = (TBranch*) inWTreeTree->GetBranch("TheWTree");
    
    WTree* wTree = 0;
    m_br->SetAddress(&wTree);
    
    int nEvent = inWTreeTree->GetEntries();
    
    TTree* inConfigTree = (TTree*) inFile->Get("configTreeWTreeFile");
    TBranch* d_br = (TBranch*) inConfigTree->GetBranch("Dataset");
    TClonesArray* tc_dataset = new TClonesArray("Dataset",0);
    d_br->SetAddress(&tc_dataset);
    
    inConfigTree->GetEvent(0);
    Dataset* dataSet = (Dataset*) tc_dataset->At(0);
    cout << "Processing DataSet " << iDataSet << " : " << dataSetName << "  containing " << nEvent << " events" << endl;
    cout << " ***************************************** " << endl;
    cout << "Before changing (Line 1158)--> Cross section = " << dataSet->Xsection() << "  intLumi = " << dataSet->EquivalentLumi() << " Normfactor = " << dataSet->NormFactor() << endl;
    if( dataSetName.find("WSystMinus_WJets") == 0 && WSystResults == true){
      dataSet->SetEquivalentLuminosity( dataSet->EquivalentLumi() / (0.7) );  //WJets Minus 30%
      //dataSet->SetEquivalentLuminosity( dataSet->EquivalentLumi() / (0.0000001) );  //WJets Minus 100%
    }
    if(dataSetName.find("WSystPlus_WJets") == 0 && WSystResults == true){
      dataSet->SetEquivalentLuminosity( dataSet->EquivalentLumi() / (1.3) ); //WJets Plus 30 %
      //dataSet->SetEquivalentLuminosity( dataSet->EquivalentLumi() / (2.) ); //WJets Plus 100 %
    }    
    if( dataSetName.find("Data") == 0 && semiElectron == true && TriCentralJet30Trigger == true){
      dataSet->SetEquivalentLuminosity( Luminosity);
    }
    cout << "After changing --> Cross section = " << dataSet->Xsection() << "  intLumi = " << dataSet->EquivalentLumi() << " Normfactor = " << dataSet->NormFactor() << endl;
    cout << " ************************************** " << endl;

    // output ascii file stuff
    mkdir("WHelResults_ASCII/",0777);
    string outFileName = "WHelResults_ASCII/WHelResults_" + dataSetName + ".txt";
    ofstream outFile(outFileName.c_str());
    
    outFile << "Output of WHelicities_Analysis.cc" << endl;
    outFile << "First some dataSet info: " << endl;
    outFile << "dataSet Name: " << dataSetName << endl;
    outFile << "dataSet Title: " << dataSet->Title() << endl;
    //outFile << "dataSet cross-section: " << dataSet->Xsection() << endl;
    //outFile << "dataSet integrated lumi: " << dataSet->EquivalentLumi() << endl << endl;
    outFile << "dataSet cross-section: " << XSection << endl;
    outFile << "dataSet integrated lumi: " << EqLumi << endl << endl;    
    outFile << "Start of event-by-event info " << endl << endl;
    
    //Value needed to study the reconstruction efficiency of the leptonic b-jet:
    int CorrectBLeptConstruction=0;
    int ReconstructedEvents=0;

    vector<float> CosThetaValuesKinFitHadr[TotalNumberbTags],CosThetaValuesKinFitHadrAndLeptWOnly[TotalNumberbTags],CosThetaValuesKinFitHadrAndLept[TotalNumberbTags];//,LumiWeightVectorKinFitHadr[TotalNumberbTags],LumiWeightVectorKinFitHadrAndLeptWOnly[TotalNumberbTags],LumiWeightVectorKinFitHadrAndLept[TotalNumberbTags];
    int FilledEntries = 0;
    vector<double> CosThGenKinFitHadr[TotalNumberbTags],CosThGenKinFitHadrAndLeptWOnly[TotalNumberbTags],CosThGenKinFitHadrAndLept[TotalNumberbTags],EventCorrectionWeightKinFitHadr[TotalNumberbTags],EventCorrectionWeightKinFitHadrAndLeptWOnly[TotalNumberbTags],EventCorrectionWeightKinFitHadrAndLept[TotalNumberbTags];
    float binEdge[CosThetaBinNumber+1];
    float binSize = (1.-(-1.))/CosThetaBinNumber;
    for(int ii=0; ii<=CosThetaBinNumber;ii++){
      binEdge[ii] = -1 + binSize*ii;
    }

    //Initialize naming of different bTag options:
    int TCHELoop,TCHPLoop,SSVHELoop,SSVHPLoop,CSVLoop,JPLoop, JBPLoop;
    TCHELoop = 1;
    TCHPLoop = 1;
    SSVHELoop = 1;
    SSVHPLoop = 1;
    CSVLoop =1;
    JPLoop = 1;
    JBPLoop = 1;
       
    int UsedTCHE[TotalNumberbTags],UsedTCHP[TotalNumberbTags],UsedSSVHE[TotalNumberbTags],UsedSSVHP[TotalNumberbTags],UsedCSV[TotalNumberbTags], UsedJP[TotalNumberbTags], UsedJBP[TotalNumberbTags];
    while(JBPLoop <= NumberJBPbTags){
      while(JPLoop <= NumberJPbTags){
	while(CSVLoop<=NumberCSVbTags){
	  while(SSVHPLoop<=NumberSSVHPbTags){  
	    while(SSVHELoop<=NumberSSVHEbTags){
	      while(TCHPLoop<=NumberTCHPbTags){
		while(TCHELoop<=NumberTCHEbTags){
		  if(iDataSet==0){
		    bTagFileOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = bTagName.NameGiving(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags,JPLoop, NumberJPbTags,JBPLoop, NumberJBPbTags);
		    PresentationOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = bTagName.NameGivingPres(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags,JPLoop, NumberJPbTags,JBPLoop, NumberJBPbTags);

		  }
		  UsedTCHE[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = TCHELoop;
		  UsedTCHP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = TCHPLoop;
		  UsedSSVHE[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = SSVHELoop;
		  UsedSSVHP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = SSVHPLoop;
		  UsedCSV[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = CSVLoop;
		  UsedJP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = JPLoop;
		  UsedJBP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = JBPLoop;
		  
		  TCHELoop++;
		  if(TCHELoop == 14){TCHPLoop=2;SSVHELoop=2;SSVHPLoop=2;CSVLoop=2;JPLoop=2;JBPLoop=2;}
		}
		if(iDataSet==0){
		  bTagFileOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = bTagName.NameGiving(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags,JPLoop, NumberJPbTags,JBPLoop, NumberJBPbTags);
		  PresentationOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = bTagName.NameGivingPres(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags,JPLoop, NumberJPbTags,JBPLoop, NumberJBPbTags);
		}
		UsedTCHE[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = TCHELoop;
		UsedTCHP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = TCHPLoop;
		UsedSSVHE[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = SSVHELoop;
		UsedSSVHP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = SSVHPLoop;
		UsedCSV[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = CSVLoop;
		UsedJP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = JPLoop;
		UsedJBP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = JBPLoop;
		
		TCHPLoop++;
	      }
	      if(iDataSet==0){
		bTagFileOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = bTagName.NameGiving(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags,JPLoop, NumberJPbTags,JBPLoop, NumberJBPbTags);
		PresentationOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = bTagName.NameGivingPres(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags,JPLoop, NumberJPbTags,JBPLoop, NumberJBPbTags);
	      }
	      UsedTCHE[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = TCHELoop;
	      UsedTCHP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = TCHPLoop;
	      UsedSSVHE[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = SSVHELoop;
	      UsedSSVHP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = SSVHPLoop;
	      UsedCSV[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = CSVLoop;
	      UsedJP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = JPLoop;
	      UsedJBP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = JBPLoop;
	      
	      SSVHELoop++;
	    }
	    if(iDataSet==0){
	      bTagFileOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = bTagName.NameGiving(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags,JPLoop, NumberJPbTags,JBPLoop, NumberJBPbTags);
	      PresentationOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = bTagName.NameGivingPres(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags,JPLoop, NumberJPbTags,JBPLoop, NumberJBPbTags);
	    }
	    UsedTCHE[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = TCHELoop;
	    UsedTCHP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = TCHPLoop;
	    UsedSSVHE[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = SSVHELoop;
	    UsedSSVHP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = SSVHPLoop;
	    UsedCSV[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = CSVLoop;
	    UsedJP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = JPLoop;
	    UsedJBP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = JBPLoop;
	    
	    SSVHPLoop++;		
	  }	
	  if(iDataSet==0){
	    bTagFileOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = bTagName.NameGiving(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags,JPLoop, NumberJPbTags,JBPLoop, NumberJBPbTags);
	    PresentationOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = bTagName.NameGivingPres(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags,JPLoop, NumberJPbTags,JBPLoop, NumberJBPbTags);
	  }
	  UsedTCHE[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = TCHELoop;
	  UsedTCHP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = TCHPLoop;
	  UsedSSVHE[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = SSVHELoop;
	  UsedSSVHP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = SSVHPLoop;
	  UsedCSV[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = CSVLoop;
	  UsedJP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = JPLoop;
	  UsedJBP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = JBPLoop;
	  
	  CSVLoop++;
	}
	if(iDataSet==0){
	  bTagFileOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = bTagName.NameGiving(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags,JPLoop, NumberJPbTags,JBPLoop, NumberJBPbTags);
	  PresentationOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = bTagName.NameGivingPres(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags,JPLoop, NumberJPbTags,JBPLoop, NumberJBPbTags);
	}
	UsedTCHE[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = TCHELoop;
	UsedTCHP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = TCHPLoop;
	UsedSSVHE[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = SSVHELoop;
	UsedSSVHP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = SSVHPLoop;
	UsedCSV[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = CSVLoop;
	UsedJP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = JPLoop;
	UsedJBP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = JBPLoop;
	
	JPLoop++;
      }
      if(iDataSet==0){
	bTagFileOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = bTagName.NameGiving(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags,JPLoop, NumberJPbTags,JBPLoop, NumberJBPbTags);
	PresentationOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = bTagName.NameGivingPres(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags,JPLoop, NumberJPbTags,JBPLoop, NumberJBPbTags);
      }
      UsedTCHE[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = TCHELoop;
      UsedTCHP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = TCHPLoop;
      UsedSSVHE[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = SSVHELoop;
      UsedSSVHP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = SSVHPLoop;
      UsedCSV[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = CSVLoop;
      UsedJP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = JPLoop;
      UsedJBP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = JBPLoop;
      
      JBPLoop++;
    }
           
    /////////////////////////////////////////
    // Loop on events
    /////////////////////////////////////////  
    for(unsigned int iEvt=0; iEvt<nEvent; iEvt++){
    //for(unsigned int iEvt=0; iEvt<1000; iEvt++){
      
      inWTreeTree->GetEvent(iEvt);
      if(iEvt%10000 == 0)
        std::cout<<"Processing the "<<iEvt<<"th event, time = "<< ((double)clock() - start) / CLOCKS_PER_SEC <<endl;
      
      // PU reweighting???
      float avPU = ( (float)wTree->nPUBXm1() + (float)wTree->nPU() + (float)wTree->nPUBXp1() ) / 3.; // average in 3 BX!!!, as recommended
      
      double lumiWeight3D = 1.0;
      if(!(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")){
	lumiWeight3D = Lumi3DWeights.weight3D(wTree->nPUBXm1(),wTree->nPU(),wTree->nPUBXp1());
      }
      histo1D["lumiWeights3D"]->Fill(lumiWeight3D);
      
      ////////////////////////////////////////////////////////////////
      //    Get variables from wTree needed for Analysis            //
      ////////////////////////////////////////////////////////////////
      //Primary vertices:
      float nPrimaryVertices = wTree->nPV();
          
      //Different bTag values:
      vector<float> btagSSVHE = wTree->bTagSSVHE();
      vector<float> btagSSVHP = wTree->bTagSSVHP();
      vector<float> btagTCHE = wTree->bTagTCHE();
      vector<float> btagTCHP = wTree->bTagTCHP();
      vector<float> btagCSV = wTree->bTagCSV();
      vector<float> btagJP = wTree->bTagJP();
      vector<float>btagJBP = wTree->bTagJBP();
      
      //Generator information:
      int CorrectQuark1 = wTree->hadrLJet1();
      int CorrectQuark2 = wTree->hadrLJet2();
      int CorrectBHadr = wTree->hadrBJet();
      int CorrectBLept = wTree->leptBJet();
      
      //--> Use this to obtain the correct Kinematic Fit index:
      int CorrectKinFitIndex=9999;
      for(int ii=0; ii<12; ii++){
	if( ( CorrectQuark1 == Quark1Index[ii] && CorrectQuark2 == Quark2Index[ii] && CorrectBHadr == BHadrIndex[ii] ) || ( CorrectQuark1 == Quark2Index[ii] && CorrectQuark2 == Quark1Index[ii] && CorrectBHadr == BHadrIndex[ii] ) ){
	  CorrectKinFitIndex = ii;
	}
      }

      //Generator particles:
      TLorentzVector genNeutrino = wTree->standardNeutrino();
      TLorentzVector genLepton = wTree->standardLepton();     
      //TLorentzVector hadrBQuark = wTree->hadrBQuark();
      //TLorentzVector hadrLQuark1 = wTree->hadrLQuark1();
      //TLorentzVector hadrLQuark2 = wTree->hadrLQuark2();
      //TLorentzVector leptBQuark = wTree->leptBQuark();
      
      //Cos theta value on generator level:
      float CosThetaGenerator = wTree->standardCosTheta();

      //Trigger variables:
      int IsoMuTrigger = wTree->isoMuTriggerBool();
      int MuonTriggerCut = wTree->muonTriggerValue();
      int ElectronTriggerCut = wTree->electronTriggerValue();
      //histo2D["IsoMuTrigger"]->Fill(wTree->eventID(),IsoMuTrigger);
      //histo2D["MuonTriggerCut"]->Fill(wTree->eventID(),MuonTriggerCut);
      //histo2D["ElectronTriggerCut"]->Fill(wTree->eventID(),ElectronTriggerCut);
      
      //Delta R variables:
      float DeltaRJetLepton = wTree->deltaRJetLepton();
      float DeltaRMuonJet = wTree->deltaRMuonJet();
      float DeltaRElectronJet = wTree->deltaRElectronJet();
 
      vector<TLorentzVector> selectedJets = wTree->selectedJets();
      TLorentzVector lepton = wTree->lepton(); 
      TLorentzVector MET = wTree->met();
      float relativeLeptonIsolation = wTree->relIso();

      ///////////////////////////////
      //  Calculate scale factors  //
      ///////////////////////////////
      float scaleFactor = 1;
      //semiElectron case !!
      if(semiElectron == true && dataSetName.find("Data") != 0){    
	float ptJet4     = selectedJets[3].Pt();
	float absEtaJet4 = TMath::Abs(selectedJets[3].Eta());
    
	float SFjet,errorUpSFjet,errorDownSFjet;
	float meanEvtSelEffSF_ = 0.974; //(+-0.004 stat)
	float meanTriggerEffSF_ = 0.993; //(+-0.001 stat)
	float factorIso_ = 0.931;
	float factorNonIso_ = 0.069; //We also do a (very small) relIso dependent correction to account for the fact that we have a non-iso ele trigger for the first 6.9% of the dataset (cf. Eq. 4 in the note)
	if(ptJet4 < 40){
	  /// search for the index of the point associated with eta value
	  int iPoint=0;
	  if(absEtaJet4 <xVec[0]-xErrLo[iPoint]) {
	    std::cout<<"WARNING!!! |Eta| point out of range of eff. SF graph. Value set to next possible point" <<std::endl;
	    iPoint=0;
	  }
	  else if(absEtaJet4 >=xVec[N-1]+xErrLo[iPoint]) {
	    std::cout<<"WARNING!!! |Eta| point out of range of eff. SF graph. Value set to next possible point" <<std::endl;
	    iPoint=N-1;
	  }
	  else {
	    for(iPoint=0; iPoint <N; iPoint++) {
	      //std::cout<<"iPoint"<<iPoint << "; xVec[iPoint]="<<xVec[iPoint] << "; xErrLo[iPoint]=" << xErrLo[iPoint] << "; yVec[iPoint]="<<yVec[iPoint] << "; yErrLo[iPoint]=" << yErrLo[iPoint] << "; absEtaJet4="<<absEtaJet4<< std::endl;
	      if(absEtaJet4 <xVec[iPoint]+xErrLo[iPoint]) break;
	    }
	  }
	  SFjet          = yVec[iPoint];
	  errorUpSFjet   = yErrHi[iPoint];
	  errorDownSFjet = yErrLo[iPoint];
	}
	float SFele = meanTriggerEffSF_ * (factorIso_ + factorNonIso_ / (0.997+0.05*relativeLeptonIsolation-3.*relativeLeptonIsolation*relativeLeptonIsolation));
	scaleFactor = 1;//meanEvtSelEffSF_*SFele*SFjet;
      }
      //SemiMuon case !!!
      if(semiMuon == true && dataSetName.find("Data") != 0){
	float SFmuon, errorUpSFmuon, errorDownSFmuon;
	
	/// search for the index of the point associated with eta value
	int iPoint=0;
	if(lepton.Eta() <xVec[0]-xErrLo[iPoint]) {
	  std::cout<<"WARNING!!! Eta point out of range of eff. SF graph. Value set to next possible point" <<std::endl;
	  iPoint=0;
	}
	else if(lepton.Eta() >=xVec[N-1]+xErrLo[iPoint]) {
	  std::cout<<"WARNING!!! Eta point out of range of eff. SF graph. Value set to next possible point" <<std::endl;
	  iPoint=N-1;
	}
	else {
	  for(iPoint=0; iPoint <N; iPoint++) {
	    //std::cout<<"iPoint"<<iPoint << "; xVec[iPoint]="<<xVec[iPoint] << "; xErrLo[iPoint]=" << xErrLo[iPoint] << "; yVec[iPoint]="<<yVec[iPoint] << "; yErrLo[iPoint]=" << yErrLo[iPoint] << "; eta="<<lepton.Eta()<< std::endl;
	    if(lepton.Eta() <xVec[iPoint]+xErrLo[iPoint]) break;
	  }
	}
	SFmuon    = yVec[iPoint];
	errorUpSFmuon   = yErrHi[iPoint];
	errorDownSFmuon = yErrLo[iPoint];	
	
	scaleFactor = SFmuon;
	if(dataSetName.find("TriggEvtSelMinus") == true){
	  scaleFactor = SFmuon-errorDownSFmuon;
	}
	else if(dataSetName.find("TriggEvtSelPlus") == true){
	  scaleFactor = SFmuon+errorUpSFmuon;
	}
      }

      //////////////////////////////////////////////////////////////////////////////////////////////////
      if(CreateMSPlots == true){
	MSPlot["DeltaRJetLepton"]->Fill(DeltaRJetLepton,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D); 
	MSPlot["DeltaRMuonJet"]->Fill(DeltaRMuonJet,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["DeltaRElectronJet"]->Fill(DeltaRElectronJet, datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
      }                      
      
      if(dataSetName.find("Data") == 0) float MassTop = 173.1;
      else float MassTop = 172.5;
      
      for(int ii = 0; ii < 12; ii++){	
	float CosThetaOrigKins = bTagCosThetaCalculation.CalcOrigKins(BLeptIndex[ii], BHadrIndex[ii],  lepton, selectedJets, 80.4, MassTop);
	if(CosThetaOrigKins != 999){
	  TLorentzVector neutrinoOrigKins = bTagCosThetaCalculation.GetNeutrino();

	  TLorentzVector WLeptOrigKins = (neutrinoOrigKins+lepton);
	  if(WLeptOrigKins.E() < 0.){
	    cout << " Considering event with negative WLept energy !!! ( " << WLeptOrigKins.E() << " ) " << endl;
	  }
	}
      }	  

      //Chi squared values for different KinFit configurations:  --> Original kinematics Chi Squared with calculation class !!
      vector<float> ChiSqKinFitHadr;
      vector<float> KinFitProbHadr;
      vector<float> ChiSqKinFitHadrAndLeptWOnly;
      vector<float> KinFitProbHadrAndLeptWOnly;
      vector<float> ChiSqKinFitHadrAndLept;
      vector<float> KinFitProbHadrAndLept;
      
      //Original kinematics initializing:
      vector<TLorentzVector> leptBJetOrigKins, hadrBJetOrigKins, light1JetOrigKins, light2JetOrigKins, leptonOrigKins, neutrinoOrigKins;  //Get this with different method??
      //Original kinematics should not change when different KinFit configuration is used to determine ChiSquared value !! 
      //  --> Only ChiSquared value should change !!

      //Changed kinematics initializing:
      vector<TLorentzVector> leptonChangedHadrAndLeptWOnly, neutrinoChangedHadrAndLeptWOnly;
      vector<TLorentzVector> leptBChangedHadrAndLept, leptonChangedHadrAndLept, neutrinoChangedHadrAndLept;
      vector<TLorentzVector> hadrBChangedHadr, light1ChangedHadr, light2ChangedHadr; //-->Only need this ones for W and top mass distribution after KinFit !!
      vector<TLorentzVector> hadrBChangedHadrAndLeptWOnly, light1ChangedHadrAndLeptWOnly, light2ChangedHadrAndLeptWOnly;
      vector<TLorentzVector> hadrBChangedHadrAndLept, light1ChangedHadrAndLept, light2ChangedHadrAndLept;

      for(unsigned int iCombi = 0; iCombi<12;iCombi++){	
	
	//Hadronic KinFit   !!Has no changed neutrino, lepton or leptonic bJet component
	ChiSqKinFitHadr.push_back(wTree->hadrKinFitChiSq(iCombi)); 	
	if(ChiSqKinFitHadr[iCombi] != 9999) KinFitProbHadr.push_back(TMath::Prob(ChiSqKinFitHadr[iCombi],2));
	else KinFitProbHadr.push_back(-1.);
	// hadrBChangedHadr.push_back(wTree->hadrKinFitHadrB(iCombi));
	// light1ChangedHadr.push_back(wTree->hadrKinFitLight1(iCombi));
	// light2ChangedHadr.push_back(wTree->hadrKinFitLight2(iCombi));
	
	//Hadronic and leptonic W only KinFit     !!Has no changed leptonic bJet component
	ChiSqKinFitHadrAndLeptWOnly.push_back(wTree->hadrAndLeptWOnlyKinFitChiSq(iCombi));
	if(ChiSqKinFitHadrAndLeptWOnly[iCombi] != 9999) KinFitProbHadrAndLeptWOnly.push_back(TMath::Prob(ChiSqKinFitHadrAndLeptWOnly[iCombi],2));
	else KinFitProbHadrAndLeptWOnly.push_back(-1.);
	leptonChangedHadrAndLeptWOnly.push_back(wTree->hadrAndLeptWOnlyKinFitLepton(iCombi));
	neutrinoChangedHadrAndLeptWOnly.push_back(wTree->hadrAndLeptWOnlyKinFitNeutrino(iCombi));
	// hadrBChangedHadrAndLeptWOnly.push_back(wTree->hadrAndLeptWOnlyKinFitHadrB(iCombi));
	// light1ChangedHadrAndLeptWOnly.push_back(wTree->hadrAndLeptWOnlyKinFitLight1(iCombi));
	// light2ChangedHadrAndLeptWOnly.push_back(wTree->hadrAndLeptWOnlyKinFitLight2(iCombi));

	//Hadronic and leptonic KinFit
	ChiSqKinFitHadrAndLept.push_back(wTree->hadrAndLeptKinFitChiSq(iCombi));
	if(ChiSqKinFitHadrAndLept[iCombi] != 9999) KinFitProbHadrAndLept.push_back(TMath::Prob(ChiSqKinFitHadr[iCombi],2));
	else KinFitProbHadrAndLept.push_back(-1.);
	leptBChangedHadrAndLept.push_back(wTree->hadrAndLeptKinFitLeptB(iCombi));
	leptonChangedHadrAndLept.push_back(wTree->hadrAndLeptKinFitLepton(iCombi));
	neutrinoChangedHadrAndLept.push_back(wTree->hadrAndLeptKinFitNeutrino(iCombi));

	if((wTree->hadrAndLeptKinFitLepton(iCombi)+wTree->hadrAndLeptKinFitNeutrino(iCombi)).E() < 0.){	 
	  //cout << " Combination " << iCombi << " gives negative leptonic energy ! Corresponding chi squared = " << wTree->hadrAndLeptKinFitChiSq(iCombi) << endl;
	}
	hadrBChangedHadrAndLept.push_back(wTree->hadrAndLeptKinFitHadrB(iCombi));
	light1ChangedHadrAndLept.push_back(wTree->hadrAndLeptKinFitLight1(iCombi));
	light2ChangedHadrAndLept.push_back(wTree->hadrAndLeptKinFitLight2(iCombi));       
		
      }//End loop for filling KinFit stuff!!	  
      
      if(CreateMSPlots == true){
	MSPlot["ChiSqHadrComb0"]->Fill(ChiSqKinFitHadr[0], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["ChiSqHadrComb1"]->Fill(ChiSqKinFitHadr[1], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["ChiSqHadrComb2"]->Fill(ChiSqKinFitHadr[2], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["ChiSqHadrComb3"]->Fill(ChiSqKinFitHadr[3], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["ChiSqHadrComb4"]->Fill(ChiSqKinFitHadr[4], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["ChiSqHadrComb5"]->Fill(ChiSqKinFitHadr[5], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["ChiSqHadrComb6"]->Fill(ChiSqKinFitHadr[6], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["ChiSqHadrComb7"]->Fill(ChiSqKinFitHadr[7], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["ChiSqHadrComb8"]->Fill(ChiSqKinFitHadr[8], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["ChiSqHadrComb9"]->Fill(ChiSqKinFitHadr[9], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["ChiSqHadrComb10"]->Fill(ChiSqKinFitHadr[10], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["ChiSqHadrComb11"]->Fill(ChiSqKinFitHadr[11], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);

	MSPlot["KinFitProbabilityHadrComb0"]->Fill(KinFitProbHadr[0], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["KinFitProbabilityHadrComb1"]->Fill(KinFitProbHadr[1], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["KinFitProbabilityHadrComb2"]->Fill(KinFitProbHadr[2], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["KinFitProbabilityHadrComb3"]->Fill(KinFitProbHadr[3], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["KinFitProbabilityHadrComb4"]->Fill(KinFitProbHadr[4], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["KinFitProbabilityHadrComb5"]->Fill(KinFitProbHadr[5], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["KinFitProbabilityHadrComb6"]->Fill(KinFitProbHadr[6], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["KinFitProbabilityHadrComb7"]->Fill(KinFitProbHadr[7], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["KinFitProbabilityHadrComb8"]->Fill(KinFitProbHadr[8], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["KinFitProbabilityHadrComb9"]->Fill(KinFitProbHadr[9], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["KinFitProbabilityHadrComb10"]->Fill(KinFitProbHadr[10], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["KinFitProbabilityHadrComb11"]->Fill(KinFitProbHadr[11], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);

	MSPlot["ChiSqHadrAndLeptWOnlyComb0"]->Fill(ChiSqKinFitHadrAndLeptWOnly[0], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["ChiSqHadrAndLeptWOnlyComb1"]->Fill(ChiSqKinFitHadrAndLeptWOnly[1], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["ChiSqHadrAndLeptWOnlyComb2"]->Fill(ChiSqKinFitHadrAndLeptWOnly[2], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["ChiSqHadrAndLeptWOnlyComb3"]->Fill(ChiSqKinFitHadrAndLeptWOnly[3], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["ChiSqHadrAndLeptWOnlyComb4"]->Fill(ChiSqKinFitHadrAndLeptWOnly[4], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["ChiSqHadrAndLeptWOnlyComb5"]->Fill(ChiSqKinFitHadrAndLeptWOnly[5], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["ChiSqHadrAndLeptWOnlyComb6"]->Fill(ChiSqKinFitHadrAndLeptWOnly[6], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["ChiSqHadrAndLeptWOnlyComb7"]->Fill(ChiSqKinFitHadrAndLeptWOnly[7], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["ChiSqHadrAndLeptWOnlyComb8"]->Fill(ChiSqKinFitHadrAndLeptWOnly[8], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["ChiSqHadrAndLeptWOnlyComb9"]->Fill(ChiSqKinFitHadrAndLeptWOnly[9], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["ChiSqHadrAndLeptWOnlyComb10"]->Fill(ChiSqKinFitHadrAndLeptWOnly[10], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["ChiSqHadrAndLeptWOnlyComb11"]->Fill(ChiSqKinFitHadrAndLeptWOnly[11], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);

	MSPlot["KinFitProbabilityHadrAndLeptWOnlyComb0"]->Fill(KinFitProbHadrAndLeptWOnly[0], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["KinFitProbabilityHadrAndLeptWOnlyComb1"]->Fill(KinFitProbHadrAndLeptWOnly[1], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["KinFitProbabilityHadrAndLeptWOnlyComb2"]->Fill(KinFitProbHadrAndLeptWOnly[2], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["KinFitProbabilityHadrAndLeptWOnlyComb3"]->Fill(KinFitProbHadrAndLeptWOnly[3], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["KinFitProbabilityHadrAndLeptWOnlyComb4"]->Fill(KinFitProbHadrAndLeptWOnly[4], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["KinFitProbabilityHadrAndLeptWOnlyComb5"]->Fill(KinFitProbHadrAndLeptWOnly[5], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["KinFitProbabilityHadrAndLeptWOnlyComb6"]->Fill(KinFitProbHadrAndLeptWOnly[6], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["KinFitProbabilityHadrAndLeptWOnlyComb7"]->Fill(KinFitProbHadrAndLeptWOnly[7], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["KinFitProbabilityHadrAndLeptWOnlyComb8"]->Fill(KinFitProbHadrAndLeptWOnly[8], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["KinFitProbabilityHadrAndLeptWOnlyComb9"]->Fill(KinFitProbHadrAndLeptWOnly[9], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["KinFitProbabilityHadrAndLeptWOnlyComb10"]->Fill(KinFitProbHadrAndLeptWOnly[10], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["KinFitProbabilityHadrAndLeptWOnlyComb11"]->Fill(KinFitProbHadrAndLeptWOnly[11], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);

	MSPlot["ChiSqHadrAndLeptComb0"]->Fill(ChiSqKinFitHadrAndLept[0], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["ChiSqHadrAndLeptComb1"]->Fill(ChiSqKinFitHadrAndLept[1], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["ChiSqHadrAndLeptComb2"]->Fill(ChiSqKinFitHadrAndLept[2], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["ChiSqHadrAndLeptComb3"]->Fill(ChiSqKinFitHadrAndLept[3], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["ChiSqHadrAndLeptComb4"]->Fill(ChiSqKinFitHadrAndLept[4], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["ChiSqHadrAndLeptComb5"]->Fill(ChiSqKinFitHadrAndLept[5], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["ChiSqHadrAndLeptComb6"]->Fill(ChiSqKinFitHadrAndLept[6], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["ChiSqHadrAndLeptComb7"]->Fill(ChiSqKinFitHadrAndLept[7], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["ChiSqHadrAndLeptComb8"]->Fill(ChiSqKinFitHadrAndLept[8], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["ChiSqHadrAndLeptComb9"]->Fill(ChiSqKinFitHadrAndLept[9], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["ChiSqHadrAndLeptComb10"]->Fill(ChiSqKinFitHadrAndLept[10], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["ChiSqHadrAndLeptComb11"]->Fill(ChiSqKinFitHadrAndLept[11], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);

	MSPlot["KinFitProbabilityHadrAndLeptComb0"]->Fill(KinFitProbHadrAndLept[0], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["KinFitProbabilityHadrAndLeptComb1"]->Fill(KinFitProbHadrAndLept[1], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["KinFitProbabilityHadrAndLeptComb2"]->Fill(KinFitProbHadrAndLept[2], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["KinFitProbabilityHadrAndLeptComb3"]->Fill(KinFitProbHadrAndLept[3], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["KinFitProbabilityHadrAndLeptComb4"]->Fill(KinFitProbHadrAndLept[4], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["KinFitProbabilityHadrAndLeptComb5"]->Fill(KinFitProbHadrAndLept[5], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["KinFitProbabilityHadrAndLeptComb6"]->Fill(KinFitProbHadrAndLept[6], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["KinFitProbabilityHadrAndLeptComb7"]->Fill(KinFitProbHadrAndLept[7], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["KinFitProbabilityHadrAndLeptComb8"]->Fill(KinFitProbHadrAndLept[8], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["KinFitProbabilityHadrAndLeptComb9"]->Fill(KinFitProbHadrAndLept[9], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["KinFitProbabilityHadrAndLeptComb10"]->Fill(KinFitProbHadrAndLept[10], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["KinFitProbabilityHadrAndLeptComb11"]->Fill(KinFitProbHadrAndLept[11], datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
      }
  
      //ooooooooOOOOOOOOOOOOooooooooooooOOOOOOOOOOOOooooooooooooOOOOO
      //ooOOooOOoo      Reading out nTuples done           ooOOooOOoo
      //ooOOooOOoo-----------------------------------------ooOOooOOoo
      //ooOOooOOoo      Start of actual analysis           ooOOooOOoo
      //ooooooooOOOOOOOOOOOOooooooooooooOOOOOOOOOOOOooooooooooooOOOOO
      float TransverseMass = sqrt(2*(abs(lepton.Pt()))*abs(MET.Pt())*(1-cos(lepton.DeltaPhi(MET))));	//Should this not be updated to correct lepton and MET??
      
      if(CreateMSPlots == true){
	MSPlot["TransverseMassBeforeCut"]->Fill(TransverseMass, datasets[iDataSet], true, Luminosity);		
	MSPlot["RelIsoLepton"]->Fill(relativeLeptonIsolation, datasets[iDataSet], true, Luminosity);
      }
      
      //----------------------------------
      //     Require some extra cuts:
      //----------------------------------
      bool eventSelected = false;

      if((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ) NumberEventsBeforeCuts++;
      if(dataSetName.find("Data") ==0) NumberDataEventsBeforeCuts++;
      
      if(CreateMSPlots == true)
	MSPlot["nPVBeforeCuts"]->Fill(nPrimaryVertices,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
      selecTableMacro.Fill(iDataSet,0,scaleFactor*lumiWeight3D);

      //Apply extra event selections:
      if(TransverseMass > 30 && DeltaRJetLepton > 0.3 && ((DeltaRMuonJet > 0.3 && semiMuon == true) || (DeltaRElectronJet > 0.3 && semiElectron == true))){
	if((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ) NumberSelectedEvents++;
	if(dataSetName.find("Data") ==0) NumberSelectedDataEvents++;
	
	if(CreateMSPlots == true){
	  MSPlot["TransverseMassAfterCut"]->Fill(TransverseMass, datasets[iDataSet], true, Luminosity);
	  MSPlot["nPVAfterTransverseMassCut"]->Fill(nPrimaryVertices,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	  MSPlot["Jet1Pt"]->Fill(selectedJets[0].Pt(), datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);    
	  MSPlot["Jet2Pt"]->Fill(selectedJets[1].Pt(), datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 
	  MSPlot["Jet3Pt"]->Fill(selectedJets[2].Pt(), datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 
	  MSPlot["Jet4Pt"]->Fill(selectedJets[3].Pt(), datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 	  
	}

	//////////////////////////////////////////////////////////////////////
	//  No bTag discriminator combined with jet-combination selection   //
	//////////////////////////////////////////////////////////////////////	
	//Select lowest chi squared value for HadrAndLept KinFit configuration:
	int SelectedJetCombinationHadrAndLept=999;
	float LowestChiSq=9999;
	for(int ii=0;ii<12;ii++){
	  //cout << ii << " ) Lowest Chi Sq : " << LowestChiSq << " versus ChiSqHadrAndLept : " << ChiSqKinFitHadrAndLept[ii] << endl;
	  if(LowestChiSq > ChiSqKinFitHadrAndLept[ii]){
	    LowestChiSq = ChiSqKinFitHadrAndLept[ii];
	    SelectedJetCombinationHadrAndLept=ii;
	  }
	}

	//Original kinematics plots for comparison:
	if(CreateMSPlots == true)
	  MSPlot["JetCombinationChiSqOnly"]->Fill(SelectedJetCombinationHadrAndLept,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	NumberEventsChiSqOnly[iDataSet]++;
	
	if(SelectedJetCombinationHadrAndLept!=999){
	  NumberEventsChiSqOnlyJetCombFound[iDataSet]++;
	  //if(CreateMSPlots == true){
 	  //MSPlot["PtLeptonicBOriginalChiSqOnly"]->Fill(leptBJetOrigKins[SelectedJetCombinationHadrAndLept].Pt(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
 	  //MSPlot["PtHadronicBOriginalChiSqOnly"]->Fill(hadrBJetOrigKins[SelectedJetCombinationHadrAndLept].Pt(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
 	  //MSPlot["PtLight1JetOriginalChiSqOnly"]->Fill(light1JetOrigKins[SelectedJetCombinationHadrAndLept].Pt(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
 	  //MSPlot["PtLight2JetOriginalChiSqOnly"]->Fill(light2JetOrigKins[SelectedJetCombinationHadrAndLept].Pt(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	  //}

	  //Calculate cos theta* (using changed kinematics after KinFit)
	  float CosThetaCalculation=bTagCosThetaCalculation.Calculation(leptonChangedHadrAndLept[SelectedJetCombinationHadrAndLept],neutrinoChangedHadrAndLept[SelectedJetCombinationHadrAndLept],leptBChangedHadrAndLept[SelectedJetCombinationHadrAndLept]);
	    
	  TLorentzVector WLeptonicChiSqOnly = (leptonChangedHadrAndLept[SelectedJetCombinationHadrAndLept]+neutrinoChangedHadrAndLept[SelectedJetCombinationHadrAndLept]);

	  if(CosThetaCalculation != 999){
	    if(WLeptonicChiSqOnly.E() < 0.){
	      cout << " Processing event : " << iEvt << " with cos theta = " << CosThetaCalculation << " (Chi Sq only case )  --   Used combination is : " << SelectedJetCombinationHadrAndLept << " with chi squared value = " << LowestChiSq << " and kinFit probability : " <<  TMath::Prob(LowestChiSq,2) << endl;
	    }
	    	    
	    NumberEventsChiSqOnlyCosThetaFound[iDataSet]++;
	    
	    if(CreateMSPlots == true){
	      //MSPlots for no b-tag case: 	   
	      MSPlot["JetCombCosThetaFoundChiSqOnlyNoBTag"]->Fill(SelectedJetCombinationHadrAndLept,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	      MSPlot["ChiSqKinFitCosThetaFoundChiSqOnlyNoBTag"]->Fill(LowestChiSq,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	      MSPlot["ProbKinFitCosThetaFoundChiSqOnlyNoBTag"]->Fill(TMath::Prob(LowestChiSq,2),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	      
	      MSPlot["nPVChiSqOnlyNoBTag"]->Fill(nPrimaryVertices,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	      MSPlot["TransverseMassChiSqOnlyNoBTag"]->Fill(TransverseMass,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	      MSPlot["TransverseMassSelectedParticlesChiSqOnlyNoBTag"]->Fill(sqrt(2*(abs(leptonChangedHadrAndLept[SelectedJetCombinationHadrAndLept].Pt()))*abs(neutrinoChangedHadrAndLept[SelectedJetCombinationHadrAndLept].Pt())*(1-cos(leptonChangedHadrAndLept[SelectedJetCombinationHadrAndLept].DeltaPhi(neutrinoChangedHadrAndLept[SelectedJetCombinationHadrAndLept])))),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	      MSPlot["CosThetaChiSqOnlyNoBTag"]->Fill(CosThetaCalculation,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	      
	      MSPlot["Jet1PtChiSqOnlyNoBTag"]->Fill(selectedJets[0].Pt(), datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);    
	      MSPlot["Jet2PtChiSqOnlyNoBTag"]->Fill(selectedJets[1].Pt(), datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 
	      MSPlot["Jet3PtChiSqOnlyNoBTag"]->Fill(selectedJets[2].Pt(), datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 
	      MSPlot["Jet4PtChiSqOnlyNoBTag"]->Fill(selectedJets[3].Pt(), datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 	  
	      
	      MSPlot["PtLeptonicBChiSqOnlyNoBTag"]->Fill(leptBChangedHadrAndLept[SelectedJetCombinationHadrAndLept].Pt(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	      MSPlot["PtHadronicBChiSqOnlyNoBTag"]->Fill(hadrBChangedHadrAndLept[SelectedJetCombinationHadrAndLept].Pt(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	      MSPlot["PtLightJet1ChiSqOnlyNoBTag"]->Fill(light1ChangedHadrAndLept[SelectedJetCombinationHadrAndLept].Pt(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	      MSPlot["PtLightJet2ChiSqOnlyNoBTag"]->Fill(light2ChangedHadrAndLept[SelectedJetCombinationHadrAndLept].Pt(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	      
	      MSPlot["leptonPhiChiSqOnlyNoBTag"]->Fill(leptonChangedHadrAndLept[SelectedJetCombinationHadrAndLept].Phi(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	      MSPlot["leptonPxChiSqOnlyNoBTag"]->Fill(leptonChangedHadrAndLept[SelectedJetCombinationHadrAndLept].Px(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	      MSPlot["leptonPyChiSqOnlyNoBTag"]->Fill(leptonChangedHadrAndLept[SelectedJetCombinationHadrAndLept].Py(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	      MSPlot["leptonPzChiSqOnlyNoBTag"]->Fill(leptonChangedHadrAndLept[SelectedJetCombinationHadrAndLept].Pz(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	      MSPlot["leptonPtChiSqOnlyNoBTag"]->Fill(leptonChangedHadrAndLept[SelectedJetCombinationHadrAndLept].Pt(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	      MSPlot["leptonMassChiSqOnlyNoBTag"]->Fill(leptonChangedHadrAndLept[SelectedJetCombinationHadrAndLept].M(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	      
	      MSPlot["METPhiChiSqOnlyNoBTag"]->Fill(neutrinoChangedHadrAndLept[SelectedJetCombinationHadrAndLept].Phi(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	      MSPlot["METPxChiSqOnlyNoBTag"]->Fill(neutrinoChangedHadrAndLept[SelectedJetCombinationHadrAndLept].Px(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	      MSPlot["METPyChiSqOnlyNoBTag"]->Fill(neutrinoChangedHadrAndLept[SelectedJetCombinationHadrAndLept].Py(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	      MSPlot["METPzChiSqOnlyNoBTag"]->Fill(neutrinoChangedHadrAndLept[SelectedJetCombinationHadrAndLept].Pz(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	      MSPlot["METPtChiSqOnlyNoBTag"]->Fill(neutrinoChangedHadrAndLept[SelectedJetCombinationHadrAndLept].Pt(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	      MSPlot["METMassChiSqOnlyNoBTag"]->Fill(neutrinoChangedHadrAndLept[SelectedJetCombinationHadrAndLept].M(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	    }
	    
	    //Consider specific b-tag case (SSVHE M lept):
	    if(btagSSVHE[BLeptIndex[SelectedJetCombinationHadrAndLept]]>SSVHEbTag[1]){
	      
	      if(CreateMSPlots == true){
		MSPlot["JetCombCosThetaFoundChiSqOnlySpecificBTag"]->Fill(SelectedJetCombinationHadrAndLept,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
		MSPlot["ChiSqKinFitCosThetaFoundChiSqOnlySpecificBTag"]->Fill(LowestChiSq,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
		MSPlot["ProbKinFitCosThetaFoundChiSqOnlySpecificBTag"]->Fill(TMath::Prob(LowestChiSq,2),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
		
		MSPlot["nPVChiSqOnlySpecificBTag"]->Fill(nPrimaryVertices,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
		MSPlot["TransverseMassChiSqOnlySpecificBTag"]->Fill(TransverseMass,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
		MSPlot["TransverseMassSelectedParticlesChiSqOnlySpecificBTag"]->Fill(sqrt(2*(abs(leptonChangedHadrAndLept[SelectedJetCombinationHadrAndLept].Pt()))*abs(neutrinoChangedHadrAndLept[SelectedJetCombinationHadrAndLept].Pt())*(1-cos(leptonChangedHadrAndLept[SelectedJetCombinationHadrAndLept].DeltaPhi(neutrinoChangedHadrAndLept[SelectedJetCombinationHadrAndLept])))),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
		MSPlot["CosThetaChiSqOnlySpecificBTag"]->Fill(CosThetaCalculation,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
		MSPlot["Jet1PtChiSqOnlySpecificBTag"]->Fill(selectedJets[0].Pt(), datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);    
		MSPlot["Jet2PtChiSqOnlySpecificBTag"]->Fill(selectedJets[1].Pt(), datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 
		MSPlot["Jet3PtChiSqOnlySpecificBTag"]->Fill(selectedJets[2].Pt(), datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 
		MSPlot["Jet4PtChiSqOnlySpecificBTag"]->Fill(selectedJets[3].Pt(), datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 	  
		
		MSPlot["PtLeptonicBChiSqOnlySpecificBTag"]->Fill(leptBChangedHadrAndLept[SelectedJetCombinationHadrAndLept].Pt(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
		MSPlot["PtHadronicBChiSqOnlySpecificBTag"]->Fill(hadrBChangedHadrAndLept[SelectedJetCombinationHadrAndLept].Pt(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
		MSPlot["PtLightJet1ChiSqOnlySpecificBTag"]->Fill(light1ChangedHadrAndLept[SelectedJetCombinationHadrAndLept].Pt(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
		MSPlot["PtLightJet2ChiSqOnlySpecificBTag"]->Fill(light2ChangedHadrAndLept[SelectedJetCombinationHadrAndLept].Pt(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
		
		MSPlot["leptonPhiChiSqOnlySpecificBTag"]->Fill(leptonChangedHadrAndLept[SelectedJetCombinationHadrAndLept].Phi(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
		MSPlot["leptonPxChiSqOnlySpecificBTag"]->Fill(leptonChangedHadrAndLept[SelectedJetCombinationHadrAndLept].Px(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
		MSPlot["leptonPyChiSqOnlySpecificBTag"]->Fill(leptonChangedHadrAndLept[SelectedJetCombinationHadrAndLept].Py(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
		MSPlot["leptonPzChiSqOnlySpecificBTag"]->Fill(leptonChangedHadrAndLept[SelectedJetCombinationHadrAndLept].Pz(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
		MSPlot["leptonPtChiSqOnlySpecificBTag"]->Fill(leptonChangedHadrAndLept[SelectedJetCombinationHadrAndLept].Pt(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
		MSPlot["leptonMassChiSqOnlySpecificBTag"]->Fill(leptonChangedHadrAndLept[SelectedJetCombinationHadrAndLept].M(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
		
		MSPlot["METPhiChiSqOnlySpecificBTag"]->Fill(neutrinoChangedHadrAndLept[SelectedJetCombinationHadrAndLept].Phi(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
		MSPlot["METPxChiSqOnlySpecificBTag"]->Fill(neutrinoChangedHadrAndLept[SelectedJetCombinationHadrAndLept].Px(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
		MSPlot["METPyChiSqOnlySpecificBTag"]->Fill(neutrinoChangedHadrAndLept[SelectedJetCombinationHadrAndLept].Py(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
		MSPlot["METPzChiSqOnlySpecificBTag"]->Fill(neutrinoChangedHadrAndLept[SelectedJetCombinationHadrAndLept].Pz(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
		MSPlot["METPtChiSqOnlySpecificBTag"]->Fill(neutrinoChangedHadrAndLept[SelectedJetCombinationHadrAndLept].Pt(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
		MSPlot["METMassChiSqOnlySpecificBTag"]->Fill(neutrinoChangedHadrAndLept[SelectedJetCombinationHadrAndLept].M(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	      }
	    }//End specific bTag loop
	    
	  }
	}
		
	// //-------------------------------------------------------------------------------
	// //Create some alternative helicites weights to obtain different ratios:
	// //-------------------------------------------------------------------------------
	// float LongitudinalFraction = 0.645167;
	// float LeftHandedFraction = 0.321369;
	// float RightHandedFraction = 0.033464;
	// float TheoreticalDistributionValue = (LongitudinalFraction*6*(1-CosThetaGenerator*CosThetaGenerator) + (1-CosThetaGenerator)*(1-CosThetaGenerator)*3*LeftHandedFraction + RightHandedFraction*3*(1+CosThetaGenerator)*(1+CosThetaGenerator))/8;
	// for(int helicityNumbers=0;helicityNumbers<SizeArray;helicityNumbers++){
	//   if((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ){ 
	//     UsedDistributionValue[helicityNumbers]=(HelicityFraction[helicityNumbers][0]*6*(1-CosThetaGenerator*CosThetaGenerator) + (1-CosThetaGenerator)*(1-CosThetaGenerator)*3*HelicityFraction[helicityNumbers][2] + HelicityFraction[helicityNumbers][1]*3*(1+CosThetaGenerator)*(1+CosThetaGenerator))/8;
	//     HelicityWeight[helicityNumbers]=UsedDistributionValue[helicityNumbers]/TheoreticalDistributionValue;
	//   }
	// }

	//-------------------------------------------------------------
	//Obtain jet combination for the different b-tag constraints:
	//------------------------------------------------------------
	
	//Initialize bTag loop variables:
	int TCHEbTagLoop, TCHPbTagLoop, SSVHEbTagLoop, SSVHPbTagLoop, CSVbTagLoop,JPbTagLoop, JBPbTagLoop;
	int ConsideredBTagger=0; //0=tche, 1 = tchp, 2 = ssvhe, 3 = ssvhp, 4 = csv, 5 = JP & 6 = JBP		
	
	for(JBPbTagLoop =1;JBPbTagLoop< (NumberJBPbTags+1); JBPbTagLoop++){	  
	  for(JPbTagLoop=1;JPbTagLoop <=(NumberJPbTags+1); JPbTagLoop++){
	    if(ConsideredBTagger == 6) JPbTagLoop =14;   //--> Do not run over all possible TCHE bTag values!	    
	    for(CSVbTagLoop =1;CSVbTagLoop<= (NumberCSVbTags+1); CSVbTagLoop++){	  
	      if(  ConsideredBTagger == 5 || ConsideredBTagger == 6) CSVbTagLoop =14;   //--> Do not run over all possible TCHE bTag values!	    
	      for(SSVHPbTagLoop=1;SSVHPbTagLoop <=(NumberSSVHPbTags+1); SSVHPbTagLoop++){
		if(ConsideredBTagger == 4 || ConsideredBTagger == 5 || ConsideredBTagger == 6) SSVHPbTagLoop =14;   //--> Do not run over all possible TCHE bTag values!	    
		for(SSVHEbTagLoop=1; SSVHEbTagLoop<=(NumberSSVHEbTags+1); SSVHEbTagLoop++){
		  if(ConsideredBTagger == 3 || ConsideredBTagger == 4 || ConsideredBTagger == 5 || ConsideredBTagger == 6) SSVHEbTagLoop = 14;	      
		  for(TCHPbTagLoop=1; TCHPbTagLoop <= (NumberTCHPbTags+1); TCHPbTagLoop++){
		    if(ConsideredBTagger == 2 || ConsideredBTagger == 3 || ConsideredBTagger == 4 || ConsideredBTagger == 5 || ConsideredBTagger == 6) TCHPbTagLoop =14; 		
		    for(TCHEbTagLoop=1; TCHEbTagLoop <= (NumberTCHEbTags+1);TCHEbTagLoop++){
		      if(ConsideredBTagger ==1 || ConsideredBTagger ==2 || ConsideredBTagger ==3 || ConsideredBTagger ==4 || ConsideredBTagger == 5 || ConsideredBTagger == 6) TCHEbTagLoop =14;
		      int SumBTag = TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1+JPbTagLoop-1+JBPbTagLoop-1;
		      if(UsedTCHE[SumBTag]==TCHEbTagLoop && UsedTCHP[SumBTag]==TCHPbTagLoop && UsedSSVHE[SumBTag]==SSVHEbTagLoop && UsedSSVHP[SumBTag]==SSVHPbTagLoop && UsedCSV[SumBTag]==CSVbTagLoop && UsedJP[SumBTag]==JPbTagLoop && UsedJBP[SumBTag]==JBPbTagLoop){
			int bTagLoop;
			if(ConsideredBTagger == 0) bTagLoop = TCHEbTagLoop;
			else if(ConsideredBTagger == 1) bTagLoop = TCHPbTagLoop;
			else if(ConsideredBTagger == 2) bTagLoop = SSVHEbTagLoop;
			else if(ConsideredBTagger == 3) bTagLoop = SSVHPbTagLoop;
			else if(ConsideredBTagger == 4) bTagLoop = CSVbTagLoop;
			else if(ConsideredBTagger == 5) bTagLoop = JPbTagLoop;
			else if(ConsideredBTagger == 6) bTagLoop = JBPbTagLoop;
			
			int JetCombHadr = bTagJetSelection.HighestProbSelection(bTagLoop,ConsideredBTagger,KinFitProbHadr,btagTCHE,btagTCHP,btagSSVHE,btagSSVHP,btagCSV,btagJP,btagJBP);
			int JetCombHadrAndLeptWOnly =bTagJetSelection.HighestProbSelection(bTagLoop,ConsideredBTagger,KinFitProbHadrAndLeptWOnly,btagTCHE,btagTCHP,btagSSVHE,btagSSVHP,btagCSV,btagJP,btagJBP);
			int JetCombHadrAndLept = bTagJetSelection.HighestProbSelection(bTagLoop,ConsideredBTagger,KinFitProbHadrAndLept,btagTCHE,btagTCHP,btagSSVHE,btagSSVHP,btagCSV,btagJP,btagJBP);
			  
			if(bTagLoop ==1){
			  
			  if(CreateMSPlots == true){
			    MSPlot["JetCombination"]->Fill(JetCombHadrAndLept,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
			    histo1D["JetCombHadr"]->Fill(JetCombHadr);
			    histo1D["JetCombHadrAndLeptWOnly"]->Fill(JetCombHadrAndLeptWOnly);
			    histo1D["JetCombHadrAndLept"]->Fill(JetCombHadrAndLept);
			  }
			}
			
			// Only hadronic KinFit configuration:
			if(JetCombHadr!=999 && HadronicOnly == true){ 	
			  
			  if(dataSetName.find("Data") == 0) float MassTop = 173.1;
			  else float MassTop = 172.5;
			  
			  float CosThetaCalcHadr = bTagCosThetaCalculation.CalcOrigKins(BLeptIndex[JetCombHadr], BHadrIndex[JetCombHadr],  lepton, selectedJets, 80.4, MassTop);
			  if(CosThetaCalcHadr != 999){			    
			    TLorentzVector neutrinoHadr = bTagCosThetaCalculation.GetNeutrino();

			    if((lepton+neutrinoHadr).E() < 0. ){
			      cout << " Processing event : " << iEvt << " with cos theta = " << CosThetaCalcHadr << " (Chi Sq and BTag case -- Hadr)  --   Used combination is : " << JetCombHadr << " with chi squared value = " << ChiSqKinFitHadr[JetCombHadr] << " and KinFit probability : " << KinFitProbHadr[JetCombHadr] << endl;
			    }
			    
			    if((dataSetName.find("TTbarJets_SemiMu") == 0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") == 0 && semiElectron == true)){
			      histo1D["LeptWMassHadr"]->Fill((neutrinoHadr+lepton).M());
			      histo1D["LeptTopMassHadr"]->Fill((neutrinoHadr+lepton+selectedJets[BLeptIndex[JetCombHadr]]).M());
			    }		      
			    
			    NumberRemainingEventsKinFitHadr[SumBTag][iDataSet]++;
			    CosThetaValuesKinFitHadr[SumBTag].push_back(CosThetaCalcHadr);
			    //LumiWeightVectorKinFitHadr[SumBTag].push_back(lumiWeight3D);
			    EventCorrectionWeightKinFitHadr[SumBTag].push_back(scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor());		      		      
			    if(((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ) && dataSetName.find("JES") != 0){
			      if(BLeptIndex[JetCombHadr] == CorrectBLept){//Count the number of events with correctly reconstructed leptonic b-jet
				NumberBLeptCorrectEventsKinFitHadr[SumBTag][iDataSet]++;
			      }
			      CosThGenKinFitHadr[SumBTag].push_back(CosThetaGenerator);
			    }//End of signal sample
			  }
			}
			
			//Hadronic and leptonic W Only KinFit configuration:
			if(JetCombHadrAndLeptWOnly!=999 && HadronicAndLeptonicWOnly == true){ 

			  float CosThetaCalcHadrAndLeptWOnly=bTagCosThetaCalculation.Calculation(leptonChangedHadrAndLeptWOnly[JetCombHadrAndLeptWOnly],neutrinoChangedHadrAndLeptWOnly[JetCombHadrAndLeptWOnly],selectedJets[BLeptIndex[JetCombHadrAndLeptWOnly]]);
			  if(CosThetaCalcHadrAndLeptWOnly != 999){
			  
			    if((leptonChangedHadrAndLeptWOnly[JetCombHadrAndLeptWOnly]+neutrinoChangedHadrAndLeptWOnly[JetCombHadrAndLeptWOnly]).E() < 0. ){
			      cout << " Processing event : " << iEvt << " with cos theta = " << CosThetaCalcHadrAndLeptWOnly << " (Chi Sq and BTag case -- HadrAndLeptWOnly)  --   Used combination is : " << JetCombHadrAndLeptWOnly << " with chi squared value = " << ChiSqKinFitHadrAndLeptWOnly[JetCombHadrAndLeptWOnly] << " and KinFit probability : " << KinFitProbHadrAndLeptWOnly[JetCombHadrAndLeptWOnly] << endl;
			    }

			    NumberRemainingEventsKinFitHadrAndLeptWOnly[SumBTag][iDataSet]++;
			    CosThetaValuesKinFitHadrAndLeptWOnly[SumBTag].push_back(CosThetaCalcHadrAndLeptWOnly);
			    //LumiWeightVectorKinFitHadrAndLeptWOnly[SumBTag].push_back(lumiWeight3D);
			    EventCorrectionWeightKinFitHadrAndLeptWOnly[SumBTag].push_back(scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor());		      		      

			    if(((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ) && dataSetName.find("JES") != 0){
			      if(BLeptIndex[JetCombHadrAndLeptWOnly] == CorrectBLept){//Count the number of events with correctly reconstructed leptonic b-jet
				NumberBLeptCorrectEventsKinFitHadrAndLeptWOnly[SumBTag][iDataSet]++;
			      }
			      CosThGenKinFitHadrAndLeptWOnly[SumBTag].push_back(CosThetaGenerator);
			    }//End of signal sample
			  }
			}
			
			//Hadronic and leptonic KinFit configuration:
			NumberEventsChiSqAndBTag[iDataSet]++;
			if(JetCombHadrAndLept!=999 && HadronicAndLeptonic == true){ 

			  NumberEventsChiSqAndBTagJetCombFound[iDataSet]++;
			  float CosThetaCalcHadrAndLept=bTagCosThetaCalculation.Calculation(leptonChangedHadrAndLept[JetCombHadrAndLept],neutrinoChangedHadrAndLept[JetCombHadrAndLept],leptBChangedHadrAndLept[JetCombHadrAndLept]);
			  if(CosThetaCalcHadrAndLept != 999){

			    if(((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ) && dataSetName.find("JES") != 0){	
			      NumberEventsSignalCosThetaFound++;
			      if(BLeptIndex[JetCombHadrAndLept] == CorrectBLept)
				NumberEventsSignalCorrectLeptBCosThetaFound++;
			      if(BHadrIndex[JetCombHadrAndLept] == CorrectBHadr)
				NumberEventsSignalCorrectHadrBCosThetaFound++;
			      if( (Quark1Index[JetCombHadrAndLept] == CorrectQuark1 && Quark2Index[JetCombHadrAndLept] == CorrectQuark2) || (Quark1Index[JetCombHadrAndLept] == CorrectQuark2 && Quark2Index[JetCombHadrAndLept] == CorrectQuark1)){
				NumberEventsSignalCorrectLightCosThetaFound++;
				if(BHadrIndex[JetCombHadrAndLept] == CorrectBHadr)
				  NumberEventsSignalCorrectHadronicCosThetaFound++;
			      }
			      if(KinFitProbHadr[JetCombHadrAndLept]<0.01){
				NumberEventsLowProbCosThetaFound++;
				if(BLeptIndex[JetCombHadrAndLept] == CorrectBLept){
				  NumberEventsLowProbCorrectLeptBCosThetaFound++;
				  if(BHadrIndex[JetCombHadrAndLept] == CorrectBHadr)
				    NumberEventsLowProbCorrectHadrBCosThetaFound++;
				  if( (Quark1Index[JetCombHadrAndLept] == CorrectQuark1 && Quark2Index[JetCombHadrAndLept] == CorrectQuark2) || (Quark1Index[JetCombHadrAndLept] == CorrectQuark2 && Quark2Index[JetCombHadrAndLept] == CorrectQuark1)){
				    NumberEventsLowProbCorrectLightCosThetaFound++;
				    if(BHadrIndex[JetCombHadrAndLept] == CorrectBHadr)
				      NumberEventsLowProbCorrectHadronicCosThetaFound++;
				  }				  
				}
			      }
			    }
			    
			    if((leptonChangedHadrAndLept[JetCombHadrAndLept]+neutrinoChangedHadrAndLept[JetCombHadrAndLept]).E() < 0.){
			      cout << " Processing event : " << iEvt << " with cos theta = " << CosThetaCalcHadrAndLept << " (Chi Sq and BTag case -- HadrAndLept)  --   Used combination is : " << JetCombHadrAndLept << " with chi squared value = " << ChiSqKinFitHadrAndLept[JetCombHadrAndLept] << " and KinFit probability : " << KinFitProbHadrAndLept[JetCombHadrAndLept] << endl;
			    }
			    
			    NumberEventsChiSqAndBTagCosThetaFound[iDataSet]++;
			    NumberRemainingEventsKinFitHadrAndLept[SumBTag][iDataSet]++;
			    CosThetaValuesKinFitHadrAndLept[SumBTag].push_back(CosThetaCalcHadrAndLept);
			    //LumiWeightVectorKinFitHadrAndLept[SumBTag].push_back(lumiWeight3D);
			    EventCorrectionWeightKinFitHadrAndLept[SumBTag].push_back(scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor());		      		      

			    if(CreateMSPlots == true){
			      //  MSPlots for no b-tag case:   
			      if(bTagLoop ==1){
				MSPlot["JetCombCosThetaFoundNoBTag"]->Fill(JetCombHadrAndLept,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
				MSPlot["ChiSqKinFitCosThetaFoundNoBTag"]->Fill(ChiSqKinFitHadrAndLept[JetCombHadrAndLept],datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);
				MSPlot["ProbKinFitCosThetaFoundNoBTag"]->Fill(KinFitProbHadrAndLept[JetCombHadrAndLept],datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);
				
				MSPlot["TransverseMassNoBTag"]->Fill(TransverseMass, datasets[iDataSet], true, Luminosity);
				MSPlot["TransverseMassSelectedParticlesNoBTag"]->Fill(sqrt(2*(abs(leptonChangedHadrAndLept[JetCombHadrAndLept].Pt()))*abs(neutrinoChangedHadrAndLept[JetCombHadrAndLept].Pt())*(1-cos(leptonChangedHadrAndLept[JetCombHadrAndLept].DeltaPhi(neutrinoChangedHadrAndLept[JetCombHadrAndLept])))),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
				MSPlot["nPVNoBTag"]->Fill(nPrimaryVertices,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
				MSPlot["CosThetaNoBTag"]->Fill(CosThetaCalcHadrAndLept,datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);
				
				MSPlot["Jet1PtNoBTag"]->Fill(selectedJets[0].Pt(), datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);    
				MSPlot["Jet2PtNoBTag"]->Fill(selectedJets[1].Pt(), datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 
				MSPlot["Jet3PtNoBTag"]->Fill(selectedJets[2].Pt(), datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 
				MSPlot["Jet4PtNoBTag"]->Fill(selectedJets[3].Pt(), datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 	  
				
				MSPlot["METPhiNoBTag"]->Fill(neutrinoChangedHadrAndLept[JetCombHadrAndLept].Phi(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 
				MSPlot["METPxNoBTag"]->Fill(neutrinoChangedHadrAndLept[JetCombHadrAndLept].Px(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 
				MSPlot["METPyNoBTag"]->Fill(neutrinoChangedHadrAndLept[JetCombHadrAndLept].Py(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 
				MSPlot["METPzNoBTag"]->Fill(neutrinoChangedHadrAndLept[JetCombHadrAndLept].Pz(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 
				MSPlot["METPtNoBTag"]->Fill(neutrinoChangedHadrAndLept[JetCombHadrAndLept].Pt(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 
				MSPlot["METMassNoBTag"]->Fill(neutrinoChangedHadrAndLept[JetCombHadrAndLept].M(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 
				
				MSPlot["leptonPhiNoBTag"]->Fill(leptonChangedHadrAndLept[JetCombHadrAndLept].Phi(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 
				MSPlot["leptonPxNoBTag"]->Fill(leptonChangedHadrAndLept[JetCombHadrAndLept].Px(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 
				MSPlot["leptonPyNoBTag"]->Fill(leptonChangedHadrAndLept[JetCombHadrAndLept].Py(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 
				MSPlot["leptonPzNoBTag"]->Fill(leptonChangedHadrAndLept[JetCombHadrAndLept].Pz(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 
				MSPlot["leptonPtNoBTag"]->Fill(leptonChangedHadrAndLept[JetCombHadrAndLept].Pt(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 
				MSPlot["leptonMassNoBTag"]->Fill(leptonChangedHadrAndLept[JetCombHadrAndLept].M(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);
				
				MSPlot["PtLeptonicBNoBTag"]->Fill(leptBChangedHadrAndLept[JetCombHadrAndLept].Pt(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);
				MSPlot["PtHadronicBNoBTag"]->Fill(hadrBChangedHadrAndLept[JetCombHadrAndLept].Pt(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);
				MSPlot["PtLightJet1NoBTag"]->Fill(light1ChangedHadrAndLept[JetCombHadrAndLept].Pt(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);
				MSPlot["PtLightJet2NoBTag"]->Fill(light2ChangedHadrAndLept[JetCombHadrAndLept].Pt(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);

				MSPlot["HadronicWMassNoBTag"]->Fill((light1ChangedHadrAndLept[JetCombHadrAndLept]+light2ChangedHadrAndLept[JetCombHadrAndLept]).M(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);
				MSPlot["HadronicTopMassNoBTag"]->Fill((hadrBChangedHadrAndLept[JetCombHadrAndLept]+light1ChangedHadrAndLept[JetCombHadrAndLept]+light2ChangedHadrAndLept[JetCombHadrAndLept]).M(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);
				MSPlot["LeptonicWMassNoBTag"]->Fill((leptonChangedHadrAndLept[JetCombHadrAndLept]+neutrinoChangedHadrAndLept[JetCombHadrAndLept]).M(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);
				MSPlot["LeptonicTopMassNoBTag"]->Fill((leptonChangedHadrAndLept[JetCombHadrAndLept]+neutrinoChangedHadrAndLept[JetCombHadrAndLept]+leptBChangedHadrAndLept[JetCombHadrAndLept]).M(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);
			      }			    
			    }

			    //Check influence of used binnumber:
			    if(ConsideredTagger == 6 && bTagLoop == 9){
			      //cout << " Considered bTagger : " << PresentationOutput[SumBTag] << endl;
			      MSPlot["CosThetaSpecificDoubleBTag50Bins"]->Fill(CosThetaCalcHadrAndLept,datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);
			      MSPlot["CosThetaSpecificDoubleBTag30Bins"]->Fill(CosThetaCalcHadrAndLept,datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);
			    }
			    
			    //  MSPlots for specific b-tag case:  
			    int OptimalTagger, OptimalbTagLoop;
			    if(semiMuon == true){
			      OptimalTagger = 6;
			      OptimalbTagLoop = 7;
			    }
			    else if(semiElectron == true){
			      OptimalTagger = 1
			      OptimalbTagLoop = 1;
			    }
			    if(ConsideredBTagger == OptimalTagger && bTagLoop == OptimalbTagLoop){
			      //cout << " Considered bTagger : " << PresentationOutput[SumBTag] << endl;
			      NumberEventsChiSqAndBTagCosThetaFound[iDataSet]++;
				
			      if(CreateMSPlots == true){
				MSPlot["TransverseMassSpecificBTag"]->Fill(TransverseMass, datasets[iDataSet], true, Luminosity);
				MSPlot["TransverseMassSelectedParticlesSpecificBTag"]->Fill(sqrt(2*(abs(leptonChangedHadrAndLept[JetCombHadrAndLept].Pt()))*abs(neutrinoChangedHadrAndLept[JetCombHadrAndLept].Pt())*(1-cos(leptonChangedHadrAndLept[JetCombHadrAndLept].DeltaPhi(neutrinoChangedHadrAndLept[JetCombHadrAndLept])))),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
				MSPlot["nPVSpecificBTag"]->Fill(nPrimaryVertices,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
				MSPlot["Jet1PtSpecificBTag"]->Fill(selectedJets[0].Pt(), datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);    
				MSPlot["Jet2PtSpecificBTag"]->Fill(selectedJets[1].Pt(), datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 
				MSPlot["Jet3PtSpecificBTag"]->Fill(selectedJets[2].Pt(), datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 
				MSPlot["Jet4PtSpecificBTag"]->Fill(selectedJets[3].Pt(), datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 
				
				//Met and lepton info:
				MSPlot["METPhiSpecificBTag"]->Fill(neutrinoChangedHadrAndLept[JetCombHadrAndLept].Phi(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 
				MSPlot["METPxSpecificBTag"]->Fill(neutrinoChangedHadrAndLept[JetCombHadrAndLept].Px(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 
				MSPlot["METPySpecificBTag"]->Fill(neutrinoChangedHadrAndLept[JetCombHadrAndLept].Py(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 
				MSPlot["METPzSpecificBTag"]->Fill(neutrinoChangedHadrAndLept[JetCombHadrAndLept].Pz(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 
				MSPlot["METPtSpecificBTag"]->Fill(neutrinoChangedHadrAndLept[JetCombHadrAndLept].Pt(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 
				MSPlot["METMassSpecificBTag"]->Fill(neutrinoChangedHadrAndLept[JetCombHadrAndLept].M(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 
				
				MSPlot["leptonPhiSpecificBTag"]->Fill(leptonChangedHadrAndLept[JetCombHadrAndLept].Phi(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 
				MSPlot["leptonPxSpecificBTag"]->Fill(leptonChangedHadrAndLept[JetCombHadrAndLept].Px(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 
				MSPlot["leptonPySpecificBTag"]->Fill(leptonChangedHadrAndLept[JetCombHadrAndLept].Py(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 
				MSPlot["leptonPzSpecificBTag"]->Fill(leptonChangedHadrAndLept[JetCombHadrAndLept].Pz(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 
				MSPlot["leptonPtSpecificBTag"]->Fill(leptonChangedHadrAndLept[JetCombHadrAndLept].Pt(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 
				MSPlot["leptonMassSpecificBTag"]->Fill(leptonChangedHadrAndLept[JetCombHadrAndLept].M(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);
			      
				//KinFit information:
				MSPlot["ProbKinFitCosThetaFoundSpecificBTag"]->Fill(KinFitProbHadrAndLept[JetCombHadrAndLept],datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);
				MSPlot["ChiSqKinFitCosThetaFoundSpecificBTag"]->Fill(ChiSqKinFitHadrAndLept[JetCombHadrAndLept],datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);
				MSPlot["JetCombCosThetaFoundSpecificBTag"]->Fill(JetCombHadrAndLept,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
				
				//Cos theta:
				MSPlot["CosThetaSpecificBTag50Bins"]->Fill(CosThetaCalcHadrAndLept,datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);
				MSPlot["CosThetaSpecificBTag30Bins"]->Fill(CosThetaCalcHadrAndLept,datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);
				
				//Particle Kinematics:
				MSPlot["PtLeptonicBSpecificBTag"]->Fill(leptBChangedHadrAndLept[JetCombHadrAndLept].Pt(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);
				MSPlot["PtHadronicBSpecificBTag"]->Fill(hadrBChangedHadrAndLept[JetCombHadrAndLept].Pt(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);
				MSPlot["PtLightJet1SpecificBTag"]->Fill(light1ChangedHadrAndLept[JetCombHadrAndLept].Pt(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);
				MSPlot["PtLightJet2SpecificBTag"]->Fill(light2ChangedHadrAndLept[JetCombHadrAndLept].Pt(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);
				
				MSPlot["HadronicWMassSpecificBTag"]->Fill((light1ChangedHadrAndLept[JetCombHadrAndLept]+light2ChangedHadrAndLept[JetCombHadrAndLept]).M(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);
				MSPlot["HadronicTopMassSpecificBTag"]->Fill((hadrBChangedHadrAndLept[JetCombHadrAndLept]+light1ChangedHadrAndLept[JetCombHadrAndLept]+light2ChangedHadrAndLept[JetCombHadrAndLept]).M(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);
				MSPlot["LeptonicWMassSpecificBTag"]->Fill((leptonChangedHadrAndLept[JetCombHadrAndLept]+neutrinoChangedHadrAndLept[JetCombHadrAndLept]).M(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);
				MSPlot["LeptonicTopMassSpecificBTag"]->Fill((leptonChangedHadrAndLept[JetCombHadrAndLept]+neutrinoChangedHadrAndLept[JetCombHadrAndLept]+leptBChangedHadrAndLept[JetCombHadrAndLept]).M(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);

			      }
			    }//End of specific bTag loop
			    
			    if(((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ) && dataSetName.find("JES") != 0){
			      if(BLeptIndex[JetCombHadrAndLept] == CorrectBLept){//Count the number of events with correctly reconstructed leptonic b-jet
				NumberBLeptCorrectEventsKinFitHadrAndLept[SumBTag][iDataSet]++;
			      }
			      CosThGenKinFitHadrAndLept[SumBTag].push_back(CosThetaGenerator);
			    }//End of signal sample
			  }
			}
		      }
		      else{
			cout << " Looking at wrong bTagging combination in calculation loops " << endl;
			cout << " Looking at : TCHE = "<<TCHEbTagLoop <<" TCHP = "<<TCHPbTagLoop<<" SSVHE = "<<SSVHEbTagLoop<<" SSVHP = "<<SSVHPbTagLoop<<" CSV = "<<CSVbTagLoop<< " JP = " << JPbTagLoop << " JBP = " << JBPbTagLoop << endl;
			cout << " Correct combination : TCHE = "<< UsedTCHE[SumBTag] << " TCHP = " << UsedTCHP[SumBTag] << " SSVHE = " << UsedSSVHE[SumBTag] << " SSVHP = " << UsedSSVHP[SumBTag] << " CSV = " << UsedCSV[SumBTag] << " JP = " << UsedJP[SumBTag] << " JBP = " << UsedJBP[SumBTag] << endl;
		      }
		  
		      if(TCHEbTagLoop == NumberTCHEbTags){
			ConsideredBTagger = 1;
			TCHPbTagLoop = 2;
			SSVHEbTagLoop = 2;
			SSVHPbTagLoop = 2;
			CSVbTagLoop = 2;
			JPbTagLoop = 2;
			JBPbTagLoop = 2;
		      }
		    }//end of TCHE		
		    if(TCHPbTagLoop == NumberTCHPbTags) ConsideredBTagger = 2;
		  }//end of TCHP      
		  if(SSVHEbTagLoop == NumberSSVHEbTags) ConsideredBTagger = 3;
		}//end of SSVHE
		if(SSVHPbTagLoop == NumberSSVHPbTags) ConsideredBTagger = 4;
	      }//end of SSVHP
	      if(CSVbTagLoop == NumberCSVbTags) ConsideredBTagger = 5;
	    }//end of CSV
	    if(JPbTagLoop == NumberJPbTags) ConsideredBTagger = 6;
	  }//end of JP
	}//end of JBP	    
	      		
      }// End of loop over selected events      
    } // end loop over events in wTrees    
    
    std::cout << "  " << endl;
    std::cout << " size of cos theta : " << endl;
    std::cout << "            - Hadronic KinFit only                : " << CosThetaValuesKinFitHadr[0].size() << endl;
    std::cout << "            - Hadronic KinFit and leptonic W only : " << CosThetaValuesKinFitHadrAndLeptWOnly[0].size() << endl;
    std::cout << "            - Hadronic KinFit and leptonic W      : " << CosThetaValuesKinFitHadrAndLept[0].size() << endl;
    std::cout << " °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°° " << endl;
    std::cout << "  " << endl;
    std::cout << " size of generator cos theta :   -->Should not be doubled when looking at JES semiMu of semiEl" << endl;
    std::cout << "            - Hadronic KinFit only                : " << CosThGenKinFitHadr[0].size() << endl;
    std::cout << "            - Hadronic KinFit and leptonic W only : " << CosThGenKinFitHadrAndLeptWOnly[0].size() << endl;
    std::cout << "            - Hadronic KinFit and leptonic W      : " << CosThGenKinFitHadrAndLept[0].size() << endl;
    std::cout << " °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°° " << endl;
   
    if(iDataSet==(datasets.size()-1)){
      FMinusTexHadr << " bTagger  &  FL result & Diff from SM & Stat uncert & JES Syst & JER Syst & W Syst & TTScale Syst & TTMatch Syst & PU Syst & UnclusEnergy Syst & Top Mass Syst & TriggEvtSel Syst  & Total Syst uncert & Total stat + syst uncert " << endl;
      FPlusTexHadr << " bTagger  &  FR result & Diff from SM & Stat uncert & JES Syst & JER Syst & W Syst & TTScale Syst & TTMatch Syst & PU Syst & UnclusEnergy Syst & Top Mass Syst & TriggEvtSel Syst & Total Syst uncert & Total stat + syst uncert " << endl;
      FZeroTexHadr << " bTagger  &  F0 result & Diff from SM & Stat uncert & JES Syst & JER Syst & W Syst & TTScale Syst & TTMatch Syst & PU Syst & UnclusEnergy Syst & Top Mass Syst & TriggEvtSel Syst & Total Syst uncert & Total stat + syst uncert " << endl;

      FMinusTexHadrAndLeptWOnly<<" bTagger  &  FL result & Diff from SM & Stat uncert & JES Syst & JER Syst & W Syst & TTScale Syst & TTMatch Syst & PU Syst & UnclusEnergy Syst & Top Mass Syst & TriggEvtSel Syst & Total Syst uncert & Total stat + syst uncert "<<endl;
      FPlusTexHadrAndLeptWOnly<<" bTagger  &  FR result & Diff from SM & Stat uncert & JES Syst & JER Syst & W Syst & TTScale Syst & TTMatch Syst & PU Syst & UnclusEnergy Syst & Top Mass Syst & TriggEvtSel Syst & Total Syst uncert & Total stat + syst uncert "<< endl;
      FZeroTexHadrAndLeptWOnly<<" bTagger  &  F0 result & Diff from SM & Stat uncert & JES Syst & JER Syst & W Syst & TTScale Syst & TTMatch Syst & PU Syst & UnclusEnergy Syst & Top Mass Syst & TriggEvtSel Syst & Total Syst uncert & Total stat + syst uncert "<< endl;

      FMinusTexHadrAndLept << " bTagger  &  FL result & Diff from SM & Stat uncert & JES Syst & JER Syst & W Syst & TTScale Syst & TTMatch Syst & PU Syst & UnclusEnergy Syst & Top Mass Syst & TriggEvtSel Syst & Total Syst uncert & Total stat + syst uncert " << endl;
      FPlusTexHadrAndLept << " bTagger  &  FR result & Diff from SM & Stat uncert & JES Syst & JER Syst & W Syst & TTScale Syst & TTMatch Syst & PU Syst & UnclusEnergy Syst & Top Mass Syst & TriggEvtSel Syst & Total Syst uncert & Total stat + syst uncert " << endl;
      FZeroTexHadrAndLept << " bTagger  &  F0 result & Diff from SM & Stat uncert & JES Syst & JER Syst & W Syst & TTScale Syst & TTMatch Syst & PU Syst & UnclusEnergy Syst & Top Mass Syst & TriggEvtSel Syst & Total Syst uncert & Total stat + syst uncert " << endl;
    }
    
    ////////////////////////////////
    //    Execute MinuitFitter:   //
    ////////////////////////////////
    if(PerformMinuit == true){

      int SumbTag, ConsideredTagger, JBP, JP, CSV, SSVHP, SSVHE, TCHP, TCHE;
      ConsideredTagger = 0;
      for(JBP=0; JBP<NumberJBPbTags;JBP++){
	cout << " JBP value : " << JBP << endl;
	for(JP=0;JP<=NumberJPbTags;JP++){
	  if(ConsideredTagger == 6) JP =13;   //--> Do not run over all possible TCHE bTag values!	    
	  for(CSV =0;CSV<= NumberCSVbTags; CSV++){
	    if(ConsideredTagger == 5 || ConsideredTagger ==6) CSV = 13;
	    for(SSVHP=0;SSVHP <=NumberSSVHPbTags; SSVHP++){
	      if(ConsideredTagger == 4 || ConsideredTagger == 5 || ConsideredTagger ==6) SSVHP =13;   	    
	      for(SSVHE=0;SSVHE<=NumberSSVHEbTags; SSVHE++){
		if(ConsideredTagger == 3 || ConsideredTagger == 4 || ConsideredTagger == 5 || ConsideredTagger ==6) SSVHE = 13;
		for(TCHP=0; TCHP <= NumberTCHPbTags; TCHP++){
		  if(ConsideredTagger == 2 || ConsideredTagger == 3 || ConsideredTagger == 4 || ConsideredTagger == 5 || ConsideredTagger ==6) TCHP =13; 
		  for(TCHE=0; TCHE<= NumberTCHEbTags;TCHE++){
		    if(ConsideredTagger ==1 || ConsideredTagger ==2 || ConsideredTagger ==3 || ConsideredTagger ==4 || ConsideredTagger == 5 || ConsideredTagger ==6) TCHE =13;
		    SumbTag = TCHE+TCHP+SSVHE+SSVHP+CSV+JP+JBP;
		    if(UsedTCHE[SumbTag]==(TCHE+1) && UsedTCHP[SumbTag]==(TCHP+1) && UsedSSVHE[SumbTag]==(SSVHE+1) && UsedSSVHP[SumbTag]==(SSVHP+1) && UsedCSV[SumbTag]==(CSV+1) && UsedJP[SumbTag] == (JP+1) && UsedJBP[SumbTag] == (JBP+1) ){

		      if((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ){//Defining of the genttbar histo
			char hisnameHadr[100];
			char hisnameHadrAndLeptWOnly[100];
			char hisnameHadrAndLept[100];
			sprintf(hisnameHadr,"CosThetaGenHadr_TCHE%s_TCHP%s_SSVHE%s_SSVHP%s_CSV%s_JP%s_JBP%s", bTag[TCHE].c_str(),bTag[TCHP].c_str(),bTag[SSVHE].c_str(),bTag[SSVHP].c_str(),bTag[CSV].c_str(),bTag[JP].c_str(),bTag[JBP].c_str());
			sprintf(hisnameHadrAndLeptWOnly,"CosThetaGenHadrAndLeptWOnly_TCHE%s_TCHP%s_SSVHE%s_SSVHP%s_CSV%s_JP%s_JBP%s", bTag[TCHE].c_str(),bTag[TCHP].c_str(),bTag[SSVHE].c_str(),bTag[SSVHP].c_str(),bTag[CSV].c_str(),bTag[JP].c_str(),bTag[JBP].c_str());
			sprintf(hisnameHadrAndLept,"CosThetaGenHadrAndLept_TCHE%s_TCHP%s_SSVHE%s_SSVHP%s_CSV%s_JP%s_JBP%s", bTag[TCHE].c_str(),bTag[TCHP].c_str(),bTag[SSVHE].c_str(),bTag[SSVHP].c_str(),bTag[CSV].c_str(),bTag[JP].c_str(),bTag[JBP].c_str());
			for (int ibinn=0; ibinn<CosThetaBinNumber; ibinn++){
			  genttbarhistoHadr[ibinn]= new TNtuple(hisnameHadr,hisnameHadr,"costhgen:evtweight");
			  genttbarhistoHadrAndLeptWOnly[ibinn]= new TNtuple(hisnameHadrAndLeptWOnly,hisnameHadrAndLeptWOnly,"costhgen:evtweight");
			  genttbarhistoHadrAndLept[ibinn]= new TNtuple(hisnameHadrAndLept,hisnameHadrAndLept,"costhgen:evtweight");
			  genttbarhistoHadr[ibinn]->SetDirectory(0);
			  genttbarhistoHadrAndLeptWOnly[ibinn]->SetDirectory(0);
			  genttbarhistoHadrAndLept[ibinn]->SetDirectory(0);
			}
		      }
		      float scaleFactor = 1.;  //Make an array of the scaleFactor such that it can be different for every event !!
		    
		      if(HadronicOnly == true){
			//Data:
			std::string CosThetaDataStringHadr = "CosThetaDataHadr_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//JES Plus:
			std::string CosThetaJESPlusStringHadr = "CosThetaJESPlusHadr_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//JES Minus:
			std::string CosThetaJESMinusStringHadr = "CosThetaJESMinusHadr_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//JER Plus:
			std::string CosThetaJERPlusStringHadr = "CosThetaJERPlusHadr_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//JER Minus:
			std::string CosThetaJERMinusStringHadr = "CosThetaJERMinusHadr_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//WSyst Plus:
			std::string CosThetaWPlusStringHadr = "CosThetaWPlusHadr_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//WSyst Minus:
			std::string CosThetaWMinusStringHadr = "CosThetaWMinusHadr_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//TTScaling Plus:
			std::string CosThetaTTScalingUpStringHadr = "CosThetaTTScalingUpHadr_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//TTScaling Minus:
			std::string CosThetaTTScalingDownStringHadr = "CosThetaTTScalingDownHadr_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//TTMatching Plus:
			std::string CosThetaTTMatchingUpStringHadr = "CosThetaTTMatchingUpHadr_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//TTMatching Minus:
			std::string CosThetaTTMatchingDownStringHadr = "CosThetaTTMatchingDownHadr_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//PU Plus:
			std::string CosThetaPUPlusStringHadr = "CosThetaPUPlusHadr_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//PU Minus:
			std::string CosThetaPUMinusStringHadr = "CosThetaPUMinusHadr_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//UnclusEnergy Plus:
			std::string CosThetaUnclusEnergyPlusStringHadr = "CosThetaUnclusEnergyPlusHadr_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//UnclusEnergy Minus:
			std::string CosThetaUnclusEnergyMinusStringHadr = "CosThetaUnclusEnergyMinusHadr_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//TopMass Plus:
			std::string CosThetaTopMassPlusStringHadr = "CosThetaTopMassPlusHadr_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//TopMass Minus:
			std::string CosThetaTopMassMinusStringHadr = "CosThetaTopMassMinusHadr_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//TriggEvtSel Plus:
			std::string CosThetaTriggEvtSelPlusStringHadr = "CosThetaTriggEvtSelPlusHadr_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//TriggEvtSel Minus:
			std::string CosThetaTriggEvtSelMinusStringHadr = "CosThetaTriggEvtSelMinusHadr_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			std::string CosThetaSignalStringHadr = "CosThetaSignalHadr_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			std::string CosThetaBckgStringHadr = "CosThetaBckgHadr_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];

			if(iDataSet == 0){
			  //Data:
			  histo1D[CosThetaDataStringHadr]=new TH1F(CosThetaDataStringHadr.c_str(),CosThetaDataStringHadr.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaDataStringHadr]->SetDirectory(0);
			  //JESPlus:
			  histo1D[CosThetaJESPlusStringHadr]=new TH1F(CosThetaJESPlusStringHadr.c_str(),CosThetaJESPlusStringHadr.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaJESPlusStringHadr]->SetDirectory(0);
			  //JESMinus:
			  histo1D[CosThetaJESMinusStringHadr]=new TH1F(CosThetaJESMinusStringHadr.c_str(),CosThetaJESMinusStringHadr.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaJESMinusStringHadr]->SetDirectory(0);
			  //JERPlus:
			  histo1D[CosThetaJERPlusStringHadr]=new TH1F(CosThetaJERPlusStringHadr.c_str(),CosThetaJERPlusStringHadr.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaJERPlusStringHadr]->SetDirectory(0);
			  //JERMinus:
			  histo1D[CosThetaJERMinusStringHadr]=new TH1F(CosThetaJERMinusStringHadr.c_str(),CosThetaJERMinusStringHadr.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaJERMinusStringHadr]->SetDirectory(0);
			  //WJets Plus:
			  histo1D[CosThetaWPlusStringHadr]=new TH1F(CosThetaWPlusStringHadr.c_str(),CosThetaWPlusStringHadr.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaWPlusStringHadr]->SetDirectory(0);
			  //WJets Minus:
			  histo1D[CosThetaWMinusStringHadr]=new TH1F(CosThetaWMinusStringHadr.c_str(),CosThetaWMinusStringHadr.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaWMinusStringHadr]->SetDirectory(0);
			  //TTScaling Up:
			  histo1D[CosThetaTTScalingUpStringHadr]=new TH1F(CosThetaTTScalingUpStringHadr.c_str(),CosThetaTTScalingUpStringHadr.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaTTScalingUpStringHadr]->SetDirectory(0);
			  //TTScaling Down:
			  histo1D[CosThetaTTScalingDownStringHadr]=new TH1F(CosThetaTTScalingDownStringHadr.c_str(),CosThetaTTScalingDownStringHadr.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaTTScalingDownStringHadr]->SetDirectory(0);
			  //TTMatching Up:
			  histo1D[CosThetaTTMatchingUpStringHadr]=new TH1F(CosThetaTTMatchingUpStringHadr.c_str(),CosThetaTTMatchingUpStringHadr.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaTTMatchingUpStringHadr]->SetDirectory(0);
			  //TTMatching Down:
			  histo1D[CosThetaTTMatchingDownStringHadr]=new TH1F(CosThetaTTMatchingDownStringHadr.c_str(),CosThetaTTMatchingDownStringHadr.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaTTMatchingDownStringHadr]->SetDirectory(0);
			  //PU Plus:
			  histo1D[CosThetaPUPlusStringHadr]=new TH1F(CosThetaPUPlusStringHadr.c_str(),CosThetaPUPlusStringHadr.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaPUPlusStringHadr]->SetDirectory(0);
			  //PU Minus:
			  histo1D[CosThetaPUMinusStringHadr]=new TH1F(CosThetaPUMinusStringHadr.c_str(),CosThetaPUMinusStringHadr.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaPUMinusStringHadr]->SetDirectory(0);
			  //UnclusEnergy Plus:
			  histo1D[CosThetaUnclusEnergyPlusStringHadr]=new TH1F(CosThetaUnclusEnergyPlusStringHadr.c_str(),CosThetaUnclusEnergyPlusStringHadr.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaUnclusEnergyPlusStringHadr]->SetDirectory(0);
			  //UnclusEnergy Minus:
			  histo1D[CosThetaUnclusEnergyMinusStringHadr]=new TH1F(CosThetaUnclusEnergyMinusStringHadr.c_str(),CosThetaUnclusEnergyMinusStringHadr.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaUnclusEnergyMinusStringHadr]->SetDirectory(0);
			  //TopMass Plus:
			  histo1D[CosThetaTopMassPlusStringHadr]=new TH1F(CosThetaTopMassPlusStringHadr.c_str(),CosThetaTopMassPlusStringHadr.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaTopMassPlusStringHadr]->SetDirectory(0);
			  //TopMass Minus:
			  histo1D[CosThetaTopMassMinusStringHadr]=new TH1F(CosThetaTopMassMinusStringHadr.c_str(),CosThetaTopMassMinusStringHadr.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaTopMassMinusStringHadr]->SetDirectory(0);
			  //TriggEvtSel Plus:
			  histo1D[CosThetaTriggEvtSelPlusStringHadr]=new TH1F(CosThetaTriggEvtSelPlusStringHadr.c_str(),CosThetaTriggEvtSelPlusStringHadr.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaTriggEvtSelPlusStringHadr]->SetDirectory(0);
			  //TriggEvtSel Minus:
			  histo1D[CosThetaTriggEvtSelMinusStringHadr]=new TH1F(CosThetaTriggEvtSelMinusStringHadr.c_str(),CosThetaTriggEvtSelMinusStringHadr.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaTriggEvtSelMinusStringHadr]->SetDirectory(0);
			  histo1D[CosThetaSignalStringHadr]=new TH1F(CosThetaSignalStringHadr.c_str(),CosThetaSignalStringHadr.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaSignalStringHadr]->SetDirectory(0);
			  histo1D[CosThetaBckgStringHadr]=new TH1F(CosThetaBckgStringHadr.c_str(),CosThetaBckgStringHadr.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaBckgStringHadr]->SetDirectory(0);
			}
		    
			//Filling of the histograms:		    
			for(int ii=0; ii<CosThetaValuesKinFitHadr[SumbTag].size();ii++){		
		      
			  if((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ){//Defining of genttbar histo
			    for(int iBin=0; iBin< CosThetaBinNumber; iBin++){//Filling of genttbar histo:		   
			      if(CosThetaValuesKinFitHadr[SumbTag][ii] >= binEdge[iBin] && CosThetaValuesKinFitHadr[SumbTag][ii] < binEdge[iBin+1]){
				genttbarhistoHadr[iBin]->Fill(CosThGenKinFitHadr[SumbTag][ii], EventCorrectionWeightKinFitHadr[SumbTag][ii]) ; 
			      }
			      else if(CosThetaValuesKinFitHadr[SumbTag][ii] ==1){ //1 is included in last bin
				genttbarhistoHadr[CosThetaBinNumber-1]->Fill(CosThGenKinFitHadr[SumbTag][ii], EventCorrectionWeightKinFitHadr[SumbTag][ii]);
			      }	      	   
			    }
			  }
		      
			  ///////////////////////////////////////////////////////////////////////////////
			  // Change data to systematics since then nominal values will be reweighted!! //
			  ///////////////////////////////////////////////////////////////////////////////
		      
			  //Data result:
			  if(dataSetName.find("Data") == 0)
			    histo1D[CosThetaDataStringHadr]->Fill(CosThetaValuesKinFitHadr[SumbTag][ii],EventCorrectionWeightKinFitHadr[SumbTag][ii]);
			  //JES Plus:
			  if(dataSetName.find("JESPlus") == 0)
			    histo1D[CosThetaJESPlusStringHadr]->Fill(CosThetaValuesKinFitHadr[SumbTag][ii],EventCorrectionWeightKinFitHadr[SumbTag][ii]);
			  //JES Minus:
			  if(dataSetName.find("JESMinus") == 0)
			    histo1D[CosThetaJESMinusStringHadr]->Fill(CosThetaValuesKinFitHadr[SumbTag][ii],EventCorrectionWeightKinFitHadr[SumbTag][ii]);
			  //JER Plus:
			  if(dataSetName.find("JERPlus") == 0)
			    histo1D[CosThetaJERPlusStringHadr]->Fill(CosThetaValuesKinFitHadr[SumbTag][ii],EventCorrectionWeightKinFitHadr[SumbTag][ii]);
			  //JER Minus:
			  if(dataSetName.find("JERMinus") == 0)
			    histo1D[CosThetaJERMinusStringHadr]->Fill(CosThetaValuesKinFitHadr[SumbTag][ii],EventCorrectionWeightKinFitHadr[SumbTag][ii]);
			  //WJets Plus:
			  if( dataSetName.find("Data") !=0 && dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSystMinus") != 0  && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatching") != 0 && dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMass") != 0 && dataSetName.find("TriggEvtSel") != 0)
			    histo1D[CosThetaWPlusStringHadr]->Fill(CosThetaValuesKinFitHadr[SumbTag][ii],EventCorrectionWeightKinFitHadr[SumbTag][ii]);
			  //WJets Minus:
			  if( dataSetName.find("Data") !=0 && dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSystPlus") != 0  && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatching") != 0 && dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMass") != 0 && dataSetName.find("TriggEvtSel") != 0)
			    histo1D[CosThetaWMinusStringHadr]->Fill(CosThetaValuesKinFitHadr[SumbTag][ii],EventCorrectionWeightKinFitHadr[SumbTag][ii]);
			  //TTScaling Up:
			  if( dataSetName.find("Data") !=0 && dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSyst") != 0 && dataSetName.find("TTScalingDown") != 0&& dataSetName.find("TTMatching") != 0 && dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMass") != 0 && dataSetName.find("TriggEvtSel") != 0 && dataSetName.find("TTbarJets") !=0)
			    histo1D[CosThetaTTScalingUpStringHadr]->Fill(CosThetaValuesKinFitHadr[SumbTag][ii],EventCorrectionWeightKinFitHadr[SumbTag][ii]);
			  //TTScaling Down:
			  if( dataSetName.find("Data") !=0 && dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSyst") != 0 && dataSetName.find("TTScalingUp") != 0 && dataSetName.find("TTMatching") != 0 && dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMass") != 0 && dataSetName.find("TriggEvtSel") != 0 && dataSetName.find("TTbarJets") !=0 )
			    histo1D[CosThetaTTScalingDownStringHadr]->Fill(CosThetaValuesKinFitHadr[SumbTag][ii],EventCorrectionWeightKinFitHadr[SumbTag][ii]);
			  //TTMatching Up:
			  if(dataSetName.find("Data") !=0&& dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSyst") != 0 && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatchingDown") != 0 && dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMass") != 0 && dataSetName.find("TriggEvtSel") != 0 && dataSetName.find("TTbarJets") !=0 )
			    histo1D[CosThetaTTMatchingUpStringHadr]->Fill(CosThetaValuesKinFitHadr[SumbTag][ii],EventCorrectionWeightKinFitHadr[SumbTag][ii]);
			  //TTMatching Down:
			  if( dataSetName.find("Data") !=0 && dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSyst") != 0 && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatchingUp") != 0 && dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMass") != 0 && dataSetName.find("TriggEvtSel") != 0 && dataSetName.find("TTbarJets") !=0 )
			    histo1D[CosThetaTTMatchingDownStringHadr]->Fill(CosThetaValuesKinFitHadr[SumbTag][ii],EventCorrectionWeightKinFitHadr[SumbTag][ii]);	
			  //PU Plus:
			  if(dataSetName.find("Data") !=0&& dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSyst") != 0 && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatching") != 0 && dataSetName.find("PUMinus") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMass") != 0 && dataSetName.find("TriggEvtSel") != 0 && dataSetName.find("TTbarJets") !=0 )
			    histo1D[CosThetaPUPlusStringHadr]->Fill(CosThetaValuesKinFitHadr[SumbTag][ii],EventCorrectionWeightKinFitHadr[SumbTag][ii]);
			  //PU Minus:
			  if(dataSetName.find("Data") !=0&& dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSyst") != 0 && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatching") != 0 && dataSetName.find("PUPlus") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMass") != 0 && dataSetName.find("TriggEvtSel") != 0 && dataSetName.find("TTbarJets") !=0 )
			    histo1D[CosThetaPUMinusStringHadr]->Fill(CosThetaValuesKinFitHadr[SumbTag][ii],EventCorrectionWeightKinFitHadr[SumbTag][ii]);
			  //UnclusEnergy Plus:
			  if(dataSetName.find("Data") !=0&& dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSyst") != 0 && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatching") != 0 && dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergyMinus") != 0 && dataSetName.find("TopMass") != 0 && dataSetName.find("TriggEvtSel") != 0 && dataSetName.find("TTbarJets") !=0 )
			    histo1D[CosThetaUnclusEnergyPlusStringHadr]->Fill(CosThetaValuesKinFitHadr[SumbTag][ii],EventCorrectionWeightKinFitHadr[SumbTag][ii]);
			  //UnclusEnergy Minus:
			  if(dataSetName.find("Data") !=0&& dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSyst") != 0 && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatching") != 0 && dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergyPlus") != 0 && dataSetName.find("TopMass") != 0 && dataSetName.find("TriggEvtSel") != 0 && dataSetName.find("TTbarJets") !=0 )
			    histo1D[CosThetaUnclusEnergyMinusStringHadr]->Fill(CosThetaValuesKinFitHadr[SumbTag][ii],EventCorrectionWeightKinFitHadr[SumbTag][ii]);
			  //Top Mass Plus:
			  if(dataSetName.find("Data") !=0&& dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSyst") != 0 && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatching") != 0 && dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMassMinus") != 0 && dataSetName.find("TriggEvtSel") != 0 && dataSetName.find("TTbarJets") !=0 )
			    histo1D[CosThetaTopMassPlusStringHadr]->Fill(CosThetaValuesKinFitHadr[SumbTag][ii],EventCorrectionWeightKinFitHadr[SumbTag][ii]);
			  //Top Mass Minus:
			  if(dataSetName.find("Data") !=0&& dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSyst") != 0 && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatching") != 0 && dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMassPlus") != 0 && dataSetName.find("TriggEvtSel") != 0 && dataSetName.find("TTbarJets") !=0 )		      
			    histo1D[CosThetaTopMassMinusStringHadr]->Fill(CosThetaValuesKinFitHadr[SumbTag][ii],EventCorrectionWeightKinFitHadr[SumbTag][ii]);
			  //TriggEvtSel Plus:
			  if(dataSetName.find("Data") !=0&& dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSyst") != 0 && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatching") != 0 && dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMass") != 0 && dataSetName.find("TriggEvtSelMinus") != 0 && dataSetName.find("TTbarJets") !=0 )
			    histo1D[CosThetaTriggEvtSelPlusStringHadr]->Fill(CosThetaValuesKinFitHadr[SumbTag][ii],EventCorrectionWeightKinFitHadr[SumbTag][ii]);
			  //TrigEvtSel Minus:
			  if(dataSetName.find("Data") !=0&& dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSyst") != 0 && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatching") != 0 && dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMass") != 0 && dataSetName.find("TriggEvtSelPlus") != 0 && dataSetName.find("TTbarJets") !=0 )		      
			    histo1D[CosThetaTriggEvtSelMinusStringHadr]->Fill(CosThetaValuesKinFitHadr[SumbTag][ii],EventCorrectionWeightKinFitHadr[SumbTag][ii]);
	      
			  if(SignalOnly == true && ((dataSetName.find("TTbarJets_SemiMu") == 0 && semiMuon==true) || (dataSetName.find("TTbarJets_SemiEl")==0 && semiElectron == true)))	    
			    histo1D[CosThetaDataStringHadr]->Fill(CosThetaValuesKinFitHadr[SumbTag][ii],EventCorrectionWeightKinFitHadr[SumbTag][ii]);
		      
			  if(((dataSetName.find("TTbarJets_SemiMu")==0 && semiMuon==true) || (dataSetName.find("TTbarJets_SemiEl")==0 && semiElectron == true) ) && dataSetName.find("JES_") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSyst") != 0 && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatching") != 0 && dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMass") != 0 && dataSetName.find("TriggEvtSel") != 0)
			    histo1D[CosThetaSignalStringHadr]->Fill(CosThetaValuesKinFitHadr[SumbTag][ii],EventCorrectionWeightKinFitHadr[SumbTag][ii]);//Luminosity*scaleFactor*(LumiWeightVectorKinFitHadr[SumbTag][ii])*NominalNormFactor);
		      
			  if(dataSetName.find("Data") !=0 && dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSyst") != 0 && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatching") != 0 && dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMass") != 0 && dataSetName.find("TriggEvtSel") != 0 && ((dataSetName.find("TTbarJets_SemiMu")!=0 && semiMuon==true) || (dataSetName.find("TTbarJets_SemiEl")!=0 && semiElectron == true) )) 
			    histo1D[CosThetaBckgStringHadr]->Fill(CosThetaValuesKinFitHadr[SumbTag][ii],EventCorrectionWeightKinFitHadr[SumbTag][ii]);//Luminosity*scaleFactor*(LumiWeightVectorKinFitHadr[SumbTag][ii])*NominalNormFactor);
			}

			//Perform Minuit Fitter:
			if(iDataSet==(datasets.size()-1)){//Go in this loop when the last datasample is active to perform MinuitFitter

			  //Data:
			  MinuitFitter minuitFitterHadr = MinuitFitter(histo1D[CosThetaDataStringHadr], histo1D[CosThetaSignalStringHadr], histo1D[CosThetaBckgStringHadr], 0.6671, 0.3325, 0.0004,genttbarhistoHadr,ndimen);
			  //JESPlus
			  MinuitFitter minuitFitterHadrJESPlus = MinuitFitter(histo1D[CosThetaJESPlusStringHadr], histo1D[CosThetaSignalStringHadr], histo1D[CosThetaBckgStringHadr], 0.6671, 0.3325, 0.0004,genttbarhistoHadr,ndimen);			  
			  //JES Minus
			  MinuitFitter minuitFitterHadrJESMinus = MinuitFitter(histo1D[CosThetaJESMinusStringHadr], histo1D[CosThetaSignalStringHadr], histo1D[CosThetaBckgStringHadr], 0.6671, 0.3325, 0.0004,genttbarhistoHadr,ndimen);			  
			  //JERPlus
			  MinuitFitter minuitFitterHadrJERPlus = MinuitFitter(histo1D[CosThetaJERPlusStringHadr], histo1D[CosThetaSignalStringHadr], histo1D[CosThetaBckgStringHadr], 0.6671, 0.3325, 0.0004,genttbarhistoHadr,ndimen);			  
			  //JER Minus
			  MinuitFitter minuitFitterHadrJERMinus = MinuitFitter(histo1D[CosThetaJERMinusStringHadr], histo1D[CosThetaSignalStringHadr], histo1D[CosThetaBckgStringHadr], 0.6671, 0.3325, 0.0004,genttbarhistoHadr,ndimen);			  
			  //W Plus
			  MinuitFitter minuitFitterHadrWPlus = MinuitFitter(histo1D[CosThetaWPlusStringHadr], histo1D[CosThetaSignalStringHadr], histo1D[CosThetaBckgStringHadr], 0.6671, 0.3325, 0.0004,genttbarhistoHadr,ndimen);			  
			  //W Minus
			  MinuitFitter minuitFitterHadrWMinus = MinuitFitter(histo1D[CosThetaWMinusStringHadr], histo1D[CosThetaSignalStringHadr], histo1D[CosThetaBckgStringHadr], 0.6671, 0.3325, 0.0004,genttbarhistoHadr,ndimen);
			  //TTScaling Up
			  MinuitFitter minuitFitterHadrTTScalingUp = MinuitFitter(histo1D[CosThetaTTScalingUpStringHadr], histo1D[CosThetaSignalStringHadr], histo1D[CosThetaBckgStringHadr], 0.6671, 0.3325, 0.0004,genttbarhistoHadr,ndimen);
			  //TTScaling Down
			  MinuitFitter minuitFitterHadrTTScalingDown = MinuitFitter(histo1D[CosThetaTTScalingDownStringHadr], histo1D[CosThetaSignalStringHadr], histo1D[CosThetaBckgStringHadr], 0.6671, 0.3325, 0.0004,genttbarhistoHadr,ndimen);
			  //TTMatching Up
			  MinuitFitter minuitFitterHadrTTMatchingUp = MinuitFitter(histo1D[CosThetaTTMatchingUpStringHadr], histo1D[CosThetaSignalStringHadr], histo1D[CosThetaBckgStringHadr], 0.6671, 0.3325, 0.0004,genttbarhistoHadr,ndimen);
			  //TTMatching Down
			  MinuitFitter minuitFitterHadrTTMatchingDown = MinuitFitter(histo1D[CosThetaTTMatchingDownStringHadr], histo1D[CosThetaSignalStringHadr], histo1D[CosThetaBckgStringHadr], 0.6671, 0.3325, 0.0004,genttbarhistoHadr,ndimen);
			  //PU Plus
			  MinuitFitter minuitFitterHadrPUPlus = MinuitFitter(histo1D[CosThetaPUPlusStringHadr], histo1D[CosThetaSignalStringHadr], histo1D[CosThetaBckgStringHadr], 0.6671, 0.3325, 0.0004,genttbarhistoHadr,ndimen);			  
			  //PU Minus
			  MinuitFitter minuitFitterHadrPUMinus = MinuitFitter(histo1D[CosThetaPUMinusStringHadr], histo1D[CosThetaSignalStringHadr], histo1D[CosThetaBckgStringHadr], 0.6671, 0.3325, 0.0004,genttbarhistoHadr,ndimen);			  
			  //UnclusEnergy Plus
			  MinuitFitter minuitFitterHadrUnclusEnergyPlus = MinuitFitter(histo1D[CosThetaUnclusEnergyPlusStringHadr], histo1D[CosThetaSignalStringHadr], histo1D[CosThetaBckgStringHadr], 0.6671, 0.3325, 0.0004,genttbarhistoHadr,ndimen);  
			  //UnclusEnergy Minus
			  MinuitFitter minuitFitterHadrUnclusEnergyMinus = MinuitFitter(histo1D[CosThetaUnclusEnergyMinusStringHadr], histo1D[CosThetaSignalStringHadr], histo1D[CosThetaBckgStringHadr], 0.6671, 0.3325, 0.0004,genttbarhistoHadr,ndimen);
			  //TopMass Plus
			  MinuitFitter minuitFitterHadrTopMassPlus = MinuitFitter(histo1D[CosThetaTopMassPlusStringHadr], histo1D[CosThetaSignalStringHadr], histo1D[CosThetaBckgStringHadr], 0.6671, 0.3325, 0.0004,genttbarhistoHadr,ndimen);			  
			  //TopMass Minus
			  MinuitFitter minuitFitterHadrTopMassMinus = MinuitFitter(histo1D[CosThetaTopMassMinusStringHadr], histo1D[CosThetaSignalStringHadr], histo1D[CosThetaBckgStringHadr], 0.6671, 0.3325, 0.0004,genttbarhistoHadr,ndimen);			  
			  //TriggEvtSel Plus
			  MinuitFitter minuitFitterHadrTriggEvtSelPlus = MinuitFitter(histo1D[CosThetaTriggEvtSelPlusStringHadr], histo1D[CosThetaSignalStringHadr], histo1D[CosThetaBckgStringHadr], 0.6671, 0.3325, 0.0004,genttbarhistoHadr,ndimen);			  
			  //TriggEvtSel Minus
			  MinuitFitter minuitFitterHadrTriggEvtSelMinus = MinuitFitter(histo1D[CosThetaTriggEvtSelMinusStringHadr], histo1D[CosThetaSignalStringHadr], histo1D[CosThetaBckgStringHadr], 0.6671, 0.3325, 0.0004,genttbarhistoHadr,ndimen);			  
			  
			  PresentationTexHadr << PresentationOutput[SumbTag] << " & ";
			  cout << " Inside minuitFitter loop : " << PresentationOutput[SumbTag] << " & " << endl;
			  for(int ii=0; ii< nameDataSet.size(); ii++){
			    if(ii < nameDataSet.size()-1){ PresentationTexHadr << NumberRemainingEventsKinFitHadr[SumbTag][ii] << " & ";}
			    else if(ii == nameDataSet.size()-1 ){ 
			      PresentationTexHadr << NumberRemainingEventsKinFitHadr[SumbTag][ii] << " & ";//For presentation helicity values still need to be included !!
			      PresentationTexHadr << NumberBLeptCorrectEventsKinFitHadr[SumbTag][ii] << " & ";
			    }
			  }	
		  
			  //FL results:
			  float FLJESSystHadr = (abs(minuitFitterHadr.GetFLResult() - minuitFitterHadrJESPlus.GetFLResult()) + abs(minuitFitterHadr.GetFLResult() - minuitFitterHadrJESMinus.GetFLResult()))/2;
			  float FLJERSystHadr = 0;
			  if(JERResults == true) FLJERSystHadr = (abs(minuitFitterHadr.GetFLResult() - minuitFitterHadrJERPlus.GetFLResult()) + abs(minuitFitterHadr.GetFLResult() - minuitFitterHadrJERMinus.GetFLResult()))/2;
			  float FLWSystHadr = 0;
			  if(WSystResults == true) FLWSystHadr = (abs(minuitFitterHadr.GetFLResult() - minuitFitterHadrWPlus.GetFLResult()) + abs(minuitFitterHadr.GetFLResult() - minuitFitterHadrWMinus.GetFLResult()))/2;
			  float FLTTScalingHadr = 0;
			  if(TTScalingResults == true) FLTTScalingHadr = (abs(minuitFitterHadr.GetFLResult() - minuitFitterHadrTTScalingUp.GetFLResult()) + abs(minuitFitterHadr.GetFLResult() - minuitFitterHadrTTScalingDown.GetFLResult()))/2;
			  float FLTTMatchingHadr = 0;
			  if(TTMatchingResults == true) FLTTMatchingHadr = (abs(minuitFitterHadr.GetFLResult() - minuitFitterHadrTTMatchingUp.GetFLResult()) + abs(minuitFitterHadr.GetFLResult() - minuitFitterHadrTTMatchingDown.GetFLResult()))/2;
			  float FLPUSystHadr = 0;
			  if(PUResults == true) FLPUSystHadr = (abs(minuitFitterHadr.GetFLResult() - minuitFitterHadrPUPlus.GetFLResult()) + abs(minuitFitterHadr.GetFLResult() - minuitFitterHadrPUMinus.GetFLResult()))/2;
			  float FLUnclusEnergySystHadr = 0;
			  if(UnclusEnergyResults == true) FLUnclusEnergySystHadr = (abs(minuitFitterHadr.GetFLResult() - minuitFitterHadrUnclusEnergyPlus.GetFLResult()) + abs(minuitFitterHadr.GetFLResult() - minuitFitterHadrUnclusEnergyMinus.GetFLResult()))/2;
			  float FLTopMassSystHadr = 0;
			  if(TopMassResults == true) FLTopMassSystHadr = (abs(minuitFitterHadr.GetFLResult() - minuitFitterHadrTopMassPlus.GetFLResult()) + abs(minuitFitterHadr.GetFLResult() - minuitFitterHadrTopMassMinus.GetFLResult()))/2;		      
			  float FLTriggEvtSelSystHadr = 0;
			  if(TriggEvtSelResults == true) FLTriggEvtSelSystHadr = (abs(minuitFitterHadr.GetFLResult() - minuitFitterHadrTriggEvtSelPlus.GetFLResult()) + abs(minuitFitterHadr.GetFLResult() - minuitFitterHadrTriggEvtSelMinus.GetFLResult()))/2;		      
			  FMinusTexHadr << "\\hline" << endl;
			  FMinusTexHadr << PresentationOutput[SumbTag] << " & " << minuitFitterHadr.GetFLResult() << " & " << minuitFitterHadr.GetFLResult()-SMfrResult <<" & " << minuitFitterHadr.GetFLError() << " & " << FLJESSystHadr << " & " << FLJERSystHadr << " & " << FLWSystHadr << " & " << FLTTScalingHadr << " & " << FLTTMatchingHadr << " & " << FLPUSystHadr << " & " << FLUnclusEnergySystHadr << " & " << FLTopMassSystHadr << " & " << FLTriggEvtSelSystHadr << " & " << sqrt(FLJESSystHadr*FLJESSystHadr + FLJERSystHadr*FLJERSystHadr + FLWSystHadr*FLWSystHadr + FLTTScalingHadr*FLTTScalingHadr + FLTTMatchingHadr*FLTTMatchingHadr + FLPUSystHadr*FLPUSystHadr + FLUnclusEnergySystHadr*FLUnclusEnergySystHadr + FLTopMassSystHadr*FLTopMassSystHadr + FLTriggEvtSelSystHadr*FLTriggEvtSelSystHadr) << " & " << sqrt(FLJESSystHadr*FLJESSystHadr + FLJERSystHadr*FLJERSystHadr + FLWSystHadr*FLWSystHadr + FLTTScalingHadr*FLTTScalingHadr + FLTTMatchingHadr*FLTTMatchingHadr + FLPUSystHadr*FLPUSystHadr + FLUnclusEnergySystHadr*FLUnclusEnergySystHadr + FLTopMassSystHadr*FLTopMassSystHadr + FLTriggEvtSelSystHadr*FLTriggEvtSelSystHadr + (minuitFitterHadr.GetFLError())*(minuitFitterHadr.GetFLError())) << "\\\\" << endl;
		  
			  //FR results:
			  float FRJESSystHadr = (abs(minuitFitterHadr.GetFRResult() - minuitFitterHadrJESPlus.GetFRResult()) + abs(minuitFitterHadr.GetFRResult() - minuitFitterHadrJESMinus.GetFRResult()))/2;
			  float FRJERSystHadr = 0;
			  if(JERResults == true) FRJERSystHadr = (abs(minuitFitterHadr.GetFRResult() - minuitFitterHadrJERPlus.GetFRResult()) + abs(minuitFitterHadr.GetFRResult() - minuitFitterHadrJERMinus.GetFRResult()))/2;
			  float FRWSystHadr = 0;
			  if(WSystResults == true) FRWSystHadr = (abs(minuitFitterHadr.GetFRResult() - minuitFitterHadrWPlus.GetFRResult()) + abs(minuitFitterHadr.GetFRResult() - minuitFitterHadrWMinus.GetFRResult()))/2;
			  float FRTTScalingHadr = 0;
			  if(TTScalingResults == true) FRTTScalingHadr = (abs(minuitFitterHadr.GetFRResult() - minuitFitterHadrTTScalingUp.GetFRResult()) + abs(minuitFitterHadr.GetFRResult() - minuitFitterHadrTTScalingDown.GetFRResult()))/2;
			  float FRTTMatchingHadr = 0;
			  if(TTMatchingResults == true) FRTTMatchingHadr = (abs(minuitFitterHadr.GetFRResult() - minuitFitterHadrTTMatchingUp.GetFRResult()) + abs(minuitFitterHadr.GetFRResult() - minuitFitterHadrTTMatchingDown.GetFRResult()))/2;		  
			  float FRPUSystHadr = 0;
			  if(PUResults == true) FRPUSystHadr = (abs(minuitFitterHadr.GetFRResult() - minuitFitterHadrPUPlus.GetFRResult()) + abs(minuitFitterHadr.GetFRResult() - minuitFitterHadrPUMinus.GetFRResult()))/2;
			  float FRUnclusEnergySystHadr = 0;
			  if(UnclusEnergyResults == true) FRUnclusEnergySystHadr = (abs(minuitFitterHadr.GetFRResult() - minuitFitterHadrUnclusEnergyPlus.GetFRResult()) + abs(minuitFitterHadr.GetFRResult() - minuitFitterHadrUnclusEnergyMinus.GetFRResult()))/2;
			  float FRTopMassSystHadr = 0;
			  if(TopMassResults == true) FRTopMassSystHadr = (abs(minuitFitterHadr.GetFRResult() - minuitFitterHadrTopMassPlus.GetFRResult()) + abs(minuitFitterHadr.GetFRResult() - minuitFitterHadrTopMassMinus.GetFRResult()))/2;
			  float FRTriggEvtSelSystHadr = 0;
			  if(TriggEvtSelResults == true) FRTriggEvtSelSystHadr = (abs(minuitFitterHadr.GetFRResult() - minuitFitterHadrTriggEvtSelPlus.GetFRResult()) + abs(minuitFitterHadr.GetFRResult() - minuitFitterHadrTriggEvtSelMinus.GetFRResult()))/2;
			  FPlusTexHadr << "\\hline" << endl;
			  FPlusTexHadr << PresentationOutput[SumbTag] << " & " << minuitFitterHadr.GetFRResult() << " & " << minuitFitterHadr.GetFRResult()-SMfrResult <<" & " << minuitFitterHadr.GetFRError() << " & " << FRJESSystHadr << " & " << FRJERSystHadr << " & " << FRWSystHadr << " & " << FRTTScalingHadr << " & " << FRTTMatchingHadr << " & " << FRPUSystHadr << " & " << FRUnclusEnergySystHadr << " & " << FRTopMassSystHadr << " & " << FRTriggEvtSelSystHadr << " & " << sqrt(FRJESSystHadr*FRJESSystHadr + FRJERSystHadr*FRJERSystHadr + FRWSystHadr*FRWSystHadr + FRTTScalingHadr*FRTTScalingHadr + FRTTMatchingHadr*FRTTMatchingHadr + FRPUSystHadr*FRPUSystHadr + FRUnclusEnergySystHadr*FRUnclusEnergySystHadr + FRTopMassSystHadr*FRTopMassSystHadr + FRTriggEvtSelSystHadr*FRTriggEvtSelSystHadr) << " & " << sqrt(FRJESSystHadr*FRJESSystHadr + FRJERSystHadr*FRJERSystHadr + FRWSystHadr*FRWSystHadr + FRTTScalingHadr*FRTTScalingHadr + FRTTMatchingHadr*FRTTMatchingHadr + FRPUSystHadr*FRPUSystHadr + FRUnclusEnergySystHadr*FRUnclusEnergySystHadr + FRTopMassSystHadr*FRTopMassSystHadr + FRTriggEvtSelSystHadr*FRTriggEvtSelSystHadr + (minuitFitterHadr.GetFRError())*(minuitFitterHadr.GetFRError())) << "\\\\" << endl;

			  //F0 results:
			  float F0JESSystHadr = (abs(minuitFitterHadr.GetF0Result() - minuitFitterHadrJESPlus.GetF0Result()) + abs(minuitFitterHadr.GetF0Result() - minuitFitterHadrJESMinus.GetF0Result()))/2;
			  float F0JERSystHadr = 0;
			  if(JERResults == true) F0JERSystHadr = (abs(minuitFitterHadr.GetF0Result() - minuitFitterHadrJERPlus.GetF0Result()) + abs(minuitFitterHadr.GetF0Result() - minuitFitterHadrJERMinus.GetF0Result()))/2;
			  float F0WSystHadr = 0;
			  if(WSystResults == true) F0WSystHadr = (abs(minuitFitterHadr.GetF0Result() - minuitFitterHadrWPlus.GetF0Result()) + abs(minuitFitterHadr.GetF0Result() - minuitFitterHadrWMinus.GetF0Result()))/2;
			  float F0TTScalingHadr = 0;
			  if(TTScalingResults == true) F0TTScalingHadr = (abs(minuitFitterHadr.GetF0Result() - minuitFitterHadrTTScalingUp.GetF0Result()) + abs(minuitFitterHadr.GetF0Result() - minuitFitterHadrTTScalingDown.GetF0Result()))/2;
			  float F0TTMatchingHadr = 0;
			  if(TTMatchingResults == true) F0TTMatchingHadr = (abs(minuitFitterHadr.GetF0Result() - minuitFitterHadrTTMatchingUp.GetF0Result()) + abs(minuitFitterHadr.GetF0Result() - minuitFitterHadrTTMatchingDown.GetF0Result()))/2;		  
			  float F0PUSystHadr = 0;
			  if(PUResults == true) F0PUSystHadr = (abs(minuitFitterHadr.GetF0Result() - minuitFitterHadrPUPlus.GetF0Result()) + abs(minuitFitterHadr.GetF0Result() - minuitFitterHadrPUMinus.GetF0Result()))/2;
			  float F0UnclusEnergySystHadr = 0;
			  if(UnclusEnergyResults == true) F0UnclusEnergySystHadr = (abs(minuitFitterHadr.GetF0Result() - minuitFitterHadrUnclusEnergyPlus.GetF0Result()) + abs(minuitFitterHadr.GetF0Result() - minuitFitterHadrUnclusEnergyMinus.GetF0Result()))/2;
			  float F0TopMassSystHadr = 0;
			  if(TopMassResults == true) F0TopMassSystHadr = (abs(minuitFitterHadr.GetF0Result() - minuitFitterHadrTopMassPlus.GetF0Result()) + abs(minuitFitterHadr.GetF0Result() - minuitFitterHadrTopMassMinus.GetF0Result()))/2;
			  float F0TriggEvtSelSystHadr = 0;
			  if(TriggEvtSelResults == true) F0TriggEvtSelSystHadr = (abs(minuitFitterHadr.GetF0Result() - minuitFitterHadrTriggEvtSelPlus.GetF0Result()) + abs(minuitFitterHadr.GetF0Result() - minuitFitterHadrTriggEvtSelMinus.GetF0Result()))/2;
			  FZeroTexHadr << " \\hline" << endl;
			  FZeroTexHadr << PresentationOutput[SumbTag] << " & " << minuitFitterHadr.GetF0Result() << " & " << minuitFitterHadr.GetF0Result()-SMfrResult <<" & " << minuitFitterHadr.GetF0Error() << " & " << F0JESSystHadr << " & " << F0JERSystHadr << " & " << F0WSystHadr << " & " << F0TTScalingHadr << " & " << F0TTMatchingHadr << " & " << F0PUSystHadr << " & " << F0UnclusEnergySystHadr << " & " << F0TopMassSystHadr << " & " << F0TriggEvtSelSystHadr << " & " << sqrt(F0JESSystHadr*F0JESSystHadr + F0JERSystHadr*F0JERSystHadr + F0WSystHadr*F0WSystHadr + F0TTScalingHadr*F0TTScalingHadr + F0TTMatchingHadr*F0TTMatchingHadr + F0PUSystHadr*F0PUSystHadr + F0UnclusEnergySystHadr*F0UnclusEnergySystHadr + F0TopMassSystHadr*F0TopMassSystHadr + F0TriggEvtSelSystHadr*F0TriggEvtSelSystHadr) << " & " << sqrt(F0JESSystHadr*F0JESSystHadr + F0JERSystHadr*F0JERSystHadr + F0WSystHadr*F0WSystHadr + F0TTScalingHadr*F0TTScalingHadr + F0TTMatchingHadr*F0TTMatchingHadr + F0PUSystHadr*F0PUSystHadr + F0UnclusEnergySystHadr*F0UnclusEnergySystHadr + F0TopMassSystHadr*F0TopMassSystHadr + F0TriggEvtSelSystHadr*F0TriggEvtSelSystHadr + (minuitFitterHadr.GetF0Error())*(minuitFitterHadr.GetF0Error())) << "\\\\" << endl;

			  //Statistical results and number of events:
			  PresentationTexHadr << minuitFitterHadr.GetFRResult() << " & " << minuitFitterHadr.GetFRResult()-SMfrResult <<" & " << minuitFitterHadr.GetFRError() << " & " << minuitFitterHadr.GetFLResult() << " & " << minuitFitterHadr.GetFLResult()-SMflResult << " & " << minuitFitterHadr.GetFLError() << " & " << minuitFitterHadr.GetF0Result() << " & " << minuitFitterHadr.GetF0Result()-SMf0Result << " & " << minuitFitterHadr.GetF0Error() << " \\\\ " << endl;		
			  PresentationTexHadr << " \\hline " << endl;
			}
		      }
		    
		      if(HadronicAndLeptonicWOnly == true){
			//Data:       
			std::string CosThetaDataStringHadrAndLeptW = "CosThetaDataHadrAndLeptWOnly_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//JES Plus:
			std::string CosThetaJESPlusStringHadrAndLeptW = "CosThetaJESPlusHadrAndLeptW_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//JES Minus:
			std::string CosThetaJESMinusStringHadrAndLeptW="CosThetaJESMinusHadrAndLeptW_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//JER Plus:
			std::string CosThetaJERPlusStringHadrAndLeptW = "CosThetaJERPlusHadrAndLeptW_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//JER Minus:
			std::string CosThetaJERMinusStringHadrAndLeptW="CosThetaJERMinusHadrAndLeptW_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//WSyst Plus:
			std::string CosThetaWPlusStringHadrAndLeptW = "CosThetaWPlusHadrAndLeptW_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//WSyst Minus:
			std::string CosThetaWMinusStringHadrAndLeptW = "CosThetaWMinusHadrAndLeptW_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//TTScaling Plus:
			std::string CosThetaTTScalingUpStringHadrAndLeptW = "CosThetaTTScalingUpHadrAndLeptW_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//TTScaling Minus:
			std::string CosThetaTTScalingDownStringHadrAndLeptW = "CosThetaTTScalingDownHadrAndLeptW_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//TTMatching Plus:
			std::string CosThetaTTMatchingUpStringHadrAndLeptW = "CosThetaTTMatchingUpHadrAndLeptW_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//TTMatching Minus:
			std::string CosThetaTTMatchingDownStringHadrAndLeptW = "CosThetaTTMatchingDownHadr_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//PU Plus:
			std::string CosThetaPUPlusStringHadrAndLeptW = "CosThetaPUPlusHadrAndLeptW_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//PU Minus:
			std::string CosThetaPUMinusStringHadrAndLeptW = "CosThetaPUMinusHadrAndLeptW_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//UnclusEnergy Plus:
			std::string CosThetaUnclusEnergyPlusStringHadrAndLeptW = "CosThetaUnclusEnergyPlusHadrAndLeptW_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//UnclusEnergy Minus:
			std::string CosThetaUnclusEnergyMinusStringHadrAndLeptW = "CosThetaUnclusEnergyMinusHadrAndLeptW_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//TopMass Plus:
			std::string CosThetaTopMassPlusStringHadrAndLeptW = "CosThetaTopMassPlusHadrAndLeptW_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//TopMass Minus:
			std::string CosThetaTopMassMinusStringHadrAndLeptW = "CosThetaTopMassMinusHadrAndLeptW_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//TriggEvtSel Plus:
			std::string CosThetaTriggEvtSelPlusStringHadrAndLeptW = "CosThetaTriggEvtSelPlusHadrAndLeptW_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//TriggEvtSel Minus:
			std::string CosThetaTriggEvtSelMinusStringHadrAndLeptW = "CosThetaTriggEvtSelMinusHadrAndLeptW_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			std::string CosThetaSignalStringHadrAndLeptW="CosThetaSignalHadrAndLeptWOnly_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			std::string CosThetaBckgStringHadrAndLeptW = "CosThetaBckgHadrAndLeptWOnly_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];

			if(iDataSet == 0){		      
			  //Data:
			  histo1D[CosThetaDataStringHadrAndLeptW]=new TH1F(CosThetaDataStringHadrAndLeptW.c_str(),CosThetaDataStringHadrAndLeptW.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaDataStringHadrAndLeptW]->SetDirectory(0);
			  //JESPlus:
			  histo1D[CosThetaJESPlusStringHadrAndLeptW]=new TH1F(CosThetaJESPlusStringHadrAndLeptW.c_str(),CosThetaJESPlusStringHadrAndLeptW.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaJESPlusStringHadrAndLeptW]->SetDirectory(0);
			  //JESMinus:
			  histo1D[CosThetaJESMinusStringHadrAndLeptW]=new TH1F(CosThetaJESMinusStringHadrAndLeptW.c_str(),CosThetaJESMinusStringHadrAndLeptW.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaJESMinusStringHadrAndLeptW]->SetDirectory(0);
			  //JERPlus:
			  histo1D[CosThetaJERPlusStringHadrAndLeptW]=new TH1F(CosThetaJERPlusStringHadrAndLeptW.c_str(),CosThetaJERPlusStringHadrAndLeptW.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaJERPlusStringHadrAndLeptW]->SetDirectory(0);
			  //JERMinus:
			  histo1D[CosThetaJERMinusStringHadrAndLeptW]=new TH1F(CosThetaJERMinusStringHadrAndLeptW.c_str(),CosThetaJERMinusStringHadrAndLeptW.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaJERMinusStringHadrAndLeptW]->SetDirectory(0);
			  //WJets Plus:
			  histo1D[CosThetaWPlusStringHadrAndLeptW]=new TH1F(CosThetaWPlusStringHadrAndLeptW.c_str(),CosThetaWPlusStringHadrAndLeptW.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaWPlusStringHadrAndLeptW]->SetDirectory(0);
			  //WJets Minus:
			  histo1D[CosThetaWMinusStringHadrAndLeptW]=new TH1F(CosThetaWMinusStringHadrAndLeptW.c_str(),CosThetaWMinusStringHadrAndLeptW.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaWMinusStringHadrAndLeptW]->SetDirectory(0);
			  //TTScaling Up:
			  histo1D[CosThetaTTScalingUpStringHadrAndLeptW]=new TH1F(CosThetaTTScalingUpStringHadrAndLeptW.c_str(),CosThetaTTScalingUpStringHadrAndLeptW.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaTTScalingUpStringHadrAndLeptW]->SetDirectory(0);
			  //TTScaling Down:
			  histo1D[CosThetaTTScalingDownStringHadrAndLeptW]=new TH1F(CosThetaTTScalingDownStringHadrAndLeptW.c_str(),CosThetaTTScalingDownStringHadrAndLeptW.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaTTScalingDownStringHadrAndLeptW]->SetDirectory(0);
			  //TTMatching Up:
			  histo1D[CosThetaTTMatchingUpStringHadrAndLeptW]=new TH1F(CosThetaTTMatchingUpStringHadrAndLeptW.c_str(),CosThetaTTMatchingUpStringHadrAndLeptW.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaTTMatchingUpStringHadrAndLeptW]->SetDirectory(0);
			  //TTMatching Down:
			  histo1D[CosThetaTTMatchingDownStringHadrAndLeptW]=new TH1F(CosThetaTTMatchingDownStringHadrAndLeptW.c_str(),CosThetaTTMatchingDownStringHadrAndLeptW.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaTTMatchingDownStringHadrAndLeptW]->SetDirectory(0);
			  //PU Plus:
			  histo1D[CosThetaPUPlusStringHadrAndLeptW]=new TH1F(CosThetaPUPlusStringHadrAndLeptW.c_str(),CosThetaPUPlusStringHadrAndLeptW.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaPUPlusStringHadrAndLeptW]->SetDirectory(0);
			  //PU Minus:
			  histo1D[CosThetaPUMinusStringHadrAndLeptW]=new TH1F(CosThetaPUMinusStringHadrAndLeptW.c_str(),CosThetaPUMinusStringHadrAndLeptW.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaPUMinusStringHadrAndLeptW]->SetDirectory(0);
			  //UnclusEnergy Plus:
			  histo1D[CosThetaUnclusEnergyPlusStringHadrAndLeptW]=new TH1F(CosThetaUnclusEnergyPlusStringHadrAndLeptW.c_str(),CosThetaUnclusEnergyPlusStringHadrAndLeptW.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaUnclusEnergyPlusStringHadrAndLeptW]->SetDirectory(0);
			  //UnclusEnergy Minus:
			  histo1D[CosThetaUnclusEnergyMinusStringHadrAndLeptW]=new TH1F(CosThetaUnclusEnergyMinusStringHadrAndLeptW.c_str(),CosThetaUnclusEnergyMinusStringHadrAndLeptW.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaUnclusEnergyMinusStringHadrAndLeptW]->SetDirectory(0);
			  //TopMass Plus:
			  histo1D[CosThetaTopMassPlusStringHadrAndLeptW]=new TH1F(CosThetaTopMassPlusStringHadrAndLeptW.c_str(),CosThetaTopMassPlusStringHadrAndLeptW.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaTopMassPlusStringHadrAndLeptW]->SetDirectory(0);
			  //TopMass Minus:
			  histo1D[CosThetaTopMassMinusStringHadrAndLeptW]=new TH1F(CosThetaTopMassMinusStringHadrAndLeptW.c_str(),CosThetaTopMassMinusStringHadrAndLeptW.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaTopMassMinusStringHadrAndLeptW]->SetDirectory(0);
			  //TriggEvtSel Plus:
			  histo1D[CosThetaTriggEvtSelPlusStringHadrAndLeptW]=new TH1F(CosThetaTriggEvtSelPlusStringHadrAndLeptW.c_str(),CosThetaTriggEvtSelPlusStringHadrAndLeptW.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaTriggEvtSelPlusStringHadrAndLeptW]->SetDirectory(0);
			  //TriggEvtSel Minus:
			  histo1D[CosThetaTriggEvtSelMinusStringHadrAndLeptW]=new TH1F(CosThetaTriggEvtSelMinusStringHadrAndLeptW.c_str(),CosThetaTriggEvtSelMinusStringHadrAndLeptW.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaTriggEvtSelMinusStringHadrAndLeptW]->SetDirectory(0);
			  histo1D[CosThetaSignalStringHadrAndLeptW]=new TH1F(CosThetaSignalStringHadrAndLeptW.c_str(),CosThetaSignalStringHadrAndLeptW.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaSignalStringHadrAndLeptW]->SetDirectory(0);
			  histo1D[CosThetaBckgStringHadrAndLeptW]=new TH1F(CosThetaBckgStringHadrAndLeptW.c_str(),CosThetaBckgStringHadrAndLeptW.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaBckgStringHadrAndLeptW]->SetDirectory(0);
			}

			//Filling of the histograms:		    
			for(int ii=0; ii<CosThetaValuesKinFitHadrAndLeptWOnly[SumbTag].size();ii++){		
		      
			  if((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ){//Defining of genttbar histo
			    for(int iBin=0; iBin< CosThetaBinNumber; iBin++){//Filling of genttbar histo:		   
			      if(CosThetaValuesKinFitHadrAndLeptWOnly[SumbTag][ii] >= binEdge[iBin] && CosThetaValuesKinFitHadrAndLeptWOnly[SumbTag][ii] < binEdge[iBin+1]){
				genttbarhistoHadrAndLeptWOnly[iBin]->Fill(CosThGenKinFitHadrAndLeptWOnly[SumbTag][ii], EventCorrectionWeightKinFitHadrAndLeptWOnly[SumbTag][ii]) ; 
			      }
			      else if(CosThetaValuesKinFitHadrAndLeptWOnly[SumbTag][ii] ==1){ //1 is included in last bin
				genttbarhistoHadrAndLeptWOnly[CosThetaBinNumber-1]->Fill(CosThGenKinFitHadrAndLeptWOnly[SumbTag][ii], EventCorrectionWeightKinFitHadrAndLeptWOnly[SumbTag][ii]);
			      }	      	   
			    }
			  }
		  
			  ///////////////////////////////////////////////////////////////////////////////
			  // Change data to systematics since then nominal values will be reweighted!! //
			  ///////////////////////////////////////////////////////////////////////////////

			  //Data result:
			  if(dataSetName.find("Data") == 0)
			    histo1D[CosThetaDataStringHadrAndLeptW]->Fill(CosThetaValuesKinFitHadrAndLeptWOnly[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLeptWOnly[SumbTag][ii]);
			  //JES Plus:
			  if(dataSetName.find("JESPlus") == 0)
			    histo1D[CosThetaJESPlusStringHadrAndLeptW]->Fill(CosThetaValuesKinFitHadrAndLeptWOnly[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLeptWOnly[SumbTag][ii]);
			  //JES Minus:
			  if(dataSetName.find("JESMinus") == 0)
			    histo1D[CosThetaJESMinusStringHadrAndLeptW]->Fill(CosThetaValuesKinFitHadrAndLeptWOnly[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLeptWOnly[SumbTag][ii]);
			  //JER Plus:
			  if(dataSetName.find("JERPlus") == 0)
			    histo1D[CosThetaJERPlusStringHadrAndLeptW]->Fill(CosThetaValuesKinFitHadrAndLeptWOnly[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLeptWOnly[SumbTag][ii]);
			  //JER Minus:
			  if(dataSetName.find("JERMinus") == 0)
			    histo1D[CosThetaJERMinusStringHadrAndLeptW]->Fill(CosThetaValuesKinFitHadrAndLeptWOnly[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLeptWOnly[SumbTag][ii]);
			  //WJets Plus:  --> Fill with WsystPlus sample and all other MC
			  if( dataSetName.find("Data") !=0 && dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSystMinus") != 0  && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatching") != 0 && dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMass") != 0 && dataSetName.find("TriggEvtSel") != 0)
			    histo1D[CosThetaWPlusStringHadrAndLeptW]->Fill(CosThetaValuesKinFitHadrAndLeptWOnly[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLeptWOnly[SumbTag][ii]);
			  //WJets Minus:  --> Fill with WsystMinus sample and all other MC
			  if( dataSetName.find("Data") !=0 && dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSystPlus") != 0  && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatching") != 0 && dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMass") != 0 && dataSetName.find("TriggEvtSel") != 0)
			    histo1D[CosThetaWMinusStringHadrAndLeptW]->Fill(CosThetaValuesKinFitHadrAndLeptWOnly[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLeptWOnly[SumbTag][ii]);
			  //TTScaling Up:
			  if( dataSetName.find("Data") !=0 && dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSyst") != 0 && dataSetName.find("TTScalingDown") != 0 && dataSetName.find("TTMatching") != 0 && dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMass") != 0 && dataSetName.find("TriggEvtSel") != 0 && dataSetName.find("TTbarJets") !=0)
			    histo1D[CosThetaTTScalingUpStringHadrAndLeptW]->Fill(CosThetaValuesKinFitHadrAndLeptWOnly[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLeptWOnly[SumbTag][ii]);
			  //TTScaling Down:
			  if( dataSetName.find("Data") !=0 && dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSyst") != 0 && dataSetName.find("TTScalingUp") != 0 && dataSetName.find("TTMatching") != 0 && dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMass") != 0 && dataSetName.find("TriggEvtSel") != 0 && dataSetName.find("TTbarJets") !=0 )
			    histo1D[CosThetaTTScalingDownStringHadrAndLeptW]->Fill(CosThetaValuesKinFitHadrAndLeptWOnly[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLeptWOnly[SumbTag][ii]);
			  //TTMatching Up:
			  if(dataSetName.find("Data") !=0 && dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSyst") != 0 && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatchingDown") != 0 && dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMass") != 0 && dataSetName.find("TriggEvtSel") != 0 && dataSetName.find("TTbarJets") !=0 )
			    histo1D[CosThetaTTMatchingUpStringHadrAndLeptW]->Fill(CosThetaValuesKinFitHadrAndLeptWOnly[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLeptWOnly[SumbTag][ii]);
			  //TTMatching Down:
			  if( dataSetName.find("Data") !=0 && dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSyst") != 0 && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatchingUp") != 0 && dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMass") != 0  && dataSetName.find("TriggEvtSel") != 0&& dataSetName.find("TTbarJets") !=0 )
			    histo1D[CosThetaTTMatchingDownStringHadrAndLeptW]->Fill(CosThetaValuesKinFitHadrAndLeptWOnly[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLeptWOnly[SumbTag][ii]);		  
			  //PU Plus:
			  if(dataSetName.find("Data") !=0&& dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSyst") != 0 && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatching") != 0 && dataSetName.find("PUMinus") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMass") != 0 && dataSetName.find("TriggEvtSel") != 0 && dataSetName.find("TTbarJets") !=0 )
			    histo1D[CosThetaPUPlusStringHadrAndLeptW]->Fill(CosThetaValuesKinFitHadrAndLeptWOnly[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLeptWOnly[SumbTag][ii]);
			  //PU Minus:
			  if(dataSetName.find("Data") !=0&& dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSyst") != 0 && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatching") != 0 && dataSetName.find("PUPlus") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMass") != 0 && dataSetName.find("TriggEvtSel") != 0 && dataSetName.find("TTbarJets") !=0 )
			    histo1D[CosThetaPUMinusStringHadrAndLeptW]->Fill(CosThetaValuesKinFitHadrAndLeptWOnly[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLeptWOnly[SumbTag][ii]);
			  //UnclusEnergy Plus:
			  if(dataSetName.find("Data") !=0&& dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSyst") != 0 && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatching") != 0 && dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergyMinus") != 0 && dataSetName.find("TopMass") != 0 && dataSetName.find("TriggEvtSel") != 0 && dataSetName.find("TTbarJets") !=0 )
			    histo1D[CosThetaUnclusEnergyPlusStringHadrAndLeptW]->Fill(CosThetaValuesKinFitHadrAndLeptWOnly[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLeptWOnly[SumbTag][ii]);
			  //UnclusEnergy Minus:
			  if(dataSetName.find("Data") !=0&& dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSyst") != 0 && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatching") != 0 && dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergyPlus") != 0 && dataSetName.find("TopMass") != 0 && dataSetName.find("TriggEvtSel") != 0 && dataSetName.find("TTbarJets") !=0 )
			    histo1D[CosThetaUnclusEnergyMinusStringHadrAndLeptW]->Fill(CosThetaValuesKinFitHadrAndLeptWOnly[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLeptWOnly[SumbTag][ii]);
			  //Top Mass Plus:
			  if(dataSetName.find("Data") !=0&& dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSyst") != 0 && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatching") != 0 && dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMassMinus") != 0 && dataSetName.find("TriggEvtSel") != 0 && dataSetName.find("TTbarJets") !=0 )
			    histo1D[CosThetaTopMassPlusStringHadrAndLeptW]->Fill(CosThetaValuesKinFitHadrAndLeptWOnly[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLeptWOnly[SumbTag][ii]);
			  //Top Mass Minus:
			  if(dataSetName.find("Data") !=0&& dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSyst") != 0 && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatching") != 0 && dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMassPlus") != 0 && dataSetName.find("TriggEvtSel") != 0 && dataSetName.find("TTbarJets") !=0 )		      
			    histo1D[CosThetaTopMassMinusStringHadrAndLeptW]->Fill(CosThetaValuesKinFitHadrAndLeptWOnly[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLeptWOnly[SumbTag][ii]);
			  //TriggEvtSel Plus:
			  if(dataSetName.find("Data") !=0&& dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSyst") != 0 && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatching") != 0 && dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMass") != 0 && dataSetName.find("TriggEvtSelMinus") != 0 && dataSetName.find("TTbarJets") !=0 )
			    histo1D[CosThetaTriggEvtSelPlusStringHadrAndLeptW]->Fill(CosThetaValuesKinFitHadrAndLeptWOnly[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLeptWOnly[SumbTag][ii]);
			  //TriggEvtSel Minus:
			  if(dataSetName.find("Data") !=0&& dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSyst") != 0 && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatching") != 0 && dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMass") != 0 && dataSetName.find("TriggEvtSelPlus") != 0 && dataSetName.find("TTbarJets") !=0 )		      
			    histo1D[CosThetaTriggEvtSelMinusStringHadrAndLeptW]->Fill(CosThetaValuesKinFitHadrAndLeptWOnly[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLeptWOnly[SumbTag][ii]);

			  if(SignalOnly == true && ((dataSetName.find("TTbarJets_SemiMu") == 0 && semiMuon==true) || (dataSetName.find("TTbarJets_SemiEl")==0 && semiElectron == true)))
			    histo1D[CosThetaDataStringHadrAndLeptW]->Fill(CosThetaValuesKinFitHadrAndLeptWOnly[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLeptWOnly[SumbTag][ii]);

			  if(((dataSetName.find("TTbarJets_SemiMu")==0 && semiMuon==true) || (dataSetName.find("TTbarJets_SemiEl")==0 && semiElectron == true) ) && dataSetName.find("JES_") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSyst") != 0 && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatching") != 0 && dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMass") != 0 && dataSetName.find("TriggEvtSel") != 0)
			    histo1D[CosThetaSignalStringHadrAndLeptW]->Fill(CosThetaValuesKinFitHadrAndLeptWOnly[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLeptWOnly[SumbTag][ii]);

			  if(dataSetName.find("Data") !=0 && dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSyst") != 0 && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatching") != 0 && dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMass") != 0 && dataSetName.find("TriggEvtSel") != 0 && ((dataSetName.find("TTbarJets_SemiMu")!=0 && semiMuon==true) || (dataSetName.find("TTbarJets_SemiEl")!=0 && semiElectron == true) )) 
			    histo1D[CosThetaBckgStringHadrAndLeptW]->Fill(CosThetaValuesKinFitHadrAndLeptWOnly[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLeptWOnly[SumbTag][ii]);	  
			}	

			//Perform MinuitFitter:
			if(iDataSet==(datasets.size()-1)){//Go in this loop when the last datasample is active to perform MinuitFitter
			  //Data:
			  MinuitFitter minuitFitterHadrAndLeptWOnly = MinuitFitter(histo1D[CosThetaDataStringHadrAndLeptW], histo1D[CosThetaSignalStringHadrAndLeptW], histo1D[CosThetaBckgStringHadrAndLeptW], 0.6671, 0.3325, 0.0004,genttbarhistoHadrAndLeptWOnly,ndimen);
			  //JESPlus
			  MinuitFitter minuitFitterHadrAndLeptWOnlyJESPlus = MinuitFitter(histo1D[CosThetaJESPlusStringHadrAndLeptW], histo1D[CosThetaSignalStringHadrAndLeptW], histo1D[CosThetaBckgStringHadrAndLeptW], 0.6671, 0.3325, 0.0004,genttbarhistoHadrAndLeptWOnly,ndimen);			  
			  //JES Minus
			  MinuitFitter minuitFitterHadrAndLeptWOnlyJESMinus = MinuitFitter(histo1D[CosThetaJESMinusStringHadrAndLeptW], histo1D[CosThetaSignalStringHadrAndLeptW], histo1D[CosThetaBckgStringHadrAndLeptW], 0.6671, 0.3325, 0.0004,genttbarhistoHadrAndLeptWOnly,ndimen);			  
			  //JERPlus
			  MinuitFitter minuitFitterHadrAndLeptWOnlyJERPlus = MinuitFitter(histo1D[CosThetaJERPlusStringHadrAndLeptW], histo1D[CosThetaSignalStringHadrAndLeptW], histo1D[CosThetaBckgStringHadrAndLeptW], 0.6671, 0.3325, 0.0004,genttbarhistoHadrAndLeptWOnly,ndimen);			  
			  //JER Minus
			  MinuitFitter minuitFitterHadrAndLeptWOnlyJERMinus = MinuitFitter(histo1D[CosThetaJERMinusStringHadrAndLeptW], histo1D[CosThetaSignalStringHadrAndLeptW], histo1D[CosThetaBckgStringHadrAndLeptW], 0.6671, 0.3325, 0.0004,genttbarhistoHadrAndLeptWOnly,ndimen);			  
			  //W Plus
			  MinuitFitter minuitFitterHadrAndLeptWOnlyWPlus = MinuitFitter(histo1D[CosThetaWPlusStringHadrAndLeptW], histo1D[CosThetaSignalStringHadrAndLeptW], histo1D[CosThetaBckgStringHadrAndLeptW], 0.6671, 0.3325, 0.0004,genttbarhistoHadrAndLeptWOnly,ndimen);			  
			  //W Minus
			  MinuitFitter minuitFitterHadrAndLeptWOnlyWMinus = MinuitFitter(histo1D[CosThetaWMinusStringHadrAndLeptW], histo1D[CosThetaSignalStringHadrAndLeptW], histo1D[CosThetaBckgStringHadrAndLeptW], 0.6671, 0.3325, 0.0004,genttbarhistoHadrAndLeptWOnly,ndimen);
			  //TTScaling Up
			  MinuitFitter minuitFitterHadrAndLeptWOnlyTTScalingUp = MinuitFitter(histo1D[CosThetaTTScalingUpStringHadrAndLeptW], histo1D[CosThetaSignalStringHadrAndLeptW], histo1D[CosThetaBckgStringHadrAndLeptW], 0.6671, 0.3325, 0.0004,genttbarhistoHadrAndLeptWOnly,ndimen);
			  //TTScaling Down
			  MinuitFitter minuitFitterHadrAndLeptWOnlyTTScalingDown = MinuitFitter(histo1D[CosThetaTTScalingDownStringHadrAndLeptW], histo1D[CosThetaSignalStringHadrAndLeptW], histo1D[CosThetaBckgStringHadrAndLeptW], 0.6671, 0.3325, 0.0004,genttbarhistoHadrAndLeptWOnly,ndimen);
			  //TTMatching Up
			  MinuitFitter minuitFitterHadrAndLeptWOnlyTTMatchingUp = MinuitFitter(histo1D[CosThetaTTMatchingUpStringHadrAndLeptW], histo1D[CosThetaSignalStringHadrAndLeptW], histo1D[CosThetaBckgStringHadrAndLeptW], 0.6671, 0.3325, 0.0004,genttbarhistoHadrAndLeptWOnly,ndimen);
			  //TTMatching Down
			  MinuitFitter minuitFitterHadrAndLeptWOnlyTTMatchingDown = MinuitFitter(histo1D[CosThetaTTMatchingDownStringHadrAndLeptW], histo1D[CosThetaSignalStringHadrAndLeptW], histo1D[CosThetaBckgStringHadrAndLeptW], 0.6671, 0.3325, 0.0004,genttbarhistoHadrAndLeptWOnly,ndimen);
			  //PU Plus
			  MinuitFitter minuitFitterHadrAndLeptWOnlyPUPlus = MinuitFitter(histo1D[CosThetaPUPlusStringHadrAndLeptW], histo1D[CosThetaSignalStringHadrAndLeptW], histo1D[CosThetaBckgStringHadrAndLeptW], 0.6671, 0.3325, 0.0004,genttbarhistoHadrAndLeptWOnly,ndimen);			  
			  //PU Minus
			  MinuitFitter minuitFitterHadrAndLeptWOnlyPUMinus = MinuitFitter(histo1D[CosThetaPUMinusStringHadrAndLeptW], histo1D[CosThetaSignalStringHadrAndLeptW], histo1D[CosThetaBckgStringHadrAndLeptW], 0.6671, 0.3325, 0.0004,genttbarhistoHadrAndLeptWOnly,ndimen);			  
			  //UnclusEnergy Plus
			  MinuitFitter minuitFitterHadrAndLeptWOnlyUnclusEnergyPlus = MinuitFitter(histo1D[CosThetaUnclusEnergyPlusStringHadrAndLeptW], histo1D[CosThetaSignalStringHadrAndLeptW], histo1D[CosThetaBckgStringHadrAndLeptW], 0.6671, 0.3325, 0.0004,genttbarhistoHadrAndLeptWOnly,ndimen);			  
			  //UnclusEnergy Minus
			  MinuitFitter minuitFitterHadrAndLeptWOnlyUnclusEnergyMinus = MinuitFitter(histo1D[CosThetaUnclusEnergyMinusStringHadrAndLeptW], histo1D[CosThetaSignalStringHadrAndLeptW], histo1D[CosThetaBckgStringHadrAndLeptW], 0.6671, 0.3325, 0.0004,genttbarhistoHadrAndLeptWOnly,ndimen);			  
			  //TopMass Plus
			  MinuitFitter minuitFitterHadrAndLeptWOnlyTopMassPlus = MinuitFitter(histo1D[CosThetaTopMassPlusStringHadrAndLeptW], histo1D[CosThetaSignalStringHadrAndLeptW], histo1D[CosThetaBckgStringHadrAndLeptW], 0.6671, 0.3325, 0.0004,genttbarhistoHadrAndLeptWOnly,ndimen);			  
			  //TopMass Minus
			  MinuitFitter minuitFitterHadrAndLeptWOnlyTopMassMinus = MinuitFitter(histo1D[CosThetaTopMassMinusStringHadrAndLeptW], histo1D[CosThetaSignalStringHadrAndLeptW], histo1D[CosThetaBckgStringHadrAndLeptW], 0.6671, 0.3325, 0.0004,genttbarhistoHadrAndLeptWOnly,ndimen);			  
			  //TriggEvtSel Plus
			  MinuitFitter minuitFitterHadrAndLeptWOnlyTriggEvtSelPlus = MinuitFitter(histo1D[CosThetaTriggEvtSelPlusStringHadrAndLeptW], histo1D[CosThetaSignalStringHadrAndLeptW], histo1D[CosThetaBckgStringHadrAndLeptW], 0.6671, 0.3325, 0.0004,genttbarhistoHadrAndLeptWOnly,ndimen);			  
			  //TriggEvtSel Minus
			  MinuitFitter minuitFitterHadrAndLeptWOnlyTriggEvtSelMinus = MinuitFitter(histo1D[CosThetaTriggEvtSelMinusStringHadrAndLeptW], histo1D[CosThetaSignalStringHadrAndLeptW], histo1D[CosThetaBckgStringHadrAndLeptW], 0.6671, 0.3325, 0.0004,genttbarhistoHadrAndLeptWOnly,ndimen);			  
			  
			  PresentationTexHadrAndLeptWOnly << PresentationOutput[SumbTag] << " & ";
			  cout << " Inside minuitFitter loop : " << PresentationOutput[SumbTag] << " & " << endl;
			  for(int ii=0; ii< nameDataSet.size(); ii++){
			    if(ii < nameDataSet.size()-1){ PresentationTexHadrAndLeptWOnly << NumberRemainingEventsKinFitHadrAndLeptWOnly[SumbTag][ii] << " & ";}
			    else if(ii == nameDataSet.size()-1 ){ 
			      PresentationTexHadrAndLeptWOnly << NumberRemainingEventsKinFitHadrAndLeptWOnly[SumbTag][ii] << " & ";//For presentation helicity values still need to be included !!
			      PresentationTexHadrAndLeptWOnly << NumberBLeptCorrectEventsKinFitHadrAndLeptWOnly[SumbTag][ii] << " & ";
			    }
			  }	
		  
			  //FL results:
			  float FLJESSystHadrAndLeptWOnly = (abs(minuitFitterHadrAndLeptWOnly.GetFLResult() - minuitFitterHadrAndLeptWOnlyJESPlus.GetFLResult()) + abs(minuitFitterHadrAndLeptWOnly.GetFLResult() - minuitFitterHadrAndLeptWOnlyJESMinus.GetFLResult()))/2;
			  float FLJERSystHadrAndLeptWOnly = 0;
			  if(JERResults == true) FLJERSystHadrAndLeptWOnly = (abs(minuitFitterHadrAndLeptWOnly.GetFLResult() - minuitFitterHadrAndLeptWOnlyJERPlus.GetFLResult()) + abs(minuitFitterHadrAndLeptWOnly.GetFLResult() - minuitFitterHadrAndLeptWOnlyJERMinus.GetFLResult()))/2;
			  float FLWSystHadrAndLeptWOnly = 0;
			  if(WSystResults == true) FLWSystHadrAndLeptWOnly = (abs(minuitFitterHadrAndLeptWOnly.GetFLResult() - minuitFitterHadrAndLeptWOnlyWPlus.GetFLResult()) + abs(minuitFitterHadrAndLeptWOnly.GetFLResult() - minuitFitterHadrAndLeptWOnlyWMinus.GetFLResult()))/2;
			  float FLTTScalingHadrAndLeptWOnly = 0;
			  if(TTScalingResults == true) FLTTScalingHadrAndLeptWOnly = (abs(minuitFitterHadrAndLeptWOnly.GetFLResult() - minuitFitterHadrAndLeptWOnlyTTScalingUp.GetFLResult()) + abs(minuitFitterHadrAndLeptWOnly.GetFLResult() - minuitFitterHadrAndLeptWOnlyTTScalingDown.GetFLResult()))/2;
			  float FLTTMatchingHadrAndLeptWOnly = 0;
			  if(TTMatchingResults == true) FLTTMatchingHadrAndLeptWOnly = (abs(minuitFitterHadrAndLeptWOnly.GetFLResult() - minuitFitterHadrAndLeptWOnlyTTMatchingUp.GetFLResult()) + abs(minuitFitterHadrAndLeptWOnly.GetFLResult() - minuitFitterHadrAndLeptWOnlyTTMatchingDown.GetFLResult()))/2;
			  float FLPUSystHadrAndLeptWOnly = 0;
			  if(PUResults == true) FLPUSystHadrAndLeptWOnly = (abs(minuitFitterHadrAndLeptWOnly.GetFLResult() - minuitFitterHadrAndLeptWOnlyPUPlus.GetFLResult()) + abs(minuitFitterHadrAndLeptWOnly.GetFLResult() - minuitFitterHadrAndLeptWOnlyPUMinus.GetFLResult()))/2;
			  float FLUnclusEnergySystHadrAndLeptWOnly = 0;
			  if(UnclusEnergyResults == true) FLUnclusEnergySystHadrAndLeptWOnly = (abs(minuitFitterHadrAndLeptWOnly.GetFLResult() - minuitFitterHadrAndLeptWOnlyUnclusEnergyPlus.GetFLResult()) + abs(minuitFitterHadrAndLeptWOnly.GetFLResult() - minuitFitterHadrAndLeptWOnlyUnclusEnergyMinus.GetFLResult()))/2;
			  float FLTopMassSystHadrAndLeptWOnly = 0;
			  if(TopMassResults == true) FLTopMassSystHadrAndLeptWOnly = (abs(minuitFitterHadrAndLeptWOnly.GetFLResult() - minuitFitterHadrAndLeptWOnlyTopMassPlus.GetFLResult()) + abs(minuitFitterHadrAndLeptWOnly.GetFLResult() - minuitFitterHadrAndLeptWOnlyTopMassMinus.GetFLResult()))/2;		      
			  float FLTriggEvtSelSystHadrAndLeptWOnly = 0;
			  if(TriggEvtSelResults == true) FLTriggEvtSelSystHadrAndLeptWOnly = (abs(minuitFitterHadrAndLeptWOnly.GetFLResult() - minuitFitterHadrAndLeptWOnlyTriggEvtSelPlus.GetFLResult()) + abs(minuitFitterHadrAndLeptWOnly.GetFLResult() - minuitFitterHadrAndLeptWOnlyTriggEvtSelMinus.GetFLResult()))/2;		      
			  FMinusTexHadrAndLeptWOnly << "\\hline" << endl;
			  FMinusTexHadrAndLeptWOnly << PresentationOutput[SumbTag] << " & " << minuitFitterHadrAndLeptWOnly.GetFLResult() << " & " << minuitFitterHadrAndLeptWOnly.GetFLResult()-SMfrResult <<" & " << minuitFitterHadrAndLeptWOnly.GetFLError() << " & " << FLJESSystHadrAndLeptWOnly << " & " << FLJERSystHadrAndLeptWOnly << " & " << FLWSystHadrAndLeptWOnly << " & " << FLTTScalingHadrAndLeptWOnly << " & " << FLTTMatchingHadrAndLeptWOnly << " & " << FLPUSystHadrAndLeptWOnly << " & " << FLUnclusEnergySystHadrAndLeptWOnly << " & " << FLTopMassSystHadrAndLeptWOnly << " & " << FLTriggEvtSelSystHadrAndLeptWOnly << " & " << sqrt(FLJESSystHadrAndLeptWOnly*FLJESSystHadrAndLeptWOnly + FLJERSystHadrAndLeptWOnly*FLJERSystHadrAndLeptWOnly + FLWSystHadrAndLeptWOnly*FLWSystHadrAndLeptWOnly + FLTTScalingHadrAndLeptWOnly*FLTTScalingHadrAndLeptWOnly + FLTTMatchingHadrAndLeptWOnly*FLTTMatchingHadrAndLeptWOnly + FLPUSystHadrAndLeptWOnly*FLPUSystHadrAndLeptWOnly + FLUnclusEnergySystHadrAndLeptWOnly*FLUnclusEnergySystHadrAndLeptWOnly + FLTopMassSystHadrAndLeptWOnly*FLTopMassSystHadrAndLeptWOnly + FLTriggEvtSelSystHadrAndLeptWOnly*FLTriggEvtSelSystHadrAndLeptWOnly) << " & " << sqrt(FLJESSystHadrAndLeptWOnly*FLJESSystHadrAndLeptWOnly + FLJERSystHadrAndLeptWOnly*FLJERSystHadrAndLeptWOnly + FLWSystHadrAndLeptWOnly*FLWSystHadrAndLeptWOnly + FLTTScalingHadrAndLeptWOnly*FLTTScalingHadrAndLeptWOnly + FLTTMatchingHadrAndLeptWOnly*FLTTMatchingHadrAndLeptWOnly + FLPUSystHadrAndLeptWOnly*FLPUSystHadrAndLeptWOnly + FLUnclusEnergySystHadrAndLeptWOnly*FLUnclusEnergySystHadrAndLeptWOnly + FLTopMassSystHadrAndLeptWOnly*FLTopMassSystHadrAndLeptWOnly + FLTriggEvtSelSystHadrAndLeptWOnly*FLTriggEvtSelSystHadrAndLeptWOnly + (minuitFitterHadrAndLeptWOnly.GetFLError())*(minuitFitterHadrAndLeptWOnly.GetFLError())) << "\\\\" << endl;
		  
			  //FR results:
			  float FRJESSystHadrAndLeptWOnly = (abs(minuitFitterHadrAndLeptWOnly.GetFRResult() - minuitFitterHadrAndLeptWOnlyJESPlus.GetFRResult()) + abs(minuitFitterHadrAndLeptWOnly.GetFRResult() - minuitFitterHadrAndLeptWOnlyJESMinus.GetFRResult()))/2;
			  float FRJERSystHadrAndLeptWOnly = 0;
			  if(JERResults == true) FRJERSystHadrAndLeptWOnly = (abs(minuitFitterHadrAndLeptWOnly.GetFRResult() - minuitFitterHadrAndLeptWOnlyJERPlus.GetFRResult()) + abs(minuitFitterHadrAndLeptWOnly.GetFRResult() - minuitFitterHadrAndLeptWOnlyJERMinus.GetFRResult()))/2;
			  float FRWSystHadrAndLeptWOnly = 0;
			  if(WSystResults == true) FRWSystHadrAndLeptWOnly = (abs(minuitFitterHadrAndLeptWOnly.GetFRResult() - minuitFitterHadrAndLeptWOnlyWPlus.GetFRResult()) + abs(minuitFitterHadrAndLeptWOnly.GetFRResult() - minuitFitterHadrAndLeptWOnlyWMinus.GetFRResult()))/2;
			  float FRTTScalingHadrAndLeptWOnly = 0;
			  if(TTScalingResults == true) FRTTScalingHadrAndLeptWOnly = (abs(minuitFitterHadrAndLeptWOnly.GetFRResult() - minuitFitterHadrAndLeptWOnlyTTScalingUp.GetFRResult()) + abs(minuitFitterHadrAndLeptWOnly.GetFRResult() - minuitFitterHadrAndLeptWOnlyTTScalingDown.GetFRResult()))/2;
			  float FRTTMatchingHadrAndLeptWOnly = 0;
			  if(TTMatchingResults == true) FRTTMatchingHadrAndLeptWOnly = (abs(minuitFitterHadrAndLeptWOnly.GetFRResult() - minuitFitterHadrAndLeptWOnlyTTMatchingUp.GetFRResult()) + abs(minuitFitterHadrAndLeptWOnly.GetFRResult() - minuitFitterHadrAndLeptWOnlyTTMatchingDown.GetFRResult()))/2;		  
			  float FRPUSystHadrAndLeptWOnly = 0;
			  if(PUResults == true) FRPUSystHadrAndLeptWOnly = (abs(minuitFitterHadrAndLeptWOnly.GetFRResult() - minuitFitterHadrAndLeptWOnlyPUPlus.GetFRResult()) + abs(minuitFitterHadrAndLeptWOnly.GetFRResult() - minuitFitterHadrAndLeptWOnlyPUMinus.GetFRResult()))/2;
			  float FRUnclusEnergySystHadrAndLeptWOnly = 0;
			  if(UnclusEnergyResults == true) FRUnclusEnergySystHadrAndLeptWOnly = (abs(minuitFitterHadrAndLeptWOnly.GetFRResult() - minuitFitterHadrAndLeptWOnlyUnclusEnergyPlus.GetFRResult()) + abs(minuitFitterHadrAndLeptWOnly.GetFRResult() - minuitFitterHadrAndLeptWOnlyUnclusEnergyMinus.GetFRResult()))/2;
			  float FRTopMassSystHadrAndLeptWOnly = 0;
			  if(TopMassResults == true) FRTopMassSystHadrAndLeptWOnly = (abs(minuitFitterHadrAndLeptWOnly.GetFRResult() - minuitFitterHadrAndLeptWOnlyTopMassPlus.GetFRResult()) + abs(minuitFitterHadrAndLeptWOnly.GetFRResult() - minuitFitterHadrAndLeptWOnlyTopMassMinus.GetFRResult()))/2;
			  float FRTriggEvtSelSystHadrAndLeptWOnly = 0;
			  if(TriggEvtSelResults == true) FRTriggEvtSelSystHadrAndLeptWOnly = (abs(minuitFitterHadrAndLeptWOnly.GetFRResult() - minuitFitterHadrAndLeptWOnlyTriggEvtSelPlus.GetFRResult()) + abs(minuitFitterHadrAndLeptWOnly.GetFRResult() - minuitFitterHadrAndLeptWOnlyTriggEvtSelMinus.GetFRResult()))/2;
			  FPlusTexHadrAndLeptWOnly << "\\hline" << endl;
			  FPlusTexHadrAndLeptWOnly << PresentationOutput[SumbTag] << " & " << minuitFitterHadrAndLeptWOnly.GetFRResult() << " & " << minuitFitterHadrAndLeptWOnly.GetFRResult()-SMfrResult <<" & " << minuitFitterHadrAndLeptWOnly.GetFRError() << " & " << FRJESSystHadrAndLeptWOnly << " & " << FRJERSystHadrAndLeptWOnly << " & " << FRWSystHadrAndLeptWOnly << " & " << FRTTScalingHadrAndLeptWOnly << " & " << FRTTMatchingHadrAndLeptWOnly << " & " << FRPUSystHadrAndLeptWOnly << " & " << FRUnclusEnergySystHadrAndLeptWOnly << " & " << FRTopMassSystHadrAndLeptWOnly << " & " << FRTriggEvtSelSystHadrAndLeptWOnly << " & " << sqrt(FRJESSystHadrAndLeptWOnly*FRJESSystHadrAndLeptWOnly + FRJERSystHadrAndLeptWOnly*FRJERSystHadrAndLeptWOnly + FRWSystHadrAndLeptWOnly*FRWSystHadrAndLeptWOnly + FRTTScalingHadrAndLeptWOnly*FRTTScalingHadrAndLeptWOnly + FRTTMatchingHadrAndLeptWOnly*FRTTMatchingHadrAndLeptWOnly + FRPUSystHadrAndLeptWOnly*FRPUSystHadrAndLeptWOnly + FRUnclusEnergySystHadrAndLeptWOnly*FRUnclusEnergySystHadrAndLeptWOnly + FRTopMassSystHadrAndLeptWOnly*FRTopMassSystHadrAndLeptWOnly + FRTriggEvtSelSystHadrAndLeptWOnly*FRTriggEvtSelSystHadrAndLeptWOnly) << " & " << sqrt(FRJESSystHadrAndLeptWOnly*FRJESSystHadrAndLeptWOnly + FRJERSystHadrAndLeptWOnly*FRJERSystHadrAndLeptWOnly + FRWSystHadrAndLeptWOnly*FRWSystHadrAndLeptWOnly + FRTTScalingHadrAndLeptWOnly*FRTTScalingHadrAndLeptWOnly + FRTTMatchingHadrAndLeptWOnly*FRTTMatchingHadrAndLeptWOnly + FRPUSystHadrAndLeptWOnly*FRPUSystHadrAndLeptWOnly + FRUnclusEnergySystHadrAndLeptWOnly*FRUnclusEnergySystHadrAndLeptWOnly + FRTopMassSystHadrAndLeptWOnly*FRTopMassSystHadrAndLeptWOnly + FRTriggEvtSelSystHadrAndLeptWOnly*FRTriggEvtSelSystHadrAndLeptWOnly + (minuitFitterHadrAndLeptWOnly.GetFRError())*(minuitFitterHadrAndLeptWOnly.GetFRError())) << "\\\\" << endl;

			  //F0 results:
			  float F0JESSystHadrAndLeptWOnly = (abs(minuitFitterHadrAndLeptWOnly.GetF0Result() - minuitFitterHadrAndLeptWOnlyJESPlus.GetF0Result()) + abs(minuitFitterHadrAndLeptWOnly.GetF0Result() - minuitFitterHadrAndLeptWOnlyJESMinus.GetF0Result()))/2;
			  float F0JERSystHadrAndLeptWOnly = 0;
			  if(JERResults == true) F0JERSystHadrAndLeptWOnly = (abs(minuitFitterHadrAndLeptWOnly.GetF0Result() - minuitFitterHadrAndLeptWOnlyJERPlus.GetF0Result()) + abs(minuitFitterHadrAndLeptWOnly.GetF0Result() - minuitFitterHadrAndLeptWOnlyJERMinus.GetF0Result()))/2;
			  float F0WSystHadrAndLeptWOnly = 0;
			  if(WSystResults == true) F0WSystHadrAndLeptWOnly = (abs(minuitFitterHadrAndLeptWOnly.GetF0Result() - minuitFitterHadrAndLeptWOnlyWPlus.GetF0Result()) + abs(minuitFitterHadrAndLeptWOnly.GetF0Result() - minuitFitterHadrAndLeptWOnlyWMinus.GetF0Result()))/2;
			  float F0TTScalingHadrAndLeptWOnly = 0;
			  if(TTScalingResults == true) F0TTScalingHadrAndLeptWOnly = (abs(minuitFitterHadrAndLeptWOnly.GetF0Result() - minuitFitterHadrAndLeptWOnlyTTScalingUp.GetF0Result()) + abs(minuitFitterHadrAndLeptWOnly.GetF0Result() - minuitFitterHadrAndLeptWOnlyTTScalingDown.GetF0Result()))/2;
			  float F0TTMatchingHadrAndLeptWOnly = 0;
			  if(TTMatchingResults == true) F0TTMatchingHadrAndLeptWOnly = (abs(minuitFitterHadrAndLeptWOnly.GetF0Result() - minuitFitterHadrAndLeptWOnlyTTMatchingUp.GetF0Result()) + abs(minuitFitterHadrAndLeptWOnly.GetF0Result() - minuitFitterHadrAndLeptWOnlyTTMatchingDown.GetF0Result()))/2;		  
			  float F0PUSystHadrAndLeptWOnly = 0;
			  if(PUResults == true) F0PUSystHadrAndLeptWOnly = (abs(minuitFitterHadrAndLeptWOnly.GetF0Result() - minuitFitterHadrAndLeptWOnlyPUPlus.GetF0Result()) + abs(minuitFitterHadrAndLeptWOnly.GetF0Result() - minuitFitterHadrAndLeptWOnlyPUMinus.GetF0Result()))/2;
			  float F0UnclusEnergySystHadrAndLeptWOnly = 0;
			  if(UnclusEnergyResults == true) F0UnclusEnergySystHadrAndLeptWOnly = (abs(minuitFitterHadrAndLeptWOnly.GetF0Result() - minuitFitterHadrAndLeptWOnlyUnclusEnergyPlus.GetF0Result()) + abs(minuitFitterHadrAndLeptWOnly.GetF0Result() - minuitFitterHadrAndLeptWOnlyUnclusEnergyMinus.GetF0Result()))/2;
			  float F0TopMassSystHadrAndLeptWOnly = 0;
			  if(TopMassResults == true) F0TopMassSystHadrAndLeptWOnly = (abs(minuitFitterHadrAndLeptWOnly.GetF0Result() - minuitFitterHadrAndLeptWOnlyTopMassPlus.GetF0Result()) + abs(minuitFitterHadrAndLeptWOnly.GetF0Result() - minuitFitterHadrAndLeptWOnlyTopMassMinus.GetF0Result()))/2;
			  float F0TriggEvtSelSystHadrAndLeptWOnly = 0;
			  if(TriggEvtSelResults == true) F0TriggEvtSelSystHadrAndLeptWOnly = (abs(minuitFitterHadrAndLeptWOnly.GetF0Result() - minuitFitterHadrAndLeptWOnlyTriggEvtSelPlus.GetF0Result()) + abs(minuitFitterHadrAndLeptWOnly.GetF0Result() - minuitFitterHadrAndLeptWOnlyTriggEvtSelMinus.GetF0Result()))/2;
			  FZeroTexHadrAndLeptWOnly << " \\hline" << endl;
			  FZeroTexHadrAndLeptWOnly << PresentationOutput[SumbTag] << " & " << minuitFitterHadrAndLeptWOnly.GetF0Result() << " & " << minuitFitterHadrAndLeptWOnly.GetF0Result()-SMfrResult <<" & " << minuitFitterHadrAndLeptWOnly.GetF0Error() << " & " << F0JESSystHadrAndLeptWOnly << " & " << F0JERSystHadrAndLeptWOnly << " & " << F0WSystHadrAndLeptWOnly << " & " << F0TTScalingHadrAndLeptWOnly << " & " << F0TTMatchingHadrAndLeptWOnly << " & " << F0PUSystHadrAndLeptWOnly << " & " << F0UnclusEnergySystHadrAndLeptWOnly << " & " << F0TopMassSystHadrAndLeptWOnly << " & " << F0TriggEvtSelSystHadrAndLeptWOnly << " & " << sqrt(F0JESSystHadrAndLeptWOnly*F0JESSystHadrAndLeptWOnly + F0JERSystHadrAndLeptWOnly*F0JERSystHadrAndLeptWOnly + F0WSystHadrAndLeptWOnly*F0WSystHadrAndLeptWOnly + F0TTScalingHadrAndLeptWOnly*F0TTScalingHadrAndLeptWOnly + F0TTMatchingHadrAndLeptWOnly*F0TTMatchingHadrAndLeptWOnly + F0PUSystHadrAndLeptWOnly*F0PUSystHadrAndLeptWOnly + F0UnclusEnergySystHadrAndLeptWOnly*F0UnclusEnergySystHadrAndLeptWOnly + F0TopMassSystHadrAndLeptWOnly*F0TopMassSystHadrAndLeptWOnly + F0TriggEvtSelSystHadrAndLeptWOnly*F0TriggEvtSelSystHadrAndLeptWOnly) << " & " << sqrt(F0JESSystHadrAndLeptWOnly*F0JESSystHadrAndLeptWOnly + F0JERSystHadrAndLeptWOnly*F0JERSystHadrAndLeptWOnly + F0WSystHadrAndLeptWOnly*F0WSystHadrAndLeptWOnly + F0TTScalingHadrAndLeptWOnly*F0TTScalingHadrAndLeptWOnly + F0TTMatchingHadrAndLeptWOnly*F0TTMatchingHadrAndLeptWOnly + F0PUSystHadrAndLeptWOnly*F0PUSystHadrAndLeptWOnly + F0UnclusEnergySystHadrAndLeptWOnly*F0UnclusEnergySystHadrAndLeptWOnly + F0TopMassSystHadrAndLeptWOnly*F0TopMassSystHadrAndLeptWOnly + F0TriggEvtSelSystHadrAndLeptWOnly*F0TriggEvtSelSystHadrAndLeptWOnly + (minuitFitterHadrAndLeptWOnly.GetF0Error())*(minuitFitterHadrAndLeptWOnly.GetF0Error())) << "\\\\" << endl;

			  //Statistical results and number of events:
			  PresentationTexHadrAndLeptWOnly << minuitFitterHadrAndLeptWOnly.GetFRResult() << " & " << minuitFitterHadrAndLeptWOnly.GetFRResult()-SMfrResult <<" & " << minuitFitterHadrAndLeptWOnly.GetFRError() << " & " << minuitFitterHadrAndLeptWOnly.GetFLResult() << " & " << minuitFitterHadrAndLeptWOnly.GetFLResult()-SMflResult << " & " << minuitFitterHadrAndLeptWOnly.GetFLError() << " & " << minuitFitterHadrAndLeptWOnly.GetF0Result() << " & " << minuitFitterHadrAndLeptWOnly.GetF0Result()-SMf0Result << " & " << minuitFitterHadrAndLeptWOnly.GetF0Error() << " \\\\ " << endl;		
			  PresentationTexHadrAndLeptWOnly << " \\hline " << endl;
			}
		      }

		      if(HadronicAndLeptonic == true){
			//Data:
			std::string CosThetaDataStringHadrAndLept = "CosThetaDataHadrAndLept_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//JES Plus:
			std::string CosThetaJESPlusStringHadrAndLept = "CosThetaJESPlusHadrAndLept_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//JES Minus:
			std::string CosThetaJESMinusStringHadrAndLept="CosThetaJESMinusHadrAndLept_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//JER Plus:
			std::string CosThetaJERPlusStringHadrAndLept = "CosThetaJERPlusHadrAndLept_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//JER Minus:
			std::string CosThetaJERMinusStringHadrAndLept="CosThetaJERMinusHadrAndLept_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//WSyst Plus:
			std::string CosThetaWPlusStringHadrAndLept = "CosThetaWPlusHadrAndLept_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//WSyst Minus:
			std::string CosThetaWMinusStringHadrAndLept = "CosThetaWMinusHadrAndLept_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//TTScaling Plus:
			std::string CosThetaTTScalingUpStringHadrAndLept="CosThetaTTScalingUpHadrAndLept_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//TTScaling Minus:
			std::string CosThetaTTScalingDownStringHadrAndLept = "CosThetaTTScalingDownHadrAndLept_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//TTMatching Plus:
			std::string CosThetaTTMatchingUpStringHadrAndLept = "CosThetaTTMatchingUpHadrAndLept_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//TTMatching Minus:
			std::string CosThetaTTMatchingDownStringHadrAndLept = "CosThetaTTMatchingDownHadrAndLept_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//PU Plus:
			std::string CosThetaPUPlusStringHadrAndLept = "CosThetaPUPlusHadrAndLept_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//PU Minus:
			std::string CosThetaPUMinusStringHadrAndLept = "CosThetaPUMinusHadrAndLept_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//UnclusEnergy Plus:
			std::string CosThetaUnclusEnergyPlusStringHadrAndLept = "CosThetaUnclusEnergyPlusHadrAndLept_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//UnclusEnergy Minus:
			std::string CosThetaUnclusEnergyMinusStringHadrAndLept = "CosThetaUnclusEnergyMinusHadrAndLept_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//TopMass Plus:
			std::string CosThetaTopMassPlusStringHadrAndLept = "CosThetaTopMassPlusHadrAndLept_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//TopMass Minus:
			std::string CosThetaTopMassMinusStringHadrAndLept = "CosThetaTopMassMinusHadrAndLept_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//TriggEvtSel Plus:
			std::string CosThetaTriggEvtSelPlusStringHadrAndLept = "CosThetaTriggEvtSelPlusHadrAndLept_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			//TriggEvtSel Minus:
			std::string CosThetaTriggEvtSelMinusStringHadrAndLept = "CosThetaTriggEvtSelMinusHadrAndLept_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			std::string CosThetaSignalStringHadrAndLept = "CosThetaSignalHadrAndLept_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
			std::string CosThetaBckgStringHadrAndLept = "CosThetaBckgHadrAndLept_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];		    
		    
			if(iDataSet == 0){		      

			  histo1D[CosThetaDataStringHadrAndLept]=new TH1F(CosThetaDataStringHadrAndLept.c_str(),CosThetaDataStringHadrAndLept.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaDataStringHadrAndLept]->SetDirectory(0);
			  //JESPlus:
			  histo1D[CosThetaJESPlusStringHadrAndLept]=new TH1F(CosThetaJESPlusStringHadrAndLept.c_str(),CosThetaJESPlusStringHadrAndLept.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaJESPlusStringHadrAndLept]->SetDirectory(0);
			  //JESMinus:
			  histo1D[CosThetaJESMinusStringHadrAndLept]=new TH1F(CosThetaJESMinusStringHadrAndLept.c_str(),CosThetaJESMinusStringHadrAndLept.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaJESMinusStringHadrAndLept]->SetDirectory(0);
			  //JERPlus:
			  histo1D[CosThetaJERPlusStringHadrAndLept]=new TH1F(CosThetaJERPlusStringHadrAndLept.c_str(),CosThetaJERPlusStringHadrAndLept.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaJERPlusStringHadrAndLept]->SetDirectory(0);
			  //JERMinus:
			  histo1D[CosThetaJERMinusStringHadrAndLept]=new TH1F(CosThetaJERMinusStringHadrAndLept.c_str(),CosThetaJERMinusStringHadrAndLept.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaJERMinusStringHadrAndLept]->SetDirectory(0);
			  //WJets Plus:
			  histo1D[CosThetaWPlusStringHadrAndLept]=new TH1F(CosThetaWPlusStringHadrAndLept.c_str(),CosThetaWPlusStringHadrAndLept.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaWPlusStringHadrAndLept]->SetDirectory(0);
			  //WJets Minus:
			  histo1D[CosThetaWMinusStringHadrAndLept]=new TH1F(CosThetaWMinusStringHadrAndLept.c_str(),CosThetaWMinusStringHadrAndLept.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaWMinusStringHadrAndLept]->SetDirectory(0);
			  //TTScaling Up:
			  histo1D[CosThetaTTScalingUpStringHadrAndLept]=new TH1F(CosThetaTTScalingUpStringHadrAndLept.c_str(),CosThetaTTScalingUpStringHadrAndLept.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaTTScalingUpStringHadrAndLept]->SetDirectory(0);
			  //TTScaling Down:
			  histo1D[CosThetaTTScalingDownStringHadrAndLept]=new TH1F(CosThetaTTScalingDownStringHadrAndLept.c_str(),CosThetaTTScalingDownStringHadrAndLept.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaTTScalingDownStringHadrAndLept]->SetDirectory(0);
			  //TTMatching Up:
			  histo1D[CosThetaTTMatchingUpStringHadrAndLept]=new TH1F(CosThetaTTMatchingUpStringHadrAndLept.c_str(),CosThetaTTMatchingUpStringHadrAndLept.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaTTMatchingUpStringHadrAndLept]->SetDirectory(0);
			  //TTMatching Down:
			  histo1D[CosThetaTTMatchingDownStringHadrAndLept]=new TH1F(CosThetaTTMatchingDownStringHadrAndLept.c_str(),CosThetaTTMatchingDownStringHadrAndLept.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaTTMatchingDownStringHadrAndLept]->SetDirectory(0);
			  //PU Plus:
			  histo1D[CosThetaPUPlusStringHadrAndLept]=new TH1F(CosThetaPUPlusStringHadrAndLept.c_str(),CosThetaPUPlusStringHadrAndLept.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaPUPlusStringHadrAndLept]->SetDirectory(0);
			  //PU Minus:
			  histo1D[CosThetaPUMinusStringHadrAndLept]=new TH1F(CosThetaPUMinusStringHadrAndLept.c_str(),CosThetaPUMinusStringHadrAndLept.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaPUMinusStringHadrAndLept]->SetDirectory(0);
			  //UnclusEnergy Plus:
			  histo1D[CosThetaUnclusEnergyPlusStringHadrAndLept]=new TH1F(CosThetaUnclusEnergyPlusStringHadrAndLept.c_str(),CosThetaUnclusEnergyPlusStringHadrAndLept.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaUnclusEnergyPlusStringHadrAndLept]->SetDirectory(0);
			  //UnclusEnergy Minus:
			  histo1D[CosThetaUnclusEnergyMinusStringHadrAndLept]=new TH1F(CosThetaUnclusEnergyMinusStringHadrAndLept.c_str(),CosThetaUnclusEnergyMinusStringHadrAndLept.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaUnclusEnergyMinusStringHadrAndLept]->SetDirectory(0);
			  //TopMass Plus:
			  histo1D[CosThetaTopMassPlusStringHadrAndLept]=new TH1F(CosThetaTopMassPlusStringHadrAndLept.c_str(),CosThetaTopMassPlusStringHadrAndLept.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaTopMassPlusStringHadrAndLept]->SetDirectory(0);
			  //TopMass Minus:
			  histo1D[CosThetaTopMassMinusStringHadrAndLept]=new TH1F(CosThetaTopMassMinusStringHadrAndLept.c_str(),CosThetaTopMassMinusStringHadrAndLept.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaTopMassMinusStringHadrAndLept]->SetDirectory(0);
			  //TriggEvtSel Plus:
			  histo1D[CosThetaTriggEvtSelPlusStringHadrAndLept]=new TH1F(CosThetaTriggEvtSelPlusStringHadrAndLept.c_str(),CosThetaTriggEvtSelPlusStringHadrAndLept.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaTriggEvtSelPlusStringHadrAndLept]->SetDirectory(0);
			  //TriggEvtSel Minus:
			  histo1D[CosThetaTriggEvtSelMinusStringHadrAndLept]=new TH1F(CosThetaTriggEvtSelMinusStringHadrAndLept.c_str(),CosThetaTriggEvtSelMinusStringHadrAndLept.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaTriggEvtSelMinusStringHadrAndLept]->SetDirectory(0);
			  histo1D[CosThetaSignalStringHadrAndLept]=new TH1F(CosThetaSignalStringHadrAndLept.c_str(),CosThetaSignalStringHadrAndLept.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaSignalStringHadrAndLept]->SetDirectory(0);
			  histo1D[CosThetaBckgStringHadrAndLept]=new TH1F(CosThetaBckgStringHadrAndLept.c_str(),CosThetaBckgStringHadrAndLept.c_str(),CosThetaBinNumber,-1,1);
			  histo1D[CosThetaBckgStringHadrAndLept]->SetDirectory(0);
			}

			//Filling of the histograms:
			for(int ii=0; ii<CosThetaValuesKinFitHadrAndLept[SumbTag].size();ii++){				      
			  if((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ){//Defining of genttbar histo
			    for(int iBin=0; iBin< CosThetaBinNumber; iBin++){//Filling of genttbar histo:		   
			      if(CosThetaValuesKinFitHadrAndLept[SumbTag][ii] >= binEdge[iBin] && CosThetaValuesKinFitHadrAndLept[SumbTag][ii] < binEdge[iBin+1]){
				genttbarhistoHadrAndLept[iBin]->Fill(CosThGenKinFitHadrAndLept[SumbTag][ii], EventCorrectionWeightKinFitHadrAndLept[SumbTag][ii]) ; 
			      }
			      else if(CosThetaValuesKinFitHadrAndLept[SumbTag][ii] ==1){ //1 is included in last bin
				genttbarhistoHadrAndLept[CosThetaBinNumber-1]->Fill(CosThGenKinFitHadrAndLept[SumbTag][ii], EventCorrectionWeightKinFitHadrAndLept[SumbTag][ii]);
			      }	      	   
			    }
			  }
			  ///////////////////////////////////////////////////////////////////////////////
			  // Change data to systematics since then nominal values will be reweighted!! //
			  ///////////////////////////////////////////////////////////////////////////////

			  //Data result:
			  if(dataSetName.find("Data") == 0)
			    histo1D[CosThetaDataStringHadrAndLept]->Fill(CosThetaValuesKinFitHadrAndLept[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLept[SumbTag][ii]);
			  //JES Plus:
			  if(dataSetName.find("JESPlus") == 0)
			    histo1D[CosThetaJESPlusStringHadrAndLept]->Fill(CosThetaValuesKinFitHadrAndLept[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLept[SumbTag][ii]);
			  //JES Minus:
			  if(dataSetName.find("JESMinus") == 0)
			    histo1D[CosThetaJESMinusStringHadrAndLept]->Fill(CosThetaValuesKinFitHadrAndLept[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLept[SumbTag][ii]);
			  //JER Plus:
			  if(dataSetName.find("JERPlus") == 0)
			    histo1D[CosThetaJERPlusStringHadrAndLept]->Fill(CosThetaValuesKinFitHadrAndLept[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLept[SumbTag][ii]);
			  //JER Minus:
			  if(dataSetName.find("JERMinus") == 0)
			    histo1D[CosThetaJERMinusStringHadrAndLept]->Fill(CosThetaValuesKinFitHadrAndLept[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLept[SumbTag][ii]);
			  //WJets Plus:  --> Fill with WsystPlus sample and all other MC
			  if( dataSetName.find("Data") !=0 && dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSystMinus") != 0  && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatching") != 0&& dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMass") != 0 && dataSetName.find("TriggEvtSel") != 0)
			    histo1D[CosThetaWPlusStringHadrAndLept]->Fill(CosThetaValuesKinFitHadrAndLept[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLept[SumbTag][ii]);
			  //WJets Minus:  --> Fill with WsystMinus sample and all other MC
			  if( dataSetName.find("Data") !=0 && dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSystPlus") != 0  && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatching") != 0&& dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMass") != 0 && dataSetName.find("TriggEvtSel") != 0 )
			    histo1D[CosThetaWMinusStringHadrAndLept]->Fill(CosThetaValuesKinFitHadrAndLept[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLept[SumbTag][ii]);
			  //TTScaling Up:
			  if( dataSetName.find("Data") !=0 && dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 &&dataSetName.find("WSyst") != 0 &&dataSetName.find("TTScalingDown") != 0 &&dataSetName.find("TTMatching") != 0 && dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMass") != 0 && dataSetName.find("TriggEvtSel") != 0 &&dataSetName.find("TTbarJets") !=0)
			    histo1D[CosThetaTTScalingUpStringHadrAndLept]->Fill(CosThetaValuesKinFitHadrAndLept[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLept[SumbTag][ii]);
			  //TTScaling Down:
			  if( dataSetName.find("Data") !=0 &&dataSetName.find("JES") != 0 &&dataSetName.find("JER") != 0 && dataSetName.find("WSyst") != 0 &&dataSetName.find("TTScalingUp") != 0 && dataSetName.find("TTMatching") != 0 && dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMass") != 0 && dataSetName.find("TriggEvtSel") != 0 && dataSetName.find("TTbarJets") !=0 )
			    histo1D[CosThetaTTScalingDownStringHadrAndLept]->Fill(CosThetaValuesKinFitHadrAndLept[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLept[SumbTag][ii]);
			  //TTMatching Up:
			  if(dataSetName.find("Data") !=0 && dataSetName.find("JES") != 0 &&dataSetName.find("JER") != 0 &&dataSetName.find("WSyst") != 0 &&dataSetName.find("TTScaling") != 0 &&dataSetName.find("TTMatchingDown") != 0 && dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMass") != 0 && dataSetName.find("TriggEvtSel") != 0 && dataSetName.find("TTbarJets") !=0 )
			    histo1D[CosThetaTTMatchingUpStringHadrAndLept]->Fill(CosThetaValuesKinFitHadrAndLept[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLept[SumbTag][ii]);
			  //TTMatching Down:
			  if( dataSetName.find("Data") !=0 &&dataSetName.find("JES") != 0 &&dataSetName.find("JER") != 0&& dataSetName.find("WSyst") != 0 && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatchingUp") != 0 && dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMass") != 0 && dataSetName.find("TriggEvtSel") != 0 && dataSetName.find("TTbarJets") !=0 )
			    histo1D[CosThetaTTMatchingDownStringHadrAndLept]->Fill(CosThetaValuesKinFitHadrAndLept[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLept[SumbTag][ii]);		  
			  //PU Plus:
			  if(dataSetName.find("Data") !=0&& dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSyst") != 0 && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatching") != 0 && dataSetName.find("PUMinus") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMass") != 0 && dataSetName.find("TriggEvtSel") != 0 && dataSetName.find("TTbarJets") !=0 )
			    histo1D[CosThetaPUPlusStringHadrAndLept]->Fill(CosThetaValuesKinFitHadrAndLept[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLept[SumbTag][ii]);
			  //PU Minus:
			  if(dataSetName.find("Data") !=0&& dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSyst") != 0 && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatching") != 0 && dataSetName.find("PUPlus") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMass") != 0 && dataSetName.find("TriggEvtSel") != 0 && dataSetName.find("TTbarJets") !=0 )
			    histo1D[CosThetaPUMinusStringHadrAndLept]->Fill(CosThetaValuesKinFitHadrAndLept[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLept[SumbTag][ii]);
			  //UnclusEnergy Plus:
			  if(dataSetName.find("Data") !=0&& dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSyst") != 0 && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatching") != 0 && dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergyMinus") != 0 && dataSetName.find("TopMass") != 0 && dataSetName.find("TriggEvtSel") != 0 && dataSetName.find("TTbarJets") !=0 )
			    histo1D[CosThetaUnclusEnergyPlusStringHadrAndLept]->Fill(CosThetaValuesKinFitHadrAndLept[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLept[SumbTag][ii]);
			  //UnclusEnergy Minus:
			  if(dataSetName.find("Data") !=0&& dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSyst") != 0 && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatching") != 0 && dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergyPlus") != 0 && dataSetName.find("TopMass") != 0  && dataSetName.find("TriggEvtSel") != 0&& dataSetName.find("TTbarJets") !=0 )
			    histo1D[CosThetaUnclusEnergyMinusStringHadrAndLept]->Fill(CosThetaValuesKinFitHadrAndLept[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLept[SumbTag][ii]);
			  //Top Mass Plus:
			  if(dataSetName.find("Data") !=0&& dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSyst") != 0 && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatching") != 0 && dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMassMinus") != 0 && dataSetName.find("TriggEvtSel") != 0 && dataSetName.find("TTbarJets") !=0 )
			    histo1D[CosThetaTopMassPlusStringHadrAndLept]->Fill(CosThetaValuesKinFitHadrAndLept[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLept[SumbTag][ii]);
			  //Top Mass Minus:
			  if(dataSetName.find("Data") !=0&& dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSyst") != 0 && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatching") != 0 && dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMassPlus") != 0 && dataSetName.find("TriggEvtSel") != 0 && dataSetName.find("TTbarJets") !=0 )		      
			    histo1D[CosThetaTopMassMinusStringHadrAndLept]->Fill(CosThetaValuesKinFitHadrAndLept[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLept[SumbTag][ii]);
			  //TriggEvtSel Plus:
			  if(dataSetName.find("Data") !=0&& dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSyst") != 0 && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatching") != 0 && dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMass") != 0 && dataSetName.find("TriggEvtSelMinus") != 0 && dataSetName.find("TTbarJets") !=0 )
			    histo1D[CosThetaTriggEvtSelPlusStringHadrAndLept]->Fill(CosThetaValuesKinFitHadrAndLept[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLept[SumbTag][ii]);
			  //TriggEvtSel Minus:
			  if(dataSetName.find("Data") !=0&& dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSyst") != 0 && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatching") != 0 && dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMass") != 0 && dataSetName.find("TriggEvtSelPlus") != 0 && dataSetName.find("TTbarJets") !=0 )		      
			    histo1D[CosThetaTriggEvtSelMinusStringHadrAndLept]->Fill(CosThetaValuesKinFitHadrAndLept[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLept[SumbTag][ii]);

			  if(SignalOnly == true && ((dataSetName.find("TTbarJets_SemiMu") == 0 && semiMuon==true) || (dataSetName.find("TTbarJets_SemiEl")==0 && semiElectron == true)))
			    histo1D[CosThetaDataStringHadrAndLept]->Fill(CosThetaValuesKinFitHadrAndLept[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLept[SumbTag][ii]);

			  if(((dataSetName.find("TTbarJets_SemiMu")==0 && semiMuon==true) || (dataSetName.find("TTbarJets_SemiEl")==0 && semiElectron == true) ) && dataSetName.find("JES_") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSyst") != 0 && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatching") != 0 && dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMass") != 0 && dataSetName.find("TriggEvtSel") != 0)
			    histo1D[CosThetaSignalStringHadrAndLept]->Fill(CosThetaValuesKinFitHadrAndLept[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLept[SumbTag][ii]);

			  if (dataSetName.find("Data") !=0 && dataSetName.find("JES") != 0 && dataSetName.find("JER") != 0 && dataSetName.find("WSyst") != 0 && dataSetName.find("TTScaling") != 0 && dataSetName.find("TTMatching") != 0 && dataSetName.find("PU") != 0 && dataSetName.find("UnclusEnergy") != 0 && dataSetName.find("TopMass") != 0 && dataSetName.find("TriggEvtSel") != 0 && ((dataSetName.find("TTbarJets_SemiMu")!=0 && semiMuon==true) || (dataSetName.find("TTbarJets_SemiEl")!=0 && semiElectron == true) ))
			    histo1D[CosThetaBckgStringHadrAndLept]->Fill(CosThetaValuesKinFitHadrAndLept[SumbTag][ii],EventCorrectionWeightKinFitHadrAndLept[SumbTag][ii]);		  
			}
		    
			//Perform MinuitFitter:
			if(iDataSet==(datasets.size()-1)){//Go in this loop when the last datasample is active to perform MinuitFitter
			  //Data:
			  MinuitFitter minuitFitterHadrAndLept = MinuitFitter(histo1D[CosThetaDataStringHadrAndLept], histo1D[CosThetaSignalStringHadrAndLept], histo1D[CosThetaBckgStringHadrAndLept], 0.6671, 0.3325, 0.0004,genttbarhistoHadrAndLept,ndimen);
			  //JESPlus
			  MinuitFitter minuitFitterHadrAndLeptJESPlus = MinuitFitter(histo1D[CosThetaJESPlusStringHadrAndLept], histo1D[CosThetaSignalStringHadrAndLept], histo1D[CosThetaBckgStringHadrAndLept], 0.6671, 0.3325, 0.0004,genttbarhistoHadrAndLept,ndimen);			  
			  //JES Minus
			  MinuitFitter minuitFitterHadrAndLeptJESMinus = MinuitFitter(histo1D[CosThetaJESMinusStringHadrAndLept], histo1D[CosThetaSignalStringHadrAndLept], histo1D[CosThetaBckgStringHadrAndLept], 0.6671, 0.3325, 0.0004,genttbarhistoHadrAndLept,ndimen);			  
			  //JERPlus
			  MinuitFitter minuitFitterHadrAndLeptJERPlus = MinuitFitter(histo1D[CosThetaJERPlusStringHadrAndLept], histo1D[CosThetaSignalStringHadrAndLept], histo1D[CosThetaBckgStringHadrAndLept], 0.6671, 0.3325, 0.0004,genttbarhistoHadrAndLept,ndimen);			  
			  //JER Minus
			  MinuitFitter minuitFitterHadrAndLeptJERMinus = MinuitFitter(histo1D[CosThetaJERMinusStringHadrAndLept], histo1D[CosThetaSignalStringHadrAndLept], histo1D[CosThetaBckgStringHadrAndLept], 0.6671, 0.3325, 0.0004,genttbarhistoHadrAndLept,ndimen);			  
			  //W Plus
			  MinuitFitter minuitFitterHadrAndLeptWPlus = MinuitFitter(histo1D[CosThetaWPlusStringHadrAndLept], histo1D[CosThetaSignalStringHadrAndLept], histo1D[CosThetaBckgStringHadrAndLept], 0.6671, 0.3325, 0.0004,genttbarhistoHadrAndLept,ndimen);			  
			  //W Minus
			  MinuitFitter minuitFitterHadrAndLeptWMinus = MinuitFitter(histo1D[CosThetaWMinusStringHadrAndLept], histo1D[CosThetaSignalStringHadrAndLept], histo1D[CosThetaBckgStringHadrAndLept], 0.6671, 0.3325, 0.0004,genttbarhistoHadrAndLept,ndimen);
			  //TTScaling Up
			  MinuitFitter minuitFitterHadrAndLeptTTScalingUp = MinuitFitter(histo1D[CosThetaTTScalingUpStringHadrAndLept], histo1D[CosThetaSignalStringHadrAndLept], histo1D[CosThetaBckgStringHadrAndLept], 0.6671, 0.3325, 0.0004,genttbarhistoHadrAndLept,ndimen);
			  //TTScaling Down
			  MinuitFitter minuitFitterHadrAndLeptTTScalingDown = MinuitFitter(histo1D[CosThetaTTScalingDownStringHadrAndLept], histo1D[CosThetaSignalStringHadrAndLept], histo1D[CosThetaBckgStringHadrAndLept], 0.6671, 0.3325, 0.0004,genttbarhistoHadrAndLept,ndimen);
			  //TTMatching Up
			  MinuitFitter minuitFitterHadrAndLeptTTMatchingUp = MinuitFitter(histo1D[CosThetaTTMatchingUpStringHadrAndLept], histo1D[CosThetaSignalStringHadrAndLept], histo1D[CosThetaBckgStringHadrAndLept], 0.6671, 0.3325, 0.0004,genttbarhistoHadrAndLept,ndimen);
			  //TTMatching Down
			  MinuitFitter minuitFitterHadrAndLeptTTMatchingDown = MinuitFitter(histo1D[CosThetaTTMatchingDownStringHadrAndLept], histo1D[CosThetaSignalStringHadrAndLept], histo1D[CosThetaBckgStringHadrAndLept], 0.6671, 0.3325, 0.0004,genttbarhistoHadrAndLept,ndimen);
			  //PU Plus
			  MinuitFitter minuitFitterHadrAndLeptPUPlus = MinuitFitter(histo1D[CosThetaPUPlusStringHadrAndLept], histo1D[CosThetaSignalStringHadrAndLept], histo1D[CosThetaBckgStringHadrAndLept], 0.6671, 0.3325, 0.0004,genttbarhistoHadrAndLept,ndimen);			  
			  //PU Minus
			  MinuitFitter minuitFitterHadrAndLeptPUMinus = MinuitFitter(histo1D[CosThetaPUMinusStringHadrAndLept], histo1D[CosThetaSignalStringHadrAndLept], histo1D[CosThetaBckgStringHadrAndLept], 0.6671, 0.3325, 0.0004,genttbarhistoHadrAndLept,ndimen);			  
			  //UnclusEnergy Plus
			  MinuitFitter minuitFitterHadrAndLeptUnclusEnergyPlus = MinuitFitter(histo1D[CosThetaUnclusEnergyPlusStringHadrAndLept], histo1D[CosThetaSignalStringHadrAndLept], histo1D[CosThetaBckgStringHadrAndLept], 0.6671, 0.3325, 0.0004,genttbarhistoHadrAndLept,ndimen);			  
			  //UnclusEnergy Minus
			  MinuitFitter minuitFitterHadrAndLeptUnclusEnergyMinus = MinuitFitter(histo1D[CosThetaUnclusEnergyMinusStringHadrAndLept], histo1D[CosThetaSignalStringHadrAndLept], histo1D[CosThetaBckgStringHadrAndLept], 0.6671, 0.3325, 0.0004,genttbarhistoHadrAndLept,ndimen);			  
			  //TopMass Plus
			  MinuitFitter minuitFitterHadrAndLeptTopMassPlus = MinuitFitter(histo1D[CosThetaTopMassPlusStringHadrAndLept], histo1D[CosThetaSignalStringHadrAndLept], histo1D[CosThetaBckgStringHadrAndLept], 0.6671, 0.3325, 0.0004,genttbarhistoHadrAndLept,ndimen);			  
			  //TopMass Minus
			  MinuitFitter minuitFitterHadrAndLeptTopMassMinus = MinuitFitter(histo1D[CosThetaTopMassMinusStringHadrAndLept], histo1D[CosThetaSignalStringHadrAndLept], histo1D[CosThetaBckgStringHadrAndLept], 0.6671, 0.3325, 0.0004,genttbarhistoHadrAndLept,ndimen);			  
			  //TriggEvtSel Plus
			  MinuitFitter minuitFitterHadrAndLeptTriggEvtSelPlus = MinuitFitter(histo1D[CosThetaTriggEvtSelPlusStringHadrAndLept], histo1D[CosThetaSignalStringHadrAndLept], histo1D[CosThetaBckgStringHadrAndLept], 0.6671, 0.3325, 0.0004,genttbarhistoHadrAndLept,ndimen);			  
			  //TriggEvtSel Minus
			  MinuitFitter minuitFitterHadrAndLeptTriggEvtSelMinus = MinuitFitter(histo1D[CosThetaTriggEvtSelMinusStringHadrAndLept], histo1D[CosThetaSignalStringHadrAndLept], histo1D[CosThetaBckgStringHadrAndLept], 0.6671, 0.3325, 0.0004,genttbarhistoHadrAndLept,ndimen);			  
			  
			  PresentationTexHadrAndLept << PresentationOutput[SumbTag] << " & ";
			  cout << " Inside minuitFitter loop : " << PresentationOutput[SumbTag] << " & " << endl;
			  for(int ii=0; ii< nameDataSet.size(); ii++){
			    if(ii < nameDataSet.size()-1){ PresentationTexHadrAndLept << NumberRemainingEventsKinFitHadrAndLept[SumbTag][ii] << " & ";}
			    else if(ii == nameDataSet.size()-1 ){ 
			      PresentationTexHadrAndLept << NumberRemainingEventsKinFitHadrAndLept[SumbTag][ii] << " & ";//For presentation helicity values still need to be included !!
			      PresentationTexHadrAndLept << NumberBLeptCorrectEventsKinFitHadrAndLept[SumbTag][ii] << " & ";
			    }
			  }	
		  
			  //FL results:
			  float FLJESSystHadrAndLept = (abs(minuitFitterHadrAndLept.GetFLResult() - minuitFitterHadrAndLeptJESPlus.GetFLResult()) + abs(minuitFitterHadrAndLept.GetFLResult() - minuitFitterHadrAndLeptJESMinus.GetFLResult()))/2;
			  float FLJERSystHadrAndLept = 0;
			  if(JERResults == true) FLJERSystHadrAndLept = (abs(minuitFitterHadrAndLept.GetFLResult() - minuitFitterHadrAndLeptJERPlus.GetFLResult()) + abs(minuitFitterHadrAndLept.GetFLResult() - minuitFitterHadrAndLeptJERMinus.GetFLResult()))/2;
			  float FLWSystHadrAndLept = 0;
			  if(WSystResults == true) FLWSystHadrAndLept = (abs(minuitFitterHadrAndLept.GetFLResult() - minuitFitterHadrAndLeptWPlus.GetFLResult()) + abs(minuitFitterHadrAndLept.GetFLResult() - minuitFitterHadrAndLeptWMinus.GetFLResult()))/2;
			  float FLTTScalingHadrAndLept = 0;
			  if(TTScalingResults == true) FLTTScalingHadrAndLept = (abs(minuitFitterHadrAndLept.GetFLResult() - minuitFitterHadrAndLeptTTScalingUp.GetFLResult()) + abs(minuitFitterHadrAndLept.GetFLResult() - minuitFitterHadrAndLeptTTScalingDown.GetFLResult()))/2;
			  float FLTTMatchingHadrAndLept = 0;
			  if(TTMatchingResults == true) FLTTMatchingHadrAndLept = (abs(minuitFitterHadrAndLept.GetFLResult() - minuitFitterHadrAndLeptTTMatchingUp.GetFLResult()) + abs(minuitFitterHadrAndLept.GetFLResult() - minuitFitterHadrAndLeptTTMatchingDown.GetFLResult()))/2;
			  float FLPUSystHadrAndLept = 0;
			  if(PUResults == true) FLPUSystHadrAndLept = (abs(minuitFitterHadrAndLept.GetFLResult() - minuitFitterHadrAndLeptPUPlus.GetFLResult()) + abs(minuitFitterHadrAndLept.GetFLResult() - minuitFitterHadrAndLeptPUMinus.GetFLResult()))/2;
			  float FLUnclusEnergySystHadrAndLept = 0;
			  if(UnclusEnergyResults == true) FLUnclusEnergySystHadrAndLept = (abs(minuitFitterHadrAndLept.GetFLResult() - minuitFitterHadrAndLeptUnclusEnergyPlus.GetFLResult()) + abs(minuitFitterHadrAndLept.GetFLResult() - minuitFitterHadrAndLeptUnclusEnergyMinus.GetFLResult()))/2;
			  float FLTopMassSystHadrAndLept = 0;
			  if(TopMassResults == true) FLTopMassSystHadrAndLept = (abs(minuitFitterHadrAndLept.GetFLResult() - minuitFitterHadrAndLeptTopMassPlus.GetFLResult()) + abs(minuitFitterHadrAndLept.GetFLResult() - minuitFitterHadrAndLeptTopMassMinus.GetFLResult()))/2;		      
			  float FLTriggEvtSelSystHadrAndLept = 0;
			  if(TriggEvtSelResults == true) FLTriggEvtSelSystHadrAndLept = (abs(minuitFitterHadrAndLept.GetFLResult() - minuitFitterHadrAndLeptTriggEvtSelPlus.GetFLResult()) + abs(minuitFitterHadrAndLept.GetFLResult() - minuitFitterHadrAndLeptTriggEvtSelMinus.GetFLResult()))/2;		      
			  FMinusTexHadrAndLept << "\\hline" << endl;
			  FMinusTexHadrAndLept << PresentationOutput[SumbTag] << " & " << minuitFitterHadrAndLept.GetFLResult() << " & " << minuitFitterHadrAndLept.GetFLResult()-SMfrResult <<" & " << minuitFitterHadrAndLept.GetFLError() << " & " << FLJESSystHadrAndLept << " & " << FLJERSystHadrAndLept << " & " << FLWSystHadrAndLept << " & " << FLTTScalingHadrAndLept << " & " << FLTTMatchingHadrAndLept << " & " << FLPUSystHadrAndLept << " & " << FLUnclusEnergySystHadrAndLept << " & " << FLTopMassSystHadrAndLept << " & " << FLTriggEvtSelSystHadrAndLept << " & " << sqrt(FLJESSystHadrAndLept*FLJESSystHadrAndLept + FLJERSystHadrAndLept*FLJERSystHadrAndLept + FLWSystHadrAndLept*FLWSystHadrAndLept + FLTTScalingHadrAndLept*FLTTScalingHadrAndLept + FLTTMatchingHadrAndLept*FLTTMatchingHadrAndLept + FLPUSystHadrAndLept*FLPUSystHadrAndLept + FLUnclusEnergySystHadrAndLept*FLUnclusEnergySystHadrAndLept + FLTopMassSystHadrAndLept*FLTopMassSystHadrAndLept + FLTriggEvtSelSystHadrAndLept*FLTriggEvtSelSystHadrAndLept) << " & " << sqrt(FLJESSystHadrAndLept*FLJESSystHadrAndLept + FLJERSystHadrAndLept*FLJERSystHadrAndLept + FLWSystHadrAndLept*FLWSystHadrAndLept + FLTTScalingHadrAndLept*FLTTScalingHadrAndLept + FLTTMatchingHadrAndLept*FLTTMatchingHadrAndLept + FLPUSystHadrAndLept*FLPUSystHadrAndLept + FLUnclusEnergySystHadrAndLept*FLUnclusEnergySystHadrAndLept + FLTopMassSystHadrAndLept*FLTopMassSystHadrAndLept + FLTriggEvtSelSystHadrAndLept*FLTriggEvtSelSystHadrAndLept + (minuitFitterHadrAndLept.GetFLError())*(minuitFitterHadrAndLept.GetFLError())) << "\\\\" << endl;
		  
			  //FR results:
			  float FRJESSystHadrAndLept = (abs(minuitFitterHadrAndLept.GetFRResult() - minuitFitterHadrAndLeptJESPlus.GetFRResult()) + abs(minuitFitterHadrAndLept.GetFRResult() - minuitFitterHadrAndLeptJESMinus.GetFRResult()))/2;
			  float FRJERSystHadrAndLept = 0;
			  if(JERResults == true) FRJERSystHadrAndLept = (abs(minuitFitterHadrAndLept.GetFRResult() - minuitFitterHadrAndLeptJERPlus.GetFRResult()) + abs(minuitFitterHadrAndLept.GetFRResult() - minuitFitterHadrAndLeptJERMinus.GetFRResult()))/2;
			  float FRWSystHadrAndLept = 0;
			  if(WSystResults == true) FRWSystHadrAndLept = (abs(minuitFitterHadrAndLept.GetFRResult() - minuitFitterHadrAndLeptWPlus.GetFRResult()) + abs(minuitFitterHadrAndLept.GetFRResult() - minuitFitterHadrAndLeptWMinus.GetFRResult()))/2;
			  float FRTTScalingHadrAndLept = 0;
			  if(TTScalingResults == true) FRTTScalingHadrAndLept = (abs(minuitFitterHadrAndLept.GetFRResult() - minuitFitterHadrAndLeptTTScalingUp.GetFRResult()) + abs(minuitFitterHadrAndLept.GetFRResult() - minuitFitterHadrAndLeptTTScalingDown.GetFRResult()))/2;
			  float FRTTMatchingHadrAndLept = 0;
			  if(TTMatchingResults == true) FRTTMatchingHadrAndLept = (abs(minuitFitterHadrAndLept.GetFRResult() - minuitFitterHadrAndLeptTTMatchingUp.GetFRResult()) + abs(minuitFitterHadrAndLept.GetFRResult() - minuitFitterHadrAndLeptTTMatchingDown.GetFRResult()))/2;		  
			  float FRPUSystHadrAndLept = 0;
			  if(PUResults == true) FRPUSystHadrAndLept = (abs(minuitFitterHadrAndLept.GetFRResult() - minuitFitterHadrAndLeptPUPlus.GetFRResult()) + abs(minuitFitterHadrAndLept.GetFRResult() - minuitFitterHadrAndLeptPUMinus.GetFRResult()))/2;
			  float FRUnclusEnergySystHadrAndLept = 0;
			  if(UnclusEnergyResults == true) FRUnclusEnergySystHadrAndLept = (abs(minuitFitterHadrAndLept.GetFRResult() - minuitFitterHadrAndLeptUnclusEnergyPlus.GetFRResult()) + abs(minuitFitterHadrAndLept.GetFRResult() - minuitFitterHadrAndLeptUnclusEnergyMinus.GetFRResult()))/2;
			  float FRTopMassSystHadrAndLept = 0;
			  if(TopMassResults == true) FRTopMassSystHadrAndLept = (abs(minuitFitterHadrAndLept.GetFRResult() - minuitFitterHadrAndLeptTopMassPlus.GetFRResult()) + abs(minuitFitterHadrAndLept.GetFRResult() - minuitFitterHadrAndLeptTopMassMinus.GetFRResult()))/2;
			  float FRTriggEvtSelSystHadrAndLept = 0;
			  if(TriggEvtSelResults == true) FRTriggEvtSelSystHadrAndLept = (abs(minuitFitterHadrAndLept.GetFRResult() - minuitFitterHadrAndLeptTriggEvtSelPlus.GetFRResult()) + abs(minuitFitterHadrAndLept.GetFRResult() - minuitFitterHadrAndLeptTriggEvtSelMinus.GetFRResult()))/2;
			  FPlusTexHadrAndLept << "\\hline" << endl;
			  FPlusTexHadrAndLept << PresentationOutput[SumbTag] << " & " << minuitFitterHadrAndLept.GetFRResult() << " & " << minuitFitterHadrAndLept.GetFRResult()-SMfrResult <<" & " << minuitFitterHadrAndLept.GetFRError() << " & " << FRJESSystHadrAndLept << " & " << FRJERSystHadrAndLept << " & " << FRWSystHadrAndLept << " & " << FRTTScalingHadrAndLept << " & " << FRTTMatchingHadrAndLept << " & " << FRPUSystHadrAndLept << " & " << FRUnclusEnergySystHadrAndLept << " & " << FRTopMassSystHadrAndLept << " & " << FRTriggEvtSelSystHadrAndLept << " & " << sqrt(FRJESSystHadrAndLept*FRJESSystHadrAndLept + FRJERSystHadrAndLept*FRJERSystHadrAndLept + FRWSystHadrAndLept*FRWSystHadrAndLept + FRTTScalingHadrAndLept*FRTTScalingHadrAndLept + FRTTMatchingHadrAndLept*FRTTMatchingHadrAndLept + FRPUSystHadrAndLept*FRPUSystHadrAndLept + FRUnclusEnergySystHadrAndLept*FRUnclusEnergySystHadrAndLept + FRTopMassSystHadrAndLept*FRTopMassSystHadrAndLept + FRTriggEvtSelSystHadrAndLept*FRTriggEvtSelSystHadrAndLept) << " & " << sqrt(FRJESSystHadrAndLept*FRJESSystHadrAndLept + FRJERSystHadrAndLept*FRJERSystHadrAndLept + FRWSystHadrAndLept*FRWSystHadrAndLept + FRTTScalingHadrAndLept*FRTTScalingHadrAndLept + FRTTMatchingHadrAndLept*FRTTMatchingHadrAndLept + FRPUSystHadrAndLept*FRPUSystHadrAndLept + FRUnclusEnergySystHadrAndLept*FRUnclusEnergySystHadrAndLept + FRTopMassSystHadrAndLept*FRTopMassSystHadrAndLept + FRTriggEvtSelSystHadrAndLept*FRTriggEvtSelSystHadrAndLept + (minuitFitterHadrAndLept.GetFRError())*(minuitFitterHadrAndLept.GetFRError())) << "\\\\" << endl;

			  //F0 results:
			  float F0JESSystHadrAndLept = (abs(minuitFitterHadrAndLept.GetF0Result() - minuitFitterHadrAndLeptJESPlus.GetF0Result()) + abs(minuitFitterHadrAndLept.GetF0Result() - minuitFitterHadrAndLeptJESMinus.GetF0Result()))/2;
			  float F0JERSystHadrAndLept = 0;
			  if(JERResults == true) F0JERSystHadrAndLept = (abs(minuitFitterHadrAndLept.GetF0Result() - minuitFitterHadrAndLeptJERPlus.GetF0Result()) + abs(minuitFitterHadrAndLept.GetF0Result() - minuitFitterHadrAndLeptJERMinus.GetF0Result()))/2;
			  float F0WSystHadrAndLept = 0;
			  if(WSystResults == true) F0WSystHadrAndLept = (abs(minuitFitterHadrAndLept.GetF0Result() - minuitFitterHadrAndLeptWPlus.GetF0Result()) + abs(minuitFitterHadrAndLept.GetF0Result() - minuitFitterHadrAndLeptWMinus.GetF0Result()))/2;
			  float F0TTScalingHadrAndLept = 0;
			  if(TTScalingResults == true) F0TTScalingHadrAndLept = (abs(minuitFitterHadrAndLept.GetF0Result() - minuitFitterHadrAndLeptTTScalingUp.GetF0Result()) + abs(minuitFitterHadrAndLept.GetF0Result() - minuitFitterHadrAndLeptTTScalingDown.GetF0Result()))/2;
			  float F0TTMatchingHadrAndLept = 0;
			  if(TTMatchingResults == true) F0TTMatchingHadrAndLept = (abs(minuitFitterHadrAndLept.GetF0Result() - minuitFitterHadrAndLeptTTMatchingUp.GetF0Result()) + abs(minuitFitterHadrAndLept.GetF0Result() - minuitFitterHadrAndLeptTTMatchingDown.GetF0Result()))/2;		  
			  float F0PUSystHadrAndLept = 0;
			  if(PUResults == true) F0PUSystHadrAndLept = (abs(minuitFitterHadrAndLept.GetF0Result() - minuitFitterHadrAndLeptPUPlus.GetF0Result()) + abs(minuitFitterHadrAndLept.GetF0Result() - minuitFitterHadrAndLeptPUMinus.GetF0Result()))/2;
			  float F0UnclusEnergySystHadrAndLept = 0;
			  if(UnclusEnergyResults == true) F0UnclusEnergySystHadrAndLept = (abs(minuitFitterHadrAndLept.GetF0Result() - minuitFitterHadrAndLeptUnclusEnergyPlus.GetF0Result()) + abs(minuitFitterHadrAndLept.GetF0Result() - minuitFitterHadrAndLeptUnclusEnergyMinus.GetF0Result()))/2;
			  float F0TopMassSystHadrAndLept = 0;
			  if(TopMassResults == true) F0TopMassSystHadrAndLept = (abs(minuitFitterHadrAndLept.GetF0Result() - minuitFitterHadrAndLeptTopMassPlus.GetF0Result()) + abs(minuitFitterHadrAndLept.GetF0Result() - minuitFitterHadrAndLeptTopMassMinus.GetF0Result()))/2;
			  float F0TriggEvtSelSystHadrAndLept = 0;
			  if(TriggEvtSelResults == true) F0TriggEvtSelSystHadrAndLept = (abs(minuitFitterHadrAndLept.GetF0Result() - minuitFitterHadrAndLeptTriggEvtSelPlus.GetF0Result()) + abs(minuitFitterHadrAndLept.GetF0Result() - minuitFitterHadrAndLeptTriggEvtSelMinus.GetF0Result()))/2;
			  FZeroTexHadrAndLept << " \\hline" << endl;
			  FZeroTexHadrAndLept << PresentationOutput[SumbTag] << " & " << minuitFitterHadrAndLept.GetF0Result() << " & " << minuitFitterHadrAndLept.GetF0Result()-SMfrResult <<" & " << minuitFitterHadrAndLept.GetF0Error() << " & " << F0JESSystHadrAndLept << " & " << F0JERSystHadrAndLept << " & " << F0WSystHadrAndLept << " & " << F0TTScalingHadrAndLept << " & " << F0TTMatchingHadrAndLept << " & " << F0PUSystHadrAndLept << " & " << F0UnclusEnergySystHadrAndLept << " & " << F0TopMassSystHadrAndLept << " & " << F0TriggEvtSelSystHadrAndLept << " & " << sqrt(F0JESSystHadrAndLept*F0JESSystHadrAndLept + F0JERSystHadrAndLept*F0JERSystHadrAndLept + F0WSystHadrAndLept*F0WSystHadrAndLept + F0TTScalingHadrAndLept*F0TTScalingHadrAndLept + F0TTMatchingHadrAndLept*F0TTMatchingHadrAndLept + F0PUSystHadrAndLept*F0PUSystHadrAndLept + F0UnclusEnergySystHadrAndLept*F0UnclusEnergySystHadrAndLept + F0TopMassSystHadrAndLept*F0TopMassSystHadrAndLept + F0TriggEvtSelSystHadrAndLept*F0TriggEvtSelSystHadrAndLept) << " & " << sqrt(F0JESSystHadrAndLept*F0JESSystHadrAndLept + F0JERSystHadrAndLept*F0JERSystHadrAndLept + F0WSystHadrAndLept*F0WSystHadrAndLept + F0TTScalingHadrAndLept*F0TTScalingHadrAndLept + F0TTMatchingHadrAndLept*F0TTMatchingHadrAndLept + F0PUSystHadrAndLept*F0PUSystHadrAndLept + F0UnclusEnergySystHadrAndLept*F0UnclusEnergySystHadrAndLept + F0TopMassSystHadrAndLept*F0TopMassSystHadrAndLept + F0TriggEvtSelSystHadrAndLept*F0TriggEvtSelSystHadrAndLept + (minuitFitterHadrAndLept.GetF0Error())*(minuitFitterHadrAndLept.GetF0Error())) << "\\\\" << endl;

			  //Statistical results and number of events:
			  PresentationTexHadrAndLept << minuitFitterHadrAndLept.GetFRResult() << " & " << minuitFitterHadrAndLept.GetFRResult()-SMfrResult <<" & " << minuitFitterHadrAndLept.GetFRError() << " & " << minuitFitterHadrAndLept.GetFLResult() << " & " << minuitFitterHadrAndLept.GetFLResult()-SMflResult << " & " << minuitFitterHadrAndLept.GetFLError() << " & " << minuitFitterHadrAndLept.GetF0Result() << " & " << minuitFitterHadrAndLept.GetF0Result()-SMf0Result << " & " << minuitFitterHadrAndLept.GetF0Error() << " \\\\ " << endl;		
			  PresentationTexHadrAndLept << " \\hline " << endl;
			}
		      }
		  		
		    }//end of wrong entry chosen failed requirement!
		    else{
		      cout << " Looking at wrong bTagging combination in fitting minuit !! " << endl;
		      cout << " Looking at : TCHE = " << TCHE+1<< " TCHP = " << TCHP+1 << " SSVHE = " << SSVHE+1 << " SSVHP = " << SSVHP+1 << " CSV = " << CSV+1 << " JP = " << JP+1 << " JBP = " << JBP+1<< endl;
		      cout << " Correct combination : TCHE = "<< UsedTCHE[SumbTag] << " TCHP = " << UsedTCHP[SumbTag] << " SSVHE = " << UsedSSVHE[SumbTag] << " SSVHP = " << UsedSSVHP[SumbTag] << " CSV = " << UsedCSV[SumbTag] << " JP = " << UsedJP[SumbTag] << " JBP = " << UsedJBP[SumbTag] << endl;
		    }
	      
		    if(TCHE==NumberTCHEbTags-1) {
		      ConsideredTagger =1;
		      TCHP=1;
		      SSVHE=1;
		      SSVHP=1;
		      CSV=1;
		      JP=1;
		      JBP=1;
		    }	      
		  }//end of TCHE
		  if(TCHP == NumberTCHPbTags-1) ConsideredTagger =2;
		}//end of TCHP
		if(SSVHE == NumberSSVHEbTags-1) ConsideredTagger = 3;
	      }//end of SSVHE
	      if(SSVHP == NumberSSVHPbTags-1) ConsideredTagger = 4;
	    }//end of SSVHP
	    if(CSV == NumberCSVbTags-1) ConsideredTagger = 5;
	  }//end of CSV
	  if(JP == NumberJPbTags-1) ConsideredTagger=6;
	}//end of JP
      }//end of JBP      

    } //End of boolean request: PerformMinuit
        
    outFile << endl << "Finished event-by-event info, end of file!" << endl;
    outFile.close();
    
    inFile->Close();
    delete inFile;
    
  } // end loop over datasets

  //Selection tables:
  selecTableMacro.TableCalculator(false,true,true,true);
  string selectiontableMacro = "SelectionTable_Analysis_Macro";
  if(argc >= 3){
    string sample = string(argv[2]);
    selectiontableMacro = selectiontableMacro+"_"+sample;
  }
  selectiontableMacro = selectiontableMacro+".tex";
  selecTableMacro.Write(selectiontableMacro.c_str(),1);

  //Write output of number of events:
  cout << " " << endl;
  if(semiMuon == true) cout << " TTbar SemiMu " << endl;
  else if(semiElectron == true) cout << " TTbar SemiEl " << endl;
  std::cout << " Comparing efficiency of selected events " << endl;
  cout << " NumberEventsBeforeCuts : " << NumberEventsBeforeCuts << endl;
  cout << " NumberSelectedEvents : " <<  NumberSelectedEvents << endl;
  cout << "                -------        " << endl;
  cout << " NumberDataEventsBeforeCuts : " << NumberDataEventsBeforeCuts << endl;
  cout << " NumberSelectedDataEvents : " << NumberSelectedDataEvents  << endl;
  std::cout << " " << endl;

  //Write output for KinFit probability cut:
  cout << " " << endl;
  cout << " KinFit probability information: " << endl;
  if(semiMuon == true) cout << " --------- TTbar SemiMu " << endl;
  else if(semiElectron == true) cout << " -------------- TTbar SemiEl " << endl;
  cout << " Number of signal events with cos theta found: " << NumberEventsSignalCosThetaFound << endl;
  cout << " Number of signal events with low probability (0.01) and cos theta found : " << NumberEventsLowProbCosThetaFound << endl;
  cout << " Number of signal events with correct leptonic b (cos theta found) : " << NumberEventsSignalCorrectLeptBCosThetaFound << endl;
  cout << " Number of signal events with correct hadronic b (cos theta found) : " << NumberEventsSignalCorrectHadrBCosThetaFound << endl;
  cout << " Number of signal events with correct light jets (cos theta found) : " << NumberEventsSignalCorrectLightCosThetaFound << endl;
  cout << " Number of signal events with correct hadronic topology (cos theta found) : " << NumberEventsSignalCorrectHadronicCosThetaFound << endl;
  cout << " Number of signal events with low probability (0.01) with correct leptonic b (cos theta found) : " << NumberEventsLowProbCorrectLeptBCosThetaFound << endl;
  cout << " Number of signal events with low probability (0.01) with correct hadronic b (cos theta found) : " << NumberEventsLowProbCorrectHadrBCosThetaFound << endl;
  cout << " Number of signal events with low probability (0.01) with correct light jets (cos theta found) : " << NumberEventsLowProbCorrectLightCosThetaFound << endl;
  cout << " Number of signal events with low probability (0.01) with correct hadronic topology (cos theta found) : " << NumberEventsLowProbCorrectHadronicCosThetaFound << endl;
  cout << " " << endl;

  cout << " " << endl;
  cout << "  Number comparison between ChiSqOnly and ChiSqAndBTag method for jet combination selection !! " << endl;
  cout << "  1) ChiSq Only: " << endl;
  cout << "     - Starting events : " << NumberEventsChiSqOnly[0];
  for(int ii = 1; ii<inputWTree.size(); ii++) cout << "  , " << NumberEventsChiSqOnly[ii];
  cout << " " << endl;
  cout << "    - Jet comb found : " << NumberEventsChiSqOnlyJetCombFound[0];
  for(int ii = 1; ii<inputWTree.size(); ii++) cout << "  , " << NumberEventsChiSqOnlyJetCombFound[ii];
  cout << " " << endl;
  cout << "    - Cos Theta found : " << NumberEventsChiSqOnlyCosThetaFound[0];
  for(int ii = 1; ii<inputWTree.size(); ii++) cout << "  , " << NumberEventsChiSqOnlyCosThetaFound[ii];
  cout << " " << endl;
  cout << "    - Specific bTag : " << NumberEventsChiSqOnlySpecificBTag[0];
  for(int ii = 1; ii<inputWTree.size(); ii++) cout << "  , " << NumberEventsChiSqOnlySpecificBTag[ii];
  cout << " " << endl;
  cout << "  2) ChiSq and BTag: " << endl;
  cout << "     - Starting events : " << NumberEventsChiSqAndBTag[0];
  for(int ii = 1; ii<inputWTree.size(); ii++) cout << "  , " << NumberEventsChiSqAndBTag[ii];
  cout << " " << endl;
  cout << "    - Jet comb found : " << NumberEventsChiSqAndBTagJetCombFound[0];
  for(int ii = 1; ii<inputWTree.size(); ii++) cout << "  , " << NumberEventsChiSqAndBTagJetCombFound[ii];
  cout << " " << endl;
  cout << "    - Cos Theta found : " << NumberEventsChiSqAndBTagCosThetaFound[0];
  for(int ii = 1; ii<inputWTree.size(); ii++) cout << "  , " << NumberEventsChiSqAndBTagCosThetaFound[ii];
  cout << " " << endl;
  cout << "    - Specific bTag : " << NumberEventsChiSqAndBTagSpecificBTag[0];
  for(int ii = 1; ii<inputWTree.size(); ii++) cout << "  , " << NumberEventsChiSqAndBTagSpecificBTag[ii];
  cout << " " << endl;
 
  //Hadronic:
  PresentationTexHadr << " \\end{tabular} " << endl;
  PresentationTexHadr << " \\end{center} " << endl;
  PresentationTexHadr << " \\renewcommand{\\arraystretch}{0.9} " << endl;
  PresentationTexHadr << " \\end{tiny} " << endl;
  PresentationTexHadr << " \\end{table} " << endl;
  PresentationTexHadr.close();

  FMinusTexHadr.close();
  FPlusTexHadr.close();
  FZeroTexHadr.close();

  //Hadronic and leptonic W only:
  PresentationTexHadrAndLeptWOnly << " \\end{tabular} " << endl;
  PresentationTexHadrAndLeptWOnly << " \\end{center} " << endl;
  PresentationTexHadrAndLeptWOnly << " \\renewcommand{\\arraystretch}{0.9} " << endl;
  PresentationTexHadrAndLeptWOnly << " \\end{tiny} " << endl;
  PresentationTexHadrAndLeptWOnly << " \\end{table} " << endl;
  PresentationTexHadrAndLeptWOnly.close();

  FMinusTexHadrAndLeptWOnly.close();
  FPlusTexHadrAndLeptWOnly.close();
  FZeroTexHadrAndLeptWOnly.close();
  
  //Hadronic and leptonic:
  PresentationTexHadrAndLept << " \\end{tabular} " << endl;
  PresentationTexHadrAndLept << " \\end{center} " << endl;
  PresentationTexHadrAndLept << " \\renewcommand{\\arraystretch}{0.9} " << endl;
  PresentationTexHadrAndLept << " \\end{tiny} " << endl;
  PresentationTexHadrAndLept << " \\end{table} " << endl;
  PresentationTexHadrAndLept.close();

  FMinusTexHadrAndLept.close();
  FPlusTexHadrAndLept.close();
  FZeroTexHadrAndLept.close();

  cout << "Finished running over all datasets..." << endl;
  fout->cd();

  ///////////////////
  // Writting
  //////////////////  
    
  cout << "Writing out..." << endl;
  fout->cd();

  TDirectory* tprofdir = fout->mkdir("TProfile_histograms");
  TDirectory* th1dir = fout->mkdir("1D_histograms");
  TDirectory* th2dir = fout->mkdir("2D_histograms_graphs");

  fout->cd();

  //Obtain TProfile histos for Reco-Gen vs Reco distribution:
  tprofdir->cd();

  th1dir->cd();

  // Write 1D histo's
  for(std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++){
    string name = it->first;
    if(name.find("CosTheta") != 0){
      TH1F *temp = it->second;
      temp->Write();
      TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
      tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
      //tempCanvas->Write();
      //		tempCanvas->SaveAs( (pathPNG+it->first+".pdf").c_str() );
    }
    else{
      TH1F *temp = it->second;
      temp->Write();
      TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
      tempCanvas->SaveAs( (pathPNG+"CosThetaPlots/"+it->first+".png").c_str() );
      //tempCanvas->Write();
    }
  }    

  fout->cd();
  
  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++){
    MultiSamplePlot *temp = it->second;
    string name = it->first;
    cout << " name : " << name << "    ";
    temp->Draw(false, name, false, false, true, true, true);
    temp->Write(fout, name, true, pathPNG+"MSPlot/");
    cout << " done " << endl;
  }

  // 2D
  th2dir->cd();
  for(std::map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++){
    TH2F *temp = it->second;
    temp->Write();
    TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
    tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
    //		tempCanvas->SaveAs( (pathPNG+it->first+".pdf").c_str() );
  }

  //Write TGraphAsymmErrors
  fout->cd();
  for(map<string,TGraphAsymmErrors*>::const_iterator it = graphAsymmErr.begin(); it != graphAsymmErr.end(); it++){
    TGraphAsymmErrors *temp = it->second;
    temp->Write();
    TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
    tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
    //		tempCanvas->SaveAs( (pathPNG+it->first+".pdf").c_str() );
  }
  
  //Write TGraphErrors
  fout->cd();
  for(map<string,TGraphErrors*>::const_iterator it = graphErr.begin(); it != graphErr.end(); it++){
    TGraphErrors *temp = it->second;
    temp->Write();
    TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
    tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
    //    tempCanvas->SaveAs( (pathPNG+it->first+".pdf").c_str() );
  }

  fout->Close();
  delete fout;
  
  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
  
  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "           hasn't crashed yet ;-)           " << endl;
  cout << "********************************************" << endl;
  
  return 0;
}
