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
//#include "../WHelicities/interface/MinuitFitter.h"
#include "../WHelicities/interface/MlbFitter.h"
#include "../WHelicities/interface/MinuitFitterPerforming.h"
//#include "../WHelicities/interface/KinematicFitClass.h"

//Includes for Minuit
//#include <TMinuitMinimizer.h>
//#include "TNtuple.h"  //A simple tree restricted to a list of float variables only. 

const int ndimen = 4;

const int CosThetaBinNumber = 35;       
const int NumberSSVHP=14;
const int NumberSSVHE=14;
const int NumberTCHP=14;
const int NumberTCHE=14;
const int NumberCSV=14;
const int NumberJP=14;
const int NumberJBP=14;
//TNtuple *genttbarhistoHadr[CosThetaBinNumber];  //This is the vector of ntuples containing the generated values of cos theta* for each cos theta* reconstructed bin
//TNtuple *genttbarhistoHadrAndLeptWOnly[CosThetaBinNumber];  //This is the vector of ntuples containing the generated values of cos theta* for each cos theta* reconstructed bin
//TNtuple *genttbarhistoHadrAndLept[CosThetaBinNumber];  //This is the vector of ntuples containing the generated values of cos theta* for each cos theta* reconstructed bin

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
  bool semiMuon = true;
  bool semiElectron = false;
  if(semiMuon == true){decayChannel = "SemiMu";}
  else if(semiElectron == true){decayChannel = "SemiEl";}

  string UsedTrigger;
  bool IsoMu172024Trigger = true;
  bool TriCentralJet30Trigger = false;
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
  bool bTagResults = true;
  cout << " Obtaining results for : Data = " << DataResults << " , JES = " << JESResults << " , JER = " << JERResults << " , WSyst = " << WSystResults << endl;

  //Different KinFit configurations studied:
  bool HadronicOnly = false;
  bool HadronicAndLeptonicWOnly = false;
  bool HadronicAndLeptonic = true;

  //Boolean to put MinuitFitter on and off:
  bool PerformMinuit = true;

  //Boolean to create MSPlots:
  bool CreateMSPlots = false;
 
  //1) Nominal samples:
  if(DataResults == true){
    inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_"+"Data_"+decayChannel+".root").c_str());
    nameDataSet.push_back("Data");
  }
  inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ST_SingleTop_tChannel_tbar_"+decayChannel+".root").c_str());
  nameDataSet.push_back("Nom_ST_SingleTop_tChannel_tbar");
  //inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ST_SingleTop_tChannel_t_"+decayChannel+".root").c_str());
  //nameDataSet.push_back("Nom_ST_SingleTop_tChannel_t");
  inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ST_SingleTop_tWChannel_tbar_"+decayChannel+".root").c_str());
  nameDataSet.push_back("Nom_ST_SingleTop_tWChannel_tbar");
  inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ST_SingleTop_tWChannel_t_"+decayChannel+".root").c_str());
  nameDataSet.push_back("Nom_ST_SingleTop_tWChannel_t");
  inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_TTbarJets_Other_"+decayChannel+".root").c_str());
  nameDataSet.push_back("Nom_TTbarJets_Other");
  if(semiMuon == true){
    inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_TTbarJets_SemiEl_"+decayChannel+".root").c_str());  //In muon channel case SemiEl is considered as background
    nameDataSet.push_back("Nom_TTbarJets_SemiEl");
  }
  else if(semiElectron == true){
    inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_TTbarJets_SemiMuon_"+decayChannel+".root").c_str());  //In electron channel case SemiMu is considered as background
    nameDataSet.push_back("Nom_TTbarJets_SemiMuon");
  }
  inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_WJets_"+decayChannel+".root").c_str());
  nameDataSet.push_back("Nom_WJets");
  inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ZJets_"+decayChannel+".root").c_str());
  nameDataSet.push_back("Nom_ZJets");  

  if(SignalOnly == true){
    inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ST_SingleTop_tChannel_tbar_"+decayChannel+".root").c_str());
    nameDataSet.push_back("McDta_ST_SingleTop_tChannel_tbar");
    //inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ST_SingleTop_tChannel_t_"+decayChannel+".root").c_str());
    //nameDataSet.push_back("McDta_ST_SingleTop_tChannel_t");
    inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ST_SingleTop_tWChannel_tbar_"+decayChannel+".root").c_str());
    nameDataSet.push_back("McDta_ST_SingleTop_tWChannel_tbar");
    inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ST_SingleTop_tWChannel_t_"+decayChannel+".root").c_str());
    nameDataSet.push_back("McDta_ST_SingleTop_tWChannel_t");
    inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_TTbarJets_Other_"+decayChannel+".root").c_str());
    nameDataSet.push_back("McDta_TTbarJets_Other");
    inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_TTbarJets_SemiEl_"+decayChannel+".root").c_str());  
    nameDataSet.push_back("McDta_TTbarJets_SemiEl");    
    inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_TTbarJets_SemiMuon_"+decayChannel+".root").c_str());  
    nameDataSet.push_back("McDta_TTbarJets_SemiMuon");
    inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_WJets_"+decayChannel+".root").c_str());
    nameDataSet.push_back("McDta_WJets");
    inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ZJets_"+decayChannel+".root").c_str());
    nameDataSet.push_back("McDta_ZJets");  
  }  
  else if(SignalOnly == false){
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
      //inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ST_SingleTop_tChannel_t_"+decayChannel+".root").c_str());
      //nameDataSet.push_back("TriggEvtSel_ST_SingleTop_tChannel_t");
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

    //11) BTag ScaleFactor Systematic
    if(bTagResults == true){
      //Consider bTagSyst Plus samples:
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ST_SingleTop_tChannel_tbar_"+decayChannel+".root").c_str());
      nameDataSet.push_back("bTagSystPlus_ST_SingleTop_tChannel_tbar");
      //      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_SingleTop_tChannel_t_"+decayChannel+".root").c_str());
      //      nameDataSet.push_back("bTagSystPlus_ST_SingleTop_tChannel_t");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ST_SingleTop_tWChannel_tbar_"+decayChannel+".root").c_str());
      nameDataSet.push_back("bTagSystPlus_ST_SingleTop_tWChannel_tbar");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ST_SingleTop_tWChannel_t_"+decayChannel+".root").c_str());
      nameDataSet.push_back("bTagSystPlus_ST_SingleTop_tWChannel_t");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_TTbarJets_Other_"+decayChannel+".root").c_str());
      nameDataSet.push_back("bTagSystPlus_TTbarJets_Other");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_TTbarJets_SemiEl_"+decayChannel+".root").c_str());
      nameDataSet.push_back("bTagSystPlus_TTbarJets_SemiEl");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_WJets_"+decayChannel+".root").c_str());    
      nameDataSet.push_back("bTagSystPlus_WJets");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ZJets_"+decayChannel+".root").c_str());
      nameDataSet.push_back("bTagSystPlus_ZJets");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_TTbarJets_SemiMuon_"+decayChannel+".root").c_str());
      nameDataSet.push_back("bTagSystPlus_TTbarJets_SemiMuon");

      //Consider bTagSyst Minus samples:
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ST_SingleTop_tChannel_tbar_"+decayChannel+".root").c_str());
      nameDataSet.push_back("bTagSystMinus_ST_SingleTop_tChannel_tbar");
      //inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ST_SingleTop_tChannel_t_"+decayChannel+".root").c_str());
      //nameDataSet.push_back("bTagSyst_ST_SingleTop_tChannel_t");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ST_SingleTop_tWChannel_tbar_"+decayChannel+".root").c_str());
      nameDataSet.push_back("bTagSystMinus_ST_SingleTop_tWChannel_tbar");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ST_SingleTop_tWChannel_t_"+decayChannel+".root").c_str());
      nameDataSet.push_back("bTagSystMinus_ST_SingleTop_tWChannel_t");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_TTbarJets_Other_"+decayChannel+".root").c_str());
      nameDataSet.push_back("bTagSystMinus_TTbarJets_Other");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_TTbarJets_SemiEl_"+decayChannel+".root").c_str());
      nameDataSet.push_back("bTagSystMinus_TTbarJets_SemiEl");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_WJets_"+decayChannel+".root").c_str());    
      nameDataSet.push_back("bTagSystMinus_WJets");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ZJets_"+decayChannel+".root").c_str());
      nameDataSet.push_back("bTagSystMinus_ZJets");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_TTbarJets_SemiMuon_"+decayChannel+".root").c_str());
      nameDataSet.push_back("bTagSystMinus_TTbarJets_SemiMuon");      
    }
  }//End of signalOnly = false loop
  
  //TTbarJets_SemiMuon sample should always be put as latest sample to avoid crash of TMinuitMinimizer !!
  if(semiMuon == true){
    inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_TTbarJets_SemiMuon_"+decayChannel+".root").c_str());
    nameDataSet.push_back("Nom_TTbarJets_SemiMuon");
  }
  else if(semiElectron == true){
    inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_TTbarJets_SemiEl_"+decayChannel+".root").c_str());
    nameDataSet.push_back("Nom_TTbarJets_SemiEl");
  }
  
  TFile *fout = new TFile (("WHelicities_Analysis_"+UsedTrigger+".root").c_str(), "RECREATE");
  fout->cd();

  //-----------------------------------------------
  //  Output .tex files for presentation/paper
  //-----------------------------------------------
  string bTagPerf = "BTagPerformanceStudy/";
  mkdir(("HadrKinFit/"+bTagPerf).c_str(),0777);
  mkdir(("HadrAndLeptWOnlyKinFit/"+bTagPerf).c_str(),0777);
  mkdir(("HadrAndLeptKinFit/"+bTagPerf).c_str(),0777);
 
  //oooooOOOOOOOOOOOoooooooooooOOOOOOOOOOOOOOOOoooooooooOOOOOOOOOOOoooooooooooooooooOOOOOOOOOOOOOOOoooooooooooooOOOOOOOOOOOoooooooo
  // Own results: No btag at event selection, top mass constraint for leptonic top in KinFit and offline muon pt cut of 27
  //oooooOOOOOOOOOOOoooooooooooOOOOOOOOOOOOOOOOoooooooooOOOOOOOOOOOoooooooooooooooooOOOOOOOOOOOOOOOoooooooooooooOOOOOOOOOOOoooooooo 

  //TriggEvtSel systematics initialization
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

  // initialize histograms
  histo1D["lumiWeights3D"]= new TH1F("lumiWeights3D","lumiWeights3D",25,0,50);  
  histo1D["JetCombHadr"] = new TH1F("JetCombHadr","JetCombHadr",14,-0.5,13.5);
  histo1D["JetCombHadrAndLeptWOnly"] = new TH1F("JetCombHadrAndLeptWOnly","JetCombHadrAndLeptWOnly",14,-0.5,13.5);
  histo1D["JetCombHadrAndLept"] = new TH1F("JetCombHaAndLeptdr","JetCombHadrAndLept",14,-0.5,13.5);
  histo1D["LeptWMassHadr"] = new TH1F("LeptWMassHadr","LeptWMassHadr",40,0,150);
  histo1D["LeptTopMassHadr"] = new TH1F("LeptTopMassHadr","LeptTopMassHadr",60,50,350);

  histo1D["CosThetaNoBTagData"] = new TH1F("CosThetaNoBTagData","CosThetaNoBTagData",CosThetaBinNumber,-1,1);
  histo1D["CosThetaNoBTagSemiMu"] = new TH1F("CosThetaNoBTagSemiMu","CosThetaNoBTagSemiMu",CosThetaBinNumber,-1,1);
  histo1D["CosThetaNoBTagSemiEl"] = new TH1F("CosThetaNoBTagSemiEl","CosThetaNoBTagSemiEl",CosThetaBinNumber,-1,1);
  histo1D["CosThetaNoBTagOther"] = new TH1F("CosThetaNoBTagOther","CosThetaNoBTagOther",CosThetaBinNumber,-1,1);
  histo1D["CosThetaNoBTagWJets"] = new TH1F("CosThetaNoBTagWJets","CosThetaNoBTagWJets",CosThetaBinNumber,-1,1);
  histo1D["CosThetaNoBTagZJets"] = new TH1F("CosThetaNoBTagZJets","CosThetaNoBTagZJets",CosThetaBinNumber,-1,1);
  histo1D["CosThetaNoBTagSTtChtBar"] = new TH1F("CosThetaNoBTagSTtChtBar","CosThetaNoBTagSTtChtBar",CosThetaBinNumber,-1,1);
  histo1D["CosThetaNoBTagSTtCht"] = new TH1F("CosThetaNoBTagSTtCht","CosThetaNoBTagSTtCht",CosThetaBinNumber,-1,1);
  histo1D["CosThetaNoBTagSTtWChtBar"] = new TH1F("CosThetaNoBTagSTtWChtBar","CosThetaNoBTagSTtWChtBar",CosThetaBinNumber,-1,1);
  histo1D["CosThetaNoBTagSTtWCht"] = new TH1F("CosThetaNoBTagSTtWCht","CosThetaNoBTagSTtWCht",CosThetaBinNumber,-1,1);

  histo1D["CosThetaOneBTagData"] = new TH1F("CosThetaOneBTagData","CosThetaOneBTagData",CosThetaBinNumber,-1,1);
  histo1D["CosThetaOneBTagSemiMu"] = new TH1F("CosThetaOneBTagSemiMu","CosThetaOneBTagSemiMu",CosThetaBinNumber,-1,1);
  histo1D["CosThetaOneBTagSemiEl"] = new TH1F("CosThetaOneBTagSemiEl","CosThetaOneBTagSemiEl",CosThetaBinNumber,-1,1);
  histo1D["CosThetaOneBTagOther"] = new TH1F("CosThetaOneBTagOther","CosThetaOneBTagOther",CosThetaBinNumber,-1,1);
  histo1D["CosThetaOneBTagWJets"] = new TH1F("CosThetaOneBTagWJets","CosThetaOneBTagWJets",CosThetaBinNumber,-1,1);
  histo1D["CosThetaOneBTagZJets"] = new TH1F("CosThetaOneBTagZJets","CosThetaOneBTagZJets",CosThetaBinNumber,-1,1);
  histo1D["CosThetaOneBTagSTtChtBar"] = new TH1F("CosThetaOneBTagSTtChtBar","CosThetaOneBTagSTtChtBar",CosThetaBinNumber,-1,1);
  histo1D["CosThetaOneBTagSTtCht"] = new TH1F("CosThetaOneBTagSTtCht","CosThetaOneBTagSTtCht",CosThetaBinNumber,-1,1);
  histo1D["CosThetaOneBTagSTtWChtBar"] = new TH1F("CosThetaOneBTagSTtWChtBar","CosThetaOneBTagSTtWChtBar",CosThetaBinNumber,-1,1);
  histo1D["CosThetaOneBTagSTtWCht"] = new TH1F("CosThetaOneBTagSTtWCht","CosThetaOneBTagSTtWCht",CosThetaBinNumber,-1,1);

  histo1D["CosThetaTwoBTagData"] = new TH1F("CosThetaTwoBTagData","CosThetaTwoBTagData",CosThetaBinNumber,-1,1);
  histo1D["CosThetaTwoBTagSemiMu"] = new TH1F("CosThetaTwoBTagSemiMu","CosThetaTwoBTagSemiMu",CosThetaBinNumber,-1,1);
  histo1D["CosThetaTwoBTagSemiEl"] = new TH1F("CosThetaTwoBTagSemiEl","CosThetaTwoBTagSemiEl",CosThetaBinNumber,-1,1);
  histo1D["CosThetaTwoBTagOther"] = new TH1F("CosThetaTwoBTagOther","CosThetaTwoBTagOther",CosThetaBinNumber,-1,1);
  histo1D["CosThetaTwoBTagWJets"] = new TH1F("CosThetaTwoBTagWJets","CosThetaTwoBTagWJets",CosThetaBinNumber,-1,1);
  histo1D["CosThetaTwoBTagZJets"] = new TH1F("CosThetaTwoBTagZJets","CosThetaTwoBTagZJets",CosThetaBinNumber,-1,1);
  histo1D["CosThetaTwoBTagSTtChtBar"] = new TH1F("CosThetaTwoBTagSTtChtBar","CosThetaTwoBTagSTtChtBar",CosThetaBinNumber,-1,1);
  histo1D["CosThetaTwoBTagSTtCht"] = new TH1F("CosThetaTwoBTagSTtCht","CosThetaTwoBTagSTtCht",CosThetaBinNumber,-1,1);
  histo1D["CosThetaTwoBTagSTtWChtBar"] = new TH1F("CosThetaTwoBTagSTtWChtBar","CosThetaTwoBTagSTtWChtBar",CosThetaBinNumber,-1,1);
  histo1D["CosThetaTwoBTagSTtWCht"] = new TH1F("CosThetaTwoBTagSTtWCht","CosThetaTwoBTagSTtWCht",CosThetaBinNumber,-1,1);

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
    if( dataSetName.find("Nom_TTbarJets_SemiMuon") == 0) color = kRed+1;
    if( dataSetName.find("Nom_TTbarJets_SemiEl") == 0) color = kRed-4;
    if( dataSetName.find("Nom_TTbarJets_Other") == 0) color = kRed-7;
    if( dataSetName.find("TopMass") == 0) color = kRed+1;
    if( dataSetName.find("Nom_WJets") == 0 ){
      dataSet->SetTitle("W#rightarrowl#nu");
      color = kGreen-3;
    }
    if( dataSetName.find("Nom_ZJets") == 0){
      dataSet->SetTitle("Z/#gamma*#rightarrowl^{+}l^{-}");
      color = kAzure-2;
    }
    if( dataSetName.find("Nom_ST") == 0) color = kMagenta;
    //if( dataSet->Name().find("QCD") == 0 && dataSet->Name().find("JES") !=0 ) color = kBlue;        
      
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

    MSPlot["CosThetaChiSqOnlyNoBTag"]= new MultiSamplePlot(datasets,"CosThetaChiSqOnlyNoBTag",CosThetaBinNumber,-1,2,"CosThetaChiSqOnlyNoBTag");
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

    MSPlot["HadronicWMassNoBTag"]= new MultiSamplePlot(datasets,"HadronicWMassNoBTag",50,78,82,"HadronicWMassNoBTag");
    MSPlot["HadronicTopMassNoBTag"]= new MultiSamplePlot(datasets,"HadronicTopMassNoBTag",50,170,174,"HadronicTopMassNoBTag");
    MSPlot["LeptonicWMassNoBTag"]= new MultiSamplePlot(datasets,"LeptonicWMassNoBTag",50,78,82,"LeptonicWMassNoBTag");
    MSPlot["LeptonicTopMassNoBTag"]= new MultiSamplePlot(datasets,"LeptonicTopMassNoBTag",50,170,174,"LeptonicTopMassNoBTag");
    
    //Specific b-tag applied
    MSPlot["JetCombCosThetaFoundSpecificBTag"]= new MultiSamplePlot(datasets,"JetCombCosThetaFoundSpecificBTag",13,-0.5,12.5,"JetCombCosThetaFoundSpecificBTag");
    MSPlot["ChiSqKinFitCosThetaFoundSpecificBTag"]= new MultiSamplePlot(datasets,"ChiSqKinFitCosThetaFoundSpecificBTag",40,0,150,"ChiSqKinFitCosThetaFoundSpecificBTag");
    MSPlot["ProbKinFitCosThetaFoundSpecificBTag"]= new MultiSamplePlot(datasets,"ProbKinFitCosThetaFoundSpecificBTag",20,0,1,"ProbKinFitCosThetaFoundSpecificBTag");

    MSPlot["CosThetaSpecificBTag35Bins"]= new MultiSamplePlot(datasets,"CosThetaSpecificBTag35Bins",35,-1,1,"CosThetaSpecificBTag35Bins");
    MSPlot["CosThetaSpecificBTag25Bins"]= new MultiSamplePlot(datasets,"CosThetaSpecificBTag25Bins",25,-1,1,"CosThetaSpecificBTag25Bins");
    MSPlot["CosThetaSpecificDoubleBTag35Bins"]= new MultiSamplePlot(datasets,"CosThetaSpecificDoubleBTag35Bins",35,-1,1,"CosThetaSpecificDoubleBTag35Bins");
    MSPlot["CosThetaSpecificDoubleBTag25Bins"]= new MultiSamplePlot(datasets,"CosThetaSpecificDoubleBTag25Bins",25,-1,1,"CosThetaSpecificDoubleBTag25Bins");

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

    MSPlot["HadronicWMassSpecificBTag"]= new MultiSamplePlot(datasets,"HadronicWMassSpecificBTag",50,78,82,"HadronicWMassSpecificBTag");
    MSPlot["HadronicTopMassSpecificBTag"]= new MultiSamplePlot(datasets,"HadronicTopMassSpecificBTag",50,170,174,"HadronicTopMassSpecificBTag");
    MSPlot["LeptonicWMassSpecificBTag"]= new MultiSamplePlot(datasets,"LeptonicWMassSpecificBTag",50,78,82,"LeptonicWMassSpecificBTag");
    MSPlot["LeptonicTopMassSpecificBTag"]= new MultiSamplePlot(datasets,"LeptonicTopMassSpecificBTag",50,170,174,"LeptonicTopMassSpecificBTag");
   
    //Check nPrimary vertices for different executed cuts !!
    MSPlot["nPVBeforeCuts"] = new MultiSamplePlot(datasets,"nPVBeforeCuts" , 20, 0, 20, "nPVBeforeCuts");
    MSPlot["nPVAfterTransverseMassCut"] = new MultiSamplePlot(datasets,"nPVAfterTransverseMassCut" , 20, 0, 20, "nPVAfterTransverseMassCut");

    MSPlot["TransverseMassBeforeCut"]=new MultiSamplePlot(datasets,"TransverseMassBeforeCut",50,0,200,"TransverseMassBeforeCut");
    MSPlot["TransverseMassBeforeCutTopMassPlus"]=new MultiSamplePlot(datasets,"TransverseMassBeforeCutTopMassPlus",50,0,200,"TransverseMassBeforeCutTopMassPlus");
    MSPlot["TransverseMassBeforeCutTopMassMinus"]=new MultiSamplePlot(datasets,"TransverseMassBeforeCutTopMassMinus",50,0,200,"TransverseMassBeforeCutTopMassMinus");
    MSPlot["TransverseMassAfterCut"]=new MultiSamplePlot(datasets,"TransverseMassAfterCut",50,0,200,"TransverseMassAfterCut"); 
    MSPlot["TransverseMassAfterCutTopMassPlus"]=new MultiSamplePlot(datasets,"TransverseMassAfterCutTopMassPlus",50,0,200,"TransverseMassAfterCutTopMassPlus");
    MSPlot["TransverseMassAfterCutTopMassMinus"]=new MultiSamplePlot(datasets,"TransverseMassAfterCutTopMassMinus",50,0,200,"TransverseMassAfterCutTopMassMinus");
  
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
  
    MSPlot["HadrTopMassSelJetsNobTag"]=new MultiSamplePlot(datasets,"HadrTopMassSelJetsNobTag",100,0,250,"HadrTopMassSelJetsNobTag");
    MSPlot["HadrTopMassSelJetsOnebTag"]=new MultiSamplePlot(datasets,"HadrTopMassSelJetsOnebTag",100,0,250,"HadrTopMassSelJetsOnebTag");
    MSPlot["HadrTopMassSelJetsTwobTag"]=new MultiSamplePlot(datasets,"HadrTopMassSelJetsTwobTag",100,0,250,"HadrTopMassSelJetsTwobTag");

    //Control plots chapter:
    MSPlot["METAfterMTCut"]=new MultiSamplePlot(datasets,"METAfterMTCut",100,0,500,"METAfterMTCut");
    MSPlot["4thJetPtBeforeMTCut"]=new MultiSamplePlot(datasets,"4thJetPtBeforeMTCut",50,0,500,"4thJetPtBeforeMTCut");
    MSPlot["METAfterMCSVbTag"]=new MultiSamplePlot(datasets,"METAfterMCSVbTag",100,0,500,"METAfterMCSVbTag");
    MSPlot["4thJetPtAfterMCSVbTag"]=new MultiSamplePlot(datasets,"4thJetPtAfterMCSVbTag",50,0,500,"4thJetPtAfterMCSVbTag");

    MSPlot["CosThetaOneMCSVBTag"]= new MultiSamplePlot(datasets,"CosThetaOneMCSVBTag",CosThetaBinNumber,-1,1,"CosThetaOneMCSVBTag");
    MSPlot["CosThetaTwoMCSVBTag"]= new MultiSamplePlot(datasets,"CosThetaTwoMCSVBTag",CosThetaBinNumber,-1,1,"CosThetaTwoMCSVBTag");
    MSPlot["4thJetPtTwoMCSVBTag"]=new MultiSamplePlot(datasets,"4thJetPtTwoMCSVBTag",50,0,500,"4thJetPtTwoMCSVBTag");

    MSPlot["HadrTopMassOrigKins"]=new MultiSamplePlot(datasets,"HadrTopMassOrigKins",40,0,400,"HadrTopMassOrigKins");
    MSPlot["HadrTopMassOrigKinsWJetsPlus"]=new MultiSamplePlot(datasets,"HadrTopMassOrigKinsWJetsPlus",40,0,400,"HadrTopMassOrigKinsWJetsPlus");
    MSPlot["HadrTopMassOrigKinsWJetsMinus"]=new MultiSamplePlot(datasets,"HadrTopMassOrigKinsWJetsMinus",40,0,400,"HadrTopMassOrigKinsWJetsMinus");
        
    MSPlot["CosThetaOrigKin"]=new MultiSamplePlot(datasets,"CosThetaOrigKin",CosThetaBinNumber,-1,1,"CosThetaOrigKin");
    MSPlot["NeutrinoPzOrigKin"]=new MultiSamplePlot(datasets,"NeutrinoPzOrigKin",50,-500,500,"NeutrinoPzOrigKin");
    MSPlot["METPtOrigKin"]=new MultiSamplePlot(datasets,"METPtOrigKin",100,0,500,"METPtOrigKin");
    MSPlot["LeptonPzOrigKin"]=new MultiSamplePlot(datasets,"LeptonPzOrigKin",50,-500,500,"LeptonPzOrigKin");
    MSPlot["LeptonPtOrigKin"]=new MultiSamplePlot(datasets,"LeptonPtOrigKin",50,0,500,"LeptonPtOrigKin");

    MSPlot["LightJet1PtOrigKin"]=new MultiSamplePlot(datasets,"LightJet1PtOrigKin",100,0,500,"LightJet1PtOrigKin");
    MSPlot["LightJet2PtOrigKin"]=new MultiSamplePlot(datasets,"LightJet2PtOrigKin",100,0,500,"LightJet2PtOrigKin");
    MSPlot["LeptBJetPtOrigKin"]=new MultiSamplePlot(datasets,"LeptBJetPtOrigKin",100,0,500,"LeptBJetPtOrigKin");
    MSPlot["HadrBJetPtOrigKin"]=new MultiSamplePlot(datasets,"HadrBJetPtOrigKin",100,0,500,"HadrBJetPtOrigKin");

    MSPlot["LeptonPtAfterMTCut"]=new MultiSamplePlot(datasets,"LeptonPtAfterMTCut",20,0,400,"LeptonPtAfterMTCut");
    MSPlot["LeptonPtMCSVbTag"]=new MultiSamplePlot(datasets,"LeptonPtMCSVbTag",20,0,400,"LeptonPtMCSVbTag");
    MSPlot["Jet1PtMCSVbTag"]=new MultiSamplePlot(datasets,"Jet1PtMCSVbTag",20,0,400,"Jet1PtMCSVbTag");

    MSPlot["MlbNoBTag"] = new MultiSamplePlot(datasets,"MlbNoBTag",100,0,160,"MlbNoBTag");

    MSPlot["MlbBeforeTransMassCut"]  = new MultiSamplePlot(datasets,"MlbBeforeTransMassCut",200,0,310,"MlbBeforeTransMassCut");
    MSPlot["MlbAfterTransMassCut"]  = new MultiSamplePlot(datasets,"MlbAfterTransMassCut",200,0,310,"MlbAfterTransMassCut");
    MSPlot["MlbAfterTransMassCutAfterMCSVbTag"]  = new MultiSamplePlot(datasets,"MlbAfterTransMassCutAfterMCSVbTag",200,0,310,"MlbAfterTransMassCutAfterMCSVbTag");
    MSPlot["MlbAfterTransMassCutAfterMCSVbTagSFAdded"] = new MultiSamplePlot(datasets,"MlbAfterTransMassCutAfterMCSVbTagSFAdded",200,0,310,"MlbAfterTransMassCutAfterMCSVbTagSFAdded");

    MSPlot["TransverseMassExistingJetComb"] = new MultiSamplePlot(datasets,"TransverseMassExistingJetComb",100,0,200,"TransverseMassExistingJetComb");
    MSPlot["TransverseMassExistingCos"] = new MultiSamplePlot(datasets,"TransverseMassExistingCos",100,0,200,"TransverseMassExistingCos");
   
    MSPlot["nPVAll"] = new MultiSamplePlot(datasets,"nPVAll",100,0,20,"nPVAll");
    MSPlot["nPVExistingJetComb"] = new MultiSamplePlot(datasets,"nPVExistingJetComb",100,0,20,"nPVExistingJetComb");
    MSPlot["nPVExistingCos"] = new MultiSamplePlot(datasets,"nPVExistingCos",100,0,20,"nPVExistingCos");
    MSPlot["nPVExistingCosbTag"] = new MultiSamplePlot(datasets,"nPVExistingCosbTag",100,0,20,"nPVExistingCosbTag");

    MSPlot["bTagCSVbLeptDistribution"] = new MultiSamplePlot(datasets,"bTagCSVbLeptDistribution",100,-1.1,1.1,"bTagCSVbLeptDistribution");
    MSPlot["bTagCSVbLeptDistributionAfterCut"] = new MultiSamplePlot(datasets,"bTagCSVbLeptDistributionAfterCut",100,-1.1,1.1,"bTagCSVbLeptDistributionAfterCut");
    MSPlot["bTagCSVbLeptDistributionClass"] = new MultiSamplePlot(datasets,"bTagCSVbLeptDistributionClass",100,-1.1,1.1,"bTagCSVbLeptDistributionClass");
    MSPlot["bTagCSVbLeptDistributionClassAfterCut"] = new MultiSamplePlot(datasets,"bTagCSVbLeptDistributionClassAfterCut",100,-1.1,1.1,"bTagCSVbLeptDistributionClassAfterCut");
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
  Lumi3DReWeighting Lumi3DWeightsMinus;
  Lumi3DReWeighting Lumi3DWeightsPlus;
  Lumi3DWeights = Lumi3DReWeighting("PileUpReweighting/pileup_MC_Fall11.root","PileUpReweighting/pileup_FineBin_2011Data_UpToRun180252.root", "pileup", "pileup");
  Lumi3DWeightsMinus = Lumi3DReWeighting("PileUpReweighting/pileup_MC_Fall11.root","PileUpReweighting/pileup_FineBin_2011Data_UpToRun180252.root", "pileup", "pileup");
  Lumi3DWeightsPlus = Lumi3DReWeighting("PileUpReweighting/pileup_MC_Fall11.root","PileUpReweighting/pileup_FineBin_2011Data_UpToRun180252.root", "pileup", "pileup");

  if(PUResults == true){
    Lumi3DWeightsMinus.weight3D_init(0.92);
    Lumi3DWeightsPlus.weight3D_init(1.08);
  }
  Lumi3DWeights.weight3D_init(1.0);

  /////////////////////////////////////////
  // Initializing used variables
  /////////////////////////////////////////
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
    if(dataSetName.find("TopMassPlus") == 0 && TopMassResults == true){
      dataSet->SetEquivalentLuminosity( 9572.5587 );
    }
    if(dataSetName.find("TopMassMinus") == 0 && TopMassResults == true){
      dataSet->SetEquivalentLuminosity( 9996.5650 );
    }
    cout << "After changing --> Cross section = " << dataSet->Xsection() << "  intLumi = " << dataSet->EquivalentLumi() << " Normfactor = " << dataSet->NormFactor() << endl;
    cout << " ************************************** " << endl;

    //BTag ScaleFactors:  --> Put here to avoid keeping substracting the uncertainties !!
    float LooseSF[7] = {0.95,0.93,1.,1.,1.01,0.97,0.98};
    float MediumSF[7] = {0.96,0.93,0.96,1.,0.97,0.95,0.94};
    float TightSF[7] = {0.95,0.93,0.95,0.95,0.97,0.9,0.89};
    float SFUnc = 0.04;
    
    //Apply bTaggings systematics:
    if(bTagResults == true){
      if(dataSetName.find("bTagSystMinus") == 0){
	for(int ii = 0; ii<7; ii++){
	  LooseSF[ii]=LooseSF[ii]-SFUnc;
	  MediumSF[ii]=MediumSF[ii]-SFUnc;
	  TightSF[ii]=TightSF[ii]-SFUnc;
	}
      }
      if(dataSetName.find("bTagSystPlus") == 0){
	for(int ii = 0; ii<7; ii++){
	  LooseSF[ii]=LooseSF[ii]+SFUnc;
	  MediumSF[ii]=MediumSF[ii]+SFUnc;
	  TightSF[ii]=TightSF[ii]+SFUnc;
	}
      }
    }
    
    //Value needed to study the reconstruction efficiency of the leptonic b-jet:
    int CorrectBLeptConstruction=0;
    int ReconstructedEvents=0;

    vector<float> CosThetaValuesKinFitHadr[TotalNumberbTags],CosThetaValuesKinFitHadrAndLeptWOnly[TotalNumberbTags],CosThetaValuesKinFitHadrAndLept[TotalNumberbTags];
    int FilledEntries = 0;
    vector<double> CosThGenKinFitHadr[TotalNumberbTags],CosThGenKinFitHadrAndLeptWOnly[TotalNumberbTags],CosThGenKinFitHadrAndLept[TotalNumberbTags],EventCorrectionWeightKinFitHadr[TotalNumberbTags],EventCorrectionWeightKinFitHadrAndLeptWOnly[TotalNumberbTags],EventCorrectionWeightKinFitHadrAndLept[TotalNumberbTags],MlbValuesKinFitHadrAndLept[TotalNumberbTags];
    float binEdge[CosThetaBinNumber+1];
    float binSize = (1.-(-1.))/CosThetaBinNumber;
    for(int ii=0; ii<=CosThetaBinNumber;ii++){
      binEdge[ii] = -1 + binSize*ii;
    }

    //Initialize naming of different bTag options:
    int TCHELoop = 1;
    int TCHPLoop = 1;
    int SSVHELoop = 1;
    int SSVHPLoop = 1;
    int CSVLoop =1;
    int JPLoop = 1;
    int JBPLoop = 1;
       
    while(JBPLoop <= NumberJBPbTags){
      while(JPLoop <= NumberJPbTags){
	while(CSVLoop<=NumberCSVbTags){
	  while(SSVHPLoop<=NumberSSVHPbTags){  
	    while(SSVHELoop<=NumberSSVHEbTags){
	      while(TCHPLoop<=NumberTCHPbTags){
		while(TCHELoop<=NumberTCHEbTags){
		  if(iDataSet==0){
		    PresentationOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = bTagName.NameGivingPres(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags,JPLoop, NumberJPbTags,JBPLoop, NumberJBPbTags);

		  }
		  TCHELoop++;
		  if(TCHELoop == 14){TCHPLoop=2;SSVHELoop=2;SSVHPLoop=2;CSVLoop=2;JPLoop=2;JBPLoop=2;}
		}
		if(iDataSet==0){
		  PresentationOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = bTagName.NameGivingPres(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags,JPLoop, NumberJPbTags,JBPLoop, NumberJBPbTags);
		}		
		TCHPLoop++;
	      }
	      if(iDataSet==0){
		PresentationOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = bTagName.NameGivingPres(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags,JPLoop, NumberJPbTags,JBPLoop, NumberJBPbTags);
	      }
	      SSVHELoop++;
	    }
	    if(iDataSet==0){
	      PresentationOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = bTagName.NameGivingPres(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags,JPLoop, NumberJPbTags,JBPLoop, NumberJBPbTags);
	    }
	    SSVHPLoop++;		
	  }	
	  if(iDataSet==0){
	    PresentationOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = bTagName.NameGivingPres(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags,JPLoop, NumberJPbTags,JBPLoop, NumberJBPbTags);
	  }	  
	  CSVLoop++;
	}
	if(iDataSet==0){
	  PresentationOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = bTagName.NameGivingPres(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags,JPLoop, NumberJPbTags,JBPLoop, NumberJBPbTags);
	}
	JPLoop++;
      }
      if(iDataSet==0){
	PresentationOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = bTagName.NameGivingPres(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags,JPLoop, NumberJPbTags,JBPLoop, NumberJBPbTags);
      }
      JBPLoop++;
    }
           
    /////////////////////////////////////////
    // Loop on events
    /////////////////////////////////////////  
    for(unsigned int iEvt=0; iEvt<nEvent; iEvt++){
      //for(unsigned int iEvt=0; iEvt<10; iEvt++){
      
      inWTreeTree->GetEvent(iEvt);
      if(iEvt%10000 == 0)
        std::cout<<"Processing the "<<iEvt<<"th event, time = "<< ((double)clock() - start) / CLOCKS_PER_SEC <<endl;
      
      // PU reweighting
      float avPU = ( (float)wTree->nPUBXm1() + (float)wTree->nPU() + (float)wTree->nPUBXp1() ) / 3.; // average in 3 BX!!!, as recommended
      
      double lumiWeight3D=1.0;
      if(dataSetName.find("PUMinus_") == 0)
	lumiWeight3D = Lumi3DWeightsMinus.weight3D(wTree->nPUBXm1(),wTree->nPU(),wTree->nPUBXp1());
      else if(dataSetName.find("PUPlus_") == 0)
	lumiWeight3D = Lumi3DWeightsPlus.weight3D(wTree->nPUBXm1(),wTree->nPU(),wTree->nPUBXp1());
      else if(!(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" || dataSetName.find("PU") == 0))
	lumiWeight3D = Lumi3DWeights.weight3D(wTree->nPUBXm1(),wTree->nPU(),wTree->nPUBXp1());    
      else if(dataSetName == "Data")
	lumiWeight3D = 1.0;
      histo1D["lumiWeights3D"]->Fill(lumiWeight3D);

      ////////////////////////////////////////////////////////////////
      //    Get variables from wTree needed for Analysis            //
      ////////////////////////////////////////////////////////////////
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
	
	if(dataSetName.find("TriggEvtSelMinus") == 0){
	  meanEvtSelEffSF_=meanEvtSelEffSF_-0.004;
	  meanTriggerEffSF_=meanTriggerEffSF_-0.001;
	  SFjet=SFjet-errorDownSFjet;
	}
	else if(dataSetName.find("TriggEvtSelPlus") == 0){
	  meanEvtSelEffSF_=meanEvtSelEffSF_+0.004;
	  meanTriggerEffSF_=meanTriggerEffSF_+0.001;
	  SFjet=SFjet-errorUpSFjet;
	}
	float SFele = meanTriggerEffSF_ * (factorIso_ + factorNonIso_ / (0.997+0.05*relativeLeptonIsolation-3.*relativeLeptonIsolation*relativeLeptonIsolation));
	scaleFactor = meanEvtSelEffSF_*SFele*SFjet;
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
	SFmuon = yVec[iPoint];
	errorUpSFmuon = yErrHi[iPoint];
	errorDownSFmuon = yErrLo[iPoint];	
	
	scaleFactor = SFmuon;
	if(dataSetName.find("TriggEvtSelMinus") == 0)
	  scaleFactor = SFmuon-errorDownSFmuon;	
	else if(dataSetName.find("TriggEvtSelPlus") == 0)
	  scaleFactor = SFmuon+errorUpSFmuon;	
      }

      //////////////////////////////////////////////////////////////////////////////////////////////////
      if(CreateMSPlots == true){
	MSPlot["DeltaRJetLepton"]->Fill(DeltaRJetLepton,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D); 
	MSPlot["DeltaRMuonJet"]->Fill(DeltaRMuonJet,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["DeltaRElectronJet"]->Fill(DeltaRElectronJet, datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
      }                      
      
      //Chi squared values for different KinFit configurations:  --> Original kinematics Chi Squared with calculation class !!
      vector<float> ChiSqKinFitHadr;
      vector<float> KinFitProbHadr;
      vector<float> ChiSqKinFitHadrAndLeptWOnly;
      vector<float> KinFitProbHadrAndLeptWOnly;
      vector<float> ChiSqKinFitHadrAndLept;
      vector<float> KinFitProbHadrAndLept;
      
      //Original kinematics initializing:
      //vector<TLorentzVector> leptBJetOrigKins, hadrBJetOrigKins, light1JetOrigKins, light2JetOrigKins, leptonOrigKins;  //Get this with different method??

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
	hadrBChangedHadrAndLept.push_back(wTree->hadrAndLeptKinFitHadrB(iCombi));
	light1ChangedHadrAndLept.push_back(wTree->hadrAndLeptKinFitLight1(iCombi));
	light2ChangedHadrAndLept.push_back(wTree->hadrAndLeptKinFitLight2(iCombi));       
		
      }//End loop for filling KinFit stuff!!

      float TransverseMass = sqrt(2*(abs(lepton.Pt()))*abs(MET.Pt())*(1-cos(lepton.DeltaPhi(MET))));	//Should this not be updated to correct lepton and MET?? --> Not really a difference of lepton and MET after KinFit
      if(CreateMSPlots == true){
	if(dataSetName.find("TopMass") != 0) MSPlot["TransverseMassBeforeCut"]->Fill(TransverseMass, datasets[iDataSet], true, Luminosity);
	if(dataSetName.find("TopMassMinus") != 0 && dataSetName.find("Nom_TTbarJets_") != 0) MSPlot["TransverseMassBeforeCutTopMassPlus"]->Fill(TransverseMass, datasets[iDataSet], true, Luminosity);
	if(dataSetName.find("TopMassPlus") != 0 && dataSetName.find("Nom_TTbarJets_") != 0) MSPlot["TransverseMassBeforeCutTopMassMinus"]->Fill(TransverseMass, datasets[iDataSet], true, Luminosity);
	MSPlot["4thJetPtBeforeMTCut"]->Fill(selectedJets[3].Pt(), datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 	
	//MSPlot["RelIsoLepton"]->Fill(relativeLeptonIsolation, datasets[iDataSet], true, Luminosity);
      }

      //Select KinFit combination with lowest chi squared using original kinematics:
      float MassTop;
      if(dataSetName.find("Data") == 0) MassTop = 173.1;
      else MassTop = 172.5;
      
      //Select lowest chi squared value for HadrAndLept KinFit configuration:
      int SelectedJetCombinationHadrAndLept=999;
      float LowestChiSq=999.;
      for(int ii=0;ii<12;ii++){
	if(LowestChiSq > ChiSqKinFitHadrAndLept[ii]){
	  LowestChiSq = ChiSqKinFitHadrAndLept[ii];
	  SelectedJetCombinationHadrAndLept=ii;
	}
      }
      
      if(SelectedJetCombinationHadrAndLept != 999){
	if(CreateMSPlots == true)
	  MSPlot["TransverseMassExistingJetComb"]->Fill(TransverseMass, datasets[iDataSet], true, Luminosity);

	float CosThetaOrigKins = bTagCosThetaCalculation.CalcOrigKins(BLeptIndex[SelectedJetCombinationHadrAndLept], BHadrIndex[SelectedJetCombinationHadrAndLept],  lepton, selectedJets, 80.4, MassTop);
	if(CosThetaOrigKins != 999){
	  if(CreateMSPlots == true)
	    MSPlot["TransverseMassExistingCos"]->Fill(TransverseMass, datasets[iDataSet], true, Luminosity);

	  TLorentzVector neutrinoOrigKins = bTagCosThetaCalculation.GetNeutrino();
	  TLorentzVector WLeptOrigKins = (neutrinoOrigKins+lepton);
      
	  if(CreateMSPlots == true){
	    MSPlot["CosThetaOrigKin"]->Fill(CosThetaOrigKins, datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	    MSPlot["METPtOrigKin"]->Fill(neutrinoOrigKins.Pt(), datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	    MSPlot["NeutrinoPzOrigKin"]->Fill(neutrinoOrigKins.Pz(), datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	    MSPlot["LeptonPzOrigKin"]->Fill(lepton.Pz(), datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	    MSPlot["LeptonPtOrigKin"]->Fill(lepton.Pt(), datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	    MSPlot["LeptBJetPtOrigKin"]->Fill(selectedJets[BLeptIndex[SelectedJetCombinationHadrAndLept]].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	    MSPlot["HadrBJetPtOrigKin"]->Fill(selectedJets[BHadrIndex[SelectedJetCombinationHadrAndLept]].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	    MSPlot["LightJet1PtOrigKin"]->Fill(selectedJets[Quark1Index[SelectedJetCombinationHadrAndLept]].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	    MSPlot["LightJet2PtOrigKin"]->Fill(selectedJets[Quark2Index[SelectedJetCombinationHadrAndLept]].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	    MSPlot["MlbBeforeTransMassCut"]->Fill((selectedJets[BLeptIndex[SelectedJetCombinationHadrAndLept]]+lepton).M(), datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	    
	    MSPlot["HadrTopMassOrigKins"]->Fill((selectedJets[BHadrIndex[SelectedJetCombinationHadrAndLept]]+selectedJets[Quark2Index[SelectedJetCombinationHadrAndLept]]+selectedJets[Quark2Index[SelectedJetCombinationHadrAndLept]]).M(), datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	    if(dataSetName.find("WSystMinus_WJets") != 0 && WSystResults == true)
	      MSPlot["HadrTopMassOrigKinsWJetsPlus"]->Fill((selectedJets[BHadrIndex[SelectedJetCombinationHadrAndLept]]+selectedJets[Quark2Index[SelectedJetCombinationHadrAndLept]]+selectedJets[Quark2Index[SelectedJetCombinationHadrAndLept]]).M(), datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	    if(dataSetName.find("WSystPlus_WJets") != 0 && WSystResults == true)
	      MSPlot["HadrTopMassOrigKinsWJetsMinus"]->Fill((selectedJets[BHadrIndex[SelectedJetCombinationHadrAndLept]]+selectedJets[Quark2Index[SelectedJetCombinationHadrAndLept]]+selectedJets[Quark2Index[SelectedJetCombinationHadrAndLept]]).M(), datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);	
	    
	    if(TransverseMass > 30&& DeltaRJetLepton > 0.3 && ((DeltaRMuonJet > 0.3 && semiMuon == true) || (DeltaRElectronJet > 0.3 && semiElectron == true)) && relativeLeptonIsolation <= 0.125){
	      MSPlot["METAfterMTCut"]->Fill(neutrinoOrigKins.Pt(), datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 
	      MSPlot["LeptonPtAfterMTCut"]->Fill(lepton.Pt(), datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);	
	      MSPlot["MlbAfterTransMassCut"]->Fill((selectedJets[BLeptIndex[SelectedJetCombinationHadrAndLept]]+lepton).M(), datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);	      
	      if(btagCSV[BLeptIndex[SelectedJetCombinationHadrAndLept]]>CSVbTag[1] || btagCSV[BHadrIndex[SelectedJetCombinationHadrAndLept]]>CSVbTag[1]){
		MSPlot["METAfterMCSVbTag"]->Fill(neutrinoOrigKins.Pt(), datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 
		MSPlot["LeptonPtMCSVbTag"]->Fill(lepton.Pt(), datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);
		MSPlot["Jet1PtMCSVbTag"]->Fill(selectedJets[0].Pt(), datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);
		MSPlot["MlbAfterTransMassCutAfterMCSVbTag"]->Fill((selectedJets[BLeptIndex[SelectedJetCombinationHadrAndLept]]+lepton).M(), datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
		if(dataSetName.find("Data") == 0) MSPlot["MlbAfterTransMassCutAfterMCSVbTagSFAdded"]->Fill((selectedJets[BLeptIndex[SelectedJetCombinationHadrAndLept]]+lepton).M(), datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D); //SF only for MC
		else MSPlot["MlbAfterTransMassCutAfterMCSVbTagSFAdded"]->Fill((selectedJets[BLeptIndex[SelectedJetCombinationHadrAndLept]]+lepton).M(), datasets[iDataSet], true, Luminosity*scaleFactor*0.97*lumiWeight3D);
	      }
	    }
	  }
	}
      }
      
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
      //----------------------------------
      //     Require some extra cuts:
      //----------------------------------
      bool eventSelected = false;

      if((dataSetName.find("Nom_TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("Nom_TTbarJets_SemiEl") ==0 && semiElectron == true) ) NumberEventsBeforeCuts++;
      if(dataSetName.find("Data") ==0) NumberDataEventsBeforeCuts++;
      
      if(CreateMSPlots == true)
	MSPlot["nPVBeforeCuts"]->Fill(nPrimaryVertices,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
      selecTableMacro.Fill(iDataSet,0,scaleFactor*lumiWeight3D);

      //Apply extra event selections:
      if(TransverseMass > 30 && DeltaRJetLepton > 0.3 && ((DeltaRMuonJet > 0.3 && semiMuon == true) || (DeltaRElectronJet > 0.3 && semiElectron == true)) && relativeLeptonIsolation <= 0.125){
	if((dataSetName.find("Nom_TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("Nom_TTbarJets_SemiEl") ==0 && semiElectron == true) ) NumberSelectedEvents++;
	if(dataSetName.find("Data") ==0) NumberSelectedDataEvents++;
	
	if(CreateMSPlots == true){
	  if(dataSetName.find("TopMass") != 0) MSPlot["TransverseMassAfterCut"]->Fill(TransverseMass, datasets[iDataSet], true, Luminosity);
	  if(dataSetName.find("TopMassMinus") != 0 && dataSetName.find("Nom_TTbarJets_") != 0) MSPlot["TransverseMassAfterCutTopMassPlus"]->Fill(TransverseMass, datasets[iDataSet], true, Luminosity);
	  if(dataSetName.find("TopMassPlus") != 0 && dataSetName.find("Nom_TTbarJets_") != 0) MSPlot["TransverseMassAfterCutTopMassMinus"]->Fill(TransverseMass, datasets[iDataSet], true, Luminosity);
	  MSPlot["RelIsoLepton"]->Fill(relativeLeptonIsolation, datasets[iDataSet], true, Luminosity);
	  MSPlot["nPVAfterTransverseMassCut"]->Fill(nPrimaryVertices,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	  MSPlot["Jet1Pt"]->Fill(selectedJets[0].Pt(), datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);    
	  MSPlot["Jet2Pt"]->Fill(selectedJets[1].Pt(), datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 
	  MSPlot["Jet3Pt"]->Fill(selectedJets[2].Pt(), datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 
	  MSPlot["Jet4Pt"]->Fill(selectedJets[3].Pt(), datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 	
	}

	if(CreateMSPlots == true && (btagCSV[BLeptIndex[SelectedJetCombinationHadrAndLept]]>CSVbTag[1] || btagCSV[BHadrIndex[SelectedJetCombinationHadrAndLept]]>CSVbTag[1])){
	  MSPlot["4thJetPtAfterMCSVbTag"]->Fill(selectedJets[3].Pt(), datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 	
	}

	//////////////////////////////////////////////////////////////////////
	//  No bTag discriminator combined with jet-combination selection   //
	//////////////////////////////////////////////////////////////////////	

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
	      float LowestChiSqHisto = LowestChiSq;
	      if(LowestChiSq > 150.){
		cout << " Large ChiSq Value for event " << iEvt << " , namely : " << LowestChiSq << endl;
		LowestChiSqHisto = 150.;
	      }
	      MSPlot["ChiSqKinFitCosThetaFoundChiSqOnlyNoBTag"]->Fill(LowestChiSqHisto,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	      MSPlot["ProbKinFitCosThetaFoundChiSqOnlyNoBTag"]->Fill(TMath::Prob(LowestChiSq,2),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	      
	      MSPlot["nPVChiSqOnlyNoBTag"]->Fill(nPrimaryVertices,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	      MSPlot["TransverseMassChiSqOnlyNoBTag"]->Fill(TransverseMass,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	      MSPlot["TransverseMassSelectedParticlesChiSqOnlyNoBTag"]->Fill(sqrt(2*(abs(leptonChangedHadrAndLept[SelectedJetCombinationHadrAndLept].Pt()))*abs(neutrinoChangedHadrAndLept[SelectedJetCombinationHadrAndLept].Pt())*(1-cos(leptonChangedHadrAndLept[SelectedJetCombinationHadrAndLept].DeltaPhi(neutrinoChangedHadrAndLept[SelectedJetCombinationHadrAndLept])))),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	      float CosThetaCalculationHisto = CosThetaCalculation;
	      if(CosThetaCalculation < -1. || CosThetaCalculation > 1.){
		cout << " Wrong Cos theta value for event " << iEvt << " , namely : " << CosThetaCalculation << endl;
		CosThetaCalculationHisto = 2.;
	      }
	      MSPlot["CosThetaChiSqOnlyNoBTag"]->Fill(CosThetaCalculationHisto,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	      
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
	    if(CreateMSPlots == true) MSPlot["bTagCSVbLeptDistribution"]->Fill(btagCSV[BLeptIndex[SelectedJetCombinationHadrAndLept]],datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	    if(btagCSV[BLeptIndex[SelectedJetCombinationHadrAndLept]]>CSVbTag[1]){
	      if(dataSetName.find("Data") != 0 && CreateMSPlots == true) MSPlot["bTagCSVbLeptDistributionAfterCut"]->Fill(btagCSV[BLeptIndex[SelectedJetCombinationHadrAndLept]],datasets[iDataSet], true, Luminosity*scaleFactor*0.97*lumiWeight3D);
	      else if(CreateMSPlots == true) MSPlot["bTagCSVbLeptDistributionAfterCut"]->Fill(btagCSV[BLeptIndex[SelectedJetCombinationHadrAndLept]],datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	      
	      if(CreateMSPlots == true){
		MSPlot["JetCombCosThetaFoundChiSqOnlySpecificBTag"]->Fill(SelectedJetCombinationHadrAndLept,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
		float LowestChiSqHisto = LowestChiSq;
		if(LowestChiSq > 150.){
		  cout << " Large ChiSq Value for event " << iEvt << " , namely : " << LowestChiSq << endl;
		  LowestChiSqHisto = 150.;
		}
		MSPlot["ChiSqKinFitCosThetaFoundChiSqOnlySpecificBTag"]->Fill(LowestChiSqHisto,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
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
	//   if((dataSetName.find("Nom_TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("Nom_TTbarJets_SemiEl") ==0 && semiElectron == true) ){ 
	//     UsedDistributionValue[helicityNumbers]=(HelicityFraction[helicityNumbers][0]*6*(1-CosThetaGenerator*CosThetaGenerator) + (1-CosThetaGenerator)*(1-CosThetaGenerator)*3*HelicityFraction[helicityNumbers][2] + HelicityFraction[helicityNumbers][1]*3*(1+CosThetaGenerator)*(1+CosThetaGenerator))/8;
	//     HelicityWeight[helicityNumbers]=UsedDistributionValue[helicityNumbers]/TheoreticalDistributionValue;
	//   }
	// }

	//-------------------------------------------------------------
	//Obtain jet combination for the different b-tag constraints:
	//------------------------------------------------------------
	float OriginalScaleFactor = scaleFactor;

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
		      //if(UsedTCHE[SumBTag]==TCHEbTagLoop && UsedTCHP[SumBTag]==TCHPbTagLoop && UsedSSVHE[SumBTag]==SSVHEbTagLoop && UsedSSVHP[SumBTag]==SSVHPbTagLoop && UsedCSV[SumBTag]==CSVbTagLoop && UsedJP[SumBTag]==JPbTagLoop && UsedJBP[SumBTag]==JBPbTagLoop){
			
			int bTagLoop;
			if(ConsideredBTagger == 0) bTagLoop = TCHEbTagLoop;		   
			else if(ConsideredBTagger == 1) bTagLoop = TCHPbTagLoop;
			else if(ConsideredBTagger == 2) bTagLoop = SSVHEbTagLoop;
			else if(ConsideredBTagger == 3) bTagLoop = SSVHPbTagLoop;
			else if(ConsideredBTagger == 4) bTagLoop = CSVbTagLoop;
			else if(ConsideredBTagger == 5) bTagLoop = JPbTagLoop;
			else if(ConsideredBTagger == 6) bTagLoop = JBPbTagLoop;

			if(dataSetName.find("Data")!=0){
			  if(bTagLoop == 1) scaleFactor = OriginalScaleFactor;
			  else if(bTagLoop == 2 || bTagLoop == 3 || bTagLoop == 4) scaleFactor = OriginalScaleFactor*LooseSF[ConsideredBTagger];
			  else if(bTagLoop == 5) scaleFactor = OriginalScaleFactor*LooseSF[ConsideredBTagger]*LooseSF[ConsideredBTagger];
			  else if(bTagLoop == 6 || bTagLoop == 7 || bTagLoop == 8) scaleFactor = OriginalScaleFactor*MediumSF[ConsideredBTagger];
			  else if(bTagLoop == 9) scaleFactor = OriginalScaleFactor*MediumSF[ConsideredBTagger]*MediumSF[ConsideredBTagger];
			  else if(bTagLoop == 10 || bTagLoop == 11 || bTagLoop == 12) scaleFactor = OriginalScaleFactor*TightSF[ConsideredBTagger];
			  else if(bTagLoop == 13) scaleFactor = OriginalScaleFactor*TightSF[ConsideredBTagger]*TightSF[ConsideredBTagger];
			}
			
			int JetCombHadr = bTagJetSelection.HighestProbSelection(bTagLoop,ConsideredBTagger,KinFitProbHadr,btagTCHE,btagTCHP,btagSSVHE,btagSSVHP,btagCSV,btagJP,btagJBP);
			int JetCombHadrAndLeptWOnly =bTagJetSelection.HighestProbSelection(bTagLoop,ConsideredBTagger,KinFitProbHadrAndLeptWOnly,btagTCHE,btagTCHP,btagSSVHE,btagSSVHP,btagCSV,btagJP,btagJBP);
			int JetCombHadrAndLept = bTagJetSelection.HighestProbSelection(bTagLoop,ConsideredBTagger,KinFitProbHadrAndLept,btagTCHE,btagTCHP,btagSSVHE,btagSSVHP,btagCSV,btagJP,btagJBP);
			  
			if(bTagLoop ==1){
			  if(JetCombHadrAndLept != 999 && CreateMSPlots == true){
			    MSPlot["HadrTopMassSelJetsNobTag"]->Fill((selectedJets[BHadrIndex[JetCombHadrAndLept]]+selectedJets[Quark1Index[JetCombHadrAndLept]]+selectedJets[Quark2Index[JetCombHadrAndLept]]).M(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
			    MSPlot["bTagCSVbLeptDistributionClass"]->Fill(btagCSV[BLeptIndex[JetCombHadrAndLept]],datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
			  }
			  
			  if(CreateMSPlots == true){
			    MSPlot["JetCombination"]->Fill(JetCombHadrAndLept,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
			    histo1D["JetCombHadr"]->Fill(JetCombHadr);
			    histo1D["JetCombHadrAndLeptWOnly"]->Fill(JetCombHadrAndLeptWOnly);
			    histo1D["JetCombHadrAndLept"]->Fill(JetCombHadrAndLept);
			  }
			}
			if(ConsideredBTagger == 4 && bTagLoop == 7 && JetCombHadrAndLept != 999 && CreateMSPlots == true){
			  //cout << " Selected jet combination with class : " << JetCombHadrAndLept << " has bTag value of : " << btagCSV[BLeptIndex[JetCombHadrAndLept]] << endl;
			  MSPlot["HadrTopMassSelJetsOnebTag"]->Fill((selectedJets[BHadrIndex[JetCombHadrAndLept]]+selectedJets[Quark1Index[JetCombHadrAndLept]]+selectedJets[Quark2Index[JetCombHadrAndLept]]).M(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
			  MSPlot["bTagCSVbLeptDistributionClassAfterCut"]->Fill(btagCSV[BLeptIndex[JetCombHadrAndLept]],datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
			  if(btagCSV[BLeptIndex[JetCombHadrAndLept]] < 0.679){
			    cout << " KinFit Probability = " << KinFitProbHadrAndLept[JetCombHadrAndLept] << endl;
			    cout << " Chosen combination = " << JetCombHadrAndLept << endl;
			  }
			}

			if(ConsideredBTagger == 4 && bTagLoop == 9 && JetCombHadrAndLept != 999 && CreateMSPlots == true){
			  MSPlot["HadrTopMassSelJetsTwobTag"]->Fill((selectedJets[BHadrIndex[JetCombHadrAndLept]]+selectedJets[Quark1Index[JetCombHadrAndLept]]+selectedJets[Quark2Index[JetCombHadrAndLept]]).M(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
			}
			
			// Only hadronic KinFit configuration:
			if(JetCombHadr!=999 && HadronicOnly == true){ 				  
			  
			  float CosThetaCalcHadr = bTagCosThetaCalculation.CalcOrigKins(BLeptIndex[JetCombHadr], BHadrIndex[JetCombHadr],  lepton, selectedJets, 80.4, MassTop);
			  if(CosThetaCalcHadr != 999){			    
			    TLorentzVector neutrinoHadr = bTagCosThetaCalculation.GetNeutrino();

			    if((lepton+neutrinoHadr).E() < 0. ){
			      cout << " Processing event : " << iEvt << " with cos theta = " << CosThetaCalcHadr << " (Chi Sq and BTag case -- Hadr)  --   Used combination is : " << JetCombHadr << " with chi squared value = " << ChiSqKinFitHadr[JetCombHadr] << " and KinFit probability : " << KinFitProbHadr[JetCombHadr] << endl;
			    }
			    
			    if((dataSetName.find("Nom_TTbarJets_SemiMu") == 0 && semiMuon == true) || (dataSetName.find("Nom_TTbarJets_SemiEl") == 0 && semiElectron == true)){
			      histo1D["LeptWMassHadr"]->Fill((neutrinoHadr+lepton).M());
			      histo1D["LeptTopMassHadr"]->Fill((neutrinoHadr+lepton+selectedJets[BLeptIndex[JetCombHadr]]).M());
			    }		      
			    
			    NumberRemainingEventsKinFitHadr[SumBTag][iDataSet]++;
			    CosThetaValuesKinFitHadr[SumBTag].push_back(CosThetaCalcHadr);
			    EventCorrectionWeightKinFitHadr[SumBTag].push_back(scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor()); 
			    if((dataSetName.find("Nom_TTbarJets_SemiMu") == 0 && semiMuon == true) || (dataSetName.find("Nom_TTbarJets_SemiEl") ==0 && semiElectron == true) ){
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
			    EventCorrectionWeightKinFitHadrAndLeptWOnly[SumBTag].push_back(scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor());

			    if((dataSetName.find("Nom_TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("Nom_TTbarJets_SemiEl") ==0 && semiElectron == true) ){
			      if(BLeptIndex[JetCombHadrAndLeptWOnly] == CorrectBLept){//Count the number of events with correctly reconstructed leptonic b-jet
				NumberBLeptCorrectEventsKinFitHadrAndLeptWOnly[SumBTag][iDataSet]++;
			      }
			      CosThGenKinFitHadrAndLeptWOnly[SumBTag].push_back(CosThetaGenerator);
			    }//End of signal sample
			  }
			}
			
			//Hadronic and leptonic KinFit configuration:
			if(bTagLoop == 1) NumberEventsChiSqAndBTag[iDataSet]++;
			if(bTagLoop == 1 && CreateMSPlots == true) MSPlot["nPVAll"]->Fill(nPrimaryVertices,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
			if(JetCombHadrAndLept!=999 && HadronicAndLeptonic == true){ 
			  if(bTagLoop == 1 && CreateMSPlots == true) MSPlot["nPVExistingJetComb"]->Fill(nPrimaryVertices,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);

			  if(bTagLoop == 1) NumberEventsChiSqAndBTagJetCombFound[iDataSet]++;
			  float CosThetaCalcHadrAndLept=bTagCosThetaCalculation.Calculation(leptonChangedHadrAndLept[JetCombHadrAndLept],neutrinoChangedHadrAndLept[JetCombHadrAndLept],leptBChangedHadrAndLept[JetCombHadrAndLept]);
			  if(CosThetaCalcHadrAndLept != 999){
			    if(bTagLoop == 1 && CreateMSPlots == true) MSPlot["nPVExistingCos"]->Fill(nPrimaryVertices,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
			    if(ConsideredBTagger == 4 && bTagLoop == 7 && CreateMSPlots == true) MSPlot["nPVExistingCosbTag"]->Fill(nPrimaryVertices,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);

			    if(bTagLoop == 1){
			      if(dataSetName.find("Data") == 0) histo1D["CosThetaNoBTagData"]->Fill(CosThetaCalcHadrAndLept,scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor() );
			      if(dataSetName.find("Nom_TTbarJets_SemiMuon")==0) histo1D["CosThetaNoBTagSemiMu"]->Fill(CosThetaCalcHadrAndLept,scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor() );
			      if(dataSetName.find("Nom_TTbarJets_SemiEl")==0) histo1D["CosThetaNoBTagSemiEl"]->Fill(CosThetaCalcHadrAndLept,scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor() );
			      if(dataSetName.find("Nom_TTbarJets_Other")==0) histo1D["CosThetaNoBTagOther"]->Fill(CosThetaCalcHadrAndLept,scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor() );
			      if(dataSetName.find("Nom_WJets")==0) histo1D["CosThetaNoBTagWJets"]->Fill(CosThetaCalcHadrAndLept,scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor() );
			      if(dataSetName.find("Nom_ZJets")==0) histo1D["CosThetaNoBTagZJets"]->Fill(CosThetaCalcHadrAndLept,scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor() );
			      if(dataSetName.find("Nom_ST_SingleTop_tChannel_tbar")==0) histo1D["CosThetaNoBTagSTtChtBar"]->Fill(CosThetaCalcHadrAndLept,scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor() );
			      if(dataSetName.find("Nom_ST_SingleTop_tChannel_t")==0) histo1D["CosThetaNoBTagSTtCht"]->Fill(CosThetaCalcHadrAndLept,scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor() );
			      if(dataSetName.find("Nom_ST_SingleTop_tWChannel_tbar")==0) histo1D["CosThetaNoBTagSTtWChtBar"]->Fill(CosThetaCalcHadrAndLept,scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor() );
			      if(dataSetName.find("Nom_ST_SingleTop_tWChannel_t")==0) histo1D["CosThetaNoBTagSTtWCht"]->Fill(CosThetaCalcHadrAndLept,scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor() );
			    }

			    if(bTagLoop == 7 && ConsideredBTagger == 4){
			      if(dataSetName.find("Data") == 0) histo1D["CosThetaOneBTagData"]->Fill(CosThetaCalcHadrAndLept,scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor() );
			      if(dataSetName.find("Nom_TTbarJets_SemiMuon")==0) histo1D["CosThetaOneBTagSemiMu"]->Fill(CosThetaCalcHadrAndLept,scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor() );
			      if(dataSetName.find("Nom_TTbarJets_SemiEl")==0) histo1D["CosThetaOneBTagSemiEl"]->Fill(CosThetaCalcHadrAndLept,scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor() );
			      if(dataSetName.find("Nom_TTbarJets_Other")==0) histo1D["CosThetaOneBTagOther"]->Fill(CosThetaCalcHadrAndLept,scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor() );
			      if(dataSetName.find("Nom_WJets")==0) histo1D["CosThetaOneBTagWJets"]->Fill(CosThetaCalcHadrAndLept,scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor() );
			      if(dataSetName.find("Nom_ZJets")==0) histo1D["CosThetaOneBTagZJets"]->Fill(CosThetaCalcHadrAndLept,scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor() );
			      if(dataSetName.find("Nom_ST_SingleTop_tChannel_tbar")==0) histo1D["CosThetaOneBTagSTtChtBar"]->Fill(CosThetaCalcHadrAndLept,scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor() );
			      if(dataSetName.find("Nom_ST_SingleTop_tChannel_t")==0) histo1D["CosThetaOneBTagSTtCht"]->Fill(CosThetaCalcHadrAndLept,scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor() );
			      if(dataSetName.find("Nom_ST_SingleTop_tWChannel_tbar")==0) histo1D["CosThetaOneBTagSTtWChtBar"]->Fill(CosThetaCalcHadrAndLept,scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor() );
			      if(dataSetName.find("Nom_ST_SingleTop_tWChannel_t")==0) histo1D["CosThetaOneBTagSTtWCht"]->Fill(CosThetaCalcHadrAndLept,scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor() );
			    }
		
			    if(bTagLoop == 9 && ConsideredBTagger == 4){
			      if(dataSetName.find("Data") == 0) histo1D["CosThetaTwoBTagData"]->Fill(CosThetaCalcHadrAndLept,scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor() );
			      if(dataSetName.find("Nom_TTbarJets_SemiMuon")==0) histo1D["CosThetaTwoBTagSemiMu"]->Fill(CosThetaCalcHadrAndLept,scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor() );
			      if(dataSetName.find("Nom_TTbarJets_SemiEl")==0) histo1D["CosThetaTwoBTagSemiEl"]->Fill(CosThetaCalcHadrAndLept,scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor() );
			      if(dataSetName.find("Nom_TTbarJets_Other")==0) histo1D["CosThetaTwoBTagOther"]->Fill(CosThetaCalcHadrAndLept,scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor() );
			      if(dataSetName.find("Nom_WJets")==0) histo1D["CosThetaTwoBTagWJets"]->Fill(CosThetaCalcHadrAndLept,scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor() );
			      if(dataSetName.find("Nom_ZJets")==0) histo1D["CosThetaTwoBTagZJets"]->Fill(CosThetaCalcHadrAndLept,scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor() );
			      if(dataSetName.find("Nom_ST_SingleTop_tChannel_tbar")==0) histo1D["CosThetaTwoBTagSTtChtBar"]->Fill(CosThetaCalcHadrAndLept,scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor() );
			      if(dataSetName.find("Nom_ST_SingleTop_tChannel_t")==0) histo1D["CosThetaTwoBTagSTtCht"]->Fill(CosThetaCalcHadrAndLept,scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor() );
			      if(dataSetName.find("Nom_ST_SingleTop_tWChannel_tbar")==0) histo1D["CosThetaTwoBTagSTtWChtBar"]->Fill(CosThetaCalcHadrAndLept,scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor() );
			      if(dataSetName.find("Nom_ST_SingleTop_tWChannel_t")==0) histo1D["CosThetaTwoBTagSTtWCht"]->Fill(CosThetaCalcHadrAndLept,scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor() );
			    }

// 			    //Create MSPlots for selecting optimal electron b-tagger:
// 			    if(dataSetName.find("Nom_TTbarJets_SemiEl") == 0 && semiElectron == true && CreateMSPlots == true && ConsideredBTagger == 4 && bTagLoop == 9)
// 			      scaleFactor = scaleFactor*(1.130 + 1.115)/2.;
			    
// 			    if(ConsideredBTagger == 4 && bTagLoop == 7 && semiMuon == true && (dataSetName.find("Nom_TTbarJets_SemiMuon") == 0 || dataSetName.find("McData_TTbarJets_SemiMuon") == 0))
// 			      scaleFactor = scaleFactor*1.0585;//(1.121+1.108)/2.;
// 			    else if(ConsideredBTagger == 4 && bTagLoop == 7 && semiMuon == true && dataSetName.find("Nom_TTbarJets_SemiMuon") != 0 && dataSetName.find("McData_TTbarJets_SemiMuon") != 0)
// 			      scaleFactor = scaleFactor*1.1885;//(1.121+1.108)/2.;
// 			    if(ConsideredBTagger == 4 && bTagLoop == 7 && semiMuon == true && dataSetName.find("Nom_TTbarJets_SemiMuon") == 0 )
// 			      scaleFactor = scaleFactor*1.0585;//(1.121+1.108)/2.;
// 			    else if(ConsideredBTagger == 4 && bTagLoop == 7 && semiMuon == true && dataSetName.find("Nom_TTbarJets_SemiMuon") != 0 && dataSetName.find("Data") != 0 && dataSetName.find("McDta") != 0)
// 			      scaleFactor = scaleFactor*1.1885;//(1.121+1.108)/2.;
			    
// 			    if(dataSetName.find("Nom_TTbarJets_SemiEl") == 0)	
// 			      scaleFactor = scaleFactor*1.1;
// 			    else if(dataSetName.find("Nom_TTbarJets_Other")==0)
// 			      scaleFactor = scaleFactor*0.7121;			  
 			    //if(dataSetName.find("Nom_TTbarJets_SemiMu") == 0)	
			    //scaleFactor = scaleFactor*0.9;
 			    //if(dataSetName.find("Nom_WJets")==0)
			    //scaleFactor = scaleFactor*1.83;			  
			    
			    if((dataSetName.find("Nom_TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("Nom_TTbarJets_SemiEl") ==0 && semiElectron == true) ){	
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
			    if(bTagLoop == 1) NumberEventsChiSqAndBTagCosThetaFound[iDataSet]++;
			    NumberRemainingEventsKinFitHadrAndLept[SumBTag][iDataSet]++;
			    CosThetaValuesKinFitHadrAndLept[SumBTag].push_back(CosThetaCalcHadrAndLept);
			    //MlbValuesKinFitHadrAndLept[SumBTag].push_back((leptonChangedHadrAndLept[JetCombHadrAndLept]+leptBChangedHadrAndLept[JetCombHadrAndLept]).M());
			    MlbValuesKinFitHadrAndLept[SumBTag].push_back(sqrt(2*(abs(leptonChangedHadrAndLept[JetCombHadrAndLept].Pt()))*abs(neutrinoChangedHadrAndLept[JetCombHadrAndLept].Pt())*(1-cos(leptonChangedHadrAndLept[JetCombHadrAndLept].DeltaPhi(neutrinoChangedHadrAndLept[JetCombHadrAndLept])))));
			    EventCorrectionWeightKinFitHadrAndLept[SumBTag].push_back(scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor());
			    if(bTagLoop == 1 && CreateMSPlots == true) MSPlot["MlbNoBTag"]->Fill((leptonChangedHadrAndLept[JetCombHadrAndLept]+leptBChangedHadrAndLept[JetCombHadrAndLept]).M(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
			    
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
			    if(ConsideredBTagger == 4 && bTagLoop == 9 && semiMuon == true && CreateMSPlots == true){
			      MSPlot["CosThetaSpecificDoubleBTag35Bins"]->Fill(CosThetaCalcHadrAndLept,datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);
			      MSPlot["CosThetaSpecificDoubleBTag25Bins"]->Fill(CosThetaCalcHadrAndLept,datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);
			    }
			    if(ConsideredBTagger == 4 && bTagLoop == 7 && semiElectron == true && CreateMSPlots == true){
			      MSPlot["CosThetaSpecificDoubleBTag35Bins"]->Fill(CosThetaCalcHadrAndLept,datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);
			      MSPlot["CosThetaSpecificDoubleBTag25Bins"]->Fill(CosThetaCalcHadrAndLept,datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);
			    }
			    
			    //  MSPlots for specific b-tag case:  
			    int OptimalTagger, OptimalbTagLoop;
			    if(semiMuon == true){
			      OptimalTagger = 4;
			      OptimalbTagLoop = 7;
			    }
			    else if(semiElectron == true){
			      OptimalTagger = 4;
			      OptimalbTagLoop = 9;
			    }
			    if(ConsideredBTagger == OptimalTagger && bTagLoop == OptimalbTagLoop){
			      NumberEventsChiSqAndBTagSpecificBTag[iDataSet]++;
				
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
				MSPlot["CosThetaSpecificBTag35Bins"]->Fill(CosThetaCalcHadrAndLept,datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);
				MSPlot["CosThetaSpecificBTag25Bins"]->Fill(CosThetaCalcHadrAndLept,datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);
				
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
			    
			    if((dataSetName.find("Nom_TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("Nom_TTbarJets_SemiEl") ==0 && semiElectron == true) ){
			      if(BLeptIndex[JetCombHadrAndLept] == CorrectBLept){//Count the number of events with correctly reconstructed leptonic b-jet
				NumberBLeptCorrectEventsKinFitHadrAndLept[SumBTag][iDataSet]++;
			      }
			      CosThGenKinFitHadrAndLept[SumBTag].push_back(CosThetaGenerator);
			    }//End of signal sample
			  }
			}
			//		      }
// 		      else{
// 			cout << " Looking at wrong bTagging combination in calculation loops " << endl;
// 			cout << " Looking at : TCHE = "<<TCHEbTagLoop <<" TCHP = "<<TCHPbTagLoop<<" SSVHE = "<<SSVHEbTagLoop<<" SSVHP = "<<SSVHPbTagLoop<<" CSV = "<<CSVbTagLoop<< " JP = " << JPbTagLoop << " JBP = " << JBPbTagLoop << endl;
// 		      }
		  
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
    //std::cout << "            - Hadronic KinFit only                : " << CosThetaValuesKinFitHadr[0].size() << endl;
    //std::cout << "            - Hadronic KinFit and leptonic W only : " << CosThetaValuesKinFitHadrAndLeptWOnly[0].size() << endl;
    std::cout << "            - Hadronic KinFit and leptonic W      : " << CosThetaValuesKinFitHadrAndLept[0].size() << endl;
    std::cout << "            - Hadronic KinFit and leptonic W (M CSV lept)     : " << CosThetaValuesKinFitHadrAndLept[13+13+13+13+6+1+1].size() << endl;
    cout << " Considered configuration : " << PresentationOutput[13+13+13+13+6+1+1] << endl;
    std::cout << " °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°° " << endl;
    
    ////////////////////////////////
    //    Execute MinuitFitter:   //
    ////////////////////////////////

    int dataSetSize = datasets.size();
    MinuitFitterPerforming minuitFitterHadronic = MinuitFitterPerforming(HadronicOnly, PerformMinuit, "Had", CosThetaBinNumber, ndimen, iDataSet, &CosThetaValuesKinFitHadr[TotalNumberbTags], &CosThGenKinFitHadr[TotalNumberbTags], &EventCorrectionWeightKinFitHadr[TotalNumberbTags], dataSetName, semiMuon, semiElectron, SignalOnly, DataResults, JESResults, JERResults, WSystResults, TTScalingResults, TTMatchingResults, PUResults, UnclusEnergyResults, TopMassResults, TriggEvtSelResults, bTagResults, dataSetSize);

    MinuitFitterPerforming minuitFitterHadronicAndLeptonicWOnly = MinuitFitterPerforming(HadronicAndLeptonicWOnly, PerformMinuit, "HadrAndLeptWOnly", CosThetaBinNumber, ndimen, iDataSet, &CosThetaValuesKinFitHadrAndLeptWOnly[TotalNumberbTags], &CosThGenKinFitHadrAndLeptWOnly[TotalNumberbTags], &EventCorrectionWeightKinFitHadrAndLeptWOnly[TotalNumberbTags], dataSetName, semiMuon, semiElectron, SignalOnly, DataResults, JESResults, JERResults, WSystResults, TTScalingResults, TTMatchingResults, PUResults, UnclusEnergyResults, TopMassResults, TriggEvtSelResults, bTagResults, dataSetSize);

    MinuitFitterPerforming minuitFitterHadronicAndLeptonic = MinuitFitterPerforming(HadronicAndLeptonic, PerformMinuit, "HadrAndLept", CosThetaBinNumber, ndimen, iDataSet, &CosThetaValuesKinFitHadrAndLept[TotalNumberbTags], &CosThGenKinFitHadrAndLept[TotalNumberbTags], &EventCorrectionWeightKinFitHadrAndLept[TotalNumberbTags], dataSetName, semiMuon, semiElectron, SignalOnly, DataResults, JESResults, JERResults, WSystResults, TTScalingResults, TTMatchingResults, PUResults, UnclusEnergyResults, TopMassResults, TriggEvtSelResults, bTagResults, dataSetSize);
            
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
