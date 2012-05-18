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

const int CosThetaBinNumber = 15;       
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
  
  //Probability cut value:
  float KinFitCut = 0.;

  //PU systematics:
  float PUSyst = 0;
  string systematic = "";
  if(PUSyst == 1) systematics = "PUPlus";
  else if(PUSyst == -1) systematics = "PUMinus";

  string decayChannel;
  bool semiMuon = true;
  bool semiElectron = false;
  if(semiMuon == true){decayChannel = "SemiMu";}
  else if(semiElectron == true){decayChannel = "SemiEl";}

  string dataSet;
  bool fullDataSet = true;
  if(fullDataSet == true) dataSet = "Full2011_";
  else dataSet ="";

  string UsedTrigger;
  bool IsoMu172024Trigger = false;
  bool TriCentralJet30Trigger = true;
  if(IsoMu172024Trigger == true){
    UsedTrigger = "IsoMu172024Trigger";
    if(fullDataSet == true){ 
      Luminosity = 4568.6810;
      UsedTrigger = UsedTrigger+"_Full2011";
    }
    else Luminosity = 2141.961; 
  }
  else if(TriCentralJet30Trigger == true){
    UsedTrigger = "TriCentralJet30Trigger";
    if(semiMuon == true){
      if(fullDataSet == true) Luminosity = 4656.959;
      else Luminosity = 2145.959; 
    }
    else if(semiElectron == true){
      if(fullDataSet == true) Luminosity = 4665.744;
      else Luminosity = 2161.744;
    }
  }
  else if(TriCentralJet30Trigger == false && IsoMu172024Trigger == false){
    UsedTrigger = "NoTrigger";
  }

  cout << "Executing the W Helicities analysis for an integrated luminosity of " << Luminosity << " pb^-1" << endl;

  //Booleans to load in different root files
  bool SignalOnly = true;
  bool DataResults = false;
  bool JESMinusResults = false;
  bool JESPlusResults = false;
  bool WSystResults = false;
  bool WSystPositive = false;
  string JESDirection;
  if(JESMinusResults == true)
    JESDirection = "JESMinus";
  else if(JESPlusResults == true)
    JESDirection = "JESPlus";
  
  cout << " Obtaining results for : Data = " << DataResults << " , JESMinus = " << JESMinusResults << " , JESPlus = " << JESPlusResults << " , WSyst = " << WSystResults << " for positive scaling = " << WSystPositive << endl;

  //Booleans for event selection cuts
  bool SSVHEbTag = false;
  bool MTCut = false;
  bool MuonPtCut = false;  //Require a muon Pt larger than 27 (to avoid turn-over of IsoMu17/20/24 triggers)  -- Ciemat uses 25
  string AppliedCuts = "";
  if(SSVHEbTag == true)
    AppliedCuts = AppliedCuts+"_MSSVHEbTag";
  if(MTCut == true)
    AppliedCuts = AppliedCuts+"_MTCut";
  if(MuonPtCut == true)
    AppliedCuts = AppliedCuts+"_MuonPtCut";

  //Booleans for KinFit options
  bool UseChangedKinematics = true;  //Boolean to differentiate between changed kinematics and original kinematics of fitted particles
  bool LeptTopMassConstrained = true;  //Boolean to differentiate between lept top mass constrained to world average and lept tom mass left free in KinFit
  bool LeptonicFit = true; //Boolean to differentiate between KinFit on hadronic side only and KinFit on hadronic+leptonic side
  
  //1) Nominal samples:
  if(SignalOnly == false){
    if(DataResults == true){
      if(IsoMu172024Trigger == true)
	inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_"+"Data_"+decayChannel+".root").c_str());
      else if(TriCentralJet30Trigger == true)
	inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_"+dataSet+"Data_"+decayChannel+".root").c_str());
      nameDataSet.push_back("Data");
    }
    inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ST_SingleTop_tChannel_tbar_"+decayChannel+".root").c_str());
    nameDataSet.push_back("ST_SingleTop_tChannel_tbar");
    inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ST_SingleTop_tChannel_t_"+decayChannel+".root").c_str());
    nameDataSet.push_back("ST_SingleTop_tChannel_t");
    inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ST_SingleTop_tWChannel_tbar_"+decayChannel+".root").c_str());
    nameDataSet.push_back("ST_SingleTop_tWChannel_tbar");
    inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ST_SingleTop_tWChannel_t_"+decayChannel+".root").c_str());
    nameDataSet.push_back("ST_SingleTop_tWChannel_t");
    inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_TTbarJets_Other_"+decayChannel+".root").c_str());
    nameDataSet.push_back("TTbarJets_Other");
    if(semiMuon == true){
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_TTbarJets_SemiEl_"+decayChannel+".root").c_str());  //In muon channel case SemiEl is considered as background
      nameDataSet.push_back("TTbarJets_SemiEl");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_QCD_Mu15_SemiMu.root").c_str());
      nameDataSet.push_back("QCD_Mu15");
    }
    else if(semiElectron == true){
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_TTbarJets_SemiMuon_"+decayChannel+".root").c_str());  //In electron channel case SemiMu is considered as background
      nameDataSet.push_back("TTbarJets_SemiMuon");
      //inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_QCD_Pt-20to30_BCtoE_SemiEl.root").c_str());  //Use UsedTrigger string to make sure that sample is only used in correct case
      // nameDataSet.push_back("QCD_Pt-20to30_BCtoE");
      // inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_QCD_Pt-20to30_EMEnriched_SemiEl.root").c_str());
      // nameDataSet.push_back("QCD_Pt-20to30_EMEnriched");
      // inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_QCD_Pt-30to80_BCtoE_SemiEl.root").c_str());
      // nameDataSet.push_back("QCD_Pt-30to80_BCtoE");
      // inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_QCD_Pt-30to80_EMEnriched_SemiEl.root").c_str());
      // nameDataSet.push_back("QCD_Pt-30to80_EMEnriched");
      // inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_QCD_Pt-80to170_BCtoE_SemiEl.root").c_str());
      // nameDataSet.push_back("QCD_Pt-80to170_BCtoE");
      // inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_QCD_Pt-80to170_EMEnriched_SemiEl.root").c_str());
      // nameDataSet.push_back("QCD_Pt-80to170_EMEnriched");
    }
    inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_WJets_"+decayChannel+".root").c_str());
    nameDataSet.push_back("WJets_Nominal");
    inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_ZJets_"+decayChannel+".root").c_str());
    nameDataSet.push_back("ZJets");  

    //2) JES Plus/Min samples:
    if(JESMinusResults == true || JESPlusResults == true){
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_"+JESDirection+"_1Sig_ST_SingleTop_tChannel_tbar_"+decayChannel+".root").c_str());
      nameDataSet.push_back("JES_ST_SingleTop_tChannel_tbar");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_"+JESDirection+"_1Sig_ST_SingleTop_tChannel_t_"+decayChannel+".root").c_str());
      nameDataSet.push_back("JES_ST_SingleTop_tChannel_t");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_"+JESDirection+"_1Sig_ST_SingleTop_tWChannel_tbar_"+decayChannel+".root").c_str());
      nameDataSet.push_back("JES_ST_SingleTop_tWChannel_tbar");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_"+JESDirection+"_1Sig_ST_SingleTop_tWChannel_t_"+decayChannel+".root").c_str());
      nameDataSet.push_back("JES_ST_SingleTop_tWChannel_t");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_"+JESDirection+"_1Sig_TTbarJets_Other_"+decayChannel+".root").c_str());
      nameDataSet.push_back("JES_TTbarJets_Other");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_"+JESDirection+"_1Sig_TTbarJets_SemiEl_"+decayChannel+".root").c_str());
      nameDataSet.push_back("JES_TTbarJets_SemiEl");
      if(semiElectron == true){
	inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_"+JESDirection+"_1Sig_QCD_Pt-20to30_BCtoE_SemiEl.root").c_str());  //Use UsedTrigger string to make sure that sample is only used in correct case
	nameDataSet.push_back("QCD_Pt-20to30_BCtoE");
	inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_"+JESDirection+"_1Sig_QCD_Pt-20to30_EMEnriched_SemiEl.root").c_str());
	nameDataSet.push_back("QCD_Pt-20to30_EMEnriched");
	inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_"+JESDirection+"_1Sig_QCD_Pt-30to80_BCtoE_SemiEl.root").c_str());
	nameDataSet.push_back("QCD_Pt-30to80_BCtoE");
	inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_"+JESDirection+"_1Sig_QCD_Pt-30to80_EMEnriched_SemiEl.root").c_str());
	nameDataSet.push_back("QCD_Pt-30to80_EMEnriched");
	inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_"+JESDirection+"_1Sig_QCD_Pt-80to170_BCtoE_SemiEl.root").c_str());
	nameDataSet.push_back("QCD_Pt-80to170_BCtoE");
	inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_"+JESDirection+"_1Sig_QCD_Pt-80to170_EMEnriched_SemiEl.root").c_str());
	nameDataSet.push_back("QCD_Pt-80to170_EMEnriched");
      }
      else if(semiMuon ==true){
	inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_"+JESDirection+"_1Sig_QCD_Mu15_"+decayChannel+".root").c_str());
	nameDataSet.push_back("JES_QCD_Mu15");
      }
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_"+JESDirection+"_1Sig_WJets_"+decayChannel+".root").c_str());    
      nameDataSet.push_back("JES_WJets");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_"+JESDirection+"_1Sig_ZJets_"+decayChannel+".root").c_str());
      nameDataSet.push_back("JES_ZJets");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_"+JESDirection+"_1Sig_TTbarJets_SemiMuon_"+decayChannel+".root").c_str());
      nameDataSet.push_back("JES_TTbarJets_SemiMuon");
    }
    
    //3) WJets systematics:  
    if(WSystResults == true){
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_WJets_"+decayChannel+".root").c_str());  
      nameDataSet.push_back("Syst_WJets");
    }
  }//End of signalOnly = false loop
  
  //TTbarJets_SemiMuon sample should always be put as latest sample to avoid crash of TMinuitMinimizer !!
  if(semiMuon == true){
    //inputWTree.push_back(("WTree/KinFit_WTree_SHORTTEST_DRLep03_TriCentralJet30Trigger_TTbarJets_SemiMuon_SemiMu.root").c_str());
    cout << " Using file : " << ("WTree/KinFit_WTree_SHORTTEST_DRLep03_"+UsedTrigger+"_TTbarJets_SemiMuon_"+decayChannel+".root").c_str() << endl;
    inputWTree.push_back(("WTree/KinFit_WTree_SHORTTEST_DRLep03_"+UsedTrigger+"_TTbarJets_SemiMuon_"+decayChannel+".root").c_str());
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
  PresentationTexTitle =("BTagPerformanceStudy/"+UsedTrigger+"_"+dataSet+decayChannel+AppliedCuts+".tex").c_str();
  if(DataResults == true) PresentationTexTitle =("BTagPerformanceStudy/"+UsedTrigger+"_"+dataSet+decayChannel+AppliedCuts+".tex").c_str();
  if(JESMinusResults == true) PresentationTexTitle =("BTagPerformanceStudy/"+UsedTrigger+"_"+dataSet+decayChannel+AppliedCuts+"_JESMin.tex").c_str();
  if(JESPlusResults == true) PresentationTexTitle = ("BTagPerformanceStudy/"+UsedTrigger+"_"+dataSet+decayChannel+AppliedCuts+"_JESPlus.tex").c_str();
  if(WSystResults == true && WSystPositive == false) PresentationTexTitle = ("BTagPerformanceStudy/"+UsedTrigger+"_"+dataSet+decayChannel+AppliedCuts+"_WMin.tex").c_str();
  if(WSystResults == true && WSystPositive == true) PresentationTexTitle = ("BTagPerformanceStudy/"+UsedTrigger+"_"+dataSet+decayChannel+AppliedCuts+"_WPlus.tex").c_str();
 
  ofstream PresentationTexHadr(("HadrKinFit/"+PresentationTexTitle).c_str());
  ofstream PresentationTexHadrAndLeptWOnly(("HadrAndLeptWOnlyKinFit/"+PresentationTexTitle).c_str());
  ofstream PresentationTexHadrAndLept(("HadrAndLeptKinFit/"+PresentationTexTitle).c_str());

  //ofstream PresentationTex(("BTagPerformanceStudy/WorkshopPresentation"+UsedTrigger+"NoBTagEvtSelTopFittedProbCut.tex").c_str());
  //ofstream PresentationTex(("BTagPerformanceStudy/WorkshopPresentation"+UsedTrigger+"NoBTagEvtSelTopFittedMonteCarlo.tex").c_str());
  //ofstream PresentationTex(("BTagPerformanceStudy/WorkshopPresentation"+UsedTrigger+"NoBTagEvtSelTopFittedJESMinProbCut.tex").c_str());
  //ofstream PresentationTex(("BTagPerformanceStudy/WorkshopPresentation"+UsedTrigger+"NoBTagEvtSelTopFittedJESPlusProbCut.tex").c_str());
  //ofstream PresentationTex(("BTagPerformanceStudy/WorkshopPresentation"+UsedTrigger+"NoBTagEvtSelTopFittedWMin100.tex").c_str());
  //ofstream PresentationTex(("BTagPerformanceStudy/WorkshopPresentation"+UsedTrigger+"NoBTagEvtSelTopFittedWMinProbCut.tex").c_str());
  //ofstream PresentationTex(("BTagPerformanceStudy/WorkshopPresentation"+UsedTrigger+"NoBTagEvtSelTopFittedWPlus100.tex").c_str());
  //ofstream PresentationTex(("BTagPerformanceStudy/WorkshopPresentation"+UsedTrigger+"NoBTagEvtSelTopFittedWPlusProbCut.tex").c_str());

  //oooooOOOOOOOOOOOoooooooooooOOOOOOOOOOOOOOOOoooooooooOOOOOOOOOOOoooooooooooooooooOOOOOOOOOOOOOOOoooooooooooooOOOOOOOOOOOoooooooo
  //2) Ciemat results: SSVHEM btag at event selection, leptonic top mass left free in KinFit and offline muon pt cut of 25
  //oooooOOOOOOOOOOOoooooooooooOOOOOOOOOOOOOOOOoooooooooOOOOOOOOOOOoooooooooooooooooOOOOOOOOOOOOOOOoooooooooooooOOOOOOOOOOOoooooooo

  //ofstream PresentationTex(("BTagPerformanceStudy/CiematTest"+UsedTrigger+".tex").c_str());
  //ofstream PresentationTex(("BTagPerformanceStudy/CiematMonteCarlo"+UsedTrigger+".tex").c_str());
  //ofstream PresentationTex(("BTagPerformanceStudy/Ciemat"+UsedTrigger+"JESMin.tex").c_str());
  //ofstream PresentationTex(("BTagPerformanceStudy/Ciemat"+UsedTrigger+"JESPlus.tex").c_str());
  //ofstream PresentationTex(("BTagPerformanceStudy/Ciemat"+UsedTrigger+"WMin100.tex").c_str());
  //ofstream PresentationTex(("BTagPerformanceStudy/Ciemat"+UsedTrigger+"WPlus100.tex").c_str());

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
  histo1D["LeptWMassHadrAndLeptWOnly"] = new TH1F("LeptWMassHadrAndLeptWOnly","LeptWMassHadrAndLeptWOnly",40,0,150);
  histo1D["LeptTopMassHadrAndLeptWOnly"] = new TH1F("LeptTopMassHadrAndLeptWOnly","LeptTopMassHadrAndLeptWOnly",60,50,350);
  histo1D["LeptWMassHadrAndLept"] = new TH1F("LeptWMassHadrAndLept","LeptWMassHadrAndLept",40,0,150);
  histo1D["LeptTopMassHadrAndLept"] = new TH1F("LeptTopMassHadrAndLept","LeptTopMassHadrAndLept",60,50,350);

  //Histograms to obtain ratio of cos theta* for alternative and SM helicities
  float HelicityWeight[3];
  int SizeArray = 3;
  float HelicityFraction[SizeArray][3];  //0:Longitudinal; 1:Righthanded; 2:Lefthanded
  float UsedDistributionValue[SizeArray];
  std::ostringstream HelicityNumbers[SizeArray][3];
  for(int helicityNumbers=0;helicityNumbers<SizeArray;helicityNumbers++){
    HelicityWeight[helicityNumbers]=1;
    UsedDistributionValue[helicityNumbers]=0;
    if(helicityNumbers == 0){
      HelicityFraction[helicityNumbers][0]=0.5;
      HelicityFraction[helicityNumbers][1]=0.;
      HelicityFraction[helicityNumbers][2]=0.5;
    }
    else if(helicityNumbers==1){
      HelicityFraction[helicityNumbers][0]=0.6;
      HelicityFraction[helicityNumbers][1]=0.2;
      HelicityFraction[helicityNumbers][2]=0.2;
    }
    else if(helicityNumbers==2){
      HelicityFraction[helicityNumbers][0]=0.65;
      HelicityFraction[helicityNumbers][1]=0.1;
      HelicityFraction[helicityNumbers][2]=0.25;
    }
    HelicityNumbers[helicityNumbers][0] << HelicityFraction[helicityNumbers][0];
    HelicityNumbers[helicityNumbers][1] << HelicityFraction[helicityNumbers][1];
    HelicityNumbers[helicityNumbers][2] << HelicityFraction[helicityNumbers][2];
    
    std::string HistoName = "CosThetaRight"+HelicityNumbers[helicityNumbers][1].str()+"Long"+HelicityNumbers[helicityNumbers][0].str()+"Left"+HelicityNumbers[helicityNumbers][2].str();
    TString THistoName = "CosThetaRight"+HelicityNumbers[helicityNumbers][1].str()+"Long"+HelicityNumbers[helicityNumbers][0].str()+"Left"+HelicityNumbers[helicityNumbers][2].str();
    std::string MSSVHEbTagHistoName = "MSSVHEbTagCosThetaRight"+HelicityNumbers[helicityNumbers][1].str()+"Long"+HelicityNumbers[helicityNumbers][0].str()+"Left"+HelicityNumbers[helicityNumbers][2].str();
    TString MSSVHEbTagTHistoName = "MSSVHEbTagCosThetaRight"+HelicityNumbers[helicityNumbers][1].str()+"Long"+HelicityNumbers[helicityNumbers][0].str()+"Left"+HelicityNumbers[helicityNumbers][2].str();
    std::string TCSVbTagHistoName = "TCSVbTagCosThetaRight"+HelicityNumbers[helicityNumbers][1].str()+"Long"+HelicityNumbers[helicityNumbers][0].str()+"Left"+HelicityNumbers[helicityNumbers][2].str();
    TString TCSVbTagTHistoName = "TCSVbTagCosThetaRight"+HelicityNumbers[helicityNumbers][1].str()+"Long"+HelicityNumbers[helicityNumbers][0].str()+"Left"+HelicityNumbers[helicityNumbers][2].str();

    histo1D[HistoName] = new TH1F(THistoName,THistoName,200,-1,1);               
    histo1D[MSSVHEbTagHistoName] = new TH1F(MSSVHEbTagTHistoName,MSSVHEbTagTHistoName,200,-1,1);               
    histo1D[TCSVbTagHistoName] = new TH1F(TCSVbTagTHistoName,TCSVbTagTHistoName,200,-1,1);               
  }

  float XSection;
  float EqLumi;
  vector<Dataset*> datasets; // needed for MSPlots
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
    
    //Nominal samples:
    if( dataSet->Name().find("TTbarJets_SemiMuon") == 0 && dataSet->Name().find("JES") != 0) color = kRed+1;
    if( dataSet->Name().find("TTbarJets_SemiEl") == 0 && dataSet->Name().find("JES") != 0) color = kRed-4;
    if( dataSet->Name().find("TTbarJets_Other") == 0 && dataSet->Name().find("JES") != 0 ) color = kRed-7;
    if( dataSet->Name().find("WJets") == 0 )
      {
	dataSet->SetTitle("W#rightarrowl#nu");
	color = kGreen-3;
      }
    if( dataSet->Name().find("ZJets") == 0 && dataSet->Name().find("JES") != 0 )
      {
	dataSet->SetTitle("Z/#gamma*#rightarrowl^{+}l^{-}");
	color = kAzure-2;
      }
    if( dataSet->Name().find("ST") == 0 && dataSet->Name().find("JES") != 0 ) color = kMagenta;

    if( dataSet->Name().find("QCD") == 0 && dataSet->Name().find("JES") !=0 ) color = kBlue;

    //Systematics samples:
    //JES:
    if( dataSet->Name().find("JES_TTbarJets_SemiMuon") == 0){ 
      color = kRed+1;
      dataSet->SetTitle("t#bar{t}+jets semi-#mu (JES)");
    }
    if( dataSet->Name().find("JES_TTbarJets_SemiEl") == 0){ 
      color = kRed-4;
      dataSet->SetTitle("t#bar{t}+jets semi-el (JES)");
    }
    if( dataSet->Name().find("JES_TTbarJets_Other") == 0){ 
      color = kRed-7;
      dataSet->SetTitle("t#bar{t}+jets other (JES)");
    }
    if( dataSet->Name().find("JES_WJets") == 0 )
      {
	dataSet->SetTitle("W#rightarrowl#nu (JES)");
	color = kGreen-3;
      }
    if( dataSet->Name().find("JES_ZJets") == 0)
      {
	dataSet->SetTitle("Z/#gamma*#rightarrowl^{+}l^{-} (JES)");
	color = kAzure-2;
      }
    if( dataSet->Name().find("JES_ST") == 0){ 
      color = kMagenta;
      dataSet->SetTitle("Single-Top (JES)");
    }
    
    //WJets scale up/down
    if( dataSet->Name().find("Syst_WJets") == 0 )
      {
	dataSet->SetTitle("W#rightarrowl#nu (WJets up/down)");
	color = kGreen-3;
      }    

    Dataset* tmpDS = new Dataset(dataSet->Name(), dataSet->Title(), dataSet->DoIt(), color, dataSet->LineStyle(), dataSet->LineWidth(), dataSet->NormFactor(), XSection);
    tmpDS->SetEquivalentLuminosity( EqLumi );
    datasets.push_back( tmpDS );
  }

  cout << " colors defined " << endl;  

  MSPlot["ChiSqHadr"] = new MultiSamplePlot(datasets,"ChiSqHadr",50,0,50,"ChiSqHadr");
  MSPlot["ChiSqHadrAndLeptWOnly"] = new MultiSamplePlot(datasets,"ChiSqHadrAndLeptWOnly",50,0,50,"ChiSqHadrAndLeptWOnly");
  MSPlot["ChiSqHadrAndLept"] = new MultiSamplePlot(datasets,"ChiSqHadrAndLept",50,0,50,"ChiSqHadrAndLept");

  MSPlot["LeptWMassHadr"] = new MultiSamplePlot("datasets,LeptWMassHadr",40,0,150,"LeptWMassHadr");
  MSPlot["LeptTopMassHadr"] = new MultiSamplePlot(datasets,"LeptTopMassHadr",60,50,350,"LeptTopMassHadrAndLept");
  MSPlot["LeptWMassHadrAndLeptWOnly"] = new MultiSamplePlot(datasets,"LeptWMassHadrAndLeptWOnly",40,0,150,"LeptWMassHadr");
  MSPlot["LeptTopMassHadrAndLeptWOnly"] = new MultiSamplePlot(datasets,"LeptTopMassHadrAndLeptWOnly",60,50,350,"LeptTopMassHadrAndLept");
  MSPlot["LeptWMassHadrAndLept"] = new MultiSamplePlot(datasets,"LeptWMassHadrAndLept",40,0,150,"LeptWMassHadr");
  MSPlot["LeptTopMassHadrAndLept"] = new MultiSamplePlot(datasets,"LeptTopMassHadrAndLept",60,50,350,"LeptTopMassHadrAndLept");

  // //MSPlots for mass distributions after Kinematic Fit:
  // MSPlot["HadronicWMass"] = new MultiSamplePlot(datasets,"HadronicWMass",50,50,100,"HadronicWMass");
  // MSPlot["HadronicTopMass"] = new MultiSamplePlot(datasets,"HadronicTopMass",50,150,200,"HadronicTopMass");
  // MSPlot["LeptonicWMass"] = new MultiSamplePlot(datasets,"LeptonicWMass",50,50,100,"LeptonicWMass");
  // MSPlot["LeptonicTopMass"] = new MultiSamplePlot(datasets,"LeptonicTopMass",50,150,200,"LeptonicTopMass");
  
  // //MSPlots for Cos theta distribution and KinFit probability distribution
  // MSPlot["CosThetaSSVHEMLept"]= new MultiSamplePlot(datasets, "CosThetaSSVHEMLept", CosThetaBinNumber,-1,1,"CosThetaSSVHEMLept");  
  // MSPlot["KinFitProbSSVHEMLept"]= new MultiSamplePlot(datasets, "KinFitProbSSVHEMLept", CosThetaBinNumber,0,1,"KinFitProbSSVHEMLept");  
  // MSPlot["CosThetaTCHEMLept"]= new MultiSamplePlot(datasets, "CosThetaTCHEMLept", CosThetaBinNumber,-1,1,"CosThetaTCHEMLept");  
  // MSPlot["KinFitProbTCHEMLept"]= new MultiSamplePlot(datasets, "KinFitProbTCHEMLept", CosThetaBinNumber,0,1,"KinFitProbTCHEMLept");  
  // MSPlot["CosThetaTCHPMLept"]= new MultiSamplePlot(datasets, "CosThetaTCHPMLept", CosThetaBinNumber,-1,1,"CosThetaTCHPMLept");  
  // MSPlot["CosThetaSSVHPMLept"]= new MultiSamplePlot(datasets, "CosThetaSSVHPMLept", CosThetaBinNumber,-1,1,"CosThetaSSVHPMLept");  
  // MSPlot["KinFitProbSSVHPMLept"]= new MultiSamplePlot(datasets, "KinFitProbSSVHPMLept", CosThetaBinNumber,0,1,"KinFitProbSSVHPMLept");  
  // MSPlot["CosThetaCSVMLept"]= new MultiSamplePlot(datasets, "CosThetaCSVMLept", CosThetaBinNumber,-1,1,"CosThetaCSVMLept");
 
  // //Check nPrimary vertices for different executed cuts !!
  MSPlot["nPrimaryVert"] = new MultiSamplePlot(datasets,"nPrimaryVert" , 20, 0, 20, "nPrimaryVert");
  MSPlot["nPVBeforeCuts"] = new MultiSamplePlot(datasets,"nPVBeforeCuts" , 20, 0, 20, "nPVBeforeCuts");
  MSPlot["nPVAfterSSVHEMbTag"] = new MultiSamplePlot(datasets,"nPVAfterSSVHEMbTag" , 20, 0, 20, "nPVAfterSSVHEMbTag");
  MSPlot["nPVAfterMuon27Cut"] = new MultiSamplePlot(datasets,"nPVAfterMuon27Cut" , 20, 0, 20, "nPVAfterMuon27Cut");
  MSPlot["nPVAfterTransverseMassCut"] = new MultiSamplePlot(datasets,"nPVAfterTransverseMassCut" , 20, 0, 20, "nPVAfterTransverseMassCut");
  MSPlot["nPVBeforeFoundJetComb"] = new MultiSamplePlot(datasets, "nPVBeforeFoundJetComb", 20, 0, 20,"nPVBeforeFoundJetComb");
  // MSPlot["nPVAfterFoundJetComb"] = new MultiSamplePlot(datasets, "nPVAfterFoundJetComb", 20, 0, 20,"nPVAfterFoundJetComb");
  // MSPlot["nPVAfterFoundJetCombbTag"] = new MultiSamplePlot(datasets, "nPVAfterFoundJetCombbTab", 20, 0, 20,"nPVAfterFoundJetCombbTag");
  // MSPlot["nPVAfterFoundCosTheta"] = new MultiSamplePlot(datasets, "nPVAfterFoundCosTheta", 20, 0, 20,"nPVAfterFoundCosTheta");
  // MSPlot["nPVAfterFoundCosThetabTag"] = new MultiSamplePlot(datasets, "nPVAfterFoundCosThetabTag", 20, 0, 20,"nPVAfterFoundCosThetabTag");

  MSPlot["TransverseMassBeforeCut"]=new MultiSamplePlot(datasets,"TransverseMassBeforeCut",50,0,200,"TransverseMassBeforeCut");
  MSPlot["TransverseMassAfterCut"]=new MultiSamplePlot(datasets,"TransverseMassAfterCut",50,0,200,"TransverseMassAfterCut"); 
  
  MSPlot["RelIsoLepton"]=new MultiSamplePlot(datasets,"RelIsoLepton",50,0,200,"RelIsoLepton");
  
  MSPlot["KinFitProbabilityHadr"]=new MultiSamplePlot(datasets,"KinFitProbabilityHadr",50,-1,1,"KinFitProbabilityHadr");  
  MSPlot["KinFitProbabilityHadrAndLeptWOnly"]=new MultiSamplePlot(datasets,"KinFitProbabilityHadrAndLeptWOnly",50,-1,1,"KinFitProbabilityHadrAndLeptWOnly");
  MSPlot["KinFitProbabilityHadrAndLept"]=new MultiSamplePlot(datasets,"KinFitProbabilityHadrAndLept",50,-1,1,"KinFitProbabilityHadrAndLept");
  MSPlot["NeutrinoMassHadr"]=new MultiSamplePlot(datasets,"NeutrinoMassHadr",50,-1,1,"NeutrinoMassHadr");  
  MSPlot["NeutrinoMassHadrAndLeptWOnly"]=new MultiSamplePlot(datasets,"NeutrinoMassHadrAndLeptWOnly",50,-1,1,"NeutrinoMassHadrAndLeptWOnly");
  MSPlot["NeutrinoMassHadrAndLept"]=new MultiSamplePlot(datasets,"NeutrinoMassHadrAndLept",50,-1,1,"NeutrinoMassHadrAndLept");
  
  // //Histograms to check differences between events with Negative and positive DSquared (And check influence of probability cut) 
  // MSPlot["CosThetaNobTag"]=new MultiSamplePlot(datasets,"CosThetaNobTag" ,CosThetaBinNumber,-1,1, "CosThetaNobTag");
  // MSPlot["CosThetaNobTagProbCut"]=new MultiSamplePlot(datasets,"CosThetaNobTagProbCut" ,CosThetaBinNumber*2,-2,2,"CosThetaNobTagProbCut" );
  // MSPlot["KinFitProbNobTag"]=new MultiSamplePlot(datasets, "KinFitProbNobTag",25,0,1,"KinFitProbNobTag");

  MSPlot["Jet1Pt"] = new MultiSamplePlot(datasets, "Jet1Pt", 100,0,500,"Jet1Pt");
  MSPlot["Jet2Pt"] = new MultiSamplePlot(datasets, "Jet2Pt", 100,0,500,"Jet2Pt");
  MSPlot["Jet3Pt"] = new MultiSamplePlot(datasets, "Jet3Pt", 100,0,500,"Jet3Pt");
  MSPlot["Jet4Pt"] = new MultiSamplePlot(datasets, "Jet4Pt", 100,0,500,"Jet4Pt");
  
  // if(UseChangedKinematics == true){
  //   MSPlot["BLeptPtAfterKinFit"]= new MultiSamplePlot(datasets,"BLeptPtAfterKinFit",50,0,250,"BLeptPtAfterKinFit");
  //   MSPlot["MetPtAfterKinFit"] = new MultiSamplePlot(datasets,"MetPtAfterKinFit",50,0,200,"MetPtAfterKinFit");
  //   MSPlot["MuonPtAfterKinFit"] = new MultiSamplePlot(datasets,"MuonPtAfterKinFit",40,0,150,"MuonPtAfterKinFit");
  //   MSPlot["WLeptPtAfterKinFit"] = new MultiSamplePlot(datasets,"WLeptPtAfterKinFit",50,0,250,"WLeptPtAfterKinFit");
  //   MSPlot["TopLeptPtAfterKinFit"] = new MultiSamplePlot(datasets,"TopLeptPtAfterKinFit",70,0,400,"TopLeptPtAfterKinFit");
  // }
  // else if(UseChangedKinematics == false){
  //   MSPlot["BLeptPtBeforeKinFit"]= new MultiSamplePlot(datasets,"BLeptPtBeforeKinFit",70,0,400,"BLeptPtBeforeKinFit");
  //   MSPlot["MetPtBeforeKinFit"] = new MultiSamplePlot(datasets,"MetPtBeforeKinFit",50,0,200,"MetPtBeforeKinFit");
  //   MSPlot["MuonPtBeforeKinFit"] = new MultiSamplePlot(datasets,"MuonPtBeforeKinFit",40,0,150,"MuonPtBeforeKinFit");
  //   MSPlot["WLeptPtBeforeKinFit"] = new MultiSamplePlot(datasets,"WLeptPtBeforeKinFit",50,0,250,"WLeptPtBeforeKinFit");
  //   MSPlot["TopLeptPtBeforeKinFit"] = new MultiSamplePlot(datasets,"TopLeptPtBeforeKinFit",70,0,400,"TopLeptPtBeforeKinFit");
  // }
  
  MSPlot["LeptonMass"]=new MultiSamplePlot(datasets,"LeptonMass",20,0,1,"LeptonMass");
  MSPlot["LeptonPx"]=new MultiSamplePlot(datasets,"LeptonPx",60,-20,20,"LeptonPx");
  MSPlot["LeptonPy"]=new MultiSamplePlot(datasets,"LeptonPy",60,-20,20,"LeptonPy");
  MSPlot["LeptonPz"]=new MultiSamplePlot(datasets,"LeptonPz",60,-20,20,"LeptonPz");
    
  // Zie Code Stijn voor alle gebruikte controle plots

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
  // Summer11 PU_S4, distribution obtained by averaging the number of interactions, taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupMCReweightingUtilities
  // in each beam crossing to estimate the true mean.  THIS IS THE RECOMMENDED ONE for reweighting.
  Double_t probdistFlat10[25] = {
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0630151648,
    0.0526654164,
    0.0402754482,
    0.0292988928,
    0.0194384503,
    0.0122016783,
    0.007207042,
    0.004003637,
    0.0020278322,
    0.0010739954,
    0.0004595759,
    0.0002229748,
    0.0001028162,
    4.58337152809607E-05
  };
  
  cout << " Defining MCLumi_f[25] " << endl;
   
  Double_t MCLumi_f[25] = {
    0.104109,
    0.0703573,
    0.0698445,
    0.0698254,
    0.0697054,
    0.0697907,
    0.0696751,
    0.0694486,
    0.0680332,
    0.0651044,
    0.0598036,
    0.0527395,
    0.0439513,
    0.0352202,
    0.0266714,
    0.019411,
    0.0133974,
    0.00898536,
    0.0057516,
    0.00351493,
    0.00212087,
    0.00122891,
    0.00070592,
    0.000384744,
    0.000219377
  };
  
  cout << " Defining TopDBDist2011Data_f[25] " << endl;
  
  Double_t TopDBDist2011Data_f[25] = {
    0.0127118660008111155,
    0.0273174253882752516,
    0.0647422373974094190,
    0.108494213975257103,
    0.140081296984992526,
    0.150411260268535935,
    0.142773479388604602,
    0.118012735306947752,
    0.0881395784021791473,
    0.0603740700218931975,
    0.0382939204454870868,
    0.0227366747939989136,
    0.0127228459417252551,
    0.00674674468025676568,
    0.00340977235841692389,
    0.00165292589154045016,
    0.000771798466244840342,
    0.000347480158040664431,
    0.000151563397272207710,
    0.0000642172483977206039,
    0.0000264962736283059724,
    0.0000106455374332742453,
    0.00000418355451211455042,
    0.00000161033109693768961,
    0.000000606815958689117662
  };

  cout << " Defining TrueDist2011_f[25] " << endl;
  
  Double_t TrueDist2011_f[25] = {
    0.019091,
    0.0293974,
    0.0667931,
    0.108859,
    0.139533,
    0.149342,
    0.138629,
    0.114582,
    0.0859364,
    0.059324,
    0.0381123,
    0.0229881,
    0.0131129,
    0.00711764,
    0.00369635,
    0.00184543,
    0.000889604,
    0.000415683,
    0.000188921,
    0.000146288,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0
  };
  
  std::cout << " Values defined " << endl;
  
  vector<float> TrueDist2011, MClumi, Spring11MClumi, TopDBDist2011Data;
  for( int i=0; i<25; ++i){
    TopDBDist2011Data.push_back(TopDBDist2011Data_f[i]);
    TrueDist2011.push_back(TrueDist2011_f[i]);
    MClumi.push_back(MCLumi_f[i]);
    Spring11MClumi.push_back(probdistFlat10[i]);
  }

  cout << " Starting LumiReWeighting stuff " << endl;

  //  LumiReWeighting LumiWeightsSpring11 = LumiReWeighting(Spring11MClumi, TrueDist2011);
  //  LumiReWeighting LumiWeights = LumiReWeighting(MClumi, TopDBDist2011Data);
  //LumiReWeighting LumiWeights = LumiReWeighting("PileUpReweighting/pileup_WJets_36bins.root", "PileUpReweighting/pileup_2011Data_UpToRun173692.root", "pileup2", "pileup");
  cout << " Lumi3D Reweighting stuff " << endl;
  Lumi3DReWeighting Lumi3DWeights;
  if(fullDataSet == true) 
    Lumi3DReWeighting("PileUpReweighting/pileup_MC_Fall11.root","PileUpReweighting/pileup_FineBin_2011Data_UpToRun180252.root", "pileup", "pileup");
  else if(fullDataSet == false) 
    Lumi3DWeights = Lumi3DReWeighting("PileUpReweighting/pileup_MC_Flat10PlusTail.root", "PileUpReweighting/pileup_FineBin_2011Data_UpToRun173692.root", "pileup", "pileup");

  Lumi3DWeights.weight3D_init(1.0);
  if(systematic == "PUMinus")
    Lumi3DWeights.weight3D_init(0.92);
  else if(systematic == "PUPlus")
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
  int TotalNumberbTags = NumberTCHEbTags + 1 + NumberTCHPbTags+1 + NumberSSVHEbTags + 1 + NumberSSVHPbTags + 1 + NumberCSVbTags + 1;
  
  std::string bTagFileOutput[TotalNumberbTags];
  int bTagBoolean[TotalNumberbTags];
  map<int, int[5]> ProcessedbTags;
  std::string PresentationOutput[TotalNumberbTags];
  int NumberRemainingEventsKinFitHadr[TotalNumberbTags][inputWTree.size()];
  int NumberRemainingEventsKinFitHadrAndLeptWOnly[TotalNumberbTags][inputWTree.size()];
  int NumberRemainingEventsKinFitHadrAndLept[TotalNumberbTags][inputWTree.size()];
  int NumberBLeptCorrectEventsKinFitHadr[TotalNumberbTags][inputWTree.size()];
  int NumberBLeptCorrectEventsKinFitHadrAndLeptWOnly[TotalNumberbTags][inputWTree.size()];
  int NumberBLeptCorrectEventsKinFitHadrAndLept[TotalNumberbTags][inputWTree.size()];
  int LumiWeightVectorKinFitHadr[TotalNumberbTags][inputWTree.size()];
  int LumiWeightVectorKinFitHadrAndLeptWOnly[TotalNumberbTags][inputWTree.size()];
  int LumiWeightVectorKinFitHadrAndLept[TotalNumberbTags][inputWTree.size()];

  //Initialize:
  for(int ii=0;ii<14;ii++){
    for(int jj=0;jj<14;jj++){
      for(int kk=0;kk<14;kk++){
	for(int ll=0;ll<14;ll++){
	  for(int nn=0;nn<14;nn++){
	    bTagFileOutput[ii+jj+kk+ll+nn]=" Wrong entry chosen";
	    bTagBoolean[ii+jj+kk+ll+nn]=0;
	    PresentationOutput[ii+jj+kk+ll+nn]=" Wrong entry chosen";
	    for(int mm=0;mm<inputWTree.size();mm++){
	      NumberRemainingEventsKinFitHadr[ii+jj+kk+ll+nn][mm]=0;
	      NumberRemainingEventsKinFitHadrAndLeptWOnly[ii+jj+kk+ll+nn][mm]=0;
	      NumberRemainingEventsKinFitHadrAndLept[ii+jj+kk+ll+nn][mm]=0;
	      NumberBLeptCorrectEventsKinFitHadr[ii+jj+kk+ll+nn][mm]=0;
	      NumberBLeptCorrectEventsKinFitHadrAndLeptWOnly[ii+jj+kk+ll+nn][mm]=0;
	      NumberBLeptCorrectEventsKinFitHadrAndLept[ii+jj+kk+ll+nn][mm]=0;
	    }
	  }
	}
      }
    }
  }

  //Call classes made :
  BTagName bTagName = BTagName();  //for bTagFileOutput name giving
  BTagJetSelection bTagJetSelection = BTagJetSelection();
  BTagCosThetaCalculation bTagCosThetaCalculation = BTagCosThetaCalculation();

  int LengthOfPresentationArray=0;

  int NumberSelectedEvents =0;
  int NumberEventsBeforeCuts = 0;
  int NumberEventsAfterbTag = 0;
  int NumberEventsAfterTransverseMass = 0;
  int NumberUsedEvents = 0;
  int NumberUsedCorrectEvents = 0;
  int NumberUsedWrongEvents = 0;
  int NumberUsedDataEvents = 0;
  int NumberSelectedDataEvents = 0;
  int NumberDataEventsBeforeCuts = 0;

  /////////////////////////////////////////
  // Loop on datasets
  /////////////////////////////////////////
  
  for(unsigned int iDataSet=0; iDataSet<inputWTree.size(); iDataSet++){
    cout << " " << endl;

    // string dataSetName = datasets[iDataSet]->Name();
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
    cout << "Before changing --> Cross section = " << dataSet->Xsection() << "  intLumi = " << dataSet->EquivalentLumi() << " Normfactor = " << dataSet->NormFactor() << endl;
    float NominalNormFactor = dataSet->NormFactor();
    if( dataSetName.find("Syst_WJets") == 0 && WSystResults == true){
      if(WSystPositive == false) dataSet->SetEquivalentLuminosity( dataSet->EquivalentLumi() / (0.7) );  //WJets Minus 30%
      //dataSet->SetEquivalentLuminosity( dataSet->EquivalentLumi() / (0.0000001) );  //WJets Minus 100%
      if(WSystPositive == true) dataSet->SetEquivalentLuminosity( dataSet->EquivalentLumi() / (1.3) ); //WJets Plus 30 %
      //dataSet->SetEquivalentLuminosity( dataSet->EquivalentLumi() / (2.) ); //WJets Plus 100 %
      //Normfactor value changes without having to change XSection value !!
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

    //Order of indices used in Kinematic Fit:
    int BLeptIndex[12]={3,2,3,1,2,1,3,0,2,0,1,0};
    int BHadrIndex[12]={2,3,1,3,1,2,0,3,0,2,0,1};
    int Quark1Index[12]={0,0,0,0,0,0,1,1,1,1,2,2};
    int Quark2Index[12]={1,1,2,2,3,3,2,2,3,3,3,3};

    outFile << " Kinematic Fit ordering arrays initialized " << endl;

    //Integer needed to represent the first event since iEvt = 0 does not pass the DCoefficient requirement for ttbar Other
    int FirstProcessedEvent=0;

    outFile << " FirstProcessedEvent initialized " << endl;

    vector<float> CosThetaValuesKinFitHadr[TotalNumberbTags];
    vector<float> CosThetaValuesKinFitHadrAndLeptWOnly[TotalNumberbTags];
    vector<float> CosThetaValuesKinFitHadrAndLept[TotalNumberbTags];
    vector<float> LumiWeightVectorKinFitHadr[TotalNumberbTags];
    vector<float> LumiWeightVectorKinFitHadrAndLeptWOnly[TotalNumberbTags];
    vector<float> LumiWeightVectorKinFitHadrAndLept[TotalNumberbTags];
    int FilledEntries = 0;
    vector<double> CosThGenKinFitHadr[TotalNumberbTags];
    vector<double> CosThGenKinFitHadrAndLeptWOnly[TotalNumberbTags];
    vector<double> CosThGenKinFitHadrAndLept[TotalNumberbTags];
    vector<double> EventCorrectionWeightKinFitHadr[TotalNumberbTags];
    vector<double> EventCorrectionWeightKinFitHadrAndLeptWOnly[TotalNumberbTags];
    vector<double> EventCorrectionWeightKinFitHadrAndLept[TotalNumberbTags];
    float binEdge[CosThetaBinNumber+1];
    float binSize = (1.-(-1.))/15.;
    for(int ii=0; ii<=CosThetaBinNumber;ii++){
      binEdge[ii] = -1 + binSize*ii;
    }

    //Initialize naming of different bTag options:
    int TCHELoop=1;
    int TCHPLoop=1;
    int SSVHELoop=1;
    int SSVHPLoop=1;
    int CSVLoop=1;
       
    int UsedTCHE[TotalNumberbTags];
    int UsedTCHP[TotalNumberbTags];
    int UsedSSVHE[TotalNumberbTags];
    int UsedSSVHP[TotalNumberbTags];
    int UsedCSV[TotalNumberbTags];
    if(iDataSet==0){
      while(CSVLoop<=NumberCSVbTags){
	while(SSVHPLoop<=NumberSSVHPbTags){  
	  while(SSVHELoop<=NumberSSVHEbTags){
	    while(TCHPLoop<=NumberTCHPbTags){
	      while(TCHELoop<=NumberTCHEbTags){
		bTagFileOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = bTagName.NameGiving(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags);
		bTagBoolean[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = 1;
		UsedTCHE[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = TCHELoop;
		UsedTCHP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = TCHPLoop;
		UsedSSVHE[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = SSVHELoop;
		UsedSSVHP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = SSVHPLoop;
		UsedCSV[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = CSVLoop;
		
		PresentationOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = bTagName.NameGivingPres(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags);
		TCHELoop++;
		if(TCHELoop == 14){TCHPLoop=2;SSVHELoop=2;SSVHPLoop=2;CSVLoop=2;}
	      }
	      bTagFileOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = bTagName.NameGiving(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags);
	      bTagBoolean[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = 1;
	      UsedTCHE[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = TCHELoop;
	      UsedTCHP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = TCHPLoop;
	      UsedSSVHE[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = SSVHELoop;
	      UsedSSVHP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = SSVHPLoop;
	      UsedCSV[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = CSVLoop;
	      
	      PresentationOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = bTagName.NameGivingPres(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags);
	      TCHPLoop++;
	    }
	    bTagFileOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = bTagName.NameGiving(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags);
	    bTagBoolean[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = 1;
	    UsedTCHE[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = TCHELoop;
	    UsedTCHP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = TCHPLoop;
	    UsedSSVHE[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = SSVHELoop;
	    UsedSSVHP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = SSVHPLoop;
	    UsedCSV[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = CSVLoop;

	    PresentationOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = bTagName.NameGivingPres(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags);
	    SSVHELoop++;
	  }
	  bTagFileOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = bTagName.NameGiving(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags);
	  bTagBoolean[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = 1;
	  UsedTCHE[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = TCHELoop;
	  UsedTCHP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = TCHPLoop;
	  UsedSSVHE[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = SSVHELoop;
	  UsedSSVHP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = SSVHPLoop;
	  UsedCSV[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = CSVLoop;
	  
	  PresentationOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = bTagName.NameGivingPres(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags);
	  SSVHPLoop++;		
	}	
	bTagFileOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = bTagName.NameGiving(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags);
	bTagBoolean[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = 1;
	UsedTCHE[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = TCHELoop;
	UsedTCHP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = TCHPLoop;
	UsedSSVHE[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = SSVHELoop;
	UsedSSVHP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = SSVHPLoop;
	UsedCSV[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = CSVLoop;
	
	PresentationOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = bTagName.NameGivingPres(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags);
	CSVLoop++;
      }
    }
    
    /////////////////////////////////////////
    // Loop on events
    /////////////////////////////////////////
    float ProbabilityPt10ToInf[10][10];
    float ProbabilityPt30ToInf[10][10];
    float TotalPt10ToInf[10];			    
    float TotalPt30ToInf[10];			    
    for(int ii=0;ii<10;ii++){
      TotalPt10ToInf[ii] = 0;
      TotalPt30ToInf[ii] = 0;
      for(int jj=0;jj<10;jj++){
	ProbabilityPt10ToInf[ii][jj]=0;
	ProbabilityPt30ToInf[ii][jj]=0;
      }
    }
    
    for(unsigned int iEvt=0; iEvt<nEvent; iEvt++){
      //for(unsigned int iEvt=0; iEvt<3000; iEvt++){

      inWTreeTree->GetEvent(iEvt);
      if(iEvt%10000 == 0)
        std::cout<<"Processing the "<<iEvt<<"th event, time = "<< ((double)clock() - start) / CLOCKS_PER_SEC <<endl;
      
      // PU reweighting???
      float avPU = ( (float)wTree->nPUBXm1() + (float)wTree->nPU() + (float)wTree->nPUBXp1() ) / 3.; // average in 3 BX!!!, as recommended
      //float lumiWeight = LumiWeights.ITweight( wTree->nPU() );
      
      double lumiWeight3D = 1.0;
      if(!(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")){
	lumiWeight3D = Lumi3DWeights.weight3D(wTree->nPUBXm1(),wTree->nPU(),wTree->nPUBXp1());
      }
      histo1D["lumiWeights3D"]->Fill(lumiWeight3D);
      
      // scale factor for the event
      float scaleFactor = 1.; 

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
      if((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ){
	histo1D["StandardCosThetaFit"]->Fill(CosThetaGenerator);  // Histogram with fit   	  
	histo1D["StandardCosTheta"]->Fill(CosThetaGenerator);  // Histogram without fit   	  
      }
                 
      vector<TLorentzVector> selectedJets = wTree->selectedJets();
      TLorentzVector lepton = wTree->lepton(); 
      TLorentzVector MET = wTree->met();

      //if(dataSetName.find("QCD") !=0){
      MSPlot["LeptonMass"]->Fill(lepton.M(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
      MSPlot["LeptonPx"]->Fill(lepton.Px(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
      MSPlot["LeptonPy"]->Fill(lepton.Py(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
      MSPlot["LeptonPz"]->Fill(lepton.Pz(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
      //}

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
      //vector<TLorentzVector> hadrBChangedHadr, light1ChangedHadr, light2ChangedHadr; //-->Only need this ones for W and top mass distribution after KinFit !!
      //vector<TLorentzVector> hadrBChangedHadrAndLeptWOnly, light1ChangedHadrAndLeptWOnly, light2ChangedHadrAndLeptWOnly;
      //vector<TLorentzVector> hadrBChangedHadrAndLept, light1ChangedHadrAndLept, light2ChangedHadrAndLept;

      for(unsigned int iCombi = 0; iCombi<12;iCombi++){	
	
	//Hadronic KinFit   !!Has no changed neutrino, lepton or leptonic bJet component
	ChiSqKinFitHadr.push_back(wTree->hadrKinFitChiSq(iCombi)); 
	MSPlot["ChiSqHadr"]->Fill(ChiSqKinFitHadr[iCombi], datasets[iDataSet], true, Luminosity*scaleFactor);
	
if(ChiSqKinFitHadr[iCombi] != 9999) KinFitProbHadr.push_back(TMath::Prob(ChiSqKinFitHadr[iCombi],2));
	else KinFitProbHadr.push_back(-1.);
	MSPlot["KinFitProbabilityHadr"]->Fill(KinFitProbHadr[iCombi], datasets[iDataSet], true, Luminosity);
	// hadrBChangedHadr.push_back(wTree->hadrKinFitHadrB(iCombi));
	// light1ChangedHadr.push_back(wTree->hadrKinFitLight1(iCombi));
	// light2ChangedHadr.push_back(wTree->hadrKinFitLight2(iCombi));
	
	//Hadronic and leptonic W only KinFit     !!Has no changed leptonic bJet component
	ChiSqKinFitHadrAndLeptWOnly.push_back(wTree->hadrAndLeptWOnlyKinFitChiSq(iCombi));
	MSPlot["ChiSqHadrAndLeptWOnly"]->Fill(ChiSqKinFitHadrAndLeptWOnly[iCombi], datasets[iDataSet], true, Luminosity*scaleFactor);
	if(ChiSqKinFitHadrAndLeptWOnly[iCombi] != 9999) KinFitProbHadrAndLeptWOnly.push_back(TMath::Prob(ChiSqKinFitHadrAndLeptWOnly[iCombi],2));
	else KinFitProbHadrAndLeptWOnly.push_back(-1.);
	MSPlot["KinFitProbabilityHadrAndLeptWOnly"]->Fill(KinFitProbHadrAndLeptWOnly[iCombi], datasets[iDataSet], true, Luminosity);
	leptonChangedHadrAndLeptWOnly.push_back(wTree->hadrAndLeptWOnlyKinFitLepton(iCombi));
	neutrinoChangedHadrAndLeptWOnly.push_back(wTree->hadrAndLeptWOnlyKinFitNeutrino(iCombi));
	// hadrBChangedHadrAndLeptWOnly.push_back(wTree->hadrAndLeptWOnlyKinFitHadrB(iCombi));
	// light1ChangedHadrAndLeptWOnly.push_back(wTree->hadrAndLeptWOnlyKinFitLight1(iCombi));
	// light2ChangedHadrAndLeptWOnly.push_back(wTree->hadrAndLeptWOnlyKinFitLight2(iCombi));

	//Hadronic and leptonic KinFit
	ChiSqKinFitHadrAndLept.push_back(wTree->hadrAndLeptKinFitChiSq(iCombi));
	MSPlot["ChiSqHadrAndLept"]->Fill(ChiSqKinFitHadrAndLept[iCombi], datasets[iDataSet], true, Luminosity*scaleFactor);
	if(ChiSqKinFitHadrAndLept[iCombi] != 9999) KinFitProbHadrAndLept.push_back(TMath::Prob(ChiSqKinFitHadr[iCombi],2));
	else KinFitProbHadrAndLept.push_back(-1.);
	MSPlot["KinFitProbabilityHadrAndLept"]->Fill(KinFitProbHadrAndLept[iCombi], datasets[iDataSet], true, Luminosity);
	leptBChangedHadrAndLept.push_back(wTree->hadrAndLeptKinFitLeptB(iCombi));
	leptonChangedHadrAndLept.push_back(wTree->hadrAndLeptKinFitLepton(iCombi));
	neutrinoChangedHadrAndLept.push_back(wTree->hadrAndLeptKinFitNeutrino(iCombi));
	// hadrBChangedHadrAndLept.push_back(wTree->hadrAndLeptKinFitHadrB(iCombi));
	// light1ChangedHadrAndLept.push_back(wTree->hadrAndLeptKinFitLight1(iCombi));
	// light2ChangedHadrAndLept.push_back(wTree->hadrAndLeptKinFitLight2(iCombi));       

      }//End loop for filling KinFit stuff!!
                                   
      //ooooooooOOOOOOOOOOOOooooooooooooOOOOOOOOOOOOooooooooooooOOOOO
      //ooOOooOOoo      Reading out nTuples done           ooOOooOOoo
      //ooOOooOOoo-----------------------------------------ooOOooOOoo
      //ooOOooOOoo      Start of actual analysis           ooOOooOOoo
      //ooooooooOOOOOOOOOOOOooooooooooooOOOOOOOOOOOOooooooooooooOOOOO
      float TransverseMass = sqrt(2*(abs(lepton.Pt()))*abs(MET.Pt())*(1-cos(lepton.DeltaPhi(MET))));	//Should this not be updated to correct lepton and MET??
      MSPlot["TransverseMassBeforeCut"]->Fill(TransverseMass, datasets[iDataSet], true, Luminosity);	

      float relativeLeptonIsolation = wTree->relIso();
      MSPlot["RelIsoLepton"]->Fill(relativeLeptonIsolation, datasets[iDataSet], true, Luminosity);
      
      //----------------------------------
      //     Require some extra cuts:
      //----------------------------------
      bool eventSelected = false;

      if((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ) NumberEventsBeforeCuts++;
      if(dataSetName.find("Data") ==0) NumberDataEventsBeforeCuts++;
      
      MSPlot["nPVBeforeCuts"]->Fill(nPrimaryVertices,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
      selecTableMacro.Fill(iDataSet,0,scaleFactor*lumiWeight3D);

      if(SSVHEbTag == false && MTCut == false && MuonPtCut == false){
	eventSelected = true;   //No Cuts applied
      }
      else{
	if( (SSVHEbTag == true &&  (btagSSVHE[0] > 1.74 || btagSSVHE[1] > 1.74 || btagSSVHE[2] > 1.74 || btagSSVHE[3] > 1.74))){  //Medium SSVHE bTag
	  selecTableMacro.Fill(iDataSet,1,scaleFactor*lumiWeight3D);
	  
	  if((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ) NumberEventsAfterbTag++;
	  
	  MSPlot["nPVAfterSSVHEMbTag"]->Fill(nPrimaryVertices,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	  eventSelected = true;
	}
	//Require offline muon cut of 27 (to avoid turn-on of IsoMu17/20/24 triggers)  -- Use 25 for CIEMAT comparison (value applied in tree)!
	//No extra offline muon cut for IsoMu17_TriCentralJet30 trigger --> Cut of 20 GeV is already applied
	if((MuonPtCut == true && lepton.Pt() >=27)){
	  selecTableMacro.Fill(iDataSet,2,scaleFactor*lumiWeight3D);
	  
	  MSPlot["nPVAfterMuon27Cut"]->Fill(nPrimaryVertices,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	  eventSelected = true;
	}
	if((MTCut == true && TransverseMass > 30)){
	  
	  selecTableMacro.Fill(iDataSet,3,scaleFactor*lumiWeight3D);	  
	  if((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ) NumberEventsAfterTransverseMass++;
	  
	  MSPlot["TransverseMassAfterCut"]->Fill(TransverseMass, datasets[iDataSet], true, Luminosity);	      	  
	  MSPlot["nPVAfterTransverseMassCut"]->Fill(nPrimaryVertices,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);	  
	  eventSelected = true;
	}
      }
      
      if(eventSelected == true){

	if((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ) NumberSelectedEvents++;
	if(dataSetName.find("Data") ==0) NumberSelectedDataEvents++;
	
	MSPlot["nPrimaryVert"]->Fill(nPrimaryVertices,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);    
	
	//if(dataSetName.find("QCD") != 0){
	  MSPlot["Jet1Pt"]->Fill(selectedJets[0].Pt(), datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);    
	  MSPlot["Jet2Pt"]->Fill(selectedJets[1].Pt(), datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 
	  MSPlot["Jet3Pt"]->Fill(selectedJets[2].Pt(), datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 
	  MSPlot["Jet4Pt"]->Fill(selectedJets[3].Pt(), datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 	  
	  //}
	  
	//-----------------------------------------------------
	//    Reconstruction of correct jet distribution
	//-----------------------------------------------------	  	
	//if(dataSetName.find("QCD") != 0){
	MSPlot["MetPhi"]->Fill(MET.Phi(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);  
	//MSPlot["NeutrinoKinFitPhi"]->Fill(neutrinoKinFit[HighestProbComb].Phi(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);  
	// //MSPlot["NeutrinoPhiDiff"]->Fill((genNeutrino.Phi()-neutrinoKinFit[HighestProbComb].Phi()),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);
	
	MSPlot["MetPx"]->Fill(MET.Px(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);         
	//MSPlot["NeutrinoKinFitPx"]->Fill(neutrinoKinFit[HighestProbComb].Px(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);  
	// //MSPlot["NeutrinoPxDiff"]->Fill((genNeutrino.Px()-neutrinoKinFit[HighestProbComb].Px()),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D); 
	
	MSPlot["MetPy"]->Fill(MET.Py(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);  
	// MSPlot["NeutrinoKinFitPy"]->Fill(neutrinoKinFit[HighestProbComb].Py(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);  
	// //MSPlot["NeutrinoPyDiff"]->Fill((genNeutrino.Py()-neutrinoKinFit[HighestProbComb].Py()),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);
	
	// MSPlot["NeutrinoKinFitPz"]->Fill(neutrinoKinFit[HighestProbComb].Pz(),datasets[iDataSet],true, Luminosity*scaleFactor*lumiWeight3D);  
	//}
	
	//-------------------------------------------------------------------------------
	//Create some alternative helicites weights to obtain different ratios:
	//-------------------------------------------------------------------------------
	float LongitudinalFraction = 0.645167;
	float LeftHandedFraction = 0.321369;
	float RightHandedFraction = 0.033464;
	float TheoreticalDistributionValue = (LongitudinalFraction*6*(1-CosThetaGenerator*CosThetaGenerator) + (1-CosThetaGenerator)*(1-CosThetaGenerator)*3*LeftHandedFraction + RightHandedFraction*3*(1+CosThetaGenerator)*(1+CosThetaGenerator))/8;
	for(int helicityNumbers=0;helicityNumbers<SizeArray;helicityNumbers++){
	  if((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ){ 
	    UsedDistributionValue[helicityNumbers]=(HelicityFraction[helicityNumbers][0]*6*(1-CosThetaGenerator*CosThetaGenerator) + (1-CosThetaGenerator)*(1-CosThetaGenerator)*3*HelicityFraction[helicityNumbers][2] + HelicityFraction[helicityNumbers][1]*3*(1+CosThetaGenerator)*(1+CosThetaGenerator))/8;
	    HelicityWeight[helicityNumbers]=UsedDistributionValue[helicityNumbers]/TheoreticalDistributionValue;
	  }
	}

	//-------------------------------------------------------------
	//Obtain jet combination for the different b-tag constraints:
	//------------------------------------------------------------  		
	//Initialize bTag loop variables:
	int TCHEbTagLoop, TCHPbTagLoop, SSVHEbTagLoop, SSVHPbTagLoop, CSVbTagLoop,ConsideredBTagger=0; //0=tche, 1 = tchp, 2 = ssvhe, 3 = ssvhp & 4 = csv	
	MSPlot["nPVBeforeFoundJetComb"]->Fill(nPrimaryVertices,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	
	for(CSVbTagLoop =1;CSVbTagLoop< (NumberCSVbTags+1); CSVbTagLoop++){	  
	  for(SSVHPbTagLoop=1;SSVHPbTagLoop <=(NumberSSVHPbTags+1); SSVHPbTagLoop++){
	    if(ConsideredBTagger == 4) SSVHPbTagLoop =14;   //--> Do not run over all possible TCHE bTag values!	    
	    for(SSVHEbTagLoop=1; SSVHEbTagLoop<=(NumberSSVHEbTags+1); SSVHEbTagLoop++){
	      if(ConsideredBTagger == 3 || ConsideredBTagger == 4) SSVHEbTagLoop = 14;	      
	      for(TCHPbTagLoop=1; TCHPbTagLoop <= (NumberTCHPbTags+1); TCHPbTagLoop++){
		if(ConsideredBTagger == 2 || ConsideredBTagger == 3 || ConsideredBTagger == 4) TCHPbTagLoop =14; 		
		for(TCHEbTagLoop=1; TCHEbTagLoop <= (NumberTCHEbTags+1);TCHEbTagLoop++){
		  if(ConsideredBTagger ==1 || ConsideredBTagger ==2 || ConsideredBTagger ==3 || ConsideredBTagger ==4) TCHEbTagLoop =14;
		  int SumBTag = TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1;
		  if(UsedTCHE[SumBTag]==TCHEbTagLoop && UsedTCHP[SumBTag]==TCHPbTagLoop && UsedSSVHE[SumBTag]==SSVHEbTagLoop && UsedSSVHP[SumBTag]==SSVHPbTagLoop && UsedCSV[SumBTag]==CSVbTagLoop ){
		    int bTagLoop;
		    if(ConsideredBTagger == 0) bTagLoop = TCHEbTagLoop;
		    else if(ConsideredBTagger == 1) bTagLoop = TCHPbTagLoop;
		    else if(ConsideredBTagger == 2) bTagLoop = SSVHEbTagLoop;
		    else if(ConsideredBTagger == 3) bTagLoop = SSVHPbTagLoop;
		    else if(ConsideredBTagger == 4) bTagLoop = CSVbTagLoop;

		    int JetCombHadr = bTagJetSelection.HighestProbSelection(bTagLoop,ConsideredBTagger,KinFitProbHadr,btagTCHE,btagTCHP,btagSSVHE,btagSSVHP,btagCSV);
		    int JetCombHadrAndLeptWOnly =bTagJetSelection.HighestProbSelection(bTagLoop,ConsideredBTagger,KinFitProbHadrAndLeptWOnly,btagTCHE,btagTCHP,btagSSVHE,btagSSVHP,btagCSV);
		    int JetCombHadrAndLept = bTagJetSelection.HighestProbSelection(bTagLoop,ConsideredBTagger,KinFitProbHadrAndLept,btagTCHE,btagTCHP,btagSSVHE,btagSSVHP,btagCSV);

		    if(bTagLoop ==1){
		      histo1D["JetCombHadr"]->Fill(JetCombHadr);
		      histo1D["JetCombHadrAndLeptWOnly"]->Fill(JetCombHadrAndLeptWOnly);
		      histo1D["JetCombHadrAndLept"]->Fill(JetCombHadrAndLept);
		    }

		    // Only hadronic KinFit configuration:
		    if(JetCombHadr!=999){
		      if(dataSetName.find("Data") == 0) float MassTop = 173.1;
		      else float MassTop = 172.5;
		      float CosThetaCalcHadr = bTagCosThetaCalculation.CalcOrigKins(BLeptIndex[JetCombHadr], BHadrIndex[JetCombHadr],  lepton, selectedJets, 80.4, MassTop);
		      if(CosThetaCalcHadr != 999){
			TLorentzVector neutrinoHadr = bTagCosThetaCalculation.GetNeutrino();
			if((dataSetName.find("TTbarJets_SemiMu") == 0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") == 0 && semiElectron == true)){
			  histo1D["LeptWMassHadr"]->Fill((neutrinoHadr+lepton).M());
			  histo1D["LeptTopMassHadr"]->Fill((neutrinoHadr+lepton+selectedJets[BLeptIndex[JetCombHadr]]).M());
			}		      

			NumberRemainingEventsKinFitHadr[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1][iDataSet]++;
			CosThetaValuesKinFitHadr[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1].push_back(CosThetaCalcHadr);
			LumiWeightVectorKinFitHadr[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1].push_back(lumiWeight3D);
			if(((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ) && dataSetName.find("JES") != 0){
			  if(BLeptIndex[JetCombHadr] == CorrectBLept){//Count the number of events with correctly reconstructed leptonic b-jet
			    NumberBLeptCorrectEventsKinFitHadr[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1][iDataSet]++;
			  }
			  CosThGenKinFitHadr[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1].push_back(CosThetaGenerator);
			  EventCorrectionWeightKinFitHadr[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1].push_back(scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor());		      		      
			}//End of semiMu sample
		      }
		    }

		    //Hadronic and leptonic W Only KinFit configuration:
		    if(JetCombHadrAndLeptWOnly!=999){
		      float CosThetaCalcHadrAndLeptWOnly=bTagCosThetaCalculation.Calculation(leptonChangedHadrAndLeptWOnly[JetCombHadrAndLeptWOnly],neutrinoChangedHadrAndLeptWOnly[JetCombHadrAndLeptWOnly],selectedJets[BLeptIndex[JetCombHadrAndLeptWOnly]]);
		      if(CosThetaCalcHadrAndLeptWOnly != 999){		      		      
			NumberRemainingEventsKinFitHadrAndLeptWOnly[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1][iDataSet]++;
			CosThetaValuesKinFitHadrAndLeptWOnly[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1].push_back(CosThetaCalcHadrAndLeptWOnly);
			LumiWeightVectorKinFitHadrAndLeptWOnly[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1].push_back(lumiWeight3D);
			if(((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ) && dataSetName.find("JES") != 0){
			  if(BLeptIndex[JetCombHadrAndLeptWOnly] == CorrectBLept){//Count the number of events with correctly reconstructed leptonic b-jet
			    NumberBLeptCorrectEventsKinFitHadrAndLeptWOnly[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1][iDataSet]++;
			  }
			  CosThGenKinFitHadrAndLeptWOnly[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1].push_back(CosThetaGenerator);
			  EventCorrectionWeightKinFitHadrAndLeptWOnly[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1].push_back(scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor());		      		      
			}//End of semiMu sample
		      }
		    }
		    
		    //Hadronic and leptonic KinFit configuration:
		    if(JetCombHadrAndLept!=999){
		      float CosThetaCalcHadrAndLept=bTagCosThetaCalculation.Calculation(leptonChangedHadrAndLept[JetCombHadrAndLept],neutrinoChangedHadrAndLept[JetCombHadrAndLept],leptBChangedHadrAndLept[JetCombHadrAndLept]);
		      if(CosThetaCalcHadrAndLept != 999){
			NumberRemainingEventsKinFitHadrAndLept[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1][iDataSet]++;
			CosThetaValuesKinFitHadrAndLept[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1].push_back(CosThetaCalcHadrAndLept);
			LumiWeightVectorKinFitHadrAndLept[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1].push_back(lumiWeight3D);
			if(((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ) && dataSetName.find("JES") != 0){
			  if(BLeptIndex[JetCombHadrAndLept] == CorrectBLept){//Count the number of events with correctly reconstructed leptonic b-jet
			    NumberBLeptCorrectEventsKinFitHadrAndLept[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1][iDataSet]++;
			  }
			  CosThGenKinFitHadrAndLept[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1].push_back(CosThetaGenerator);
			  EventCorrectionWeightKinFitHadrAndLept[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1].push_back(scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor());		      		      
			}//End of semiMu sample
		      }
		    }
		  }
		  else{
		    cout << " Looking at wrong bTagging combination in calculation loops " << endl;
		  }
		  
		  if(TCHEbTagLoop == NumberTCHEbTags){
		    ConsideredBTagger = 1;
		    TCHPbTagLoop = 2;
		    SSVHEbTagLoop = 2;
		    SSVHPbTagLoop = 2;
		    CSVbTagLoop = 2;
		  }
		}//end of TCHE		
		if(TCHPbTagLoop == NumberTCHPbTags) ConsideredBTagger = 2;
	      }//end of TCHP      
	      if(SSVHEbTagLoop == NumberSSVHEbTags) ConsideredBTagger = 3;
	    }//end of SSVHE
	    if(SSVHPbTagLoop == NumberSSVHPbTags) ConsideredBTagger = 4;
	  }//end of SSVHP
	}//end of CSV
	      		
      }// End of loop over selected events
      
    } // end loop over events in wTrees    
    
    std::cout << "  " << endl;
    std::cout << " size of cos theta : " << endl;
    std::cout << "            - Hadronic KinFit only                : " << CosThetaValuesKinFitHadr[0].size() << endl;
    std::cout << "            - Hadronic KinFit and leptonic W only : " << CosThetaValuesKinFitHadrAndLeptWOnly[0].size() << endl;
    std::cout << "            - Hadronic KinFit and leptonic W      : " << CosThetaValuesKinFitHadrAndLept[0].size() << endl;
    std::cout << " °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°° " << endl;

    if((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ){
      cout << " " << endl;
      cout << " -----------------------------------------------------------------------------------------------------------------------" << endl;
      cout << " Performing helicity Generator fit : " << endl;
      cout << " ------------------------------------" << endl;
      TF1 *helicityFit = new TF1("helicityFit","[0]*((((1-[1]-[2])*3*(1+x)*(1+x))+([1]*3*(1-x)*(1-x))+([2]*6*(1-x*x)))/8)",-1,1);
      histo1D["StandardCosThetaFit"]->Fit("helicityFit","Q");
      std::cout << " fit values (before event selection) : Norm =" <<helicityFit->GetParameter(0) << " , Left = " << helicityFit->GetParameter(1) << " Long = " << helicityFit->GetParameter(2) << " ==> Right = " << 1-(helicityFit->GetParameter(1))-(helicityFit->GetParameter(2))<< std::endl;
      std::cout << " fit values error (before event selection) : " << helicityFit->GetParError(0) << " " << helicityFit->GetParError(1) << " " << helicityFit->GetParError(2) << std::endl;
      cout << "                      ------------------------------------" << endl;
      histo1D["StandardCosThetaFit"]->Scale(100./(histo1D["StandardCosThetaFit"]->Integral()));
      histo1D["StandardCosThetaFit"]->SetMinimum(0);
      histo1D["StandardCosThetaFit"]->SetMaximum(0.8);
      TF1 *helicityFit2 = new TF1("helicityFit2","((([0]*3*(1+x)*(1+x))+([1]*3*(1-x)*(1-x))+([2]*6*(1-x*x)))/8)",-1,1);
      histo1D["StandardCosThetaFit"]->Fit("helicityFit2","Q");
      std::cout << " fit values 2 (before event selection) : " << helicityFit2->GetParameter(0) << " " << helicityFit2->GetParameter(1) << " " << helicityFit2->GetParameter(2) << std::endl;
      std::cout << " fit values error 2 (before event selection) : " << helicityFit2->GetParError(0) << " " << helicityFit2->GetParError(1) << " " << helicityFit2->GetParError(2) << std::endl;
      cout << " -----------------------------------------------------------------------------------------------------------------------" << endl;    
    }

    //Execute MinuitFitter:
    int SumbTag, ConsideredTagger, CSV, SSVHP, SSVHE, TCHP, TCHE;
    for(CSV =0;CSV< NumberCSVbTags; CSV++){
      for(SSVHP=0;SSVHP <=NumberSSVHPbTags; SSVHP++){
	if(ConsideredTagger == 4) SSVHP =13;   //--> Do not run over all possible TCHE bTag values!	    
	for(SSVHE=0;SSVHE<=NumberSSVHEbTags; SSVHE++){
	  if(ConsideredTagger == 3 || ConsideredTagger == 4) SSVHE = 13;
	  for(TCHP=0; TCHP <= NumberTCHPbTags; TCHP++){
	    if(ConsideredTagger == 2 || ConsideredTagger == 3 || ConsideredTagger == 4) TCHP =13; 
	    for(TCHE=0; TCHE<= NumberTCHEbTags;TCHE++){
	      if(ConsideredTagger ==1 || ConsideredTagger ==2 || ConsideredTagger ==3 || ConsideredTagger ==4) TCHE =13;
	      SumbTag = TCHE+TCHP+SSVHE+SSVHP+CSV;
	      if(UsedTCHE[SumbTag]==(TCHE+1) && UsedTCHP[SumbTag]==(TCHP+1) && UsedSSVHE[SumbTag]==(SSVHE+1) && UsedSSVHP[SumbTag]==(SSVHP+1) && UsedCSV[SumbTag]==(CSV+1) ){		
		
		std::string CosThetaDataStringHadr = "CosThetaDataHadr_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV];
		std::string CosThetaSignalStringHadr = "CosThetaSignalHadr_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV];
		std::string CosThetaBckgStringHadr = "CosThetaBckgHadr_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV];
		
		std::string CosThetaDataStringHadrAndLeptW = "CosThetaDataHadrAndLeptWOnly_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV];
		std::string CosThetaSignalStringHadrAndLeptW="CosThetaSignalHadrAndLeptWOnly_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV];
		std::string CosThetaBckgStringHadrAndLeptW = "CosThetaBckgHadrAndLeptWOnly_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV];
		
		std::string CosThetaDataStringHadrAndLept = "CosThetaDataHadrAndLept_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV];
		std::string CosThetaSignalStringHadrAndLept = "CosThetaSignalHadrAndLept_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV];
		std::string CosThetaBckgStringHadrAndLept = "CosThetaBckgHadrAndLept_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV];
	      
		if(iDataSet == 0){
		  histo1D[CosThetaDataStringHadr]=new TH1F(CosThetaDataStringHadr.c_str(),CosThetaDataStringHadr.c_str(),CosThetaBinNumber,-1,1);
		  histo1D[CosThetaDataStringHadr]->SetDirectory(0);
		  histo1D[CosThetaSignalStringHadr]=new TH1F(CosThetaSignalStringHadr.c_str(),CosThetaSignalStringHadr.c_str(),CosThetaBinNumber,-1,1);
		  histo1D[CosThetaSignalStringHadr]->SetDirectory(0);
		  histo1D[CosThetaBckgStringHadr]=new TH1F(CosThetaBckgStringHadr.c_str(),CosThetaBckgStringHadr.c_str(),CosThetaBinNumber,-1,1);
		  histo1D[CosThetaBckgStringHadr]->SetDirectory(0);

		  histo1D[CosThetaDataStringHadrAndLeptW]=new TH1F(CosThetaDataStringHadrAndLeptW.c_str(),CosThetaDataStringHadrAndLeptW.c_str(),CosThetaBinNumber,-1,1);
		  histo1D[CosThetaDataStringHadrAndLeptW]->SetDirectory(0);
		  histo1D[CosThetaSignalStringHadrAndLeptW]=new TH1F(CosThetaSignalStringHadrAndLeptW.c_str(),CosThetaSignalStringHadrAndLeptW.c_str(),CosThetaBinNumber,-1,1);
		  histo1D[CosThetaSignalStringHadrAndLeptW]->SetDirectory(0);
		  histo1D[CosThetaBckgStringHadrAndLeptW]=new TH1F(CosThetaBckgStringHadrAndLeptW.c_str(),CosThetaBckgStringHadrAndLeptW.c_str(),CosThetaBinNumber,-1,1);
		  histo1D[CosThetaBckgStringHadrAndLeptW]->SetDirectory(0);

		  histo1D[CosThetaDataStringHadrAndLept]=new TH1F(CosThetaDataStringHadrAndLept.c_str(),CosThetaDataStringHadrAndLept.c_str(),CosThetaBinNumber,-1,1);
		  histo1D[CosThetaDataStringHadrAndLept]->SetDirectory(0);
		  histo1D[CosThetaSignalStringHadrAndLept]=new TH1F(CosThetaSignalStringHadrAndLept.c_str(),CosThetaSignalStringHadrAndLept.c_str(),CosThetaBinNumber,-1,1);
		  histo1D[CosThetaSignalStringHadrAndLept]->SetDirectory(0);
		  histo1D[CosThetaBckgStringHadrAndLept]=new TH1F(CosThetaBckgStringHadrAndLept.c_str(),CosThetaBckgStringHadrAndLept.c_str(),CosThetaBinNumber,-1,1);
		  histo1D[CosThetaBckgStringHadrAndLept]->SetDirectory(0);
		}
	      
		if((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ){//Defining of the genttbar histo
		  char hisnameHadr[100];
		  char hisnameHadrAndLeptWOnly[100];
		  char hisnameHadrAndLept[100];
		  sprintf(hisnameHadr,"CosThetaGenHadr_TCHE%s_TCHP%s_SSVHE%s_SSVHP%s_CSV%s", bTag[TCHE].c_str(),bTag[TCHP].c_str(),bTag[SSVHE].c_str(),bTag[SSVHP].c_str(),bTag[CSV].c_str());
		  sprintf(hisnameHadrAndLeptWOnly,"CosThetaGenHadrAndLeptWOnly_TCHE%s_TCHP%s_SSVHE%s_SSVHP%s_CSV%s", bTag[TCHE].c_str(),bTag[TCHP].c_str(),bTag[SSVHE].c_str(),bTag[SSVHP].c_str(),bTag[CSV].c_str());
		  sprintf(hisnameHadrAndLept,"CosThetaGenHadrAndLept_TCHE%s_TCHP%s_SSVHE%s_SSVHP%s_CSV%s", bTag[TCHE].c_str(),bTag[TCHP].c_str(),bTag[SSVHE].c_str(),bTag[SSVHP].c_str(),bTag[CSV].c_str());
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

		//Hadronic KinFit configuration: Filling of histograms		  		
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
		  
		  //Change data to systematics since then nominal values will be reweighted!!
		  if(WSystResults == true && dataSetName.find("WJets_Nominal") != 0)  //WJets syst
		    histo1D[CosThetaDataStringHadr]->Fill(CosThetaValuesKinFitHadr[SumbTag][ii],scaleFactor*(LumiWeightVectorKinFitHadr[SumbTag][ii])*Luminosity*dataSet->NormFactor());
		  if(WSystResults == false && (dataSetName.find("Data") == 0 || dataSetName.find("JES") == 0)) //Nominal and JES systematics
		    histo1D[CosThetaDataStringHadr]->Fill(CosThetaValuesKinFitHadr[SumbTag][ii],scaleFactor*(LumiWeightVectorKinFitHadr[SumbTag][ii])*Luminosity*dataSet->NormFactor());
		  if(SignalOnly == true && (dataSetName.find("TTbarJets_SemiMu") == 0))
		    histo1D[CosThetaDataStringHadr]->Fill(CosThetaValuesKinFitHadr[SumbTag][ii],scaleFactor*(LumiWeightVectorKinFitHadr[SumbTag][ii])*Luminosity*dataSet->NormFactor());  
		  
		  if(((dataSetName.find("TTbarJets_SemiMu")==0 && semiMuon==true) || (dataSetName.find("TTbarJets_SemiEl")==0 && semiElectron == true) ) && dataSetName.find("JES_") != 0)
		    histo1D[CosThetaSignalStringHadr]->Fill(CosThetaValuesKinFitHadr[SumbTag][ii],Luminosity*scaleFactor*(LumiWeightVectorKinFitHadr[SumbTag][ii])*NominalNormFactor);
		  else if(dataSetName.find("Data") !=0 && dataSetName.find("JES") != 0 && dataSetName.find("Syst") != 0) 
		    histo1D[CosThetaBckgStringHadr]->Fill(CosThetaValuesKinFitHadr[SumbTag][ii],Luminosity*scaleFactor*(LumiWeightVectorKinFitHadr[SumbTag][ii])*NominalNormFactor);
		}

		//Hadronic and Leptonic W Only KinFit Configuration: Filling of histograms
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
		  
		  //Change data to systematics since then nominal values will be reweighted!!
		  if(WSystResults == true && dataSetName.find("WJets_Nominal") != 0)  //WJets syst
		    histo1D[CosThetaDataStringHadrAndLeptW]->Fill(CosThetaValuesKinFitHadrAndLeptWOnly[SumbTag][ii],scaleFactor*(LumiWeightVectorKinFitHadrAndLeptWOnly[SumbTag][ii])*Luminosity*dataSet->NormFactor());
		  if(WSystResults == false && (dataSetName.find("Data") == 0 || dataSetName.find("JES") == 0)) //Nominal and JES systematics
		    histo1D[CosThetaDataStringHadrAndLeptW]->Fill(CosThetaValuesKinFitHadrAndLeptWOnly[SumbTag][ii],scaleFactor*(LumiWeightVectorKinFitHadrAndLeptWOnly[SumbTag][ii])*Luminosity*dataSet->NormFactor());
		  if(SignalOnly == true && (dataSetName.find("TTbarJets_SemiMu") == 0))
		    histo1D[CosThetaDataStringHadrAndLeptW]->Fill(CosThetaValuesKinFitHadrAndLeptWOnly[SumbTag][ii],scaleFactor*(LumiWeightVectorKinFitHadrAndLeptWOnly[SumbTag][ii])*Luminosity*dataSet->NormFactor());  
		  
		  if(((dataSetName.find("TTbarJets_SemiMu")==0 && semiMuon==true) || (dataSetName.find("TTbarJets_SemiEl")==0 && semiElectron == true) ) && dataSetName.find("JES_") != 0)
		    histo1D[CosThetaSignalStringHadrAndLeptW]->Fill(CosThetaValuesKinFitHadrAndLeptWOnly[SumbTag][ii],Luminosity*scaleFactor*(LumiWeightVectorKinFitHadrAndLeptWOnly[SumbTag][ii])*NominalNormFactor);
		  else if(dataSetName.find("Data") !=0 && dataSetName.find("JES") != 0 && dataSetName.find("Syst") != 0) 
		    histo1D[CosThetaBckgStringHadrAndLeptW]->Fill(CosThetaValuesKinFitHadrAndLeptWOnly[SumbTag][ii],Luminosity*scaleFactor*(LumiWeightVectorKinFitHadrAndLeptWOnly[SumbTag][ii])*NominalNormFactor);
		}

		//Hadronic and Leptonic KinFit Configuration: Filling of histograms
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
		  
		  //Change data to systematics since then nominal values will be reweighted!!
		  if(WSystResults == true && dataSetName.find("WJets_Nominal") != 0)  //WJets syst
		    histo1D[CosThetaDataStringHadrAndLept]->Fill(CosThetaValuesKinFitHadrAndLept[SumbTag][ii],scaleFactor*(LumiWeightVectorKinFitHadrAndLept[SumbTag][ii])*Luminosity*dataSet->NormFactor());
		  if(WSystResults == false && (dataSetName.find("Data") == 0 || dataSetName.find("JES") == 0)) //Nominal and JES systematics
		    histo1D[CosThetaDataStringHadrAndLept]->Fill(CosThetaValuesKinFitHadrAndLept[SumbTag][ii],scaleFactor*(LumiWeightVectorKinFitHadrAndLept[SumbTag][ii])*Luminosity*dataSet->NormFactor());
		  if(SignalOnly == true && (dataSetName.find("TTbarJets_SemiMu") == 0))
		    histo1D[CosThetaDataStringHadrAndLept]->Fill(CosThetaValuesKinFitHadrAndLept[SumbTag][ii],scaleFactor*(LumiWeightVectorKinFitHadrAndLept[SumbTag][ii])*Luminosity*dataSet->NormFactor());  
		  
		  if(((dataSetName.find("TTbarJets_SemiMu")==0 && semiMuon==true) || (dataSetName.find("TTbarJets_SemiEl")==0 && semiElectron == true) ) && dataSetName.find("JES_") != 0)
		    histo1D[CosThetaSignalStringHadrAndLept]->Fill(CosThetaValuesKinFitHadrAndLept[SumbTag][ii],Luminosity*scaleFactor*(LumiWeightVectorKinFitHadrAndLept[SumbTag][ii])*NominalNormFactor);
		  else if(dataSetName.find("Data") !=0 && dataSetName.find("JES") != 0 && dataSetName.find("Syst") != 0) 
		    histo1D[CosThetaBckgStringHadrAndLept]->Fill(CosThetaValuesKinFitHadrAndLept[SumbTag][ii],Luminosity*scaleFactor*(LumiWeightVectorKinFitHadrAndLept[SumbTag][ii])*NominalNormFactor);
		}


		if(iDataSet==(datasets.size()-1)){//Go in this loop when the last datasample is active to perform MinuitFitter
		  
		  //Hadronic KinFit configuration: MinuitFitter		  		
		  MinuitFitter minuitFitterHadr = MinuitFitter(histo1D[CosThetaDataStringHadr], histo1D[CosThetaSignalStringHadr], histo1D[CosThetaBckgStringHadr], 0.6671, 0.3325, 0.0004,genttbarhistoHadr,ndimen);			  
		  PresentationTexHadr << PresentationOutput[SumbTag] << " & ";
		  cout << " Inside minuitFitter loop : " << PresentationOutput[SumbTag] << " & " << endl;
		  for(int ii=0; ii< nameDataSet.size(); ii++){
		    //if(nameDataSet[ii].find("JES") != 0 && nameDataSet[ii].find("Syst") != 0 && ii < nameDataSet.size()-1){ //Only output for samples with no systematics
		    if(ii < nameDataSet.size()-1){ PresentationTexHadr << NumberRemainingEventsKinFitHadr[SumbTag][ii] << " & ";}
		    else if(ii == nameDataSet.size()-1 ){ 
		      PresentationTexHadr << NumberRemainingEventsKinFitHadr[SumbTag][ii] << " & ";//For presentation helicity values still need to be included !!
		      PresentationTexHadr << NumberBLeptCorrectEventsKinFitHadr[SumbTag][ii] << " & ";
		    }
		  }		  
		  
		  PresentationTexHadr << minuitFitterHadr.GetFRResult() << " & " << minuitFitterHadr.GetFRResult()-SMfrResult <<" & " << minuitFitterHadr.GetFRError() << " & " << minuitFitterHadr.GetFLResult() << " & " << minuitFitterHadr.GetFLResult()-SMflResult << " & " << minuitFitterHadr.GetFLError() << " & " << minuitFitterHadr.GetF0Result() << " & " << minuitFitterHadr.GetF0Result()-SMf0Result << " & " << minuitFitterHadr.GetF0Error() << " \\\\ " << endl;		
		  PresentationTexHadr << " \\hline " << endl;

		  //Hadronic and Leptonic W Only KinFit configuration: MinuitFitter		  		
		  MinuitFitter minuitFitterHadrAndLeptWOnly = MinuitFitter(histo1D[CosThetaDataStringHadrAndLeptW], histo1D[CosThetaSignalStringHadrAndLeptW], histo1D[CosThetaBckgStringHadrAndLeptW], 0.6671, 0.3325, 0.0004,genttbarhistoHadrAndLeptWOnly,ndimen);			  
		  PresentationTexHadrAndLeptWOnly << PresentationOutput[SumbTag] << " & ";
		  for(int ii=0; ii< nameDataSet.size(); ii++){
		    //if(nameDataSet[ii].find("JES") != 0 && nameDataSet[ii].find("Syst") != 0 && ii < nameDataSet.size()-1){ //Only output for samples with no systematics
		    if(ii < nameDataSet.size()-1){ PresentationTexHadrAndLeptWOnly << NumberRemainingEventsKinFitHadrAndLeptWOnly[SumbTag][ii] << " & ";}
		    else if(ii == nameDataSet.size()-1 ){ 
		      PresentationTexHadrAndLeptWOnly << NumberRemainingEventsKinFitHadrAndLeptWOnly[SumbTag][ii] << " & ";//For presentation helicity values still need to be included !!
		      PresentationTexHadrAndLeptWOnly << NumberBLeptCorrectEventsKinFitHadrAndLeptWOnly[SumbTag][ii] << " & ";
		    }
		  }		  
		  
		  PresentationTexHadrAndLeptWOnly << minuitFitterHadrAndLeptWOnly.GetFRResult() << " & " << minuitFitterHadrAndLeptWOnly.GetFRResult()-SMfrResult <<" & " << minuitFitterHadrAndLeptWOnly.GetFRError() << " & " << minuitFitterHadrAndLeptWOnly.GetFLResult() << " & " << minuitFitterHadrAndLeptWOnly.GetFLResult()-SMflResult << " & " << minuitFitterHadrAndLeptWOnly.GetFLError() << " & " << minuitFitterHadrAndLeptWOnly.GetF0Result() << " & " << minuitFitterHadrAndLeptWOnly.GetF0Result()-SMf0Result << " & " << minuitFitterHadrAndLeptWOnly.GetF0Error() << " \\\\ " << endl;		
		  PresentationTexHadrAndLeptWOnly << " \\hline " << endl;

		  //Hadronic and Leptonic KinFit configuration: MinuitFitter		  		
		  MinuitFitter minuitFitterHadrAndLept = MinuitFitter(histo1D[CosThetaDataStringHadrAndLept], histo1D[CosThetaSignalStringHadrAndLept], histo1D[CosThetaBckgStringHadrAndLept], 0.6671, 0.3325, 0.0004,genttbarhistoHadrAndLept,ndimen);			  
		  PresentationTexHadrAndLept << PresentationOutput[SumbTag] << " & ";
		  for(int ii=0; ii< nameDataSet.size(); ii++){
		    //if(nameDataSet[ii].find("JES") != 0 && nameDataSet[ii].find("Syst") != 0 && ii < nameDataSet.size()-1){ //Only output for samples with no systematics
		    if(ii < nameDataSet.size()-1){ PresentationTexHadrAndLept << NumberRemainingEventsKinFitHadrAndLept[SumbTag][ii] << " & ";}
		    else if(ii == nameDataSet.size()-1 ){ 
		      PresentationTexHadrAndLept << NumberRemainingEventsKinFitHadrAndLept[SumbTag][ii] << " & ";//For presentation helicity values still need to be included !!
		      PresentationTexHadrAndLept << NumberBLeptCorrectEventsKinFitHadrAndLept[SumbTag][ii] << " & ";
		    }
		  }		  
		  
		  PresentationTexHadrAndLept << minuitFitterHadrAndLept.GetFRResult() << " & " << minuitFitterHadrAndLept.GetFRResult()-SMfrResult <<" & " << minuitFitterHadrAndLept.GetFRError() << " & " << minuitFitterHadrAndLept.GetFLResult() << " & " << minuitFitterHadrAndLept.GetFLResult()-SMflResult << " & " << minuitFitterHadrAndLept.GetFLError() << " & " << minuitFitterHadrAndLept.GetF0Result() << " & " << minuitFitterHadrAndLept.GetF0Result()-SMf0Result << " & " << minuitFitterHadrAndLept.GetF0Error() << " \\\\ " << endl;		
		  PresentationTexHadrAndLept << " \\hline " << endl;

		}//End of minuitFitter
		
	      }//end of wrong entry chosen failed requirement!
	      else cout << " Looking at wrong bTagging combination in fitting minuit !! " << endl;
	      
	      if(TCHE==NumberTCHEbTags-1) {
		ConsideredTagger =1;
		TCHP=1;
		SSVHE=1;
		SSVHP=1;
		CSV=1;
	      }	      
	    }//end of TCHE
	    if(TCHP == NumberTCHPbTags-1) ConsideredTagger =2;
	  }//end of TCHP
	  if(SSVHE == NumberSSVHEbTags-1) ConsideredTagger = 3;
	}//end of SSVHE
	if(SSVHP == NumberSSVHPbTags-1) ConsideredTagger = 4;
      }//end of SSVHP
    }//end of CSV
        
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
  cout << " NumberEventsAfterbTag : " <<  NumberEventsAfterbTag << endl;
  cout << " NumberEventsAfterTransverseMass : " << NumberEventsAfterTransverseMass  <<endl;
  cout << " NumberSelectedEvents : " <<  NumberSelectedEvents << endl;
  cout << " NumberUsedEvents : " <<  NumberUsedEvents << endl;
  cout << " NumberUsedCorrectEvents : " << NumberUsedCorrectEvents << endl;
  cout << " NumberUsedWrongEvents : " <<  NumberUsedWrongEvents  << endl;
  cout << "                -------        " << endl;
  cout << " NumberDataEventsBeforeCuts : " << NumberDataEventsBeforeCuts << endl;
  cout << " NumberSelectedDataEvents : " << NumberSelectedDataEvents  << endl;
  cout << " NumberUsedDataEvents : " << NumberUsedDataEvents  << endl;
  std::cout << " " << endl;

  //Hadronic:
  PresentationTexHadr << " \\end{tabular} " << endl;
  PresentationTexHadr << " \\end{center} " << endl;
  PresentationTexHadr << " \\renewcommand{\\arraystretch}{0.9} " << endl;
  PresentationTexHadr << " \\end{tiny} " << endl;
  PresentationTexHadr << " \\end{table} " << endl;
  PresentationTexHadr.close();

  //Hadronic and leptonic W only:
  PresentationTexHadrAndLeptWOnly << " \\end{tabular} " << endl;
  PresentationTexHadrAndLeptWOnly << " \\end{center} " << endl;
  PresentationTexHadrAndLeptWOnly << " \\renewcommand{\\arraystretch}{0.9} " << endl;
  PresentationTexHadrAndLeptWOnly << " \\end{tiny} " << endl;
  PresentationTexHadrAndLeptWOnly << " \\end{table} " << endl;
  PresentationTexHadrAndLeptWOnly.close();
  
  //Hadronic and leptonic:
  PresentationTexHadrAndLept << " \\end{tabular} " << endl;
  PresentationTexHadrAndLept << " \\end{center} " << endl;
  PresentationTexHadrAndLept << " \\renewcommand{\\arraystretch}{0.9} " << endl;
  PresentationTexHadrAndLept << " \\end{tiny} " << endl;
  PresentationTexHadrAndLept << " \\end{table} " << endl;
  PresentationTexHadrAndLept.close();

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
