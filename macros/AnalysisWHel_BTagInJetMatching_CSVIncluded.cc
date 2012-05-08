
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
#include "../WHelicities/interface/BTagName_CSV.h"
#include "../WHelicities/interface/BTagJetSelection_CSV.h"
#include "../WHelicities/interface/BTagCosThetaCalculation.h"
#include "../WHelicities/interface/MinuitFitter.h"

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
TNtuple *genttbarhisto[CosThetaBinNumber];  //This is the vector of ntuples containing the generated values of cos theta* for each cos theta* reconstructed bin

//new version received from Mara
class MyChi2 : public ROOT::Math::IMultiGenFunction {
public:
  // Mandatory functions 
  
  virtual MyChi2 * Clone() const {return new MyChi2;};
  virtual unsigned int NDim() const {return ndimen;}
  
  // Constructors
  MyChi2(){};
  MyChi2(TH1F* & datah, TH1F* &signalh, TH1F* & bkgh,double ff0, double ffl, double ffr) : datah_(datah),signalh_(signalh),bkgh_(bkgh),ff0_(ff0),ffl_(ffl),ffr_(ffr){};  

  double DoEval(const double* x) const {


    double f0Fit = x[0];
    double flFit = x[1];
    double frFit = 0;
    double nn=1.;
    if (ndimen==3) { 
      frFit = 1-x[0]-x[1]; 
      nn = x[2];
    }
    else if (ndimen==2){
      flFit = 1.-x[0];
      nn = x[1];
    }

    double ChiSquaredAllBins =0.;
    double chi2 = 0.;
    
    int nncells = signalh_->GetSize()-2;  // -1 and not -2, to get the last bin

    // initialization
    double nmcaux[nncells];
    for (int ibin=0; ibin< nncells; ibin++){
      nmcaux[ibin]=0.;
    }
  
    // loop over ntuples containing gen-level costh info relative to each costhRec bin
    Int_t nGENentries[nncells];
    for (int ihgen=0; ihgen<nncells; ihgen++){ 
 
      // Here I am assuming that the input array only contains entries that 
      // have costheta* reconstructed
      // If no costhetarec for this igen, one should skip this entry in the sum
      // or send it to an underflow/overflow bin
      
      Float_t costhgen, evtweight;
      nGENentries[ihgen] = genttbarhisto[ihgen]->GetEntries();
      
      genttbarhisto[ihgen]->SetBranchAddress("costhgen",&costhgen);
      genttbarhisto[ihgen]->SetBranchAddress("evtweight",&evtweight);
    
      if (nGENentries[ihgen] ==0) {
	nmcaux[ihgen] =0;
	//unrwgt[ihgen]=0; //Double which has the function of summing up the number of generated events after normalization, but before reweighting in the fitter.
	// Double is necessary because cos theta * only needs to be reweighted for generated Muon leptonic decays
	// --> Not needed in my code since we can split in ttbarSemiMu ??
      }
      else{
        for (Int_t igen=0; igen<nGENentries[ihgen]; ++igen) {
          genttbarhisto[ihgen]->GetEntry(igen);
          double xx = costhgen;
          double thisevweight = evtweight;

          double SM = 0.3325*(1-xx)*(1-xx)*3/8+0.6671*(1-xx*xx)*3/4 + 0.0004*(1+xx)*(1+xx)*3/8 ;
          double newmodel = flFit*(1-xx)*(1-xx)*3/8+f0Fit*(1-xx*xx)*3/4 + frFit*(1+xx)*(1+xx)*3/8 ;
          
          nmcaux[ihgen] +=thisevweight*newmodel/SM;
        }
      }	  
    }
    for (int ibin=0; ibin< nncells; ibin++){
      double ndata_i = datah_->GetBinContent(ibin+1);       
      double nbkg_i = bkgh_->GetBinContent(ibin+1);     
      double  nmc_i =  nbkg_i +  nn*nmcaux[ibin] ;   // nn--> free normalization found by the Fit
      
      if (nmc_i>0){
	chi2 += ( nmc_i - ndata_i * TMath::Log(nmc_i) ); 
	ChiSquaredAllBins=ChiSquaredAllBins+(((nmc_i-ndata_i)/(sqrt(ndata_i)))*((nmc_i-ndata_i)/(sqrt(ndata_i))));//nmc_i = all samples ndata_ = semiMu	
      }      
    }    
    return ChiSquaredAllBins;
  }
private:
  
  TH1F* datah_;
  TH1F* signalh_;
  TH1F* bkgh_;
  double ff0_;
  double ffl_;
  double ffr_;
};

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
  
  vector<string> inputWTree;
  vector<string> nameDataSet;

  //--------------------------------------------------------------
  //     DataSamples needed for 2011 Data for muon channel:     --
  //--------------------------------------------------------------
  
  //Probability cut value:
  float KinFitCut = 0.;

  string decayChannel;
  bool semiMuon = false;
  bool semiElectron = true;
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
  bool SignalOnly = false;
  bool DataResults = true;
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
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_QCD_Pt-20to30_BCtoE_SemiEl.root").c_str());  //Use UsedTrigger string to make sure that sample is only used in correct case
      nameDataSet.push_back("QCD_Pt-20to30_BCtoE");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_QCD_Pt-20to30_EMEnriched_SemiEl.root").c_str());
      nameDataSet.push_back("QCD_Pt-20to30_EMEnriched");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_QCD_Pt-30to80_BCtoE_SemiEl.root").c_str());
      nameDataSet.push_back("QCD_Pt-30to80_BCtoE");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_QCD_Pt-30to80_EMEnriched_SemiEl.root").c_str());
      nameDataSet.push_back("QCD_Pt-30to80_EMEnriched");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_QCD_Pt-80to170_BCtoE_SemiEl.root").c_str());
      nameDataSet.push_back("QCD_Pt-80to170_BCtoE");
      inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_QCD_Pt-80to170_EMEnriched_SemiEl.root").c_str());
      nameDataSet.push_back("QCD_Pt-80to170_EMEnriched");
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
    inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_TTbarJets_SemiMuon_"+decayChannel+".root").c_str());
    nameDataSet.push_back("TTbarJets_SemiMuon");
  }
  else if(semiElectron == true){
    inputWTree.push_back(("WTree/KinFit_WTree_"+UsedTrigger+"_TTbarJets_SemiEl_"+decayChannel+".root").c_str());
    nameDataSet.push_back("TTbarJets_SemiEl");
  }
  
  TFile *fout = new TFile (("WHelicities_Analysis_"+UsedTrigger+".root").c_str(), "RECREATE");
  fout->cd();

  //Create .txt file to write away all information obtained for bTag study
  mkdir("BTagPerformanceStudy_CSV/",0777);
  //ofstream bTagFile("BTagPerformanceStudy/BTagInMatchingOutput.txt");
  //ofstream bTagFileTex("../../../../Documents/Vub/PhD/BTagPerformanceStudy/bTagInMatching.tex");

  //-----------------------------------------------
  //  Output .tex files for presentation/paper
  //-----------------------------------------------
  //PresentationTex.precision(2);

  //oooooOOOOOOOOOOOoooooooooooOOOOOOOOOOOOOOOOoooooooooOOOOOOOOOOOoooooooooooooooooOOOOOOOOOOOOOOOoooooooooooooOOOOOOOOOOOoooooooo
  //1) Own results: No btag at event selection, top mass constraint for leptonic top in KinFit and offline muon pt cut of 27
  //oooooOOOOOOOOOOOoooooooooooOOOOOOOOOOOOOOOOoooooooooOOOOOOOOOOOoooooooooooooooooOOOOOOOOOOOOOOOoooooooooooooOOOOOOOOOOOoooooooo

  string PresentationTexTitle;
  PresentationTexTitle =("BTagPerformanceStudy_CSV/"+UsedTrigger+"_"+dataSet+decayChannel+AppliedCuts+".tex").c_str();
  if(DataResults == true) PresentationTexTitle =("BTagPerformanceStudy_CSV/"+UsedTrigger+"_"+dataSet+decayChannel+AppliedCuts+".tex").c_str();
  if(JESMinusResults == true) PresentationTexTitle =("BTagPerformanceStudy_CSV/"+UsedTrigger+"_"+dataSet+decayChannel+AppliedCuts+"_JESMin.tex").c_str();
  if(JESPlusResults == true) PresentationTexTitle = ("BTagPerformanceStudy_CSV/"+UsedTrigger+"_"+dataSet+decayChannel+AppliedCuts+"_JESPlus.tex").c_str();
  if(WSystResults == true && WSystPositive == false) PresentationTexTitle = ("BTagPerformanceStudy_CSV/"+UsedTrigger+"_"+dataSet+decayChannel+AppliedCuts+"_WMin.tex").c_str();
  if(WSystResults == true && WSystPositive == true) PresentationTexTitle = ("BTagPerformanceStudy_CSV/"+UsedTrigger+"_"+dataSet+decayChannel+AppliedCuts+"_WPlus.tex").c_str();
  
  ofstream PresentationTex(PresentationTexTitle.c_str());

  //ofstream PresentationTex(("BTagPerformanceStudy_CSV/WorkshopPresentation"+UsedTrigger+"NoBTagEvtSelTopFittedProbCut.tex").c_str());
  //ofstream PresentationTex(("BTagPerformanceStudy_CSV/WorkshopPresentation"+UsedTrigger+"NoBTagEvtSelTopFittedMonteCarlo.tex").c_str());
  //ofstream PresentationTex(("BTagPerformanceStudy_CSV/WorkshopPresentation"+UsedTrigger+"NoBTagEvtSelTopFittedJESMinProbCut.tex").c_str());
  //ofstream PresentationTex(("BTagPerformanceStudy_CSV/WorkshopPresentation"+UsedTrigger+"NoBTagEvtSelTopFittedJESPlusProbCut.tex").c_str());
  //ofstream PresentationTex(("BTagPerformanceStudy_CSV/WorkshopPresentation"+UsedTrigger+"NoBTagEvtSelTopFittedWMin100.tex").c_str());
  //ofstream PresentationTex(("BTagPerformanceStudy_CSV/WorkshopPresentation"+UsedTrigger+"NoBTagEvtSelTopFittedWMinProbCut.tex").c_str());
  //ofstream PresentationTex(("BTagPerformanceStudy_CSV/WorkshopPresentation"+UsedTrigger+"NoBTagEvtSelTopFittedWPlus100.tex").c_str());
  //ofstream PresentationTex(("BTagPerformanceStudy_CSV/WorkshopPresentation"+UsedTrigger+"NoBTagEvtSelTopFittedWPlusProbCut.tex").c_str());

  //oooooOOOOOOOOOOOoooooooooooOOOOOOOOOOOOOOOOoooooooooOOOOOOOOOOOoooooooooooooooooOOOOOOOOOOOOOOOoooooooooooooOOOOOOOOOOOoooooooo
  //2) Ciemat results: SSVHEM btag at event selection, leptonic top mass left free in KinFit and offline muon pt cut of 25
  //oooooOOOOOOOOOOOoooooooooooOOOOOOOOOOOOOOOOoooooooooOOOOOOOOOOOoooooooooooooooooOOOOOOOOOOOOOOOoooooooooooooOOOOOOOOOOOoooooooo

  //ofstream PresentationTex(("BTagPerformanceStudy_CSV/CiematTest"+UsedTrigger+".tex").c_str());
  //ofstream PresentationTex(("BTagPerformanceStudy_CSV/CiematMonteCarlo"+UsedTrigger+".tex").c_str());
  //ofstream PresentationTex(("BTagPerformanceStudy_CSV/Ciemat"+UsedTrigger+"JESMin.tex").c_str());
  //ofstream PresentationTex(("BTagPerformanceStudy_CSV/Ciemat"+UsedTrigger+"JESPlus.tex").c_str());
  //ofstream PresentationTex(("BTagPerformanceStudy_CSV/Ciemat"+UsedTrigger+"WMin100.tex").c_str());
  //ofstream PresentationTex(("BTagPerformanceStudy_CSV/Ciemat"+UsedTrigger+"WPlus100.tex").c_str());

  //-----------------------------------------
  //  Start of filling of .tex files !!
  //-----------------------------------------
  PresentationTex << " \\begin{table} " << endl;
  PresentationTex << " \\begin{tiny} " << endl;
  PresentationTex << " \\renewcommand{\\arraystretch}{1.2} " << endl;
  PresentationTex << " \\begin{center} " << endl;
  for(int ii=0; ii<nameDataSet.size(); ii++){
    if(ii==0){PresentationTex << "  \\begin{tabular}{|c|";}
    else if(ii < (nameDataSet.size()-1)){PresentationTex << "c|";}
    else{PresentationTex << "c|c|c|c|c|c|c|c|c|c|c|} " << endl;}
  }
  PresentationTex << " \\hline " << endl;
  for(int ii=0; ii < nameDataSet.size();ii++){
    PresentationTex << " & " << nameDataSet[ii];
  }
  PresentationTex << " & bLept correct & F+ & F+ - SM & $\\delta$ F+ & F- & F- - SM & $\\delta$ F- & F0 & F0 - SM & $\\delta$ F0 \\\\ " << endl;
  PresentationTex << " \\hline " << endl;

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

  histo1D["HadronicWMass"] = new TH1F("HadronicWMass","HadronicWMass",50,50,100);
  histo1D["HadronicTopMass"] = new TH1F("HadronicTopMass","HadronicTopMass",50,150,200);
  histo1D["LeptonicWMass"] = new TH1F("LeptonicWMass","LeptonicWMass",50,50,100);
  histo1D["LeptonicTopMass"] = new TH1F("LeptonicTopMass","LeptonicTopMass",50,150,200); 

  histo1D["KinFitProbCorrectBLept"] = new TH1F("KinFitProbCorrectBLept","KinFitProbCorrectBLept",25,0,1);
  histo1D["KinFitProbWrongBLept"] = new TH1F("KinFitProbWrongBLept","KinFitProbWrongBLept",25,0,1);

  histo1D["KinFitProbCorrectBLeptSSVHEM"] = new TH1F("KinFitProbCorrectBLeptSSVHEM","KinFitProbCorrectBLeptSSVHEM",25,0,1);
  histo1D["KinFitProbWrongBLeptSSVHEM"] = new TH1F("KinFitProbWrongBLeptSSVHEM","KinFitProbWrongBLeptSSVHEM",25,0,1);

  histo1D["CosThetaCorrectBLept"]=new TH1F("CosThetaCorrectBLept","CosThetaCorrectBLept",CosThetaBinNumber,-1,1);
  histo1D["CosThetaWrongBLept"]=new TH1F("CosThetaWrongBLept","CosThetaWrongBLept",CosThetaBinNumber,-1,1);
  histo1D["CosThetaCorrectBLeptSSVHEM"]=new TH1F("CosThetaCorrectBLeptSSVHEM","CosThetaCorrectBLeptSSVHEM",CosThetaBinNumber,-1,1);
  histo1D["CosThetaWrongBLeptSSVHEM"]=new TH1F("CosThetaWrongBLeptSSVHEM","CosThetaWrongBLeptSSVHEM",CosThetaBinNumber,-1,1);
  
  histo1D["CosThetaRes"] = new TH1F("CosThetaRes","CosThetaRes",200,-2,2);
  histo1D["CosThetaDistributionSignal"] = new TH1F("CosThetaDistributionSignal","CosThetaDistributionSignal",100,-1,1);
  histo1D["CosThetaDistributionSignalRatio"] = new TH1F("CosThetaDistributionSignalRatio","CosThetaDistributionSignalRatio",100,-1,1);
  histo1D["CosThetaDistributionSignalGen"] = new TH1F("CosThetaDistributionSignalGen","CosThetaDistributionSignalGen",100,-1,1);
  histo1D["MuonPtSignal"] = new TH1F("MuonPtSignal","MuonPtSignal",100,0,200);
  TProfile *CosThetaDiffVSRecoSignal = new TProfile("CosThetaDiffVSRecoSignal","TProfile of Cos Theta* (reco-gen) vs Cos Theta* (reco) for muon Pt 20 to inf GeV",10,-1,1,-2,2);

  //Histos to determine optimal muon Pt cut:
  histo1D["CosThetaDifferencePt10ToInf"] = new TH1F("CosThetaDifferencePt10ToInf","CosThetaDifferencePt10ToInf",200,-2,2);
  histo1D["CosThetaRecoPt10ToInf"] = new TH1F("CosThetaRecoPt10ToInf","CosThetaRecoPt10ToInf",100,-1,1);
  histo1D["CosThetaRecoPt10ToInfForRatio"] = new TH1F("CosThetaRecoPt10ToInfForRatio","CosThetaRecoPt10ToInfForRatio",100,-1,1);
  histo1D["CosThetaGenPt10ToInf"] = new TH1F("CosThetaGenPt10ToInf","CosThetaGenPt10ToInf",100,-1,1);
  TProfile *CosThetaDiffVSRecoPt10ToInf = new TProfile("CosThetaDiffVSRecoPt10ToInf","TProfile of Cos Theta* (reco-gen) vs Cos Theta* (reco) for muon Pt 10 to inf GeV",10,-1,1,-2,2);
  histo2D["CosThetaGenVsRecoPt10ToInf"] = new TH2F("CosThetaGenVsRecoPt10ToInf","CosThetaGenVsRecoPt10ToInf",50,-1,1,20,-1,1);
  histo1D["CosThetaDiffPt10ToInfGenFirstBin"] = new TH1F("CosThetaDiffPt10ToInfGenFirstBin","CosThetaDiffPt10ToInfGenFirstBin",100,-2,2);
  histo1D["CosThetaDiffPt10ToInfGenSecondBin"] = new TH1F("CosThetaDiffPt10ToInfGenSecondBin","CosThetaDiffPt10ToInfGenSecondBin",100,-2,2);
  histo1D["CosThetaDiffPt10ToInfGenThirdBin"] = new TH1F("CosThetaDiffPt10ToInfGenThirdBin","CosThetaDiffPt10ToInfGenThirdBin",100,-2,2);
  histo1D["CosThetaDiffPt10ToInfGenFourthBin"] = new TH1F("CosThetaDiffPt10ToInfGenFourthBin","CosThetaDiffPt10ToInfGenFourthBin",100,-2,2);
  histo1D["CosThetaDiffPt10ToInfGenFifthBin"] = new TH1F("CosThetaDiffPt10ToInfGenFifthBin","CosThetaDiffPt10ToInfGenFifthBin",100,-2,2);
  histo1D["CosThetaDiffPt10ToInfGenSixthBin"] = new TH1F("CosThetaDiffPt10ToInfGenSixthBin","CosThetaDiffPt10ToInfGenSixthBin",100,-2,2);
  histo1D["CosThetaDiffPt10ToInfGenSeventhBin"] = new TH1F("CosThetaDiffPt10ToInfGenSeventhBin","CosThetaDiffPt10ToInfGenSeventhBin",100,-2,2);
  histo1D["CosThetaDiffPt10ToInfGenEighthBin"] = new TH1F("CosThetaDiffPt10ToInfGenEighthBin","CosThetaDiffPt10ToInfGenEighthBin",100,-2,2);
  histo1D["CosThetaDiffPt10ToInfGenNinthBin"] = new TH1F("CosThetaDiffPt10ToInfGenNinthBin","CosThetaDiffPt10ToInfGenNinthBin",100,-2,2);
  histo1D["CosThetaDiffPt10ToInfGenTenthBin"] = new TH1F("CosThetaDiffPt10ToInfGenTenthBin","CosThetaDiffPt10ToInfGenTenthBin",100,-2,2);
  histo1D["CosThetaRecoPt10ToInfGenFirstBin"] = new TH1F("CosThetaRecoPt10ToInfGenFirstBin","CosThetaRecoPt10ToInfGenFirstBin",100,-2,2);
  histo1D["CosThetaRecoPt10ToInfGenSecondBin"] = new TH1F("CosThetaRecoPt10ToInfGenSecondBin","CosThetaRecoPt10ToInfGenSecondBin",100,-2,2);
  histo1D["CosThetaRecoPt10ToInfGenThirdBin"] = new TH1F("CosThetaRecoPt10ToInfGenThirdBin","CosThetaRecoPt10ToInfGenThirdBin",100,-2,2);
  histo1D["CosThetaRecoPt10ToInfGenFourthBin"] = new TH1F("CosThetaRecoPt10ToInfGenFourthBin","CosThetaRecoPt10ToInfGenFourthBin",100,-2,2);
  histo1D["CosThetaRecoPt10ToInfGenFifthBin"] = new TH1F("CosThetaRecoPt10ToInfGenFifthBin","CosThetaRecoPt10ToInfGenFifthBin",100,-2,2);
  histo1D["CosThetaRecoPt10ToInfGenSixthBin"] = new TH1F("CosThetaRecoPt10ToInfGenSixthBin","CosThetaRecoPt10ToInfGenSixthBin",100,-2,2);
  histo1D["CosThetaRecoPt10ToInfGenSeventhBin"] = new TH1F("CosThetaRecoPt10ToInfGenSeventhBin","CosThetaRecoPt10ToInfGenSeventhBin",100,-2,2);
  histo1D["CosThetaRecoPt10ToInfGenEighthBin"] = new TH1F("CosThetaRecoPt10ToInfGenEighthBin","CosThetaRecoPt10ToInfGenEighthBin",100,-2,2);
  histo1D["CosThetaRecoPt10ToInfGenNinthBin"] = new TH1F("CosThetaRecoPt10ToInfGenNinthBin","CosThetaRecoPt10ToInfGenNinthBin",100,-2,2);
  histo1D["CosThetaRecoPt10ToInfGenTenthBin"] = new TH1F("CosThetaRecoPt10ToInfGenTenthBin","CosThetaRecoPt10ToInfGenTenthBin",100,-2,2);
  histo1D["CosThetaRecoPt10ToInfGenFirstBin"] = new TH1F("CosThetaRecoPt10ToInfGenFirstBin","CosThetaRecoPt10ToInfGenFirstBin",50,-1,1);

  histo1D["CosThetaDifferencePt15ToInf"] = new TH1F("CosThetaDifferencePt15ToInf","CosThetaDifferencePt15ToInf",200,-2,2);
  histo1D["CosThetaRecoPt15ToInf"] = new TH1F("CosThetaRecoPt15ToInf","CosThetaRecoPt15ToInf",100,-1,1);
  histo1D["CosThetaRecoPt15ToInfForRatio"] = new TH1F("CosThetaRecoPt15ToInfForRatio","CosThetaRecoPt15ToInfForRatio",100,-1,1);
  histo1D["CosThetaGenPt15ToInf"] = new TH1F("CosThetaGenPt15ToInf","CosThetaGenPt15ToInf",100,-1,1);
  TProfile *CosThetaDiffVSRecoPt15ToInf = new TProfile("CosThetaDiffVSRecoPt15ToInf","TProfile of Cos Theta* (reco-gen) vs Cos Theta* (reco) for muon Pt 15 to inf GeV",10,-1,1,-2,2);
  histo2D["CosThetaGenVsRecoPt15ToInf"] = new TH2F("CosThetaGenVsRecoPt15ToInf","CosThetaGenVsRecoPt15ToInf",50,-1,1,20,-1,1);
  histo1D["CosThetaDiffPt15ToInfGenFirstBin"] = new TH1F("CosThetaDiffPt15ToInfGenFirstBin","CosThetaDiffPt15ToInfGenFirstBin",100,-2,2);
  histo1D["CosThetaDiffPt15ToInfGenSecondBin"] = new TH1F("CosThetaDiffPt15ToInfGenSecondBin","CosThetaDiffPt15ToInfGenSecondBin",100,-2,2);
  histo1D["CosThetaDiffPt15ToInfGenThirdBin"] = new TH1F("CosThetaDiffPt15ToInfGenThirdBin","CosThetaDiffPt15ToInfGenThirdBin",100,-2,2);
  histo1D["CosThetaDiffPt15ToInfGenFourthBin"] = new TH1F("CosThetaDiffPt15ToInfGenFourthBin","CosThetaDiffPt15ToInfGenFourthBin",100,-2,2);
  histo1D["CosThetaDiffPt15ToInfGenFifthBin"] = new TH1F("CosThetaDiffPt15ToInfGenFifthBin","CosThetaDiffPt15ToInfGenFifthBin",100,-2,2);
  histo1D["CosThetaDiffPt15ToInfGenSixthBin"] = new TH1F("CosThetaDiffPt15ToInfGenSixthBin","CosThetaDiffPt15ToInfGenSixthBin",100,-2,2);
  histo1D["CosThetaDiffPt15ToInfGenSeventhBin"] = new TH1F("CosThetaDiffPt15ToInfGenSeventhBin","CosThetaDiffPt15ToInfGenSeventhBin",100,-2,2);
  histo1D["CosThetaDiffPt15ToInfGenEighthBin"] = new TH1F("CosThetaDiffPt15ToInfGenEighthBin","CosThetaDiffPt15ToInfGenEighthBin",100,-2,2);
  histo1D["CosThetaDiffPt15ToInfGenNinthBin"] = new TH1F("CosThetaDiffPt15ToInfGenNinthBin","CosThetaDiffPt15ToInfGenNinthBin",100,-2,2);
  histo1D["CosThetaDiffPt15ToInfGenTenthBin"] = new TH1F("CosThetaDiffPt15ToInfGenTenthBin","CosThetaDiffPt15ToInfGenTenthBin",100,-2,2);
  histo1D["CosThetaRecoPt15ToInfGenFirstBin"] = new TH1F("CosThetaRecoPt15ToInfGenFirstBin","CosThetaRecoPt15ToInfGenFirstBin",100,-2,2);
  histo1D["CosThetaRecoPt15ToInfGenSecondBin"] = new TH1F("CosThetaRecoPt15ToInfGenSecondBin","CosThetaRecoPt15ToInfGenSecondBin",100,-2,2);
  histo1D["CosThetaRecoPt15ToInfGenThirdBin"] = new TH1F("CosThetaRecoPt15ToInfGenThirdBin","CosThetaRecoPt15ToInfGenThirdBin",100,-2,2);
  histo1D["CosThetaRecoPt15ToInfGenFourthBin"] = new TH1F("CosThetaRecoPt15ToInfGenFourthBin","CosThetaRecoPt15ToInfGenFourthBin",100,-2,2);
  histo1D["CosThetaRecoPt15ToInfGenFifthBin"] = new TH1F("CosThetaRecoPt15ToInfGenFifthBin","CosThetaRecoPt15ToInfGenFifthBin",100,-2,2);
  histo1D["CosThetaRecoPt15ToInfGenSixthBin"] = new TH1F("CosThetaRecoPt15ToInfGenSixthBin","CosThetaRecoPt15ToInfGenSixthBin",100,-2,2);
  histo1D["CosThetaRecoPt15ToInfGenSeventhBin"] = new TH1F("CosThetaRecoPt15ToInfGenSeventhBin","CosThetaRecoPt15ToInfGenSeventhBin",100,-2,2);
  histo1D["CosThetaRecoPt15ToInfGenEighthBin"] = new TH1F("CosThetaRecoPt15ToInfGenEighthBin","CosThetaRecoPt15ToInfGenEighthBin",100,-2,2);
  histo1D["CosThetaRecoPt15ToInfGenNinthBin"] = new TH1F("CosThetaRecoPt15ToInfGenNinthBin","CosThetaRecoPt15ToInfGenNinthBin",100,-2,2);
  histo1D["CosThetaRecoPt15ToInfGenTenthBin"] = new TH1F("CosThetaRecoPt15ToInfGenTenthBin","CosThetaRecoPt15ToInfGenTenthBin",100,-2,2);
  histo1D["CosThetaRecoPt15ToInfGenFirstBin"] = new TH1F("CosThetaRecoPt15ToInfGenFirstBin","CosThetaRecoPt15ToInfGenFirstBin",50,-1,1);

  histo1D["CosThetaDifferencePt20ToInf"] = new TH1F("CosThetaDifferencePt20ToInf","CosThetaDifferencePt20ToInf",200,-2,2);
  histo1D["CosThetaRecoPt20ToInf"] = new TH1F("CosThetaRecoPt20ToInf","CosThetaRecoPt20ToInf",100,-1,1);
  histo1D["CosThetaRecoPt20ToInfForRatio"] = new TH1F("CosThetaRecoPt20ToInfForRatio","CosThetaRecoPt20ToInfForRatio",100,-1,1);
  histo1D["CosThetaGenPt20ToInf"] = new TH1F("CosThetaGenPt20ToInf","CosThetaGenPt20ToInf",100,-1,1);
  TProfile *CosThetaDiffVSRecoPt20ToInf = new TProfile("CosThetaDiffVSRecoPt20ToInf","TProfile of Cos Theta* (reco-gen) vs Cos Theta* (reco) for muon Pt 20 to inf GeV",10,-1,1,-2,2);
  histo2D["CosThetaGenVsRecoPt20ToInf"] = new TH2F("CosThetaGenVsRecoPt20ToInf","CosThetaGenVsRecoPt20ToInf",50,-1,1,20,-1,1);
  histo1D["CosThetaDiffPt20ToInfGenFirstBin"] = new TH1F("CosThetaDiffPt20ToInfGenFirstBin","CosThetaDiffPt20ToInfGenFirstBin",100,-2,2);
  histo1D["CosThetaDiffPt20ToInfGenSecondBin"] = new TH1F("CosThetaDiffPt20ToInfGenSecondBin","CosThetaDiffPt20ToInfGenSecondBin",100,-2,2);
  histo1D["CosThetaDiffPt20ToInfGenThirdBin"] = new TH1F("CosThetaDiffPt20ToInfGenThirdBin","CosThetaDiffPt20ToInfGenThirdBin",100,-2,2);
  histo1D["CosThetaDiffPt20ToInfGenFourthBin"] = new TH1F("CosThetaDiffPt20ToInfGenFourthBin","CosThetaDiffPt20ToInfGenFourthBin",100,-2,2);
  histo1D["CosThetaDiffPt20ToInfGenFifthBin"] = new TH1F("CosThetaDiffPt20ToInfGenFifthBin","CosThetaDiffPt20ToInfGenFifthBin",100,-2,2);
  histo1D["CosThetaDiffPt20ToInfGenSixthBin"] = new TH1F("CosThetaDiffPt20ToInfGenSixthBin","CosThetaDiffPt20ToInfGenSixthBin",100,-2,2);
  histo1D["CosThetaDiffPt20ToInfGenSeventhBin"] = new TH1F("CosThetaDiffPt20ToInfGenSeventhBin","CosThetaDiffPt20ToInfGenSeventhBin",100,-2,2);
  histo1D["CosThetaDiffPt20ToInfGenEighthBin"] = new TH1F("CosThetaDiffPt20ToInfGenEighthBin","CosThetaDiffPt20ToInfGenEighthBin",100,-2,2);
  histo1D["CosThetaDiffPt20ToInfGenNinthBin"] = new TH1F("CosThetaDiffPt20ToInfGenNinthBin","CosThetaDiffPt20ToInfGenNinthBin",100,-2,2);
  histo1D["CosThetaDiffPt20ToInfGenTenthBin"] = new TH1F("CosThetaDiffPt20ToInfGenTenthBin","CosThetaDiffPt20ToInfGenTenthBin",100,-2,2);
  histo1D["CosThetaRecoPt20ToInfGenFirstBin"] = new TH1F("CosThetaRecoPt20ToInfGenFirstBin","CosThetaRecoPt20ToInfGenFirstBin",100,-2,2);
  histo1D["CosThetaRecoPt20ToInfGenSecondBin"] = new TH1F("CosThetaRecoPt20ToInfGenSecondBin","CosThetaRecoPt20ToInfGenSecondBin",100,-2,2);
  histo1D["CosThetaRecoPt20ToInfGenThirdBin"] = new TH1F("CosThetaRecoPt20ToInfGenThirdBin","CosThetaRecoPt20ToInfGenThirdBin",100,-2,2);
  histo1D["CosThetaRecoPt20ToInfGenFourthBin"] = new TH1F("CosThetaRecoPt20ToInfGenFourthBin","CosThetaRecoPt20ToInfGenFourthBin",100,-2,2);
  histo1D["CosThetaRecoPt20ToInfGenFifthBin"] = new TH1F("CosThetaRecoPt20ToInfGenFifthBin","CosThetaRecoPt20ToInfGenFifthBin",100,-2,2);
  histo1D["CosThetaRecoPt20ToInfGenSixthBin"] = new TH1F("CosThetaRecoPt20ToInfGenSixthBin","CosThetaRecoPt20ToInfGenSixthBin",100,-2,2);
  histo1D["CosThetaRecoPt20ToInfGenSeventhBin"] = new TH1F("CosThetaRecoPt20ToInfGenSeventhBin","CosThetaRecoPt20ToInfGenSeventhBin",100,-2,2);
  histo1D["CosThetaRecoPt20ToInfGenEighthBin"] = new TH1F("CosThetaRecoPt20ToInfGenEighthBin","CosThetaRecoPt20ToInfGenEighthBin",100,-2,2);
  histo1D["CosThetaRecoPt20ToInfGenNinthBin"] = new TH1F("CosThetaRecoPt20ToInfGenNinthBin","CosThetaRecoPt20ToInfGenNinthBin",100,-2,2);
  histo1D["CosThetaRecoPt20ToInfGenTenthBin"] = new TH1F("CosThetaRecoPt20ToInfGenTenthBin","CosThetaRecoPt20ToInfGenTenthBin",100,-2,2);
  histo1D["CosThetaRecoPt20ToInfGenFirstBin"] = new TH1F("CosThetaRecoPt20ToInfGenFirstBin","CosThetaRecoPt20ToInfGenFirstBin",50,-1,1);

  histo1D["CosThetaDifferencePt25ToInf"] = new TH1F("CosThetaDifferencePt25ToInf","CosThetaDifferencePt25ToInf",200,-2,2);
  histo1D["CosThetaRecoPt25ToInf"] = new TH1F("CosThetaRecoPt25ToInf","CosThetaRecoPt25ToInf",100,-1,1);
  histo1D["CosThetaRecoPt25ToInfForRatio"] = new TH1F("CosThetaRecoPt25ToInfForRatio","CosThetaRecoPt25ToInfForRatio",100,-1,1);
  histo1D["CosThetaGenPt25ToInf"] = new TH1F("CosThetaGenPt25ToInf","CosThetaGenPt25ToInf",100,-1,1);
  TProfile *CosThetaDiffVSRecoPt25ToInf = new TProfile("CosThetaDiffVSRecoPt25ToInf","TProfile of Cos Theta* (reco-gen) vs Cos Theta* (reco) for muon Pt 25 to inf GeV",10,-1,1,-2,2);
  histo2D["CosThetaGenVsRecoPt25ToInf"] = new TH2F("CosThetaGenVsRecoPt25ToInf","CosThetaGenVsRecoPt25ToInf",50,-1,1,20,-1,1);
  histo1D["CosThetaDiffPt25ToInfGenFirstBin"] = new TH1F("CosThetaDiffPt25ToInfGenFirstBin","CosThetaDiffPt25ToInfGenFirstBin",100,-2,2);
  histo1D["CosThetaDiffPt25ToInfGenSecondBin"] = new TH1F("CosThetaDiffPt25ToInfGenSecondBin","CosThetaDiffPt25ToInfGenSecondBin",100,-2,2);
  histo1D["CosThetaDiffPt25ToInfGenThirdBin"] = new TH1F("CosThetaDiffPt25ToInfGenThirdBin","CosThetaDiffPt25ToInfGenThirdBin",100,-2,2);
  histo1D["CosThetaDiffPt25ToInfGenFourthBin"] = new TH1F("CosThetaDiffPt25ToInfGenFourthBin","CosThetaDiffPt25ToInfGenFourthBin",100,-2,2);
  histo1D["CosThetaDiffPt25ToInfGenFifthBin"] = new TH1F("CosThetaDiffPt25ToInfGenFifthBin","CosThetaDiffPt25ToInfGenFifthBin",100,-2,2);
  histo1D["CosThetaDiffPt25ToInfGenSixthBin"] = new TH1F("CosThetaDiffPt25ToInfGenSixthBin","CosThetaDiffPt25ToInfGenSixthBin",100,-2,2);
  histo1D["CosThetaDiffPt25ToInfGenSeventhBin"] = new TH1F("CosThetaDiffPt25ToInfGenSeventhBin","CosThetaDiffPt25ToInfGenSeventhBin",100,-2,2);
  histo1D["CosThetaDiffPt25ToInfGenEighthBin"] = new TH1F("CosThetaDiffPt25ToInfGenEighthBin","CosThetaDiffPt25ToInfGenEighthBin",100,-2,2);
  histo1D["CosThetaDiffPt25ToInfGenNinthBin"] = new TH1F("CosThetaDiffPt25ToInfGenNinthBin","CosThetaDiffPt25ToInfGenNinthBin",100,-2,2);
  histo1D["CosThetaDiffPt25ToInfGenTenthBin"] = new TH1F("CosThetaDiffPt25ToInfGenTenthBin","CosThetaDiffPt25ToInfGenTenthBin",100,-2,2);
  histo1D["CosThetaRecoPt25ToInfGenFirstBin"] = new TH1F("CosThetaRecoPt25ToInfGenFirstBin","CosThetaRecoPt25ToInfGenFirstBin",100,-2,2);
  histo1D["CosThetaRecoPt25ToInfGenSecondBin"] = new TH1F("CosThetaRecoPt25ToInfGenSecondBin","CosThetaRecoPt25ToInfGenSecondBin",100,-2,2);
  histo1D["CosThetaRecoPt25ToInfGenThirdBin"] = new TH1F("CosThetaRecoPt25ToInfGenThirdBin","CosThetaRecoPt25ToInfGenThirdBin",100,-2,2);
  histo1D["CosThetaRecoPt25ToInfGenFourthBin"] = new TH1F("CosThetaRecoPt25ToInfGenFourthBin","CosThetaRecoPt25ToInfGenFourthBin",100,-2,2);
  histo1D["CosThetaRecoPt25ToInfGenFifthBin"] = new TH1F("CosThetaRecoPt25ToInfGenFifthBin","CosThetaRecoPt25ToInfGenFifthBin",100,-2,2);
  histo1D["CosThetaRecoPt25ToInfGenSixthBin"] = new TH1F("CosThetaRecoPt25ToInfGenSixthBin","CosThetaRecoPt25ToInfGenSixthBin",100,-2,2);
  histo1D["CosThetaRecoPt25ToInfGenSeventhBin"] = new TH1F("CosThetaRecoPt25ToInfGenSeventhBin","CosThetaRecoPt25ToInfGenSeventhBin",100,-2,2);
  histo1D["CosThetaRecoPt25ToInfGenEighthBin"] = new TH1F("CosThetaRecoPt25ToInfGenEighthBin","CosThetaRecoPt25ToInfGenEighthBin",100,-2,2);
  histo1D["CosThetaRecoPt25ToInfGenNinthBin"] = new TH1F("CosThetaRecoPt25ToInfGenNinthBin","CosThetaRecoPt25ToInfGenNinthBin",100,-2,2);
  histo1D["CosThetaRecoPt25ToInfGenTenthBin"] = new TH1F("CosThetaRecoPt25ToInfGenTenthBin","CosThetaRecoPt25ToInfGenTenthBin",100,-2,2);
  histo1D["CosThetaRecoPt25ToInfGenFirstBin"] = new TH1F("CosThetaRecoPt25ToInfGenFirstBin","CosThetaRecoPt25ToInfGenFirstBin",50,-1,1);

  histo1D["CosThetaDifferencePt30ToInf"] = new TH1F("CosThetaDifferencePt30ToInf","CosThetaDifferencePt30ToInf",200,-2,2);
  histo1D["CosThetaRecoPt30ToInf"] = new TH1F("CosThetaRecoPt30ToInf","CosThetaRecoPt30ToInf",100,-1,1);
  histo1D["CosThetaRecoPt30ToInfForRatio"] = new TH1F("CosThetaRecoPt30ToInfForRatio","CosThetaRecoPt30ToInfForRatio",100,-1,1);
  histo1D["CosThetaGenPt30ToInf"] = new TH1F("CosThetaGenPt30ToInf","CosThetaGenPt30ToInf",100,-1,1);
  TProfile *CosThetaDiffVSRecoPt30ToInf = new TProfile("CosThetaDiffVSRecoPt30ToInf","TProfile of Cos Theta* (reco-gen) vs Cos Theta* (reco) for muon Pt 30 to inf GeV",10,-1,1,-2,2);
  histo2D["CosThetaGenVsRecoPt30ToInf"] = new TH2F("CosThetaGenVsRecoPt30ToInf","CosThetaGenVsRecoPt30ToInf",50,-1,1,20,-1,1);
  histo1D["CosThetaDiffPt30ToInfGenFirstBin"] = new TH1F("CosThetaDiffPt30ToInfGenFirstBin","CosThetaDiffPt30ToInfGenFirstBin",100,-2,2);
  histo1D["CosThetaDiffPt30ToInfGenSecondBin"] = new TH1F("CosThetaDiffPt30ToInfGenSecondBin","CosThetaDiffPt30ToInfGenSecondBin",100,-2,2);
  histo1D["CosThetaDiffPt30ToInfGenThirdBin"] = new TH1F("CosThetaDiffPt30ToInfGenThirdBin","CosThetaDiffPt30ToInfGenThirdBin",100,-2,2);
  histo1D["CosThetaDiffPt30ToInfGenFourthBin"] = new TH1F("CosThetaDiffPt30ToInfGenFourthBin","CosThetaDiffPt30ToInfGenFourthBin",100,-2,2);
  histo1D["CosThetaDiffPt30ToInfGenFifthBin"] = new TH1F("CosThetaDiffPt30ToInfGenFifthBin","CosThetaDiffPt30ToInfGenFifthBin",100,-2,2);
  histo1D["CosThetaDiffPt30ToInfGenSixthBin"] = new TH1F("CosThetaDiffPt30ToInfGenSixthBin","CosThetaDiffPt30ToInfGenSixthBin",100,-2,2);
  histo1D["CosThetaDiffPt30ToInfGenSeventhBin"] = new TH1F("CosThetaDiffPt30ToInfGenSeventhBin","CosThetaDiffPt30ToInfGenSeventhBin",100,-2,2);
  histo1D["CosThetaDiffPt30ToInfGenEighthBin"] = new TH1F("CosThetaDiffPt30ToInfGenEighthBin","CosThetaDiffPt30ToInfGenEighthBin",100,-2,2);
  histo1D["CosThetaDiffPt30ToInfGenNinthBin"] = new TH1F("CosThetaDiffPt30ToInfGenNinthBin","CosThetaDiffPt30ToInfGenNinthBin",100,-2,2);
  histo1D["CosThetaDiffPt30ToInfGenTenthBin"] = new TH1F("CosThetaDiffPt30ToInfGenTenthBin","CosThetaDiffPt30ToInfGenTenthBin",100,-2,2);
  histo1D["CosThetaRecoPt30ToInfGenFirstBin"] = new TH1F("CosThetaRecoPt30ToInfGenFirstBin","CosThetaRecoPt30ToInfGenFirstBin",100,-2,2);
  histo1D["CosThetaRecoPt30ToInfGenSecondBin"] = new TH1F("CosThetaRecoPt30ToInfGenSecondBin","CosThetaRecoPt30ToInfGenSecondBin",100,-2,2);
  histo1D["CosThetaRecoPt30ToInfGenThirdBin"] = new TH1F("CosThetaRecoPt30ToInfGenThirdBin","CosThetaRecoPt30ToInfGenThirdBin",100,-2,2);
  histo1D["CosThetaRecoPt30ToInfGenFourthBin"] = new TH1F("CosThetaRecoPt30ToInfGenFourthBin","CosThetaRecoPt30ToInfGenFourthBin",100,-2,2);
  histo1D["CosThetaRecoPt30ToInfGenFifthBin"] = new TH1F("CosThetaRecoPt30ToInfGenFifthBin","CosThetaRecoPt30ToInfGenFifthBin",100,-2,2);
  histo1D["CosThetaRecoPt30ToInfGenSixthBin"] = new TH1F("CosThetaRecoPt30ToInfGenSixthBin","CosThetaRecoPt30ToInfGenSixthBin",100,-2,2);
  histo1D["CosThetaRecoPt30ToInfGenSeventhBin"] = new TH1F("CosThetaRecoPt30ToInfGenSeventhBin","CosThetaRecoPt30ToInfGenSeventhBin",100,-2,2);
  histo1D["CosThetaRecoPt30ToInfGenEighthBin"] = new TH1F("CosThetaRecoPt30ToInfGenEighthBin","CosThetaRecoPt30ToInfGenEighthBin",100,-2,2);
  histo1D["CosThetaRecoPt30ToInfGenNinthBin"] = new TH1F("CosThetaRecoPt30ToInfGenNinthBin","CosThetaRecoPt30ToInfGenNinthBin",100,-2,2);
  histo1D["CosThetaRecoPt30ToInfGenTenthBin"] = new TH1F("CosThetaRecoPt30ToInfGenTenthBin","CosThetaRecoPt30ToInfGenTenthBin",100,-2,2);
  histo1D["CosThetaRecoPt30ToInfGenFirstBin"] = new TH1F("CosThetaRecoPt30ToInfGenFirstBin","CosThetaRecoPt30ToInfGenFirstBin",50,-1,1);
  
  histo1D["CosThetaDifferencePt35ToInf"] = new TH1F("CosThetaDifferencePt35ToInf","CosThetaDifferencePt35ToInf",200,-2,2);
  histo1D["CosThetaRecoPt35ToInf"] = new TH1F("CosThetaRecoPt35ToInf","CosThetaRecoPt35ToInf",100,-1,1);
  histo1D["CosThetaRecoPt35ToInfForRatio"] = new TH1F("CosThetaRecoPt35ToInfForRatio","CosThetaRecoPt35ToInfForRatio",100,-1,1);
  histo1D["CosThetaGenPt35ToInf"] = new TH1F("CosThetaGenPt35ToInf","CosThetaGenPt35ToInf",100,-1,1);
  TProfile *CosThetaDiffVSRecoPt35ToInf = new TProfile("CosThetaDiffVSRecoPt35ToInf","TProfile of Cos Theta* (reco-gen) vs Cos Theta* (reco) for muon Pt 35 to inf GeV",10,-1,1,-2,2);
  histo2D["CosThetaGenVsRecoPt35ToInf"] = new TH2F("CosThetaGenVsRecoPt35ToInf","CosThetaGenVsRecoPt35ToInf",50,-1,1,20,-1,1);
  histo1D["CosThetaDiffPt35ToInfGenFirstBin"] = new TH1F("CosThetaDiffPt35ToInfGenFirstBin","CosThetaDiffPt35ToInfGenFirstBin",100,-2,2);
  histo1D["CosThetaDiffPt35ToInfGenSecondBin"] = new TH1F("CosThetaDiffPt35ToInfGenSecondBin","CosThetaDiffPt35ToInfGenSecondBin",100,-2,2);
  histo1D["CosThetaDiffPt35ToInfGenThirdBin"] = new TH1F("CosThetaDiffPt35ToInfGenThirdBin","CosThetaDiffPt35ToInfGenThirdBin",100,-2,2);
  histo1D["CosThetaDiffPt35ToInfGenFourthBin"] = new TH1F("CosThetaDiffPt35ToInfGenFourthBin","CosThetaDiffPt35ToInfGenFourthBin",100,-2,2);
  histo1D["CosThetaDiffPt35ToInfGenFifthBin"] = new TH1F("CosThetaDiffPt35ToInfGenFifthBin","CosThetaDiffPt35ToInfGenFifthBin",100,-2,2);
  histo1D["CosThetaDiffPt35ToInfGenSixthBin"] = new TH1F("CosThetaDiffPt35ToInfGenSixthBin","CosThetaDiffPt35ToInfGenSixthBin",100,-2,2);
  histo1D["CosThetaDiffPt35ToInfGenSeventhBin"] = new TH1F("CosThetaDiffPt35ToInfGenSeventhBin","CosThetaDiffPt35ToInfGenSeventhBin",100,-2,2);
  histo1D["CosThetaDiffPt35ToInfGenEighthBin"] = new TH1F("CosThetaDiffPt35ToInfGenEighthBin","CosThetaDiffPt35ToInfGenEighthBin",100,-2,2);
  histo1D["CosThetaDiffPt35ToInfGenNinthBin"] = new TH1F("CosThetaDiffPt35ToInfGenNinthBin","CosThetaDiffPt35ToInfGenNinthBin",100,-2,2);
  histo1D["CosThetaDiffPt35ToInfGenTenthBin"] = new TH1F("CosThetaDiffPt35ToInfGenTenthBin","CosThetaDiffPt35ToInfGenTenthBin",100,-2,2);
  histo1D["CosThetaRecoPt35ToInfGenFirstBin"] = new TH1F("CosThetaRecoPt35ToInfGenFirstBin","CosThetaRecoPt35ToInfGenFirstBin",100,-2,2);
  histo1D["CosThetaRecoPt35ToInfGenSecondBin"] = new TH1F("CosThetaRecoPt35ToInfGenSecondBin","CosThetaRecoPt35ToInfGenSecondBin",100,-2,2);
  histo1D["CosThetaRecoPt35ToInfGenThirdBin"] = new TH1F("CosThetaRecoPt35ToInfGenThirdBin","CosThetaRecoPt35ToInfGenThirdBin",100,-2,2);
  histo1D["CosThetaRecoPt35ToInfGenFourthBin"] = new TH1F("CosThetaRecoPt35ToInfGenFourthBin","CosThetaRecoPt35ToInfGenFourthBin",100,-2,2);
  histo1D["CosThetaRecoPt35ToInfGenFifthBin"] = new TH1F("CosThetaRecoPt35ToInfGenFifthBin","CosThetaRecoPt35ToInfGenFifthBin",100,-2,2);
  histo1D["CosThetaRecoPt35ToInfGenSixthBin"] = new TH1F("CosThetaRecoPt35ToInfGenSixthBin","CosThetaRecoPt35ToInfGenSixthBin",100,-2,2);
  histo1D["CosThetaRecoPt35ToInfGenSeventhBin"] = new TH1F("CosThetaRecoPt35ToInfGenSeventhBin","CosThetaRecoPt35ToInfGenSeventhBin",100,-2,2);
  histo1D["CosThetaRecoPt35ToInfGenEighthBin"] = new TH1F("CosThetaRecoPt35ToInfGenEighthBin","CosThetaRecoPt35ToInfGenEighthBin",100,-2,2);
  histo1D["CosThetaRecoPt35ToInfGenNinthBin"] = new TH1F("CosThetaRecoPt35ToInfGenNinthBin","CosThetaRecoPt35ToInfGenNinthBin",100,-2,2);
  histo1D["CosThetaRecoPt35ToInfGenTenthBin"] = new TH1F("CosThetaRecoPt35ToInfGenTenthBin","CosThetaRecoPt35ToInfGenTenthBin",100,-2,2);
  histo1D["CosThetaRecoPt35ToInfGenFirstBin"] = new TH1F("CosThetaRecoPt35ToInfGenFirstBin","CosThetaRecoPt35ToInfGenFirstBin",50,-1,1);

  histo1D["CosThetaResLastBin"] = new TH1F("CosThetaResLastBin","CosThetaResLastBin",200,-2,2);
  histo1D["CosThetaResSecondNegBin"] = new TH1F("CosThetaResSecondNegBin","CosThetaResSecondNegBin",200,-2,2);
  histo1D["CosThetaResThirdNegBin"] = new TH1F("CosThetaResThirdNegBin","CosThetaResThirdNegBin",200,-2,2);
  histo1D["CosThetaResFourthNegBin"] = new TH1F("CosThetaResFourthNegBin","CosThetaResFourthNegBin",200,-2,2);
  histo1D["CosThetaResFifthNegBin"] = new TH1F("CosThetaResFifthNegBin","CosThetaResFifthNegBin",200,-2,2);
  //histo1D["CosThetaResSixtNegBin"] = new TH1F("CosThetaResSixtNegBin","CosThetaResSixtNegBin",200,-2,2);
  //histo1D["CosThetaResSeventhNegBin"] = new TH1F("CosThetaResSeventhNegBin","CosThetaResSeventhNegBin",200,-2,2);
  //histo1D["CosThetaResEighthNegBin"] = new TH1F("CosThetaResEighthNegBin","CosThetaResEighthNegBin",200,-2,2);
  //histo1D["CosThetaResNinthNegBin"] = new TH1F("CosThetaResNinthNegBin","CosThetaResNinthNegBin",200,-2,2);
  //histo1D["CosThetaResTenthNegBin"] = new TH1F("CosThetaResTenthNegBin","CosThetaResTenthNegBin",200,-2,2);
  //histo1D["CosThetaResTenthBin"] = new TH1F("CosThetaResTenthBin","CosThetaResTenthBin",200,-2,2);
  //histo1D["CosThetaResNinthBin"] = new TH1F("CosThetaResNinthBin","CosThetaResNinthBin",200,-2,2);
  //histo1D["CosThetaResEighthBin"] = new TH1F("CosThetaResEighthBin","CosThetaResEighthBin",200,-2,2);
  //histo1D["CosThetaResSeventhBin"] = new TH1F("CosThetaResSeventhBin","CosThetaResSeventhBin",200,-2,2);
  //histo1D["CosThetaResSixtBin"] = new TH1F("CosThetaResSixtBin","CosThetaResSixtBin",200,-2,2);
  histo1D["CosThetaResFifthBin"] = new TH1F("CosThetaResFifthBin","CosThetaResFifthBin",200,-2,2);
  histo1D["CosThetaResFourthBin"] = new TH1F("CosThetaResFourthBin","CosThetaResFourthBin",200,-2,2);
  histo1D["CosThetaResThirdBin"] = new TH1F("CosThetaResThirdBin","CosThetaResThirdBin",200,-2,2);
  histo1D["CosThetaResSecondBin"] = new TH1F("CosThetaResSecondBin","CosThetaResSecondBin",200,-2,2);
  histo1D["CosThetaResFirstBin"] = new TH1F("CosThetaResFirstBin","CosThetaResFirstBin",200,-2,2);

  histo1D["CosThetaResTCHELLept"] = new TH1F("CosThetaResTCHELLept","CosThetaResTCHELLept",200,-2,2);
  histo1D["CosThetaResTCHEMLept"] = new TH1F("CosThetaResTCHEMLept","CosThetaResTCHEMLept",200,-2,2);
  histo1D["CosThetaResTCHPMLept"] = new TH1F("CosThetaResTCHPMLept","CosThetaResTCHPMLept",200,-2,2);
  histo1D["CosThetaResTCHPTLept"] = new TH1F("CosThetaResTCHPTLept","CosThetaResTCHPTLept",200,-2,2);
  histo1D["CosThetaResSSVHEMLept"] = new TH1F("CosThetaResSSVHEMLept","CosThetaResSSVHEMLept",200,-2,2);
  histo1D["CosThetaResSSVHPTLept"] = new TH1F("CosThetaResSSVHPTLept","CosThetaResSSVHPTLept",200,-2,2);
  histo1D["CosThetaResCSVLLept"] = new TH1F("CosThetaResCSVLLept","CosThetaResCSVLLept",200,-2,2);
  histo1D["CosThetaResCSVMLept"] = new TH1F("CosThetaResCSVMLept","CosThetaResCSVMLept",200,-2,2);
  histo1D["CosThetaResCSVTLept"] = new TH1F("CosThetaResCSVMTLept","CosThetaResCSVMTLept",200,-2,2);

  //Histograms to obtain ratio of cos theta* for alternative and SM helicities
  float HelicityWeight[3];
  int SizeArray = 3;
  float HelicityFraction[SizeArray][3];  //0:Longitudinal; 1:Righthanded; 2:Lefthanded
  float UsedDistributionValue[SizeArray];
  std::ostringstream HelicityNumbers[SizeArray][3];
  cout << " size of array : " << SizeArray << endl;
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

  //MSPlots for mass distributions after Kinematic Fit:
  MSPlot["HadronicWMass"] = new MultiSamplePlot(datasets,"HadronicWMass",50,50,100,"HadronicWMass");
  MSPlot["HadronicTopMass"] = new MultiSamplePlot(datasets,"HadronicTopMass",50,150,200,"HadronicTopMass");
  MSPlot["LeptonicWMass"] = new MultiSamplePlot(datasets,"LeptonicWMass",50,50,100,"LeptonicWMass");
  MSPlot["LeptonicTopMass"] = new MultiSamplePlot(datasets,"LeptonicTopMass",50,150,200,"LeptonicTopMass");
  
  //MSPlots for Cos theta distribution and KinFit probability distribution
  MSPlot["CosThetaSSVHEMLept"]= new MultiSamplePlot(datasets, "CosThetaSSVHEMLept", CosThetaBinNumber,-1,1,"CosThetaSSVHEMLept");  
  MSPlot["KinFitProbSSVHEMLept"]= new MultiSamplePlot(datasets, "KinFitProbSSVHEMLept", CosThetaBinNumber,0,1,"KinFitProbSSVHEMLept");  
  MSPlot["CosThetaTCHEMLept"]= new MultiSamplePlot(datasets, "CosThetaTCHEMLept", CosThetaBinNumber,-1,1,"CosThetaTCHEMLept");  
  MSPlot["KinFitProbTCHEMLept"]= new MultiSamplePlot(datasets, "KinFitProbTCHEMLept", CosThetaBinNumber,0,1,"KinFitProbTCHEMLept");  
  MSPlot["CosThetaTCHPMLept"]= new MultiSamplePlot(datasets, "CosThetaTCHPMLept", CosThetaBinNumber,-1,1,"CosThetaTCHPMLept");  
  MSPlot["CosThetaSSVHPMLept"]= new MultiSamplePlot(datasets, "CosThetaSSVHPMLept", CosThetaBinNumber,-1,1,"CosThetaSSVHPMLept");  
  MSPlot["KinFitProbSSVHPMLept"]= new MultiSamplePlot(datasets, "KinFitProbSSVHPMLept", CosThetaBinNumber,0,1,"KinFitProbSSVHPMLept");  
  MSPlot["CosThetaCSVMLept"]= new MultiSamplePlot(datasets, "CosThetaCSVMLept", CosThetaBinNumber,-1,1,"CosThetaCSVMLept");
 
  //Check nPrimary vertices for different executed cuts !!
  MSPlot["nPrimaryVert"] = new MultiSamplePlot(datasets,"nPrimaryVert" , 20, 0, 20, "nPrimaryVert");
  MSPlot["nPVBeforeCuts"] = new MultiSamplePlot(datasets,"nPVBeforeCuts" , 20, 0, 20, "nPVBeforeCuts");
  MSPlot["nPVAfterSSVHEMbTag"] = new MultiSamplePlot(datasets,"nPVAfterSSVHEMbTag" , 20, 0, 20, "nPVAfterSSVHEMbTag");
  MSPlot["nPVAfterMuon27Cut"] = new MultiSamplePlot(datasets,"nPVAfterMuon27Cut" , 20, 0, 20, "nPVAfterMuon27Cut");
  MSPlot["nPVAfterTransverseMassCut"] = new MultiSamplePlot(datasets,"nPVAfterTransverseMassCut" , 20, 0, 20, "nPVAfterTransverseMassCut");
  MSPlot["nPVBeforeFoundJetComb"] = new MultiSamplePlot(datasets, "nPVBeforeFoundJetComb", 20, 0, 20,"nPVBeforeFoundJetComb");
  MSPlot["nPVAfterFoundJetComb"] = new MultiSamplePlot(datasets, "nPVAfterFoundJetComb", 20, 0, 20,"nPVAfterFoundJetComb");
  MSPlot["nPVAfterFoundJetCombbTag"] = new MultiSamplePlot(datasets, "nPVAfterFoundJetCombbTab", 20, 0, 20,"nPVAfterFoundJetCombbTag");
  MSPlot["nPVAfterFoundCosTheta"] = new MultiSamplePlot(datasets, "nPVAfterFoundCosTheta", 20, 0, 20,"nPVAfterFoundCosTheta");
  MSPlot["nPVAfterFoundCosThetabTag"] = new MultiSamplePlot(datasets, "nPVAfterFoundCosThetabTag", 20, 0, 20,"nPVAfterFoundCosThetabTag");

  MSPlot["TransverseMassBeforeCut"]=new MultiSamplePlot(datasets,"TransverseMassBeforeCut",50,0,200,"TransverseMassBeforeCut");
  MSPlot["TransverseMassAfterCut"]=new MultiSamplePlot(datasets,"TransverseMassAfterCut",50,0,200,"TransverseMassAfterCut"); 
  
  //Histograms to check differences between events with Negative and positive DSquared (And check influence of probability cut) 
  MSPlot["CosThetaNobTag"]=new MultiSamplePlot(datasets,"CosThetaNobTag" ,CosThetaBinNumber,-1,1, "CosThetaNobTag");
  MSPlot["CosThetaNobTagProbCut"]=new MultiSamplePlot(datasets,"CosThetaNobTagProbCut" ,CosThetaBinNumber*2,-2,2,"CosThetaNobTagProbCut" );
  MSPlot["KinFitProbNobTag"]=new MultiSamplePlot(datasets, "KinFitProbNobTag",25,0,1,"KinFitProbNobTag");

  MSPlot["Jet1Pt"] = new MultiSamplePlot(datasets, "Jet1Pt", 100,0,500,"Jet1Pt");
  MSPlot["Jet2Pt"] = new MultiSamplePlot(datasets, "Jet2Pt", 100,0,500,"Jet2Pt");
  MSPlot["Jet3Pt"] = new MultiSamplePlot(datasets, "Jet3Pt", 100,0,500,"Jet3Pt");
  MSPlot["Jet4Pt"] = new MultiSamplePlot(datasets, "Jet4Pt", 100,0,500,"Jet4Pt");
  
  MSPlot["MetPhi"]= new MultiSamplePlot(datasets, "MetPhi",100,-3.3,3.3,"MetPhi");
  // MSPlot["NeutrinoKinFitPhi"]= new MultiSamplePlot(datasets,"NeutrinoKinFitPhi" ,50,-3.3,3.3,"NeutrinoKinFitPhi");
  // //MSPlot["NeutrinoPhiDiff"]= new MultiSamplePlot(datasets, "NeutrinoPhiDiff",100,-7,7,"NeutrinoPhiDiff");

  MSPlot["MetPx"]= new MultiSamplePlot(datasets, "MetPx",200,-400,400,"MetPx");
  // MSPlot["NeutrinoKinFitPx"]= new MultiSamplePlot(datasets,"NeutrinoKinFitPx" ,50,-200,200,"NeutrinoKinFitPx");
  // //MSPlot["NeutrinoPxDiff"]= new MultiSamplePlot(datasets, "NeutrinoPxDiff",100,-400,400,"NeutrinoPxDiff");
  
  MSPlot["MetPy"]= new MultiSamplePlot(datasets, "MetPy",200,-400,400,"MetPy");
  // MSPlot["NeutrinoKinFitPy"]= new MultiSamplePlot(datasets,"NeutrinoKinFitPy" ,50,-200,200,"NeutrinoKinFitPy");
  // //MSPlot["NeutrinoPyDiff"]= new MultiSamplePlot(datasets, "NeutrinoPyDiff",100,-400,400,"NeutrinoPyDiff");  

  if(UseChangedKinematics == true){
    MSPlot["BLeptPtAfterKinFit"]= new MultiSamplePlot(datasets,"BLeptPtAfterKinFit",50,0,250,"BLeptPtAfterKinFit");
    MSPlot["MetPtAfterKinFit"] = new MultiSamplePlot(datasets,"MetPtAfterKinFit",50,0,200,"MetPtAfterKinFit");
    MSPlot["MuonPtAfterKinFit"] = new MultiSamplePlot(datasets,"MuonPtAfterKinFit",40,0,150,"MuonPtAfterKinFit");
    MSPlot["WLeptPtAfterKinFit"] = new MultiSamplePlot(datasets,"WLeptPtAfterKinFit",50,0,250,"WLeptPtAfterKinFit");
    MSPlot["TopLeptPtAfterKinFit"] = new MultiSamplePlot(datasets,"TopLeptPtAfterKinFit",70,0,400,"TopLeptPtAfterKinFit");
  }
  else if(UseChangedKinematics == false){
    MSPlot["BLeptPtBeforeKinFit"]= new MultiSamplePlot(datasets,"BLeptPtBeforeKinFit",70,0,400,"BLeptPtBeforeKinFit");
    MSPlot["MetPtBeforeKinFit"] = new MultiSamplePlot(datasets,"MetPtBeforeKinFit",50,0,200,"MetPtBeforeKinFit");
    MSPlot["MuonPtBeforeKinFit"] = new MultiSamplePlot(datasets,"MuonPtBeforeKinFit",40,0,150,"MuonPtBeforeKinFit");
    MSPlot["WLeptPtBeforeKinFit"] = new MultiSamplePlot(datasets,"WLeptPtBeforeKinFit",50,0,250,"WLeptPtBeforeKinFit");
    MSPlot["TopLeptPtBeforeKinFit"] = new MultiSamplePlot(datasets,"TopLeptPtBeforeKinFit",70,0,400,"TopLeptPtBeforeKinFit");
  }

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
  if(fullDataSet == true) Lumi3DWeights = Lumi3DReWeighting("PileUpReweighting/pileup_MC_Flat10PlusTail.root", "PileUpReweighting/pileup_FineBin_2011Data_UpToRun180252.root", "pileup", "pileup");
  else if(fullDataSet == false) Lumi3DWeights = Lumi3DReWeighting("PileUpReweighting/pileup_MC_Flat10PlusTail.root", "PileUpReweighting/pileup_FineBin_2011Data_UpToRun173692.root", "pileup", "pileup");

  Lumi3DWeights.weight3D_init(1.0);

  PoissonMeanShifter PShiftUp = PoissonMeanShifter(0.6); // PU-systematic
  PoissonMeanShifter PShiftDown = PoissonMeanShifter(-0.6); // PU-systematic

  /////////////////////////////////////////
  // Initializing used variables
  /////////////////////////////////////////
  float MassW=83.6103;
  float MassTop = 172.956;
  float SigmaW=11.1534;  //Obtained from gaussian fit on Top and W distribution with simulated information
  float SigmaTop=18.232;

  std::string bTagValues[14];
  for(int ii=1;ii<=14;ii++){
    std::stringstream out;
    out << ii;
    bTagValues[ii-1]  =  out.str();
  }

  int NumberTCHEbTags = 13;
  int NumberTCHPbTags = 13;
  int NumberSSVHEbTags = 13;
  int NumberSSVHPbTags=13;
  int NumberCSVbTags=13;
  int TotalNumberbTags = NumberTCHEbTags + 1 + NumberTCHPbTags+1 + NumberSSVHEbTags + 1 + NumberSSVHPbTags + 1 + NumberCSVbTags + 1;
  
  std::string bTagFileOutput[TotalNumberbTags];
  std::string PresentationOutput[TotalNumberbTags];
  int NumberRemainingEvents[TotalNumberbTags][inputWTree.size()];
  int NumberRemainingEventsOrigKins[TotalNumberbTags][inputWTree.size()];
  int NumberBLeptCorrectEvents[TotalNumberbTags][inputWTree.size()];
  int NumberBLeptCorrectEventsOrigKins[TotalNumberbTags][inputWTree.size()];
  //Initialize:
  for(int ii=0;ii<14;ii++){
    for(int jj=0;jj<14;jj++){
      for(int kk=0;kk<14;kk++){
	for(int ll=0;ll<14;ll++){
	  for(int nn=0;nn<14;nn++){
	    bTagFileOutput[ii+jj+kk+ll+nn]=" Wrong entry chosen";
	    PresentationOutput[ii+jj+kk+ll+nn]=" Wrong entry chosen";
	    for(int mm=0;mm<inputWTree.size();mm++){
	      NumberRemainingEvents[ii+jj+kk+ll+nn][mm]=0;
	      NumberRemainingEventsOrigKins[ii+jj+kk+ll+nn][mm]=0;
	      NumberBLeptCorrectEvents[ii+jj+kk+ll+nn][mm]=0;
	      NumberBLeptCorrectEventsOrigKins[ii+jj+kk+ll+nn][mm]=0;
	    }
	  }
	}
      }
    }
  }

  //Call classes made :
  BTagName_CSV bTagName = BTagName_CSV();  //for bTagFileOutput name giving
  BTagJetSelection_CSV bTagJetSelection = BTagJetSelection_CSV();
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

    vector<float> CosThetaValues[TotalNumberbTags];
    vector<float> CosThetaValuesOrigKins[TotalNumberbTags];
    vector<float> LumiWeightVector[TotalNumberbTags];
    vector<float> LumiWeightVectorOrigKins[TotalNumberbTags];
    int FilledEntries = 0;
    vector<double> CosThGen[TotalNumberbTags];
    vector<double> CosThGenOrigKins[TotalNumberbTags];
    vector<double> EventCorrectionWeight[TotalNumberbTags];
    vector<double> EventCorrWeightOrigKins[TotalNumberbTags];
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
       
    if(iDataSet==0){
      while(CSVLoop<=NumberCSVbTags){
	while(SSVHPLoop<=NumberSSVHPbTags){  
	  while(SSVHELoop<=NumberSSVHEbTags){
	    while(TCHPLoop<=NumberTCHPbTags){
	      while(TCHELoop<=NumberTCHEbTags){
		bTagFileOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = bTagName.NameGiving(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags);
		PresentationOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = bTagName.NameGivingPres(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags);
		TCHELoop++;
		if(TCHELoop == 14){TCHPLoop=2;SSVHELoop=2;SSVHPLoop=2;CSVLoop=2;}
	      }
	      bTagFileOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = bTagName.NameGiving(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags);
	      PresentationOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = bTagName.NameGivingPres(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags);
	      TCHPLoop++;
	    }
	    bTagFileOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = bTagName.NameGiving(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags);
	    PresentationOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = bTagName.NameGivingPres(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags);
	    SSVHELoop++;
	  }
	  bTagFileOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = bTagName.NameGiving(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags);
	  PresentationOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = bTagName.NameGivingPres(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags);
	  SSVHPLoop++;		
	}	
	bTagFileOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1] = bTagName.NameGiving(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags);
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

      //    for(unsigned int iEvt=0; iEvt<10000; iEvt++){ nEvent = 10000; //nEvent and end of iEvt loop needs to be the same for correctly performing the Minuit Fitter

      inWTreeTree->GetEvent(iEvt);
      if(iEvt%10000 == 0)
        std::cout<<"Processing the "<<iEvt<<"th event, time = "<< ((double)clock() - start) / CLOCKS_PER_SEC <<endl;
      
      // PU reweighting???
      float avPU = ( (float)wTree->nPUBXm1() + (float)wTree->nPU() + (float)wTree->nPUBXp1() ) / 3.; // average in 3 BX!!!, as recommended
      //float lumiWeight = LumiWeights.ITweight( wTree->nPU() );
      
      double lumiWeight3D;
      if(dataSetName.find("Data") != 0){
	lumiWeight3D = Lumi3DWeights.weight3D(wTree->nPUBXm1(),wTree->nPU(),wTree->nPUBXp1());
      }

      //if(dataSetName.find("Data") == 0) lumiWeight = 1;
      if(dataSetName.find("Data") == 0) lumiWeight3D = 1;
      //if( dataSetName.find("Fall10") != string::npos ) lumiWeight = 1; //no PU in Fall10!
      //histo1D["lumiWeights"]->Fill(lumiWeight);
      histo1D["lumiWeights3D"]->Fill(lumiWeight3D);
      
      // scale factor for the event
      float scaleFactor = 1.; 

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
      
      // ChiSquared values and corresponding Leptonic b-jet index from Kinematic Fit (on hadronic side only)
      //      vector<float> ChiSquaredValue; 
      
      /////////////////////////////////////////////////////////////////////////////
      //Read in Chi squared and particles after Hadronic and Leptonic KinFit:
      /////////////////////////////////////////////////////////////////////////////
      
      vector<float> ChiSqKinFit;
      
      vector<TLorentzVector> selectedJets = wTree->selectedJets();
      TLorentzVector lepton = wTree->muon(); //Needed for offline muon cut!  
      //if(dataSetName.find("QCD") !=0){
	MSPlot["LeptonMass"]->Fill(lepton.M(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["LeptonPx"]->Fill(lepton.Px(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["LeptonPy"]->Fill(lepton.Py(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	MSPlot["LeptonPz"]->Fill(lepton.Pz(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	//}

      TLorentzVector MET = wTree->met();
                  
      vector<TLorentzVector> leptBJetKinFit;
      vector<TLorentzVector> leptonKinFit;
      vector<TLorentzVector> neutrinoKinFit;
      vector<TLorentzVector> hadrBJetKinFit;
      vector<TLorentzVector> light1JetKinFit;
      vector<TLorentzVector> light2JetKinFit;
      
      //Use Original kinematics together with changed kinematics to understand the origin of the data-mc discrepancy seen in cos theta* distribution
      vector<TLorentzVector> leptBJetOrig;
      vector<TLorentzVector> leptonOrig;
      vector<TLorentzVector> neutrinoOrig;
      vector<float> ChiSqKinFitOrig;
      for(unsigned int iCombi=0; iCombi<12; iCombi++){      	
	ChiSqKinFitOrig.push_back(wTree->chi2KinFit(iCombi));  //Chi squared obtained by performing a KinFit on hadronic side of the event only. Use original kinematics to reconstruct the event to avoid too much influence from the Kinematic Fit.
      }

      //--------------------------------
      //  Use original kinematics : 
      //--------------------------------
      if(UseChangedKinematics == false){
	//ooooooOOOOOOOOOOOooooooooooOOOOOOOOOoooooooo
	// Apply KinFit on hadronic + leptonic side: 
	//ooooooOOOOOOOOOOOooooooooooOOOOOOOOOoooooooo
	if(LeptonicFit == true){
	  if(iEvt == 0) std::cout << " -- Performing analysis with original kinematics before Hadronic + Leptonic KinFit " << endl;
	  //selectedJets = wTree->selectedJets();
	  //lepton = wTree->muon();
	  //MET = wTree->met();
	  
	  if(LeptTopMassConstrained == true){
	    if(iEvt == 0) std::cout << "   --  Performing analysis with original kinematics before KinFit on hadronic + leptonic side with theoretical mass constraint on leptonic top !!!! " << std::endl;
	    for(unsigned int iCombi=0; iCombi<12; iCombi++){      	
	      ChiSqKinFit.push_back(wTree->chi2FullKinFit(iCombi));   
	    }
	  }
	  else{
	    if(iEvt == 0) std::cout << "   --  Performing analysis with original kinematics before KinFit on hadronic + leptonic side with leptonic top mass left free !!!! " << std::endl;
	    for(unsigned int iCombi=0; iCombi<12; iCombi++){      	
	      ChiSqKinFit.push_back(wTree->chi2FullKinFitMassFit(iCombi));   
	    }
	  }
	}
	//ooooooOOOOOOOOOOOooooooooooOOOOOOOOOOO
	// Apply KinFit on hadronic side only: 
	//ooooooOOOOOOOOOOOooooooooooOOOOOOOOOOO
	else if(LeptonicFit == false){
	  if(iEvt == 0) std::cout << " -- Performing analysis with original kinematics before Hadronic KinFit " << endl;
	  //selectedJets = wTree->selectedJets();
	  //lepton = wTree->muon();
	  //MET = wTree->met();
	  
	  if(iEvt == 0) std::cout << "   --  Performing analysis with original kinematics before KinFit on hadronic side with theoretical mass constraint !!!! " << std::endl;
	  for(unsigned int iCombi=0; iCombi<12; iCombi++){      	
	    ChiSqKinFit.push_back(wTree->chi2KinFit(iCombi));   
	  }	  	  
	}
      }
      //--------------------------------
      //  Use changed kinematics : 
      //--------------------------------
      else if(UseChangedKinematics == true){   
	//ooooooOOOOOOOOOOOooooooooooOOOOOOOOOoooooooo
	// Apply KinFit on hadronic + leptonic side: 
	//ooooooOOOOOOOOOOOooooooooooOOOOOOOOOoooooooo
	if(LeptonicFit == true){
	  if(LeptTopMassConstrained == true){
	    if(iEvt == 0) std::cout << "   --  Performing analysis with kinematics changed after KinFit on hadronic + leptonic side with theoretical mass constraint on leptonic top !!!! " << std::endl;
	    for(unsigned int iCombi=0; iCombi<12; iCombi++){      	
	      ChiSqKinFit.push_back(wTree->chi2FullKinFit(iCombi));   
	      leptBJetKinFit.push_back(wTree->fittedFullLeptB(iCombi));
	      leptonKinFit.push_back(wTree->fittedFullLepton(iCombi));
	      neutrinoKinFit.push_back(wTree->fittedFullNeutrino(iCombi));
	      hadrBJetKinFit.push_back(wTree->fittedFullHadrB(iCombi));
	      light1JetKinFit.push_back(wTree->fittedFullLight1(iCombi));
	      light2JetKinFit.push_back(wTree->fittedFullLight2(iCombi));
	    }
	  }
	  else{
	    if(iEvt == 0) std::cout << "   --  Performing analysis with kinematics changed after KinFit on hadronic + leptonic side with leptonic top mass left free !!!! " << std::endl;
	    for(unsigned int iCombi=0; iCombi<12; iCombi++){      		    
	      ChiSqKinFit.push_back(wTree->chi2FullKinFitMassFit(iCombi));   
	      leptBJetKinFit.push_back(wTree->fittedFullLeptBMassFit(iCombi));
	      leptonKinFit.push_back(wTree->fittedFullLeptonMassFit(iCombi));
	      neutrinoKinFit.push_back(wTree->fittedFullNeutrinoMassFit(iCombi));   
	      hadrBJetKinFit.push_back(wTree->fittedFullHadrB(iCombi));
	      light1JetKinFit.push_back(wTree->fittedFullLight1(iCombi));
	      light2JetKinFit.push_back(wTree->fittedFullLight2(iCombi));
	    }	  
	  }	
	}
	//ooooooOOOOOOOOOOOooooooooooOOOOOOOOOOO
	// Apply KinFit on hadronic side only: 
	//ooooooOOOOOOOOOOOooooooooooOOOOOOOOOOO	
	else if(LeptonicFit == false){
	  if(LeptTopMassConstrained == true){
	    if(iEvt == 0) std::cout << "   --  Performing analysis with kinematics changed after KinFit on hadronic side with theoretical mass constraint !!!! " << std::endl;
	    for(unsigned int iCombi=0; iCombi<12; iCombi++){      	
	      ChiSqKinFit.push_back(wTree->chi2KinFit(iCombi));   
	      leptBJetKinFit.push_back(wTree->fittedLeptB(iCombi));
	      leptonKinFit.push_back(wTree->fittedLepton(iCombi));
	      neutrinoKinFit.push_back(wTree->fittedNeutrino(iCombi));   
	      hadrBJetKinFit.push_back(wTree->fittedFullHadrB(iCombi));
	      light1JetKinFit.push_back(wTree->fittedFullLight1(iCombi));
	      light2JetKinFit.push_back(wTree->fittedFullLight2(iCombi));
	    }
	  }
	}
      }
      
      //ooooooooOOOOOOOOOOOOooooooooooooOOOOOOOOOOOOooooooooooooOOOOO
      //ooOOooOOoo      Reading out nTuples done           ooOOooOOoo
      //ooOOooOOoo-----------------------------------------ooOOooOOoo
      //ooOOooOOoo      Start of actual analysis           ooOOooOOoo
      //ooooooooOOOOOOOOOOOOooooooooooooOOOOOOOOOOOOooooooooooooOOOOO
      float TransverseMass = sqrt(2*(abs(lepton.Pt()))*abs(MET.Pt())*(1-cos(lepton.DeltaPhi(MET))));	//Should this not be updated to correct lepton and MET??
      MSPlot["TransverseMassBeforeCut"]->Fill(TransverseMass, datasets[iDataSet], true, Luminosity);	

      //float reliso = (lepton.chargedHadronIso()+lepton.neutralHadronIso()+lepton.photonIso())/lepton.Pt();
      //cout << " relIso = " << reliso << endl;

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
	  
	  if(TransverseMass < 30) cout << " Event should be removed by MTCut " << endl;		    
	  
	  selecTableMacro.Fill(iDataSet,3,scaleFactor*lumiWeight3D);
	  
	  if((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ) NumberEventsAfterTransverseMass++;
	  
	  MSPlot["TransverseMassAfterCut"]->Fill(TransverseMass, datasets[iDataSet], true, Luminosity);	      	  
	  MSPlot["nPVAfterTransverseMassCut"]->Fill(nPrimaryVertices,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	  
	  eventSelected = true;
	}
      }
      
      if(eventSelected == true){

	//Check whether mass constraints applied in Kinematic Fit are recovered in mass distributions:
	for(int ii=0;ii<12;ii++){
	  if(ChiSqKinFit[ii] != 9999){  //Only study converging events!
	    if((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ){
	      histo1D["HadronicWMass"]->Fill((light1JetKinFit[ii]+light2JetKinFit[ii]).M());
	      histo1D["HadronicTopMass"]->Fill((light1JetKinFit[ii]+light2JetKinFit[ii]+hadrBJetKinFit[ii]).M());
	      histo1D["LeptonicWMass"]->Fill((leptonKinFit[ii]+neutrinoKinFit[ii]).M());
	      histo1D["LeptonicTopMass"]->Fill((leptonKinFit[ii]+neutrinoKinFit[ii]+leptBJetKinFit[ii]).M());
	    }
	    MSPlot["HadronicWMass"]->Fill((light1JetKinFit[ii]+light2JetKinFit[ii]).M(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	    MSPlot["HadronicTopMass"]->Fill((light1JetKinFit[ii]+light2JetKinFit[ii]+hadrBJetKinFit[ii]).M(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	    MSPlot["LeptonicWMass"]->Fill((leptonKinFit[ii]+neutrinoKinFit[ii]).M(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	    MSPlot["LeptonicTopMass"]->Fill((leptonKinFit[ii]+neutrinoKinFit[ii]+leptBJetKinFit[ii]).M(),datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
	  }
	}

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
	  
	//Obtain probabilities:
	vector<float> ProbabilityOfKinFit;	 
	vector<float> KinFitProbOrigKins;
	vector<float> ProbSOverSB;                 //Need to change WTree classes such that Mlb method is removed !!!!
	for(int jj = 0; jj < 12 ; jj++){ 
	  if(ChiSqKinFit[jj]!=9999){ ProbabilityOfKinFit.push_back(TMath::Prob(ChiSqKinFit[jj],2));}
	  else{ProbabilityOfKinFit.push_back(-1.);}

	  if(ChiSqKinFitOrig[jj] != 9999){ KinFitProbOrigKins.push_back(TMath::Prob(ChiSqKinFitOrig[jj],2));}
	  else{KinFitProbOrigKins.push_back(-1.);}

	  ProbSOverSB.push_back( 1.);  	  
	}	
	
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

	//-------------------------------------------------------------
	//Obtain jet combination for the different b-tag constraints:
	//------------------------------------------------------------      		
	
	//Initialize bTag loop variables:
	int TCHEbTagLoop=1;
	int TCHPbTagLoop=1; 
	int SSVHEbTagLoop=1;
	int SSVHPbTagLoop=1;
	int CSVbTagLoop=1;
      
	int ConsideredBTagger; //0=tche, 1 = tchp, 2 = ssvhe, 3 = ssvhp & 4 = csv

	//Create some alternative helicites weights to obtain different ratios:
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
	
	while(CSVbTagLoop <= NumberCSVbTags){
	  while(SSVHPbTagLoop <=NumberSSVHPbTags ){
	    while(SSVHEbTagLoop <=NumberSSVHEbTags ){
	      while(TCHPbTagLoop <= NumberTCHPbTags){
		while(TCHEbTagLoop <= NumberTCHEbTags){	   
		  ConsideredBTagger=0;
		  
		  int JetCombination = bTagJetSelection.HighestProbSelection(TCHEbTagLoop,ConsideredBTagger,ProbabilityOfKinFit,ProbSOverSB,btagTCHE,btagTCHP,btagSSVHE,btagSSVHP,btagCSV);

		  if(TCHEbTagLoop == 1){
		    MSPlot["nPVBeforeFoundJetComb"]->Fill(nPrimaryVertices,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
		  }
		  
		  int BLeptonicIndex = BLeptIndex[JetCombination];
		  int BHadronicIndex = BHadrIndex[JetCombination];
		    
		  if(JetCombination != 999){//&& JetCombination == CorrectKinFitIndex){
		    
		    if(TCHEbTagLoop == 1){MSPlot["nPVAfterFoundJetComb"]->Fill(nPrimaryVertices,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);}
		    if(TCHEbTagLoop == 7){MSPlot["nPVAfterFoundJetCombbTag"]->Fill(nPrimaryVertices,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);}
		    
		    float CosThetaCalculation;
		    if(UseChangedKinematics == true) CosThetaCalculation = bTagCosThetaCalculation.Calculation(leptonKinFit[JetCombination],neutrinoKinFit[JetCombination],leptBJetKinFit[JetCombination]);
		    else if(UseChangedKinematics == false) CosThetaCalculation = bTagCosThetaCalculation.CalcOrigKins(BLeptonicIndex,BHadronicIndex,lepton,selectedJets,MassW, MassTop);
		
		    if(CosThetaCalculation != 999){	

		      //Resolution histograms:
		      if((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ){
			if(TCHEbTagLoop == 3){ histo1D["CosThetaResTCHELLept"]->Fill(CosThetaCalculation - CosThetaGenerator);}
			if(TCHEbTagLoop == 7){ histo1D["CosThetaResTCHEMLept"]->Fill(CosThetaCalculation - CosThetaGenerator);}
		      }

		      //No btag case: 
		      if(TCHEbTagLoop == 1){
			MSPlot["nPVAfterFoundCosTheta"]->Fill(nPrimaryVertices,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);

			//Obtain cos theta* histo's for alternative helcities to investigate ratio with StandardModel:
			if(dataSetName.find("Data") != 0){	
			  for(int ii=0;ii<SizeArray;ii++){
			    histo1D["CosThetaRight"+HelicityNumbers[ii][1].str()+"Long"+HelicityNumbers[ii][0].str()+"Left"+HelicityNumbers[ii][2].str()]->Fill(CosThetaCalculation,(HelicityWeight[ii]*Luminosity*scaleFactor*lumiWeight3D*datasets[iDataSet]->NormFactor()));  
			  }
			}
			
			//Information about signal dataset:
			if((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ){
			  NumberUsedEvents++;
			  
			  //Histos using all signal events:
			  histo1D["CosThetaRes"]->Fill(CosThetaCalculation - CosThetaGenerator); 
			  histo1D["CosThetaDistributionSignal"]->Fill(CosThetaCalculation);
			  histo1D["CosThetaDistributionSignalRatio"]->Fill(CosThetaCalculation); //Extra histo to avoid double naming when calculating ratio!!
			  histo1D["CosThetaDistributionSignalGen"]->Fill(CosThetaGenerator);
			  histo1D["MuonPtSignal"]->Fill(lepton.Pt());
			  CosThetaDiffVSRecoSignal->Fill(CosThetaCalculation, (CosThetaCalculation - CosThetaGenerator) );
			  
			  if(CosThetaCalculation <= - 0.8) histo1D["CosThetaResLastBin"]->Fill(CosThetaCalculation - CosThetaGenerator);			
			  if(CosThetaCalculation > -0.8 && CosThetaCalculation <= -0.6) histo1D["CosThetaResSecondNegBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			  if(CosThetaCalculation > -0.6 && CosThetaCalculation <= -0.4) histo1D["CosThetaResThirdNegBin"]->Fill(CosThetaCalculation - CosThetaGenerator); 
			  if(CosThetaCalculation > -0.4 && CosThetaCalculation <= -0.2) histo1D["CosThetaResFourthNegBin"]->Fill(CosThetaCalculation - CosThetaGenerator); 
			  if(CosThetaCalculation > -0.2 && CosThetaCalculation <= -0.) histo1D["CosThetaResFifthNegBin"]->Fill(CosThetaCalculation - CosThetaGenerator); 
			  //if(CosThetaCalculation > -0.5 && CosThetaCalculation <= -0.4) histo1D["CosThetaResSixtNegBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			  //if(CosThetaCalculation > -0.4 && CosThetaCalculation <= -0.3) histo1D["CosThetaResSeventhNegBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			  //if(CosThetaCalculation > -0.3 && CosThetaCalculation <= -0.2) histo1D["CosThetaResEighthNegBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			  //if(CosThetaCalculation > -0.2 && CosThetaCalculation <= -0.1) histo1D["CosThetaResNinthNegBin"]->Fill(CosThetaCalculation - CosThetaGenerator); 
			  //if(CosThetaCalculation > -0.1 && CosThetaCalculation <= 0.) histo1D["CosThetaResTenthNegBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			  //if(CosThetaCalculation > 0. && CosThetaCalculation <= 0.1) histo1D["CosThetaResTenthBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			  //if(CosThetaCalculation > 0.1 && CosThetaCalculation <= 0.2) histo1D["CosThetaResNinthBin"]->Fill(CosThetaCalculation - CosThetaGenerator); 
			  //if(CosThetaCalculation > 0.2 && CosThetaCalculation <= 0.3) histo1D["CosThetaResEighthBin"]->Fill(CosThetaCalculation - CosThetaGenerator); 
			  //if(CosThetaCalculation > 0.3 && CosThetaCalculation <= 0.4) histo1D["CosThetaResSeventhBin"]->Fill(CosThetaCalculation - CosThetaGenerator); 
			  //if(CosThetaCalculation > 0.4 && CosThetaCalculation <= 0.5) histo1D["CosThetaResSixtBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			  if(CosThetaCalculation > 0. && CosThetaCalculation <= 0.2) histo1D["CosThetaResFifthBin"]->Fill(CosThetaCalculation - CosThetaGenerator); 
			  if(CosThetaCalculation > 0.2 && CosThetaCalculation <= 0.4) histo1D["CosThetaResFourthBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			  if(CosThetaCalculation > 0.4 && CosThetaCalculation <= 0.6) histo1D["CosThetaResThirdBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			  if(CosThetaCalculation > 0.6 && CosThetaCalculation <= 0.8) histo1D["CosThetaResSecondBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			  if(CosThetaCalculation > 0.8) histo1D["CosThetaResFirstBin"]->Fill(CosThetaCalculation - CosThetaGenerator);

			  //Histos using different muon Pt to determine optimal muon pt cut:
			  if(leptonKinFit[JetCombination].Pt()>=10){ //No real difference in using leptonKinFit or muon since muon doesn't change much with KinFit!!
			    histo1D["CosThetaDifferencePt10ToInf"]->Fill(CosThetaCalculation - CosThetaGenerator); 
			    histo1D["CosThetaRecoPt10ToInf"]->Fill(CosThetaCalculation);
			    histo1D["CosThetaRecoPt10ToInfForRatio"]->Fill(CosThetaCalculation);
			    histo1D["CosThetaGenPt10ToInf"]->Fill(CosThetaGenerator);
			    CosThetaDiffVSRecoPt10ToInf->Fill(CosThetaCalculation, (CosThetaCalculation - CosThetaGenerator) );
			    histo2D["CosThetaGenVsRecoPt10ToInf"]->Fill(CosThetaCalculation,CosThetaGenerator);

			    if(CosThetaGenerator <= -0.8){
			      TotalPt10ToInf[9]++;
			      histo1D["CosThetaDiffPt10ToInfGenTenthBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt10ToInfGenTenthBin"]->Fill(CosThetaCalculation);
			      if(CosThetaCalculation <= -0.8) ProbabilityPt10ToInf[9][9]++;
			      if(CosThetaCalculation > -0.8 && CosThetaCalculation <= -0.6) ProbabilityPt10ToInf[9][8]++;
			      if(CosThetaCalculation > -0.6 && CosThetaCalculation <= -0.4) ProbabilityPt10ToInf[9][7]++;
			      if(CosThetaCalculation > -0.4 && CosThetaCalculation <= -0.2) ProbabilityPt10ToInf[9][6]++;
			      if(CosThetaCalculation > -0.2 && CosThetaCalculation <= 0.) ProbabilityPt10ToInf[9][5]++;
			      if(CosThetaCalculation > 0. && CosThetaCalculation <= 0.2) ProbabilityPt10ToInf[9][4]++;
			      if(CosThetaCalculation > 0.2 && CosThetaCalculation <= 0.4) ProbabilityPt10ToInf[9][3]++;
			      if(CosThetaCalculation > 0.4 && CosThetaCalculation <= 0.6) ProbabilityPt10ToInf[9][2]++;
			      if(CosThetaCalculation > 0.6 && CosThetaCalculation <= 0.8) ProbabilityPt10ToInf[9][1]++;
			      if(CosThetaCalculation > 0.8 && CosThetaCalculation <= 1.) ProbabilityPt10ToInf[9][0]++;			      
			    }	
			    else if(CosThetaGenerator > -0.8 && CosThetaGenerator <= -0.6){
			      TotalPt10ToInf[8]++;
			      histo1D["CosThetaDiffPt10ToInfGenNinthBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt10ToInfGenNinthBin"]->Fill(CosThetaCalculation);
			      if(CosThetaCalculation <= -0.8) ProbabilityPt10ToInf[8][9]++;
			      if(CosThetaCalculation > -0.8 && CosThetaCalculation <= -0.6) ProbabilityPt10ToInf[8][8]++;
			      if(CosThetaCalculation > -0.6 && CosThetaCalculation <= -0.4) ProbabilityPt10ToInf[8][7]++;
			      if(CosThetaCalculation > -0.4 && CosThetaCalculation <= -0.2) ProbabilityPt10ToInf[8][6]++;
			      if(CosThetaCalculation > -0.2 && CosThetaCalculation <= 0.) ProbabilityPt10ToInf[8][5]++;
			      if(CosThetaCalculation > 0. && CosThetaCalculation <= 0.2) ProbabilityPt10ToInf[8][4]++;
			      if(CosThetaCalculation > 0.2 && CosThetaCalculation <= 0.4) ProbabilityPt10ToInf[8][3]++;
			      if(CosThetaCalculation > 0.4 && CosThetaCalculation <= 0.6) ProbabilityPt10ToInf[8][2]++;
			      if(CosThetaCalculation > 0.6 && CosThetaCalculation <= 0.8) ProbabilityPt10ToInf[8][1]++;
			      if(CosThetaCalculation > 0.8 && CosThetaCalculation <= 1.) ProbabilityPt10ToInf[8][0]++;
			    }	
			    else if(CosThetaGenerator > -0.6 && CosThetaGenerator <= -0.4){
			      TotalPt10ToInf[7]++;
			      histo1D["CosThetaDiffPt10ToInfGenEighthBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt10ToInfGenEighthBin"]->Fill(CosThetaCalculation);
			      if(CosThetaCalculation <= -0.8) ProbabilityPt10ToInf[7][9]++;
			      if(CosThetaCalculation > -0.8 && CosThetaCalculation <= -0.6) ProbabilityPt10ToInf[7][8]++;
			      if(CosThetaCalculation > -0.6 && CosThetaCalculation <= -0.4) ProbabilityPt10ToInf[7][7]++;
			      if(CosThetaCalculation > -0.4 && CosThetaCalculation <= -0.2) ProbabilityPt10ToInf[7][6]++;
			      if(CosThetaCalculation > -0.2 && CosThetaCalculation <= 0.) ProbabilityPt10ToInf[7][5]++;
			      if(CosThetaCalculation > 0. && CosThetaCalculation <= 0.2) ProbabilityPt10ToInf[7][4]++;
			      if(CosThetaCalculation > 0.2 && CosThetaCalculation <= 0.4) ProbabilityPt10ToInf[7][3]++;
			      if(CosThetaCalculation > 0.4 && CosThetaCalculation <= 0.6) ProbabilityPt10ToInf[7][2]++;
			      if(CosThetaCalculation > 0.6 && CosThetaCalculation <= 0.8) ProbabilityPt10ToInf[7][1]++;
			      if(CosThetaCalculation > 0.8 && CosThetaCalculation <= 1.) ProbabilityPt10ToInf[7][0]++;
			    }	
			    else if(CosThetaGenerator > -0.4 && CosThetaGenerator <= -0.2){
			      TotalPt10ToInf[6]++;
			      histo1D["CosThetaDiffPt10ToInfGenSeventhBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt10ToInfGenSeventhBin"]->Fill(CosThetaCalculation);
			      if(CosThetaCalculation <= -0.8) ProbabilityPt10ToInf[6][9]++;
			      if(CosThetaCalculation > -0.8 && CosThetaCalculation <= -0.6) ProbabilityPt10ToInf[6][8]++;
			      if(CosThetaCalculation > -0.6 && CosThetaCalculation <= -0.4) ProbabilityPt10ToInf[6][7]++;
			      if(CosThetaCalculation > -0.4 && CosThetaCalculation <= -0.2) ProbabilityPt10ToInf[6][6]++;
			      if(CosThetaCalculation > -0.2 && CosThetaCalculation <= 0.) ProbabilityPt10ToInf[6][5]++;
			      if(CosThetaCalculation > 0. && CosThetaCalculation <= 0.2) ProbabilityPt10ToInf[6][4]++;
			      if(CosThetaCalculation > 0.2 && CosThetaCalculation <= 0.4) ProbabilityPt10ToInf[6][3]++;
			      if(CosThetaCalculation > 0.4 && CosThetaCalculation <= 0.6) ProbabilityPt10ToInf[6][2]++;
			      if(CosThetaCalculation > 0.6 && CosThetaCalculation <= 0.8) ProbabilityPt10ToInf[6][1]++;
			      if(CosThetaCalculation > 0.8 && CosThetaCalculation <= 1.) ProbabilityPt10ToInf[6][0]++;
			    }	
			    else if(CosThetaGenerator > -0.2 && CosThetaGenerator <= 0.){
			      TotalPt10ToInf[5]++;
			      histo1D["CosThetaDiffPt10ToInfGenSixthBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt10ToInfGenSixthBin"]->Fill(CosThetaCalculation);
			      if(CosThetaCalculation <= -0.8) ProbabilityPt10ToInf[5][9]++;
			      if(CosThetaCalculation > -0.8 && CosThetaCalculation <= -0.6) ProbabilityPt10ToInf[5][8]++;
			      if(CosThetaCalculation > -0.6 && CosThetaCalculation <= -0.4) ProbabilityPt10ToInf[5][7]++;
			      if(CosThetaCalculation > -0.4 && CosThetaCalculation <= -0.2) ProbabilityPt10ToInf[5][6]++;
			      if(CosThetaCalculation > -0.2 && CosThetaCalculation <= 0.) ProbabilityPt10ToInf[5][5]++;
			      if(CosThetaCalculation > 0. && CosThetaCalculation <= 0.2) ProbabilityPt10ToInf[5][4]++;
			      if(CosThetaCalculation > 0.2 && CosThetaCalculation <= 0.4) ProbabilityPt10ToInf[5][3]++;
			      if(CosThetaCalculation > 0.4 && CosThetaCalculation <= 0.6) ProbabilityPt10ToInf[5][2]++;
			      if(CosThetaCalculation > 0.6 && CosThetaCalculation <= 0.8) ProbabilityPt10ToInf[5][1]++;
			      if(CosThetaCalculation > 0.8 && CosThetaCalculation <= 1.) ProbabilityPt10ToInf[5][0]++;
			    }	
			    else if(CosThetaGenerator > 0. && CosThetaGenerator <= 0.2){
			      TotalPt10ToInf[4]++;
			      histo1D["CosThetaDiffPt10ToInfGenFifthBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt10ToInfGenFifthBin"]->Fill(CosThetaCalculation);
			      if(CosThetaCalculation <= -0.8) ProbabilityPt10ToInf[4][9]++;
			      if(CosThetaCalculation > -0.8 && CosThetaCalculation <= -0.6) ProbabilityPt10ToInf[4][8]++;
			      if(CosThetaCalculation > -0.6 && CosThetaCalculation <= -0.4) ProbabilityPt10ToInf[4][7]++;
			      if(CosThetaCalculation > -0.4 && CosThetaCalculation <= -0.2) ProbabilityPt10ToInf[4][6]++;
			      if(CosThetaCalculation > -0.2 && CosThetaCalculation <= 0.) ProbabilityPt10ToInf[4][5]++;
			      if(CosThetaCalculation > 0. && CosThetaCalculation <= 0.2) ProbabilityPt10ToInf[4][4]++;
			      if(CosThetaCalculation > 0.2 && CosThetaCalculation <= 0.4) ProbabilityPt10ToInf[4][3]++;
			      if(CosThetaCalculation > 0.4 && CosThetaCalculation <= 0.6) ProbabilityPt10ToInf[4][2]++;
			      if(CosThetaCalculation > 0.6 && CosThetaCalculation <= 0.8) ProbabilityPt10ToInf[4][1]++;
			      if(CosThetaCalculation > 0.8 && CosThetaCalculation <= 1.) ProbabilityPt10ToInf[4][0]++;
			    }	
			    else if(CosThetaGenerator > 0.2 && CosThetaGenerator <= 0.4){
			      TotalPt10ToInf[3]++;
			      histo1D["CosThetaDiffPt10ToInfGenFourthBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt10ToInfGenFourthBin"]->Fill(CosThetaCalculation);
			      if(CosThetaCalculation <= -0.8) ProbabilityPt10ToInf[3][9]++;
			      if(CosThetaCalculation > -0.8 && CosThetaCalculation <= -0.6) ProbabilityPt10ToInf[3][8]++;
			      if(CosThetaCalculation > -0.6 && CosThetaCalculation <= -0.4) ProbabilityPt10ToInf[3][7]++;
			      if(CosThetaCalculation > -0.4 && CosThetaCalculation <= -0.2) ProbabilityPt10ToInf[3][6]++;
			      if(CosThetaCalculation > -0.2 && CosThetaCalculation <= 0.) ProbabilityPt10ToInf[3][5]++;
			      if(CosThetaCalculation > 0. && CosThetaCalculation <= 0.2) ProbabilityPt10ToInf[3][4]++;
			      if(CosThetaCalculation > 0.2 && CosThetaCalculation <= 0.4) ProbabilityPt10ToInf[3][3]++;
			      if(CosThetaCalculation > 0.4 && CosThetaCalculation <= 0.6) ProbabilityPt10ToInf[3][2]++;
			      if(CosThetaCalculation > 0.6 && CosThetaCalculation <= 0.8) ProbabilityPt10ToInf[3][1]++;
			      if(CosThetaCalculation > 0.8 && CosThetaCalculation <= 1.) ProbabilityPt10ToInf[3][0]++;
			    }	
			    else if(CosThetaGenerator > 0.4 && CosThetaGenerator <= 0.6){
			      TotalPt10ToInf[2]++;
			      histo1D["CosThetaDiffPt10ToInfGenThirdBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt10ToInfGenThirdBin"]->Fill(CosThetaCalculation);
			      if(CosThetaCalculation <= -0.8) ProbabilityPt10ToInf[2][9]++;
			      if(CosThetaCalculation > -0.8 && CosThetaCalculation <= -0.6) ProbabilityPt10ToInf[2][8]++;
			      if(CosThetaCalculation > -0.6 && CosThetaCalculation <= -0.4) ProbabilityPt10ToInf[2][7]++;
			      if(CosThetaCalculation > -0.4 && CosThetaCalculation <= -0.2) ProbabilityPt10ToInf[2][6]++;
			      if(CosThetaCalculation > -0.2 && CosThetaCalculation <= 0.) ProbabilityPt10ToInf[2][5]++;
			      if(CosThetaCalculation > 0. && CosThetaCalculation <= 0.2) ProbabilityPt10ToInf[2][4]++;
			      if(CosThetaCalculation > 0.2 && CosThetaCalculation <= 0.4) ProbabilityPt10ToInf[2][3]++;
			      if(CosThetaCalculation > 0.4 && CosThetaCalculation <= 0.6) ProbabilityPt10ToInf[2][2]++;
			      if(CosThetaCalculation > 0.6 && CosThetaCalculation <= 0.8) ProbabilityPt10ToInf[2][1]++;
			      if(CosThetaCalculation > 0.8 && CosThetaCalculation <= 1.) ProbabilityPt10ToInf[2][0]++;
			    }	
			    else if(CosThetaGenerator > 0.6 && CosThetaGenerator <= 0.8){
			      TotalPt10ToInf[1]++;
			      histo1D["CosThetaDiffPt10ToInfGenSecondBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt10ToInfGenSecondBin"]->Fill(CosThetaCalculation);
			      if(CosThetaCalculation <= -0.8) ProbabilityPt10ToInf[1][9]++;
			      if(CosThetaCalculation > -0.8 && CosThetaCalculation <= -0.6) ProbabilityPt10ToInf[1][8]++;
			      if(CosThetaCalculation > -0.6 && CosThetaCalculation <= -0.4) ProbabilityPt10ToInf[1][7]++;
			      if(CosThetaCalculation > -0.4 && CosThetaCalculation <= -0.2) ProbabilityPt10ToInf[1][6]++;
			      if(CosThetaCalculation > -0.2 && CosThetaCalculation <= 0.) ProbabilityPt10ToInf[1][5]++;
			      if(CosThetaCalculation > 0. && CosThetaCalculation <= 0.2) ProbabilityPt10ToInf[1][4]++;
			      if(CosThetaCalculation > 0.2 && CosThetaCalculation <= 0.4) ProbabilityPt10ToInf[1][3]++;
			      if(CosThetaCalculation > 0.4 && CosThetaCalculation <= 0.6) ProbabilityPt10ToInf[1][2]++;
			      if(CosThetaCalculation > 0.6 && CosThetaCalculation <= 0.8) ProbabilityPt10ToInf[1][1]++;
			      if(CosThetaCalculation > 0.8 && CosThetaCalculation <= 1.) ProbabilityPt10ToInf[1][0]++;
			    }	
			    else if(CosThetaGenerator > 0.8 && CosThetaGenerator <= 1.){
			      TotalPt10ToInf[0]++;
			      histo1D["CosThetaDiffPt10ToInfGenFirstBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt10ToInfGenFirstBin"]->Fill(CosThetaCalculation);
			      if(CosThetaCalculation <= -0.8) ProbabilityPt10ToInf[0][9]++;
			      if(CosThetaCalculation > -0.8 && CosThetaCalculation <= -0.6) ProbabilityPt10ToInf[0][8]++;
			      if(CosThetaCalculation > -0.6 && CosThetaCalculation <= -0.4) ProbabilityPt10ToInf[0][7]++;
			      if(CosThetaCalculation > -0.4 && CosThetaCalculation <= -0.2) ProbabilityPt10ToInf[0][6]++;
			      if(CosThetaCalculation > -0.2 && CosThetaCalculation <= 0.) ProbabilityPt10ToInf[0][5]++;
			      if(CosThetaCalculation > 0. && CosThetaCalculation <= 0.2) ProbabilityPt10ToInf[0][4]++;
			      if(CosThetaCalculation > 0.2 && CosThetaCalculation <= 0.4) ProbabilityPt10ToInf[0][3]++;
			      if(CosThetaCalculation > 0.4 && CosThetaCalculation <= 0.6) ProbabilityPt10ToInf[0][2]++;
			      if(CosThetaCalculation > 0.6 && CosThetaCalculation <= 0.8) ProbabilityPt10ToInf[0][1]++;
			      if(CosThetaCalculation > 0.8 && CosThetaCalculation <= 1.) ProbabilityPt10ToInf[0][0]++;
			    }	
			      
			  }
			  if(leptonKinFit[JetCombination].Pt()>=15){
			    histo1D["CosThetaDifferencePt15ToInf"]->Fill(CosThetaCalculation - CosThetaGenerator); 
			    histo1D["CosThetaRecoPt15ToInf"]->Fill(CosThetaCalculation);
			    histo1D["CosThetaRecoPt15ToInfForRatio"]->Fill(CosThetaCalculation);
			    histo1D["CosThetaGenPt15ToInf"]->Fill(CosThetaGenerator);
			    CosThetaDiffVSRecoPt15ToInf->Fill(CosThetaCalculation, (CosThetaCalculation - CosThetaGenerator) );
			    histo2D["CosThetaGenVsRecoPt15ToInf"]->Fill(CosThetaCalculation,CosThetaGenerator);

			    if(CosThetaGenerator <= -0.8){
			      histo1D["CosThetaDiffPt15ToInfGenTenthBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt15ToInfGenTenthBin"]->Fill(CosThetaCalculation);
			    }	
			    else if(CosThetaGenerator > -0.8 && CosThetaGenerator <= -0.6){
			      histo1D["CosThetaDiffPt15ToInfGenNinthBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt15ToInfGenNinthBin"]->Fill(CosThetaCalculation);
			    }	
			    else if(CosThetaGenerator > -0.6 && CosThetaGenerator <= -0.4){
			      histo1D["CosThetaDiffPt15ToInfGenEighthBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt15ToInfGenEighthBin"]->Fill(CosThetaCalculation);
			    }	
			    else if(CosThetaGenerator > -0.4 && CosThetaGenerator <= -0.2){
			      histo1D["CosThetaDiffPt15ToInfGenSeventhBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt15ToInfGenSeventhBin"]->Fill(CosThetaCalculation);
			    }	
			    else if(CosThetaGenerator > -0.2 && CosThetaGenerator <= 0.){
			      histo1D["CosThetaDiffPt15ToInfGenSixthBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt15ToInfGenSixthBin"]->Fill(CosThetaCalculation);
			    }	
			    else if(CosThetaGenerator > 0. && CosThetaGenerator <= 0.2){
			      histo1D["CosThetaDiffPt15ToInfGenFifthBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt15ToInfGenFifthBin"]->Fill(CosThetaCalculation);
			    }	
			    else if(CosThetaGenerator > 0.2 && CosThetaGenerator <= 0.4){
			      histo1D["CosThetaDiffPt15ToInfGenFourthBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt15ToInfGenFourthBin"]->Fill(CosThetaCalculation);
			    }	
			    else if(CosThetaGenerator > 0.4 && CosThetaGenerator <= 0.6){
			      histo1D["CosThetaDiffPt15ToInfGenThirdBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt15ToInfGenThirdBin"]->Fill(CosThetaCalculation);
			    }	
			    else if(CosThetaGenerator > 0.6 && CosThetaGenerator <= 0.8){
			      histo1D["CosThetaDiffPt15ToInfGenSecondBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt15ToInfGenSecondBin"]->Fill(CosThetaCalculation);
			    }	
			    else if(CosThetaGenerator > 0.8 && CosThetaGenerator <= 1.){
			      histo1D["CosThetaDiffPt15ToInfGenFirstBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt15ToInfGenFirstBin"]->Fill(CosThetaCalculation);
			    }	

			  }
			  if(leptonKinFit[JetCombination].Pt()>=20){
			    histo1D["CosThetaDifferencePt20ToInf"]->Fill(CosThetaCalculation - CosThetaGenerator); 
			    histo1D["CosThetaRecoPt20ToInf"]->Fill(CosThetaCalculation);
			    histo1D["CosThetaRecoPt20ToInfForRatio"]->Fill(CosThetaCalculation);
			    histo1D["CosThetaGenPt20ToInf"]->Fill(CosThetaGenerator);
			    CosThetaDiffVSRecoPt20ToInf->Fill(CosThetaCalculation, (CosThetaCalculation - CosThetaGenerator) );
			    histo2D["CosThetaGenVsRecoPt20ToInf"]->Fill(CosThetaCalculation,CosThetaGenerator);

			    if(CosThetaGenerator <= -0.8){
			      histo1D["CosThetaDiffPt20ToInfGenTenthBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt20ToInfGenTenthBin"]->Fill(CosThetaCalculation);
			    }	
			    else if(CosThetaGenerator > -0.8 && CosThetaGenerator <= -0.6){
			      histo1D["CosThetaDiffPt20ToInfGenNinthBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt20ToInfGenNinthBin"]->Fill(CosThetaCalculation);
			    }	
			    else if(CosThetaGenerator > -0.6 && CosThetaGenerator <= -0.4){
			      histo1D["CosThetaDiffPt20ToInfGenEighthBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt20ToInfGenEighthBin"]->Fill(CosThetaCalculation);
			    }	
			    else if(CosThetaGenerator > -0.4 && CosThetaGenerator <= -0.2){
			      histo1D["CosThetaDiffPt20ToInfGenSeventhBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt20ToInfGenSeventhBin"]->Fill(CosThetaCalculation);
			    }	
			    else if(CosThetaGenerator > -0.2 && CosThetaGenerator <= 0.){
			      histo1D["CosThetaDiffPt20ToInfGenSixthBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt20ToInfGenSixthBin"]->Fill(CosThetaCalculation);
			    }	
			    else if(CosThetaGenerator > 0. && CosThetaGenerator <= 0.2){
			      histo1D["CosThetaDiffPt20ToInfGenFifthBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt20ToInfGenFifthBin"]->Fill(CosThetaCalculation);
			    }	
			    else if(CosThetaGenerator > 0.2 && CosThetaGenerator <= 0.4){
			      histo1D["CosThetaDiffPt20ToInfGenFourthBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt20ToInfGenFourthBin"]->Fill(CosThetaCalculation);
			    }	
			    else if(CosThetaGenerator > 0.4 && CosThetaGenerator <= 0.6){
			      histo1D["CosThetaDiffPt20ToInfGenThirdBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt20ToInfGenThirdBin"]->Fill(CosThetaCalculation);
			    }	
			    else if(CosThetaGenerator > 0.6 && CosThetaGenerator <= 0.8){
			      histo1D["CosThetaDiffPt20ToInfGenSecondBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt20ToInfGenSecondBin"]->Fill(CosThetaCalculation);
			    }	
			    else if(CosThetaGenerator > 0.8 && CosThetaGenerator <= 1.){
			      histo1D["CosThetaDiffPt20ToInfGenFirstBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt20ToInfGenFirstBin"]->Fill(CosThetaCalculation);
			    }	

			  }
			  if(leptonKinFit[JetCombination].Pt()>=25){
			    histo1D["CosThetaDifferencePt25ToInf"]->Fill(CosThetaCalculation - CosThetaGenerator); 
			    histo1D["CosThetaRecoPt25ToInf"]->Fill(CosThetaCalculation);
			    histo1D["CosThetaRecoPt25ToInfForRatio"]->Fill(CosThetaCalculation);
			    histo1D["CosThetaGenPt25ToInf"]->Fill(CosThetaGenerator);
			    CosThetaDiffVSRecoPt25ToInf->Fill(CosThetaCalculation, (CosThetaCalculation - CosThetaGenerator) );
			    histo2D["CosThetaGenVsRecoPt25ToInf"]->Fill(CosThetaCalculation,CosThetaGenerator);

			    if(CosThetaGenerator <= -0.8){
			      histo1D["CosThetaDiffPt25ToInfGenTenthBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt25ToInfGenTenthBin"]->Fill(CosThetaCalculation);
			    }	
			    else if(CosThetaGenerator > -0.8 && CosThetaGenerator <= -0.6){
			      histo1D["CosThetaDiffPt25ToInfGenNinthBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt25ToInfGenNinthBin"]->Fill(CosThetaCalculation);
			    }	
			    else if(CosThetaGenerator > -0.6 && CosThetaGenerator <= -0.4){
			      histo1D["CosThetaDiffPt25ToInfGenEighthBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt25ToInfGenEighthBin"]->Fill(CosThetaCalculation);
			    }	
			    else if(CosThetaGenerator > -0.4 && CosThetaGenerator <= -0.2){
			      histo1D["CosThetaDiffPt25ToInfGenSeventhBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt25ToInfGenSeventhBin"]->Fill(CosThetaCalculation);
			    }	
			    else if(CosThetaGenerator > -0.2 && CosThetaGenerator <= 0.){
			      histo1D["CosThetaDiffPt25ToInfGenSixthBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt25ToInfGenSixthBin"]->Fill(CosThetaCalculation);
			    }	
			    else if(CosThetaGenerator > 0. && CosThetaGenerator <= 0.2){
			      histo1D["CosThetaDiffPt25ToInfGenFifthBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt25ToInfGenFifthBin"]->Fill(CosThetaCalculation);
			    }	
			    else if(CosThetaGenerator > 0.2 && CosThetaGenerator <= 0.4){
			      histo1D["CosThetaDiffPt25ToInfGenFourthBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt25ToInfGenFourthBin"]->Fill(CosThetaCalculation);
			    }	
			    else if(CosThetaGenerator > 0.4 && CosThetaGenerator <= 0.6){
			      histo1D["CosThetaDiffPt25ToInfGenThirdBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt25ToInfGenThirdBin"]->Fill(CosThetaCalculation);
			    }	
			    else if(CosThetaGenerator > 0.6 && CosThetaGenerator <= 0.8){
			      histo1D["CosThetaDiffPt25ToInfGenSecondBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt25ToInfGenSecondBin"]->Fill(CosThetaCalculation);
			    }	
			    else if(CosThetaGenerator > 0.8 && CosThetaGenerator <= 1.){
			      histo1D["CosThetaDiffPt25ToInfGenFirstBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt25ToInfGenFirstBin"]->Fill(CosThetaCalculation);
			    }	

			  }
			  if(leptonKinFit[JetCombination].Pt()>=30){
			    histo1D["CosThetaDifferencePt30ToInf"]->Fill(CosThetaCalculation - CosThetaGenerator); 
			    histo1D["CosThetaRecoPt30ToInf"]->Fill(CosThetaCalculation);
			    histo1D["CosThetaRecoPt30ToInfForRatio"]->Fill(CosThetaCalculation);
			    histo1D["CosThetaGenPt30ToInf"]->Fill(CosThetaGenerator);
			    CosThetaDiffVSRecoPt30ToInf->Fill(CosThetaCalculation, (CosThetaCalculation - CosThetaGenerator) );
			    histo2D["CosThetaGenVsRecoPt30ToInf"]->Fill(CosThetaCalculation,CosThetaGenerator);

			    if(CosThetaGenerator <= -0.8){
			      TotalPt30ToInf[9]++;
			      histo1D["CosThetaDiffPt30ToInfGenTenthBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt30ToInfGenTenthBin"]->Fill(CosThetaCalculation);
			      if(CosThetaCalculation <= -0.8) ProbabilityPt30ToInf[9][9]++;
			      if(CosThetaCalculation > -0.8 && CosThetaCalculation <= -0.6) ProbabilityPt30ToInf[9][8]++;
			      if(CosThetaCalculation > -0.6 && CosThetaCalculation <= -0.4) ProbabilityPt30ToInf[9][7]++;
			      if(CosThetaCalculation > -0.4 && CosThetaCalculation <= -0.2) ProbabilityPt30ToInf[9][6]++;
			      if(CosThetaCalculation > -0.2 && CosThetaCalculation <= 0.) ProbabilityPt30ToInf[9][5]++;
			      if(CosThetaCalculation > 0. && CosThetaCalculation <= 0.2) ProbabilityPt30ToInf[9][4]++;
			      if(CosThetaCalculation > 0.2 && CosThetaCalculation <= 0.4) ProbabilityPt30ToInf[9][3]++;
			      if(CosThetaCalculation > 0.4 && CosThetaCalculation <= 0.6) ProbabilityPt30ToInf[9][2]++;
			      if(CosThetaCalculation > 0.6 && CosThetaCalculation <= 0.8) ProbabilityPt30ToInf[9][1]++;
			      if(CosThetaCalculation > 0.8 && CosThetaCalculation <= 1.) ProbabilityPt30ToInf[9][0]++;
			    }	
			    else if(CosThetaGenerator > -0.8 && CosThetaGenerator <= -0.6){
			      TotalPt30ToInf[8]++;
			      histo1D["CosThetaDiffPt30ToInfGenNinthBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt30ToInfGenNinthBin"]->Fill(CosThetaCalculation);
			      if(CosThetaCalculation <= -0.8) ProbabilityPt30ToInf[8][9]++;
			      if(CosThetaCalculation > -0.8 && CosThetaCalculation <= -0.6) ProbabilityPt30ToInf[8][8]++;
			      if(CosThetaCalculation > -0.6 && CosThetaCalculation <= -0.4) ProbabilityPt30ToInf[8][7]++;
			      if(CosThetaCalculation > -0.4 && CosThetaCalculation <= -0.2) ProbabilityPt30ToInf[8][6]++;
			      if(CosThetaCalculation > -0.2 && CosThetaCalculation <= 0.) ProbabilityPt30ToInf[8][5]++;
			      if(CosThetaCalculation > 0. && CosThetaCalculation <= 0.2) ProbabilityPt30ToInf[8][4]++;
			      if(CosThetaCalculation > 0.2 && CosThetaCalculation <= 0.4) ProbabilityPt30ToInf[8][3]++;
			      if(CosThetaCalculation > 0.4 && CosThetaCalculation <= 0.6) ProbabilityPt30ToInf[8][2]++;
			      if(CosThetaCalculation > 0.6 && CosThetaCalculation <= 0.8) ProbabilityPt30ToInf[8][1]++;
			      if(CosThetaCalculation > 0.8 && CosThetaCalculation <= 1.) ProbabilityPt30ToInf[8][0]++;
			    }	
			    else if(CosThetaGenerator > -0.6 && CosThetaGenerator <= -0.4){
			      TotalPt30ToInf[7]++;
			      histo1D["CosThetaDiffPt30ToInfGenEighthBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt30ToInfGenEighthBin"]->Fill(CosThetaCalculation);
			      if(CosThetaCalculation <= -0.8) ProbabilityPt30ToInf[7][9]++;
			      if(CosThetaCalculation > -0.8 && CosThetaCalculation <= -0.6) ProbabilityPt30ToInf[7][8]++;
			      if(CosThetaCalculation > -0.6 && CosThetaCalculation <= -0.4) ProbabilityPt30ToInf[7][7]++;
			      if(CosThetaCalculation > -0.4 && CosThetaCalculation <= -0.2) ProbabilityPt30ToInf[7][6]++;
			      if(CosThetaCalculation > -0.2 && CosThetaCalculation <= 0.) ProbabilityPt30ToInf[7][5]++;
			      if(CosThetaCalculation > 0. && CosThetaCalculation <= 0.2) ProbabilityPt30ToInf[7][4]++;
			      if(CosThetaCalculation > 0.2 && CosThetaCalculation <= 0.4) ProbabilityPt30ToInf[7][3]++;
			      if(CosThetaCalculation > 0.4 && CosThetaCalculation <= 0.6) ProbabilityPt30ToInf[7][2]++;
			      if(CosThetaCalculation > 0.6 && CosThetaCalculation <= 0.8) ProbabilityPt30ToInf[7][1]++;
			      if(CosThetaCalculation > 0.8 && CosThetaCalculation <= 1.) ProbabilityPt30ToInf[7][0]++;
			    }	
			    else if(CosThetaGenerator > -0.4 && CosThetaGenerator <= -0.2){
			      TotalPt30ToInf[6]++;
			      histo1D["CosThetaDiffPt30ToInfGenSeventhBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt30ToInfGenSeventhBin"]->Fill(CosThetaCalculation);
			      if(CosThetaCalculation <= -0.8) ProbabilityPt30ToInf[6][9]++;
			      if(CosThetaCalculation > -0.8 && CosThetaCalculation <= -0.6) ProbabilityPt30ToInf[6][8]++;
			      if(CosThetaCalculation > -0.6 && CosThetaCalculation <= -0.4) ProbabilityPt30ToInf[6][7]++;
			      if(CosThetaCalculation > -0.4 && CosThetaCalculation <= -0.2) ProbabilityPt30ToInf[6][6]++;
			      if(CosThetaCalculation > -0.2 && CosThetaCalculation <= 0.) ProbabilityPt30ToInf[6][5]++;
			      if(CosThetaCalculation > 0. && CosThetaCalculation <= 0.2) ProbabilityPt30ToInf[6][4]++;
			      if(CosThetaCalculation > 0.2 && CosThetaCalculation <= 0.4) ProbabilityPt30ToInf[6][3]++;
			      if(CosThetaCalculation > 0.4 && CosThetaCalculation <= 0.6) ProbabilityPt30ToInf[6][2]++;
			      if(CosThetaCalculation > 0.6 && CosThetaCalculation <= 0.8) ProbabilityPt30ToInf[6][1]++;
			      if(CosThetaCalculation > 0.8 && CosThetaCalculation <= 1.) ProbabilityPt30ToInf[6][0]++;
			    }	
			    else if(CosThetaGenerator > -0.2 && CosThetaGenerator <= 0.){
			      TotalPt30ToInf[5]++;
			      histo1D["CosThetaDiffPt30ToInfGenSixthBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt30ToInfGenSixthBin"]->Fill(CosThetaCalculation);
			      if(CosThetaCalculation <= -0.8) ProbabilityPt30ToInf[5][9]++;
			      if(CosThetaCalculation > -0.8 && CosThetaCalculation <= -0.6) ProbabilityPt30ToInf[5][8]++;
			      if(CosThetaCalculation > -0.6 && CosThetaCalculation <= -0.4) ProbabilityPt30ToInf[5][7]++;
			      if(CosThetaCalculation > -0.4 && CosThetaCalculation <= -0.2) ProbabilityPt30ToInf[5][6]++;
			      if(CosThetaCalculation > -0.2 && CosThetaCalculation <= 0.) ProbabilityPt30ToInf[5][5]++;
			      if(CosThetaCalculation > 0. && CosThetaCalculation <= 0.2) ProbabilityPt30ToInf[5][4]++;
			      if(CosThetaCalculation > 0.2 && CosThetaCalculation <= 0.4) ProbabilityPt30ToInf[5][3]++;
			      if(CosThetaCalculation > 0.4 && CosThetaCalculation <= 0.6) ProbabilityPt30ToInf[5][2]++;
			      if(CosThetaCalculation > 0.6 && CosThetaCalculation <= 0.8) ProbabilityPt30ToInf[5][1]++;
			      if(CosThetaCalculation > 0.8 && CosThetaCalculation <= 1.) ProbabilityPt30ToInf[5][0]++;
			    }	
			    else if(CosThetaGenerator > 0. && CosThetaGenerator <= 0.2){
			      TotalPt30ToInf[4]++;
			      histo1D["CosThetaDiffPt30ToInfGenFifthBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt30ToInfGenFifthBin"]->Fill(CosThetaCalculation);
			      if(CosThetaCalculation <= -0.8) ProbabilityPt30ToInf[4][9]++;
			      if(CosThetaCalculation > -0.8 && CosThetaCalculation <= -0.6) ProbabilityPt30ToInf[4][8]++;
			      if(CosThetaCalculation > -0.6 && CosThetaCalculation <= -0.4) ProbabilityPt30ToInf[4][7]++;
			      if(CosThetaCalculation > -0.4 && CosThetaCalculation <= -0.2) ProbabilityPt30ToInf[4][6]++;
			      if(CosThetaCalculation > -0.2 && CosThetaCalculation <= 0.) ProbabilityPt30ToInf[4][5]++;
			      if(CosThetaCalculation > 0. && CosThetaCalculation <= 0.2) ProbabilityPt30ToInf[4][4]++;
			      if(CosThetaCalculation > 0.2 && CosThetaCalculation <= 0.4) ProbabilityPt30ToInf[4][3]++;
			      if(CosThetaCalculation > 0.4 && CosThetaCalculation <= 0.6) ProbabilityPt30ToInf[4][2]++;
			      if(CosThetaCalculation > 0.6 && CosThetaCalculation <= 0.8) ProbabilityPt30ToInf[4][1]++;
			      if(CosThetaCalculation > 0.8 && CosThetaCalculation <= 1.) ProbabilityPt30ToInf[4][0]++;
			    }	
			    else if(CosThetaGenerator > 0.2 && CosThetaGenerator <= 0.4){
			      TotalPt30ToInf[3]++;
			      histo1D["CosThetaDiffPt30ToInfGenFourthBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt30ToInfGenFourthBin"]->Fill(CosThetaCalculation);
			      if(CosThetaCalculation <= -0.8) ProbabilityPt30ToInf[3][9]++;
			      if(CosThetaCalculation > -0.8 && CosThetaCalculation <= -0.6) ProbabilityPt30ToInf[3][8]++;
			      if(CosThetaCalculation > -0.6 && CosThetaCalculation <= -0.4) ProbabilityPt30ToInf[3][7]++;
			      if(CosThetaCalculation > -0.4 && CosThetaCalculation <= -0.2) ProbabilityPt30ToInf[3][6]++;
			      if(CosThetaCalculation > -0.2 && CosThetaCalculation <= 0.) ProbabilityPt30ToInf[3][5]++;
			      if(CosThetaCalculation > 0. && CosThetaCalculation <= 0.2) ProbabilityPt30ToInf[3][4]++;
			      if(CosThetaCalculation > 0.2 && CosThetaCalculation <= 0.4) ProbabilityPt30ToInf[3][3]++;
			      if(CosThetaCalculation > 0.4 && CosThetaCalculation <= 0.6) ProbabilityPt30ToInf[3][2]++;
			      if(CosThetaCalculation > 0.6 && CosThetaCalculation <= 0.8) ProbabilityPt30ToInf[3][1]++;
			      if(CosThetaCalculation > 0.8 && CosThetaCalculation <= 1.) ProbabilityPt30ToInf[3][0]++;
			    }	
			    else if(CosThetaGenerator > 0.4 && CosThetaGenerator <= 0.6){
			      TotalPt30ToInf[2]++;
			      histo1D["CosThetaDiffPt30ToInfGenThirdBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt30ToInfGenThirdBin"]->Fill(CosThetaCalculation);
			      if(CosThetaCalculation <= -0.8) ProbabilityPt30ToInf[2][9]++;
			      if(CosThetaCalculation > -0.8 && CosThetaCalculation <= -0.6) ProbabilityPt30ToInf[2][8]++;
			      if(CosThetaCalculation > -0.6 && CosThetaCalculation <= -0.4) ProbabilityPt30ToInf[2][7]++;
			      if(CosThetaCalculation > -0.4 && CosThetaCalculation <= -0.2) ProbabilityPt30ToInf[2][6]++;
			      if(CosThetaCalculation > -0.2 && CosThetaCalculation <= 0.) ProbabilityPt30ToInf[2][5]++;
			      if(CosThetaCalculation > 0. && CosThetaCalculation <= 0.2) ProbabilityPt30ToInf[2][4]++;
			      if(CosThetaCalculation > 0.2 && CosThetaCalculation <= 0.4) ProbabilityPt30ToInf[2][3]++;
			      if(CosThetaCalculation > 0.4 && CosThetaCalculation <= 0.6) ProbabilityPt30ToInf[2][2]++;
			      if(CosThetaCalculation > 0.6 && CosThetaCalculation <= 0.8) ProbabilityPt30ToInf[2][1]++;
			      if(CosThetaCalculation > 0.8 && CosThetaCalculation <= 1.) ProbabilityPt30ToInf[2][0]++;
			    }	
			    else if(CosThetaGenerator > 0.6 && CosThetaGenerator <= 0.8){
			      TotalPt30ToInf[1]++;
			      histo1D["CosThetaDiffPt30ToInfGenSecondBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt30ToInfGenSecondBin"]->Fill(CosThetaCalculation);
			      if(CosThetaCalculation <= -0.8) ProbabilityPt30ToInf[1][9]++;
			      if(CosThetaCalculation > -0.8 && CosThetaCalculation <= -0.6) ProbabilityPt30ToInf[1][8]++;
			      if(CosThetaCalculation > -0.6 && CosThetaCalculation <= -0.4) ProbabilityPt30ToInf[1][7]++;
			      if(CosThetaCalculation > -0.4 && CosThetaCalculation <= -0.2) ProbabilityPt30ToInf[1][6]++;
			      if(CosThetaCalculation > -0.2 && CosThetaCalculation <= 0.) ProbabilityPt30ToInf[1][5]++;
			      if(CosThetaCalculation > 0. && CosThetaCalculation <= 0.2) ProbabilityPt30ToInf[1][4]++;
			      if(CosThetaCalculation > 0.2 && CosThetaCalculation <= 0.4) ProbabilityPt30ToInf[1][3]++;
			      if(CosThetaCalculation > 0.4 && CosThetaCalculation <= 0.6) ProbabilityPt30ToInf[1][2]++;
			      if(CosThetaCalculation > 0.6 && CosThetaCalculation <= 0.8) ProbabilityPt30ToInf[1][1]++;
			      if(CosThetaCalculation > 0.8 && CosThetaCalculation <= 1.) ProbabilityPt30ToInf[1][0]++;
			    }	
			    else if(CosThetaGenerator > 0.8 && CosThetaGenerator <= 1.){
			      TotalPt30ToInf[0]++;
			      histo1D["CosThetaDiffPt30ToInfGenFirstBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt30ToInfGenFirstBin"]->Fill(CosThetaCalculation);
			      if(CosThetaCalculation <= -0.8) ProbabilityPt30ToInf[0][9]++;
			      if(CosThetaCalculation > -0.8 && CosThetaCalculation <= -0.6) ProbabilityPt30ToInf[0][8]++;
			      if(CosThetaCalculation > -0.6 && CosThetaCalculation <= -0.4) ProbabilityPt30ToInf[0][7]++;
			      if(CosThetaCalculation > -0.4 && CosThetaCalculation <= -0.2) ProbabilityPt30ToInf[0][6]++;
			      if(CosThetaCalculation > -0.2 && CosThetaCalculation <= 0.) ProbabilityPt30ToInf[0][5]++;
			      if(CosThetaCalculation > 0. && CosThetaCalculation <= 0.2) ProbabilityPt30ToInf[0][4]++;
			      if(CosThetaCalculation > 0.2 && CosThetaCalculation <= 0.4) ProbabilityPt30ToInf[0][3]++;
			      if(CosThetaCalculation > 0.4 && CosThetaCalculation <= 0.6) ProbabilityPt30ToInf[0][2]++;
			      if(CosThetaCalculation > 0.6 && CosThetaCalculation <= 0.8) ProbabilityPt30ToInf[0][1]++;
			      if(CosThetaCalculation > 0.8 && CosThetaCalculation <= 1.) ProbabilityPt30ToInf[0][0]++;
			    }	

			  }
			  if(leptonKinFit[JetCombination].Pt()>=35){
			    histo1D["CosThetaDifferencePt35ToInf"]->Fill(CosThetaCalculation - CosThetaGenerator); 
			    histo1D["CosThetaRecoPt35ToInf"]->Fill(CosThetaCalculation);
			    histo1D["CosThetaRecoPt35ToInfForRatio"]->Fill(CosThetaCalculation);
			    histo1D["CosThetaGenPt35ToInf"]->Fill(CosThetaGenerator);
			    CosThetaDiffVSRecoPt35ToInf->Fill(CosThetaCalculation, (CosThetaCalculation - CosThetaGenerator) );
			    histo2D["CosThetaGenVsRecoPt35ToInf"]->Fill(CosThetaCalculation,CosThetaGenerator);

			    if(CosThetaGenerator <= -0.8){
			      histo1D["CosThetaDiffPt35ToInfGenTenthBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt35ToInfGenTenthBin"]->Fill(CosThetaCalculation);
			    }	
			    else if(CosThetaGenerator > -0.8 && CosThetaGenerator <= -0.6){
			      histo1D["CosThetaDiffPt35ToInfGenNinthBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt35ToInfGenNinthBin"]->Fill(CosThetaCalculation);
			    }	
			    else if(CosThetaGenerator > -0.6 && CosThetaGenerator <= -0.4){
			      histo1D["CosThetaDiffPt35ToInfGenEighthBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt35ToInfGenEighthBin"]->Fill(CosThetaCalculation);
			    }	
			    else if(CosThetaGenerator > -0.4 && CosThetaGenerator <= -0.2){
			      histo1D["CosThetaDiffPt35ToInfGenSeventhBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt35ToInfGenSeventhBin"]->Fill(CosThetaCalculation);
			    }	
			    else if(CosThetaGenerator > -0.2 && CosThetaGenerator <= 0.){
			      histo1D["CosThetaDiffPt35ToInfGenSixthBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt35ToInfGenSixthBin"]->Fill(CosThetaCalculation);
			    }	
			    else if(CosThetaGenerator > 0. && CosThetaGenerator <= 0.2){
			      histo1D["CosThetaDiffPt35ToInfGenFifthBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt35ToInfGenFifthBin"]->Fill(CosThetaCalculation);
			    }	
			    else if(CosThetaGenerator > 0.2 && CosThetaGenerator <= 0.4){
			      histo1D["CosThetaDiffPt35ToInfGenFourthBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt35ToInfGenFourthBin"]->Fill(CosThetaCalculation);
			    }	
			    else if(CosThetaGenerator > 0.4 && CosThetaGenerator <= 0.6){
			      histo1D["CosThetaDiffPt35ToInfGenThirdBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt35ToInfGenThirdBin"]->Fill(CosThetaCalculation);
			    }	
			    else if(CosThetaGenerator > 0.6 && CosThetaGenerator <= 0.8){
			      histo1D["CosThetaDiffPt35ToInfGenSecondBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt35ToInfGenSecondBin"]->Fill(CosThetaCalculation);
			    }	
			    else if(CosThetaGenerator > 0.8 && CosThetaGenerator <= 1.){
			      histo1D["CosThetaDiffPt35ToInfGenFirstBin"]->Fill(CosThetaCalculation - CosThetaGenerator);
			      histo1D["CosThetaRecoPt35ToInfGenFirstBin"]->Fill(CosThetaCalculation);
			    }	
			    
			  }

			  //Leptonic side of event correct/wrong reconstructed:
			  if(BLeptonicIndex == CorrectBLept){
			    NumberUsedCorrectEvents++;
			    histo1D["KinFitProbCorrectBLept"]->Fill(ProbabilityOfKinFit[JetCombination]);
			    histo1D["CosThetaCorrectBLept"]->Fill(CosThetaCalculation);
			  }
			  else if(BLeptonicIndex != CorrectBLept){
			    NumberUsedWrongEvents++;
			    histo1D["KinFitProbWrongBLept"]->Fill(ProbabilityOfKinFit[JetCombination]);
			    histo1D["CosThetaWrongBLept"]->Fill(CosThetaCalculation);
			  }
			}//end of signal dataset
			if(dataSetName.find("Data") ==0) NumberUsedDataEvents++;

			//Kinematic histograms (before and after KinFit):
			if(UseChangedKinematics == true){
			  MSPlot["BLeptPtAfterKinFit"]->Fill(leptBJetKinFit[JetCombination].Pt(),datasets[iDataSet],true,Luminosity*scaleFactor*lumiWeight3D);
			  MSPlot["MetPtAfterKinFit"]->Fill(neutrinoKinFit[JetCombination].Pt(),datasets[iDataSet],true,Luminosity*scaleFactor*lumiWeight3D);
			  MSPlot["MuonPtAfterKinFit"]->Fill(leptonKinFit[JetCombination].Pt(),datasets[iDataSet],true,Luminosity*scaleFactor*lumiWeight3D);
			  MSPlot["WLeptPtAfterKinFit"]->Fill((neutrinoKinFit[JetCombination]+leptonKinFit[JetCombination]).Pt(),datasets[iDataSet],true,Luminosity*scaleFactor*lumiWeight3D);
			  MSPlot["TopLeptPtAfterKinFit"]->Fill((neutrinoKinFit[JetCombination]+leptonKinFit[JetCombination]+leptBJetKinFit[JetCombination]).Pt(),datasets[iDataSet],true,Luminosity*scaleFactor*lumiWeight3D);
			}
			else{
			  MSPlot["BLeptPtBeforeKinFit"]->Fill((selectedJets[BLeptonicIndex]).Pt(),datasets[iDataSet],true,Luminosity*scaleFactor*lumiWeight3D);
			  MSPlot["MetPtBeforeKinFit"]->Fill(MET.Pt(),datasets[iDataSet],true,Luminosity*scaleFactor*lumiWeight3D);
			  MSPlot["MuonPtBeforeKinFit"]->Fill(lepton.Pt(),datasets[iDataSet],true,Luminosity*scaleFactor*lumiWeight3D);
			  MSPlot["WLeptPtBeforeKinFit"]->Fill((MET+lepton).Pt(),datasets[iDataSet],true,Luminosity*scaleFactor*lumiWeight3D);
			  MSPlot["TopLeptPtBeforeKinFit"]->Fill((selectedJets[BLeptonicIndex]+MET+lepton).Pt(),datasets[iDataSet],true,Luminosity*scaleFactor*lumiWeight3D);
			}
			MSPlot["CosThetaNobTag"]->Fill(CosThetaCalculation,datasets[iDataSet],true,Luminosity*scaleFactor*lumiWeight3D);
			MSPlot["KinFitProbNobTag"]->Fill(ProbabilityOfKinFit[JetCombination],datasets[iDataSet],true,Luminosity*scaleFactor*lumiWeight3D);
		      }	

		      
		      //Filling of CosThetaValues array for all bTag cases:
		      if(ProbabilityOfKinFit[JetCombination] >= KinFitCut){
		      
			NumberRemainingEvents[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1][iDataSet]++;
			CosThetaValues[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1].push_back(CosThetaCalculation);
			LumiWeightVector[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1].push_back(lumiWeight3D);
			if(((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ) && dataSetName.find("JES") != 0){
			  if(BLeptonicIndex == CorrectBLept){//Count the number of events with correctly reconstructed leptonic b-jet
			    NumberBLeptCorrectEvents[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1][iDataSet]++;
			  }
			  CosThGen[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1].push_back(CosThetaGenerator);
			  EventCorrectionWeight[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1].push_back(scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor());		      		      
			}//End of semiMu sample
		      }
		      
		      if(TCHEbTagLoop == 7){ //optimal bTag case
			MSPlot["CosThetaTCHEMLept"]->Fill(CosThetaCalculation,datasets[iDataSet],true,Luminosity*scaleFactor*lumiWeight3D);
			MSPlot["KinFitProbTCHEMLept"]->Fill(ProbabilityOfKinFit[JetCombination],datasets[iDataSet],true,Luminosity*scaleFactor*lumiWeight3D);
			MSPlot["nPVAfterFoundCosThetabTag"]->Fill(nPrimaryVertices,datasets[iDataSet], true, Luminosity*scaleFactor*lumiWeight3D);
		      }		
		    }
		  }	      
		  
		  TCHEbTagLoop++;
		  if(TCHEbTagLoop == 14){
		    TCHPbTagLoop=2;
		    SSVHEbTagLoop=2;
		    SSVHPbTagLoop=2;
		    CSVbTagLoop=2;
		  }
		}
		ConsideredBTagger=1;
	      
		int JetCombination = bTagJetSelection.HighestProbSelection(TCHPbTagLoop,ConsideredBTagger,ProbabilityOfKinFit,ProbSOverSB,btagTCHE,btagTCHP,btagSSVHE,btagSSVHP,btagCSV);
		  
		int BLeptonicIndex = BLeptIndex[JetCombination];
		int BHadronicIndex = BHadrIndex[JetCombination];

		if(JetCombination != 999){//&& JetCombination == CorrectKinFitIndex){
		  float CosThetaCalculation;
		  if(UseChangedKinematics == true) CosThetaCalculation = bTagCosThetaCalculation.Calculation(leptonKinFit[JetCombination],neutrinoKinFit[JetCombination],leptBJetKinFit[JetCombination]);
		  else if(UseChangedKinematics == false) CosThetaCalculation = bTagCosThetaCalculation.CalcOrigKins(BLeptonicIndex,BHadronicIndex,lepton,selectedJets,MassW, MassTop);
		  
		  if(CosThetaCalculation != 999 && ProbabilityOfKinFit[JetCombination] >= KinFitCut){

		    if((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ){
		      if(TCHPbTagLoop == 7){ histo1D["CosThetaResTCHPMLept"]->Fill(CosThetaCalculation - CosThetaGenerator);}
		      if(TCHPbTagLoop == 11){ histo1D["CosThetaResTCHPTLept"]->Fill(CosThetaCalculation - CosThetaGenerator);}
		    }

		    if(TCHPbTagLoop == 7){
		      MSPlot["CosThetaTCHPMLept"]->Fill(CosThetaCalculation,datasets[iDataSet],true,Luminosity*scaleFactor*lumiWeight3D);
		    }
	    
		    NumberRemainingEvents[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1][iDataSet]++;
		    CosThetaValues[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1].push_back(CosThetaCalculation);
		    LumiWeightVector[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1].push_back(lumiWeight3D);
		    if(((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ) && dataSetName.find("JES") != 0){
		      if(BLeptonicIndex == CorrectBLept){//Count the number of events with correctly reconstructed leptonic b-jet
			NumberBLeptCorrectEvents[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1][iDataSet]++;
		      }
		      CosThGen[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1].push_back(CosThetaGenerator);
		      EventCorrectionWeight[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1].push_back(scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor());		      		      
		    }//End of semiMu sample
		  }
		}
				
		TCHPbTagLoop++;
	      }
	      ConsideredBTagger=2;

	      int JetCombination = bTagJetSelection.HighestProbSelection(SSVHEbTagLoop,ConsideredBTagger,ProbabilityOfKinFit,ProbSOverSB,btagTCHE,btagTCHP,btagSSVHE,btagSSVHP,btagCSV);
		  
	      int BLeptonicIndex = BLeptIndex[JetCombination];
	      int BHadronicIndex = BHadrIndex[JetCombination];
		    
	      if(JetCombination != 999){//&& JetCombination == CorrectKinFitIndex){
		float CosThetaCalculation;
		if(UseChangedKinematics == true) CosThetaCalculation = bTagCosThetaCalculation.Calculation(leptonKinFit[JetCombination],neutrinoKinFit[JetCombination],leptBJetKinFit[JetCombination]);
		else if(UseChangedKinematics == false) CosThetaCalculation = bTagCosThetaCalculation.CalcOrigKins(BLeptonicIndex,BHadronicIndex,lepton,selectedJets,MassW, MassTop);
		
		if(CosThetaCalculation != 999 && ProbabilityOfKinFit[JetCombination] >= KinFitCut){

		  if((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ){
		    if(SSVHEbTagLoop == 7){ histo1D["CosThetaResSSVHEMLept"]->Fill(CosThetaCalculation - CosThetaGenerator);}
		  }

		  if(SSVHEbTagLoop == 7){
		    
		    //Obtain cos theta* histo's for alternative helcities to investigate ratio with StandardModel:
		    if(dataSetName.find("Data") != 0){	
		      for(int ii=0;ii<SizeArray;ii++){
			histo1D["MSSVHEbTagCosThetaRight"+HelicityNumbers[ii][1].str()+"Long"+HelicityNumbers[ii][0].str()+"Left"+HelicityNumbers[ii][2].str()]->Fill(CosThetaCalculation,(HelicityWeight[ii]*Luminosity*scaleFactor*lumiWeight3D*datasets[iDataSet]->NormFactor()));  
		      }
		    }

		    MSPlot["CosThetaSSVHEMLept"]->Fill(CosThetaCalculation,datasets[iDataSet],true,Luminosity*scaleFactor*lumiWeight3D);
		    MSPlot["KinFitProbSSVHEMLept"]->Fill(ProbabilityOfKinFit[JetCombination],datasets[iDataSet],true,Luminosity*scaleFactor*lumiWeight3D);

		    if(((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ) && dataSetName.find("JES") != 0){
		      histo1D["CosThetaResSSVHEMLept"]->Fill(CosThetaCalculation - CosThetaGenerator);
		      if(BLeptonicIndex == CorrectBLept){
			histo1D["KinFitProbCorrectBLeptSSVHEM"]->Fill(ProbabilityOfKinFit[JetCombination]);
			histo1D["CosThetaCorrectBLeptSSVHEM"]->Fill(CosThetaCalculation);
		      }
		      else if(BLeptonicIndex != CorrectBLept){
			histo1D["KinFitProbWrongBLeptSSVHEM"]->Fill(ProbabilityOfKinFit[JetCombination]);
			histo1D["CosThetaWrongBLeptSSVHEM"]->Fill(CosThetaCalculation);
		      }

		    }
		  }

		  NumberRemainingEvents[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1][iDataSet]++;
		  CosThetaValues[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1].push_back(CosThetaCalculation);
		  LumiWeightVector[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1].push_back(lumiWeight3D);
		  if(((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ) && dataSetName.find("JES") != 0){
		    if(BLeptonicIndex == CorrectBLept){//Count the number of events with correctly reconstructed leptonic b-jet
		      NumberBLeptCorrectEvents[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1][iDataSet]++;
		    }
		    CosThGen[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1].push_back(CosThetaGenerator);
		    EventCorrectionWeight[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1].push_back(scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor()); 
		  }//End of semiMu sample
		  
		}
	      }	  
	      SSVHEbTagLoop++;
	    }
	    ConsideredBTagger=3;

	    int JetCombination = bTagJetSelection.HighestProbSelection(SSVHPbTagLoop,ConsideredBTagger,ProbabilityOfKinFit,ProbSOverSB,btagTCHE,btagTCHP,btagSSVHE,btagSSVHP,btagCSV);
		  
	    int BLeptonicIndex = BLeptIndex[JetCombination];
	    int BHadronicIndex = BHadrIndex[JetCombination];
		    
	    if(JetCombination != 999){//&& JetCombination == CorrectKinFitIndex){		
	      float CosThetaCalculation;
	      if(UseChangedKinematics == true) CosThetaCalculation = bTagCosThetaCalculation.Calculation(leptonKinFit[JetCombination],neutrinoKinFit[JetCombination],leptBJetKinFit[JetCombination]);	      
	      else if(UseChangedKinematics == false) CosThetaCalculation = bTagCosThetaCalculation.CalcOrigKins(BLeptonicIndex,BHadronicIndex,lepton,selectedJets,MassW, MassTop);

	      if(CosThetaCalculation != 999 && ProbabilityOfKinFit[JetCombination] >= KinFitCut){

		if((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ){
		  if(SSVHPbTagLoop == 11){ histo1D["CosThetaResSSVHPTLept"]->Fill(CosThetaCalculation - CosThetaGenerator); }
		}

		if(SSVHPbTagLoop == 7){
		  MSPlot["CosThetaSSVHPMLept"]->Fill(CosThetaCalculation,datasets[iDataSet],true,Luminosity*scaleFactor*lumiWeight3D);
		  MSPlot["KinFitProbSSVHPMLept"]->Fill(ProbabilityOfKinFit[JetCombination],datasets[iDataSet],true,Luminosity*scaleFactor*lumiWeight3D);
		}
		    
		NumberRemainingEvents[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1][iDataSet]++;
		CosThetaValues[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1].push_back(CosThetaCalculation);
		LumiWeightVector[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1].push_back(lumiWeight3D);
		if(((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ) && dataSetName.find("JES") != 0){
		  if(BLeptonicIndex == CorrectBLept){//Count the number of events with correctly reconstructed leptonic b-jet
		    NumberBLeptCorrectEvents[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1][iDataSet]++;
		  }
		  CosThGen[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1].push_back(CosThetaGenerator);
		  EventCorrectionWeight[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1].push_back(scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor());    
		}//End of semiMu sample	      
	      }	
	    }	
	    SSVHPbTagLoop++;
	  }
	  ConsideredBTagger=4;

	  int JetCombination = bTagJetSelection.HighestProbSelection(CSVbTagLoop,ConsideredBTagger,ProbabilityOfKinFit,ProbSOverSB,btagTCHE,btagTCHP,btagSSVHE,btagSSVHP,btagCSV);
		  
	  int BLeptonicIndex = BLeptIndex[JetCombination];
	  int BHadronicIndex = BHadrIndex[JetCombination];
		    
	  if(JetCombination != 999 ){//&& JetCombination == CorrectKinFitIndex){
	    float CosThetaCalculation;
	    if(UseChangedKinematics == true) CosThetaCalculation = bTagCosThetaCalculation.Calculation(leptonKinFit[JetCombination],neutrinoKinFit[JetCombination],leptBJetKinFit[JetCombination]);
	    else if(UseChangedKinematics == false) CosThetaCalculation = bTagCosThetaCalculation.CalcOrigKins(BLeptonicIndex,BHadronicIndex,lepton,selectedJets,MassW, MassTop);
	    
	    if(CosThetaCalculation != 999 && ProbabilityOfKinFit[JetCombination] >= KinFitCut){

	      if((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ){
		if(CSVbTagLoop == 3){ histo1D["CosThetaResCSVLLept"]->Fill(CosThetaCalculation - CosThetaGenerator);}
		if(CSVbTagLoop == 7){ histo1D["CosThetaResCSVMLept"]->Fill(CosThetaCalculation - CosThetaGenerator);}
		if(CSVbTagLoop == 11){ histo1D["CosThetaResCSVTLept"]->Fill(CosThetaCalculation - CosThetaGenerator);}
	      }
	      
	      if(CSVbTagLoop == 7){
		MSPlot["CosThetaCSVMLept"]->Fill(CosThetaCalculation,datasets[iDataSet],true,Luminosity*scaleFactor*lumiWeight3D);
	      }

	      if(CSVbTagLoop == 11){

		//Obtain cos theta* histo's for alternative helcities to investigate ratio with StandardModel:
		if(dataSetName.find("Data") != 0){	
		  for(int ii=0;ii<SizeArray;ii++){
		    histo1D["TCSVbTagCosThetaRight"+HelicityNumbers[ii][1].str()+"Long"+HelicityNumbers[ii][0].str()+"Left"+HelicityNumbers[ii][2].str()]->Fill(CosThetaCalculation,(HelicityWeight[ii]*Luminosity*scaleFactor*lumiWeight3D*datasets[iDataSet]->NormFactor()));  
		  }
		}
		
	      }
		    
	      NumberRemainingEvents[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1][iDataSet]++;
	      CosThetaValues[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1].push_back(CosThetaCalculation);
	      LumiWeightVector[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1].push_back(lumiWeight3D);
	      if(((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ) && dataSetName.find("JES") != 0){
		if(BLeptonicIndex == CorrectBLept){//Count the number of events with correctly reconstructed leptonic b-jet
		  NumberBLeptCorrectEvents[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1][iDataSet]++;
		}
		CosThGen[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1].push_back(CosThetaGenerator);
		EventCorrectionWeight[TCHEbTagLoop-1+TCHPbTagLoop-1+SSVHEbTagLoop-1+SSVHPbTagLoop-1+CSVbTagLoop-1].push_back(scaleFactor*Luminosity*lumiWeight3D*dataSet->NormFactor());     
	      }//End of semiMu sample	    
	    }	
	  }	
	  CSVbTagLoop++;
	}
		
      }// End of loop over selected events
      
    } // end loop over events in wTrees    
    
    std::cout << "  " << endl;
    std::cout << " size of cos theta : " << CosThetaValues[0].size() << endl;
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

    cout << " " << endl;
    cout << " oooooooooooOOOOOOOOOOOOoooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooo " << endl;
    cout << "                                               Probability matrix Reco vs Gen (Pt 10 To Inf)" << endl;
    cout << "                                              ----------------------------------- " << endl;
    cout << "    Gen " << endl;
    cout << "    1    " << endl;
    cout << "         " <<  ProbabilityPt10ToInf[0][9]/TotalPt10ToInf[0] << "    " << ProbabilityPt10ToInf[0][8]/TotalPt10ToInf[0] << "    " << ProbabilityPt10ToInf[0][7]/TotalPt10ToInf[0] << "    " << ProbabilityPt10ToInf[0][6]/TotalPt10ToInf[0] << "    " << ProbabilityPt10ToInf[0][5]/TotalPt10ToInf[0] << "    " << ProbabilityPt10ToInf[0][4]/TotalPt10ToInf[0] << "    " << ProbabilityPt10ToInf[0][3]/TotalPt10ToInf[0] << "    " << ProbabilityPt10ToInf[0][2]/TotalPt10ToInf[0] << "    " << ProbabilityPt10ToInf[0][1]/TotalPt10ToInf[0] << "    " << ProbabilityPt10ToInf[0][0]/TotalPt10ToInf[0] << endl;
    cout << "   0.8   " << endl;
    cout << "         " <<  ProbabilityPt10ToInf[1][9]/TotalPt10ToInf[1] << "    " << ProbabilityPt10ToInf[1][8]/TotalPt10ToInf[1] << "    " << ProbabilityPt10ToInf[1][7]/TotalPt10ToInf[1] << "    " << ProbabilityPt10ToInf[1][6]/TotalPt10ToInf[1] << "    " << ProbabilityPt10ToInf[1][5]/TotalPt10ToInf[1] << "    " << ProbabilityPt10ToInf[1][4]/TotalPt10ToInf[1] << "    " << ProbabilityPt10ToInf[1][3]/TotalPt10ToInf[1] << "    " << ProbabilityPt10ToInf[1][2]/TotalPt10ToInf[1] << "  " << ProbabilityPt10ToInf[1][1]/TotalPt10ToInf[1] << "    " << ProbabilityPt10ToInf[1][0]/TotalPt10ToInf[1] << endl;
    cout << "   0.6   " << endl;
    cout << "         " <<  ProbabilityPt10ToInf[2][9]/TotalPt10ToInf[2] << "    " << ProbabilityPt10ToInf[2][8]/TotalPt10ToInf[2] << "    " << ProbabilityPt10ToInf[2][7]/TotalPt10ToInf[2] << "    " << ProbabilityPt10ToInf[2][6]/TotalPt10ToInf[2] << "    " << ProbabilityPt10ToInf[2][5]/TotalPt10ToInf[2] << "    " << ProbabilityPt10ToInf[2][4]/TotalPt10ToInf[2] << "    " << ProbabilityPt10ToInf[2][3]/TotalPt10ToInf[2] << "    " << ProbabilityPt10ToInf[2][2]/TotalPt10ToInf[2] << "    " << ProbabilityPt10ToInf[2][1]/TotalPt10ToInf[2] << "    " << ProbabilityPt10ToInf[2][0]/TotalPt10ToInf[2] << endl;
    cout << "   0.4   " << endl;
    cout << "         " <<  ProbabilityPt10ToInf[3][9]/TotalPt10ToInf[3] << "    " << ProbabilityPt10ToInf[3][8]/TotalPt10ToInf[3] << "    " << ProbabilityPt10ToInf[3][7]/TotalPt10ToInf[3] << "    " << ProbabilityPt10ToInf[3][6]/TotalPt10ToInf[3] << "    " << ProbabilityPt10ToInf[3][5]/TotalPt10ToInf[3] << "    " << ProbabilityPt10ToInf[3][4]/TotalPt10ToInf[3] << "    " << ProbabilityPt10ToInf[3][3]/TotalPt10ToInf[3] << "    " << ProbabilityPt10ToInf[3][2]/TotalPt10ToInf[3] << "    " << ProbabilityPt10ToInf[3][1]/TotalPt10ToInf[3] << "    " << ProbabilityPt10ToInf[3][0]/TotalPt10ToInf[3] << endl;
    cout << "   0.2   " << endl;
    cout << "         " <<  ProbabilityPt10ToInf[4][9]/TotalPt10ToInf[4] << "    " << ProbabilityPt10ToInf[4][8]/TotalPt10ToInf[4] << "    " << ProbabilityPt10ToInf[4][7]/TotalPt10ToInf[4] << "    " << ProbabilityPt10ToInf[4][6]/TotalPt10ToInf[4] << "    " << ProbabilityPt10ToInf[4][5]/TotalPt10ToInf[4] << "    " << ProbabilityPt10ToInf[4][4]/TotalPt10ToInf[4] << "    " << ProbabilityPt10ToInf[4][3]/TotalPt10ToInf[4] << "    " << ProbabilityPt10ToInf[4][2]/TotalPt10ToInf[4] << "    " << ProbabilityPt10ToInf[4][1]/TotalPt10ToInf[4] << "    " << ProbabilityPt10ToInf[4][0]/TotalPt10ToInf[4] << endl;
    cout << "   0.   " << endl;
    cout << "         " <<  ProbabilityPt10ToInf[5][9]/TotalPt10ToInf[5] << "    " << ProbabilityPt10ToInf[5][8]/TotalPt10ToInf[5] << "    " << ProbabilityPt10ToInf[5][7]/TotalPt10ToInf[5] << "    " << ProbabilityPt10ToInf[5][6]/TotalPt10ToInf[5] << "    " << ProbabilityPt10ToInf[5][5]/TotalPt10ToInf[5] << "    " << ProbabilityPt10ToInf[5][4]/TotalPt10ToInf[5] << "    " << ProbabilityPt10ToInf[5][3]/TotalPt10ToInf[5] << "    " << ProbabilityPt10ToInf[5][2]/TotalPt10ToInf[5] << "    " << ProbabilityPt10ToInf[5][1]/TotalPt10ToInf[5] << "    " << ProbabilityPt10ToInf[5][0]/TotalPt10ToInf[5] << endl;
    cout << "   -0.2   " << endl;
    cout << "         " <<  ProbabilityPt10ToInf[6][9]/TotalPt10ToInf[6] << "    " << ProbabilityPt10ToInf[6][8]/TotalPt10ToInf[6] << "    " << ProbabilityPt10ToInf[6][7]/TotalPt10ToInf[6] << "    " << ProbabilityPt10ToInf[6][6]/TotalPt10ToInf[6] << "    " << ProbabilityPt10ToInf[6][5]/TotalPt10ToInf[6] << "    " << ProbabilityPt10ToInf[6][4]/TotalPt10ToInf[6] << "    " << ProbabilityPt10ToInf[6][3]/TotalPt10ToInf[6] << "    " << ProbabilityPt10ToInf[6][2]/TotalPt10ToInf[6] << "    " << ProbabilityPt10ToInf[6][1]/TotalPt10ToInf[6] << "    " << ProbabilityPt10ToInf[6][0]/TotalPt10ToInf[6] << endl;
    cout << "   -0.4   " << endl;
    cout << "         " <<  ProbabilityPt10ToInf[7][9]/TotalPt10ToInf[7] << "    " << ProbabilityPt10ToInf[7][8]/TotalPt10ToInf[7] << "    " << ProbabilityPt10ToInf[7][7]/TotalPt10ToInf[7] << "    " << ProbabilityPt10ToInf[7][6]/TotalPt10ToInf[7] << "    " << ProbabilityPt10ToInf[7][5]/TotalPt10ToInf[7] << "    " << ProbabilityPt10ToInf[7][4]/TotalPt10ToInf[7] << "    " << ProbabilityPt10ToInf[7][3]/TotalPt10ToInf[7] << "    " << ProbabilityPt10ToInf[7][2]/TotalPt10ToInf[7] << "    " << ProbabilityPt10ToInf[7][1]/TotalPt10ToInf[7] << "    " << ProbabilityPt10ToInf[7][0]/TotalPt10ToInf[7] << endl;
    cout << "   -0.6   " << endl;
    cout << "         " <<  ProbabilityPt10ToInf[8][9]/TotalPt10ToInf[8] << "    " << ProbabilityPt10ToInf[8][8]/TotalPt10ToInf[8] << "    " << ProbabilityPt10ToInf[8][7]/TotalPt10ToInf[8] << "    " << ProbabilityPt10ToInf[8][6]/TotalPt10ToInf[8] << "    " << ProbabilityPt10ToInf[8][5]/TotalPt10ToInf[8] << "    " << ProbabilityPt10ToInf[8][4]/TotalPt10ToInf[8] << "    " << ProbabilityPt10ToInf[8][3]/TotalPt10ToInf[8] << "    " << ProbabilityPt10ToInf[8][2]/TotalPt10ToInf[8] << "    " << ProbabilityPt10ToInf[8][1]/TotalPt10ToInf[8] << "    " << ProbabilityPt10ToInf[8][0]/TotalPt10ToInf[8] << endl;
    cout << "   -0.8   " << endl;
    cout << "         " <<  ProbabilityPt10ToInf[9][9]/TotalPt10ToInf[9] << "    " << ProbabilityPt10ToInf[9][8]/TotalPt10ToInf[9] << "    " << ProbabilityPt10ToInf[9][7]/TotalPt10ToInf[9] << "    " << ProbabilityPt10ToInf[9][6]/TotalPt10ToInf[9] << "    " << ProbabilityPt10ToInf[9][5]/TotalPt10ToInf[9] << "    " << ProbabilityPt10ToInf[9][4]/TotalPt10ToInf[9] << "    " << ProbabilityPt10ToInf[9][3]/TotalPt10ToInf[9] << "    " << ProbabilityPt10ToInf[9][2]/TotalPt10ToInf[9] << "    " << ProbabilityPt10ToInf[9][1]/TotalPt10ToInf[9] << "    " << ProbabilityPt10ToInf[9][0]/TotalPt10ToInf[9] << endl;
    cout << "   -1.  -1         -0.8         -0.6        -0.4         -0.2       0.        0.2           0.4          0.6           0.8          1      Reco" << endl;
    cout << " " << endl;
    cout << "                                               Probability matrix Reco vs Gen (Pt 30 To Inf)" << endl;
    cout << "                                              ----------------------------------- " << endl;
    cout << "    Gen " << endl;
    cout << "    1    " << endl;
    cout << "         " <<  ProbabilityPt30ToInf[0][9]/TotalPt30ToInf[0] << "    " << ProbabilityPt30ToInf[0][8]/TotalPt30ToInf[0] << "    " << ProbabilityPt30ToInf[0][7]/TotalPt30ToInf[0] << "    " << ProbabilityPt30ToInf[0][6]/TotalPt30ToInf[0] << "    " << ProbabilityPt30ToInf[0][5]/TotalPt30ToInf[0] << "    " << ProbabilityPt30ToInf[0][4]/TotalPt30ToInf[0] << "    " << ProbabilityPt30ToInf[0][3]/TotalPt30ToInf[0] << "    " << ProbabilityPt30ToInf[0][2]/TotalPt30ToInf[0] << "    " << ProbabilityPt30ToInf[0][1]/TotalPt30ToInf[0] << "    " << ProbabilityPt30ToInf[0][0]/TotalPt30ToInf[0] << endl;
    cout << "   0.8   " << endl;
    cout << "         " <<  ProbabilityPt30ToInf[1][9]/TotalPt30ToInf[1] << "    " << ProbabilityPt30ToInf[1][8]/TotalPt30ToInf[1] << "    " << ProbabilityPt30ToInf[1][7]/TotalPt30ToInf[1] << "    " << ProbabilityPt30ToInf[1][6]/TotalPt30ToInf[1] << "    " << ProbabilityPt30ToInf[1][5]/TotalPt30ToInf[1] << "    " << ProbabilityPt30ToInf[1][4]/TotalPt30ToInf[1] << "    " << ProbabilityPt30ToInf[1][3]/TotalPt30ToInf[1] << "    " << ProbabilityPt30ToInf[1][2]/TotalPt30ToInf[1] << "  " << ProbabilityPt30ToInf[1][1]/TotalPt30ToInf[1] << "    " << ProbabilityPt30ToInf[1][0]/TotalPt30ToInf[1] << endl;
    cout << "   0.6   " << endl;
    cout << "         " <<  ProbabilityPt30ToInf[2][9]/TotalPt30ToInf[2] << "    " << ProbabilityPt30ToInf[2][8]/TotalPt30ToInf[2] << "    " << ProbabilityPt30ToInf[2][7]/TotalPt30ToInf[2] << "    " << ProbabilityPt30ToInf[2][6]/TotalPt30ToInf[2] << "    " << ProbabilityPt30ToInf[2][5]/TotalPt30ToInf[2] << "    " << ProbabilityPt30ToInf[2][4]/TotalPt30ToInf[2] << "    " << ProbabilityPt30ToInf[2][3]/TotalPt30ToInf[2] << "    " << ProbabilityPt30ToInf[2][2]/TotalPt30ToInf[2] << "    " << ProbabilityPt30ToInf[2][1]/TotalPt30ToInf[2] << "    " << ProbabilityPt30ToInf[2][0]/TotalPt30ToInf[2] << endl;
    cout << "   0.4   " << endl;
    cout << "         " <<  ProbabilityPt30ToInf[3][9]/TotalPt30ToInf[3] << "    " << ProbabilityPt30ToInf[3][8]/TotalPt30ToInf[3] << "    " << ProbabilityPt30ToInf[3][7]/TotalPt30ToInf[3] << "    " << ProbabilityPt30ToInf[3][6]/TotalPt30ToInf[3] << "    " << ProbabilityPt30ToInf[3][5]/TotalPt30ToInf[3] << "    " << ProbabilityPt30ToInf[3][4]/TotalPt30ToInf[3] << "    " << ProbabilityPt30ToInf[3][3]/TotalPt30ToInf[3] << "    " << ProbabilityPt30ToInf[3][2]/TotalPt30ToInf[3] << "    " << ProbabilityPt30ToInf[3][1]/TotalPt30ToInf[3] << "    " << ProbabilityPt30ToInf[3][0]/TotalPt30ToInf[3] << endl;
    cout << "   0.2   " << endl;
    cout << "         " <<  ProbabilityPt30ToInf[4][9]/TotalPt30ToInf[4] << "    " << ProbabilityPt30ToInf[4][8]/TotalPt30ToInf[4] << "    " << ProbabilityPt30ToInf[4][7]/TotalPt30ToInf[4] << "    " << ProbabilityPt30ToInf[4][6]/TotalPt30ToInf[4] << "    " << ProbabilityPt30ToInf[4][5]/TotalPt30ToInf[4] << "    " << ProbabilityPt30ToInf[4][4]/TotalPt30ToInf[4] << "    " << ProbabilityPt30ToInf[4][3]/TotalPt30ToInf[4] << "    " << ProbabilityPt30ToInf[4][2]/TotalPt30ToInf[4] << "    " << ProbabilityPt30ToInf[4][1]/TotalPt30ToInf[4] << "    " << ProbabilityPt30ToInf[4][0]/TotalPt30ToInf[4] << endl;
    cout << "   0.   " << endl;
    cout << "         " <<  ProbabilityPt30ToInf[5][9]/TotalPt30ToInf[5] << "    " << ProbabilityPt30ToInf[5][8]/TotalPt30ToInf[5] << "    " << ProbabilityPt30ToInf[5][7]/TotalPt30ToInf[5] << "    " << ProbabilityPt30ToInf[5][6]/TotalPt30ToInf[5] << "    " << ProbabilityPt30ToInf[5][5]/TotalPt30ToInf[5] << "    " << ProbabilityPt30ToInf[5][4]/TotalPt30ToInf[5] << "    " << ProbabilityPt30ToInf[5][3]/TotalPt30ToInf[5] << "    " << ProbabilityPt30ToInf[5][2]/TotalPt30ToInf[5] << "    " << ProbabilityPt30ToInf[5][1]/TotalPt30ToInf[5] << "    " << ProbabilityPt30ToInf[5][0]/TotalPt30ToInf[5] << endl;
    cout << "   -0.2   " << endl;
    cout << "         " <<  ProbabilityPt30ToInf[6][9]/TotalPt30ToInf[6] << "    " << ProbabilityPt30ToInf[6][8]/TotalPt30ToInf[6] << "    " << ProbabilityPt30ToInf[6][7]/TotalPt30ToInf[6] << "    " << ProbabilityPt30ToInf[6][6]/TotalPt30ToInf[6] << "    " << ProbabilityPt30ToInf[6][5]/TotalPt30ToInf[6] << "    " << ProbabilityPt30ToInf[6][4]/TotalPt30ToInf[6] << "    " << ProbabilityPt30ToInf[6][3]/TotalPt30ToInf[6] << "    " << ProbabilityPt30ToInf[6][2]/TotalPt30ToInf[6] << "    " << ProbabilityPt30ToInf[6][1]/TotalPt30ToInf[6] << "    " << ProbabilityPt30ToInf[6][0]/TotalPt30ToInf[6] << endl;
    cout << "   -0.4   " << endl;
    cout << "         " <<  ProbabilityPt30ToInf[7][9]/TotalPt30ToInf[7] << "    " << ProbabilityPt30ToInf[7][8]/TotalPt30ToInf[7] << "    " << ProbabilityPt30ToInf[7][7]/TotalPt30ToInf[7] << "    " << ProbabilityPt30ToInf[7][6]/TotalPt30ToInf[7] << "    " << ProbabilityPt30ToInf[7][5]/TotalPt30ToInf[7] << "    " << ProbabilityPt30ToInf[7][4]/TotalPt30ToInf[7] << "    " << ProbabilityPt30ToInf[7][3]/TotalPt30ToInf[7] << "    " << ProbabilityPt30ToInf[7][2]/TotalPt30ToInf[7] << "    " << ProbabilityPt30ToInf[7][1]/TotalPt30ToInf[7] << "    " << ProbabilityPt30ToInf[7][0]/TotalPt30ToInf[7] << endl;
    cout << "   -0.6   " << endl;
    cout << "         " <<  ProbabilityPt30ToInf[8][9]/TotalPt30ToInf[8] << "    " << ProbabilityPt30ToInf[8][8]/TotalPt30ToInf[8] << "    " << ProbabilityPt30ToInf[8][7]/TotalPt30ToInf[8] << "    " << ProbabilityPt30ToInf[8][6]/TotalPt30ToInf[8] << "    " << ProbabilityPt30ToInf[8][5]/TotalPt30ToInf[8] << "    " << ProbabilityPt30ToInf[8][4]/TotalPt30ToInf[8] << "    " << ProbabilityPt30ToInf[8][3]/TotalPt30ToInf[8] << "    " << ProbabilityPt30ToInf[8][2]/TotalPt30ToInf[8] << "    " << ProbabilityPt30ToInf[8][1]/TotalPt30ToInf[8] << "    " << ProbabilityPt30ToInf[8][0]/TotalPt30ToInf[8] << endl;
    cout << "   -0.8   " << endl;
    cout << "         " <<  ProbabilityPt30ToInf[9][9]/TotalPt30ToInf[9] << "    " << ProbabilityPt30ToInf[9][8]/TotalPt30ToInf[9] << "    " << ProbabilityPt30ToInf[9][7]/TotalPt30ToInf[9] << "    " << ProbabilityPt30ToInf[9][6]/TotalPt30ToInf[9] << "    " << ProbabilityPt30ToInf[9][5]/TotalPt30ToInf[9] << "    " << ProbabilityPt30ToInf[9][4]/TotalPt30ToInf[9] << "    " << ProbabilityPt30ToInf[9][3]/TotalPt30ToInf[9] << "    " << ProbabilityPt30ToInf[9][2]/TotalPt30ToInf[9] << "    " << ProbabilityPt30ToInf[9][1]/TotalPt30ToInf[9] << "    " << ProbabilityPt30ToInf[9][0]/TotalPt30ToInf[9] << endl;
    cout << "   -1.  -1         -0.8         -0.6        -0.4         -0.2       0.        0.2           0.4          0.6           0.8          1      Reco" << endl;
 
    
    int TCHE=0;
    int TCHP=0;
    int SSVHE=0;
    int SSVHP=0;
    int CSV=0;
    
    int ndimen=3;
    
    while(CSV<13){
      while(SSVHP<13){ //Try to improve code by making a Minuit fitter class !!
	while(SSVHE<13){
	  while(TCHP<13){
	    while(TCHE<13){
	      
	      std::string CosThetaDataString = "CosThetaData_TCHE"+bTagValues[TCHE]+"_TCHP"+bTagValues[TCHP]+"_SSVHE"+bTagValues[SSVHE]+"_SSVHP"+bTagValues[SSVHP]+"_CSV"+bTagValues[CSV];
	      std::string CosThetaSignalString = "CosThetaSignal_TCHE"+bTagValues[TCHE]+"_TCHP"+bTagValues[TCHP]+"_SSVHE"+bTagValues[SSVHE]+"_SSVHP"+bTagValues[SSVHP]+"_CSV"+bTagValues[CSV];
	      std::string CosThetaBckgString = "CosThetaBckg_TCHE"+bTagValues[TCHE]+"_TCHP"+bTagValues[TCHP]+"_SSVHE"+bTagValues[SSVHE]+"_SSVHP"+bTagValues[SSVHP]+"_CSV"+bTagValues[CSV];
	      std::string CosThetaGeneratorString = "CosThetaGenerator_TCHE"+bTagValues[TCHE]+"_TCHP"+bTagValues[TCHP]+"_SSVHE"+bTagValues[SSVHE]+"_SSVHP"+bTagValues[SSVHP]+"_CSV"+bTagValues[CSV];
	      
	      if(iDataSet == 0){
		histo1D[CosThetaDataString]=new TH1F(CosThetaDataString.c_str(),CosThetaDataString.c_str(),CosThetaBinNumber,-1,1);
		histo1D[CosThetaDataString]->SetDirectory(0);
		histo1D[CosThetaSignalString]=new TH1F(CosThetaSignalString.c_str(),CosThetaSignalString.c_str(),CosThetaBinNumber,-1,1);
		histo1D[CosThetaSignalString]->SetDirectory(0);
		histo1D[CosThetaBckgString]=new TH1F(CosThetaBckgString.c_str(),CosThetaBckgString.c_str(),CosThetaBinNumber,-1,1);
		histo1D[CosThetaBckgString]->SetDirectory(0);
	      }
	      
	      if((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ){   //Filling of the genttbar histo
		char hisname[100];
		sprintf(hisname,"CosThetaGenerator_TCHE%s_TCHP%s_SSVHE%s_SSVHP%s_CSV%s", bTagValues[TCHE].c_str(),bTagValues[TCHP].c_str(),bTagValues[SSVHE].c_str(),bTagValues[SSVHP].c_str(),bTagValues[CSV].c_str());
		for (int ibinn=0; ibinn<CosThetaBinNumber; ibinn++){
		  genttbarhisto[ibinn]= new TNtuple(hisname,hisname,"costhgen:evtweight");
		  genttbarhisto[ibinn]->SetDirectory(0);
		}
		for(int ii=0; ii<CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV].size();ii++){		
		  for(int iBin=0; iBin< CosThetaBinNumber; iBin++){
		    
		    if(CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV][ii] >= binEdge[iBin] && CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV][ii] < binEdge[iBin+1]){
		      double costhgood = CosThGen[TCHE+TCHP+SSVHE+SSVHP+CSV][ii]; 
		      double thisWeight = EventCorrectionWeight[TCHE+TCHP+SSVHE+SSVHP+CSV][ii];
		      genttbarhisto[iBin]->Fill(costhgood , thisWeight) ; 
		    }
		    else if(CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV][ii] ==1){ //1 is included in last bin
		      genttbarhisto[CosThetaBinNumber-1]->Fill(CosThGen[TCHE+TCHP+SSVHE+SSVHP+CSV][ii], EventCorrectionWeight[TCHE+TCHP+SSVHE+SSVHP+CSV][ii]);
		    }	      	   
		  }
		}
	      }
	      
	      float scaleFactor = 1.;  //Make an array of the scaleFactor such that it can be different for every event !!
	      
	      for(int ii=0;ii<CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV].size();ii++){
		//Change data to systematics since then nominal values will be reweighted!!
		if(WSystResults == true && dataSetName.find("WJets_Nominal") != 0)  //WJets syst
		  histo1D[CosThetaDataString]->Fill(CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV][ii],scaleFactor*(LumiWeightVector[TCHE+TCHP+SSVHE+SSVHP+CSV][ii])*Luminosity*dataSet->NormFactor());    				
		if(WSystResults == false && (dataSetName.find("Data") == 0 || dataSetName.find("JES") == 0)) //Nominal and JES systematics
		  histo1D[CosThetaDataString]->Fill(CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV][ii],scaleFactor*(LumiWeightVector[TCHE+TCHP+SSVHE+SSVHP+CSV][ii])*Luminosity*dataSet->NormFactor());    				
		if(SignalOnly == true && (dataSetName.find("TTbarJets_SemiMu") == 0)){
		  histo1D[CosThetaDataString]->Fill(CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV][ii],scaleFactor*(LumiWeightVector[TCHE+TCHP+SSVHE+SSVHP+CSV][ii])*Luminosity*dataSet->NormFactor());    				
		}		
		if(((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ) && dataSetName.find("JES_") != 0)
		  histo1D[CosThetaSignalString]->Fill(CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV][ii],Luminosity*scaleFactor*(LumiWeightVector[TCHE+TCHP+SSVHE+SSVHP+CSV][ii])*NominalNormFactor);
		else if(dataSetName.find("Data") !=0 && dataSetName.find("JES") != 0 && dataSetName.find("Syst") != 0) 
		  histo1D[CosThetaBckgString]->Fill(CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV][ii],Luminosity*scaleFactor*(LumiWeightVector[TCHE+TCHP+SSVHE+SSVHP+CSV][ii])*NominalNormFactor);
	      }
	      
	      if(iDataSet==(datasets.size()-1)){//Go in this loop when the last datasample is active

		MinuitFitter minuitFitter = MinuitFitter(histo1D[CosThetaDataString], histo1D[CosThetaSignalString], histo1D[CosThetaBckgString], 0.6671, 0.3325, 0.0004,genttbarhisto,ndimen);
		
		MyChi2 myChi2 (histo1D[CosThetaDataString], histo1D[CosThetaSignalString], histo1D[CosThetaBckgString], 0.6671, 0.3325, 0.0004); 		

		TMinuitMinimizer minimizer;
		minimizer.SetFunction(myChi2);
		minimizer.SetErrorDef(1.0);   // 1.0 for chi2, 0.5 for -logL
		minimizer.SetVariable(0, "F0", 0.6671, 1.e-4);
		minimizer.SetVariable(1, "FL", 0.3325, 1.e-4);
		if (ndimen==3)  minimizer.SetVariable(2, "Normal", 1., 1.e-4); 

		bool isValid = minimizer.Minimize();
		double f0result, frresult,flresult;
		double ef0result, efrresult,eflresult;
		if (isValid) {
		  //minimizer.PrintResults();
		  f0result =minimizer.X()[0]; ef0result= minimizer.Errors()[0];
		  flresult= minimizer.X()[1]; eflresult= minimizer.Errors()[1];
		  
		  PresentationTex << PresentationOutput[TCHE+TCHP+SSVHE+SSVHP+CSV] << " & ";
		  for(int ii=0; ii< nameDataSet.size(); ii++){
		    //if(nameDataSet[ii].find("JES") != 0 && nameDataSet[ii].find("Syst") != 0 && ii < nameDataSet.size()-1){ //Only output for samples with no systematics
		    if(ii < nameDataSet.size()-1){ 
		      PresentationTex << NumberRemainingEvents[TCHE+TCHP+SSVHE+SSVHP+CSV][ii] << " & ";
		    }
		    else if(ii == nameDataSet.size()-1 ){ 
		      PresentationTex << NumberRemainingEvents[TCHE+TCHP+SSVHE+SSVHP+CSV][ii] << " & "; //For presentation this is not the end of the table, helicity values still need to be included !!
		      PresentationTex << NumberBLeptCorrectEvents[TCHE+TCHP+SSVHE+SSVHP+CSV][ii] << " & ";
		    }
		  }		  
		  if (ndimen==3){
		    frresult = 1.-minimizer.X()[0]-minimizer.X()[1];
		    double er0=minimizer.Errors()[0];
		    double er1=minimizer.Errors()[1];
		    double cov01 = minimizer.CovMatrix(0,1);		  
		    efrresult = sqrt( er0*er0 + er1*er1 + 2.*cov01);
		    PresentationTex << frresult << " & " << frresult-SMfrResult <<" & " << efrresult << " & " << flresult << " & " << flresult-SMflResult << " & " << eflresult << " & " << f0result << " & " << f0result-SMf0Result << " & " << ef0result << " \\\\ " << endl;

		    cout << " frresult with analysis code : " << frresult << " +- " << efrresult << endl;
		    cout << " frresult with class         : " << minuitFitter.GetFRResult() << " +- " << minuitFitter.GetFRError() << endl;

		  } 
		  else {
		    frresult = 1.-minimizer.X()[0]-minimizer.X()[1];
		    double er0=minimizer.Errors()[0];
		    double er1=minimizer.Errors()[1];
		    double cov01 = minimizer.CovMatrix(0,1);
		    
		    efrresult = sqrt( er0*er0 + er1*er1 + 2.*cov01);
		    PresentationTex << frresult << " & " << frresult-SMfrResult <<" & " << efrresult << " & " << flresult << " & " << flresult-SMflResult << " & " << eflresult << " & " << f0result << " & " << f0result-SMf0Result << " & " << ef0result << " \\\\ " << endl;
		  }

		} 
		PresentationTex << " \\hline " << endl;
		
		LengthOfPresentationArray++;
		for (int ibinn=0; ibinn<CosThetaBinNumber; ibinn++){
		  delete genttbarhisto[ibinn];
		}
	      }//End of Minuit fitter loop
	      
	      TCHE++;
	      if(TCHE==13) {
		TCHP=1;
		SSVHE=1;
		SSVHP=1;
		CSV=1;
	      }	      

	    }//end of tche while loop
	    TCHE=13;	    	    
	    
	    std::string CosThetaDataString = "CosThetaData_TCHE"+bTagValues[TCHE]+"_TCHP"+bTagValues[TCHP]+"_SSVHE"+bTagValues[SSVHE]+"_SSVHP"+bTagValues[SSVHP]+"_CSV"+bTagValues[CSV];
	    std::string CosThetaSignalString = "CosThetaSignal_TCHE"+bTagValues[TCHE]+"_TCHP"+bTagValues[TCHP]+"_SSVHE"+bTagValues[SSVHE]+"_SSVHP"+bTagValues[SSVHP]+"_CSV"+bTagValues[CSV];
	    std::string CosThetaBckgString = "CosThetaBckg_TCHE"+bTagValues[TCHE]+"_TCHP"+bTagValues[TCHP]+"_SSVHE"+bTagValues[SSVHE]+"_SSVHP"+bTagValues[SSVHP]+"_CSV"+bTagValues[CSV];
	    std::string CosThetaGeneratorString = "CosThetaGenerator_TCHE"+bTagValues[TCHE]+"_TCHP"+bTagValues[TCHP]+"_SSVHE"+bTagValues[SSVHE]+"_SSVHP"+bTagValues[SSVHP]+"_CSV"+bTagValues[CSV];
	    
	    if(iDataSet == 0){
	      histo1D[CosThetaDataString]=new TH1F(CosThetaDataString.c_str(),CosThetaDataString.c_str(),CosThetaBinNumber,-1,1);
	      histo1D[CosThetaDataString]->SetDirectory(0);
	      histo1D[CosThetaSignalString]=new TH1F(CosThetaSignalString.c_str(),CosThetaSignalString.c_str(),CosThetaBinNumber,-1,1);
	      histo1D[CosThetaSignalString]->SetDirectory(0);
	      histo1D[CosThetaBckgString]=new TH1F(CosThetaBckgString.c_str(),CosThetaBckgString.c_str(),CosThetaBinNumber,-1,1);
	      histo1D[CosThetaBckgString]->SetDirectory(0);
	    }
	    
	    if((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ){   //Filling of the genttbar histo
	      char hisname[100];
	      sprintf(hisname,"CosThetaGenerator_TCHE%s_TCHP%s_SSVHE%s_SSVHP%s_CSV%s", bTagValues[TCHE].c_str(),bTagValues[TCHP].c_str(),bTagValues[SSVHE].c_str(),bTagValues[SSVHP].c_str(),bTagValues[CSV].c_str());
	      for (int ibinn=0; ibinn<CosThetaBinNumber; ibinn++){
		genttbarhisto[ibinn]= new TNtuple(hisname,hisname,"costhgen:evtweight");
		genttbarhisto[ibinn]->SetDirectory(0);
	      }
	      for(int ii=0; ii<CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV].size();ii++){		
		for(int iBin=0; iBin< CosThetaBinNumber; iBin++){
		  
		  if(CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV][ii] >= binEdge[iBin] && CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV][ii] < binEdge[iBin+1]){
		    double costhgood = CosThGen[TCHE+TCHP+SSVHE+SSVHP+CSV][ii]; 
		    double thisWeight = EventCorrectionWeight[TCHE+TCHP+SSVHE+SSVHP+CSV][ii];
		    genttbarhisto[iBin]->Fill(costhgood , thisWeight) ; 
		  }
		  else if(CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV][ii] ==1){ //1 is included in last bin
		    genttbarhisto[CosThetaBinNumber-1]->Fill(CosThGen[TCHE+TCHP+SSVHE+SSVHP+CSV][ii], EventCorrectionWeight[TCHE+TCHP+SSVHE+SSVHP+CSV][ii]);
		  }	      	   
		}
	      }
	    }
	    
	    float scaleFactor = 1.;  //Make an array of the scaleFactor such that it can be different for every event !!
	    	    
	    for(int ii=0;ii<CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV].size();ii++){
	      //Change data to systematics since then nominal values will be reweighted!!
	      if(WSystResults == true && dataSetName.find("WJets_Nominal") != 0)  //WJets syst
		histo1D[CosThetaDataString]->Fill(CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV][ii],scaleFactor*(LumiWeightVector[TCHE+TCHP+SSVHE+SSVHP+CSV][ii])*Luminosity*dataSet->NormFactor());    
	      if(WSystResults == false && (dataSetName.find("Data") == 0 || dataSetName.find("JES") == 0))//Nominal and JES systematics
		histo1D[CosThetaDataString]->Fill(CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV][ii],scaleFactor*(LumiWeightVector[TCHE+TCHP+SSVHE+SSVHP+CSV][ii])*Luminosity*dataSet->NormFactor());    
	      
	      if(((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ) && dataSetName.find("JES_") != 0)
		histo1D[CosThetaSignalString]->Fill(CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV][ii],Luminosity*scaleFactor*(LumiWeightVector[TCHE+TCHP+SSVHE+SSVHP+CSV][ii])*NominalNormFactor);
	      else if(dataSetName.find("Data") !=0 && dataSetName.find("JES") != 0 && dataSetName.find("Syst") != 0) 
		histo1D[CosThetaBckgString]->Fill(CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV][ii],Luminosity*scaleFactor*(LumiWeightVector[TCHE+TCHP+SSVHE+SSVHP+CSV][ii])*NominalNormFactor); 
	    }
	    
	    if(iDataSet==(datasets.size()-1)){//Go in this loop when the last datasample is active
	      
	      MyChi2 myChi2 (histo1D[CosThetaDataString], histo1D[CosThetaSignalString], histo1D[CosThetaBckgString], 0.6671, 0.3325, 0.0004); 
	      
	      TMinuitMinimizer minimizer;
	      minimizer.SetFunction(myChi2);
	      minimizer.SetErrorDef(1.0);   // 1.0 for chi2, 0.5 for -logL
	      minimizer.SetVariable(0, "F0", 0.6671, 1.e-4);
	      minimizer.SetVariable(1, "FL", 0.3325, 1.e-4);
	      if (ndimen==3)  minimizer.SetVariable(2, "Normal", 1., 1.e-4); 
	      
	      bool isValid = minimizer.Minimize();
	      double f0result, frresult,flresult;
	      double ef0result, efrresult,eflresult;
	      if (isValid) {
		//minimizer.PrintResults();
		f0result =minimizer.X()[0]; ef0result= minimizer.Errors()[0];
		flresult= minimizer.X()[1]; eflresult= minimizer.Errors()[1];
		
		PresentationTex << PresentationOutput[TCHE+TCHP+SSVHE+SSVHP+CSV] << " & ";
		for(int ii=0; ii< nameDataSet.size(); ii++){
		  if(ii < nameDataSet.size()-1){ 
		    PresentationTex << NumberRemainingEvents[TCHE+TCHP+SSVHE+SSVHP+CSV][ii] << " & ";
		  }
		  else if(ii == nameDataSet.size()-1 ){ 
		    PresentationTex << NumberRemainingEvents[TCHE+TCHP+SSVHE+SSVHP+CSV][ii] << " & "; //For presentation this is not the end of the table, helicity values still need to be included !!
		    PresentationTex << NumberBLeptCorrectEvents[TCHE+TCHP+SSVHE+SSVHP+CSV][ii] << " & ";
		  }
		}
		
		if (ndimen==3){
		  frresult = 1.-minimizer.X()[0]-minimizer.X()[1];
		  double er0=minimizer.Errors()[0];
		  double er1=minimizer.Errors()[1];
		  double cov01 = minimizer.CovMatrix(0,1);		  
		  efrresult = sqrt( er0*er0 + er1*er1 + 2.*cov01);

		  PresentationTex <<frresult<< " & " << frresult-SMfrResult <<" & " << efrresult << " & " << flresult << " & " << flresult-SMflResult << " & " << eflresult << " & " << f0result << " & " << f0result-SMf0Result << " & " << ef0result << " \\\\ " << endl;
		} 
		else {
		  frresult = 1.-minimizer.X()[0]-minimizer.X()[1];
		  double er0=minimizer.Errors()[0];
		  double er1=minimizer.Errors()[1];
		  double cov01 = minimizer.CovMatrix(0,1);		  
		  efrresult = sqrt( er0*er0 + er1*er1 + 2.*cov01);

		  PresentationTex <<frresult<< " & " << frresult-SMfrResult <<" & " << efrresult << " & " << flresult << " & " << flresult-SMflResult << " & " << eflresult << " & " << f0result << " & " << f0result-SMf0Result << " & " << ef0result << " \\\\ " << endl;
		}
	      } 
	      PresentationTex << " \\hline " << endl;
	      
	      for (int ibinn=0; ibinn<CosThetaBinNumber; ibinn++){
		delete genttbarhisto[ibinn];
	      }

	    }//End of Minuit fitter loop	      
	    
	    TCHP++;
	  }
	  TCHP=13;	  	  
	  
	  std::string CosThetaDataString = "CosThetaData_TCHE"+bTagValues[TCHE]+"_TCHP"+bTagValues[TCHP]+"_SSVHE"+bTagValues[SSVHE]+"_SSVHP"+bTagValues[SSVHP]+"_CSV"+bTagValues[CSV];
	  std::string CosThetaSignalString = "CosThetaSignal_TCHE"+bTagValues[TCHE]+"_TCHP"+bTagValues[TCHP]+"_SSVHE"+bTagValues[SSVHE]+"_SSVHP"+bTagValues[SSVHP]+"_CSV"+bTagValues[CSV];
	  std::string CosThetaBckgString = "CosThetaBckg_TCHE"+bTagValues[TCHE]+"_TCHP"+bTagValues[TCHP]+"_SSVHE"+bTagValues[SSVHE]+"_SSVHP"+bTagValues[SSVHP]+"_CSV"+bTagValues[CSV];
	  std::string CosThetaGeneratorString = "CosThetaGenerator_TCHE"+bTagValues[TCHE]+"_TCHP"+bTagValues[TCHP]+"_SSVHE"+bTagValues[SSVHE]+"_SSVHP"+bTagValues[SSVHP]+"_CSV"+bTagValues[CSV];
	  
	  if(iDataSet == 0){
	    histo1D[CosThetaDataString]=new TH1F(CosThetaDataString.c_str(),CosThetaDataString.c_str(),CosThetaBinNumber,-1,1);
	    histo1D[CosThetaDataString]->SetDirectory(0);
	    histo1D[CosThetaSignalString]=new TH1F(CosThetaSignalString.c_str(),CosThetaSignalString.c_str(),CosThetaBinNumber,-1,1);
	    histo1D[CosThetaSignalString]->SetDirectory(0);
	    histo1D[CosThetaBckgString]=new TH1F(CosThetaBckgString.c_str(),CosThetaBckgString.c_str(),CosThetaBinNumber,-1,1);
	    histo1D[CosThetaBckgString]->SetDirectory(0);
	  }
	  
	  if((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ){   //Filling of the genttbar histo
	    char hisname[100];
	    sprintf(hisname,"CosThetaGenerator_TCHE%s_TCHP%s_SSVHE%s_SSVHP%s_CSV%s", bTagValues[TCHE].c_str(),bTagValues[TCHP].c_str(),bTagValues[SSVHE].c_str(),bTagValues[SSVHP].c_str(),bTagValues[CSV].c_str());
	    for (int ibinn=0; ibinn<CosThetaBinNumber; ibinn++){
	      genttbarhisto[ibinn]= new TNtuple(hisname,hisname,"costhgen:evtweight");
	      genttbarhisto[ibinn]->SetDirectory(0);
	    }
	    for(int ii=0; ii<CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV].size();ii++){		
	      for(int iBin=0; iBin< CosThetaBinNumber; iBin++){
		
		if(CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV][ii] >= binEdge[iBin] && CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV][ii] < binEdge[iBin+1]){
		  double costhgood = CosThGen[TCHE+TCHP+SSVHE+SSVHP+CSV][ii]; 
		  double thisWeight = EventCorrectionWeight[TCHE+TCHP+SSVHE+SSVHP+CSV][ii];
		  genttbarhisto[iBin]->Fill(costhgood , thisWeight) ; 
		}
		else if(CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV][ii] ==1){ //1 is included in last bin
		  genttbarhisto[CosThetaBinNumber-1]->Fill(CosThGen[TCHE+TCHP+SSVHE+SSVHP+CSV][ii], EventCorrectionWeight[TCHE+TCHP+SSVHE+SSVHP+CSV][ii]);
		}	      	   
	      }
	    }
	  }
	  
	  float scaleFactor = 1.;  //Make an array of the scaleFactor such that it can be different for every event !!
	  	  
	  for(int ii=0;ii<CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV].size();ii++){
	    //Change data to systematics since then nominal values will be reweighted!!
	    if(WSystResults == true && dataSetName.find("WJets_Nominal") != 0)  //WJets syst
	      histo1D[CosThetaDataString]->Fill(CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV][ii],scaleFactor*LumiWeightVector[TCHE+TCHP+SSVHE+SSVHP+CSV][ii]*Luminosity*dataSet->NormFactor());    
	    if(WSystResults == false && (dataSetName.find("Data") == 0 || dataSetName.find("JES") == 0))//Nominal and JES systematics	    
	      histo1D[CosThetaDataString]->Fill(CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV][ii],scaleFactor*LumiWeightVector[TCHE+TCHP+SSVHE+SSVHP+CSV][ii]*Luminosity*dataSet->NormFactor());    
	    
	    if(((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true))  && dataSetName.find("JES_") != 0)   
	      histo1D[CosThetaSignalString]->Fill(CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV][ii],Luminosity*scaleFactor*LumiWeightVector[TCHE+TCHP+SSVHE+SSVHP+CSV][ii]*NominalNormFactor);
	    else if(dataSetName.find("WJets_Nominal")==0 || dataSetName.find("TTbarJets_Other") ==0) 
	      histo1D[CosThetaBckgString]->Fill(CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV][ii],Luminosity*scaleFactor*LumiWeightVector[TCHE+TCHP+SSVHE+SSVHP+CSV][ii]*NominalNormFactor); 
	  }
	  
	  if(iDataSet==(datasets.size()-1)){//Go in this loop when the last datasample is active

	    
	    MyChi2 myChi2 (histo1D[CosThetaDataString], histo1D[CosThetaSignalString], histo1D[CosThetaBckgString], 0.6671, 0.3325, 0.0004); 
	    
	    TMinuitMinimizer minimizer;
	    minimizer.SetFunction(myChi2);
	    minimizer.SetErrorDef(1.0);   // 1.0 for chi2, 0.5 for -logL
	    minimizer.SetVariable(0, "F0", 0.6671, 1.e-4);
	    minimizer.SetVariable(1, "FL", 0.3325, 1.e-4);
	    if (ndimen==3)  minimizer.SetVariable(2, "Normal", 1., 1.e-4); 
	    
	    bool isValid = minimizer.Minimize();
	    double f0result, frresult,flresult;
	    double ef0result, efrresult,eflresult;
	    if (isValid) {
	      //minimizer.PrintResults();
	      f0result =minimizer.X()[0]; ef0result= minimizer.Errors()[0];
	      flresult= minimizer.X()[1]; eflresult= minimizer.Errors()[1];
	      
	      
	      PresentationTex << PresentationOutput[TCHE+TCHP+SSVHE+SSVHP+CSV] << " & ";
	      for(int ii=0; ii< nameDataSet.size(); ii++){
		if(ii < nameDataSet.size()-1){ 
		  PresentationTex << NumberRemainingEvents[TCHE+TCHP+SSVHE+SSVHP+CSV][ii] << " & ";
		}
		else if(ii == nameDataSet.size()-1 ){ 
		  PresentationTex << NumberRemainingEvents[TCHE+TCHP+SSVHE+SSVHP+CSV][ii] << " & "; //For presentation this is not the end of the table, helicity values still need to be included !!
		  PresentationTex << NumberBLeptCorrectEvents[TCHE+TCHP+SSVHE+SSVHP+CSV][ii] << " & ";
		}
	      }
	      
	      if (ndimen==3){
		frresult = 1.-minimizer.X()[0]-minimizer.X()[1];
		double er0=minimizer.Errors()[0];
		double er1=minimizer.Errors()[1];
		double cov01 = minimizer.CovMatrix(0,1);		  
		efrresult = sqrt( er0*er0 + er1*er1 + 2.*cov01);
		
		PresentationTex << frresult << " & " << frresult-SMfrResult <<" & " << efrresult << " & " << flresult << " & " << flresult-SMflResult << " & " << eflresult << " & " << f0result << " & " << f0result-SMf0Result << " & " << ef0result << " \\\\ " << endl;
	      } 
	      else {
		frresult = 1.-minimizer.X()[0]-minimizer.X()[1];
		double er0=minimizer.Errors()[0];
		double er1=minimizer.Errors()[1];
		double cov01 = minimizer.CovMatrix(0,1);		
		efrresult = sqrt( er0*er0 + er1*er1 + 2.*cov01);

		PresentationTex << frresult << " & " << frresult-SMfrResult <<" & " << efrresult << " & " << flresult << " & " << flresult-SMflResult << " & " << eflresult << " & " << f0result << " & " << f0result-SMf0Result << " & " << ef0result << " \\\\ " << endl;
	      }
	    } 
	    PresentationTex << " \\hline " << endl;
	    
	    for (int ibinn=0; ibinn<CosThetaBinNumber; ibinn++){
	      delete genttbarhisto[ibinn];
	    }

	  }//End of Minuit fitter loop	      
	  
	  SSVHE++;
	}
	SSVHE=13;		
	
	std::string CosThetaDataString = "CosThetaData_TCHE"+bTagValues[TCHE]+"_TCHP"+bTagValues[TCHP]+"_SSVHE"+bTagValues[SSVHE]+"_SSVHP"+bTagValues[SSVHP]+"_CSV"+bTagValues[CSV];
	std::string CosThetaSignalString = "CosThetaSignal_TCHE"+bTagValues[TCHE]+"_TCHP"+bTagValues[TCHP]+"_SSVHE"+bTagValues[SSVHE]+"_SSVHP"+bTagValues[SSVHP]+"_CSV"+bTagValues[CSV];
	std::string CosThetaBckgString = "CosThetaBckg_TCHE"+bTagValues[TCHE]+"_TCHP"+bTagValues[TCHP]+"_SSVHE"+bTagValues[SSVHE]+"_SSVHP"+bTagValues[SSVHP]+"_CSV"+bTagValues[CSV];
	std::string CosThetaGeneratorString = "CosThetaGenerator_TCHE"+bTagValues[TCHE]+"_TCHP"+bTagValues[TCHP]+"_SSVHE"+bTagValues[SSVHE]+"_SSVHP"+bTagValues[SSVHP]+"_CSV"+bTagValues[CSV];
	
	if(iDataSet == 0){
	  histo1D[CosThetaDataString]=new TH1F(CosThetaDataString.c_str(),CosThetaDataString.c_str(),CosThetaBinNumber,-1,1);
	  histo1D[CosThetaDataString]->SetDirectory(0);
	  histo1D[CosThetaSignalString]=new TH1F(CosThetaSignalString.c_str(),CosThetaSignalString.c_str(),CosThetaBinNumber,-1,1);
	  histo1D[CosThetaSignalString]->SetDirectory(0);
	  histo1D[CosThetaBckgString]=new TH1F(CosThetaBckgString.c_str(),CosThetaBckgString.c_str(),CosThetaBinNumber,-1,1);
	  histo1D[CosThetaBckgString]->SetDirectory(0);
	}
	
	if((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ){   //Filling of the genttbar histo
	  char hisname[100];
	  sprintf(hisname,"CosThetaGenerator_TCHE%s_TCHP%s_SSVHE%s_SSVHP%s_CSV%s", bTagValues[TCHE].c_str(),bTagValues[TCHP].c_str(),bTagValues[SSVHE].c_str(),bTagValues[SSVHP].c_str(),bTagValues[CSV].c_str());
	  for (int ibinn=0; ibinn<CosThetaBinNumber; ibinn++){
	    genttbarhisto[ibinn]= new TNtuple(hisname,hisname,"costhgen:evtweight");
	    genttbarhisto[ibinn]->SetDirectory(0);
	  }
	  for(int ii=0; ii<CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV].size();ii++){		
	    for(int iBin=0; iBin< CosThetaBinNumber; iBin++){
	      
	      if(CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV][ii] >= binEdge[iBin] && CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV][ii] < binEdge[iBin+1]){
		double costhgood = CosThGen[TCHE+TCHP+SSVHE+SSVHP+CSV][ii]; 
		double thisWeight = EventCorrectionWeight[TCHE+TCHP+SSVHE+SSVHP+CSV][ii];
		genttbarhisto[iBin]->Fill(costhgood , thisWeight) ; 
	      }
	      else if(CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV][ii] ==1){ //1 is included in last bin
		genttbarhisto[CosThetaBinNumber-1]->Fill(CosThGen[TCHE+TCHP+SSVHE+SSVHP+CSV][ii], EventCorrectionWeight[TCHE+TCHP+SSVHE+SSVHP+CSV][ii]);
	      }	      	   
	    }
	  }
	}
	
	float scaleFactor = 1.;  //Make an array of the scaleFactor such that it can be different for every event !!
	
	for(int ii=0;ii<CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV].size();ii++){
	  //Change data to systematics since then nominal values will be reweighted!!
	  if(WSystResults == true && dataSetName.find("WJets_Nominal") != 0)  //WJets syst
	    histo1D[CosThetaDataString]->Fill(CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV][ii],scaleFactor*LumiWeightVector[TCHE+TCHP+SSVHE+SSVHP+CSV][ii]*Luminosity*dataSet->NormFactor());    
	  if(WSystResults == false && (dataSetName.find("Data") == 0 || dataSetName.find("JES") == 0))//Nominal and JES systematics
	    histo1D[CosThetaDataString]->Fill(CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV][ii],scaleFactor*LumiWeightVector[TCHE+TCHP+SSVHE+SSVHP+CSV][ii]*Luminosity*dataSet->NormFactor());    
	  
	    if(((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ) && dataSetName.find("JES_") != 0)
	    histo1D[CosThetaSignalString]->Fill(CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV][ii],Luminosity*scaleFactor*LumiWeightVector[TCHE+TCHP+SSVHE+SSVHP+CSV][ii]*NominalNormFactor);
	  else if(dataSetName.find("WJets_Nominal")==0 || dataSetName.find("TTbarJets_Other") ==0) //Syst
	    histo1D[CosThetaBckgString]->Fill(CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV][ii],Luminosity*scaleFactor*LumiWeightVector[TCHE+TCHP+SSVHE+SSVHP+CSV][ii]*NominalNormFactor); 
	}
	
	if(iDataSet==(datasets.size()-1)){//Go in this loop when the last datasample is active
	  
	  MyChi2 myChi2 (histo1D[CosThetaDataString], histo1D[CosThetaSignalString], histo1D[CosThetaBckgString], 0.6671, 0.3325, 0.0004); 
	  
	  TMinuitMinimizer minimizer;
	  minimizer.SetFunction(myChi2);
	  minimizer.SetErrorDef(1.0);   // 1.0 for chi2, 0.5 for -logL
	  minimizer.SetVariable(0, "F0", 0.6671, 1.e-4);
	  minimizer.SetVariable(1, "FL", 0.3325, 1.e-4);
	  if (ndimen==3)  minimizer.SetVariable(2, "Normal", 1., 1.e-4); 
	  
	  bool isValid = minimizer.Minimize();
	  double f0result, frresult,flresult;
	  double ef0result, efrresult,eflresult;
	  if (isValid) {
	    //minimizer.PrintResults();
	    f0result =minimizer.X()[0]; ef0result= minimizer.Errors()[0];
	    flresult= minimizer.X()[1]; eflresult= minimizer.Errors()[1];
	    
	    PresentationTex << PresentationOutput[TCHE+TCHP+SSVHE+SSVHP+CSV] << " & ";
	    for(int ii=0; ii< nameDataSet.size(); ii++){
	      if(ii < nameDataSet.size()-1){ 
		PresentationTex << NumberRemainingEvents[TCHE+TCHP+SSVHE+SSVHP+CSV][ii] << " & ";
	      }
	      else if(ii == nameDataSet.size()-1 ){ 
		PresentationTex << NumberRemainingEvents[TCHE+TCHP+SSVHE+SSVHP+CSV][ii] << " & "; //For presentation this is not the end of the table, helicity values still need to be included !!
		PresentationTex << NumberBLeptCorrectEvents[TCHE+TCHP+SSVHE+SSVHP+CSV][ii] << " & ";
	      }
	    }
	    
	    if (ndimen==3){
	      frresult = 1.-minimizer.X()[0]-minimizer.X()[1];
	      double er0=minimizer.Errors()[0];
	      double er1=minimizer.Errors()[1];
	      double cov01 = minimizer.CovMatrix(0,1);		  
	      efrresult = sqrt( er0*er0 + er1*er1 + 2.*cov01);
	      
	      PresentationTex << frresult << " & " << frresult-SMfrResult <<" & " << efrresult << " & " << flresult << " & " << flresult-SMflResult << " & " << eflresult << " & " << f0result << " & " << f0result-SMf0Result << " & " << ef0result << " \\\\ " << endl;
	    } 
	    else {
	      frresult = 1.-minimizer.X()[0]-minimizer.X()[1];
	      double er0=minimizer.Errors()[0];
	      double er1=minimizer.Errors()[1];
	      double cov01 = minimizer.CovMatrix(0,1);	      
	      efrresult = sqrt( er0*er0 + er1*er1 + 2.*cov01);
	      
	      PresentationTex << frresult << " & " << frresult-SMfrResult <<" & " << efrresult << " & " << flresult << " & " << flresult-SMflResult << " & " << eflresult << " & " << f0result << " & " << f0result-SMf0Result << " & " << ef0result << " \\\\ " << endl;
	    }
	  } 
	  PresentationTex << " \\hline " << endl;
	  
	  for (int ibinn=0; ibinn<CosThetaBinNumber; ibinn++){
	    delete genttbarhisto[ibinn];
	  }

	}//End of Minuit fitter loop	
	
	SSVHP++;
      }    
      SSVHP=13;            
      
      std::string CosThetaDataString = "CosThetaData_TCHE"+bTagValues[TCHE]+"_TCHP"+bTagValues[TCHP]+"_SSVHE"+bTagValues[SSVHE]+"_SSVHP"+bTagValues[SSVHP]+"_CSV"+bTagValues[CSV];
      std::string CosThetaSignalString = "CosThetaSignal_TCHE"+bTagValues[TCHE]+"_TCHP"+bTagValues[TCHP]+"_SSVHE"+bTagValues[SSVHE]+"_SSVHP"+bTagValues[SSVHP]+"_CSV"+bTagValues[CSV];
      std::string CosThetaBckgString = "CosThetaBckg_TCHE"+bTagValues[TCHE]+"_TCHP"+bTagValues[TCHP]+"_SSVHE"+bTagValues[SSVHE]+"_SSVHP"+bTagValues[SSVHP]+"_CSV"+bTagValues[CSV];
      std::string CosThetaGeneratorString = "CosThetaGenerator_TCHE"+bTagValues[TCHE]+"_TCHP"+bTagValues[TCHP]+"_SSVHE"+bTagValues[SSVHE]+"_SSVHP"+bTagValues[SSVHP]+"_CSV"+bTagValues[CSV];
	
      if(iDataSet == 0){
	histo1D[CosThetaDataString]=new TH1F(CosThetaDataString.c_str(),CosThetaDataString.c_str(),CosThetaBinNumber,-1,1);
	histo1D[CosThetaDataString]->SetDirectory(0);
	histo1D[CosThetaSignalString]=new TH1F(CosThetaSignalString.c_str(),CosThetaSignalString.c_str(),CosThetaBinNumber,-1,1);
	histo1D[CosThetaSignalString]->SetDirectory(0);
	histo1D[CosThetaBckgString]=new TH1F(CosThetaBckgString.c_str(),CosThetaBckgString.c_str(),CosThetaBinNumber,-1,1);
	histo1D[CosThetaBckgString]->SetDirectory(0);
      }
	
      if((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ){   //Filling of the genttbar histo
	char hisname[100];
	sprintf(hisname,"CosThetaGenerator_TCHE%s_TCHP%s_SSVHE%s_SSVHP%s_CSV%s", bTagValues[TCHE].c_str(),bTagValues[TCHP].c_str(),bTagValues[SSVHE].c_str(),bTagValues[SSVHP].c_str(),bTagValues[CSV].c_str());
	for (int ibinn=0; ibinn<CosThetaBinNumber; ibinn++){
	  genttbarhisto[ibinn]= new TNtuple(hisname,hisname,"costhgen:evtweight");
	  genttbarhisto[ibinn]->SetDirectory(0);
	}
	for(int ii=0; ii<CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV].size();ii++){		
	  for(int iBin=0; iBin< CosThetaBinNumber; iBin++){
	      
	    if(CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV][ii] >= binEdge[iBin] && CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV][ii] < binEdge[iBin+1]){
	      double costhgood = CosThGen[TCHE+TCHP+SSVHE+SSVHP+CSV][ii]; 
	      double thisWeight = EventCorrectionWeight[TCHE+TCHP+SSVHE+SSVHP+CSV][ii];
	      genttbarhisto[iBin]->Fill(costhgood , thisWeight) ; 
	    }
	    else if(CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV][ii] ==1){ //1 is included in last bin
	      genttbarhisto[CosThetaBinNumber-1]->Fill(CosThGen[TCHE+TCHP+SSVHE+SSVHP+CSV][ii], EventCorrectionWeight[TCHE+TCHP+SSVHE+SSVHP+CSV][ii]);
	    }	      	   
	  }
	}
      }
	
      float scaleFactor = 1.;  //Make an array of the scaleFactor such that it can be different for every event !!
	
      for(int ii=0;ii<CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV].size();ii++){
	//Change data to systematics since then nominal values will be reweighted!!
	if(WSystResults == true && dataSetName.find("WJets_Nominal") != 0)  //WJets syst
	  histo1D[CosThetaDataString]->Fill(CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV][ii],scaleFactor*LumiWeightVector[TCHE+TCHP+SSVHE+SSVHP+CSV][ii]*Luminosity*dataSet->NormFactor());    
	if(WSystResults == false && (dataSetName.find("Data") == 0 || dataSetName.find("JES") == 0))//Nominal and JES systematics
	  histo1D[CosThetaDataString]->Fill(CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV][ii],scaleFactor*LumiWeightVector[TCHE+TCHP+SSVHE+SSVHP+CSV][ii]*Luminosity*dataSet->NormFactor());    
	  
	  if(((dataSetName.find("TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("TTbarJets_SemiEl") ==0 && semiElectron == true) ) && dataSetName.find("JES_") != 0)
	  histo1D[CosThetaSignalString]->Fill(CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV][ii],Luminosity*scaleFactor*LumiWeightVector[TCHE+TCHP+SSVHE+SSVHP+CSV][ii]*NominalNormFactor);
	else if(dataSetName.find("WJets_Nominal")==0 || dataSetName.find("TTbarJets_Other") ==0) //Syst
	  histo1D[CosThetaBckgString]->Fill(CosThetaValues[TCHE+TCHP+SSVHE+SSVHP+CSV][ii],Luminosity*scaleFactor*LumiWeightVector[TCHE+TCHP+SSVHE+SSVHP+CSV][ii]*NominalNormFactor); 
      }
	
      if(iDataSet==(datasets.size()-1)){//Go in this loop when the last datasample is active

	  
	MyChi2 myChi2 (histo1D[CosThetaDataString], histo1D[CosThetaSignalString], histo1D[CosThetaBckgString], 0.6671, 0.3325, 0.0004); 
	  
	TMinuitMinimizer minimizer;
	minimizer.SetFunction(myChi2);
	minimizer.SetErrorDef(1.0);   // 1.0 for chi2, 0.5 for -logL
	minimizer.SetVariable(0, "F0", 0.6671, 1.e-4);
	minimizer.SetVariable(1, "FL", 0.3325, 1.e-4);
	if (ndimen==3)  minimizer.SetVariable(2, "Normal", 1., 1.e-4); 
	  
	bool isValid = minimizer.Minimize();
	double f0result, frresult,flresult;
	double ef0result, efrresult,eflresult;
	if (isValid) {
	  //minimizer.PrintResults();
	  f0result =minimizer.X()[0]; ef0result= minimizer.Errors()[0];
	  flresult= minimizer.X()[1]; eflresult= minimizer.Errors()[1];
	    
	  PresentationTex << PresentationOutput[TCHE+TCHP+SSVHE+SSVHP+CSV] << " & ";
	  for(int ii=0; ii< nameDataSet.size(); ii++){
	    if(ii < nameDataSet.size()-1){ 
	      PresentationTex << NumberRemainingEvents[TCHE+TCHP+SSVHE+SSVHP+CSV][ii] << " & ";
	    }
	    else if(ii == nameDataSet.size()-1 ){ 
	      PresentationTex << NumberRemainingEvents[TCHE+TCHP+SSVHE+SSVHP+CSV][ii] << " & "; //For presentation this is not the end of the table, helicity values still need to be included !!
	      PresentationTex << NumberBLeptCorrectEvents[TCHE+TCHP+SSVHE+SSVHP+CSV][ii] << " & ";
	    }
	  }
	    
	  if (ndimen==3){
	    frresult = 1.-minimizer.X()[0]-minimizer.X()[1];
	    double er0=minimizer.Errors()[0];
	    double er1=minimizer.Errors()[1];
	    double cov01 = minimizer.CovMatrix(0,1);		  
	    efrresult = sqrt( er0*er0 + er1*er1 + 2.*cov01);
	      
	    PresentationTex << frresult << " & " << frresult-SMfrResult <<" & " << efrresult << " & " << flresult << " & " << flresult-SMflResult << " & " << eflresult << " & " << f0result << " & " << f0result-SMf0Result << " & " << ef0result << " \\\\ " << endl;
	  } 
	  else {
	    frresult = 1.-minimizer.X()[0]-minimizer.X()[1];
	    double er0=minimizer.Errors()[0];
	    double er1=minimizer.Errors()[1];
	    double cov01 = minimizer.CovMatrix(0,1);	      
	    efrresult = sqrt( er0*er0 + er1*er1 + 2.*cov01);
	      
	    PresentationTex << frresult << " & " << frresult-SMfrResult <<" & " << efrresult << " & " << flresult << " & " << flresult-SMflResult << " & " << eflresult << " & " << f0result << " & " << f0result-SMf0Result << " & " << ef0result << " \\\\ " << endl;
	  }
	} 
	PresentationTex << " \\hline " << endl;
	  
	for (int ibinn=0; ibinn<CosThetaBinNumber; ibinn++){
	  delete genttbarhisto[ibinn];
	}

      }//End of Minuit fitter loop	
      
      CSV++;
    }
        
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

  PresentationTex << " \\end{tabular} " << endl;
  PresentationTex << " \\end{center} " << endl;
  PresentationTex << " \\renewcommand{\\arraystretch}{0.9} " << endl;
  PresentationTex << " \\end{tiny} " << endl;
  PresentationTex << " \\end{table} " << endl;
  PresentationTex.close();		

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
  //TCanvas *CosThetaDiffVSRecoCanvasSignal = new TCanvas("CosThetaDiffVSRecoCanvasSignal","TProfile distribution of Cos #theta (reco-gen) vs Cos #theta (reco) for muon Pt 20 to inf GeV");
  TCanvas *CosThetaDiffVSRecoCanvasPt10ToInf = new TCanvas("CosThetaDiffVSRecoCanvasPt10ToInf","TProfile distribution of Cos #theta (reco-gen) vs Cos #theta (reco) for muon Pt 10 to inf GeV");
  TCanvas *CosThetaDiffVSRecoCanvasPt15ToInf = new TCanvas("CosThetaDiffVSRecoCanvasPt15ToInf","TProfile distribution of Cos #theta (reco-gen) vs Cos #theta (reco) for muon Pt 15 to inf GeV");
  TCanvas *CosThetaDiffVSRecoCanvasPt20ToInf = new TCanvas("CosThetaDiffVSRecoCanvasPt20ToInf","TProfile distribution of Cos #theta (reco-gen) vs Cos #theta (reco) for muon Pt 20 to inf GeV");
  TCanvas *CosThetaDiffVSRecoCanvasPt25ToInf = new TCanvas("CosThetaDiffVSRecoCanvasPt25ToInf","TProfile distribution of Cos #theta (reco-gen) vs Cos #theta (reco) for muon Pt 25 to inf GeV");
  TCanvas *CosThetaDiffVSRecoCanvasPt30ToInf = new TCanvas("CosThetaDiffVSRecoCanvasPt30ToInf","TProfile distribution of Cos #theta (reco-gen) vs Cos #theta (reco) for muon Pt 30 to inf GeV");
  TCanvas *CosThetaDiffVSRecoCanvasPt35ToInf = new TCanvas("CosThetaDiffVSRecoCanvasPt35ToInf","TProfile distribution of Cos #theta (reco-gen) vs Cos #theta (reco) for muon Pt 35 to inf GeV");
//   CosThetaDiffVSRecoCanvasSignal->cd();
//   CosThetaDiffVSRecoSignal->Draw();
//   CosThetaDiffVSRecoSignal->Write();
//   CosThetaDiffVSRecoCanvasSignal->Update();
//   CosThetaDiffVSRecoCanvasSignal->Write();
//   CosThetaDiffVSRecoCanvasSignal->SaveAs("CosThetaDiffVSRecoSignalProfile.png");
//   CosThetaDiffVSRecoCanvasSignal->SaveAs("CosThetaDiffVSRecoSignalProfile.root");

  CosThetaDiffVSRecoCanvasPt10ToInf->cd();
  CosThetaDiffVSRecoPt10ToInf->Draw();
  CosThetaDiffVSRecoPt10ToInf->Write();
  CosThetaDiffVSRecoCanvasPt10ToInf->Update();
  CosThetaDiffVSRecoCanvasPt10ToInf->Write();
  CosThetaDiffVSRecoCanvasPt10ToInf->SaveAs("CosThetaDiffVSRecoPt10ToInfProfile.png");
  CosThetaDiffVSRecoCanvasPt10ToInf->SaveAs("CosThetaDiffVSRecoPt10ToInfProfile.root");

  CosThetaDiffVSRecoCanvasPt15ToInf->cd();
  CosThetaDiffVSRecoPt15ToInf->Draw();
  CosThetaDiffVSRecoPt15ToInf->Write();
  CosThetaDiffVSRecoCanvasPt15ToInf->Update();
  CosThetaDiffVSRecoCanvasPt15ToInf->Write();
  CosThetaDiffVSRecoCanvasPt15ToInf->SaveAs("CosThetaDiffVSRecoPt15ToInfProfile.png");
  CosThetaDiffVSRecoCanvasPt15ToInf->SaveAs("CosThetaDiffVSRecoPt15ToInfProfile.root");

  CosThetaDiffVSRecoCanvasPt20ToInf->cd();
  CosThetaDiffVSRecoPt20ToInf->Draw();
  CosThetaDiffVSRecoPt20ToInf->Write();
  CosThetaDiffVSRecoCanvasPt20ToInf->Update();
  CosThetaDiffVSRecoCanvasPt20ToInf->Write();
  CosThetaDiffVSRecoCanvasPt20ToInf->SaveAs("CosThetaDiffVSRecoPt20ToInfProfile.png");
  CosThetaDiffVSRecoCanvasPt20ToInf->SaveAs("CosThetaDiffVSRecoPt20ToInfProfile.root");

  CosThetaDiffVSRecoCanvasPt25ToInf->cd();
  CosThetaDiffVSRecoPt25ToInf->Draw();
  CosThetaDiffVSRecoPt25ToInf->Write();
  CosThetaDiffVSRecoCanvasPt25ToInf->Update();
  CosThetaDiffVSRecoCanvasPt25ToInf->Write();
  CosThetaDiffVSRecoCanvasPt25ToInf->SaveAs("CosThetaDiffVSRecoPt25ToInfProfile.png");
  CosThetaDiffVSRecoCanvasPt25ToInf->SaveAs("CosThetaDiffVSRecoPt25ToInfProfile.root");

  CosThetaDiffVSRecoCanvasPt30ToInf->cd();
  CosThetaDiffVSRecoPt30ToInf->Draw();
  CosThetaDiffVSRecoPt30ToInf->Write();
  CosThetaDiffVSRecoCanvasPt30ToInf->Update();
  CosThetaDiffVSRecoCanvasPt30ToInf->Write();
  CosThetaDiffVSRecoCanvasPt30ToInf->SaveAs("CosThetaDiffVSRecoPt30ToInfProfile.png");
  CosThetaDiffVSRecoCanvasPt30ToInf->SaveAs("CosThetaDiffVSRecoPt30ToInfProfile.root");

  CosThetaDiffVSRecoCanvasPt35ToInf->cd();
  CosThetaDiffVSRecoPt35ToInf->Draw();
  CosThetaDiffVSRecoPt35ToInf->Write();
  CosThetaDiffVSRecoCanvasPt35ToInf->Update();
  CosThetaDiffVSRecoCanvasPt35ToInf->Write();
  CosThetaDiffVSRecoCanvasPt35ToInf->SaveAs("CosThetaDiffVSRecoPt35ToInfProfile.png");
  CosThetaDiffVSRecoCanvasPt35ToInf->SaveAs("CosThetaDiffVSRecoPt35ToInfProfile.root");

  th1dir->cd();

  //Calculate ratio histos for cos theta distribution reco/gen for different muon pt cuts:
//   histo1D["CosThetaRatioMuon20To30"] = (TH1F*) histo1D["CosThetaDistributionPtSmaller30Ratio"]->Clone();
//   histo1D["CosThetaRatioMuon20To30"]->Divide(histo1D["CosThetaDistributionPtSmaller30Gen"]);
//   histo1D["CosThetaRatioMuon20To30"]->SaveAs("CosThetaRatioMuon20To30.root");
//   histo1D["CosThetaRatioMuon20To30"]->SaveAs("CosThetaRatioMuon20To30.png");
//   histo1D["CosThetaRatioMuon20To30"]->Write();

  histo1D["CosThetaRatioMuon10ToInf"] = (TH1F*) histo1D["CosThetaRecoPt10ToInfForRatio"]->Clone();
  histo1D["CosThetaRatioMuon10ToInf"]->Divide(histo1D["CosThetaGenPt10ToInf"]);
  histo1D["CosThetaRatioMuon10ToInf"]->SaveAs("CosThetaRatioMuon10ToInf.root");
  histo1D["CosThetaRatioMuon10ToInf"]->SaveAs("CosThetaRatioMuon10ToInf.png");
  histo1D["CosThetaRatioMuon10ToInf"]->Write();

  histo1D["CosThetaRatioMuon15ToInf"] = (TH1F*) histo1D["CosThetaRecoPt15ToInfForRatio"]->Clone();
  histo1D["CosThetaRatioMuon15ToInf"]->Divide(histo1D["CosThetaGenPt15ToInf"]);
  histo1D["CosThetaRatioMuon15ToInf"]->SaveAs("CosThetaRatioMuon15ToInf.root");
  histo1D["CosThetaRatioMuon15ToInf"]->SaveAs("CosThetaRatioMuon15ToInf.png");
  histo1D["CosThetaRatioMuon15ToInf"]->Write();

  histo1D["CosThetaRatioMuon20ToInf"] = (TH1F*) histo1D["CosThetaRecoPt20ToInfForRatio"]->Clone();
  histo1D["CosThetaRatioMuon20ToInf"]->Divide(histo1D["CosThetaGenPt20ToInf"]);
  histo1D["CosThetaRatioMuon20ToInf"]->SaveAs("CosThetaRatioMuon20ToInf.root");
  histo1D["CosThetaRatioMuon20ToInf"]->SaveAs("CosThetaRatioMuon20ToInf.png");
  histo1D["CosThetaRatioMuon20ToInf"]->Write();

  histo1D["CosThetaRatioMuon25ToInf"] = (TH1F*) histo1D["CosThetaRecoPt25ToInfForRatio"]->Clone();
  histo1D["CosThetaRatioMuon25ToInf"]->Divide(histo1D["CosThetaGenPt25ToInf"]);
  histo1D["CosThetaRatioMuon25ToInf"]->SaveAs("CosThetaRatioMuon25ToInf.root");
  histo1D["CosThetaRatioMuon25ToInf"]->SaveAs("CosThetaRatioMuon25ToInf.png");
  histo1D["CosThetaRatioMuon25ToInf"]->Write();

  histo1D["CosThetaRatioMuon30ToInf"] = (TH1F*) histo1D["CosThetaRecoPt30ToInfForRatio"]->Clone();
  histo1D["CosThetaRatioMuon30ToInf"]->Divide(histo1D["CosThetaGenPt30ToInf"]);
  histo1D["CosThetaRatioMuon30ToInf"]->SaveAs("CosThetaRatioMuon30ToInf.root");
  histo1D["CosThetaRatioMuon30ToInf"]->SaveAs("CosThetaRatioMuon30ToInf.png");
  histo1D["CosThetaRatioMuon30ToInf"]->Write();

  histo1D["CosThetaRatioMuon35ToInf"] = (TH1F*) histo1D["CosThetaRecoPt35ToInfForRatio"]->Clone();
  histo1D["CosThetaRatioMuon35ToInf"]->Divide(histo1D["CosThetaGenPt35ToInf"]);
  histo1D["CosThetaRatioMuon35ToInf"]->SaveAs("CosThetaRatioMuon35ToInf.root");
  histo1D["CosThetaRatioMuon35ToInf"]->SaveAs("CosThetaRatioMuon35ToInf.png");
  histo1D["CosThetaRatioMuon35ToInf"]->Write();
  
  //Write resolution histos for separate bins in one canvas:
  TCanvas *ResolutionBinsNegCanvas = new TCanvas("ResolutionBinsNegCanvas","Resolution histos for separate bins (neg)");
  TCanvas *ResolutionBinsCanvas = new TCanvas("ResolutionBinsCanvas","Resolution histos for separate bins");
  ResolutionBinsNegCanvas->Divide(2,3);
  ResolutionBinsCanvas->Divide(2,3);

  ResolutionBinsNegCanvas->cd(1);
  histo1D["CosThetaResLastBin"]->Draw();
  ResolutionBinsNegCanvas->Update();
  ResolutionBinsNegCanvas->cd(2);
  histo1D["CosThetaResSecondNegBin"]->Draw();
  ResolutionBinsNegCanvas->Update();
  ResolutionBinsNegCanvas->cd(3);
  histo1D["CosThetaResThirdNegBin"]->Draw();
  ResolutionBinsNegCanvas->Update();
  ResolutionBinsNegCanvas->cd(4);
  histo1D["CosThetaResFourthNegBin"]->Draw();
  ResolutionBinsNegCanvas->Update();
  ResolutionBinsNegCanvas->cd(5);
  histo1D["CosThetaResFifthNegBin"]->Draw();
  ResolutionBinsNegCanvas->Update();
  ResolutionBinsNegCanvas->Write();
  ResolutionBinsNegCanvas->SaveAs("ResolutionBinsNegCanvas.png");
  
  ResolutionBinsCanvas->cd(1);
  histo1D["CosThetaResFifthBin"]->Draw();
  ResolutionBinsCanvas->Update();
  ResolutionBinsCanvas->cd(2);
  histo1D["CosThetaResFourthBin"]->Draw();
  ResolutionBinsCanvas->Update();
  ResolutionBinsCanvas->cd(3);
  histo1D["CosThetaResThirdBin"]->Draw();
  ResolutionBinsCanvas->Update();
  ResolutionBinsCanvas->cd(4);
  histo1D["CosThetaResSecondBin"]->Draw();
  ResolutionBinsCanvas->Update();
  ResolutionBinsCanvas->cd(5);
  histo1D["CosThetaResFirstBin"]->Draw();
  ResolutionBinsCanvas->Update();
  ResolutionBinsCanvas->Write();
  ResolutionBinsCanvas->SaveAs("ResolutionBinsCanvas.png");
  

//   //Calculate ratio for cut on probability of Kinematic Fit:
//   histo1D["CosThetaNobTagProbCut1D"]->Scale(100./(histo1D["CosThetaNobTagProbCut1D"]->Integral()));
//   histo1D["CosThetaNobTag1D"]->Scale(100./(histo1D["CosThetaNobTag1D"]->Integral()));
  
//   histo1D["ratioHisto"] = (TH1F*) histo1D["CosThetaNobTagProbCut1D"]->Clone();
//   //ratioHisto->SetName("Ratio");
//   //ratioHisto->Reset();
//   histo1D["CosThetaNobTagClone"] = (TH1F*) histo1D["CosThetaNobTag1D"]->Clone();
//   histo1D["ratioHisto"]->Divide(histo1D["CosThetaNobTagClone"]);
//   histo1D["ratioHisto"]->SaveAs("ratioHisto.root");
//   histo1D["ratioHisto"]->SaveAs("ratioHisto.png");

//   histo1D["KinFitRatio"] = (TH1F*) histo1D["KinFitProbCorrectBLept"]->Clone();
//   histo1D["KinFitProbWrongBLeptClone"] = (TH1F*) histo1D["KinFitProbWrongBLept"]->Clone();
//   histo1D["KinFitRatio"]->Divide(histo1D["KinFitProbWrongBLeptClone"]);
//   histo1D["KinFitRatio"]->SaveAs("KinFitRatioHisto.root");
//   histo1D["KinFitRatio"]->SaveAs("KinFitRatioHisto.png"); 

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
      tempCanvas->SaveAs( (pathPNG+"/CosThetaPlots/"+it->first+".png").c_str() );
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
