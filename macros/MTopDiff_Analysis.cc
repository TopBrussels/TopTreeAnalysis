#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <algorithm>

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

// RooFit librairies
#include "RooDataHist.h"
#include "RooHist.h"
#include "RooRealVar.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooGaussian.h"
#include "RooBifurGauss.h"
#include "RooLandau.h"
#include "RooCBShape.h"

#include "../Tools/interface/PlottingTools.h"
#include "../Tools/interface/MultiSamplePlot.h"
#include "../Content/interface/Dataset.h"
#include "../JESMeasurement/interface/LightMonster.h"
#include "../JESMeasurement/interface/MonsterTools.h"
#include "../JESMeasurement/IdeogramTools/interface/MassFitInfoFiller.h"
#include "../MCInformation/interface/LumiReWeighting.h"
#include "../MCInformation/interface/Lumi3DReWeighting.h"
#include "../MCInformation/interface/ResolutionFit.h"

#include "../../TopQuarkAnalysis/TopMassIdeogram/interface/MassFitInfoHeader.h"
#include "../../TopQuarkAnalysis/TopMassIdeogram/interface/MassFitInfoEvent.h"
#include "../../TopQuarkAnalysis/TopMassIdeogram/interface/MassLikelihoodTool.h"
#include "../../TopQuarkAnalysis/TopMassIdeogram/interface/MassLikelihoodWriter.h"

#include "Style.C"

using namespace std;
using namespace reweight;
using namespace RooFit;

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
  cout << " Beginning of the program for Top Mass Diff Analysis ! " << endl;
  cout << "*******************************************************" << endl;
 
//  setTDRStyle();
  setMyStyle();
  
  string pathPNG = "PlotsTopMassDiff/";
  mkdir(pathPNG.c_str(),0777);
  mkdir((pathPNG+"MSPlot/").c_str(),0777);
  
//  float LuminosityMu = 4656, LuminosityEl = 4666;
  float LuminosityMu = 4953, LuminosityEl = 4961;
  
  bool doAllMSPlots = false;
  
  cout << "Executing the Top Mass difference analysis for an integrated luminosity of semi-mu " << LuminosityMu << " pb^-1 and semi-el " << LuminosityEl << " pb^-1" << endl;
  
  vector<string> inputMonsters;
  if (argc >= 2)
		inputMonsters.push_back( string(argv[1]) );
  else
  {
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_TTbarJets_powheg_Nominal_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_TtbarJets_Chamonix_Nominal_SemiLep.root");
  
    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_Data_ElectronHad_4p7fb_Nominal_SemiLep.root");
    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_Data_MuHad_4p7fb_Nominal_SemiLep.root");
  
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_Data_ElectronHad_4p7fb_InvertedIso_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_Data_MuHad_4p7fb_InvertedIso_SemiLep.root");
  
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_ST_SingleTop_tChannel_tbar_Nominal_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_ST_SingleTop_tChannel_t_Nominal_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_ST_SingleTop_tWChannel_tbar_Nominal_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_ST_SingleTop_tWChannel_t_Nominal_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_ZJets_Nominal_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_WJets_Nominal_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_TTbarJets_Nominal_SemiLep.root");
  
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_QCD_Mu15_Nominal_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_QCD_Pt-20to30_BCtoE_Nominal_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_QCD_Pt-30to80_BCtoE_Nominal_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_QCD_Pt-80to170_BCtoE_Nominal_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_QCD_Pt-20to30_EMEnriched_Nominal_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_QCD_Pt-30to80_EMEnriched_Nominal_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_QCD_Pt-80to170_EMEnriched_Nominal_SemiLep.root");
  
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_TTbarJets_mass161_5_Nominal_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_TTbarJets_mass163_5_Nominal_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_TTbarJets_mass166_5_Nominal_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_TTbarJets_mass169_5_Nominal_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_TTbarJets_mass175_5_Nominal_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_TTbarJets_mass178_5_Nominal_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_TTbarJets_mass181_5_Nominal_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_TTbarJets_mass184_5_Nominal_SemiLep.root");
  
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_ST_SingleTop_tChannel_tbar_JERMinus_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_ST_SingleTop_tChannel_tbar_JERPlus_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_ST_SingleTop_tChannel_tbar_JESMinus_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_ST_SingleTop_tChannel_tbar_JESPlus_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_ST_SingleTop_tChannel_tbar_AlignMinus_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_ST_SingleTop_tChannel_tbar_AlignPlus_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_ST_SingleTop_tChannel_t_JERMinus_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_ST_SingleTop_tChannel_t_JERPlus_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_ST_SingleTop_tChannel_t_JESMinus_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_ST_SingleTop_tChannel_t_JESPlus_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_ST_SingleTop_tChannel_t_AlignMinus_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_ST_SingleTop_tChannel_t_AlignPlus_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_ST_SingleTop_tWChannel_tbar_JERMinus_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_ST_SingleTop_tWChannel_tbar_JERPlus_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_ST_SingleTop_tWChannel_tbar_JESMinus_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_ST_SingleTop_tWChannel_tbar_JESPlus_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_ST_SingleTop_tWChannel_tbar_AlignMinus_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_ST_SingleTop_tWChannel_tbar_AlignPlus_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_ST_SingleTop_tWChannel_t_JERMinus_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_ST_SingleTop_tWChannel_t_JERPlus_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_ST_SingleTop_tWChannel_t_JESMinus_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_ST_SingleTop_tWChannel_t_JESPlus_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_ST_SingleTop_tWChannel_t_AlignMinus_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_ST_SingleTop_tWChannel_t_AlignPlus_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_TTbarJets_JERMinus_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_TTbarJets_JERPlus_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_TTbarJets_JESMinus_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_TTbarJets_JESPlus_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_TTbarJets_AlignMinus_SemiLep.root"); broken file
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_TTbarJets_AlignPlus_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_WJets_JERMinus_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_WJets_JERPlus_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_WJets_JESMinus_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_WJets_JESPlus_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_WJets_AlignMinus_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_WJets_AlignPlus_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_ZJets_JERMinus_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_ZJets_JERPlus_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_ZJets_JESMinus_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_ZJets_JESPlus_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_ZJets_AlignMinus_SemiLep.root");
//    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_ZJets_AlignPlus_SemiLep.root");
    
/*    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_TTbarJets_matchingdown_Nominal_SemiLep.root");
    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_TTbarJets_matchingup_Nominal_SemiLep.root");
    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_TTbarJets_scaledown_Nominal_SemiLep.root");
    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_TTbarJets_scaleup_Nominal_SemiLep.root");
    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_WJets_matchingdown_Nominal_SemiLep.root");
    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_WJets_matchingup_Nominal_SemiLep.root");
    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_WJets_scaledown_Nominal_SemiLep.root");
    inputMonsters.push_back("Monsters/KinFit_LightMonsters_TopMassDiff_WJets_scaleup_Nominal_SemiLep.root");*/
  }
  
  TFile *fout = new TFile ("MTopDiff_Analysis.root", "RECREATE");
  fout->cd();
  
  // initialize histograms
  cout << "Initializing histograms" << endl;
  
 	//control plots
 	histo1D["lumiWeights"] = new TH1F("lumiWeights","lumiWeights;lumiWeight;#events",100,0,4);
 	histo1D["lumiWeightsUp"] = new TH1F("lumiWeightsUp","lumiWeightsUp;lumiWeightUp;#events",100,0,4);
 	histo1D["lumiWeightsDown"] = new TH1F("lumiWeightsDown","lumiWeightsDown;lumiWeightDown;#events",100,0,4);
	histo1D["mTop_gen"] = new TH1F("mTop_gen","mTop_gen;mTop;#events",500,150,200);
	histo1D["mTopDiff_gen"] = new TH1F("mTopDiff_gen","mTopDiff_gen;mTop;#events",100,-5,5);
	
	histo1D["RunNr_Data_mu"] = new TH1F("RunNr_Data_mu","RunNr_Data_mu;RunNr;#events",180400-160400,160400,180400);
	histo1D["RunNr_Data_el"] = new TH1F("RunNr_Data_el","RunNr_Data_el;RunNr;#events",180400-160400,160400,180400);
	histo1D["RunNr_DataDrivenQCD_mu"] = new TH1F("RunNr_DataDrivenQCD_mu","RunNr_DataDrivenQCD_mu;RunNr;#events",180400-160400,160400,180400);
	histo1D["RunNr_DataDrivenQCD_el"] = new TH1F("RunNr_DataDrivenQCD_el","RunNr_DataDrivenQCD_el;RunNr;#events",180400-160400,160400,180400);
  histo2D["RunNr_Data_VS_Lepton_pt_mu"] = new TH2F("RunNr_Data_VS_Lepton_pt_mu","RunNr_Data_VS_Lepton_pt_mu;RunNr;Lepton Pt",180400-160400/100,160400,180400,50,0,250);
  histo2D["RunNr_Data_VS_Lepton_pt_el"] = new TH2F("RunNr_Data_VS_Lepton_pt_el","RunNr_Data_VS_Lepton_pt_el;RunNr;Lepton Pt",180400-160400/100,160400,180400,50,0,250);
  histo2D["RunNr_DataDrivenQCD_VS_Lepton_pt_mu"] = new TH2F("RunNr_DataDrivenQCD_VS_Lepton_pt_mu","RunNr_DataDrivenQCD_VS_Lepton_pt_mu;RunNr;Lepton Pt",180400-160400/100,160400,180400,50,0,250);
  histo2D["RunNr_DataDrivenQCD_VS_Lepton_pt_el"] = new TH2F("RunNr_DataDrivenQCD_VS_Lepton_pt_el","RunNr_DataDrivenQCD_VS_Lepton_pt_el;RunNr;Lepton Pt",180400-160400/100,160400,180400,50,0,250);
  histo2D["RunNr_Data_VS_Lepton_eta_mu"] = new TH2F("RunNr_Data_VS_Lepton_eta_mu","RunNr_Data_VS_Lepton_eta_mu;RunNr;Lepton Eta",180400-160400/100,160400,180400,50,-2.5,2.5);
  histo2D["RunNr_Data_VS_Lepton_eta_el"] = new TH2F("RunNr_Data_VS_Lepton_eta_el","RunNr_Data_VS_Lepton_eta_el;RunNr;Lepton Eta",180400-160400/100,160400,180400,50,-2.5,2.5);
  histo2D["RunNr_DataDrivenQCD_VS_Lepton_eta_mu"] = new TH2F("RunNr_DataDrivenQCD_VS_Lepton_eta_mu","RunNr_DataDrivenQCD_VS_Lepton_eta_mu;RunNr;Lepton Eta",180400-160400/100,160400,180400,50,-2.5,2.5);
  histo2D["RunNr_DataDrivenQCD_VS_Lepton_eta_el"] = new TH2F("RunNr_DataDrivenQCD_VS_Lepton_eta_el","RunNr_DataDrivenQCD_VS_Lepton_eta_el;RunNr;Lepton Eta",180400-160400/100,160400,180400,50,-2.5,2.5);
  
	histo1D["mTop_Fit"] = new TH1F("mTop_Fit","mTop_Fit;mTop;#jet-combis",125,0,500);
	histo1D["sigmaMtop_Fit"] = new TH1F("sigmaMtop_Fit","sigmaMtop_Fit;sigmaMtop;#jet-combis",100,0,50);
	histo1D["MinChi2ndf_Fit"] = new TH1F("MinChi2ndf_Fit","MinChi2ndf_Fit;sigmaMtop;#jet-combis",80,0,20);
	histo1D["Prob_MinChi2_Fit"] = new TH1F("Prob_MinChi2_Fit","Prob_MinChi2_Fit;Prob_MinChi2;#jet-combis",50,0,1);
  
  histo1D["mTop_NoFit_goodcombi"] = new TH1F("mTop_NoFit_goodcombi","mTop_NoFit_goodcombi;#events",100,50,250);
  histo1D["mTopL7_NoFit_goodcombi"] = new TH1F("mTopL7_NoFit_goodcombi","mTopL7_NoFit_goodcombi;#events",100,50,250);
	histo1D["mTop_Fit_goodcombi"] = new TH1F("mTop_Fit_goodcombi","mTop_Fit_goodcombi;mTop_goodcombi;#events",100,50,250);
	histo1D["sigmaMtop_Fit_goodcombi"] = new TH1F("sigmaMtop_Fit_goodcombi","sigmaMtop_Fit_goodcombi;sigmaMtop;#events",50,0,25);
	histo1D["MinChi2ndf_Fit_goodcombi"] = new TH1F("MinChi2ndf_Fit_goodcombi","MinChi2ndf_Fit_goodcombi;MinChi2ndf;#events",80,0,20);
	histo1D["Prob_MinChi2_Fit_goodcombi"] = new TH1F("Prob_MinChi2_Fit_goodcombi","Prob_MinChi2_Fit_goodcombi;Prob_MinChi2_goodcombi;#events",50,0,1);
	
	histo1D["nUsedCombis_TTJets"] = new TH1F("nUsedCombis_TTJets","nUsedCombis_TTJets",13,-0.5,12.5);
	histo1D["nUsedCombis_WJets"] = new TH1F("nUsedCombis_WJets","nUsedCombis_WJets",13,-0.5,12.5);
	histo1D["maxIndexHadrJets"] = new TH1F("maxIndexHadrJets","maxIndexHadrJets",10,-0.5,9.5);
  
  histo1D["bjets_SSVHE"] = new TH1F("bjets_SSVHE","bjets_SSVHE;bjets_SSVHE;#b-jets",60,0,6);
  
  histo1D["lightJets_pt"] = new TH1F("lightJets_pt","lightJets_pt;lightJets_pt;#jets",100,20,270);
  histo1D["lightJets_eta"] = new TH1F("lightJets_eta","lightJets_eta;lightJets_eta;#jets",100,-2.5,2.5);
  histo1D["bJets_pt"] = new TH1F("bJets_pt","bJets_pt;bJets_pt;#jets",100,20,270);
  histo1D["bJets_eta"] = new TH1F("bJets_eta","bJets_eta;bJets_eta;#jets",100,-2.5,2.5);
  
  histo1D["EtaMostForwardJet"] = new TH1F("EtaMostForwardJet","EtaMostForwardJet",50,-2.5,2.5);
  
  // Ideogram shapes
  histo1D["wJets_mTopFitted"] = new TH1F("wJets_mTopFitted","wJets_mTopFitted;Fitted top quark mass;#jet-combis",100,0,500);
  histo1D["wJets_mTopFitted_noChi2weight"] = new TH1F("wJets_mTopFitted_noChi2weight","wJets_mTopFitted_noChi2weight;Fitted top quark mass;#jet-combis",100,0,500);
  histo1D["wJets_mTopFitted_muPlus"] = new TH1F("wJets_mTopFitted_muPlus","wJets_mTopFitted_muPlus;Fitted top quark mass;#jet-combis",100,0,500);
  histo1D["wJets_mTopFitted_muMinus"] = new TH1F("wJets_mTopFitted_muMinus","wJets_mTopFitted_muMinus;Fitted top quark mass;#jet-combis",100,0,500);
  
  histo1D["ttJets_mTopFitted_badcombi"] = new TH1F("ttJets_mTopFitted_badcombi","ttJets_mTopFitted_badcomb;Fitted top quark mass;#jet-combis",200,0,500);
  histo1D["ttJets_mTopFitted_badcombi_noChi2weight"] = new TH1F("ttJets_mTopFitted_badcombi_noChi2weight","ttJets_mTopFitted_badcomb_noChi2weight;Fitted top quark mass;#jet-combis",200,0,500);
  histo1D["ttJets_mTopFitted_badcombi_muPlus"] = new TH1F("ttJets_mTopFitted_badcombi_muPlus","ttJets_mTopFitted_badcomb_muPlus;Fitted top quark mass;#jet-combis",200,0,500);
  histo1D["ttJets_mTopFitted_badcombi_muMinus"] = new TH1F("ttJets_mTopFitted_badcombi_muMinus","ttJets_mTopFitted_badcomb_muMinus;Fitted top quark mass;#jet-combis",200,0,500);
  
  histo1D["combiWeight"] = new TH1F("combiWeight","combiWeight",100,0,1.00001);
  histo1D["wJets_combiWeight"] = new TH1F("wJets_combiWeight","wJets_combiWeight",100,0,1.00001);
  histo2D["wJets_mTopFitted_VS_combiWeight"] = new TH2F("wJets_mTopFitted_VS_combiWeight","wJets_mTopFitted_VS_combiWeight",100,0,500,100,0,1.00001);
  
  cout << "Initializing MSPlots" << endl;
  vector<Dataset*> datasetsSemiMu, datasetsSemiEl, dataSetsAll, dataSetsOneData; // needed for MSPlots
  bool addedData = false;
  for(unsigned int iDataSet=0; iDataSet<inputMonsters.size(); iDataSet++)
  {
    TFile* inFile = new TFile(inputMonsters[iDataSet].c_str(),"READ");
    TTree* inConfigTree = (TTree*) inFile->Get("configTreeMonsterFile");
    TBranch* d_br = (TBranch*) inConfigTree->GetBranch("Dataset");
    TClonesArray* tc_dataset = new TClonesArray("Dataset",0);
    d_br->SetAddress(&tc_dataset);
    inConfigTree->GetEvent(0);
    Dataset* dataSet = (Dataset*) tc_dataset->At(0);
    int color = 0;
    if( dataSet->Name().find("Data_ElectronHad") == 0 )
    {
      if( dataSet->Name().find("InvertedIso") != string::npos )
      {
        dataSet->SetEquivalentLuminosity(LuminosityEl * 15133.7/2200.);
        dataSet->SetName("QCD_El_DAtaDriven");
        dataSet->SetTitle("QCD (DAta-driven)");
      }
      else
      {
        dataSet->SetTitle("Data_El");
        dataSet->SetName("Data");
        dataSet->SetEquivalentLuminosity(LuminosityEl);
      }
    }
    else if( dataSet->Name().find("Data_MuHad") == 0 )
    {
      if( dataSet->Name().find("InvertedIso") != string::npos )
      {
        dataSet->SetEquivalentLuminosity(LuminosityMu * 41106.4/675.6);
        dataSet->SetName("QCD_Mu_DAtaDriven");
        dataSet->SetTitle("QCD (DAta-driven)");
      }
      else
      {
        dataSet->SetTitle("Data_Mu");
        dataSet->SetName("Data");
        dataSet->SetEquivalentLuminosity(LuminosityMu);
      }
    }
    if( dataSet->Name().find("QCD") == 0 ) color = kYellow;
    if( dataSet->Name().find("TT") == 0 )
    {
      color = kRed+1;
      dataSet->SetTitle("t#bar{t}");
    }
    if( dataSet->Name().find("TTbarJets_Other") == 0 ) color = kRed-7;
    if( dataSet->Name().find("Chamonix") == 0 )
    {
      color = kRed-7;
      dataSet->SetTitle("TTJets_Chamonix");
    }
    if( dataSet->Name().find("WJets") == 0 )
    {
      dataSet->SetTitle("W#rightarrowl#nu");
      color = kGreen-3;
    }
    if( dataSet->Name().find("ZJets") == 0 )
    {
      dataSet->SetTitle("Z/#gamma*#rightarrowl^{+}l^{-}");
      color = kAzure-2;
      dataSet->SetEquivalentLuminosity(11229.3432);
    }
    if( dataSet->Name().find("ST") == 0 || dataSet->Name().find("SingleTop") ==0 )
      color = kMagenta;
//    if( dataSet->Name().find("ST_tChannel") == 0 ) color = kMagenta+2;
    Dataset* tmpDS = new Dataset(dataSet->Name(), dataSet->Title(), dataSet->DoIt(), color, dataSet->LineStyle(), dataSet->LineWidth(), dataSet->NormFactor(), dataSet->Xsection());
    tmpDS->SetEquivalentLuminosity( dataSet->EquivalentLumi() );
    if( ! ( dataSet->Name().find("Data") == 0 ) || addedData == false )
    {
      if(dataSet->Name().find("Data") == 0) addedData = true;
      datasetsSemiMu.push_back( tmpDS );
      datasetsSemiEl.push_back( tmpDS );
      dataSetsOneData.push_back( tmpDS );
    }
    dataSetsAll.push_back( tmpDS );
  }
  
  //MSPlots
  MSPlot["nPV_beforechi2_mu"] = new MultiSamplePlot(datasetsSemiMu,"nPV_beforechi2_mu",40,-0.5,39.5,"Number of PV, before chi2 cuts","Nr. of events");
  MSPlot["nPV_beforechi2_el"] = new MultiSamplePlot(datasetsSemiEl,"nPV_beforechi2_el",40,-0.5,39.5,"Number of PV, before chi2 cuts","Nr. of events");
  MSPlot["nPV_mu"] = new MultiSamplePlot(datasetsSemiMu,"nPV_mu",40,-0.5,39.5,"Number of PV","Nr. of events");
  MSPlot["nPV_el"] = new MultiSamplePlot(datasetsSemiEl,"nPV_el",40,-0.5,39.5,"Number of PV","Nr. of events");
  MSPlot["RelPFISo_mu"] = new MultiSamplePlot(datasetsSemiMu,"RelPFIso_mu",100,0,1,"RelPFIso","Nr. of events");
  MSPlot["RelPFISo_el"] = new MultiSamplePlot(datasetsSemiMu,"RelPFIso_el",100,0,1,"RelPFIso","Nr. of events");
  
  MSPlot["M3"] = new MultiSamplePlot(datasetsSemiMu,"M3",160,0,800,"M3","Nr. of events");
  MSPlot["M3_PtTop300"] = new MultiSamplePlot(datasetsSemiMu,"M3_PtTop300",80,0,800,"M3_PtTop300","Nr. of events");  
  MSPlot["M3_PtTTbar200"] = new MultiSamplePlot(datasetsSemiMu,"M3_PtTTbar200",80,0,800,"M3_PtTTbar200","Nr. of events");  
  MSPlot["PtTop"] = new MultiSamplePlot(datasetsSemiMu,"PtTop",150,0,450,"PtTop","Nr. of events");
  MSPlot["PtTTbar"] = new MultiSamplePlot(datasetsSemiMu,"PtTTbar",150,0,450,"PtTTbar","Nr. of events");
  
  MSPlot["AllJets_pt_leptonPlus"] = new MultiSamplePlot(dataSetsOneData,"AllJets_pt_leptonPlus",50,30,530,"Jet p_{T} (GeV)","Jets / 10 GeV");
  MSPlot["AllJets_pt_leptonPlus"]->setMaxY(70000.);
  MSPlot["AllJets_pt_leptonMinus"] = new MultiSamplePlot(dataSetsOneData,"AllJets_pt_leptonMinus",50,30,530,"Jet p_{T} (GeV)","Jets / 10 GeV");
  MSPlot["AllJets_pt_leptonMinus"]->setMaxY(70000.);
  MSPlot["Leading_jet_pt_leptonPlus"] = new MultiSamplePlot(dataSetsOneData,"Leading_jet_pt_leptonPlus",50,20,520,"Leading Jet p_{T} (GeV)","Nr. of events / 10 GeV");
  MSPlot["Leading_jet_pt_leptonMinus"] = new MultiSamplePlot(dataSetsOneData,"Leading_jet_pt_leptonMinus",50,20,520,"Leading Jet p_{T} (GeV)","Nr. of events / 10 GeV");
  MSPlot["4th_jet_pt_leptonPlus"] = new MultiSamplePlot(dataSetsOneData,"4th_jet_pt_leptonPlus",50,20,120,"4th Jet p_{T} (GeV)","Nr. of events / 2 GeV");
  MSPlot["4th_jet_pt_leptonMinus"] = new MultiSamplePlot(dataSetsOneData,"4th_jet_pt_leptonMinus",50,20,120,"4th Jet p_{T} (GeV)","Nr. of events / 2 GeV");
  MSPlot["MET_leptonPlus"] = new MultiSamplePlot(dataSetsOneData,"MET_leptonPlus",100,0,250,"MET (GeV)","Nr. of events");
  MSPlot["MET_leptonMinus"] = new MultiSamplePlot(dataSetsOneData,"MET_leptonMinus",100,0,250,"MET (GeV)","Nr. of events");
  MSPlot["Nr_of_jets"] = new MultiSamplePlot(datasetsSemiEl,"Nr_of_jets",5,3.5,8.5,"Nr. of jets","Nr. of events");
  MSPlot["Nr_of_jets_leptonPlus"] = new MultiSamplePlot(dataSetsOneData,"Nr_of_jets_leptonPlus",5,3.5,8.5,"No. of jets","Events");
  MSPlot["Nr_of_jets_leptonPlus"]->setMaxY(50000.);
  MSPlot["Nr_of_jets_leptonMinus"] = new MultiSamplePlot(dataSetsOneData,"Nr_of_jets_leptonMinus",5,3.5,8.5,"No. of jets","Events");
  MSPlot["Nr_of_jets_leptonMinus"]->setMaxY(50000.);
  MSPlot["MinChi2ndf_Fit_leptonPlus"] = new MultiSamplePlot(dataSetsOneData,"MinChi2ndf_Fit_leptonPlus",50,0,10,"Min #chi^{2}","Events");
  MSPlot["MinChi2ndf_Fit_leptonPlus"]->setMaxY(18000.);
  MSPlot["MinChi2ndf_Fit_leptonMinus"] = new MultiSamplePlot(dataSetsOneData,"MinChi2ndf_Fit_leptonMinus",50,0,10,"Min #chi^{2}","Events");
  MSPlot["MinChi2ndf_Fit_leptonMinus"]->setMaxY(28000.);
  MSPlot["mTop_Fit_leptonPlus"] = new MultiSamplePlot(dataSetsOneData,"mTop_Fit_leptonPlus",50,0,1000,"Fitted Top Mass (GeV)","Events / 20 GeV");
  MSPlot["mTop_Fit_leptonPlus"]->setMaxY(6000.);
  MSPlot["mTop_Fit_leptonMinus"] = new MultiSamplePlot(dataSetsOneData,"mTop_Fit_leptonMinus",50,0,1000,"Fitted Top Mass (GeV)","Events / 20 GeV");
  MSPlot["mTop_Fit_leptonMinus"]->setMaxY(6000.);
  
  MSPlot["mTop_Fit_AllCombi_muPlus"] = new MultiSamplePlot(datasetsSemiMu,"mTop_Fit_AllCombi_muPlus",100,0,1000,"Fitted Top Mass (GeV)","Nr. of events / 10 GeV");
  MSPlot["mTop_Fit_AllCombi_muMinus"] = new MultiSamplePlot(datasetsSemiMu,"mTop_Fit_AllCombi_muMinus",100,0,1000,"Fitted Top Mass (GeV)","Nr. of events / 10 GeV");
  MSPlot["mTop_Fit_AllCombi_elPlus"] = new MultiSamplePlot(datasetsSemiEl,"mTop_Fit_AllCombi_elPlus",100,0,1000,"Fitted Top Mass (GeV)","Nr. of events / 10 GeV");
  MSPlot["mTop_Fit_AllCombi_elMinus"] = new MultiSamplePlot(datasetsSemiEl,"mTop_Fit_AllCombi_elMinus",100,0,1000,"Fitted Top Mass (GeV)","Nr. of events / 10 GeV");
  MSPlot["sigmaMtop_Fit_AllCombi_muPlus"] = new MultiSamplePlot(datasetsSemiMu,"sigmaMtop_Fit_AllCombi_muPlus",120,0,60,"Fitted Top Mass Sigma (GeV)","Nr. of events / 0.5 GeV");
  MSPlot["sigmaMtop_Fit_AllCombi_muMinus"] = new MultiSamplePlot(datasetsSemiMu,"sigmaMtop_Fit_AllCombi_muMinus",120,0,60,"Fitted Top Mass Sigma (GeV)","Nr. of events / 0.5 GeV");
  MSPlot["sigmaMtop_Fit_AllCombi_elPlus"] = new MultiSamplePlot(datasetsSemiEl,"sigmaMtop_Fit_AllCombi_elPlus",120,0,60,"Fitted Top Mass Sigma (GeV)","Nr. of events / 0.5 GeV");
  MSPlot["sigmaMtop_Fit_AllCombi_elMinus"] = new MultiSamplePlot(datasetsSemiEl,"sigmaMtop_Fit_AllCombi_elMinus",120,0,60,"Fitted Top Mass Sigma (GeV)","Nr. of events / 0.5 GeV");
  MSPlot["MinChi2ndf_Fit_AllCombi_muPlus"] = new MultiSamplePlot(datasetsSemiMu,"MinChi2ndf_Fit_AllCombi_muPlus",80,0,20,"Min #chi^{2}/ndf","Nr. of events");
  MSPlot["MinChi2ndf_Fit_AllCombi_muMinus"] = new MultiSamplePlot(datasetsSemiMu,"MinChi2ndf_Fit_AllCombi_muMinus",80,0,20,"Min #chi^{2}/ndf","Nr. of events");
  MSPlot["MinChi2ndf_Fit_AllCombi_elPlus"] = new MultiSamplePlot(datasetsSemiEl,"MinChi2ndf_Fit_AllCombi_elPlus",80,0,20,"Min #chi^{2}/ndf","Nr. of events");
  MSPlot["MinChi2ndf_Fit_AllCombi_elMinus"] = new MultiSamplePlot(datasetsSemiEl,"MinChi2ndf_Fit_AllCombi_elMinus",80,0,20,"Min #chi^{2}/ndf","Nr. of events");
  
  if(doAllMSPlots)
  {
    MSPlot["Lepton_pt_muPlus"] = new MultiSamplePlot(datasetsSemiMu,"Lepton_pt_muPlus",50,0,250,"Positive Muon p_{T} (GeV)","Nr. of events / 5 GeV");
    MSPlot["Lepton_pt_muMinus"] = new MultiSamplePlot(datasetsSemiMu,"Lepton_pt_muMinus",50,0,250,"Negative Muon p_{T} (GeV)","Nr. of events / 5 GeV");
    MSPlot["Lepton_pt_elPlus"] = new MultiSamplePlot(datasetsSemiEl,"Lepton_pt_elPlus",50,0,250,"Positive Electron p_{T} (GeV)","Nr. of events / 5 GeV");
    MSPlot["Lepton_pt_elMinus"] = new MultiSamplePlot(datasetsSemiEl,"Lepton_pt_elMinus",50,0,250,"Negative Electron p_{T} (GeV)","Nr. of events / 5 GeV");
    MSPlot["Lepton_et_muPlus"] = new MultiSamplePlot(datasetsSemiMu,"Lepton_et_muPlus",50,0,250,"Positive Muon E_{T} (GeV)","Nr. of events / 5 GeV");
    MSPlot["Lepton_et_muMinus"] = new MultiSamplePlot(datasetsSemiMu,"Lepton_et_muMinus",50,0,250,"Negative Muon E_{T} (GeV)","Nr. of events / 5 GeV");
    MSPlot["Lepton_et_elPlus"] = new MultiSamplePlot(datasetsSemiEl,"Lepton_et_elPlus",50,0,250,"Positive Electron E_{T} (GeV)","Nr. of events / 5 GeV");
    MSPlot["Lepton_et_elMinus"] = new MultiSamplePlot(datasetsSemiEl,"Lepton_et_elMinus",50,0,250,"Negative Electron E_{T} (GeV)","Nr. of events / 5 GeV");
    MSPlot["Lepton_eta_muPlus"] = new MultiSamplePlot(datasetsSemiMu,"Lepton_eta_muPlus",50,-2.5,2.5,"Positive Muon #eta","Nr. of events");
    MSPlot["Lepton_eta_muMinus"] = new MultiSamplePlot(datasetsSemiMu,"Lepton_eta_muMinus",50,-2.5,2.5,"Negative Muon #eta","Nr. of events");
    MSPlot["Lepton_eta_elPlus"] = new MultiSamplePlot(datasetsSemiEl,"Lepton_eta_elPlus",50,-2.5,2.5,"Positive Electron #eta","Nr. of events");
    MSPlot["Lepton_eta_elMinus"] = new MultiSamplePlot(datasetsSemiEl,"Lepton_eta_elMinus",50,-2.5,2.5,"Negative Electron #eta","Nr. of events");
    MSPlot["Lepton_phi_muPlus"] = new MultiSamplePlot(datasetsSemiMu,"Lepton_phi_muPlus",35,-3.5,3.5,"Positive Muon #phi","Nr. of events");
    MSPlot["Lepton_phi_muMinus"] = new MultiSamplePlot(datasetsSemiMu,"Lepton_phi_muMinus",35,-3.5,3.5,"Negative Muon #phi","Nr. of events");
    MSPlot["Lepton_phi_elPlus"] = new MultiSamplePlot(datasetsSemiEl,"Lepton_phi_elPlus",35,-3.5,3.5,"Positive Electron #phi","Nr. of events");
    MSPlot["Lepton_phi_elMinus"] = new MultiSamplePlot(datasetsSemiEl,"Lepton_phi_elMinus",35,-3.5,3.5,"Negative Electron #phi","Nr. of events");
    MSPlot["AllJets_pt_muPlus"] = new MultiSamplePlot(datasetsSemiMu,"AllJets_pt_muPlus",60,20,320,"All jets p_{T} (GeV)","Nr. of jets / 5 GeV");
    MSPlot["AllJets_pt_muMinus"] = new MultiSamplePlot(datasetsSemiMu,"AllJets_pt_muMinus",60,20,320,"All jets p_{T} (GeV)","Nr. of jets / 5 GeV");
    MSPlot["AllJets_pt_elPlus"] = new MultiSamplePlot(datasetsSemiEl,"AllJets_pt_elPlus",60,20,320,"All jets p_{T} (GeV)","Nr. of jets / 5 GeV");
    MSPlot["AllJets_pt_elMinus"] = new MultiSamplePlot(datasetsSemiEl,"AllJets_pt_elMinus",60,20,320,"All jets p_{T} (GeV)","Nr. of jets / 5 GeV");
    MSPlot["AllJets_eta_muPlus"] = new MultiSamplePlot(datasetsSemiMu,"AllJets_eta_muPlus",50,-2.5,2.5,"All jets #eta","Nr. of jets");
    MSPlot["AllJets_eta_muMinus"] = new MultiSamplePlot(datasetsSemiMu,"AllJets_eta_muMinus",50,-2.5,2.5,"All jets #eta","Nr. of jets");
    MSPlot["AllJets_eta_elPlus"] = new MultiSamplePlot(datasetsSemiEl,"AllJets_eta_elPlus",50,-2.5,2.5,"All jets #eta","Nr. of jets");
    MSPlot["AllJets_eta_elMinus"] = new MultiSamplePlot(datasetsSemiEl,"AllJets_eta_elMinus",50,-2.5,2.5,"All jets #eta","Nr. of jets");
    MSPlot["AllJets_phi_muPlus"] = new MultiSamplePlot(datasetsSemiMu,"AllJets_phi_muPlus",35,-3.5,3.5,"All jets #phi","Nr. of jets");
    MSPlot["AllJets_phi_muMinus"] = new MultiSamplePlot(datasetsSemiMu,"AllJets_phi_muMinus",35,-3.5,3.5,"All jets #phi","Nr. of jets");
    MSPlot["AllJets_phi_elPlus"] = new MultiSamplePlot(datasetsSemiEl,"AllJets_phi_elPlus",35,-3.5,3.5,"All jets #phi","Nr. of jets");
    MSPlot["AllJets_phi_elMinus"] = new MultiSamplePlot(datasetsSemiEl,"AllJets_phi_elMinus",35,-3.5,3.5,"All jets #phi","Nr. of jets");
	  MSPlot["Leading_jet_pt_muPlus"] = new MultiSamplePlot(datasetsSemiMu,"Leading_jet_pt_muPlus",50,20,520,"Leading Jet p_{T} (GeV)","Nr. of events / 10 GeV");
	  MSPlot["Leading_jet_pt_muMinus"] = new MultiSamplePlot(datasetsSemiMu,"Leading_jet_pt_muMinus",50,20,520,"Leading Jet p_{T} (GeV)","Nr. of events / 10 GeV");
	  MSPlot["Leading_jet_pt_elPlus"] = new MultiSamplePlot(datasetsSemiEl,"Leading_jet_pt_elPlus",50,20,520,"Leading Jet p_{T} (GeV)","Nr. of events / 10 GeV");
	  MSPlot["Leading_jet_pt_elMinus"] = new MultiSamplePlot(datasetsSemiEl,"Leading_jet_pt_elMinus",50,20,520,"Leading Jet p_{T} (GeV)","Nr. of events / 10 GeV");
	  MSPlot["Leading_jet_eta_muPlus"] = new MultiSamplePlot(datasetsSemiMu,"Leading_jet_eta_muPlus",50,-2.5,2.5,"Leading Jet #eta","Nr. of events");
	  MSPlot["Leading_jet_eta_muMinus"] = new MultiSamplePlot(datasetsSemiMu,"Leading_jet_eta_muMinus",50,-2.5,2.5,"Leading Jet #eta","Nr. of events");
	  MSPlot["Leading_jet_eta_elPlus"] = new MultiSamplePlot(datasetsSemiEl,"Leading_jet_eta_elPlus",50,-2.5,2.5,"Leading Jet #eta","Nr. of events");
	  MSPlot["Leading_jet_eta_elMinus"] = new MultiSamplePlot(datasetsSemiEl,"Leading_jet_eta_elMinus",50,-2.5,2.5,"Leading Jet #eta","Nr. of events");
	  MSPlot["Leading_jet_phi_muPlus"] = new MultiSamplePlot(datasetsSemiMu,"Leading_jet_phi_muPlus",35,-3.5,3.5,"Leading Jet #phi","Nr. of events");
	  MSPlot["Leading_jet_phi_muMinus"] = new MultiSamplePlot(datasetsSemiMu,"Leading_jet_phi_muMinus",35,-3.5,3.5,"Leading Jet #phi","Nr. of events");
	  MSPlot["Leading_jet_phi_elPlus"] = new MultiSamplePlot(datasetsSemiEl,"Leading_jet_phi_elPlus",35,-3.5,3.5,"Leading Jet #phi","Nr. of events");
	  MSPlot["Leading_jet_phi_elMinus"] = new MultiSamplePlot(datasetsSemiEl,"Leading_jet_phi_elMinus",35,-3.5,3.5,"Leading Jet #phi","Nr. of events");
	  MSPlot["2nd_jet_pt_muPlus"] = new MultiSamplePlot(datasetsSemiMu,"2nd_jet_pt_muPlus",60,20,320,"2nd Jet p_{T} (GeV)","Nr. of events / 5 GeV");
	  MSPlot["2nd_jet_pt_muMinus"] = new MultiSamplePlot(datasetsSemiMu,"2nd_jet_pt_muMinus",60,20,320,"2nd Jet p_{T} (GeV)","Nr. of events / 5 GeV");
	  MSPlot["2nd_jet_pt_elPlus"] = new MultiSamplePlot(datasetsSemiEl,"2nd_jet_pt_elPlus",60,20,320,"2nd Jet p_{T} (GeV)","Nr. of events / 5 GeV");
	  MSPlot["2nd_jet_pt_elMinus"] = new MultiSamplePlot(datasetsSemiEl,"2nd_jet_pt_elMinus",60,20,320,"2nd Jet p_{T} (GeV)","Nr. of events / 5 GeV");
	  MSPlot["2nd_jet_eta_muPlus"] = new MultiSamplePlot(datasetsSemiMu,"2nd_jet_eta_muPlus",50,-2.5,2.5,"2nd Jet #eta","Nr. of events");
	  MSPlot["2nd_jet_eta_muMinus"] = new MultiSamplePlot(datasetsSemiMu,"2nd_jet_eta_muMinus",50,-2.5,2.5,"2nd Jet #eta","Nr. of events");
	  MSPlot["2nd_jet_eta_elPlus"] = new MultiSamplePlot(datasetsSemiEl,"2nd_jet_eta_elPlus",50,-2.5,2.5,"2nd Jet #eta","Nr. of events");
	  MSPlot["2nd_jet_eta_elMinus"] = new MultiSamplePlot(datasetsSemiEl,"2nd_jet_eta_elMinus",50,-2.5,2.5,"2nd Jet #eta","Nr. of events");
	  MSPlot["2nd_jet_phi_muPlus"] = new MultiSamplePlot(datasetsSemiMu,"2nd_jet_phi_muPlus",35,-3.5,3.5,"2nd Jet #phi","Nr. of events");
	  MSPlot["2nd_jet_phi_muMinus"] = new MultiSamplePlot(datasetsSemiMu,"2nd_jet_phi_muMinus",35,-3.5,3.5,"2nd Jet #phi","Nr. of events");
	  MSPlot["2nd_jet_phi_elPlus"] = new MultiSamplePlot(datasetsSemiEl,"2nd_jet_phi_elPlus",35,-3.5,3.5,"2nd Jet #phi","Nr. of events");
	  MSPlot["2nd_jet_phi_elMinus"] = new MultiSamplePlot(datasetsSemiEl,"2nd_jet_phi_elMinus",35,-3.5,3.5,"2nd Jet #phi","Nr. of events");
	  MSPlot["3rd_jet_pt_muPlus"] = new MultiSamplePlot(datasetsSemiMu,"3rd_jet_pt_muPlus",50,20,220,"3rd Jet p_{T} (GeV)","Nr. of events / 4 GeV");
	  MSPlot["3rd_jet_pt_muMinus"] = new MultiSamplePlot(datasetsSemiMu,"3rd_jet_pt_muMinus",50,20,220,"3rd Jet p_{T} (GeV)","Nr. of events / 4 GeV");
	  MSPlot["3rd_jet_pt_elPlus"] = new MultiSamplePlot(datasetsSemiEl,"3rd_jet_pt_elPlus",50,20,220,"3rd Jet p_{T} (GeV)","Nr. of events / 4 GeV");
	  MSPlot["3rd_jet_pt_elMinus"] = new MultiSamplePlot(datasetsSemiEl,"3rd_jet_pt_elMinus",50,20,220,"3rd Jet p_{T} (GeV)","Nr. of events / 4 GeV");
	  MSPlot["3rd_jet_eta_muPlus"] = new MultiSamplePlot(datasetsSemiMu,"3rd_jet_eta_muPlus",50,-2.5,2.5,"3rd Jet #eta","Nr. of events");
	  MSPlot["3rd_jet_eta_muMinus"] = new MultiSamplePlot(datasetsSemiMu,"3rd_jet_eta_muMinus",50,-2.5,2.5,"3rd Jet #eta","Nr. of events");
	  MSPlot["3rd_jet_eta_elPlus"] = new MultiSamplePlot(datasetsSemiEl,"3rd_jet_eta_elPlus",50,-2.5,2.5,"3rd Jet #eta","Nr. of events");
	  MSPlot["3rd_jet_eta_elMinus"] = new MultiSamplePlot(datasetsSemiEl,"3rd_jet_eta_elMinus",50,-2.5,2.5,"3rd Jet #eta","Nr. of events");
	  MSPlot["3rd_jet_phi_muPlus"] = new MultiSamplePlot(datasetsSemiMu,"3rd_jet_phi_muPlus",35,-3.5,3.5,"3rd Jet #phi","Nr. of events");
	  MSPlot["3rd_jet_phi_muMinus"] = new MultiSamplePlot(datasetsSemiMu,"3rd_jet_phi_muMinus",35,-3.5,3.5,"3rd Jet #phi","Nr. of events");
	  MSPlot["3rd_jet_phi_elPlus"] = new MultiSamplePlot(datasetsSemiEl,"3rd_jet_phi_elPlus",35,-3.5,3.5,"3rd Jet #phi","Nr. of events");
	  MSPlot["3rd_jet_phi_elMinus"] = new MultiSamplePlot(datasetsSemiEl,"3rd_jet_phi_elMinus",35,-3.5,3.5,"3rd Jet #phi","Nr. of events");
	  MSPlot["4th_jet_pt_muPlus"] = new MultiSamplePlot(datasetsSemiMu,"4th_jet_pt_muPlus",50,20,120,"4th Jet p_{T} (GeV)","Nr. of events / 2 GeV");
	  MSPlot["4th_jet_pt_muMinus"] = new MultiSamplePlot(datasetsSemiMu,"4th_jet_pt_muMinus",50,20,120,"4th Jet p_{T} (GeV)","Nr. of events / 2 GeV");
	  MSPlot["4th_jet_pt_elPlus"] = new MultiSamplePlot(datasetsSemiEl,"4th_jet_pt_elPlus",50,20,120,"4th Jet p_{T} (GeV)","Nr. of events / 2 GeV");
	  MSPlot["4th_jet_pt_elMinus"] = new MultiSamplePlot(datasetsSemiEl,"4th_jet_pt_elMinus",50,20,120,"4th Jet p_{T} (GeV)","Nr. of events / 2 GeV");
	  MSPlot["4th_jet_eta_muPlus"] = new MultiSamplePlot(datasetsSemiMu,"4th_jet_eta_muPlus",50,-2.5,2.5,"4th Jet #eta","Nr. of events");
	  MSPlot["4th_jet_eta_muMinus"] = new MultiSamplePlot(datasetsSemiMu,"4th_jet_eta_muMinus",50,-2.5,2.5,"4th Jet #eta","Nr. of events");
	  MSPlot["4th_jet_eta_elPlus"] = new MultiSamplePlot(datasetsSemiEl,"4th_jet_eta_elPlus",50,-2.5,2.5,"4th Jet #eta","Nr. of events");
	  MSPlot["4th_jet_eta_elMinus"] = new MultiSamplePlot(datasetsSemiEl,"4th_jet_eta_elMinus",50,-2.5,2.5,"4th Jet #eta","Nr. of events");
	  MSPlot["4th_jet_phi_muPlus"] = new MultiSamplePlot(datasetsSemiMu,"4th_jet_phi_muPlus",35,-3.5,3.5,"4th Jet #phi","Nr. of events");
    MSPlot["4th_jet_phi_muMinus"] = new MultiSamplePlot(datasetsSemiMu,"4th_jet_phi_muMinus",35,-3.5,3.5,"4th Jet #phi","Nr. of events");
   	MSPlot["4th_jet_phi_elPlus"] = new MultiSamplePlot(datasetsSemiEl,"4th_jet_phi_elPlus",35,-3.5,3.5,"4th Jet #phi","Nr. of events");
    MSPlot["4th_jet_phi_elMinus"] = new MultiSamplePlot(datasetsSemiEl,"4th_jet_phi_elMinus",35,-3.5,3.5,"4th Jet #phi","Nr. of events");
   
    MSPlot["Nr_of_jets_muPlus"] = new MultiSamplePlot(datasetsSemiMu,"Nr_of_jets_muPlus",5,3.5,8.5,"Nr. of jets","Nr. of events");
    MSPlot["Nr_of_jets_muMinus"] = new MultiSamplePlot(datasetsSemiMu,"Nr_of_jets_muMinus",5,3.5,8.5,"Nr. of jets","Nr. of events");
    MSPlot["Nr_of_jets_elPlus"] = new MultiSamplePlot(datasetsSemiEl,"Nr_of_jets_elPlus",5,3.5,8.5,"Nr. of jets","Nr. of events");
    MSPlot["Nr_of_jets_elMinus"] = new MultiSamplePlot(datasetsSemiEl,"Nr_of_jets_elMinus",5,3.5,8.5,"Nr. of jets","Nr. of events");
    MSPlot["nBtags_muPlus"] = new MultiSamplePlot(datasetsSemiMu,"nBtags_muPlus",6,-0.5,5.5,"Nr. of b-tagged jets","Nr. of events");
    MSPlot["nBtags_muMinus"] = new MultiSamplePlot(datasetsSemiMu,"nBtags_muMinus",6,-0.5,5.5,"Nr. of b-tagged jets","Nr. of events");
    MSPlot["nBtags_elPlus"] = new MultiSamplePlot(datasetsSemiEl,"nBtags_elPlus",6,-0.5,5.5,"Nr. of b-tagged jets","Nr. of events");
    MSPlot["nBtags_elMinus"] = new MultiSamplePlot(datasetsSemiEl,"nBtags_elMinus",6,-0.5,5.5,"Nr. of b-tagged jets","Nr. of events");
    MSPlot["M3_muPlus"] = new MultiSamplePlot(datasetsSemiMu,"M3_muPlus",80,0,800,"M3 (GeV)","Nr. of events / 10 GeV");
    MSPlot["M3_muMinus"] = new MultiSamplePlot(datasetsSemiMu,"M3_muMinus",80,0,800,"M3 (GeV)","Nr. of events / 10 GeV");
    MSPlot["M3_elPlus"] = new MultiSamplePlot(datasetsSemiEl,"M3_elPlus",80,0,800,"M3 (GeV)","Nr. of events / 10 GeV");
    MSPlot["M3_elMinus"] = new MultiSamplePlot(datasetsSemiEl,"M3_elMinus",80,0,800,"M3 (GeV)","Nr. of events / 10 GeV");
    MSPlot["Angle_1stjet_2ndjet_muPlus"] = new MultiSamplePlot(datasetsSemiMu,"Angle_1stjet_2ndjet_muPlus",32,0,3.2,"Angle(leading jet, 2nd jet)","Nr. of events");
    MSPlot["Angle_1stjet_2ndjet_muMinus"] = new MultiSamplePlot(datasetsSemiMu,"Angle_1stjet_2ndjet_muMinus",32,0,3.2,"Angle(leading jet, 2nd jet)","Nr. of events");
    MSPlot["Angle_1stjet_2ndjet_elPlus"] = new MultiSamplePlot(datasetsSemiEl,"Angle_1stjet_2ndjet_elPlus",32,0,3.2,"Angle(leading jet, 2nd jet)","Nr. of events");
    MSPlot["Angle_1stjet_2ndjet_elMinus"] = new MultiSamplePlot(datasetsSemiEl,"Angle_1stjet_2ndjet_elMinus",32,0,3.2,"Angle(leading jet, 2nd jet)","Nr. of events");
    MSPlot["Angle_1stjet_3rdjet_muPlus"] = new MultiSamplePlot(datasetsSemiMu,"Angle_1stjet_3rdjet_muPlus",32,0,3.2,"Angle(leading jet, 3rd jet)","Nr. of events");
    MSPlot["Angle_1stjet_3rdjet_muMinus"] = new MultiSamplePlot(datasetsSemiMu,"Angle_1stjet_3rdjet_muMinus",32,0,3.2,"Angle(leading jet, 3rd jet)","Nr. of events");
    MSPlot["Angle_1stjet_3rdjet_elPlus"] = new MultiSamplePlot(datasetsSemiEl,"Angle_1stjet_3rdjet_elPlus",32,0,3.2,"Angle(leading jet, 3rd jet)","Nr. of events");
    MSPlot["Angle_1stjet_3rdjet_elMinus"] = new MultiSamplePlot(datasetsSemiEl,"Angle_1stjet_3rdjet_elMinus",32,0,3.2,"Angle(leading jet, 3rd jet)","Nr. of events");
    MSPlot["MET_muPlus"] = new MultiSamplePlot(datasetsSemiMu,"MET_muPlus",100,0,250,"MET (GeV)","Nr. of events");
    MSPlot["MET_muMinus"] = new MultiSamplePlot(datasetsSemiMu,"MET_muMinus",100,0,250,"MET (GeV)","Nr. of events");
    MSPlot["MET_elPlus"] = new MultiSamplePlot(datasetsSemiEl,"MET_elPlus",100,0,250,"MET (GeV)","Nr. of events");
    MSPlot["MET_elMinus"] = new MultiSamplePlot(datasetsSemiEl,"MET_elMinus",100,0,250,"MET (GeV)","Nr. of events");
    
    MSPlot["mW_UnFitted_muPlus"] = new MultiSamplePlot(datasetsSemiMu,"mW_UnFitted_muPlus",80,0,200,"Unfitted W Mass (GeV)","Nr. of events / 2.5 GeV");
    MSPlot["mW_UnFitted_muMinus"] = new MultiSamplePlot(datasetsSemiMu,"mW_UnFitted_muMinus",80,0,200,"Unfitted W Mass (GeV)","Nr. of events / 2.5 GeV");
    MSPlot["mW_UnFitted_elPlus"] = new MultiSamplePlot(datasetsSemiEl,"mW_UnFitted_elPlus",80,0,200,"Unfitted W Mass (GeV)","Nr. of events / 2.5 GeV");
    MSPlot["mW_UnFitted_elMinus"] = new MultiSamplePlot(datasetsSemiEl,"mW_UnFitted_elMinus",80,0,200,"Unfitted W Mass (GeV)","Nr. of events / 2.5 GeV");
    MSPlot["mTop_UnFitted_muPlus"] = new MultiSamplePlot(datasetsSemiMu,"mTop_UnFitted_muPlus",100,0,1000,"Unfitted Top Mass (GeV)","Nr. of events / 10 GeV");
    MSPlot["mTop_UnFitted_muMinus"] = new MultiSamplePlot(datasetsSemiMu,"mTop_UnFitted_muMinus",100,0,1000,"Unfitted Top Mass (GeV)","Nr. of events / 10 GeV");
    MSPlot["mTop_UnFitted_elPlus"] = new MultiSamplePlot(datasetsSemiEl,"mTop_UnFitted_elPlus",100,0,1000,"Unfitted Top Mass (GeV)","Nr. of events / 10 GeV");
    MSPlot["mTop_UnFitted_elMinus"] = new MultiSamplePlot(datasetsSemiEl,"mTop_UnFitted_elMinus",100,0,1000,"Unfitted Top Mass (GeV)","Nr. of events / 10 GeV");
    MSPlot["mTop_Fit_muPlus"] = new MultiSamplePlot(datasetsSemiMu,"mTop_Fit_muPlus",100,0,1000,"Fitted Top Mass (GeV)","Nr. of events / 10 GeV");
    MSPlot["mTop_Fit_muMinus"] = new MultiSamplePlot(datasetsSemiMu,"mTop_Fit_muMinus",100,0,1000,"Fitted Top Mass (GeV)","Nr. of events / 10 GeV");
    MSPlot["mTop_Fit_elPlus"] = new MultiSamplePlot(datasetsSemiEl,"mTop_Fit_elPlus",100,0,1000,"Fitted Top Mass (GeV)","Nr. of events / 10 GeV");
    MSPlot["mTop_Fit_elMinus"] = new MultiSamplePlot(datasetsSemiEl,"mTop_Fit_elMinus",100,0,1000,"Fitted Top Mass (GeV)","Nr. of events / 10 GeV");
    MSPlot["sigmaMtop_Fit_muPlus"] = new MultiSamplePlot(datasetsSemiMu,"sigmaMtop_Fit_muPlus",120,0,60,"Fitted Top Mass Sigma (GeV)","Nr. of events / 0.5 GeV");
    MSPlot["sigmaMtop_Fit_muMinus"] = new MultiSamplePlot(datasetsSemiMu,"sigmaMtop_Fit_muMinus",120,0,60,"Fitted Top Mass Sigma (GeV)","Nr. of events / 0.5 GeV");
    MSPlot["sigmaMtop_Fit_elPlus"] = new MultiSamplePlot(datasetsSemiEl,"sigmaMtop_Fit_elPlus",120,0,60,"Fitted Top Mass Sigma (GeV)","Nr. of events / 0.5 GeV");
    MSPlot["sigmaMtop_Fit_elMinus"] = new MultiSamplePlot(datasetsSemiEl,"sigmaMtop_Fit_elMinus",120,0,60,"Fitted Top Mass Sigma (GeV)","Nr. of events / 0.5 GeV");
    MSPlot["MinChi2ndf_Fit_muPlus"] = new MultiSamplePlot(datasetsSemiMu,"MinChi2ndf_Fit_muPlus",80,0,20,"Min #chi^{2}/ndf","Nr. of events");
    MSPlot["MinChi2ndf_Fit_muMinus"] = new MultiSamplePlot(datasetsSemiMu,"MinChi2ndf_Fit_muMinus",80,0,20,"Min #chi^{2}/ndf","Nr. of events");
    MSPlot["MinChi2ndf_Fit_elPlus"] = new MultiSamplePlot(datasetsSemiEl,"MinChi2ndf_Fit_elPlus",80,0,20,"Min #chi^{2}/ndf","Nr. of events");
    MSPlot["MinChi2ndf_Fit_elMinus"] = new MultiSamplePlot(datasetsSemiEl,"MinChi2ndf_Fit_elMinus",80,0,20,"Min #chi^{2}/ndf","Nr. of events");
    MSPlot["Prob_MinChi2_Fit_muPlus"] = new MultiSamplePlot(datasetsSemiMu,"Prob_MinChi2_Fit_muPlus",110,0,1,"Prob(Min #chi^{2})","Nr. of events");
    MSPlot["Prob_MinChi2_Fit_muMinus"] = new MultiSamplePlot(datasetsSemiMu,"Prob_MinChi2_Fit_muMinus",110,0,1,"Prob(Min #chi^{2})","Nr. of events");
    MSPlot["Prob_MinChi2_Fit_elPlus"] = new MultiSamplePlot(datasetsSemiEl,"Prob_MinChi2_Fit_elPlus",110,0,1,"Prob(Min #chi^{2})","Nr. of events");
    MSPlot["Prob_MinChi2_Fit_elMinus"] = new MultiSamplePlot(datasetsSemiEl,"Prob_MinChi2_Fit_elMinus",110,0,1,"Prob(Min #chi^{2})","Nr. of events");
    MSPlot["nUsedCombis_muPlus"] = new MultiSamplePlot(datasetsSemiMu,"nUsedCombis_muPlus",13,-0.5,12.5,"Nr. of used jet combinations","Nr. of events");
    MSPlot["nUsedCombis_muMinus"] = new MultiSamplePlot(datasetsSemiMu,"nUsedCombis_muMinus",13,-0.5,12.5,"Nr. of used jet combinations","Nr. of events");
    MSPlot["nUsedCombis_elPlus"] = new MultiSamplePlot(datasetsSemiEl,"nUsedCombis_elPlus",13,-0.5,12.5,"Nr. of used jet combinations","Nr. of events");
    MSPlot["nUsedCombis_elMinus"] = new MultiSamplePlot(datasetsSemiEl,"nUsedCombis_elMinus",13,-0.5,12.5,"Nr. of used jet combinations","Nr. of events");
  }
  
  //setLumi off the plots (then scaling to lumi is done at the end, for proper error calculation)
  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
  {
    string name = it->first;
    if(name.find("mu") < name.size()) it->second->setDataLumi( LuminosityMu );
    else if(name.find("el") < name.size()) it->second->setDataLumi( LuminosityEl );
    else it->second->setDataLumi( ( LuminosityMu + LuminosityEl ) / 2. );
  }
  
  // initialize LumiReWeighting stuff
  cout << "Initializing LumiReWeighting stuff" << endl;
//  LumiReWeighting LumiWeightsMu = LumiReWeighting("PileUpReweighting/pileup_WJets_36bins.root", "PileUpReweighting/pileup_2011Data_UpToRun173692.root", "pileup2", "pileup");
//  LumiReWeighting LumiWeights = LumiReWeighting("PileUpReweighting/pileup_WJets_36bins.root", "PileUpReweighting/pileup_2011Data_UpToRun173692.root", "pileup2", "pileup");
  LumiReWeighting LumiWeights = LumiReWeighting("PileUpReweighting/pileup_WJets_36bins.root", "PileUpReweighting/pileup_2011Data_UpToRun180252.root", "pileup2", "pileup");
  PoissonMeanShifter PShiftUp = PoissonMeanShifter(0.6); // PU-systematic
  PoissonMeanShifter PShiftDown = PoissonMeanShifter(-0.6); // PU-systematic
  
  // 3D-reweighting
  Lumi3DReWeighting Lumi3DWeights = Lumi3DReWeighting("PileUpReweighting/pileup_MC_Flat10PlusTail.root", "PileUpReweighting/pileup_FineBin_2011Data_UpToRun180252.root", "pileup", "pileup");
//  Lumi3DReWeighting Lumi3DWeights = Lumi3DReWeighting("PileUpReweighting/pileup_MC_Flat10PlusTail.root", "PileUpReweighting/pileup_FineBin_2011Data_UpToRun173692.root", "pileup", "pileup");
  Lumi3DWeights.weight3D_init(1.);
  Lumi3DReWeighting Lumi3DWeightsUp = Lumi3DReWeighting("PileUpReweighting/pileup_MC_Flat10PlusTail.root", "PileUpReweighting/pileup_FineBin_2011Data_UpToRun180252.root", "pileup", "pileup");
  Lumi3DWeightsUp.weight3D_init(1.08);
  Lumi3DReWeighting Lumi3DWeightsDown = Lumi3DReWeighting("PileUpReweighting/pileup_MC_Flat10PlusTail.root", "PileUpReweighting/pileup_FineBin_2011Data_UpToRun180252.root", "pileup", "pileup");
  Lumi3DWeightsDown.weight3D_init(0.92);
  
  cout << " Initialized LumiReWeighting stuff" << endl;
  
  // load L7Corrections
  ResolutionFit* resFitLightJets = new ResolutionFit("LightJet");
  resFitLightJets->LoadResolutions("/home/stijn/JES+MTop/Presentatie/2011/11_07_2011_L7Corrections/Summer11/lightJetReso.root");
  ResolutionFit* resFitBJets = new ResolutionFit("BJet");
  resFitBJets->LoadResolutions("/home/stijn/JES+MTop/Presentatie/2011/11_07_2011_L7Corrections/Summer11/bJetReso.root");
  
  // load ideogram cfg stuff
  string ideogramParameterFilename("../config/ideogramConfig");
  ideogram::MassLikelihoodWriter writer(ideogramParameterFilename);
  Double_t maxChi2 = writer.initialParameter->GetValue("maxChi2",(Double_t) -1.0);
  Double_t bTagCut = writer.initialParameter->GetValue("bTagCut",(Double_t) -1.0);
  Double_t bTagEff = writer.initialParameter->GetValue("bTagEff",(Double_t) -1.0);
  Double_t misTagRate = writer.initialParameter->GetValue("misTagRate",(Double_t) -1.0);
  
  cout << "Loaded settings from ideogram config file:  " << ideogramParameterFilename << endl;
  cout << " --> maxChi2 : " << maxChi2 << endl;
  cout << " --> bTagCut : " << bTagCut << endl;
  cout << " --> bTagEff : " << bTagEff << endl;
  cout << " --> misTagRate : " << misTagRate << endl;
  
  for(unsigned int iDataSet=0; iDataSet<inputMonsters.size(); iDataSet++)
  {
    TFile* inFile = new TFile(inputMonsters[iDataSet].c_str(),"READ");
    
    TTree* inMonstersTree = (TTree*) inFile->Get("MonsterTree");
    TBranch* m_br = (TBranch*) inMonstersTree->GetBranch("TheLightMonster");
    LightMonster* monster = 0;
    m_br->SetAddress(&monster);
    int nEvent = inMonstersTree->GetEntries();
    
    Dataset* dataSet = dataSetsAll[iDataSet];//(Dataset*) tc_dataset->At(0);
    string dataSetName = dataSet->Name();
    cout << "Processing DataSet: " << dataSetName << "  containing " << nEvent << " events" << endl;
    cout << "Cross section = " << dataSet->Xsection() << "  intLumi = " << dataSet->EquivalentLumi() << "  NormFactor = " << dataSet->NormFactor() << endl;
    
    ideogram::MassFitInfoHeader headerEl, headerMu;
    ideogram::MassFitInfoEvent eventEl, eventMu;
    
    fillMassFitInfoHeader(dataSet, headerEl, true, false);
    fillMassFitInfoHeader(dataSet, headerMu, true, false);
    
    // output ascii file stuff
    mkdir("FitResults_ASCII/",0777);
    string outFileNameSemiMu = "FitResults_ASCII/FitResults_" + dataSetName + "_SemiMu.txt";
    string outFileNameSemiEl = "FitResults_ASCII/FitResults_" + dataSetName + "_SemiEl.txt";
    if( dataSet->Title().find("Data_El") == 0 )
    {
      outFileNameSemiMu = "FitResults_ASCII/FitResults_" + dataSetName + "_El_SemiMu.txt";
      outFileNameSemiEl = "FitResults_ASCII/FitResults_" + dataSetName + "_El_SemiEl.txt";
    }
    else if( dataSet->Title().find("Data_Mu") == 0 )
    {
      outFileNameSemiMu = "FitResults_ASCII/FitResults_" + dataSetName + "_Mu_SemiMu.txt";
      outFileNameSemiEl = "FitResults_ASCII/FitResults_" + dataSetName + "_Mu_SemiEl.txt";
    }
    ofstream outFileSemiMu(outFileNameSemiMu.c_str());
    ofstream outFileSemiEl(outFileNameSemiEl.c_str());
    
    outFileSemiMu << "Output of MTopDiff_Analysis.cc" << endl;
    outFileSemiMu << "BeginHeader" << endl;
    outFileSemiMu << "channel: 1" << endl;
    if( dataSetName.find("Data") == 0 || dataSetName.find("DAtaDriven") != string::npos ) outFileSemiMu << "isMC: 0" << endl;
    else outFileSemiMu << "isMC: 1" << endl;
    outFileSemiMu << "First some dataSet info: " << endl;
    outFileSemiMu << "datasetName: " << dataSetName << endl;
    outFileSemiMu << "datasetTitle: " << dataSet->Title() << endl;
    outFileSemiMu << "datasetCrossSection: " << dataSet->Xsection() << endl;
    outFileSemiMu << "datasetIntegratedLuminosity: " << dataSet->EquivalentLumi() << endl;
    outFileSemiMu << "EndHeader" << endl << endl;
    outFileSemiMu << "Start of event-by-event info " << endl << endl;
    
    outFileSemiEl << "Output of MTopDiff_Analysis.cc" << endl;
    outFileSemiEl << "BeginHeader" << endl;
    outFileSemiEl << "channel: 0" << endl;
    if( dataSetName.find("Data") == 0 || dataSetName.find("DAtaDriven") != string::npos ) outFileSemiEl << "isMC: 0" << endl;
    else outFileSemiEl << "isMC: 1" << endl;
    outFileSemiEl << "First some dataSet info: " << endl;
    outFileSemiEl << "datasetName: " << dataSetName << endl;
    outFileSemiEl << "datasetTitle: " << dataSet->Title() << endl;
    outFileSemiEl << "datasetCrossSection: " << dataSet->Xsection() << endl;
    outFileSemiEl << "datasetIntegratedLuminosity: " << dataSet->EquivalentLumi() << endl;
    outFileSemiEl << "EndHeader" << endl << endl;
    outFileSemiEl << "Start of event-by-event info " << endl << endl;
    
    // output ascii file after likelihood calculation
    mkdir("LikelihoodResults_ASCII/",0777);
    string outFileNameLikelihoodSemiMu = "LikelihoodResults_ASCII/LikelihoodResults_" + dataSetName + "_SemiMu.txt";
    string outFileNameLikelihoodSemiEl = "LikelihoodResults_ASCII/LikelihoodResults_" + dataSetName + "_SemiEl.txt";
    if( dataSet->Title().find("Data_El") == 0 )
    {
      outFileNameLikelihoodSemiMu = "LikelihoodResults_ASCII/LikelihoodResults_" + dataSetName + "_El_SemiMu.txt";
      outFileNameLikelihoodSemiEl = "LikelihoodResults_ASCII/LikelihoodResults_" + dataSetName + "_El_SemiEl.txt";
    }
    else if( dataSet->Title().find("Data_Mu") == 0 )
    {
      outFileNameLikelihoodSemiMu = "LikelihoodResults_ASCII/LikelihoodResults_" + dataSetName + "_Mu_SemiMu.txt";
      outFileNameLikelihoodSemiEl = "LikelihoodResults_ASCII/LikelihoodResults_" + dataSetName + "_Mu_SemiEl.txt";
    }
    ofstream outFileLikelihoodSemiMu(outFileNameLikelihoodSemiMu.c_str());
    ofstream outFileLikelihoodSemiEl(outFileNameLikelihoodSemiEl.c_str());
    
    outFileLikelihoodSemiMu << headerMu << endl;
    outFileLikelihoodSemiMu << "Start of event-by-event info " << endl << endl;
    outFileLikelihoodSemiEl << headerEl << endl;
    outFileLikelihoodSemiEl << "Start of event-by-event info " << endl << endl;
    
    int nSemiMu = 0, nSemiEl = 0;
    int nSemiMu_TTSemiMu = 0, nSemiEl_TTSemiEl = 0;
    int nSemiMuInASCII = 0, nSemiElInASCII = 0;
    float nMuPlus = 0, nElPlus = 0, nMuMinus = 0, nElMinus = 0;
    int nRealBjets = 0, nRealBjetsBtag = 0;
    
    for(unsigned int iEvt=0; iEvt<nEvent; iEvt++)
//    for(unsigned int iEvt=0; iEvt<1000; iEvt++)
    {
      inMonstersTree->GetEvent(iEvt);
//      cout << "event: " << iEvt << endl;
      if(iEvt%10000 == 0)
        std::cout<<"Processing the "<<iEvt<<"th event, time = "<< ((double)clock() - start) / CLOCKS_PER_SEC <<endl;
      
//      if( monster->runID() > 174000 ) continue;
      
      // PU reweighting???
//      float avPU = ( (float)monster->nPUBXm1() + (float)monster->nPU() + (float)monster->nPUBXp1() ) / 3.; // average in 3 BX!!!, as recommended
//      float lumiWeight = LumiWeights.ITweight( (float) monster->nPU() );
//      lumiWeight = 1;
      float lumiWeight = 1., lumiWeightUp = 1., lumiWeightDown = 1.;
      if( ! (dataSetName.find("Data") == 0 || dataSetName.find("DAtaDriven") != string::npos) )
      {
        lumiWeight = Lumi3DWeights.weight3D(monster->nPUBXm1(), monster->nPU(), monster->nPUBXp1());
        lumiWeightUp = Lumi3DWeightsUp.weight3D(monster->nPUBXm1(), monster->nPU(), monster->nPUBXp1());
        lumiWeightDown = Lumi3DWeightsDown.weight3D(monster->nPUBXm1(), monster->nPU(), monster->nPUBXp1());
//        lumiWeight = LumiWeights.ITweight( monster->nPU() );
//        lumiWeightUp = lumiWeight * PShiftUp.ShiftWeight( (float) monster->nPU() );
//        lumiWeightDown = lumiWeight * PShiftDown.ShiftWeight( (float) monster->nPU() );
      }
      
      // muoncharge-string
      string leptonCharge, leptonDecay, leptonChargeIncl;
      float Luminosity;
      if(monster->selectedSemiMu())
      {
        Luminosity = LuminosityMu;
        leptonCharge = "_mu";
        leptonDecay = "_mu";
        if(monster->semiMuDecay()) nSemiMu_TTSemiMu++;
        else nSemiMu++;
      }
      else
      {
        Luminosity = LuminosityEl;
        leptonCharge = "_el";
        leptonDecay = "_el";
        if(monster->semiElDecay()) nSemiEl_TTSemiEl++;
        else nSemiEl++;
      }
      if(monster->leptonCharge() == 1)
      {
        leptonCharge += "Plus";
        leptonChargeIncl = "_leptonPlus";
      }
      else if(monster->leptonCharge() == -1)
      {
        leptonCharge += "Minus";
        leptonChargeIncl = "_leptonMinus";
      }
      else cout << "leptonCharge = " << monster->leptonCharge() << endl;
      
      // the selected jets!
      vector<TLorentzVector> selectedJets = monster->selectedJets();
      TLorentzVector lepton = monster->lepton();
      TLorentzVector MET = monster->met();
      vector<float> btagSSVHE = monster->bTagSSVHE();
      
    	if(dataSetName.find("DAtaDriven") != string::npos)
    	{
    	  if(monster->selectedSemiMu())
    	  {
    	    if(monster->runID() >= 160431 && monster->runID() <= 165633) monster->setEventWeight( monster->eventWeight() * 0.0814781/0.717103 );
    	    else if(monster->runID() >= 165970 && monster->runID() <= 173692) monster->setEventWeight( monster->eventWeight() * 0.39043/0.11016 );
    	    else if(monster->runID() >= 175860 && monster->runID() <= 180252) monster->setEventWeight( monster->eventWeight() * 0.528092/0.172738 );
    	    else cout << "Unknown run for Data-driven QCD reweighting:  " << monster->runID() << endl;
    	  }
    	  else
    	  {
     	    if(monster->runID() >= 160431 && monster->runID() <= 165633) monster->setEventWeight( monster->eventWeight() * 0.0849205/0.248866 );
    	    else if(monster->runID() >= 165970 && monster->runID() <= 173692) monster->setEventWeight( monster->eventWeight() * 0.448948/0.34194 );
    	    else if(monster->runID() >= 175860 && monster->runID() <= 180252) monster->setEventWeight( monster->eventWeight() * 0.466132/0.409194 );
    	    else cout << "Unknown run for Data-driven QCD reweighting:  " << monster->runID() << endl;
    	  }
    	  histo1D["RunNr_DataDrivenQCD"+leptonDecay]->Fill(monster->runID(), monster->eventWeight());
    	  histo2D["RunNr_DataDrivenQCD_VS_Lepton_pt"+leptonDecay]->Fill(monster->runID(),lepton.Pt());
    	  histo2D["RunNr_DataDrivenQCD_VS_Lepton_eta"+leptonDecay]->Fill(monster->runID(),lepton.Eta());
      }
      if(dataSetName.find("Data") == 0)
      {
        histo1D["RunNr_Data"+leptonDecay]->Fill(monster->runID(), monster->eventWeight());
        histo2D["RunNr_Data_VS_Lepton_pt"+leptonDecay]->Fill(monster->runID(),lepton.Pt());
        histo2D["RunNr_Data_VS_Lepton_eta"+leptonDecay]->Fill(monster->runID(),lepton.Eta());
      }
      
      histo1D["lumiWeights"]->Fill(lumiWeight);
      histo1D["lumiWeightsUp"]->Fill(lumiWeightUp);
      histo1D["lumiWeightsDown"]->Fill(lumiWeightDown);
      if(dataSetName.find("Data") == 0 || dataSetName.find("DAtaDriven") != string::npos)
        lumiWeight = 1;
      if( dataSetName.find("Fall10") != string::npos )
        lumiWeight = 1; //no PU in Fall10!
      
      int nBtags = 0;
      for(unsigned int i=0; i<selectedJets.size(); i++)
        if(btagSSVHE[i] > bTagCut) nBtags++;
      
      if( dataSetName.find("TTbarJets") == 0 || dataSetName.find("TT_") == 0 )
      {
        histo1D["mTop_gen"]->Fill(monster->topMass());
        histo1D["mTop_gen"]->Fill(monster->antiTopMass());
        histo1D["mTopDiff_gen"]->Fill( monster->topMass() - monster->antiTopMass() );
        float etaMostForwardJet = selectedJets[0].Eta();
        for(unsigned int i=1; i< selectedJets.size(); i++)
          if( fabs(selectedJets[i].Eta()) > fabs(etaMostForwardJet) )
            etaMostForwardJet = selectedJets[i].Eta();
        histo1D["EtaMostForwardJet"]->Fill(etaMostForwardJet);
        if( monster->hadrBJet() < 9999 && monster->hadrLJet1() < 9999 && monster->hadrLJet2() < 9999 )
        {
          int maxIndex = max(monster->hadrLJet1(), monster->hadrLJet2());
          maxIndex = max(maxIndex, monster->hadrBJet());
          histo1D["maxIndexHadrJets"]->Fill(maxIndex);
        }
      }
//      cout << "Luminosity: " << Luminosity << "  monster->eventWeight(): " << monster->eventWeight() << "  lumiWeight: " << lumiWeight << endl;
      MSPlot["nPV_beforechi2"+leptonDecay]->Fill(monster->nPV(), dataSet, true, monster->eventWeight()*lumiWeight);
      MSPlot["RelPFISo"+leptonDecay]->Fill(monster->leptonPFRelIso(), dataSet, true, monster->eventWeight()*lumiWeight);
      
      // Extra mTop and sigma(mTop) from fitresults
      vector<float> mTop, sigmaMtop, minChi2;
      vector<int> hadrBindex, light1index, light2index;
      int nFinalCombis=0; // final means after chi2<10 cut
      int mcCombiIndex = -1;
      for(unsigned int iCombi=0; iCombi<12; iCombi++)
      {
        unsigned int* combi = monster->mvaResult(iCombi);
        bool goodCombi = hadrJetsMVAMatched(monster, iCombi);
        if(goodCombi) mcCombiIndex = nFinalCombis;
        
//        cout << "goodCombi: " << goodCombi << endl;
        
        float mTopFit = monster->mTopFit(iCombi);
        float sigmaMTopFit = monster->sigmaMTopFit(iCombi);
        float chi2MTopFit = monster->chi2MTopFit(iCombi);
        
//        if( chi2MTopFit < 99999 ) // 99999 means no mTopFit value was calculated, see FullKinFit.cc
        if( chi2MTopFit < maxChi2 )
        {
          nFinalCombis++;
          light1index.push_back(combi[0]);
          light2index.push_back(combi[1]);
          hadrBindex.push_back(combi[2]);
          mTop.push_back(mTopFit);
          sigmaMtop.push_back(sigmaMTopFit);
          minChi2.push_back(chi2MTopFit);
          
          // control plots
          if( dataSetName.find("TTbarJets") == 0 )
          {
            histo1D["mTop_Fit"]->Fill(mTopFit);
            histo1D["sigmaMtop_Fit"]->Fill(sigmaMTopFit);
            histo1D["MinChi2ndf_Fit"]->Fill(chi2MTopFit);
            histo1D["Prob_MinChi2_Fit"]->Fill(TMath::Prob(chi2MTopFit, 1));
            
            if( goodCombi )
            {
              float mTopNoFit = ( selectedJets[combi[0]] + selectedJets[combi[1]] + selectedJets[combi[2]] ).M();
              float lightCorr1 = 1 - resFitLightJets->EtCorrection(&selectedJets[combi[0]]);
              float lightCorr2 = 1 - resFitLightJets->EtCorrection(&selectedJets[combi[1]]);
              float bCorr = 1 - resFitBJets->EtCorrection(&selectedJets[combi[2]]);
              float mTopNoFitL7 = ( selectedJets[combi[0]]*lightCorr1 + selectedJets[combi[1]]*lightCorr2 + selectedJets[combi[2]]*bCorr ).M();
              
              histo1D["mTop_NoFit_goodcombi"]->Fill(mTopNoFit);
              histo1D["mTopL7_NoFit_goodcombi"]->Fill(mTopNoFitL7);
              histo1D["mTop_Fit_goodcombi"]->Fill(mTopFit);
              histo1D["sigmaMtop_Fit_goodcombi"]->Fill(sigmaMTopFit);
            	histo1D["MinChi2ndf_Fit_goodcombi"]->Fill(chi2MTopFit);
            	histo1D["Prob_MinChi2_Fit_goodcombi"]->Fill(TMath::Prob(chi2MTopFit, 1));
            }
          }
        }
      }
      
      int bestCombiIndex = 99;
      float minChi2BestCombi = 9999;
      for(unsigned int iCombi=0; iCombi<mTop.size(); iCombi++)
      {
        if(minChi2[iCombi] < minChi2BestCombi)
        {
          bestCombiIndex = iCombi;
          minChi2BestCombi = minChi2[iCombi];
        }
      }
      
    	if( dataSet->Name().find("TT") == 0 )
    	  histo1D["nUsedCombis_TTJets"]->Fill(mTop.size());
    	else if( dataSet->Name().find("WJets") == 0 )
    	  histo1D["nUsedCombis_WJets"]->Fill(mTop.size());
      if(doAllMSPlots) MSPlot["nUsedCombis"+leptonCharge]->Fill(mTop.size(), dataSet, true, monster->eventWeight()*lumiWeight);
      
      if( mTop.size() > 0 ) // event still has at least one good kinFitted jet-combi
      {
        if(monster->hadrBJet() < 9999)
        {
          nRealBjets++;
          if(btagSSVHE[monster->hadrBJet()] > bTagCut) nRealBjetsBtag++;
        }
        
        if(monster->leptBJet() < 9999)
        {
          nRealBjets++;
          if(btagSSVHE[monster->leptBJet()] > bTagCut) nRealBjetsBtag++;
        }
        
        if(leptonCharge == "_muPlus") nMuPlus += (monster->eventWeight()*lumiWeight);
        else if(leptonCharge == "_elPlus") nElPlus += (monster->eventWeight()*lumiWeight);
        else if(leptonCharge == "_muMinus") nMuMinus += (monster->eventWeight()*lumiWeight);
        else if(leptonCharge == "_elMinus") nElMinus += (monster->eventWeight()*lumiWeight);
        else cout << "unknown leptonCharge:  " << leptonCharge << endl;
        
        MSPlot["nPV"+leptonDecay]->Fill(monster->nPV(), dataSet, true, monster->eventWeight()*lumiWeight);
        MSPlot["Nr_of_jets"]->Fill(selectedJets.size(), dataSet, true, monster->eventWeight()*lumiWeight);
        MSPlot["Nr_of_jets"+leptonChargeIncl]->Fill(selectedJets.size(), dataSet, true, monster->eventWeight()*lumiWeight);
        MSPlot["MinChi2ndf_Fit"+leptonChargeIncl]->Fill(minChi2[bestCombiIndex], dataSet, true, monster->eventWeight()*lumiWeight);
        MSPlot["mTop_Fit"+leptonChargeIncl]->Fill(mTop[bestCombiIndex], dataSet, true, monster->eventWeight()*lumiWeight);
        MSPlot["Leading_jet_pt"+leptonChargeIncl]->Fill(selectedJets[0].Pt(), dataSet, true, monster->eventWeight()*lumiWeight);
        MSPlot["4th_jet_pt"+leptonChargeIncl]->Fill(selectedJets[3].Pt(), dataSet, true, monster->eventWeight()*lumiWeight);
        MSPlot["MET"+leptonChargeIncl]->Fill(MET.Pt(), dataSet, true, monster->eventWeight()*lumiWeight);
        for(int i=0; i<mTop.size(); i++)
        {
          MSPlot["mTop_Fit_AllCombi"+leptonCharge]->Fill(mTop[i], dataSet, true, monster->eventWeight()*lumiWeight);
          MSPlot["sigmaMtop_Fit_AllCombi"+leptonCharge]->Fill(sigmaMtop[i], dataSet, true, monster->eventWeight()*lumiWeight);
          MSPlot["MinChi2ndf_Fit_AllCombi"+leptonCharge]->Fill(minChi2[i], dataSet, true, monster->eventWeight()*lumiWeight);
        }
        
        if(doAllMSPlots)
        {
          MSPlot["Lepton_pt"+leptonCharge]->Fill(lepton.Pt(), dataSet, true, monster->eventWeight()*lumiWeight);
          MSPlot["Lepton_et"+leptonCharge]->Fill(lepton.Et(), dataSet, true, monster->eventWeight()*lumiWeight);
          MSPlot["Lepton_eta"+leptonCharge]->Fill(lepton.Eta(), dataSet, true, monster->eventWeight()*lumiWeight);
          MSPlot["Lepton_phi"+leptonCharge]->Fill(lepton.Phi(), dataSet, true, monster->eventWeight()*lumiWeight);
          MSPlot["Leading_jet_pt"+leptonCharge]->Fill(selectedJets[0].Pt(), dataSet, true, monster->eventWeight()*lumiWeight);
	        MSPlot["Leading_jet_eta"+leptonCharge]->Fill(selectedJets[0].Eta(), dataSet, true, monster->eventWeight()*lumiWeight);
	        MSPlot["Leading_jet_phi"+leptonCharge]->Fill(selectedJets[0].Phi(), dataSet, true, monster->eventWeight()*lumiWeight);
	        MSPlot["2nd_jet_pt"+leptonCharge]->Fill(selectedJets[1].Pt(), dataSet, true, monster->eventWeight()*lumiWeight);
	        MSPlot["2nd_jet_eta"+leptonCharge]->Fill(selectedJets[1].Eta(), dataSet, true, monster->eventWeight()*lumiWeight);
	        MSPlot["2nd_jet_phi"+leptonCharge]->Fill(selectedJets[1].Phi(), dataSet, true, monster->eventWeight()*lumiWeight);
	        MSPlot["3rd_jet_pt"+leptonCharge]->Fill(selectedJets[2].Pt(), dataSet, true, monster->eventWeight()*lumiWeight);
	        MSPlot["3rd_jet_eta"+leptonCharge]->Fill(selectedJets[2].Eta(), dataSet, true, monster->eventWeight()*lumiWeight);
	        MSPlot["3rd_jet_phi"+leptonCharge]->Fill(selectedJets[2].Phi(), dataSet, true, monster->eventWeight()*lumiWeight);
	        MSPlot["4th_jet_pt"+leptonCharge]->Fill(selectedJets[3].Pt(), dataSet, true, monster->eventWeight()*lumiWeight);
	        MSPlot["4th_jet_eta"+leptonCharge]->Fill(selectedJets[3].Eta(), dataSet, true, monster->eventWeight()*lumiWeight);
	        MSPlot["4th_jet_phi"+leptonCharge]->Fill(selectedJets[3].Phi(), dataSet, true, monster->eventWeight()*lumiWeight);
          MSPlot["Nr_of_jets"+leptonCharge]->Fill(selectedJets.size(), dataSet, true, monster->eventWeight()*lumiWeight);
          MSPlot["nBtags"+leptonCharge]->Fill(nBtags, dataSet, true, monster->eventWeight()*lumiWeight);
          MSPlot["Angle_1stjet_2ndjet"+leptonCharge]->Fill(selectedJets[0].Angle(selectedJets[1].Vect()), dataSet, true, monster->eventWeight()*lumiWeight);
          MSPlot["Angle_1stjet_3rdjet"+leptonCharge]->Fill(selectedJets[0].Angle(selectedJets[2].Vect()), dataSet, true, monster->eventWeight()*lumiWeight);
          MSPlot["MET"+leptonCharge]->Fill(MET.Pt(), dataSet, true, monster->eventWeight()*lumiWeight);
          
          float mTopUnfitted = (selectedJets[light1index[bestCombiIndex]]+selectedJets[light2index[bestCombiIndex]]+selectedJets[hadrBindex[bestCombiIndex]]).M();
          float mWUnfitted = (selectedJets[light1index[bestCombiIndex]] + selectedJets[light2index[bestCombiIndex]]).M();
          MSPlot["mW_UnFitted"+leptonCharge]->Fill(mWUnfitted, dataSet, true, monster->eventWeight()*lumiWeight);
          MSPlot["mTop_UnFitted"+leptonCharge]->Fill(mTopUnfitted, dataSet, true, monster->eventWeight()*lumiWeight);
          MSPlot["mTop_Fit"+leptonCharge]->Fill(mTop[bestCombiIndex], dataSet, true, monster->eventWeight()*lumiWeight);
          MSPlot["sigmaMtop_Fit"+leptonCharge]->Fill(sigmaMtop[bestCombiIndex], dataSet, true, monster->eventWeight()*lumiWeight);
          MSPlot["MinChi2ndf_Fit"+leptonCharge]->Fill(minChi2[bestCombiIndex], dataSet, true, monster->eventWeight()*lumiWeight);
          MSPlot["Prob_MinChi2_Fit"+leptonCharge]->Fill(TMath::Prob(minChi2[bestCombiIndex], 1), dataSet, true, monster->eventWeight()*lumiWeight);
        }
        
        float M3 = -1, maxPt = -1, PtTTbar = 0;
        for(int i=0;i<selectedJets.size();i++)
        {
          MSPlot["AllJets_pt"+leptonChargeIncl]->Fill(selectedJets[i].Pt(), dataSet, true, monster->eventWeight()*lumiWeight);
          if(doAllMSPlots)
          {
            MSPlot["AllJets_pt"+leptonCharge]->Fill(selectedJets[i].Pt(), dataSet, true, monster->eventWeight()*lumiWeight);
            MSPlot["AllJets_eta"+leptonCharge]->Fill(selectedJets[i].Eta(), dataSet, true, monster->eventWeight()*lumiWeight);
            MSPlot["AllJets_phi"+leptonCharge]->Fill(selectedJets[i].Phi(), dataSet, true, monster->eventWeight()*lumiWeight);
          }  
          for(int j=0;j<i;j++)
          {
            for(int k=0;k<j;k++)
            {
              float combinedPt = (selectedJets[i]+selectedJets[j]+selectedJets[k]).Pt();
              if(combinedPt > maxPt)
              {
                maxPt = combinedPt;
                M3 = (selectedJets[i]+selectedJets[j]+selectedJets[k]).M();
                PtTTbar = (selectedJets[i]+selectedJets[j]+selectedJets[k]+MET+lepton).Pt();
              }
            }
          }
        }
        MSPlot["M3"]->Fill(M3, dataSet, true, monster->eventWeight()*lumiWeight);
        MSPlot["PtTop"]->Fill(maxPt, dataSet, true, monster->eventWeight()*lumiWeight);
        MSPlot["PtTTbar"]->Fill(PtTTbar, dataSet, true, monster->eventWeight()*lumiWeight);
        if(PtTTbar>200) MSPlot["M3_PtTTbar200"]->Fill(M3, dataSet, true, monster->eventWeight()*lumiWeight);
        if(maxPt>300) MSPlot["M3_PtTop300"]->Fill(M3, dataSet, true, monster->eventWeight()*lumiWeight);
        if(doAllMSPlots) MSPlot["M3"+leptonCharge]->Fill(M3, dataSet, true, monster->eventWeight()*lumiWeight);
        
        // Write everything out to the ASCII file
        if(monster->selectedSemiMu())
        {
          fillMassFitInfoEvent(monster, headerMu, eventMu, lumiWeight, lumiWeightUp, lumiWeightDown);
          if(applyChi2Cut(eventMu, maxChi2))
          {
            setBtagWeights(eventMu, bTagCut, bTagEff, misTagRate);
//            cout << eventMu << endl;
            writer(outFileLikelihoodSemiMu,eventMu);
          }
          
          nSemiMuInASCII++;
          outFileSemiMu << "------------------------------------------" << endl;
          outFileSemiMu << "BeginEvent" << endl;
          outFileSemiMu << "RunNr: " << monster->runID() << endl;
          outFileSemiMu << "LumiSection: " << monster->lumiBlockID() << endl;
          outFileSemiMu << "EventNr: " << monster->eventID() << endl;
          outFileSemiMu << "EventWeight: " << monster->eventWeight() << "  " << lumiWeight << "  " << lumiWeightDown << "  " << lumiWeightUp << endl;
          outFileSemiMu << selectedJets.size() << "  " << monster->nPV() << "  " << monster->nPU() << "  " << monster->topMass() << "  " << monster->antiTopMass() << "  " << monster->topDecayedLept() << endl;
          outFileSemiMu << lepton.Phi() << "  " << lepton.Eta() << "  " << lepton.Pt() << "  " << monster->leptonCharge() << "  " << MET.Px() << "  " << MET.Py() << endl;
          outFileSemiMu << selectedJets[0].Pt() << "  " << selectedJets[0].Eta() << "  " << selectedJets[0].Phi() << "  " << btagSSVHE[0] << "  "
            << selectedJets[1].Pt() << "  " << selectedJets[1].Eta() << "  " << selectedJets[1].Phi() << "  " << btagSSVHE[1] << "  "
            << selectedJets[2].Pt() << "  " << selectedJets[2].Eta() << "  " << selectedJets[2].Phi() << "  " << btagSSVHE[2] << "  " 
            << selectedJets[3].Pt() << "  " << selectedJets[3].Eta() << "  " << selectedJets[3].Phi() << "  " << btagSSVHE[3] << endl;
          outFileSemiMu << "goodCombi: " << monster->hadrLJet1() << "  " << monster->hadrLJet2() << "  " << monster->hadrBJet() << "  index: " << mcCombiIndex << endl;
          for(unsigned int iCombi=0; iCombi<mTop.size(); iCombi++)
            outFileSemiMu << mTop[iCombi] << "  " << sigmaMtop[iCombi] << "  " << minChi2[iCombi] << "  " << light1index[iCombi] << "  " << light2index[iCombi] << "  " << hadrBindex[iCombi] << endl;
          outFileSemiMu << "FinishedEvent" << endl;
          outFileSemiMu << "------------------------------------------" << endl;
        }
        else // selected semiEl
        {
          fillMassFitInfoEvent(monster, headerEl, eventEl, lumiWeight, lumiWeightUp, lumiWeightDown);
          if(applyChi2Cut(eventEl, maxChi2))
          {
            setBtagWeights(eventEl, bTagCut, bTagEff, misTagRate);
//            cout << eventEl << endl;
            writer(outFileLikelihoodSemiEl,eventEl);
          }
          
          nSemiElInASCII++;
          outFileSemiEl << "------------------------------------------" << endl;
          outFileSemiEl << "BeginEvent" << endl;
          outFileSemiEl << "RunNr: " << monster->runID() << endl;
          outFileSemiEl << "LumiSection: " << monster->lumiBlockID() << endl;
          outFileSemiEl << "EventNr: " << monster->eventID() << endl;
          outFileSemiEl << "EventWeight: " << monster->eventWeight() << "  " << lumiWeight << "  " << lumiWeightDown << "  " << lumiWeightUp << endl;
          outFileSemiEl << selectedJets.size() << "  " << monster->nPV() << "  " << monster->nPU() << "  " << monster->topMass() << "  " << monster->antiTopMass() << "  " << monster->topDecayedLept() << endl;
          outFileSemiEl << lepton.Phi() << "  " << lepton.Eta() << "  " << lepton.Pt() << "  " << monster->leptonCharge() << "  " << MET.Px() << "  " << MET.Py() << endl;
          outFileSemiEl << selectedJets[0].Pt() << "  " << selectedJets[0].Eta() << "  " << selectedJets[0].Phi() << "  " << btagSSVHE[0] << "  "
            << selectedJets[1].Pt() << "  " << selectedJets[1].Eta() << "  " << selectedJets[1].Phi() << "  " << btagSSVHE[1] << "  "
            << selectedJets[2].Pt() << "  " << selectedJets[2].Eta() << "  " << selectedJets[2].Phi() << "  " << btagSSVHE[2] << "  " 
            << selectedJets[3].Pt() << "  " << selectedJets[3].Eta() << "  " << selectedJets[3].Phi() << "  " << btagSSVHE[3] << endl;
          outFileSemiEl << "goodCombi: " << monster->hadrLJet1() << "  " << monster->hadrLJet2() << "  " << monster->hadrBJet() << "  index: " << mcCombiIndex << endl;
          for(unsigned int iCombi=0; iCombi<mTop.size(); iCombi++)
            outFileSemiEl << mTop[iCombi] << "  " << sigmaMtop[iCombi] << "  " << minChi2[iCombi] << "  " << light1index[iCombi] << "  " << light2index[iCombi] << "  " << hadrBindex[iCombi] << endl;
          outFileSemiEl << "FinishedEvent" << endl;
          outFileSemiEl << "------------------------------------------" << endl;
        }
        // normalize chi2-weights to 1
/*        float totalWeight = 0;
        for(unsigned int iCombi=0; iCombi<mTop.size(); iCombi++) totalWeight += TMath::Exp(-0.5*minChi2[iCombi]);
        
        // shape study
        for(unsigned int iCombi=0; iCombi<mTop.size(); iCombi++)
        {
//          if( minChi2[iCombi] < maxChi2 )
          {
            histo1D["combiWeight"]->Fill(TMath::Exp(-0.5*minChi2[iCombi]));
            if( dataSetName.find("WJets") == 0 )
            {
              histo1D["wJets_mTopFitted"]->Fill(mTop[iCombi], monster->eventWeight()*lumiWeight*TMath::Exp(-0.5*minChi2[iCombi]));
              histo1D["wJets_mTopFitted_noChi2weight"]->Fill(mTop[iCombi], monster->eventWeight()*lumiWeight);
              histo1D["wJets_mTopFitted"+leptonCharge]->Fill(mTop[iCombi], monster->eventWeight()*lumiWeight*TMath::Exp(-0.5*minChi2[iCombi]));
              histo1D["wJets_combiWeight"]->Fill(TMath::Exp(-0.5*minChi2[iCombi]));
              histo2D["wJets_mTopFitted_VS_combiWeight"]->Fill(mTop[iCombi],TMath::Exp(-0.5*minChi2[iCombi]));
            }
            if( dataSetName.find("TTbarJets") == 0 && mcCombiIndex != iCombi )
            {
              histo1D["ttJets_mTopFitted_badcombi"]->Fill(mTop[iCombi], monster->eventWeight()*lumiWeight*TMath::Exp(-0.5*minChi2[iCombi]));
              histo1D["ttJets_mTopFitted_badcombi_noChi2weight"]->Fill(mTop[iCombi], monster->eventWeight()*lumiWeight);
              histo1D["ttJets_mTopFitted_badcombi"+leptonCharge]->Fill(mTop[iCombi], monster->eventWeight()*lumiWeight*TMath::Exp(-0.5*minChi2[iCombi]));
            }
          }
        }*/
      }
    } // end loop over monsters
    
    outFileSemiMu << endl << "Finished event-by-event info, end of file!" << endl;
    outFileSemiMu.close();
    outFileSemiEl << endl << "Finished event-by-event info, end of file!" << endl;
    outFileSemiEl.close();
    
    outFileLikelihoodSemiMu << endl << "Finished event-by-event info, end of file!" << endl;
    outFileLikelihoodSemiMu.close();
    outFileLikelihoodSemiEl << endl << "Finished event-by-event info, end of file!" << endl;
    outFileLikelihoodSemiEl.close();
    
    cout << dataSetName << " nSemiMuInASCII = " << nSemiMuInASCII << "  nSemiElInASCII = " << nSemiElInASCII << endl;
    
    cout << "nEvents selected before lumi-rescaling:" << endl;
    cout << "muPlus: " << nMuPlus << " +- " << sqrt(nMuPlus) << "  muMinus: " << nMuMinus << " +- " << sqrt(nMuMinus) << "  elPlus: " << nElPlus << " +- " << sqrt(nElPlus) << "  elMinus: " << nElMinus << " +- " << sqrt(nElMinus) << endl;
    float factorMu = LuminosityMu*dataSet->NormFactor();
    float factorEl = LuminosityEl*dataSet->NormFactor();
    cout << "nEvents selected after lumi-rescaling to:  Mu: " << LuminosityMu << "  El: " << LuminosityEl << "  pb-1" << endl;
    cout << "muPlus: " << nMuPlus*factorMu << " +- " << sqrt(nMuPlus)*factorMu << "  muMinus: " << nMuMinus*factorMu << " +- " << sqrt(nMuMinus)*factorMu << "  elPlus: " << nElPlus*factorEl << " +- " << sqrt(nElPlus)*factorEl << "  elMinus: " << nElMinus*factorEl << " +- " << sqrt(nElMinus)*factorEl << endl;
    
    if( dataSetName.find("TTbarJets") == 0 )
    {
      cout << " nSemiLepton_NoSel = " << dataSet->Xsection() * dataSet->EquivalentLumi() * 4/27 << endl;
      cout << " nSemiMu_TTSemiMu = " << nSemiMu_TTSemiMu << "  selEff = " << (float) nSemiMu_TTSemiMu / ( dataSet->Xsection() * dataSet->EquivalentLumi() * 4/27 ) << endl;
      cout << " nSemiEl_TTSemiEl = " << nSemiEl_TTSemiEl << "  selEff = " << (float) nSemiEl_TTSemiEl / ( dataSet->Xsection() * dataSet->EquivalentLumi() * 4/27 ) << endl;
      cout << " nOther_NoSel = " << dataSet->Xsection() * dataSet->EquivalentLumi() * 23/27 << endl;
      cout << " nSemiMu = " << nSemiMu << "  selEff = " << (float) nSemiMu / ( dataSet->Xsection() * dataSet->EquivalentLumi() * 23/27 ) << endl;
      cout << " nSemiEl = " << nSemiEl << "  selEff = " << (float) nSemiEl / ( dataSet->Xsection() * dataSet->EquivalentLumi() * 23/27 ) << endl;
      
      cout << " b-tag eff = " <<  ( (float) nRealBjetsBtag ) / ( (float) nRealBjets ) *100 << " +- " << sqrt( (nRealBjetsBtag/pow((float)nRealBjets,2)) + (pow((float)nRealBjetsBtag,2)/pow((float)nRealBjets,3)) ) *100 << " %" << endl;
    }
    else
    {
      cout << " n_NoSel = " << dataSet->Xsection() * dataSet->EquivalentLumi() << endl;
      cout << " nSemiMu = " << nSemiMu << "  selEff = " << (float) nSemiMu / ( dataSet->Xsection() * dataSet->EquivalentLumi() ) << endl;
      cout << " nSemiEl = " << nSemiEl << "  selEff = " << (float) nSemiEl / ( dataSet->Xsection() * dataSet->EquivalentLumi() ) << endl;
    }
    inFile->Close();
    delete inFile;
  } // end loop over datasets
  
  cout << "Finished running over all datasets..." << endl;
  fout->cd();
  
  float nSemiMu1 = histo1D["RunNr_Data_mu"]->Integral(histo1D["RunNr_Data_mu"]->GetXaxis()->FindBin(160577), histo1D["RunNr_Data_mu"]->GetXaxis()->FindBin(165633));
  float nSemiMu2 = histo1D["RunNr_Data_mu"]->Integral(histo1D["RunNr_Data_mu"]->GetXaxis()->FindBin(165970), histo1D["RunNr_Data_mu"]->GetXaxis()->FindBin(173692));
  float nSemiMu3 = histo1D["RunNr_Data_mu"]->Integral(histo1D["RunNr_Data_mu"]->GetXaxis()->FindBin(175860), histo1D["RunNr_Data_mu"]->GetXaxis()->FindBin(180252));
  cout << "nSemiMu in Data: " << nSemiMu1 << " " << nSemiMu2 << " " << nSemiMu3 << endl;
  cout << "Fractions:       " << nSemiMu1/(nSemiMu1+nSemiMu2+nSemiMu3) << " " << nSemiMu2/(nSemiMu1+nSemiMu2+nSemiMu3) << " " << nSemiMu3/(nSemiMu1+nSemiMu2+nSemiMu3) << endl;
  
  float nSemiMuQCD1 = histo1D["RunNr_DataDrivenQCD_mu"]->Integral(histo1D["RunNr_DataDrivenQCD_mu"]->GetXaxis()->FindBin(160577), histo1D["RunNr_DataDrivenQCD_mu"]->GetXaxis()->FindBin(165633));
  float nSemiMuQCD2 = histo1D["RunNr_DataDrivenQCD_mu"]->Integral(histo1D["RunNr_DataDrivenQCD_mu"]->GetXaxis()->FindBin(165970), histo1D["RunNr_DataDrivenQCD_mu"]->GetXaxis()->FindBin(173692));
  float nSemiMuQCD3 = histo1D["RunNr_DataDrivenQCD_mu"]->Integral(histo1D["RunNr_DataDrivenQCD_mu"]->GetXaxis()->FindBin(175860), histo1D["RunNr_DataDrivenQCD_mu"]->GetXaxis()->FindBin(180252));
	cout << "nSemiMuQCD in Data: " << nSemiMuQCD1 << " " << nSemiMuQCD2 << " " << nSemiMuQCD3 << endl;
  cout << "Fractions:       " << nSemiMuQCD1/(nSemiMuQCD1+nSemiMuQCD2+nSemiMuQCD3) << " " << nSemiMuQCD2/(nSemiMuQCD1+nSemiMuQCD2+nSemiMuQCD3) << " " << nSemiMuQCD3/(nSemiMuQCD1+nSemiMuQCD2+nSemiMuQCD3) << endl;
  
  float nSemiEl1 = histo1D["RunNr_Data_el"]->Integral(histo1D["RunNr_Data_el"]->GetXaxis()->FindBin(160577), histo1D["RunNr_Data_el"]->GetXaxis()->FindBin(165633));
  float nSemiEl2 = histo1D["RunNr_Data_el"]->Integral(histo1D["RunNr_Data_el"]->GetXaxis()->FindBin(165970), histo1D["RunNr_Data_el"]->GetXaxis()->FindBin(173692));
  float nSemiEl3 = histo1D["RunNr_Data_el"]->Integral(histo1D["RunNr_Data_el"]->GetXaxis()->FindBin(175860), histo1D["RunNr_Data_el"]->GetXaxis()->FindBin(180252));
	cout << "nSemiEl in Data: " << nSemiEl1 << " " << nSemiEl2 << " " << nSemiEl3 << endl;
  cout << "Fractions:       " << nSemiEl1/(nSemiEl1+nSemiEl2+nSemiEl3) << " " << nSemiEl2/(nSemiEl1+nSemiEl2+nSemiEl3) << " " << nSemiEl3/(nSemiEl1+nSemiEl2+nSemiEl3) << endl;
  
  float nSemiElQCD1 = histo1D["RunNr_DataDrivenQCD_el"]->Integral(histo1D["RunNr_DataDrivenQCD_el"]->GetXaxis()->FindBin(160577), histo1D["RunNr_DataDrivenQCD_el"]->GetXaxis()->FindBin(165633));
  float nSemiElQCD2 = histo1D["RunNr_DataDrivenQCD_el"]->Integral(histo1D["RunNr_DataDrivenQCD_el"]->GetXaxis()->FindBin(165970), histo1D["RunNr_DataDrivenQCD_el"]->GetXaxis()->FindBin(173692));
  float nSemiElQCD3 = histo1D["RunNr_DataDrivenQCD_el"]->Integral(histo1D["RunNr_DataDrivenQCD_el"]->GetXaxis()->FindBin(175860), histo1D["RunNr_DataDrivenQCD_el"]->GetXaxis()->FindBin(180252));
	cout << "nSemiElQCD in Data: " << nSemiElQCD1 << " " << nSemiElQCD2 << " " << nSemiElQCD3 << endl;
  cout << "Fractions:       " << nSemiElQCD1/(nSemiElQCD1+nSemiElQCD2+nSemiElQCD3) << " " << nSemiElQCD2/(nSemiElQCD1+nSemiElQCD2+nSemiElQCD3) << " " << nSemiElQCD3/(nSemiElQCD1+nSemiElQCD2+nSemiElQCD3) << endl;
  
  histo1D["mTop_NoFit_goodcombi"]->Fit("gaus","QR","",155,190);
  histo1D["mTopL7_NoFit_goodcombi"]->Fit("gaus","QR","",155,190);
  histo1D["mTop_Fit_goodcombi"]->Fit("gaus","QR","",160,185);
  
  // Fit some histo's
  if( histo1D["wJets_mTopFitted"]->GetEntries() > 2 )
  {
	  RooRealVar mTopFitWjets("mTopFitWjets","mTopFitWjets",0.,500.);
	  RooDataHist dataHistWjets("dataHistWjets","dataHistWjets",mTopFitWjets,histo1D["wJets_mTopFitted"]);
	
	  RooRealVar mean("mean","mean",0.,500.);
	  RooRealVar sigma("sigma","sigma",0.,500.);
    RooLandau wJetsShapeLandau("wJetsShapeLandau","wJetsShapeLandau",mTopFitWjets,mean,sigma);
    wJetsShapeLandau.fitTo(dataHistWjets, SumW2Error(true), PrintLevel(-3), Verbose(false), Extended(false));
    
	  RooRealVar mean2010("mean2010","mean2010",169.092); // 2010 values
	  RooRealVar sigma2010("sigma2010","sigma2010",24.825); // 2010 values
    RooLandau wJetsShape2010("wJetsShape2010","wJetsShape2010",mTopFitWjets,mean2010,sigma2010);
    
    RooPlot* plot = mTopFitWjets.frame();
    dataHistWjets.plotOn(plot);
    wJetsShapeLandau.plotOn(plot);
//    wJetsShape2010.plotOn(plot, LineColor(6));
    
    TCanvas* cFit = new TCanvas("mTopFitWjets","mTopFitWjets");
    cFit->cd();
    plot->Draw();
    cFit->SaveAs( (pathPNG+"mTopFitWjets.png").c_str() );
    cFit->Write();
    
    RooPlot* plotMuPlus = mTopFitWjets.frame();
    RooDataHist dataHistWjetsMuPlus("dataHistWjetsMuPlus","dataHistWjetsMuPlus",mTopFitWjets,histo1D["wJets_mTopFitted_muPlus"]);
    dataHistWjetsMuPlus.plotOn(plotMuPlus);
    wJetsShapeLandau.plotOn(plotMuPlus);
    
    TCanvas* cFitMuPlus = new TCanvas("mTopFitWjetsMuPlus","mTopFitWjetsMuPlus");
    cFitMuPlus->cd();
    plotMuPlus->Draw();
    cFitMuPlus->SaveAs( (pathPNG+"mTopFitWjets_MuPlus.png").c_str() );
    cFitMuPlus->Write();
    
    RooPlot* plotMuMinus = mTopFitWjets.frame();
    RooDataHist dataHistWjetsMuMinus("dataHistWjetsMuMinus","dataHistWjetsMuMinus",mTopFitWjets,histo1D["wJets_mTopFitted_muMinus"]);
    dataHistWjetsMuMinus.plotOn(plotMuMinus);
    wJetsShapeLandau.plotOn(plotMuMinus);
    
    TCanvas* cFitMuMinus = new TCanvas("mTopFitWjetsMuMinus","mTopFitWjetsMuMinus");
    cFitMuMinus->cd();
    plotMuMinus->Draw();
    cFitMuMinus->SaveAs( (pathPNG+"mTopFitWjets_MuMinus.png").c_str() );
    cFitMuMinus->Write();
    
    cout << "Shape-fit results for WJets: " << endl;
    cout << "Landau function with:  mean = " << mean.getVal() << " +- " << mean.getError() << "  sigma = " << sigma.getVal() << " +- " << sigma.getError() << endl;
  }
  
  if( histo1D["ttJets_mTopFitted_badcombi"]->GetEntries() > 2 )
  {
	  RooRealVar mTopFitTTjets("mTopFitTTjets","mTopFitTTjets",0.,500.);
	  RooDataHist dataHistTTjets("dataHistTTjets","dataHistTTjets",mTopFitTTjets,histo1D["ttJets_mTopFitted_badcombi"]);
	  
	  RooRealVar mean("mean","mean",154.,1.,500.);
	  RooRealVar sigma("sigma","sigma",20.,2.,500.);
	  RooRealVar alpha("alpha","alpha",-0.35,-10.,10.);
	  RooRealVar N("N","N",5);
	  RooCBShape TTJetsShapeCB("TTJetsShapeCB","TTJetsShapeCB",mTopFitTTjets,mean,sigma,alpha,N);
	  
	  RooRealVar meanLandau("meanLandau","meanLandau",125.,1.,500.);
	  RooRealVar sigmaLandau("sigmaLandau","sigmaLandau",10.,2.,500.);
    RooLandau TTJetsShapeLandau("TTJetsShapeLandau","TTJetsShapeLandau",mTopFitTTjets,meanLandau,sigmaLandau);
    
    // --- Construct signal+background PDF ---
    RooRealVar nsig("nsig","#signal events",0.,100000000);
    RooRealVar nbkg("nbkg","#background events",0.,100000000);
    RooAddPdf TTJetsShapeTotal("TTJetsShapeTotal","TTJetsShapeTotal",RooArgList(TTJetsShapeCB,TTJetsShapeLandau),RooArgList(nsig,nbkg));
    TTJetsShapeTotal.fitTo(dataHistTTjets, SumW2Error(false), PrintLevel(-3), Verbose(false), Extended(true));
    
    float truemass = 172.5;
    float cbmean = 15.3679 + 1.01314*truemass; // 2010 value
    float cbsigma = 152.705 - 0.426815*truemass; // 2010 value
//    float cbalpha = 250.0 + 0.0*truemass; // 2010 value, fixed???
    float cbalpha = -3308.87 + 18.208*truemass; // 2010 value
	  RooRealVar mean2010("mean2010","mean2010",cbmean); // 2010 value
	  RooRealVar sigma2010("sigma2010","sigma2010",cbsigma); // 201 value
	  RooRealVar alpha2010("alpha2010","alpha2010",cbalpha); // 2010 value
	  RooRealVar N2010("N2010","N2010",5); // 2010 value
    RooCBShape TTJetsShape2010("TTJetsShape2010","TTJetsShape2010",mTopFitTTjets,mean2010,sigma2010,alpha2010,N2010); // 2010 shape
    
    RooRealVar mean2011("mean2011","mean2011",166.567);
	  RooRealVar sigma2011("sigma2011","sigma2011",21.9164);
	  RooRealVar alpha2011("alpha2011","alpha2011",-0.397863);
	  RooCBShape TTJetsShapeCB2011("TTJetsShapeCB2011","TTJetsShapeCB2011",mTopFitTTjets,mean2011,sigma2011,alpha2011,N);
	  
	  RooRealVar meanLandau2011("meanLandau2011","meanLandau2011",138.145);
	  RooRealVar sigmaLandau2011("sigmaLandau2011","sigmaLandau2011",15.918);
    RooLandau TTJetsShapeLandau2011("TTJetsShapeLandau2011","TTJetsShapeLandau2011",mTopFitTTjets,meanLandau2011,sigmaLandau2011);
    
    // --- Construct signal+background PDF ---
    RooRealVar nsig2011("nsig2011","#signal events",54735.7);
    RooRealVar nbkg2011("nbkg2011","#background events",34110);
    RooAddPdf TTJetsShapeTotal2011("TTJetsShapeTotal2011","TTJetsShapeTotal2011",RooArgList(TTJetsShapeCB2011,TTJetsShapeLandau2011),RooArgList(nsig2011,nbkg2011));

    
    RooPlot* plot = mTopFitTTjets.frame();
    dataHistTTjets.plotOn(plot);
    TTJetsShapeTotal.plotOn(plot);
    TTJetsShapeTotal.plotOn(plot, Components(TTJetsShapeCB), LineStyle(kDashed), LineColor(2));
    TTJetsShapeTotal.plotOn(plot, Components(TTJetsShapeLandau), LineStyle(kDashed), LineColor(3));
//    TTJetsShapeTotal2011.plotOn(plot, LineStyle(kDashed), LineColor(1));
//    TTJetsShape2010.plotOn(plot, LineColor(6));
    
    TCanvas* cFit = new TCanvas("mTopFitTTjets","mTopFitTTjets");
    cFit->cd();
    plot->Draw();
    cFit->SaveAs( (pathPNG+"mTopFitTTjets.png").c_str() );
    cFit->Write();
    
    RooPlot* plotMuPlus = mTopFitTTjets.frame();
    RooDataHist dataHistTTjetsMuPlus("dataHistTTjetsMuPlus","dataHistTTjetsMuPlus",mTopFitTTjets,histo1D["ttJets_mTopFitted_badcombi_muPlus"]);
    dataHistTTjetsMuPlus.plotOn(plotMuPlus);
    TTJetsShapeTotal.plotOn(plotMuPlus);
    
    TCanvas* cFitMuPlus = new TCanvas("mTopFitTTjetsMuPlus","mTopFitTTjetsMuPlus");
    cFitMuPlus->cd();
    plotMuPlus->Draw();
    cFitMuPlus->SaveAs( (pathPNG+"mTopFitTTjets_MuPlus.png").c_str() );
    cFitMuPlus->Write();
    
    RooPlot* plotMuMinus = mTopFitTTjets.frame();
    RooDataHist dataHistTTjetsMuMinus("dataHistTTjetsMuMinus","dataHistTTjetsMuMinus",mTopFitTTjets,histo1D["ttJets_mTopFitted_badcombi_muMinus"]);
    dataHistTTjetsMuMinus.plotOn(plotMuMinus);
    TTJetsShapeTotal.plotOn(plotMuMinus);
    
    TCanvas* cFitMuMinus = new TCanvas("mTopFitTTjetsMuMinus","mTopFitTTjetsMuMinus");
    cFitMuMinus->cd();
    plotMuMinus->Draw();
    cFitMuMinus->SaveAs( (pathPNG+"mTopFitTTjets_MuMinus.png").c_str() );
    cFitMuMinus->Write();
    
    cout << "Shape-fit results for TTJets:   f * Crystal-Ball + (1 - f) * Landau " << endl;
    TTJetsShapeTotal.Print("t");
//    cout << "Landau function with:  mean = " << meanLandau.getVal() << " +- " << meanLandau.getError() << "  sigma = " << sigmaLandau.getVal() << " +- " << sigmaLandau.getError() << endl;
//    cout << "CB function with:  mean = " << mean.getVal() << " +- " << mean.getError() << "  sigma = " << sigma.getVal() << " +- " << sigma.getError() <<
//      "  alpha = " << alpha.getVal() << " +- " << alpha.getError() << "  N = " << N.getVal() << " +- " << N.getError() << endl;
//    cout << " N1 = " << nbkg.getVal() << " +- " << nbkg.getError() << "  N2 = " << nsig.getVal() << " +- " << nsig.getError() << endl;
  }
  
  cout << "Writing out..." << endl;
  fout->cd();
  
  TDirectory* th1dir = fout->mkdir("1D_histograms");
  th1dir->cd();
  // Write 1D histo's
  for(std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
	{
		TH1F *temp = it->second;
		temp->Write();
		TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
		tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
		tempCanvas->Write();
//		tempCanvas->SaveAs( (pathPNG+it->first+".pdf").c_str() );
	}
	
	fout->cd();
  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
  {
    MultiSamplePlot *temp = it->second;
    string name = it->first;
    temp->showNumberEntries(true);
    if(name.find("muPlus") < name.size()) temp->addText("#mu^{+}+jets");
    else if(name.find("muMinus") < name.size()) temp->addText("#mu^{-}+jets");
    else if(name.find("elPlus") < name.size()) temp->addText("e^{+}+jets");
    else if(name.find("elMinus") < name.size()) temp->addText("e^{-}+jets");
    if(name.find("leptonPlus") < name.size()) temp->addText("l^{+}+jets");
    else if(name.find("leptonMinus") < name.size()) temp->addText("l^{-}+jets");
    temp->Draw(false, name, false, true, true, true, true, 1, true);
//    temp->Draw(false, name, false, false, false, false, false);
    temp->Write(fout, name, true, pathPNG+"MSPlot/","png");
//    temp->Write(fout, name, true, pathPNG+"MSPlot/","pdf");
  }
	
	// 2D
  TDirectory* th2dir = fout->mkdir("2D_histograms_graphs");
  th2dir->cd();
	for(std::map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++)
	{
		TH2F *temp = it->second;
		temp->Write();
		TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
		tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
//		tempCanvas->SaveAs( (pathPNG+it->first+".pdf").c_str() );
	}

	//Write TGraphAsymmErrors
	fout->cd();
	for(map<string,TGraphAsymmErrors*>::const_iterator it = graphAsymmErr.begin(); it != graphAsymmErr.end(); it++)
	{
	  TGraphAsymmErrors *temp = it->second;
	  temp->Write();
    TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
		tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
//		tempCanvas->SaveAs( (pathPNG+it->first+".pdf").c_str() );
	}
	
  //Write TGraphErrors
  fout->cd();
  for(map<string,TGraphErrors*>::const_iterator it = graphErr.begin(); it != graphErr.end(); it++)
  {
    TGraphErrors *temp = it->second;
    temp->Write();
    TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
    tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
//    tempCanvas->SaveAs( (pathPNG+it->first+".pdf").c_str() );
  }
  
  fout->Close();
  
  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;

  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "           hasn't crashed yet ;-)           " << endl;
  cout << "********************************************" << endl;
  
  return 0;
}

