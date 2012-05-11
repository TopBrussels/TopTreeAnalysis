// general
#include "TH1.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TF1.h"
#include "TKey.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TLine.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TFormula.h"
#include "TAxis.h"
#include "TRandom3.h"
#include "TMath.h"
#include "THStack.h"
#include "TString.h"

#include <cmath>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <string>
//#include "tdrstyle.C"

//user code
#include "../Tools/interface/PlottingTools.h"
#include "../Tools/interface/MultiSamplePlot.h"
#include "../Tools/interface/TTreeLoader.h"
#include "../Tools/interface/AnalysisEnvironmentLoader.h"
#include "../Content/interface/AnalysisEnvironment.h"
#include "../Content/interface/Dataset.h"
#include "../MCInformation/interface/MCWeighter.h"

#include "Style.C"

using namespace std;

/// MultiSamplePlot
//map<string,vector<pair<TH1F*,Dataset*> > > MSPlot_input;
map<string,MultiSamplePlot*> MSPlot;

TH1F*** create2D_TH1F_Array(int d1, int d2) {
    TH1F*** array = new TH1F**[d1];
    for (int i = 0; i < d1; i++) {
        array[i] = new TH1F*[d2];
    }
    return array;
}

void destroy2D_TH1F_Array(TH1F*** array, int d1){
    for (int i = 0; i < d1; i++) {
        delete [] array[i];
    }
    delete [] array;
}

int main(unsigned int argc, char *argv[])
{
	setTDRStyle();
  //gStyle->SetOptFit(0);
  //gStyle->SetOptStat(0);
	////////////////////////////////////
	// define and open the inputfiles //
	////////////////////////////////////
	string Outputpath = "Plots_FittedModel_PASPlotsTotalUncert_12May2012/";
	gSystem->mkdir(Outputpath.c_str());
	
	//string inputfiledir = "CalculateLimits/Apr2012_Systematics/INPUT_2012_04_16/";
	string inputfiledir = "/user/pvmulder/NewEraOfDataAnalysis/TopTree/CMSSW_4_2_8_patch7/src/TopBrussels/TopTreeAnalysis/macros/OutputFiles_InclFourthGenTreeAnalyzer_8May_WmassMu/";
	bool mergeSignal = true;
	
	cout<<"The arguments passed to the executable are: "<<endl;
  for(unsigned int i=1;i<argc;i++)
	{
		cout<<argv[i]<<endl;
	}
	
	bool semiLepton = true; //if true, semiElectron and semiMuon will be combined
	string channelpostfix = "";
	bool semiElectron = false; // use semiElectron channel?
  bool semiMuon = false; // use semiMuon channel?
	
	if(semiLepton)
	{
		semiMuon = true; //first start with semiMuon
		semiElectron = false;
	}
	
	if (argc >= 2)
	{	
	  semiMuon = atoi(argv[1]);
		semiElectron = !semiMuon;
	}
	
	vector<vector<pair<TH1F*,Dataset*> > > vec_AllBoxes_semiMu;
	vector<vector<pair<TH1F*,Dataset*> > > vec_AllBoxes_semiEl;
	vector<vector<pair<TH1F*,Dataset*> > > vec_AllBoxes_semiLep;
	vector<vector<pair<TH1F*,Dataset*> > > vec_AllDistributions_Nonfitted_semiMu;
	vector<vector<pair<TH1F*,Dataset*> > > vec_AllDistributions_Nonfitted_semiEl;
	vector<vector<pair<TH1F*,Dataset*> > > vec_AllDistributions_Nonfitted_semiLep;
	
	
	vector<vector<float> > QuadUncertSumPlusvec_AllBoxes_semiMu;
	vector<vector<float> > QuadUncertSumMinusvec_AllBoxes_semiMu;	
	vector<vector<float> > QuadUncertSumPlusvec_AllBoxes_semiEl;
	vector<vector<float> > QuadUncertSumMinusvec_AllBoxes_semiEl;
	vector<vector<float> > QuadUncertSumPlusvec_AllBoxes_semiLep;
	vector<vector<float> > QuadUncertSumMinusvec_AllBoxes_semiLep;
	
	string basedir = "1D_histograms/";
	const unsigned int nplots = 9;
	string boxdistributions[nplots] = {"HT_1B_2W","HT_1B_3W","nEvents_1B_4W","HT_2B_1W","HT_2B_2W","HT_2B_3W","nEvents_2B_4W","Mtop_MVA_1B_2W","Mtop_MVA_2B_2W"}; 
  string boxText[nplots] = {"1B_2W","1B_3W","1B_4W","2B_1W","2B_2W","2B_3W","2B_4W","1B_2W","2B_2W"}; 
	string basedir_Nonfitted = "MultiSamplePlot_allDiJetMasses/";
	const unsigned int nplots_Nonfitted = 1;
	string distributions_Nonfitted[nplots_Nonfitted] = {"allDiJetMasses_"};

	//make sure the 3 single top processes are the 6 last in the array!!! -> outdated?	
	//string datasetnamesSM[ndatasetnamesSM] = {TTbarJets_SemiLepton_datasetname,"TTbarJets_Other","WJets_HF","WJets_LF","ZJets","WW","WZ","ZZ","ttW","ttZ","samesignWWjj","ST_tChannel_t","ST_tWChannel_t","ST_sChannel_t","ST_tChannel_tbar","ST_tWChannel_tbar","ST_sChannel_tbar"};

	string axisnames[nplots] = {"S_{T} (GeV)","S_{T} (GeV)","","S_{T} (GeV)","S_{T} (GeV)","S_{T} (GeV)","","m_{bW} (GeV)","m_{bW} (GeV)"};
	string axisnames_Nonfitted[nplots_Nonfitted] = {"Dijet mass (GeV)"};
	
	const unsigned int nsystfiles = 11; //9
	string syst_filespostfix[nsystfiles] = {"Nominal","JESPlus","JESMinus","JERPlus","JERMinus","bTagPlus","bTagMinus","misTagPlus","misTagMinus","PUPlus","PUMinus"}; //"PUPlus","PUMinus"
//	string syst_filespostfix[nsystfiles] = {"Nominal"}; //"PUPlus","PUMinus"

	
	//string nominal_fileName = inputfiledir+"InclFourthGenSearch_TreeAnalyzer"+channelpostfix+"_Nominal"+".root";
	///cout << " Nominal filename = " << nominal_fileName << endl;
	//TFile* nominal_file = TFile::Open(nominal_fileName.c_str());

   vector < Dataset* > datasets;
	 vector < Dataset* > datasets_AfterMerging;
	 
while(semiElectron==true || semiMuon==true){
  	
	string TTbarJets_SemiLepton_datasetname = "";
	
  if(semiMuon){
       cout << " --> Using the semiMuon channel..." << endl;
       channelpostfix = "_semiMu";
			 TTbarJets_SemiLepton_datasetname = "TTbarJets_SemiMuon";
  }
  else if(semiElectron){
       cout << " --> Using the semiElectron channel..." << endl;
       channelpostfix = "_semiEl";
			 TTbarJets_SemiLepton_datasetname = "TTbarJets_SemiElectron";
  }

	
	//xml file
  string xmlFileName = "";
	if(semiElectron) xmlFileName = "../config/myFourthGenconfig_Electron_Fall11.xml";
  else if(semiMuon) xmlFileName = "../config/myFourthGenconfig_Muon_Fall11.xml";	
  const char *xmlfile = xmlFileName.c_str();
  cout << "used config file: " << xmlfile << endl;
	
	////////////////////////////////////
  /// AnalysisEnvironment  
  ////////////////////////////////////

  AnalysisEnvironment anaEnv;
  cout<<" - Loading environment ..."<<endl;
  AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);
  int verbose = anaEnv.Verbose;
  float anaEnvLuminosity = anaEnv.Luminosity;	// in 1/pb 
  cout << "analysis environment luminosity for rescaling "<< anaEnvLuminosity << endl;
	
	/////////////////////
  // Load Datasets
  /////////////////////

  TTreeLoader treeLoader;
  datasets.clear();
	datasets_AfterMerging.clear();
  cout << " - Loading datasets ..." << endl;
  treeLoader.LoadDatasets (datasets, xmlfile);
	
	
	vector<string> syst_fileNames;
	vector<TFile*> syst_files;
	for(unsigned int i=0;i<nsystfiles;i++)
	{
	  string syst_fileName = inputfiledir+"InclFourthGenSearch_TreeAnalyzer"+channelpostfix+"_"+syst_filespostfix[i]+".root";
		syst_fileNames.push_back(syst_fileName);
		cout << " Systematic filename = " << syst_fileName << endl;
		syst_files.push_back(TFile::Open(syst_fileName.c_str()));
	}

	//float nuisanceparameters = {Lumi,alpha_btag,alpha_jer,alpha_jes,};
	
//source: /user/gvonsem/CMSSW/CMSSW_4_2_8_patch4/src/42X_CommonCVS/TopBrussels/TopTreeAnalysis/macros/CalculateLimits/Apr2012_Systematics/RooStatsLimit/CLs_PBS/logs/SignalFix_combinedmuonelectron_floating_w_CKMA0_mass550_16042012_234658.txt
/*	float Lumi =  0.993358;
	float alpha_btag =  0.231355;
	float alpha_jer =  -0.363060;
	float alpha_jes   = 0.394299;
	float alpha_lepton_eff_el  = -0.735652;
	float alpha_lepton_eff_mu =  -1.39080;
	float alpha_mistag = -0.446058;
	float alpha_norm_ssww_syst =  0.260625;
	float alpha_norm_st_syst =  1.47604;
	float alpha_norm_top_syst =  0.144836;
	float alpha_norm_ttw_syst = -0.0564080;
	float alpha_norm_ttz_syst = -0.128374;
	float alpha_norm_ww_syst = -0.00984092;
	float alpha_norm_wz_syst =  0.0379692;
	float alpha_norm_z_syst =  0.0141265;
	float alpha_norm_zz_syst =  0.00286032;
	float norm_w_hf  =  1.0000;//note; W 'floating'!
	float norm_w_lf =   0.950000;//note; W 'floating'!
	*/
/*//source: /user/gvonsem/CMSSW/CMSSW_4_2_8_patch4/src/42X_CommonCVS/TopBrussels/TopTreeAnalysis/macros/CalculateLimits/Apr2012_Systematics/RooStatsLimit/2012_04_16_SignalFix_combinedmuonelectron_floating_w/CKMA0/mass550/hf_tprime.xml
	float rel_unc_mu_eff 	= 0.03;
	float rel_unc_el_eff 	= 0.05;
	float rel_unc_top = 0.08;
	float rel_unc_w_hf = 0.50;
  float rel_unc_w_lf = 0.05;
	float rel_unc_z = 0.05;
	float rel_unc_ww = 0.15;
	float rel_unc_wz = 0.17;
	float rel_unc_zz = 0.42;
	float rel_unc_ttw = 0.23;
	float rel_unc_ttz = 0.39;
	float rel_unc_ssww = 0.95;
	float rel_unc_st = 0.30;

//source: /user/gvonsem/CMSSW/CMSSW_4_2_8_patch4/src/42X_CommonCVS/TopBrussels/TopTreeAnalysis/macros/CalculateLimits/Apr2012_Systematics/RooStatsLimit/CLs_PBS/logs/combined_NuisParamTest_CKMA0_mass550_25042012_183212.txt
  float Lumi  =    0.966712;
  float alpha_btag  = 0.326854;
  float alpha_jer =  -0.252486;
	float alpha_jes =   0.420719;
	float alpha_lepton_eff_el = -1.20028;
	float alpha_lepton_eff_mu = -1.38267;
	float alpha_mistag = -1.21184;
  float alpha_norm_ssww_syst = -0.122225;
  float alpha_norm_st_syst =  1.61569;
  float alpha_norm_top_syst =  0.766222;
  float alpha_norm_ttw_syst =  0.134449;
  float alpha_norm_ttz_syst = -0.0574255;
  float alpha_norm_ww_syst = -0.0186962;
  float alpha_norm_wz_syst =  0.129456;
  float alpha_norm_z_syst =  0.000369715;
  float alpha_norm_zz_syst =  0.00123987;
  float alpha_unc_ch_e =  0.0522277;
  float alpha_unc_fl_e =  0.107090;
  float alpha_unc_fl_mu =  0.0522324;
  float norm_ttz  =   1.00000;
  float norm_w_hf  =  1.00000;
  float norm_w_lf  =  0.950000;
*/

//source: /user/pvmulder/NewEraOfDataAnalysis/TopTree/CMSSW_4_2_8_patch7/src/TopBrussels/TopTreeAnalysis/macros/CalculateLimits/RooStatsLimit/2012_05_10_combined/CKMA0/mass450/hf_tprime.xml
	float rel_unc_lumi = 0.022; //2.2%
	float rel_unc_mu_eff 	= 0.03;
	float rel_unc_el_eff 	= 0.05;
	float rel_unc_top = 0.12;
	float rel_unc_w = 0.50;
	float rel_unc_z = 0.05;
	float rel_unc_ww = 0.35;
	float rel_unc_wz = 0.42;
	float rel_unc_zz = 0.27;
	float rel_unc_ttw = 0.19;
	float rel_unc_ttz = 0.28;
	float rel_unc_ssww = 0.49;
	float rel_unc_st = 0.30;
//source: /user/pvmulder/NewEraOfDataAnalysis/TopTree/CMSSW_4_2_8_patch7/src/TopBrussels/TopTreeAnalysis/macros/CalculateLimits/RooStatsLimit/CLs_PBS/logs/combined_CKMA0_mass550_10052012_173351.txt
	float Lumi  =   1.0045;
  float alpha_btag  = 0.35;
  float alpha_jer =  -0.22;
	float alpha_jes =  0.44 ;
	float alpha_lepton_eff_el = -0.18;
	float alpha_lepton_eff_mu = 0.37;
	float alpha_mistag = -0.85;
  float alpha_norm_ssww_syst = 0.13;
  float alpha_norm_st_syst = 2.5;
  float alpha_norm_top_syst =  -0.33;
  float alpha_norm_ttw_syst =  0.12;
  float alpha_norm_ttz_syst = -0.04;
  float alpha_norm_ww_syst = -0.03;
  float alpha_norm_wz_syst = 0.03;
  float alpha_norm_z_syst =  0.003;
  float alpha_norm_zz_syst =  -0.002;
  float alpha_unc_ch_e =  0.04;
  float alpha_unc_fl_e =  0.07;
  float alpha_unc_fl_mu = 0.02;
  float alpha_norm_w  =  -0.18; //the level norm factor is actually norm_w: 0.91, so 9% down, so a 0.18sigma effect because the uncertainty is 50%
	float alpha_pu = 0.79;

  map<string,pair<float,float> > alpha_XS_systs;
	//hardcoded, and correspond to the SM datasetnames in the configs (they don't have to be put on add="1" btw), do not change this without a reason
	alpha_XS_systs["TTbarJets_SemiMuon"] = make_pair(rel_unc_top,alpha_norm_top_syst);	
	alpha_XS_systs["TTbarJets_SemiElectron"] = make_pair(rel_unc_top,alpha_norm_top_syst);
	alpha_XS_systs["TTbarJets_Other"] = make_pair(rel_unc_top,alpha_norm_top_syst);
	alpha_XS_systs["WJets1Jet"] = make_pair(rel_unc_w,alpha_norm_w); //note; W 'floating'!
	alpha_XS_systs["WJets2Jets"] = make_pair(rel_unc_w,alpha_norm_w); //note; W 'floating'!
	alpha_XS_systs["WJets3Jet"] = make_pair(rel_unc_w,alpha_norm_w); //note; W 'floating'!
	alpha_XS_systs["WJets4Jets"] = make_pair(rel_unc_w,alpha_norm_w); //note; W 'floating'!
	//alpha_XS_systs["WJets_HF"] = make_pair(rel_unc_w_hf,norm_w_hf); //note; W 'floating'!
	//alpha_XS_systs["WJets_LF"] = make_pair(rel_unc_w_lf,norm_w_lf); //note; W 'floating'!
	alpha_XS_systs["ZJets"] = make_pair(rel_unc_z,alpha_norm_z_syst);
	alpha_XS_systs["ST_tChannel_t"] = make_pair(rel_unc_st,alpha_norm_st_syst);
	alpha_XS_systs["ST_tChannel_tbar"] = make_pair(rel_unc_st,alpha_norm_st_syst);
	alpha_XS_systs["ST_sChannel_t"] = make_pair(rel_unc_st,alpha_norm_st_syst);
	alpha_XS_systs["ST_sChannel_tbar"] = make_pair(rel_unc_st,alpha_norm_st_syst);
	alpha_XS_systs["ST_tWChannel_t"] = make_pair(rel_unc_st,alpha_norm_st_syst);
	alpha_XS_systs["ST_tWChannel_tbar"] = make_pair(rel_unc_st,alpha_norm_st_syst);
	alpha_XS_systs["WW"] = make_pair(rel_unc_ww,alpha_norm_ww_syst);
	alpha_XS_systs["WZ"] = make_pair(rel_unc_wz,alpha_norm_wz_syst);
	alpha_XS_systs["ZZ"] = make_pair(rel_unc_zz,alpha_norm_zz_syst);
	alpha_XS_systs["ttW"] = make_pair(rel_unc_ttw,alpha_norm_ttw_syst);
	alpha_XS_systs["ttZ"] = make_pair(rel_unc_ttz,alpha_norm_ttz_syst);
	alpha_XS_systs["samesignWWjj"] = make_pair(rel_unc_ssww,alpha_norm_ssww_syst);
	
	
	///////////////////////////////////////
	// get histograms for all plots //
	///////////////////////////////////////	
	
	map<string,TH1F*** > boxes_processes_map; //second entry should become a matrix of TH1F*'s
	for(unsigned int s = 0; s < nsystfiles; s++)
	{
		//TH1F* boxes_processes[nplots][datasets.size()];
		TH1F*** boxes_processes = create2D_TH1F_Array(nplots,datasets.size());
		for(unsigned int box_i = 0; box_i < nplots; box_i++)
		{
			for(unsigned int d = 0; d<datasets.size(); d++)
			{
				string dataSetName = datasets[d]->Name();
				string histoname = basedir + boxdistributions[box_i].c_str() + dataSetName;
				
				if(dataSetName=="Data")
					boxes_processes[box_i][d] = (TH1F*) syst_files[0]->Get(histoname.c_str()); //element of matrix is filled; 0 = nominal (not very well coded...)
				else
				{
				  boxes_processes[box_i][d] = (TH1F*) syst_files[s]->Get(histoname.c_str()); //element of matrix is filled
					//cout<<(boxes_processes[box_i][d])->GetName()<<endl;
				}
			}
		}
		boxes_processes_map[syst_filespostfix[s]] = boxes_processes;
	}	
	cout << "retrieved histograms for processes for all plots (towards fitted plots)" << endl;
	
	
	//map<string,TH1F*** > distributions_Nonfitted_processes_map; //second entry should become a matrix of TH1F*'s
	//TH1F*** distributions_Nonfitted_processes = create2D_TH1F_Array(nplots_Nonfitted,datasets.size());
	for(unsigned int i = 0; i < nplots_Nonfitted; i++)
	{
		  vector<pair<TH1F*,Dataset*> > vec;
			TH1F* htemp_SignalA1 = 0;
		  TH1F* htemp_SignalA08 = 0;
		  Dataset* dataSet_SignalA1 = new Dataset("NP_overlay_SignalA1","Signal for A=1",true,1,1,2,1.,1.); //color and style overwritten in MultiSamplePlot
		  Dataset* dataSet_SignalA08 = new Dataset("NP_overlay_SignalA08","Signal for A=0.8",true,1,1,2,1.,1.); //color and style overwritten in MultiSamplePlot
		  pair<TH1F*,Dataset*> histodataset_pair_SignalA1;
		  pair<TH1F*,Dataset*> histodataset_pair_SignalA08;
		  float A;
		  bool firstoverlaysignalA1 = true, firstoverlaysignalA08 = true;
			for(unsigned int d = 0; d<datasets.size(); d++)
			{
				string dataSetName = datasets[d]->Name();
				string histoname = basedir_Nonfitted + distributions_Nonfitted[i].c_str() + dataSetName;				
				TH1F* htemp = (TH1F*) syst_files[0]->Get(histoname.c_str());
				if(dataSetName.find("NP_overlay")==0 && mergeSignal)
				{
			 	//cout<<dataSetName<<endl; 
				  A = 1;								
					if(dataSetName.find("NP_overlay_Tprime")<=0 && firstoverlaysignalA1)
					{ //first of the signal processes in the config...
					  htemp_SignalA1 = (TH1F*) htemp->Clone(); //clone is necessary! (otherwise htemp_SignalA1 and htemp_SignalA08 the same pointer eventually...)
						htemp_SignalA1->Scale(A);
						firstoverlaysignalA1 = false;
					}
			 		else
					{
					  if(!firstoverlaysignalA1)
						{
					  	if(dataSetName.find("NP_overlay_STprime")<=0)
							  htemp_SignalA1->Add(htemp,1.-A);
							if(dataSetName.find("NP_overlay_Bprime")<=0)
							  htemp_SignalA1->Add(htemp,A);
							if(dataSetName.find("NP_overlay_SBprime")<=0)
							  htemp_SignalA1->Add(htemp,1.-A);
							if(dataSetName.find("NP_overlay_BprimeTprime")<=0)
							  htemp_SignalA1->Add(htemp,A);
						}						
					}
					
				  A = 0.8;							
					if(dataSetName.find("NP_overlay_Tprime")<=0 && firstoverlaysignalA08)
					{ //first of the signal processes in the config...
					  htemp_SignalA08 = (TH1F*) htemp->Clone();//clone is necessary! (otherwise htemp_SignalA1 and htemp_SignalA08 the same pointer eventually...)
						htemp_SignalA1->Scale(A);
						firstoverlaysignalA08 = false;
					}
			 		else
					{
					  if(!firstoverlaysignalA08)
						{
					  	if(dataSetName.find("NP_overlay_STprime")<=0)
							  htemp_SignalA08->Add(htemp,1.-A);
							if(dataSetName.find("NP_overlay_Bprime")<=0)
							  htemp_SignalA08->Add(htemp,A);
							if(dataSetName.find("NP_overlay_SBprime")<=0)
							  htemp_SignalA08->Add(htemp,1.-A);
							if(dataSetName.find("NP_overlay_BprimeTprime")<=0)
							  htemp_SignalA08->Add(htemp,A);
						}						
					}					
				}
				else
				{
					pair<TH1F*,Dataset*> histodataset_pair;	
					histodataset_pair = make_pair((TH1F*) htemp->Clone(),datasets[d]);
					vec.push_back(histodataset_pair);
				}	
			}
			
			if(mergeSignal)	
			{
					histodataset_pair_SignalA1 = make_pair(htemp_SignalA1,dataSet_SignalA1);
					vec.push_back(histodataset_pair_SignalA1);
					datasets_AfterMerging.push_back(dataSet_SignalA1);
					histodataset_pair_SignalA08 = make_pair(htemp_SignalA08,dataSet_SignalA08);
					vec.push_back(histodataset_pair_SignalA08);
					datasets_AfterMerging.push_back(dataSet_SignalA08);
			}
	
			if(!semiLepton)
		  	MSPlot[distributions_Nonfitted[i]] = new MultiSamplePlot(vec,axisnames_Nonfitted[i],"#Events");

			if(semiLepton)
			{		
				if(semiMuon)
			  	vec_AllDistributions_Nonfitted_semiMu.push_back(vec);
				else if(semiElectron)
			  	vec_AllDistributions_Nonfitted_semiEl.push_back(vec);
			}
	}
  //distributions_Nonfitted_processes_map["Nominal"] = distributions_Nonfitted_processes;
	cout << "retrieved histograms for processes for all plots (towards non-fitted plots)" << endl;

	
	
	//loop over plots: remake each plot using fitted nuisance parameters
	for(unsigned int box_i=0; box_i<nplots; box_i++)
	{
		cout << "remake the " << boxdistributions[box_i] << " histogram" << endl;
	
		
		TH1F ** contribution_lumi = new TH1F*[datasets.size()];
		TH1F ** contribution_mueff = new TH1F*[datasets.size()];
		TH1F ** contribution_eleff = new TH1F*[datasets.size()];
		TH1F ** contribution_XS = new TH1F*[datasets.size()];
		TH1F ** contribution_JES = new TH1F*[datasets.size()];
		TH1F ** contribution_JER = new TH1F*[datasets.size()];
		TH1F ** contribution_misTag = new TH1F*[datasets.size()];
		TH1F ** contribution_bTag = new TH1F*[datasets.size()];
		TH1F ** contribution_PU = new TH1F*[datasets.size()];
		TH1F ** contribution_All = new TH1F*[datasets.size()];
		//
		TH1F ** contributionSigmaPlus_lumi = new TH1F*[datasets.size()];
		TH1F ** contributionSigmaPlus_mueff = new TH1F*[datasets.size()];
		TH1F ** contributionSigmaPlus_eleff = new TH1F*[datasets.size()];
		TH1F ** contributionSigmaPlus_XS = new TH1F*[datasets.size()];
		TH1F ** contributionSigmaPlus_JES = new TH1F*[datasets.size()];
		TH1F ** contributionSigmaPlus_JER = new TH1F*[datasets.size()];
		TH1F ** contributionSigmaPlus_misTag = new TH1F*[datasets.size()];
		TH1F ** contributionSigmaPlus_bTag = new TH1F*[datasets.size()];
		TH1F ** contributionSigmaPlus_PU = new TH1F*[datasets.size()];
		TH1F ** contributionSigmaPlus_All = new TH1F*[datasets.size()];
		//
		TH1F ** contributionSigmaMinus_lumi = new TH1F*[datasets.size()];
		TH1F ** contributionSigmaMinus_mueff = new TH1F*[datasets.size()];
		TH1F ** contributionSigmaMinus_eleff = new TH1F*[datasets.size()];
		TH1F ** contributionSigmaMinus_XS = new TH1F*[datasets.size()];
		TH1F ** contributionSigmaMinus_JES = new TH1F*[datasets.size()];
		TH1F ** contributionSigmaMinus_JER = new TH1F*[datasets.size()];
		TH1F ** contributionSigmaMinus_misTag = new TH1F*[datasets.size()];
		TH1F ** contributionSigmaMinus_bTag = new TH1F*[datasets.size()];
		TH1F ** contributionSigmaMinus_PU = new TH1F*[datasets.size()];
		TH1F ** contributionSigmaMinus_All = new TH1F*[datasets.size()];
		
		for (unsigned int d = 0; d < datasets.size(); d++) {
        contribution_lumi[d] = new TH1F();
				contribution_mueff[d] = new TH1F();
				contribution_eleff[d] = new TH1F();
				contribution_XS[d] = new TH1F();
				contribution_JES[d] = new TH1F();
				contribution_JER[d] = new TH1F();
				contribution_misTag[d] = new TH1F();
				contribution_bTag[d] = new TH1F();
				contribution_PU[d] = new TH1F();
				contribution_All[d] = new TH1F();
				//
				contributionSigmaPlus_lumi[d] = new TH1F();
				contributionSigmaPlus_mueff[d] = new TH1F();
				contributionSigmaPlus_eleff[d] = new TH1F();
				contributionSigmaPlus_XS[d] = new TH1F();
				contributionSigmaPlus_JES[d] = new TH1F();
				contributionSigmaPlus_JER[d] = new TH1F();
				contributionSigmaPlus_misTag[d] = new TH1F();
				contributionSigmaPlus_bTag[d] = new TH1F();
				contributionSigmaPlus_PU[d] = new TH1F();
				contributionSigmaPlus_All[d] = new TH1F();
				//
				contributionSigmaMinus_lumi[d] = new TH1F();
				contributionSigmaMinus_mueff[d] = new TH1F();
				contributionSigmaMinus_eleff[d] = new TH1F();
				contributionSigmaMinus_XS[d] = new TH1F();
				contributionSigmaMinus_JES[d] = new TH1F();
				contributionSigmaMinus_JER[d] = new TH1F();
				contributionSigmaMinus_misTag[d] = new TH1F();
				contributionSigmaMinus_bTag[d] = new TH1F();
				contributionSigmaMinus_PU[d] = new TH1F();
				contributionSigmaMinus_All[d] = new TH1F();							
    }
		
		vector<pair<TH1F*,Dataset*> > vec;
		TH1F* htemp_SignalA1 = 0;
		TH1F* htemp_SignalA08 = 0;
		Dataset* dataSet_SignalA1 = new Dataset("NP_overlay_SignalA1","Signal for A=1",true,1,1,2,1.,1.); //color and style overwritten in MultiSamplePlot
		Dataset* dataSet_SignalA08 = new Dataset("NP_overlay_SignalA08","Signal for A=0.8",true,1,1,2,1.,1.); //color and style overwritten in MultiSamplePlot
		pair<TH1F*,Dataset*> histodataset_pair_SignalA1;
		pair<TH1F*,Dataset*> histodataset_pair_SignalA08;
		float A;
		bool firstoverlaysignalA1 = true, firstoverlaysignalA08 = true;
		datasets_AfterMerging.clear();
		vector<float> QuadUncertSumPlus,QuadUncertSumMinus; //for each bin of the current plot an entry
		for(int bini=0; bini<((boxes_processes_map["Nominal"])[box_i][0])->GetNbinsX()+1; bini++)
		{	
				QuadUncertSumPlus.push_back(0);
				QuadUncertSumMinus.push_back(0);
		}
		
		for(unsigned int d = 0; d<datasets.size(); d++)
		{
			string dataSetName = datasets[d]->Name();	
			pair<TH1F*,Dataset*> histodataset_pair;
			if(dataSetName=="Data")
			{
			  //cout << "Data..." << endl;
			  TH1F * htemp_Data = (boxes_processes_map["Nominal"])[box_i][d];	
				histodataset_pair = make_pair((TH1F*) htemp_Data->Clone(),datasets[d]);
				vec.push_back(histodataset_pair);
				datasets_AfterMerging.push_back(datasets[d]);
				continue;			
			}
			
			if(dataSetName.find("NP")<=0 && !(dataSetName.find("NP_overlay")<=0))
				continue; //no use in looping over signal that is not overlayed
			
			TH1F * htemp_Nominal = (boxes_processes_map["Nominal"])[box_i][d];		
			
			
			//cout << "lumi contribution" << endl;				
			//cout<<htemp->GetName()<<endl;
			contribution_lumi[d] = (TH1F*) htemp_Nominal->Clone();
			contribution_lumi[d]->Scale(Lumi-1); //if lumi > 1 ===> contribution is positive //else contribution is negative												
			
			if(!(dataSetName.find("NP")<=0))
			{
				contributionSigmaPlus_lumi[d] = (TH1F*) htemp_Nominal->Clone();
				contributionSigmaPlus_lumi[d]->Scale(rel_unc_lumi);
				contributionSigmaMinus_lumi[d] = (TH1F*) htemp_Nominal->Clone();
				contributionSigmaMinus_lumi[d]->Scale(rel_unc_lumi);
				for(int bini=0; bini<htemp_Nominal->GetNbinsX()+1; bini++)
				{	
				  float BinErrorPlus = contributionSigmaPlus_lumi[d]->GetBinContent(bini);
					float BinErrorMinus = contributionSigmaMinus_lumi[d]->GetBinContent(bini);
					QuadUncertSumPlus[bini] = QuadUncertSumPlus[bini] + pow(BinErrorPlus,2);
					QuadUncertSumMinus[bini] = QuadUncertSumMinus[bini] + pow(BinErrorMinus,2);
				}
			}

			
			
			
			//cout << "lepton eff contribution" << endl;
			contribution_mueff[d] = (TH1F*) htemp_Nominal->Clone();
			contribution_mueff[d]->Scale(rel_unc_mu_eff*alpha_lepton_eff_mu);
			contribution_eleff[d] = (TH1F*) htemp_Nominal->Clone();
			contribution_eleff[d]->Scale(rel_unc_el_eff*alpha_lepton_eff_el);
			
			if(!(dataSetName.find("NP")<=0))
			{
				contributionSigmaPlus_mueff[d] = (TH1F*) htemp_Nominal->Clone();
				contributionSigmaPlus_mueff[d]->Scale(rel_unc_mu_eff);
				contributionSigmaMinus_mueff[d] = (TH1F*) htemp_Nominal->Clone();
				contributionSigmaMinus_mueff[d]->Scale(rel_unc_mu_eff);
				for(int bini=0; bini<htemp_Nominal->GetNbinsX()+1; bini++)
				{	
			  	float BinErrorPlus = 0;
					float BinErrorMinus = 0;
					if(semiMuon)
					{
						BinErrorPlus = contributionSigmaPlus_mueff[d]->GetBinContent(bini);
						BinErrorMinus = contributionSigmaMinus_mueff[d]->GetBinContent(bini);
					}
					if(semiElectron)
					{
						BinErrorPlus = contributionSigmaPlus_eleff[d]->GetBinContent(bini);
						BinErrorMinus = contributionSigmaMinus_eleff[d]->GetBinContent(bini);
					}
					QuadUncertSumPlus[bini] = QuadUncertSumPlus[bini] + pow(BinErrorPlus,2);
					QuadUncertSumMinus[bini] = QuadUncertSumMinus[bini] + pow(BinErrorMinus,2);
				}
			}
			
			
			
			if(!(dataSetName.find("NP")<=0))
			{
				//cout << "cross section contributions" << endl;
				contribution_XS[d] = (TH1F*) htemp_Nominal->Clone();
				contribution_XS[d]->Scale(alpha_XS_systs[dataSetName].first * alpha_XS_systs[dataSetName].second);
			
				contributionSigmaPlus_XS[d] = (TH1F*) htemp_Nominal->Clone();
				contributionSigmaPlus_XS[d]->Scale(alpha_XS_systs[dataSetName].first);
				contributionSigmaMinus_XS[d] = (TH1F*) htemp_Nominal->Clone();
				contributionSigmaMinus_XS[d]->Scale(alpha_XS_systs[dataSetName].first);
				for(int bini=0; bini<htemp_Nominal->GetNbinsX()+1; bini++)
				{	
			  	float BinErrorPlus = contributionSigmaPlus_XS[d]->GetBinContent(bini);
					float BinErrorMinus = contributionSigmaMinus_XS[d]->GetBinContent(bini);
					QuadUncertSumPlus[bini] = QuadUncertSumPlus[bini] + pow(BinErrorPlus,2);
					QuadUncertSumMinus[bini] = QuadUncertSumMinus[bini] + pow(BinErrorMinus,2);
				}
			
			}



			//cout << "JES contribution" << endl;
			TH1F* htemp_JES = 0;
			if(alpha_jes>0)
			{
				htemp_JES = (boxes_processes_map["JESPlus"])[box_i][d]; //not flexible enough, contradictory with the aim of some things in the code above...	
			}
			else
			{
				htemp_JES = (boxes_processes_map["JESMinus"])[box_i][d];
			}
			contribution_JES[d] = (TH1F*) htemp_JES->Clone();	
			contribution_JES[d]->Add(htemp_Nominal,-1);
			contribution_JES[d]->Scale(fabs(alpha_jes));
			
			if(!(dataSetName.find("NP")<=0))
			{
				contributionSigmaPlus_JES[d] = (TH1F*) ((boxes_processes_map["JESPlus"])[box_i][d])->Clone();
				contributionSigmaPlus_JES[d]->Add(htemp_Nominal,-1);
				contributionSigmaMinus_JES[d] = (TH1F*) ((boxes_processes_map["JESMinus"])[box_i][d])->Clone();
				contributionSigmaMinus_JES[d]->Add(htemp_Nominal,-1);
				for(int bini=0; bini<htemp_Nominal->GetNbinsX()+1; bini++)
				{	
			  	float BinErrorPlus = contributionSigmaPlus_JES[d]->GetBinContent(bini);
					float BinErrorMinus = contributionSigmaMinus_JES[d]->GetBinContent(bini);
					QuadUncertSumPlus[bini] = QuadUncertSumPlus[bini] + pow(BinErrorPlus,2);
					QuadUncertSumMinus[bini] = QuadUncertSumMinus[bini] + pow(BinErrorMinus,2);
				}
			}
			
			
			
			//cout << "JER contribution" << endl;
			TH1F* htemp_JER = 0;
			if(alpha_jer>0)
			{
				htemp_JER = (boxes_processes_map["JERPlus"])[box_i][d];
			}
			else
			{
				htemp_JER = (boxes_processes_map["JERMinus"])[box_i][d];
			}
			contribution_JER[d] = (TH1F*) htemp_JER->Clone();	
			contribution_JER[d]->Add(htemp_Nominal,-1);
			contribution_JER[d]->Scale(fabs(alpha_jer));
			
			if(!(dataSetName.find("NP")<=0))
			{
				contributionSigmaPlus_JER[d] = (TH1F*) ((boxes_processes_map["JERPlus"])[box_i][d])->Clone();
				contributionSigmaPlus_JER[d]->Add(htemp_Nominal,-1);
				contributionSigmaMinus_JER[d] = (TH1F*) ((boxes_processes_map["JERMinus"])[box_i][d])->Clone();
				contributionSigmaMinus_JER[d]->Add(htemp_Nominal,-1);
				for(int bini=0; bini<htemp_Nominal->GetNbinsX()+1; bini++)
				{	
				  float BinErrorPlus = contributionSigmaPlus_JER[d]->GetBinContent(bini);
					float BinErrorMinus = contributionSigmaMinus_JER[d]->GetBinContent(bini);
					QuadUncertSumPlus[bini] = QuadUncertSumPlus[bini] + pow(BinErrorPlus,2);
					QuadUncertSumMinus[bini] = QuadUncertSumMinus[bini] + pow(BinErrorMinus,2);
				}
			}		
			
			
				
			
			//cout << "misTag contribution" << endl;
			TH1F* htemp_misTag = 0;
			if(alpha_mistag>0)
			{
				htemp_misTag = (boxes_processes_map["misTagPlus"])[box_i][d];
			}
			else
			{
				htemp_misTag = (boxes_processes_map["misTagMinus"])[box_i][d];
			}
			contribution_misTag[d] = (TH1F*) htemp_misTag->Clone();	
			contribution_misTag[d]->Add(htemp_Nominal,-1);
			contribution_misTag[d]->Scale(fabs(alpha_mistag));
			
			if(!(dataSetName.find("NP")<=0))
			{
				contributionSigmaPlus_misTag[d] = (TH1F*) ((boxes_processes_map["misTagPlus"])[box_i][d])->Clone();
				contributionSigmaPlus_misTag[d]->Add(htemp_Nominal,-1);
				contributionSigmaMinus_misTag[d] = (TH1F*) ((boxes_processes_map["misTagMinus"])[box_i][d])->Clone();
				contributionSigmaMinus_misTag[d]->Add(htemp_Nominal,-1);
				for(int bini=0; bini<htemp_Nominal->GetNbinsX()+1; bini++)
				{	
			  	float BinErrorPlus = contributionSigmaPlus_misTag[d]->GetBinContent(bini);
					float BinErrorMinus = contributionSigmaMinus_misTag[d]->GetBinContent(bini);
					QuadUncertSumPlus[bini] = QuadUncertSumPlus[bini] + pow(BinErrorPlus,2);
					QuadUncertSumMinus[bini] = QuadUncertSumMinus[bini] + pow(BinErrorMinus,2);
				}
			}
			
			
			
			
			//cout << "bTag contribution" << endl;
			TH1F* htemp_bTag = 0;
			if(alpha_btag>0)
			{
				htemp_bTag = (boxes_processes_map["bTagPlus"])[box_i][d];
			}
			else
			{
				htemp_bTag = (boxes_processes_map["bTagMinus"])[box_i][d];
			}
			contribution_bTag[d] = (TH1F*) htemp_bTag->Clone();	
			contribution_bTag[d]->Add(htemp_Nominal,-1);
			contribution_bTag[d]->Scale(fabs(alpha_btag));
			
			if(!(dataSetName.find("NP")<=0))
			{
				contributionSigmaPlus_bTag[d] = (TH1F*) ((boxes_processes_map["bTagPlus"])[box_i][d])->Clone();
				contributionSigmaPlus_bTag[d]->Add(htemp_Nominal,-1);
				contributionSigmaMinus_bTag[d] = (TH1F*) ((boxes_processes_map["bTagMinus"])[box_i][d])->Clone();
				contributionSigmaMinus_bTag[d]->Add(htemp_Nominal,-1);
				for(int bini=0; bini<htemp_Nominal->GetNbinsX()+1; bini++)
				{	
				  float BinErrorPlus = contributionSigmaPlus_bTag[d]->GetBinContent(bini);
					float BinErrorMinus = contributionSigmaMinus_bTag[d]->GetBinContent(bini);
					QuadUncertSumPlus[bini] = QuadUncertSumPlus[bini] + pow(BinErrorPlus,2);
					QuadUncertSumMinus[bini] = QuadUncertSumMinus[bini] + pow(BinErrorMinus,2);
				}
			}
			
			
			
			
			//cout << "PU contribution" << endl;
			TH1F* htemp_PU = 0;
			if(alpha_pu>0)
			{
				htemp_PU = (boxes_processes_map["PUPlus"])[box_i][d];
			}
			else
			{
				htemp_PU = (boxes_processes_map["PUMinus"])[box_i][d];
			}
			contribution_PU[d] = (TH1F*) htemp_PU->Clone();	
			contribution_PU[d]->Add(htemp_Nominal,-1);
			contribution_PU[d]->Scale(fabs(alpha_pu));
			
			if(!(dataSetName.find("NP")<=0))
			{
				contributionSigmaPlus_PU[d] = (TH1F*) ((boxes_processes_map["PUPlus"])[box_i][d])->Clone();
				contributionSigmaPlus_PU[d]->Add(htemp_Nominal,-1);
				contributionSigmaMinus_PU[d] = (TH1F*) ((boxes_processes_map["PUMinus"])[box_i][d])->Clone();
				contributionSigmaMinus_PU[d]->Add(htemp_Nominal,-1);
				for(int bini=0; bini<htemp_Nominal->GetNbinsX()+1; bini++)
				{	
				  float BinErrorPlus = contributionSigmaPlus_PU[d]->GetBinContent(bini);
					float BinErrorMinus = contributionSigmaMinus_PU[d]->GetBinContent(bini);
					QuadUncertSumPlus[bini] = QuadUncertSumPlus[bini] + pow(BinErrorPlus,2);
					QuadUncertSumMinus[bini] = QuadUncertSumMinus[bini] + pow(BinErrorMinus,2);
				}
			}

			
			
			//cout << "add the contributions of all systematic effects to the nominal histogram" << endl;
			contribution_All[d] = (TH1F*) htemp_Nominal->Clone();
			contribution_All[d]->Add(contribution_lumi[d]);
			if(!(dataSetName.find("NP")<=0)) contribution_All[d]->Add(contribution_XS[d]);
			contribution_All[d]->Add(contribution_JES[d]);
			contribution_All[d]->Add(contribution_JER[d]);
			contribution_All[d]->Add(contribution_misTag[d]);
			contribution_All[d]->Add(contribution_bTag[d]);
			contribution_All[d]->Add(contribution_PU[d]);
			if(channelpostfix == "_semiMu") contribution_All[d]->Add(contribution_mueff[d]);
			if(channelpostfix == "_semiEl") contribution_All[d]->Add(contribution_eleff[d]);
			
			
			
			//now put them in a format for the MultiSamplePlot	
						
			if(dataSetName.find("NP_overlay")==0 && mergeSignal)
			{
			 	//cout<<dataSetName<<endl; 
				  A = 1;								
					if(dataSetName.find("NP_overlay_Tprime")<=0 && firstoverlaysignalA1)
					{ //first of the signal processes in the config...
					  htemp_SignalA1 = (TH1F*) contribution_All[d]->Clone(); //clone is necessary! (otherwise htemp_SignalA1 and htemp_SignalA08 the same pointer eventually...)
						htemp_SignalA1->Scale(A);
						firstoverlaysignalA1 = false;
					}
			 		else
					{
					  if(!firstoverlaysignalA1)
						{
					  	if(dataSetName.find("NP_overlay_STprime")<=0)
							  htemp_SignalA1->Add(contribution_All[d],1.-A);
							if(dataSetName.find("NP_overlay_Bprime")<=0)
							  htemp_SignalA1->Add(contribution_All[d],A);
							if(dataSetName.find("NP_overlay_SBprime")<=0)
							  htemp_SignalA1->Add(contribution_All[d],1.-A);
							if(dataSetName.find("NP_overlay_BprimeTprime")<=0)
							  htemp_SignalA1->Add(contribution_All[d],A);
						}						
					}
					
				  A = 0.8;							
					if(dataSetName.find("NP_overlay_Tprime")<=0 && firstoverlaysignalA08)
					{ //first of the signal processes in the config...
					  htemp_SignalA08 = (TH1F*) contribution_All[d]->Clone();//clone is necessary! (otherwise htemp_SignalA1 and htemp_SignalA08 the same pointer eventually...)
						htemp_SignalA1->Scale(A);
						firstoverlaysignalA08 = false;
					}
			 		else
					{
					  if(!firstoverlaysignalA08)
						{
					  	if(dataSetName.find("NP_overlay_STprime")<=0)
							  htemp_SignalA08->Add(contribution_All[d],1.-A);
							if(dataSetName.find("NP_overlay_Bprime")<=0)
							  htemp_SignalA08->Add(contribution_All[d],A);
							if(dataSetName.find("NP_overlay_SBprime")<=0)
							  htemp_SignalA08->Add(contribution_All[d],1.-A);
							if(dataSetName.find("NP_overlay_BprimeTprime")<=0)
							  htemp_SignalA08->Add(contribution_All[d],A);
						}						
					}					
			}
			else
			{		
				histodataset_pair = make_pair(contribution_All[d],datasets[d]);
				vec.push_back(histodataset_pair);
				datasets_AfterMerging.push_back(datasets[d]);
			}
		
		}
		
		for(int bini=0; bini<((boxes_processes_map["Nominal"])[box_i][0])->GetNbinsX()+1; bini++)
		{
				QuadUncertSumPlus[bini] = sqrt(QuadUncertSumPlus[bini]);
				QuadUncertSumMinus[bini] = sqrt(QuadUncertSumMinus[bini]);
			  
		}
		//now, for the current channel and box, and after the dataset loop, we should have two vectors (of the bin errors up en down, the vector entries are for different bins of the current histo) to be turned in the error band...
	
		
		if(mergeSignal)	
		{
			histodataset_pair_SignalA1 = make_pair(htemp_SignalA1,dataSet_SignalA1);
			vec.push_back(histodataset_pair_SignalA1);
			datasets_AfterMerging.push_back(dataSet_SignalA1);
			histodataset_pair_SignalA08 = make_pair(htemp_SignalA08,dataSet_SignalA08);
			vec.push_back(histodataset_pair_SignalA08);
			datasets_AfterMerging.push_back(dataSet_SignalA08);
		}
		
		if(!semiLepton)
		  MSPlot[boxdistributions[box_i]] = new MultiSamplePlot(vec,axisnames[box_i],"#Events","Fitted nuisance parameters applied"); //last argument is text which is written on the canvas

		if(semiLepton)
		{		
			if(semiMuon)
			{
			  vec_AllBoxes_semiMu.push_back(vec);
				QuadUncertSumPlusvec_AllBoxes_semiMu.push_back(QuadUncertSumPlus); //as much entries as boxes
				QuadUncertSumMinusvec_AllBoxes_semiMu.push_back(QuadUncertSumMinus);
			}
			else if(semiElectron)
			{
			  vec_AllBoxes_semiEl.push_back(vec);
				QuadUncertSumPlusvec_AllBoxes_semiEl.push_back(QuadUncertSumPlus); //as much entries as boxes
				QuadUncertSumMinusvec_AllBoxes_semiEl.push_back(QuadUncertSumMinus);
			}
		}
	}
	
	
	/*//just to test
	if(semiMuon)
	{
		cout<<"for muon channel:"<<endl;
		for(int m=0;m<QuadUncertSumPlusvec_AllBoxes_semiMu.size();m++)
		{
			cout<<"* Histogram "<<m<<endl;
			for(int k=0;k<QuadUncertSumPlusvec_AllBoxes_semiMu[m].size();k++)
			{
				cout<<"  --> error for bin k = "<<k<<": "<<QuadUncertSumPlusvec_AllBoxes_semiMu[m][k]<<endl;
			}
		}
	}*/	
	
		
	
	if(!semiLepton)
	   break;
	else
	{	   
		if(semiElectron==false)
		{
			semiMuon = false;
		  semiElectron = true;
		}
		else
		{
			semiMuon = false;
		  semiElectron = false;
		}
	}
}

	if(semiLepton)
	{
		for(int m=0;m<QuadUncertSumPlusvec_AllBoxes_semiMu.size();m++)
		{
		  
		  cout<<"* Histogram "<<m<<endl;
			vector<float> dummyvect;
			for(int k=0;k<QuadUncertSumPlusvec_AllBoxes_semiMu[m].size();k++)
			{
				dummyvect.push_back(sqrt(pow(QuadUncertSumPlusvec_AllBoxes_semiMu[m][k],2) + pow(QuadUncertSumPlusvec_AllBoxes_semiEl[m][k],2)));
				cout<<"  --> total systematic error for bin k = "<<k<<": "<<dummyvect[k]<<endl;
			}
			QuadUncertSumPlusvec_AllBoxes_semiLep.push_back(dummyvect);
		}
		
	}


 if(semiLepton)
 {  
	for(unsigned int box_i=0; box_i<nplots; box_i++)
	{
		//cout<<"box_i = "<<box_i<<endl;
	  vector<pair<TH1F*,Dataset*> > vec;
		
		//data; do some dirty stuff to order the datasets differently... the data still has to go first, but then I would like to revers the order
		string dataSetName = datasets_AfterMerging[0]->Name();
		//cout<<"  dataSetName "<<0<<" = "<<dataSetName<<endl;
		pair<TH1F*,Dataset*> histodataset_pair;
		TH1F* htemp = (TH1F*) (((vec_AllBoxes_semiMu[box_i])[0]).first)->Clone();
		TH1F* htemp_semiEl = (TH1F*) (((vec_AllBoxes_semiEl[box_i])[0]).first)->Clone();
		htemp->Add(htemp_semiEl);
		//cout<<"  htemp->Integral() = "<<htemp->Integral()<<endl;
		histodataset_pair = make_pair(htemp,datasets_AfterMerging[0]);
		vec.push_back(histodataset_pair);
	  
		for(unsigned int d = datasets_AfterMerging.size()-1; d>0; d--)
		{
			string dataSetName = datasets_AfterMerging[d]->Name();
			//cout<<"  dataSetName "<<d<<" = "<<dataSetName<<endl;
			pair<TH1F*,Dataset*> histodataset_pair;
			TH1F* htemp = (TH1F*) (((vec_AllBoxes_semiMu[box_i])[d]).first)->Clone();
			TH1F* htemp_semiEl = (TH1F*) (((vec_AllBoxes_semiEl[box_i])[d]).first)->Clone();
			htemp->Add(htemp_semiEl);
			//cout<<"  htemp->Integral() = "<<htemp->Integral()<<endl;
			histodataset_pair = make_pair(htemp,datasets_AfterMerging[d]);
			vec.push_back(histodataset_pair);
		}

		
		//MSPlot[boxdistributions[box_i]] = new MultiSamplePlot(vec,axisnames[box_i],"Events","Fitted nuisance parameters applied"); //last argument is text which is written on the canvas
		string text = "lepton+jets, " + boxText[box_i];
		MSPlot[boxdistributions[box_i]] = new MultiSamplePlot(vec,axisnames[box_i],"Events",text,false); //second last argument is text which is written on the canvas
	}
	
	for(unsigned int i=0; i<nplots_Nonfitted; i++)
	{
	  vector<pair<TH1F*,Dataset*> > vec;
		
		//data; do some dirty stuff to order the datasets differently... the data still has to go first, but then I would like to revers the order
		string dataSetName = datasets_AfterMerging[0]->Name();
		//cout<<"  dataSetName "<<0<<" = "<<dataSetName<<endl;
		pair<TH1F*,Dataset*> histodataset_pair;
		TH1F* htemp = (TH1F*) (((vec_AllDistributions_Nonfitted_semiMu[i])[0]).first)->Clone();
		TH1F* htemp_semiEl = (TH1F*) (((vec_AllDistributions_Nonfitted_semiEl[i])[0]).first)->Clone();
		htemp->Add(htemp_semiEl);
		//cout<<"  htemp->Integral() = "<<htemp->Integral()<<endl;
		histodataset_pair = make_pair(htemp,datasets_AfterMerging[0]);
		vec.push_back(histodataset_pair);
	  for(unsigned int d = datasets_AfterMerging.size()-1; d>0; d--)
		{
			string dataSetName = datasets_AfterMerging[d]->Name();
			//cout<<"dataSetName "<<d<<" = "<<dataSetName<<endl;
			pair<TH1F*,Dataset*> histodataset_pair;
			TH1F* htemp = (TH1F*) (((vec_AllDistributions_Nonfitted_semiMu[i])[d]).first)->Clone();
			TH1F* htemp_semiEl = (TH1F*) (((vec_AllDistributions_Nonfitted_semiEl[i])[d]).first)->Clone();
			htemp->Add(htemp_semiEl);
			histodataset_pair = make_pair(htemp,datasets_AfterMerging[d]);
			vec.push_back(histodataset_pair);
		}
		//do some dirty stuff to order the datasets differently..
		
		string text = "lepton+jets";
		MSPlot[distributions_Nonfitted[i]] = new MultiSamplePlot(vec,axisnames_Nonfitted[i],"Events",text,false); //second last argument is text which is written on the canvas
	}
	
	channelpostfix = "_semiLep";
 }
	 
	//do some dirty stuff to order the datasets differently...
	
	 
	 
	//Output ROOT file
  string rootFileName (Outputpath+"FittedModel"+channelpostfix+".root");
  TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");
	fout->cd();
	bool savePDF = true;
	string pathPDF (Outputpath+"FittedModelPDF"+channelpostfix);
  pathPDF = pathPDF +"/"; 		
  if(savePDF) mkdir(pathPDF.c_str(),0777);
	for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
  {
        MultiSamplePlot *temp = it->second;
				string name;
        if(it->first != "allDiJetMasses_") name = it->first + "_Fitted" +channelpostfix; //dirty hardcoded
				else name = it->first +channelpostfix;
        //temp->Draw(false, name, true, true, true, true, true,5,false, true, true); //old way
				if( (it->first).find("Mtop") <= (it->first).size() )
				  temp->Draw(false, name, true, false, true, false, true,8,false, false, false, false, true);//(bool addRandomPseudoData, string label, bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST,int scaleNPsignal, bool addRatio, bool mergeVV, bool mergeTTV, bool mergeOtherTTJetsWJetsST) // previously: temp->Draw(false, name, true, true, true, true, true,5,false, true,true)
  			else
				  temp->Draw(false, name, true, false, true, false, true,8,false, false, false, true, false);
				temp->Write(fout, name, savePDF, pathPDF,"pdf");//bool savePDF     , pathPDF+"MSPlot/"
	}

}

