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
	
	cout<<"The arguments passed to the executable are: "<<endl;
  for(unsigned int i=1;i<argc;i++)
	{
		cout<<argv[i]<<endl;
	}
	
	string channelpostfix = "";
	bool semiElectron = false; // use semiElectron channel?
  bool semiMuon = true; // use semiMuon channel?
	string TTbarJets_SemiLepton_datasetname = "";
	if (argc >= 2)
	{	
	  semiMuon = atoi(argv[1]);
		semiElectron = !semiMuon;
	}
	
	if(semiElectron && semiMuon)
  {
     cout << "  --> Using both semiMuon and semiElectron channel? Choose only one!" << endl;
     exit(1);
  }
  else
  {
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
  }
	
	////////////////////////////////////
	// define and open the inputfiles //
	////////////////////////////////////
	string Outputpath = "Plots_FittedModel/";
	gSystem->mkdir(Outputpath.c_str());
	
	string inputfiledir = "CalculateLimits/Apr2012_Systematics/INPUT_2012_04_16/";
	
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
  vector < Dataset* > datasets;
  cout << " - Loading datasets ..." << endl;
  treeLoader.LoadDatasets (datasets, xmlfile);
	
	
	
	string basedir = "1D_histograms/";
	const unsigned int nplots = 10;
	string boxdistributions[nplots] = {"HT_1B_1W","HT_1B_2W","HT_1B_3W","nEvents_1B_4W","HT_2B_1W","HT_2B_2W","HT_2B_3W","nEvents_2B_4W","Mtop_MVA_1B_2W","Mtop_MVA_2B_2W"};

	//make sure the 3 single top processes are the 6 last in the array!!! -> outdated?	
	//string datasetnamesSM[ndatasetnamesSM] = {TTbarJets_SemiLepton_datasetname,"TTbarJets_Other","WJets_HF","WJets_LF","ZJets","WW","WZ","ZZ","ttW","ttZ","samesignWWjj","ST_tChannel_t","ST_tWChannel_t","ST_sChannel_t","ST_tChannel_tbar","ST_tWChannel_tbar","ST_sChannel_tbar"};

	string axisnames[nplots] = {"H_{T} (GeV)","H_{T} (GeV)","H_{T} (GeV)","","H_{T} (GeV)","H_{T} (GeV)","H_{T} (GeV)","","m_{bW} (GeV)","m_{bW} (GeV)"};
	
	const unsigned int nsystfiles = 9;
	string syst_filespostfix[nsystfiles] = {"Nominal","JESPlus","JESMinus","JERPlus","JERMinus","bTagPlus","bTagMinus","misTagPlus","misTagMinus"}; //"PUPlus","PUMinus"
	
	//string nominal_fileName = inputfiledir+"InclFourthGenSearch_TreeAnalyzer"+channelpostfix+"_Nominal"+".root";
	///cout << " Nominal filename = " << nominal_fileName << endl;
	//TFile* nominal_file = TFile::Open(nominal_fileName.c_str());
	
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
	float Lumi =  0.993358;
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
//source: /user/gvonsem/CMSSW/CMSSW_4_2_8_patch4/src/42X_CommonCVS/TopBrussels/TopTreeAnalysis/macros/CalculateLimits/Apr2012_Systematics/RooStatsLimit/2012_04_16_SignalFix_combinedmuonelectron_floating_w/CKMA0/mass550/hf_tprime.xml
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

	
  map<string,pair<float,float> > alpha_XS_systs;
	//hardcoded, and correspond to the SM datasetnames in the configs (they don't have to be put on add="1" btw), do not change this without a reason
	alpha_XS_systs["TTbarJets_SemiMuon"] = make_pair(rel_unc_top,alpha_norm_top_syst);	
	alpha_XS_systs["TTbarJets_SemiElectron"] = make_pair(rel_unc_top,alpha_norm_top_syst);
	alpha_XS_systs["TTbarJets_Other"] = make_pair(rel_unc_top,alpha_norm_top_syst);
	alpha_XS_systs["WJets_HF"] = make_pair(rel_unc_w_hf,norm_w_hf); //note; W 'floating'!
	alpha_XS_systs["WJets_LF"] = make_pair(rel_unc_w_lf,norm_w_lf); //note; W 'floating'!
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
	cout << "retrieved histograms for processes for all plots" << endl;
	
	
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
		TH1F ** contribution_All = new TH1F*[datasets.size()];
		for (unsigned int d = 0; d < datasets.size(); d++) {
        contribution_lumi[d] = new TH1F();
				contribution_mueff[d] = new TH1F();
				contribution_eleff[d] = new TH1F();
				contribution_XS[d] = new TH1F();
				contribution_JES[d] = new TH1F();
				contribution_JER[d] = new TH1F();
				contribution_misTag[d] = new TH1F();
				contribution_bTag[d] = new TH1F();
				contribution_All[d] = new TH1F();								
    }
		
		vector<pair<TH1F*,Dataset*> > vec;
		for(unsigned int d = 0; d<datasets.size(); d++)
		{
			string dataSetName = datasets[d]->Name();	
			pair<TH1F*,Dataset*> histodataset_pair;
			if(dataSetName=="Data")
			{
			  cout << "Data..." << endl;
			  TH1F * htemp_Data = (boxes_processes_map["Nominal"])[box_i][d];	
				histodataset_pair = make_pair((TH1F*) htemp_Data->Clone(),datasets[d]);
				vec.push_back(histodataset_pair);
				continue;			
			}
			
			if(dataSetName.find("NP")<=0 && !(dataSetName.find("NP_overlay")<=0))
				continue; //no use in looping over signal that is not overlayed
			
			cout << "lumi contribution" << endl;
			TH1F * htemp_Nominal = (boxes_processes_map["Nominal"])[box_i][d];					
			//cout<<htemp->GetName()<<endl;
			contribution_lumi[d] = (TH1F*) htemp_Nominal->Clone();
			contribution_lumi[d]->Scale(Lumi-1); //if lumi > 1 ===> contribution is positive //else contribution is negative
			
			cout << "lepton eff contribution" << endl;
			contribution_mueff[d] = (TH1F*) htemp_Nominal->Clone();
			contribution_mueff[d]->Scale(rel_unc_mu_eff*alpha_lepton_eff_mu);
			contribution_eleff[d] = (TH1F*) htemp_Nominal->Clone();
			contribution_eleff[d]->Scale(rel_unc_el_eff*alpha_lepton_eff_el);
			
			if(!(dataSetName.find("NP")<=0))
			{
				cout << "cross section contributions" << endl;
				contribution_XS[d] = (TH1F*) htemp_Nominal->Clone();
				contribution_XS[d]->Scale(alpha_XS_systs[dataSetName].first * alpha_XS_systs[dataSetName].second);
			}
			
			//cout << "PU contribution" << endl;
			
			cout << "JES contribution" << endl;
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
			
			cout << "JER contribution" << endl;
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
			
			cout << "misTag contribution" << endl;
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
			
			cout << "bTag contribution" << endl;
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
			
			cout << "add the contributions of all systematic effects to the nominal histogram" << endl;
			contribution_All[d] = (TH1F*) htemp_Nominal->Clone();
			contribution_All[d]->Add(contribution_lumi[d]);
			if(!(dataSetName.find("NP")<=0)) contribution_All[d]->Add(contribution_XS[d]);
			contribution_All[d]->Add(contribution_JES[d]);
			contribution_All[d]->Add(contribution_JER[d]);
			contribution_All[d]->Add(contribution_misTag[d]);
			contribution_All[d]->Add(contribution_bTag[d]);
			if(channelpostfix == "_semiMu") contribution_All[d]->Add(contribution_mueff[d]);
			if(channelpostfix == "_semiEl") contribution_All[d]->Add(contribution_eleff[d]);
			
			//no put them in a format for the MultiSamplePlot
			histodataset_pair = make_pair(contribution_All[d],datasets[d]);
			vec.push_back(histodataset_pair);
		}
		
		MSPlot[boxdistributions[box_i]] = new MultiSamplePlot(vec,axisnames[box_i],"#Events","Fitted nuisance parameters applied"); //last argument is text which is written on the canvas


	}
	
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
        string name = it->first + "_Fitted"+channelpostfix;
        temp->Draw(false, name, true, true, true, true, true,5,false, true, true);//(bool addRandomPseudoData, string label, bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST,int scaleNPsignal, bool addRatio, bool mergeVV, bool mergeTTV)
        temp->Write(fout, name, savePDF, pathPDF,"pdf");//bool savePDF     , pathPDF+"MSPlot/"
  }

}

