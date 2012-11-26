#include "TStyle.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
/*
// Root headers
#include "TArrow.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TMarker.h"
#include "TPaveStats.h"
#include "TRandom3.h"
*/
//user code
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "../../Content/interface/Dataset.h"
#include "../../Tools/interface/MultiSamplePlot.h"
#include "../interface/TopFCNC_Evt.h"
#include "../interface/TopFCNC_KinFit.h"

#include "../../macros/Style.C"

using namespace std;
using namespace TopTree;

// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;
// MultiSamplePlot
map<string,MultiSamplePlot*> MSPlot;

struct HighestCVSBtag{
    bool operator()( TRootJet j1, TRootJet j2 ) const{
    	return j1.btag_combinedSecondaryVertexBJetTags() > j2.btag_combinedSecondaryVertexBJetTags();
    }
};
/*https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagPerformanceOP
  Tagger name  	            WP name WP Discr cut
  TrackCountingHighPur 	    TCHPT 	3.41
  JetProbability 	          JPL 	  0.275
  JetProbability 	          JPM 	  0.545
  JetProbability 	          JPT 	  0.790
  CombinedSecondaryVertex 	CSVL 	  0.244
  CombinedSecondaryVertex 	CSVM 	  0.679
  CombinedSecondaryVertex 	CSVT 	  0.898
*/

int main (int argc, char *argv[])
{

  clock_t start = clock();

  cout << "*************************************************************" << endl;
  cout << " Beginning of the program for the FCNC analysis ! " << endl;
  cout << "*************************************************************" << endl;

  //SetStyle if needed
  //setTDRStyle();
  setGregStyle();
  //setMyStyle();


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////// Configuration ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Choose leptonic channel
  Bool_t UseMuChannel = true;

  //Output ROOT file
  string postfix = "_Analysis";
  string channelpostfix = (UseMuChannel ? "_DiMuTrigger" : "_DiElecTrigger");
  string comments = "_Run2012A";
  string rootFileName ("TopFCNC"+postfix+channelpostfix+comments+".root");
  string treepath = "$HOME/AnalysisCode/CMSSW_53X/TopBrussels/TopTreeAnalysis/TopFCNC/macros/TopFCNC_EventSelection_DiMuTrigger_Run2012A_Trees/";
  string resopath = "$HOME/AnalysisCode/CMSSW_53X/TopBrussels/TopTreeAnalysis/TopFCNC/macros/ResolutionFiles/";
  string resoprefix = "TTbar_FCNC_";

  Float_t Luminosity = -1.;
  Float_t EventWeight = 1.;

  vector<Dataset*> dataSets; // needed for MSPlots
  vector<string>   inputFiles;
  
  inputFiles.push_back(treepath+"TopFCNC_EventSelection"+channelpostfix+comments+"_TTree_Data.root");
/*
  inputFiles.push_back(treepath+"TopFCNC_EventSelection"+channelpostfix+comments+"_TTree_ST_tbar_tWch_DR.root");
  inputFiles.push_back(treepath+"TopFCNC_EventSelection"+channelpostfix+comments+"_TTree_ST_t_tWch_DR.root");
  inputFiles.push_back(treepath+"TopFCNC_EventSelection"+channelpostfix+comments+"_TTree_ST_tbar_tch.root");
  inputFiles.push_back(treepath+"TopFCNC_EventSelection"+channelpostfix+comments+"_TTree_ST_t_tch.root");
  inputFiles.push_back(treepath+"TopFCNC_EventSelection"+channelpostfix+comments+"_TTree_WWJetsTo2L2Nu.root");
  inputFiles.push_back(treepath+"TopFCNC_EventSelection"+channelpostfix+comments+"_TTree_WZJetsTo2L2Q.root");
  inputFiles.push_back(treepath+"TopFCNC_EventSelection"+channelpostfix+comments+"_TTree_WZJetsTo3LNu.root");
  inputFiles.push_back(treepath+"TopFCNC_EventSelection"+channelpostfix+comments+"_TTree_ZZJetsTo2L2Nu.root");
  inputFiles.push_back(treepath+"TopFCNC_EventSelection"+channelpostfix+comments+"_TTree_ZZJetsTo2L2Q.root");
  inputFiles.push_back(treepath+"TopFCNC_EventSelection"+channelpostfix+comments+"_TTree_ZZJetsTo4L.root");
*/
  inputFiles.push_back(treepath+"TopFCNC_EventSelection"+channelpostfix+comments+"_TTree_ttjets.root");
/*
  inputFiles.push_back(treepath+"TopFCNC_EventSelection"+channelpostfix+comments+"_TTree_zjets.root");
//  inputFiles.push_back(treepath+"TopFCNC_EventSelection"+channelpostfix+comments+"_TTree_ttbar_fcnc.root");
*/
  for(unsigned int iDataSet=0; iDataSet<inputFiles.size(); iDataSet++)
  {
    TFile* inputFile = new TFile(inputFiles[iDataSet].c_str(),"READ");
    TTree* inConfigTree = (TTree*) inputFile->Get("configTree");
    TBranch* d_br = (TBranch*) inConfigTree->GetBranch("Dataset");
    TClonesArray* tc_dataset = new TClonesArray("Dataset",0);
    d_br->SetAddress(&tc_dataset);
    inConfigTree->GetEvent(0);
    Dataset* dataSet = (Dataset*) tc_dataset->At(0);
    dataSets.push_back( (Dataset*) dataSet->Clone() );
    if(dataSet->Name().find("Data") == 0 || dataSet->Name().find("data") == 0 || dataSet->Name().find("DATA") == 0 )
      Luminosity = dataSet->EquivalentLumi();
    delete tc_dataset;
    inputFile->Close();
    delete inputFile;
  }
  
  if( Luminosity<0 ) Luminosity = 5000.;
  cout<<"Executing analysis for an integrated luminosity of " << Luminosity << " pb^-1" << endl;
    

  TFile      *fout  = new TFile (rootFileName.c_str(), "RECREATE");
  TDirectory *myDir = 0;

  ResolutionFit *resFitLeptons = 0;
  if(UseMuChannel){
    resFitLeptons = new ResolutionFit("Muon");
    resFitLeptons->LoadResolutions(resopath+resoprefix+"muonReso.root");
  }
  else{
    resFitLeptons = new ResolutionFit("Electron");
    resFitLeptons->LoadResolutions(resopath+resoprefix+"electronReso.root");
  }
  
  ResolutionFit *resFitBJets = new ResolutionFit("BJet");
  resFitBJets->LoadResolutions(resopath+resoprefix+"bJetReso.root");

  ResolutionFit *resFitQJets = new ResolutionFit("QJet");
  resFitQJets->LoadResolutions(resopath+resoprefix+"qJetReso.root");

  ResolutionFit *resFitLightJets = new ResolutionFit("LightJet");
  resFitLightJets->LoadResolutions(resopath+resoprefix+"lightJetReso.root");
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////// Histograms /////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////
  ///////////////////// MS plots /////////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  MSPlot["KinFit_Prob"]        = new MultiSamplePlot(dataSets,"KinFit_Prob",100,0,1,"Prob.");
  MSPlot["KinFit_Chi2"]        = new MultiSamplePlot(dataSets,"KinFit_Chi2",500,0,100,"#chi^{2}");
  MSPlot["KinFit_ReducedChi2"] = new MultiSamplePlot(dataSets,"KinFit_ReducedChi2",100,0,1,"#chi^{2}/Ndf");

  MSPlot["KinFit_HadWMass"]    = new MultiSamplePlot(dataSets,"KinFit_HadWMass",120,50,110,"m_{W} [Gev/c^{2}]");
  MSPlot["KinFit_HadTopMass"]  = new MultiSamplePlot(dataSets,"KinFit_HadTopMass",280,100,240,"m^{SM}_{top} [Gev/c^{2}]");
  MSPlot["KinFit_FcncTopMass"] = new MultiSamplePlot(dataSets,"KinFit_FcncTopMass",280,100,240,"m^{FCNC}_{top} [Gev/c^{2}]");

  MSPlot["KinFit_Prob_AtLeast1Btag_CSVM"]        = new MultiSamplePlot(dataSets,"KinFit_Prob_AtLeast1Btag_CSVM",100,0,1,"Prob.");
  MSPlot["KinFit_Chi2_AtLeast1Btag_CSVM"]        = new MultiSamplePlot(dataSets,"KinFit_Chi2_AtLeast1Btag_CSVM",500,0,100,"#chi^{2}");
  MSPlot["KinFit_ReducedChi2_AtLeast1Btag_CSVM"] = new MultiSamplePlot(dataSets,"KinFit_ReducedChi2_AtLeast1Btag_CSVM",100,0,1,"#chi^{2}/Ndf");

  MSPlot["KinFit_HadWMass_AtLeast1Btag_CSVM"]    = new MultiSamplePlot(dataSets,"KinFit_HadWMass_AtLeast1Btag_CSVM",120,50,110,"m_{W} [Gev/c^{2}]");
  MSPlot["KinFit_HadTopMass_AtLeast1Btag_CSVM"]  = new MultiSamplePlot(dataSets,"KinFit_HadTopMass_AtLeast1Btag_CSVM",280,100,240,"m^{SM}_{top} [Gev/c^{2}]");
  MSPlot["KinFit_FcncTopMass_AtLeast1Btag_CSVM"] = new MultiSamplePlot(dataSets,"KinFit_FcncTopMass_AtLeast1Btag_CSVM",280,100,240,"m^{FCNC}_{top} [Gev/c^{2}]");

  MSPlot["KinFit_Prob_NoBtag_CVSL"]        = new MultiSamplePlot(dataSets,"KinFit_Prob_NoBtag_CVSL",100,0,1,"Prob.");
  MSPlot["KinFit_Chi2_NoBtag_CVSL"]        = new MultiSamplePlot(dataSets,"KinFit_Chi2_NoBtag_CVSL",500,0,100,"#chi^{2}");
  MSPlot["KinFit_ReducedChi2_NoBtag_CVSL"] = new MultiSamplePlot(dataSets,"KinFit_ReducedChi2_NoBtag_CVSL",100,0,1,"#chi^{2}/Ndf");

  MSPlot["KinFit_HadWMass_NoBtag_CVSL"]    = new MultiSamplePlot(dataSets,"KinFit_HadWMass_NoBtag_CVSL",120,50,110,"m_{W} [Gev/c^{2}]");
  MSPlot["KinFit_HadTopMass_NoBtag_CVSL"]  = new MultiSamplePlot(dataSets,"KinFit_HadTopMass_NoBtag_CVSL",280,100,240,"m^{SM}_{top} [Gev/c^{2}]");
  MSPlot["KinFit_FcncTopMass_NoBtag_CVSL"] = new MultiSamplePlot(dataSets,"KinFit_FcncTopMass_NoBtag_CVSL",280,100,240,"m^{FCNC}_{top} [Gev/c^{2}]");

  ////////////////////////////////////////////////////////////////////
  ////////////////// 1D histograms  //////////////////////////////////
  ////////////////////////////////////////////////////////////////////


  cout << " - Declared histograms ..." <<  endl;
	
  ////////////////////////////////////////////////////////////////////
  ////////////////// Plots  //////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  string pathPNG = "TopFCNC"+postfix+channelpostfix+comments;
  pathPNG += "_MSPlots/"; 	
//  pathPNG = pathPNG +"/"; 	
  mkdir(pathPNG.c_str(),0777);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////// Analysis //////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  for(unsigned int iDataSet=0; iDataSet<inputFiles.size(); iDataSet++)
  {
    string dataSetName = dataSets[iDataSet]->Name();
    cout << "Processing DataSet: " << dataSetName << endl;

    TFile* inFile = TFile::Open(inputFiles[iDataSet].c_str());
    //TFile* inFile = new TFile(inputFiles[iDataSet].c_str(),"READ");
    
    TTree* inTree = (TTree*)   inFile->Get("Tree");
    TBranch* m_br = (TBranch*) inTree->GetBranch("TheTopFCNC_Evt");
    
    int nEvent = inTree->GetEntries();
    cout<< " - number of entries: "<<nEvent<<endl;
    
    TTree* inConfigTree = (TTree*) inFile->Get("configTree");
    TBranch* d_br = (TBranch*) inConfigTree->GetBranch("Dataset");
    TClonesArray* tc_dataset = new TClonesArray("Dataset",0);
    d_br->SetAddress(&tc_dataset);
    
    inConfigTree->GetEvent(0);
    Dataset* dataSet = (Dataset*) tc_dataset->At(0);
    
    TopFCNC_Evt* topFCNC_Evt = 0;
    m_br->SetAddress(&topFCNC_Evt);
    m_br->SetAutoDelete(kTRUE);
    
    TopFCNC_KinFit *topFCNC_KinFit = new TopFCNC_KinFit(dataSet, resFitLeptons, resFitBJets, resFitQJets, resFitLightJets);
    topFCNC_KinFit->SetMaxNbIter(60);
    topFCNC_KinFit->SetVerbosity(false);
    topFCNC_KinFit->SetFitVerbosity(false);

    double topMass = 172.5;
    //if(dataSet->Name().find("Data") == 0 || dataSet->Name().find("data") == 0 || dataSet->Name().find("DATA") == 0 )
    //  topMass = 173.3;
    // Top quark mass LHC average = 172.6 GeV/c²
    double kin_prob        = -1.;
    double kin_chi2        = -1.;
    double kin_chi2ByNdf   = -1.;
    double kin_hadWmass    = -1.;
    double kin_hadtopmass  = -1.;
    double kin_fcnctopmass = -1.;

    for(int iEvt=0; iEvt<nEvent; iEvt++)
    {
      inTree->GetEvent(iEvt);
//      cout << "event: " << iEvt << endl;
      if(iEvt%1000 == 0)
		    std::cout<<"Processing the "<<iEvt<<"th event, time = "<< ((double)clock()-start)/CLOCKS_PER_SEC << " ("<<100*(iEvt)/(nEvent)<<"%)"<<flush<<"\r";

      if(!topFCNC_Evt->isDiLeptonic()) continue;
//      cout<<"Nb of selected jets: "<<topFCNC_Evt->selectedJets().size()<<endl;
      
      // float WMass = 80.4, float Zmass = 91.2, float topMass = 172.5/173.3 MC vs data);
      topFCNC_KinFit->FitEvent(topFCNC_Evt,80.4,91.2,topMass);

      topFCNC_Evt->ReconstructEvt();

      EventWeight = Luminosity*topFCNC_Evt->eventWeight();

      kin_prob        = topFCNC_KinFit->GetProb();
      kin_chi2        = topFCNC_KinFit->GetChi2();
      if(topFCNC_KinFit->GetNdof()!=0)
        kin_chi2ByNdf   = kin_chi2/topFCNC_KinFit->GetNdof();
/*
      cout<<"HadQ1 properties :"<<topFCNC_Evt->quark1FromW()<<endl;
      cout<<"HadQ2 properties :"<<topFCNC_Evt->quark2FromW()<<endl;
      cout<<"Had W properties :"<<topFCNC_Evt->W()<<endl;
*/
      kin_hadWmass    = topFCNC_Evt->W().M();
      kin_hadtopmass  = topFCNC_Evt->smDecayTop().M();
      kin_fcnctopmass = topFCNC_Evt->fcncDecayTop().M();
      
      MSPlot["KinFit_Prob"]->Fill(kin_prob, dataSet, true, EventWeight);
      MSPlot["KinFit_Chi2"]->Fill(kin_chi2, dataSet, true, EventWeight);
      if(topFCNC_KinFit->GetNdof()!=0)
        MSPlot["KinFit_ReducedChi2"]->Fill(kin_chi2ByNdf, dataSet, true, EventWeight);

      MSPlot["KinFit_HadWMass"]   ->Fill(kin_hadWmass, dataSet, true, EventWeight);
      MSPlot["KinFit_HadTopMass"] ->Fill(kin_hadtopmass, dataSet, true, EventWeight);
      MSPlot["KinFit_FcncTopMass"]->Fill(kin_fcnctopmass, dataSet, true, EventWeight);

      vector<TRootJet> selectedJets = topFCNC_Evt->selectedJets();
      sort(selectedJets.begin(),selectedJets.end(),HighestCVSBtag());
      double bdisc = selectedJets[0].btag_combinedSecondaryVertexBJetTags();

      if(bdisc>0.679){ // DO NOT FORGET THE B-TAGGING SF WHEN MC !!!!!!
        MSPlot["KinFit_Prob_AtLeast1Btag_CSVM"]->Fill(kin_prob, dataSet, true, EventWeight);
        MSPlot["KinFit_Chi2_AtLeast1Btag_CSVM"]->Fill(kin_chi2, dataSet, true, EventWeight);
        if(topFCNC_KinFit->GetNdof()!=0)
          MSPlot["KinFit_ReducedChi2_AtLeast1Btag_CSVM"]->Fill(kin_chi2ByNdf, dataSet, true, EventWeight);

        MSPlot["KinFit_HadWMass_AtLeast1Btag_CSVM"]   ->Fill(kin_hadWmass, dataSet, true, EventWeight);
        MSPlot["KinFit_HadTopMass_AtLeast1Btag_CSVM"] ->Fill(kin_hadtopmass, dataSet, true, EventWeight);
        MSPlot["KinFit_FcncTopMass_AtLeast1Btag_CSVM"]->Fill(kin_fcnctopmass, dataSet, true, EventWeight);
      }
      else if(bdisc<0.244){ // DO NOT FORGET THE B-TAGGING SF WHEN MC !!!!!!!!!!
        MSPlot["KinFit_Prob_NoBtag_CVSL"]->Fill(kin_prob, dataSet, true, EventWeight);
        MSPlot["KinFit_Chi2_NoBtag_CVSL"]->Fill(kin_chi2, dataSet, true, EventWeight);
        if(topFCNC_KinFit->GetNdof()!=0)
          MSPlot["KinFit_ReducedChi2_NoBtag_CVSL"]->Fill(kin_chi2ByNdf, dataSet, true, EventWeight);

        MSPlot["KinFit_HadWMass_NoBtag_CVSL"]   ->Fill(kin_hadWmass, dataSet, true, EventWeight);
        MSPlot["KinFit_HadTopMass_NoBtag_CVSL"] ->Fill(kin_hadtopmass, dataSet, true, EventWeight);
        MSPlot["KinFit_FcncTopMass_NoBtag_CVSL"]->Fill(kin_fcnctopmass, dataSet, true, EventWeight);
      }

    } // loop on events
    
    delete topFCNC_KinFit;
  } // loop on datasets

  fout->cd();
  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
  {
	  MultiSamplePlot *temp = it->second;
	  //temp->addText("CMS preliminary");
	  string name = it->first;
	  cout<<"Booking MS :"<<name<<endl;
	  temp->Draw(false, name, true, true, true, true, true,1,true); // merge TT/QCD/W/Z/ST/
	  //Draw(bool addRandomPseudoData = false, string label = string("CMSPlot"), bool mergeTT = false, bool mergeQCD = false, bool mergeW = false, bool mergeZ = false, bool mergeST = false, int scaleNPSignal = 1, bool addRatio = false, bool mergeVV = false, bool mergeTTV = false);
	  temp->Write(fout, name, true, pathPNG, "pdf");
  }
  
  //delete  
  delete fout;

  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;

  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "********************************************" << endl;

  return 0;
}

