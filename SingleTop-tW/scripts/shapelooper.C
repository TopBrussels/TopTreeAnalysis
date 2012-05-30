#include "shapelooper.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "inputs.h"
#include "Riostream.h"
#include <vector>
#include <string>
#include "TFile.h"
#include "TChain.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"

void shapelooper::Loop(){
  // running default loop:
  myLoop(0,0,0,0,0);
}
void shapelooper::myLoop(int nsel, int mode, int syst, bool up, bool silent)
{
  
  char plotName[300];
  sprintf(plotName,"test");
  
  if (nsel == 0)                	{sprintf(plotName,"tt");}
  else if (nsel == 1)   		{sprintf(plotName,"twdr");}
  else if (nsel == -1)   		{sprintf(plotName,"twds");}
  else if (nsel == 2)   		{sprintf(plotName,"zjets");}
  else if (nsel == 3)   		{sprintf(plotName,"di");}
  else if (nsel == 4)			{sprintf(plotName, "tschan");}
  else if (nsel == 5)   		{sprintf(plotName,"wjets");}
  else if (nsel == 6)   		{sprintf(plotName,"qcd_mu");}
  else if (nsel == 7)                	{sprintf(plotName,"others");}
  else if (nsel == 8)                	{sprintf(plotName,"tt_scaleup");}
  else if (nsel == 9)                	{sprintf(plotName,"tt_scaledown");}
  else if (nsel == 10)                	{sprintf(plotName,"tw_sdo");}
  else if (nsel == 11)                	{sprintf(plotName,"tw_sup");}
  else if (nsel == 12)                	{sprintf(plotName,"tt_matchingup");}
  else if (nsel == 13)                	{sprintf(plotName,"tt_matchingdown");}
  
  else if (nsel == 555)                	{sprintf(plotName,"mc");}
  
  else if (nsel == 666)                	{sprintf(plotName,"data");}
  else if (nsel == 6661)                {sprintf(plotName,"data1");}
  else if (nsel == 6662)                {sprintf(plotName,"data2");}
  
  bool nosf = false;
  
  char newRootFile[300];
  double lumi = luminosity; 
  if (mode == 0 )        lumi = 4904.338;
  else if ( mode == 1)   lumi = 4919.924;
  else if ( mode == 2)   lumi = 4895.249;
  else if ( mode == 3)   lumi = 4.9;
  sprintf(newRootFile,"results/Histos_cutbased_CR.root");
  
  TFile f_var(newRootFile, "UPDATE");
  
  if(!silent){
    std::cout << "[Info:] results root file " << newRootFile << std::endl;
  }

  bool SFplus = false;
  bool SFminus = false;
  
  char modeName[300];
  if (mode == 0) sprintf(modeName,"emu");
  else if (mode == 1) sprintf(modeName,"mumu");
  else sprintf(modeName,"ee");
  
  char longPlotName[300];
  if (syst == 0) sprintf(longPlotName,"%s", plotName);
  else if (syst == 1 && up) sprintf(longPlotName,"%s__PU__plus", plotName);
  else if (syst == 1 && !up) sprintf(longPlotName,"%s__PU__minus", plotName);
  else if (syst == 2 && up) sprintf(longPlotName,"%s__JES__plus", plotName);
  else if (syst == 2 && !up) sprintf(longPlotName,"%s__JES__minus", plotName);
  else if (syst == 3 && up) sprintf(longPlotName,"%s__JER__plus", plotName);
  else if (syst == 3 && !up) sprintf(longPlotName,"%s__JER__minus", plotName);
  else if (syst == 4 && up) sprintf(longPlotName,"%s__MET__plus", plotName);
  else if (syst == 4 && !up) sprintf(longPlotName,"%s__MET__minus", plotName);
  else if (syst == 5 && up) {
    sprintf(longPlotName,"%s__btag__plus", plotName);
    SFplus = true;
  }  
  else if (syst == 5 && !up) {
    sprintf(longPlotName,"%s__btag__minus", plotName);
    SFminus = true;
  }  

  
  //////////
  char title[300];
  
  /// Classic plotmaker plots
  sprintf(title,"%s1j1t__%s",modeName, longPlotName);
  TH1F* histo_1j1t = new TH1F( title, " ", 1, 1, 2);
  histo_1j1t->Sumw2();
  
  sprintf(title,"%s2j1t__%s",modeName, longPlotName);
  TH1F* histo_2j1t = new TH1F( title, " ", 1, 1, 2);
  histo_2j1t->Sumw2();
  
  sprintf(title,"%s2j2t__%s",modeName, longPlotName);
  TH1F* histo_2j2t = new TH1F( title, " ", 1, 1, 2);
  histo_2j2t->Sumw2();
  
  sprintf(title,"%sall__%s",modeName, longPlotName);
  TH1F* histo_all = new TH1F( title, " ", 3, 1, 4);
  histo_all->Sumw2();
  
  // tt control regions
  sprintf(title,"R_%s_%s",modeName, longPlotName);
  TH1F* histo_R = new TH1F( title, " ", 40,  0, 40 );
  histo_R->Sumw2();
  
  
  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    
    if (lumi != lum && nsel != 666 && mode !=3){
      if (jentry == 0)std::cout << "[Warning:] This tree was made with a different luminosity (" << lum << ") than " << lumi << std::endl;
      //xlWeight*=(lumi/lum);
    }
    
    if(ptLepton->size() != 2){
      std::cout << "[Warning:] Something is wrong, your Tree is not correctly filled" << std::endl;
      break;
    } else {
      
      TLorentzVector lepton0(pxLepton->at(0),pyLepton->at(0), pzLepton->at(0), eLepton->at(0));
      TLorentzVector lepton1(pxLepton->at(1),pyLepton->at(1), pzLepton->at(1), eLepton->at(1)); 
      TLorentzVector pair = lepton0+lepton1;
      
      double phipairmet_t = 0;
      double pi_m = 3.1416/2;
      phipairmet_t = pi_m;
      
      TVector3 vmet(metPx, metPy, 0);
      
      double promet = metPt*sin(phipairmet_t);
      TVector3 m0(pxLepton->at(0),pyLepton->at(0), pzLepton->at(0));
      TVector3 m1(pxLepton->at(1),pyLepton->at(1), pzLepton->at(1)); 
      if (fabs(m0.DeltaPhi(vmet)) < phipairmet_t) phipairmet_t = fabs(m0.DeltaPhi(vmet));
      if (fabs(m1.DeltaPhi(vmet)) < phipairmet_t) phipairmet_t = fabs(m1.DeltaPhi(vmet));
      
      if (phipairmet_t == pi_m) promet = metPt;
      
      if (pair.M() > 20){
	
	double SFval, SFerror;
	
	if ( nsel == 666 || nosf){
	  SFval = 1;
	  SFerror = 0;
	} else if (nsel == 0){
	  SFval = 0.956;
	  SFerror = 0.030;
	} else {
	  SFval = 0.96;
	  SFerror = 0.04;
	}
	
	int nJetsBT = 0;
	int nTightJetsBT = 0;
	int nJets = 0;
	bool bTagged = false;
	int iJet = -5;
	int iSF;
	double tempSF = SFval;
	if (SFminus) 	tempSF = SFval - SFerror;
	if (SFplus) 	tempSF = SFval + SFerror;
	int SFvalue = int(tempSF*100);
	
	for (unsigned int i =0; i < ptJet->size(); i ++){ 
	  TLorentzVector tempJet(pxJet->at(i),pyJet->at(i), pzJet->at(i), eJet->at(i));
	  if (ptJet->at(i) > 30 && TMath::Min(fabs(lepton0.DeltaR(tempJet)), fabs(lepton1.DeltaR(tempJet))) > 0.3) {
	    nJets++;
	    iJet = i;
	    if (btSSVHEJet->at(i) > 1.74){
	      iSF = rand() % 100;
	      if (iSF < SFvalue ){
		bTagged = true;
		nJetsBT++;
		nTightJetsBT++;
	      } 
	    } 
	  } else if (btSSVHEJet->at(i) > 1.74){
	    iSF = rand() % 100;
	    if (iSF < SFvalue ) nJetsBT++;
	  }
	}
	
	bool invMass = false;
	if      (mode == 0) invMass = true;
	else if (mode == 1  && (pair.M() > invMax || pair.M() < invMin)) invMass = true;
	else if (mode == 2 && (pair.M() > invMax || pair.M() < invMin)) invMass = true;
	
	if (invMass) {
	  if (TMath::Min(metPt, tmetPt) >= metCut || mode == 0){
	    if(nJets !=0){
	  
	      TLorentzVector jet(pxJet->at(iJet),pyJet->at(iJet), pzJet->at(iJet), eJet->at(iJet));
	      
	      double ptSysPx = lepton0.Px() + lepton1.Px() + jet.Px() + metPx;
	      double ptSysPy = lepton0.Py() + lepton1.Py() + jet.Py() + metPy;
	      double ptSystem = sqrt(ptSysPx*ptSysPx + ptSysPy*ptSysPy);
	      double ht = lepton0.Pt() + lepton1.Pt() + jet.Pt() + metPt; 
	      
	      //if (ptSystem <= ptsysCut){
	      if (mode != 999){	
		//if (ht > htMin || mode !=0){
		if (mode  != 999){  
		 
		  //DY normalization
		  if (nsel == 2 && mode == 1){
		    if (metPt <= 5) xlWeight*=0.800093;
		    else if (metPt <= 10) xlWeight*=0.822401;
		    else if (metPt <= 15) xlWeight*=0.879812;
		    else if (metPt <= 20) xlWeight*=0.948094;
		    else if (metPt <= 25) xlWeight*=1.03923;
		    else if (metPt <= 30) xlWeight*=1.14394;
		    else if (metPt <= 35) xlWeight*=1.25926;
		    else if (metPt <= 40) xlWeight*=1.39992;
		    else if (metPt <= 45) xlWeight*=1.49694;
		    else if (metPt <= 50) xlWeight*=1.60845;
		    else if (metPt <= 60) xlWeight*=1.78741;
		    else  xlWeight*=1.4131;
		  } else if (nsel == 2 && mode == 2){
		    if (metPt <= 5) xlWeight*=0.791871;
		    else if (metPt <= 10) xlWeight*=0.823337;
		    else if (metPt <= 15) xlWeight*=0.893007;
		    else if (metPt <= 20) xlWeight*=0.965893;
		    else if (metPt <= 25) xlWeight*=1.07273;
		    else if (metPt <= 30) xlWeight*=1.17567;
		    else if (metPt <= 35) xlWeight*=1.28749;
		    else if (metPt <= 40) xlWeight*=1.43283;
		    else if (metPt <= 45) xlWeight*=1.51069;
		    else if (metPt <= 50) xlWeight*=1.80464;
		    else if (metPt <= 60) xlWeight*=1.86129;
		    else  xlWeight*= 1.68509;
		  }else if (nsel == 2){
		    if (metPt <= 5) xlWeight*=0.8;
		    else if (metPt <= 10) xlWeight*=0.88;
		    else if (metPt <= 15) xlWeight*=0.88;
		    else if (metPt <= 20) xlWeight*=0.95;
		    else if (metPt <= 25) xlWeight*=1.05;
		    else if (metPt <= 30) xlWeight*=1.15;
		    else if (metPt <= 35) xlWeight*=1.27;
		    else if (metPt <= 40) xlWeight*=1.41;
		    else if (metPt <= 45) xlWeight*=1.50;
		    else if (metPt <= 50) xlWeight*=1.70;
		    else if (metPt <= 60) xlWeight*=1.80;
		    else  xlWeight*= 1.5;
		  }
		  
		  //Histos by channel
		  if (nJets == 1 && nTightJetsBT == 1 && bTagged && nJetsBT == 1 && ptSystem <= ptsysCut && (ht > htMin || mode !=0))histo_1j1t->Fill(1.5, xlWeight);
		  if (nJets == 2 && nTightJetsBT == 1 && bTagged && nJetsBT == 1 )histo_2j1t->Fill(1.5, xlWeight);
		  if (nJets == 2 && nTightJetsBT == 2 && bTagged && nJetsBT == 2 )histo_2j2t->Fill(1.5, xlWeight);
		  
		  //Global histo
		  if (nJets == 1 && nTightJetsBT == 1 && bTagged && nJetsBT == 1 && ptSystem <= ptsysCut && (ht > htMin || mode !=0))histo_all->Fill(1.5, xlWeight);
		  if (nJets == 2 && nTightJetsBT == 1 && bTagged && nJetsBT == 1 )histo_all->Fill(2.5, xlWeight);
		  if (nJets == 2 && nTightJetsBT == 2 && bTagged && nJetsBT == 2 )histo_all->Fill(3.5, xlWeight);
		  
		  //All possible regions
		  if (nJets == 1 && nTightJetsBT == 1 && bTagged && nJetsBT == 1 && ptSystem <= ptsysCut && (ht > htMin || mode !=0))histo_R->Fill(1, xlWeight); //signal
		  if (nJets == 1 && nTightJetsBT == 2)  histo_R->Fill(2, xlWeight);
		  if (nJets == 1 && nTightJetsBT > 0)  histo_R->Fill(3, xlWeight);
		  if (nJets == 1 && nTightJetsBT > 1)  histo_R->Fill(4, xlWeight);
		  if (nJets == 2 && nTightJetsBT == 0)  histo_R->Fill(5, xlWeight);
		  if (nJets == 2 && nTightJetsBT == 1)  histo_R->Fill(6, xlWeight); //CR1 no ht no ptsys
		  if (nJets == 2 && nTightJetsBT == 2)  histo_R->Fill(7, xlWeight); //CR2 no ht no ptsys
		  if (nJets == 2 && nTightJetsBT > 0)  histo_R->Fill(8, xlWeight);
		  if (nJets == 2 && nTightJetsBT > 1)  histo_R->Fill(9, xlWeight);
		  if (nJets > 1 && nTightJetsBT == 0)  histo_R->Fill(10, xlWeight);
		  if (nJets > 1 && nTightJetsBT == 1)  histo_R->Fill(11, xlWeight);
		  if (nJets > 1 && nTightJetsBT == 2)  histo_R->Fill(12, xlWeight);
		  if (nJets > 1 && nTightJetsBT !=0 )  histo_R->Fill(13, xlWeight);
		  if (nJets > 1 && nTightJetsBT > 1 )  histo_R->Fill(14, xlWeight);
		  if (nJets == 3 && nTightJetsBT ==3 )  histo_R->Fill(15, xlWeight);
		  if (nJets == 1 && nTightJetsBT ==1 && bTagged && nJetsBT == 1)  histo_R->Fill(16, xlWeight);
		  if (nJets == 2 && nTightJetsBT == 1 && ptSystem <= ptsysCut && (ht > htMin || mode !=0))  histo_R->Fill(17, xlWeight); //CR 1 regular
		  if (nJets == 2 && nTightJetsBT == 2 && ptSystem <= ptsysCut && (ht > htMin || mode !=0))  histo_R->Fill(18, xlWeight); //CR 2 regular
		  if (nJets == 2 && nTightJetsBT == 1 && nJetsBT == 1 && ptSystem <= ptsysCut && (ht > htMin || mode !=0))  histo_R->Fill(19, xlWeight);
		  if (nJets == 2 && nTightJetsBT == 2 && nJetsBT == 2 && ptSystem <= ptsysCut && (ht > htMin || mode !=0))  histo_R->Fill(20, xlWeight);
		  if (nJets == 2 && nTightJetsBT == 1 && nJetsBT == 1)  histo_R->Fill(21, xlWeight); //CR1 no ht no ptsys tighter
		  if (nJets == 2 && nTightJetsBT == 2 && nJetsBT == 2)  histo_R->Fill(22, xlWeight); //CR2 no ht no ptsys tighter
		  
		  
		}
	      }
	    }
	  }
	}
      } //mll pre
    } // 2 leptons
  }// event loop.
  
  
  
  if (!silent){ 
    cout << "------------------------------------------" << endl;
    cout << "[1j1t:] " << histo_R->GetBinContent(2) << "\t\t[2j1t:] " << histo_R->GetBinContent(18) << "\t[2j2t:] " << histo_R->GetBinContent(19) << endl;
    cout << "Removing Ht and Pt sys\t" << "\t[2j1t:] " << histo_R->GetBinContent(7) << "\t[2j2t:] " << histo_R->GetBinContent(8) << endl;
    cout << "------------------------------------------" << endl; 
  }
  f_var.Write();
  f_var.Close();
}
