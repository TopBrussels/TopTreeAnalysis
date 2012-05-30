#include "metlooper.h"
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

void metlooper::Loop(){
  // running default loop:
  myLoop(0,0,0);
}
void metlooper::myLoop(int nsel, int mode, bool silent)
{

  char plotName[300];
  sprintf(plotName,"test");
  
  if (nsel == 0)                	{sprintf(plotName,"tt");}
  else if (nsel == 1)   		{sprintf(plotName,"twdr");}
  else if (nsel == -1)   		{sprintf(plotName,"twds");}
  else if (nsel == 2)   		{sprintf(plotName,"zjets");}
  else if (nsel == 3)   		{sprintf(plotName,"di");}
  else if (nsel == 4)			{sprintf(plotName, "st");}
  else if (nsel == 5)   		{sprintf(plotName,"wjets");}
  else if (nsel == 6)   		{sprintf(plotName,"qcd_mu");}
  else if (nsel == 7)                	{sprintf(plotName,"others");}
  
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
  sprintf(newRootFile,"results/met_%dpb_%d.root", (int)lumi, mode);
 
  TFile f_var(newRootFile, "UPDATE");
  
  if(!silent){
    std::cout << "[Info:] results root file " << newRootFile << std::endl;
  }
  
  
  //////////
  char title[300];
  
  /// Classic plotmaker plots
  Double_t xbins[13] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 200} ;
  sprintf(title,"met_%s",plotName);
  TH1F* histo_met = new TH1F( title, " ", 12, xbins );
  histo_met->Sumw2();
  
  Double_t xbinsd[7] = {0, 10, 20, 30, 40, 50, 200} ;
  sprintf(title,"metd_%s",plotName);
  TH1F* histo_metd = new TH1F( title, " ", 6, xbinsd );
  histo_metd->Sumw2();
  
  

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
	

	

	

	/*
	    if (nJets == 1){
	 
	      TLorentzVector jet(pxJet->at(iJet),pyJet->at(iJet), pzJet->at(iJet), eJet->at(iJet));
		
	      double ptSysPx = lepton0.Px() + lepton1.Px() + jet.Px() + metPx;
	      double ptSysPy = lepton0.Py() + lepton1.Py() + jet.Py() + metPy;
	      double ptSystem = sqrt(ptSysPx*ptSysPx + ptSysPy*ptSysPy);
	      double ht = lepton0.Pt() + lepton1.Pt() + jet.Pt() + metPt; 
	   
	      if (nJets == 1 && nTightJetsBT == 1 && bTagged && nJetsBT == 1){
		
		*/
		//if (ptSystem <= ptsysCut){
		
		  //if (ht > htMin || mode !=0){
		    
		    if (pair.M() <= invMax && pair.M() >= invMin){
		    
		     /*
    if (nsel == 2 && mode == 1){
      if (metPt <= 5) xlWeight*=0.917087;
      else if (metPt <= 10) xlWeight*=0.94228;
      else if (metPt <= 15) xlWeight*=1.02343;
      else if (metPt <= 20) xlWeight*=1.11051;
      else if (metPt <= 25) xlWeight*=1.20101;
      else if (metPt <= 30) xlWeight*=1.26738;
      else if (metPt <= 35) xlWeight*=1.46525;
      else if (metPt <= 40) xlWeight*=1.47895;
      else if (metPt <= 45) xlWeight*=1.46037;
      else if (metPt <= 50) xlWeight*=1.28589;
      else if (metPt <= 60) xlWeight*=1.15589;
      else  xlWeight*=0.689844;
    } else if (nsel == 2 && mode == 2){
      if (metPt <= 5) xlWeight*=0.917655;
      else if (metPt <= 10) xlWeight*=1.04012;
      else if (metPt <= 15) xlWeight*=0.986518;
      else if (metPt <= 20) xlWeight*=1.11664;
      else if (metPt <= 25) xlWeight*=1.16103;
      else if (metPt <= 30) xlWeight*=1.22194;
      else if (metPt <= 35) xlWeight*=1.43669;
      else if (metPt <= 40) xlWeight*=1.79543;
      else if (metPt <= 45) xlWeight*=1.23258;
      else if (metPt <= 50) xlWeight*=1.69726;
      else if (metPt <= 60) xlWeight*=1.89196;
      else  xlWeight*= 0.743117;
    }

		    */
		    //dilepton
		    /*
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
    }*/
		    
		     histo_met->Fill(metPt, xlWeight);
		     histo_metd->Fill(metPt, xlWeight);
		    //}
		  //}
		}
	     // }
	   // }
	  
	  

      } //mll pre
    } // 2 leptons
  }// event loop.
  
  
  
  if (!silent){ 
    cout << "------------------------------------------" << endl;
    cout << "[You are done with:] " << plotName <<  endl;
    cout << "------------------------------------------" << endl; 
  }
  f_var.Write();
  f_var.Close();
}
