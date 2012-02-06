// rebeca@cern.ch
#include "TLorentzVector.h"
#include "TVector3.h"
#include "inputs.h"


void chain(int nsel = 0, int mode = 0, bool silent = false){  
  
 
  bool SFplus = false;
  bool SFminus = false;
  // samples used
  double x_sec = 0.;
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
  
  // if (mode != 0 &&  mode !=1 && mode !=2) mode = 0;
  if (!silent){
    if      (mode == 0) 	cout << " Electron-Muon Mixed channel " << endl;
    else if (mode == 1) 	cout << " Di-Muon channel " << endl;
    else if (mode == 2) 	cout << " Di-Electron channel " << endl;
    cout << "*************************************************" << endl;
  }
  
  bool nosf = true;
  char myRootFile[300];
  sprintf(myRootFile,"outputs/out_%d_%s.root", mode, plotName);
  TFile *input = TFile::Open(myRootFile);
  
  // tree variables
  /////
  double xlWeight; 
  double lum;
  
  double metPt;
  double metPx;
  double metPy;
  
  int npu;
  int nvertex;
      
  std::vector<double> *ptLepton;
  std::vector<double> *pxLepton;
  std::vector<double> *pyLepton;
  std::vector<double> *pzLepton;
  std::vector<double> *eLepton;
  std::vector<double> *qLepton;
  
  std::vector<double> *ptJet;
  std::vector<double> *pxJet;
  std::vector<double> *pyJet;
  std::vector<double> *pzJet;
  std::vector<double> *eJet;
  std::vector<double> *qJet;
  std::vector<double> *btSSVHEJet;
  std::vector<double> *btSSVHPJet;
  std::vector<double> *btTCHEJet;
  std::vector<double> *btTCHPJet;
  
  TTree*Tree = (TTree*) gROOT->FindObject("myTree");
  
  Tree->SetBranchAddress("xlWeight", &xlWeight);
  Tree->SetBranchAddress("lum", &lum);
  
  Tree->SetBranchAddress("npu", &npu);
  Tree->SetBranchAddress("nvertex", &nvertex);
  
  Tree->SetBranchAddress("metPt", &metPt);
  Tree->SetBranchAddress("metPx", &metPx);
  Tree->SetBranchAddress("metPy", &metPy);
  
  Tree->SetBranchAddress("ptLepton",&ptLepton);
  Tree->SetBranchAddress("pxLepton",&pxLepton);
  Tree->SetBranchAddress("pyLepton",&pyLepton);
  Tree->SetBranchAddress("pzLepton",&pzLepton);
  Tree->SetBranchAddress("eLepton",&eLepton);
  Tree->SetBranchAddress("qLepton",&qLepton);
  
  Tree->SetBranchAddress("ptJet",&ptJet);
  Tree->SetBranchAddress("pxJet",&pxJet);
  Tree->SetBranchAddress("pyJet",&pyJet);
  Tree->SetBranchAddress("pzJet",&pzJet);
  Tree->SetBranchAddress("eJet",&eJet);
  Tree->SetBranchAddress("qJet",&qJet);
  Tree->SetBranchAddress("btSSVHEJet",&btSSVHEJet);
  Tree->SetBranchAddress("btSSVHPJet",&btSSVHPJet);
  Tree->SetBranchAddress("btTCHEJet",&btTCHEJet);
  Tree->SetBranchAddress("btTCHPJet",&btTCHPJet);
  
  int nEvents = Tree->GetEntries();
  if(!silent){
    cout << endl;
    cout << "******************************************" << endl;
    cout << "------------------------------------------" << endl;
    cout << "Starting the analysis: " << plotName <<  endl;
    cout << "------------------------------------------" << endl;
    cout << "Number of Raw events: " <<  nEvents << endl;
    cout << "------------------------------------------" << endl;
    cout << "******************************************" << endl;
  }
  
  char newRootFile[300];
  double lumi = luminosity;
  if (mode == 0 )        lumi = 4626.297;
  else if ( mode == 1)   lumi = 4534.871;
  else if ( mode == 2)   lumi = 4593.348;
  else if ( mode == 3)   lumi = 4.5;
  sprintf(newRootFile,"results/an_%dpb_%d.root", lumi, mode);
  
  TFile f_var(newRootFile, "UPDATE");
  
  //////////
  char title[300];
  sprintf(title,"cuts_%s",plotName);
  TH1F* histo = new TH1F( title, " ", 10,  0, 10 );
  histo->Sumw2();
  
  sprintf(title,"met_%s",plotName);
  TH1F* histo_met = new TH1F( title, " ", 100,  0, 200 );
  histo_met->Sumw2();
  
  sprintf(title,"met_high_%s",plotName);
  TH1F* histo_met_high = new TH1F( title, " ", 100,  0, 200 );
  histo_met_high->Sumw2();
  
  sprintf(title,"met_low_%s",plotName);
  TH1F* histo_met_low = new TH1F( title, " ", 100,  0, 200 );
  histo_met_low->Sumw2();
  
  sprintf(title,"promet_%s",plotName);
  TH1F* histo_promet = new TH1F( title, " ", 100,  0, 200 );
  histo_promet->Sumw2();
  
  sprintf(title,"met_cut_%s",plotName);
  TH1F* histo_met_cut = new TH1F( title, " ", 100,  0, 200 );
  histo_met_cut->Sumw2();
  
  sprintf(title,"met_bt_%s",plotName);
  TH1F* histo_met_bt = new TH1F( title, " ", 100,  0, 200 );
  histo_met_bt->Sumw2();
  
  sprintf(title,"mll_%s",plotName);
  TH1F* histo_mll = new TH1F( title, " ", 100,  0, 200 );
  histo_mll->Sumw2();
  
  sprintf(title,"mll_after_%s",plotName);
  TH1F* histo_mll_after = new TH1F( title, " ", 100,  0, 200 );
  histo_mll_after->Sumw2();
  
  sprintf(title,"njets_%s",plotName);
  TH1F* histo_njets = new TH1F( title, " ", 10,  0, 10 );
  histo_njets->Sumw2();
  
  sprintf(title,"njets_cut_%s",plotName);
  TH1F* histo_njets_cut = new TH1F( title, " ", 10,  -0.5, 9.5 );
  histo_njets_cut->Sumw2();
  
  sprintf(title,"njetsbt_%s",plotName);
  TH1F* histo_njetsbt = new TH1F( title, " ", 10,  -0.5, 9.5 );
  histo_njetsbt->Sumw2();
  
  sprintf(title,"njetsbt_cut_%s",plotName);
  TH1F* histo_njetsbt_cut = new TH1F( title, " ", 10,   -0.5, 9.5 );
  histo_njetsbt_cut->Sumw2();
  
  sprintf(title,"njetsbt_high_%s",plotName);
  TH1F* histo_njetsbt_high = new TH1F( title, " ", 10,   -0.5, 9.5 );
  histo_njetsbt_high->Sumw2();

  sprintf(title,"njetsbt_low_%s",plotName);
  TH1F* histo_njetsbt_low = new TH1F( title, " ", 10,   -0.5, 9.5 );
  histo_njetsbt_low->Sumw2();
  
  sprintf(title,"njets_high_%s",plotName);
  TH1F* histo_njets_high = new TH1F( title, " ", 10,   -0.5, 9.5 );
  histo_njets_high->Sumw2();

  sprintf(title,"njets_low_%s",plotName);
  TH1F* histo_njets_low = new TH1F( title, " ", 10,   -0.5, 9.5 );
  histo_njets_low->Sumw2();
  
  sprintf(title,"ptsys_%s",plotName);
  TH1F* histo_ptsys = new TH1F( title, " ", 100,  0, 200 );
  histo_ptsys->Sumw2();
  
  sprintf(title,"ptsys_high_%s",plotName);
  TH1F* histo_ptsys_high = new TH1F( title, " ", 100,  0, 200 );
  histo_ptsys_high->Sumw2();
  
  sprintf(title,"ptsys_low_%s",plotName);
  TH1F* histo_ptsys_low = new TH1F( title, " ", 100,  0, 200 );
  histo_ptsys_low->Sumw2();
  
  sprintf(title,"ht_%s",plotName);
  TH1F* histo_ht = new TH1F( title, " ", 300,  0, 600 );
  histo_ht->Sumw2();
  
  sprintf(title,"ht_high_%s",plotName);
  TH1F* histo_ht_high = new TH1F( title, " ", 300,  0, 600 );
  histo_ht_high->Sumw2();
  
  sprintf(title,"ht_low_%s",plotName);
  TH1F* histo_ht_low = new TH1F( title, " ", 300,  0, 600 );
  histo_ht_low->Sumw2();
  
  sprintf(title,"ht_cut_%s",plotName);
  TH1F* histo_ht_cut = new TH1F( title, " ", 300,  0, 600 );
  histo_ht_cut->Sumw2();
  
  sprintf(title,"pt_max_%s",plotName);
  TH1F* histo_pt_max = new TH1F( title, " ", 100,  0, 200 );
  histo_pt_max->Sumw2();
  
  sprintf(title,"pt_min_%s",plotName);
  TH1F* histo_pt_min = new TH1F( title, " ", 100,  0, 200 );
  histo_pt_min->Sumw2();
  
  sprintf(title,"pt_leading_%s",plotName);
  TH1F* histo_pt_leading = new TH1F( title, " ", 100,  0, 200 );
  histo_pt_leading->Sumw2();
  
  sprintf(title,"btagHE_%s",plotName);
  TH2F* histo_btagHE = new TH2F( title, " ", 300,  -200, 100, 100, -2, 7);
  histo_btagHE->Sumw2();
  
  sprintf(title,"btagHP_%s",plotName);
  TH2F* histo_btagHP = new TH2F( title, " ", 300,  -200, 100, 100, -2, 7);
  histo_btagHP->Sumw2();
  
  sprintf(title,"etalepton_%s",plotName);
  TH1F* histo_etalepton = new TH1F( title, " ", 101,  -3, 3);
  histo_etalepton->Sumw2();
  
  sprintf(title,"ht_bf_%s",plotName);
  TH1F* histo_ht_bf = new TH1F( title, " ", 300,  0, 600 );
  histo_ht_bf->Sumw2();
  
  sprintf(title,"ptsys_bf_%s",plotName);
  TH1F* histo_ptsys_bf = new TH1F( title, " ", 100,  0, 200 );
  histo_ptsys_bf->Sumw2();
  
  
  
  sprintf(title,"nvertex_%s",plotName);
  TH1F* histo_nvertex = new TH1F( title, " ", 30,   -0.5, 29.5 );
  histo_nvertex->Sumw2();
  
  sprintf(title,"npu_%s",plotName);
  TH1F* histo_npu = new TH1F( title, " ", 30,   -0.5, 29.5 );
  histo_npu->Sumw2();
  
  sprintf(title,"R_%s",plotName);
  TH1F* histo_R = new TH1F( title, " ", 40,  0, 40 );
  histo_R->Sumw2();
  
  
  //////////
  for(int event = 0; event<nEvents; event++){
    
    Tree->GetEntry(event);
    
    if (lumi != lum && nsel != 666 && mode !=3){
      if (event == 0) cout << "Warning: This tree was made with a different luminosity (" << lum << ") than " << lumi << endl;
      xlWeight*=(lumi/lum);
    }
    
    if(ptLepton->size() != 2) cout << "Something is wrong, your Tree is not correctly filled" << endl;
    else {
      
      
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
	histo->Fill(1, xlWeight);
	
	// double SFval = 0.95;  //Summer11 version
	double SFval, SFerror;
	if ( nsel == 666 || !nosf){
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

	for (int i =0; i < ptJet->size(); i ++){ 
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
	//
      
	
	
	histo_pt_max->Fill(TMath::Max(lepton0.Pt(), lepton1.Pt()), xlWeight);
	histo_pt_min->Fill(TMath::Min(lepton0.Pt(), lepton1.Pt()), xlWeight);
	histo_njets->Fill(nJets,  xlWeight);
	histo_njetsbt->Fill(nJetsBT,  xlWeight);
	histo_mll->Fill(pair.M(),  xlWeight);
	histo_met->Fill(metPt,  xlWeight);
	histo_promet->Fill(promet, xlWeight);
	
	if (nvertex > 5){
	  histo_met_high->Fill(metPt,  xlWeight);
	  histo_njets_high->Fill(nJets,  xlWeight);
	  histo_njetsbt_high->Fill(nJetsBT,  xlWeight);
	} else {
	  histo_met_low->Fill(metPt,  xlWeight);
	  histo_njets_low->Fill(nJets,  xlWeight);
	  histo_njetsbt_low->Fill(nJetsBT,  xlWeight);  
	}
	
	if (nJets) histo_pt_leading->Fill(ptJet->at(0), xlWeight);
	if (nJets == 1){
	  histo_etalepton->Fill(lepton0.Eta(), xlWeight);
	  TLorentzVector jet(pxJet->at(iJet),pyJet->at(iJet), pzJet->at(iJet), eJet->at(iJet));
	   
	  double ptSysPx1 = lepton0.Px() + lepton1.Px() + jet.Px() + metPx;
	  double ptSysPy1 = lepton0.Py() + lepton1.Py() + jet.Py() + metPy;
	  double ptSystem1 = sqrt(ptSysPx1*ptSysPx1 + ptSysPy1*ptSysPy1);
	  double ht1 = lepton0.Pt() + lepton1.Pt() + jet.Pt() + metPt; 
	  histo_ptsys_bf->Fill(ptSystem1, xlWeight);
	  histo_ht_bf->Fill(ht1, xlWeight);
	 }
	bool invMass = false;
	if      (mode == 0) invMass = true;
	else if (mode == 1  && (pair.M() > invMax || pair.M() < invMin)) invMass = true;
	else if (mode == 2 && (pair.M() > invMax || pair.M() < invMin)) invMass = true;
	
	if (invMass){
	  histo->Fill(2, xlWeight);
	  histo_mll_after->Fill(pair.M(),  xlWeight);
	  histo_met_cut->Fill(metPt,  xlWeight);
	  
	  if (metPt >= metCut || mode ==0){
	    //if (promet >= metCut || mode ==0){
	    histo->Fill(3, xlWeight);
	    histo_njets_cut->Fill(nJets, xlWeight);
	    if (nJets == 1){
	       histo->Fill(4, xlWeight);
	       histo_njetsbt_cut->Fill(nJetsBT, xlWeight);
	       
	       if (nJets == 1 && nTightJetsBT == 1 && bTagged && nJetsBT == 1){
		 histo->Fill(5, xlWeight);
		 
		 double ptSysPx = lepton0.Px() + lepton1.Px() + jet.Px() + metPx;
		 double ptSysPy = lepton0.Py() + lepton1.Py() + jet.Py() + metPy;
		 double ptSystem = sqrt(ptSysPx*ptSysPx + ptSysPy*ptSysPy);
		 double ht = lepton0.Pt() + lepton1.Pt() + jet.Pt() + metPt; 
		 
		 histo_ptsys->Fill(ptSystem, xlWeight);
		 histo_ht->Fill(ht, xlWeight);
		 histo_btagHE->Fill(btTCHEJet->at(iJet), btSSVHEJet->at(iJet), xlWeight);
		 histo_btagHP->Fill(btTCHPJet->at(iJet), btSSVHPJet->at(iJet), xlWeight);
		 histo_met_bt->Fill(metPt, xlWeight);
		 
		 histo_nvertex->Fill(nvertex, xlWeight);
		 histo_npu->Fill(npu, xlWeight);
		 
		 if (nvertex > 5) {
		   histo_ptsys_high->Fill(ptSystem, xlWeight);
		   histo_ht_high->Fill(ht, xlWeight);
		 } else {
		   histo_ptsys_low->Fill(ptSystem, xlWeight);
		   histo_ht_low->Fill(ht, xlWeight);
		   
		 }
		 
		 if (ptSystem <= ptsysCut){
		   histo->Fill(6, xlWeight);
		   histo_ht_cut->Fill(ht, xlWeight);
		   if (ht > htMin || mode !=0){
		     histo->Fill(7, xlWeight);
		     
		     
		   }
		 }
	       }
	     }
	   }
	   
	   
	   //tt control region from here
	   if (metPt >= metCut || mode ==0){
	     
	     if (nJets != 0){
	     
	       TLorentzVector jet(pxJet->at(iJet),pyJet->at(iJet), pzJet->at(iJet), eJet->at(iJet));
	       
	       double ptSysPx = lepton0.Px() + lepton1.Px() + jet.Px() + metPx;
	       double ptSysPy = lepton0.Py() + lepton1.Py() + jet.Py() + metPy;
	       double ptSystem = sqrt(ptSysPx*ptSysPx + ptSysPy*ptSysPy);
	       double ht = lepton0.Pt() + lepton1.Pt() + jet.Pt() + metPt; 
	       
	       if (ptSystem <= ptsysCut){
		 
		 if (ht > htMin || mode !=0){
		   
		   if (nJets == 1 && nTightJetsBT == 1 && bTagged && nJetsBT == 1)histo_R->Fill(1, xlWeight);
		   if (nJets == 1 && nTightJetsBT == 2)  histo_R->Fill(2, xlWeight);
		   if (nJets == 1 && nTightJetsBT > 0)  histo_R->Fill(3, xlWeight);
		   if (nJets == 1 && nTightJetsBT > 1)  histo_R->Fill(4, xlWeight);
		   if (nJets == 2 && nTightJetsBT == 0)  histo_R->Fill(5, xlWeight);
		   if (nJets == 2 && nTightJetsBT == 1)  histo_R->Fill(6, xlWeight);
		   if (nJets == 2 && nTightJetsBT == 2)  histo_R->Fill(7, xlWeight);
		   if (nJets == 2 && nTightJetsBT > 0)  histo_R->Fill(8, xlWeight);
		   if (nJets == 2 && nTightJetsBT > 1)  histo_R->Fill(9, xlWeight);
		   if (nJets > 1 && nTightJetsBT == 0)  histo_R->Fill(10, xlWeight);
		   if (nJets > 1 && nTightJetsBT == 1)  histo_R->Fill(11, xlWeight);
		   if (nJets > 1 && nTightJetsBT == 2)  histo_R->Fill(12, xlWeight);
		   if (nJets > 1 && nTightJetsBT !=0 )  histo_R->Fill(13, xlWeight);
		   if (nJets > 1 && nTightJetsBT > 1 )  histo_R->Fill(14, xlWeight);
		   if (nJets == 3 && nTightJetsBT ==3 )  histo_R->Fill(15, xlWeight);
		 }
	       }
	     }
	   } //tt CR
	 } // mll
      } //mll pre
    } // 2 leptons
  } // events
  
  
  if (!silent){ 
    cout << "------------------------------------------" << endl;
    cout << "Results: " << plotName <<  endl;
    cout << "------------------------------------------" << endl;  
    for (int i = 2; i < 9; i++){
      if (i == 2) cout << " leptons: " <<  histo->GetBinContent(i) << " +/-  " <<  histo->GetBinError(i)  << endl;
      if (i == 3) cout << " inv. mass: " <<  histo->GetBinContent(i) << " +/-  " <<  histo->GetBinError(i)  << endl;
      if (i == 4) cout << " met: " <<  histo->GetBinContent(i) << " +/-  " <<  histo->GetBinError(i)  << endl;
      if (i == 5) cout << " jet: " <<  histo->GetBinContent(i) << " +/-  " <<  histo->GetBinError(i)  << endl;
      if (i == 6) cout << " jet_bt: " <<  histo->GetBinContent(i) << " +/-  " <<  histo->GetBinError(i)  << endl;
      if (i == 7) cout << " pt system: " <<  histo->GetBinContent(i) << " +/-  " <<  histo->GetBinError(i)  << endl;
      if (i == 8) cout << " ht: " <<  histo->GetBinContent(i) << " +/-  " <<  histo->GetBinError(i)  << endl;
    }
    cout << "------------------------------------------" << endl; 
  }
  f_var.Write();
  f_var.Close();
  
  
}


