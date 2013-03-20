#include "isis_looper.h"
#include <TH2.h>
#include <TH3.h>
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


void isis_looper::Loop(){
  // running default loop:
  myLoop(0,0,0);
}

void isis_looper::myLoop(int nsel, int mode, bool silent)
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
  
  /*
  //JER  
  else if (nsel == -10)                   {sprintf(plotName,"tt");}
  else if (nsel ==  10)                   {sprintf(plotName,"tt");}
  else if (nsel == -20)                   {sprintf(plotName,"twdr");}
  else if (nsel ==  20)                   {sprintf(plotName,"twdr");}
  
  //JES
  else if (nsel == -11)                   {sprintf(plotName,"tt");}
  else if (nsel ==  11)                   {sprintf(plotName,"tt");}
  else if (nsel == -21)                   {sprintf(plotName,"twdr");}
  else if (nsel ==  21)                   {sprintf(plotName,"twdr");}
  
    //PU
  else if (nsel == -12)                   {sprintf(plotName,"tt");}
  else if (nsel ==  12)                   {sprintf(plotName,"tt");}
  else if (nsel == -22)                   {sprintf(plotName,"twdr");}
  else if (nsel ==  22)                   {sprintf(plotName,"twdr");}
  
  //SF
  else if (nsel == -13)                   {sprintf(plotName,"tt");}
  else if (nsel ==  13)                   {sprintf(plotName,"tt");}
  else if (nsel == -23)                   {sprintf(plotName,"twdr");}
  else if (nsel ==  23)                   {sprintf(plotName,"twdr");}

  //MET
  else if (nsel == -14)                   {sprintf(plotName,"tt");}
  else if (nsel ==  14)                   {sprintf(plotName,"tt");}
  else if (nsel == -24)                   {sprintf(plotName,"twdr");}
  else if (nsel ==  24)                   {sprintf(plotName,"twdr");}
  */
  
 
  
  
  char newRootFile[300];
  double lumi = luminosity; 
    if      (mode == 0)	 lumi = 11966.617;  	
    else if (mode == 1) lumi = 12067.294;  	
    else if (mode == 2) lumi = 12093.792;  	
  
  
  if(nsel == -10 || nsel == -20){
  sprintf(newRootFile,"results/JERsysDown_an_%dpb_%d.root", (int)lumi, mode);
  }
  else if(nsel == 10 || nsel == 20){
   sprintf(newRootFile,"results/JERsysUp_an_%dpb_%d.root", (int)lumi, mode);
  }
  else if(nsel == -11 || nsel == -21 ){
   sprintf(newRootFile,"results/JESsysDown_an_%dpb_%d.root", (int)lumi, mode);
  }
  else if(nsel == 11 || nsel == 21 ){
   sprintf(newRootFile,"results/JESsysUp_an_%dpb_%d.root", (int)lumi, mode);
  }
  else if(nsel == -12 || nsel == -22 ){
   sprintf(newRootFile,"results/PUsysDown_an_%dpb_%d.root", (int)lumi, mode);
  }
  else if(nsel == 12 || nsel == 22 ){
   sprintf(newRootFile,"results/PUsysUp_an_%dpb_%d.root", (int)lumi, mode);
  }
  else if(nsel == -13 || nsel == -23 ){
   sprintf(newRootFile,"results/SFsysDown_an_%dpb_%d.root", (int)lumi, mode);
  }
  else if(nsel == 13 || nsel == 23 ){
   sprintf(newRootFile,"results/SFsysUp_an_%dpb_%d.root", (int)lumi, mode);
  }
  else if(nsel == -14 || nsel == -24 ){
   sprintf(newRootFile,"results/METsysDown_an_%dpb_%d.root", (int)lumi, mode);
  }
  else if(nsel == 14 || nsel == 24 ){
   sprintf(newRootFile,"results/METsysUp_an_%dpb_%d.root", (int)lumi, mode);
  }
  else{
  
    sprintf(newRootFile,"results/an_%dpb_%d.root", (int)lumi, mode);
  }
 
 
 
  TFile f_var(newRootFile, "UPDATE");
  
  
  if(!silent){
    std::cout << "[Info:] results root file " << newRootFile << std::endl;
  }
  
  
  //////////
  char title[300];
  
  sprintf(title,"3d_btagged_tightjets_%s",plotName);
  TH3F* histo_3d_btagged_tightjets = new TH3F(title, " Btagged  jets",  10,  -0.5, 9.5,   10,  -0.5, 9.5,10,  -0.5, 9.5  ); 
  histo_3d_btagged_tightjets->Sumw2();
 
  sprintf(title,"2d_btagged_tightjets_noHt_%s",plotName);
  	TH2F* histo_2d_btagged_tightjets_noHt = new TH2F(title, " Btagged  jets - no Ht cut",  10,  -0.5, 9.5,10,  -0.5, 9.5  );
  histo_2d_btagged_tightjets_noHt->Sumw2();
  
  sprintf(title,"cuts_%s",plotName);
  TH1F* histo = new TH1F( title, " ", 10,  0, 10 );
  histo->Sumw2();
  
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
  
  sprintf(title,"mll_after_%s",plotName);
  TH1F* histo_mll_after = new TH1F( title, " ", 100,  0, 200 );
  histo_mll_after->Sumw2();
  
  sprintf(title,"njets_cut_%s",plotName);
  TH1F* histo_njets_cut = new TH1F( title, " ", 10,  -0.5, 9.5 );
  histo_njets_cut->Sumw2();
  
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
  
  sprintf(title,"ptsys_high_%s",plotName);
  TH1F* histo_ptsys_high = new TH1F( title, " ", 100,  0, 200 );
  histo_ptsys_high->Sumw2();
  
  sprintf(title,"ptsys_low_%s",plotName);
  TH1F* histo_ptsys_low = new TH1F( title, " ", 100,  0, 200 );
  histo_ptsys_low->Sumw2();
  
  sprintf(title,"ht_high_%s",plotName);
  TH1F* histo_ht_high = new TH1F( title, " ", 300,  0, 600 );
  histo_ht_high->Sumw2();
  
  sprintf(title,"ht_low_%s",plotName);
  TH1F* histo_ht_low = new TH1F( title, " ", 300,  0, 600 );
  histo_ht_low->Sumw2();
  
  sprintf(title,"ht_cut_%s",plotName);
  TH1F* histo_ht_cut = new TH1F( title, " ", 300,  0, 600 );
  histo_ht_cut->Sumw2();
  
   sprintf(title,"ht_nomet_high_%s",plotName);
  TH1F* histo_ht_nomet_high = new TH1F( title, " ", 300,  0, 600 );
  histo_ht_nomet_high->Sumw2();
  
  sprintf(title,"ht_nomet_low_%s",plotName);
  TH1F* histo_ht_nomet_low = new TH1F( title, " ", 300,  0, 600 );
  histo_ht_nomet_low->Sumw2();
  
  sprintf(title,"ht_nomet_cut_%s",plotName);
  TH1F* histo_ht_nomet_cut = new TH1F( title, " ", 300,  0, 600 );
  histo_ht_nomet_cut->Sumw2();
  
  sprintf(title,"pt_max_%s",plotName);
  TH1F* histo_pt_max = new TH1F( title, " ", 100,  0, 200 );
  histo_pt_max->Sumw2();
  
  sprintf(title,"pt_min_%s",plotName);
  TH1F* histo_pt_min = new TH1F( title, " ", 100,  0, 200 );
  histo_pt_min->Sumw2();
  
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
  
   sprintf(title,"ht_nomet_bf_%s",plotName);
  TH1F* histo_ht_nomet_bf = new TH1F( title, " ", 300,  0, 600 );
  histo_ht_nomet_bf->Sumw2();
  
  sprintf(title,"ptsys_bf_%s",plotName);
  TH1F* histo_ptsys_bf = new TH1F( title, " ", 100,  0, 200 );
  histo_ptsys_bf->Sumw2();
  
  sprintf(title,"npu_%s",plotName);
  TH1F* histo_npu = new TH1F( title, " ", 30,   -0.5, 29.5 );
  histo_npu->Sumw2();
  
  
  /// Classic plotmaker plots
  sprintf(title,"met_%s",plotName);
  TH1F* histo_met = new TH1F( title, " ", 100,  0, 200 );
  histo_met->Sumw2();
  
  sprintf(title,"mll_%s",plotName);
  TH1F* histo_mll = new TH1F( title, " ", 100,  0, 200 );
  histo_mll->Sumw2();
  
  sprintf(title,"njets_%s",plotName);
  TH1F* histo_njets = new TH1F( title, " ", 10,  0, 10 );
  histo_njets->Sumw2();
  
  sprintf(title,"njetsbt_%s",plotName);
  TH1F* histo_njetsbt = new TH1F( title, " ", 10,  -0.5, 9.5 );
  histo_njetsbt->Sumw2();
  
  sprintf(title,"ptsys_%s",plotName);
  TH1F* histo_ptsys = new TH1F( title, " ", 100,  0, 200 );
  histo_ptsys->Sumw2();
  
  sprintf(title,"ht_%s",plotName);
  TH1F* histo_ht = new TH1F( title, " ", 300,  0, 600 );
  histo_ht->Sumw2();
  
    
  sprintf(title,"ht_nomet_%s",plotName);
  TH1F* histo_ht_nomet = new TH1F( title, " ", 300,  0, 600 );
  histo_ht_nomet->Sumw2();
  
  sprintf(title,"pt_leading_%s",plotName);
  TH1F* histo_pt_leading = new TH1F( title, " ", 100,  0, 200 );
  histo_pt_leading->Sumw2();
  
  
  sprintf(title,"eta_leading_%s",plotName);
  TH1F* histo_eta_leading = new TH1F( title, " ", 101,  -3, 3);
  histo_eta_leading->Sumw2();
  
  sprintf(title,"nvertex_%s",plotName);
  TH1F* histo_nvertex = new TH1F( title, " ", 70,   -0.5, 69.5 );
  histo_nvertex->Sumw2();
 
  
  // 1 jet level
  /// Classic plotmaker plots
  sprintf(title,"met_1j_%s",plotName);
  TH1F* histo_met_1j = new TH1F( title, " ", 100,  0, 200 );
  histo_met_1j->Sumw2();
  
  sprintf(title,"mll_1j_%s",plotName);
  TH1F* histo_mll_1j = new TH1F( title, " ", 100,  0, 200 );
  histo_mll_1j->Sumw2();
  
  sprintf(title,"ptsys_1j_%s",plotName);
  TH1F* histo_ptsys_1j = new TH1F( title, " ", 100,  0, 200 );
  histo_ptsys_1j->Sumw2();
  
  sprintf(title,"ht_1j_%s",plotName);
  TH1F* histo_ht_1j = new TH1F( title, " ", 300,  0, 600 );
  histo_ht_1j->Sumw2();
  
   sprintf(title,"ht_nomet_1j_%s",plotName);
  TH1F* histo_ht_nomet_1j = new TH1F( title, " ", 300,  0, 600 );
  histo_ht_nomet_1j->Sumw2();
  
  sprintf(title,"pt_leading_1j_%s",plotName);
  TH1F* histo_pt_leading_1j = new TH1F( title, " ", 100,  0, 200 );
  histo_pt_leading_1j->Sumw2();
  
  sprintf(title,"nvertex_1j_%s",plotName);
  TH1F* histo_nvertex_1j = new TH1F( title, " ", 30,   -0.5, 29.5 );
  histo_nvertex_1j->Sumw2();
  
  // 1 jet level
  /// Classic plotmaker plots
  sprintf(title,"met_1j1t_%s",plotName);
  TH1F* histo_met_1j1t = new TH1F( title, " ", 100,  0, 200 );
  histo_met_1j1t->Sumw2();
  
  sprintf(title,"mll_1j1t_%s",plotName);
  TH1F* histo_mll_1j1t = new TH1F( title, " ", 100,  0, 200 );
  histo_mll_1j1t->Sumw2();
  
  sprintf(title,"ptsys_1j1t_%s",plotName);
  TH1F* histo_ptsys_1j1t = new TH1F( title, " ", 100,  0, 200 );
  histo_ptsys_1j1t->Sumw2();
  
  sprintf(title,"ht_1j1t_%s",plotName);
  TH1F* histo_ht_1j1t = new TH1F( title, " ", 300,  0, 600 );
  histo_ht_1j1t->Sumw2();
  
   sprintf(title,"ht_nomet_1j1t_%s",plotName);
  TH1F* histo_ht_nomet_1j1t = new TH1F( title, " ", 300,  0, 600 );
  histo_ht_nomet_1j1t->Sumw2();
  
  sprintf(title,"pt_leading_1j1t_%s",plotName);
  TH1F* histo_pt_leading_1j1t = new TH1F( title, " ", 100,  0, 200 );
  histo_pt_leading_1j1t->Sumw2();
  
  sprintf(title,"nvertex_1j1t_%s",plotName);
  TH1F* histo_nvertex_1j1t = new TH1F( title, " ", 30,   -0.5, 29.5 );
  histo_nvertex_1j1t->Sumw2();
  
  
  // 2j1t
  /// Classic plotmaker plots
  sprintf(title,"met_2j1t_%s",plotName);
  TH1F* histo_met_2j1t = new TH1F( title, " ", 100,  0, 200 );
  histo_met_2j1t->Sumw2();
  
  sprintf(title,"mll_2j1t_%s",plotName);
  TH1F* histo_mll_2j1t = new TH1F( title, " ", 100,  0, 200 );
  histo_mll_2j1t->Sumw2();
  
  sprintf(title,"ptsys_2j1t_%s",plotName);
  TH1F* histo_ptsys_2j1t = new TH1F( title, " ", 100,  0, 200 );
  histo_ptsys_2j1t->Sumw2();
  
  sprintf(title,"ht_2j1t_%s",plotName);
  TH1F* histo_ht_2j1t = new TH1F( title, " ", 300,  0, 600 );
  histo_ht_2j1t->Sumw2();
  
    sprintf(title,"ht_nomet_2j1t_%s",plotName);
  TH1F* histo_ht_nomet_2j1t = new TH1F( title, " ", 300,  0, 600 );
  histo_ht_nomet_2j1t->Sumw2();
  
  sprintf(title,"pt_leading_2j1t_%s",plotName);
  TH1F* histo_pt_leading_2j1t = new TH1F( title, " ", 100,  0, 200 );
  histo_pt_leading_2j1t->Sumw2();
  
  sprintf(title,"nvertex_2j1t_%s",plotName);
  TH1F* histo_nvertex_2j1t = new TH1F( title, " ", 30,   -0.5, 29.5 );
  histo_nvertex_2j1t->Sumw2();
  
  
  // 2j2t
  /// Classic plotmaker plots
  sprintf(title,"met_2j2t_%s",plotName);
  TH1F* histo_met_2j2t = new TH1F( title, " ", 100,  0, 200 );
  histo_met_2j2t->Sumw2();
  
  sprintf(title,"mll_2j2t_%s",plotName);
  TH1F* histo_mll_2j2t = new TH1F( title, " ", 100,  0, 200 );
  histo_mll_2j2t->Sumw2();
  
  sprintf(title,"ptsys_2j2t_%s",plotName);
  TH1F* histo_ptsys_2j2t = new TH1F( title, " ", 100,  0, 200 );
  histo_ptsys_2j2t->Sumw2();
  
  sprintf(title,"ht_2j2t_%s",plotName);
  TH1F* histo_ht_2j2t = new TH1F( title, " ", 300,  0, 600 );
  histo_ht_2j2t->Sumw2();
  
   sprintf(title,"ht_nomet_2j2t_%s",plotName);
  TH1F* histo_ht_nomet_2j2t = new TH1F( title, " ", 300,  0, 600 );
  histo_ht_nomet_2j2t->Sumw2();
  
  sprintf(title,"pt_leading_2j2t_%s",plotName);
  TH1F* histo_pt_leading_2j2t = new TH1F( title, " ", 100,  0, 200 );
  histo_pt_leading_2j2t->Sumw2();
  
  sprintf(title,"nvertex_2j2t_%s",plotName);
  TH1F* histo_nvertex_2j2t = new TH1F( title, " ", 30,   -0.5, 29.5 );
  histo_nvertex_2j2t->Sumw2();

  
  // all regions
  sprintf(title,"R_%s",plotName);
  TH1F* histo_R = new TH1F( title, " ", 40,  0.5, 40.5 );
  histo_R->Sumw2();
  
  
  // checking pu reweighting
  sprintf(title,"nvertex_final_%s",plotName);
  TH1F* histo_nvertex_final = new TH1F( title, " ", 70,   -0.5, 69.5 );
  histo_nvertex_final->Sumw2();
  
  sprintf(title,"nvertex_final_3D_%s",plotName);
  TH1F* histo_nvertex_final_3D = new TH1F( title, " ", 70,   -0.5, 69.5 );
  histo_nvertex_final_3D->Sumw2();
  
  sprintf(title,"nvertex_final_purw_%s",plotName);
  TH1F* histo_nvertex_final_purw = new TH1F( title, " ", 70,   -0.5, 69.5 );
  histo_nvertex_final_purw->Sumw2();
  


  sprintf(title,"nvertex_2lep_%s",plotName);
  TH1F* histo_nvertex_2lep = new TH1F( title, " ", 70,   -0.5, 69.5 );
  histo_nvertex_2lep->Sumw2();
//<<<<<<< isis_looper.C
  
  //_________________________________________________________________________________
  //added by Isis
  
  // -- Leading jet --- 
    sprintf(title,"eta_jet_%s",plotName);
  TH1F* histo_eta_jet = new TH1F( title, " ", 50,-5, 5 );  //#bins - begin - end
  histo_eta_jet->Sumw2();
  
 
     sprintf(title,"eta_jet_110_%s",plotName);
  TH1F* histo_eta_jet_110 = new TH1F( title, " ", 50,-5, 5 );  //#bins - begin - end
  histo_eta_jet_110->Sumw2();
  
      sprintf(title,"eta_jet_90_%s",plotName);
  TH1F* histo_eta_jet_90 = new TH1F( title, " ", 50,-5, 5 );  //#bins - begin - end
  histo_eta_jet_90->Sumw2();
  
  
      sprintf(title,"eta_jet_70_%s",plotName);
  TH1F* histo_eta_jet_70 = new TH1F( title, " ", 50,-5, 5 );  //#bins - begin - end
  histo_eta_jet_70->Sumw2();
  
  
      sprintf(title,"eta_jet_50_%s",plotName);
  TH1F* histo_eta_jet_50 = new TH1F( title, " ", 50,-5, 5 );  //#bins - begin - end
  histo_eta_jet_50->Sumw2();
  
      sprintf(title,"eta_jet_30_%s",plotName);
  TH1F* histo_eta_jet_30 = new TH1F( title, " ", 50,-5, 5 );  //#bins - begin - end
  histo_eta_jet_30->Sumw2();
 
 
   // -- Second Leading jet --- 
    sprintf(title,"eta_jet1_%s",plotName);
  TH1F* histo_eta_jet1 = new TH1F( title, " ", 50,-5, 5 );  //#bins - begin - end
  histo_eta_jet1->Sumw2();
  
 
     sprintf(title,"eta_jet1_110_%s",plotName);
  TH1F* histo_eta_jet1_110 = new TH1F( title, " ", 50,-5, 5 );  //#bins - begin - end
  histo_eta_jet1_110->Sumw2();
  
      sprintf(title,"eta_jet1_90_%s",plotName);
  TH1F* histo_eta_jet1_90 = new TH1F( title, " ", 50,-5, 5 );  //#bins - begin - end
  histo_eta_jet1_90->Sumw2();
  
  
      sprintf(title,"eta_jet1_70_%s",plotName);
  TH1F* histo_eta_jet1_70 = new TH1F( title, " ", 50,-5, 5 );  //#bins - begin - end
  histo_eta_jet1_70->Sumw2();
  
  
      sprintf(title,"eta_jet1_50_%s",plotName);
  TH1F* histo_eta_jet1_50 = new TH1F( title, " ", 50,-5, 5 );  //#bins - begin - end
  histo_eta_jet1_50->Sumw2();
  
      sprintf(title,"eta_jet1_30_%s",plotName);
  TH1F* histo_eta_jet1_30 = new TH1F( title, " ", 50,-5, 5 );  //#bins - begin - end
  histo_eta_jet1_30->Sumw2();
  
  // ---- Z/gamma controlregion 
 
    sprintf(title,"mll_zgamma_%s",plotName);
  TH1F* histo_mll_zgamma = new TH1F( title, " ", 50,  0, 200 );
  histo_mll_zgamma->Sumw2();
  
  
    sprintf(title,"met_zgamma_%s",plotName);
  TH1F* histo_met_zgamma= new TH1F( title, " ", 50,  0, 200 );
  histo_met_zgamma->Sumw2();
  
      sprintf(title,"mll_outside_%s",plotName);
  TH1F* histo_mll_outside = new TH1F( title, " ", 50,  0, 200 );
  histo_mll_outside->Sumw2();
  
  
    sprintf(title,"met_outside_%s",plotName);
  TH1F* histo_met_outside= new TH1F( title, " ", 50,  0, 200 );
  histo_met_outside->Sumw2();
  
//=======
//>>>>>>> 1.1.2.3



  

//__________________________________________________ END HISTO DEF _________________________________________________________


  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {  // start eventloop
    
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    
    if (lumi != lum && nsel != 666 && mode !=3){
      if (jentry
       == 0)std::cout << "[Warning:] This tree was made with a different luminosity (" << lum << ") than " << lumi << std::endl;
      //xlWeight*=(lumi/lum);
    }
    
  //  if(puweight > 5){ continue ;}
    double ht_temp = 0; 
    double ptsysX_temp = 0; 
    double ptsysY_temp = 0; 
    
    if(ptLepton->size() != 2){
      std::cout << "[Warning:] Something is wrong, your Tree is not correctly filled" << std::endl;
      break;
    } else {
      histo->Fill(0.,xlWeight);
      histo_nvertex_2lep->Fill(nvertex,xlWeight);
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
        
       

	int nJetsBT = 0;
	int nTightJetsBT = 0;
	int nJets = 0;	
	bool bTagged = false;;
	
	if( ptJet->size() != Btagjet->size()){
		cout << "ERROR: something went wrong with btagging " << endl; 
	}
	
	for (unsigned int i =0; i < ptJet->size(); i ++){ 
	  TLorentzVector tempJet(pxJet->at(i),pyJet->at(i), pzJet->at(i), eJet->at(i));
	  bool btag = Btagjet->at(i);
	  
	 
	 // cout << "jet: " << i << "btagjet " << btag << endl; 
	  
          	
	  if (ptJet->at(i) > 30 && fabs(etaJet->at(i)) < 2.4 && TMath::Min(fabs(lepton0.DeltaR(tempJet)), fabs(lepton1.DeltaR(tempJet))) > 0.3) {
		ht_temp += ptJet->at(i);
		ptsysX_temp += pxJet->at(i);
		ptsysY_temp += pyJet->at(i); 
		
	    	nJets++;

	   	 if (btag){
	        	nTightJetsBT++;
	                nJetsBT++;
	    	} // end  if 
	  } // end ptJet->at(i) > 30 && TMath::Min(fabs(lepton0.DeltaR(tempJet)), fabs(lepton1.DeltaR(tempJet))) > 0.3
	  else if (btag && ptJet->at(i) > 20 && fabs(etaJet->at(i)) < 2.4){
	  	bTagged = true;
	  	nJetsBT++;
	  
	  }
	} // end for (unsigned int i =0; i < ptJet->size(); i ++)
	
	//cout << " number of btagged jets: " << nJetsBT << " tight: " << nTightJetsBT << endl; 


	histo_pt_max->Fill(TMath::Max(lepton0.Pt(), lepton1.Pt()), xlWeight);
	histo_pt_min->Fill(TMath::Min(lepton0.Pt(), lepton1.Pt()), xlWeight);
	histo_njets->Fill(nJets,  xlWeight);
	histo_njetsbt->Fill(nJetsBT,  xlWeight);
	histo_mll->Fill(pair.M(),  xlWeight);
	histo_met->Fill(metPt,  xlWeight);
	histo_promet->Fill(promet, xlWeight);
	
	
	if (nJets) histo_pt_leading->Fill(ptJet->at(0), xlWeight);
	
	if (nJets){
	  TLorentzVector jet_aux(pxJet->at(0),pyJet->at(0), pzJet->at(0), eJet->at(0));
	  histo_eta_leading->Fill(jet_aux.Eta(), xlWeight);
	}
	
	if (nJets == 1){
	  histo_etalepton->Fill(lepton0.Eta(), xlWeight);
	} // end (nJets == 1)
	
	bool invMass = false;
	if      (mode == 0) invMass = true;
	else if (mode == 1  && (pair.M() > invMax || pair.M() < invMin)) invMass = true;
	else if (mode == 2 && (pair.M() > invMax || pair.M() < invMin)) invMass = true;
	
	//____________________________
        //Z/gamma control region
        
	//if (pair.M() > 81 || pair.M() < 101)
	if(!invMass || (mode == 0 && (pair.M() > 81 || pair.M() < 101) )) 
	{	  
	         
		  histo_mll_zgamma->Fill(pair.M(),  xlWeight);
	          histo_met_zgamma->Fill(metPt,  xlWeight);
	          
	}else{
		histo_mll_outside->Fill(pair.M(),  xlWeight);
	          histo_met_outside->Fill(metPt,  xlWeight);
	
	}	 

	//______________________________________


           
	  
         


	if (invMass){
	  histo->Fill(2, xlWeight);
	  histo_mll_after->Fill(pair.M(),  xlWeight);
	  histo_met_cut->Fill(metPt,  xlWeight);
	  
	  
	  
	  
	  if (metPt >= 50 || mode ==0){
	   
	    histo->Fill(3, xlWeight);
	    histo_njets_cut->Fill(nJets, xlWeight);
	    
	    double ptSysPx = lepton0.Px() + lepton1.Px() + ptsysX_temp + metPx;
	    double ptSysPy = lepton0.Py() + lepton1.Py() + ptsysY_temp + metPy;
	    double ptSystem = sqrt(ptSysPx *ptSysPx + ptSysPy *ptSysPy);
	    double ht = lepton0.Pt() + lepton1.Pt() + ht_temp + metPt;
	    double ht_nomet = lepton0.Pt() + lepton1.Pt() + ht_temp ; 
	    
	    if(nJets==2 && nTightJetsBT == 1 && nJetsBT == 1 &&  (ht > 160 || mode !=0)) histo_R->Fill(2, xlWeight); // after all cuts 2j1t veto on loose
	    if(nJets==2 && nTightJetsBT ==2 && nJetsBT == 2 &&  (ht > 160 || mode !=0)) histo_R->Fill(3, xlWeight);  // after all cuts 2j2t veto on loose
	    
	    if(nJets==2 && nTightJetsBT == 1 && nJetsBT == 1) histo_R->Fill(5, xlWeight); // before 1jet cut and before ht> 160 cut: 2j1t veto on loose
	    if(nJets==2 && nTightJetsBT ==2 && nJetsBT == 2) histo_R->Fill(6, xlWeight); // before 1jet cut and before ht> 160 cut: 2j2t veto on loose
	    if(nJets==2 && nJetsBT == 1) histo_R->Fill(8, xlWeight);// before 1jet cut and before ht> 160 cut: 2j1t 
	    if(nJets==2 && nJetsBT == 2) histo_R->Fill(9, xlWeight);  // before 1jet cut and before ht> 160 cut: 2j2t 
	    
	    if(ht > 160 || mode !=0){
	    	histo_3d_btagged_tightjets->Fill(nJetsBT,nTightJetsBT,nJets, xlWeight);
	    }
	    histo_2d_btagged_tightjets_noHt->Fill(nTightJetsBT,nJets, xlWeight);
	    
	    
	    if (nJets == 1){
	      histo->Fill(4, xlWeight);
	      histo_njetsbt_cut->Fill(nJetsBT, xlWeight);
	      
	      TLorentzVector jet(pxJet->at(0),pyJet->at(0), pzJet->at(0), eJet->at(0));
		
	      
	      
	      if(nTightJetsBT == 1) {histo_R->Fill(7, xlWeight); }
	      
	      
	      
	      if (nJets == 1 && nTightJetsBT == 1 && nJetsBT == 1 ){
		histo->Fill(5, xlWeight);
	
	      
		histo_ptsys->Fill(ptSystem, xlWeight);
		histo_ht->Fill(ht, xlWeight);
		histo_ht_nomet->Fill(ht_nomet, xlWeight);
	
		
		histo_nvertex->Fill(nvertex, xlWeight);
		histo_npu->Fill(npu, xlWeight);
	
		
		histo_R->Fill(4, xlWeight);   // 1j1t1L no ht cut 
		
		if (ht > 160 || mode !=0){
		  histo->Fill(6, xlWeight);
		  histo_ht_cut->Fill(ht, xlWeight);
		  histo_ht_nomet_cut->Fill(ht_nomet, xlWeight);
		  
		  //Example to access the pu reweighting!
		  histo_nvertex_final->Fill(nvertex, rawWeight);
		  histo_nvertex_final_purw->Fill(nvertex, rawWeight*puweight);
		  
		  
		  
		  
		  histo_R->Fill(1, xlWeight); //signal passing all the cuts
		  
		} // end if (ht > htMin || mode !=0)
	      } // end  if (nJets == 1 && nTightJetsBT == 1 && bTagged && nJetsBT == 1)
	    } // end if (nJets == 1
	  } // end if (met cut)
	  

	} // mll
      } //mll pre
    } // 2 leptons
  }// event loop.
  
  
  
  if (!silent){ 
    cout << "------------------------------------------" << endl;
    cout << "[Results:] " << plotName <<  endl;
    cout << "------------------------------------------" << endl;  
    for (int i = 2; i < 9; i++){
      if (i == 2) cout << " leptons: " <<  histo->GetBinContent(i) << " +/-  " <<  histo->GetBinError(i)  << endl;
      if (i == 3) cout << " inv. mass: " <<  histo->GetBinContent(i) << " +/-  " <<  histo->GetBinError(i)  << endl;
      if (i == 4) cout << " met: " <<  histo->GetBinContent(i) << " +/-  " <<  histo->GetBinError(i)  << endl;
      if (i == 5) cout << " jet: " <<  histo->GetBinContent(i) << " +/-  " <<  histo->GetBinError(i)  << endl;
      if (i == 6) cout << " jet_bt: " <<  histo->GetBinContent(i) << " +/-  " <<  histo->GetBinError(i)  << endl;
      if (i == 7) cout << " ht: " <<  histo->GetBinContent(i) << " +/-  " <<  histo->GetBinError(i)  << endl;
    }
    
     cout << "------------------------------------------" << endl;
    cout << "[Results: regions - ] " << plotName <<  endl;
    cout << "------------------------------------------" << endl;  
     cout << " tight - loose jets are vetoed" << endl; 
      cout << " 1jet 1tag: " <<  histo_R->GetBinContent(1) << " +/-  " <<  histo_R->GetBinError(1)  << endl;
       cout << "2jets 1tag: " <<  histo_R->GetBinContent(2) << " +/-  " <<  histo_R->GetBinError(2)  << endl;
       cout << "2jets 2tags: " <<  histo_R->GetBinContent(3) << " +/-  " <<  histo_R->GetBinError(3)  << endl;
       cout << endl; 
    cout << "tight - veto on loose jets - no ht cut" << endl; 
     cout << " 1jet 1tag: " <<  histo_R->GetBinContent(4) << " +/-  " <<  histo_R->GetBinError(4)  << endl;
       cout << "2jets 1tag: " <<  histo_R->GetBinContent(5) << " +/-  " <<  histo_R->GetBinError(5)  << endl;
       cout << "2jets 2tags: " <<  histo_R->GetBinContent(6) << " +/-  " <<  histo_R->GetBinError(6)  << endl;
       cout << endl; 
    cout << "tight - no on loose jets - no ht cut" << endl; 
     cout << " 1jet 1tag: " <<  histo_R->GetBinContent(7) << " +/-  " <<  histo_R->GetBinError(7)  << endl;
       cout << "2jets 1tag: " <<  histo_R->GetBinContent(8) << " +/-  " <<  histo_R->GetBinError(8)  << endl;
       cout << "2jets 2tags: " <<  histo_R->GetBinContent(9) << " +/-  " <<  histo_R->GetBinError(9)  << endl;
       
 	cout << "--------------------------------------------------" << endl;
	cout << " btag check (no logistics) " << endl; 
	cout << "--------------------------------------------------" << endl;
	
	cout << "VETOING LOOSE JETS: "<< endl;
	cout << "1jet 1tag: " << histo_3d_btagged_tightjets->GetBinContent(2,2,2) << "+/-" << histo_3d_btagged_tightjets->GetBinError(2,2,2) << endl; 
	cout << "2jet 1tag: " << histo_3d_btagged_tightjets->GetBinContent(2,2,3) << "+/-" << histo_3d_btagged_tightjets->GetBinError(2,2,3) << endl;
	cout << "2jet 2tag: " << histo_3d_btagged_tightjets->GetBinContent(3,3,3) << "+/-" << histo_3d_btagged_tightjets->GetBinError(3,3,3) << endl;  
	cout << endl;  
        cout << " NO VETO ON LOOSE JETS - NO HT CUT " << endl; 
	cout << "1jet 1tag: " << histo_2d_btagged_tightjets_noHt->GetBinContent(2,2) << "+/-" << histo_2d_btagged_tightjets_noHt->GetBinError(2,2) << endl; 
	cout << "2jet 1tag: " << histo_2d_btagged_tightjets_noHt->GetBinContent(2,3) << "+/-" << histo_2d_btagged_tightjets_noHt->GetBinError(2,3) << endl;
	cout << "2jet 2tag: " << histo_2d_btagged_tightjets_noHt->GetBinContent(3,3) << "+/-" << histo_2d_btagged_tightjets_noHt->GetBinError(3,3) << endl;
  /*  
    cout << "------------------------------------------" << endl; 
    cout << "[eta values:]" << plotName << endl;
    cout << "------------------------------------------" << endl; 
    for (int j =1 ; j<7; j++){
      if(j == 1) {
         double amount = 0; 
	 //double amounterror = 0; 
      
         for(int k = 0; k< 50; k++) {   
	     amount = amount + histo_eta_jet_30->GetBinContent(k); 
	  //   amounterror = sqrt(amounterror^2 + (histo_eta_jet_30->GetBinError(k))^2);
           //cout << "pt higher then 30: " << histo_eta_jet_30->GetBinContent(k) << " +/- " << histo_eta_jet_30->GetBinError(k) << endl; 
         }
	 
	 cout << "pt higher then 30: " << amount <<  endl;
      
      }
        if(j == 2) {
         double amount = 0; 
	 //double amounterror = 0; 
      
         for(int k = 0; k< 50; k++) {   
	     amount = amount + histo_eta_jet_50->GetBinContent(k); 
	  //   amounterror = sqrt(amounterror^2 + (histo_eta_jet_30->GetBinError(k))^2);
           //cout << "pt higher then 50: " << histo_eta_jet_30->GetBinContent(k) << " +/- " << histo_eta_jet_30->GetBinError(k) << endl; 
         }
	 
	 cout << "pt higher then 50: " << amount <<  endl;
      
      }
      
            if(j == 3) {
         double amount = 0; 
	 //double amounterror = 0; 
      
         for(int k = 0; k< 50; k++) {   
	     amount = amount + histo_eta_jet_70->GetBinContent(k); 
	  //   amounterror = sqrt(amounterror^2 + (histo_eta_jet_30->GetBinError(k))^2);
           //cout << "pt higher then 30: " << histo_eta_jet_30->GetBinContent(k) << " +/- " << histo_eta_jet_30->GetBinError(k) << endl; 
         }
	 
	 cout << "pt higher then 70: " << amount <<  endl;
      
      }
            if(j == 4) {
         double amount = 0; 
	 //double amounterror = 0; 
      
         for(int k = 0; k< 50; k++) {   
	     amount = amount + histo_eta_jet_70->GetBinContent(k); 
	  //   amounterror = sqrt(amounterror^2 + (histo_eta_jet_30->GetBinError(k))^2);
           //cout << "pt higher then 30: " << histo_eta_jet_30->GetBinContent(k) << " +/- " << histo_eta_jet_30->GetBinError(k) << endl; 
         }
	 
	 cout << "pt higher then 70: " << amount <<  endl;
      
      }
            if(j == 5) {
         double amount = 0; 
	 //double amounterror = 0; 
      
         for(int k = 0; k< 50; k++) {   
	     amount = amount + histo_eta_jet_90->GetBinContent(k); 
	  //   amounterror = sqrt(amounterror^2 + (histo_eta_jet_30->GetBinError(k))^2);
           //cout << "pt higher then 30: " << histo_eta_jet_30->GetBinContent(k) << " +/- " << histo_eta_jet_30->GetBinError(k) << endl; 
         }
	 
	 cout << "pt higher then 90: " << amount <<  endl;
      
      }
            if(j == 6) {
         double amount = 0; 
	 //double amounterror = 0; 
      
         for(int k = 0; k< 50; k++) {   
	     amount = amount + histo_eta_jet_110->GetBinContent(k); 
	  //   amounterror = sqrt(amounterror^2 + (histo_eta_jet_30->GetBinError(k))^2);
           //cout << "pt higher then 30: " << histo_eta_jet_30->GetBinContent(k) << " +/- " << histo_eta_jet_30->GetBinError(k) << endl; 
         }
	 
	 cout << "pt higher then 110: " << amount <<  endl;
      
      }
    
    }
    
  */  
    
    
  }
  f_var.Write(newRootFile,TFile::kOverwrite);
  f_var.Close();
}
