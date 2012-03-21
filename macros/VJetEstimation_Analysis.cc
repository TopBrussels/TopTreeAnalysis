#include "TStyle.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

// Root headers
#include "TArrow.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include  "TGraphErrors.h"
#include "TMath.h"
#include "TMarker.h"
#include "TPaveStats.h"
#include "TRandom3.h"

// RooFit headers
#include "RooAddition.h"
#include "RooConstVar.h"
#include "RooCategory.h"
#include "RooCatType.h"
#include "RooCBShape.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFormula.h"
#include "RooFormulaVar.h"
#include "RooHistPdf.h"
#include "RooKeysPdf.h"
#include "RooPolynomial.h"
#include "RooRealVar.h"
#include "RooSimultaneous.h"

#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooProdPdf.h"

#include "RooMCStudy.h"
#include "RooMinuit.h"
#include "RooMinimizer.h"
#include "RooNLLVar.h"

#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooTable.h"
#include "Roo1DTable.h"

#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/ConfInterval.h"
#include "RooStats/PointSetInterval.h"
#include "RooStats/ConfidenceBelt.h"
#include "RooStats/FeldmanCousins.h"
#include "RooStats/ModelConfig.h"

#include "RooStats/HypoTestInverterOriginal.h"
#include "RooStats/HypoTestInverterResult.h"
#include "RooStats/HypoTestInverterPlot.h"
#include "RooStats/HybridCalculatorOriginal.h"

#include "Style.C"

using namespace RooFit;
using namespace RooStats ;
using namespace std;
//using namespace TopTree;

int main (int argc, char *argv[])
{

  clock_t start = clock();

  cout << "*************************************************************" << endl;
  cout << " Beginning of the program for the VJetEstimation analysis ! " << endl;
  cout << "*************************************************************" << endl;

  //SetStyle if needed
  //setTDRStyle();
  setGregStyle();
  //setMyStyle();
  
  unsigned int btag_wp_idx = 0; // B-tagging working point index
  unsigned int njets = 4;       // Nb of selected jets
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////Configuration ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //Input ROOT file
  char inputname[100];
	sprintf(inputname,"VJetEstimation_RooFit_WS_%d_wp_%d_jets.root",btag_wp_idx,njets);
  TFile *fin = TFile::Open(inputname);

  //Retrieve workspace from file
  sprintf(inputname,"w_%d_wp_%d_jets",btag_wp_idx,njets);
  RooWorkspace* w = (RooWorkspace*) fin->Get(inputname) ;
  w->Print("t") ;
  //Output ROOT file
  string channelpostfix = "_SemiMuon";
  string postfix = "_Analysis";
  string rootFileName ("VJetEstimation"+postfix+channelpostfix+".root");
  TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");
  TDirectory *myDir = 0;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////// Analysis ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // R e t r i e v e   p d f ,   d a t a   f r o m   w o r k s p a c e
  // -----------------------------------------------------------------

  // Retrieve x,model and data from workspace
	RooRealVar* Ntt   = w->var("Ntt") ; float Ntt_true = Ntt->getVal(); //Ntt->setVal(Ntt_true-1000);
	RooRealVar* Nv    = w->var("Nv") ;  float Nv_true = Nv->getVal();   //Nv->setVal(Nv_true+1000);
	RooRealVar* eb    = w->var("eb") ;    eb->setConstant(kFALSE);
	//RooRealVar* eudsc = w->var("eudsc") ; eudsc->setConstant(kTRUE);
	//RooRealVar* euds  = w->var("euds") ;  euds->setConstant(kTRUE);
	RooCategory* nbjets = (RooCategory*)w->cat("nbjets") ;
  if(verbose){
  	cout<<" MC Values for : "<<endl;
	  cout<<" Ntt = "<<Ntt_true<<endl;
	  cout<<" Nv  = "<<Nv_true<<endl;
  }

	sprintf(inputname,"model_%d_wp_%d_jets",btag_wp_idx,njets);
  RooAbsPdf* model = w->pdf(inputname) ;
  RooAbsData* data = w->data("data") ;
  //data->table(*nbjets)->Print("v");
  // Print structure of composite p.d.f.
  //model->Print("t") ;

  // Contrainsts on free parameters

  //RooGaussian eb_const("eb_const","eb_const",eb,RooConst(btag_eff),RooConst(btag_err));
  //RooGaussian eudsc_const("eudsc_const","eudsc_const",eudsc,RooConst(mistag_eff),RooConst(mistag_err));
  //RooGaussian euds_const("euds_const","euds_const",euds,RooConst(mistag_p_eff),RooConst(mistag_p_err));

  //RooProdPdf model_constraint("model_constraint","model with constraint",RooArgSet(model,eb_const)) ;
  //RooProdPdf model_constraint("model_constraint","model with constraint",RooArgSet(model,eudsc_const)) ;
  //RooProdPdf model_constraint("model_constraint","model with constraint",RooArgSet(model,euds_const)) ;
  //RooProdPdf model_constraint("model_constraint","model with constraint",RooArgSet(model,eb_const,eudsc_const,euds_const)) ;

  // F i t   m o d e l   t o   d a t a ,   p l o t   m o d e l 
  // ---------------------------------------------------------

  RooMCStudy* mcstudy = new RooMCStudy(*model,*nbjets,Binned(kTRUE),Silence(),Extended(kTRUE),FitOptions(Save(0),Optimize(0),Minos(1),Extended(1),Strategy(2),PrintEvalErrors(-1)));
  //RooMCStudy* mcstudy = new RooMCStudy(model_constraint,nbjets,Constrain(RooArgSet(eb_const,eudsc_const,euds_const)),Binned(kTRUE),Silence(),Extended(kTRUE),FitOptions(Save(0),Minos(1),Extended(1),Strategy(2),PrintEvalErrors(-1))) ;
  //RooMCStudy* mcstudy = new RooMCStudy(model,nbjets,ExternalConstraints(eb_const),Binned(kTRUE),Silence(),Extended(kTRUE),FitOptions(Save(0),Minos(1),Extended(1),Strategy(2),PrintEvalErrors(-1))) ;
  //RooMCStudy* mcstudy = new RooMCStudy(model_constraint,nbjets,ExternalConstraints(RooArgSet(eb_const,eudsc_const,euds_const)),Binned(kTRUE),Silence(),Extended(kTRUE),FitOptions(Save(0),Minos(1),Extended(1),Strategy(2),PrintEvalErrors(-1))) ;

  // G e n e r a t e   a n d   f i t   e v e n t s
  // ---------------------------------------------

  unsigned int nbOfPE = 10;                    // Nb of pseudo-exp
  unsigned int nExpected = data->numEntries(); // Sum of Nttlike and Nvlike
  // Generate and fit nbOfPE samples of Poisson(nExpected) events
  mcstudy->generateAndFit(nbOfPE,nExpected) ;

  // P l o t   f i t   r e s u l t s 
  // ---------------------------------------------------------------------

  unsigned int pullNbOfBins = 100;
  
  fout->cd();
  myDir = fout->mkdir("Pulls");

  TCanvas* c1_1 = new TCanvas("myCanva_Nttpull","",600,600);
  c1_1->cd();
  TH1* hNttpull = mcstudy->fitParDataSet().createHistogram("Nttpull",pullNbOfBins) ;
  hNttpull->Draw("E1");
  hNttpull->SetMarkerStyle(20);
  hNttpull->Fit("gaus","IEM");
  TF1* fit_Nttpull = (TF1*)hNttpull->FindObject("gaus");
  fit_Nttpull->SetLineWidth(3);
  fit_Nttpull->SetLineColor(kBlue);
  hNttpull->GetXaxis()->SetTitle("Pull");
  hNttpull->GetYaxis()->SetTitle("Number of events");
  hNttpull->GetYaxis()->SetTitleOffset(1.35);
  gPad->Update();
  TPaveStats* stat_Nttpull = (TPaveStats*)hNttpull->FindObject("stats");
  stat_Nttpull->SetX1NDC(0.60);
  stat_Nttpull->SetX2NDC(0.98);
  stat_Nttpull->SetY1NDC(0.81);
  stat_Nttpull->SetY2NDC(0.98);
  stat_Nttpull->Draw("same");
  c1_1->Write();

  // Closing files
  fin->Close();
  fout->Close();

  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;

  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "********************************************" << endl;

  return 0;
}

