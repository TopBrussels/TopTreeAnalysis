#include "../interface/ExpCorrCalculator.h"

ExpCorrCalculator::ExpCorrCalculator(string name)
{
  name_ = name;
  nPtBins_ = 4;
  wasCalculated_ = false;
  
  diffFunction_ = 0;
  dEbFunction_ = 0;
  
  Float_t jetPtBinning[] = {30.,60.,75.,100.,3000.};
  for(unsigned int i=0; i < nPtBins_+1; i++)
    jetPtBinning_[i] = jetPtBinning[i];
  
  string ptLightvsptBTitle = "PtLightVSPtB_" + name_;
  histo2D_[ptLightvsptBTitle] = new TH2F(ptLightvsptBTitle.c_str(),ptLightvsptBTitle.c_str(),nPtBins_,jetPtBinning_,nPtBins_,jetPtBinning_);
  
  string btitle = "ExpectedDEb_" + name_;
	string ltitle = "ExpectedDEl_" + name_;
	histo1D_[btitle] = new TH1F(btitle.c_str(),"Expected energy correction for bjets;#Delta E_{b};#events",80,-1,1);
  histo1D_[ltitle] = new TH1F(ltitle.c_str(),"Expected energy correction for light jets;#Delta E_{l};#events",80,-1,1);
  string ptBjetTitle = "PtBjet_" + name_;
  string ptLightjetTitle = "PtLightjet_" + name_;
  for(unsigned int nBin=0; nBin<nPtBins_; nBin++)
  {
    stringstream ss1, ss2;
    ss1 << jetPtBinning_[nBin];
    ss2 << jetPtBinning_[nBin+1];
    string tmpBtitle = btitle + "_PtB" + ss1.str() + "-" + ss2.str();
    string tmpLtitle = ltitle + "_PtB" + ss1.str() + "-" + ss2.str();
    string tmpPtBjetTitle = ptBjetTitle + "_PtB" + ss1.str() + "-" + ss2.str();
    string tmpBtitleLight = btitle + "_PtLight" + ss1.str() + "-" + ss2.str();
    string tmpLtitleLight = ltitle + "_PtLight" + ss1.str() + "-" + ss2.str();
    string tmpPtLightjetTitle = ptLightjetTitle + "_PtLight" + ss1.str() + "-" + ss2.str();
  	histo1D_[tmpBtitle] = new TH1F(tmpBtitle.c_str(),"Expected energy correction for bjets;#Delta E_{b};#events",80,-1,1);
    histo1D_[tmpLtitle] = new TH1F(tmpLtitle.c_str(),"Expected energy correction for light jets;#Delta E_{l};#events",80,-1,1);
    histo1D_[tmpPtBjetTitle] = new TH1F(tmpPtBjetTitle.c_str(),"Pt of b-jets",1000000,jetPtBinning_[nBin]-1,jetPtBinning_[nBin+1]+1);
    histo1D_[tmpBtitleLight] = new TH1F(tmpBtitleLight.c_str(),"Expected energy correction for bjets;#Delta E_{b};#events",80,-1,1);
    histo1D_[tmpLtitleLight] = new TH1F(tmpLtitleLight.c_str(),"Expected energy correction for light jets;#Delta E_{l};#events",80,-1,1);
    histo1D_[tmpPtLightjetTitle] = new TH1F(tmpPtLightjetTitle.c_str(),"Pt of light-jets",1000000,jetPtBinning_[nBin]-1,jetPtBinning_[nBin+1]+1);
    for(unsigned int nBin2=0; nBin2<nPtBins_; nBin2++)
    {
      stringstream ss3, ss4;
      ss3 << jetPtBinning_[nBin2];
      ss4 << jetPtBinning_[nBin2+1];
      string tmpLtitle2D = ltitle + "_PtB" + ss1.str() + "-" + ss2.str() + "_PtLight" + ss3.str() + "-" + ss4.str();
      string tmpBtitle2D = btitle + "_PtB" + ss1.str() + "-" + ss2.str() + "_PtLight" + ss3.str() + "-" + ss4.str();
      histo1D_[tmpLtitle2D] = new TH1F(tmpLtitle2D.c_str(),"Expected energy correction for light jets;#Delta E_{l};#events",40,-1,1);
      histo1D_[tmpBtitle2D] = new TH1F(tmpBtitle2D.c_str(),"Expected energy correction for bjets;#Delta E_{b};#events",40,-1,1);
    }
  }
  
  string btitlePtbins = btitle + "_PtBjetBins";
  string ltitlePtbins = ltitle + "_PtBjetBins";
  histo1D_[btitlePtbins] = new TH1F(btitlePtbins.c_str(),btitlePtbins.c_str(),nPtBins_,jetPtBinning_);
  histo1D_[ltitlePtbins] = new TH1F(ltitlePtbins.c_str(),ltitlePtbins.c_str(),nPtBins_,jetPtBinning_);
}

ExpCorrCalculator::~ExpCorrCalculator()
{
  if(diffFunction_) delete diffFunction_;
  if(dEbFunction_) delete dEbFunction_;
}

void ExpCorrCalculator::FillLightJet(float DEl, float PtLight, float PtB)
{
	string ltitle = "ExpectedDEl_" + name_;
	string ptLightjetTitle = "PtLightjet_" + name_;
  histo1D_[ltitle]->Fill(DEl);
  for(unsigned int nBin=0; nBin<nPtBins_; nBin++)
  {
    if(PtB >= jetPtBinning_[nBin] && PtB < jetPtBinning_[nBin+1])
    {
      stringstream ss1, ss2;
      ss1 << jetPtBinning_[nBin];
      ss2 << jetPtBinning_[nBin+1];
      string tmpLtitle = ltitle + "_PtB" + ss1.str() + "-" + ss2.str();
      histo1D_[tmpLtitle]->Fill(DEl);
    }
    if(PtLight >= jetPtBinning_[nBin] && PtLight < jetPtBinning_[nBin+1])
    {
      stringstream ss1, ss2;
      ss1 << jetPtBinning_[nBin];
      ss2 << jetPtBinning_[nBin+1];
      string tmpLtitleLight = ltitle + "_PtLight" + ss1.str() + "-" + ss2.str();
      string tmpPtLightjetTitle = ptLightjetTitle + "_PtLight" + ss1.str() + "-" + ss2.str();
      histo1D_[tmpLtitleLight]->Fill(DEl);
      histo1D_[tmpPtLightjetTitle]->Fill(PtLight);
      for(unsigned int nBin2=0; nBin2<nPtBins_; nBin2++)
      {
        if(PtB >= jetPtBinning_[nBin2] && PtB < jetPtBinning_[nBin2+1])
        {
          stringstream ss3, ss4;
          ss3 << jetPtBinning_[nBin2];
          ss4 << jetPtBinning_[nBin2+1];
          string tmpLtitle2D = ltitle + "_PtB" + ss1.str() + "-" + ss2.str() + "_PtLight" + ss3.str() + "-" + ss4.str();
          histo1D_[tmpLtitle2D]->Fill(DEl);
        }
      }
    }
  }
}

void ExpCorrCalculator::FillBJet(float DEb, float PtLight, float PtB)
{
	string btitle = "ExpectedDEb_" + name_;
	string ptBjetTitle = "PtBjet_" + name_;
  histo1D_[btitle]->Fill(DEb);
  
  string ptLightvsptBTitle = "PtLightVSPtB_" + name_;
  histo2D_[ptLightvsptBTitle]->Fill(PtLight,PtB);
  for(unsigned int nBin=0; nBin<nPtBins_; nBin++)
  {
    if(PtB >= jetPtBinning_[nBin] && PtB < jetPtBinning_[nBin+1])
    {
      stringstream ss1, ss2;
      ss1 << jetPtBinning_[nBin];
      ss2 << jetPtBinning_[nBin+1];
      string tmpBtitle = btitle + "_PtB" + ss1.str() + "-" + ss2.str();
      string tmpPtBjetTitle = ptBjetTitle + "_PtB" + ss1.str() + "-" + ss2.str();
      histo1D_[tmpBtitle]->Fill(DEb);
      histo1D_[tmpPtBjetTitle]->Fill(PtB);
    }
    if(PtLight >= jetPtBinning_[nBin] && PtLight < jetPtBinning_[nBin+1])
    {
      stringstream ss1, ss2;
      ss1 << jetPtBinning_[nBin];
      ss2 << jetPtBinning_[nBin+1];
      string tmpBtitleLight = btitle + "_PtLight" + ss1.str() + "-" + ss2.str();
      histo1D_[tmpBtitleLight]->Fill(DEb);
      for(unsigned int nBin2=0; nBin2<nPtBins_; nBin2++)
      {
        if(PtB >= jetPtBinning_[nBin2] && PtB < jetPtBinning_[nBin2+1])
        {
          stringstream ss3, ss4;
          ss3 << jetPtBinning_[nBin2];
          ss4 << jetPtBinning_[nBin2+1];
          string tmpBtitle2D = btitle + "_PtB" + ss1.str() + "-" + ss2.str() + "_PtLight" + ss3.str() + "-" + ss4.str();
          histo1D_[tmpBtitle2D]->Fill(DEb);
        }
      }
    }
  }
}

void ExpCorrCalculator::Fill3Jets(vector<TLorentzVector> TLV)
{
  jetTLVs_.push_back(TLV);
}

void ExpCorrCalculator::Fill3MCParticles(vector<TLorentzVector> TLV)
{
  mcParticlesTLVs_.push_back(TLV);
}

void ExpCorrCalculator::FillMuons(vector<TLorentzVector> TLV)
{
  muonTLVs_.push_back(TLV);
}

void ExpCorrCalculator::FillElectrons(vector<TLorentzVector> TLV)
{
  electronTLVs_.push_back(TLV);
}

void ExpCorrCalculator::FillMaxMVA(float maxMVA)
{
  maxMVAs_.push_back(maxMVA);
}

pair<float,float> ExpCorrCalculator::GetFittedMass(TH1F* histo)
{
  pair<float,float> values;
  string func_title = string(histo->GetName())+"_Fitted";
  
  TF1* fit_func = new TF1(func_title.c_str(),"gaus");
  double left = histo->GetMean() - histo->GetRMS();
  double right = histo->GetMean() + histo->GetRMS();
  fit_func->SetRange(left,right);
  histo->Fit(fit_func,"RQ");
  values.first = fit_func->GetParameter(1);
  values.second = fit_func->GetParError(1);
  
  return values;
}

void ExpCorrCalculator::Calculate(bool printToScreen)
{
  string ltitle = "ExpectedDEl_" + name_;
  string btitle = "ExpectedDEb_" + name_;
  string btitlePtbins = btitle + "_PtBjetBins";
  string ltitlePtbins = ltitle + "_PtBjetBins";
  string lMethod = "doubleGaus";
  string bMethod = "mean";
  pair<float,pair<float,float> > inclLCorr = GetExpectation(histo1D_[ltitle], lMethod);
  pair<float,pair<float,float> > inclBCorr = GetExpectation(histo1D_[btitle], bMethod);

//  inclLCorr.first = -9.47001/100; // Estimated values of JES analysis
//  inclBCorr.first = -0.9923/100; // Estimated values of JES analysis
//  inclBCorr.first = 0.027532; // Expected from Mass Extrapolation (see Stijn's notes)
//  inclLCorr.first = -0.103231; // Expected from Mass Extrapolation (see Stijn's notes)

  if(printToScreen)
  {
    cout << "ExpCorrCalculator Results for " << name_ << endl;
    cout << ltitle << " = " << inclLCorr.first*100 << " +- " << inclLCorr.second.first*100 << " % (stat)" << " +- " << inclLCorr.second.second*100 << " % (syst)" << endl;
    cout << btitle << " = " << inclBCorr.first*100 << " +- " << inclBCorr.second.first*100 << " % (stat)" << " +- " << inclBCorr.second.second*100 << " % (syst)" << endl;
  }
  
  string ltitle2D = ltitle + "_2DPtBins";
  string btitle2D = btitle + "_2DPtBins";
  string difftitle2D = "ExpectedDEdiff_" + name_ + "_2DPtBins";
  histo2D_[ltitle2D] = new TH2F(ltitle2D.c_str(),ltitle2D.c_str(),nPtBins_,jetPtBinning_,nPtBins_,jetPtBinning_);
  histo2D_[ltitle2D]->GetXaxis()->SetTitle("Pt Light Jet");
  histo2D_[ltitle2D]->GetXaxis()->SetLimits(0,250);
  histo2D_[ltitle2D]->GetYaxis()->SetTitle("Pt B Jet");
  histo2D_[ltitle2D]->GetYaxis()->SetRangeUser(0,250);
  
  histo2D_[btitle2D] = new TH2F(btitle2D.c_str(),btitle2D.c_str(),nPtBins_,jetPtBinning_,nPtBins_,jetPtBinning_);
  histo2D_[btitle2D]->GetXaxis()->SetTitle("Pt Light Jet");
  histo2D_[btitle2D]->GetXaxis()->SetLimits(0,250);
  histo2D_[btitle2D]->GetYaxis()->SetTitle("Pt B Jet");
  histo2D_[btitle2D]->GetYaxis()->SetRangeUser(0,250);
  
  histo2D_[difftitle2D] = new TH2F(difftitle2D.c_str(),difftitle2D.c_str(),nPtBins_,jetPtBinning_,nPtBins_,jetPtBinning_);
  histo2D_[difftitle2D]->GetXaxis()->SetTitle("Pt Light Jet");
  histo2D_[difftitle2D]->GetXaxis()->SetLimits(0,250);
  histo2D_[difftitle2D]->GetYaxis()->SetTitle("Pt B Jet");
  histo2D_[difftitle2D]->GetYaxis()->SetRangeUser(0,250);
  
  string ptBjetTitle = "PtBjet_" + name_;
  string ptLightjetTitle = "PtLightjet_" + name_;
  Float_t PtBJets[nPtBins_], PtBJetsErr[nPtBins_], DElPtB[nPtBins_], ErrDElPtB[nPtBins_], DEbPtB[nPtBins_], ErrDEbPtB[nPtBins_],
    DEdiffPtB[nPtBins_], ErrDEdiffPtB[nPtBins_], PtLightJets[nPtBins_], PtLightJetsErr[nPtBins_], DElPtLight[nPtBins_], ErrDElPtLight[nPtBins_], 
    DEbPtLight[nPtBins_], ErrDEbPtLight[nPtBins_], DEdiffPtLight[nPtBins_], ErrDEdiffPtLight[nPtBins_];
  for(unsigned int nBin=0; nBin<nPtBins_; nBin++)
  {
    stringstream ss1, ss2;
    ss1 << jetPtBinning_[nBin];
    ss2 << jetPtBinning_[nBin+1];
    string tmpLtitle = ltitle + "_PtB" + ss1.str() + "-" + ss2.str();
    string tmpBtitle = btitle + "_PtB" + ss1.str() + "-" + ss2.str();
    string tmpPtBjetTitle = ptBjetTitle + "_PtB" + ss1.str() + "-" + ss2.str();
    string tmpBtitleLight = btitle + "_PtLight" + ss1.str() + "-" + ss2.str();
    string tmpLtitleLight = ltitle + "_PtLight" + ss1.str() + "-" + ss2.str();
    string tmpPtLightjetTitle = ptLightjetTitle + "_PtLight" + ss1.str() + "-" + ss2.str();
    
    PtBJets[nBin] = histo1D_[tmpPtBjetTitle]->GetMean();
    PtBJetsErr[nBin] = histo1D_[tmpPtBjetTitle]->GetMeanError();
    PtLightJets[nBin] = histo1D_[tmpPtLightjetTitle]->GetMean();
    PtLightJetsErr[nBin] = histo1D_[tmpPtLightjetTitle]->GetMeanError();
    
    pair<float,pair<float,float> > tmpLCorrPtB = GetExpectation(histo1D_[tmpLtitle], lMethod);
    DElPtB[nBin] = tmpLCorrPtB.first;
    ErrDElPtB[nBin] = tmpLCorrPtB.second.first;
//    DElPtB[nBin] = inclLCorr.first; // inclusive corrections!
//    ErrDElPtB[nBin] = inclLCorr.second.first; //inclusive corrections!

    pair<float,pair<float,float> > tmpBCorrPtB = GetExpectation(histo1D_[tmpBtitle], bMethod);
    DEbPtB[nBin] = tmpBCorrPtB.first;
    ErrDEbPtB[nBin] = tmpBCorrPtB.second.first;
//    DEbPtB[nBin] = inclBCorr.first;
//    ErrDEbPtB[nBin] = inclBCorr.second.first;
    
    DEdiffPtB[nBin] = DEbPtB[nBin] - DElPtB[nBin];
    ErrDEdiffPtB[nBin] = sqrt(pow(ErrDEbPtB[nBin],2) + pow(ErrDElPtB[nBin],2));
    
    histo1D_[ltitlePtbins]->Fill(PtBJets[nBin],DElPtB[nBin]);
    histo1D_[btitlePtbins]->Fill(PtBJets[nBin],DEbPtB[nBin]);
    
    pair<float,pair<float,float> > tmpLCorrPtLight = GetExpectation(histo1D_[tmpLtitleLight], lMethod);
    DElPtLight[nBin] = tmpLCorrPtLight.first;
    ErrDElPtLight[nBin] = tmpLCorrPtLight.second.first;
    
    pair<float,pair<float,float> > tmpBCorrPtLight = GetExpectation(histo1D_[tmpBtitleLight], bMethod);
    DEbPtLight[nBin] = tmpBCorrPtLight.first;
    ErrDEbPtLight[nBin] = tmpBCorrPtLight.second.first;
   
    DEdiffPtLight[nBin] = DEbPtLight[nBin] - DElPtLight[nBin];
    ErrDEdiffPtLight[nBin] = sqrt(pow(ErrDEbPtLight[nBin],2) + pow(ErrDElPtLight[nBin],2));
    
    for(unsigned int nBin2=0; nBin2<nPtBins_; nBin2++)
    {
      stringstream ss3, ss4;
      ss3 << jetPtBinning_[nBin2];
      ss4 << jetPtBinning_[nBin2+1];
      string tmpLtitle2D = ltitle + "_PtB" + ss1.str() + "-" + ss2.str() + "_PtLight" + ss3.str() + "-" + ss4.str();
      string tmpBtitle2D = btitle + "_PtB" + ss1.str() + "-" + ss2.str() + "_PtLight" + ss3.str() + "-" + ss4.str();
      
      pair<float,pair<float,float> > tmpLCorr2D = GetExpectation(histo1D_[tmpLtitle2D], lMethod);
      histo2D_[ltitle2D]->Fill((jetPtBinning_[nBin]+jetPtBinning_[nBin+1])/2,(jetPtBinning_[nBin2]+jetPtBinning_[nBin2+1])/2,tmpLCorr2D.first);
      
      pair<float,pair<float,float> > tmpBCorr2D = GetExpectation(histo1D_[tmpBtitle2D], bMethod);
      histo2D_[btitle2D]->Fill((jetPtBinning_[nBin]+jetPtBinning_[nBin+1])/2,(jetPtBinning_[nBin2]+jetPtBinning_[nBin2+1])/2,tmpBCorr2D.first);
      
      histo2D_[difftitle2D]->Fill((jetPtBinning_[nBin]+jetPtBinning_[nBin+1])/2,(jetPtBinning_[nBin2]+jetPtBinning_[nBin2+1])/2,tmpBCorr2D.first-tmpLCorr2D.first);
    }
  }
  
  string lGraphTitle = ltitle + "_PtBJetBins";
  string bGraphTitle = btitle + "_PtBJetBins";
  string diffGraphTitle = "ExpectedDEdiff_" + name_ + "_PtBJetBins";
  graphErr_[lGraphTitle] = new TGraphErrors(nPtBins_,PtBJets,DElPtB,PtBJetsErr,ErrDElPtB);
  graphErr_[lGraphTitle]->SetNameTitle(lGraphTitle.c_str(),lGraphTitle.c_str());
  graphErr_[lGraphTitle]->GetXaxis()->SetTitle("Pt B Jet");
  graphErr_[lGraphTitle]->GetYaxis()->SetTitle("#DeltaEl (%)");
  
  graphErr_[bGraphTitle] = new TGraphErrors(nPtBins_,PtBJets,DEbPtB,PtBJetsErr,ErrDEbPtB);
  graphErr_[bGraphTitle]->SetNameTitle(bGraphTitle.c_str(),bGraphTitle.c_str());
  graphErr_[bGraphTitle]->GetXaxis()->SetTitle("Pt B Jet");
  graphErr_[bGraphTitle]->GetYaxis()->SetTitle("#DeltaEb (%)");
  
  dEbFunction_ = new TF1((bGraphTitle+"_fit").c_str(),"[0]+[1]/x+[2]/(x*x)",10,1800.);
  graphErr_[bGraphTitle]->Fit((bGraphTitle+"_fit").c_str(),"MQR");
  
  graphErr_[diffGraphTitle] = new TGraphErrors(nPtBins_,PtBJets,DEdiffPtB,PtBJetsErr,ErrDEdiffPtB);
  graphErr_[diffGraphTitle]->SetNameTitle(diffGraphTitle.c_str(),diffGraphTitle.c_str());
  graphErr_[diffGraphTitle]->GetXaxis()->SetTitle("Pt B Jet");
  graphErr_[diffGraphTitle]->GetYaxis()->SetTitle("#DeltaEb - #DeltaEl (%)");
  
  diffFunction_ = new TF1((diffGraphTitle+"_fit").c_str(),"[0]+[1]/x+[2]/(x*x)",10,1800.); //"sqrt([0]*[0]+[1]*[1]/(x*x))",10,1800.);
  graphErr_[diffGraphTitle]->Fit((diffGraphTitle+"_fit").c_str(),"MQR");
  
  string lGraphTitlePtLight = ltitle + "_PtLightJetBins";
  string bGraphTitlePtLight = btitle + "_PtLightJetBins";
  string diffGraphTitlePtLight = "ExpectedDEdiff_" + name_ + "_PtLightJetBins";
  graphErr_[lGraphTitlePtLight] = new TGraphErrors(nPtBins_,PtLightJets,DElPtLight,PtLightJetsErr,ErrDElPtLight);
  graphErr_[lGraphTitlePtLight]->SetNameTitle(lGraphTitlePtLight.c_str(),lGraphTitlePtLight.c_str());
  graphErr_[lGraphTitlePtLight]->GetXaxis()->SetTitle("Pt Light Jet");
  graphErr_[lGraphTitlePtLight]->GetYaxis()->SetTitle("#DeltaEl (%)");
  
  graphErr_[bGraphTitlePtLight] = new TGraphErrors(nPtBins_,PtLightJets,DEbPtLight,PtLightJetsErr,ErrDEbPtLight);
  graphErr_[bGraphTitlePtLight]->SetNameTitle(bGraphTitlePtLight.c_str(),bGraphTitlePtLight.c_str());
  graphErr_[bGraphTitlePtLight]->GetXaxis()->SetTitle("Pt Light Jet");
  graphErr_[bGraphTitlePtLight]->GetYaxis()->SetTitle("#DeltaEb (%)");
  
  graphErr_[diffGraphTitlePtLight] = new TGraphErrors(nPtBins_,PtLightJets,DEdiffPtLight,PtLightJetsErr,ErrDEdiffPtLight);
  graphErr_[diffGraphTitlePtLight]->SetNameTitle(diffGraphTitlePtLight.c_str(),diffGraphTitlePtLight.c_str());
  graphErr_[diffGraphTitlePtLight]->GetXaxis()->SetTitle("Pt Light Jet");
  graphErr_[diffGraphTitlePtLight]->GetYaxis()->SetTitle("#DeltaEb - #DeltaEl (%)");
	
	// Top Mass stuff
	if( mcParticlesTLVs_.size() != jetTLVs_.size() )
	  cout << "mcParticlesTLVs_.size() != jetTLVs_.size()" << endl;
	
	string DEbTitleDEdiffCorr = "ExpectedDEb_DEdiffCorrected_" + name_;
	string DEallTitleDEdiffCorr = "ExpectedDEall_DEdiffCorrected_" + name_;
	histo1D_[DEbTitleDEdiffCorr] = new TH1F(DEbTitleDEdiffCorr.c_str(),DEbTitleDEdiffCorr.c_str(),80,-1,1);
	histo1D_[DEallTitleDEdiffCorr] = new TH1F(DEallTitleDEdiffCorr.c_str(),DEallTitleDEdiffCorr.c_str(),80,-1,1);
	
  string unCorrTitleW = "unCorr_" + name_ + "_mW";
	string corrTitleW = "Corrected_" + name_ + "_mW";
	string corrPtBinsTitleW = "Corrected_PtBins_" + name_ + "_mW";
	string corr2DPtBinsTitleW = "Corrected_2DPtBins_" + name_ + "_mW";
	string corrDEdiffTitleW = "Corrected_DEdiff_" + name_ + "_mW";
	string unCorrTitleTop = "unCorr_" + name_ + "_mTop";
	string corrTitleTop = "Corrected_" + name_ + "_mTop";
	string corrPtBinsTitleTop = "Corrected_PtBins_" + name_ + "_mTop";
	string corr2DPtBinsTitleTop = "Corrected_2DPtBins_" + name_ + "_mTop";
	string corrDEdiffTitleTop = "Corrected_DEdiff_" + name_ + "_mTop";
	histo1D_[unCorrTitleW] = new TH1F(unCorrTitleW.c_str(),unCorrTitleW.c_str(),50,0,150);
	histo1D_[corrTitleW] = new TH1F(corrTitleW.c_str(),corrTitleW.c_str(),50,0,150);
	histo1D_[corrPtBinsTitleW] = new TH1F(corrPtBinsTitleW.c_str(),corrPtBinsTitleW.c_str(),50,0,150);
	histo1D_[corr2DPtBinsTitleW] = new TH1F(corr2DPtBinsTitleW.c_str(),corr2DPtBinsTitleW.c_str(),50,0,150);
	histo1D_[corrDEdiffTitleW] = new TH1F(corrDEdiffTitleW.c_str(),corrDEdiffTitleW.c_str(),50,0,150);
	histo1D_[unCorrTitleTop] = new TH1F(unCorrTitleTop.c_str(),unCorrTitleTop.c_str(),50,100,250);
	histo1D_[corrTitleTop] = new TH1F(corrTitleTop.c_str(),corrTitleTop.c_str(),50,100,250);
	histo1D_[corrPtBinsTitleTop] = new TH1F(corrPtBinsTitleTop.c_str(),corrPtBinsTitleTop.c_str(),50,100,250);
	histo1D_[corr2DPtBinsTitleTop] = new TH1F(corr2DPtBinsTitleTop.c_str(),corr2DPtBinsTitleTop.c_str(),50,100,250);
	histo1D_[corrDEdiffTitleTop] = new TH1F(corrDEdiffTitleTop.c_str(),corrDEdiffTitleTop.c_str(),50,100,250);

	string angleQQreco = "SpaceAngleQQ_RECO_"+name_;
	string angleWBreco = "SpaceAngleWB_RECO_"+name_;
	string angleQQgen = "SpaceAngleQQ_GEN_"+name_;
	string angleWBgen = "SpaceAngleWB_GEN_"+name_;
	histo1D_[angleQQreco] = new TH1F(angleQQreco.c_str(),angleQQreco.c_str(),64,0,3.2);
	histo1D_[angleQQgen] = new TH1F(angleQQgen.c_str(),angleQQgen.c_str(),64,0,3.2);
	histo1D_[angleWBreco] = new TH1F(angleWBreco.c_str(),angleWBreco.c_str(),64,0,3.2);
	histo1D_[angleWBgen] = new TH1F(angleWBgen.c_str(),angleWBgen.c_str(),64,0,3.2);
	
	string angleQQdiff = "SpaceAngleQQ_diff_"+name_;
	string angleWBdiff = "SpaceAngleWB_diff_"+name_;
	histo1D_[angleQQdiff] = new TH1F(angleQQdiff.c_str(),angleQQdiff.c_str(),50,-1,1);
	histo1D_[angleWBdiff] = new TH1F(angleWBdiff.c_str(),angleWBdiff.c_str(),50,-1,1);
	
	string angleQQdiffVSdRQQ = "SpaceAngleQQ_diffVSdRQQ_RECO_"+name_;
	string angleWBdiffVSdRWB = "SpaceAngleWB_diffVSdRWB_RECO_"+name_;
	string angleQQdiffVSangleQQ = "SpaceAngleQQ_diffVSAngleQQ_RECO_"+name_;
	string angleWBdiffVSangleWB = "SpaceAngleWB_diffVSAngleWB_RECO_"+name_;
	string angleQQdiffVSdEtaQQ = "SpaceAngleQQ_diffVSdEtaQQ_RECO_"+name_;
	string angleWBdiffVSdEtaWB = "SpaceAngleWB_diffVSdEtaWB_RECO_"+name_;
	
	histo2D_[angleQQdiffVSdRQQ] = new TH2F(angleQQdiffVSdRQQ.c_str(),angleQQdiffVSdRQQ.c_str(),32,0,6.4,200,-1,1);
	histo2D_[angleWBdiffVSdRWB] = new TH2F(angleWBdiffVSdRWB.c_str(),angleWBdiffVSdRWB.c_str(),32,0,6.4,200,-1,1);
	histo2D_[angleQQdiffVSangleQQ] = new TH2F(angleQQdiffVSangleQQ.c_str(),angleQQdiffVSangleQQ.c_str(),32,0,6.4,200,-1,1);
	histo2D_[angleWBdiffVSangleWB] = new TH2F(angleWBdiffVSangleWB.c_str(),angleWBdiffVSangleWB.c_str(),32,0,6.4,200,-1,1);
	histo2D_[angleQQdiffVSdEtaQQ] = new TH2F(angleQQdiffVSdEtaQQ.c_str(),angleQQdiffVSdEtaQQ.c_str(),32,0,6.4,200,-1,1);
	histo2D_[angleWBdiffVSdEtaWB] = new TH2F(angleWBdiffVSdEtaWB.c_str(),angleWBdiffVSdEtaWB.c_str(),32,0,6.4,200,-1,1);
	
  Float_t AngleBinning[17] = {0.0,0.4,0.8,1.2,1.6,2.0,2.4,2.8,3.2,3.6,4.0,4.4,4.8,5.2,5.6,6.0,6.4};
  Float_t AngleDiffBinning[17] = {-0.5,-0.4,-0.3,-0.2,-0.15,-0.1,-0.05,-0.025,0.,0.025,0.05,0.1,0.15,0.2,0.3,0.4,0.5};

  string corrEdiffL = "corrETdiffL_" + name_;
  string corrEdiffB = "corrETdiffB_" + name_;
  
  histo1D_[corrEdiffL] = new TH1F(corrEdiffL.c_str(),corrEdiffL.c_str(),80,-1,1);
  histo1D_[corrEdiffB] = new TH1F(corrEdiffB.c_str(),corrEdiffB.c_str(),80,-1,1);
  
  string corr2DETdiffL = "corr2DETdiffL_" + name_;
  string corr2DETdiffB = "corr2DETdiffB_" + name_;
  
  histo1D_[corr2DETdiffL] = new TH1F(corr2DETdiffL.c_str(),corr2DETdiffL.c_str(),80,-1,1);
  histo1D_[corr2DETdiffB] = new TH1F(corr2DETdiffB.c_str(),corr2DETdiffB.c_str(),80,-1,1);
  
  string PtBVSrelEdiffB = "PtBVSrelEdiffB_" + name_;
  histo2D_[PtBVSrelEdiffB] = new TH2F(PtBVSrelEdiffB.c_str(),PtBVSrelEdiffB.c_str(),50,0,200,120,-1,2);
  
  string EtaBVSrelEdiffB = "EtaBVSrelEdiffB_" + name_;
  histo2D_[EtaBVSrelEdiffB] = new TH2F(EtaBVSrelEdiffB.c_str(),EtaBVSrelEdiffB.c_str(),50,-2.5,2.5,120,-1,2);
  
  string angleBquarkjetVSrelEdiffB = "angleBquarkjetVSrelEdiffB_" + name_;
  histo2D_[angleBquarkjetVSrelEdiffB] = new TH2F(angleBquarkjetVSrelEdiffB.c_str(),angleBquarkjetVSrelEdiffB.c_str(),50,0,0.5,120,-1,2);
  
  string MinDRBmuVSrelEdiffB = "MinDRBmuVSrelEdiffB_" + name_;
  histo2D_[MinDRBmuVSrelEdiffB] = new TH2F(MinDRBmuVSrelEdiffB.c_str(),MinDRBmuVSrelEdiffB.c_str(),120,0,6,120,-1,2);
      
  string MinDRBelVSrelEdiffB = "MinDRBelVSrelEdiffB_" + name_;
  histo2D_[MinDRBelVSrelEdiffB] = new TH2F(MinDRBelVSrelEdiffB.c_str(),MinDRBelVSrelEdiffB.c_str(),120,0,6,120,-1,2);
      
  string MinDRBlepVSrelEdiffB = "MinDRBlepVSrelEdiffB_" + name_;
  histo2D_[MinDRBlepVSrelEdiffB] = new TH2F(MinDRBlepVSrelEdiffB.c_str(),MinDRBlepVSrelEdiffB.c_str(),120,0,6,120,-1,2);
  
  string MaxMVAVSrelEdiffB = "MaxMVAVSrelEdiffB_" + name_;
  histo2D_[MaxMVAVSrelEdiffB] = new TH2F(MaxMVAVSrelEdiffB.c_str(),MaxMVAVSrelEdiffB.c_str(),50,0,1,120,-1,2);
  
  string MinDRJetsVSrelEdiffB = "MinDRJetsVSrelEdiffB_" + name_;
  histo2D_[MinDRJetsVSrelEdiffB] = new TH2F(MinDRJetsVSrelEdiffB.c_str(),MinDRJetsVSrelEdiffB.c_str(),60,0,3,120,-1,2);
      
//  float systBcorr = -inclBCorr.first;
  float systBcorr = 0.;
  
  bool histosInitialized = false; // are the histos initialized?
  
	for(unsigned int i=0; i<jetTLVs_.size(); i++)
	{
 	  float lightCorrFactor2D = histo2D_[ltitle2D]->GetBinContent( histo2D_[ltitle2D]->FindBin( ( jetTLVs_[i][0].Pt() + jetTLVs_[i][1].Pt() )/2 , jetTLVs_[i][2].Pt() ) );
	  float bCorrFactor2D = histo2D_[btitle2D]->GetBinContent( histo2D_[btitle2D]->FindBin( ( jetTLVs_[i][0].Pt() + jetTLVs_[i][1].Pt() )/2 , jetTLVs_[i][2].Pt() ) );
	  	  	  
	  float relEdiffB = (mcParticlesTLVs_[i][2].Et() - jetTLVs_[i][2].Et()) / (jetTLVs_[i][2].Et());

    float corrRelEdiffL1 = (mcParticlesTLVs_[i][0].Et() - jetTLVs_[i][0].Et()*(1+inclLCorr.first)) / (jetTLVs_[i][0].Et()*(1+inclLCorr.first));
    float corrRelEdiffL2 = (mcParticlesTLVs_[i][1].Et() - jetTLVs_[i][1].Et()*(1+inclLCorr.first)) / (jetTLVs_[i][1].Et()*(1+inclLCorr.first));
    float corrRelEdiffB = (mcParticlesTLVs_[i][2].Et() - jetTLVs_[i][2].Et()*(1+inclBCorr.first)) / (jetTLVs_[i][2].Et()*(1+inclBCorr.first));

    float corr2DRelEdiffL1 = (mcParticlesTLVs_[i][0].Et() - jetTLVs_[i][0].Et()*(1+lightCorrFactor2D)) / (jetTLVs_[i][0].Et()*(1+lightCorrFactor2D));
    float corr2DRelEdiffL2 = (mcParticlesTLVs_[i][1].Et() - jetTLVs_[i][1].Et()*(1+lightCorrFactor2D)) / (jetTLVs_[i][1].Et()*(1+lightCorrFactor2D));
    float corr2DRelEdiffB = (mcParticlesTLVs_[i][2].Et() - jetTLVs_[i][2].Et()*(1+bCorrFactor2D)) / (jetTLVs_[i][2].Et()*(1+bCorrFactor2D));

//    if( fabs(corrRelEdiffB) > 0.3 ) continue;
//    if( jetTLVs_[i][2].Pt() < 60 ) continue;
    
    histo1D_[corrEdiffL]->Fill(corrRelEdiffL1);
    histo1D_[corrEdiffL]->Fill(corrRelEdiffL2);
    histo1D_[corrEdiffB]->Fill(corrRelEdiffB);
    
    histo1D_[corr2DETdiffL]->Fill(corr2DRelEdiffL1);
    histo1D_[corr2DETdiffL]->Fill(corr2DRelEdiffL2);
    histo1D_[corr2DETdiffB]->Fill(corr2DRelEdiffB);
    
    float minDRMuBJet = 9999.;
    for(unsigned int j=0; j<muonTLVs_[i].size(); j++)
    {
      if(muonTLVs_[i][j].DeltaR(jetTLVs_[i][2]) < minDRMuBJet)
        minDRMuBJet = muonTLVs_[i][j].DeltaR(jetTLVs_[i][2]);
    }
    
    float minDRElBJet = 9999.;
    for(unsigned int j=0; j<electronTLVs_[i].size(); j++)
    {
      if(electronTLVs_[i][j].DeltaR(jetTLVs_[i][2]) < minDRElBJet)
        minDRElBJet = electronTLVs_[i][j].DeltaR(jetTLVs_[i][2]);
    }
    
    float minDRLeptonBJet = minDRMuBJet;
    if(minDRLeptonBJet > minDRElBJet)
      minDRLeptonBJet = minDRElBJet;
    
    float minDRJets = jetTLVs_[i][0].DeltaR(jetTLVs_[i][1]);
    if( jetTLVs_[i][2].DeltaR(jetTLVs_[i][1]) < minDRJets ) minDRJets = jetTLVs_[i][2].DeltaR(jetTLVs_[i][1]);
    if( jetTLVs_[i][2].DeltaR(jetTLVs_[i][0]) < minDRJets ) minDRJets = jetTLVs_[i][2].DeltaR(jetTLVs_[i][0]);
    
    histo2D_[PtBVSrelEdiffB]->Fill(jetTLVs_[i][2].Pt(),relEdiffB);
    histo2D_[EtaBVSrelEdiffB]->Fill(jetTLVs_[i][2].Eta(),relEdiffB);
    histo2D_[angleBquarkjetVSrelEdiffB]->Fill( (jetTLVs_[i][2]).Angle( mcParticlesTLVs_[i][2].Vect() ) ,relEdiffB);
    histo2D_[MinDRBmuVSrelEdiffB]->Fill(minDRMuBJet,relEdiffB);
    histo2D_[MinDRBelVSrelEdiffB]->Fill(minDRElBJet,relEdiffB);
    histo2D_[MinDRBlepVSrelEdiffB]->Fill(minDRLeptonBJet,relEdiffB);
    histo2D_[MaxMVAVSrelEdiffB]->Fill(maxMVAs_[i],relEdiffB);
    histo2D_[MinDRJetsVSrelEdiffB]->Fill(minDRJets,relEdiffB);

	  float mWunCorr = ( (jetTLVs_[i][0]) + (jetTLVs_[i][1]) ).M();
	  float mTopunCorr = ( (jetTLVs_[i][0]) + (jetTLVs_[i][1]) + (jetTLVs_[i][2]) ).M();
	  histo1D_[unCorrTitleW.c_str()]->Fill(mWunCorr);
	  histo1D_[unCorrTitleTop.c_str()]->Fill(mTopunCorr);

	  float mTopCorr = ( (jetTLVs_[i][0])*(1+inclLCorr.first) + (jetTLVs_[i][1])*(1+inclLCorr.first) + (jetTLVs_[i][2])*(1+inclBCorr.first + systBcorr) ).M();
	  float mWCorr = ( (jetTLVs_[i][0])*(1+inclLCorr.first) + (jetTLVs_[i][1])*(1+inclLCorr.first) ).M();
//	  float mTopCorr = ( (jetTLVs_[i][0])*(1+inclLCorr.first) + (jetTLVs_[i][1])*(1+inclLCorr.first) + (jetTLVs_[i][2])*(1 + dEbFunction_->Eval(jetTLVs_[i][2].Pt()) + systBcorr) ).M();
	  histo1D_[corrTitleW]->Fill(mWCorr);
	  histo1D_[corrTitleTop]->Fill(mTopCorr);
    
	  float mWCorr2D = ( (jetTLVs_[i][0])*(1+lightCorrFactor2D) + (jetTLVs_[i][1])*(1+lightCorrFactor2D) ).M();
	  float mTopCorr2D = ( (jetTLVs_[i][0])*(1+lightCorrFactor2D) + (jetTLVs_[i][1])*(1+lightCorrFactor2D) + (jetTLVs_[i][2])*(1+bCorrFactor2D) ).M();
	  histo1D_[corr2DPtBinsTitleW]->Fill(mWCorr2D);
	  histo1D_[corr2DPtBinsTitleTop]->Fill(mTopCorr2D);

	  float spaceAngleQQ_RECO = (jetTLVs_[i][0]).Angle( jetTLVs_[i][1].Vect() );
	  float spaceAngleWB_RECO = (jetTLVs_[i][0] + jetTLVs_[i][1]).Angle( jetTLVs_[i][2].Vect() );
	  float dRQQ_RECO = (jetTLVs_[i][0]).DeltaR( jetTLVs_[i][1] );
	  float dRWB_RECO = (jetTLVs_[i][0] + jetTLVs_[i][1]).DeltaR( jetTLVs_[i][2] );
	  histo1D_[angleQQreco]->Fill(spaceAngleQQ_RECO);
	  histo1D_[angleWBreco]->Fill(spaceAngleWB_RECO);
	  
	  // indices for jetTLVs_ and mcParticlesTLVs_ are the same by construction...
	  float spaceAngleQQ_GEN = (mcParticlesTLVs_[i][0]).Angle( mcParticlesTLVs_[i][1].Vect() );
	  float spaceAngleWB_GEN = (mcParticlesTLVs_[i][0] + mcParticlesTLVs_[i][1]).Angle( mcParticlesTLVs_[i][2].Vect() );
	  histo1D_[angleQQgen]->Fill(spaceAngleQQ_GEN);
	  histo1D_[angleWBgen]->Fill(spaceAngleWB_GEN);
	  histo1D_[angleQQdiff]->Fill(spaceAngleQQ_RECO - spaceAngleQQ_GEN);
	  histo1D_[angleWBdiff]->Fill(spaceAngleWB_RECO - spaceAngleWB_GEN);
	  
	  float spaceAngleQQdiff = spaceAngleQQ_RECO - spaceAngleQQ_GEN;
	  float spaceAngleWBdiff = spaceAngleWB_RECO - spaceAngleWB_GEN;
	  histo2D_[angleQQdiffVSdRQQ]->Fill(dRQQ_RECO, spaceAngleQQdiff);
	  histo2D_[angleWBdiffVSdRWB]->Fill(dRWB_RECO, spaceAngleWBdiff);
	  histo2D_[angleQQdiffVSangleQQ]->Fill(spaceAngleQQ_RECO, spaceAngleQQdiff);
	  histo2D_[angleWBdiffVSangleWB]->Fill(spaceAngleWB_RECO, spaceAngleWBdiff);
	  histo2D_[angleQQdiffVSdEtaQQ]->Fill(fabs(jetTLVs_[i][0].Eta() - jetTLVs_[i][1].Eta()), spaceAngleQQdiff);
	  histo2D_[angleWBdiffVSdEtaWB]->Fill(fabs((jetTLVs_[i][0] + jetTLVs_[i][1]).Eta() - jetTLVs_[i][2].Eta()), spaceAngleWBdiff);
	  
	  for(unsigned int j=0; j<nPtBins_; j++)
	  {
	    stringstream ss1, ss2;
      ss1 << jetPtBinning_[j];
      ss2 << jetPtBinning_[j+1];
      string PtBinTitle = "_Pt" + ss1.str() + "-" + ss2.str();

      string unCorrTitleW = "unCorr_" + name_ + "_mW" + PtBinTitle;
	    string corrTitleW = "Corrected_" + name_ + "_mW" + PtBinTitle;
	    string unCorrTitleTop = "unCorr_" + name_ + "_mTop" + PtBinTitle;
	    string corrTitleTop = "Corrected_" + name_ + "_mTop" + PtBinTitle;
	    
	    if( ! histosInitialized ) // initialize the histos
	    {
  	    histo1D_[unCorrTitleW] = new TH1F(unCorrTitleW.c_str(),unCorrTitleW.c_str(),50,0,150);
	      histo1D_[corrTitleW] = new TH1F(corrTitleW.c_str(),corrTitleW.c_str(),50,0,150);
	      histo1D_[unCorrTitleTop] = new TH1F(unCorrTitleTop.c_str(),unCorrTitleTop.c_str(),50,100,250);
	      histo1D_[corrTitleTop] = new TH1F(corrTitleTop.c_str(),corrTitleTop.c_str(),50,100,250);
	    }
      
      if( jetTLVs_[i][2].Pt() >= jetPtBinning_[j] && jetTLVs_[i][2].Pt() < jetPtBinning_[j+1] )
      {
        float relEdiffL1 = ( mcParticlesTLVs_[i][0].Et() - jetTLVs_[i][0].Et() ) / (jetTLVs_[i][0].Et());
        float relEdiffL2 = ( mcParticlesTLVs_[i][1].Et() - jetTLVs_[i][1].Et() ) / (jetTLVs_[i][1].Et());
        float relEdiffBDEdiffCorr = ( mcParticlesTLVs_[i][2].Et() - (jetTLVs_[i][2].Et())*(1+DEdiffPtB[j]) ) / ( (jetTLVs_[i][2].Et())*(1+DEdiffPtB[j]) );
	  
        histo1D_[DEbTitleDEdiffCorr]->Fill(relEdiffBDEdiffCorr);
       	histo1D_[DEallTitleDEdiffCorr]->Fill(relEdiffL1);
      	histo1D_[DEallTitleDEdiffCorr]->Fill(relEdiffL2);
      	histo1D_[DEallTitleDEdiffCorr]->Fill(relEdiffBDEdiffCorr);

        float mTopCorrPtBin = ( (jetTLVs_[i][0])*(1+DElPtB[j]) + (jetTLVs_[i][1])*(1+DElPtB[j]) + (jetTLVs_[i][2])*(1+DEbPtB[j]) ).M();
	      float mWCorrPtBin = ( (jetTLVs_[i][0])*(1+DElPtB[j]) + (jetTLVs_[i][1])*(1+DElPtB[j]) ).M();
        
        histo1D_[unCorrTitleW]->Fill(mWunCorr);
	      histo1D_[corrTitleW]->Fill(mWCorrPtBin);
	      histo1D_[corrPtBinsTitleW]->Fill(mWCorrPtBin);
	      histo1D_[unCorrTitleTop]->Fill(mTopunCorr);
	      histo1D_[corrTitleTop]->Fill(mTopCorrPtBin);
	      histo1D_[corrPtBinsTitleTop]->Fill(mTopCorrPtBin);
      }
	  }
	  
    for(unsigned int j=0; j<16; j++)
    {
      stringstream ss1, ss2;
	    ss1 << AngleBinning[j];
	    ss2 << AngleBinning[j+1];

      string mWvsdRQQ = "mW_dRQQ" + ss1.str() + "-" + ss2.str() + "_" + name_;
      string mWvsAngleQQ = "mW_AngleQQ" + ss1.str() + "-" + ss2.str() + "_" + name_;
      string mTopvsdRWB = "mTop_dRWB" + ss1.str() + "-" + ss2.str() + "_" + name_;
      string mTopvsAngleWB = "mTop_AngleWB" + ss1.str() + "-" + ss2.str() + "_" + name_;
      
      stringstream ss3, ss4;
      ss3 << AngleDiffBinning[j];
      ss4 << AngleDiffBinning[j+1];
      
      string mWvsAngleQQdiff = "mW_AngleQQdiff" + ss3.str() + "-" + ss4.str() + "_" + name_;
      string mTopvsAngleWBdiff = "mTop_AngleWBdiff" + ss3.str() + "-" + ss4.str() + "_" + name_;
      
      if( ! histosInitialized ) // initialize the histos
      {
        histo1D_[mWvsdRQQ] = new TH1F(mWvsdRQQ.c_str(),mWvsdRQQ.c_str(),50,0,200);
        histo1D_[mWvsAngleQQ] = new TH1F(mWvsAngleQQ.c_str(),mWvsAngleQQ.c_str(),50,0,200);
        histo1D_[mWvsAngleQQdiff] = new TH1F(mWvsAngleQQdiff.c_str(),mWvsAngleQQdiff.c_str(),50,0,200);
        histo1D_[mTopvsdRWB] = new TH1F(mTopvsdRWB.c_str(),mTopvsdRWB.c_str(),100,0,400);
        histo1D_[mTopvsAngleWB] = new TH1F(mTopvsAngleWB.c_str(),mTopvsAngleWB.c_str(),100,0,400);
        histo1D_[mTopvsAngleWBdiff] = new TH1F(mTopvsAngleWBdiff.c_str(),mTopvsAngleWBdiff.c_str(),100,0,400);
        if(j == 15) histosInitialized = true; // after the last loop, they are all initialized
      }
      
      if( dRQQ_RECO >= AngleBinning[j] && dRQQ_RECO < AngleBinning[j+1] )
        histo1D_[mWvsdRQQ]->Fill(mWCorr);
      if( spaceAngleQQ_RECO >= AngleBinning[j] && spaceAngleQQ_RECO < AngleBinning[j+1] )
        histo1D_[mWvsAngleQQ]->Fill(mWCorr);
      if( spaceAngleQQdiff >= AngleDiffBinning[j] && spaceAngleQQdiff < AngleDiffBinning[j+1] )
        histo1D_[mWvsAngleQQdiff]->Fill(mWCorr);
      if( dRWB_RECO >= AngleBinning[j] && dRWB_RECO < AngleBinning[j+1] )
        histo1D_[mTopvsdRWB]->Fill(mTopCorr);
      if( spaceAngleWB_RECO >= AngleBinning[j] && spaceAngleWB_RECO < AngleBinning[j+1] )
        histo1D_[mTopvsAngleWB]->Fill(mTopCorr);
      if( spaceAngleWBdiff >= AngleDiffBinning[j] && spaceAngleWBdiff < AngleDiffBinning[j+1] )
        histo1D_[mTopvsAngleWBdiff]->Fill(mTopCorr);
    }
	}
	
  pair<float,pair<float,float> > inclCorrDEdiff = GetExpectation(histo1D_[DEallTitleDEdiffCorr], "gaus");
  for(unsigned int i=0; i<jetTLVs_.size(); i++)
	{
	  float mTopCorrDEdiff = ( (jetTLVs_[i][0])*(1+inclCorrDEdiff.first) + (jetTLVs_[i][1])*(1+inclCorrDEdiff.first) + (jetTLVs_[i][2])*(1+inclCorrDEdiff.first) ).M();
    float mWCorrDEdiff = ( (jetTLVs_[i][0])*(1+inclCorrDEdiff.first) + (jetTLVs_[i][1])*(1+inclCorrDEdiff.first) ).M();
    histo1D_[corrDEdiffTitleTop]->Fill(mTopCorrDEdiff);
    histo1D_[corrDEdiffTitleW]->Fill(mWCorrDEdiff);
  }
	
	Float_t UnCorrMW[nPtBins_], UnCorrMWErr[nPtBins_], CorrMW[nPtBins_], CorrMWErr[nPtBins_], UnCorrMTop[nPtBins_], UnCorrMTopErr[nPtBins_], CorrMTop[nPtBins_], CorrMTopErr[nPtBins_];
	
  for(unsigned int j=0; j<nPtBins_; j++)
  {
    UnCorrMW[j] = UnCorrMWErr[j] = CorrMW[j] = CorrMWErr[j] = UnCorrMTop[j] = UnCorrMTopErr[j] = CorrMTop[j] = CorrMTopErr[j] = 0;
    
    stringstream ss1, ss2;
    ss1 << jetPtBinning_[j];
    ss2 << jetPtBinning_[j+1];
    string PtBinTitle = "_Pt" + ss1.str() + "-" + ss2.str();

    string unCorrTitleW = "unCorr_" + name_ + "_mW" + PtBinTitle;
    string corrTitleW = "Corrected_" + name_ + "_mW" + PtBinTitle;
    string unCorrTitleTop = "unCorr_" + name_ + "_mTop" + PtBinTitle;
    string corrTitleTop = "Corrected_" + name_ + "_mTop" + PtBinTitle;
    
    if( histo1D_[unCorrTitleW]->GetEntries() > 1 )
    {
      pair<float,float> fitUnCorrMW = GetFittedMass(histo1D_[unCorrTitleW]);
      UnCorrMW[j] = fitUnCorrMW.first;
      UnCorrMWErr[j] = fitUnCorrMW.second;
      
      pair<float,float> fitCorrMW = GetFittedMass(histo1D_[corrTitleW]);
      CorrMW[j] = fitCorrMW.first;
      CorrMWErr[j] = fitCorrMW.second;
      
      pair<float,float> fitUnCorrMTop = GetFittedMass(histo1D_[unCorrTitleTop]);
      UnCorrMTop[j] = fitUnCorrMTop.first;
      UnCorrMTopErr[j] = fitUnCorrMTop.second;
      
      pair<float,float> fitCorrMTop = GetFittedMass(histo1D_[corrTitleTop]);
      CorrMTop[j] = fitCorrMTop.first;
      CorrMTopErr[j] = fitCorrMTop.second;
    }
  }
  
  string unCorrMWvsPtBJet = "unCorrMW_VS_PtBJet_" + name_;
  string corrMWvsPtBJet = "corrMW_VS_PtBJet_" + name_;
  string unCorrMTopvsPtBJet = "unCorrMTop_VS_PtBJet_" + name_;
  string corrMTopvsPtBJet = "corrMTop_VS_PtBJet_" + name_;
  
  graphErr_[unCorrMWvsPtBJet] = new TGraphErrors(32,PtBJets,UnCorrMW,PtBJetsErr,UnCorrMWErr);
  graphErr_[unCorrMWvsPtBJet]->SetNameTitle(unCorrMWvsPtBJet.c_str(),unCorrMWvsPtBJet.c_str());
  graphErr_[unCorrMWvsPtBJet]->GetXaxis()->SetTitle("P_{T} B-jet");
  graphErr_[unCorrMWvsPtBJet]->GetXaxis()->SetLimits(0,250);
  graphErr_[unCorrMWvsPtBJet]->GetYaxis()->SetTitle("M_{W}^{Fit, UnCorr}");
  graphErr_[unCorrMWvsPtBJet]->GetYaxis()->SetRangeUser(85,95);

  graphErr_[corrMWvsPtBJet] = new TGraphErrors(32,PtBJets,CorrMW,PtBJetsErr,CorrMWErr);
  graphErr_[corrMWvsPtBJet]->SetNameTitle(corrMWvsPtBJet.c_str(),corrMWvsPtBJet.c_str());
  graphErr_[corrMWvsPtBJet]->GetXaxis()->SetTitle("P_{T} B-jet");
  graphErr_[corrMWvsPtBJet]->GetXaxis()->SetLimits(0,250);
  graphErr_[corrMWvsPtBJet]->GetYaxis()->SetTitle("M_{W}^{Fit, ExpCorr}");
  graphErr_[corrMWvsPtBJet]->GetYaxis()->SetRangeUser(75,85);

  graphErr_[unCorrMTopvsPtBJet] = new TGraphErrors(32,PtBJets,UnCorrMTop,PtBJetsErr,UnCorrMTopErr);
  graphErr_[unCorrMTopvsPtBJet]->SetNameTitle(unCorrMTopvsPtBJet.c_str(),unCorrMTopvsPtBJet.c_str());
  graphErr_[unCorrMTopvsPtBJet]->GetXaxis()->SetTitle("P_{T} B-jet");
  graphErr_[unCorrMTopvsPtBJet]->GetXaxis()->SetLimits(0,250);
  graphErr_[unCorrMTopvsPtBJet]->GetYaxis()->SetTitle("M_{Top}^{Fit, UnCorr}");
  graphErr_[unCorrMTopvsPtBJet]->GetYaxis()->SetRangeUser(160,195);

  graphErr_[corrMTopvsPtBJet] = new TGraphErrors(32,PtBJets,CorrMTop,PtBJetsErr,CorrMTopErr);
  graphErr_[corrMTopvsPtBJet]->SetNameTitle(corrMTopvsPtBJet.c_str(),corrMTopvsPtBJet.c_str());
  graphErr_[corrMTopvsPtBJet]->GetXaxis()->SetTitle("P_{T} B-jet");
  graphErr_[corrMTopvsPtBJet]->GetXaxis()->SetLimits(0,250);
  graphErr_[corrMTopvsPtBJet]->GetYaxis()->SetTitle("M_{Top}^{Fit, ExpCorr}");
  graphErr_[corrMTopvsPtBJet]->GetYaxis()->SetRangeUser(155,180);

	profile_[(angleQQdiffVSdRQQ+"_profile").c_str()] = histo2D_[angleQQdiffVSdRQQ.c_str()]->ProfileX((angleQQdiffVSdRQQ+"_profile").c_str());
	profile_[(angleQQdiffVSdRQQ+"_profile").c_str()]->GetXaxis()->SetTitle("#DeltaR(Q,Q) RECO");
	profile_[(angleQQdiffVSdRQQ+"_profile").c_str()]->GetYaxis()->SetTitle("#Omega(Q,Q)_{RECO} - #Omega(Q,Q)_{GEN}");
	
	profile_[(angleWBdiffVSdRWB+"_profile").c_str()] = histo2D_[angleWBdiffVSdRWB.c_str()]->ProfileX((angleWBdiffVSdRWB+"_profile").c_str());
	profile_[(angleWBdiffVSdRWB+"_profile").c_str()]->GetXaxis()->SetTitle("#DeltaR(W,B) RECO");
	profile_[(angleWBdiffVSdRWB+"_profile").c_str()]->GetYaxis()->SetTitle("#Omega(W,B)_{RECO} - #Omega(W,B)_{GEN}");
	
	profile_[(angleQQdiffVSangleQQ+"_profile").c_str()] = histo2D_[angleQQdiffVSangleQQ.c_str()]->ProfileX((angleQQdiffVSangleQQ+"_profile").c_str());
	profile_[(angleQQdiffVSangleQQ+"_profile").c_str()]->GetXaxis()->SetTitle("#Omega(Q,Q) RECO");
	profile_[(angleQQdiffVSangleQQ+"_profile").c_str()]->GetYaxis()->SetTitle("#Omega(Q,Q)_{RECO} - #Omega(Q,Q)_{GEN}");
	
	profile_[(angleWBdiffVSangleWB+"_profile").c_str()] = histo2D_[angleWBdiffVSangleWB.c_str()]->ProfileX((angleWBdiffVSangleWB+"_profile").c_str());
	profile_[(angleWBdiffVSangleWB+"_profile").c_str()]->GetXaxis()->SetTitle("#Omega(W,B) RECO");
	profile_[(angleWBdiffVSangleWB+"_profile").c_str()]->GetYaxis()->SetTitle("#Omega(W,B)_{RECO} - #Omega(W,B)_{GEN}");
	
	profile_[(angleQQdiffVSdEtaQQ+"_profile").c_str()] = histo2D_[angleQQdiffVSdEtaQQ.c_str()]->ProfileX((angleQQdiffVSdEtaQQ+"_profile").c_str());
	profile_[(angleQQdiffVSdEtaQQ+"_profile").c_str()]->GetXaxis()->SetTitle("#Delta#eta(Q,Q) RECO");
	profile_[(angleQQdiffVSdEtaQQ+"_profile").c_str()]->GetYaxis()->SetTitle("#Omega(Q,Q)_{RECO} - #Omega(Q,Q)_{GEN}");
	
	profile_[(angleWBdiffVSdEtaWB+"_profile").c_str()] = histo2D_[angleWBdiffVSdEtaWB.c_str()]->ProfileX((angleWBdiffVSdEtaWB+"_profile").c_str());
	profile_[(angleWBdiffVSdEtaWB+"_profile").c_str()]->GetXaxis()->SetTitle("#Delta#eta(W,B) RECO");
	profile_[(angleWBdiffVSdEtaWB+"_profile").c_str()]->GetYaxis()->SetTitle("#Omega(W,B)_{RECO} - #Omega(W,B)_{GEN}");
	
  Float_t angle[16], angleErr[16], mWdR[16], mWdRErr[16], mWAngle[16], mWAngleErr[16], mTopdR[16], mTopdRErr[16], mTopAngle[16], mTopAngleErr[16];
  Float_t angleDiff[16], mWAngleDiff[16], mWAngleDiffErr[16], mTopAngleDiff[16], mTopAngleDiffErr[16];
  
	for(unsigned int j=0; j<16; j++)
  {
    stringstream ss1, ss2;
    ss1 << AngleBinning[j];
    ss2 << AngleBinning[j+1];

    string mWvsdRQQ = "mW_dRQQ" + ss1.str() + "-" + ss2.str() + "_" + name_;
    string mWvsAngleQQ = "mW_AngleQQ" + ss1.str() + "-" + ss2.str() + "_" + name_;
    string mTopvsdRWB = "mTop_dRWB" + ss1.str() + "-" + ss2.str() + "_" + name_;
    string mTopvsAngleWB = "mTop_AngleWB" + ss1.str() + "-" + ss2.str() + "_" + name_;
    
    angle[j] = (AngleBinning[j+1] + AngleBinning[j]) / 2;
    angleErr[j] = mWAngle[j] = mWAngleErr[j] = mWdR[j] = mWdRErr[j] = mTopdR[j] = mTopdRErr[j] = mTopAngle[j] = mTopAngleErr[j] = 0;
    
    stringstream ss3, ss4;
    ss3 << AngleDiffBinning[j];
    ss4 << AngleDiffBinning[j+1];
    
    string mWvsAngleQQdiff = "mW_AngleQQdiff" + ss3.str() + "-" + ss4.str() + "_" + name_;
    string mTopvsAngleWBdiff = "mTop_AngleWBdiff" + ss3.str() + "-" + ss4.str() + "_" + name_;
    
    angleDiff[j] = (AngleDiffBinning[j+1] + AngleDiffBinning[j]) / 2;
    mWAngleDiff[j] = mWAngleDiffErr[j] = mTopAngleDiff[j] = mTopAngleDiffErr[j] = 0;
    
    if( histo1D_[mWvsdRQQ]->GetEntries() > 1 )
    {
      pair<float,float> fitWdR = GetFittedMass(histo1D_[mWvsdRQQ]);
      mWAngle[j] = fitWdR.first;
      mWAngleErr[j] = fitWdR.second;
    }
    
    if( histo1D_[mWvsAngleQQ]->GetEntries() > 1 )
    {  
      pair<float,float> fitWAngle = GetFittedMass(histo1D_[mWvsAngleQQ]);
      mWdR[j] = fitWAngle.first;
      mWdRErr[j] = fitWAngle.second;
    }
    
    if( histo1D_[mWvsAngleQQdiff]->GetEntries() > 1 )
    { 
      pair<float,float> fitWAngleDiff = GetFittedMass(histo1D_[mWvsAngleQQdiff]);
      mWAngleDiff[j] = fitWAngleDiff.first;
      mWAngleDiffErr[j] = fitWAngleDiff.second;
    }
    
    if( histo1D_[mTopvsdRWB]->GetEntries() > 1 )
    {
      pair<float,float> fitTopdR = GetFittedMass(histo1D_[mTopvsdRWB]);
      mTopdR[j] = fitTopdR.first;
      mTopdRErr[j] = fitTopdR.second;
    }
    
    if( histo1D_[mTopvsAngleWB]->GetEntries() > 1 )
    {
      pair<float,float> fitTopAngle = GetFittedMass(histo1D_[mTopvsAngleWB]);
      mTopAngle[j] = fitTopAngle.first;
      mTopAngleErr[j] = fitTopAngle.second;
    }

    if( histo1D_[mTopvsAngleWBdiff]->GetEntries() > 1 )
    {
      pair<float,float> fitTopAngleDiff = GetFittedMass(histo1D_[mTopvsAngleWBdiff]);
      mTopAngleDiff[j] = fitTopAngleDiff.first;
      mTopAngleDiffErr[j] = fitTopAngleDiff.second;
    }
  }

  string mWvsdRQQ = "mW_vs_dRQQ_" + name_;
  string mWvsAngleQQ = "mW_vs_AngleQQ_" + name_;
  string mWvsAngleQQdiff = "mW_vs_AngleQQdiff_" + name_;
  string mTopvsdRWB = "mTop_vs_dRWB_" + name_;
  string mTopvsAngleWB = "mTop_vs_AngleWB_" + name_;
  string mTopvsAngleWBdiff = "mTop_vs_AngleWBdiff_" + name_;
  
  graphErr_[mWvsdRQQ] = new TGraphErrors(32,angle,mWdR,angleErr,mWdRErr);
  graphErr_[mWvsdRQQ]->SetNameTitle(mWvsdRQQ.c_str(),mWvsdRQQ.c_str());
  graphErr_[mWvsdRQQ]->GetXaxis()->SetTitle("#DeltaR(Q,Q)");
  graphErr_[mWvsdRQQ]->GetXaxis()->SetLimits(0.,5.);
  graphErr_[mWvsdRQQ]->GetYaxis()->SetTitle("M_{W}^{Fit, ExpCorr}");
  graphErr_[mWvsdRQQ]->GetYaxis()->SetRangeUser(70,95);
  
  graphErr_[mWvsAngleQQ] = new TGraphErrors(32,angle,mWAngle,angleErr,mWAngleErr);
  graphErr_[mWvsAngleQQ]->SetNameTitle(mWvsAngleQQ.c_str(),mWvsAngleQQ.c_str());
  graphErr_[mWvsAngleQQ]->GetXaxis()->SetTitle("#Omega(Q,Q)"); 
  graphErr_[mWvsAngleQQ]->GetXaxis()->SetLimits(0.,5.);
  graphErr_[mWvsAngleQQ]->GetYaxis()->SetTitle("M_{W}^{Fit, ExpCorr}");
  graphErr_[mWvsAngleQQ]->GetYaxis()->SetRangeUser(70,95);
  
  graphErr_[mWvsAngleQQdiff] = new TGraphErrors(32,angleDiff,mWAngleDiff,angleErr,mWAngleDiffErr);
  graphErr_[mWvsAngleQQdiff]->SetNameTitle(mWvsAngleQQdiff.c_str(),mWvsAngleQQdiff.c_str());
  graphErr_[mWvsAngleQQdiff]->GetXaxis()->SetTitle("#Omega(Q,Q)_{RECO} - #Omega(Q,Q)_{GEN}"); 
  graphErr_[mWvsAngleQQdiff]->GetXaxis()->SetLimits(-0.5,0.5);
  graphErr_[mWvsAngleQQdiff]->GetYaxis()->SetTitle("M_{W}^{Fit, ExpCorr}");
  graphErr_[mWvsAngleQQdiff]->GetYaxis()->SetRangeUser(70,95);
  
  graphErr_[mTopvsdRWB] = new TGraphErrors(32,angle,mTopdR,angleErr,mTopdRErr);
  graphErr_[mTopvsdRWB]->SetNameTitle(mTopvsdRWB.c_str(),mTopvsdRWB.c_str());
  graphErr_[mTopvsdRWB]->GetXaxis()->SetTitle("#DeltaR(W,B)"); 
  graphErr_[mTopvsdRWB]->GetXaxis()->SetLimits(0.,5.);
  graphErr_[mTopvsdRWB]->GetYaxis()->SetTitle("M_{Top}^{Fit, ExpCorr}");
  graphErr_[mTopvsdRWB]->GetYaxis()->SetRangeUser(160,180);
  
  graphErr_[mTopvsAngleWB] = new TGraphErrors(32,angle,mTopAngle,angleErr,mTopAngleErr);
  graphErr_[mTopvsAngleWB]->SetNameTitle(mTopvsAngleWB.c_str(),mTopvsAngleWB.c_str());
  graphErr_[mTopvsAngleWB]->GetXaxis()->SetTitle("#Omega(W,B)"); 
  graphErr_[mTopvsAngleWB]->GetXaxis()->SetLimits(0.,5.);
  graphErr_[mTopvsAngleWB]->GetYaxis()->SetTitle("M_{Top}^{Fit, ExpCorr}");
  graphErr_[mTopvsAngleWB]->GetYaxis()->SetRangeUser(160,180);
  
  graphErr_[mTopvsAngleWBdiff] = new TGraphErrors(32,angleDiff,mTopAngleDiff,angleErr,mTopAngleDiffErr);
  graphErr_[mTopvsAngleWBdiff]->SetNameTitle(mTopvsAngleWBdiff.c_str(),mTopvsAngleWBdiff.c_str());
  graphErr_[mTopvsAngleWBdiff]->GetXaxis()->SetTitle("#Omega(W,B)_{RECO} - #Omega(W,B)_{GEN}"); 
  graphErr_[mTopvsAngleWBdiff]->GetXaxis()->SetLimits(-0.5,0.5);
  graphErr_[mTopvsAngleWBdiff]->GetYaxis()->SetTitle("M_{Top}^{Fit, ExpCorr}");
  graphErr_[mTopvsAngleWBdiff]->GetYaxis()->SetRangeUser(160,180);
  
  GetExpectation(histo1D_[corrEdiffL], lMethod);
  GetExpectation(histo1D_[corrEdiffB], bMethod);
  GetExpectation(histo1D_[corr2DETdiffL], lMethod);
  GetExpectation(histo1D_[corr2DETdiffB], bMethod);
  
  pair<float,float> fitUnCorrW = GetFittedMass(histo1D_[unCorrTitleW]);
  pair<float,float> fitCorrW = GetFittedMass(histo1D_[corrTitleW]);
  pair<float,float> fitCorrPtBinsW = GetFittedMass(histo1D_[corrPtBinsTitleW]);
  pair<float,float> fitCorr2DPtBinsW = GetFittedMass(histo1D_[corr2DPtBinsTitleW]);
  pair<float,float> fitCorrDEdiffW = GetFittedMass(histo1D_[corrDEdiffTitleW]);
  pair<float,float> fitUnCorrTop = GetFittedMass(histo1D_[unCorrTitleTop]);
  pair<float,float> fitCorrTop = GetFittedMass(histo1D_[corrTitleTop]);
  pair<float,float> fitCorrPtBinsTop = GetFittedMass(histo1D_[corrPtBinsTitleTop]);
  pair<float,float> fitCorr2DPtBinsTop = GetFittedMass(histo1D_[corr2DPtBinsTitleTop]);
  pair<float,float> fitCorrDEdiffTop = GetFittedMass(histo1D_[corrDEdiffTitleTop]);
	
  if(printToScreen)
  {
    cout << "unCorr mW = " << fitUnCorrW.first << " +- " << fitUnCorrW.second << " GeV" << endl;
    cout << "Corr mW = " << fitCorrW.first << " +- " << fitCorrW.second << " GeV" << endl;
    cout << "Corr PtBins mW = " << fitCorrPtBinsW.first << " +- " << fitCorrPtBinsW.second << " GeV" << endl;
    cout << "Corr 2DPtBins mW = " << fitCorr2DPtBinsW.first << " +- " << fitCorr2DPtBinsW.second << " GeV" <<  endl;
    cout << "Corr DEdiff mW = " << fitCorrDEdiffW.first << " +- " << fitCorrDEdiffW.second << " GeV" << endl;
    cout << "unCorr mTop = " << fitUnCorrTop.first << " +- " << fitUnCorrTop.second << " GeV" << endl;
    cout << "Corr mTop = " << fitCorrTop.first << " +- " << fitCorrTop.second << " GeV" << endl;
    cout << "Corr PtBins mTop = " << fitCorrPtBinsTop.first << " +- " << fitCorrPtBinsTop.second << " GeV" << endl;
    cout << "Corr 2DPtBins mTop = " << fitCorr2DPtBinsTop.first << " +- " << fitCorr2DPtBinsTop.second << " GeV" << endl;
    cout << "Corr DEdiff mTop = " << fitCorrDEdiffTop.first << " +- " << fitCorrDEdiffTop.second << " GeV" << endl;
  }

	wasCalculated_ = true;
}

void ExpCorrCalculator::Write(TFile* fout, bool savePNG, string pathPNG)
{
  if(jetTLVs_.size() < 2)
  {
    cout << "ExpCorrCalculator " << name_ << ":  Not enough entries, will not produce plots!" << endl;
    return;
  }
  if( ! wasCalculated_ ) Calculate(true);

  string newPathPNG = pathPNG+name_+"/";
  mkdir(newPathPNG.c_str(),0777);
  
  fout->cd();
  string dirname = "ExpCorrCalculator_" + name_;
	if(fout->Get(dirname.c_str())==0)
		fout->mkdir(dirname.c_str());
	
	fout->cd(dirname.c_str());

	for(std::map<std::string,TH1F*>::const_iterator it = histo1D_.begin(); it != histo1D_.end(); it++)
	{
		TH1F *temp = it->second;
		temp->Write();
		TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
		tempCanvas->SaveAs( (newPathPNG+it->first+".png").c_str() );
	}
  
  for(map<string,TGraphErrors*>::const_iterator it = graphErr_.begin(); it != graphErr_.end(); it++)
  {
    TGraphErrors *temp = it->second;
    temp->Write();
    TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
    tempCanvas->SaveAs( (newPathPNG+it->first+".png").c_str() );
  }
  
  for(map<string,TH2F*>::const_iterator it = histo2D_.begin(); it != histo2D_.end(); it++)
  {
    TH2F *temp = it->second;
    temp->Write();
    TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
    tempCanvas->SaveAs( (newPathPNG+it->first+".png").c_str() );
  }
  
  for(map<string,TProfile*>::const_iterator it = profile_.begin(); it != profile_.end(); it++)
  {
    TProfile *temp = it->second;
    temp->Write();
    TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
    tempCanvas->SaveAs( (newPathPNG+it->first+".png").c_str() );
  }
 
  TFile* tmp = new TFile("ExpectedDEdiff_new.root","UPDATE");
  tmp->cd();
  
  string lGraphTitle = "ExpectedDEl_" + name_ + "_PtBjetBins";
  string bGraphTitle = "ExpectedDEb_" + name_ + "_PtBjetBins";
  string ltitle2D = "ExpectedDEl_" + name_ + "_2DPtBins";
  string btitle2D = "ExpectedDEb_" + name_ + "_2DPtBins";
  
  if( tmp->Get(diffFunction_->GetName()) == 0 || tmp->Get(histo1D_[lGraphTitle]->GetName()) == 0 || tmp->Get(histo1D_[bGraphTitle]->GetName()) == 0
    || tmp->Get(histo2D_[ltitle2D]->GetName()) == 0 || tmp->Get(histo2D_[btitle2D]->GetName()) == 0 )
  {
	  diffFunction_->Write();
	  histo1D_[lGraphTitle]->Write();
	  histo1D_[bGraphTitle]->Write();
	  histo2D_[ltitle2D]->Write();
	  histo2D_[btitle2D]->Write();
	}
	else
	  cout << "Found already items with name: " << name_ << endl;
	
  tmp->Close();
  
 	string angleQQreco = "SpaceAngleQQ_RECO_"+name_;
	string angleWBreco = "SpaceAngleWB_RECO_"+name_;
	string angleQQgen = "SpaceAngleQQ_GEN_"+name_;
	string angleWBgen = "SpaceAngleWB_GEN_"+name_;

 	vector< pair<TH1F*,string> > listOfHistos1;
	listOfHistos1.push_back( pair<TH1F*,string>(histo1D_[angleQQreco], angleQQreco) );
	listOfHistos1.push_back( pair<TH1F*,string>(histo1D_[angleQQgen], angleQQgen) );
	TCanvas* tmpCanv1 = TCanvasCreator(listOfHistos1, "l", ("SpaceAngleQQ_"+name_).c_str() );
	tmpCanv1->SaveAs( (newPathPNG+"SpaceAngleQQ_"+name_+".png").c_str() );
	delete tmpCanv1;
	
 	vector< pair<TH1F*,string> > listOfHistos2;
	listOfHistos2.push_back( pair<TH1F*,string>(histo1D_[angleWBreco], angleWBreco) );
	listOfHistos2.push_back( pair<TH1F*,string>(histo1D_[angleWBgen], angleWBgen) );
	TCanvas* tmpCanv2 = TCanvasCreator(listOfHistos1, "l", ("SpaceAngleWB_"+name_).c_str() );
	tmpCanv2->SaveAs( (newPathPNG+"SpaceAngleWB_"+name_+".png").c_str() );
	delete tmpCanv2;
}

pair<float,float> ExpCorrCalculator::GetExpCorrL()
{
  if( ! wasCalculated_ ) Calculate(false);
  return inclLCorr_;
}

pair<float,float> ExpCorrCalculator::GetExpCorrB()
{
  if( ! wasCalculated_ ) Calculate(false);
  return inclBCorr_;
}


