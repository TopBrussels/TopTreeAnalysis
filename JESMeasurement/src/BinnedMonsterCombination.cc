#include "../interface/BinnedMonsterCombination.h"

BinnedMonsterCombination::BinnedMonsterCombination(string title, bool topMass, vector<float> binning, bool doMvaBinned)
{
  title_ = title;
  topMass_ = topMass;
  doMvaBinned_ = doMvaBinned;
  vector<float> maxMVAbins;
  maxMVAbins.push_back(0);
  maxMVAbins.push_back(0.1);
  maxMVAbins.push_back(0.2);
  maxMVAbins.push_back(0.3);
  maxMVAbins.push_back(0.4);
  maxMVAbins.push_back(0.5);
  maxMVAbins.push_back(0.6);
  maxMVAbins.push_back(0.7);
  maxMVAbins.push_back(0.8);
  maxMVAbins.push_back(0.9);
  maxMVAbins.push_back(1.);
  binning_ = binning;
  
  vector<TH1F*> mWbinnedHistos, mTopbinnedHistos;
  for(unsigned int i=0; i<binning_.size()-1; i++)
  {
    stringstream ss1, ss2;
    ss1 << binning_[i];
    ss2 << binning_[i+1];
    string binTitle = title + "_Bin_" + ss1.str() + "_" + ss2.str();
    monsterCombis_.push_back( new MonsterCombination(binTitle, topMass_) );
    if(doMvaBinned_)
      monsterCombisMVAbinned_.push_back( new MVABinnedMonsterCombination("mvaBinned_"+binTitle, topMass_, maxMVAbins) );
    xBinHistos_.push_back( new TH1F((binTitle+"_xBin").c_str(),(binTitle+"_xBin").c_str(),10000,binning_[i]-0.1,binning[i+1]+0.1) );
    mWbinnedHistos.push_back( new TH1F((binTitle+"_mWbin").c_str(),(binTitle+"_mWbin").c_str(),75,0,150) );
    mTopbinnedHistos.push_back( new TH1F((binTitle+"_mTopbin").c_str(),(binTitle+"_mTopbin").c_str(),75,100,250) );
  }
  binnedVars_["mWbin"] = mWbinnedHistos;
  binnedVars_["mTopbin"] = mTopbinnedHistos;
}

BinnedMonsterCombination::~BinnedMonsterCombination()
{
  
}

void BinnedMonsterCombination::Fill(Monster* monster, float binValue, float weight, int jetCombi)
{
  unsigned int* MVAres = monster->mvaResult(jetCombi);
  for(unsigned int i=0; i<monsterCombis_.size(); i++)
  {
    if(binValue > binning_[i] && binValue <= binning_[i+1])
    {
      monsterCombis_[i]->addMonster(monster, weight, jetCombi);
      monsterCombis_[i]->addExpected(monster, weight, jetCombi);
      if( doMvaBinned_ ) monsterCombisMVAbinned_[i]->Fill(monster, weight);
      xBinHistos_[i]->Fill(binValue);
      
      binnedVars_["mWbin"][i]->Fill( ( monster->selectedJet(MVAres[0]) + monster->selectedJet(MVAres[1]) ).M() );
      binnedVars_["mTopbin"][i]->Fill( ( monster->selectedJet(MVAres[0]) + monster->selectedJet(MVAres[1]) + monster->selectedJet(MVAres[2]) ).M() );
    }
  }
}

void BinnedMonsterCombination::FillExpected(Monster* monster, float binValue, float weights)
{
  for(unsigned int i=0; i<monsterCombis_.size(); i++)
    if(binValue > binning_[i] && binValue <= binning_[i+1])
      monsterCombis_[i]->addExpected(monster, weights);
}

TH2F* BinnedMonsterCombination::AverageMonster(float binValue)
{
  for(unsigned int i=0; i<monsterCombis_.size(); i++)
    if(binValue > binning_[i] && binValue <= binning_[i+1])
      return monsterCombis_[i]->AverageMonster();
  return 0;
}

void BinnedMonsterCombination::SubTractAvMonster(TH2F* avMonster, float binValue, float nBadCombis)
{
  for(unsigned int i=0; i<monsterCombis_.size(); i++)
    if(binValue > binning_[i] && binValue <= binning_[i+1])
      monsterCombis_[i]->SubTractAvMonster(avMonster, nBadCombis);
}

void BinnedMonsterCombination::Write(TFile* f, string pathforPNG, string xTitle, bool silent)
{
  f->cd();
	string dirname = "BinnedMonsterCombination_" + title_;
	if(f->Get(dirname.c_str())==0)
		f->mkdir(dirname.c_str());
	f->cd(dirname.c_str());
	
	string masstechdirname = dirname + "_massBinnedTechStuff";
	if(f->Get(masstechdirname.c_str())==0)
		f->mkdir(masstechdirname.c_str());
	
	if(!silent)
	  cout << dirname << endl;
	
	string pathTechnicalPNG = pathforPNG + "/TechnicalPlots/";
  mkdir(pathTechnicalPNG.c_str(),0777);
	
	string pathPNGMVABinnedResults = pathforPNG + "MVABinnedResults/";
  mkdir(pathPNGMVABinnedResults.c_str(),0777);
  
  string pathPNGMassBinnedPlots = pathforPNG + "MassBinnedPlots/";
  mkdir(pathPNGMassBinnedPlots.c_str(),0777);
  
  string pathPNGMassBinnedTechPlots = pathPNGMassBinnedPlots + "TechnicalPlots/";
  mkdir(pathPNGMassBinnedTechPlots.c_str(),0777);
  
  if( doMvaBinned_ )
    for(unsigned int i=0; i<monsterCombisMVAbinned_.size(); i++)
      monsterCombisMVAbinned_[i]->Write(f, pathPNGMVABinnedResults, false, 0, 0);
  
  float DEl[20] = {-9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9};
  float DElError[20] = {0.};
  float mTop[20] = {-9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9};
  float mTopError[20] = {0.};
  float xBin[20] = {-9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9};
  float xBinError[20] = {0.};
  
  float DElExpected[20] = {-9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9};
  float DElErrorExpected[20] = {0.};
  float DEbExpected[20] = {-9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9};
  float DEbErrorExpected[20] = {0.};
  float xBinExpected[20] = {-9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9};
  float xBinErrorExpected[20] = {0.};
  
  // mW and mTop Gauss-fit binned
  float mWgaus[20] = {-9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9};
  float mWErrorgaus[20] = {0.};
  float mTopgaus[20] = {-9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9};
  float mTopErrorgaus[20] = {0.};
  
  bool foundEntries = false;
  
  // variables to calculate the range of the plots
  float minDEl = 9999;
  float maxDEl = -9999;
  float minDEb = 9999;
  float maxDEb = -9999;
  float minmWgaus = 9999;
  float maxmWgaus = -9999;
  float minmTopgaus = 9999;
  float maxmTopgaus = -9999;
  
  for(unsigned int i=0; i<monsterCombis_.size(); i++)
  {
    f->cd(dirname.c_str());
    
    string calibCurveName = "";
    if(title_.find("chargeBins") == 0)
    {
      if(i == 0)
        calibCurveName = "MuMinus";
      else if(i == 1)
        calibCurveName = "MuPlus";
    }
    vector< vector<float> > result = monsterCombis_[i]->GetEstimation(pathTechnicalPNG, f, silent, calibCurveName);
    if(result.size() == 0)
    {
      cout << "BinnedMonsterCombination: " << title_ << ":  no input entries found" << endl;
      continue;
    }
    foundEntries = true;
    
    DEl[i] = result[0][0];
    DElError[i] = result[0][1];
    mTop[i] = result[1][0];
    mTopError[i] = result[1][1];
    xBin[i] = xBinHistos_[i]->GetMean();
    xBinError[i] = xBinHistos_[i]->GetMeanError();
    
    if(DEl[i]-DElError[i] < minDEl) minDEl = DEl[i]-DElError[i];
    if(DEl[i]+DElError[i] > maxDEl) maxDEl = DEl[i]+DElError[i];
    if(mTop[i]-mTopError[i] < minDEb) minDEb = mTop[i]-mTopError[i];
    if(mTop[i]+mTopError[i] > maxDEb) maxDEb = mTop[i]+mTopError[i];
    
    // expected stuff
    xBinExpected[i] = (binning_[i]+binning_[i+1]) / 2;
    xBinErrorExpected[i] = fabs(xBinExpected[i] - binning_[i]);
    
    pair<float,std::pair<float,float> > expL = monsterCombis_[i]->GetExpectedDEl();
    if(expL.first > -999.)
    {  
      DElExpected[i] = expL.first*100;
      DElErrorExpected[i] = expL.second.first*100;
    }
    
    pair<float,std::pair<float,float> > expB = monsterCombis_[i]->GetExpectedDEb();
    if(expB.first > -999.)
    {      
      DEbExpected[i] = expB.first*100;
      DEbErrorExpected[i] = expB.second.first*100;
    }
    
    if( DElExpected[i]-DElErrorExpected[i] < minDEl && DElExpected[i] != -9 && DElErrorExpected[i] < 4 ) minDEl = DElExpected[i]-DElErrorExpected[i];
    if( DElExpected[i]+DElErrorExpected[i] > maxDEl && DElExpected[i] != -9 && DElErrorExpected[i] < 4 ) maxDEl = DElExpected[i]+DElErrorExpected[i];
    if( ! topMass_ )
    {
      if( DEbExpected[i]-DEbErrorExpected[i] < minDEb && DEbExpected[i] != -9 && DEbErrorExpected[i] < 4 ) minDEb = DEbExpected[i]-DEbErrorExpected[i];
      if( DEbExpected[i]+DEbErrorExpected[i] > maxDEb && DEbExpected[i] != -9 && DEbErrorExpected[i] < 4 ) maxDEb = DEbExpected[i]+DEbErrorExpected[i];
    }

    if(!silent)
    {
      cout << " -  Expected DEl = " << expL.first*100 << " +- " << expL.second.first*100 << endl;
      cout << " -  Expected DEb = " << expB.first*100 << " +- " << expB.second.first*100 << endl;
    }
    
    // fit mW and mTop with gaussians
   	f->cd(masstechdirname.c_str());
   	
    TF1* fitfuncW = new TF1("fit_W","gaus");
    TH1F* histoW = binnedVars_["mWbin"][i];
    double rms = histoW->GetRMS();
    double maxbin =  histoW->GetBinCenter(histoW->GetMaximumBin());
    fitfuncW->SetRange(maxbin-1.2*rms,maxbin+1.2*rms);
    histoW->Fit(fitfuncW,"RQ");
    TCanvas* canvasW = TCanvasCreator(histoW, histoW->GetTitle());
    canvasW->SaveAs( (pathPNGMassBinnedTechPlots+"/"+histoW->GetTitle()+".png").c_str() );
    histoW->Write();
    
    TF1* fitfuncTop = new TF1("fit_Top","gaus");
    TH1F* histoTop = binnedVars_["mTopbin"][i];
    rms = histoTop->GetRMS();
    maxbin =  histoTop->GetBinCenter(histoTop->GetMaximumBin());
    fitfuncTop->SetRange(maxbin-1.2*rms,maxbin+1.2*rms);
    histoTop->Fit(fitfuncTop,"RQ");
    TCanvas* canvasTop = TCanvasCreator(histoTop, histoTop->GetTitle());
    canvasTop->SaveAs( (pathPNGMassBinnedTechPlots+"/"+histoTop->GetTitle()+".png").c_str() );
    histoTop->Write();
    
    mWgaus[i] = fitfuncW->GetParameter(1);
    mWErrorgaus[i] = fitfuncW->GetParError(1);
    mTopgaus[i] = fitfuncTop->GetParameter(1);
    mTopErrorgaus[i] = fitfuncTop->GetParError(1);
    
    if( mWgaus[i] + mWErrorgaus[i] > maxmWgaus && mWErrorgaus[i] < 20. && mWgaus[i] < 250. ) maxmWgaus = mWgaus[i] + mWErrorgaus[i];
    if( mWgaus[i] - mWErrorgaus[i] < minmWgaus && mWErrorgaus[i] < 20. && mWgaus[i] < 250. ) minmWgaus = mWgaus[i] - mWErrorgaus[i];
    
    if( mTopgaus[i] + mTopErrorgaus[i] > maxmTopgaus && mTopErrorgaus[i] < 20. && mTopgaus[i] < 350. ) maxmTopgaus = mTopgaus[i] + mTopErrorgaus[i];
    if( mTopgaus[i] - mTopErrorgaus[i] < minmTopgaus && mTopErrorgaus[i] < 20. && mTopgaus[i] < 350. ) minmTopgaus = mTopgaus[i] - mTopErrorgaus[i];
    delete fitfuncW;
    delete fitfuncTop;
  }
  if(!foundEntries) return;
  
  f->cd(dirname.c_str());
	TGraphErrors* DElExpectedBinnedGraph = new TGraphErrors(20,xBinExpected,DElExpected,xBinErrorExpected,DElErrorExpected);
  DElExpectedBinnedGraph->SetNameTitle((title_+"_ExpDElBinned").c_str(),(title_+"_ExpDElBinned").c_str());
  DElExpectedBinnedGraph->GetXaxis()->SetLimits(binning_[0],binning_[binning_.size()-1]);
  DElExpectedBinnedGraph->GetXaxis()->SetTitle(xTitle.c_str());
  DElExpectedBinnedGraph->GetYaxis()->SetRangeUser(minDEl-0.1,maxDEl+0.1);
  DElExpectedBinnedGraph->GetYaxis()->SetTitle("#DeltaE_{l} (%)");
	DElExpectedBinnedGraph->SetMarkerStyle(1);
	DElExpectedBinnedGraph->SetFillColor(3);
  
  TGraphErrors* DElBinnedGraph = new TGraphErrors(20,xBin,DEl,xBinError,DElError);
  DElBinnedGraph->GetXaxis()->SetLimits(binning_[0],binning_[binning_.size()-1]);
  DElBinnedGraph->SetNameTitle((title_+"_EstDElBinned").c_str(),(title_+"_EstDElBinned").c_str());
  TCanvas* canvasDEl = TCanvasCreator(DElExpectedBinnedGraph, title_+"_DElBinned", "AP2");
  canvasDEl->cd();
  DElBinnedGraph->Draw("P");
  canvasDEl->SaveAs( (pathforPNG+"/"+title_+"_DElBinned.png").c_str() );
  canvasDEl->Write();
  DElExpectedBinnedGraph->Write();
  DElBinnedGraph->Write();

  TGraphErrors* mTopBinnedGraph = new TGraphErrors(20,xBin,mTop,xBinError,mTopError);
  mTopBinnedGraph->GetXaxis()->SetTitle(xTitle.c_str());
  mTopBinnedGraph->GetXaxis()->SetLimits(binning_[0],binning_[binning_.size()-1]);
  mTopBinnedGraph->GetYaxis()->SetRangeUser(minDEb-0.1,maxDEb+0.1);
  if(topMass_)
  {
    mTopBinnedGraph->SetNameTitle((title_+"_EstmTopBinned").c_str(),(title_+"_EstmTopBinned").c_str());
    mTopBinnedGraph->GetYaxis()->SetTitle("m_{top}");
    TCanvas* canvasMtop = TCanvasCreator(mTopBinnedGraph, title_+"_mTopBinned", "AP");
    canvasMtop->SaveAs( (pathforPNG+"/"+title_+"_mTopBinned.png").c_str() );
    canvasMtop->Write();
  }
  else
  {
   	TGraphErrors* DEbExpectedBinnedGraph = new TGraphErrors(20,xBinExpected,DEbExpected,xBinErrorExpected,DEbErrorExpected);
   	DEbExpectedBinnedGraph->SetNameTitle((title_+"_ExpDEbBinned").c_str(),(title_+"_ExpDEbBinned").c_str());
    DEbExpectedBinnedGraph->GetXaxis()->SetTitle(xTitle.c_str());
    DEbExpectedBinnedGraph->GetXaxis()->SetLimits(binning_[0],binning_[binning_.size()-1]);
	  DEbExpectedBinnedGraph->SetMarkerStyle(1);
  	DEbExpectedBinnedGraph->SetFillColor(3);
  	DEbExpectedBinnedGraph->GetYaxis()->SetRangeUser(minDEb-0.1,maxDEb+0.1);

    mTopBinnedGraph->SetNameTitle((title_+"_EstDEbBinned").c_str(),(title_+"_EstDEbBinned").c_str());
    DEbExpectedBinnedGraph->GetYaxis()->SetTitle("#DeltaE_{b} (%)");
    TCanvas* canvasDEb = TCanvasCreator(DEbExpectedBinnedGraph, title_+"_DEbBinned", "AP2");
    canvasDEb->cd();
    mTopBinnedGraph->Draw("P");
    canvasDEb->SaveAs( (pathforPNG+"/"+title_+"_DEbBinned.png").c_str() );
    canvasDEb->Write();
    DEbExpectedBinnedGraph->Write();
  }
  mTopBinnedGraph->Write();
  
  TGraphErrors* mWgausBinnedGraph = new TGraphErrors(20,xBinExpected,mWgaus,xBinErrorExpected,mWErrorgaus);
  mWgausBinnedGraph->SetNameTitle((title_+"_mWgausBinned").c_str(),(title_+"_mWgausBinned").c_str());
  mWgausBinnedGraph->GetXaxis()->SetTitle(xTitle.c_str());
  mWgausBinnedGraph->GetYaxis()->SetTitle("m_{W}^{Gaus-fit}");
  mWgausBinnedGraph->GetXaxis()->SetLimits(binning_[0],binning_[binning_.size()-1]);
  mWgausBinnedGraph->GetYaxis()->SetRangeUser(minmWgaus-1.,maxmWgaus+1.);
  TCanvas* canvasMWgaus = TCanvasCreator(mWgausBinnedGraph, title_+"_mWgausBinned", "AP");
  canvasMWgaus->SaveAs( (pathPNGMassBinnedPlots+"/"+title_+"_mWgausBinned.png").c_str() );
  canvasMWgaus->Write();
  mWgausBinnedGraph->Write();
  
  TGraphErrors* mTopgausBinnedGraph = new TGraphErrors(20,xBinExpected,mTopgaus,xBinErrorExpected,mTopErrorgaus);
  mTopgausBinnedGraph->SetNameTitle((title_+"_mTopgausBinned").c_str(),(title_+"_mTopgausBinned").c_str());
  mTopgausBinnedGraph->GetXaxis()->SetTitle(xTitle.c_str());
  mTopgausBinnedGraph->GetYaxis()->SetTitle("m_{Top}^{Gaus-fit}");
  mTopgausBinnedGraph->GetXaxis()->SetLimits(binning_[0],binning_[binning_.size()-1]);
  mTopgausBinnedGraph->GetYaxis()->SetRangeUser(minmTopgaus-1.,maxmTopgaus+1.);
  TCanvas* canvasMTopgaus = TCanvasCreator(mTopgausBinnedGraph, title_+"_mTopgausBinned", "AP");
  canvasMTopgaus->SaveAs( (pathPNGMassBinnedPlots+"/"+title_+"_mTopgausBinned.png").c_str() );
  canvasMTopgaus->Write();
  mTopgausBinnedGraph->Write();
}

