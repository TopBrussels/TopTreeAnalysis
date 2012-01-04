#include "../interface/MVABinnedMonsterCombination.h"

MVABinnedMonsterCombination::MVABinnedMonsterCombination(string title, bool topMass, vector<float> binning)
{
  title_ = title;
  topMass_ = topMass;
  binning_ = binning;
  
  maxMVA_ = new TH1F((title+"_maxMVA").c_str(),(title+"_maxMVA").c_str(),50,0,1);
  
  expDElHisto_ = new TH1F((title_+"_Expected_DEl").c_str(),(title_+"_Expected_DEl").c_str(),80,-1,1);
  expDEbHisto_ = new TH1F((title_+"_Expected_DEb").c_str(),(title_+"_Expected_DEb").c_str(),80,-1,1);
  
  vector<TH1F*> ptLbinnedHistos, ptBbinnedHistos;
  for(unsigned int i=0; i<binning_.size()-1; i++)
  {
    stringstream ss1, ss2;
    ss1 << binning_[i];
    ss2 << binning_[i+1];
    string binTitle = title + "_Bin_" + ss1.str() + "_" + ss2.str();
    monsterCombis_.push_back( new MonsterCombination(binTitle, topMass_) );
    xBinHistos_.push_back( new TH1F((binTitle+"_xBin").c_str(),(binTitle+"_xBin").c_str(),10000,binning_[i]-0.1,binning[i+1]+0.1) );
    ptLbinnedHistos.push_back( new TH1F((binTitle+"_PtLbin").c_str(), (binTitle+"_PtLbin").c_str(), 10000, 0, 3500) );
    ptBbinnedHistos.push_back( new TH1F((binTitle+"_PtBbin").c_str(), (binTitle+"_PtBbin").c_str(), 10000, 0, 3500) );
  }
  mvaBinnedVars_["PtLbin"] = ptLbinnedHistos;
  mvaBinnedVars_["PtBbin"] = ptBbinnedHistos;
}

MVABinnedMonsterCombination::~MVABinnedMonsterCombination()
{
  delete maxMVA_;
  delete expDElHisto_;
  delete expDEbHisto_;
  for(unsigned int i=0; i<xBinHistos_.size(); i++)
    delete xBinHistos_[i];
}

void MVABinnedMonsterCombination::Fill(Monster* monster, float weight, int jetCombi)
{
  unsigned int* MVAres = monster->mvaResult(jetCombi);
  float binValue = monster->mvaVal(jetCombi);
  for(unsigned int i=0; i<monsterCombis_.size(); i++)
  {
    if(binValue > binning_[i] && binValue <= binning_[i+1])
    {
      monsterCombis_[i]->addMonster(monster, weight, jetCombi);
      xBinHistos_[i]->Fill(binValue);
      monsterCombis_[i]->addExpected(monster, weight, jetCombi);
      
      mvaBinnedVars_["PtLbin"][i]->Fill(monster->selectedJet(MVAres[0]).Pt());
      mvaBinnedVars_["PtLbin"][i]->Fill(monster->selectedJet(MVAres[1]).Pt());
      mvaBinnedVars_["PtBbin"][i]->Fill(monster->selectedJet(MVAres[2]).Pt());
    }
  }
  
  if( hadrJetsMVAMatched(monster, jetCombi) )
  {
    expDElHisto_->Fill((monster->hadrLQuark1().Et() - monster->selectedJet(MVAres[0]).Et())/monster->selectedJet(MVAres[0]).Et());
    expDElHisto_->Fill((monster->hadrLQuark2().Et() - monster->selectedJet(MVAres[1]).Et())/monster->selectedJet(MVAres[1]).Et());
    expDEbHisto_->Fill((monster->hadrBQuark().Et() - monster->selectedJet(MVAres[2]).Et())/monster->selectedJet(MVAres[2]).Et());
  }
  
  maxMVA_->Fill(binValue, weight);
}

void MVABinnedMonsterCombination::FillExpected(Monster* monster, float binValue, float weights)
{
  for(unsigned int i=0; i<monsterCombis_.size(); i++)
    if(binValue > binning_[i] && binValue <= binning_[i+1])
      monsterCombis_[i]->addExpected(monster, weights);
}

void MVABinnedMonsterCombination::Write(TFile* f, string pathforPNG, bool templateFit, TH1F* goodMVAtemplate, TH1F* badMVAtemplate, string xTitle, bool silent)
{
  f->cd();
	string dirname = "maxMVABinned_MonsterCombination_" + title_;
	if(f->Get(dirname.c_str())==0)
		f->mkdir(dirname.c_str());
	f->cd(dirname.c_str());
	
	string pathTechnicalPNG = pathforPNG + "/TechnicalPlots/";
  mkdir(pathTechnicalPNG.c_str(),0777);
	
  float DEl[20] = {-9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9};
  float DElError[20] = {0.};
  float mTop[20] = {-9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9};
  float mTopError[20] = {0.};
  float xBin[20] = {-9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9};
  float xBinError[20] = {0.};
  
  bool foundEntries = false;
  
  // template fitting to max MVA distribution
  float** resultsTemplateFit = {0};
  if(templateFit)
    resultsTemplateFit = estimateGoodCombisMVA(binning_, goodMVAtemplate, badMVAtemplate, maxMVA_, pathTechnicalPNG, f, "maxMVA_templateFit_"+title_+"_");
  
  // variables to calculate the range of the plots
  float minDEl = 9999;
  float maxDEl = -9999;
  float minDEb = 9999;
  float maxDEb = -9999;
  
  for(unsigned int i=0; i<monsterCombis_.size(); i++)
  {
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
      cout << "MVABinnedMonsterCombination: " << title_ << ":  no input entries found" << endl;
      continue;
    }
    foundEntries = true;
    
    DEl[i] = result[0][0];
    DElError[i] = result[0][1];
    mTop[i] = result[1][0];
    mTopError[i] = result[1][1];
    if(templateFit)
    {
      xBin[i] = resultsTemplateFit[i][0]; // xBinHistos_[i]->GetMean();
      xBinError[i] = resultsTemplateFit[i][1]; // xBinHistos_[i]->GetMeanError();
    }
    else
    {
      xBin[i] = 0; // xBinHistos_[i]->GetMean();
      xBinError[i] = 1; // xBinHistos_[i]->GetMeanError();
    }
    if(DEl[i]-DElError[i] < minDEl) minDEl = DEl[i]-DElError[i];
    if(DEl[i]+DElError[i] > maxDEl) maxDEl = DEl[i]+DElError[i];
    if(mTop[i]-mTopError[i] < minDEb) minDEb = mTop[i]-mTopError[i];
    if(mTop[i]+mTopError[i] > maxDEb) maxDEb = mTop[i]+mTopError[i];
  }
  
  if(!foundEntries) return;
  
  // Expected stuff
  float expDEl[1] = {-9};
  float expDElError[1] = {-9};
  float expDEb[1] = {-9};
  float expDEbError[1] = {-9};
  float xExpected[1] = {0.5};
  float xExpectedError[1] = {0.5};
  
  if(expDElHisto_->Integral(1, 80) > 0.1)
  {
    pair<float,std::pair<float,float> > expL = GetExpectation(expDElHisto_, "gaus");
    
    expDEl[0] = expL.first*100;
    expDElError[0] = expL.second.first*100;
  }
  TCanvas* tmpExpLCanv = TCanvasCreator(expDElHisto_,expDElHisto_->GetName());
  tmpExpLCanv->SaveAs( (pathTechnicalPNG+"/"+expDElHisto_->GetName()+".png").c_str() );
  tmpExpLCanv->Write();

  if(expDEbHisto_->Integral(1, 80) > 0.1)
  {
    pair<float,std::pair<float,float> > expB = GetExpectation(expDEbHisto_, "gaus");
    
    expDEb[0] = expB.first*100;
    expDEbError[0] = expB.second.first*100;
  }
  TCanvas* tmpExpBCanv = TCanvasCreator(expDEbHisto_,expDEbHisto_->GetName());
  tmpExpBCanv->SaveAs( (pathTechnicalPNG+"/"+expDEbHisto_->GetName()+".png").c_str() );
  tmpExpBCanv->Write();
  
  if(expDEl[0]-expDElError[0] < minDEl) minDEl = expDEl[0]-expDElError[0];
  if(expDEl[0]+expDElError[0] > maxDEl) maxDEl = expDEl[0]+expDElError[0];
  if( ! topMass_ )
  {
    if(expDEb[0]-expDEbError[0] < minDEb) minDEb = expDEb[0]-expDEbError[0];
    if(expDEb[0]+expDEbError[0] > maxDEb) maxDEb = expDEb[0]+expDEbError[0];
  }
  
  TGraphErrors* DElBinnedGraph = new TGraphErrors(20,xBin,DEl,xBinError,DElError);
  TGraphErrors* DElExpectedGraph = new TGraphErrors(2,xExpected,expDEl,xExpectedError,expDElError);
  DElExpectedGraph->GetXaxis()->SetTitle(xTitle.c_str());
  DElExpectedGraph->GetXaxis()->SetLimits(0.,1.1);
  DElExpectedGraph->GetYaxis()->SetRangeUser(minDEl-1,maxDEl+1);
  DElExpectedGraph->GetYaxis()->SetTitle("#DeltaE_{l} (%)");
 	DElExpectedGraph->SetMarkerStyle(1);
	DElExpectedGraph->SetFillColor(3);
  
  TCanvas* canvasDEl = TCanvasCreator(DElExpectedGraph, title_+"_DElBinned", "AP2");
  canvasDEl->cd();
  DElBinnedGraph->Draw("P");
  canvasDEl->SaveAs( (pathforPNG+"/"+title_+"_DElBinned.png").c_str() );
  canvasDEl->Write();
  delete canvasDEl;
  delete DElBinnedGraph;
  delete DElExpectedGraph;
  
  TGraphErrors* mTopBinnedGraph = new TGraphErrors(20,xBin,mTop,xBinError,mTopError);
  if(topMass_)
  {
    mTopBinnedGraph->GetXaxis()->SetTitle(xTitle.c_str());
    mTopBinnedGraph->GetXaxis()->SetLimits(0.,1.1);
    mTopBinnedGraph->GetYaxis()->SetRangeUser(minDEb-1,maxDEb+1);
    mTopBinnedGraph->SetNameTitle((title_+"_mTopBinned").c_str(),(title_+"_mTopBinned").c_str());
    mTopBinnedGraph->GetYaxis()->SetTitle("m_{top}");
    TCanvas* canvasMtop = TCanvasCreator(mTopBinnedGraph, title_+"_mTopBinned");
    canvasMtop->SaveAs( (pathforPNG+"/"+title_+"_mTopBinned.png").c_str() );
    canvasMtop->Write();
    delete canvasMtop;
  }
  else
  {
    TGraphErrors* DEbExpectedGraph = new TGraphErrors(2,xExpected,expDEb,xExpectedError,expDEbError);
    DEbExpectedGraph->GetXaxis()->SetTitle(xTitle.c_str());
    DEbExpectedGraph->GetXaxis()->SetLimits(0.,1.1);
    DEbExpectedGraph->GetYaxis()->SetRangeUser(minDEb-1,maxDEb+1);
    DEbExpectedGraph->GetYaxis()->SetTitle("#DeltaE_{b} (%)");
    DEbExpectedGraph->SetMarkerStyle(1);
    DEbExpectedGraph->SetFillColor(3);
    
    TCanvas* canvasDEb = TCanvasCreator(DEbExpectedGraph, title_+"_DEbBinned", "AP2");
    canvasDEb->cd();
    mTopBinnedGraph->Draw("P");
    canvasDEb->SaveAs( (pathforPNG+"/"+title_+"_DEbBinned.png").c_str() );
    canvasDEb->Write();
    delete canvasDEb;
    delete DEbExpectedGraph;
  }
  delete mTopBinnedGraph;
  
  // Plot results in bins of the MVA
  float xBinMVA[20] = {-9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9};
  float xBinMVAError[20] = {-9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9};
  
  for(unsigned int i=0; i<xBinHistos_.size(); i++)
  {
    xBinMVA[i] = xBinHistos_[i]->GetMean();
    xBinMVAError[i] = xBinHistos_[i]->GetMeanError();
  }
  
  TGraphErrors* DElMVAbinned = new TGraphErrors(20,xBinMVA,DEl,xBinMVAError,DElError);
  DElMVAbinned->GetXaxis()->SetLimits(0.,1.1);
  DElMVAbinned->GetXaxis()->SetTitle("maxMVA");
  DElMVAbinned->GetYaxis()->SetTitle("#DeltaE_{l}");
  
  TCanvas* canvasDElMVAbinned = TCanvasCreator(DElMVAbinned, title_+"_DEl_maxMVABinned", "AP");
  canvasDElMVAbinned->SaveAs( (pathforPNG+"/"+title_+"_DEl_maxMVABinned.png").c_str() );
  canvasDElMVAbinned->Write();
  delete DElMVAbinned;
  delete canvasDElMVAbinned;
  
  TGraphErrors* mTopMVAbinned = new TGraphErrors(20,xBinMVA,mTop,xBinMVAError,mTopError);
  mTopMVAbinned->GetXaxis()->SetLimits(0.,1.1);
  mTopMVAbinned->GetXaxis()->SetTitle("maxMVA");
  if(topMass_)
  {
    mTopMVAbinned->GetYaxis()->SetTitle("m_{top}");
    
    TCanvas* canvasmTopMVAbinned = TCanvasCreator(mTopMVAbinned, title_+"_mTop_maxMVABinned", "AP");
    canvasmTopMVAbinned->SaveAs( (pathforPNG+"/"+title_+"_mTop_maxMVABinned.png").c_str() );
    canvasmTopMVAbinned->Write();
    delete canvasmTopMVAbinned;
  }
  else
  {
    mTopMVAbinned->GetYaxis()->SetTitle("#DeltaE_{b}");
    
    TCanvas* canvasDEbMVAbinned = TCanvasCreator(mTopMVAbinned, title_+"_DEb_maxMVABinned", "AP");
    canvasDEbMVAbinned->SaveAs( (pathforPNG+"/"+title_+"_DEb_maxMVABinned.png").c_str() );
    canvasDEbMVAbinned->Write();
    delete canvasDEbMVAbinned;
  }
  delete mTopMVAbinned;
  
  // properties in maxMVA-bins  
  for(map<string, vector<TH1F*> >::const_iterator it = mvaBinnedVars_.begin(); it != mvaBinnedVars_.end(); it++)
  {
    float yBin[20] = {-9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9};
    float yBinError[20] = {-9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9};
    float yMin = 99999, yMax = -99999;
    
    for(unsigned int i=0; i<it->second.size(); i++)
    {
      yBin[i] = it->second[i]->GetMean();
      yBinError[i] = it->second[i]->GetMeanError();
      if(yBin[i]+yBinError[i] > yMax) yMax = yBin[i]+yBinError[i];
      if(yBin[i]-yBinError[i] < yMin) yMin = yBin[i]-yBinError[i];
    }
    
    TGraphErrors* graph = new TGraphErrors(20,xBinMVA,yBin,xBinMVAError,yBinError);
    graph->GetXaxis()->SetTitle("maxMVA");
    graph->GetYaxis()->SetTitle( (it->first).c_str() );
    graph->GetXaxis()->SetLimits(0.,1.1);
    graph->GetYaxis()->SetRangeUser(yMin-1, yMax+1);
    
    string graphTitle = title_+"_"+it->first+"_maxMVABinned";
    TCanvas* canvasGraph = TCanvasCreator(graph, graphTitle, "AP");
    canvasGraph->SaveAs( (pathforPNG+"/"+graphTitle+".png").c_str() );
    canvasGraph->Write();
    delete canvasGraph;
    delete graph;
  }
}

