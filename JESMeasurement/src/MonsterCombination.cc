#include "../interface/MonsterCombination.h"

MonsterCombination::MonsterCombination(string title, bool topMass)
{
  FullKinFit tmpKinfit = FullKinFit(0, 0, 0, topMass, false, true);
  superMonster_ = tmpKinfit.DummyMonster();
  superMonster_->SetNameTitle(("superMonster_"+title).c_str(),("superMonster_"+title).c_str());
  title_ = title;
  topMass_ = topMass;
  nEventsWeighted_ = 0;
  
  expDElHisto_ = TH1F((title_+"_ExpDEl").c_str(), (title_+"_ExpDEl").c_str(),80,-1,1);
  expDEbHisto_ = TH1F((title_+"_ExpDEb").c_str(), (title_+"_ExpDEb").c_str(),80,-1,1);
}

MonsterCombination::~MonsterCombination()
{
  if(superMonster_) delete superMonster_;
}

void MonsterCombination::addMonster(Monster* monster, float weight, int jetCombi)
{
//  hadrLJet1_.push_back(monster->hadrLJet1());
//  hadrLJet2_.push_back(monster->hadrLJet2());
//  hadrBJet_.push_back(monster->hadrBJet());
  nEventsWeighted_ += weight;
  
	for (int nBinX = 1; nBinX < superMonster_->GetNbinsX()+1; nBinX++)
	{
		for (int nBinY = 1; nBinY < superMonster_->GetNbinsY()+1; nBinY++)
		{
		  float kinFitResult = monster->fitResult(jetCombi, nBinX-1, nBinY-1);
		  int globalBinNr = superMonster_->GetBin(nBinX, nBinY);
 			if (kinFitResult > 0)
   			superMonster_->AddBinContent(globalBinNr,-2*log(kinFitResult)*weight);
 			else
 			  superMonster_->AddBinContent(globalBinNr,-2*log(pow(10.,-16.))*weight);
		}
	}
}

void MonsterCombination::addExpected(Monster* monster, float weight, int jetCombi)
{
  // Fill the expected corrections
  if( hadrJetsMVAMatched(monster, jetCombi) )
  {
    unsigned int* MVAres = monster->mvaResult(jetCombi);
    expDElHisto_.Fill((monster->hadrLQuark1().Et() - (monster->selectedJet(MVAres[0])).Et())/(monster->selectedJet(MVAres[0])).Et());
    expDElHisto_.Fill((monster->hadrLQuark2().Et() - (monster->selectedJet(MVAres[1])).Et())/(monster->selectedJet(MVAres[1])).Et());
    expDEbHisto_.Fill((monster->hadrBQuark().Et() - (monster->selectedJet(MVAres[2])).Et())/(monster->selectedJet(MVAres[2])).Et());
  }
}

TH2F* MonsterCombination::AverageMonster()
{
  TH2F* avMonster = (TH2F*)superMonster_->Clone();
  for (int nBinX = 1; nBinX < superMonster_->GetNbinsX()+1; nBinX++)
	{
		for (int nBinY = 1; nBinY < superMonster_->GetNbinsY()+1; nBinY++)
		{
		  int globalBinNr = superMonster_->GetBin(nBinX, nBinY);
		  if(globalBinNr == 100)
  		  cout << "avmonster before reweighting: " << avMonster->GetBinContent(globalBinNr) << "   " << superMonster_->GetBinContent(globalBinNr) << endl;
		  avMonster->SetBinContent(globalBinNr, superMonster_->GetBinContent(globalBinNr) / nEventsWeighted_ );
		  if(globalBinNr == 100)
		    cout << "avmonster after reweighting: " << avMonster->GetBinContent(globalBinNr) << endl;
		}
  }
  return avMonster;
}

void MonsterCombination::SubTractAvMonster(TH2F* avMonster, float nBadCombis)
{
  for (int nBinX = 1; nBinX < superMonster_->GetNbinsX()+1; nBinX++)
	{
		for (int nBinY = 1; nBinY < superMonster_->GetNbinsY()+1; nBinY++)
		{
		  int globalBinNr = superMonster_->GetBin(nBinX, nBinY);
		  if(globalBinNr == 100)
  		  cout << "before subtraction: " << superMonster_->GetBinContent(globalBinNr) << endl;
		  superMonster_->SetBinContent(globalBinNr, superMonster_->GetBinContent(globalBinNr) - (nBadCombis * avMonster->GetBinContent(globalBinNr) ) );
		  if(globalBinNr == 100)
  		  cout << "after subtraction: " << superMonster_->GetBinContent(globalBinNr) << endl;
		}
	}
}

pair<float, pair<float,float> > MonsterCombination::GetExpectedDEl()
{
  if(expDElHisto_.Integral(1, 80) > 0.1)
    return GetExpectation(&expDElHisto_, "gaus");
  pair<float,float> dummyPair = pair<float,float>(-9999., -9999.);
  return pair<float, pair<float,float> > (-9999., dummyPair);
}

pair<float, pair<float,float> > MonsterCombination::GetExpectedDEb()
{
  if(expDEbHisto_.Integral(1, 80) > 0.1)
    return GetExpectation(&expDEbHisto_, "gaus");
  pair<float,float> dummyPair = pair<float,float>(-9999., -9999.);
  return pair<float, pair<float,float> > (-9999., dummyPair);
}

vector< vector<float> > MonsterCombination::GetEstimation(string pathforPNG, TFile* f, bool silent, string calibCurveName) //get the estimated corrections
{
	if(superMonster_ == NULL)
  {
    cout << "superMonster_ for " << title_ << " does not exist, returning..." << endl;
    vector< vector<float> > empty;
	  return empty;
  }
	if(superMonster_->GetBinContent(20, 20) < 1 && superMonster_->GetEntries() < 1)
	{
	  cout << "Not enough entries in the SuperMonster: " << title_ << " : "<< superMonster_->GetEntries() << " So, not calculating the final results for this..." << endl;
	  vector< vector<float> > empty;
	  return empty;
	}
	TH2F* superMonsterNoSubtr = (TH2F*) superMonster_->Clone();
  vector< vector<float> > result = FitSuperMonster(superMonsterNoSubtr, pathforPNG, f, title_, topMass_);
  float corrL = result[0][0];
  float uncCorrL = result[0][1];
  float corrB = result[1][0];
  float uncCorrB = result[1][1];
  
	f->cd();
	string dirname = "PlotsMonsterCombination_" + title_;
	if(f->Get(dirname.c_str())==0)
		f->mkdir(dirname.c_str());
	f->cd(dirname.c_str());
  
  TCanvas* tmpExpLCanv = TCanvasCreator(&expDElHisto_,expDElHisto_.GetName());
  tmpExpLCanv->SaveAs( (pathforPNG+"/"+expDElHisto_.GetName()+".png").c_str() );
  tmpExpLCanv->Write();
  
  TCanvas* tmpExpBCanv = TCanvasCreator(&expDEbHisto_,expDEbHisto_.GetName());
  tmpExpBCanv->SaveAs( (pathforPNG+"/"+expDEbHisto_.GetName()+".png").c_str() );
  tmpExpBCanv->Write();
  
/*  if(topMass_)
  {
    float* corrUnc = correctMTopErrors(uncCorrB, uncCorrL); // in this case, corrB is the top mass ;-)
    uncCorrL = corrUnc[1];
    uncCorrB = corrUnc[0];
    
    float** calibResults = calibrateMTop(corrB, uncCorrB, corrL, uncCorrL, calibCurveName); // apply calibrationCurve
    corrB = calibResults[0][0];
    uncCorrB = calibResults[0][1];
    corrL = calibResults[1][0];
    uncCorrL = calibResults[1][1];
  }
  else
  {
    float* corrUnc = correctMTopErrors(uncCorrL, uncCorrB);
    uncCorrL = corrUnc[1];
    uncCorrB = corrUnc[0];
  }*/
  
  if( ! silent )
  {
    cout << "Results for " << title_ << endl;
    cout << " -  Estimated light Correction = " << corrL << "+-" << uncCorrL << " %" << endl;
    if(topMass_)
    {
      cout << " -  Estimated  Top  Mass  = " << corrB << "+-" << uncCorrB << " GeV"<< endl;
    }
    else
    {
      cout << " -  Estimated  b  Correction = " << corrB << "+-" << uncCorrB << " %"<< endl;
    }
  }
  
/*  // Apply the corrections
  map<string,TH1F*> histo1D;
  
	histo1D["mW_noCorr"] = new TH1F("mW_noCorr","W boson mass noCorr;m_{W};#events",100,0,200);
	histo1D["mtop_noCorr"] = new TH1F("mtop_noCorr","top quark mass noCorr;m_{t};#events",100,0,400);
	
	histo1D["mW_Corr"] = new TH1F("mW_Corr","W boson mass Corrected;m_{W};#events",100,0,200);
 	histo1D["mtop_Corr"] = new TH1F("mtop_Corr","top quark mass Corrected;m_{t};#events",100,0,400);
  
  // test iterative correction factors
//  corrL = -13.14;
//  corrB = -10.46;
//  corrL = -2.95;
//  corrB = -0.93;
  float scale = 1.;
  
  for(unsigned int i=0; i<hadrLJet1_.size(); i++)
  {
    TLorentzVector unCorrLJet1 = (hadrLJet1_[i]*scale);
    TLorentzVector unCorrLJet2 = (hadrLJet2_[i]*scale);
    TLorentzVector unCorrBJet = (hadrBJet_[i]*scale);
    
    float mWnoCorr = (unCorrLJet1+unCorrLJet2).M();
    float mTopNoCorr = (unCorrLJet1+unCorrLJet2+unCorrBJet).M();
//    cout << "corrL = " << corrL << "  unCorr Pt = " << (hadrLJet1_[i]).Pt() << "  scaledPt = " << (unCorrLJet1).Pt() << "  corr Pt = " << (unCorrLJet1*(1.+0.01*corrL)).Pt() << endl;
    float mWcorr = (unCorrLJet1*(1.+0.01*corrL) + unCorrLJet2*(1.+0.01*corrL)).M();
    float mTopCorr = 0;
    if(topMass_) mTopCorr = (unCorrLJet1*(1.+0.01*corrL) + unCorrLJet2*(1.+0.01*corrL) + unCorrBJet).M();
    else mTopCorr = (unCorrLJet1*(1.+0.01*corrL) + unCorrLJet2*(1.+0.01*corrL) + unCorrBJet*(1.+0.01*corrB)).M();
//    cout << "mWunCorr = " << (hadrLJet1_[i]+hadrLJet2_[i]).M() << "  mWscaled = " << mWnoCorr << "  mWcorr = " << mWcorr << endl;
    
    histo1D["mW_noCorr"]->Fill(mWnoCorr);
    histo1D["mtop_noCorr"]->Fill(mTopNoCorr);
    histo1D["mW_Corr"]->Fill(mWcorr);
    histo1D["mtop_Corr"]->Fill(mTopCorr);
  }
  
  float maxMWnoCorr = histo1D["mW_noCorr"]->GetXaxis()->GetBinCenter(histo1D["mW_noCorr"]->GetMaximumBin());
  float maxMTopNoCorr = histo1D["mtop_noCorr"]->GetXaxis()->GetBinCenter(histo1D["mtop_noCorr"]->GetMaximumBin());
  float maxMWcorr = histo1D["mW_Corr"]->GetXaxis()->GetBinCenter(histo1D["mW_Corr"]->GetMaximumBin());
  float maxMTopCorr = histo1D["mtop_Corr"]->GetXaxis()->GetBinCenter(histo1D["mtop_Corr"]->GetMaximumBin());
  
  histo1D["mW_noCorr"]->Fit("gaus","Q","",maxMWnoCorr-10., maxMWnoCorr+10.);
	histo1D["mtop_noCorr"]->Fit("gaus","Q","",maxMTopNoCorr-15., maxMTopNoCorr+15.);
	
	histo1D["mW_Corr"]->Fit("gaus","Q","",maxMWcorr-10., maxMWcorr+10.);
 	histo1D["mtop_Corr"]->Fit("gaus","Q","",maxMTopCorr-15., maxMTopCorr+15.);
  
  // Write 1D histo's
  for(std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
	{
		TH1F *temp = it->second;
		temp->Write();
		TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
		tempCanvas->SaveAs( (pathforPNG+"/"+it->first+"_"+title_+".png").c_str() );
		delete it->second;
		delete tempCanvas;
	}*/
  
  return result;
}

