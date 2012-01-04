#include "../interface/MonsterTools.h"

float maxProb(Monster* monster, int JetCombi)
{
  float max = 0;
  for(unsigned iJESl=0; iJESl<41; iJESl++)
  {
    for(unsigned iJESb=0; iJESb<41; iJESb++)
    {
      float prob = monster->fitResult(JetCombi, iJESl, iJESb);
      if(prob > max) max = prob;
    }
  }
  return max;
}

int nHoles(Monster* monster, bool topMass, int JetCombi)
{
  Int_t nHoles = 0, maxNb = 41;
  if(topMass) maxNb = 26;
	for (int iJESl=0; iJESl < 41; iJESl++)
 		for (int iJESb=0; iJESb < maxNb; iJESb++)
   		if (monster->fitResult(JetCombi, iJESl, iJESb) == 0)
				nHoles++;
  return nHoles;
}

float minDisHoleMax(TH2F* monster)
{
  float minDist = 9999.;
	Int_t maxBinx = 0; Int_t maxBiny = 0; Int_t maxBinz = 0;
  Int_t maxBin = monster->GetMaximumBin(); monster->GetBinXYZ(maxBin,maxBinx,maxBiny,maxBinz);
	for (int nBinX = 1; nBinX < monster->GetNbinsX()+1; nBinX++)
 	{
 		for (int nBinY = 1; nBinY < monster->GetNbinsY()+1; nBinY++)
 		{
   		if (monster->GetBinContent(nBinX,nBinY) == 0)
   		{
   		  float distHoleMax = sqrt( pow( (float) nBinX - maxBinx ,2) + pow( (float) nBinY - maxBiny ,2) );
				if( minDist > distHoleMax ) minDist = distHoleMax;
   		}
   	}
  }
  return minDist;
}

float probNoCorr(Monster* monster, TH2F* dummyMonster, bool topMass, int jetCombi)
{
  int binx = dummyMonster->GetXaxis()->FindBin(1.);
  int biny = -1;
  
  unsigned int* MVAres = monster->mvaResult(jetCombi);
  float mTop = ( monster->selectedJet(MVAres[0]) + monster->selectedJet(MVAres[1]) + monster->selectedJet(MVAres[2]) ).M();
  if(topMass)
    biny = dummyMonster->GetYaxis()->FindBin(mTop);
//    biny = histo.GetYaxis()->FindBin(172.5);
  else
    biny = dummyMonster->GetYaxis()->FindBin(1.);
  return monster->fitResult(jetCombi, binx-1, biny-1);
}

bool hadrJetsMVAMatched(LightMonster* monster, int jetCombi)
{
  bool goodCombi = false;
  unsigned int* MVAres = monster->mvaResult(jetCombi);
  int hadrBJetIndex = MVAres[2], lightJet1Index = MVAres[0], lightJet2Index = MVAres[1];
  if( hadrBJetIndex == monster->hadrBJet() && ( ( lightJet1Index == monster->hadrLJet1() && lightJet2Index == monster->hadrLJet2() )
    || ( lightJet2Index == monster->hadrLJet1() && lightJet1Index == monster->hadrLJet2() ) ) )
    goodCombi = true;
  return goodCombi;
}

bool hadrJetsMVAMatched(Monster* monster, int jetCombi)
{
  bool goodCombi = false;
  unsigned int* MVAres = monster->mvaResult(jetCombi);
  int hadrBJetIndex = MVAres[2], lightJet1Index = MVAres[0], lightJet2Index = MVAres[1];
  if( hadrBJetIndex == monster->hadrBJet() && ( ( lightJet1Index == monster->hadrLJet1() && lightJet2Index == monster->hadrLJet2() )
    || ( lightJet2Index == monster->hadrLJet1() && lightJet1Index == monster->hadrLJet2() ) ) )
    goodCombi = true;
  return goodCombi;
}

bool wJetsMVAMatched(Monster* monster, int jetCombi)
{
  bool goodCombi = false;
  unsigned int* MVAres = monster->mvaResult(jetCombi);
  int lightJet1Index = MVAres[0], lightJet2Index = MVAres[1];
  if( ( lightJet1Index == monster->hadrLJet1() && lightJet2Index == monster->hadrLJet2() )
    || ( lightJet2Index == monster->hadrLJet1() && lightJet1Index == monster->hadrLJet2() ) )
    goodCombi = true;
  return goodCombi;
}

int nWJetsMVAMatched(Monster* monster, int jetCombi)
{
  int nWJets=0;
  unsigned int* MVAres = monster->mvaResult(jetCombi);
  int lightJet1Index = MVAres[0], lightJet2Index = MVAres[1];
  if( lightJet1Index == monster->hadrLJet1() || lightJet1Index == monster->hadrLJet2() ) nWJets++;
  if( lightJet2Index == monster->hadrLJet1() || lightJet2Index == monster->hadrLJet2() ) nWJets++;
  return nWJets;
}

int nHadrBJetsMVAMatched(Monster* monster, int jetCombi)
{
  int nBJets=0;
  unsigned int* MVAres = monster->mvaResult(jetCombi);
  int hadrBJetIndex = MVAres[2];
  if( hadrBJetIndex == monster->hadrBJet() ) nBJets++;
  return nBJets;
}

TH2F* histoMonster(Monster* monster, TH2F* dummyMonster, int jetCombi)
{
//  TH2F* returnMonster = (TH2F*) dummyMonster->Clone();
  TH2F* returnMonster = (TH2F*) dummyMonster;
  for (int nBinX = 1; nBinX < returnMonster->GetNbinsX()+1; nBinX++)
	{
		for (int nBinY = 1; nBinY < returnMonster->GetNbinsY()+1; nBinY++)
		{
		  float kinFitResult = monster->fitResult(jetCombi, nBinX-1, nBinY-1);
 			returnMonster->SetBinContent(nBinX, nBinY, kinFitResult);
		}
	}
	return returnMonster;
}

vector<double> CalculateParabola(TH1D* parabola, int size, bool fit)
{ 
  vector<double> getParabolaValues;
  Int_t minBin = parabola->GetMinimumBin();

  //cout << " DChi2 value of minbin-1 " << parabola->GetBinContent(minBin-1) << endl;
  //cout << " DChi2 value of minbin " << parabola->GetBinContent(minBin) << endl;
  //cout << " DChi2 value of minbin+1 " << parabola->GetBinContent(minBin+1) << endl;
  
  if (!fit)
  {
    double x1 = parabola->GetBinCenter(minBin);
    double x2 = parabola->GetBinCenter(minBin+size);
    double x3 = parabola->GetBinCenter(minBin-size);
    double y1 = parabola->GetBinContent(minBin);
    double y2 = parabola->GetBinContent(minBin+size);
    double y3 = parabola->GetBinContent(minBin-size);
    
    //y = ax^2 + bx + c
    double a = (y3-y1)/((x3-x1)*(x3-x2))-(y2-y1)/((x2-x1)*(x3-x2));
    double b = (y2-y1)/(x2-x1)-(y3-y1)*(x2+x1)/((x3-x1)*(x3-x2))+(y2-y1)*(x2+x1)/((x2-x1)*(x3-x2));
    double c = y1-(y3-y1)*(x1*x1)/((x3-x1)*(x3-x2))+(y2-y1)*(x1*x1)/((x2-x1)*(x3-x2))-(y2-y1)*(x1)/(x2-x1)+(y3-y1)*(x2+x1)*(x1)/((x3-x1)*(x3-x2))-(y2-y1)*(x2+x1)*(x1)/((x2-x1)*(x3-x2));
    
    double minim = - b/(2*a); //minimum gives JES correction 
    double minimchi2 = a*minim*minim + b*minim +c;
    double minplus = (-b + sqrt(b*b-4*a*(c-minimchi2-1)))/(2*a); //upper bound right
    double minmin  = (-b - sqrt(b*b-4*a*(c-minimchi2-1)))/(2*a); //lower bound left
    
    double uncert1 = minplus-minim;
    double uncert2 = minim-minmin;
    double uncert=0.;
    if(uncert1 > uncert2) { uncert = uncert1; } else { uncert = uncert2; }
    
    getParabolaValues.push_back(a);
    getParabolaValues.push_back(b);
    getParabolaValues.push_back(c);
    getParabolaValues.push_back(minim);
    getParabolaValues.push_back(uncert);

    //cout << "Analytic min " << minim << " minmin " << minmin << " minplus " << minplus << " uncert " << uncert << endl;
  }
  else
  {
    double lowLim = parabola->GetBinCenter(minBin) - parabola->GetBinCenter(minBin - size); //(size*0.02);
    double highLim = parabola->GetBinCenter(minBin) + parabola->GetBinCenter(minBin + size); //(size*0.02);

    string func_title = "Fitted";
    TF1* func = new TF1(func_title.c_str(),"pol2");

    func->SetRange(lowLim,highLim);
    //func->FixParameter(0,0);
    
    parabola->Fit(func,"RQ0");

    double a = func->GetParameter(2);
    double b = func->GetParameter(1);
    double c = func->GetParameter(0);

    double minim = - b/(2*a); //minimum gives JES correction 
    double minimchi2 = a*minim*minim + b*minim +c;
    
    double minplus = (-b + sqrt(b*b-4*a*(c-minimchi2-1)))/(2*a); //upper bound right
    double minmin  = (-b - sqrt(b*b-4*a*(c-minimchi2-1)))/(2*a); //lower bound left
    
    double uncert1 = minplus-minim;
    double uncert2 = minim-minmin;
    double uncert=0.;
    if(uncert1 > uncert2) { uncert = uncert1; } else { uncert = uncert2; }

    getParabolaValues.push_back(a);
    getParabolaValues.push_back(b);
    getParabolaValues.push_back(c);
    getParabolaValues.push_back(minim);
    getParabolaValues.push_back(uncert);

//    cout << "a " << a << "  b " << b << "  c " << c << endl;
//    cout << "Fitted min " << minim << "  minmin " << minmin << "  minplus " << minplus << "  uncert " << uncert << endl;
  }

  return vector<double>(getParabolaValues);
}

float mTopMax(TH2F* monster)
{
	Int_t maxBinx = 0; Int_t maxBiny = 0; Int_t maxBinz = 0;
  Int_t maxBin = monster->GetMaximumBin(); monster->GetBinXYZ(maxBin,maxBinx,maxBiny,maxBinz);
  return monster->GetYaxis()->GetBinCenter(maxBiny);
}

vector< vector<float> > FitSuperMonster(TH2F* superMonster, string pathforPNG, TFile* f, string title, bool topMass)
{
  f->cd();
	string dirname = "PlotsMonsterCombination_" + title;
	if(f->Get(dirname.c_str())==0)
		f->mkdir(dirname.c_str());
	f->cd(dirname.c_str());
	superMonster->Write();
	
  float corrL = 0;
  float uncCorrL = 0;
  float corrB = 0;
  float uncCorrB = 0;
  float chi2 = -1;
  int ndf = -1;
	
  TH2F* histo = (TH2F*) superMonster->Clone();
  Int_t binx = 0; Int_t biny = 0; Int_t binz = 0;
  Int_t bin = histo->GetMinimumBin(); histo->GetBinXYZ(bin,binx,biny,binz);
  float min = histo->GetBinContent(binx,biny);
  
  // Shift bin content
  for(int i=0; i<=histo->GetNbinsX(); i++)
  {
    for(int j=0; j<=histo->GetNbinsY(); j++)
    {
//      histo->SetBinError(i, j, 1000000000000000.);
      histo->SetBinContent(i, j, histo->GetBinContent(i, j) - min );
    }
  }
  
	if(topMass)
	{
    int nBins = 9; // 2D-fit on a nBins x nBins grid
//    nBins = 18; // temporary for mass difference between top and anti-top
    bool goLeft = ( histo->GetBinContent(binx-1,biny) < histo->GetBinContent(binx+1,biny) );
    bool goDown = ( histo->GetBinContent(binx,biny-1) < histo->GetBinContent(binx,biny+1) );
    int nBinsLeft = -1, nBinsRight = -1, nBinsBottom = -1, nBinsTop = -1;
    
    if( (nBins % 2) != 0 )
    {
      nBinsLeft = nBins/2;
      nBinsRight = nBins/2;
      nBinsBottom = nBins/2;
      nBinsTop = nBins/2;
    }
    else
    {
      if(goLeft)
      {
        nBinsLeft = nBins/2;
        nBinsRight = nBins/2 - 1;
      }
      else
      {
        nBinsLeft = nBins/2 - 1;
        nBinsRight = nBins/2;
      }
      if(goDown)
      {
        nBinsBottom = nBins/2;
        nBinsTop = nBins/2 - 1;
      }
      else
      {
        nBinsBottom = nBins/2 - 1;
        nBinsTop = nBins/2;
      }
    }
    
    float stepXLeft = histo->GetXaxis()->GetBinCenter(binx-nBinsLeft) - histo->GetXaxis()->GetBinWidth(binx)/3;
    float stepXRight = histo->GetXaxis()->GetBinCenter(binx+nBinsRight) + histo->GetXaxis()->GetBinWidth(binx)/3;
    float stepYBottom = histo->GetYaxis()->GetBinCenter(biny-nBinsBottom) - histo->GetYaxis()->GetBinWidth(biny)/3;
    float stepYTop = histo->GetYaxis()->GetBinCenter(biny+nBinsTop) + histo->GetYaxis()->GetBinWidth(biny)/3;
	  
	  string func_title = title + "_2DFitted";
//	  TF2* func = new TF2(func_title.c_str(),"[0]+[1]*((x-[2])*cos([5])+(y-[4])*sin([5]))^2+[3]*(-(x-[2])*sin([5])+(y-[4])*cos([5]))^2");
//	  TF2* func = new TF2(func_title.c_str(),"[0]+[1]*(x*cos([5])+y*sin([5])-[2])^2+[3]*(-x*sin([5])+y*cos([5])-[4])^2"); extraction of the actual values is wrong for this function...
	  TF2* func = new TF2(func_title.c_str(),"[0]+[1]*(x-[2])^2+[3]*(y-[4])^2+[5]*(x-[2])*(y-[4])");
	  func->SetRange(stepXLeft, stepYBottom, stepXRight, stepYTop);
//	  func->SetParLimits(1,0,999999);
//	  func->SetParameter(1,1000);
//    func->SetParLimits(3,0,999999);
	  func->SetParLimits(2,0.4,1.5);
	  func->SetParameter(2,0.95);
	  func->SetParLimits(4,140,195);
	  func->SetParameter(4,172.5);
//	  histo->Fit(func,"RQ0 WW");
	  histo->Fit(func,"RQ0");
	  
	  histo->GetXaxis()->SetRangeUser(stepXLeft, stepXRight);
	  histo->GetYaxis()->SetRangeUser(stepYBottom, stepYTop);
	  
	  TH2F* histoRatio = PlotRatio(histo, func, func_title+"_Ratio");
	  
	  // create canvasses
    TCanvas* tempCanvas = TCanvasCreator(histo, func, func_title, "lego2"); //TCanvasCreator(pLight, func, func_title);
    tempCanvas->SaveAs( (pathforPNG+func_title+".png").c_str() );
    tempCanvas->Write();
    
    TCanvas* tempCanvas2 = TCanvasCreator(histo, title+"_SmallRange");
    tempCanvas2->SaveAs( (pathforPNG+title+"_SmallRange.png").c_str() );
    tempCanvas2->Write();
    
    TCanvas* tempCanvasRatio = TCanvasCreator(histoRatio, func_title+"_Ratio");
    tempCanvasRatio->SaveAs( (pathforPNG+func_title+"_Ratio.png").c_str() );
    tempCanvasRatio->Write();
    
    // calculate the errors
    float p0 = func->GetParameter(0), p1 = func->GetParameter(1), p2 = func->GetParameter(2), p3 = func->GetParameter(3), p4 = func->GetParameter(4);
    float minChi2 = func->Eval(p2,p4);
    
    float Dx = -4*p1*( p0 - ( minChi2 + 1 ) );
    float xErrPlus = fabs( ( 2*p1*p2 + sqrt(Dx) ) / (2*p1) - p2 );
    float xErrMin = fabs( ( 2*p1*p2 - sqrt(Dx) ) / (2*p1) - p2 );
    float xErr = xErrPlus;
    if(xErrMin > xErrPlus) xErr = xErrMin;
        
    float Dy = -4*p3*( p0 - ( minChi2 + 1 ) );
    float yErrPlus = fabs( ( 2*p3*p4 + sqrt(Dy) ) / (2*p3) - p4 );
    float yErrMin = fabs( ( 2*p3*p4 - sqrt(Dy) ) / (2*p3) - p4 );
    float yErr = yErrPlus;
    if(yErrMin > yErrPlus) yErr = yErrMin;
    
	  // store the analysis results
    corrL = -(1-p2)*100;
    uncCorrL = xErr*100;
    corrB = p4;
    uncCorrB = yErr;
    chi2 = func->GetChisquare();
    ndf = func->GetNDF();

    //clean up
    delete func;
	}
	else
	{
    // project to light flavour correction and shift to 0
    TH1D* pLight = histo->ProjectionX((title+"_projx").c_str(),biny,biny);
    pLight->GetXaxis()->SetTitle("#Delta E_{l}");
    pLight->GetYaxis()->SetTitle("#Delta #chi^{2}");
        
    double shift = pLight->GetBinContent(pLight->GetMinimumBin());

    for (int i=1; i<pLight->GetNbinsX()+1; i++)
    	pLight->SetBinContent(i,pLight->GetBinContent(i)-shift);
	  pLight->Write();

    // fit the parabola with a "pol2"
    vector<double> resultsF_l = CalculateParabola(pLight,5,true);
    double aF_l = resultsF_l[0];
    double bF_l = resultsF_l[1];
    double cF_l = resultsF_l[2];
        
    string func_title = string(pLight->GetName())+"_Fitted";

    TF1* func = new TF1(func_title.c_str(),"pol2",pLight->GetBinCenter(0)-(pLight->GetBinWidth(0)/2),pLight->GetBinCenter(pLight->GetNbinsX())+(pLight->GetBinWidth(0)/2));
    func->SetParameter(2,aF_l);
    func->SetParameter(1,bF_l);
    func->SetParameter(0,cF_l);
    
//    for(unsigned int i=0; i<50; i++)
//      cout << "bincenter[" << i << "] = " << histo->GetXaxis()->GetBinCenter(i) << endl;
    
    // set binx to DEl = 0
//    binx = histo->GetXaxis()->FindBin(1.);
//    cout << "new binx: " << binx << endl;
    
    // project to b flavour correction and shift to 0
    TH1D* pB = histo->ProjectionY((title+"_projy").c_str(),binx,binx);
    pB->GetXaxis()->SetTitle("#Delta E_{b}");
    pB->GetYaxis()->SetTitle("#Delta #chi^{2}");

    shift = pB->GetBinContent(pB->GetMinimumBin());
    for (int i=1; i<pB->GetNbinsX()+1; i++)
     	pB->SetBinContent(i,pB->GetBinContent(i)-shift);
	  pB->Write();

    // fit the parabola with a "pol2"
    vector<double> resultsF_b = CalculateParabola(pB,5,true);
    double aF_b = resultsF_b[0];
    double bF_b = resultsF_b[1];
    double cF_b = resultsF_b[2];

    string func_titleb = string(pB->GetName())+"_Fitted";

    TF1* funcb = new TF1(func_titleb.c_str(),"pol2",pB->GetBinCenter(0)-(pB->GetBinWidth(0)/2),pB->GetBinCenter(pB->GetNbinsX())+(pB->GetBinWidth(0)/2));
    funcb->SetParameter(2,aF_b);
    funcb->SetParameter(1,bF_b);
    funcb->SetParameter(0,cF_b);

    // create canvasses
    TCanvas* tempCanvasLight = TCanvasCreator(pLight, func, func_title);
    tempCanvasLight->SaveAs( (pathforPNG+func_title+".png").c_str() );
    tempCanvasLight->Write();
    TCanvas* tempCanvasB = TCanvasCreator(pB, funcb, func_titleb);
    tempCanvasB->SaveAs( (pathforPNG+func_titleb+".png").c_str() );
    tempCanvasB->Write();

    //clean up
    delete pLight;
    delete pB;
    delete func;
    delete funcb;
      
	  // show the analysis results :-)
    corrL = -(1-resultsF_l[3])*100;
    uncCorrL = resultsF_l[4]*100;
    if(topMass)
    {
      corrB = resultsF_b[3];
      uncCorrB = resultsF_b[4];
    }
    else
    {
      corrB = -(1-resultsF_b[3])*100;
      uncCorrB = resultsF_b[4]*100;
    }
  }
  
  vector< vector<float> > result;
  vector<float> lResult, bResult;
  lResult.push_back(corrL);
  lResult.push_back(uncCorrL);
  bResult.push_back(corrB);
  bResult.push_back(uncCorrB);
  result.push_back(lResult);
  result.push_back(bResult);
  return result;
}

void fillExpCorr(Monster *monster, ExpCorrCalculator *expCorr)
{
  if( hadrJetsMVAMatched(monster) )
  {
    float WJet1CorrFactor = 1;
    float WJet2CorrFactor = 1;
    float BJetCorrFactor = 1;

//    WJet1CorrFactor = selectedLooseJets_[hadronicWJet1_.first]->getJetCorrFactor("L1L2L3L4L5_uds") / selectedLooseJets_[hadronicWJet1_.first]->getJetCorrFactor("L1L2L3");
//    WJet2CorrFactor = selectedLooseJets_[hadronicWJet2_.first]->getJetCorrFactor("L1L2L3L4L5_uds") / selectedLooseJets_[hadronicWJet2_.first]->getJetCorrFactor("L1L2L3");
//    BJetCorrFactor = selectedLooseJets_[hadronicBJet_.first]->getJetCorrFactor("L1L2L3L4L5_b") / selectedLooseJets_[hadronicBJet_.first]->getJetCorrFactor("L1L2L3");

    vector<TLorentzVector> vTLV;
    vTLV.push_back( monster->hadrLJet1() * WJet1CorrFactor );
    vTLV.push_back( monster->hadrLJet2() * WJet2CorrFactor );
    vTLV.push_back( monster->hadrBJet() * BJetCorrFactor );
    
    vector<TLorentzVector> vTLVmc;
    vTLVmc.push_back( monster->hadrLQuark1() );
    vTLVmc.push_back( monster->hadrLQuark2() );
    vTLVmc.push_back( monster->hadrBQuark() );
    
    expCorr->Fill3Jets(vTLV);
    expCorr->Fill3MCParticles(vTLVmc);
    expCorr->FillMaxMVA(monster->mvaVal(0));
    
    // fill with dummy stuff
    vector<TLorentzVector> dummyVTLV;
    dummyVTLV.push_back(TLorentzVector(0.1,0.1,999,0.1));
    expCorr->FillMuons(dummyVTLV);
    expCorr->FillElectrons(dummyVTLV);
    
    float PtLight = (vTLV[0].Pt() + vTLV[1].Pt()) / 2;
    expCorr->FillLightJet( (vTLVmc[0].Et() - vTLV[0].Et()) / vTLV[0].Et(), PtLight, vTLV[2].Pt() );
    expCorr->FillLightJet( (vTLVmc[1].Et() - vTLV[1].Et()) / vTLV[1].Et(), PtLight, vTLV[2].Pt() );
    expCorr->FillBJet( (vTLVmc[2].Et() - vTLV[2].Et()) / vTLV[2].Et(), PtLight, vTLV[2].Pt() );
  }
}

