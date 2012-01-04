#include "../interface/JESTools.h"

using namespace std;

pair<float,pair<float,float> > GetExpectation(TH1F* histo, string method)
{
  std::map<double, std::pair<float, float> > values;
  string func_title = string(histo->GetName())+"_Fitted";
  
  if(method.compare("gaus") == 0)
  {
    for (unsigned int i=0; i<3; i++)
    {
//      unsigned int i=1;
      double nRMS = 0;
      if (i==0) nRMS=1; if (i==1) nRMS=1.5; if (i==2) nRMS=2.;
//      if (i==0) nRMS=0.5; if (i==1) nRMS=1; if (i==2) nRMS=1.5;

      TF1 *fitfunc_1, *fitfunc_2, *fitfunc_3;
      double rms = histo->GetRMS();
      double maxbin =  histo->GetBinCenter(histo->GetMaximumBin());
      fitfunc_1 = new TF1("fit_1","gaus");
      fitfunc_1->SetRange(maxbin-nRMS*rms,maxbin+nRMS*rms);
      histo->Fit(fitfunc_1,"RQN");
//      cout << "fit1:  mean: " << fitfunc_1->GetParameter(1) << "  rms: " << fitfunc_1->GetParameter(2) << endl;
      fitfunc_2 = new TF1("fit_2","gaus");
      fitfunc_2->SetRange(fitfunc_1->GetParameter(1)-nRMS*fitfunc_1->GetParameter(2),fitfunc_1->GetParameter(1)+nRMS*fitfunc_1->GetParameter(2));
      histo->Fit(fitfunc_2,"RQN");
//      cout << "fit2:  mean: " << fitfunc_2->GetParameter(1) << "  rms: " << fitfunc_2->GetParameter(2) << endl;
      fitfunc_3 = new TF1("fit_3","gaus");
      fitfunc_3->SetRange(fitfunc_2->GetParameter(1)-nRMS*fitfunc_2->GetParameter(2),fitfunc_2->GetParameter(1)+nRMS*fitfunc_2->GetParameter(2));
      histo->Fit(fitfunc_3,"RQN");
//      cout << "fit3:  mean: " << fitfunc_3->GetParameter(1) << "  rms: " << fitfunc_3->GetParameter(2) << endl;
      
      TF1 *fitfunc = 0;
      fitfunc = new TF1(func_title.c_str(),"gaus");
      fitfunc->SetRange(fitfunc_3->GetParameter(1)-nRMS*fitfunc_3->GetParameter(2),fitfunc_3->GetParameter(1)+nRMS*fitfunc_3->GetParameter(2));
      if(i == 1)
        histo->Fit(fitfunc,"RQ");
      else
        histo->Fit(fitfunc,"RQN");
//      cout << "finalFit:  mean: " << fitfunc->GetParameter(1) << "  rms: " << fitfunc->GetParameter(2) << endl;
            
      float par1 = fitfunc->GetParameter(1);
      float Errpar1 = fitfunc->GetParError(1);
      
      delete fitfunc;
      delete fitfunc_1;
      delete fitfunc_2;
      delete fitfunc_3;
      
      values[nRMS-0.5] = pair<float,float>(par1,Errpar1);
    }
  }
  else if(method.compare("doubleGaus") == 0)
  {
    double rms = histo->GetRMS();
    double maxbin =  histo->GetBinCenter(histo->GetMaximumBin());
    
    TF1 *fitfunc = 0;
    fitfunc = new TF1(func_title.c_str(),"gaus(0)+gaus(3)");
    fitfunc->SetParLimits(0,1,9000);
    fitfunc->SetParameter(0,140);
    fitfunc->SetParLimits(1,-0.4,0.6);
    fitfunc->SetParameter(1,-0.02);
    fitfunc->SetParLimits(2,0.02,2);
    fitfunc->SetParameter(2,0.1);
    fitfunc->SetParLimits(3,1,9000);
    fitfunc->SetParameter(3,50);
    fitfunc->SetParLimits(4,0.05,1);
    fitfunc->SetParameter(4,0.1);
    fitfunc->SetParLimits(5,0.02,3);
    fitfunc->SetParameter(5,0.2);
//    fitfunc->SetRange(maxbin-2*rms,0.6);
    fitfunc->SetRange(-0.3,0.8);

    histo->Fit(fitfunc,"RQ");
    
    TF1 *backgr = new TF1((func_title+"_BG").c_str(),"gaus(0)");
    backgr->SetParameter(0,fitfunc->GetParameter(3));
    backgr->SetParameter(1,fitfunc->GetParameter(4));
    backgr->SetParameter(2,fitfunc->GetParameter(5));
    backgr->SetRange(maxbin-2*rms,0.6);
    histo->GetListOfFunctions()->Add(backgr);
      
    int parToUse = -1;
    if( fitfunc->GetParameter(1) < fitfunc->GetParameter(4) )
      parToUse = 1;
    else
      parToUse = 4;
    float par1 = fitfunc->GetParameter(parToUse);
    float Errpar1 = fitfunc->GetParError(parToUse);
    
    values[0.5] = pair<float,float>(par1,0);
    values[1] = pair<float,float>(par1,Errpar1);
    values[1.5] = pair<float,float>(par1,0);
  }
  else if(method.compare("mean") == 0)
  {
    float par1 = histo->GetMean();
    float Errpar1 = histo->GetMeanError();
    values[0.5] = pair<float,float>(par1,0);
    values[1] = pair<float,float>(par1,Errpar1);
    values[1.5] = pair<float,float>(par1,0);
  }
  else cout << "Unknown method in ExpCorrCalculator::GetExpectation :  " << method << endl;
  
  //calculating systematic error
  float syst = 0;
  if (fabs(values[0.5].first-values[1].first) > fabs(values[1].first-values[1.5].first))
    syst = fabs(values[0.5].first-values[1].first);
  else
    syst = fabs(values[1].first-values[1.5].first);
  
  //cout << "Systematic error " << syst << endl;
  return pair<float,pair<float, float> >(values[1].first,pair<float,float>(values[1].second,syst));
}

float* correctMTopErrors(float mTopError, float DElError) // for non unity pull distribution
{
  static float rescaled[2];
  rescaled[0] = mTopError * 1.326; // calculated on March 30th
  rescaled[1] = DElError * 1.733; // calculated on March 30th
  return rescaled;
}

float* correctJESErrors(float DElError, float DEbError) // for non unity pull distribution
{
  static float rescaled[2];
  rescaled[0] = DElError * 1.824; // Calculated on April 4th
  rescaled[1] = DEbError * 1.549; // Calculated on April 4th
  return rescaled;
}

float** calibrateMTop(float mTop, float mTopError, float DEl, float DElError, string calibCurveName) // apply calibrationCurve
{
  static float ** rescaled = new float*[2];
  for(int i=0; i<2; i++)
    rescaled[i] = new float[2];
  
  if(calibCurveName == "MuPlus")
  {
    rescaled[0][0] = ( mTop - 126.3 ) / 0.255; // calculated for 2D-parabola fit
    rescaled[0][1] = fabs( ( ( mTop - 126.3 ) / 0.255 ) - ( ( mTop + mTopError - 126.3 ) / 0.255 ) ); // calculated for 2D-parabola fit
//    rescaled[0][0] = ( mTop - 125.7 ) / 0.261; // Calculated for DEl = 0
//    rescaled[0][1] = fabs( ( ( mTop - 125.7 ) / 0.261 ) - ( ( mTop + mTopError - 125.7 ) / 0.261 ) ); // Calculated for DEl = 0
    rescaled[1][0] = DEl; // no calibration curve calculated yet
    rescaled[1][1] = DElError; // no calibration curve calculated yet
  }
  else if(calibCurveName == "MuMinus")
  {
    rescaled[0][0] = ( mTop - 121.2 ) / 0.289; // calculated for 2D-parabola fit
    rescaled[0][1] = fabs( ( ( mTop - 121.2 ) / 0.289 ) - ( ( mTop + mTopError - 121.2 ) / 0.289 ) ); // calculated for 2D-parabola fit
//    rescaled[0][0] = ( mTop - 120.4 ) / 0.296; // Calculated for DEl = 0
//    rescaled[0][1] = fabs( ( ( mTop - 120.4 ) / 0.296 ) - ( ( mTop + mTopError - 120.4 ) / 0.296 ) ); // Calculated for DEl = 0
    rescaled[1][0] = DEl; // no calibration curve calculated yet
    rescaled[1][1] = DElError; // no calibration curve calculated yet
  }
  else
  {
//    rescaled[0][0] = ( mTop - 123.9 ) / 0.2709; // calculated for 2D-parabola fit
//    rescaled[0][1] = fabs( ( ( mTop - 123.9 ) / 0.2709 ) - ( ( mTop + mTopError - 123.9 ) / 0.2709 ) ); // calculated for 2D-parabola fit
    cout << "No seperate MuPlus MuMinus Calicration Curves   calibCurveName = " << calibCurveName << endl;
    rescaled[0][0] = mTop;
    rescaled[0][1] = mTopError;
    rescaled[1][0] = DEl; // no calibration curve calculated yet
    rescaled[1][1] = DElError; // no calibration curve calculated yet
  }
  return rescaled;
}

float** estimateGoodCombisMVA(vector<float> binning, TH1F* goodTemplate, TH1F* badTemplate, TH1F* totalMVA, string pathforPNG, TFile* f, string title) // estimate the number of good and bad jet combi's as a function of the MVA value
{
  using namespace RooFit;
  
  f->cd();
	string dirname = "PlotsMvaTemplateFit_" + title;
	if(f->Get(dirname.c_str())==0)
		f->mkdir(dirname.c_str());
	f->cd(dirname.c_str());
  
  static float ** result = new float*[10];
  for(int i=0; i<10; i++)
    result[i] = new float[2];
  
  // do the template fit on the maxMVA distribution
  cout << "Performing template fit on MaxMVA distribution" << endl;
	RooRealVar maxMVA("maxMVA","maxMVA",0.,1.);
	
  RooDataHist goodCombiHist("goodCombiHist","goodCombiHist",maxMVA,goodTemplate);
  RooHistPdf goodCombiPDF("goodCombiPDF","goodCombiPDF",maxMVA,goodCombiHist);

  RooDataHist badCombiHist("badCombiHist","badCombiHist",maxMVA,badTemplate);
  RooHistPdf badCombiPDF("badCombiPDF","badCombiPDF",maxMVA,badCombiHist);

  RooRealVar NgoodCombi("NgoodCombi","NgoodCombi",100.,0.,10000.);
  RooRealVar NbadCombi("NbadCombi","NbadCombi",100.,0.,10000.);

  RooAddPdf model("model","model",RooArgList(goodCombiPDF,badCombiPDF),RooArgList(NgoodCombi,NbadCombi));

  RooDataHist dataHist("dataHist","dataHist",maxMVA,totalMVA);
  
  model.fitTo(dataHist, SumW2Error(true), PrintLevel(-3), Verbose(false), Extended(false));
    
  RooPlot* plot = maxMVA.frame();
  dataHist.plotOn(plot);
  model.plotOn(plot);
  model.plotOn(plot, Components(badCombiPDF), LineStyle(kDashed), LineColor(2));
    
  TCanvas* cFit = new TCanvas(("maxMVA_fit"+title).c_str(),("maxMVA_fit_"+title).c_str());
  cFit->cd();
  plot->Draw();
  cFit->SaveAs( (pathforPNG+"maxMVA_fit_"+title+".png").c_str() );
  cFit->Write();
    
  cout << "Full range" << endl;
  cout << "Estimation of NgoodCombi :  " << NgoodCombi.getVal() << " +- " << NgoodCombi.getError() << "      " << NgoodCombi.getErrorHi() << "      " << NgoodCombi.getErrorLo() << endl;
  cout << "Estimation of NbadCombi  :  " << NbadCombi.getVal() << " +- " << NbadCombi.getError() << "      " << NbadCombi.getErrorHi() << "      " << NbadCombi.getErrorLo() << endl;
  cout << "NgoodCombi + NbadCombi : " << NgoodCombi.getVal() + NbadCombi.getVal() << endl;
  
  for(unsigned int j=0; j<binning.size()-1; j++)
  {
    float min = binning[j], max = binning[j+1];
    
    maxMVA.setRange("window",min,max);
    RooRealVar NgoodCombiWindow("NgoodCombiWindow","NgoodCombiWindow",100.,0.,10000.);
    RooRealVar NbadCombiWindow("NbadCombiWindow","NbadCombiWindow",100.,0.,10000.);
    
    RooExtendPdf eGoodCombiPDF("eGoodCombiPDF","eGoodCombiPDF",goodCombiPDF,NgoodCombiWindow,"window");
    RooExtendPdf eBadCombiPDF("eBadCombiPDF","eBadCombiPDF",badCombiPDF,NbadCombiWindow,"window");
    
    RooAddPdf modelWindow("modelWindow","modelWindow",RooArgList(eGoodCombiPDF,eBadCombiPDF));
    
    modelWindow.fitTo(dataHist, SumW2Error(false), PrintLevel(-1), Verbose(false), Extended(true));
    
    RooPlot* plotTmp = maxMVA.frame();
    dataHist.plotOn(plotTmp);
    modelWindow.plotOn(plotTmp);
    modelWindow.plotOn(plotTmp, Components(eBadCombiPDF), LineStyle(kDashed), LineColor(2));
    
    stringstream ss1, ss2;
    ss1 << min;
    ss2 << max;
    string plotTitle = "maxMVA_fit_range_"+ss1.str()+"-"+ss2.str()+"_"+title;
    
    TCanvas* cFitTmp = new TCanvas(plotTitle.c_str(),plotTitle.c_str());
    cFitTmp->cd();
    plotTmp->Draw();
    cFitTmp->SaveAs( (pathforPNG+plotTitle+".png").c_str() );
    cFitTmp->Write();
    
    float SoverSB = NgoodCombiWindow.getVal() / ( NbadCombiWindow.getVal() + NgoodCombiWindow.getVal() );
    float ErrorSoverSB = sqrt( pow(NgoodCombiWindow.getError() * NbadCombiWindow.getVal() ,2) + pow(NbadCombiWindow.getError() * NgoodCombiWindow.getVal() ,2) )
      / pow( NgoodCombiWindow.getVal() + NbadCombiWindow.getVal(), 2);
    
    cout << "Range: " << min << " - " << max << endl;
    cout << "Estimation of NgoodCombi :  " << NgoodCombiWindow.getVal() << " +- " << NgoodCombiWindow.getError() << "  Estimation of NbadCombi  :  " << NbadCombiWindow.getVal() << " +- " << NbadCombiWindow.getError() << endl;
    cout << "NgoodCombi + NbadCombi : " << NgoodCombiWindow.getVal() + NbadCombiWindow.getVal() << endl;
    cout << "NgoodCombi / ( NbadCombi + NgoodCombi ) : " << SoverSB << " +- " << ErrorSoverSB << endl;
    
    int lowBin = goodTemplate->GetXaxis()->FindBin(min);
    int highBin = goodTemplate->GetXaxis()->FindBin(max-0.001);
    
    cout << "lowBin: " << lowBin << " value: " << goodTemplate->GetXaxis()->GetBinLowEdge(lowBin) << "   highBin: " << highBin << " value: " << goodTemplate->GetXaxis()->GetBinUpEdge(highBin) << endl;
    cout << "Expected Ngood / ( Nbad + Ngood ) : " << goodTemplate->Integral(lowBin, highBin) / ( badTemplate->Integral(lowBin, highBin) + goodTemplate->Integral(lowBin, highBin) ) << endl << endl;
    
    result[j][0] = SoverSB;
    result[j][1] = ErrorSoverSB;
    
    delete plotTmp;
    delete cFitTmp;
  }
  
  delete plot;
  delete cFit;
  
  return result;
}

