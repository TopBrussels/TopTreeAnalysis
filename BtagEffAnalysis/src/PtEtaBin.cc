#include "../interface/PtEtaBin.h"

//constructor
PtEtaBin::PtEtaBin(int debug, int nVar1, int nVar0, int nBdisc, double ptbinlow, double etabinlow, double ptbinup, double etabinup, bool varBinSize, bool doShift){
  
  debug_=debug;
  nVar1_=nVar1;
  nVar0_=nVar0;
  nBdisc_=nBdisc;
  ptbinlow_=ptbinlow;
  etabinlow_=etabinlow;  
  ptbinup_=ptbinup;
  etabinup_=etabinup;  
  varBinSize_=varBinSize;
  lowerBin_=0;
  upperBin_=-1;
  centralLowBin_=0;  
  centralUpBin_=0;  
  doShift_=doShift;
  
  TString strGenericName="nDisc_"; strGenericName+=nBdisc_; strGenericName+="_ptbinlow_"; strGenericName+= (int) ptbinlow_; strGenericName+="_etabinlow_"; strGenericName+=(int) round(etabinlow_*10); strGenericName+="_";

  genericName = new TString();
  genericName->Append(strGenericName);
  
  //if(debug_>1) cout << "+----> genericName " << *genericName << endl;
  //  cout << "+----> genericName " << *genericName << endl;

  //initialize all histograms to NULL, all member histograms should be in here. Potential memory leaks might be due to forgotten histograms and/or histograms which were created in member functions 'on-the-fly' not being member histograms.

  // the 1DX and 1DY plots of TH2(Bkg,Sng,Data)_(Var12,ControlVar12) are not defined null!

  TH1Sng_chisq = NULL;
  TH1Bkg_chisq = NULL;
  TH1Data_chisq = NULL;
  TH1Sng_chisqControl = NULL;
  TH1Bkg_chisqControl = NULL;
  TH1Data_chisqControl = NULL;

  TH2Sng = NULL;
  TH2Bkg = NULL;
  TH2Data = NULL;   
  TH2Sng_ProfileX = NULL;
  TH2Bkg_ProfileX = NULL;
  TH2Data_ProfileX = NULL; 
  TH2Sng_Left = NULL;
  TH2Bkg_Left = NULL;
  TH2Data_Left = NULL;
  TH2Sng_Right = NULL;
  TH2Bkg_Right = NULL;
  TH2Data_Right = NULL;
  TH2SngVar12_Left = NULL;
  TH2BkgVar12_Left = NULL;
  TH2DataVar12_Left = NULL;
  TH2SngVar12_Right = NULL;
  TH2BkgVar12_Right = NULL;
  TH2DataVar12_Right = NULL;
  TH1Sng_Var0 = NULL;
  TH1Bkg_Var0 = NULL;
  TH1Bkg_R_Var0 = NULL;
  TH1Bkg_W_Var0 = NULL;
    TH1Data_Var0 = NULL;
    TH1Data_Var0_XS = NULL;
  TH1SoverSB_Var0 = NULL;
  TH1SoverSB_ControlVar = NULL;
  TH1Sng_BtagAll = NULL;
  TH1Bkg_BtagAll = NULL;
  TH1Data_BtagAll = NULL;  
  TH1Sng_BtagEffAll = NULL;
  TH1Bkg_BtagEffAll = NULL;
  TH1Data_BtagEffAll = NULL;  
  TH2Sng_Var12 = NULL;
  TH2Bkg_Var12 = NULL;
  TH2Bkg_W_Var12 = NULL;
  TH2Bkg_R_Var12 = NULL;
  TH2Data_Var12 = NULL;
  TH1Sng_ControlVar = NULL;
  TH1Bkg_ControlVar = NULL;
  TH1Bkg_W_ControlVar = NULL;
  TH1Bkg_R_ControlVar = NULL;
  TH1Data_ControlVar = NULL;   
  TH1Sng_ControlVarReweigh = NULL;
  TH1Bkg_ControlVarReweigh = NULL;
  TH1Bkg_W_ControlVarReweigh = NULL;
  TH1Bkg_R_ControlVarReweigh = NULL;
  TH1Data_ControlVarReweigh = NULL; 
  TH2Data_SignalControlVar12 = NULL;
  TH2Data_SignalControlVar12_1DX = NULL;
  TH2Data_SignalControlVar12_1DY = NULL;
  
  TH2Sng_ControlVar12 = NULL;
  TH2Bkg_ControlVar12 = NULL;
  TH2Bkg_W_ControlVar12 = NULL;
  TH2Bkg_R_ControlVar12 = NULL;
  TH2Data_ControlVar12 = NULL;
  TH1Sng_BtagEffMC = NULL;
  TH1Bkg_BtagEffMC = NULL;
  TH1Data_BtagEffMC = NULL;
  TH2Sng_Left1DX = NULL;
  TH2Bkg_Left1DX = NULL;
  TH2Data_Left1DX = NULL;
  TH2Sng_Right1DX = NULL;
  TH2Bkg_Right1DX = NULL;
  TH2Data_Right1DX = NULL;
  TH2SngVar12_Left1DY = NULL;
  TH2BkgVar12_Left1DY = NULL;
  TH2DataVar12_Left1DY = NULL;
  TH2SngVar12_Right1DY = NULL;
  TH2BkgVar12_Right1DY = NULL;
  TH2DataVar12_Right1DY = NULL;
  TH2SngVar12_Right1DYReweigh = NULL;
  TH2BkgVar12_Right1DYReweigh = NULL;
  TH2DataVar12_Right1DYReweigh = NULL;
  TH2SngVar12_Right1DYReweighRatio = NULL;
  TH2BkgVar12_Right1DYReweighRatio = NULL;
  TH2DataVar12_Right1DYReweighRatio = NULL;
  TH2Data_RightLeft1DX = NULL;
  TH2Data_LeftRight1DX = NULL;
  TFData_LeftRight1DXFit = NULL;
  TH2DataVar12_LeftRight = NULL;
  TH2Data_Left1DY = NULL;
  TH2Sng_Left1DY = NULL;
  TH2Bkg_Left1DY = NULL;
  TH2Data_Right1DY = NULL;
  TH2Sng_Right1DY = NULL;
  TH2Bkg_Right1DY = NULL;

  TH2Bkg_Right1DYReweigh = NULL;
  TH2Sng_Right1DYReweigh = NULL;
  TH2Data_Right1DYReweigh = NULL;
  TH2Bkg_Left1DYReweigh = NULL;
  TH2Sng_Left1DYReweigh = NULL;
  TH2Data_Left1DYReweigh = NULL;

  TH2Bkg_Right1DYReweighRatio = NULL;
  TH2Sng_Right1DYReweighRatio = NULL;
  TH2Data_Right1DYReweighRatio = NULL;
  TH2Bkg_Left1DYReweighRatio = NULL;
  TH2Sng_Left1DYReweighRatio = NULL;
  TH2Data_Left1DYReweighRatio = NULL;

  TH2Bkg_Right1DY_LeftRightRatio = NULL;
  TH2Sng_Right1DY_LeftRightRatio = NULL;
  TH2Data_Right1DY_LeftRightRatio = NULL;

  TH2Bkg_Right1DYReweigh_LeftRightRatio = NULL;
  TH2Sng_Right1DYReweigh_LeftRightRatio = NULL;
  TH2Data_Right1DYReweigh_LeftRightRatio = NULL;
  TH2Bkg_Left1DYReweigh_LeftRightRatio = NULL;
  TH2Sng_Left1DYReweigh_LeftRightRatio = NULL;
  TH2Data_Left1DYReweigh_LeftRightRatio = NULL;

  TH1Data_BtagMeasured = NULL;
  TH1Data_BtagMCMeasured = NULL;
  TH1Data_BtagEffMeasured = NULL;
  TH1Data_BtagEffMCMeasured = NULL;
  TH1Data_BtagEffMeasuredDiff = NULL;
  TH1Data_BtagEffMCMeasuredDiff = NULL;
  TH1Data_BtagMeasuredLR = NULL;
  TH1Data_BtagMCMeasuredLR = NULL;
  TH1Data_BtagEffMeasuredLR = NULL;
  TH1Data_BtagEffMCMeasuredLR = NULL;
  TH1Data_BtagEffMeasuredLRDiff = NULL;
  TH1Data_BtagEffMCMeasuredLRDiff = NULL;
  TH1Data_BtagMeasuredRR = NULL;
  TH1Data_BtagMCMeasuredRR = NULL;
  TH1Data_BtagEffMeasuredRR = NULL;
  TH1Data_BtagEffMCMeasuredRR = NULL;
  TH1Data_BtagEffMeasuredRRDiff = NULL;
  TH1Data_BtagEffMCMeasuredRRDiff = NULL;

  //cout << "done constructot PtEtaBin"<< endl;
}

//destructor
PtEtaBin::~PtEtaBin(){

  //cout << "start deleting bin" << endl;  
  delete TH1Sng_chisq;
  delete TH1Bkg_chisq;
  delete TH1Data_chisq;
  delete TH1Sng_chisqControl;
  delete TH1Bkg_chisqControl;
  delete TH1Data_chisqControl;

  delete genericName;
  delete TH2Sng;
  delete TH2Bkg;
  delete TH2Data;    
  delete TH2Sng_ProfileX;
  delete TH2Bkg_ProfileX;
  delete TH2Data_ProfileX;  
  delete TH2Sng_Left;
  delete TH2Bkg_Left;
  delete TH2Data_Left;
  delete TH2Sng_Right;
  delete TH2Bkg_Right;
  delete TH2Data_Right;
  delete TH2SngVar12_Left;
  delete TH2BkgVar12_Left;
  delete TH2DataVar12_Left;
  delete TH2SngVar12_Right;
  delete TH2BkgVar12_Right;
  delete TH2DataVar12_Right;
  delete TH1Sng_Var0;
  delete TH1Bkg_Var0;
  delete TH1Bkg_W_Var0;
  delete TH1Bkg_R_Var0;
    delete TH1Data_Var0;
    delete TH1Data_Var0_XS;
  delete TH1SoverSB_Var0;
  delete TH1SoverSB_ControlVar;
  delete TH1Sng_BtagAll;
  delete TH1Bkg_BtagAll;
  delete TH1Data_BtagAll; 
  delete TH1Sng_BtagEffAll;
  delete TH1Bkg_BtagEffAll;
  delete TH1Data_BtagEffAll; 
  delete TH2Sng_Var12;
  delete TH2Bkg_Var12;
  delete TH2Bkg_W_Var12;
  delete TH2Bkg_R_Var12;
  delete TH2Data_Var12;
  delete TH1Sng_ControlVar;
  delete TH1Bkg_ControlVar;
  delete TH1Bkg_W_ControlVar;
  delete TH1Bkg_R_ControlVar;
  delete TH1Data_ControlVar;  
  delete TH1Sng_ControlVarReweigh;
  delete TH1Bkg_ControlVarReweigh;
  delete TH1Bkg_W_ControlVarReweigh;
  delete TH1Bkg_R_ControlVarReweigh;
  delete TH1Data_ControlVarReweigh; 
  delete TH2Sng_ControlVar12;
  delete TH2Bkg_ControlVar12;
  delete TH2Bkg_W_ControlVar12;
  delete TH2Bkg_R_ControlVar12;
  delete TH2Data_ControlVar12;
  delete TH2Data_SignalControlVar12;
  delete TH2Data_SignalControlVar12_1DX;
  delete TH2Data_SignalControlVar12_1DY;
  
  delete TH1Sng_BtagEffMC;
  delete TH1Bkg_BtagEffMC;
  delete TH1Data_BtagEffMC;
  delete TH2Sng_Left1DX;
  delete TH2Bkg_Left1DX;
  delete TH2Data_Left1DX;
  delete TH2Sng_Right1DX;
  delete TH2Bkg_Right1DX;
  delete TH2Data_Right1DX;
  delete TH2SngVar12_Left1DY;
  delete TH2BkgVar12_Left1DY;
  delete TH2DataVar12_Left1DY;
  delete TH2SngVar12_Right1DY;
  delete TH2BkgVar12_Right1DY;
  delete TH2DataVar12_Right1DY;
  delete TH2SngVar12_Right1DYReweigh;
  delete TH2BkgVar12_Right1DYReweigh;
  delete TH2DataVar12_Right1DYReweigh;
  delete TH2SngVar12_Right1DYReweighRatio;
  delete TH2BkgVar12_Right1DYReweighRatio;
  delete TH2DataVar12_Right1DYReweighRatio;
  delete TH2Data_RightLeft1DX;
  delete TH2Data_LeftRight1DX;
  delete TFData_LeftRight1DXFit;
  delete TH2DataVar12_LeftRight;
  delete TH2Data_Left1DY;
  delete TH2Sng_Left1DY;
  delete TH2Bkg_Left1DY;
  delete TH2Data_Right1DY;
  delete TH2Sng_Right1DY;
  delete TH2Bkg_Right1DY;
  //delete TH2Data_Left1DYReweigh;

  delete TH2Bkg_Right1DYReweigh;
  delete TH2Sng_Right1DYReweigh;
  delete TH2Data_Right1DYReweigh;
  delete TH2Bkg_Left1DYReweigh;
  delete TH2Sng_Left1DYReweigh;
  delete TH2Data_Left1DYReweigh;

  delete TH2Bkg_Right1DYReweighRatio;
  delete TH2Sng_Right1DYReweighRatio;
  delete TH2Data_Right1DYReweighRatio;
  delete TH2Bkg_Left1DYReweighRatio;
  delete TH2Sng_Left1DYReweighRatio;
  delete TH2Data_Left1DYReweighRatio;

  delete TH2Bkg_Right1DY_LeftRightRatio;
  delete TH2Sng_Right1DY_LeftRightRatio;
  delete TH2Data_Right1DY_LeftRightRatio;

  delete TH2Bkg_Right1DYReweigh_LeftRightRatio;
  delete TH2Sng_Right1DYReweigh_LeftRightRatio;
  delete TH2Data_Right1DYReweigh_LeftRightRatio;
  delete TH2Bkg_Left1DYReweigh_LeftRightRatio;
  delete TH2Sng_Left1DYReweigh_LeftRightRatio;
  delete TH2Data_Left1DYReweigh_LeftRightRatio;

  delete TH1Data_BtagMeasured;
  delete TH1Data_BtagMCMeasured;
  delete TH1Data_BtagEffMeasured;
  delete TH1Data_BtagEffMCMeasured;
  delete TH1Data_BtagEffMeasuredDiff;
  delete TH1Data_BtagEffMCMeasuredDiff;
  delete TH1Data_BtagMeasuredLR;
  delete TH1Data_BtagMCMeasuredLR;
  delete TH1Data_BtagEffMeasuredLR;
  delete TH1Data_BtagEffMCMeasuredLR;
  delete TH1Data_BtagEffMeasuredLRDiff;
  delete TH1Data_BtagEffMCMeasuredLRDiff;
  delete TH1Data_BtagMeasuredRR;
  delete TH1Data_BtagMCMeasuredRR;
  delete TH1Data_BtagEffMeasuredRR;
  delete TH1Data_BtagEffMCMeasuredRR;
  delete TH1Data_BtagEffMeasuredRRDiff;
  delete TH1Data_BtagEffMCMeasuredRRDiff;
	
	delete TH1Data_Pt_Left;
	delete TH1Sng_Pt_Left;
	delete TH1Data_Pt_Right;
	delete TH1Data_Pt_Control;
	
	
	for (std::map<TString,TH1D*>::const_iterator it=histo1D.begin(); it != histo1D.end(); ++it) {
		
	//cout << "blaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaat " << it->first << endl;
		//it->second->Delete();
		delete histo1D[it->first];
	}
	//exit(1);
}

void PtEtaBin::SetVarBins(vector<float> rangesbTag) {

	for (unsigned int i=0; i<rangesbTag.size(); i++)
		rangesbTag_[i]=rangesbTag[i];
	
}

void PtEtaBin::DefineSignalSamplePlots(int nBdiscrAlgos,int nBinsVar1, double lowRangVar1, double upRangeVar1, int nBinsVar2, double lowRangVar2, double upRangeVar2, int nBinsVar1SC, double lowRangVar1SC, double upRangeVar1SC, int nBinsVar2SC, double lowRangVar2SC, double upRangeVar2SC, int nBinsVar0, double lowRangeVar0, double upRangeVar0, int nBinsBtag[], double lowRangeBtag[], double upRangeBtag[]){

  //if(debug_>2) cout << "PtEtaBin::DefineSignalSamplePlots::START" << endl;
  
  //+if(debug_>2) cout << "PtEtaBin::DefineSignalSamplePlots " << *genericName << endl;

  GiveName(&titleSng_); titleSng_+="TH2Sng"; 
  GiveName(&titleBkg_); titleBkg_+="TH2Bkg";
  GiveName(&titleData_); titleData_+="TH2Data"; 

  //if(debug_>2) cout << "PtEtaBin::DefineSignalSamplePlots " << titleSng_ << endl;
  
  if(!varBinSize_){ /*would be needed if making a 2D histo of pt vs btag
    TH2Sng = new TH2D(titleSng_,titleSng_,nBinsVar1,lowRangVar1,upRangeVar1,nBinsBtag[nBdiscrAlgos],lowRangeBtag[nBdiscrAlgos],upRangeBtag[nBdiscrAlgos]);
    TH2Bkg = new TH2D(titleBkg_,titleBkg_,nBinsVar1,lowRangVar1,upRangeVar1,nBinsBtag[nBdiscrAlgos],lowRangeBtag[nBdiscrAlgos],upRangeBtag[nBdiscrAlgos]);
    TH2Data = new TH2D(titleData_,titleData_,nBinsVar1,lowRangVar1,upRangeVar1,nBinsBtag[nBdiscrAlgos],lowRangeBtag[nBdiscrAlgos],upRangeBtag[nBdiscrAlgos]);*/
    TH2Sng = new TH2D(titleSng_,titleSng_,nBinsVar0,lowRangeVar0,upRangeVar0,nBinsBtag[nBdiscrAlgos],lowRangeBtag[nBdiscrAlgos],upRangeBtag[nBdiscrAlgos]);
    TH2Bkg = new TH2D(titleBkg_,titleBkg_,nBinsVar0,lowRangeVar0,upRangeVar0,nBinsBtag[nBdiscrAlgos],lowRangeBtag[nBdiscrAlgos],upRangeBtag[nBdiscrAlgos]);
    TH2Data = new TH2D(titleData_,titleData_,nBinsVar0,lowRangeVar0,upRangeVar0,nBinsBtag[nBdiscrAlgos],lowRangeBtag[nBdiscrAlgos],upRangeBtag[nBdiscrAlgos]);
  }

  GiveName(&titleSng_Left_); titleSng_Left_+="TH2Sng_Left"; 
  GiveName(&titleBkg_Left_); titleBkg_Left_+="TH2Bkg_Left";
  GiveName(&titleData_Left_); titleData_Left_+="TH2Data_Left"; 
 
  if(!varBinSize_){
    TH2Sng_Left = new TH2D(titleSng_Left_,titleSng_Left_,nBinsVar1,lowRangVar1,upRangeVar1,nBinsBtag[nBdiscrAlgos],lowRangeBtag[nBdiscrAlgos],upRangeBtag[nBdiscrAlgos]);
    TH2Bkg_Left = new TH2D(titleBkg_Left_,titleBkg_Left_,nBinsVar1,lowRangVar1,upRangeVar1,nBinsBtag[nBdiscrAlgos],lowRangeBtag[nBdiscrAlgos],upRangeBtag[nBdiscrAlgos]);
    TH2Data_Left = new TH2D(titleData_Left_,titleData_Left_,nBinsVar1,lowRangVar1,upRangeVar1,nBinsBtag[nBdiscrAlgos],lowRangeBtag[nBdiscrAlgos],upRangeBtag[nBdiscrAlgos]);
  }
  
  //int nvarBinsVar1=20;

  if(varBinSize_){
    //This is dirty and annoying but I have the array of binranges defined here
  
    /*would be needed if making a 2D histo of pt vs btag
    TH2Sng = new TH2D(titleSng_,titleSng_,nBinsVar1,lowRangVar1,upRangeVar1,nBinsBtag[nBdiscrAlgos],rangesbTag_);
    TH2Bkg = new TH2D(titleBkg_,titleBkg_,nBinsVar1,lowRangVar1,upRangeVar1,nBinsBtag[nBdiscrAlgos],rangesbTag_);
    TH2Data = new TH2D(titleData_,titleData_,nBinsVar1,lowRangVar1,upRangeVar1,nBinsBtag[nBdiscrAlgos],rangesbTag_);  */
    TH2Sng = new TH2D(titleSng_,titleSng_,nBinsVar0,lowRangeVar0,upRangeVar0,nBinsBtag[nBdiscrAlgos],rangesbTag_);
    TH2Bkg = new TH2D(titleBkg_,titleBkg_,nBinsVar0,lowRangeVar0,upRangeVar0,nBinsBtag[nBdiscrAlgos],rangesbTag_);
    TH2Data = new TH2D(titleData_,titleData_,nBinsVar0,lowRangeVar0,upRangeVar0,nBinsBtag[nBdiscrAlgos],rangesbTag_);  

    /*TH2Sng_Left = new TH2D(titleSng_Left_,titleSng_Left_,nvarBinsVar1,rangesVar1_,nBinsBtag,rangesbTag_);
    TH2Bkg_Left = new TH2D(titleBkg_Left_,titleBkg_Left_,nvarBinsVar1,rangesVar1_,nBinsBtag,rangesbTag_);
    TH2Data_Left = new TH2D(titleData_Left_,titleData_Left_,nvarBinsVar1,rangesVar1_,nBinsBtag,rangesbTag_);*/
    TH2Sng_Left = new TH2D(titleSng_Left_,titleSng_Left_,nBinsVar1,lowRangVar1,upRangeVar1,nBinsBtag[nBdiscrAlgos],rangesbTag_);
    TH2Bkg_Left = new TH2D(titleBkg_Left_,titleBkg_Left_,nBinsVar1,lowRangVar1,upRangeVar1,nBinsBtag[nBdiscrAlgos],rangesbTag_);
    TH2Data_Left = new TH2D(titleData_Left_,titleData_Left_,nBinsVar1,lowRangVar1,upRangeVar1,nBinsBtag[nBdiscrAlgos],rangesbTag_);
  } 

  GiveName(&titleSng_Right_); titleSng_Right_+="TH2Sng_Right"; 
  GiveName(&titleBkg_Right_); titleBkg_Right_+="TH2Bkg_Right"; 
  GiveName(&titleData_Right_); titleData_Right_+="TH2Data_Right";
    
  if(!varBinSize_){
    TH2Sng_Right = new TH2D(titleSng_Right_,titleSng_Right_,nBinsVar1,lowRangVar1,upRangeVar1,nBinsBtag[nBdiscrAlgos],lowRangeBtag[nBdiscrAlgos],upRangeBtag[nBdiscrAlgos]);
    TH2Bkg_Right = new TH2D(titleBkg_Right_,titleBkg_Right_,nBinsVar1,lowRangVar1,upRangeVar1,nBinsBtag[nBdiscrAlgos],lowRangeBtag[nBdiscrAlgos],upRangeBtag[nBdiscrAlgos]);
    TH2Data_Right = new TH2D(titleData_Right_,titleData_Right_,nBinsVar1,lowRangVar1,upRangeVar1,nBinsBtag[nBdiscrAlgos],lowRangeBtag[nBdiscrAlgos],upRangeBtag[nBdiscrAlgos]);
  }

  if(varBinSize_){
    /*TH2Sng_Right = new TH2D(titleSng_Right_,titleSng_Right_,nvarBinsVar1,rangesVar1_,nBinsBtag,rangesbTag_);
    TH2Bkg_Right = new TH2D(titleBkg_Right_,titleBkg_Right_,nvarBinsVar1,rangesVar1_,nBinsBtag,rangesbTag_);
    TH2Data_Right = new TH2D(titleData_Right_,titleData_Right_,nvarBinsVar1,rangesVar1_,nBinsBtag,rangesbTag_);*/
    TH2Sng_Right = new TH2D(titleSng_Right_,titleSng_Right_,nBinsVar1,lowRangVar1,upRangeVar1,nBinsBtag[nBdiscrAlgos],rangesbTag_);
    TH2Bkg_Right = new TH2D(titleBkg_Right_,titleBkg_Right_,nBinsVar1,lowRangVar1,upRangeVar1,nBinsBtag[nBdiscrAlgos],rangesbTag_);
    TH2Data_Right = new TH2D(titleData_Right_,titleData_Right_,nBinsVar1,lowRangVar1,upRangeVar1,nBinsBtag[nBdiscrAlgos],rangesbTag_);
  }

  //add also 2D histograms containing the pt and eta

  GiveName(&titleSngVar12_Left_); titleSngVar12_Left_+="TH2SngVar12_Left"; 
  GiveName(&titleBkgVar12_Left_); titleBkgVar12_Left_+="TH2BkgVar12_Left";
  GiveName(&titleDataVar12_Left_); titleDataVar12_Left_+="TH2DataVar12_Left"; 
  
  TH2SngVar12_Left = new TH2D(titleSngVar12_Left_,titleSngVar12_Left_,nBinsVar1,lowRangVar1,upRangeVar1,nBinsVar2,lowRangVar2,upRangeVar2);
  TH2BkgVar12_Left = new TH2D(titleBkgVar12_Left_,titleBkgVar12_Left_,nBinsVar1,lowRangVar1,upRangeVar1,nBinsVar2,lowRangVar2,upRangeVar2);
  TH2DataVar12_Left = new TH2D(titleDataVar12_Left_,titleDataVar12_Left_,nBinsVar1,lowRangVar1,upRangeVar1,nBinsVar2,lowRangVar2,upRangeVar2);
  
  GiveName(&titleSngVar12_Right_); titleSngVar12_Right_+="TH2SngVar12_Right"; 
  GiveName(&titleBkgVar12_Right_); titleBkgVar12_Right_+="TH2BkgVar12_Right"; 
  GiveName(&titleDataVar12_Right_); titleDataVar12_Right_+="TH2DataVar12_Right";
    
  TH2SngVar12_Right = new TH2D(titleSngVar12_Right_,titleSngVar12_Right_,nBinsVar1,lowRangVar1,upRangeVar1,nBinsVar2,lowRangVar2,upRangeVar2);
  TH2BkgVar12_Right = new TH2D(titleBkgVar12_Right_,titleBkgVar12_Right_,nBinsVar1,lowRangVar1,upRangeVar1,nBinsVar2,lowRangVar2,upRangeVar2);
  TH2DataVar12_Right = new TH2D(titleDataVar12_Right_,titleDataVar12_Right_,nBinsVar1,lowRangVar1,upRangeVar1,nBinsVar2,lowRangVar2,upRangeVar2);

  //define 1D histograms
  TH1Sng_chisq = NULL;
  TH1Bkg_chisq = NULL;
  TH1Data_chisq = NULL;
  TH1Sng_chisqControl = NULL;
  TH1Bkg_chisqControl = NULL;
  TH1Data_chisqControl = NULL;

  GiveName(&titleSng_chisq_);  titleSng_chisq_+="TH1Sng_chisq";
  GiveName(&titleBkg_chisq_);  titleBkg_chisq_+="TH1Bkg_chisq";
  GiveName(&titleData_chisq_);  titleData_chisq_+="TH1Data_chisq";
  GiveName(&titleSng_chisqControl_);  titleSng_chisqControl_+="TH1Sng_chisqControl";
  GiveName(&titleBkg_chisqControl_);  titleBkg_chisqControl_+="TH1Bkg_chisqControl";
  GiveName(&titleData_chisqControl_);  titleData_chisqControl_+="TH1Data_chisqControl";
  TH1Sng_chisq         = new TH1D(titleSng_chisq_,titleSng_chisq_,2000,0,100);
  TH1Bkg_chisq         = new TH1D(titleBkg_chisq_,titleBkg_chisq_,2000,0,100);;
  TH1Data_chisq        = new TH1D(titleData_chisq_,titleData_chisq_,2000,0,100);;
  TH1Sng_chisqControl  = new TH1D(titleSng_chisqControl_,titleSng_chisqControl_,2000,0,100);;
  TH1Bkg_chisqControl  = new TH1D(titleBkg_chisqControl_,titleBkg_chisqControl_,2000,0,100);;
  TH1Data_chisqControl = new TH1D(titleData_chisqControl_,titleData_chisqControl_,2000,0,100);;


  GiveName(&titleSng_Var0_); titleSng_Var0_+="TH1Sng_Var0"; 
  GiveName(&titleBkg_Var0_); titleBkg_Var0_+="TH1Bkg_Var0"; 
  GiveName(&titleBkg_W_Var0_); titleBkg_W_Var0_+="TH1Bkg_W_Var0"; 
  GiveName(&titleBkg_R_Var0_); titleBkg_R_Var0_+="TH1Bkg_R_Var0"; 
    GiveName(&titleData_Var0_); titleData_Var0_+="TH1Data_Var0"; 
    GiveName(&titleData_Var0_XS_); titleData_Var0_XS_+="TH1Data_Var0_XS"; 
  
  TH1Sng_Var0 = new TH1D(titleSng_Var0_,titleSng_Var0_,nBinsVar0,lowRangeVar0,upRangeVar0);
  TH1Bkg_Var0 = new TH1D(titleBkg_Var0_,titleBkg_Var0_,nBinsVar0,lowRangeVar0,upRangeVar0);
  TH1Bkg_W_Var0 = new TH1D(titleBkg_W_Var0_,titleBkg_W_Var0_,nBinsVar0,lowRangeVar0,upRangeVar0);
  TH1Bkg_R_Var0 = new TH1D(titleBkg_R_Var0_,titleBkg_R_Var0_,nBinsVar0,lowRangeVar0,upRangeVar0);
    TH1Data_Var0 = new TH1D(titleData_Var0_,titleData_Var0_,nBinsVar0,lowRangeVar0,upRangeVar0);
    TH1Data_Var0_XS = new TH1D(titleData_Var0_XS_,titleData_Var0_XS_,nBinsVar0,lowRangeVar0,upRangeVar0);
		
  GiveName(&titleSoverSB_Var0_); titleSoverSB_Var0_+="TH1SoverSB_Var0"; 
  TH1SoverSB_Var0 = new TH1D(titleSoverSB_Var0_,titleSoverSB_Var0_,nBinsVar0,lowRangeVar0,upRangeVar0);
   
  /*  GiveName(&titleSng_Var1_); titleSng_Var1_+="TH1Sng_Var1"; //obsolete
  GiveName(&titleBkg_Var1_); titleBkg_Var1_+="TH1Bkg_Var1"; 
  GiveName(&titleData_Var1_); titleData_Var1_+="TH1Data_Var1"; 
  
  TH1Sng_Var1 = new TH1D(titleSng_Var1_,titleSng_Var1_,nBinsVar1,lowRangVar1,upRangeVar1);
  TH1Bkg_Var1 = new TH1D(titleBkg_Var1_,titleBkg_Var1_,nBinsVar1,lowRangVar1,upRangeVar1);
  TH1Data_Var1 = new TH1D(titleData_Var1_,titleData_Var1_,nBinsVar1,lowRangVar1,upRangeVar1);*/
  
  GiveName(&titleSng_Var12_); titleSng_Var12_+="TH2Sng_Var12"; 
  GiveName(&titleBkg_Var12_); titleBkg_Var12_+="TH2Bkg_Var12"; 
  GiveName(&titleBkg_W_Var12_); titleBkg_W_Var12_+="TH2Bkg_W_Var12"; 
  GiveName(&titleBkg_R_Var12_); titleBkg_R_Var12_+="TH2Bkg_R_Var12"; 
  GiveName(&titleData_Var12_); titleData_Var12_+="TH2Data_Var12"; 
  
  /*TH2Sng_Var12 = new TH2D(titleSng_Var12_,titleSng_Var12_,nBinsVar1SC,lowRangVar1SC,upRangeVar1SC,nBinsVar2SC,lowRangVar2SC,upRangeVar2SC);
  TH2Bkg_Var12 = new TH2D(titleBkg_Var12_,titleBkg_Var12_,nBinsVar1SC,lowRangVar1SC,upRangeVar1SC,nBinsVar2SC,lowRangVar2SC,upRangeVar2SC);
  TH2Bkg_W_Var12 = new TH2D(titleBkg_W_Var12_,titleBkg_W_Var12_,nBinsVar1SC,lowRangVar1SC,upRangeVar1SC,nBinsVar2SC,lowRangVar2SC,upRangeVar2SC);
  TH2Bkg_R_Var12 = new TH2D(titleBkg_R_Var12_,titleBkg_R_Var12_,nBinsVar1SC,lowRangVar1SC,upRangeVar1SC,nBinsVar2SC,lowRangVar2SC,upRangeVar2SC);
  TH2Data_Var12 = new TH2D(titleData_Var12_,titleData_Var12_,nBinsVar1SC,lowRangVar1SC,upRangeVar1SC,nBinsVar2SC,lowRangVar2SC,upRangeVar2SC);
  */

  
  //  double rangesPtSC[21] = {30,33.9517,38.0413,42.3827,47.1889,52.4274,58.0206,63.9226,70.6819,77.5725,85.4028,94.2801,104.116,115.42,128.936,145.192,165.937,194.182,237.506,316.667,500};
  
  double rangesPtSC[19] = {30,33.9517,38.0413,42.3827,47.1889,52.4274,58.0206,63.9226,70.6819,77.5725,85.4028,94.2801,104.116,115.42,128.936,145.192,165.937,250,9999.};

  //double rangesPtSC[51] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400,410,420,430,440,450,460,470,480,490,500};

  //double rangesEtaSC[11] = {0,0.24,0.48,0.72,0.96,1.20,1.44,1.68,1.92,2.16,2.4};
  double rangesEtaSC[21] = {0,0.12,0.24,0.36,0.48,0.60,0.72,0.84,0.96,1.08,1.20,1.32,1.44,1.56,1.68,1.80,1.92,2.04,2.16,2.28,2.4};
  //double rangesEtaSC[21] = {0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1};
  //double rangesEtaSC[21] = {-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
  //double rangesEtaSC[21] = {0,0.175,0.175*2,0.175*3,0.175*4,0.175*5,0.175*6,0.175*7,0.175*8,0.175*9,0.175*10,0.175*11,0.175*12,0.175*13,0.175*14,0.175*15,0.175*16,0.175*17,0.175*18,0.175*19,0.175*20};
  
  TH2Sng_Var12 = new TH2D(titleSng_Var12_,titleSng_Var12_,nBinsVar1SC,rangesPtSC,nBinsVar2SC,rangesEtaSC);
  TH2Bkg_Var12 = new TH2D(titleBkg_Var12_,titleBkg_Var12_,nBinsVar1SC,rangesPtSC,nBinsVar2SC,rangesEtaSC);
  TH2Bkg_W_Var12 = new TH2D(titleBkg_W_Var12_,titleBkg_W_Var12_,nBinsVar1SC,rangesPtSC,nBinsVar2SC,rangesEtaSC);
  TH2Bkg_R_Var12 = new TH2D(titleBkg_R_Var12_,titleBkg_R_Var12_,nBinsVar1SC,rangesPtSC,nBinsVar2SC,rangesEtaSC);
  TH2Data_Var12 = new TH2D(titleData_Var12_,titleData_Var12_,nBinsVar1SC,rangesPtSC,nBinsVar2SC,rangesEtaSC);
  
  
  GiveName(&titleSng_BtagAll_); titleSng_BtagAll_+="TH1Sng_BtagAll"; 
  GiveName(&titleBkg_BtagAll_); titleBkg_BtagAll_+="TH1Bkg_BtagAll"; 
  GiveName(&titleData_BtagAll_); titleData_BtagAll_+="TH1Data_BtagAll"; 
  
  if(!varBinSize_){
    TH1Sng_BtagAll = new TH1D(titleSng_BtagAll_,titleSng_BtagAll_,nBinsBtag[nBdiscrAlgos],lowRangeBtag[nBdiscrAlgos],upRangeBtag[nBdiscrAlgos]);
    TH1Bkg_BtagAll = new TH1D(titleBkg_BtagAll_,titleBkg_BtagAll_,nBinsBtag[nBdiscrAlgos],lowRangeBtag[nBdiscrAlgos],upRangeBtag[nBdiscrAlgos]);
    TH1Data_BtagAll = new TH1D(titleData_BtagAll_,titleData_BtagAll_,nBinsBtag[nBdiscrAlgos],lowRangeBtag[nBdiscrAlgos],upRangeBtag[nBdiscrAlgos]);
  }
  
  if(varBinSize_){
    TH1Sng_BtagAll = new TH1D(titleSng_BtagAll_,titleSng_BtagAll_,nBinsBtag[nBdiscrAlgos],rangesbTag_);
    TH1Bkg_BtagAll = new TH1D(titleBkg_BtagAll_,titleBkg_BtagAll_,nBinsBtag[nBdiscrAlgos],rangesbTag_);
    TH1Data_BtagAll = new TH1D(titleData_BtagAll_,titleData_BtagAll_,nBinsBtag[nBdiscrAlgos],rangesbTag_);
  }
 

  GiveName(&titleSng_BtagEffAll_); titleSng_BtagEffAll_+="TH1Sng_BtagEffAll"; 
  GiveName(&titleBkg_BtagEffAll_); titleBkg_BtagEffAll_+="TH1Bkg_BtagEffAll"; 
  GiveName(&titleData_BtagEffAll_); titleData_BtagEffAll_+="TH1Data_BtagEffAll"; 
  
  if(!varBinSize_){
    TH1Sng_BtagEffAll = new TH1D(titleSng_BtagEffAll_,titleSng_BtagEffAll_,nBinsBtag[nBdiscrAlgos],lowRangeBtag[nBdiscrAlgos],upRangeBtag[nBdiscrAlgos]);
    TH1Bkg_BtagEffAll = new TH1D(titleBkg_BtagEffAll_,titleBkg_BtagEffAll_,nBinsBtag[nBdiscrAlgos],lowRangeBtag[nBdiscrAlgos],upRangeBtag[nBdiscrAlgos]);
    TH1Data_BtagEffAll = new TH1D(titleData_BtagEffAll_,titleData_BtagEffAll_,nBinsBtag[nBdiscrAlgos],lowRangeBtag[nBdiscrAlgos],upRangeBtag[nBdiscrAlgos]);
  }
  
  if(varBinSize_){
    TH1Sng_BtagEffAll = new TH1D(titleSng_BtagEffAll_,titleSng_BtagEffAll_,nBinsBtag[nBdiscrAlgos],rangesbTag_);
    TH1Bkg_BtagEffAll = new TH1D(titleBkg_BtagEffAll_,titleBkg_BtagEffAll_,nBinsBtag[nBdiscrAlgos],rangesbTag_);
    TH1Data_BtagEffAll = new TH1D(titleData_BtagEffAll_,titleData_BtagEffAll_,nBinsBtag[nBdiscrAlgos],rangesbTag_);
  }


  TH2Sng_Left->Sumw2(); 
  TH2Bkg_Left->Sumw2(); 
  TH2Data_Left->Sumw2(); 
  TH2Sng_Right->Sumw2(); 
  TH2Bkg_Right->Sumw2(); 
  TH2Data_Right->Sumw2(); 
  TH1Sng_Var0->Sumw2(); 
  TH1Bkg_Var0->Sumw2(); 
  TH1Bkg_W_Var0->Sumw2(); 
  TH1Bkg_R_Var0->Sumw2(); 
    TH1Data_Var0->Sumw2();   
    TH1Data_Var0_XS->Sumw2();   
  TH1SoverSB_Var0->Sumw2(); 

  //  TH1Sng_Var1->Sumw2(); 
  //TH1Bkg_Var1->Sumw2(); 
  //  TH1Data_Var1->Sumw2(); 
  TH2Sng_Var12->Sumw2();
  TH2Bkg_Var12->Sumw2();
  TH2Bkg_W_Var12->Sumw2();
  TH2Bkg_R_Var12->Sumw2();
  TH2Data_Var12->Sumw2();
 
  TH1Sng_BtagAll->Sumw2(); 
  TH1Bkg_BtagAll->Sumw2(); 
  TH1Data_BtagAll->Sumw2(); 

  TH2Sng->Sumw2(); 
  TH2Bkg->Sumw2(); 
  TH2Data->Sumw2();  

  TH2SngVar12_Left->Sumw2(); 
  TH2BkgVar12_Left->Sumw2(); 
  TH2DataVar12_Left->Sumw2(); 
  TH2SngVar12_Right->Sumw2(); 
  TH2BkgVar12_Right->Sumw2(); 
  TH2DataVar12_Right->Sumw2(); 
	
	GiveName(&titleData_Pt_Left); titleData_Pt_Left+="TH1Data_Pt_Left"; 
	TH1Data_Pt_Left = new TH1D(titleData_Pt_Left,titleData_Pt_Left,100,0,500);
    
    GiveName(&titleSng_Pt_Left); titleSng_Pt_Left+="TH1Sng_Pt_Left"; 
	TH1Sng_Pt_Left = new TH1D(titleSng_Pt_Left,titleSng_Pt_Left,100,0,500);
	
	GiveName(&titleData_Pt_Right); titleData_Pt_Right+="TH1Data_Pt_Right"; 
	TH1Data_Pt_Right = new TH1D(titleData_Pt_Right,titleData_Pt_Right,100,0,500);

	GiveName(&titleData_Pt_Control); titleData_Pt_Control+="TH1Data_Pt_Control"; 
	TH1Data_Pt_Control = new TH1D(titleData_Pt_Control,titleData_Pt_Control,100,0,500);
	
	TH1Data_Pt_Left->Sumw2(); 
	TH1Sng_Pt_Left->Sumw2(); 
	TH1Data_Pt_Right->Sumw2(); 
	TH1Data_Pt_Control->Sumw2(); 
    	
	//nBinsVar0=70;
	//lowRangeVar0=0;
	//upRangeVar0=700;	
	
	TString name = ""; GiveName(&name); name+="TH1Data_Var0_bTagL"; 
	histo1D["TH1Data_Var0_bTagL"] = new TH1D(name,name,nBinsVar0,lowRangeVar0,upRangeVar0);
	
	name = ""; GiveName(&name); name+="TH1Sng_Var0_bTagL"; 
	histo1D["TH1Sng_Var0_bTagL"] = new TH1D(name,name,nBinsVar0,lowRangeVar0,upRangeVar0);
	
	name = ""; GiveName(&name); name+="TH1Bkg_Var0_bTagL"; 
	histo1D["TH1Bkg_Var0_bTagL"] = new TH1D(name,name,nBinsVar0,lowRangeVar0,upRangeVar0);
	
	name = ""; GiveName(&name); name+="TH1Data_Var0_bTagM"; 
	histo1D["TH1Data_Var0_bTagM"] = new TH1D(name,name,nBinsVar0,lowRangeVar0,upRangeVar0);
	
	name = ""; GiveName(&name); name+="TH1Sng_Var0_bTagM"; 
	histo1D["TH1Sng_Var0_bTagM"] = new TH1D(name,name,nBinsVar0,lowRangeVar0,upRangeVar0);
	
	name = ""; GiveName(&name); name+="TH1Bkg_Var0_bTagM"; 
	histo1D["TH1Bkg_Var0_bTagM"] = new TH1D(name,name,nBinsVar0,lowRangeVar0,upRangeVar0);
	
	
	
	name = ""; GiveName(&name); name+="TH1Data_Var0_bTagT"; 
	histo1D["TH1Data_Var0_bTagT"] = new TH1D(name,name,nBinsVar0,lowRangeVar0,upRangeVar0);
	
	name = ""; GiveName(&name); name+="TH1Sng_Var0_bTagT"; 
	histo1D["TH1Sng_Var0_bTagT"] = new TH1D(name,name,nBinsVar0,lowRangeVar0,upRangeVar0);
	
	name = ""; GiveName(&name); name+="TH1Bkg_Var0_bTagT"; 
	histo1D["TH1Bkg_Var0_bTagT"] = new TH1D(name,name,nBinsVar0,lowRangeVar0,upRangeVar0);
	
	
	
	name = ""; GiveName(&name); name+="TH1Data_Var0_VVMC"; 
	histo1D["TH1Data_Var0_VVMC"] = new TH1D(name,name,nBinsVar0,lowRangeVar0,upRangeVar0);
	
	name = ""; GiveName(&name); name+="TH1Sng_Var0_VVMC"; 
	histo1D["TH1Sng_Var0_VVMC"] = new TH1D(name,name,nBinsVar0,lowRangeVar0,upRangeVar0);
	
	name = ""; GiveName(&name); name+="TH1Bkg_Var0_VVMC"; 
	histo1D["TH1Bkg_Var0_VVMC"] = new TH1D(name,name,nBinsVar0,lowRangeVar0,upRangeVar0);
	
	name = ""; GiveName(&name); name+="TH1Data_Var0_VVMC_bTagL"; 
	histo1D["TH1Data_Var0_VVMC_bTagL"] = new TH1D(name,name,nBinsVar0,lowRangeVar0,upRangeVar0);
	
	name = ""; GiveName(&name); name+="TH1Sng_Var0_VVMC_bTagL"; 
	histo1D["TH1Sng_Var0_VVMC_bTagL"] = new TH1D(name,name,nBinsVar0,lowRangeVar0,upRangeVar0);
	
	name = ""; GiveName(&name); name+="TH1Bkg_Var0_VVMC_bTagL"; 
	histo1D["TH1Bkg_Var0_VVMC_bTagL"] = new TH1D(name,name,nBinsVar0,lowRangeVar0,upRangeVar0);
	
	name = ""; GiveName(&name); name+="TH1Data_Var0_VVMC_bTagM"; 
	histo1D["TH1Data_Var0_VVMC_bTagM"] = new TH1D(name,name,nBinsVar0,lowRangeVar0,upRangeVar0);
	
	name = ""; GiveName(&name); name+="TH1Sng_Var0_VVMC_bTagM"; 
	histo1D["TH1Sng_Var0_VVMC_bTagM"] = new TH1D(name,name,nBinsVar0,lowRangeVar0,upRangeVar0);
	
	name = ""; GiveName(&name); name+="TH1Bkg_Var0_VVMC_bTagM"; 
	histo1D["TH1Bkg_Var0_VVMC_bTagM"] = new TH1D(name,name,nBinsVar0,lowRangeVar0,upRangeVar0);	
	
	name = ""; GiveName(&name); name+="TH1Data_Var0_VVMC_bTagT"; 
	histo1D["TH1Data_Var0_VVMC_bTagT"] = new TH1D(name,name,nBinsVar0,lowRangeVar0,upRangeVar0);
	
	name = ""; GiveName(&name); name+="TH1Sng_Var0_VVMC_bTagT"; 
	histo1D["TH1Sng_Var0_VVMC_bTagT"] = new TH1D(name,name,nBinsVar0,lowRangeVar0,upRangeVar0);
	
	name = ""; GiveName(&name); name+="TH1Bkg_Var0_VVMC_bTagT"; 
	histo1D["TH1Bkg_Var0_VVMC_bTagT"] = new TH1D(name,name,nBinsVar0,lowRangeVar0,upRangeVar0);
	
	name = ""; GiveName(&name); name+="TH1Data_Var0_VVTTotherMC"; 
	histo1D["TH1Data_Var0_VVTTotherMC"] = new TH1D(name,name,nBinsVar0,lowRangeVar0,upRangeVar0);
	
	name = ""; GiveName(&name); name+="TH1Data_Var0_TTbar"; 
	histo1D["TH1Data_Var0_TTbar"] = new TH1D(name,name,nBinsVar0,lowRangeVar0,upRangeVar0);
	
	name = ""; GiveName(&name); name+="TH1Sng_Var0_TTbar"; 
	histo1D["TH1Sng_Var0_TTbar"] = new TH1D(name,name,nBinsVar0,lowRangeVar0,upRangeVar0);
	
	name = ""; GiveName(&name); name+="TH1Bkg_Var0_TTbar"; 
	histo1D["TH1Bkg_Var0_TTbar"] = new TH1D(name,name,nBinsVar0,lowRangeVar0,upRangeVar0);
	
	name = ""; GiveName(&name); name+="TH1Data_Var0_TTbar_bTagL"; 
	histo1D["TH1Data_Var0_TTbar_bTagL"] = new TH1D(name,name,nBinsVar0,lowRangeVar0,upRangeVar0);
	
	name = ""; GiveName(&name); name+="TH1Sng_Var0_TTbar_bTagL"; 
	histo1D["TH1Sng_Var0_TTbar_bTagL"] = new TH1D(name,name,nBinsVar0,lowRangeVar0,upRangeVar0);
	
	name = ""; GiveName(&name); name+="TH1Bkg_Var0_TTbar_bTagL"; 
	histo1D["TH1Bkg_Var0_TTbar_bTagL"] = new TH1D(name,name,nBinsVar0,lowRangeVar0,upRangeVar0);	
	
	name = ""; GiveName(&name); name+="TH1Data_Var0_TTbar_bTagM"; 
	histo1D["TH1Data_Var0_TTbar_bTagM"] = new TH1D(name,name,nBinsVar0,lowRangeVar0,upRangeVar0);
	
	name = ""; GiveName(&name); name+="TH1Sng_Var0_TTbar_bTagM"; 
	histo1D["TH1Sng_Var0_TTbar_bTagM"] = new TH1D(name,name,nBinsVar0,lowRangeVar0,upRangeVar0);
	
	name = ""; GiveName(&name); name+="TH1Bkg_Var0_TTbar_bTagM"; 
	histo1D["TH1Bkg_Var0_TTbar_bTagM"] = new TH1D(name,name,nBinsVar0,lowRangeVar0,upRangeVar0);	
	
	name = ""; GiveName(&name); name+="TH1Data_Var0_TTbar_bTagT"; 
	histo1D["TH1Data_Var0_TTbar_bTagT"] = new TH1D(name,name,nBinsVar0,lowRangeVar0,upRangeVar0);
	
	name = ""; GiveName(&name); name+="TH1Sng_Var0_TTbar_bTagT"; 
	histo1D["TH1Sng_Var0_TTbar_bTagT"] = new TH1D(name,name,nBinsVar0,lowRangeVar0,upRangeVar0);
	
	name = ""; GiveName(&name); name+="TH1Bkg_Var0_TTbar_bTagT"; 
	histo1D["TH1Bkg_Var0_TTbar_bTagT"] = new TH1D(name,name,nBinsVar0,lowRangeVar0,upRangeVar0);	
	
	
	name = ""; GiveName(&name); name+="TH1Data_Var0_TTbarSemiMu"; 
	histo1D["TH1Data_Var0_TTbarSemiMu"] = new TH1D(name,name,nBinsVar0,lowRangeVar0,upRangeVar0);
	
	name = ""; GiveName(&name); name+="TH1Data_Var0_TTbarOther"; 
	histo1D["TH1Data_Var0_TTbarOther"] = new TH1D(name,name,nBinsVar0,lowRangeVar0,upRangeVar0);
	
    //for (std::map<TString,TH1D*>::iterator it=histo1D.begin(); it != histo1D.end(); it++)
    //    it->second->Sumw2();
	
}

void PtEtaBin::DefineControlSamplePlots(int nBinsControlVar1, double lowRangControlVar1, double upRangeControlVar1, int nBinsControlVar2, double lowRangControlVar2, double upRangeControlVar2, int nBinsControlVar, double lowRangeControlVar, double upRangeControlVar){

  GiveName(&titleSng_ControlVar_); titleSng_ControlVar_+="TH1Sng_ControlVar"; 
  GiveName(&titleBkg_ControlVar_); titleBkg_ControlVar_+="TH1Bkg_ControlVar"; 
  GiveName(&titleBkg_W_ControlVar_); titleBkg_W_ControlVar_+="TH1Bkg_W_ControlVar"; 
  GiveName(&titleBkg_R_ControlVar_); titleBkg_R_ControlVar_+="TH1Bkg_R_ControlVar"; 
  GiveName(&titleData_ControlVar_); titleData_ControlVar_+="TH1Data_ControlVar"; 
 
  TH1Sng_ControlVar = new TH1D(titleSng_ControlVar_,titleSng_ControlVar_,nBinsControlVar,lowRangeControlVar,upRangeControlVar);
  TH1Bkg_ControlVar = new TH1D(titleBkg_ControlVar_,titleBkg_ControlVar_,nBinsControlVar,lowRangeControlVar,upRangeControlVar);
  TH1Bkg_W_ControlVar = new TH1D(titleBkg_W_ControlVar_,titleBkg_W_ControlVar_,nBinsControlVar,lowRangeControlVar,upRangeControlVar);
  TH1Bkg_R_ControlVar = new TH1D(titleBkg_R_ControlVar_,titleBkg_R_ControlVar_,nBinsControlVar,lowRangeControlVar,upRangeControlVar);
  TH1Data_ControlVar = new TH1D(titleData_ControlVar_,titleData_ControlVar_,nBinsControlVar,lowRangeControlVar,upRangeControlVar);
  
  TH1Sng_ControlVar->Sumw2(); 
  TH1Bkg_ControlVar->Sumw2(); 
  TH1Bkg_W_ControlVar->Sumw2(); 
  TH1Bkg_R_ControlVar->Sumw2();
  TH1Data_ControlVar->Sumw2();  

  GiveName(&titleSoverSB_ControlVar_); titleSoverSB_ControlVar_+="TH1SoverSB_ControlVar"; 
  TH1SoverSB_ControlVar = new TH1D(titleSoverSB_ControlVar_,titleSoverSB_ControlVar_,nBinsControlVar,lowRangeControlVar,upRangeControlVar);

  TH1SoverSB_ControlVar->Sumw2();
	
  // templates for Xsection 
	
	TString name = "";
	
	//name = ""; GiveName(&name); name+="TH1Data_Var0_VVData"; 
	//histo1D["TH1Data_Var0_VVData"] = new TH1D(name,name,nBinsControlVar,lowRangeControlVar,upRangeControlVar);
	
	name = ""; GiveName(&name); name+="TH1Data_Var0_VVData"; 
	histo1D["TH1Data_Var0_VVData"] = new TH1D(name,name,nBinsControlVar,lowRangeControlVar,upRangeControlVar);

	//name = ""; GiveName(&name); name+="TH1Sng_Var0"; 
	//histo1D["TH1Sng_Var0"] = new TH1D(name,name,nBinsControlVar,lowRangeControlVar,upRangeControlVar);
	
	//name = ""; GiveName(&name); name+="TH1Sng_Var0_bTagL"; 
	//histo1D["TH1Sng_Var0_bTagL"] = new TH1D(name,name,nBinsControlVar,lowRangeControlVar,upRangeControlVar);

	
	for (std::map<TString,TH1D*>::const_iterator it=histo1D.begin(); it != histo1D.end(); ++it)
		
		it->second->Sumw2();


  GiveName(&titleSng_ControlVarReweigh_); titleSng_ControlVarReweigh_+="TH1Sng_ControlVarReweigh"; 
  GiveName(&titleBkg_ControlVarReweigh_); titleBkg_ControlVarReweigh_+="TH1Bkg_ControlVarReweigh"; 
  GiveName(&titleBkg_W_ControlVarReweigh_); titleBkg_W_ControlVarReweigh_+="TH1Bkg_W_ControlVarReweigh"; 
  GiveName(&titleBkg_R_ControlVarReweigh_); titleBkg_R_ControlVarReweigh_+="TH1Bkg_R_ControlVarReweigh"; 
  GiveName(&titleData_ControlVarReweigh_); titleData_ControlVarReweigh_+="TH1Data_ControlVarReweigh"; 
 
  TH1Sng_ControlVarReweigh = new TH1D(titleSng_ControlVarReweigh_,titleSng_ControlVarReweigh_,nBinsControlVar,lowRangeControlVar,upRangeControlVar);
  TH1Bkg_ControlVarReweigh = new TH1D(titleBkg_ControlVarReweigh_,titleBkg_ControlVarReweigh_,nBinsControlVar,lowRangeControlVar,upRangeControlVar);
  TH1Bkg_W_ControlVarReweigh = new TH1D(titleBkg_W_ControlVarReweigh_,titleBkg_W_ControlVarReweigh_,nBinsControlVar,lowRangeControlVar,upRangeControlVar);
  TH1Bkg_R_ControlVarReweigh = new TH1D(titleBkg_R_ControlVarReweigh_,titleBkg_R_ControlVarReweigh_,nBinsControlVar,lowRangeControlVar,upRangeControlVar);
  TH1Data_ControlVarReweigh = new TH1D(titleData_ControlVarReweigh_,titleData_ControlVarReweigh_,nBinsControlVar,lowRangeControlVar,upRangeControlVar);

  //nBinsControlVar=20;
  //TH1Sng_ControlVarReweigh = new TH1D(titleSng_ControlVarReweigh_,titleSng_ControlVarReweigh_,nBinsControlVar,rangesVar1_);
  //TH1Bkg_ControlVarReweigh = new TH1D(titleBkg_ControlVarReweigh_,titleBkg_ControlVarReweigh_,nBinsControlVar,rangesVar1_);
  //TH1Data_ControlVarReweigh = new TH1D(titleData_ControlVarReweigh_,titleData_ControlVarReweigh_,nBinsControlVar,rangesVar1_);

  TH1Sng_ControlVarReweigh->Sumw2(); 
  TH1Bkg_ControlVarReweigh->Sumw2(); 
  TH1Bkg_W_ControlVarReweigh->Sumw2(); 
  TH1Bkg_R_ControlVarReweigh->Sumw2(); 
  TH1Data_ControlVarReweigh->Sumw2();  

  /*GiveName(&titleSng_ControlVar1_); titleSng_ControlVar1_+="TH1Sng_ControlVar1"; 
  GiveName(&titleBkg_ControlVar1_); titleBkg_ControlVar1_+="TH1Bkg_ControlVar1"; 
  GiveName(&titleData_ControlVar1_); titleData_ControlVar1_+="TH1Data_ControlVar1"; 
 
  TH1Sng_ControlVar1 = new TH1D(titleSng_ControlVar1_,titleSng_ControlVar1_,nBinsControlVar1,lowRangControlVar1,upRangeControlVar1);
  TH1Bkg_ControlVar1 = new TH1D(titleBkg_ControlVar1_,titleBkg_ControlVar1_,nBinsControlVar1,lowRangControlVar1,upRangeControlVar1);
  TH1Data_ControlVar1 = new TH1D(titleData_ControlVar1_,titleData_ControlVar1_,nBinsControlVar1,lowRangControlVar1,upRangeControlVar1);
  
  TH1Sng_ControlVar1->Sumw2(); 
  TH1Bkg_ControlVar1->Sumw2(); 
  TH1Data_ControlVar1->Sumw2();  */

  GiveName(&titleSng_ControlVar12_); titleSng_ControlVar12_+="TH2Sng_ControlVar12"; 
  GiveName(&titleBkg_ControlVar12_); titleBkg_ControlVar12_+="TH2Bkg_ControlVar12"; 
  GiveName(&titleBkg_W_ControlVar12_); titleBkg_W_ControlVar12_+="TH2Bkg_W_ControlVar12"; 
  GiveName(&titleBkg_R_ControlVar12_); titleBkg_R_ControlVar12_+="TH2Bkg_R_ControlVar12"; 
  GiveName(&titleData_ControlVar12_); titleData_ControlVar12_+="TH2Data_ControlVar12"; 
  
  /*TH2Sng_ControlVar12 = new TH2D(titleSng_ControlVar12_,titleSng_ControlVar12_,nBinsControlVar1,lowRangControlVar1,upRangeControlVar1,nBinsControlVar2,lowRangControlVar2,upRangeControlVar2);
  TH2Bkg_ControlVar12 = new TH2D(titleBkg_ControlVar12_,titleBkg_ControlVar12_,nBinsControlVar1,lowRangControlVar1,upRangeControlVar1,nBinsControlVar2,lowRangControlVar2,upRangeControlVar2);
  TH2Bkg_W_ControlVar12 = new TH2D(titleBkg_W_ControlVar12_,titleBkg_W_ControlVar12_,nBinsControlVar1,lowRangControlVar1,upRangeControlVar1,nBinsControlVar2,lowRangControlVar2,upRangeControlVar2);
  TH2Bkg_R_ControlVar12 = new TH2D(titleBkg_R_ControlVar12_,titleBkg_R_ControlVar12_,nBinsControlVar1,lowRangControlVar1,upRangeControlVar1,nBinsControlVar2,lowRangControlVar2,upRangeControlVar2);
  TH2Data_ControlVar12 = new TH2D(titleData_ControlVar12_,titleData_ControlVar12_,nBinsControlVar1,lowRangControlVar1,upRangeControlVar1,nBinsControlVar2,lowRangControlVar2,upRangeControlVar2); */
  
  
  //  double rangesPtSC[21] = {30,33.9517,38.0413,42.3827,47.1889,52.4274,58.0206,63.9226,70.6819,77.5725,85.4028,94.2801,104.116,115.42,128.936,145.192,165.937,194.182,237.506,316.667,500};
  double rangesPtSC[19] = {30,33.9517,38.0413,42.3827,47.1889,52.4274,58.0206,63.9226,70.6819,77.5725,85.4028,94.2801,104.116,115.42,128.936,145.192,165.937,250,9999.};
  
  //double rangesPtSC[19] = {30,50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350,370,400};
  //double rangesPtSC[51] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400,410,420,430,440,450,460,470,480,490,500};
  
  //double rangesEtaSC[11] = {0,0.24,0.48,0.72,0.96,1.20,1.44,1.68,1.92,2.16,2.4};
  double rangesEtaSC[21] = {0,0.12,0.24,0.36,0.48,0.60,0.72,0.84,0.96,1.08,1.20,1.32,1.44,1.56,1.68,1.80,1.92,2.04,2.16,2.28,2.4};
  //double rangesEtaSC[21] = {0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1};
  //double rangesEtaSC[21] = {-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
  //double rangesEtaSC[21] = {0,0.175,0.175*2,0.175*3,0.175*4,0.175*5,0.175*6,0.175*7,0.175*8,0.175*9,0.175*10,0.175*11,0.175*12,0.175*13,0.175*14,0.175*15,0.175*16,0.175*17,0.175*18,0.175*19,0.175*20};
  
  
  TH2Sng_ControlVar12 = new TH2D(titleSng_ControlVar12_,titleSng_ControlVar12_,nBinsControlVar1,rangesPtSC,nBinsControlVar2,rangesEtaSC);
  TH2Bkg_ControlVar12 = new TH2D(titleBkg_ControlVar12_,titleBkg_ControlVar12_,nBinsControlVar1,rangesPtSC,nBinsControlVar2,rangesEtaSC);
  TH2Bkg_W_ControlVar12 = new TH2D(titleBkg_W_ControlVar12_,titleBkg_W_ControlVar12_,nBinsControlVar1,rangesPtSC,nBinsControlVar2,rangesEtaSC);
  TH2Bkg_R_ControlVar12 = new TH2D(titleBkg_R_ControlVar12_,titleBkg_R_ControlVar12_,nBinsControlVar1,rangesPtSC,nBinsControlVar2,rangesEtaSC);
  TH2Data_ControlVar12 = new TH2D(titleData_ControlVar12_,titleData_ControlVar12_,nBinsControlVar1,rangesPtSC,nBinsControlVar2,rangesEtaSC);
  

  TH2Sng_ControlVar12->Sumw2();
  TH2Bkg_ControlVar12->Sumw2();
  TH2Bkg_W_ControlVar12->Sumw2();
  TH2Bkg_R_ControlVar12->Sumw2();
  TH2Data_ControlVar12->Sumw2();
 
}

void PtEtaBin::SetError(TH1D* histo){

  for(int i=0; i<=histo->GetXaxis()->GetNbins()+1; i++){
    histo->SetBinError(i,sqrt(histo->GetBinContent(i)));
  }

}

void PtEtaBin::SetErrorsSignalSamples(){

  SetError(TH1Sng_BtagAll);   
  SetError(TH1Bkg_BtagAll);
  SetError(TH1Data_BtagAll); 
  SetError(TH1Sng_Var0);
  SetError(TH1Bkg_Var0);
  SetError(TH1Bkg_W_Var0);
  SetError(TH1Bkg_R_Var0);
  SetError(TH1Data_Var0);

}

void PtEtaBin::SetErrorsControlSamples(){

  SetError(TH1Sng_ControlVar);   
  SetError(TH1Bkg_ControlVar);
  SetError(TH1Bkg_W_ControlVar);
  SetError(TH1Bkg_R_ControlVar);
  SetError(TH1Data_ControlVar);

}

//void PtEtaBin::FillSignalSamplePlots(double weight, int partonFlavour, double bTag, double var1, double var2, double var0, int partonFlavourControl, double bTagControl, double controlVar1, double controlVar2, double controlVar0, double lowCutVar0, double centralCutVar0, double upCutVar0){//I should replace the name var1 with var1
void PtEtaBin::FillSignalSamplePlots(double weight, double weight_nonrew, int partonFlavour, bool isW, bool isR, double chisq, double bTag, double* bTagCuts, double var1, double var2, double var0, double lowCutVar0, double centralLowCutVar0, double centralUpCutVar0, double upCutVar0, double defaultbTag){

	if(var1>ptbinlow_ && var1<ptbinup_){
		
		
		if(var2>etabinlow_ && var2<etabinup_){ 
			
			if((var0>lowCutVar0 && var0<centralLowCutVar0) || (var0>=centralUpCutVar0 && var0<upCutVar0)){
				if(fabs(partonFlavour)==5){
					TH2Sng->Fill(var0,bTag,weight);
					TH1Sng_chisq->Fill(chisq);
				} else {
					TH2Bkg->Fill(var0,bTag,weight);
					TH1Bkg_chisq->Fill(chisq);
				}
				TH2Data->Fill(var0,bTag,weight);
				TH1Data_chisq->Fill(chisq);
			}
			
			if(var0>lowCutVar0 && var0<centralLowCutVar0){
				if(fabs(partonFlavour)==5){
					TH2Sng_Left->Fill(var1,bTag,weight);
					TH2SngVar12_Left->Fill(var1,var2,weight);
				} else {
					TH2Bkg_Left->Fill(var1,bTag,weight);
					TH2BkgVar12_Left->Fill(var1,var2,weight);
				}
				TH2Data_Left->Fill(var1,bTag,weight);
				TH2DataVar12_Left->Fill(var1,var2,weight);
				TH1Data_Pt_Left->Fill(var1,weight);
                //if(fabs(partonFlavour)==5)
                TH1Sng_Pt_Left->Fill(var1,weight_nonrew);
                //if (var1 < 0) cout << var0 << " " << var1 << endl; exit(1);
			}
			
			//original right plot
			if(var0>=centralUpCutVar0 && var0<upCutVar0){
				if(fabs(partonFlavour)==5){
					TH2Sng_Right->Fill(var1,bTag,weight);
					TH2SngVar12_Right->Fill(var1,var2,weight);
				} else {
					TH2Bkg_Right->Fill(var1,bTag,weight);
					TH2BkgVar12_Right->Fill(var1,var2,weight);
				}
				TH2Data_Right->Fill(var1,bTag,weight);
				TH2DataVar12_Right->Fill(var1,var2,weight);
				TH1Data_Pt_Right->Fill(var1,weight);
				
				//cout << "right weight " << weight << " -- " << bTag << endl;
				
				
			}
		}
	} else {
		
		//cout << ptbinlow_ << " " << ptbinup_ << "  - " << var1 << endl;
		
	}
	
	//double SCweight=(2735.27/3992.1)*2281.25/3095.5;//dirty hard-coded
	//this one is to use the right region from control sample
	//  double SCweight=(2735.27/3992.1);//dirty hard-coded
	//if(controlVar0>=centralCutVar0 && controlVar0<upCutVar0){
	
	//this one is to use the left region of the control sample
	//double SCweight=(2735.27/9927.5);//dirty hard-coded
	//double SCweight=1;//dirty hard-coded
	//double SCweight=2294./6721.;//dirty hard-coded -> rescaling the non-b content 90 160 210
	//double SCweight=637.389/12265.; 
	/*if(controlVar0>lowCutVar0 && controlVar0<centralCutVar0){
	 if(fabs(partonFlavourControl)==5){
	 TH2Sng_Right->Fill(controlVar1,bTagControl,weight*SCweight);
	 TH2SngVar12_Right->Fill(controlVar1,controlVar2,weight*SCweight);
	 } else {
	 TH2Bkg_Right->Fill(controlVar1,bTagControl,weight*SCweight);
	 TH2BkgVar12_Right->Fill(controlVar1,controlVar2,weight*SCweight);
	 }
	 TH2Data_Right->Fill(controlVar1,bTagControl,weight*SCweight);
	 TH2DataVar12_Right->Fill(controlVar1,controlVar2,weight*SCweight);
	 }*/
	
	//    if(var0>lowCutVar0 && var0<centralLowCutVar0){
	
	if(var1>ptbinlow_ && var1<ptbinup_){
		if(var2>etabinlow_ && var2<etabinup_){ 
			
			if(fabs(partonFlavour)==5){
				//cout << var0 << endl;
				TH1Sng_Var0->Fill(var0,weight);
				TH1Sng_BtagAll->Fill(bTag,weight_nonrew);
				
				
			} else { 
				if(isW) TH1Bkg_W_Var0->Fill(var0,weight);
				if(isR) TH1Bkg_R_Var0->Fill(var0,weight);
				TH1Bkg_Var0->Fill(var0,weight); 
				TH1Bkg_BtagAll->Fill(bTag,weight);
			}
			
			TH1Data_Var0->Fill(var0,weight); 
			TH1Data_Var0_XS->Fill(var0,weight_nonrew); 
			TH1Data_BtagAll->Fill(bTag,weight);
			
			if (bTag > bTagCuts[0]) {
				histo1D["TH1Data_Var0_bTagL"]->Fill(var0,weight_nonrew);
				if(fabs(partonFlavour)==5) {
					histo1D["TH1Sng_Var0_bTagL"]->Fill(var0,weight_nonrew);
				} else {
					histo1D["TH1Bkg_Var0_bTagL"]->Fill(var0,weight_nonrew);
				}
			}
			
			if (bTag > bTagCuts[1]) {
				histo1D["TH1Data_Var0_bTagM"]->Fill(var0,weight_nonrew);
				if(fabs(partonFlavour)==5) {
					histo1D["TH1Sng_Var0_bTagM"]->Fill(var0,weight_nonrew);
				} else {
					histo1D["TH1Bkg_Var0_bTagM"]->Fill(var0,weight_nonrew);
				}
			}				
			if (bTag > bTagCuts[2]) {
				histo1D["TH1Data_Var0_bTagT"]->Fill(var0,weight_nonrew);
				if(fabs(partonFlavour)==5) {
					histo1D["TH1Sng_Var0_bTagT"]->Fill(var0,weight_nonrew);
				} else {
					histo1D["TH1Bkg_Var0_bTagT"]->Fill(var0,weight_nonrew);
				}
				
			}
		}
	}
	//}  
	
	//if((var0>lowCutVar0 && var0<centralLowCutVar0) || (var0>=centralUpCutVar0 && var0<upCutVar0)) {
	//Start::filling of Var12 reweighing plots (full mlb range and L+R mlb range)
	if(fabs(partonFlavour)==5){
		if(defaultbTag<3) TH2Sng_Var12->Fill(var1,var2,weight);
	} else {
		if(defaultbTag<3) TH2Bkg_Var12->Fill(var1,var2,weight);
		if(defaultbTag<3) if(isW) TH2Bkg_W_Var12->Fill(var1,var2,weight);
		if(defaultbTag<3) if(isR) TH2Bkg_R_Var12->Fill(var1,var2,weight);
	}
	
	//if(isW) if(defaultbTag<3) TH2Data_Var12->Fill(var1,var2,weight); //in the case I want to do the reweighting for a certain component only (eg W or R or W&R)
	if(defaultbTag<3) TH2Data_Var12->Fill(var1,var2,weight); 
	//}
	
	
	//End::filling of Var12 reweighing plots (full mlb range and L+R mlb range)
	
	//Start::Filling the plots for the reweighing between different chi-square cuts
	//if(fabs(partonFlavour)==5){
	/*if(fabs(partonFlavour)!=5){
	 TH2Sng_Var12->Fill(var1,var2,weight);
	 if(chisq<3) TH2Bkg_Var12->Fill(var1,var2,weight);
	 //if(chisq<3) TH2Sng_Var12->Fill(var1,var2,weight);
	 //TH2Bkg_Var12->Fill(var1,var2,weight);
	 }*/
	
	//define the bin where the cut is being applied (Important is to not rebin afterwards!)
	for(int i=0; i<=TH1Data_Var0->GetXaxis()->GetNbins()+1; i++){
		//cout << "low " << TH2Data_Fix->GetXaxis()->GetBinLowEdge(i) << " up " << TH2Data_Fix->GetXaxis()->GetBinUpEdge(i) << endl;
		if(lowCutVar0>=TH1Data_Var0->GetXaxis()->GetBinLowEdge(i) && lowCutVar0<TH1Data_Var0->GetXaxis()->GetBinUpEdge(i)) lowerBin_=i;
		if(centralLowCutVar0>TH1Data_Var0->GetXaxis()->GetBinLowEdge(i) && centralLowCutVar0<=TH1Data_Var0->GetXaxis()->GetBinUpEdge(i)){
			centralLowBin_=i;
		}
		if(centralUpCutVar0>TH1Data_Var0->GetXaxis()->GetBinLowEdge(i) && centralUpCutVar0<=TH1Data_Var0->GetXaxis()->GetBinUpEdge(i)){
			centralUpBin_=i;
		}
		if(upCutVar0>TH1Data_Var0->GetXaxis()->GetBinLowEdge(i) && upCutVar0<=TH1Data_Var0->GetXaxis()->GetBinUpEdge(i)) upperBin_=i;
	}
	
	
	
	if(lowCutVar0!=TH1Data_Var0->GetXaxis()->GetBinLowEdge(lowerBin_)) cout << "a ERROR PtEtaBin::FillSignalSamplePlots lowCutVar0 does is not agree with binning " << lowCutVar0 << " " << TH1Data_Var0->GetXaxis()->GetBinLowEdge(lowerBin_) <<  endl;
	if(centralLowCutVar0!=TH1Data_Var0->GetXaxis()->GetBinUpEdge(centralLowBin_)) cout << "b ERROR PtEtaBin::FillSignalSamplePlots lowCutVar0 does is not agree with binning " << centralLowCutVar0 << " " << TH1Data_Var0->GetXaxis()->GetBinLowEdge(centralLowBin_) << endl;
	if(centralUpCutVar0!=TH1Data_Var0->GetXaxis()->GetBinUpEdge(centralUpBin_)) cout << "c ERROR PtEtaBin::FillSignalSamplePlots lowCutVar0 does is not agree with binning " << centralUpCutVar0 << " " << TH1Data_Var0->GetXaxis()->GetBinLowEdge(centralUpBin_) << endl;
	if(upCutVar0!=TH1Data_Var0->GetXaxis()->GetBinUpEdge(upperBin_)) cout << "d ERROR PtEtaBin::FillSignalSamplePlots lowCutVar0 does is not agree with binning " << upCutVar0 << " " << TH1Data_Var0->GetXaxis()->GetBinUpEdge(upperBin_) << endl;
	
	//   cout << "+----> Bin : " << lowerBin_ << " " << centralLowBin_ << " " << centralUpBin_ << " " << upperBin_ << " " << endl;
	
	
	
}

void PtEtaBin::FillControlSamplePlots(double weight, int partonFlavour, bool isW, bool isR,  double chisq, double controlVar1, double controlVar2, double controlVar, double lowCutVar0, double centralLowCutVar0, double centralUpCutVar0, double upCutVar0){
 
	if(controlVar1>ptbinlow_ && controlVar1<ptbinup_){
		if(controlVar2>etabinlow_ && controlVar2<etabinup_){
			if(fabs(partonFlavour)==5){
				TH1Sng_ControlVar->Fill(controlVar,weight);
			} else {
				TH1Bkg_ControlVar->Fill(controlVar,weight);
				if(isW) TH1Bkg_W_ControlVar->Fill(controlVar,weight);
				if(isR) TH1Bkg_R_ControlVar->Fill(controlVar,weight);
			}  
			TH1Data_ControlVar->Fill(controlVar,weight);
			TH1Data_Pt_Control->Fill(controlVar1,weight);
									
		}
	}  
	
	//if((controlVar>lowCutVar0 && controlVar<centralLowCutVar0) || (controlVar>=centralUpCutVar0 && controlVar<upCutVar0)){
	
	if(controlVar1>ptbinlow_ && controlVar1<ptbinup_){
		if(controlVar2>etabinlow_ && controlVar2<etabinup_){
			if((controlVar>lowCutVar0 && controlVar<centralLowCutVar0) || (controlVar>=centralUpCutVar0 && controlVar<upCutVar0)){
				
				if(fabs(partonFlavour)==5){
					TH1Sng_chisqControl->Fill(chisq);
				} else {
					TH1Bkg_chisqControl->Fill(chisq);
				}
				TH1Data_chisqControl->Fill(chisq);
			}
		}
	}
	
	
	if(fabs(partonFlavour)==5){
		TH2Sng_ControlVar12->Fill(controlVar1,controlVar2,weight);
	} else {
		TH2Bkg_ControlVar12->Fill(controlVar1,controlVar2,weight);
		if(isW) TH2Bkg_W_ControlVar12->Fill(controlVar1,controlVar2,weight);
		if(isR) TH2Bkg_R_ControlVar12->Fill(controlVar1,controlVar2,weight);
	}
	
	//if(isW) TH2Data_ControlVar12->Fill(controlVar1,controlVar2,weight);//in the case I want to do the reweighting for a certain component only (eg W or R or W&R)
	TH2Data_ControlVar12->Fill(controlVar1,controlVar2,weight);
    //cout << controlVar1 << " " << controlVar2 << " " << weight << endl << endl << endl;
	//}
	
}

void PtEtaBin::FillXStemplates(double weight, string dataSetName, int partonFlavour, double btag, double *btagCuts, double controlVar0, double lowCutVar0, double centralLowCutVar0, double centralUpCutVar0, double upCutVar0) {

	//cout << weight << endl;
	
	/*cout << "PtEtaBin::FillXStemplates::bTag " << btag << endl;
	
	cout << "PtEtaBin::FillXStemplates::bTagCuts[0] " << btagCuts[0] << endl;
	cout << "PtEtaBin::FillXStemplates::bTagCuts[1] " << btagCuts[1] << endl;
	cout << "PtEtaBin::FillXStemplates::bTagCuts[2] " << btagCuts[2] << endl;*/
	
	/*for (unsigned int i=0; i<3; i++) { 
		
		//cout << i << endl;
		
		stringstream cut; cut << btagCuts[i];
		
		string title = "TH1Data_Var0_btagCut_"+cut.str()+"_VVMC";
		
		if (histo1D.find(title) == histo1D.end()) {
			cout << "uhu?" << endl;
			TString name = ""; GiveName(&name); name+="TH1Data_Var0_VVMC_btagCut"+cut.str(); 
			histo1D[title] = new TH1D(name,name,histo1D["TH1Data_Var0_VVMC_btagCut"]->GetNbinsX(),histo1D["TH1Data_Var0_VVMC_btagCut"]->GetBinLowEdge(1),histo1D["TH1Data_Var0_VVMC_btagCut"]->GetBinLowEdge(histo1D["TH1Data_Var0_VVMC_btagCut"]->GetNbinsX())); 
			//histo1D[title]->SetNameTitle(name,name);
		}
													   }*/
	
	//exit(1);
		
	if(dataSetName != "Data" && dataSetName != "data" && dataSetName != "DATA") {

		if(dataSetName.find("TTbarJets") != 0) {
			histo1D["TH1Data_Var0_VVMC"]->Fill(controlVar0,weight);
			if (fabs(partonFlavour)==5)
				histo1D["TH1Sng_Var0_VVMC"]->Fill(controlVar0,weight);
			else 
				histo1D["TH1Bkg_Var0_VVMC"]->Fill(controlVar0,weight);
			//cout << "VVJETS " << dataSetName << endl;
			
		}
		
		if(dataSetName.find("TTbarJets_Semi") != 0) {
			histo1D["TH1Data_Var0_VVTTotherMC"]->Fill(controlVar0,weight);	
		}
		
		if(dataSetName.find("TTbarJets") == 0) {
			
			nTTbar_+=weight;

			//cout << "pom " << controlVar0 << " " << weight << endl;
			histo1D["TH1Data_Var0_TTbar"]->Fill(controlVar0,weight);
			if (fabs(partonFlavour)==5)
				histo1D["TH1Sng_Var0_TTbar"]->Fill(controlVar0,weight);
			else 
				histo1D["TH1Bkg_Var0_TTbar"]->Fill(controlVar0,weight);

		}
		
		if(dataSetName.find("TTbarJets_Semi") == 0) {
			histo1D["TH1Data_Var0_TTbarSemiMu"]->Fill(controlVar0,weight);
		}
		
		if(dataSetName.find("TTbarJets_Other") == 0) {
			histo1D["TH1Data_Var0_TTbarOther"]->Fill(controlVar0,weight);
		}
		
		if (btag > btagCuts[0]) {
			
			//cout << btag << " Using cut " << btagCuts[0] << endl;
			 
			if(dataSetName.find("TTbarJets") == 0) {
				histo1D["TH1Data_Var0_TTbar_bTagL"]->Fill(controlVar0,weight);
				if (fabs(partonFlavour)==5)
					histo1D["TH1Sng_Var0_TTbar_bTagL"]->Fill(controlVar0,weight);
				else 
					histo1D["TH1Bkg_Var0_TTbar_bTagL"]->Fill(controlVar0,weight);				
			}
			else {
				histo1D["TH1Data_Var0_VVMC_bTagL"]->Fill(controlVar0,weight);
				if (fabs(partonFlavour)==5)
					histo1D["TH1Sng_Var0_VVMC_bTagL"]->Fill(controlVar0,weight);
				else 
					histo1D["TH1Bkg_Var0_VVMC_bTagL"]->Fill(controlVar0,weight);
			}

		}
		
		if (btag > btagCuts[1]) {
			
			if(dataSetName.find("TTbarJets") == 0) {
				histo1D["TH1Data_Var0_TTbar_bTagM"]->Fill(controlVar0,weight);
				if (fabs(partonFlavour)==5)
					histo1D["TH1Sng_Var0_TTbar_bTagM"]->Fill(controlVar0,weight);
				else 
					histo1D["TH1Bkg_Var0_TTbar_bTagM"]->Fill(controlVar0,weight);				
			}
			else {
				histo1D["TH1Data_Var0_VVMC_bTagM"]->Fill(controlVar0,weight);
				if (fabs(partonFlavour)==5)
					histo1D["TH1Sng_Var0_VVMC_bTagM"]->Fill(controlVar0,weight);
				else 
					histo1D["TH1Bkg_Var0_VVMC_bTagM"]->Fill(controlVar0,weight);
			}
		}
		
		if (btag > btagCuts[2]) {
						
			if(dataSetName.find("TTbarJets") == 0) {
				histo1D["TH1Data_Var0_TTbar_bTagT"]->Fill(controlVar0,weight);
				if (fabs(partonFlavour)==5)
					histo1D["TH1Sng_Var0_TTbar_bTagT"]->Fill(controlVar0,weight);
				else 
					histo1D["TH1Bkg_Var0_TTbar_bTagT"]->Fill(controlVar0,weight);				
			}
			else {
				histo1D["TH1Data_Var0_VVMC_bTagT"]->Fill(controlVar0,weight);
			if (fabs(partonFlavour)==5)
				histo1D["TH1Sng_Var0_VVMC_bTagL"]->Fill(controlVar0,weight);
			else 
				histo1D["TH1Bkg_Var0_VVMC_bTagL"]->Fill(controlVar0,weight);
			}
		}
		
	}
	
}

void PtEtaBin::MakeSoverSBPlots(){
 
  TH1D *TH1temp = new TH1D("temp","temp",TH1Sng_Var0->GetXaxis()->GetNbins(),TH1Sng_Var0->GetXaxis()->GetBinLowEdge(1),TH1Sng_Var0->GetXaxis()->GetBinUpEdge(TH1Sng_Var0->GetXaxis()->GetNbins())); 
     
  TH1temp->Add(TH1Sng_Var0,TH1Bkg_Var0);
  TH1SoverSB_Var0->Divide(TH1Sng_Var0,TH1temp);

  TH1D *TH1tempControl = new TH1D("tempcontrol","tempcontrol",TH1Sng_ControlVar->GetXaxis()->GetNbins(),TH1Sng_ControlVar->GetXaxis()->GetBinLowEdge(1),TH1Sng_ControlVar->GetXaxis()->GetBinUpEdge(TH1Sng_ControlVar->GetXaxis()->GetNbins())); 

  TH1tempControl->Add(TH1Sng_ControlVar,TH1Bkg_ControlVar);
  TH1SoverSB_ControlVar->Divide(TH1Sng_ControlVar,TH1tempControl);

  delete TH1temp;
  delete TH1tempControl;

}

void PtEtaBin::MakeMCEffPlots(){
	
	GiveName(&titleSng_BtagEffMC_); titleSng_BtagEffMC_+="TH1Sng_BtagEffMC"; 
	GiveName(&titleBkg_BtagEffMC_); titleBkg_BtagEffMC_+="TH1Bkg_BtagEffMC"; 
	GiveName(&titleData_BtagEffMC_); titleData_BtagEffMC_+="TH1Data_BtagEffMC"; 
	
	if(!varBinSize_){
		TH1Sng_BtagEffMC = new TH1D(titleSng_BtagEffMC_,titleSng_BtagEffMC_,TH1Sng_BtagAll->GetXaxis()->GetNbins(),TH1Sng_BtagAll->GetXaxis()->GetBinLowEdge(1),TH1Sng_BtagAll->GetXaxis()->GetBinUpEdge(TH1Sng_BtagAll->GetXaxis()->GetNbins()));
		TH1Bkg_BtagEffMC = new TH1D(titleBkg_BtagEffMC_,titleBkg_BtagEffMC_,TH1Bkg_BtagAll->GetXaxis()->GetNbins(),TH1Bkg_BtagAll->GetXaxis()->GetBinLowEdge(1),TH1Bkg_BtagAll->GetXaxis()->GetBinUpEdge(TH1Bkg_BtagAll->GetXaxis()->GetNbins()));
		TH1Data_BtagEffMC = new TH1D(titleData_BtagEffMC_,titleData_BtagEffMC_,TH1Data_BtagAll->GetXaxis()->GetNbins(),TH1Data_BtagAll->GetXaxis()->GetBinLowEdge(1),TH1Data_BtagAll->GetXaxis()->GetBinUpEdge(TH1Data_BtagAll->GetXaxis()->GetNbins()));
	}
	
	if(varBinSize_){
		TH1Sng_BtagEffMC = new TH1D(titleSng_BtagEffMC_,titleSng_BtagEffMC_,TH2Sng_Left->GetYaxis()->GetNbins(),rangesbTag_);
		TH1Bkg_BtagEffMC = new TH1D(titleBkg_BtagEffMC_,titleBkg_BtagEffMC_,TH2Bkg_Left->GetYaxis()->GetNbins(),rangesbTag_);
		TH1Data_BtagEffMC = new TH1D(titleData_BtagEffMC_,titleData_BtagEffMC_,TH2Data_Left->GetYaxis()->GetNbins(),rangesbTag_);
	}
		
	TH1Sng_BtagEffMC->Sumw2();
	TH1Bkg_BtagEffMC->Sumw2();
	TH1Data_BtagEffMC->Sumw2();
	
	double integral=0; integral=TH2Sng_Left1DY->Integral(0,TH2Sng_Left1DY->GetXaxis()->GetNbins()+1);
	
	//effcalculation
	
	double min=0;
	double mine=0;
	double plus=0;
	double pluse=0;
	for(int i=0; i<=TH2Sng_Left1DY->GetXaxis()->GetNbins()+1; i++){
		plus+=TH2Sng_Left1DY->GetBinContent(i);
		pluse+=pow(TH2Sng_Left1DY->GetBinError(i),2);
	}
	
	TH1Sng_BtagEffMC->SetBinContent(0,1);
	for(int i=0; i<=TH2Sng_Left1DY->GetXaxis()->GetNbins(); i++){
		double Eff=TH2Sng_Left1DY->Integral(i,TH2Sng_Left1DY->GetXaxis()->GetNbins()+1)/integral;//The corresponding cut is the low edge of the bin
		//TH1Sng_BtagEffMC->SetBinContent(i,Eff); 
		
		plus-=TH2Sng_Left1DY->GetBinContent(i);
		pluse-=pow(TH2Sng_Left1DY->GetBinError(i),2);
		if(plus<1){plus=0; pluse=0;}
		min+=TH2Sng_Left1DY->GetBinContent(i);
		mine+=pow(TH2Sng_Left1DY->GetBinError(i),2);
		
		double tmpcontent=plus/(plus+min);
		
		double tmperror=sqrt(
							 pow(sqrt(pluse) * ( 1/(plus+min) - plus/( (plus+min)*(plus+min) )  ) ,2) + 
							 pow(sqrt(mine)  *                  plus/( (plus+min)*(plus+min) )    ,2)
							 );
		
		//cout << " MC : " << titleSng_BtagEffMC_ << " " << tmpcontent << " " << Eff << endl;
		
		TH1Sng_BtagEffMC->SetBinContent(i+1,tmpcontent); 
		TH1Sng_BtagEffMC->SetBinError(i+1,tmperror); 
	}
	
	integral=0;
	integral=TH2Bkg_Left1DY->Integral(0,TH2Bkg_Left1DY->GetXaxis()->GetNbins()+1);
	for(int i=0; i<=TH2Bkg_Left1DY->GetXaxis()->GetNbins()+1; i++){
		double Eff=TH2Bkg_Left1DY->Integral(i,TH2Bkg_Left1DY->GetXaxis()->GetNbins()+1)/integral;
		TH1Bkg_BtagEffMC->SetBinContent(i,Eff); 
		TH1Bkg_BtagEffMC->SetBinError(i,0); 
		//The corresponding cut is the low edge of the bin
	}
	integral=0;
	integral=TH2Data_Left1DY->Integral(0,TH2Data_Left1DY->GetXaxis()->GetNbins()+1);
	for(int i=0; i<=TH2Data_Left1DY->GetXaxis()->GetNbins()+1; i++){
		double Eff=TH2Data_Left1DY->Integral(i,TH2Data_Left1DY->GetXaxis()->GetNbins()+1)/integral;
		TH1Data_BtagEffMC->SetBinContent(i,Eff); 
		TH1Data_BtagEffMC->SetBinError(i,0); 
		//The corresponding cut is the low edge of the bin
	}
	integral=0;
	
	integral=0; 
	integral=TH1Sng_BtagAll->Integral(0,TH1Sng_BtagAll->GetXaxis()->GetNbins()+1); //the plots here are for the distributions of the whole signal sample  
	
	double minA=0;
	double minAe=0;
	double plusA=0;
	double plusAe=0;  
	for(int i=0; i<=TH1Sng_BtagAll->GetXaxis()->GetNbins()+1; i++){
		plusA+=TH1Sng_BtagAll->GetBinContent(i);
		plusAe+=pow(TH1Sng_BtagAll->GetBinError(i),2);
	} 
	TH1Sng_BtagEffAll->SetBinContent(0,1);
	for(int i=0; i<=TH1Sng_BtagAll->GetXaxis()->GetNbins(); i++){
		double Eff=TH1Sng_BtagAll->Integral(i,TH1Sng_BtagAll->GetXaxis()->GetNbins()+1)/integral;    //The corresponding cut is the low edge of the bin
		//TH1Sng_BtagEffAll->SetBinContent(i,Eff); 
		
		plusA-=TH1Sng_BtagAll->GetBinContent(i);
		plusAe-=pow(TH1Sng_BtagAll->GetBinError(i),2);
		if(plusA<1){plusA=0; plusAe=0;}
		minA+=TH1Sng_BtagAll->GetBinContent(i);
		minAe+=pow(TH1Sng_BtagAll->GetBinError(i),2);
		
		double tmpcontent=plusA/(plusA+minA);
		
		double tmperror=sqrt(
							 pow(sqrt(plusAe) * ( 1/(plusA+minA) - plusA/( (plusA+minA)*(plusA+minA) )  ) ,2) + 
							 pow(sqrt(minAe)  *                  plusA/( (plusA+minA)*(plusA+minA) )    ,2)
							 );
		
		TH1Sng_BtagEffAll->SetBinContent(i+1,tmpcontent); 
		TH1Sng_BtagEffAll->SetBinError(i+1,tmperror); 
		
		//cout << " All : " << titleSng_BtagEffMC_ << " " << tmpcontent << " " << Eff << endl;
		
	}
	
	integral=0;
	integral=TH1Bkg_BtagAll->Integral(0,TH1Bkg_BtagAll->GetXaxis()->GetNbins()+1);
	for(int i=0; i<=TH1Bkg_BtagAll->GetXaxis()->GetNbins()+1; i++){
		double Eff=TH1Bkg_BtagAll->Integral(i,TH1Bkg_BtagAll->GetXaxis()->GetNbins()+1)/integral;
		TH1Bkg_BtagEffAll->SetBinContent(i,Eff); 
		TH1Bkg_BtagEffAll->SetBinError(i,0); 
		//The corresponding cut is the low edge of the bin
	}
	integral=0;  
	integral=TH1Data_BtagAll->Integral(0,TH1Data_BtagAll->GetXaxis()->GetNbins()+1);
	for(int i=0; i<=TH1Data_BtagAll->GetXaxis()->GetNbins()+1; i++){
		double Eff=TH1Data_BtagAll->Integral(i,TH1Data_BtagAll->GetXaxis()->GetNbins()+1)/integral;
		TH1Data_BtagEffAll->SetBinContent(i,Eff); 
		TH1Data_BtagEffAll->SetBinError(i,0); 
		//The corresponding cut is the low edge of the bin
	}
	integral=0;
}

void PtEtaBin::MakeProfileXplots(){
  
  GiveName(&titleSng_ProfileX_); titleSng_ProfileX_+="TH2Sng_ProfileX"; 
  GiveName(&titleBkg_ProfileX_); titleBkg_ProfileX_+="TH2Bkg_ProfileX"; 
  GiveName(&titleData_ProfileX_); titleData_ProfileX_+="TH2Data_ProfileX"; 
  //if(debug_>0) cout << "WARNING::PtEtaBin::MakeProfileXplots() - this gives only a reasonable result for tagger nr.0 and for fixed bins where one bin has edge 0!" << endl;
  int lowBin=0;//this is only useful for the 0th tagger 
  int highBin=0;//this is only useful for the 0th tagger 
  for(int i=0; i<TH2Sng->GetYaxis()->GetNbins()+1; i++){
    if(TH2Sng->GetYaxis()->GetBinLowEdge(i)==0) {
      //cout << "PtEtaBin::MakeProfileXplots() Found bin low edge " << i << endl;
      lowBin=i-1;
    }
    if(TH2Sng->GetYaxis()->GetBinLowEdge(i)==10) {
      cout << "PtEtaBin::MakeProfileXplots() Found bin low edge " << i << endl;
      highBin=i-1;
    }
  }

  TH2Sng_ProfileX = TH2Sng->ProfileX(titleSng_ProfileX_,lowBin,TH2Sng->GetYaxis()->GetNbins()+1,"e");
  TH2Bkg_ProfileX = TH2Bkg->ProfileX(titleBkg_ProfileX_,lowBin,TH2Bkg->GetYaxis()->GetNbins()+1,"e");
  TH2Data_ProfileX = TH2Data->ProfileX(titleData_ProfileX_,lowBin,TH2Data->GetYaxis()->GetNbins()+1,"e");

  //TH2Sng_ProfileX = TH2Sng->ProfileX(titleSng_ProfileX_,lowBin,highBin,"e");
  //TH2Bkg_ProfileX = TH2Bkg->ProfileX(titleBkg_ProfileX_,lowBin,highBin,"e");
  //TH2Data_ProfileX = TH2Data->ProfileX(titleData_ProfileX_,lowBin,highBin,"e");

}

void PtEtaBin::Make1DXplots(){

  GiveName(&titleSng_Left1DX_); titleSng_Left1DX_+="TH2Sng_Left1DX"; 
  GiveName(&titleBkg_Left1DX_); titleBkg_Left1DX_+="TH2Bkg_Left1DX"; 
  GiveName(&titleData_Left1DX_); titleData_Left1DX_+="TH2Data_Left1DX"; 
  
  TH2Sng_Left1DX = TH2Sng_Left->ProjectionX(titleSng_Left1DX_,0,TH2Sng_Left->GetYaxis()->GetNbins()+1,"e");
  TH2Bkg_Left1DX = TH2Bkg_Left->ProjectionX(titleBkg_Left1DX_,0,TH2Bkg_Left->GetYaxis()->GetNbins()+1,"e");
  TH2Data_Left1DX = TH2Data_Left->ProjectionX(titleData_Left1DX_,0,TH2Data_Left->GetYaxis()->GetNbins()+1,"e");

  GiveName(&titleSng_Right1DX_); titleSng_Right1DX_+="TH2Sng_Right1DX"; 
  GiveName(&titleBkg_Right1DX_); titleBkg_Right1DX_+="TH2Bkg_Right1DX"; 
  GiveName(&titleData_Right1DX_); titleData_Right1DX_+="TH2Data_Right1DX"; 
    
  TH2Sng_Right1DX = TH2Sng_Right->ProjectionX(titleSng_Right1DX_,0,TH2Sng_Right->GetYaxis()->GetNbins()+1,"e");
  TH2Bkg_Right1DX = TH2Bkg_Right->ProjectionX(titleBkg_Right1DX_,0,TH2Bkg_Right->GetYaxis()->GetNbins()+1,"e");
  TH2Data_Right1DX = TH2Data_Right->ProjectionX(titleData_Right1DX_,0,TH2Data_Right->GetYaxis()->GetNbins()+1,"e");

}

void PtEtaBin::MakeSCprojectionPlots(){

  GiveName(&titleSng_Var12_1DX_); titleSng_Var12_1DX_+="TH2Sng_Var12_1DX_"; 
  GiveName(&titleBkg_Var12_1DX_); titleBkg_Var12_1DX_+="TH2Bkg_Var12_1DX_"; 
  GiveName(&titleBkg_W_Var12_1DX_); titleBkg_W_Var12_1DX_+="TH2Bkg_W_Var12_1DX_"; 
  GiveName(&titleBkg_R_Var12_1DX_); titleBkg_R_Var12_1DX_+="TH2Bkg_R_Var12_1DX_"; 
  GiveName(&titleData_Var12_1DX_); titleData_Var12_1DX_+="TH2Data_Var12_1DX_"; 
  GiveName(&titleSng_Var12_1DY_); titleSng_Var12_1DY_+="TH2Sng_Var12_1DY_"; 
  GiveName(&titleBkg_Var12_1DY_); titleBkg_Var12_1DY_+="TH2Bkg_Var12_1DY_"; 
  GiveName(&titleBkg_W_Var12_1DY_); titleBkg_W_Var12_1DY_+="TH2Bkg_W_Var12_1DY_"; 
  GiveName(&titleBkg_R_Var12_1DY_); titleBkg_R_Var12_1DY_+="TH2Bkg_R_Var12_1DY_"; 
  GiveName(&titleData_Var12_1DY_); titleData_Var12_1DY_+="TH2Data_Var12_1DY_"; 
  GiveName(&titleSng_ControlVar12_1DX_); titleSng_ControlVar12_1DX_+="TH2Sng_ControlVar12_1DX_"; 
  GiveName(&titleBkg_ControlVar12_1DX_); titleBkg_ControlVar12_1DX_+="TH2Bkg_ControlVar12_1DX_"; 
  GiveName(&titleBkg_W_ControlVar12_1DX_); titleBkg_W_ControlVar12_1DX_+="TH2Bkg_W_ControlVar12_1DX_"; 
  GiveName(&titleBkg_R_ControlVar12_1DX_); titleBkg_R_ControlVar12_1DX_+="TH2Bkg_R_ControlVar12_1DX_"; 
  GiveName(&titleData_ControlVar12_1DX_); titleData_ControlVar12_1DX_+="TH2Data_ControlVar12_1DX_"; 
  GiveName(&titleSng_ControlVar12_1DY_); titleSng_ControlVar12_1DY_+="TH2Sng_ControlVar12_1DY_"; 
  GiveName(&titleBkg_ControlVar12_1DY_); titleBkg_ControlVar12_1DY_+="TH2Bkg_ControlVar12_1DY_"; 
  GiveName(&titleBkg_W_ControlVar12_1DY_); titleBkg_W_ControlVar12_1DY_+="TH2Bkg_W_ControlVar12_1DY_"; 
  GiveName(&titleBkg_R_ControlVar12_1DY_); titleBkg_R_ControlVar12_1DY_+="TH2Bkg_R_ControlVar12_1DY_"; 
  GiveName(&titleData_ControlVar12_1DY_); titleData_ControlVar12_1DY_+="TH2Data_ControlVar12_1DY_"; 

  TH2Sng_Var12_1DX = TH2Sng_Var12->ProjectionX(titleSng_Var12_1DX_,0,TH2Sng_Var12->GetYaxis()->GetNbins()+1,"e");
  TH2Bkg_Var12_1DX = TH2Bkg_Var12->ProjectionX(titleBkg_Var12_1DX_,0,TH2Bkg_Var12->GetYaxis()->GetNbins()+1,"e");
  TH2Bkg_W_Var12_1DX = TH2Bkg_W_Var12->ProjectionX(titleBkg_W_Var12_1DX_,0,TH2Bkg_W_Var12->GetYaxis()->GetNbins()+1,"e");
  TH2Bkg_R_Var12_1DX = TH2Bkg_R_Var12->ProjectionX(titleBkg_R_Var12_1DX_,0,TH2Bkg_R_Var12->GetYaxis()->GetNbins()+1,"e");
  TH2Data_Var12_1DX = TH2Data_Var12->ProjectionX(titleData_Var12_1DX_,0,TH2Data_Var12->GetYaxis()->GetNbins()+1,"e");
 
  TH2Sng_Var12_1DY = TH2Sng_Var12->ProjectionY(titleSng_Var12_1DY_,0,TH2Sng_Var12->GetXaxis()->GetNbins()+1,"e");
  TH2Bkg_Var12_1DY = TH2Bkg_Var12->ProjectionY(titleBkg_Var12_1DY_,0,TH2Bkg_Var12->GetXaxis()->GetNbins()+1,"e");
  TH2Bkg_W_Var12_1DY = TH2Bkg_W_Var12->ProjectionY(titleBkg_W_Var12_1DY_,0,TH2Bkg_W_Var12->GetXaxis()->GetNbins()+1,"e");
  TH2Bkg_R_Var12_1DY = TH2Bkg_R_Var12->ProjectionY(titleBkg_R_Var12_1DY_,0,TH2Bkg_R_Var12->GetXaxis()->GetNbins()+1,"e");
  TH2Data_Var12_1DY = TH2Data_Var12->ProjectionY(titleData_Var12_1DY_,0,TH2Data_Var12->GetXaxis()->GetNbins()+1,"e");

  TH2Sng_ControlVar12_1DX = TH2Sng_ControlVar12->ProjectionX(titleSng_ControlVar12_1DX_,0,TH2Sng_ControlVar12->GetYaxis()->GetNbins()+1,"e");
  TH2Bkg_ControlVar12_1DX = TH2Bkg_ControlVar12->ProjectionX(titleBkg_ControlVar12_1DX_,0,TH2Bkg_ControlVar12->GetYaxis()->GetNbins()+1,"e");
  TH2Bkg_W_ControlVar12_1DX = TH2Bkg_W_ControlVar12->ProjectionX(titleBkg_W_ControlVar12_1DX_,0,TH2Bkg_W_ControlVar12->GetYaxis()->GetNbins()+1,"e");
  TH2Bkg_R_ControlVar12_1DX = TH2Bkg_R_ControlVar12->ProjectionX(titleBkg_R_ControlVar12_1DX_,0,TH2Bkg_R_ControlVar12->GetYaxis()->GetNbins()+1,"e");
  TH2Data_ControlVar12_1DX = TH2Data_ControlVar12->ProjectionX(titleData_ControlVar12_1DX_,0,TH2Data_ControlVar12->GetYaxis()->GetNbins()+1,"e");
 
  TH2Sng_ControlVar12_1DY = TH2Sng_ControlVar12->ProjectionY(titleSng_ControlVar12_1DY_,0,TH2Sng_ControlVar12->GetXaxis()->GetNbins()+1,"e");
  TH2Bkg_ControlVar12_1DY = TH2Bkg_ControlVar12->ProjectionY(titleBkg_ControlVar12_1DY_,0,TH2Bkg_ControlVar12->GetXaxis()->GetNbins()+1,"e");
  TH2Bkg_W_ControlVar12_1DY = TH2Bkg_W_ControlVar12->ProjectionY(titleBkg_W_ControlVar12_1DY_,0,TH2Bkg_W_ControlVar12->GetXaxis()->GetNbins()+1,"e");
  TH2Bkg_R_ControlVar12_1DY = TH2Bkg_R_ControlVar12->ProjectionY(titleBkg_R_ControlVar12_1DY_,0,TH2Bkg_R_ControlVar12->GetXaxis()->GetNbins()+1,"e");
  TH2Data_ControlVar12_1DY = TH2Data_ControlVar12->ProjectionY(titleData_ControlVar12_1DY_,0,TH2Data_ControlVar12->GetXaxis()->GetNbins()+1,"e");

  //TH2Sng_Var12_1DX->Sumw2();//not needed, inherited from 2D histo

}

void PtEtaBin::SetError1DXplots(){

  SetError(TH2Sng_Left1DX);
  SetError(TH2Bkg_Left1DX);
  SetError(TH2Data_Left1DX);
  SetError(TH2Sng_Right1DX);
  SetError(TH2Bkg_Right1DX);
  SetError(TH2Data_Right1DX);

}

void PtEtaBin::Make2DRatioPlot(){

  double rightScale = TH2DataVar12_Right->Integral(0,TH2DataVar12_Right->GetXaxis()->GetNbins()+1,0,TH2DataVar12_Right->GetYaxis()->GetNbins()+1);
  double leftScale  = TH2DataVar12_Left->Integral(0,TH2DataVar12_Left->GetXaxis()->GetNbins()+1,0,TH2DataVar12_Left->GetYaxis()->GetNbins()+1);
  
  TH2DataVar12_Right->Scale(1/rightScale); // Todo: check explicitely that the errors are propagated correctly.
  TH2DataVar12_Left->Scale(1/leftScale);

  GiveName(&titleDataVar12_LeftRight_); titleDataVar12_LeftRight_+="TH2DataVar12_LeftRight"; 
  TH2DataVar12_LeftRight = new TH2D(titleDataVar12_LeftRight_,titleDataVar12_LeftRight_,TH2DataVar12_Right->GetXaxis()->GetNbins(),TH2DataVar12_Right->GetXaxis()->GetBinLowEdge(1),TH2DataVar12_Right->GetXaxis()->GetBinUpEdge(TH2DataVar12_Right->GetXaxis()->GetNbins()),TH2DataVar12_Right->GetYaxis()->GetNbins(),TH2DataVar12_Right->GetYaxis()->GetBinLowEdge(1),TH2DataVar12_Right->GetYaxis()->GetBinUpEdge(TH2DataVar12_Right->GetYaxis()->GetNbins()));

  for(int i=0; i<=TH2DataVar12_Right->GetXaxis()->GetNbins()+1; i++){ 
    for(int j=0; j<=TH2DataVar12_Right->GetYaxis()->GetNbins()+1; j++){ 
      if(TH2DataVar12_Left->GetBinContent(i,j)!=0 && TH2DataVar12_Right->GetBinContent(i,j)!=0){
	TH2DataVar12_LeftRight->SetBinContent(i,j,TH2DataVar12_Left->GetBinContent(i,j)/TH2DataVar12_Right->GetBinContent(i,j));
	
	//still to be implemented
	//      TH2Data_RightLeft1DX->SetBinError(i,sqrt(pow((TH2Data_Right1DX->GetBinError(i)/TH2Data_Left1DX->GetBinContent(i)),2)+pow((TH2Data_Right1DX->GetBinContent(i)/TH2Data_Left1DX->GetBinError(i)*TH2Data_Left1DX->GetBinError(i)),2)));
      } else {   
	TH2DataVar12_LeftRight->SetBinContent(i,j,-0.01);
      }
    }
  }
  
 
  TH2DataVar12_Right->Scale(rightScale);
  TH2DataVar12_Left->Scale(leftScale);
  
  /*double rightScale = TH2BkgVar12_Right->Integral(0,TH2BkgVar12_Right->GetXaxis()->GetNbins()+1,0,TH2BkgVar12_Right->GetYaxis()->GetNbins()+1);
  double leftScale  = TH2BkgVar12_Left->Integral(0,TH2BkgVar12_Left->GetXaxis()->GetNbins()+1,0,TH2BkgVar12_Left->GetYaxis()->GetNbins()+1);
  
  TH2BkgVar12_Right->Scale(1/rightScale); // Todo: check explicitely that the errors are propagated correctly.
  TH2BkgVar12_Left->Scale(1/leftScale);

  GiveName(&titleDataVar12_LeftRight_); titleDataVar12_LeftRight_+="TH2DataVar12_LeftRight"; 
  TH2DataVar12_LeftRight = new TH2D(titleDataVar12_LeftRight_,titleDataVar12_LeftRight_,TH2BkgVar12_Right->GetXaxis()->GetNbins(),TH2BkgVar12_Right->GetXaxis()->GetBinLowEdge(1),TH2BkgVar12_Right->GetXaxis()->GetBinUpEdge(TH2BkgVar12_Right->GetXaxis()->GetNbins()),TH2BkgVar12_Right->GetYaxis()->GetNbins(),TH2BkgVar12_Right->GetYaxis()->GetBinLowEdge(1),TH2BkgVar12_Right->GetYaxis()->GetBinUpEdge(TH2BkgVar12_Right->GetYaxis()->GetNbins()));

  for(int i=0; i<=TH2BkgVar12_Right->GetXaxis()->GetNbins()+1; i++){ 
    for(int j=0; j<=TH2BkgVar12_Right->GetYaxis()->GetNbins()+1; j++){ 
      if(TH2BkgVar12_Left->GetBinContent(i,j)!=0 && TH2BkgVar12_Right->GetBinContent(i,j)!=0){
	TH2DataVar12_LeftRight->SetBinContent(i,j,TH2BkgVar12_Left->GetBinContent(i,j)/TH2BkgVar12_Right->GetBinContent(i,j));
	
	//still to be implemented
	//      TH2Data_RightLeft1DX->SetBinError(i,sqrt(pow((TH2Bkg_Right1DX->GetBinError(i)/TH2Bkg_Left1DX->GetBinContent(i)),2)+pow((TH2Bkg_Right1DX->GetBinContent(i)/TH2Bkg_Left1DX->GetBinError(i)*TH2Bkg_Left1DX->GetBinError(i)),2)));
      } else {   
	TH2DataVar12_LeftRight->SetBinContent(i,j,-0.01);
      }
    }
  }
  
  TH2BkgVar12_Right->Scale(rightScale);
  TH2BkgVar12_Left->Scale(leftScale);*/


}

void PtEtaBin::MakeSCVar12RatioPlot(){
 

  if(TH2Data_Var12->GetXaxis()->GetNbins() !=  TH2Data_ControlVar12->GetXaxis()->GetNbins()) cout << "WARNING::PtEtaBin::MakeSCVar12RatioPlot WARNING:histograms have different number of bins" << endl;

  GiveName(&titleData_SignalControlVar12_); titleData_SignalControlVar12_+="TH2Data_SignalControlVar12_";
  GiveName(&titleData_SignalControlVar12_1DX_); titleData_SignalControlVar12_1DX_+="TH2Data_SignalControlVar12_1DX_";
  GiveName(&titleData_SignalControlVar12_1DY_); titleData_SignalControlVar12_1DY_+="TH2Data_SignalControlVar12_1DY_";
    
  /*TH2Data_SignalControlVar12 = new TH2D(titleData_SignalControlVar12_,titleData_SignalControlVar12_,
					TH2Data_Var12->GetXaxis()->GetNbins(),
					TH2Data_Var12->GetXaxis()->GetBinLowEdge(1),
					TH2Data_Var12->GetXaxis()->GetBinUpEdge(TH2Data_Var12->GetXaxis()->GetNbins()),
					TH2Data_Var12->GetYaxis()->GetNbins(),
					TH2Data_Var12->GetYaxis()->GetBinLowEdge(1),
					TH2Data_Var12->GetYaxis()->GetBinUpEdge(TH2Data_Var12->GetYaxis()->GetNbins()));
  
  TH2Data_SignalControlVar12_1DX = new TH1D(titleData_SignalControlVar12_1DX_,titleData_SignalControlVar12_1DX_,
					    TH2Data_Var12->GetXaxis()->GetNbins(),
					    TH2Data_Var12->GetXaxis()->GetBinLowEdge(1),
					    TH2Data_Var12->GetXaxis()->GetBinUpEdge(TH2Data_Var12->GetXaxis()->GetNbins()));

  TH2Data_SignalControlVar12_1DY = new TH1D(titleData_SignalControlVar12_1DY_,titleData_SignalControlVar12_1DY_,
					    TH2Data_Var12->GetYaxis()->GetNbins(),
					    TH2Data_Var12->GetYaxis()->GetBinLowEdge(1),
					    TH2Data_Var12->GetYaxis()->GetBinUpEdge(TH2Data_Var12->GetYaxis()->GetNbins()));
  */

  
  //  double rangesPtSC[21] = {30,33.9517,38.0413,42.3827,47.1889,52.4274,58.0206,63.9226,70.6819,77.5725,85.4028,94.2801,104.116,115.42,128.936,145.192,165.937,194.182,237.506,316.667,500};
  double rangesPtSC[19] = {30,33.9517,38.0413,42.3827,47.1889,52.4274,58.0206,63.9226,70.6819,77.5725,85.4028,94.2801,104.116,115.42,128.936,145.192,165.937,250,9999.};
  
  //double rangesPtSC[19] = {30,50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350,370,400};
  //double rangesPtSC[51] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400,410,420,430,440,450,460,470,480,490,500};
  
  
  //double rangesEtaSC[11] = {0,0.24,0.48,0.72,0.96,1.20,1.44,1.68,1.92,2.16,2.4};
  double rangesEtaSC[21] = {0,0.12,0.24,0.36,0.48,0.60,0.72,0.84,0.96,1.08,1.20,1.32,1.44,1.56,1.68,1.80,1.92,2.04,2.16,2.28,2.4};
  //double rangesEtaSC[21] = {0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1};
  //double rangesEtaSC[21] = {-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
  //double rangesEtaSC[21] = {0,0.175,0.175*2,0.175*3,0.175*4,0.175*5,0.175*6,0.175*7,0.175*8,0.175*9,0.175*10,0.175*11,0.175*12,0.175*13,0.175*14,0.175*15,0.175*16,0.175*17,0.175*18,0.175*19,0.175*20};
  
  
  TH2Data_SignalControlVar12 = new TH2D(titleData_SignalControlVar12_,titleData_SignalControlVar12_,
					TH2Data_Var12->GetXaxis()->GetNbins(),
					rangesPtSC,
					TH2Data_Var12->GetYaxis()->GetNbins(),
					rangesEtaSC);
  
  
  TH2Data_SignalControlVar12_1DX = new TH1D(titleData_SignalControlVar12_1DX_,titleData_SignalControlVar12_1DX_,
					    TH2Data_Var12->GetXaxis()->GetNbins(),
					    rangesPtSC);
  
  TH2Data_SignalControlVar12_1DY = new TH1D(titleData_SignalControlVar12_1DY_,titleData_SignalControlVar12_1DY_,
					    TH2Data_Var12->GetYaxis()->GetNbins(),
					    rangesEtaSC);
  

  TH2Data_SignalControlVar12->Sumw2();
  TH2Data_SignalControlVar12_1DX->Sumw2();
  TH2Data_SignalControlVar12_1DY->Sumw2();

  //Start:: This part is for the data driven reweigh
  double signalScale=TH2Data_Var12->Integral(0,TH2Data_Var12->GetXaxis()->GetNbins()+1,0,TH2Data_Var12->GetYaxis()->GetNbins()+1);
  double controlScale=TH2Data_ControlVar12->Integral(0,TH2Data_ControlVar12->GetXaxis()->GetNbins()+1,0,TH2Data_ControlVar12->GetXaxis()->GetNbins()+1);
    
  double signalScaleY=TH2Data_Var12_1DY->Integral(0,TH2Data_Var12_1DY->GetXaxis()->GetNbins()+1);
  double controlScaleY=TH2Data_ControlVar12_1DY->Integral(0,TH2Data_ControlVar12_1DY->GetXaxis()->GetNbins()+1);
    
  double signalScaleX=TH2Data_Var12_1DX->Integral(0,TH2Data_Var12_1DX->GetXaxis()->GetNbins()+1);
  double controlScaleX=TH2Data_ControlVar12_1DX->Integral(0,TH2Data_ControlVar12_1DX->GetXaxis()->GetNbins()+1);
    
  TH2Data_Var12->Scale(1/signalScale);
  TH2Data_ControlVar12->Scale(1/controlScale);

  TH2Data_Var12_1DX->Scale(1/signalScaleX);
  TH2Data_ControlVar12_1DX->Scale(1/controlScaleX);

  TH2Data_Var12_1DY->Scale(1/signalScaleY);
  TH2Data_ControlVar12_1DY->Scale(1/controlScaleY);


  //cout << genericName << "  " << controlScaleX << "  " << controlScaleY << endl;

  TH2Data_SignalControlVar12->Divide(TH2Data_Var12,TH2Data_ControlVar12);
  TH2Data_SignalControlVar12_1DX->Divide(TH2Data_Var12_1DX,TH2Data_ControlVar12_1DX);
  TH2Data_SignalControlVar12_1DY->Divide(TH2Data_Var12_1DY,TH2Data_ControlVar12_1DY);

  /*  int ifix=-1;
  
  for(int i=0; i<=TH2Data_SignalControlVar12->GetXaxis()->GetNbins()+1; i++){ 
    for(int j=0; j<=TH2Data_SignalControlVar12->GetYaxis()->GetNbins()+1; j++){ 
      
      if(TH2Data_SignalControlVar12->GetXaxis()->GetBinUpEdge(i)>200) {
	if(ifix<0) ifix=i-1;
	TH2Data_SignalControlVar12->SetBinContent(i,j,TH2Data_SignalControlVar12->GetBinContent(ifix,j));
      }
      
    }
    }*/

  TH2Data_Var12->Scale(signalScale);
  TH2Data_ControlVar12->Scale(controlScale);
 
  TH2Data_Var12_1DX->Scale(signalScaleX);
  TH2Data_ControlVar12_1DX->Scale(controlScaleX);
  
  TH2Data_Var12_1DY->Scale(signalScaleY);
  TH2Data_ControlVar12_1DY->Scale(controlScaleY);
  //End:: This part is for the data driven reweigh

  //-----------

  //Start:: This part is for the reweigh of the signal between different Chi Square cuts
  /*double nominalScale=TH2Sng_Var12->Integral(0,TH2Sng_Var12->GetXaxis()->GetNbins()+1,0,TH2Sng_Var12->GetYaxis()->GetNbins()+1);
  double cuttedScale=TH2Bkg_Var12->Integral(0,TH2Bkg_Var12->GetXaxis()->GetNbins()+1,0,TH2Bkg_Var12->GetYaxis()->GetNbins()+1);
    
  double nominalScaleY=TH2Sng_Var12_1DY->Integral(0,TH2Sng_Var12_1DY->GetXaxis()->GetNbins()+1);
  double cuttedScaleY=TH2Bkg_Var12_1DY->Integral(0,TH2Bkg_Var12_1DY->GetXaxis()->GetNbins()+1);
    
  double nominalScaleX=TH2Sng_Var12_1DX->Integral(0,TH2Sng_Var12_1DX->GetXaxis()->GetNbins()+1);
  double cuttedScaleX=TH2Bkg_Var12_1DX->Integral(0,TH2Bkg_Var12_1DX->GetXaxis()->GetNbins()+1);
    
  TH2Sng_Var12->Scale(1/nominalScale);
  TH2Bkg_Var12->Scale(1/cuttedScale);

  TH2Sng_Var12_1DX->Scale(1/nominalScaleX);
  TH2Bkg_Var12_1DX->Scale(1/cuttedScaleX);

  TH2Sng_Var12_1DY->Scale(1/nominalScaleY);
  TH2Bkg_Var12_1DY->Scale(1/cuttedScaleY);
  
  for(int i=0; i<=TH2Sng_Var12->GetXaxis()->GetNbins(); i++){ 
      for(int j=0; j<=TH2Sng_Var12->GetYaxis()->GetNbins(); j++){ 
        if(TH2Sng_Var12->GetBinContent(i,j)!=0 && TH2Bkg_Var12->GetBinContent(i,j)!=0){
	  TH2Data_SignalControlVar12->SetBinContent(i,j,TH2Sng_Var12->GetBinContent(i,j)/TH2Bkg_Var12->GetBinContent(i,j));
	} else {
	TH2Data_SignalControlVar12->SetBinContent(i,j,-0.1);
      }
    }
  }
  
  for(int j=0; j<=TH2Sng_Var12->GetYaxis()->GetNbins(); j++){ 
    if(TH2Sng_Var12_1DY->GetBinContent(j)!=0 && TH2Bkg_Var12_1DY->GetBinContent(j)!=0){
      TH2Data_SignalControlVar12_1DY->SetBinContent(j,TH2Sng_Var12_1DY->GetBinContent(j)/TH2Bkg_Var12_1DY->GetBinContent(j));
    } else {
      TH2Data_SignalControlVar12_1DY->SetBinContent(j,-0.1);
    }
  }

  for(int i=0; i<=TH2Sng_Var12->GetXaxis()->GetNbins(); i++){ 
    if(TH2Sng_Var12_1DX->GetBinContent(i)!=0 && TH2Bkg_Var12_1DX->GetBinContent(i)!=0){
      TH2Data_SignalControlVar12_1DX->SetBinContent(i,TH2Sng_Var12_1DX->GetBinContent(i)/TH2Bkg_Var12_1DX->GetBinContent(i));
    } else {
      TH2Data_SignalControlVar12_1DX->SetBinContent(i,-0.1);
    }
  }
  
  TH2Sng_Var12->Scale(nominalScale);
  TH2Bkg_Var12->Scale(cuttedScale);

  TH2Sng_Var12_1DX->Scale(nominalScaleX);
  TH2Bkg_Var12_1DX->Scale(cuttedScaleX);

  TH2Sng_Var12_1DY->Scale(nominalScaleY);
  TH2Bkg_Var12_1DY->Scale(cuttedScaleY);*/
  //End:: This part is for the reweigh of the signal between different Chi Square cuts

}

void PtEtaBin::MakeXRatioPlot(bool useFit){
  //X stands for Pt

  //I first have to normalize the pt plots, I rescale them to the original ones at the end of this method
  double rightScale=TH2Data_Right1DX->Integral(0,TH2Data_Right1DX->GetXaxis()->GetNbins()+1);
  double leftScale=TH2Data_Left1DX->Integral(0,TH2Data_Left1DX->GetXaxis()->GetNbins()+1);
  
  TH2Data_Right1DX->Scale(1/rightScale); // I've set sumw2 so that the error gets propagated correctly (I've checked this explicitely)
  TH2Data_Left1DX->Scale(1/leftScale);

  GiveName(&titleData_RightLeft1DX_); titleData_RightLeft1DX_+="TH2Data_RightLeft1DX"; 
 
  TH2Data_RightLeft1DX = new TH1D(titleData_RightLeft1DX_,titleData_RightLeft1DX_,TH2Data_Right1DX->GetXaxis()->GetNbins(),TH2Data_Right1DX->GetXaxis()->GetBinLowEdge(1),TH2Data_Right1DX->GetXaxis()->GetBinUpEdge(TH2Data_Right1DX->GetXaxis()->GetNbins()));

  for(int i=0; i<=TH2Data_Right1DX->GetXaxis()->GetNbins()+1; i++){ 
    if(TH2Data_Left1DX->GetBinContent(i)!=0 && TH2Data_Right1DX->GetBinContent(i)!=0){
      TH2Data_RightLeft1DX->SetBinContent(i,TH2Data_Right1DX->GetBinContent(i)/TH2Data_Left1DX->GetBinContent(i));
      TH2Data_RightLeft1DX->SetBinError(i,sqrt(pow((TH2Data_Right1DX->GetBinError(i)/TH2Data_Left1DX->GetBinContent(i)),2)+pow((TH2Data_Right1DX->GetBinContent(i)/TH2Data_Left1DX->GetBinError(i)*TH2Data_Left1DX->GetBinError(i)),2)));
    } else {   
      TH2Data_RightLeft1DX->SetBinContent(i,-0.01);
    }
  }
 
  GiveName(&titleData_LeftRight1DX_); titleData_LeftRight1DX_+="TH2Data_LeftRight1DX"; 
  //if(!varBinSize_){
    TH2Data_LeftRight1DX = new TH1D(titleData_LeftRight1DX_,titleData_LeftRight1DX_,TH2Data_Left1DX->GetXaxis()->GetNbins(),TH2Data_Left1DX->GetXaxis()->GetBinLowEdge(1),TH2Data_Left1DX->GetXaxis()->GetBinUpEdge(TH2Data_Left1DX->GetXaxis()->GetNbins()));
    //}
    //if(varBinSize_){
    //TH2Data_LeftRight1DX = new TH1D(titleData_LeftRight1DX_,titleData_LeftRight1DX_,TH2Data_Left1DX->GetXaxis()->GetNbins(),rangesVar1_);
    //}

  for(int i=0; i<=TH2Data_Left1DX->GetXaxis()->GetNbins()+1; i++){
    if(TH2Data_Left1DX->GetBinContent(i)!=0 && TH2Data_Right1DX->GetBinContent(i)!=0){
      //cout << TH2Data_Left1DX->GetBinContent(i) << " " << TH2Data_Right1DX->GetBinContent(i) << " " << TH2Data_Left1DX->GetBinContent(i)/TH2Data_Right1DX->GetBinContent(i)<< endl;
      TH2Data_LeftRight1DX->SetBinContent(i,TH2Data_Left1DX->GetBinContent(i)/TH2Data_Right1DX->GetBinContent(i)); 
      TH2Data_LeftRight1DX->SetBinError(i,sqrt(pow((TH2Data_Left1DX->GetBinError(i)/TH2Data_Right1DX->GetBinContent(i)),2)+pow((TH2Data_Left1DX->GetBinContent(i)/TH2Data_Right1DX->GetBinError(i)*TH2Data_Right1DX->GetBinError(i)),2)));
    } else {
      TH2Data_LeftRight1DX->SetBinContent(i,-0.01);
    }
  }
  
  TH2Data_Right1DX->Scale(rightScale);
  TH2Data_Left1DX->Scale(leftScale);
 
  //here I'm using the background pt distributions left and right
  /*  double rightScale=TH2Bkg_Right1DX->Integral(0,TH2Bkg_Right1DX->GetXaxis()->GetNbins()+1);
  double leftScale=TH2Bkg_Left1DX->Integral(0,TH2Bkg_Left1DX->GetXaxis()->GetNbins()+1);
  
  TH2Bkg_Right1DX->Scale(1/rightScale); // I've set sumw2 so that the error gets propagated correctly (I've checked this explicitely)
  TH2Bkg_Left1DX->Scale(1/leftScale);

  GiveName(&titleData_RightLeft1DX_); titleData_RightLeft1DX_+="TH2Data_RightLeft1DX"; 
 
  TH2Data_RightLeft1DX = new TH1D(titleData_RightLeft1DX_,titleData_RightLeft1DX_,TH2Bkg_Right1DX->GetXaxis()->GetNbins(),TH2Bkg_Right1DX->GetXaxis()->GetBinLowEdge(1),TH2Bkg_Right1DX->GetXaxis()->GetBinUpEdge(TH2Bkg_Right1DX->GetXaxis()->GetNbins()));

  for(int i=0; i<=TH2Bkg_Right1DX->GetXaxis()->GetNbins()+1; i++){ 
    if(TH2Bkg_Left1DX->GetBinContent(i)!=0 && TH2Bkg_Right1DX->GetBinContent(i)!=0){
      //TH2Data_RightLeft1DX->SetBinContent(i,TH2Bkg_Right1DX->GetBinContent(i)/TH2Bkg_Left1DX->GetBinContent(i));
      TH2Data_RightLeft1DX->SetBinContent(i,0);
      //TH2Data_RightLeft1DX->SetBinError(i,sqrt(pow((TH2Bkg_Right1DX->GetBinError(i)/TH2Bkg_Left1DX->GetBinContent(i)),2)+pow((TH2Bkg_Right1DX->GetBinContent(i)/TH2Bkg_Left1DX->GetBinError(i)*TH2Bkg_Left1DX->GetBinError(i)),2)));
      TH2Data_RightLeft1DX->SetBinError(i,0);
    } else {   
      TH2Data_RightLeft1DX->SetBinContent(i,-0.01);
    }
  }
 
  GiveName(&titleData_LeftRight1DX_); titleData_LeftRight1DX_+="TH2Data_LeftRight1DX"; 
 
  TH2Data_LeftRight1DX = new TH1D(titleData_LeftRight1DX_,titleData_LeftRight1DX_,TH2Bkg_Left1DX->GetXaxis()->GetNbins(),TH2Bkg_Left1DX->GetXaxis()->GetBinLowEdge(1),TH2Bkg_Left1DX->GetXaxis()->GetBinUpEdge(TH2Bkg_Left1DX->GetXaxis()->GetNbins()));

  for(int i=0; i<=TH2Bkg_Left1DX->GetXaxis()->GetNbins()+1; i++){
    if(TH2Bkg_Left1DX->GetBinContent(i)!=0 && TH2Bkg_Right1DX->GetBinContent(i)!=0){
      TH2Data_LeftRight1DX->SetBinContent(i,TH2Bkg_Left1DX->GetBinContent(i)/TH2Bkg_Right1DX->GetBinContent(i)); 
      TH2Data_LeftRight1DX->SetBinError(i,sqrt(pow((TH2Bkg_Left1DX->GetBinError(i)/TH2Bkg_Right1DX->GetBinContent(i)),2)+pow((TH2Bkg_Left1DX->GetBinContent(i)/TH2Bkg_Right1DX->GetBinError(i)*TH2Bkg_Right1DX->GetBinError(i)),2)));
    } else {
      TH2Data_LeftRight1DX->SetBinContent(i,-0.01);
    }
  }

  TH2Bkg_Right1DX->Scale(rightScale);
  TH2Bkg_Left1DX->Scale(leftScale);*/

  //fit the leftright with an exponential function and use the fitted points
 
  if(useFit){
    //TF1 * TFData_LeftRight1DXFit;
    //TString titleData_LeftRight1DX_;
    GiveName(&titleData_LeftRight1DXFit_); titleData_LeftRight1DXFit_+="TH2Data_LeftRight1DXFit"; 
    
    //    TFData_LeftRight1DXFit = new TF1(titleData_LeftRight1DXFit_,"expo",30,300);
    //TFData_LeftRight1DXFit = new TF1(titleData_LeftRight1DXFit_,"exp([0]+[1]*x)+[2]*x",30,300);
    
	TFData_LeftRight1DXFit = new TF1(titleData_LeftRight1DXFit_,"exp([0]+[1]*x)+[2]",31,300);
    TFData_LeftRight1DXFit->SetParameter(0,1.016);
    TFData_LeftRight1DXFit->SetParameter(1,-0.014);
    TFData_LeftRight1DXFit->SetParameter(2,0.25);
	
	//TFData_LeftRight1DXFit = new TF1(titleData_LeftRight1DXFit_,"pol3",30,300);
   
    TH2Data_LeftRight1DX->Fit(TFData_LeftRight1DXFit,"RQ");
    
    /*    cout << "value par0: " << TFData_LeftRight1DXFit->GetParameter(0) << endl;
    cout << "value par1: " << TFData_LeftRight1DXFit->GetParameter(1) << endl;
    cout << "value par2: " << TFData_LeftRight1DXFit->GetParameter(2) << endl;
    */

    //cout << "error par0: " << TFData_LeftRight1DXFit->GetParError(0) << " par1: " << TFData_LeftRight1DXFit->GetParError(1) << endl;
  
    //cout << "2 sigma range par0: " << TFData_LeftRight1DXFit->GetParameter(0)-2*TFData_LeftRight1DXFit->GetParError(0) << " - " << TFData_LeftRight1DXFit->GetParameter(0)+2*TFData_LeftRight1DXFit->GetParError(0) << endl;
    //cout << "2 sigma range par1: " << TFData_LeftRight1DXFit->GetParameter(1)-2*TFData_LeftRight1DXFit->GetParError(1) << " - " << TFData_LeftRight1DXFit->GetParameter(1)+2*TFData_LeftRight1DXFit->GetParError(1) << endl;
  
    /*for(int i=0; i<=TH2Data_LeftRight1DX->GetXaxis()->GetNbins()+1; i++){
      TH2Data_LeftRight1DX->SetBinContent(i,TFData_LeftRight1DXFit->Eval(TH2Data_LeftRight1DX->GetBinCenter(i)));
      }*/
	  
	  //saving reweighting function
	  
	  /*
	   cout << "saving reweigting function" << endl;
	   TFile* r = new TFile("ReweightFunction.root","RECREATE");
	  
	  r->cd();
	  
	  TH2Data_LeftRight1DX->Write();
	  
	  TFData_LeftRight1DXFit->Write();
	  
	  r->Close();*/
	  
	  //loading reweighting function
	  
	  /*for (int i=0; i<100;i++)
	  cout << "loading reweigting function" << endl;

	  TFile* r = new TFile("ReweightFunction.root","READ");
	   
	  r->cd();
		TH2Data_LeftRight1DX = (TH1D*) r->Get("nDisc_10_ptbinlow_0_etabinlow_-9990_TH2Data_LeftRight1DX");
	  	TFData_LeftRight1DXFit = (TF1*) r->Get("nDisc_10_ptbinlow_0_etabinlow_-9990_TH2Data_LeftRight1DXFit");*/
	   
	   //r->Close();
	  
	  
  }
}

void PtEtaBin::GetLeftRightPars(double *par0, double *par1, double *par2){

  // double partmp=TFData_LeftRight1DXFit->GetParameter(0);

  *par0=TFData_LeftRight1DXFit->GetParameter(0);
  *par1=TFData_LeftRight1DXFit->GetParameter(1);
  *par2=TFData_LeftRight1DXFit->GetParameter(2);
  

 
  /*cout << "PtEtaBin::GetLeftRightPars par0: " << *par0 << endl;
  cout << "PtEtaBin::GetLeftRightPars par1: " << *par1 << endl;
  cout << "PtEtaBin::GetLeftRightPars par2: " << *par2 << endl;*/

}

void PtEtaBin::SetLeftRightPars(double *par0, double *par1, double *par2){

  /*cout << "PtEtaBin::SetLeftRightPars par0: " << *par0 << endl;
  cout << "PtEtaBin::SetLeftRightPars par1: " << *par1 << endl;
  cout << "PtEtaBin::SetLeftRightPars par2: " << *par2 << endl;*/

  TFData_LeftRight1DXFit->SetParameter(0,*par0);
  TFData_LeftRight1DXFit->SetParameter(1,*par1);
  TFData_LeftRight1DXFit->SetParameter(2,*par2);

}

void PtEtaBin::Make1DYplots(){
 
  GiveName(&titleSng_Left1DY_); titleSng_Left1DY_+="TH2Sng_Left1DY"; 
  GiveName(&titleBkg_Left1DY_); titleBkg_Left1DY_+="TH2Bkg_Left1DY"; 
  GiveName(&titleData_Left1DY_); titleData_Left1DY_+="TH2Data_Left1DY"; 

  GiveName(&titleSng_Right1DY_); titleSng_Right1DY_+="TH2Sng_Right1DY"; 
  GiveName(&titleBkg_Right1DY_); titleBkg_Right1DY_+="TH2Bkg_Right1DY"; 
  GiveName(&titleData_Right1DY_); titleData_Right1DY_+="TH2Data_Right1DY"; 

  TH2Sng_Right1DY = TH2Sng_Right->ProjectionY(titleSng_Right1DY_,0,TH2Sng_Right->GetXaxis()->GetNbins()+1,"e");
  TH2Bkg_Right1DY = TH2Bkg_Right->ProjectionY(titleBkg_Right1DY_,0,TH2Bkg_Right->GetXaxis()->GetNbins()+1,"e");
  TH2Data_Right1DY = TH2Data_Right->ProjectionY(titleData_Right1DY_,0,TH2Data_Right->GetXaxis()->GetNbins()+1,"e");

  TH2Sng_Left1DY = TH2Sng_Left->ProjectionY(titleSng_Left1DY_,0,TH2Sng_Left->GetXaxis()->GetNbins()+1,"e");
  TH2Bkg_Left1DY = TH2Bkg_Left->ProjectionY(titleBkg_Left1DY_,0,TH2Bkg_Left->GetXaxis()->GetNbins()+1,"e");
  TH2Data_Left1DY = TH2Data_Left->ProjectionY(titleData_Left1DY_,0,TH2Data_Left->GetXaxis()->GetNbins()+1,"e");

  /*  cout << "Sng left " << TH2Sng_Left1DY->Integral(0,TH2Sng_Left1DY->GetXaxis()->GetNbins()+1) << endl;
  cout << "Bkg left " << TH2Bkg_Left1DY->Integral(0,TH2Bkg_Left1DY->GetXaxis()->GetNbins()+1) << endl;
  cout << "Data left " << TH2Data_Left1DY->Integral(0,TH2Data_Left1DY->GetXaxis()->GetNbins()+1) << endl;

  cout << "Sng right " << TH2Sng_Right1DY->Integral(0,TH2Sng_Right1DY->GetXaxis()->GetNbins()+1) << endl;
  cout << "Bkg right " << TH2Bkg_Right1DY->Integral(0,TH2Bkg_Right1DY->GetXaxis()->GetNbins()+1) << endl;
  cout << "Data right " << TH2Data_Right1DY->Integral(0,TH2Data_Right1DY->GetXaxis()->GetNbins()+1) << endl;*/

  //TEMP HACK:
  /*  for(int i=0; i<=TH2Sng_Left1DY->GetXaxis()->GetNbins()+1; i++){
    TH2Sng_Left1DY->SetBinContent(i,0);
  }
  for(int i=0; i<=TH2Bkg_Left1DY->GetXaxis()->GetNbins()+1; i++){
    TH2Bkg_Left1DY->SetBinContent(i,0);
  } 
  for(int i=0; i<=TH2Data_Left1DY->GetXaxis()->GetNbins()+1; i++){
    TH2Data_Left1DY->SetBinContent(i,0);
    }*/
}

void PtEtaBin::SetError1DYplots(){

  SetError(TH2Data_Right1DY);
  SetError(TH2Sng_Right1DY);
  SetError(TH2Bkg_Right1DY);
  SetError(TH2Data_Left1DY);
  SetError(TH2Sng_Left1DY);
  SetError(TH2Bkg_Left1DY);
  
}

void PtEtaBin::Make1DYVar2plots(){

  GiveName(&titleSngVar12_Left1DY_); titleSngVar12_Left1DY_+="TH2SngVar12_Left1DY"; 
  GiveName(&titleBkgVar12_Left1DY_); titleBkgVar12_Left1DY_+="TH2BkgVar12_Left1DY"; 
  GiveName(&titleDataVar12_Left1DY_); titleDataVar12_Left1DY_+="TH2DataVar12_Left1DY"; 
  
  TH2SngVar12_Left1DY = TH2SngVar12_Left->ProjectionY(titleSngVar12_Left1DY_,0,TH2SngVar12_Left->GetYaxis()->GetNbins()+1,"e");
  TH2BkgVar12_Left1DY = TH2BkgVar12_Left->ProjectionY(titleBkgVar12_Left1DY_,0,TH2BkgVar12_Left->GetYaxis()->GetNbins()+1,"e");
  TH2DataVar12_Left1DY = TH2DataVar12_Left->ProjectionY(titleDataVar12_Left1DY_,0,TH2DataVar12_Left->GetYaxis()->GetNbins()+1,"e");

  GiveName(&titleSngVar12_Right1DY_); titleSngVar12_Right1DY_+="TH2SngVar12_Right1DY"; 
  GiveName(&titleBkgVar12_Right1DY_); titleBkgVar12_Right1DY_+="TH2BkgVar12_Right1DY"; 
  GiveName(&titleDataVar12_Right1DY_); titleDataVar12_Right1DY_+="TH2DataVar12_Right1DY"; 
    
  TH2SngVar12_Right1DY = TH2SngVar12_Right->ProjectionY(titleSngVar12_Right1DY_,0,TH2SngVar12_Right->GetYaxis()->GetNbins()+1,"e");

  TH2BkgVar12_Right1DY = TH2BkgVar12_Right->ProjectionY(titleBkgVar12_Right1DY_,0,TH2BkgVar12_Right->GetYaxis()->GetNbins()+1,"e");
  TH2DataVar12_Right1DY = TH2DataVar12_Right->ProjectionY(titleDataVar12_Right1DY_,0,TH2DataVar12_Right->GetYaxis()->GetNbins()+1,"e");
  
} // maybe obsolete

void PtEtaBin::SetError1DYVar2plots(){

  SetError(TH2SngVar12_Left1DY);
  SetError(TH2BkgVar12_Left1DY);
  SetError(TH2DataVar12_Left1DY);
  SetError(TH2SngVar12_Right1DY);
  SetError(TH2BkgVar12_Right1DY);
  SetError(TH2DataVar12_Right1DY);
 
}

void PtEtaBin::ReweighLeft(){//obsolete
  return;
  //does make the 1DY plots +reweight one of them
  /*  if(debug_>2) cout << "+----> PtEtaBin::ReweighLeft - start method" << endl; 

  GiveName(&titleData_Left1DYReweigh_); titleData_Left1DYReweigh_+="TH2Data_Left1DYReweigh"; 

  TH2Data_Left1DYReweigh = new TH1D(titleData_Left1DYReweigh_,titleData_Left1DYReweigh_,TH2Data_Left->GetYaxis()->GetNbins(),TH2Data_Left->GetYaxis()->GetBinLowEdge(1),TH2Data_Left->GetYaxis()->GetBinUpEdge(TH2Data_Left->GetYaxis()->GetNbins()));

  if(debug_>2) cout << "+----> PtEtaBin::ReweighLeft - finish declaring histogram names" << endl; 

  CalculateReweigh(TH2Data_Left, TH2Data_RightLeft1DX, TH2Data_Left1DYReweigh);
  */
} // certainly obsolete


void PtEtaBin::ReweighRight(){
  //the useFit here should be the same as the one of the makeXRatioPlot method, maybe make a global variable for this
  //does make the 1DY plots +reweight one of them
  if(debug_>2) cout << "+----> PtEtaBin::ReweighRight - start method" << endl; 

  GiveName(&titleSng_Right1DYReweigh_); titleSng_Right1DYReweigh_+="TH2Sng_Right1DYReweigh"; 
  GiveName(&titleBkg_Right1DYReweigh_); titleBkg_Right1DYReweigh_+="TH2Bkg_Right1DYReweigh"; 
  GiveName(&titleData_Right1DYReweigh_); titleData_Right1DYReweigh_+="TH2Data_Right1DYReweigh"; 

  if(!varBinSize_){
    TH2Sng_Right1DYReweigh = new TH1D(titleSng_Right1DYReweigh_,titleSng_Right1DYReweigh_,TH2Sng_Right->GetYaxis()->GetNbins(),TH2Sng_Right->GetYaxis()->GetBinLowEdge(1),TH2Sng_Right->GetYaxis()->GetBinUpEdge(TH2Sng_Right->GetYaxis()->GetNbins()));
    TH2Bkg_Right1DYReweigh = new TH1D(titleBkg_Right1DYReweigh_,titleBkg_Right1DYReweigh_,TH2Bkg_Right->GetYaxis()->GetNbins(),TH2Bkg_Right->GetYaxis()->GetBinLowEdge(1),TH2Bkg_Right->GetYaxis()->GetBinUpEdge(TH2Bkg_Right->GetYaxis()->GetNbins()));
    TH2Data_Right1DYReweigh = new TH1D(titleData_Right1DYReweigh_,titleData_Right1DYReweigh_,TH2Data_Right->GetYaxis()->GetNbins(),TH2Data_Right->GetYaxis()->GetBinLowEdge(1),TH2Data_Right->GetYaxis()->GetBinUpEdge(TH2Data_Right->GetYaxis()->GetNbins()));
  }
  if(varBinSize_){
    TH2Sng_Right1DYReweigh = new TH1D(titleSng_Right1DYReweigh_,titleSng_Right1DYReweigh_,TH2Sng_Right->GetYaxis()->GetNbins(),rangesbTag_);
    TH2Bkg_Right1DYReweigh = new TH1D(titleBkg_Right1DYReweigh_,titleBkg_Right1DYReweigh_,TH2Bkg_Right->GetYaxis()->GetNbins(),rangesbTag_);
    TH2Data_Right1DYReweigh = new TH1D(titleData_Right1DYReweigh_,titleData_Right1DYReweigh_,TH2Data_Right->GetYaxis()->GetNbins(),rangesbTag_);
  }    

  GiveName(&titleSng_Left1DYReweigh_); titleSng_Left1DYReweigh_+="TH2Sng_Left1DYReweigh"; 
  GiveName(&titleBkg_Left1DYReweigh_); titleBkg_Left1DYReweigh_+="TH2Bkg_Left1DYReweigh"; 
  GiveName(&titleData_Left1DYReweigh_); titleData_Left1DYReweigh_+="TH2Data_Left1DYReweigh"; 

  if(!varBinSize_){
    TH2Sng_Left1DYReweigh = new TH1D(titleSng_Left1DYReweigh_,titleSng_Left1DYReweigh_,TH2Sng_Left->GetYaxis()->GetNbins(),TH2Sng_Left->GetYaxis()->GetBinLowEdge(1),TH2Sng_Left->GetYaxis()->GetBinUpEdge(TH2Sng_Left->GetYaxis()->GetNbins()));
    TH2Bkg_Left1DYReweigh = new TH1D(titleBkg_Left1DYReweigh_,titleBkg_Left1DYReweigh_,TH2Bkg_Left->GetYaxis()->GetNbins(),TH2Bkg_Left->GetYaxis()->GetBinLowEdge(1),TH2Bkg_Left->GetYaxis()->GetBinUpEdge(TH2Bkg_Left->GetYaxis()->GetNbins()));
    TH2Data_Left1DYReweigh = new TH1D(titleData_Left1DYReweigh_,titleData_Left1DYReweigh_,TH2Data_Left->GetYaxis()->GetNbins(),TH2Data_Left->GetYaxis()->GetBinLowEdge(1),TH2Data_Left->GetYaxis()->GetBinUpEdge(TH2Data_Left->GetYaxis()->GetNbins()));
  }
  if(varBinSize_){
    TH2Sng_Left1DYReweigh = new TH1D(titleSng_Left1DYReweigh_,titleSng_Left1DYReweigh_,TH2Sng_Left->GetYaxis()->GetNbins(),rangesbTag_);
    TH2Bkg_Left1DYReweigh = new TH1D(titleBkg_Left1DYReweigh_,titleBkg_Left1DYReweigh_,TH2Bkg_Left->GetYaxis()->GetNbins(),rangesbTag_);
    TH2Data_Left1DYReweigh = new TH1D(titleData_Left1DYReweigh_,titleData_Left1DYReweigh_,TH2Data_Left->GetYaxis()->GetNbins(),rangesbTag_);
  }    

  if(debug_>2) cout << "+----> PtEtaBin::ReweighRight - finish declaring histogram names" << endl; 

  //if(!useFit){CalculateReweigh(TH2Sng_Right, TH2Data_LeftRight1DX, TH2Sng_Right1DYReweigh);}//obsolete
  //if(!useFit){CalculateReweigh(TH2Bkg_Right, TH2Data_LeftRight1DX, TH2Bkg_Right1DYReweigh);}
  //if(!useFit){CalculateReweigh(TH2Data_Right, TH2Data_LeftRight1DX, TH2Data_Right1DYReweigh);}
  //else { cout << "The useFit boolean is set to true, now the TH2Data_Right1DYReweigh histogram has been made but still needs to be fille by an extra pass over the data" << endl; }

  //start of Right1DYReweighRatio
  GiveName(&titleSng_Right1DYReweighRatio_); titleSng_Right1DYReweighRatio_+="TH2Sng_Right1DYReweighRatio"; 
  GiveName(&titleBkg_Right1DYReweighRatio_); titleBkg_Right1DYReweighRatio_+="TH2Bkg_Right1DYReweighRatio"; 
  GiveName(&titleData_Right1DYReweighRatio_); titleData_Right1DYReweighRatio_+="TH2Data_Right1DYReweighRatio"; 

  if(!varBinSize_){
    TH2Sng_Right1DYReweighRatio = new TH1D(titleSng_Right1DYReweighRatio_,titleSng_Right1DYReweighRatio_,TH2Sng_Right->GetYaxis()->GetNbins(),TH2Sng_Right->GetYaxis()->GetBinLowEdge(1),TH2Sng_Right->GetYaxis()->GetBinUpEdge(TH2Sng_Right->GetYaxis()->GetNbins()));
    TH2Bkg_Right1DYReweighRatio = new TH1D(titleBkg_Right1DYReweighRatio_,titleBkg_Right1DYReweighRatio_,TH2Bkg_Right->GetYaxis()->GetNbins(),TH2Bkg_Right->GetYaxis()->GetBinLowEdge(1),TH2Bkg_Right->GetYaxis()->GetBinUpEdge(TH2Bkg_Right->GetYaxis()->GetNbins()));
    TH2Data_Right1DYReweighRatio = new TH1D(titleData_Right1DYReweighRatio_,titleData_Right1DYReweighRatio_,TH2Data_Right->GetYaxis()->GetNbins(),TH2Data_Right->GetYaxis()->GetBinLowEdge(1),TH2Data_Right->GetYaxis()->GetBinUpEdge(TH2Data_Right->GetYaxis()->GetNbins()));
  }
  if(varBinSize_){
    TH2Sng_Right1DYReweighRatio = new TH1D(titleSng_Right1DYReweighRatio_,titleSng_Right1DYReweighRatio_,TH2Sng_Right->GetYaxis()->GetNbins(),rangesbTag_);
    TH2Bkg_Right1DYReweighRatio = new TH1D(titleBkg_Right1DYReweighRatio_,titleBkg_Right1DYReweighRatio_,TH2Bkg_Right->GetYaxis()->GetNbins(),rangesbTag_);
    TH2Data_Right1DYReweighRatio = new TH1D(titleData_Right1DYReweighRatio_,titleData_Right1DYReweighRatio_,TH2Data_Right->GetYaxis()->GetNbins(),rangesbTag_);
  }
  //end of Right1DYReweighRatio

  //start of Left1DYReweighRatio
  GiveName(&titleSng_Left1DYReweighRatio_); titleSng_Left1DYReweighRatio_+="TH2Sng_Left1DYReweighRatio"; 
  GiveName(&titleBkg_Left1DYReweighRatio_); titleBkg_Left1DYReweighRatio_+="TH2Bkg_Left1DYReweighRatio"; 
  GiveName(&titleData_Left1DYReweighRatio_); titleData_Left1DYReweighRatio_+="TH2Data_Left1DYReweighRatio"; 

  if(!varBinSize_){
    TH2Sng_Left1DYReweighRatio = new TH1D(titleSng_Left1DYReweighRatio_,titleSng_Left1DYReweighRatio_,TH2Sng_Left->GetYaxis()->GetNbins(),TH2Sng_Left->GetYaxis()->GetBinLowEdge(1),TH2Sng_Left->GetYaxis()->GetBinUpEdge(TH2Sng_Left->GetYaxis()->GetNbins()));
    TH2Bkg_Left1DYReweighRatio = new TH1D(titleBkg_Left1DYReweighRatio_,titleBkg_Left1DYReweighRatio_,TH2Bkg_Left->GetYaxis()->GetNbins(),TH2Bkg_Left->GetYaxis()->GetBinLowEdge(1),TH2Bkg_Left->GetYaxis()->GetBinUpEdge(TH2Bkg_Left->GetYaxis()->GetNbins()));
    TH2Data_Left1DYReweighRatio = new TH1D(titleData_Left1DYReweighRatio_,titleData_Left1DYReweighRatio_,TH2Data_Left->GetYaxis()->GetNbins(),TH2Data_Left->GetYaxis()->GetBinLowEdge(1),TH2Data_Left->GetYaxis()->GetBinUpEdge(TH2Data_Left->GetYaxis()->GetNbins()));
  }
  if(varBinSize_){
    TH2Sng_Left1DYReweighRatio = new TH1D(titleSng_Left1DYReweighRatio_,titleSng_Left1DYReweighRatio_,TH2Sng_Left->GetYaxis()->GetNbins(),rangesbTag_);
    TH2Bkg_Left1DYReweighRatio = new TH1D(titleBkg_Left1DYReweighRatio_,titleBkg_Left1DYReweighRatio_,TH2Bkg_Left->GetYaxis()->GetNbins(),rangesbTag_);
    TH2Data_Left1DYReweighRatio = new TH1D(titleData_Left1DYReweighRatio_,titleData_Left1DYReweighRatio_,TH2Data_Left->GetYaxis()->GetNbins(),rangesbTag_);
  }
  //end of Left1DYReweighRatio

  //start of Right1DY_LeftRightRatio
  GiveName(&titleSng_Right1DY_LeftRightRatio_); titleSng_Right1DY_LeftRightRatio_+="TH2Sng_Right1DY_LeftRightRatio"; 
  GiveName(&titleBkg_Right1DY_LeftRightRatio_); titleBkg_Right1DY_LeftRightRatio_+="TH2Bkg_Right1DY_LeftRightRatio"; 
  GiveName(&titleData_Right1DY_LeftRightRatio_); titleData_Right1DY_LeftRightRatio_+="TH2Data_Right1DY_LeftRightRatio"; 

  if(!varBinSize_){
    TH2Sng_Right1DY_LeftRightRatio = new TH1D(titleSng_Right1DY_LeftRightRatio_,titleSng_Right1DY_LeftRightRatio_,TH2Sng_Right->GetYaxis()->GetNbins(),TH2Sng_Right->GetYaxis()->GetBinLowEdge(1),TH2Sng_Right->GetYaxis()->GetBinUpEdge(TH2Sng_Right->GetYaxis()->GetNbins()));
    TH2Bkg_Right1DY_LeftRightRatio = new TH1D(titleBkg_Right1DY_LeftRightRatio_,titleBkg_Right1DY_LeftRightRatio_,TH2Bkg_Right->GetYaxis()->GetNbins(),TH2Bkg_Right->GetYaxis()->GetBinLowEdge(1),TH2Bkg_Right->GetYaxis()->GetBinUpEdge(TH2Bkg_Right->GetYaxis()->GetNbins()));
    TH2Data_Right1DY_LeftRightRatio = new TH1D(titleData_Right1DY_LeftRightRatio_,titleData_Right1DY_LeftRightRatio_,TH2Data_Right->GetYaxis()->GetNbins(),TH2Data_Right->GetYaxis()->GetBinLowEdge(1),TH2Data_Right->GetYaxis()->GetBinUpEdge(TH2Data_Right->GetYaxis()->GetNbins()));
  }
  if(varBinSize_){
    TH2Sng_Right1DY_LeftRightRatio = new TH1D(titleSng_Right1DY_LeftRightRatio_,titleSng_Right1DY_LeftRightRatio_,TH2Sng_Right->GetYaxis()->GetNbins(),rangesbTag_);
    TH2Bkg_Right1DY_LeftRightRatio = new TH1D(titleBkg_Right1DY_LeftRightRatio_,titleBkg_Right1DY_LeftRightRatio_,TH2Bkg_Right->GetYaxis()->GetNbins(),rangesbTag_);
    TH2Data_Right1DY_LeftRightRatio = new TH1D(titleData_Right1DY_LeftRightRatio_,titleData_Right1DY_LeftRightRatio_,TH2Data_Right->GetYaxis()->GetNbins(),rangesbTag_);
  }
  //end of Right1DY_LeftRightRatio

  //start of Right1DYReweigh_LeftRightRatio
  GiveName(&titleSng_Right1DYReweigh_LeftRightRatio_); titleSng_Right1DYReweigh_LeftRightRatio_+="TH2Sng_Right1DYReweigh_LeftRightRatio"; 
  GiveName(&titleBkg_Right1DYReweigh_LeftRightRatio_); titleBkg_Right1DYReweigh_LeftRightRatio_+="TH2Bkg_Right1DYReweigh_LeftRightRatio"; 
  GiveName(&titleData_Right1DYReweigh_LeftRightRatio_); titleData_Right1DYReweigh_LeftRightRatio_+="TH2Data_Right1DYReweigh_LeftRightRatio"; 

  if(!varBinSize_){
    TH2Sng_Right1DYReweigh_LeftRightRatio = new TH1D(titleSng_Right1DYReweigh_LeftRightRatio_,titleSng_Right1DYReweigh_LeftRightRatio_,TH2Sng_Right->GetYaxis()->GetNbins(),TH2Sng_Right->GetYaxis()->GetBinLowEdge(1),TH2Sng_Right->GetYaxis()->GetBinUpEdge(TH2Sng_Right->GetYaxis()->GetNbins()));
    TH2Bkg_Right1DYReweigh_LeftRightRatio = new TH1D(titleBkg_Right1DYReweigh_LeftRightRatio_,titleBkg_Right1DYReweigh_LeftRightRatio_,TH2Bkg_Right->GetYaxis()->GetNbins(),TH2Bkg_Right->GetYaxis()->GetBinLowEdge(1),TH2Bkg_Right->GetYaxis()->GetBinUpEdge(TH2Bkg_Right->GetYaxis()->GetNbins()));
    TH2Data_Right1DYReweigh_LeftRightRatio = new TH1D(titleData_Right1DYReweigh_LeftRightRatio_,titleData_Right1DYReweigh_LeftRightRatio_,TH2Data_Right->GetYaxis()->GetNbins(),TH2Data_Right->GetYaxis()->GetBinLowEdge(1),TH2Data_Right->GetYaxis()->GetBinUpEdge(TH2Data_Right->GetYaxis()->GetNbins()));
  }
  if(varBinSize_){
    TH2Sng_Right1DYReweigh_LeftRightRatio = new TH1D(titleSng_Right1DYReweigh_LeftRightRatio_,titleSng_Right1DYReweigh_LeftRightRatio_,TH2Sng_Right->GetYaxis()->GetNbins(),rangesbTag_);
    TH2Bkg_Right1DYReweigh_LeftRightRatio = new TH1D(titleBkg_Right1DYReweigh_LeftRightRatio_,titleBkg_Right1DYReweigh_LeftRightRatio_,TH2Bkg_Right->GetYaxis()->GetNbins(),rangesbTag_);
    TH2Data_Right1DYReweigh_LeftRightRatio = new TH1D(titleData_Right1DYReweigh_LeftRightRatio_,titleData_Right1DYReweigh_LeftRightRatio_,TH2Data_Right->GetYaxis()->GetNbins(),rangesbTag_);
  }
  //end of Right1DYReweigh_LeftRightRatio
 
  //start of Left1DYReweigh_LeftRightRatio
  GiveName(&titleSng_Left1DYReweigh_LeftRightRatio_); titleSng_Left1DYReweigh_LeftRightRatio_+="TH2Sng_Left1DYReweigh_LeftRightRatio"; 
  GiveName(&titleBkg_Left1DYReweigh_LeftRightRatio_); titleBkg_Left1DYReweigh_LeftRightRatio_+="TH2Bkg_Left1DYReweigh_LeftRightRatio"; 
  GiveName(&titleData_Left1DYReweigh_LeftRightRatio_); titleData_Left1DYReweigh_LeftRightRatio_+="TH2Data_Left1DYReweigh_LeftRightRatio"; 

  if(!varBinSize_){
    TH2Sng_Left1DYReweigh_LeftRightRatio = new TH1D(titleSng_Left1DYReweigh_LeftRightRatio_,titleSng_Left1DYReweigh_LeftRightRatio_,TH2Sng_Left->GetYaxis()->GetNbins(),TH2Sng_Left->GetYaxis()->GetBinLowEdge(1),TH2Sng_Left->GetYaxis()->GetBinUpEdge(TH2Sng_Left->GetYaxis()->GetNbins()));
    TH2Bkg_Left1DYReweigh_LeftRightRatio = new TH1D(titleBkg_Left1DYReweigh_LeftRightRatio_,titleBkg_Left1DYReweigh_LeftRightRatio_,TH2Bkg_Left->GetYaxis()->GetNbins(),TH2Bkg_Left->GetYaxis()->GetBinLowEdge(1),TH2Bkg_Left->GetYaxis()->GetBinUpEdge(TH2Bkg_Left->GetYaxis()->GetNbins()));
    TH2Data_Left1DYReweigh_LeftRightRatio = new TH1D(titleData_Left1DYReweigh_LeftRightRatio_,titleData_Left1DYReweigh_LeftRightRatio_,TH2Data_Left->GetYaxis()->GetNbins(),TH2Data_Left->GetYaxis()->GetBinLowEdge(1),TH2Data_Left->GetYaxis()->GetBinUpEdge(TH2Data_Left->GetYaxis()->GetNbins()));
  }
  if(varBinSize_){
    TH2Sng_Left1DYReweigh_LeftRightRatio = new TH1D(titleSng_Left1DYReweigh_LeftRightRatio_,titleSng_Left1DYReweigh_LeftRightRatio_,TH2Sng_Left->GetYaxis()->GetNbins(),rangesbTag_);
    TH2Bkg_Left1DYReweigh_LeftRightRatio = new TH1D(titleBkg_Left1DYReweigh_LeftRightRatio_,titleBkg_Left1DYReweigh_LeftRightRatio_,TH2Bkg_Left->GetYaxis()->GetNbins(),rangesbTag_);
    TH2Data_Left1DYReweigh_LeftRightRatio = new TH1D(titleData_Left1DYReweigh_LeftRightRatio_,titleData_Left1DYReweigh_LeftRightRatio_,TH2Data_Left->GetYaxis()->GetNbins(),rangesbTag_);
  }
  //end of Left1DYReweigh_LeftRightRatio

  //also reweigh the var2 plot
  GiveName(&titleSngVar12_Right1DYReweigh_); titleSngVar12_Right1DYReweigh_+="TH2SngVar12_Right1DYReweigh"; 
  GiveName(&titleBkgVar12_Right1DYReweigh_); titleBkgVar12_Right1DYReweigh_+="TH2BkgVar12_Right1DYReweigh"; 
  GiveName(&titleDataVar12_Right1DYReweigh_); titleDataVar12_Right1DYReweigh_+="TH2DataVar12_Right1DYReweigh"; 
 
  TH2SngVar12_Right1DYReweigh = new TH1D(titleSngVar12_Right1DYReweigh_,titleSngVar12_Right1DYReweigh_,TH2SngVar12_Right->GetYaxis()->GetNbins(),TH2SngVar12_Right->GetYaxis()->GetBinLowEdge(1),TH2SngVar12_Right->GetYaxis()->GetBinUpEdge(TH2SngVar12_Right->GetYaxis()->GetNbins()));
  TH2BkgVar12_Right1DYReweigh = new TH1D(titleBkgVar12_Right1DYReweigh_,titleBkgVar12_Right1DYReweigh_,TH2BkgVar12_Right->GetYaxis()->GetNbins(),TH2BkgVar12_Right->GetYaxis()->GetBinLowEdge(1),TH2BkgVar12_Right->GetYaxis()->GetBinUpEdge(TH2BkgVar12_Right->GetYaxis()->GetNbins()));
  TH2DataVar12_Right1DYReweigh = new TH1D(titleDataVar12_Right1DYReweigh_,titleDataVar12_Right1DYReweigh_,TH2DataVar12_Right->GetYaxis()->GetNbins(),TH2DataVar12_Right->GetYaxis()->GetBinLowEdge(1),TH2DataVar12_Right->GetYaxis()->GetBinUpEdge(TH2DataVar12_Right->GetYaxis()->GetNbins()));
 
  //if(!useFit){CalculateReweigh(TH2DataVar12_Right, TH2Data_LeftRight1DX, TH2DataVar12_Right1DYReweigh);}
  //if(!useFit){CalculateReweigh(TH2SngVar12_Right, TH2Data_LeftRight1DX, TH2SngVar12_Right1DYReweigh);}
  //if(!useFit){CalculateReweigh(TH2BkgVar12_Right, TH2Data_LeftRight1DX, TH2BkgVar12_Right1DYReweigh);}
  //else { cout << "The useFit boolean is set to true, now the  histogram has been made but still needs to be fille by an extra pass over the data" << endl; }

  GiveName(&titleSngVar12_Right1DYReweighRatio_); titleSngVar12_Right1DYReweighRatio_+="TH2SngVar12_Right1DYReweighRatio"; 
  GiveName(&titleBkgVar12_Right1DYReweighRatio_); titleBkgVar12_Right1DYReweighRatio_+="TH2BkgVar12_Right1DYReweighRatio"; 
  GiveName(&titleDataVar12_Right1DYReweighRatio_); titleDataVar12_Right1DYReweighRatio_+="TH2DataVar12_Right1DYReweighRatio"; 
 
  TH2SngVar12_Right1DYReweighRatio = new TH1D(titleSngVar12_Right1DYReweighRatio_,titleSngVar12_Right1DYReweighRatio_,TH2SngVar12_Right->GetYaxis()->GetNbins(),TH2SngVar12_Right->GetYaxis()->GetBinLowEdge(1),TH2SngVar12_Right->GetYaxis()->GetBinUpEdge(TH2SngVar12_Right->GetYaxis()->GetNbins()));
  TH2BkgVar12_Right1DYReweighRatio = new TH1D(titleBkgVar12_Right1DYReweighRatio_,titleBkgVar12_Right1DYReweighRatio_,TH2BkgVar12_Right->GetYaxis()->GetNbins(),TH2BkgVar12_Right->GetYaxis()->GetBinLowEdge(1),TH2BkgVar12_Right->GetYaxis()->GetBinUpEdge(TH2BkgVar12_Right->GetYaxis()->GetNbins()));
  TH2DataVar12_Right1DYReweighRatio = new TH1D(titleDataVar12_Right1DYReweighRatio_,titleDataVar12_Right1DYReweighRatio_,TH2DataVar12_Right->GetYaxis()->GetNbins(),TH2DataVar12_Right->GetYaxis()->GetBinLowEdge(1),TH2DataVar12_Right->GetYaxis()->GetBinUpEdge(TH2DataVar12_Right->GetYaxis()->GetNbins()));
 
  TH2Sng_Right1DYReweigh->Sumw2();
  TH2Bkg_Right1DYReweigh->Sumw2();
  TH2Data_Right1DYReweigh->Sumw2();
  
  TH2Sng_Left1DYReweigh->Sumw2();
  TH2Bkg_Left1DYReweigh->Sumw2();
  TH2Data_Left1DYReweigh->Sumw2();
  
  TH2Sng_Right1DYReweighRatio->Sumw2();
  TH2Bkg_Right1DYReweighRatio->Sumw2();
  TH2Data_Right1DYReweighRatio->Sumw2();
  
  TH2Sng_Left1DYReweighRatio->Sumw2(); 
  TH2Bkg_Left1DYReweighRatio->Sumw2(); 
  TH2Data_Left1DYReweighRatio->Sumw2();
 
  TH2Sng_Right1DY_LeftRightRatio->Sumw2(); 
  TH2Bkg_Right1DY_LeftRightRatio->Sumw2(); 
  TH2Data_Right1DY_LeftRightRatio->Sumw2(); 
  
  TH2Sng_Right1DYReweigh_LeftRightRatio->Sumw2(); 
  TH2Bkg_Right1DYReweigh_LeftRightRatio->Sumw2(); 
  TH2Data_Right1DYReweigh_LeftRightRatio->Sumw2(); 
   
  TH2Sng_Left1DYReweigh_LeftRightRatio->Sumw2();  
  TH2Bkg_Left1DYReweigh_LeftRightRatio->Sumw2(); 
  TH2Data_Left1DYReweigh_LeftRightRatio->Sumw2(); 
  
  TH2SngVar12_Right1DYReweigh->Sumw2();
  TH2BkgVar12_Right1DYReweigh->Sumw2(); 
  TH2DataVar12_Right1DYReweigh->Sumw2(); 
  
  TH2SngVar12_Right1DYReweighRatio->Sumw2(); 
  TH2BkgVar12_Right1DYReweighRatio->Sumw2(); 
  TH2DataVar12_Right1DYReweighRatio->Sumw2(); 
}

/*obsolete void PtEtaBin::ReweighRightChangeFitParams(double par0, double par1){

  TFData_LeftRight1DXFit->SetParameter(0,par0);
  TFData_LeftRight1DXFit->SetParameter(1,par1);

  }*/

void PtEtaBin::FillReweighControl(double btag, double* btagCuts, bool do2D, double weight, int partonFlavour, bool isW, bool isR, double controlVar1, double controlVar2, double controlVar0,  double lowCutVar0, double centralUpCutVar0, double centralLowCutVar0, double upCutVar0, double chisq){

  bool BinFound = false;
  
  double SignalControlWeight=0;
  SignalControlWeight=weight;
  double ptlo=0;
  double pthi=0;
  double etalo=0;
  double etahi=0;

  //cout << controlVar1 << "   " << controlVar2 << endl;

  if(do2D) {
    for(int i=0; i<=TH2Data_SignalControlVar12->GetXaxis()->GetNbins()+1; i++){
      for(int j=0; j<=TH2Data_SignalControlVar12->GetYaxis()->GetNbins()+1; j++){
	//cout << i << " " << TH2Data_SignalControlVar12->GetXaxis()->GetBinLowEdge(i) << " " << TH2Data_SignalControlVar12->GetXaxis()->GetBinUpEdge(i) << endl;
	//cout << j << " " << TH2Data_SignalControlVar12->GetYaxis()->GetBinLowEdge(j) << " " << TH2Data_SignalControlVar12->GetYaxis()->GetBinUpEdge(j) << endl;
	
	ptlo=TH2Data_SignalControlVar12->GetXaxis()->GetBinLowEdge(i);
	pthi=TH2Data_SignalControlVar12->GetXaxis()->GetBinUpEdge(i);
	etalo=TH2Data_SignalControlVar12->GetYaxis()->GetBinLowEdge(j);
	etahi=TH2Data_SignalControlVar12->GetYaxis()->GetBinUpEdge(j);
	if(i==TH2Data_SignalControlVar12->GetXaxis()->GetNbins()+1) pthi=999999.;
	if(j==TH2Data_SignalControlVar12->GetYaxis()->GetNbins()+1) etahi=999999.;
	
	if(controlVar1>=ptlo && controlVar1<pthi){
	  if(controlVar2>=etalo && controlVar2<etahi){
	    
	    BinFound=true;
	    //cout << " var1 " << controlVar1 << " var 2 " << controlVar2 << "   TH2Data_SignalControlVar12->GetBinContent(i,j) " << TH2Data_SignalControlVar12->GetBinContent(i,j) << endl;
	    if(TH2Data_SignalControlVar12->GetBinContent(i,j)>0){
			SignalControlWeight = SignalControlWeight*TH2Data_SignalControlVar12->GetBinContent(i,j); 
	      //cout << "SignalControlWeight :" << SignalControlWeight << endl;
	      //if(controlVar0<100) cout << "SignalControlWeight: " << TH2Data_SignalControlVar12->GetBinContent(i,j) << endl;

	    }
	  }
	}
      }
    }
    if(!BinFound && controlVar1!=-999 && controlVar2!=999) cout << "FillReweighControl::No Bin Found containing " << controlVar1 << " " << controlVar2 << endl;
    //if(!BinFound) cout << "FillReweighControl::No Bin Found containing " << controlVar1 << " " << controlVar2 << endl;
  } else {
  
    for(int i=0; i<=TH2Data_SignalControlVar12_1DX->GetXaxis()->GetNbins()+1; i++){
      ptlo=TH2Data_SignalControlVar12_1DX->GetXaxis()->GetBinLowEdge(i);
      pthi=TH2Data_SignalControlVar12_1DX->GetXaxis()->GetBinUpEdge(i);
      if(i==TH2Data_SignalControlVar12_1DX->GetXaxis()->GetNbins()+1) pthi=999999.;
      
      if(controlVar1>=ptlo && controlVar1<pthi){
	BinFound=true;
	if(TH2Data_SignalControlVar12_1DX->GetBinContent(i)>0){
	  SignalControlWeight = SignalControlWeight*TH2Data_SignalControlVar12_1DX->GetBinContent(i);
	}
      }
    }
    //if(!BinFound) cout << "FillReweighControl::No Bin Found containing " << controlVar1 << endl;
    

    //to do the 1D reweigh only in 2nd var eta (also tried deltaR)
    /*for(int i=0; i<=TH2Data_SignalControlVar12_1DY->GetXaxis()->GetNbins()+1; i++){
      ptlo=TH2Data_SignalControlVar12_1DY->GetXaxis()->GetBinLowEdge(i);
      pthi=TH2Data_SignalControlVar12_1DY->GetXaxis()->GetBinUpEdge(i);
      if(i==TH2Data_SignalControlVar12_1DY->GetXaxis()->GetNbins()+1) pthi=999999.;
      
      if(controlVar2>=ptlo && controlVar2<pthi){
	BinFound=true;
	if(TH2Data_SignalControlVar12_1DY->GetBinContent(i)>0){
	  SignalControlWeight = SignalControlWeight*TH2Data_SignalControlVar12_1DY->GetBinContent(i);
	}
      }
    }
    if(!BinFound) cout << "FillReweighControl::No Bin Found containing " << controlVar2 << endl;*/
    
  }
	
	
	//if (SignalControlWeight > 40) SignalControlWeight=0;
  
   
  if(controlVar1>ptbinlow_ && controlVar1<ptbinup_){
    if(controlVar2>etabinlow_ && controlVar2<etabinup_){  //if(fabs(partonFlavour)!=5) if(chisq<3) TH1Data_ControlVarReweigh->Fill(controlVar0,SignalControlWeight); //abuse for chisq-reweigh
      //if(fabs(partonFlavour)!=5) if(chisq<10) TH1Data_ControlVarReweigh->Fill(controlVar0,1); //abuse for chisq-reweigh
      
      TH1Data_ControlVarReweigh->Fill(controlVar0,SignalControlWeight); 
        
        //cout << controlVar0 << " -> " << SignalControlWeight << endl << endl << endl; exit(1);
	
		histo1D["TH1Data_Var0_VVData"]->Fill(controlVar0,SignalControlWeight); 
				
		
      //TH1Data_ControlVarReweigh->Fill(controlVar0,1); 
      if(fabs(partonFlavour)==5) TH1Sng_ControlVarReweigh->Fill(controlVar0,SignalControlWeight);
      if(fabs(partonFlavour)!=5) TH1Bkg_ControlVarReweigh->Fill(controlVar0,SignalControlWeight);
      if(fabs(partonFlavour)!=5) if(isW) TH1Bkg_W_ControlVarReweigh->Fill(controlVar0,SignalControlWeight);
      if(fabs(partonFlavour)!=5) if(isR) TH1Bkg_R_ControlVarReweigh->Fill(controlVar0,SignalControlWeight);
    }
  }
}

//void PtEtaBin::FillReweighRight(bool do2D, bool useFit, double weight, int partonFlavour, double bTag, double var1, double var2, double var0, int partonFlavourControl, double bTagControl, double controlVar1, double controlVar2, double controlVar0,  double lowCutVar0, double centralCutVar0, double upCutVar0){
void PtEtaBin::FillReweighRight(bool do2D, bool useFit, double weight, int partonFlavour, double bTag, double var1, double var2, double var0, double lowCutVar0, double centralLowCutVar0, double centralUpCutVar0, double upCutVar0, double ptbinlo, double ptbinhi, double etabinlo, double etabinhi){

  //if(debug_>2) cout << "+----> PtEtaBin::FillReweighRight - start method" << endl; 
  
  //if(debug_>2) cout << "weight " << weight << " partonFlavour " << partonFlavour << " bTag " << bTag << " var1 " << var1 << " var0 " << var0 << " lowCutVar0 " << lowCutVar0 << " centralCutVar0 " << centralCutVar0 << " upCutVar0 " << upCutVar0 << endl; 
  //warning, the lowCutVar0 are inputs of this function, they should be IDENTICAL to the ones of the fill Signal samples method!! Maybe make a global variable for this
  
  //doShift stands fro filling the left and right reweigh plot with the shifted value
  
  //cout << "weight- before " << weight << endl;
  if(!doShift_){
    
    bool BinFound = false;
    
    //************************************************************** 
    //This part is for the reweighing of the btag distribution right
    
    //if(controlVar0>lowCutVar0 && controlVar0<centralCutVar0){ //to reweigh left region, not used
    if(var0>=centralUpCutVar0 && var0<upCutVar0){
      
      //if(debug_>2) cout << "TFData_LeftRight1DXFit->Eval(var1) " << TFData_LeftRight1DXFit->Eval(var1) << endl;
      //following line is to do the reweighing using the fit, if commented and NTupleAnalyzer is used with fit I can mis-use the code for the shifting
      if(!do2D){

	if(ptbinhi>0 && ptbinlo>0){
	  //weight = weight * ((ptbinhi-ptbinlo)/TFData_LeftRight1DXFit->Integral(ptbinlo,ptbinhi));
	  //cout << ptbinlo << " " << ptbinhi << " integral " << TFData_LeftRight1DXFit->Integral(ptbinlo,ptbinhi) << endl;
	}

	//cout << "fillrew  " << weight << endl;

	//cout << TFData_LeftRight1DXFit->GetParameter(0) << " " << TFData_LeftRight1DXFit->GetParameter(1) << endl;
		  //cout << "Old weight " << weight << endl;

	if(useFit) weight = weight * TFData_LeftRight1DXFit->Eval(var1);//was used for the reweighting with the fitted pt ratio histo
		  
		  //if (useFit)cout << "New weight " << weight << endl;
		  //cout << "POM WEIGHT " << weight << endl;
		  
		  //exit(0);
	//cout << "fillrew  " << TFData_LeftRight1DXFit->Eval(var1) << " " << weight << endl;
	

	if(!useFit){
	  //cout << "NOT USING THE FIT --> please debug the (var1>300)" << endl;
	  for(int err=0; err<100; err++){
	    //cout << "err: NOT USING THE FIT" << endl;
	  }
	  // 1D way (only var1 = pt)
	  for(int i=0; i<=TH2Data_LeftRight1DX->GetXaxis()->GetNbins()+1; i++){//was used to reweigh with the discrete pt ratio histogram
	    if(var1>=TH2Data_LeftRight1DX->GetXaxis()->GetBinLowEdge(i) && var1<TH2Data_LeftRight1DX->GetXaxis()->GetBinUpEdge(i)){
	      
	      BinFound=true;
	      //cout << " var1 " << var1 << " TH2Data_LeftRight1DX->GetBinContent(i) " << TH2Data_LeftRight1DX->GetBinContent(i) << endl;
	      if(TH2Data_LeftRight1DX->GetBinContent(i)>0){
		if(var1>300){
		  //weight = weight*0.3; 
		} else if(var1<50) {
		  //weight = weight*TH2Data_LeftRight1DX->GetBinContent(i)*0.95;		
		} else {
		  //weight = weight*TH2Data_LeftRight1DX->GetBinContent(i);		
		}
	      }
	      
	    }
	  }
	}	
      }
      if(do2D){
	// 2D way (var1 and var2, pt and eta)  
	for(int i=0; i<=TH2DataVar12_LeftRight->GetXaxis()->GetNbins()+1; i++){//was used to reweigh with the discrete pt ratio histogram
	  for(int j=0; j<=TH2DataVar12_LeftRight->GetYaxis()->GetNbins()+1; j++){//was used to reweigh with the discrete pt ratio histogram
	    if(var1>=TH2DataVar12_LeftRight->GetXaxis()->GetBinLowEdge(i) && var1<TH2DataVar12_LeftRight->GetXaxis()->GetBinUpEdge(i)){
	      if(var2>=TH2DataVar12_LeftRight->GetYaxis()->GetBinLowEdge(j) && var2<TH2DataVar12_LeftRight->GetYaxis()->GetBinUpEdge(j)){
		BinFound=true;
		//cout << " var1 " << var1 << " TH2Data_LeftRight1DX->GetBinContent(i) " << TH2Data_LeftRight1DX->GetBinContent(i) << endl;
		if(TH2DataVar12_LeftRight->GetBinContent(i,j)>0){
		  weight = weight*TH2DataVar12_LeftRight->GetBinContent(i,j);
		}
	      }
	    }
	  }
	}
	if(!BinFound) cout << "FillReweighRight::No Bin Found containing " << var1 << " " << var2 << endl;
      }
      
      if(var1>ptbinlow_ && var1<ptbinup_){
	if(var2>etabinlow_ && var2<etabinup_){ 
	  TH2Data_Right1DYReweigh->Fill(bTag,weight);
	  if(fabs(partonFlavour)==5) TH2Sng_Right1DYReweigh->Fill(bTag,weight);
	  if(fabs(partonFlavour)!=5) TH2Bkg_Right1DYReweigh->Fill(bTag,weight);
	} else {
	  cout << etabinlow_ << " " << etabinup_ << "  - " << var2 << endl;
	}
      }
      
      //cout << "rightreweigh weight " << weight << " -- " << bTag << endl;

    }
  }
  //cout << "weight- after " << weight <<  endl;
  
  double shift=0.;
  if(doShift_){
    //************************************************************** 
    //This part is for the left and right shifting
    //double shift=-0.2;//fix shift
    //double a=0.001538; // full range fit
    //double a=0.001942; // full range fit
    double a=0.001823; // full range fit
    //double a=0.001676; // 50 to 300 fit
    double b=-130*a;//130 is the reference mlb where everything is shifted to
    shift=-(a*var0+b);
    //cout << "shift at " << var0 << " is " << shift << endl;
    if(var0>=centralUpCutVar0 && var0<upCutVar0){
      TH2Data_Right1DYReweigh->Fill(bTag+shift,weight);
      if(fabs(partonFlavour)==5) TH2Sng_Right1DYReweigh->Fill(bTag+shift,weight);
      if(fabs(partonFlavour)!=5) TH2Bkg_Right1DYReweigh->Fill(bTag+shift,weight);
    }
  }
  
 
  TH2DataVar12_Right1DYReweigh->Fill(var2,weight);
  if(fabs(partonFlavour)==5) TH2SngVar12_Right1DYReweigh->Fill(var2,weight);
  if(fabs(partonFlavour)!=5) TH2BkgVar12_Right1DYReweigh->Fill(var2,weight);
  
  // }
  if(doShift_){
    if(var0>lowCutVar0 && var0<centralLowCutVar0){
      
      TH2Data_Left1DYReweigh->Fill(bTag+shift,weight);
      if(fabs(partonFlavour)==5) TH2Sng_Left1DYReweigh->Fill(bTag+shift,weight);
      if(fabs(partonFlavour)!=5) TH2Bkg_Left1DYReweigh->Fill(bTag+shift,weight);
      
    }
  }
  
  //if(debug_>2) cout << "+----> PtEtaBin::FillReweighRight - finish method" << endl; 
  
}

void PtEtaBin::MakeReweighRatio(){
  
  /* TH2Sng_Right1DY_LeftRightRatio->Divide(TH2Sng_Left1DY,TH2Sng_Right1DY);
  TH2Bkg_Right1DY_LeftRightRatio->Divide(TH2Bkg_Left1DY,TH2Bkg_Right1DY);
  TH2Data_Right1DY_LeftRightRatio->Divide(TH2Data_Left1DY,TH2Data_Right1DY);
  
  TH2Sng_Right1DYReweigh_LeftRightRatio->Divide(TH2Sng_Left1DY,TH2Sng_Right1DYReweigh);
  TH2Bkg_Right1DYReweigh_LeftRightRatio->Divide(TH2Bkg_Left1DY,TH2Bkg_Right1DYReweigh);
  TH2Data_Right1DYReweigh_LeftRightRatio->Divide(TH2Data_Left1DY,TH2Data_Right1DYReweigh);*/

  TH2Sng_Right1DY_LeftRightRatio->Add(TH2Sng_Left1DY,TH2Sng_Right1DY,1/TH2Sng_Left1DY->Integral(0,TH2Sng_Left1DY->GetXaxis()->GetNbins()+1),-1/TH2Sng_Right1DY->Integral(0,TH2Sng_Right1DY->GetXaxis()->GetNbins()+1));
  TH2Bkg_Right1DY_LeftRightRatio->Add(TH2Bkg_Left1DY,TH2Bkg_Right1DY,1/TH2Bkg_Left1DY->Integral(0,TH2Bkg_Left1DY->GetXaxis()->GetNbins()+1),-1/TH2Bkg_Right1DY->Integral(0,TH2Bkg_Right1DY->GetXaxis()->GetNbins()+1));
  TH2Data_Right1DY_LeftRightRatio->Add(TH2Data_Left1DY,TH2Data_Right1DY,1/TH2Data_Left1DY->Integral(0,TH2Data_Left1DY->GetXaxis()->GetNbins()+1),-1/TH2Data_Right1DY->Integral(0,TH2Data_Right1DY->GetXaxis()->GetNbins()+1));
  
  TH2Sng_Right1DYReweigh_LeftRightRatio->Add(TH2Sng_Left1DY,TH2Sng_Right1DYReweigh,1/TH2Sng_Left1DY->Integral(0,TH2Sng_Left1DY->GetXaxis()->GetNbins()+1),-1/TH2Sng_Right1DYReweigh->Integral(0,TH2Sng_Right1DYReweigh->GetXaxis()->GetNbins()+1));
  TH2Bkg_Right1DYReweigh_LeftRightRatio->Add(TH2Bkg_Left1DY,TH2Bkg_Right1DYReweigh,1/TH2Bkg_Left1DY->Integral(0,TH2Bkg_Left1DY->GetXaxis()->GetNbins()+1),-1/TH2Bkg_Right1DYReweigh->Integral(0,TH2Bkg_Right1DYReweigh->GetXaxis()->GetNbins()+1));
  TH2Data_Right1DYReweigh_LeftRightRatio->Add(TH2Data_Left1DY,TH2Data_Right1DYReweigh,1/TH2Data_Left1DY->Integral(0,TH2Data_Left1DY->GetXaxis()->GetNbins()+1),-1/TH2Data_Right1DYReweigh->Integral(0,TH2Data_Right1DYReweigh->GetXaxis()->GetNbins()+1));

  
  TH2Sng_Left1DYReweigh_LeftRightRatio->Add(TH2Sng_Left1DYReweigh,TH2Sng_Right1DYReweigh,1/TH2Sng_Left1DYReweigh->Integral(0,TH2Sng_Left1DYReweigh->GetXaxis()->GetNbins()+1),-1/TH2Sng_Right1DYReweigh->Integral(0,TH2Sng_Right1DYReweigh->GetXaxis()->GetNbins()+1));
  TH2Bkg_Left1DYReweigh_LeftRightRatio->Add(TH2Bkg_Left1DYReweigh,TH2Bkg_Right1DYReweigh,1/TH2Bkg_Left1DYReweigh->Integral(0,TH2Bkg_Left1DYReweigh->GetXaxis()->GetNbins()+1),-1/TH2Bkg_Right1DYReweigh->Integral(0,TH2Bkg_Right1DYReweigh->GetXaxis()->GetNbins()+1));
  TH2Data_Left1DYReweigh_LeftRightRatio->Add(TH2Data_Left1DYReweigh,TH2Data_Right1DYReweigh,1/TH2Data_Left1DYReweigh->Integral(0,TH2Data_Left1DY->GetXaxis()->GetNbins()+1),-1/TH2Data_Right1DYReweigh->Integral(0,TH2Data_Right1DYReweigh->GetXaxis()->GetNbins()+1));


  TH2Sng_Right1DYReweighRatio->Divide(TH2Sng_Right1DYReweigh,TH2Sng_Right1DY);
  TH2Bkg_Right1DYReweighRatio->Divide(TH2Bkg_Right1DYReweigh,TH2Bkg_Right1DY);
  TH2Data_Right1DYReweighRatio->Divide(TH2Data_Right1DYReweigh,TH2Data_Right1DY);

  TH2Sng_Left1DYReweighRatio->Divide(TH2Sng_Left1DYReweigh,TH2Sng_Left1DY);
  TH2Bkg_Left1DYReweighRatio->Divide(TH2Bkg_Left1DYReweigh,TH2Bkg_Left1DY);
  TH2Data_Left1DYReweighRatio->Divide(TH2Data_Left1DYReweigh,TH2Data_Left1DY);

  TH2SngVar12_Right1DYReweighRatio->Divide(TH2SngVar12_Right1DYReweigh,TH2SngVar12_Right1DY);
  TH2BkgVar12_Right1DYReweighRatio->Divide(TH2BkgVar12_Right1DYReweigh,TH2BkgVar12_Right1DY);
  TH2DataVar12_Right1DYReweighRatio->Divide(TH2DataVar12_Right1DYReweigh,TH2DataVar12_Right1DY);

}

void PtEtaBin::CalculateReweigh(TH2D* THData_,TH1D *THratio,TH1D *THReweigh){//obsolete
  return;

  double adder=0;
  for(int j=0; j<=THData_->GetYaxis()->GetNbins()+1; j++){
    for(int i=0; i<=THData_->GetXaxis()->GetNbins()+1; i++){
      //if(debug_>2) cout << "+----> PtEtaBin::ReweighLeft - TH2Data_RightLeft1DX->GetBinContent(" << i << ") = " << TH2Data_RightLeft1DX->GetBinContent(i) << endl;     
      //if(debug_>2) cout << "+----> PtEtaBin::ReweighLeft - TH2Data_->GetBinContent(" << i << "," << j << ") = " << TH2Data_->GetBinContent(i,j) << endl; 
      if(THratio->GetBinContent(i)>0){ //safety measure to prevent reading out -999 bins
	adder+=(THData_->GetBinContent(i,j)*THratio->GetBinContent(i));
      } else {
	if (THData_->GetBinContent(i,j)>0) cout << "!! ERROR: PtEtaBin::ReweighLeft  TH2Data_->GetBinContent(i,j)>0 but not taken into account - bin : " << i << " " << j << " " << THData_->GetBinContent(i,j) << endl;
      }
    }
    THReweigh->SetBinContent(j,adder);
    adder=0;
  }

}

void PtEtaBin::GetPercentiles(double left, double center, double right, int *returnval){

  /*cout << "percentile boundaries: " << "0 - " << lowerBin_ << " - "  << centralLowBin_ << " - " << upperBin_ << " - inf" << endl;

  double first=TH1Bkg_Var0->Integral(0,lowerBin_-1)/TH1Bkg_Var0->Integral(0, TH1Bkg_Var0->GetXaxis()->GetNbins()+1);
  double second=TH1Bkg_Var0->Integral(lowerBin_,centralLowBin_)/TH1Bkg_Var0->Integral(0, TH1Bkg_Var0->GetXaxis()->GetNbins()+1);
  double third=TH1Bkg_Var0->Integral(centralUpBin_+1,upperBin_)/TH1Bkg_Var0->Integral(0, TH1Bkg_Var0->GetXaxis()->GetNbins()+1);
  double fourth=TH1Bkg_Var0->Integral(upperBin_+1, TH1Bkg_Var0->GetXaxis()->GetNbins()+1)/TH1Bkg_Var0->Integral(0, TH1Bkg_Var0->GetXaxis()->GetNbins()+1);

  //cout << "percentile fractions:  " << first << " - " << second << " - " << third << " - " << fourth << endl;
  cout << "percentile fractions of current choice of LR:  " << first << " - " << first+second << " - " << first+second+third << endl;
  cout << "percentile fractions wanted LR (input)      :  " << left << " - " << center << " - " << right << endl;*/
  
  /*returnval[0]=first;
  returnval[1]=first+second;
  returnval[2]=first+second+third;*/ // was commented
  
  /*double binSum=0;
  bool low=true;
  bool up=true;
  bool middle=true;
  for(int bin=0; bin<TH1Bkg_Var0->GetXaxis()->GetNbins()+1; bin++){
    binSum+=TH1Bkg_Var0->GetBinContent(bin)/TH1Bkg_Var0->Integral(0, TH1Bkg_Var0->GetXaxis()->GetNbins()+1);
    if(binSum>left+0.0001 && low) {
      cout << "boundary-bin : " << bin << endl;
      cout << "boundary     : " << TH1Bkg_Var0->GetXaxis()->GetBinLowEdge(bin) << endl; 
      low=false;
      returnval[0]=(int) TH1Bkg_Var0->GetXaxis()->GetBinLowEdge(bin);
    }
    if(binSum>center+0.0001 && middle) {
      cout << "boundary-bin : " << bin << endl; 
      cout << "boundary     : " << TH1Bkg_Var0->GetXaxis()->GetBinLowEdge(bin) << endl; 
      middle=false;
      returnval[1]=(int) TH1Bkg_Var0->GetXaxis()->GetBinLowEdge(bin); 
    }
    if(binSum>right+0.0001 && up) {
      cout << "boundary-bin : " << bin << endl; 
      cout << "boundary     : " << TH1Bkg_Var0->GetXaxis()->GetBinUpEdge(bin) << endl; 
      up=false;
      returnval[2]=(int) TH1Bkg_Var0->GetXaxis()->GetBinLowEdge(bin); 
    }
  }*/
  
 

}; // very obsolete

void PtEtaBin::GetLRratio(bool dofitprint, bool doPrint, double FMCBias, bool doNewF, double newF){

  double leftAdder=0;
  double rightAdder=0;
  
  for(int i=lowerBin_; i<=centralLowBin_; i++){
    leftAdder+=TH1Data_ControlVar->GetBinContent(i);
  }
  for(int i=centralUpBin_+1; i<=upperBin_; i++){   
    rightAdder+=TH1Data_ControlVar->GetBinContent(i);
  }

  leftNonBe_      = 0;
  rightNonBe_     = 0;
  leftNonBMCe_      = 0;
  rightNonBMCe_     = 0;
  leftNonBReweighe_      = 0;
  rightNonBReweighe_     = 0;

  leftNonB_      = TH1Data_ControlVar->Integral(lowerBin_,centralLowBin_);
  rightNonB_     = TH1Data_ControlVar->Integral(centralUpBin_+1,upperBin_);
  leftNonBMC_      = TH1Bkg_Var0->Integral(lowerBin_,centralLowBin_);
 
  //rightNonBMC_     = TH1Bkg_Var0->Integral(centralUpBin_+1,upperBin_);
	rightNonBMC_     = TH2Bkg_Right1DYReweigh->Integral(0,TH2Bkg_Right1DYReweigh->GetNbinsX()+1);
	//rightNonBMC_     = TH2Bkg_Right1DY->Integral(0,TH2Bkg_Right1DY->GetNbinsX()+1);
	
  leftNonBReweigh_      = TH1Data_ControlVarReweigh->Integral(lowerBin_,centralLowBin_);
  rightNonBReweigh_     = TH1Data_ControlVarReweigh->Integral(centralUpBin_+1,upperBin_);
	
  //cout << "~~~~~~~~~~~~~~~~~~~~~~~ ControlVarReweigh BinContent " << TH1Data_ControlVarReweigh->GetBinContent(20) << "  -->Fix this difference!!" << endl;
  //cout << "~~~~~~~~~~~~~~~~~~~~~~~ ControlVarReweigh BinContent " << TH1Data_ControlVarReweigh->GetBinContent(60) << "  -->Fix this difference!!" << endl;
  //cout << "~~~~~~~~~~~~~~~~~~~~~~~ ControlVarReweigh BinContent " << TH1Data_ControlVarReweigh->GetBinContent(90) << "  -->Fix this difference!!" << endl;

  /*calculation of the error if we want to have it dependent on the weights
    for(int i=lowerBin_; i<=centralBin_; i++){
    leftNonBe_+=TH1Data_ControlVar->GetBinError(i)*TH1Data_ControlVar->GetBinError(i);
    leftNonBMCe_+=TH1Bkg_Var0->GetBinError(i)*TH1Bkg_Var0->GetBinError(i);
    leftNonBReweighe_+=TH1Data_ControlVarReweigh->GetBinError(i)*TH1Data_ControlVarReweigh->GetBinError(i);
  }
  leftNonBe_=sqrt(leftNonBe_);
  leftNonBMCe_=sqrt(leftNonBMCe_);
  leftNonBReweighe_=sqrt(leftNonBReweighe_);

  for(int i=centralBin_+1; i<=upperBin_; i++){
    rightNonBe_+=TH1Data_ControlVar->GetBinError(i)*TH1Data_ControlVar->GetBinError(i);
    rightNonBMCe_+=TH1Bkg_Var0->GetBinError(i)*TH1Bkg_Var0->GetBinError(i);
    rightNonBReweighe_+=TH1Data_ControlVarReweigh->GetBinError(i)*TH1Data_ControlVarReweigh->GetBinError(i);
  }
  rightNonBe_=sqrt(rightNonBe_);
  rightNonBMCe_=sqrt(rightNonBMCe_);
  rightNonBReweighe_=sqrt(rightNonBReweighe_);*/

  leftNonBe_             = sqrt(leftNonB_);
  rightNonBe_            = sqrt(rightNonB_);
  leftNonBMCe_           = sqrt(leftNonBMC_);
  rightNonBMCe_          = sqrt(rightNonBMC_);
  leftNonBReweighe_      = sqrt(leftNonBReweigh_);
  rightNonBReweighe_     = sqrt(rightNonBReweigh_);

  F_        = leftNonB_/rightNonB_;
  FMC_      = leftNonBMC_/rightNonBMC_;
  FReweigh_ = leftNonBReweigh_/rightNonBReweigh_;
	
	//FReweigh_ = 1.08;

  //cout << "########  " << FMC_ << " " << leftNonBMC_ <<  " " << rightNonBMC_ <<  " " << leftNonBMC_/rightNonBMC_ << endl;
  //cout << "######## FMC=" << FMC_ << " FReweigh="  << FReweigh_ << "   rel. bias=" << (FReweigh_-FMC_)/FMC_ << endl;

  //FMC_=FMC_+FMCBias;
  FMC_=FMC_*FMCBias;

  if(doNewF){
    F_      = newF;
    FReweigh_ = newF;

  }

  Fe_        = sqrt(pow(leftNonBe_       /rightNonB_       ,2)+pow(leftNonB_       *rightNonBe_       /(rightNonB_       *rightNonB_)       ,2));
  FMCe_      = sqrt(pow(leftNonBMCe_     /rightNonBMC_     ,2)+pow(leftNonBMC_     *rightNonBMCe_     /(rightNonBMC_     *rightNonBMC_)     ,2));
  FReweighe_ = sqrt(pow(leftNonBReweighe_/rightNonBReweigh_,2)+pow(leftNonBReweigh_*rightNonBReweighe_/(rightNonBReweigh_*rightNonBReweigh_),2));

  if(doNewF){
    Fe_      = 0;
    FReweighe_ = 0;

  }


  if(debug_>1) cout << "+--> ====================================================================+ "<< endl;
  
  if(doPrint){
    if(leftAdder!=leftNonB_ || rightAdder!=rightNonB_)  if(debug_>1) cout << "+--> ERROR!! PtEtaBin::GetLRratio - From Control sample leftAdder " << leftAdder << " rightAdder "  << rightAdder << "     - ratio:  " << leftAdder/rightAdder << endl;
    
    if(debug_>1) cout << "+--> --------------------------------------------------------------------+ "<< endl;
    if(debug_>1) cout << "+--> From Control sample " << endl;
    if(debug_>1) cout << "+--> leftb     " << TH1Sng_ControlVar->Integral(lowerBin_,centralLowBin_) << " rightb    "  << TH1Sng_ControlVar->Integral(centralUpBin_+1,upperBin_) << endl;
    if(debug_>1) cout << "+--> leftNonb  " << TH1Bkg_ControlVar->Integral(lowerBin_,centralLowBin_) << " rightNonb "  << TH1Bkg_ControlVar->Integral(centralUpBin_+1,upperBin_) << endl;
    
    if(debug_>1) cout << "+--> From ControlReweigh sample  " << endl;
    if(debug_>1) cout << "+--> leftb     " << TH1Sng_ControlVarReweigh->Integral(lowerBin_,centralLowBin_) << " rightb    "  << TH1Sng_ControlVarReweigh->Integral(centralUpBin_+1,upperBin_) << endl;
    if(debug_>1) cout << "+--> leftNonb  " << TH1Bkg_ControlVarReweigh->Integral(lowerBin_,centralLowBin_) << " rightNonb "  << TH1Bkg_ControlVarReweigh->Integral(centralUpBin_+1,upperBin_) << endl;
    
    if(debug_>1) cout << "+--> From Control sample          left b + Non b " << leftNonB_ << "+-" << leftNonBe_               << " right b + Non b  "  << rightNonB_ << "+-" << rightNonBe_ << "     - ratio:  " << leftNonB_/rightNonB_ <<  "+-" << Fe_ << endl;  
    
    if(debug_>1) cout << "+--> From ControlReweigh sample   left b + Non b " << leftNonBReweigh_ << "+-" << leftNonBReweighe_ << " right b + Non b  "  << rightNonBReweigh_ << "+-" << rightNonBReweighe_ << "     - ratio:  " << leftNonBReweigh_/rightNonBReweigh_ <<  "+-" << FReweighe_ << endl;
    
    
    if(debug_>1) cout << "+--> From Signal sample           leftNonb        " << leftNonBMC_ << "+-" << leftNonBMCe_          << " rightNonb        "  << rightNonBMC_ << "+-" << rightNonBMCe_ << "     - ratio:  " << leftNonBMC_/rightNonBMC_  << "+-" << FMCe_ << endl;

    if(debug_>1) cout << "From Control sample          - ratio:  " << leftNonB_/rightNonB_ <<  " $\\pm$ " << Fe_ << "   " << Fe_/(leftNonB_/rightNonB_) << endl;  
    if(debug_>1) cout << "From ControlReweigh sample   - ratio:  " << leftNonBReweigh_/rightNonBReweigh_ <<  " $\\pm$ " << FReweighe_ << "   " << FReweighe_/(leftNonBReweigh_/rightNonBReweigh_) << endl;
    if(debug_>1) cout << "From Signal sample           - ratio:  " << leftNonBMC_/rightNonBMC_  << " $\\pm$ " << FMCe_ << "   " << FMCe_/(leftNonBMC_/rightNonBMC_) << endl;
    if(debug_>1) cout << "(ControlReweig-Signal)/Signal:         " << ((leftNonBReweigh_/rightNonBReweigh_)-(leftNonBMC_/rightNonBMC_))/(leftNonBMC_/rightNonBMC_)  << " $\\pm$ " << sqrt( pow( FReweighe_/(leftNonBMC_/rightNonBMC_) ,2) + pow( (leftNonBReweigh_/rightNonBReweigh_)/((leftNonBMC_/rightNonBMC_)*(leftNonBMC_/rightNonBMC_))*FMCe_ ,2) ) << endl;

    if(debug_>1) cout << "From Signal sample (bis)          - ratio:  " << FMC_  << " $\\pm$ " << FMCe_ << "   " << FMCe_/FMC_ << endl;
    if(debug_>1) cout << "(ControlReweig-Signal)/Signal:  (bis)       " <<  (FReweigh_-FMC_)/FMC_<< " $\\pm$ " << sqrt( pow((1/FMC_)*FReweighe_,2) +pow((FReweigh_/(FMC_*FMC_))*FMCe_,2) ) << endl;



    if(debug_>1) cout << "+--> ------------------------------------------------" << endl;
	  if(debug_>1) cout << "+--> From Signal sample leftb  " << TH1Sng_Var0->Integral(lowerBin_,centralLowBin_) << " rightb "  << TH1Sng_Var0->Integral(centralUpBin_+1,upperBin_) << endl;
	  if(debug_>1) cout << "+--> From Signal sample b fraction left : " << TH1Sng_Var0->Integral(lowerBin_,centralLowBin_)/(leftNonBMC_+TH1Sng_Var0->Integral(lowerBin_,centralLowBin_)) << endl;
	  if(debug_>1) cout << "+--> From Signal sample b fraction right: "  << TH1Sng_Var0->Integral(centralUpBin_+1,upperBin_)/(TH1Sng_Var0->Integral(centralUpBin_+1,upperBin_)+rightNonBMC_) << endl;
	  if(debug_>1) cout << "+--> From Signal sample leftnb  " << TH1Bkg_Var0->Integral(lowerBin_,centralLowBin_) << " rightnb "  << TH1Bkg_Var0->Integral(centralUpBin_+1,upperBin_) << endl;
	  
    TH1Bkg_Var0->Integral(lowerBin_,centralLowBin_);
	  
	  if (debug_>1)cout << "RightRew " << TH2Bkg_Right1DYReweigh->Integral(0,TH2Bkg_Right1DYReweigh->GetNbinsX()+1) << endl;
	  if (debug_>1)cout << "Right " << TH2Bkg_Right1DY->Integral(0,TH2Bkg_Right1DY->GetNbinsX()+1) << endl;  
	  
	  if (debug_>1)cout << "Data RightRew " << TH2Data_Right1DYReweigh->Integral(0,TH2Data_Right1DYReweigh->GetNbinsX()+1) << endl;
	  if (debug_>1)cout << "Data Right " << TH2Data_Right1DY->Integral(0,TH2Data_Right1DY->GetNbinsX()+1) << endl;
	  
	  rightNonBMC_=TH2Bkg_Right1DYReweigh->Integral(0,TH2Bkg_Right1DYReweigh->GetNbinsX()+1);
	  
	  double left = TH2Data_Left1DY->Integral(0,TH2Data_Left1DY->GetNbinsX()+1); //TH1Data_ControlVar->Integral(lowerBin_,centralLowBin_);
	  double right = TH2Data_Right1DY->Integral(0,TH2Data_Right1DY->GetNbinsX()+1); //TH1Data_ControlVar->Integral(centralUpBin_+1,upperBin_);
	  
	  double Fprime = TH2Bkg_Left1DY->Integral(0,TH2Bkg_Left1DY->GetNbinsX()+1)/TH2Bkg_Right1DY->Integral(0,TH2Bkg_Right1DY->GetNbinsX()+1);
	  
	  if (debug_>1)cout << " L (" << left << ") - F (" << FMC_ << ") * R (" << right << ") = " << left-(FMC_*right) << endl; 
	  if (debug_>1)cout << " L (" << left << ") - F (" << Fprime << ") * R (" << right << ") = " << left-(Fprime*right) << endl;    	  
	  
	  if (debug_>1)cout << " L (" << left << ") - F (" << FMC_ << ") * R (" << TH1Bkg_Var0->Integral(centralUpBin_+1,upperBin_) << ") = " << left-(FMC_*TH1Bkg_Var0->Integral(centralUpBin_+1,upperBin_)) << endl; 

	  if(debug_>1) cout << "+--> --------------------------------------------------------------------+ "<< endl;
      
  }

}

vector<double> PtEtaBin::doMLJTemplateFit(string chi2cut) {
    
    //cout << lumi_ << endl; exit(1);
	
	vector<double> fitResults;
	///---------------------------------------------------------------------
	///---------------------------------------------------------------------
	//This section will use roofit to do estimate the ttbar cross section
	///---------------------------------------------------------------------
	///---------------------------------------------------------------------
	
	/*for (unsigned int i=0; i<TH1Data_ControlVarReweigh->GetNbinsX()+1; i++) {
	 
	 TH1Data_Var0_VVData_->SetBinContent(i,TH1Data_ControlVarReweigh->GetBinContent(i));
	 
	 }*/
	
	// load ttbar and vv templates from MC
	
	TString filename=""; 
	GiveName(&filename); 
    TString filename2=""; 
	GiveName(&filename2); 
	
	filename="FitTemplates/"+filename+"_Chi2Cut_"+chi2cut+"_MLJTemplates.root";
	filename2="FitTemplates/"+filename2+"_Chi2Cut_"+chi2cut+"_MLJTemplates_lum.txt";
	
	//filename="FitTemplates/"+filename+"_MLJTemplates.root";
	
	if (debug_ > 0) cout << "doMLJTemplateFit:: Loading Mlj templates from " << filename << endl;
	
	TFile* ftmp = new TFile(filename,"READ");
	
	ftmp->cd();
	
	TString tmpname = ""; GiveName(&tmpname); tmpname+="TH1Data_Var0_TTbar";
	TH1D* ttbar = (TH1D*) ftmp->Get(tmpname);
	
	tmpname = ""; GiveName(&tmpname); tmpname+="TH1Data_Var0_VVMC";
	TH1D* vv = (TH1D*) ftmp->Get(tmpname);
	
	if (!ttbar || !vv) {
		
		mkdir("FitTemplates",0777);
		
		cout << "doMLJTemplateFit:: Mlj Templates not found or corrupt, recreating " << filename << "!  Please restart the code!" << endl;
		
		ftmp->Close();
		
		TFile* ftmp = new TFile(filename,"RECREATE");
		
		histo1D["TH1Data_Var0_TTbar"]->Write();
		histo1D["TH1Sng_Var0_TTbar"]->Write();
		histo1D["TH1Bkg_Var0_TTbar"]->Write();

		histo1D["TH1Sng_Var0_VVMC"]->Write();
		histo1D["TH1Bkg_Var0_VVMC"]->Write();
		histo1D["TH1Data_Var0_VVMC"]->Write();
		
		histo1D["TH1Data_Var0_TTbar_bTagL"]->Write();
		histo1D["TH1Sng_Var0_TTbar_bTagL"]->Write();
		histo1D["TH1Bkg_Var0_TTbar_bTagL"]->Write();
		
		histo1D["TH1Data_Var0_TTbar_bTagM"]->Write();	
		histo1D["TH1Sng_Var0_TTbar_bTagM"]->Write();
		histo1D["TH1Bkg_Var0_TTbar_bTagM"]->Write();
		
		histo1D["TH1Data_Var0_TTbar_bTagT"]->Write();
		histo1D["TH1Sng_Var0_TTbar_bTagT"]->Write();
		histo1D["TH1Bkg_Var0_TTbar_bTagT"]->Write();
		
		histo1D["TH1Data_Var0_VVMC_bTagL"]->Write();
		histo1D["TH1Sng_Var0_VVMC_bTagL"]->Write();
		histo1D["TH1Bkg_Var0_VVMC_bTagL"]->Write();

		histo1D["TH1Data_Var0_VVMC_bTagM"]->Write();
		histo1D["TH1Sng_Var0_VVMC_bTagM"]->Write();
		histo1D["TH1Bkg_Var0_VVMC_bTagM"]->Write();

		histo1D["TH1Data_Var0_VVMC_bTagT"]->Write();
		histo1D["TH1Sng_Var0_VVMC_bTagT"]->Write();
		histo1D["TH1Bkg_Var0_VVMC_bTagT"]->Write();
				
		ftmp->Close();
        
        fstream lum(filename2, ios::out | ios::trunc);
        
        lum << lumi_;
        
        lum.close();
        
        doMLJTemplateFit(chi2cut);
		
	} else if (ttbar && vv) {
		
		if (debug_ > 0)cout << "doMLJTemplateFit:: Mlj Templates loaded from " << filename << "!" << endl;
		
		
		histo1D["TH1Data_Var0_TTbar"]->Delete();
		histo1D["TH1Sng_Var0_TTbar"]->Delete();
		histo1D["TH1Bkg_Var0_TTbar"]->Delete();
		
		histo1D["TH1Sng_Var0_VVMC"]->Delete();
		histo1D["TH1Bkg_Var0_VVMC"]->Delete();
		histo1D["TH1Data_Var0_VVMC"]->Delete();
		
		histo1D["TH1Data_Var0_TTbar_bTagL"]->Delete();
		histo1D["TH1Sng_Var0_TTbar_bTagL"]->Delete();
		histo1D["TH1Bkg_Var0_TTbar_bTagL"]->Delete();
		
		histo1D["TH1Data_Var0_TTbar_bTagM"]->Delete();
		histo1D["TH1Sng_Var0_TTbar_bTagM"]->Delete();
		histo1D["TH1Bkg_Var0_TTbar_bTagM"]->Delete();
		
		histo1D["TH1Data_Var0_TTbar_bTagT"]->Delete();
		histo1D["TH1Sng_Var0_TTbar_bTagT"]->Delete();
		histo1D["TH1Bkg_Var0_TTbar_bTagT"]->Delete();

		histo1D["TH1Data_Var0_VVMC_bTagL"]->Delete();
		histo1D["TH1Data_Var0_VVMC_bTagM"]->Delete();
		histo1D["TH1Data_Var0_VVMC_bTagT"]->Delete();

		histo1D["TH1Data_Var0_TTbar"] = (TH1D*) ttbar->Clone();
		
		TString tmpname = ""; GiveName(&tmpname); tmpname+="TH1Sng_Var0_TTbar";
		histo1D["TH1Sng_Var0_TTbar"] = (TH1D*) ftmp->Get(tmpname)->Clone();
		
		tmpname = ""; GiveName(&tmpname); tmpname+="TH1Bkg_Var0_TTbar";
		histo1D["TH1Bkg_Var0_TTbar"] = (TH1D*) ftmp->Get(tmpname)->Clone();

		tmpname = ""; GiveName(&tmpname); tmpname+="TH1Data_Var0_VVMC";
		histo1D["TH1Data_Var0_VVMC"] = (TH1D*) vv->Clone();

		tmpname = ""; GiveName(&tmpname); tmpname+="TH1Sng_Var0_VVMC";
		histo1D["TH1Sng_Var0_VVMC"] = (TH1D*) vv->Clone();

		tmpname = ""; GiveName(&tmpname); tmpname+="TH1Bkg_Var0_VVMC";
		histo1D["TH1Bkg_Var0_VVMC"] = (TH1D*) vv->Clone();

		tmpname = ""; GiveName(&tmpname); tmpname+="TH1Data_Var0_TTbar_bTagL";
		histo1D["TH1Data_Var0_TTbar_bTagL"] = (TH1D*) ftmp->Get(tmpname)->Clone();
		tmpname = ""; GiveName(&tmpname); tmpname+="TH1Sng_Var0_TTbar_bTagL";
		histo1D["TH1Sng_Var0_TTbar_bTagL"] = (TH1D*) ftmp->Get(tmpname)->Clone();
		tmpname = ""; GiveName(&tmpname); tmpname+="TH1Bkg_Var0_TTbar_bTagL";
		histo1D["TH1Bkg_Var0_TTbar_bTagL"] = (TH1D*) ftmp->Get(tmpname)->Clone();
		
		tmpname = ""; GiveName(&tmpname); tmpname+="TH1Data_Var0_TTbar_bTagM";
		histo1D["TH1Data_Var0_TTbar_bTagM"] = (TH1D*) ftmp->Get(tmpname)->Clone();		
		tmpname = ""; GiveName(&tmpname); tmpname+="TH1Sng_Var0_TTbar_bTagM";
		histo1D["TH1Sng_Var0_TTbar_bTagM"] = (TH1D*) ftmp->Get(tmpname)->Clone();		
		tmpname = ""; GiveName(&tmpname); tmpname+="TH1Bkg_Var0_TTbar_bTagM";
		histo1D["TH1Bkg_Var0_TTbar_bTagM"] = (TH1D*) ftmp->Get(tmpname)->Clone();
		
		tmpname = ""; GiveName(&tmpname); tmpname+="TH1Data_Var0_TTbar_bTagT";
		histo1D["TH1Data_Var0_TTbar_bTagT"] = (TH1D*) ftmp->Get(tmpname)->Clone();		
		tmpname = ""; GiveName(&tmpname); tmpname+="TH1Sng_Var0_TTbar_bTagT";
		histo1D["TH1Sng_Var0_TTbar_bTagT"] = (TH1D*) ftmp->Get(tmpname)->Clone();		
		tmpname = ""; GiveName(&tmpname); tmpname+="TH1Bkg_Var0_TTbar_bTagT";
		histo1D["TH1Bkg_Var0_TTbar_bTagT"] = (TH1D*) ftmp->Get(tmpname)->Clone();
		
		tmpname = ""; GiveName(&tmpname); tmpname+="TH1Data_Var0_VVMC_bTagL";
		histo1D["TH1Data_Var0_VVMC_bTagL"] = (TH1D*) ftmp->Get(tmpname)->Clone();		
		tmpname = ""; GiveName(&tmpname); tmpname+="TH1Data_Var0_VVMC_bTagM";
		histo1D["TH1Data_Var0_VVMC_bTagM"] = (TH1D*) ftmp->Get(tmpname)->Clone();		
		tmpname = ""; GiveName(&tmpname); tmpname+="TH1Data_Var0_VVMC_bTagT";
		histo1D["TH1Data_Var0_VVMC_bTagT"] = (TH1D*) ftmp->Get(tmpname)->Clone();
		//delete ttbar;
		
		//delete vv;
		//ftmp->Close();
        
        
        fstream lum(filename2, ios::in );
        
        string line;
        
        while (! lum.eof())
            lum >> line;
        
        templateLumi_ = (double)atof(line.c_str());
        
        lum.close();

	}
	
	int b=histo1D["TH1Sng_Var0_TTbar"]->Integral(0,histo1D["TH1Sng_Var0_TTbar"]->GetNbinsX()+1);
	int nb=histo1D["TH1Bkg_Var0_TTbar"]->Integral(0,histo1D["TH1Bkg_Var0_TTbar"]->GetNbinsX()+1);
	
	float fractionBtotal=(float)b/((float)b+(float)nb);
	
	if (debug_ > 0) cout << "!!! b " << b << " nb " << nb << endl;
	
	int nb_btagL=histo1D["TH1Bkg_Var0_TTbar_bTagL"]->Integral(0,histo1D["TH1Bkg_Var0_TTbar_bTagL"]->GetNbinsX()+1);
	int nb_btagM=histo1D["TH1Bkg_Var0_TTbar_bTagM"]->Integral(0,histo1D["TH1Bkg_Var0_TTbar_bTagM"]->GetNbinsX()+1);
	int nb_btagT=histo1D["TH1Bkg_Var0_TTbar_bTagT"]->Integral(0,histo1D["TH1Bkg_Var0_TTbar_bTagT"]->GetNbinsX()+1);
	
	float misTagRateL=(float)nb_btagL/(float)nb;
	float misTagRateM=(float)nb_btagM/(float)nb;
	float misTagRateT=(float)nb_btagT/(float)nb;
		
	if (debug_ > 0) cout << "!!! nb_btagL " << nb_btagL << " nb " << nb << endl;
	
	// fit without using btag
	
	if (debug_ > 0) cout << "doMLJTemplateFit: * Performing NOMINAL fits without BTAG cut" << endl;
	
	vector<float> tmp = doTemplateFit(histo1D["TH1Data_Var0_TTbar"],histo1D["TH1Data_Var0_VVMC"],histo1D["TH1Data_Var0_VVData"],TH1Data_Var0_XS,(TString)"Fit_NOCUT");
	for (unsigned int t=0;t<tmp.size();t++) fitResults.push_back(tmp[t]);
	  
	fitResults.push_back(0); // no mistagrate

	if (debug_ > 0) cout << "doMLJTemplateFit: * Performing BTAGCUT fits with BTAG cut LOOSE" << endl;
	
	//histo1D["TH1Data_Var0_bTagL"])
	
	if (debug_ > 0) cout << "bTagL Mistag rate: " << misTagRateL << " b/(nb+b): " << fractionBtotal << endl;
	if (debug_ > 0)cout << "bTagL eff TTbar: " << histo1D["TH1Sng_Var0_TTbar_bTagL"]->Integral()/histo1D["TH1Sng_Var0_TTbar"]->Integral() << endl;
	if (debug_ > 0) cout << "bTagL eff Other: " << histo1D["TH1Sng_Var0_VVMC_bTagL"]->Integral()/histo1D["TH1Sng_Var0_VVMC"]->Integral() << endl;
	
	if (debug_ > 0) cout << "#ttbar " << histo1D["TH1Data_Var0_TTbar_bTagL"]->Integral() << " #non-ttbar " << histo1D["TH1Data_Var0_VVMC_bTagL"]->Integral() << " data " << histo1D["TH1Data_Var0_bTagL"]->Integral() << endl;
	
	tmp.clear();tmp = doTemplateFit(histo1D["TH1Data_Var0_TTbar_bTagL"],histo1D["TH1Data_Var0_VVMC_bTagL"],histo1D["TH1Data_Var0_VVData"],histo1D["TH1Data_Var0_bTagL"],(TString)"Fit_BtagLCut");
	//tmp.clear(); tmp = doTemplateFit(histo1D["TH1Data_Var0_TTbar_bTagL"],histo1D["TH1Data_Var0_VVMC_bTagL"],histo1D["TH1Data_Var0_VVData"],histo1D["TH1Data_Var0_bTagL"]);
	for (unsigned int t=0;t<tmp.size();t++) fitResults.push_back(tmp[t]);
	
	fitResults.push_back(misTagRateL);
	
	if (debug_ > 0) cout << "doMLJTemplateFit: * Performing BTAGCUT fits with BTAG cut MEDIUM" << endl;
	
	if (debug_ > 0) cout << "#ttbar " << histo1D["TH1Data_Var0_TTbar_bTagM"]->Integral() << " #non-ttbar " << histo1D["TH1Data_Var0_VVMC_bTagM"]->Integral() << " data " << histo1D["TH1Data_Var0_bTagM"]->Integral() << endl;
	
	if (debug_ > 0) cout << "bTagM Mistag rate: " << misTagRateM << " b/(nb+b): " << fractionBtotal << endl;

	tmp.clear(); tmp = doTemplateFit(histo1D["TH1Data_Var0_TTbar_bTagM"],histo1D["TH1Data_Var0_VVMC_bTagM"],histo1D["TH1Data_Var0_VVData"],histo1D["TH1Data_Var0_bTagM"],(TString)"Fit_BtagMCut");
	for (unsigned int t=0;t<tmp.size();t++) fitResults.push_back(tmp[t]);
	
	fitResults.push_back(misTagRateM);

	if (debug_ > 0) cout << "doMLJTemplateFit: * Performing BTAGCUT fits with BTAG cut TIGHT" << endl;
	
	if (debug_ > 0) cout << "#ttbar " << histo1D["TH1Data_Var0_TTbar_bTagM"]->Integral() << " #non-ttbar " << histo1D["TH1Data_Var0_VVMC_bTagM"]->Integral() << " data " << histo1D["TH1Data_Var0_bTagM"]->Integral() << endl;
	
	if (debug_ > 0) cout << "bTagT Mistag rate: " << misTagRateT << " b/(nb+b): " << fractionBtotal << endl;

	tmp.clear(); tmp = doTemplateFit(histo1D["TH1Data_Var0_TTbar_bTagT"],histo1D["TH1Data_Var0_VVMC_bTagT"],histo1D["TH1Data_Var0_VVData"],histo1D["TH1Data_Var0_bTagT"],(TString)"Fit_BtagTCut");
	for (unsigned int t=0;t<tmp.size();t++) fitResults.push_back(tmp[t]); 

	fitResults.push_back(misTagRateT);

	fitResults.push_back(fractionBtotal);

	return fitResults;
}

void PtEtaBin::MeasureEff(bool doSCreweigh){

  GiveName(&titleData_BtagEffMeasured_); titleData_BtagEffMeasured_+="TH1Data_BtagEffMeasured"; 
  GiveName(&titleData_BtagMeasured_); titleData_BtagMeasured_+="TH1Data_BtagMeasured"; 
  GiveName(&titleData_BtagEffMCMeasured_); titleData_BtagEffMCMeasured_+="TH1Data_BtagEffMCMeasured"; 
  GiveName(&titleData_BtagMCMeasured_); titleData_BtagMCMeasured_+="TH1Data_BtagMCMeasured"; 
 
  if(!varBinSize_){
    TH1Data_BtagEffMeasured = new TH1D(titleData_BtagEffMeasured_,titleData_BtagEffMeasured_,TH2Data_Left->GetYaxis()->GetNbins(),TH2Data_Left->GetYaxis()->GetBinLowEdge(1),TH2Data_Left->GetYaxis()->GetBinUpEdge(TH2Data_Left->GetYaxis()->GetNbins()));
    TH1Data_BtagMeasured = new TH1D(titleData_BtagMeasured_,titleData_BtagMeasured_,TH2Data_Left->GetYaxis()->GetNbins(),TH2Data_Left->GetYaxis()->GetBinLowEdge(1),TH2Data_Left->GetYaxis()->GetBinUpEdge(TH2Data_Left->GetYaxis()->GetNbins()));
    TH1Data_BtagEffMCMeasured = new TH1D(titleData_BtagEffMCMeasured_,titleData_BtagEffMCMeasured_,TH2Data_Left->GetYaxis()->GetNbins(),TH2Data_Left->GetYaxis()->GetBinLowEdge(1),TH2Data_Left->GetYaxis()->GetBinUpEdge(TH2Data_Left->GetYaxis()->GetNbins()));
    TH1Data_BtagMCMeasured = new TH1D(titleData_BtagMCMeasured_,titleData_BtagMCMeasured_,TH2Data_Left->GetYaxis()->GetNbins(),TH2Data_Left->GetYaxis()->GetBinLowEdge(1),TH2Data_Left->GetYaxis()->GetBinUpEdge(TH2Data_Left->GetYaxis()->GetNbins()));
  }      

  if(varBinSize_){
    TH1Data_BtagEffMeasured = new TH1D(titleData_BtagEffMeasured_,titleData_BtagEffMeasured_,TH2Data_Left->GetYaxis()->GetNbins(),rangesbTag_);
    TH1Data_BtagMeasured = new TH1D(titleData_BtagMeasured_,titleData_BtagMeasured_,TH2Data_Left->GetYaxis()->GetNbins(),rangesbTag_);
    TH1Data_BtagEffMCMeasured = new TH1D(titleData_BtagEffMCMeasured_,titleData_BtagEffMCMeasured_,TH2Data_Left->GetYaxis()->GetNbins(),rangesbTag_);
    TH1Data_BtagMCMeasured = new TH1D(titleData_BtagMCMeasured_,titleData_BtagMCMeasured_,TH2Data_Left->GetYaxis()->GetNbins(),rangesbTag_);
  }
  
  //diff plots
  GiveName(&titleData_BtagEffMeasuredDiff_); titleData_BtagEffMeasuredDiff_+="TH1Data_BtagEffMeasuredDiff"; 
  GiveName(&titleData_BtagEffMCMeasuredDiff_); titleData_BtagEffMCMeasuredDiff_+="TH1Data_BtagEffMCMeasuredDiff"; 
  
  if(!varBinSize_){
    TH1Data_BtagEffMeasuredDiff = new TH1D(titleData_BtagEffMeasuredDiff_,titleData_BtagEffMeasuredDiff_,TH2Data_Left->GetYaxis()->GetNbins(),TH2Data_Left->GetYaxis()->GetBinLowEdge(1),TH2Data_Left->GetYaxis()->GetBinUpEdge(TH2Data_Left->GetYaxis()->GetNbins()));
    TH1Data_BtagEffMCMeasuredDiff = new TH1D(titleData_BtagEffMCMeasuredDiff_,titleData_BtagEffMCMeasuredDiff_,TH2Data_Left->GetYaxis()->GetNbins(),TH2Data_Left->GetYaxis()->GetBinLowEdge(1),TH2Data_Left->GetYaxis()->GetBinUpEdge(TH2Data_Left->GetYaxis()->GetNbins()));
  }
  if(varBinSize_){
    TH1Data_BtagEffMeasuredDiff = new TH1D(titleData_BtagEffMeasuredDiff_,titleData_BtagEffMeasuredDiff_,TH2Data_Left->GetYaxis()->GetNbins(),rangesbTag_);
    TH1Data_BtagEffMCMeasuredDiff = new TH1D(titleData_BtagEffMCMeasuredDiff_,titleData_BtagEffMCMeasuredDiff_,TH2Data_Left->GetYaxis()->GetNbins(),rangesbTag_);
  }

  TH1Data_BtagEffMeasured->Sumw2();
  TH1Data_BtagMeasured->Sumw2();
  TH1Data_BtagEffMCMeasured->Sumw2();
  TH1Data_BtagMCMeasured->Sumw2();

  TH1Data_BtagEffMeasuredDiff->Sumw2();
  TH1Data_BtagEffMCMeasuredDiff->Sumw2();

  EffCalculation(doSCreweigh, TH2Data_Left1DY,TH2Data_Right1DY,TH1Data_BtagMeasured,TH1Data_BtagEffMeasured,TH1Data_BtagMCMeasured,TH1Data_BtagEffMCMeasured);

  DiffCalculation(TH1Sng_BtagEffMC,TH1Data_BtagEffMeasured,TH1Data_BtagEffMeasuredDiff);
  DiffCalculation(TH1Sng_BtagEffMC,TH1Data_BtagEffMCMeasured,TH1Data_BtagEffMCMeasuredDiff);

}

void PtEtaBin::MeasureEffLR(bool doSCreweigh){//obsolete (did not implement the variable bin size)

  GiveName(&titleData_BtagEffMeasuredLR_); titleData_BtagEffMeasuredLR_+="TH1Data_BtagEffMeasuredLR";
  GiveName(&titleData_BtagMeasuredLR_); titleData_BtagMeasuredLR_+="TH1Data_BtagMeasuredLR";
  GiveName(&titleData_BtagEffMCMeasuredLR_); titleData_BtagEffMCMeasuredLR_+="TH1Data_BtagEffMCMeasuredLR";
  GiveName(&titleData_BtagMCMeasuredLR_); titleData_BtagMCMeasuredLR_+="TH1Data_BtagMCMeasuredLR";

  TH1Data_BtagEffMeasuredLR = new TH1D(titleData_BtagEffMeasuredLR_,titleData_BtagEffMeasuredLR_,TH1Data_BtagAll->GetXaxis()->GetNbins(),TH1Data_BtagAll->GetXaxis()->GetBinLowEdge(1),TH1Data_BtagAll->GetXaxis()->GetBinUpEdge(TH1Data_BtagAll->GetXaxis()->GetNbins()));
  TH1Data_BtagMeasuredLR = new TH1D(titleData_BtagMeasuredLR_,titleData_BtagMeasuredLR_,TH1Data_BtagAll->GetXaxis()->GetNbins(),TH1Data_BtagAll->GetXaxis()->GetBinLowEdge(1),TH1Data_BtagAll->GetXaxis()->GetBinUpEdge(TH1Data_BtagAll->GetXaxis()->GetNbins()));
  TH1Data_BtagEffMCMeasuredLR = new TH1D(titleData_BtagEffMCMeasuredLR_,titleData_BtagEffMCMeasuredLR_,TH1Data_BtagAll->GetXaxis()->GetNbins(),TH1Data_BtagAll->GetXaxis()->GetBinLowEdge(1),TH1Data_BtagAll->GetXaxis()->GetBinUpEdge(TH1Data_BtagAll->GetXaxis()->GetNbins()));
  TH1Data_BtagMCMeasuredLR = new TH1D(titleData_BtagMCMeasuredLR_,titleData_BtagMCMeasuredLR_,TH1Data_BtagAll->GetXaxis()->GetNbins(),TH1Data_BtagAll->GetXaxis()->GetBinLowEdge(1),TH1Data_BtagAll->GetXaxis()->GetBinUpEdge(TH1Data_BtagAll->GetXaxis()->GetNbins()));

  //diff plots
  GiveName(&titleData_BtagEffMeasuredLRDiff_); titleData_BtagEffMeasuredLRDiff_+="TH1Data_BtagEffMeasuredLRDiff";
  GiveName(&titleData_BtagEffMCMeasuredLRDiff_); titleData_BtagEffMCMeasuredLRDiff_+="TH1Data_BtagEffMCMeasuredLRDiff";
  
  TH1Data_BtagEffMeasuredLRDiff = new TH1D(titleData_BtagEffMeasuredLRDiff_,titleData_BtagEffMeasuredLRDiff_,TH1Data_BtagAll->GetXaxis()->GetNbins(),TH1Data_BtagAll->GetXaxis()->GetBinLowEdge(1),TH1Data_BtagAll->GetXaxis()->GetBinUpEdge(TH1Data_BtagAll->GetXaxis()->GetNbins()));
  TH1Data_BtagEffMCMeasuredLRDiff = new TH1D(titleData_BtagEffMCMeasuredLRDiff_,titleData_BtagEffMCMeasuredLRDiff_,TH1Data_BtagAll->GetXaxis()->GetNbins(),TH1Data_BtagAll->GetXaxis()->GetBinLowEdge(1),TH1Data_BtagAll->GetXaxis()->GetBinUpEdge(TH1Data_BtagAll->GetXaxis()->GetNbins()));

  TH1Data_BtagEffMeasuredLR->Sumw2();
  TH1Data_BtagMeasuredLR->Sumw2();
  TH1Data_BtagEffMCMeasuredLR->Sumw2();
  TH1Data_BtagMCMeasuredLR->Sumw2();
  TH1Data_BtagEffMeasuredLRDiff->Sumw2();
  TH1Data_BtagEffMCMeasuredLRDiff->Sumw2();
  
  //  EffCalculation(doSCreweigh, TH2Data_Left1DYReweigh,TH2Data_Right1DY,TH1Data_BtagMeasuredLR,TH1Data_BtagEffMeasuredLR,TH1Data_BtagMCMeasuredLR,TH1Data_BtagEffMCMeasuredLR);

  DiffCalculation(TH1Sng_BtagEffMC,TH1Data_BtagEffMeasuredLR,TH1Data_BtagEffMeasuredLRDiff);
  DiffCalculation(TH1Sng_BtagEffMC,TH1Data_BtagEffMCMeasuredLR,TH1Data_BtagEffMCMeasuredLRDiff);


} // obsolete

void PtEtaBin::MeasureEffRR(bool doSCreweigh){

  GiveName(&titleData_BtagEffMeasuredRR_); titleData_BtagEffMeasuredRR_+="TH1Data_BtagEffMeasuredRR";
  GiveName(&titleData_BtagMeasuredRR_); titleData_BtagMeasuredRR_+="TH1Data_BtagMeasuredRR";
  GiveName(&titleData_BtagEffMCMeasuredRR_); titleData_BtagEffMCMeasuredRR_+="TH1Data_BtagEffMCMeasuredRR";
  GiveName(&titleData_BtagMCMeasuredRR_); titleData_BtagMCMeasuredRR_+="TH1Data_BtagMCMeasuredRR";
  
  if(!varBinSize_){
    TH1Data_BtagEffMeasuredRR = new TH1D(titleData_BtagEffMeasuredRR_,titleData_BtagEffMeasuredRR_,TH1Data_BtagAll->GetXaxis()->GetNbins(),TH1Data_BtagAll->GetXaxis()->GetBinLowEdge(1),TH1Data_BtagAll->GetXaxis()->GetBinUpEdge(TH1Data_BtagAll->GetXaxis()->GetNbins()));
    TH1Data_BtagMeasuredRR = new TH1D(titleData_BtagMeasuredRR_,titleData_BtagMeasuredRR_,TH1Data_BtagAll->GetXaxis()->GetNbins(),TH1Data_BtagAll->GetXaxis()->GetBinLowEdge(1),TH1Data_BtagAll->GetXaxis()->GetBinUpEdge(TH1Data_BtagAll->GetXaxis()->GetNbins()));
    TH1Data_BtagEffMCMeasuredRR = new TH1D(titleData_BtagEffMCMeasuredRR_,titleData_BtagEffMCMeasuredRR_,TH1Data_BtagAll->GetXaxis()->GetNbins(),TH1Data_BtagAll->GetXaxis()->GetBinLowEdge(1),TH1Data_BtagAll->GetXaxis()->GetBinUpEdge(TH1Data_BtagAll->GetXaxis()->GetNbins()));
    TH1Data_BtagMCMeasuredRR = new TH1D(titleData_BtagMCMeasuredRR_,titleData_BtagMCMeasuredRR_,TH1Data_BtagAll->GetXaxis()->GetNbins(),TH1Data_BtagAll->GetXaxis()->GetBinLowEdge(1),TH1Data_BtagAll->GetXaxis()->GetBinUpEdge(TH1Data_BtagAll->GetXaxis()->GetNbins()));
  }
  if(varBinSize_){
    TH1Data_BtagEffMeasuredRR = new TH1D(titleData_BtagEffMeasuredRR_,titleData_BtagEffMeasuredRR_,TH1Data_BtagAll->GetXaxis()->GetNbins(),rangesbTag_);
    TH1Data_BtagMeasuredRR = new TH1D(titleData_BtagMeasuredRR_,titleData_BtagMeasuredRR_,TH1Data_BtagAll->GetXaxis()->GetNbins(),rangesbTag_);
    TH1Data_BtagEffMCMeasuredRR = new TH1D(titleData_BtagEffMCMeasuredRR_,titleData_BtagEffMCMeasuredRR_,TH1Data_BtagAll->GetXaxis()->GetNbins(),rangesbTag_);
    TH1Data_BtagMCMeasuredRR = new TH1D(titleData_BtagMCMeasuredRR_,titleData_BtagMCMeasuredRR_,TH1Data_BtagAll->GetXaxis()->GetNbins(),rangesbTag_);
  }
  
  //diff plots
  GiveName(&titleData_BtagEffMeasuredRRDiff_); titleData_BtagEffMeasuredRRDiff_+="TH1Data_BtagEffMeasuredRRDiff";
  GiveName(&titleData_BtagEffMCMeasuredRRDiff_); titleData_BtagEffMCMeasuredRRDiff_+="TH1Data_BtagEffMCMeasuredRRDiff";
  
  if(!varBinSize_){
    TH1Data_BtagEffMeasuredRRDiff = new TH1D(titleData_BtagEffMeasuredRRDiff_,titleData_BtagEffMeasuredRRDiff_,TH1Data_BtagAll->GetXaxis()->GetNbins(),TH1Data_BtagAll->GetXaxis()->GetBinLowEdge(1),TH1Data_BtagAll->GetXaxis()->GetBinUpEdge(TH1Data_BtagAll->GetXaxis()->GetNbins()));
    TH1Data_BtagEffMCMeasuredRRDiff = new TH1D(titleData_BtagEffMCMeasuredRRDiff_,titleData_BtagEffMCMeasuredRRDiff_,TH1Data_BtagAll->GetXaxis()->GetNbins(),TH1Data_BtagAll->GetXaxis()->GetBinLowEdge(1),TH1Data_BtagAll->GetXaxis()->GetBinUpEdge(TH1Data_BtagAll->GetXaxis()->GetNbins()));
  }
  if(varBinSize_){
    TH1Data_BtagEffMeasuredRRDiff = new TH1D(titleData_BtagEffMeasuredRRDiff_,titleData_BtagEffMeasuredRRDiff_,TH1Data_BtagAll->GetXaxis()->GetNbins(),rangesbTag_);
    TH1Data_BtagEffMCMeasuredRRDiff = new TH1D(titleData_BtagEffMCMeasuredRRDiff_,titleData_BtagEffMCMeasuredRRDiff_,TH1Data_BtagAll->GetXaxis()->GetNbins(),rangesbTag_);
  }

  TH1Data_BtagEffMeasuredRR->Sumw2();
  TH1Data_BtagMeasuredRR->Sumw2();
  TH1Data_BtagEffMCMeasuredRR->Sumw2();
  TH1Data_BtagMCMeasuredRR->Sumw2();
  TH1Data_BtagEffMeasuredRRDiff->Sumw2();
  TH1Data_BtagEffMCMeasuredRRDiff->Sumw2();
  
  if(!doShift_){
    EffCalculation(doSCreweigh, TH2Data_Left1DY,TH2Data_Right1DYReweigh,TH1Data_BtagMeasuredRR,TH1Data_BtagEffMeasuredRR,TH1Data_BtagMCMeasuredRR,TH1Data_BtagEffMCMeasuredRR);
  }

  if(doShift_){
    //I'm now making the plot for shifting both left and right
    EffCalculation(doSCreweigh, TH2Data_Left1DYReweigh,TH2Data_Right1DYReweigh,TH1Data_BtagMeasuredRR,TH1Data_BtagEffMeasuredRR,TH1Data_BtagMCMeasuredRR,TH1Data_BtagEffMCMeasuredRR);
  }
	
	//cout << EffCalcDetails_["F"] << endl;
	
	/*for (std::map<string,float>::const_iterator it = EffCalcDetails_.begin(); it != EffCalcDetails_.end(); ++it)
		cout << it->first << " -> " << it->second << endl;
	
	exit(0);
*/
  DiffCalculation(TH1Sng_BtagEffMC,TH1Data_BtagEffMeasuredRR,TH1Data_BtagEffMeasuredRRDiff);
  DiffCalculation(TH1Sng_BtagEffMC,TH1Data_BtagEffMCMeasuredRR,TH1Data_BtagEffMCMeasuredRRDiff);

}

void PtEtaBin::CoutWPEff(bool RRincluded, bool doSCreweigh, double wp, double* WParray, bool pseudo, bool more, int ptbin, int etabin, double runNb){
  
  bool wpFound=false;
  
  //cout << "wp: " << wp << endl;

  for(int i=0; i<TH1Data_BtagEffMeasuredRR->GetXaxis()->GetNbins(); i++){
    //if(fabs(wp-TH1Data_BtagEffMeasuredRR->GetXaxis()->GetBinUpEdge(i))<0.0001){
    if(fabs(wp-TH1Data_BtagEffMeasuredRR->GetXaxis()->GetBinLowEdge(i))<0.0001){
      wpFound=true; 

      double bused=TH1Data_BtagMeasured->Integral(0,TH1Data_BtagMeasured->GetNbinsX()+1);

      cout << runNb << ") " << ptbin << " " << etabin << " " << wp << " FMC_ " << FMC_ << " \\pm " << FMCe_ << endl;
      cout << runNb << ") " << ptbin << " " << etabin << " " << wp << " F_ " << F_ << " \\pm " << Fe_ << endl;
      cout << runNb << ") " << ptbin << " " << etabin << " " << wp << " #b's used " << bused << endl;
      cout << runNb << ") " << ptbin << " " << etabin << " " << wp << " FReweigh_ " << FReweigh_ << " \\pm " << FReweighe_ << endl;
      cout << runNb << ") " << ptbin << " " << etabin << " " << wp << " b-tag eff (All): " << TH1Sng_BtagEffAll->GetBinContent(i) << " \\pm " << TH1Sng_BtagEffAll->GetBinError(i) << endl;
      cout << runNb << ") " << ptbin << " " << etabin << " " << wp << " b-tag eff (MC): " << TH1Sng_BtagEffMC->GetBinContent(i) << " \\pm " << TH1Sng_BtagEffMC->GetBinError(i) << endl;
      cout << runNb << ") " << ptbin << " " << etabin << " " << wp << " b-tag eff (FMC) : " << TH1Data_BtagEffMCMeasured->GetBinContent(i) << " \\pm " << TH1Data_BtagEffMCMeasured->GetBinError(i) << endl;
      cout << runNb << ") " << ptbin << " " << etabin << " " << wp << " b-tag eff (Fdata) : " << TH1Data_BtagEffMeasured->GetBinContent(i) << " \\pm " << TH1Data_BtagEffMeasured->GetBinError(i) << endl;
      if(RRincluded) cout << runNb << ") " << ptbin << " " << etabin << " " << wp << " b-tag eff (RR+FMC) : " << TH1Data_BtagEffMCMeasuredRR->GetBinContent(i) << " \\pm " << TH1Data_BtagEffMCMeasuredRR->GetBinError(i) << endl;
      if(RRincluded) cout << runNb << ") " << ptbin << " " << etabin << " " << wp << " b-tag eff (RR+Fdata) : " << TH1Data_BtagEffMeasuredRR->GetBinContent(i) << " \\pm " << TH1Data_BtagEffMeasuredRR->GetBinError(i) << endl;
   
    }
  }
}

void PtEtaBin::GetWPEff(bool RRincluded, bool doSCreweigh, double wp, double* WParray, bool pseudo, bool more, double runNb){

  bool wpFound=false;
  
  double Fused=0;
  double Fusede=0;

  if(doSCreweigh){
    Fused=FReweigh_;
    Fusede=FReweighe_;
  } else {
    Fused=F_;
    Fusede=Fe_;
  }
  
  for(int i=0; i<22; i++){
    WParray[i]=-1;
  }

  for(int i=1; i<TH1Data_BtagEffMeasuredRR->GetXaxis()->GetNbins()+1; i++){
	  //cout << "BinLowEdge " << TH1Data_BtagEffMeasuredRR->GetXaxis()->GetBinLowEdge(i) << " wp " << wp << " diff " << wp-TH1Data_BtagEffMeasuredRR->GetXaxis()->GetBinLowEdge(i) << endl;
    //if(fabs(wp-TH1Data_BtagEffMeasuredRR->GetXaxis()->GetBinUpEdge(i))<0.0001){
    if(fabs(wp-TH1Data_BtagEffMeasuredRR->GetXaxis()->GetBinLowEdge(i))<0.0001){
      wpFound=true;
        //WParray[0] = TH1Sng_BtagEffMC->GetBinContent(i);
        //WParray[1] = TH1Sng_BtagEffMC->GetBinError(i);
        WParray[0] = TH1Sng_BtagEffAll->GetBinContent(i);
        WParray[1] = TH1Sng_BtagEffAll->GetBinError(i);
      WParray[2] = TH1Data_BtagEffMCMeasured->GetBinContent(i);
      WParray[3] = TH1Data_BtagEffMCMeasured->GetBinError(i);
      WParray[4] = TH1Data_BtagEffMeasured->GetBinContent(i);
      WParray[5] = TH1Data_BtagEffMeasured->GetBinError(i);
      if(RRincluded) WParray[6] = TH1Data_BtagEffMCMeasuredRR->GetBinContent(i);
      if(RRincluded) WParray[7] = TH1Data_BtagEffMCMeasuredRR->GetBinError(i);
      if(RRincluded) WParray[8] = TH1Data_BtagEffMeasuredRR->GetBinContent(i);
      if(RRincluded) WParray[9] = TH1Data_BtagEffMeasuredRR->GetBinError(i);
      WParray[10] = FMC_;
      WParray[11] = FMCe_;
      WParray[12] = F_;
      WParray[13] = Fe_;
      WParray[14] = FReweigh_;
      WParray[15] = FReweighe_;
      WParray[16] = 0;
      WParray[17] = 0;
      WParray[18] = TFData_LeftRight1DXFit->GetParameter(0);
      WParray[19] = 0;
      WParray[20] = TFData_LeftRight1DXFit->GetParameter(1);
      WParray[21] = 0;
      //extend with parameters of SCreweigh-fit (or do this unfitted)
		
		//cout << "-------------++++++++++++++++++ " << WParray[8] << endl;
      
      //cout << TH1Sng_Btag->Integral(i,TH1Sng_Btag->GetNbinsX()+1) << " " << TH1Sng_Btag->Integral(0,TH1Sng_Btag->GetNbinsX()+1) << " " << TH1Sng_Btag->Integral(i,TH1Sng_Btag->GetNbinsX()+1)/TH1Sng_Btag->Integral(0,TH1Sng_Btag->GetNbinsX()+1) << endl;
      //cout << TH1Sng_Btag->GetNbinsX()+1 << " " << TH1Data_BtagMCMeasuredRR->GetNbinsX()+1 << endl;
     
      /*if(more){ //temporary
	cout << runNb << ") #jets  " << TH1Data_BtagMCMeasured->Integral(i,TH1Data_BtagMCMeasured->GetNbinsX()+1) << " " << TH1Data_BtagMCMeasured->Integral(0,TH1Data_BtagMCMeasured->GetNbinsX()+1) << " " << TH1Data_BtagMCMeasured->Integral(i,TH1Data_BtagMCMeasured->GetNbinsX()+1)/TH1Data_BtagMCMeasured->Integral(0,TH1Data_BtagMCMeasured->GetNbinsX()+1) << endl;
	cout << runNb << ") #jets  " << TH1Data_BtagMeasured->Integral(i,TH1Data_BtagMeasured->GetNbinsX()+1) << " " << TH1Data_BtagMeasured->Integral(0,TH1Data_BtagMeasured->GetNbinsX()+1) << " " << TH1Data_BtagMeasured->Integral(i,TH1Data_BtagMeasured->GetNbinsX()+1)/TH1Data_BtagMeasured->Integral(0,TH1Data_BtagMeasured->GetNbinsX()+1) << endl;
	cout << runNb << ") #jets  " << TH1Data_BtagMCMeasuredRR->Integral(i,TH1Data_BtagMCMeasuredRR->GetNbinsX()+1) << " " << TH1Data_BtagMCMeasuredRR->Integral(0,TH1Data_BtagMCMeasuredRR->GetNbinsX()+1) << " " << TH1Data_BtagMCMeasuredRR->Integral(i,TH1Data_BtagMCMeasuredRR->GetNbinsX()+1)/TH1Data_BtagMCMeasuredRR->Integral(0,TH1Data_BtagMCMeasuredRR->GetNbinsX()+1) << endl;
	cout << runNb << ") #jets  " << TH1Data_BtagMeasuredRR->Integral(i,TH1Data_BtagMeasuredRR->GetNbinsX()+1) << " " << TH1Data_BtagMeasuredRR->Integral(0,TH1Data_BtagMeasuredRR->GetNbinsX()+1) << " " << TH1Data_BtagMeasuredRR->Integral(i,TH1Data_BtagMeasuredRR->GetNbinsX()+1)/TH1Data_BtagMeasuredRR->Integral(0,TH1Data_BtagMeasuredRR->GetNbinsX()+1) << endl;
      }*/

    }
  }



  if(!wpFound) cout << runNb << ") WARNING::PtEtaBin::GetWPEffRR : The working point has not been found, the returned array is corrupt" << endl;
  if(!pseudo) {
    if(more){
      cout << runNb << ") Frw: " << WParray[14] << " +- " << WParray[15] << "     F(MC): " << WParray[10] << " +- " << WParray[11] << endl;
      cout << runNb << ") par0: " << WParray[18] << " +- " << WParray[19] << "  par1: " << WParray[20] << " +- " <<  WParray[21]<< endl;
    }
    if(doShift_) cout << "warning, the measured values are for the MC shift" << endl;
    cout << runNb << ") wp: " << wp << endl;
    cout << runNb << ") Eff_b(MC)      : " <<  WParray[0] << " $\\pm$ " << WParray[1] << endl;
    cout << runNb << ") NR - Eff_b(FMC): " << WParray[2] << " $\\pm$ " << WParray[3] << "  " << WParray[3]/WParray[2] <<  " -- " << (WParray[2]-WParray[0])/WParray[0] <<  endl;
    cout << runNb << ") NR - Eff_b     : " << WParray[4] << " $\\pm$ " << WParray[5] << "  " << WParray[5]/WParray[4] <<  " -- " << (WParray[4]-WParray[0])/WParray[0] <<  endl;
    cout << runNb << ") RR - Eff_b(FMC): " << WParray[6] << " $\\pm$ " << WParray[7] << "  " << WParray[7]/WParray[6] <<  " -- " << (WParray[6]-WParray[0])/WParray[0] << endl;
    cout << runNb<< ") RR - Eff_b     : " << WParray[8] << " $\\pm$ " << WParray[9] << "  " << WParray[9]/WParray[8] <<  " -- " << (WParray[8]-WParray[0])/WParray[0] << endl;
  }

  //dirty saturday morning writing!
  //cout << "PtEtaBin::GetWPEffRR::start" << endl;
  /* for(int i=0; i<TH1Data_BtagEffMeasuredRR->GetXaxis()->GetNbins(); i++){
    if(fabs(wp-TH1Data_BtagEffMeasuredRR->GetXaxis()->GetBinUpEdge(i))<0.0001){
      //cout << "PtEtaBin::GetWPEffRR " << *eff << endl;
      *eff = TH1Data_BtagEffMeasuredRR->GetBinContent(i);//Is this correct???? shoudl it be upbin or lowbin???
      //cout << "PtEtaBin::GetWPEffRR " << *eff << endl;
      //err=0; //to be implemeted
      *MCeff = TH1Data_BtagEffMCMeasuredRR->GetBinContent(i);
      *FRatio = leftNonB_/rightNonB_;
      *MCFRatio = leftNonBMC_/rightNonBMC_;
      *FitPar0 = TFData_LeftRight1DXFit->GetParameter(0);
      *FitPar1 = TFData_LeftRight1DXFit->GetParameter(1);
    }
    }*/
}

void::PtEtaBin::DiffCalculation(TH1D* Estimate, TH1D* Orig, TH1D* result){
  //Orig stands for the MC plot
  //Estimate stands for the plot containing the btag efficiency estimator

  result->Add(Orig,Estimate,1,-1);
  result->Divide(Orig);

}

void PtEtaBin::EffCalculation(bool doSCreweigh, TH1D * THLeft, TH1D * THRight, TH1D* THBtag, TH1D *THEff, TH1D *THBtagMC, TH1D *THEffMC){
	
	if(debug_>2) cout << "+----> PtEtaBin::EffCalculation - THEff->GetXaxis()->GetNbins() " << THEff->GetXaxis()->GetNbins() << endl;  
	
	EffCalcDetails_.clear();
	
	double Fused=0;
	double Fusede=0;
	
	if(doSCreweigh){
		Fused=FReweigh_;
		Fusede=FReweighe_;
	} else {
		Fused=F_;
		Fusede=Fe_;
    }
	
	EffCalcDetails_["F"] = Fused;
	EffCalcDetails_["Fe"] = Fusede;
	EffCalcDetails_["Fmc"] = FMC_;
	EffCalcDetails_["Fmce"] = FMCe_;
	
	//HACK FOR PSEUDOEXP
	/*  if(doSCreweigh){
	 //    Fused=FReweigh_;
	 //    Fusede=FReweighe_;
	 Fused=F_;
	 Fusede=Fe_;
	 } else {
	 Fused=F_;
	 Fusede=Fe_;
	 }
	 */
	
	//FOR pseudo-experiments, keep FMC fixed to the overall F value (pseudoexp)
	//FMC_=0.959747;
	//FMCe_=0;
	
	//cout << *genericName << endl;
	
	//EffCalcDetails_["NLeft"]=0;
	for(int i=0; i<=THLeft->GetXaxis()->GetNbins()+1; i++)
		EffCalcDetails_["NLeft"] += THLeft->GetBinContent(i);
	
	//cout << "***********************************************************************************************************" << EffCalcDetails_["NLeft"]<<endl;

	//EffCalcDetails_["NRight"]=0;
	for(int i=0; i<=THRight->GetXaxis()->GetNbins()+1; i++)
		EffCalcDetails_["NRight"] += THRight->GetBinContent(i);
	
	for(int i=0; i<=THEff->GetXaxis()->GetNbins()+1; i++){
		
		THBtag->SetBinContent(i,THLeft->GetBinContent(i)-Fused*THRight->GetBinContent(i));
		THBtag->SetBinError(i,
							sqrt(
								 THLeft->GetBinContent(i)+
								 Fusede*Fusede*THRight->GetBinContent(i)*THRight->GetBinContent(i)+
								 Fused*Fused*THRight->GetBinContent(i)
								 ));
		
		THBtagMC->SetBinContent(i,THLeft->GetBinContent(i)-FMC_*THRight->GetBinContent(i));
		THBtagMC->SetBinError(i,sqrt(
									 THLeft->GetBinContent(i)+
									 FMCe_*FMCe_*THRight->GetBinContent(i)*THRight->GetBinContent(i)+
									 FMC_*FMC_*THRight->GetBinContent(i)
									 ));  
		
		
		 //  cout << THLeft->GetBinContent(i)-Fused*THRight->GetBinContent(i) << " " << THLeft->GetBinContent(i) << " " << Fused<< " " << THRight->GetBinContent(i) << endl;
		/*cout << "A) " << sqrt(
		 THLeft->GetBinContent(i)+
		 FMCe_*FMCe_*THRight->GetBinContent(i)*THRight->GetBinContent(i)+
		 FMC_*FMC_*THRight->GetBinContent(i)
		 ) << "    " <<
		 THLeft->GetBinContent(i) << "    " <<
		 FMCe_*FMCe_*THRight->GetBinContent(i)*THRight->GetBinContent(i) << "    " <<
		 FMC_*FMC_*THRight->GetBinContent(i) << "    "
		 << endl;*/
		/*cout << "B) " << sqrt(
		 THLeft->GetBinContent(i)+
		 Fusede*Fusede*THRight->GetBinContent(i)*THRight->GetBinContent(i)+
		 Fused*Fused*THRight->GetBinContent(i)
		 ) << "    " <<
		 THLeft->GetBinContent(i) << "    " <<
		 Fusede*Fusede*THRight->GetBinContent(i)*THRight->GetBinContent(i) << "    " <<
		 Fused*Fused*THRight->GetBinContent(i) << "    "
		 << endl;*/
		
	}
	
	//cout << "for efficiencies" << endl;
	
	double minMC=0;
	double minMCe=0;
	double plusMC=0;
	double plusMCe=0;
	double min=0;
	double mine=0;
	double plus=0;
	double pluse=0;
	for(int i=0; i<=THBtagMC->GetNbinsX()+1; i++){
		plusMC+=THBtagMC->GetBinContent(i);
		plusMCe+=pow(THBtagMC->GetBinError(i),2);
		plus+=THBtag->GetBinContent(i);
		pluse+=pow(THBtag->GetBinError(i),2);
		
		//cout << "Btag Bin " << i << " bin low " << THBtag->GetBinLowEdge(i) << " Content " << THBtag->GetBinContent(i) <<  " error: " << THBtag->GetBinError(i) << " ratio : " << THBtag->GetBinError(i)/THBtag->GetBinContent(i) << endl;
 	}
	
	THEff->SetBinContent(0,1);
	THEffMC->SetBinContent(0,1);
	
	for(int i=0; i<=THBtagMC->GetNbinsX(); i++){
		plusMC-=THBtagMC->GetBinContent(i);
		plusMCe-=pow(THBtagMC->GetBinError(i),2);
		plus-=THBtag->GetBinContent(i);
		pluse-=pow(THBtag->GetBinError(i),2);
		if(plusMC<1){plusMC=0; plusMCe=0;}
		if(plus<1){plus=0; pluse=0;}
		minMC+=THBtagMC->GetBinContent(i);
		minMCe+=pow(THBtagMC->GetBinError(i),2);
		min+=THBtag->GetBinContent(i);
		mine+=pow(THBtag->GetBinError(i),2);
		
		double tmpcontentMC=plusMC/(plusMC+minMC);
		
		double tmperrorMC=sqrt(
							   pow(sqrt(plusMCe) * ( 1/(plusMC+minMC) - plusMC/( (plusMC+minMC)*(plusMC+minMC) )  ) ,2) + 
							   pow(sqrt(minMCe)  *                  plusMC/( (plusMC+minMC)*(plusMC+minMC) )    ,2)
							   );
		double tmpcontent=plus/(plus+min);
		
		double tmperror=sqrt(
							 pow(sqrt(pluse) * ( 1/(plus+min) - plus/( (plus+min)*(plus+min) )  ) ,2) + 
							 pow(sqrt(mine)  *                  plus/( (plus+min)*(plus+min) )    ,2)
							 );
		
		/*cout << sqrt(
		 pow(sqrt(plusMCe) * ( 1/(plusMC+minMC) - plusMC/( (plusMC+minMC)*(plusMC+minMC) )  ) ,2) + 
		 pow(sqrt(minMCe)  *                  plusMC/( (plusMC+minMC)*(plusMC+minMC) )    ,2)
		 ) << "   " << 
		 plusMCe << "   " << 
		 minMCe << "   "
		 << endl;
		 cout << sqrt(
		 pow(sqrt(pluse) * ( 1/(plus+min) - plus/( (plus+min)*(plus+min) )  ) ,2) + 
		 pow(sqrt(mine)  *                  plus/( (plus+min)*(plus+min) )    ,2)
		 )<< "   " << 
		 pluse << "   " << 
		 mine << "   " 
		 <<endl;
		 */
		
		//cout << "Eff in bin " << i << " = " << tmpcontent << "+-" << tmperror << endl;
		
		THEff->SetBinContent(i+1,tmpcontent);
		THEff->SetBinError(i+1,tmperror);
		THEffMC->SetBinContent(i+1,tmpcontentMC);
		THEffMC->SetBinError(i+1,tmperrorMC);
		
		stringstream a; a << THBtag->GetBinLowEdge(i);
		
		string WP = "WP_"+a.str();
		
		EffCalcDetails_[WP+"_plus"]=plus;
		EffCalcDetails_[WP+"_plusMC"]=plusMC;
		EffCalcDetails_[WP+"_min"]=min;
		EffCalcDetails_[WP+"_minMC"]=minMC;

		if (THEff->GetBinLowEdge(i+1) == 1.7 || THEff->GetBinLowEdge(i+1) == 3.3 || THEff->GetBinLowEdge(i+1) == 10.2) {
		
			cout << "++> WP: " << THEff->GetBinLowEdge(i+1) << endl;
			cout << "min: " << min << " +- " << sqrt(mine) << " ratio(abs) " << sqrt(mine)/min << endl;
			cout << "plus " << plus << " +- " << sqrt(pluse) << " ratio " << sqrt(pluse)/plus << endl;
			cout << "Eff " << tmpcontent << " +- " << tmperror << endl;
			cout << endl;
		}
		
	}
	
	//old, with bug :)
	
	/*
	 double bEff=0;
	 double bEffMC=0;
	 double ebEff=0;
	 double ebEffMC=0;
	 double plusbEff=0;
	 double plusbEffMC=0;
	 double eplusbEff=0;
	 double eplusbEffMC=0;
	 double minbEff=0;
	 double minbEffMC=0;
	 double eminbEff=0;
	 double eminbEffMC=0;
	 
	 
	 for(int i=0; i<=THEff->GetXaxis()->GetNbins()+1; i++){
	 
	 
	 
	 double plus=THBtagMC->Integral(i,THBtagMC->GetNbinsX()+1);
	 double min=0;
	 int iprime=i-1;
	 if(iprime>=0){
	 //if(i2>0){
	 min=THBtagMC->Integral(0,iprime);
	 }
	 // double Eff=plus/(plus+min);
	 
	 //cout << "A) " << plus/(plus+min) << " = " << plus << " " << (plus+min) << endl;
	 //cout << THEffMC->GetBinContent(i) << endl;
	 
	 double plusRegionIntegralLeft=THLeft->Integral(i,THLeft->GetNbinsX()+1);
	 double plusRegionIntegralRight=THRight->Integral(i,THRight->GetNbinsX()+1); 
	 double eplusRegionIntegralLeft=sqrt(THLeft->Integral(i,THLeft->GetNbinsX()+1));
	 double eplusRegionIntegralRight=sqrt(THRight->Integral(i,THRight->GetNbinsX()+1));
	 
	 double minRegionIntegralLeft=0;
	 double minRegionIntegralRight=0;
	 double eminRegionIntegralLeft=0;
	 double eminRegionIntegralRight=0;
	 
	 int i2=i-1;
	 //if(i2<0) i2=0;
	 if(i2>=0){
	 minRegionIntegralLeft=THLeft->Integral(0,i2);
	 minRegionIntegralRight=THRight->Integral(0,i2);
	 eminRegionIntegralLeft=sqrt(THLeft->Integral(0,i2));
	 eminRegionIntegralRight=sqrt(THRight->Integral(0,i2));
	 }
	 
	 //cout << plusRegionIntegralLeft << " + " << minRegionIntegralLeft << " = " << plusRegionIntegralLeft+minRegionIntegralLeft << " ?= " << THLeft->Integral(0,THLeft->GetNbinsX()+1) << endl;    
	 plusbEff = plusRegionIntegralLeft-Fused*plusRegionIntegralRight;
	 plusbEffMC = plusRegionIntegralLeft-FMC_*plusRegionIntegralRight;
	 minbEff = minRegionIntegralLeft-Fused*minRegionIntegralRight;
	 minbEffMC = minRegionIntegralLeft-FMC_*minRegionIntegralRight;
	 
	 bEff=plusbEff/(plusbEff+minbEff);
	 bEffMC=plusbEffMC/(plusbEffMC+minbEffMC);
	 
	 //cout << "B) " << plusbEffMC/(plusbEffMC+minbEffMC) << " = " << plusbEffMC << " " << (plusbEffMC+minbEffMC) << endl;
	 
	 eplusbEff = sqrt(eplusRegionIntegralLeft*eplusRegionIntegralLeft + Fusede*plusRegionIntegralRight*Fusede*plusRegionIntegralRight + Fused*eplusRegionIntegralRight*Fused*eplusRegionIntegralRight);
	 
	 //cout << eplusbEff << " ?= " << sqrt(plusbEff) << endl;
	 
	 eplusbEffMC = sqrt(eplusRegionIntegralLeft*eplusRegionIntegralLeft + FMCe_*plusRegionIntegralRight*FMCe_*plusRegionIntegralRight + FMC_*eplusRegionIntegralRight*FMC_*eplusRegionIntegralRight);
	 
	 eminbEff = sqrt(eminRegionIntegralLeft*eminRegionIntegralLeft + Fusede*minRegionIntegralRight*Fusede*minRegionIntegralRight + Fused*eminRegionIntegralRight*Fused*eminRegionIntegralRight);
	 
	 eminbEffMC = sqrt(eminRegionIntegralLeft*eminRegionIntegralLeft + FMCe_*minRegionIntegralRight*FMCe_*minRegionIntegralRight + FMC_*eminRegionIntegralRight*FMC_*eminRegionIntegralRight);
	 
	 
	 //    ebEff=sqrt(eplusbEff/totalbEff * eplusbEff/totalbEff + plusbEff/(totalbEff*totalbEff)*etotalbEff * plusbEff/(totalbEff*totalbEff)*etotalbEff);
	 ebEff=sqrt((pow((minbEff*eplusbEff),2)+pow((plusbEff*eminbEff),2))/pow((plusbEff+minbEff),4));
	 
	 //    ebEffMC=sqrt(eplusbEffMC/totalbEffMC * eplusbEffMC/totalbEffMC + plusbEffMC/(totalbEffMC*totalbEffMC)*etotalbEff * plusbEffMC/(totalbEffMC*totalbEffMC)*etotalbEffMC);
	 ebEffMC=sqrt((pow((minbEffMC*eplusbEffMC),2)+pow((plusbEffMC*eminbEffMC),2))/pow((plusbEffMC+minbEffMC),4));
	 
	 THEff->SetBinContent(i,bEff);
	 THEff->SetBinError(i,ebEff);
	 THEffMC->SetBinContent(i,bEffMC);
	 THEffMC->SetBinError(i,ebEffMC);
	 
	 
	 cout << i << "   " << bEffMC << "    " << ebEffMC << endl;    
	 
	 }*/
	
	
	
	//old
	/*  double plusbEff=0;
	 double plusbEffMC=0;
	 double totalbEff=0;
	 double totalbEffMC=0;
	 double bEff=0;
	 double bEffMC=0;
	 
	 double eplusbEff=0;
	 double eplusbEffMC=0;
	 double etotalbEff=0;
	 double etotalbEffMC=0;
	 double ebEff=0;
	 double ebEffMC=0;
	 
	 double totalRegionIntegralLeft=THLeft->Integral(0,THLeft->GetNbinsX()+1);
	 double totalRegionIntegralRight=THRight->Integral(0,THRight->GetNbinsX()+1); 
	 double etotalRegionIntegralLeft=sqrt(THLeft->Integral(0,THLeft->GetNbinsX()+1));
	 double etotalRegionIntegralRight=sqrt(THRight->Integral(0,THRight->GetNbinsX()+1));
	 
	 if(debug_>1) cout << "+--> From Signal sample totalRegionIntegralLeft " << totalRegionIntegralLeft << " totalRegionIntegralRight " << totalRegionIntegralRight << endl; 
	 
	 double Fused=0;
	 double Fusede=0;
	 
	 if(doSCreweigh){
	 Fused=FReweigh_;
	 Fusede=FReweighe_;
	 } else {
	 Fused=F_;
	 Fusede=Fe_;
	 }
	 
	 //totalbEff = totalRegionIntegralLeft-(leftNonB_/rightNonB_)*totalRegionIntegralRight;
	 //totalbEffMC = totalRegionIntegralLeft-(leftNonBMC_/rightNonBMC_)*totalRegionIntegralRight;
	 
	 totalbEff = totalRegionIntegralLeft-Fused*totalRegionIntegralRight;
	 totalbEffMC = totalRegionIntegralLeft-FMC_*totalRegionIntegralRight;
	 
	 etotalbEff = sqrt(etotalRegionIntegralLeft*etotalRegionIntegralLeft + Fusede*totalRegionIntegralRight*Fusede*totalRegionIntegralRight + Fused*etotalRegionIntegralRight*Fused*etotalRegionIntegralRight);
	 etotalbEffMC = sqrt(etotalRegionIntegralLeft*etotalRegionIntegralLeft + FMCe_*totalRegionIntegralRight*FMCe_*totalRegionIntegralRight + FMC_*etotalRegionIntegralRight*FMC_*etotalRegionIntegralRight);
	 
	 //etotalbEff = sqrt(etotalRegionIntegralLeft*etotalRegionIntegralLeft + ((leftNonB_/rightNonB_)*etotalRegionIntegralRight*(leftNonB_/rightNonB_)*etotalRegionIntegralRight) + ((sqrt(leftNonB_)/rightNonB_)*totalRegionIntegralRight*(sqrt(leftNonB_)/rightNonB_)*totalRegionIntegralRight) + ((leftNonB_/(rightNonB_*rightNonB_))*sqrt(rightNonB_)*totalRegionIntegralRight*(leftNonB_/(rightNonB_*rightNonB_))*sqrt(rightNonB_)*totalRegionIntegralRight));
	 
	 if(debug_>2) cout << "+----> PtEtaBin::EffCalculation - THEff->GetXaxis()->GetNbins() " << THEff->GetXaxis()->GetNbins() << endl;
	 
	 for(int i=0; i<=THEff->GetXaxis()->GetNbins()+1; i++){
	 //THBtag->SetBinContent(i,THLeft->GetBinContent(i)-(leftNonB_/rightNonB_)*THRight->GetBinContent(i));
	 //THBtag->SetBinError(i,sqrt(THLeft->GetBinContent(i)-(leftNonB_/rightNonB_)*THRight->GetBinContent(i)));
     
	 THBtag->SetBinContent(i,THLeft->GetBinContent(i)-Fused*THRight->GetBinContent(i));
	 THBtag->SetBinError(i,sqrt(THLeft->GetBinContent(i)+Fusede*Fusede*THRight->GetBinContent(i)*THRight->GetBinContent(i)+Fused*Fused*THRight->GetBinContent(i)));
	 
	 //todo: double check; the sum of the errors here should be the same as etotalbEff
	 
	 THBtagMC->SetBinContent(i,THLeft->GetBinContent(i)-FMC_*THRight->GetBinContent(i));
	 THBtag->SetBinError(i,sqrt(THLeft->GetBinContent(i)+FMCe_*FMCe_*THRight->GetBinContent(i)*THRight->GetBinContent(i)+FMC_*FMC_*THRight->GetBinContent(i)));  
	 
	 double plusRegionIntegralLeft=THLeft->Integral(i,THLeft->GetNbinsX()+1);
	 double plusRegionIntegralRight=THRight->Integral(i,THRight->GetNbinsX()+1); 
	 double eplusRegionIntegralLeft=sqrt(THLeft->Integral(i,THLeft->GetNbinsX()+1));
	 double eplusRegionIntegralRight=sqrt(THRight->Integral(i,THRight->GetNbinsX()+1));
	 
	 double minRegionIntegralLeft=THLeft->Integral(0,i);
	 edit -> added this, but first clarify what is the appropriate way to handle the errors
	 double minRegionIntegralRight=THRight->Integral(0,i);
	 double eminRegionIntegralLeft=sqrt(THLeft->Integral(0,i));
	 double eminRegionIntegralRight=sqrt(THRight->Integral(0,i));
	 
	 plusbEff = plusRegionIntegralLeft-Fused*plusRegionIntegralRight;
	 plusbEffMC = plusRegionIntegralLeft-FMC_*plusRegionIntegralRight;
	 
	 bEff=plusbEff/totalbEff;
	 bEffMC=plusbEffMC/totalbEffMC;
	 
	 eplusbEff = sqrt(eplusRegionIntegralLeft*eplusRegionIntegralLeft + Fusede*plusRegionIntegralRight*Fusede*plusRegionIntegralRight + Fused*eplusRegionIntegralRight*Fused*eplusRegionIntegralRight);
	 
	 //eplusbEff = sqrt(eplusRegionIntegralLeft*eplusRegionIntegralLeft + ((leftNonB_/rightNonB_)*eplusRegionIntegralRight*(leftNonB_/rightNonB_)*eplusRegionIntegralRight) + ((sqrt(leftNonB_)/rightNonB_)*plusRegionIntegralRight*(sqrt(leftNonB_)/rightNonB_)*plusRegionIntegralRight) + ((leftNonB_/(rightNonB_*rightNonB_))*sqrt(rightNonB_)*plusRegionIntegralRight*(leftNonB_/(rightNonB_*rightNonB_))*sqrt(rightNonB_)*plusRegionIntegralRight));
	 ebEff=sqrt(eplusbEff/totalbEff * eplusbEff/totalbEff + plusbEff/(totalbEff*totalbEff)*etotalbEff * plusbEff/(totalbEff*totalbEff)*etotalbEff);
	 
	 eplusbEffMC = sqrt(eplusRegionIntegralLeft*eplusRegionIntegralLeft + FMCe_*plusRegionIntegralRight*FMCe_*plusRegionIntegralRight + FMC_*eplusRegionIntegralRight*FMC_*eplusRegionIntegralRight);
	 ebEffMC=sqrt(eplusbEffMC/totalbEffMC * eplusbEffMC/totalbEffMC + plusbEffMC/(totalbEffMC*totalbEffMC)*etotalbEff * plusbEffMC/(totalbEffMC*totalbEffMC)*etotalbEffMC);
	 
	 
	 THEff->SetBinContent(i,bEff);
	 THEff->SetBinError(i,ebEff);
	 THEffMC->SetBinContent(i,bEffMC);
	 THEffMC->SetBinError(i,ebEffMC);
	 
	 }*/
	//end of old
	
	/* for(int i=0; i<=THEff->GetXaxis()->GetNbins()+1; i++){
	 THBtag->SetBinContent(i,THLeft->GetBinContent(i)-(leftNonB_/rightNonB_)*THRight->GetBinContent(i));
	 THBtag->SetBinError(i,sqrt(THLeft->GetBinContent(i)-(leftNonB_/rightNonB_)*THRight->GetBinContent(i)));
	 THBtagMC->SetBinContent(i,THLeft->GetBinContent(i)-(leftNonBMC_/rightNonBMC_)*THRight->GetBinContent(i));
	 
	 if(THBtag->GetBinContent(i)<0) THBtag->SetBinContent(i,0);
	 if(THBtagMC->GetBinContent(i)<0) THBtagMC->SetBinContent(i,0);
	 
	 }
	 
	 for(int i=0; i<=THEff->GetXaxis()->GetNbins()+1; i++){
	 
	 double totalbtag=0;
	 double etotalbtag=0;
	 double plusbtag=0;
	 double eplusbtag=0;
	 double eEff=0;
	 totalbtag = THBtag->Integral(0,THBtag->GetNbinsX()+1);
	 etotalbtag = sqrt(THBtag->Integral(0,THBtag->GetNbinsX()+1));
	 plusbtag = THBtag->Integral(i,THBtag->GetNbinsX()+1);
	 eplusbtag = sqrt(THBtag->Integral(i,THBtag->GetNbinsX()+1));
	 eEff=sqrt(pow(eplusbtag/totalbtag,2)+pow(plusbtag*etotalbtag/(totalbtag*totalbtag),2)); //this error is not correct since the number of contents in each bin is the difference of two other values which on themselves have a poisonian error.
	 
	 double totalbtagMC=0;
	 double plusbtagMC=0;
	 totalbtagMC = THBtagMC->Integral(0,THBtag->GetNbinsX()+1);
	 plusbtagMC = THBtagMC->Integral(i,THBtag->GetNbinsX()+1);
	 
	 THEff->SetBinContent(i,plusbtag/totalbtag);
	 //THEff->SetBinError(i,eEff);
	 
	 THEffMC->SetBinContent(i,plusbtagMC/totalbtagMC);
	 
	 //if(i==0) cout << "THEff->GetBinContent(0) " << THEff->GetBinContent(i) << endl;
	 //if(i==0) cout << "THEffMC->GetBinContent(0) " << THEffMC->GetBinContent(i) << endl;
	 
	 }*/
	
	
}

void PtEtaBin::WriteSignalSamplePlots(){
  if(debug_>4) cout << "+--> PtEtaBin::WriteSignalSamplePlots" << endl;
  //if(debug_>1) cout << "Signal sample mlj, mean:    " << TH1Bkg_Var0->GetMean() << endl;
  //if(debug_>1) cout << "Signal sample mlj, spread:  " << TH1Bkg_Var0->GetRMS() << endl;

  TH2Sng_Left->Write();
  TH2Bkg_Left->Write();
  TH2Data_Left->Write();

  TH2Sng->Write();
  TH2Bkg->Write();
  TH2Data->Write();

  //cout << "correlation between mlj and btag SNG: " << TH2Sng->GetCorrelationFactor() << endl;
  //cout << "correlation between mlj and btag BKG: " << TH2Bkg->GetCorrelationFactor() << endl;
  //cout << "correlation between mlj and btag DATA: " << TH2Data->GetCorrelationFactor() << endl;

  TH2Sng_Right->Write();
  TH2Bkg_Right->Write();
  TH2Data_Right->Write();

  TH2SngVar12_Left->Write();
  TH2BkgVar12_Left->Write();
  TH2DataVar12_Left->Write();
  
  TH2SngVar12_Right->Write();
  TH2BkgVar12_Right->Write();
  TH2DataVar12_Right->Write();

  TH1Sng_chisq->Write();
  TH1Bkg_chisq->Write();
  TH1Data_chisq->Write();
  TH1Sng_chisqControl->Write();
  TH1Bkg_chisqControl->Write();
  TH1Data_chisqControl->Write();

  TH1Sng_Var0->Write();
  TH1Bkg_Var0->Write();
  TH1Bkg_W_Var0->Write();
  TH1Bkg_R_Var0->Write();
    TH1Data_Var0->Write();
    TH1Data_Var0_XS->Write();
	
	for (std::map<TString,TH1D*>::const_iterator it=histo1D.begin(); it != histo1D.end(); ++it)
		
		it->second->Write();
	
	

  TH2Sng_Var12->Write();
  TH2Bkg_Var12->Write();
  TH2Bkg_W_Var12->Write();
  TH2Bkg_R_Var12->Write();
  TH2Data_Var12->Write();

  TH1Sng_BtagAll->Write();
  TH1Bkg_BtagAll->Write();
  TH1Data_BtagAll->Write();
	
	TH1Data_Pt_Left->Write();
	TH1Sng_Pt_Left->Write();
    TH1Data_Pt_Right->Write();
	TH1Data_Pt_Control->Write();
}

void PtEtaBin::WriteSoverSBPlot(){
  if(debug_>4) cout << "+--> PtEtaBin::WriteSoverSBPlot" << endl;
  TH1SoverSB_Var0->Write();  
  TH1SoverSB_ControlVar->Write();

}


void PtEtaBin::WriteControlSamplePlots(){
  if(debug_>4) cout << "+--> PtEtaBin::WriteControlSamplePlots" << endl;
  //if(debug_>1) cout << "Signal sample mlj, mean:    " << TH1Data_ControlVarReweigh->GetMean() << endl;
  //if(debug_>1) cout << "Signal sample mlj, spread:  " << TH1Data_ControlVarReweigh->GetRMS() << endl;

  TH1Sng_ControlVar->Write();
  TH1Bkg_ControlVar->Write();
  TH1Bkg_W_ControlVar->Write();
  TH1Bkg_R_ControlVar->Write();
  TH1Data_ControlVar->Write(); 
 
  TH1Sng_ControlVarReweigh->Write();
  TH1Bkg_ControlVarReweigh->Write();
  TH1Bkg_W_ControlVarReweigh->Write();
  TH1Bkg_R_ControlVarReweigh->Write();
  TH1Data_ControlVarReweigh->Write(); 
  
  TH2Sng_ControlVar12->Write();
  TH2Bkg_ControlVar12->Write();
  TH2Bkg_W_ControlVar12->Write();
  TH2Bkg_R_ControlVar12->Write();
  TH2Data_ControlVar12->Write();

}

void PtEtaBin::WriteMCEffPlots(bool doCout){
  if(debug_>4) cout << "+--> PtEtaBin::MakeMCEffPlots" << endl;
  TH1Sng_BtagEffMC->Write();
  TH1Sng_BtagEffAll->Write();
  //TH1Bkg_BtagEff->Write();
  //TH1Data_BtagEff->Write();
  if(doCout){
    cout << "################## define working points ##################" << endl;
    
    bool low=false;
    bool medium=false;
    bool tight=false;
    for(int i=0; i<TH1Sng_BtagEffAll->GetNbinsX(); i++){
      //    cout << TH1Sng_BtagEffAll->GetBinContent(i) << " " << TH1Sng_BtagEffAll->GetBinLowEdge(i) << endl;
      if(!low) if(TH1Sng_BtagEffAll->GetBinContent(i)<0.75) {cout << TH1Sng_BtagEffAll->GetBinLowEdge(i) << endl; low=true;}
      //if(TH1Sng_BtagEffAll->GetBinContent(i)<0.75) cout << TH1Sng_BtagEffAll->GetBinLowEdge(i) << endl; low=true;
      if(!medium) if(TH1Sng_BtagEffAll->GetBinContent(i)<0.5) {cout << TH1Sng_BtagEffAll->GetBinLowEdge(i) << endl; medium=true;}
      if(!tight) if(TH1Sng_BtagEffAll->GetBinContent(i)<0.25) {cout << TH1Sng_BtagEffAll->GetBinLowEdge(i) << endl; tight=true;}
      
    }  
  }
  
}

void PtEtaBin::WriteProfileXplots(){
  if(debug_>4) cout << "+--> PtEtaBin::WriteProfileXplots" << endl;
  TH2Sng_ProfileX->Write();
  TH2Bkg_ProfileX->Write();
  TH2Data_ProfileX->Write(); 
}

void PtEtaBin::Write1DXplots(bool doCout){
  if(debug_>4) cout << "+--> PtEtaBin::Write1DXplots" << endl;
  TH2Sng_Left1DX->Write();

  doCout=false;
  if(doCout){

    double integral = TH2Sng_Left1DX->Integral(0,TH2Sng_Left1DX->GetXaxis()->GetNbins()+1);

    for(int i=0; i<TH2Sng_Left1DX->GetXaxis()->GetNbins()+1; i++){
      cout << "bin:  [ " << TH2Sng_Left1DX->GetXaxis()->GetBinLowEdge(i) << " - " << TH2Sng_Left1DX->GetXaxis()->GetBinUpEdge(i) << " ] : " << TH2Sng_Left1DX->GetBinContent(i)/integral << endl;
    }
 

    cout << "{" << endl;

    for(int i=0; i<TH2Sng_Left1DX->GetXaxis()->GetNbins()+1; i++){
      cout << TH2Sng_Left1DX->GetBinContent(i)/integral << ",";
    }
    cout << "}" << endl;

  }

  TH2Bkg_Left1DX->Write();
  TH2Data_Left1DX->Write(); 
  TH2Sng_Right1DX->Write();
  TH2Bkg_Right1DX->Write();
  TH2Data_Right1DX->Write();
}

void PtEtaBin::Write1DYplots(){
  if(debug_>4) cout << "+--> PtEtaBin::Write1DYplots" << endl;
   TH2Data_Left1DY->Write();
   TH2Sng_Left1DY->Write();
   TH2Bkg_Left1DY->Write();
   
   TH2Sng_Right1DY->Write();
   TH2Bkg_Right1DY->Write();
   TH2Data_Right1DY->Write();
}

void PtEtaBin::Write1DYVar2plots(){
  if(debug_>4) cout << "+--> PtEtaBin::Write1DYVar2plots" << endl;
  TH2SngVar12_Left1DY->Write();
  TH2BkgVar12_Left1DY->Write();
  TH2DataVar12_Left1DY->Write(); 
  TH2SngVar12_Right1DY->Write();
  TH2BkgVar12_Right1DY->Write();
  TH2DataVar12_Right1DY->Write();
}

void PtEtaBin::WriteEffPlots(){  
  if(debug_>4) cout << "+--> PtEtaBin::WriteEffPlots" << endl;
  TH1Data_BtagMeasured->Write();
  TH1Data_BtagMCMeasured->Write();
  TH1Data_BtagEffMeasured->Write();
  TH1Data_BtagEffMCMeasured->Write(); 
  TH1Data_BtagEffMeasuredDiff->Write();
  TH1Data_BtagEffMCMeasuredDiff->Write();
}

void PtEtaBin::WriteXRatioPlot(){
  if(debug_>4) cout << "+--> PtEtaBin::WriteXRatioPlot" << endl;
  TH2Data_RightLeft1DX->Write();  
  TH2Data_LeftRight1DX->Write();
	TFData_LeftRight1DXFit->Write();

}

void PtEtaBin::WriteSCPlot(){
  if(debug_>4) cout << "+--> PtEtaBin::WriteSCPlot" << endl;
  
  TH2Sng_Var12_1DX->Write();
  TH2Bkg_Var12_1DX->Write();
  TH2Bkg_W_Var12_1DX->Write();
  TH2Bkg_R_Var12_1DX->Write();
  TH2Data_Var12_1DX->Write();
  
  TH2Sng_Var12_1DY->Write();
  TH2Bkg_Var12_1DY->Write();
  TH2Bkg_W_Var12_1DY->Write();
  TH2Bkg_R_Var12_1DY->Write();
  TH2Data_Var12_1DY->Write();
  TH2Sng_ControlVar12_1DX->Write();
  TH2Bkg_ControlVar12_1DX->Write();
  TH2Bkg_W_ControlVar12_1DX->Write();
  TH2Bkg_R_ControlVar12_1DX->Write();
  TH2Data_ControlVar12_1DX->Write();
  TH2Sng_ControlVar12_1DY->Write();
  TH2Bkg_ControlVar12_1DY->Write();
  TH2Bkg_W_ControlVar12_1DY->Write();
  TH2Bkg_R_ControlVar12_1DY->Write();
  TH2Data_ControlVar12_1DY->Write();

  TH2Data_SignalControlVar12->Write();
  TH2Data_SignalControlVar12_1DX->Write();
  TH2Data_SignalControlVar12_1DY->Write();
  
}

void PtEtaBin::Write2DRatioPlot(){
  if(debug_>4) cout << "+--> PtEtaBin::Write2DRatioPlot" << endl;
  TH2DataVar12_LeftRight->Write();
}

void PtEtaBin::WriteEffPlotsLR(){
  if(debug_>4) cout << "+--> PtEtaBin::WriteEffPlotsLR" << endl;
  //  TH2Data_Left1DYReweigh->Write();
  TH1Data_BtagMeasuredLR->Write();
  TH1Data_BtagMCMeasuredLR->Write();
  TH1Data_BtagEffMeasuredLR->Write();
  TH1Data_BtagEffMCMeasuredLR->Write();	
  TH1Data_BtagEffMeasuredLRDiff->Write();
  TH1Data_BtagEffMCMeasuredLRDiff->Write();
}

void PtEtaBin::WriteEffPlotsRR(){ 
  if(debug_>4) cout << "+--> PtEtaBin::WriteEffPlotsRR" << endl; 

  TH2SngVar12_Right1DYReweigh->Write();
  TH2BkgVar12_Right1DYReweigh->Write();
  TH2DataVar12_Right1DYReweigh->Write(); //this is not really the proper place to do this! should write a new function for it
  
  TH2SngVar12_Right1DYReweighRatio->Write();
  TH2BkgVar12_Right1DYReweighRatio->Write();
  TH2DataVar12_Right1DYReweighRatio->Write(); //this is not really the proper place to do this! should write a new function for it

  TH2Data_Right1DYReweigh->Write(); 
  TH2Sng_Right1DYReweigh->Write(); 
  TH2Bkg_Right1DYReweigh->Write(); 
  TH2Data_Left1DYReweigh->Write(); 
  TH2Sng_Left1DYReweigh->Write(); 
  TH2Bkg_Left1DYReweigh->Write(); 

  TH2Data_Right1DY_LeftRightRatio->Write(); 
  TH2Sng_Right1DY_LeftRightRatio->Write(); 
  TH2Bkg_Right1DY_LeftRightRatio->Write(); 

  TH2Data_Right1DYReweigh_LeftRightRatio->Write(); 
  TH2Sng_Right1DYReweigh_LeftRightRatio->Write(); 
  TH2Bkg_Right1DYReweigh_LeftRightRatio->Write(); 
  TH2Data_Left1DYReweigh_LeftRightRatio->Write(); 
  TH2Sng_Left1DYReweigh_LeftRightRatio->Write(); 
  TH2Bkg_Left1DYReweigh_LeftRightRatio->Write(); 

  TH2Data_Right1DYReweighRatio->Write(); 
  TH2Sng_Right1DYReweighRatio->Write(); 
  TH2Bkg_Right1DYReweighRatio->Write(); 
  TH2Data_Left1DYReweighRatio->Write(); 
  TH2Sng_Left1DYReweighRatio->Write(); 
  TH2Bkg_Left1DYReweighRatio->Write(); 

  TH1Data_BtagMeasuredRR->Write();
  TH1Data_BtagMCMeasuredRR->Write();
  TH1Data_BtagEffMeasuredRR->Write();
  TH1Data_BtagEffMCMeasuredRR->Write();	
  TH1Data_BtagEffMeasuredRRDiff->Write();
  TH1Data_BtagEffMCMeasuredRRDiff->Write();
}

void PtEtaBin::WriteHistosToRootFile(bool w1, bool w2, bool w3, bool w4, bool w5, bool w6, bool w7, bool w8, bool w9, bool w10, bool w11, bool dow5cout){
  if(w1) WriteSignalSamplePlots();
  if(w2) WriteControlSamplePlots();
  if(w3) WriteSoverSBPlot();
  if(w4) WriteMCEffPlots(dow5cout);
  if(w5) Write1DXplots(dow5cout);
  if(w5) Write1DYVar2plots();
  if(w6) Write1DYplots();
  if(w7) WriteProfileXplots();
  if(w8) WriteEffPlots();
  if(w9) WriteXRatioPlot();
  if(w9) WriteSCPlot();
  if(w9) Write2DRatioPlot();
  if(w10) WriteEffPlotsLR();
  if(w11) WriteEffPlotsRR();
}

void PtEtaBin::WriteHistoToPSFile(TString *pathname, bool all){

  DrawEffNR(pathname);
  DrawbtagNR(pathname);
  DrawEffRR(pathname);
  DrawEffDiffNR(pathname);
  DrawEffDiffRR(pathname);
  if(all){
    DrawControlRWXPlot(pathname);
    DrawControlRWYPlot(pathname);
    DrawBasicFPlot(pathname);
    DrawControlFPlot(pathname);
    DrawControlRWFPlot(pathname);
    DrawControlOnlyRPlot(pathname);
    DrawControlOnlyWPlot(pathname);
    DrawSignalFPlot(pathname);
    DrawPtDistro(pathname);
    DrawSoverSB(pathname);
    DrawControlWSignalControlXPlot(pathname);
    DrawControlWSignalControlYPlot(pathname);
    DrawControlRSignalControlXPlot(pathname);
    DrawControlRSignalControlYPlot(pathname);
  }
}

void PtEtaBin::DrawSoverSB(TString *pathname){ //in fact I should change the name var1 with var1

  TString fname = *pathname;
  fname += "SoverSB";
  //fname +="_nDisc_"; fname+=nBdisc_; fname +="_pt_"; fname += (int) pt_; fname +="_eta_"; fname +=(int) round(eta_*10); 
  //  fname = GiveName(fname);
  GiveName(&fname);
  fname +=".eps";
  TCanvas c("c","",1);
  c.cd();
  c.SetGrid();

  TH1SoverSB_Var0->SetAxisRange(0,0.7,"Y");
  TH1SoverSB_Var0->SetTitle("");
  TH1SoverSB_Var0->GetXaxis()->SetTitle("m_{lj} (GeV)");
  TH1SoverSB_Var0->GetYaxis()->SetTitle("b flavor purity");
  TH1SoverSB_Var0->Draw("");
 
  //  TLegend leg(0.65,0.8,0.9,0.9);
  //leg.AddEntry(TH1SoverSB_Var0,"S over (S+B)","");

  //leg.SetFillColor(0); 
  //leg.Draw();

  c.Print(fname,"Landscape"); 

  TString fname2 = *pathname;
  fname2 += "SandB";
  GiveName(&fname2);
  fname2 +=".eps";
  TCanvas c2("c2","",1);
  c2.cd();
  c2.SetGrid();

  TH1Sng_Var0->SetTitle("");
  TH1Sng_Var0->GetXaxis()->SetTitle("m_{lj} (GeV)");
  TH1Sng_Var0->GetYaxis()->SetTitle("# jets/fb^{-1}");
  TH1Sng_Var0->SetLineWidth(2);
  TH1Sng_Var0->Draw("histo");
  TH1Bkg_Var0->SetTitle("");
  TH1Bkg_Var0->Draw("samehisto");
  TH1Bkg_Var0->SetLineStyle(2);
  TH1Bkg_Var0->SetLineWidth(2);

  TLegend leg2(0.65,0.8,0.9,0.9);
  leg2.AddEntry(TH1Sng_Var0,"b jets","l");
  leg2.AddEntry(TH1Bkg_Var0,"non-b jets","l");

  leg2.SetFillColor(0); 
  leg2.Draw();

  c2.Print(fname2,"Landscape"); 

}

void PtEtaBin::DrawPtDistro(TString *pathname){ //in fact I should change the name var1 with var1

  TString fname = *pathname;
  fname += "PtPlot";
  //fname +="_nDisc_"; fname+=nBdisc_; fname +="_pt_"; fname += (int) pt_; fname +="_eta_"; fname +=(int) round(eta_*10); 
  //  fname = GiveName(fname);
  GiveName(&fname);
  fname +=".eps";
  TCanvas c("c","",1);
  c.cd();

  double rightScale=TH2Data_Right1DX->Integral(0,TH2Data_Right1DX->GetXaxis()->GetNbins()+1);
  double leftScale=TH2Data_Left1DX->Integral(0,TH2Data_Left1DX->GetXaxis()->GetNbins()+1);
  
  TH2Data_Right1DX->Scale(1/rightScale);
  TH2Data_Left1DX->Scale(1/leftScale);

  TH2Data_Left1DX->SetTitle("");
  TH2Data_Left1DX->GetXaxis()->SetTitle("p_{T} (GeV)");
  TH2Data_Left1DX->GetYaxis()->SetTitle("a.u.");
  TH2Data_Left1DX->Draw("histo");
  TH2Data_Left1DX->SetLineColor(3);
  TH2Data_Right1DX->Draw("samehisto");
  TH2Data_Right1DX->SetLineColor(4);

  TLegend leg(0.65,0.8,0.9,0.9);
  leg.AddEntry(TH2Data_Left1DX,"P_{T,left}","l");
  leg.AddEntry(TH2Data_Right1DX,"P_{T,right}","l");
  leg.SetFillColor(0); 
  leg.Draw();

  c.Print(fname,"Landscape"); 

  TH2Data_Right1DX->Scale(rightScale);
  TH2Data_Left1DX->Scale(leftScale);


  TString fname2 = *pathname;

  fname2 += "PtRatioPlot";
  GiveName(&fname2);
  fname2 +=".eps";
  TCanvas c2("c2","",1);
  c2.cd();
  TH2Data_LeftRight1DX->SetTitle("");
  TH2Data_LeftRight1DX->GetXaxis()->SetTitle("p_{T} (GeV)");
  TH2Data_LeftRight1DX->GetYaxis()->SetTitle("a.u.");
  TH2Data_LeftRight1DX->Draw();

  TLegend leg2(0.65,0.8,0.9,0.9);
  leg2.AddEntry(TH2Data_LeftRight1DX,"P_{T,left} / P_{T,right}  ","p");
  leg2.SetFillColor(0); 
  leg2.Draw();

  c2.Print(fname2,"Landscape"); 
  
}

void PtEtaBin::DrawControlRWXPlot(TString *pathname){

  //These plots reflect the reweighting which will be done when using data. There is no distinction between b or non b, that's all jets folks!

  TString fname = *pathname;
  fname += "ControlRWXPlot_";
  GiveName(&fname);
  fname +=".eps";
  TCanvas c("c","",1);
  c.cd();
  //c.SetGrid();

  double scale=TH2Data_Var12_1DX->Integral(1,TH2Data_Var12_1DX->GetXaxis()->GetNbins()+1);
  double scaleControl=TH2Data_ControlVar12_1DX->Integral(1,TH2Data_ControlVar12_1DX->GetXaxis()->GetNbins()+1);

  TH2Data_Var12_1DX->Scale(1/scale);
  TH2Data_ControlVar12_1DX->Scale(1/scaleControl);

  TH2Data_ControlVar12_1DX->SetTitle("");
  TH2Data_ControlVar12_1DX->GetXaxis()->SetTitle("p_{T} (GeV)");
  TH2Data_ControlVar12_1DX->GetYaxis()->SetTitle("a.u.");
  TH2Data_ControlVar12_1DX->Draw("histoe");  
  //TH2Data_ControlVar12_1DX->Draw("PE");  
  TH2Data_ControlVar12_1DX->SetLineColor(2);
  //  TH2Data_Var12_1DX->SetLineColor(3);

  TH2Data_Var12_1DX->SetLineStyle(2);
  TH2Data_Var12_1DX->SetLineWidth(2);
  TH2Data_Var12_1DX->Draw("samehistoe");
  //TH2Data_Var12_1DX->Draw("samepe");
  //  TH2Data_ControlVar12_1DX->SetLineColor(4);
  //TH2Data_Var12_1DX->SetLineColor(2);
  
  //TH2Data_SignalControlVar12_1DX
  //TH2Data_Var12_1DX->GetBinContent(i)/TH2Data_ControlVar12_1DX->GetBinContent(i)

  TLegend leg(0.65,0.8,0.9,0.9);
  leg.AddEntry(TH2Data_Var12_1DX,"signal sample","l");
  leg.AddEntry(TH2Data_ControlVar12_1DX,"control sample","l");
  leg.SetFillColor(0); 
  leg.Draw();

  c.Print(fname,"Landscape");  
  TH2Data_Var12_1DX->Scale(scale);
  TH2Data_ControlVar12_1DX->Scale(scaleControl);
 
  TString fname2 = *pathname;
  fname2 += "ControlRWXPlot_Ratio_";
  GiveName(&fname2);
  fname2 +=".eps";
  TCanvas c2("c2","",1);
  c2.cd();
  //c2.SetGrid();

  TH2Data_SignalControlVar12_1DX->SetTitle("");
  TH2Data_SignalControlVar12_1DX->GetXaxis()->SetTitle("p_{T} (GeV)");
  TH2Data_SignalControlVar12_1DX->GetYaxis()->SetRangeUser(-4,30);
  
  TH2Data_SignalControlVar12_1DX->GetYaxis()->SetTitle("weight");
  TH2Data_SignalControlVar12_1DX->Draw("e");

  TLegend leg2(0.65,0.8,0.9,0.9);
  leg2.AddEntry(TH2Data_SignalControlVar12_1DX,"signal/control","l");
  leg2.SetFillColor(0); 
  leg2.Draw();

  c2.Print(fname2,"Landscape");  

}

void PtEtaBin::DrawControlRWYPlot(TString *pathname){

  TString fname = *pathname;
  fname += "ControlRWYPlot_";
  GiveName(&fname);
  fname +=".eps";
  TCanvas c("c","",1);
  c.cd();
  //c.SetGrid();

  double scale=TH2Data_Var12_1DY->Integral(1,TH2Data_Var12_1DY->GetXaxis()->GetNbins()+1);
  double scaleControl=TH2Data_ControlVar12_1DY->Integral(1,TH2Data_ControlVar12_1DY->GetXaxis()->GetNbins()+1);

  TH2Data_Var12_1DY->Scale(1/scale);
  TH2Data_ControlVar12_1DY->Scale(1/scaleControl);

  TH2Data_ControlVar12_1DY->SetTitle("");
  TH2Data_ControlVar12_1DY->GetXaxis()->SetTitle("#eta");
  TH2Data_ControlVar12_1DY->GetYaxis()->SetTitle("a.u.");
  TH2Data_ControlVar12_1DY->Draw("histoe"); 
  TH2Data_ControlVar12_1DY->SetLineColor(2);

  //  TH2Data_Var12_1DY->SetLineColor(3);
  TH2Data_Var12_1DY->SetLineStyle(2);
  TH2Data_Var12_1DY->SetLineWidth(2);
  TH2Data_Var12_1DY->Draw("samehistoe");
  //  TH2Data_ControlVar12_1DY->SetLineColor(4);
  //  TH2Data_Var12_1DY->SetLineColor(2);

  //TH2Data_SignalControlVar12_1DY
  //TH2Data_Var12_1DY->GetBinContent(i)/TH2Data_ControlVar12_1DY->GetBinContent(i)

  TLegend leg(0.65,0.8,0.9,0.9);
  leg.AddEntry(TH2Data_Var12_1DY,"signal sample","l");
  leg.AddEntry(TH2Data_ControlVar12_1DY,"control sample","l");
  leg.SetFillColor(0); 
  leg.Draw();

  c.Print(fname,"Landscape");  
  TH2Data_Var12_1DY->Scale(scale);
  TH2Data_ControlVar12_1DY->Scale(scaleControl);
 
  TString fname2 = *pathname;
  fname2 += "ControlRWYPlot_Ratio_";
  GiveName(&fname2);
  fname2 +=".eps";
  TCanvas c2("c2","",1);
  c2.cd();
  //c2.SetGrid();

  TH2Data_SignalControlVar12_1DY->SetTitle("");
  TH2Data_SignalControlVar12_1DY->GetXaxis()->SetTitle("#eta");
  //  TH2Data_SignalControlVar12_1DY->GetYaxis()->SetRangeUser(-4,30);
  
  TH2Data_SignalControlVar12_1DY->GetYaxis()->SetTitle("weight");
  TH2Data_SignalControlVar12_1DY->Draw();

  TLegend leg2(0.65,0.8,0.9,0.9);
  leg2.AddEntry(TH2Data_SignalControlVar12_1DY,"signal/control","l");
  leg2.SetFillColor(0); 
  leg2.Draw();

  c2.Print(fname2,"Landscape");  

}

void PtEtaBin::DrawControlWSignalControlXPlot(TString *pathname){

  TString fname = *pathname;
  fname += "ControlWSignalControlXPlot_";
  GiveName(&fname);
  fname +=".eps";
  TCanvas c("c","",1);
  c.cd();
  //c.SetGrid();

  double scale=TH2Bkg_W_Var12_1DX->Integral(1,TH2Bkg_W_Var12_1DX->GetXaxis()->GetNbins()+1);
  double scaleControl=TH2Bkg_W_ControlVar12_1DX->Integral(1,TH2Bkg_W_ControlVar12_1DX->GetXaxis()->GetNbins()+1);

  TH2Bkg_W_Var12_1DX->Scale(1/scale);
  TH2Bkg_W_ControlVar12_1DX->Scale(1/scaleControl);

  TH2Bkg_W_ControlVar12_1DX->SetTitle("");
  TH2Bkg_W_ControlVar12_1DX->GetXaxis()->SetTitle("p_{T} GeV");
  TH2Bkg_W_ControlVar12_1DX->GetYaxis()->SetTitle("a.u.");
  TH2Bkg_W_ControlVar12_1DX->Draw(""); 
  TH2Bkg_W_ControlVar12_1DX->SetMarkerColor(kRed);

  //TH2Bkg_W_Var12_1DX->SetMarkerColor(kRed);
  //  TH2Bkg_W_ControlVar12_1DX->SetLineWidth(2);
  TH2Bkg_W_Var12_1DX->Draw("histosame");
 
  TLegend leg(0.65,0.8,0.9,0.9);
  leg.AddEntry(TH2Bkg_W_Var12_1DX,"signal sample","l");
  leg.AddEntry(TH2Bkg_W_ControlVar12_1DX,"control sample","p");
  leg.SetFillColor(0); 
  leg.Draw();

  c.Print(fname,"Landscape");  

  TH2Bkg_W_Var12_1DX->Scale(scale);
  TH2Bkg_W_ControlVar12_1DX->Scale(scaleControl);
 
}

void PtEtaBin::DrawControlWSignalControlYPlot(TString *pathname){

  TString fname = *pathname;
  fname += "ControlWSignalControlYPlot_";
  GiveName(&fname);
  fname +=".eps";
  TCanvas c("c","",1);
  c.cd();
  //c.SetGrid();

  double scale=TH2Bkg_W_Var12_1DY->Integral(1,TH2Bkg_W_Var12_1DY->GetXaxis()->GetNbins()+1);
  double scaleControl=TH2Bkg_W_ControlVar12_1DY->Integral(1,TH2Bkg_W_ControlVar12_1DY->GetXaxis()->GetNbins()+1);

  TH2Bkg_W_Var12_1DY->Scale(1/scale);
  TH2Bkg_W_ControlVar12_1DY->Scale(1/scaleControl);

  TH2Bkg_W_ControlVar12_1DY->SetTitle("");
  TH2Bkg_W_ControlVar12_1DY->GetXaxis()->SetTitle("#eta");
  TH2Bkg_W_ControlVar12_1DY->GetYaxis()->SetTitle("a.u.");
  TH2Bkg_W_ControlVar12_1DY->Draw(""); 
  TH2Bkg_W_ControlVar12_1DY->SetMarkerColor(kRed);

  //TH2Bkg_W_Var12_1DY->SetMarkerColor(kRed);
  //  TH2Bkg_W_ControlVar12_1DY->SetLineWidth(2);
  TH2Bkg_W_Var12_1DY->Draw("histosame");
 
  TLegend leg(0.65,0.8,0.9,0.9);
  leg.AddEntry(TH2Bkg_W_Var12_1DY,"signal sample","l");
  leg.AddEntry(TH2Bkg_W_ControlVar12_1DY,"control sample","p");
  leg.SetFillColor(0); 
  leg.Draw();

  c.Print(fname,"Landscape");  

  TH2Bkg_W_Var12_1DY->Scale(scale);
  TH2Bkg_W_ControlVar12_1DY->Scale(scaleControl);
 
}

void PtEtaBin::DrawControlRSignalControlXPlot(TString *pathname){

  TString fname = *pathname;
  fname += "ControlRSignalControlXPlot_";
  GiveName(&fname);
  fname +=".eps";
  TCanvas c("c","",1);
  c.cd();
  //c.SetGrid();

  double scale=TH2Bkg_R_Var12_1DX->Integral(1,TH2Bkg_R_Var12_1DX->GetXaxis()->GetNbins()+1);
  double scaleControl=TH2Bkg_R_ControlVar12_1DX->Integral(1,TH2Bkg_R_ControlVar12_1DX->GetXaxis()->GetNbins()+1);

  TH2Bkg_R_Var12_1DX->Scale(1/scale);
  TH2Bkg_R_ControlVar12_1DX->Scale(1/scaleControl);

  TH2Bkg_R_ControlVar12_1DX->SetTitle("");
  TH2Bkg_R_ControlVar12_1DX->GetXaxis()->SetTitle("p_{T} GeV");
  TH2Bkg_R_ControlVar12_1DX->GetYaxis()->SetTitle("a.u.");
  TH2Bkg_R_ControlVar12_1DX->Draw(""); 
  TH2Bkg_R_ControlVar12_1DX->SetMarkerColor(kRed);

  //TH2Bkg_R_Var12_1DX->SetMarkerColor(kRed);
  //  TH2Bkg_R_ControlVar12_1DX->SetLineWidth(2);
  TH2Bkg_R_Var12_1DX->Draw("histosame");
 
  TLegend leg(0.65,0.8,0.9,0.9);
  leg.AddEntry(TH2Bkg_R_Var12_1DX,"signal sample","l");
  leg.AddEntry(TH2Bkg_R_ControlVar12_1DX,"control sample","p");
  leg.SetFillColor(0); 
  leg.Draw();

  c.Print(fname,"Landscape");  

  TH2Bkg_R_Var12_1DX->Scale(scale);
  TH2Bkg_R_ControlVar12_1DX->Scale(scaleControl);
 
}

void PtEtaBin::DrawControlRSignalControlYPlot(TString *pathname){

  TString fname = *pathname;
  fname += "ControlRSignalControlYPlot_";
  GiveName(&fname);
  fname +=".eps";
  TCanvas c("c","",1);
  c.cd();
  //c.SetGrid();

  double scale=TH2Bkg_R_Var12_1DY->Integral(1,TH2Bkg_R_Var12_1DY->GetXaxis()->GetNbins()+1);
  double scaleControl=TH2Bkg_R_ControlVar12_1DY->Integral(1,TH2Bkg_R_ControlVar12_1DY->GetXaxis()->GetNbins()+1);

  TH2Bkg_R_Var12_1DY->Scale(1/scale);
  TH2Bkg_R_ControlVar12_1DY->Scale(1/scaleControl);

  TH2Bkg_R_ControlVar12_1DY->SetTitle("");
  TH2Bkg_R_ControlVar12_1DY->GetXaxis()->SetTitle("#eta");
  TH2Bkg_R_ControlVar12_1DY->GetYaxis()->SetTitle("a.u.");
  TH2Bkg_R_ControlVar12_1DY->Draw(""); 
  TH2Bkg_R_ControlVar12_1DY->SetMarkerColor(kRed);

  //TH2Bkg_R_Var12_1DY->SetMarkerColor(kRed);
  //  TH2Bkg_R_ControlVar12_1DY->SetLineWidth(2);
  TH2Bkg_R_Var12_1DY->Draw("histosame");
 
  TLegend leg(0.65,0.8,0.9,0.9);
  leg.AddEntry(TH2Bkg_R_Var12_1DY,"signal sample","l");
  leg.AddEntry(TH2Bkg_R_ControlVar12_1DY,"control sample","p");
  leg.SetFillColor(0); 
  leg.Draw();

  c.Print(fname,"Landscape");  

  TH2Bkg_R_Var12_1DY->Scale(scale);
  TH2Bkg_R_ControlVar12_1DY->Scale(scaleControl);
 
}

void PtEtaBin::DrawBasicFPlot(TString *pathname){

  TString fname = *pathname;
  fname += "Control_BasicF_plot";
  //fname += "ControlPlot";
  GiveName(&fname);
  fname +=".eps";
  TCanvas c("c","",1);
  c.cd();
  c.SetGrid();

  /*  TH1Data_ControlVar->Rebin(5);
  TH1Data_ControlVarReweigh->Rebin(5);
  TH1Bkg_Var0->Rebin(5);
  */
 
  cout << "PtEtaBin::DrawControlPlot::WARNING : not using underflow to normalize " << endl;
  //double scaleControl=TH1Data_ControlVar->Integral(0,TH1Data_ControlVar->GetXaxis()->GetNbins()+1);
  double scaleControl=TH1Data_ControlVar->Integral(1,TH1Data_ControlVar->GetXaxis()->GetNbins()+1);
  //  double scaleControlReweigh=TH1Data_ControlVarReweigh->Integral(0,TH1Data_ControlVarReweigh->GetXaxis()->GetNbins()+1);
  double scaleControlReweigh=TH1Data_ControlVarReweigh->Integral(1,TH1Data_ControlVarReweigh->GetXaxis()->GetNbins()+1);//*1.1;
  double scaleVar0=TH1Bkg_Var0->Integral(0,TH1Bkg_Var0->GetXaxis()->GetNbins()+1);

  TH1Data_ControlVar->Scale(1/scaleControl);
  TH1Data_ControlVarReweigh->Scale(1/scaleControlReweigh);
  TH1Bkg_Var0->Scale(1/scaleVar0);

  /*  double scaleControlReweigh_W=TH1Bkg_W_ControlVarReweigh->Integral(1,TH1Bkg_W_ControlVarReweigh->GetXaxis()->GetNbins()+1);
  double scaleVar0_W=TH1Bkg_W_Var0->Integral(0,TH1Bkg_W_Var0->GetXaxis()->GetNbins()+1);
  TH1Bkg_W_ControlVarReweigh->Scale((2/(3*scaleControlReweigh_W)));
  TH1Bkg_W_Var0->Scale((2/(3*scaleVar0_W)));

  double scaleControlReweigh_R=TH1Bkg_R_ControlVarReweigh->Integral(1,TH1Bkg_R_ControlVarReweigh->GetXaxis()->GetNbins()+1);
  double scaleVar0_R=TH1Bkg_R_Var0->Integral(0,TH1Bkg_R_Var0->GetXaxis()->GetNbins()+1);
  TH1Bkg_R_ControlVarReweigh->Scale((1/(3*scaleControlReweigh_R)));
  TH1Bkg_R_Var0->Scale((1/(3*scaleVar0_R)));*/

  double scaleControlReweigh_W=TH1Bkg_W_ControlVar->Integral(1,TH1Bkg_W_ControlVar->GetXaxis()->GetNbins()+1);//*1.1;
  double scaleVar0_W=TH1Bkg_W_Var0->Integral(0,TH1Bkg_W_Var0->GetXaxis()->GetNbins()+1);
  TH1Bkg_W_ControlVar->Scale((2/(3*scaleControlReweigh_W)));
  TH1Bkg_W_Var0->Scale((2/(3*scaleVar0_W)));

  double scaleControlReweigh_R=TH1Bkg_R_ControlVar->Integral(1,TH1Bkg_R_ControlVar->GetXaxis()->GetNbins()+1);//*1.1;
  double scaleVar0_R=TH1Bkg_R_Var0->Integral(0,TH1Bkg_R_Var0->GetXaxis()->GetNbins()+1);
  TH1Bkg_R_ControlVar->Scale((1/(3*scaleControlReweigh_R)));
  TH1Bkg_R_Var0->Scale((1/(3*scaleVar0_R)));

  //  TH1Data_ControlVar->Rebin(10);
  //TH1Bkg_Var0->Rebin(10);


 
  TH1Data_ControlVar->SetMarkerColor(kRed);
  TH1Data_ControlVar->SetMarkerStyle(4);
  TH1Data_ControlVar->SetTitle("");
  TH1Data_ControlVar->GetXaxis()->SetTitle("m_{lj}(GeV)");
  TH1Data_ControlVar->GetYaxis()->SetTitle("a.u.");
  TH1Data_ControlVar->Draw();

  TH1Data_ControlVarReweigh->SetMarkerColor(kRed);
  // TH1Data_ControlVarReweigh->SetTitle("");
  //TH1Data_ControlVarReweigh->GetXaxis()->SetTitle("m_{lj}(GeV)");
  //TH1Data_ControlVarReweigh->GetYaxis()->SetTitle("a.u.");
  //TH1Data_ControlVarReweigh->Draw();
  TH1Data_ControlVarReweigh->Draw("same");

  /*  TH1Bkg_W_ControlVarReweigh->SetMarkerColor(kBlue);
  TH1Bkg_W_ControlVarReweigh->Draw("same");

  TH1Bkg_R_ControlVarReweigh->SetMarkerColor(kGreen);
  TH1Bkg_R_ControlVarReweigh->Draw("same");*/
 
  /* TH1Bkg_W_ControlVar->SetMarkerColor(kBlue);
  TH1Bkg_W_ControlVar->Draw("same");

  TH1Bkg_R_ControlVar->SetMarkerColor(kGreen);
  TH1Bkg_R_ControlVar->Draw("same");*/

  TH1Bkg_Var0->SetLineStyle(2);
  TH1Bkg_Var0->SetLineWidth(2);
  TH1Bkg_Var0->Draw("samehistoe");

  /* TH1Bkg_W_Var0->SetLineStyle(2);
  TH1Bkg_W_Var0->SetLineWidth(2);
  TH1Bkg_W_Var0->SetLineColor(kBlue);
  TH1Bkg_W_Var0->Draw("samehisto");

  TH1Bkg_R_Var0->SetLineStyle(2);
  TH1Bkg_R_Var0->SetLineWidth(2);
  TH1Bkg_R_Var0->SetLineColor(kGreen);
  TH1Bkg_R_Var0->Draw("samehisto");*/

  TLegend leg(0.55,0.8,0.9,0.9);
  leg.AddEntry(TH1Data_ControlVar,"Control sample","p");
  leg.AddEntry(TH1Data_ControlVarReweigh,"Control (reweigh)","p");
  leg.AddEntry(TH1Bkg_Var0,"Signal","l");
  //  leg.AddEntry(TH1Bkg_W_ControlVarReweigh,"W, Control (reweigh), not to scale","p");
  //leg.AddEntry(TH1Bkg_W_Var0,"W, Signal, not to scale","l");
  //leg.AddEntry(TH1Bkg_R_ControlVarReweigh,"R, Control (reweigh), not to scale","p");
  //leg.AddEntry(TH1Bkg_R_Var0,"R, Signal, not to scale","l");
  leg.SetFillColor(0); 
  leg.Draw();

  c.Print(fname,"Landscape");  

  TH1Data_ControlVar->Scale(scaleControl);
  TH1Data_ControlVarReweigh->Scale(scaleControlReweigh);
  TH1Bkg_Var0->Scale(scaleVar0);
 
  //TH1Bkg_W_ControlVarReweigh->Scale((2*scaleControlReweigh_W));
  //TH1Bkg_W_Var0->Scale((2*scaleVar0_W));
  //TH1Bkg_R_ControlVarReweigh->Scale((2*scaleControlReweigh_R));
  //TH1Bkg_R_Var0->Scale((2*scaleVar0_R));
  
  TH1Bkg_W_ControlVar->Scale(1/(2/(3*scaleControlReweigh_W)));
  TH1Bkg_W_Var0->Scale(1/(2/(3*scaleVar0_W)));

  TH1Bkg_R_ControlVar->Scale(1/(1/(3*scaleControlReweigh_R)));
  TH1Bkg_R_Var0->Scale(1/(1/(3*scaleVar0_R)));


}

void PtEtaBin::DrawControlOnlyRPlot(TString *pathname){
 
  TString fname = *pathname;
  fname += "ControlOnlyRPlot";
  //fname += "ControlPlot";
  GiveName(&fname);
  fname +=".eps";
  TCanvas c("c","",1);
  c.cd();
  c.SetGrid();
 
  double scale=0;
  double scaleRW=0;
  double scaleSignal=0;

  scale=TH1Bkg_R_ControlVar->Integral(1,TH1Bkg_R_ControlVar->GetXaxis()->GetNbins()+1);
  scaleRW=TH1Bkg_R_ControlVarReweigh->Integral(1,TH1Bkg_R_ControlVarReweigh->GetXaxis()->GetNbins()+1);
  scaleSignal=TH1Bkg_R_Var0->Integral(1,TH1Bkg_R_Var0->GetXaxis()->GetNbins()+1);

  TH1Bkg_R_ControlVar->Scale(1/scale);
  TH1Bkg_R_ControlVarReweigh->Scale(1/scaleRW);
  TH1Bkg_R_Var0->Scale(1/scaleSignal);

  //  TH1Sng_ControlVar->SetMarkerColor(kRed);
  TH1Bkg_R_ControlVar->SetTitle("");
  TH1Bkg_R_ControlVar->GetXaxis()->SetTitle("m_{lj}(GeV)");
  TH1Bkg_R_ControlVar->GetYaxis()->SetTitle("a.u.");
  TH1Bkg_R_ControlVar->Draw("histo");
  TH1Bkg_R_ControlVar->SetLineColor(7);


  TH1Bkg_R_ControlVarReweigh->SetLineColor(6);
  TH1Bkg_R_ControlVarReweigh->Draw("samehisto");
 
  TH1Bkg_R_Var0->Draw("samehisto");

  TLegend leg(0.55,0.8,0.9,0.9);
  leg.AddEntry(TH1Bkg_R_Var0,"R, signal sample","l");
  leg.AddEntry(TH1Bkg_R_ControlVar,"R, control sample","l");
  leg.AddEntry(TH1Bkg_R_ControlVarReweigh,"R, control (reweigh)","l");
  leg.SetFillColor(0); 
  leg.Draw();

  c.Print(fname,"Landscape");  

  TH1Bkg_R_ControlVar->Scale(scale);
  TH1Bkg_R_ControlVarReweigh->Scale(scaleRW);
  TH1Bkg_R_Var0->Scale(scaleSignal);

}

void PtEtaBin::DrawControlOnlyWPlot(TString *pathname){
 
  TString fname = *pathname;
  fname += "ControlOnlyWPlot";
  //fname += "ControlPlot";
  GiveName(&fname);
  fname +=".eps";
  TCanvas c("c","",1);
  c.cd();
  c.SetGrid();
 
  double scale=0;
  double scaleRW=0;
  double scaleSignal=0;

  scale=TH1Bkg_W_ControlVar->Integral(1,TH1Bkg_W_ControlVar->GetXaxis()->GetNbins()+1);
  scaleRW=TH1Bkg_W_ControlVarReweigh->Integral(1,TH1Bkg_W_ControlVarReweigh->GetXaxis()->GetNbins()+1);
  scaleSignal=TH1Bkg_W_Var0->Integral(1,TH1Bkg_W_Var0->GetXaxis()->GetNbins()+1);

  TH1Bkg_W_ControlVar->Scale(1/scale);
  TH1Bkg_W_ControlVarReweigh->Scale(1/scaleRW);
  TH1Bkg_W_Var0->Scale(1/scaleSignal);

  //  TH1Sng_ControlVar->SetMarkerColor(kRed);
  TH1Bkg_W_ControlVar->SetTitle("");
  TH1Bkg_W_ControlVar->GetXaxis()->SetTitle("m_{lj}(GeV)");
  TH1Bkg_W_ControlVar->GetYaxis()->SetTitle("a.u.");
  TH1Bkg_W_ControlVar->Draw("histo");
  TH1Bkg_W_ControlVar->SetLineColor(7);


  TH1Bkg_W_ControlVarReweigh->SetLineColor(6);
  TH1Bkg_W_ControlVarReweigh->Draw("samehisto");
 
  TH1Bkg_W_Var0->Draw("samehisto");

  TLegend leg(0.55,0.8,0.9,0.9);
  leg.AddEntry(TH1Bkg_W_Var0,"W, signal sample","l");
  leg.AddEntry(TH1Bkg_W_ControlVar,"W, control sample","l");
  leg.AddEntry(TH1Bkg_W_ControlVarReweigh,"W, control (reweigh)","l");
  leg.SetFillColor(0); 
  leg.Draw();

  c.Print(fname,"Landscape");   

  TH1Bkg_W_ControlVar->Scale(scale);
  TH1Bkg_W_ControlVarReweigh->Scale(scaleRW);
  TH1Bkg_W_Var0->Scale(scaleSignal);

}

void PtEtaBin::DrawControlFPlot(TString *pathname){
 
  TString fname = *pathname;
  fname += "ControlFPlot";
  //fname += "ControlPlot";
  GiveName(&fname);
  fname +=".eps";
  TCanvas c("c","",1);
  c.cd();
  c.SetGrid();
 
  double scaleb=0;
  double scalenonb=0;
  double scaleR=0;
  double scaleW=0;

  scaleb=TH1Sng_ControlVar->Integral(1,TH1Sng_ControlVar->GetXaxis()->GetNbins()+1);
  scalenonb=TH1Bkg_ControlVar->Integral(1,TH1Bkg_ControlVar->GetXaxis()->GetNbins()+1);
  scaleR=TH1Bkg_R_ControlVar->Integral(1,TH1Bkg_R_ControlVar->GetXaxis()->GetNbins()+1);
  scaleW=TH1Bkg_W_ControlVar->Integral(1,TH1Bkg_W_ControlVar->GetXaxis()->GetNbins()+1);

  TH1Sng_ControlVar->Scale(1/scaleb);
  TH1Bkg_ControlVar->Scale(1/scalenonb);
  TH1Bkg_R_ControlVar->Scale(1/scaleR);
  TH1Bkg_W_ControlVar->Scale(1/scaleW);

  //  TH1Sng_ControlVar->SetMarkerColor(kRed);
  TH1Sng_ControlVar->SetTitle("");
  TH1Sng_ControlVar->GetXaxis()->SetTitle("m_{lj}(GeV)");
  TH1Sng_ControlVar->GetYaxis()->SetTitle("a.u.");
  TH1Sng_ControlVar->Draw("histo");

  TH1Bkg_ControlVar->SetLineColor(kRed);
  TH1Bkg_ControlVar->Draw("samehisto");

  TH1Bkg_W_ControlVar->SetLineColor(kGreen);
  TH1Bkg_W_ControlVar->Draw("samehisto");

  TH1Bkg_R_ControlVar->SetLineColor(kBlue);
  TH1Bkg_R_ControlVar->Draw("samehisto");

  TLegend leg(0.55,0.8,0.9,0.9);
  leg.AddEntry(TH1Sng_ControlVar,"b, control sample","l");
  leg.AddEntry(TH1Bkg_ControlVar,"non b, control sample","l");
  leg.AddEntry(TH1Bkg_R_ControlVar,"R, control sample","l");
  leg.AddEntry(TH1Bkg_W_ControlVar,"W, control sample","l");
  leg.SetFillColor(0); 
  leg.Draw();

  c.Print(fname,"Landscape");  

  TH1Sng_ControlVar->Scale(scaleb);
  TH1Bkg_ControlVar->Scale(scalenonb);
  TH1Bkg_R_ControlVar->Scale(scaleR);
  TH1Bkg_W_ControlVar->Scale(scaleW);

}

void PtEtaBin::DrawControlRWFPlot(TString *pathname){
 
  TString fname = *pathname;
  fname += "ControlRWFPlot";
  //fname += "ControlPlot";
  GiveName(&fname);
  fname +=".eps";
  TCanvas c("c","",1);
  c.cd();
  c.SetGrid();
 
  double scaleb=0;
  double scalenonb=0;
  double scaleR=0;
  double scaleW=0;

  scaleb=TH1Sng_ControlVarReweigh->Integral(1,TH1Sng_ControlVarReweigh->GetXaxis()->GetNbins()+1);
  scalenonb=TH1Bkg_ControlVarReweigh->Integral(1,TH1Bkg_ControlVarReweigh->GetXaxis()->GetNbins()+1);
  scaleR=TH1Bkg_R_ControlVarReweigh->Integral(1,TH1Bkg_R_ControlVarReweigh->GetXaxis()->GetNbins()+1);
  scaleW=TH1Bkg_W_ControlVarReweigh->Integral(1,TH1Bkg_W_ControlVarReweigh->GetXaxis()->GetNbins()+1);

  TH1Sng_ControlVarReweigh->Scale(1/scaleb);
  TH1Bkg_ControlVarReweigh->Scale(1/scalenonb);
  TH1Bkg_R_ControlVarReweigh->Scale(1/scaleR);
  TH1Bkg_W_ControlVarReweigh->Scale(1/scaleW);

  //  TH1Sng_ControlVarReweigh->SetMarkerColor(kRed);
  TH1Sng_ControlVarReweigh->SetTitle("");
  TH1Sng_ControlVarReweigh->GetXaxis()->SetTitle("m_{lj}(GeV)");
  TH1Sng_ControlVarReweigh->GetYaxis()->SetTitle("a.u.");
  TH1Sng_ControlVarReweigh->Draw("histo");

  TH1Bkg_ControlVarReweigh->SetLineColor(kRed);
  TH1Bkg_ControlVarReweigh->Draw("samehisto");

  TH1Bkg_W_ControlVarReweigh->SetLineColor(kGreen);
  TH1Bkg_W_ControlVarReweigh->Draw("samehisto");

  TH1Bkg_R_ControlVarReweigh->SetLineColor(kBlue);
  TH1Bkg_R_ControlVarReweigh->Draw("samehisto");

  TLegend leg(0.55,0.8,0.9,0.9);
  leg.AddEntry(TH1Sng_ControlVarReweigh,"b, control (reweigh)","l");
  leg.AddEntry(TH1Bkg_ControlVarReweigh,"non b, control (reweigh)","l");
  leg.AddEntry(TH1Bkg_R_ControlVarReweigh,"R, control (reweigh)","l");
  leg.AddEntry(TH1Bkg_W_ControlVarReweigh,"W, control (reweigh)","l");
  leg.SetFillColor(0); 
  leg.Draw();

  c.Print(fname,"Landscape"); 

  TH1Sng_ControlVarReweigh->Scale(scaleb);
  TH1Bkg_ControlVarReweigh->Scale(scalenonb);
  TH1Bkg_R_ControlVarReweigh->Scale(scaleR);
  TH1Bkg_W_ControlVarReweigh->Scale(scaleW);

}

void PtEtaBin::DrawSignalFPlot(TString *pathname){
 
  TString fname = *pathname;
  fname += "ControlSignalFPlot";
  //fname += "ControlPlot";
  GiveName(&fname);
  fname +=".eps";
  TCanvas c("c","",1);
  c.cd();
  c.SetGrid();
 
  double scaleb=0;
  double scalenonb=0;
  double scaleR=0;
  double scaleW=0;

  scaleb=TH1Sng_Var0->Integral(1,TH1Sng_Var0->GetXaxis()->GetNbins()+1);
  scalenonb=TH1Bkg_Var0->Integral(1,TH1Bkg_Var0->GetXaxis()->GetNbins()+1);
  scaleR=TH1Bkg_R_Var0->Integral(1,TH1Bkg_R_Var0->GetXaxis()->GetNbins()+1);
  scaleW=TH1Bkg_W_Var0->Integral(1,TH1Bkg_W_Var0->GetXaxis()->GetNbins()+1);

  TH1Sng_Var0->Scale(1/scaleb);
  TH1Bkg_Var0->Scale(1/scalenonb);
  TH1Bkg_R_Var0->Scale(1/scaleR);
  TH1Bkg_W_Var0->Scale(1/scaleW);

  //  TH1Sng_Var0->SetMarkerColor(kRed);
  /*  TH1Sng_Var0->SetTitle("");
  TH1Sng_Var0->GetXaxis()->SetTitle("m_{lj}(GeV)");
  TH1Sng_Var0->GetYaxis()->SetTitle("a.u.");
  TH1Sng_Var0->Draw("histo");*/

  TH1Bkg_Var0->SetLineColor(kRed);
  TH1Bkg_Var0->Draw("histo");
  TH1Bkg_Var0->SetTitle("");
  TH1Bkg_Var0->GetXaxis()->SetTitle("m_{lj}(GeV)");
  TH1Bkg_Var0->GetYaxis()->SetTitle("a.u.");

  TH1Bkg_W_Var0->SetLineColor(kGreen);
  TH1Bkg_W_Var0->Draw("samehisto");

  TH1Bkg_R_Var0->SetLineColor(kBlue);
  TH1Bkg_R_Var0->Draw("samehisto");

  TLegend leg(0.55,0.8,0.9,0.9);
  leg.AddEntry(TH1Sng_Var0,"b, signal sample","l");
  leg.AddEntry(TH1Bkg_Var0,"non b, signal sample","l");
  leg.AddEntry(TH1Bkg_R_Var0,"R, signal sample","l");
  leg.AddEntry(TH1Bkg_W_Var0,"W, signal sample","l");
  leg.SetFillColor(0); 
  leg.Draw();

  c.Print(fname,"Landscape");   

  TH1Sng_Var0->Scale(scaleb);
  TH1Bkg_Var0->Scale(scalenonb);
  TH1Bkg_R_Var0->Scale(scaleR);
  TH1Bkg_W_Var0->Scale(scaleW);

}

/*void PtEtaBin::DrawControlPlot(){
 
  TString fname = *pathname;
  fname += "ControlRWPlot";
  //fname += "ControlPlot";
  GiveName(&fname);
  fname +=".eps";
  TCanvas c("c","",1);
  c.cd();
  c.SetGrid();
   
  TLegend leg(0.55,0.8,0.9,0.9);
  leg.AddEntry(TH1Bkg_R_Var0,"R, Signal, not to scale","l");
  leg.SetFillColor(0); 
  leg.Draw();

  c.Print(fname,"Landscape");  
  }*/

void PtEtaBin::DrawbtagNR(TString *pathname){

  TString fname = *pathname;
  fname += "BtagLR_NR";
  GiveName(&fname);
  fname +=".eps";
  TCanvas c("c","",1);
  c.cd();
  c.SetGrid();
  //c.SetLogy();

  cout << "Underflow left  : " << TH2Data_Left1DY->GetBinContent(0) << endl;
  cout << "Underflow right : " << TH2Data_Right1DY->GetBinContent(0) << endl;

  //double scaleLeft=TH2Data_Left1DY->Integral(0,TH2Data_Left1DY->GetXaxis()->GetNbins()+1);
  //double scaleRight=TH2Data_Right1DY->Integral(0,TH2Data_Right1DY->GetXaxis()->GetNbins()+1);

  //TH2Data_Left1DY->Scale(1/scaleLeft);
  //TH2Data_Right1DY->Scale(1/scaleRight);

  TH2Data_Left1DY->SetTitle("");
  TH2Data_Left1DY->SetMarkerColor(3);
  TH2Data_Right1DY->SetMarkerColor(4);
  //TH2Data_Left1DY->SetMarkerStyle(1);
  //TH2Data_Right1DY->SetMarkerStyle(1);
  TH2Data_Left1DY->Draw("");  
  TH2Data_Left1DY->GetXaxis()->SetTitle("b-tag discriminant");
  TH2Data_Left1DY->GetYaxis()->SetTitle("# jets/fb^{-1}");
  TH2Data_Right1DY->Draw("same");
  
  TLegend leg(0.60,0.8,0.9,0.9);
  leg.AddEntry(TH2Data_Left1DY,"Signal sample - Left","p");
  leg.AddEntry(TH2Data_Right1DY,"Signal sample Right","p");
  leg.SetFillColor(0); 
  leg.Draw();

  /* double scaleLeft=TH2Bkg_Left1DY->Integral(0,TH2Bkg_Left1DY->GetXaxis()->GetNbins()+1);
  double scaleRight=TH2Bkg_Right1DY->Integral(0,TH2Bkg_Right1DY->GetXaxis()->GetNbins()+1);

  TH2Bkg_Left1DY->Scale(1/scaleLeft);
  TH2Bkg_Right1DY->Scale(1/scaleRight);

  TH2Bkg_Left1DY->SetTitle("");
  TH2Bkg_Left1DY->SetMarkerColor(3);
  TH2Bkg_Right1DY->SetMarkerColor(4);
  //TH2Bkg_Left1DY->SetMarkerStyle(1);
  //TH2Bkg_Right1DY->SetMarkerStyle(1);
  TH2Bkg_Left1DY->GetXaxis()->SetRangeUser(-2,4);  
  TH2Bkg_Left1DY->Draw("");  
  TH2Bkg_Left1DY->GetXaxis()->SetTitle("b-tag discriminant");
  TH2Bkg_Left1DY->GetYaxis()->SetTitle("# jets/fb^{-1}");
  TH2Bkg_Right1DY->Draw("same");
  
  TLegend leg(0.05,0.8,0.35,0.9);
  leg.AddEntry(TH2Bkg_Left1DY,"Signal sample (non b) - Left","p");
  leg.AddEntry(TH2Bkg_Right1DY,"Signal sample (non b) Right","p");
  leg.SetFillColor(0); 
  leg.Draw();*/

  c.Print(fname,"Landscape");  
  //TH2Bkg_Left1DY->Scale(scaleLeft);
  //TH2Bkg_Right1DY->Scale(scaleRight);

  //TH2Data_Left1DY->Scale(scaleLeft);
  //TH2Data_Right1DY->Scale(scaleRight);

 
  TString fname2 = *pathname;
  fname2 += "Btag_NR";
  GiveName(&fname2);
  fname2 +=".eps";
  TCanvas c2("c2","",1);
  c2.cd();
  c2.SetGrid();
  c2.SetLogy();

  cout << "Underflow MCmeasured : " << TH1Data_BtagMCMeasured->GetBinContent(0) << endl;
 
  double scaleMeasured=TH1Data_BtagMCMeasured->Integral(0,TH1Data_BtagMCMeasured->GetXaxis()->GetNbins()+1);
  double scaleMC=TH1Sng_BtagAll->Integral(0,TH1Sng_BtagAll->GetXaxis()->GetNbins()+1);
  TH1Sng_BtagAll->Scale(1/scaleMC);
  TH1Data_BtagMCMeasured->Scale(1/scaleMeasured);

  cout << "Underflow MCmeasured (norm) : " << TH1Data_BtagMCMeasured->GetBinContent(0) << endl;
  cout << "Underflow expected (norm)   : " << TH1Sng_BtagAll->GetBinContent(0) << endl;

  TH1Sng_BtagAll->SetTitle("");
  TH1Sng_BtagAll->SetMarkerColor(1);
  TH1Data_BtagMCMeasured->SetMarkerColor(2);
  TH1Sng_BtagAll->GetXaxis()->SetTitle("b-tag discriminant");
  TH1Sng_BtagAll->GetYaxis()->SetTitle("a.u.");
  TH1Sng_BtagAll->Draw("histo");
  TH1Data_BtagMCMeasured->Draw("samePE");
  
  TLegend leg2(0.65,0.8,0.9,0.9);
  leg2.AddEntry(TH1Sng_BtagAll,"MC","l");
  leg2.AddEntry(TH1Data_BtagMCMeasured,"Method(F_{MC})","p");
  leg2.SetFillColor(0); 
  leg2.Draw();

  c2.Print(fname2,"Landscape");  

  TH1Sng_BtagAll->Scale(scaleMC);
  TH1Data_BtagMCMeasured->Scale(scaleMeasured);


}

void PtEtaBin::DrawEffNR(TString *pathname){

  TString fname = *pathname;
  fname += "Btag_Eff_NR";
  GiveName(&fname);
  fname +=".eps";
  TCanvas c("c","",1);
  c.cd();
  c.SetGrid();
 
  TH1Sng_BtagEffMC->SetTitle("");
  TH1Sng_BtagEffMC->GetXaxis()->SetTitle("b-tag cut");
  TH1Sng_BtagEffMC->GetYaxis()->SetTitle("#epsilon_{b,taggable}");
  //TH1Sng_BtagEffMC->SetLineWidth(2);
  TH1Sng_BtagEffMC->Draw("histo");
  
  TH1Data_BtagEffMeasured->Draw("PEsame");
  TH1Data_BtagEffMeasured->SetMarkerColor(6);

  //TH1Data_BtagEffMCMeasured->Draw("samehisto");
  //TH1Data_BtagEffMCMeasured->SetLineColor(kRed);
  TH1Data_BtagEffMCMeasured->Draw("samePE");
  TH1Data_BtagEffMCMeasured->SetMarkerColor(2);
    
  TLegend leg(0.65,0.8,0.9,0.9);
  leg.AddEntry(TH1Sng_BtagEffMC,"MC","l");
  //leg.AddEntry(TH1Data_BtagEffMCMeasured,"Method(F_{MC})","l");
  leg.AddEntry(TH1Data_BtagEffMCMeasured,"Method(F_{MC})","p");
  //leg.AddEntry(TH1Data_BtagEffMeasured,"Method","l");
  leg.SetFillColor(0); 
  leg.Draw();

  c.Print(fname,"Landscape");
}

void PtEtaBin::DrawEffDiffNR(TString *pathname){

  TString fname = *pathname;
  fname += "Btag_Eff_Diff_NR";
  GiveName(&fname);
  fname +=".eps";
  TCanvas c("c","",1);
  c.cd();
  c.SetGrid();

  TH1Data_BtagEffMCMeasuredDiff->SetTitle(""); 
  TH1Data_BtagEffMCMeasuredDiff->GetXaxis()->SetTitle("b-tag cut");
  TH1Data_BtagEffMCMeasuredDiff->GetYaxis()->SetTitle("(#epsilon_{b}-#epsilon_{b,MC})_{taggable}");
  //  TH1Data_BtagEffMCMeasuredDiff->Draw("histo");
  //TH1Data_BtagEffMCMeasuredDiff->SetLineColor(kRed);
  TH1Data_BtagEffMCMeasuredDiff->Draw("PE");
  TH1Data_BtagEffMCMeasuredDiff->SetMarkerColor(kRed);

  TH1Data_BtagEffMeasuredDiff->Draw("PEsame");
  TH1Data_BtagEffMeasuredDiff->SetMarkerColor(6);
  
  TLegend leg(0.65,0.8,0.9,0.9);
  //leg.AddEntry(TH1Data_BtagEffMCMeasuredDiff,"Method(F_{MC})","l");
  leg.AddEntry(TH1Data_BtagEffMCMeasuredDiff,"Method(F_{MC})","p");
  //leg.AddEntry(TH1Data_BtagEffMeasuredDiff,"b-tag Efficiency measured","l");
  leg.SetFillColor(0); 
  leg.Draw();

  c.Print(fname,"Landscape");
}

void PtEtaBin::DrawEffLR(TString *pathname){

  TString fname = *pathname;
  fname += "Btag_Eff_LR";
  GiveName(&fname);
  fname +=".eps";
  TCanvas c("c","",1);
  c.cd();
  
  TH1Sng_BtagEffMC->SetTitle("");
  TH1Sng_BtagEffMC->Draw("histo");
  
  TH1Data_BtagEffMeasuredLR->Draw("same");
  TH1Data_BtagEffMeasuredLR->SetLineColor(6);

  TH1Data_BtagEffMCMeasuredLR->Draw("samehisto");
  TH1Data_BtagEffMCMeasuredLR->SetLineColor(kRed);
  
  TLegend leg(0.65,0.8,0.9,0.9);
  leg.AddEntry(TH1Sng_BtagEffMC,"b-tag Efficiency full MC","l");
  leg.AddEntry(TH1Data_BtagEffMCMeasuredLR,"b-tag Efficiency F from MC","l");
  leg.AddEntry(TH1Data_BtagEffMeasuredLR,"b-tag Efficiency measured","l");
  leg.SetFillColor(0); 
  leg.Draw();

  c.Print(fname,"Landscape");
}

void PtEtaBin::DrawEffDiffLR(TString *pathname){

  TString fname = *pathname;
  fname += "Btag_Eff_Diff_LR";
  GiveName(&fname);
  fname +=".eps";
  TCanvas c("c","",1);
  c.cd();
 
  TH1Data_BtagEffMeasuredLRDiff->SetTitle("");
  TH1Data_BtagEffMeasuredLRDiff->Draw();
  TH1Data_BtagEffMeasuredLRDiff->SetLineColor(6);

  TH1Data_BtagEffMCMeasuredLRDiff->Draw("samehisto");
  TH1Data_BtagEffMCMeasuredLRDiff->SetLineColor(kRed);
  
  TLegend leg(0.65,0.8,0.9,0.9);
  leg.AddEntry(TH1Data_BtagEffMCMeasuredLRDiff,"b-tag Efficiency F from MC","l");
  leg.AddEntry(TH1Data_BtagEffMeasuredLRDiff,"b-tag Efficiency measured","l");
  leg.SetFillColor(0); 
  leg.Draw();

  c.Print(fname,"Landscape");
}

void PtEtaBin::DrawEffRR(TString *pathname){
  TString fname = *pathname;
  fname += "Btag_Eff_RR";
  GiveName(&fname);
  fname +=".eps";
  TCanvas c("c","",1);
  c.cd();
  c.SetGrid();
 
  TH1Sng_BtagEffMC->SetTitle("");
  TH1Sng_BtagEffMC->GetXaxis()->SetTitle("b-tag cut");
  TH1Sng_BtagEffMC->GetYaxis()->SetTitle("#epsilon_{b}");
  //TH1Sng_BtagEffMC->SetLineWidth(2);
  TH1Sng_BtagEffMC->Draw("histo");
  
  TH1Data_BtagEffMeasuredRR->Draw("samePE");
  TH1Data_BtagEffMeasuredRR->SetMarkerColor(6);

  //TH1Data_BtagEffMCMeasured->Draw("samehisto");
  //TH1Data_BtagEffMCMeasured->SetLineColor(kRed);
  TH1Data_BtagEffMCMeasuredRR->Draw("samePE");
  TH1Data_BtagEffMCMeasuredRR->SetMarkerColor(2);
    
  TLegend leg(0.65,0.8,0.9,0.9);
  leg.AddEntry(TH1Sng_BtagEffMC,"MC","l");
  //leg.AddEntry(TH1Data_BtagEffMCMeasuredRR,"Method(F_{MC})","l");
  leg.AddEntry(TH1Data_BtagEffMCMeasuredRR,"Method(F_{MC})","p");
  //leg.AddEntry(TH1Data_BtagEffMeasuredRR,"Method","l");
  leg.AddEntry(TH1Data_BtagEffMeasuredRR,"full method","p");
  leg.SetFillColor(0); 
  leg.Draw();

  c.Print(fname,"Landscape");


}

void PtEtaBin::DrawEffDiffRR(TString *pathname){
  TString fname = *pathname;
  fname += "Btag_Eff_Diff_RR";
  GiveName(&fname);
  fname +=".eps";
  TCanvas c("c","",1);
  c.cd();
  c.SetGrid();

  TH1Data_BtagEffMCMeasuredRRDiff->SetTitle(""); 
  TH1Data_BtagEffMCMeasuredRRDiff->GetXaxis()->SetTitle("b-tag cut");
  TH1Data_BtagEffMCMeasuredRRDiff->GetYaxis()->SetTitle("#epsilon_{b}-#epsilon_{b,MC}");
  //  TH1Data_BtagEffMCMeasuredRRDiff->Draw("histo");
  //TH1Data_BtagEffMCMeasuredRRDiff->SetLineColor(kRed);

  TH1Data_BtagEffMeasuredRRDiff->Draw("PE");
  TH1Data_BtagEffMeasuredRRDiff->SetTitle("");
  TH1Data_BtagEffMeasuredRRDiff->GetXaxis()->SetTitle("b-tag cut");
  TH1Data_BtagEffMeasuredRRDiff->GetYaxis()->SetTitle("epsilon_{b}-#epsilon_{b,MC}");
  TH1Data_BtagEffMeasuredRRDiff->SetMarkerColor(6);

  TH1Data_BtagEffMCMeasuredRRDiff->Draw("PEsame");
  TH1Data_BtagEffMCMeasuredRRDiff->SetMarkerColor(kRed);

  
  TLegend leg(0.65,0.8,0.9,0.9);
  //leg.AddEntry(TH1Data_BtagEffMCMeasuredRRDiff,"Method(F_{MC})","l");
  leg.AddEntry(TH1Data_BtagEffMCMeasuredRRDiff,"Method(F_{MC})","p");
  leg.AddEntry(TH1Data_BtagEffMeasuredRRDiff,"full method","p");
  leg.SetFillColor(0); 
  leg.Draw();

  c.Print(fname,"Landscape");
 

}

void PtEtaBin::GiveName(TString *name){
  
  name->Append(*genericName);

}

double PtEtaBin::getmljMeanVal(){return TH1Data_ControlVarReweigh->GetMean();}
double PtEtaBin::getmljMeannoRWVal(){return TH1Data_ControlVar->GetMean();}
double PtEtaBin::getmljMeanMCVal(){return TH1Bkg_Var0->GetMean();}
double PtEtaBin::getmljSigmaVal(){return TH1Data_ControlVarReweigh->GetRMS();}
double PtEtaBin::getmljSigmanoRWVal(){return TH1Data_ControlVar->GetRMS();}
double PtEtaBin::getmljSigmaMCVal(){return TH1Bkg_Var0->GetRMS();}

double PtEtaBin::getmlj_W_MeanVal(){return TH1Bkg_W_ControlVarReweigh->GetMean();}
double PtEtaBin::getmlj_W_MeannoRWVal(){return TH1Bkg_W_ControlVar->GetMean();}
double PtEtaBin::getmlj_W_MeanMCVal(){return TH1Bkg_W_Var0->GetMean();}
double PtEtaBin::getmlj_W_SigmaVal(){return TH1Bkg_W_ControlVarReweigh->GetRMS();}
double PtEtaBin::getmlj_W_SigmanoRWVal(){return TH1Bkg_W_ControlVar->GetRMS();}
double PtEtaBin::getmlj_W_SigmaMCVal(){return TH1Bkg_W_Var0->GetRMS();}

double PtEtaBin::getmlj_R_MeanVal(){return TH1Bkg_R_ControlVarReweigh->GetMean();}
double PtEtaBin::getmlj_R_MeannoRWVal(){return TH1Bkg_R_ControlVar->GetMean();}
double PtEtaBin::getmlj_R_MeanMCVal(){return TH1Bkg_R_Var0->GetMean();}
double PtEtaBin::getmlj_R_SigmaVal(){return TH1Bkg_R_ControlVarReweigh->GetRMS();}
double PtEtaBin::getmlj_R_SigmanoRWVal(){return TH1Bkg_R_ControlVar->GetRMS();}
double PtEtaBin::getmlj_R_SigmaMCVal(){return TH1Bkg_R_Var0->GetRMS();}

vector<float> PtEtaBin::doTemplateFit (TH1D* ttbar, TH1D* vvmc, TH1D* vvdata,TH1D* data,TString PrefixPlot) {
	
	if (debug_>1) cout << "nTTbar = " << nTTbar_ << " WAS FILLED" << endl;
    
    if (debug_>1) cout << "Doing template fit for L=" << lumi_ << ", the templates where created with L=" << templateLumi_ << " (ttbar: " << ttbar->Integral() << ", VV: " << vvmc->Integral() << ")" << endl;
	
	float ScaleTTbar = ttbar->Integral();
	float ScaleVV = vvmc->Integral();
	float ScaleVVD = vvdata->Integral();
    
	//cout << ScaleTTbar << " " << ScaleVV << endl;
	
    // start values for the fit
    double min = 0.5;
    double max = 1.5;
    double start_ttb = ScaleTTbar*(lumi_/templateLumi_);
    double start_VV = ScaleVV*(lumi_/templateLumi_);
    
    if (debug_>1) cout << "ttbar Fit parameters -> start: " << start_ttb << " min: " << start_ttb*min << " max: " << start_ttb*max << endl;
    if (debug_>1) cout << "VV Fit parameters -> start: " << start_VV << " min: " << start_VV*min << " max: " << start_VV*max << endl;
    
	ttbar->Scale(1/ScaleTTbar);
	vvmc->Scale(1./ScaleVV);
	vvdata->Scale(1./ScaleVVD);
	
	RooRealVar mlb("mlb","m_{#muj} (GeV/c^{2})",20,500);
	RooDataHist ttDataHist("ttDataHist","ttDataHist",mlb,ttbar);
	//RooDataHist ttDataHist("ttDataHist","ttDataHist",mlb,ttTemplateHist);
	RooHistPdf ttPDF("ttPDF","ttPDF",mlb,ttDataHist);
	
	RooDataHist VDataHist("VDataHist","VDataHist",mlb,vvmc);
	//RooDataHist VDataHist("VDataHist","VDataHist",mlb,VTemplateHist);
	RooHistPdf VPDF("VPDF","VPDF",mlb,VDataHist);
	
	RooDataHist VDataHist_control("VDataHist_control","VDataHist_control",mlb,vvdata);
	RooHistPdf VPDF_control("VPDF_control","VPDF_control",mlb,VDataHist_control);
	
	RooRealVar Ntt("Ntt","N_{t#bar{t}-like}",start_ttb,start_ttb*min,start_ttb*max);
    RooRealVar NV( "NV", "N_{V-like}",start_VV,start_VV*min,start_VV*max);

    //RooRealVar Ntt("Ntt","N_{t#bar{t}-like}",1000,0,10000);
	//RooRealVar NV( "NV", "N_{V-like}",1000,0,10000);
	
	RooRealVar Ntt_control("Ntt_control","N_{t#bar{t}-like}_control",start_ttb,start_ttb*min,start_ttb*max);
    RooRealVar NV_control( "NV_control", "N_{V-like}_control",start_VV,start_VV*min,start_VV*max);

    //RooRealVar Ntt_control("Ntt_control","N_{t#bar{t}-like}_control",1000,0,10000);
	//RooRealVar NV_control( "NV_control", "N_{V-like}_control",1000,0,10000);
	
	RooAddPdf model("model","model",RooArgList(ttPDF,VPDF),RooArgList(Ntt,NV));
	RooAddPdf model_control("model","model",RooArgList(ttPDF,VPDF_control),RooArgList(Ntt_control,NV_control));
	
	RooDataHist dataHist("dataHist","dataHist",mlb,data);
		
	model.fitTo(dataHist,RooFit::SumW2Error(false), RooFit::PrintLevel(-3), RooFit::Verbose(false),RooFit::Extended(true));

	//model_control.fitTo(dataHist, RooFit::SumW2Error(false), RooFit::PrintLevel(-3), RooFit::Verbose(false),RooFit::Extended(true));
	
	if (debug_ > 0) cout << "Estimation of Ntt  :  " << Ntt.getVal() << " +- " << Ntt.getError() << "      " << Ntt.getErrorHi() << "      " << Ntt.getErrorLo() << endl;
	if (debug_ > 0) cout << "Estimation of NV   :  " << NV.getVal() << " +- " << NV.getError() << "      " << NV.getErrorHi() << "      " << NV.getErrorLo() << endl;
	//cout << "Estimation of Ntt_control  :  " << Ntt_control.getVal() << " +- " << Ntt_control.getError() << "      " << Ntt_control.getErrorHi() << "      " << Ntt_control.getErrorLo() << endl;
	//cout << "Estimation of NV_control   :  " << NV_control.getVal() << " +- " << NV_control.getError() << "      " << NV_control.getErrorHi() << "      " << NV_control.getErrorLo() << endl;
	
	vector<float> fitResults;
	fitResults.push_back(Ntt.getVal());
	fitResults.push_back(Ntt.getError());
	//fitResults.push_back(Ntt.getErrorHi());
	//fitResults.push_back(Ntt.getErrorLo());
	
	fitResults.push_back(NV.getVal());
	fitResults.push_back(NV.getError());
	//fitResults.push_back(NV.getErrorHi());
	//fitResults.push_back(NV.getErrorLo());
	
	/*fitResults.push_back(Ntt_control.getVal());
	fitResults.push_back(Ntt_control.getError());
	fitResults.push_back(Ntt_control.getErrorHi());
	fitResults.push_back(Ntt_control.getErrorLo());
	
	fitResults.push_back(NV_control.getVal());
	fitResults.push_back(NV_control.getError());
	fitResults.push_back(NV_control.getErrorHi());
	fitResults.push_back(NV_control.getErrorLo());*/
	
	ttbar->Scale(ScaleTTbar);
	vvmc->Scale(ScaleVV);
	vvdata->Scale(ScaleVVD);
	
    if (debug_ > 0) cout << "Data integral " << data->Integral() << endl;
	if (debug_ > 0) cout << "Ratio Fit (VMC template)/TTbar MC template Integral" << Ntt.getVal() << "/" << ttbar->Integral() << "=" << Ntt.getVal()/ttbar->Integral() << endl;
	if (debug_ > 0) cout << "Ratio Fit (VMC template)/V MC template Integral " << NV.getVal() << "/" << vvmc->Integral() << "=" << NV.getVal()/vvmc->Integral() << endl;
	
	//cout << "Ratio Fit (Control sample)/TTbar MC template Integral" << Ntt_control.getVal() << "/" << ttbar->Integral() << "=" << Ntt_control.getVal()/ttbar->Integral() << endl;
	//cout << "Ratio Fit (Control sample)/V MC template Integral " << NV_control.getVal() << "/" << vvmc->Integral() << "=" << NV_control.getVal()/vvmc->Integral() << endl;
	
	//create dirs to store plots
	
	if (debug_ > 0) {
				
		TString dir = ""; GiveName(&dir);
		
		mkdir("FitOutput",0777);
		mkdir(("FitOutput/"+(string)dir).c_str(),0777);
		
		string savePath = ("FitOutput/"+(string)dir+"/");
		
		// give the plot a name
		
		GiveName(&PrefixPlot);	
		cout << PrefixPlot << endl;
		
		// save the MC plot
		
		string MCPlot = "MCTemplate_"+(string)PrefixPlot;
		
		RooPlot* plot = mlb.frame();
		dataHist.plotOn(plot);
		model.plotOn(plot);
		model.plotOn(plot, RooFit::Components(VPDF), RooFit::LineStyle(kDashed),RooFit::LineColor(2));
		model.plotOn(plot, RooFit::Components(ttPDF), RooFit::LineStyle(kDashed),RooFit::LineColor(3));
		
		TCanvas* cFit = new TCanvas(MCPlot.c_str(),MCPlot.c_str());
		cFit->cd();
		
		plot->SetTitle("");
		
		plot->Draw();
		
		TLegend *leg1 = new TLegend(0.65,0.73,0.96,0.97);
		//leg1->SetFillColor(kWhite);
		//leg1->SetLineColor(kWhite);
		
		//leg1->SetHeader(PrefixPlot);
		
		TH1F* dummy = new TH1F("dummy","dummy",1,0,1); 
		dummy->SetLineColor(kBlack); dummy->SetMarkerColor(kBlack); 
		
		leg1->AddEntry(dummy,"Data", "P"); 
		
		TH1F* dummy2 = new TH1F("dummy2","dummy2",1,0,1); 
		dummy2->SetLineColor(kBlue); 
		dummy2->SetMarkerColor(kBlue); 
		
		leg1->AddEntry(dummy2,"t#bar{t} + Background","LP"); 
		
		TH1F* dummy3 = new TH1F("dummy3","dummy3",1,0,1); 
		dummy3->SetLineColor(kGreen); 
		dummy3->SetMarkerColor(kGreen);
		dummy3->SetLineStyle(kDashed);
		
		leg1->AddEntry(dummy3,"t#bar{t}", "LP"); 
		
		TH1F* dummy4 = new TH1F("dummy4","dummy4",1,0,1); 
		dummy4->SetLineColor(kRed); 
		dummy4->SetMarkerColor(kRed); 
		dummy4->SetLineStyle(kDashed);
		
		leg1->AddEntry(dummy4,"Background", "LP"); 
		
		
		leg1->Draw();
		
		
		//TPaveText* text1 = new TPaveText(0.5,0.5,0.8,0.6,"NDC");
		//text1->AddText("CMS Preliminary");
		//text1->Draw();	
		
		cFit->SaveAs((savePath+MCPlot+".png").c_str());
		cFit->SaveAs((savePath+MCPlot+".C").c_str());
		
		TLatex* text = new TLatex(0.13,0.955,"CMS Simulation");
		
		text->SetTextSize(0.05);
		
		text->SetNDC();
		
		text->Draw();

		cFit->SaveAs((savePath+MCPlot+"_RESTYLED_.png").c_str());
		//cFit->Write();
		
		delete plot;
		delete cFit;
		delete leg1;
		delete dummy;
		delete dummy2;
		delete dummy3;
		delete dummy4;
		
		FitPlotPaths.push_back(savePath+MCPlot+".png");
		
		// save the MC plot
		
		/*string CSPlot = "CSTemplate_"+(string)PrefixPlot;
		 
		 RooPlot* plotCS = mlb.frame();
		 dataHist.plotOn(plotCS);
		 model_control.plotOn(plotCS);
		 model_control.plotOn(plotCS, RooFit::Components(VPDF_control), RooFit::LineStyle(kDashed),RooFit::LineColor(2));
		 model_control.plotOn(plotCS, RooFit::Components(ttPDF), RooFit::LineStyle(kDashed),RooFit::LineColor(3));
		 
		 TCanvas* cFitCS = new TCanvas(CSPlot.c_str(),CSPlot.c_str());
		 cFitCS->cd();
		 plotCS->Draw();
		 cFitCS->SaveAs((savePath+CSPlot+".png").c_str());
		 //cFitCS->Write();
		 
		 delete plotCS;
		 delete cFitCS;*/
		
	}
		
	return fitResults;

}
