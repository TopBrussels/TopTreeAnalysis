 #include "../interface/MinuitFitter.h" 

// Constructors
MinuitFitter::MinuitFitter(){
}
  
MinuitFitter::MinuitFitter(TH1F* datah, TH1F* signalh, TH1F* bkgh,double ff0, double ffl, double ffr, TNtuple *genttbarhisto[CosThetaNbins], int ndimen){

  //MinuitFitter();

  myChi2_ = Annik::MyChi2(datah, signalh, bkgh, ff0, ffl, ffr); 
  myChi2_.SetNDim(ndimen);
  myChi2_.SetNTuple(genttbarhisto);
  
  minimizer_ = new TMinuitMinimizer("Migrad");
  minimizer_->SetFunction(myChi2_);
  minimizer_->SetErrorDef(1.0);   // 1.0 for chi2, 0.5 for -logL
  minimizer_->SetVariable(0, "F0", ff0, 1.e-4);
  minimizer_->SetVariable(1, "FL", ffl, 1.e-4);
  if (ndimen==3)  minimizer_->SetVariable(2, "Normal", 1., 1.e-4); 
  
  bool isValid = minimizer_->Minimize();
  if (isValid) {
    //minimizer_->PrintResults();
    f0result_ =minimizer_->X()[0]; ef0result_= minimizer_->Errors()[0];
    flresult_= minimizer_->X()[1]; eflresult_= minimizer_->Errors()[1];
    
    //PresentationTex << PresentationOutput[TCHE+TCHP+SSVHE+SSVHP+CSV] << " & ";
    //for(int ii=0; ii< nameDataSet.size(); ii++){
    //if(ii < nameDataSet.size()-1){ 
    //PresentationTex << NumberRemainingEvents[TCHE+TCHP+SSVHE+SSVHP+CSV][ii] << " & ";
    //}
    //else if(ii == nameDataSet.size()-1 ){ 
    //PresentationTex << NumberRemainingEvents[TCHE+TCHP+SSVHE+SSVHP+CSV][ii] << " & "; //For presentation this is not the end of the table, helicity values still need to be included !!
    //PresentationTex << NumberBLeptCorrectEvents[TCHE+TCHP+SSVHE+SSVHP+CSV][ii] << " & ";
    //}
    //}
    
    if (ndimen==3){
      frresult_ = 1.-minimizer_->X()[0]-minimizer_->X()[1];
      double er0=minimizer_->Errors()[0];
      double er1=minimizer_->Errors()[1];
      double cov01 = minimizer_->CovMatrix(0,1);		  
      efrresult_ = sqrt( er0*er0 + er1*er1 + 2.*cov01);
      
      //PresentationTex <<frresult<< " & " << frresult-SMfrResult <<" & " << efrresult << " & " << flresult << " & " << flresult-SMflResult << " & " << eflresult << " & " << f0result << " & " << f0result-SMf0Result << " & " << ef0result << " \\\\ " << endl;
    } 
    else {
      frresult_ = 1.-minimizer_->X()[0]-minimizer_->X()[1];
      double er0=minimizer_->Errors()[0];
      double er1=minimizer_->Errors()[1];
      double cov01 = minimizer_->CovMatrix(0,1);		  
      efrresult_ = sqrt( er0*er0 + er1*er1 + 2.*cov01);
      
      //PresentationTex <<frresult<< " & " << frresult-SMfrResult <<" & " << efrresult << " & " << flresult << " & " << flresult-SMflResult << " & " << eflresult << " & " << f0result << " & " << f0result-SMf0Result << " & " << ef0result << " \\\\ " << endl;
    }
  }     

  //PresentationTex << " \\hline " << endl;
  
  //for (int ibinn=0; ibinn<CosThetaNbins; ibinn++){  //Waarvoor nodig???  --> Buiten klasse zetten!
  //delete genttbarhisto[ibinn];
  //}
}

MinuitFitter::~MinuitFitter()
{
  if(minimizer_) delete minimizer_;
}

