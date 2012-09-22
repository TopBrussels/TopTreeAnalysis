#include "../interface/MinuitFitterPerforming.h"

//Constructors
MinuitFitterPerforming::MinuitFitterPerforming(){}

MinuitFitterPerforming::MinuitFitterPerforming(string Name, int CosThetaBinNumber, int ndimen, int iDataSet, vector<float> CosThetaValuesKinFit[], vector<float> CosThGenKinFit[], vector<float> EventCorrectionWeightKinFit[], string dataSetName, bool semiMuon, bool semiElectron, bool SignalOnly, bool DataResults, bool JESResults, bool JERResults, bool WSystResults, bool TTScalingResults, bool TTMatchingResults, bool PUResults, bool UnclusEnergyResults, bool TopMassResults, bool TriggEvtSelResults, bool bTagResults, int dataSetsSize){

  BTagName bTagName = BTagName();  //for bTagFileOutput name giving
  TNtuple *genttbarhisto[CosThetaBinNumber];  //This is the vector of ntuples containing the generated values of cos theta* for each cos theta* reconstructed bin
  map<string,TH1F*> histo1D;

  ofstream FMinusTex((Name+"KinFit/FMinus.tex").c_str());
  ofstream FPlusTex((Name+"KinFit/FPlus.tex").c_str());
  ofstream FZeroTex((Name+"KinFit/FZero.tex").c_str());
  ofstream NormTex((Name+"KinFit/Norm.tex").c_str());
  ofstream NormBckgTex((Name+"KinFit/NormBckg.tex").c_str());

  int NumberTCHEbTags = 13;
  int NumberTCHPbTags = 13;
  int NumberSSVHEbTags = 13;
  int NumberSSVHPbTags=13;
  int NumberCSVbTags=13;
  int NumberJPbTags=13;
  int NumberJBPbTags=13;
  int TotalNumberbTags = NumberTCHEbTags + 1 + NumberTCHPbTags+1 + NumberSSVHEbTags + 1 + NumberSSVHPbTags + 1 + NumberCSVbTags + 1 + NumberJPbTags + 1 + NumberJBPbTags + 1;
  std::string PresentationOutput[TotalNumberbTags];

  if(iDataSet==0){
    FMinusTex<<" bTagger & FL result & Stat uncert & JES Syst & JER Syst & W Syst & TTScale Syst & TTMatch Syst & TTMatch Up & PU Syst & UnclusEnergy Syst & Top Mass Syst & TriggEvtSel Syst & BTagSyst & Total Syst uncert & Total uncert & Total Syst uncert (Up) & Total uncert (Up) "<<endl;
    FPlusTex<<" bTagger & FR result & Stat uncert & JES Syst & JER Syst & W Syst & TTScale Syst & TTMatch Syst & TTMatch Up & PU Syst & UnclusEnergy Syst & Top Mass Syst & TriggEvtSel Syst & BTagSyst & Total Syst uncert & Total uncert & Total Syst uncert (Up) & Total uncert (Up) "<<endl;
    FZeroTex<<" bTagger & F0 result & Stat uncert & JES Syst & JER Syst & W Syst & TTScale Syst & TTMatch Syst & TTMatch Up & PU Syst & UnclusEnergy Syst & Top Mass Syst & TriggEvtSel Syst & BTagSyst & Total Syst uncert & Total uncert & Total Syst uncert (Up) & Total uncert (Up) "<<endl;
    NormTex<<" bTagger  & Normalisation & Stat uncert & JES Syst & JER Syst & W Syst & TTScale Syst & TTMatch Syst & TTMatch Up & PU Syst & UnclusEnergy Syst & Top Mass Syst & TriggEvtSel Syst & BTagSyst & Total Syst uncert & Total uncert & Total Syst uncert (Up) & Total uncert (Up) "<< endl;
    NormBckgTex<<" bTagger  & Normalisation (Bckg) & Stat uncert & JES Syst & JER Syst & W Syst & TTScale Syst & TTMatch Syst & TTMatch Up & PU Syst & UnclusEnergy Syst & Top Mass Syst & TriggEvtSel Syst & BTagSyst & Total Syst uncert & Total uncert & Total Syst uncert (Up) & Total uncert (Up)"<< endl;
  }  

  //Initialize naming of different bTag options:
  int TCHELoop = 1;
  int TCHPLoop = 1;
  int SSVHELoop = 1;
  int SSVHPLoop = 1;
  int CSVLoop =1;
  int JPLoop = 1;
  int JBPLoop = 1;
    
  int UsedTCHE[TotalNumberbTags],UsedTCHP[TotalNumberbTags],UsedSSVHE[TotalNumberbTags],UsedSSVHP[TotalNumberbTags],UsedCSV[TotalNumberbTags], UsedJP[TotalNumberbTags], UsedJBP[TotalNumberbTags];
  while(JBPLoop <= NumberJBPbTags){
    while(JPLoop <= NumberJPbTags){
      while(CSVLoop<=NumberCSVbTags){
	while(SSVHPLoop<=NumberSSVHPbTags){  
	  while(SSVHELoop<=NumberSSVHEbTags){
	    while(TCHPLoop<=NumberTCHPbTags){
	      while(TCHELoop<=NumberTCHEbTags){
		if(iDataSet==0){
		  PresentationOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = bTagName.NameGivingPres(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags,JPLoop, NumberJPbTags,JBPLoop, NumberJBPbTags);
		  
		}
		UsedTCHE[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = TCHELoop;
		UsedTCHP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = TCHPLoop;
		UsedSSVHE[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = SSVHELoop;
		UsedSSVHP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = SSVHPLoop;
		UsedCSV[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = CSVLoop;
		UsedJP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = JPLoop;
		UsedJBP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = JBPLoop;
		
		TCHELoop++;
		if(TCHELoop == 14){TCHPLoop=2;SSVHELoop=2;SSVHPLoop=2;CSVLoop=2;JPLoop=2;JBPLoop=2;}
	      }
	      if(iDataSet==0){
		PresentationOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = bTagName.NameGivingPres(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags,JPLoop, NumberJPbTags,JBPLoop, NumberJBPbTags);
	      }
	      UsedTCHE[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = TCHELoop;
	      UsedTCHP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = TCHPLoop;
	      UsedSSVHE[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = SSVHELoop;
	      UsedSSVHP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = SSVHPLoop;
	      UsedCSV[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = CSVLoop;
	      UsedJP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = JPLoop;
	      UsedJBP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = JBPLoop;
	      
	      TCHPLoop++;
	    }
	    if(iDataSet==0){
	      PresentationOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = bTagName.NameGivingPres(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags,JPLoop, NumberJPbTags,JBPLoop, NumberJBPbTags);
	    }
	    UsedTCHE[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = TCHELoop;
	    UsedTCHP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = TCHPLoop;
	    UsedSSVHE[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = SSVHELoop;
	    UsedSSVHP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = SSVHPLoop;
	    UsedCSV[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = CSVLoop;
	    UsedJP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = JPLoop;
	    UsedJBP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = JBPLoop;
	      
	    SSVHELoop++;
	  }
	  if(iDataSet==0){
	    PresentationOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = bTagName.NameGivingPres(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags,JPLoop, NumberJPbTags,JBPLoop, NumberJBPbTags);
	  }
	  UsedTCHE[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = TCHELoop;
	  UsedTCHP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = TCHPLoop;
	  UsedSSVHE[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = SSVHELoop;
	  UsedSSVHP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = SSVHPLoop;
	  UsedCSV[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = CSVLoop;
	  UsedJP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = JPLoop;
	  UsedJBP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = JBPLoop;
	    
	  SSVHPLoop++;		
	}	
	if(iDataSet==0){
	  PresentationOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = bTagName.NameGivingPres(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags,JPLoop, NumberJPbTags,JBPLoop, NumberJBPbTags);
	}
	UsedTCHE[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = TCHELoop;
	UsedTCHP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = TCHPLoop;
	UsedSSVHE[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = SSVHELoop;
	UsedSSVHP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = SSVHPLoop;
	UsedCSV[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = CSVLoop;
	UsedJP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = JPLoop;
	UsedJBP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = JBPLoop;
	  
	CSVLoop++;
      }
      if(iDataSet==0){
	PresentationOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = bTagName.NameGivingPres(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags,JPLoop, NumberJPbTags,JBPLoop, NumberJBPbTags);
      }
      UsedTCHE[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = TCHELoop;
      UsedTCHP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = TCHPLoop;
      UsedSSVHE[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = SSVHELoop;
      UsedSSVHP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = SSVHPLoop;
      UsedCSV[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = CSVLoop;
      UsedJP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = JPLoop;
      UsedJBP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = JBPLoop;
	
      JPLoop++;
    }
    if(iDataSet==0){
      PresentationOutput[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = bTagName.NameGivingPres(TCHELoop,NumberTCHEbTags,TCHPLoop,NumberTCHPbTags,SSVHELoop,NumberSSVHEbTags,SSVHPLoop,NumberSSVHPbTags,CSVLoop,NumberCSVbTags,JPLoop, NumberJPbTags,JBPLoop, NumberJBPbTags);
    }
    UsedTCHE[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = TCHELoop;
    UsedTCHP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = TCHPLoop;
    UsedSSVHE[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = SSVHELoop;
    UsedSSVHP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = SSVHPLoop;
    UsedCSV[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = CSVLoop;
    UsedJP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = JPLoop;
    UsedJBP[TCHELoop-1+TCHPLoop-1+SSVHELoop-1+SSVHPLoop-1+CSVLoop-1+JPLoop-1+JBPLoop-1] = JBPLoop;
      
    JBPLoop++;
  }

  std::string bTag[14];
  for(int ii=1;ii<=14;ii++){
    std::stringstream out;
    out << ii;
    bTag[ii-1]  =  out.str();
  }

  float binEdge[CosThetaBinNumber+1];
  float binSize = (1.-(-1.))/CosThetaBinNumber;
  for(int ii=0; ii<=CosThetaBinNumber;ii++){
    binEdge[ii] = -1 + binSize*ii;
  }

  int SumbTag, ConsideredTagger, JBP, JP, CSV, SSVHP, SSVHE, TCHP, TCHE;
  ConsideredTagger = 0;
  for(JBP=0; JBP<NumberJBPbTags;JBP++){
    for(JP=0;JP<=NumberJPbTags;JP++){	  
      if(ConsideredTagger == 6) JP =13;   //--> Do not run over all possible TCHE bTag values!
      for(CSV =0;CSV<= NumberCSVbTags; CSV++){
	if(ConsideredTagger == 5 || ConsideredTagger ==6) CSV = 13;
	for(SSVHP=0;SSVHP <=NumberSSVHPbTags; SSVHP++){
	  if(ConsideredTagger == 4 || ConsideredTagger == 5 || ConsideredTagger ==6) SSVHP =13;   	    
	  for(SSVHE=0;SSVHE<=NumberSSVHEbTags; SSVHE++){
	    if(ConsideredTagger == 3 || ConsideredTagger == 4 || ConsideredTagger == 5 || ConsideredTagger ==6) SSVHE = 13;
	    for(TCHP=0; TCHP <= NumberTCHPbTags; TCHP++){
	      if(ConsideredTagger == 2 || ConsideredTagger == 3 || ConsideredTagger == 4 || ConsideredTagger == 5 || ConsideredTagger ==6) TCHP =13; 
	      for(TCHE=0; TCHE<= NumberTCHEbTags;TCHE++){
		if(ConsideredTagger ==1 || ConsideredTagger ==2 || ConsideredTagger ==3 || ConsideredTagger ==4 || ConsideredTagger == 5 || ConsideredTagger ==6) TCHE =13;
		SumbTag = TCHE+TCHP+SSVHE+SSVHP+CSV+JP+JBP;
		if(UsedTCHE[SumbTag]==(TCHE+1) && UsedTCHP[SumbTag]==(TCHP+1) && UsedSSVHE[SumbTag]==(SSVHE+1) && UsedSSVHP[SumbTag]==(SSVHP+1) && UsedCSV[SumbTag]==(CSV+1) && UsedJP[SumbTag] == (JP+1) && UsedJBP[SumbTag] == (JBP+1) ){
		  
		  if((dataSetName.find("Nom_TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("Nom_TTbarJets_SemiEl") ==0 && semiElectron == true)){//Defining of the genttbar histo 	 
		    char hisname[100]; 	 
		    sprintf(hisname,"CosThetaGen_TCHE%s_TCHP%s_SSVHE%s_SSVHP%s_CSV%s_JP%s_JBP%s", bTag[TCHE].c_str(),bTag[TCHP].c_str(),bTag[SSVHE].c_str(),bTag[SSVHP].c_str(),bTag[CSV].c_str(),bTag[JP].c_str(),bTag[JBP].c_str()); 	 
		    for (int ibinn=0; ibinn<CosThetaBinNumber; ibinn++){ 	 
		      genttbarhisto[ibinn]= new TNtuple(hisname,hisname,"costhgen:evtweight"); 	 
		      genttbarhisto[ibinn]->SetDirectory(0); 	 
		    } 	 
		  } 	 		 
		  
		  if(iDataSet == 0){		      
		    
		    //Nominal:
		    CosThetaString = "CosTheta_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
		    //Data:
		    CosThetaDataString = "CosThetaData_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
		    //JES :
		    CosThetaJESPlusString = "CosThetaJESPlus_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
		    CosThetaJESMinusString="CosThetaJESMinus_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
		    //JER :
		    CosThetaJERPlusString = "CosThetaJERPlus_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
		    CosThetaJERMinusString="CosThetaJERMinus_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
		    //WSyst :
		    CosThetaWPlusString = "CosThetaWPlus_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
		    CosThetaWMinusString = "CosThetaWMinus_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
		    //TTScaling :
		    CosThetaTTScalingUpString="CosThetaTTScalingUp_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
		    CosThetaTTScalingDownString = "CosThetaTTScalingDown_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
		    //TTMatching :
		    CosThetaTTMatchingUpString = "CosThetaTTMatchingUp_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
		    CosThetaTTMatchingDownString = "CosThetaTTMatchingDown_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
		    //PU :
		    CosThetaPUPlusString = "CosThetaPUPlus_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
		    CosThetaPUMinusString = "CosThetaPUMinus_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
		    //UnclusEnergy :
		    CosThetaUnclusEnergyPlusString = "CosThetaUnclusEnergyPlus_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
		    CosThetaUnclusEnergyMinusString = "CosThetaUnclusEnergyMinus_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
		    //TopMass :
		    CosThetaTopMassPlusString = "CosThetaTopMassPlus_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
		    CosThetaTopMassMinusString = "CosThetaTopMassMinus_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
		    //TriggEvtSel :
		    CosThetaTriggEvtSelPlusString = "CosThetaTriggEvtSelPlus_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
		    CosThetaTriggEvtSelMinusString = "CosThetaTriggEvtSelMinus_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
		    //BTagSyst :
		    CosThetaBTagSystPlusString = "CosThetaBTagSystPlus_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
		    CosThetaBTagSystMinusString = "CosThetaBTagSystMinus_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
	
		    CosThetaSignalString = "CosThetaSignal_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
		    //		    MlbSignalString = "MlbSignal_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
		    CosThetaBckgString = "CosThetaBckg_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];		    
		    //		    MlbBckgString = "MlbBckg_TCHE"+bTag[TCHE]+"_TCHP"+bTag[TCHP]+"_SSVHE"+bTag[SSVHE]+"_SSVHP"+bTag[SSVHP]+"_CSV"+bTag[CSV]+"_JP"+bTag[JP]+"_JBP"+bTag[JBP];
		    
		    //Nominal
		    histo1D[CosThetaString]=new TH1F(CosThetaString.c_str(),CosThetaString.c_str(),CosThetaBinNumber,-1,1);
		    histo1D[CosThetaString]->SetDirectory(0);
		    //		    histo1D[MlbString]=new TH1F(MlbString.c_str(),MlbString.c_str(),100,0,200);
		    //		    histo1D[MlbString]->SetDirectory(0);
		    //Data
		    histo1D[CosThetaDataString]=new TH1F(CosThetaDataString.c_str(),CosThetaDataString.c_str(),CosThetaBinNumber,-1,1);
		    histo1D[CosThetaDataString]->SetDirectory(0);
		    //		    histo1D[MlbDataString]=new TH1F(MlbDataString.c_str(),MlbDataString.c_str(),100,0,200);
		    //		    histo1D[MlbDataString]->SetDirectory(0);
		    //JES :
		    histo1D[CosThetaJESPlusString]=new TH1F(CosThetaJESPlusString.c_str(),CosThetaJESPlusString.c_str(),CosThetaBinNumber,-1,1);
		    histo1D[CosThetaJESPlusString]->SetDirectory(0);
		    histo1D[CosThetaJESMinusString]=new TH1F(CosThetaJESMinusString.c_str(),CosThetaJESMinusString.c_str(),CosThetaBinNumber,-1,1);
		    histo1D[CosThetaJESMinusString]->SetDirectory(0);
		    //JER :
		    histo1D[CosThetaJERPlusString]=new TH1F(CosThetaJERPlusString.c_str(),CosThetaJERPlusString.c_str(),CosThetaBinNumber,-1,1);
		    histo1D[CosThetaJERPlusString]->SetDirectory(0);
		    histo1D[CosThetaJERMinusString]=new TH1F(CosThetaJERMinusString.c_str(),CosThetaJERMinusString.c_str(),CosThetaBinNumber,-1,1);
		    histo1D[CosThetaJERMinusString]->SetDirectory(0);
		    //WJets :
		    histo1D[CosThetaWPlusString]=new TH1F(CosThetaWPlusString.c_str(),CosThetaWPlusString.c_str(),CosThetaBinNumber,-1,1);
		    histo1D[CosThetaWPlusString]->SetDirectory(0);
		    histo1D[CosThetaWMinusString]=new TH1F(CosThetaWMinusString.c_str(),CosThetaWMinusString.c_str(),CosThetaBinNumber,-1,1);
		    histo1D[CosThetaWMinusString]->SetDirectory(0);
		    //TTScaling :
		    histo1D[CosThetaTTScalingUpString]=new TH1F(CosThetaTTScalingUpString.c_str(),CosThetaTTScalingUpString.c_str(),CosThetaBinNumber,-1,1);
		    histo1D[CosThetaTTScalingUpString]->SetDirectory(0);
		    histo1D[CosThetaTTScalingDownString]=new TH1F(CosThetaTTScalingDownString.c_str(),CosThetaTTScalingDownString.c_str(),CosThetaBinNumber,-1,1);
		    histo1D[CosThetaTTScalingDownString]->SetDirectory(0);
		    //TTMatching :
		    histo1D[CosThetaTTMatchingUpString]=new TH1F(CosThetaTTMatchingUpString.c_str(),CosThetaTTMatchingUpString.c_str(),CosThetaBinNumber,-1,1);
		    histo1D[CosThetaTTMatchingUpString]->SetDirectory(0);
		    histo1D[CosThetaTTMatchingDownString]=new TH1F(CosThetaTTMatchingDownString.c_str(),CosThetaTTMatchingDownString.c_str(),CosThetaBinNumber,-1,1);
		    histo1D[CosThetaTTMatchingDownString]->SetDirectory(0);
		    //PU :
		    histo1D[CosThetaPUPlusString]=new TH1F(CosThetaPUPlusString.c_str(),CosThetaPUPlusString.c_str(),CosThetaBinNumber,-1,1);
		    histo1D[CosThetaPUPlusString]->SetDirectory(0);
		    histo1D[CosThetaPUMinusString]=new TH1F(CosThetaPUMinusString.c_str(),CosThetaPUMinusString.c_str(),CosThetaBinNumber,-1,1);
		    histo1D[CosThetaPUMinusString]->SetDirectory(0);
		    //UnclusEnergy :
		    histo1D[CosThetaUnclusEnergyPlusString]=new TH1F(CosThetaUnclusEnergyPlusString.c_str(),CosThetaUnclusEnergyPlusString.c_str(),CosThetaBinNumber,-1,1);
		    histo1D[CosThetaUnclusEnergyPlusString]->SetDirectory(0);
		    histo1D[CosThetaUnclusEnergyMinusString]=new TH1F(CosThetaUnclusEnergyMinusString.c_str(),CosThetaUnclusEnergyMinusString.c_str(),CosThetaBinNumber,-1,1);
		    histo1D[CosThetaUnclusEnergyMinusString]->SetDirectory(0);
		    //TopMass :
		    histo1D[CosThetaTopMassPlusString]=new TH1F(CosThetaTopMassPlusString.c_str(),CosThetaTopMassPlusString.c_str(),CosThetaBinNumber,-1,1);
		    histo1D[CosThetaTopMassPlusString]->SetDirectory(0);
		    histo1D[CosThetaTopMassMinusString]=new TH1F(CosThetaTopMassMinusString.c_str(),CosThetaTopMassMinusString.c_str(),CosThetaBinNumber,-1,1);
		    histo1D[CosThetaTopMassMinusString]->SetDirectory(0);
		    //TriggEvtSel :
		    histo1D[CosThetaTriggEvtSelPlusString]=new TH1F(CosThetaTriggEvtSelPlusString.c_str(),CosThetaTriggEvtSelPlusString.c_str(),CosThetaBinNumber,-1,1);
		    histo1D[CosThetaTriggEvtSelPlusString]->SetDirectory(0);
		    histo1D[CosThetaTriggEvtSelMinusString]=new TH1F(CosThetaTriggEvtSelMinusString.c_str(),CosThetaTriggEvtSelMinusString.c_str(),CosThetaBinNumber,-1,1);
		    histo1D[CosThetaTriggEvtSelMinusString]->SetDirectory(0);
		    //bTagSyst :
		    histo1D[CosThetaBTagSystPlusString]=new TH1F(CosThetaBTagSystPlusString.c_str(),CosThetaBTagSystPlusString.c_str(),CosThetaBinNumber,-1,1);
		    histo1D[CosThetaBTagSystPlusString]->SetDirectory(0);
		    histo1D[CosThetaBTagSystMinusString]=new TH1F(CosThetaBTagSystMinusString.c_str(),CosThetaBTagSystMinusString.c_str(),CosThetaBinNumber,-1,1);
		    histo1D[CosThetaBTagSystMinusString]->SetDirectory(0);
			  
		    histo1D[CosThetaSignalString]=new TH1F(CosThetaSignalString.c_str(),CosThetaSignalString.c_str(),CosThetaBinNumber,-1,1);
		    histo1D[CosThetaSignalString]->SetDirectory(0);
		    //		    histo1D[MlbSignalString]=new TH1F(MlbSignalString.c_str(),MlbSignalString.c_str(),100,0,200);
		    //		    histo1D[MlbSignalString]->SetDirectory(0);

		    histo1D[CosThetaBckgString]=new TH1F(CosThetaBckgString.c_str(),CosThetaBckgString.c_str(),CosThetaBinNumber,-1,1);
		    histo1D[CosThetaBckgString]->SetDirectory(0);
		    //		    histo1D[MlbBckgString]=new TH1F(MlbBckgString.c_str(),MlbBckgString.c_str(),100,0,200);
		    //		    histo1D[MlbBckgString]->SetDirectory(0);
		  }

		  //Filling of the histograms:
		  for(int ii=0; ii<CosThetaValuesKinFit[SumbTag].size();ii++){	

		    if((dataSetName.find("Nom_TTbarJets_SemiMu") ==0 && semiMuon == true) || (dataSetName.find("Nom_TTbarJets_SemiEl") ==0 && semiElectron == true) ){//Defining of genttbar histo
		      for(int iBin=0; iBin< CosThetaBinNumber; iBin++){//Filling of genttbar histo:		   
			if(CosThetaValuesKinFit[SumbTag][ii] >= binEdge[iBin] && CosThetaValuesKinFit[SumbTag][ii] < binEdge[iBin+1]){
			  genttbarhisto[iBin]->Fill(CosThGenKinFit[SumbTag][ii], EventCorrectionWeightKinFit[SumbTag][ii]);
			}
			else if(CosThetaValuesKinFit[SumbTag][ii] ==1){ //1 is included in last bin
			  genttbarhisto[CosThetaBinNumber-1]->Fill(CosThGenKinFit[SumbTag][ii], EventCorrectionWeightKinFit[SumbTag][ii]);
			}	      	   
		      }
		    }
		    ///////////////////////////////////////////////////////////////////////////////
		    // Change data to systematics since then nominal values will be reweighted!! //
		    ///////////////////////////////////////////////////////////////////////////////

		    //Nominal result:
		    if(dataSetName.find("Nom_") == 0){
		      histo1D[CosThetaString]->Fill(CosThetaValuesKinFit[SumbTag][ii],EventCorrectionWeightKinFit[SumbTag][ii]);
		      //		      histo1D[MlbString]->Fill(MlbValuesKinFit[SumbTag][ii],EventCorrectionWeightKinFit[SumbTag][ii]);
		    }
		    //Data result:
		    if(dataSetName.find("Data") == 0){
		      histo1D[CosThetaDataString]->Fill(CosThetaValuesKinFit[SumbTag][ii],EventCorrectionWeightKinFit[SumbTag][ii]);
		      //		      histo1D[MlbDataString]->Fill(MlbValuesKinFit[SumbTag][ii],EventCorrectionWeightKinFit[SumbTag][ii]);
		    }
		    //JES:
		    if(dataSetName.find("JESPlus") == 0)
		      histo1D[CosThetaJESPlusString]->Fill(CosThetaValuesKinFit[SumbTag][ii],EventCorrectionWeightKinFit[SumbTag][ii]);
		    if(dataSetName.find("JESMinus") == 0)
		      histo1D[CosThetaJESMinusString]->Fill(CosThetaValuesKinFit[SumbTag][ii],EventCorrectionWeightKinFit[SumbTag][ii]);
		    //JER :
		    if(dataSetName.find("JERPlus") == 0)
		      histo1D[CosThetaJERPlusString]->Fill(CosThetaValuesKinFit[SumbTag][ii],EventCorrectionWeightKinFit[SumbTag][ii]);
		    if(dataSetName.find("JERMinus") == 0)
		      histo1D[CosThetaJERMinusString]->Fill(CosThetaValuesKinFit[SumbTag][ii],EventCorrectionWeightKinFit[SumbTag][ii]);
		    //WJets :  --> Fill with WsystPlus sample and all other MC
		    if((dataSetName.find("Nom_") == 0 && dataSetName.find("Nom_WJets") != 0 && WSystResults == true) || dataSetName.find("WSystPlus") == 0 )
		      histo1D[CosThetaWPlusString]->Fill(CosThetaValuesKinFit[SumbTag][ii],EventCorrectionWeightKinFit[SumbTag][ii]);
		    if((dataSetName.find("Nom_") == 0 && dataSetName.find("Nom_WJets") != 0 && WSystResults == true) || dataSetName.find("WSystMinus") == 0 )
		      histo1D[CosThetaWMinusString]->Fill(CosThetaValuesKinFit[SumbTag][ii],EventCorrectionWeightKinFit[SumbTag][ii]);
		    //TTScaling :
		    if((dataSetName.find("Nom_") == 0 && dataSetName.find("Nom_TTbarJets") != 0 && TTScalingResults == true) || dataSetName.find("TTScalingUp") == 0)
		      histo1D[CosThetaTTScalingUpString]->Fill(CosThetaValuesKinFit[SumbTag][ii],EventCorrectionWeightKinFit[SumbTag][ii]);
		    if((dataSetName.find("Nom_") == 0 && dataSetName.find("Nom_TTbarJets") != 0 && TTScalingResults == true) || dataSetName.find("TTScalingDown") == 0)
		      histo1D[CosThetaTTScalingDownString]->Fill(CosThetaValuesKinFit[SumbTag][ii],EventCorrectionWeightKinFit[SumbTag][ii]);
		    //TTMatching :
		    if((dataSetName.find("Nom_") == 0 && dataSetName.find("Nom_TTbarJets") != 0 && TTMatchingResults == true) || dataSetName.find("TTMatchingUp") == 0 )
		      histo1D[CosThetaTTMatchingUpString]->Fill(CosThetaValuesKinFit[SumbTag][ii],EventCorrectionWeightKinFit[SumbTag][ii]);
		    if((dataSetName.find("Nom_") == 0 && dataSetName.find("Nom_TTbarJets") != 0 && TTMatchingResults == true) || dataSetName.find("TTMatchingDown") == 0 )
		      histo1D[CosThetaTTMatchingDownString]->Fill(CosThetaValuesKinFit[SumbTag][ii],EventCorrectionWeightKinFit[SumbTag][ii]);
		    //PU :
		    if(dataSetName.find("PUPlus_") == 0)
		      histo1D[CosThetaPUPlusString]->Fill(CosThetaValuesKinFit[SumbTag][ii],EventCorrectionWeightKinFit[SumbTag][ii]);
		    if(dataSetName.find("PUMinus_") == 0)
		      histo1D[CosThetaPUMinusString]->Fill(CosThetaValuesKinFit[SumbTag][ii],EventCorrectionWeightKinFit[SumbTag][ii]);
		    //UnclusEnergy :
		    if(dataSetName.find("UnclusEnergyPlus_") == 0)
		      histo1D[CosThetaUnclusEnergyPlusString]->Fill(CosThetaValuesKinFit[SumbTag][ii],EventCorrectionWeightKinFit[SumbTag][ii]);
		    if(dataSetName.find("UnclusEnergyMinus_") == 0)
		      histo1D[CosThetaUnclusEnergyMinusString]->Fill(CosThetaValuesKinFit[SumbTag][ii],EventCorrectionWeightKinFit[SumbTag][ii]);
		    //Top Mass :
		    if((dataSetName.find("Nom_") == 0 && dataSetName.find("Nom_TTbarJets") != 0 && TopMassResults == true) ||dataSetName.find("TopMassPlus") == 0)
		      histo1D[CosThetaTopMassPlusString]->Fill(CosThetaValuesKinFit[SumbTag][ii],EventCorrectionWeightKinFit[SumbTag][ii]);
		    if((dataSetName.find("Nom_") == 0 && dataSetName.find("Nom_TTbarJets") != 0 && TopMassResults == true) ||dataSetName.find("TopMassMinus") == 0)
		      histo1D[CosThetaTopMassMinusString]->Fill(CosThetaValuesKinFit[SumbTag][ii],EventCorrectionWeightKinFit[SumbTag][ii]);
		    //TriggEvtSel :
		    if(dataSetName.find("TriggEvtSelPlus_") == 0)
		      histo1D[CosThetaTriggEvtSelPlusString]->Fill(CosThetaValuesKinFit[SumbTag][ii],EventCorrectionWeightKinFit[SumbTag][ii]);
		    if(dataSetName.find("TriggEvtSelMinus_") == 0)
		      histo1D[CosThetaTriggEvtSelMinusString]->Fill(CosThetaValuesKinFit[SumbTag][ii],EventCorrectionWeightKinFit[SumbTag][ii]);
		    //bTag Syst :
		    if(dataSetName.find("bTagSystPlus_") == 0)
		      histo1D[CosThetaBTagSystPlusString]->Fill(CosThetaValuesKinFit[SumbTag][ii],EventCorrectionWeightKinFit[SumbTag][ii]);
		    if(dataSetName.find("bTagSystMinus_") == 0)
		      histo1D[CosThetaBTagSystMinusString]->Fill(CosThetaValuesKinFit[SumbTag][ii],EventCorrectionWeightKinFit[SumbTag][ii]);

		    if(SignalOnly == true && dataSetName.find("McDta") == 0){
		      histo1D[CosThetaDataString]->Fill(CosThetaValuesKinFit[SumbTag][ii],EventCorrectionWeightKinFit[SumbTag][ii]);
		      //		      histo1D[MlbDataString]->Fill(MlbValuesKinFit[SumbTag][ii],EventCorrectionWeightKinFit[SumbTag][ii]);
		    }
		    if((dataSetName.find("Nom_TTbarJets_SemiMu")==0 && semiMuon==true) || (dataSetName.find("Nom_TTbarJets_SemiEl")==0 && semiElectron == true) ){
		      histo1D[CosThetaSignalString]->Fill(CosThetaValuesKinFit[SumbTag][ii],EventCorrectionWeightKinFit[SumbTag][ii]);
		      //		      histo1D[MlbSignalString]->Fill(MlbValuesKinFit[SumbTag][ii],EventCorrectionWeightKinFit[SumbTag][ii]);
		    }
		    if(dataSetName.find("Nom_") == 0 && ((dataSetName.find("Nom_TTbarJets_SemiMu")!=0 && semiMuon==true) || (dataSetName.find("Nom_TTbarJets_SemiEl")!=0 && semiElectron == true))){
		      histo1D[CosThetaBckgString]->Fill(CosThetaValuesKinFit[SumbTag][ii],EventCorrectionWeightKinFit[SumbTag][ii]);			
		      //		      histo1D[MlbBckgString]->Fill(MlbValuesKinFit[SumbTag][ii],EventCorrectionWeightKinFit[SumbTag][ii]);
		    }
		  }
		    
		  //Perform MinuitFitter:
		  if(iDataSet==(dataSetsSize-1)){//Go in this loop when the last datasample is active to perform MinuitFitter			  
			  
		    if(CSV == 6 || CSV == 7){
		      cout << " Size of 'Data' histo for Nominal minuitFitter : " << histo1D[CosThetaString]->GetEntries() << " | " << histo1D[CosThetaString]->GetEffectiveEntries() << endl;
		      cout << " Size of 'Data' histo                          : " << histo1D[CosThetaDataString]->GetEntries() << " | " << histo1D[CosThetaDataString]->GetEffectiveEntries() << endl;
		      cout << " Size of Signal histo                          : " << histo1D[CosThetaSignalString]->GetEntries() << " | " << histo1D[CosThetaSignalString]->GetEffectiveEntries() << endl;
		      cout << " Size of Bckg histo                            : " << histo1D[CosThetaBckgString]->GetEntries() << " | " << histo1D[CosThetaBckgString]->GetEffectiveEntries() << endl;
		    }
		    
		    //		    MlbFitter MlbFitter = MlbFitter(histo1D[MlbString], histo1D[MlbSignalString], histo1D[MlbBckgString]);
		    //		    MlbFitter MlbFitterData = MlbFitter(histo1D[MlbDataString], histo1D[MlbSignalString], histo1D[MlbBckgString]);
		    //		    MlbNormTex << " \\hline \n ";
		    //		    MlbNormTex << PresentationOutput[SumbTag] << " & " << MlbFitterData.GetNorm() << " & " << MlbFitterData.GetNormError() << " & " << " & " << MlbFitterData.GetNormBckg() << " & " << MlbFitterData.GetNormBckgError() << " & " << MlbFitter.GetNorm() << " & " << MlbFitter.GetNormError() << " & " << MlbFitter.GetNormBckg() << " & " << MlbFitter.GetNormBckgError() << endl;

		    //Nominal:
		    MinuitFitter minuitFitter = MinuitFitter(histo1D[CosThetaString], histo1D[CosThetaSignalString], histo1D[CosThetaBckgString], 0.6671, 0.3325, 0.0004,genttbarhisto,ndimen);
		    //Data:
		    MinuitFitter minuitFitterData = MinuitFitter(histo1D[CosThetaDataString], histo1D[CosThetaSignalString], histo1D[CosThetaBckgString], 0.6671, 0.3325, 0.0004,genttbarhisto,ndimen);
		    //JES
		    MinuitFitter minuitFitterJESPlus = MinuitFitter(histo1D[CosThetaJESPlusString], histo1D[CosThetaSignalString],histo1D[CosThetaBckgString], 0.6671, 0.3325, 0.0004,genttbarhisto,ndimen);
		    MinuitFitter minuitFitterJESMinus = MinuitFitter(histo1D[CosThetaJESMinusString],histo1D[CosThetaSignalString],histo1D[CosThetaBckgString],0.6671, 0.3325, 0.0004,genttbarhisto,ndimen);
		    //JER
		    MinuitFitter minuitFitterJERPlus = MinuitFitter(histo1D[CosThetaJERPlusString],histo1D[CosThetaSignalString], histo1D[CosThetaBckgString], 0.6671, 0.3325, 0.0004,genttbarhisto,ndimen);
		    MinuitFitter minuitFitterJERMinus = MinuitFitter(histo1D[CosThetaJERMinusString],histo1D[CosThetaSignalString],histo1D[CosThetaBckgString], 0.6671,0.3325, 0.0004,genttbarhisto,ndimen);
		    //WSyst
		    MinuitFitter minuitFitterWPlus = MinuitFitter(histo1D[CosThetaWPlusString], histo1D[CosThetaSignalString],histo1D[CosThetaBckgString], 0.6671, 0.3325, 0.0004,genttbarhisto,ndimen);
		    MinuitFitter minuitFitterWMinus = MinuitFitter(histo1D[CosThetaWMinusString], histo1D[CosThetaSignalString], histo1D[CosThetaBckgString], 0.6671, 0.3325, 0.0004,genttbarhisto,ndimen);
		    //TTScaling
		    MinuitFitter minuitFitterTTScalingUp=MinuitFitter(histo1D[CosThetaTTScalingUpString],histo1D[CosThetaSignalString],histo1D[CosThetaBckgString],0.6671,0.3325,0.0004,genttbarhisto,ndimen);
		    MinuitFitter minuitFitterTTScalingDown = MinuitFitter(histo1D[CosThetaTTScalingDownString], histo1D[CosThetaSignalString], histo1D[CosThetaBckgString], 0.6671, 0.3325, 0.0004,genttbarhisto,ndimen);
		    //TTMatching
		    MinuitFitter minuitFitterTTMatchingUp = MinuitFitter(histo1D[CosThetaTTMatchingUpString], histo1D[CosThetaSignalString], histo1D[CosThetaBckgString], 0.6671, 0.3325, 0.0004,genttbarhisto,ndimen);
		    MinuitFitter minuitFitterTTMatchingDown = MinuitFitter(histo1D[CosThetaTTMatchingDownString], histo1D[CosThetaSignalString], histo1D[CosThetaBckgString],0.6671,0.3325, 0.0004,genttbarhisto,ndimen);
		    //PU
		    MinuitFitter minuitFitterPUPlus = MinuitFitter(histo1D[CosThetaPUPlusString], histo1D[CosThetaSignalString], histo1D[CosThetaBckgString], 0.6671, 0.3325, 0.0004,genttbarhisto,ndimen);
		    MinuitFitter minuitFitterPUMinus = MinuitFitter(histo1D[CosThetaPUMinusString], histo1D[CosThetaSignalString],histo1D[CosThetaBckgString], 0.6671, 0.3325, 0.0004,genttbarhisto,ndimen);
		    //UnclusEnergy
		    MinuitFitter minuitFitterUnclusEnergyPlus=MinuitFitter(histo1D[CosThetaUnclusEnergyPlusString],histo1D[CosThetaSignalString],histo1D[CosThetaBckgString],0.6671,0.3325,0.0004,genttbarhisto,ndimen);
		    MinuitFitter minuitFitterUnclusEnergyMinus=MinuitFitter(histo1D[CosThetaUnclusEnergyMinusString],histo1D[CosThetaSignalString],histo1D[CosThetaBckgString],0.6671,0.3325,0.0004,genttbarhisto,ndimen);
		    //TopMass
		    MinuitFitter minuitFitterTopMassPlus=MinuitFitter(histo1D[CosThetaTopMassPlusString],histo1D[CosThetaSignalString],histo1D[CosThetaBckgString],0.6671,0.3325,0.0004,genttbarhisto,ndimen);
		    MinuitFitter minuitFitterTopMassMinus=MinuitFitter(histo1D[CosThetaTopMassMinusString],histo1D[CosThetaSignalString],histo1D[CosThetaBckgString],0.6671,0.3325,0.0004,genttbarhisto,ndimen);
		    //TriggEvtSel
		    MinuitFitter minuitFitterTriggEvtSelPlus = MinuitFitter(histo1D[CosThetaTriggEvtSelPlusString],histo1D[CosThetaSignalString],histo1D[CosThetaBckgString],0.6671,0.3325,0.0004,genttbarhisto,ndimen);
		    MinuitFitter minuitFitterTriggEvtSelMinus =MinuitFitter(histo1D[CosThetaTriggEvtSelMinusString],histo1D[CosThetaSignalString],histo1D[CosThetaBckgString],0.6671,0.3325,0.0004,genttbarhisto,ndimen);
		    //BTag
		    MinuitFitter minuitFitterBTagSystPlus=MinuitFitter(histo1D[CosThetaBTagSystPlusString],histo1D[CosThetaSignalString],histo1D[CosThetaBckgString],0.6671,0.3325,0.0004,genttbarhisto,ndimen);
		    MinuitFitter minuitFitterBTagSystMinus=MinuitFitter(histo1D[CosThetaBTagSystMinusString],histo1D[CosThetaSignalString],histo1D[CosThetaBckgString],0.6671,0.3325,0.0004,genttbarhisto,ndimen);
		  
		    //FL results:
		    float FLJESSyst = 0;
		    if(JESResults == true) FLJESSyst = (abs(minuitFitter.GetFLResult() - minuitFitterJESPlus.GetFLResult()) + abs(minuitFitter.GetFLResult() - minuitFitterJESMinus.GetFLResult()))/2;
		    float FLJERSyst = 0;
		    if(JERResults == true) FLJERSyst = (abs(minuitFitter.GetFLResult() - minuitFitterJERPlus.GetFLResult()) + abs(minuitFitter.GetFLResult() - minuitFitterJERMinus.GetFLResult()))/2;
		    float FLWSyst = 0;
		    if(WSystResults == true) FLWSyst = (abs(minuitFitter.GetFLResult() - minuitFitterWPlus.GetFLResult()) + abs(minuitFitter.GetFLResult() - minuitFitterWMinus.GetFLResult()))/2;
		    float FLTTScaling = 0;
		    if(TTScalingResults == true) FLTTScaling = (abs(minuitFitter.GetFLResult() - minuitFitterTTScalingUp.GetFLResult()) + abs(minuitFitter.GetFLResult() - minuitFitterTTScalingDown.GetFLResult()))/2;
		    float FLTTMatching = 0;
		    float FLTTMatchingUp = 0;
		    if(TTMatchingResults == true){
		      FLTTMatching = (abs(minuitFitter.GetFLResult() - minuitFitterTTMatchingUp.GetFLResult()) + abs(minuitFitter.GetFLResult() - minuitFitterTTMatchingDown.GetFLResult()))/2;
		      FLTTMatchingUp = abs(minuitFitter.GetFLResult() - minuitFitterTTMatchingUp.GetFLResult());
		    }
		    float FLPUSyst = 0;
		    if(PUResults == true) FLPUSyst = (abs(minuitFitter.GetFLResult() - minuitFitterPUPlus.GetFLResult()) + abs(minuitFitter.GetFLResult() - minuitFitterPUMinus.GetFLResult()))/2;
		    float FLUnclusEnergySyst = 0;
		    if(UnclusEnergyResults == true) FLUnclusEnergySyst = (abs(minuitFitter.GetFLResult() - minuitFitterUnclusEnergyPlus.GetFLResult()) + abs(minuitFitter.GetFLResult() - minuitFitterUnclusEnergyMinus.GetFLResult()))/2;
		    float FLTopMassSyst = 0;
		    if(TopMassResults == true) FLTopMassSyst = (abs(minuitFitter.GetFLResult() - minuitFitterTopMassPlus.GetFLResult()) + abs(minuitFitter.GetFLResult() - minuitFitterTopMassMinus.GetFLResult()))/6;  //Shift to 1 GeV
		    float FLTriggEvtSelSyst = 0;
		    if(TriggEvtSelResults == true) FLTriggEvtSelSyst = (abs(minuitFitter.GetFLResult() - minuitFitterTriggEvtSelPlus.GetFLResult()) + abs(minuitFitter.GetFLResult() - minuitFitterTriggEvtSelMinus.GetFLResult()))/2;		      
		    float FLBTagSyst = 0;
		    if(bTagResults == true) FLBTagSyst = (abs(minuitFitter.GetFLResult() - minuitFitterBTagSystPlus.GetFLResult()) + abs(minuitFitter.GetFLResult() - minuitFitterBTagSystMinus.GetFLResult()))/2;		      
		    FMinusTex << "\\hline" << endl;
		    FMinusTex << PresentationOutput[SumbTag] << " & " << minuitFitterData.GetFLResult() <<" & " << minuitFitterData.GetFLError() << " & " << FLJESSyst << " & " << FLJERSyst << " & " << FLWSyst << " & " << FLTTScaling << " & " << FLTTMatching << " & " << FLTTMatchingUp << " & " << FLPUSyst << " & " << FLUnclusEnergySyst << " & " << FLTopMassSyst << " & " << FLTriggEvtSelSyst << " & " << FLBTagSyst << " & " << sqrt(FLJESSyst*FLJESSyst + FLJERSyst*FLJERSyst + FLWSyst*FLWSyst + FLTTScaling*FLTTScaling + FLTTMatching*FLTTMatching + FLPUSyst*FLPUSyst + FLUnclusEnergySyst*FLUnclusEnergySyst + FLTopMassSyst*FLTopMassSyst + FLTriggEvtSelSyst*FLTriggEvtSelSyst + FLBTagSyst*FLBTagSyst) << " & " << sqrt(FLJESSyst*FLJESSyst + FLJERSyst*FLJERSyst + FLWSyst*FLWSyst + FLTTScaling*FLTTScaling + FLTTMatching*FLTTMatching + FLPUSyst*FLPUSyst + FLUnclusEnergySyst*FLUnclusEnergySyst + FLTopMassSyst*FLTopMassSyst + FLTriggEvtSelSyst*FLTriggEvtSelSyst + FLBTagSyst*FLBTagSyst + (minuitFitterData.GetFLError())*(minuitFitterData.GetFLError())) << " & " << sqrt(FLJESSyst*FLJESSyst + FLJERSyst*FLJERSyst + FLWSyst*FLWSyst + FLTTScaling*FLTTScaling + FLTTMatchingUp*FLTTMatchingUp + FLPUSyst*FLPUSyst + FLUnclusEnergySyst*FLUnclusEnergySyst + FLTopMassSyst*FLTopMassSyst + FLTriggEvtSelSyst*FLTriggEvtSelSyst + FLBTagSyst*FLBTagSyst) << " & " << sqrt(FLJESSyst*FLJESSyst + FLJERSyst*FLJERSyst + FLWSyst*FLWSyst + FLTTScaling*FLTTScaling + FLTTMatchingUp*FLTTMatchingUp + FLPUSyst*FLPUSyst + FLUnclusEnergySyst*FLUnclusEnergySyst + FLTopMassSyst*FLTopMassSyst + FLTriggEvtSelSyst*FLTriggEvtSelSyst + FLBTagSyst*FLBTagSyst + (minuitFitterData.GetFLError())*(minuitFitterData.GetFLError())) << endl;
		  
		    //FR results:
		    float FRJESSyst = 0;
		    if(JESResults == true) FRJESSyst = (abs(minuitFitter.GetFRResult() - minuitFitterJESPlus.GetFRResult()) + abs(minuitFitter.GetFRResult() - minuitFitterJESMinus.GetFRResult()))/2;
		    float FRJERSyst = 0;
		    if(JERResults == true) FRJERSyst = (abs(minuitFitter.GetFRResult() - minuitFitterJERPlus.GetFRResult()) + abs(minuitFitter.GetFRResult() - minuitFitterJERMinus.GetFRResult()))/2;
		    float FRWSyst = 0;
		    if(WSystResults == true) FRWSyst = (abs(minuitFitter.GetFRResult() - minuitFitterWPlus.GetFRResult()) + abs(minuitFitter.GetFRResult() - minuitFitterWMinus.GetFRResult()))/2;
		    float FRTTScaling = 0;
		    if(TTScalingResults == true) FRTTScaling = (abs(minuitFitter.GetFRResult() - minuitFitterTTScalingUp.GetFRResult()) + abs(minuitFitter.GetFRResult() - minuitFitterTTScalingDown.GetFRResult()))/2;
		    float FRTTMatching = 0;
		    float FRTTMatchingUp = 0;
		    if(TTMatchingResults == true){
		      FRTTMatching = (abs(minuitFitter.GetFRResult() - minuitFitterTTMatchingUp.GetFRResult()) + abs(minuitFitter.GetFRResult() - minuitFitterTTMatchingDown.GetFRResult()))/2;		  
		      FRTTMatchingUp = abs(minuitFitter.GetFRResult() - minuitFitterTTMatchingUp.GetFRResult());		  
		    }
		    float FRPUSyst = 0;
		    if(PUResults == true) FRPUSyst = (abs(minuitFitter.GetFRResult() - minuitFitterPUPlus.GetFRResult()) + abs(minuitFitter.GetFRResult() - minuitFitterPUMinus.GetFRResult()))/2;
		    float FRUnclusEnergySyst = 0;
		    if(UnclusEnergyResults == true) FRUnclusEnergySyst = (abs(minuitFitter.GetFRResult() - minuitFitterUnclusEnergyPlus.GetFRResult()) + abs(minuitFitter.GetFRResult() - minuitFitterUnclusEnergyMinus.GetFRResult()))/2;
		    float FRTopMassSyst = 0;
		    if(TopMassResults == true) FRTopMassSyst = (abs(minuitFitter.GetFRResult() - minuitFitterTopMassPlus.GetFRResult()) + abs(minuitFitter.GetFRResult() - minuitFitterTopMassMinus.GetFRResult()))/6;  //Shift to 1 GeV
		    float FRTriggEvtSelSyst = 0;
		    if(TriggEvtSelResults == true) FRTriggEvtSelSyst = (abs(minuitFitter.GetFRResult() - minuitFitterTriggEvtSelPlus.GetFRResult()) + abs(minuitFitter.GetFRResult() - minuitFitterTriggEvtSelMinus.GetFRResult()))/2;
		    float FRBTagSyst = 0;
		    if(bTagResults == true) FRBTagSyst = (abs(minuitFitter.GetFRResult() - minuitFitterBTagSystPlus.GetFRResult()) + abs(minuitFitter.GetFRResult() - minuitFitterBTagSystMinus.GetFRResult()))/2;		      
		    FPlusTex << "\\hline" << endl;
		    FPlusTex << PresentationOutput[SumbTag] << " & " << minuitFitterData.GetFRResult() <<" & " << minuitFitterData.GetFRError() << " & " << FRJESSyst << " & " << FRJERSyst << " & " << FRWSyst << " & " << FRTTScaling << " & " << FRTTMatching << " & " << FRTTMatchingUp << " & " << FRPUSyst << " & " << FRUnclusEnergySyst << " & " << FRTopMassSyst << " & " << FRTriggEvtSelSyst << " & " << FRBTagSyst << " & " << sqrt(FRJESSyst*FRJESSyst + FRJERSyst*FRJERSyst + FRWSyst*FRWSyst + FRTTScaling*FRTTScaling + FRTTMatching*FRTTMatching + FRPUSyst*FRPUSyst + FRUnclusEnergySyst*FRUnclusEnergySyst + FRTopMassSyst*FRTopMassSyst + FRTriggEvtSelSyst*FRTriggEvtSelSyst + FRBTagSyst*FRBTagSyst) << " & " << sqrt(FRJESSyst*FRJESSyst + FRJERSyst*FRJERSyst + FRWSyst*FRWSyst + FRTTScaling*FRTTScaling + FRTTMatching*FRTTMatching + FRPUSyst*FRPUSyst + FRUnclusEnergySyst*FRUnclusEnergySyst + FRTopMassSyst*FRTopMassSyst + FRTriggEvtSelSyst*FRTriggEvtSelSyst + FRBTagSyst*FRBTagSyst + (minuitFitterData.GetFRError())*(minuitFitterData.GetFRError())) << " & " << sqrt(FRJESSyst*FRJESSyst + FRJERSyst*FRJERSyst + FRWSyst*FRWSyst + FRTTScaling*FRTTScaling + FRTTMatchingUp*FRTTMatchingUp + FRPUSyst*FRPUSyst + FRUnclusEnergySyst*FRUnclusEnergySyst + FRTopMassSyst*FRTopMassSyst + FRTriggEvtSelSyst*FRTriggEvtSelSyst + FRBTagSyst*FRBTagSyst) << " & " << sqrt(FRJESSyst*FRJESSyst + FRJERSyst*FRJERSyst + FRWSyst*FRWSyst + FRTTScaling*FRTTScaling + FRTTMatchingUp*FRTTMatchingUp + FRPUSyst*FRPUSyst + FRUnclusEnergySyst*FRUnclusEnergySyst + FRTopMassSyst*FRTopMassSyst + FRTriggEvtSelSyst*FRTriggEvtSelSyst + FRBTagSyst*FRBTagSyst + (minuitFitterData.GetFRError())*(minuitFitterData.GetFRError())) << endl;

		    //F0 results:
		    float F0JESSyst = 0;
		    if(JESResults == true) F0JESSyst = (abs(minuitFitter.GetF0Result() - minuitFitterJESPlus.GetF0Result()) + abs(minuitFitter.GetF0Result() - minuitFitterJESMinus.GetF0Result()))/2;
		    float F0JERSyst = 0;
		    if(JERResults == true) F0JERSyst = (abs(minuitFitter.GetF0Result() - minuitFitterJERPlus.GetF0Result()) + abs(minuitFitter.GetF0Result() - minuitFitterJERMinus.GetF0Result()))/2;
		    float F0WSyst = 0;
		    if(WSystResults == true) F0WSyst = (abs(minuitFitter.GetF0Result() - minuitFitterWPlus.GetF0Result()) + abs(minuitFitter.GetF0Result() - minuitFitterWMinus.GetF0Result()))/2;
		    float F0TTScaling = 0;
		    if(TTScalingResults == true) F0TTScaling = (abs(minuitFitter.GetF0Result() - minuitFitterTTScalingUp.GetF0Result()) + abs(minuitFitter.GetF0Result() - minuitFitterTTScalingDown.GetF0Result()))/2;
		    float F0TTMatching = 0;
		    float F0TTMatchingUp = 0;
		    if(TTMatchingResults == true){
		      F0TTMatching = (abs(minuitFitter.GetF0Result() - minuitFitterTTMatchingUp.GetF0Result()) + abs(minuitFitter.GetF0Result() - minuitFitterTTMatchingDown.GetF0Result()))/2;		  
		      F0TTMatchingUp = abs(minuitFitter.GetF0Result() - minuitFitterTTMatchingUp.GetF0Result());		  
		    }
		    float F0PUSyst = 0;
		    if(PUResults == true) F0PUSyst = (abs(minuitFitter.GetF0Result() - minuitFitterPUPlus.GetF0Result()) + abs(minuitFitter.GetF0Result() - minuitFitterPUMinus.GetF0Result()))/2;
		    float F0UnclusEnergySyst = 0;
		    if(UnclusEnergyResults == true) F0UnclusEnergySyst = (abs(minuitFitter.GetF0Result() - minuitFitterUnclusEnergyPlus.GetF0Result()) + abs(minuitFitter.GetF0Result() - minuitFitterUnclusEnergyMinus.GetF0Result()))/2;
		    float F0TopMassSyst = 0;
		    if(TopMassResults == true) F0TopMassSyst = (abs(minuitFitter.GetF0Result() - minuitFitterTopMassPlus.GetF0Result()) + abs(minuitFitter.GetF0Result() - minuitFitterTopMassMinus.GetF0Result()))/6;  //Shift to 1 GeV
		    float F0TriggEvtSelSyst = 0;
		    if(TriggEvtSelResults == true) F0TriggEvtSelSyst = (abs(minuitFitter.GetF0Result() - minuitFitterTriggEvtSelPlus.GetF0Result()) + abs(minuitFitter.GetF0Result() - minuitFitterTriggEvtSelMinus.GetF0Result()))/2;
		    float F0BTagSyst = 0;
		    if(bTagResults == true) F0BTagSyst = (abs(minuitFitter.GetF0Result() - minuitFitterBTagSystPlus.GetF0Result()) + abs(minuitFitter.GetF0Result() - minuitFitterBTagSystMinus.GetF0Result()))/2;		      
		    FZeroTex << " \\hline" << endl;
		    FZeroTex << PresentationOutput[SumbTag] << " & " << minuitFitterData.GetF0Result() <<" & " << minuitFitterData.GetF0Error() << " & " << F0JESSyst << " & " << F0JERSyst << " & " << F0WSyst << " & " << F0TTScaling << " & " << F0TTMatching << " & " << F0TTMatchingUp << " & " << F0PUSyst << " & " << F0UnclusEnergySyst << " & " << F0TopMassSyst << " & " << F0TriggEvtSelSyst << " & " << F0BTagSyst << " & " << sqrt(F0JESSyst*F0JESSyst + F0JERSyst*F0JERSyst + F0WSyst*F0WSyst + F0TTScaling*F0TTScaling + F0TTMatching*F0TTMatching + F0PUSyst*F0PUSyst + F0UnclusEnergySyst*F0UnclusEnergySyst + F0TopMassSyst*F0TopMassSyst + F0TriggEvtSelSyst*F0TriggEvtSelSyst + F0BTagSyst*F0BTagSyst) << " & " << sqrt(F0JESSyst*F0JESSyst + F0JERSyst*F0JERSyst + F0WSyst*F0WSyst + F0TTScaling*F0TTScaling + F0TTMatching*F0TTMatching + F0PUSyst*F0PUSyst + F0UnclusEnergySyst*F0UnclusEnergySyst + F0TopMassSyst*F0TopMassSyst + F0TriggEvtSelSyst*F0TriggEvtSelSyst + F0BTagSyst*F0BTagSyst + (minuitFitterData.GetF0Error())*(minuitFitterData.GetF0Error())) << " & " << sqrt(F0JESSyst*F0JESSyst + F0JERSyst*F0JERSyst + F0WSyst*F0WSyst + F0TTScaling*F0TTScaling + F0TTMatchingUp*F0TTMatchingUp + F0PUSyst*F0PUSyst + F0UnclusEnergySyst*F0UnclusEnergySyst + F0TopMassSyst*F0TopMassSyst + F0TriggEvtSelSyst*F0TriggEvtSelSyst + F0BTagSyst*F0BTagSyst) << " & " << sqrt(F0JESSyst*F0JESSyst + F0JERSyst*F0JERSyst + F0WSyst*F0WSyst + F0TTScaling*F0TTScaling + F0TTMatchingUp*F0TTMatchingUp + F0PUSyst*F0PUSyst + F0UnclusEnergySyst*F0UnclusEnergySyst + F0TopMassSyst*F0TopMassSyst + F0TriggEvtSelSyst*F0TriggEvtSelSyst + F0BTagSyst*F0BTagSyst + (minuitFitterData.GetF0Error())*(minuitFitterData.GetF0Error())) << endl;

		    //Normalisation results:
		    float NormJESSyst = 0;
		    if(JESResults == true) NormJESSyst = (abs(minuitFitter.GetNorm() - minuitFitterJESPlus.GetNorm()) + abs(minuitFitter.GetNorm() - minuitFitterJESMinus.GetNorm()))/2;
		    float NormJERSyst = 0;
		    if(JERResults == true) NormJERSyst = (abs(minuitFitter.GetNorm() - minuitFitterJERPlus.GetNorm()) + abs(minuitFitter.GetNorm() - minuitFitterJERMinus.GetNorm()))/2;
		    float NormWSyst = 0;
		    if(WSystResults == true) NormWSyst = (abs(minuitFitter.GetNorm() - minuitFitterWPlus.GetNorm()) + abs(minuitFitter.GetNorm() - minuitFitterWMinus.GetNorm()))/2;
		    float NormTTScaling = 0;
		    if(TTScalingResults == true) NormTTScaling = (abs(minuitFitter.GetNorm() - minuitFitterTTScalingUp.GetNorm()) + abs(minuitFitter.GetNorm() - minuitFitterTTScalingDown.GetNorm()))/2;
		    float NormTTMatching = 0;
		    float NormTTMatchingUp = 0;
		    if(TTMatchingResults == true){
		      NormTTMatching = (abs(minuitFitter.GetNorm() - minuitFitterTTMatchingUp.GetNorm()) + abs(minuitFitter.GetNorm() - minuitFitterTTMatchingDown.GetNorm()))/2;
		      NormTTMatchingUp = abs(minuitFitter.GetNorm() - minuitFitterTTMatchingUp.GetNorm());
		    }
		    float NormPUSyst = 0;
		    if(PUResults == true) NormPUSyst = (abs(minuitFitter.GetNorm() - minuitFitterPUPlus.GetNorm()) + abs(minuitFitter.GetNorm() - minuitFitterPUMinus.GetNorm()))/2;
		    float NormUnclusEnergySyst = 0;
		    if(UnclusEnergyResults == true) NormUnclusEnergySyst = (abs(minuitFitter.GetNorm() - minuitFitterUnclusEnergyPlus.GetNorm()) + abs(minuitFitter.GetNorm() - minuitFitterUnclusEnergyMinus.GetNorm()))/2;
		    float NormTopMassSyst = 0;
		    if(TopMassResults==true) NormTopMassSyst=(abs(minuitFitter.GetNorm()-minuitFitterTopMassPlus.GetNorm())+abs(minuitFitter.GetNorm()-minuitFitterTopMassMinus.GetNorm()))/6;//Shift to 1 GeV
		    float NormTriggEvtSelSyst = 0;
		    if(TriggEvtSelResults == true) NormTriggEvtSelSyst = (abs(minuitFitter.GetNorm() - minuitFitterTriggEvtSelPlus.GetNorm()) + abs(minuitFitter.GetNorm() - minuitFitterTriggEvtSelMinus.GetNorm()))/2;
		    float NormBTagSyst = 0;
		    if(bTagResults == true) NormBTagSyst = (abs(minuitFitter.GetNorm() - minuitFitterBTagSystPlus.GetNorm()) + abs(minuitFitter.GetNorm() - minuitFitterBTagSystMinus.GetNorm()))/2;

		    if(CSV == 6 || CSV == 7){
		      cout << " Norm result (Data) : " << minuitFitterData.GetNorm() << endl;
		      cout << " Norm result (Nominal) : " << minuitFitter.GetNorm() << endl;
		      cout << " Configuration : " << PresentationOutput[SumbTag] << endl;
		    }

		    NormTex << " \\hline" << endl;
		    NormTex << PresentationOutput[SumbTag] << " & " << minuitFitterData.GetNorm() <<" & " << minuitFitterData.GetNormError() << " & " << NormJESSyst << " & " << NormJERSyst << " & " << NormWSyst << " & " << NormTTScaling << " & " << NormTTMatching << " & " << NormTTMatchingUp << " & " << NormPUSyst << " & " << NormUnclusEnergySyst << " & " << NormTopMassSyst << " & " << NormTriggEvtSelSyst << " & " << NormBTagSyst << " & " << sqrt(NormJESSyst*NormJESSyst + NormJERSyst*NormJERSyst + NormWSyst*NormWSyst + NormTTScaling*NormTTScaling + NormTTMatching*NormTTMatching + NormPUSyst*NormPUSyst + NormUnclusEnergySyst*NormUnclusEnergySyst + NormTopMassSyst*NormTopMassSyst + NormTriggEvtSelSyst*NormTriggEvtSelSyst + NormBTagSyst*NormBTagSyst) << " & " << sqrt(NormJESSyst*NormJESSyst + NormJERSyst*NormJERSyst + NormWSyst*NormWSyst + NormTTScaling*NormTTScaling + NormTTMatching*NormTTMatching + NormPUSyst*NormPUSyst + NormUnclusEnergySyst*NormUnclusEnergySyst + NormTopMassSyst*NormTopMassSyst + NormTriggEvtSelSyst*NormTriggEvtSelSyst + NormBTagSyst*NormBTagSyst + (minuitFitterData.GetNormError())*(minuitFitterData.GetNormError())) << " & " << sqrt(NormJESSyst*NormJESSyst + NormJERSyst*NormJERSyst + NormWSyst*NormWSyst + NormTTScaling*NormTTScaling + NormTTMatchingUp*NormTTMatchingUp + NormPUSyst*NormPUSyst + NormUnclusEnergySyst*NormUnclusEnergySyst + NormTopMassSyst*NormTopMassSyst + NormTriggEvtSelSyst*NormTriggEvtSelSyst + NormBTagSyst*NormBTagSyst) << " & " << sqrt(NormJESSyst*NormJESSyst + NormJERSyst*NormJERSyst + NormWSyst*NormWSyst + NormTTScaling*NormTTScaling + NormTTMatchingUp*NormTTMatchingUp + NormPUSyst*NormPUSyst + NormUnclusEnergySyst*NormUnclusEnergySyst + NormTopMassSyst*NormTopMassSyst + NormTriggEvtSelSyst*NormTriggEvtSelSyst + NormBTagSyst*NormBTagSyst + (minuitFitterData.GetNormError())*(minuitFitterData.GetNormError())) << endl;

		    //Normalisation results:
		    float NormBckgJESSyst = 0;
		    if(JESResults == true) NormBckgJESSyst = (abs(minuitFitter.GetNormBckg() - minuitFitterJESPlus.GetNormBckg()) + abs(minuitFitter.GetNormBckg() - minuitFitterJESMinus.GetNormBckg()))/2;
		    float NormBckgJERSyst = 0;
		    if(JERResults == true) NormBckgJERSyst = (abs(minuitFitter.GetNormBckg() - minuitFitterJERPlus.GetNormBckg()) + abs(minuitFitter.GetNormBckg() - minuitFitterJERMinus.GetNormBckg()))/2;
		    float NormBckgWSyst = 0;
		    if(WSystResults == true) NormBckgWSyst = (abs(minuitFitter.GetNormBckg() - minuitFitterWPlus.GetNormBckg()) + abs(minuitFitter.GetNormBckg() - minuitFitterWMinus.GetNormBckg()))/2;
		    float NormBckgTTScaling = 0;
		    if(TTScalingResults == true) NormBckgTTScaling = (abs(minuitFitter.GetNormBckg() - minuitFitterTTScalingUp.GetNormBckg()) + abs(minuitFitter.GetNormBckg() - minuitFitterTTScalingDown.GetNormBckg()))/2;
		    float NormBckgTTMatching = 0;
		    float NormBckgTTMatchingUp = 0;
		    if(TTMatchingResults == true){
		      NormBckgTTMatching = (abs(minuitFitter.GetNormBckg() - minuitFitterTTMatchingUp.GetNormBckg()) + abs(minuitFitter.GetNormBckg() - minuitFitterTTMatchingDown.GetNormBckg()))/2;
		      NormBckgTTMatchingUp = abs(minuitFitter.GetNormBckg() - minuitFitterTTMatchingUp.GetNormBckg());
		    }
		    float NormBckgPUSyst = 0;
		    if(PUResults == true) NormBckgPUSyst = (abs(minuitFitter.GetNormBckg() - minuitFitterPUPlus.GetNormBckg()) + abs(minuitFitter.GetNormBckg() - minuitFitterPUMinus.GetNormBckg()))/2;
		    float NormBckgUnclusEnergySyst = 0;
		    if(UnclusEnergyResults == true) NormBckgUnclusEnergySyst = (abs(minuitFitter.GetNormBckg() - minuitFitterUnclusEnergyPlus.GetNormBckg()) + abs(minuitFitter.GetNormBckg() - minuitFitterUnclusEnergyMinus.GetNormBckg()))/2;
		    float NormBckgTopMassSyst = 0;
		    if(TopMassResults==true) NormBckgTopMassSyst=(abs(minuitFitter.GetNormBckg()-minuitFitterTopMassPlus.GetNormBckg())+abs(minuitFitter.GetNormBckg()-minuitFitterTopMassMinus.GetNormBckg()))/6;//Shift to 1 GeV
		    float NormBckgTriggEvtSelSyst = 0;
		    if(TriggEvtSelResults == true) NormBckgTriggEvtSelSyst = (abs(minuitFitter.GetNormBckg() - minuitFitterTriggEvtSelPlus.GetNormBckg()) + abs(minuitFitter.GetNormBckg() - minuitFitterTriggEvtSelMinus.GetNormBckg()))/2;
		    float NormBckgBTagSyst = 0;
		    if(bTagResults == true) NormBckgBTagSyst = (abs(minuitFitter.GetNormBckg() - minuitFitterBTagSystPlus.GetNormBckg()) + abs(minuitFitter.GetNormBckg() - minuitFitterBTagSystMinus.GetNormBckg()))/2;

		    if(CSV == 6 || CSV == 7){
		      cout << " NormBckg result (Data) : " << minuitFitterData.GetNormBckg() << endl;
		      cout << " NormBckg result (Nominal) : " << minuitFitter.GetNormBckg() << endl;
		      cout << " Configuration : " << PresentationOutput[SumbTag] << endl;
		    }

		    NormBckgTex << " \\hline" << endl;
		    NormBckgTex << PresentationOutput[SumbTag] << " & " << minuitFitterData.GetNormBckg() <<" & " << minuitFitterData.GetNormBckgError() << " & " << NormBckgJESSyst << " & " << NormBckgJERSyst << " & " << NormBckgWSyst << " & " << NormBckgTTScaling << " & " << NormBckgTTMatching << " & " << NormBckgTTMatchingUp << " & " << NormBckgPUSyst << " & " << NormBckgUnclusEnergySyst << " & " << NormBckgTopMassSyst << " & " << NormBckgTriggEvtSelSyst << " & " << NormBckgBTagSyst << " & " << sqrt(NormBckgJESSyst*NormBckgJESSyst + NormBckgJERSyst*NormBckgJERSyst + NormBckgWSyst*NormBckgWSyst + NormBckgTTScaling*NormBckgTTScaling + NormBckgTTMatching*NormBckgTTMatching + NormBckgPUSyst*NormBckgPUSyst + NormBckgUnclusEnergySyst*NormBckgUnclusEnergySyst + NormBckgTopMassSyst*NormBckgTopMassSyst + NormBckgTriggEvtSelSyst*NormBckgTriggEvtSelSyst + NormBckgBTagSyst*NormBckgBTagSyst) << " & " << sqrt(NormBckgJESSyst*NormBckgJESSyst + NormBckgJERSyst*NormBckgJERSyst + NormBckgWSyst*NormBckgWSyst + NormBckgTTScaling*NormBckgTTScaling + NormBckgTTMatching*NormBckgTTMatching + NormBckgPUSyst*NormBckgPUSyst + NormBckgUnclusEnergySyst*NormBckgUnclusEnergySyst + NormBckgTopMassSyst*NormBckgTopMassSyst + NormBckgTriggEvtSelSyst*NormBckgTriggEvtSelSyst + NormBckgBTagSyst*NormBckgBTagSyst + (minuitFitterData.GetNormBckgError())*(minuitFitterData.GetNormBckgError())) << " & " << sqrt(NormBckgJESSyst*NormBckgJESSyst + NormBckgJERSyst*NormBckgJERSyst + NormBckgWSyst*NormBckgWSyst + NormBckgTTScaling*NormBckgTTScaling + NormBckgTTMatchingUp*NormBckgTTMatchingUp + NormBckgPUSyst*NormBckgPUSyst + NormBckgUnclusEnergySyst*NormBckgUnclusEnergySyst + NormBckgTopMassSyst*NormBckgTopMassSyst + NormBckgTriggEvtSelSyst*NormBckgTriggEvtSelSyst + NormBckgBTagSyst*NormBckgBTagSyst) << " & " << sqrt(NormBckgJESSyst*NormBckgJESSyst + NormBckgJERSyst*NormBckgJERSyst + NormBckgWSyst*NormBckgWSyst + NormBckgTTScaling*NormBckgTTScaling + NormBckgTTMatchingUp*NormBckgTTMatchingUp + NormBckgPUSyst*NormBckgPUSyst + NormBckgUnclusEnergySyst*NormBckgUnclusEnergySyst + NormBckgTopMassSyst*NormBckgTopMassSyst + NormBckgTriggEvtSelSyst*NormBckgTriggEvtSelSyst + NormBckgBTagSyst*NormBckgBTagSyst + (minuitFitterData.GetNormBckgError())*(minuitFitterData.GetNormBckgError())) << endl;
			  

		  }		

		}//end of wrong entry chosen failed requirement!
		else{
		  cout << " Looking at wrong bTagging combination in fitting minuit !! " << endl;
		  cout << " Looking at : TCHE = " << TCHE+1<< " TCHP = " << TCHP+1 << " SSVHE = " << SSVHE+1 << " SSVHP = " << SSVHP+1 << " CSV = " << CSV+1 << " JP = " << JP+1 << " JBP = " << JBP+1<< endl;
		  cout << " Correct combination : TCHE = "<< UsedTCHE[SumbTag] << " TCHP = " << UsedTCHP[SumbTag] << " SSVHE = " << UsedSSVHE[SumbTag] << " SSVHP = " << UsedSSVHP[SumbTag] << " CSV = " << UsedCSV[SumbTag] << " JP = " << UsedJP[SumbTag] << " JBP = " << UsedJBP[SumbTag] << endl;
		}
		
		if(TCHE==NumberTCHEbTags-1) {
		  ConsideredTagger =1;
		  TCHP=1;
		  SSVHE=1;
		  SSVHP=1;
		  CSV=1;
		  JP=1;
		  JBP=1;
		}		    
	      }//end of TCHE
	      if(TCHP == NumberTCHPbTags-1) ConsideredTagger =2;
	    }//end of TCHP
	    if(SSVHE == NumberSSVHEbTags-1) ConsideredTagger = 3;
	  }//end of SSVHE
	  if(SSVHP == NumberSSVHPbTags-1) ConsideredTagger = 4;
	}//end of SSVHP
	if(CSV == NumberCSVbTags-1) ConsideredTagger = 5;
      }//end of CSV
      if(JP == NumberJPbTags-1) ConsideredTagger=6;
    }//end of JP
  }//end of JBP      
  
  //} //End of boolean request: PerformMinuit
    
  //Closing of the texfiles!
  FMinusTex.close();
  FPlusTex.close();
  FZeroTex.close();
  NormTex.close();
  NormBckgTex.close();
}
