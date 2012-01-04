#include "../interface/WReconstruction.h"

//-------- Constructor ----------//
WReconstruction::WReconstruction(){

  WMassKinFit = 80.4;

  //  resFitLightJets_ = new ResolutionFit("LightJet");
  //  resFitLightJets_->LoadResolutions("lightJetReso.root");
  //  resFitBJets_ = new ResolutionFit("BJet");
  //  resFitBJets_->LoadResolutions("bJetReso.root");  

  WMassKinFit = 80.4;
  MassW=83.6103;
  MassTop = 172.956;
  SigmaW=11.1534;  
  SigmaTop=18.232;
}

//-------- Destructor ----------//
WReconstruction::~WReconstruction(){

}

//-------- Initializing variables ------------//
void WReconstruction::Initializing(bool Simulation,std::vector<float> jetsPx,std::vector<float> jetsPy,std::vector<float> jetsPz,std::vector<float> jetsE, std::vector<float> TCHEbTagValues,std::vector<float> SSVHEbTagValues,std::vector<float> EtResolutionLightJets,std::vector<float> ThetaResolutionLightJets,std::vector<float> PhiResolutionLightJets,std::vector<float> EtResolutionBJets,std::vector<float> ThetaResolutionBJets,std::vector<float> PhiResolutionBJets){

  jetsPx_.clear();
  jetsPy_.clear();
  jetsPz_.clear();
  jetsE_.clear();
  TCHEbTags_.clear();
  SSVHEbTags_.clear();
  JetNumbers.clear();
  bTagBoolean.clear();
  EtResLight_.clear();
  ThetaResLight_.clear();
  PhiResLight_.clear();
  EtResB_.clear();
  ThetaResB_.clear();
  PhiResB_.clear();

  NumberCombinations = 0;
  BTagValueMVA=0;
  ChiSquaredFit = 100000;
  JetSize=0;
  TCHEbTagSize=0;
  SSVHEbTagSize=0;
  bTagSelected=0;

  if(!Simulation){TopMassKinFit = 173.1;} //Top mass value in data
  else{TopMassKinFit = 172.5;} //Top mass value for simulation samples

  if(jetsPx.size() == jetsPy.size() && jetsPx.size() == jetsPz.size() && jetsPx.size() == jetsE.size() && jetsPx.size() == EtResolutionLightJets.size() && jetsPx.size() == ThetaResolutionLightJets.size() && jetsPx.size() == PhiResolutionLightJets.size()&& jetsPx.size() == EtResolutionBJets.size()&& jetsPx.size() == ThetaResolutionBJets.size() && jetsPx.size() == PhiResolutionBJets.size() && jetsPx.size() == TCHEbTagValues.size() && jetsPx.size() == SSVHEbTagValues.size()){

    for(unsigned int ii=0;ii<jetsPx.size();ii++){
      jetsPx_.push_back(jetsPx[ii]);
      jetsPy_.push_back(jetsPy[ii]);
      jetsPz_.push_back(jetsPz[ii]);
      jetsE_.push_back(jetsE[ii]);
      EtResLight_.push_back(EtResolutionLightJets[ii]);
      ThetaResLight_.push_back(ThetaResolutionLightJets[ii]);
      PhiResLight_.push_back(PhiResolutionLightJets[ii]);
      EtResB_.push_back(EtResolutionBJets[ii]);
      ThetaResB_.push_back(ThetaResolutionBJets[ii]);
      PhiResB_.push_back(PhiResolutionBJets[ii]);      
      TCHEbTags_.push_back(TCHEbTagValues[ii]);
      SSVHEbTags_.push_back(SSVHEbTagValues[ii]);
    }
  }
  else{cerr << " Sizes of jet information (px, py, pz, E, resolutions or bTags) are not equal !! " << std::endl;}

}

//-------- Control checks for correct working of class ----------//
bool WReconstruction::ControlChecks(){
  JetSize = jetsPx_.size();
  TCHEbTagSize = TCHEbTags_.size();
  SSVHEbTagSize = SSVHEbTags_.size();

  //Check if array sizes of jets corresponds to array sizes of the bTag values associated to the jets
  if(JetSize != TCHEbTagSize){
    std::cerr << " Vector of jets and TCHE bTag values does not have the same length !! " << std::endl;
    return false;
  }
  if(JetSize != SSVHEbTagSize){
    std::cerr << " Vector of jets and SSVHE bTag values does not have the same lengt !! " << std::endl;
    return false;
  }
  
  //Check if array size length is at least larger than four
  if(JetSize<4){
    std::cerr << " Length of vector of jets is smaller than four!! " << std::endl;
    return false;
  }

  return true;
}

//-------- Calculate jet topology with bTag analysis before application of Kinematic Fit ---------//
std::vector<int> WReconstruction::BTagAnalysis(){
  
  if(!ControlChecks()){
    std::cerr << " Error discovered in ControlChecks() !! " << std::endl;
  }

  //Check for all twelve possible combinations which two jets have the highest TCHE bTag value.
  //These two possiblilities are considered as the two bJets and will be passed on to Kinematic Fit.
  //
  //Candidates: i&j = light jets ; k = hadronic b jet ; l = leptonic b jet
  for(int i=0;i<3;i++){
    for(int j=i+1;j<4;j++){
      for(int k=0;k<4;k++){
	if(k!=i && k!=j){
	  for(int l = 0; l < 4 ; l++){ 
	    if(l!=i && l!=j && l!=k){

	      //TCHE bTag values for the four highest pT jets
	      float btag_i = TCHEbTags_[i];
	      if(btag_i < -90) btag_i = 0;
	      float btag_j = TCHEbTags_[j];
	      if(btag_j < -90) btag_j = 0;
	      float btag_k = TCHEbTags_[k];
	      if(btag_k < -90) btag_k = 0;
	      float btag_l = TCHEbTags_[l];
	      if(btag_l < -90) btag_l = 0;		      
	      
	      //Apply bTag constraint on the TCHE bTag value of the two bJet candidates
	      bTagBoolean.push_back(applicationBTag(k,l));

	      //Calculate bTag value to discriminate between light jets and b-quark jets
	      btag[NumberCombinations] = pow(btag_k,2) + pow(btag_l,2);
	      btag[NumberCombinations] = btag[NumberCombinations] / ( pow(btag_i,2) + pow(btag_j,2) + pow(btag_k,2) + pow(btag_l,2) );

	      //Leptonic b jet defining:
	      BLeptonicIndexMVA[NumberCombinations]=l;
	      
	      //Treat events which get no bTag values specially!	      
	      //Here these events are defined to have the topology in which the two b-jets have the highest pT of the four present jets.
	      ProblemEvent = false;
	      if(btag_i ==0 && btag_j ==0 && btag_k ==0 && btag_k ==0){
		btag[NumberCombinations]=0;
		ProblemEvent=true;
	      }
	    }
	  }
	  //Number values for light jets and hadronic b:
	  QuarkOneIndex[NumberCombinations]=i;
	  QuarkTwoIndex[NumberCombinations]=j;
	  BHadronicIndex[NumberCombinations]=k;		
	  
	  NumberCombinations++;
	}
      }
    }
  }
  //From the twelve combinations the combinations with the highest BTag value are selected
  //Two combinations are selected since the hadronic and leptonic b jet can be interchanged
  for(int ii=0;ii<12;ii++){  
    if(BTagValueMVA<btag[ii] && !(fabs(BTagValueMVA-btag[ii])<0.000001)){
      BTagValueMVA=btag[ii];
      MVACombination[0]=ii;
      MVACombination[1]=ii+1;	           
    }	
    if(fabs(BTagValueMVA-btag[ii])<0.000001){
      MVACombination[1]=ii;
    }
  }	 
  //Check whether for this event the required bTag value was obtained!!
  bTagSelected=bTagBoolean[MVACombination[0]];

  //Consider events for which no bTag value was obtained
  //In this configuration the b-jets have the highest pT which is more reasonable
  if(ProblemEvent==true){
    MVACombination[0]=10;
    MVACombination[1]=11;
  }
  
  //Put the two obtained combinations from the bTag analysis into a kinematic fit 
  //Kinematic fit is applied onto hadronic top only to select the most possible hadronic b jet from the two possibilities
  for(int jj=0;jj<2; jj++){    
    MVAComb=MVACombination[jj];  

    ChiSquared[jj]=KinFitCalculation(QuarkOneIndex[MVAComb],QuarkTwoIndex[MVAComb],BHadronicIndex[MVAComb]);
  }
  //Calculate the lowest chi squared value for the two combinations
  //For this combination the masses match the best the expected ones for the hadronic top topology
  CombinationOfMVA=CalculateLowestChiSquared(2);

  //Obtain the correct numbering for the selected jet combination from the twelve possiblities
  UsedCombination=MVACombination[CombinationOfMVA];
  BLeptonicIndex=BLeptonicIndexMVA[UsedCombination];  

  //Obtain the numbering in arranged Pt value for respectively the hadronic b, leptonic b and two light jets
  //Also whether the event has passed the required bTag constraint is represented by this vector.
  CorrectJets=QuarkDistribution(UsedCombination);

  return CorrectJets;
}

float WReconstruction::KinFitCalculation(int NumberLight1, int NumberLight2, int NumberBHadr){

  // prepare everything for the Kinematic Fit
  //  TMatrixD Ml1(3,3), Ml2(3,3), Mb(3,3);
  // Ml1.Zero(); Ml2.Zero(); Mb.Zero();
  // Ml1(0,0) = pow(resFitLightJets_->EtResolution(&Light1), 2);
  // Ml1(1,1) = pow(resFitLightJets_->ThetaResolution(&Light1), 2);
  // Ml1(2,2) = pow(resFitLightJets_->PhiResolution(&Light1), 2);
  // Ml2(0,0) = pow(resFitLightJets_->EtResolution(&Light2), 2);
  // Ml2(1,1) = pow(resFitLightJets_->ThetaResolution(&Light2), 2);
  // Ml2(2,2) = pow(resFitLightJets_->PhiResolution(&Light2), 2);
  // Mb(0,0) = pow(resFitBJets_->EtResolution(&BJet), 2);
  // Mb(1,1) = pow(resFitBJets_->ThetaResolution(&BJet), 2);
  // Mb(2,2) = pow(resFitBJets_->PhiResolution(&BJet), 2);

  TLorentzVector lightJet1;
  TLorentzVector lightJet2;
  TLorentzVector bJet;
  lightJet1.SetPxPyPzE(jetsPx_[NumberLight1],jetsPy_[NumberLight1],jetsPz_[NumberLight1],jetsE_[NumberLight1]);
  lightJet2.SetPxPyPzE(jetsPx_[NumberLight2],jetsPy_[NumberLight2], jetsPz_[NumberLight2], jetsE_[NumberLight2]);
  bJet.SetPxPyPzE(jetsPx_[NumberBHadr],jetsPy_[NumberBHadr],jetsPz_[NumberBHadr],jetsE_[NumberBHadr]);
    
  TMatrixD Ml1(3,3), Ml2(3,3), Mb(3,3);
  Ml1.Zero(); Ml2.Zero(); Mb.Zero();
  Ml1(0,0) = EtResLight_[NumberLight1];
  Ml1(1,1) = ThetaResLight_[NumberLight1];
  Ml1(2,2) = PhiResLight_[NumberLight1];
  Ml2(0,0) = EtResLight_[NumberLight2];
  Ml2(1,1) = ThetaResLight_[NumberLight2];
  Ml2(2,2) = PhiResLight_[NumberLight2];
  Mb(0,0) = EtResB_[NumberBHadr];
  Mb(1,1) = ThetaResB_[NumberBHadr];
  Mb(2,2) = PhiResB_[NumberBHadr];

  TKinFitter *theFitter = new TKinFitter("hadtopFit", "hadtopFit");  
  TFitParticleEtThetaPhi *fitLight1 = new TFitParticleEtThetaPhi("lightJet1", "lightJet1", &lightJet1, &Ml1);
  TFitParticleEtThetaPhi *fitLight2 = new TFitParticleEtThetaPhi("lightJet2", "lightJet2", &lightJet2, &Ml2);
  TFitParticleEtThetaPhi *fitB = new TFitParticleEtThetaPhi("bJet", "bJet", &bJet, &Mb);
  theFitter->setVerbosity(0);
  theFitter->addMeasParticles(fitLight1,fitLight2,fitB);
  
  TFitConstraintM *consW = new TFitConstraintM("WBosonMass", "MassConstraint", 0, 0, WMassKinFit);
  TFitConstraintM *consTop = new TFitConstraintM("TopQuarkMass", "MassConstraint", 0, 0, TopMassKinFit );//Different mass for MC and Data!!
  consW->addParticles1(fitLight1,fitLight2);
  consTop->addParticles1(fitB,fitLight1,fitLight2);
  
  theFitter->addConstraint(consW);
  theFitter->addConstraint(consTop);
  theFitter->setMaxNbIter(30);
  theFitter->setMaxDeltaS(5e-5);
  theFitter->setMaxF(1e-4);
  
  //do the fit!
  theFitter->fit();
  if (theFitter->getStatus() == 0) 
    ChiSquaredFit=theFitter->getS();
  //else
  //cout << "FIT NOT CONVERGED" << endl;
  
  delete theFitter;
  delete fitLight1;
  delete fitLight2;
  delete fitB;
  delete consW;
  delete consTop;

  return ChiSquaredFit;
}

float WReconstruction::applicationBTag(int bJet1, int bJet2){

  //if(TCHEbTags_[bJet1] >1.7 && TCHEbTags_[bJet2] > 1.7) return 1;
  //else return 0;

  if(SSVHEbTags_[bJet1] > 1.74 || SSVHEbTags_[bJet2] > 1.74) return 1;
  else return 0;

  //return 1; //In case no bTag constraint should be applied
}

int WReconstruction::CalculateLowestChiSquared(int sizeLoop){

  ChiSquaredValue=ChiSquared[0];
  int Combination=0;
  for(int ii=0;ii<sizeLoop;ii++){
    if(ChiSquaredValue>ChiSquared[ii]){
      ChiSquaredValue=ChiSquared[ii];
      Combination=ii;        
    } 
  }
  return Combination;
}

std::vector<int> WReconstruction::QuarkDistribution(int QuarkDistributionValue){

  JetNumbers.push_back(BHadronicIndex[QuarkDistributionValue]);
  JetNumbers.push_back(BLeptonicIndex);
  JetNumbers.push_back(QuarkOneIndex[QuarkDistributionValue]);
  JetNumbers.push_back(QuarkTwoIndex[QuarkDistributionValue]);

  //Put in vector whether event has passed the applied bTag constraints
  JetNumbers.push_back(bTagSelected);

  return JetNumbers;
}
