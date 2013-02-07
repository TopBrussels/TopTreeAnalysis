#include "../interface/BTagCosThetaCalculation.h"

BTagCosThetaCalculation::BTagCosThetaCalculation(){

}

BTagCosThetaCalculation::~BTagCosThetaCalculation(){

}

float BTagCosThetaCalculation::Calculation(TLorentzVector lepton, TLorentzVector Neutrino, TLorentzVector leptonicBJet){

  float CosTheta = 999;	  
  
  //----------------------------------------------
  //  Calculating cos theta value
  //----------------------------------------------
  
  TRootMCParticle WLeptonic = (Neutrino+lepton);  
  TRootMCParticle TopLeptonic = (Neutrino+lepton+leptonicBJet);
  TLorentzVector TopWRF = (Neutrino+lepton+leptonicBJet);
  TLorentzVector leptWRF = lepton;

  //Angle between Top in WRF and lepton in WRF
  TopWRF.Boost(-WLeptonic.BoostVector());
  leptWRF.Boost(-WLeptonic.BoostVector());
 
  //Calculating cos:	      
  float ThetaTevatron = ROOT::Math::VectorUtil::Angle( TopWRF, leptWRF );
  CosTheta = -(TMath::Cos(ThetaTevatron));
  //Cos theta is defined as the angle between the lepton and the reversed direction of the top quark, both boosted to the W-boson rest frame.
  //Still reversed direction doesn't need to be calculated since angle between lepton and top and between lepton and reversed top is proportional to theta and Pi-theta.
  //For these the angles the following relation holds: cos(theta) = - cos(Pi-theta) 
  // --> Hence the need of the minus sign in the CosTheta definition!!

  if(WLeptonic.E() < 0.){
    cout << " Event with negative WLept energy!!! (BTagCosThetaCalculation class) Cos theta = " << CosTheta << endl;
  }
  
  return CosTheta;	
}        

float BTagCosThetaCalculation::CalcOrigKins(int BLeptonicIndex, int BHadronicIndex,  TLorentzVector muon, vector<TLorentzVector> selectedJets, float MassW, float MassTop){
  
  float CosThetaOrigKins = 999; 
  //---------------------------------------------------------------
  //     Calculating MET_Pz() (equation ax² + bx + c = 0 ):     
  //---------------------------------------------------------------
  //
  // MET_Px() and MET_Py() is known since there is no Px() and Py() component before the collision.
  // The Pz() component before the collision is not known for proton-proton collisions.
  //
  float NeutrinoPx = -(muon +selectedJets[0]+ selectedJets[1] + selectedJets[2] + selectedJets[3] ).Px();
  float NeutrinoPy = -(muon +selectedJets[0]+ selectedJets[1] + selectedJets[2] + selectedJets[3] ).Py();	
  float NeutrinoPz=999;
  float NeutrinoE=999;
  Neutrino.SetPxPyPzE(NeutrinoPx,NeutrinoPy,NeutrinoPz,NeutrinoE);  //Initialize in case no neutrino can be reconstructed!!	  
      
  //Calculating solutions for quadratic equation  
  float aCoefficient = 4*pow(muon.E(),2)-4*pow(muon.Pz(),2);
  float bCoefficient = 4*(muon.Pz())*(pow(muon.M(),2)-pow(MassW,2)-2*(muon.Px())*NeutrinoPx-2*(muon.Py())*NeutrinoPy);
  float cCoefficient = -pow(muon.M(),4)-pow(MassW,4)-4*pow(muon.Px(),2)*pow(NeutrinoPx,2)-4*pow(muon.Py(),2)*pow(NeutrinoPy,2)+4*(pow(muon.M(),2)-pow(MassW,2))*((muon.Px())*NeutrinoPx+(muon.Py())*NeutrinoPy)-8*(muon.Px())*NeutrinoPx*(muon.Py())*NeutrinoPy+4*pow(muon.E(),2)*pow(NeutrinoPx,2)+4*pow(muon.E(),2)*pow(NeutrinoPy,2)+2*(pow(MassW,2))*(pow(muon.M(),2));
      
  float DCoefficient = pow(bCoefficient,2)-4*aCoefficient*cCoefficient;           
  if(DCoefficient>=0){	  //Still need to find solution for D<0 !!!	      
    float NeutrinoPzOne = ((-bCoefficient + sqrt(DCoefficient))/(aCoefficient*2));
    float NeutrinoPzTwo = ((-bCoefficient - sqrt(DCoefficient))/(aCoefficient*2));	  
	
    float NeutrinoEOne = sqrt(NeutrinoPx*NeutrinoPx+NeutrinoPy*NeutrinoPy+NeutrinoPzOne*NeutrinoPzOne);
    float NeutrinoETwo = sqrt(NeutrinoPx*NeutrinoPx+NeutrinoPy*NeutrinoPy+NeutrinoPzTwo*NeutrinoPzTwo);
    TLorentzVector NeutrinoOne;
    NeutrinoOne.SetPxPyPzE(NeutrinoPx,NeutrinoPy,NeutrinoPzOne,NeutrinoEOne);  //Has to be initialized like this !
    TLorentzVector NeutrinoTwo;
    NeutrinoTwo.SetPxPyPzE(NeutrinoPx,NeutrinoPy,NeutrinoPzTwo,NeutrinoETwo);
	
    float TopMassDiffOne = (MassTop -((NeutrinoOne+muon+selectedJets[BLeptonicIndex]).M()) );	  
    float TopMassDiffTwo = (MassTop -((NeutrinoTwo+muon+selectedJets[BLeptonicIndex]).M()) );	   
	
    //-----   Selected neutrino has the smallest deviation from the top mass   -----
    //
    if(fabs(TopMassDiffOne)<fabs(TopMassDiffTwo)){
      NeutrinoPz=NeutrinoPzOne;
      NeutrinoE = NeutrinoEOne;
    }
    else{
      NeutrinoPz=NeutrinoPzTwo;
      NeutrinoE = NeutrinoETwo;
    }	  	  
    Neutrino.SetPxPyPzE(NeutrinoPx,NeutrinoPy,NeutrinoPz,NeutrinoE);

    //----------------------------------------------
    //  Calculating cos theta value
    //----------------------------------------------

    TLorentzVector WLeptonic = (Neutrino+muon);
  
    TLorentzVector TopWRF = (Neutrino+muon+selectedJets[BLeptonicIndex]);
    TLorentzVector leptWRF = muon;

    TopWRF.Boost(-WLeptonic.BoostVector());
    leptWRF.Boost (-WLeptonic.BoostVector());
  
    //Calculating cos:	      
    float ThetaTevatron = ROOT::Math::VectorUtil::Angle( TopWRF, leptWRF );
    CosThetaOrigKins = -(TMath::Cos(ThetaTevatron));
		
		if(WLeptonic.E() < 0.){
    cout << " Event with negative WLept energy!!! (BTagCosThetaCalculation class) Cos theta = " << CosThetaOrigKins << endl;
  	}       	
 
  }  //end of D>=0 loop   
	//std::cout << "Cos theta Check!" << std::endl;
  return CosThetaOrigKins;
	
}        
