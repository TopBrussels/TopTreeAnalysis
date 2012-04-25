#include "../interface/TopFCNC_GenEvt.h"

//ClassImp(TopFCNC_GenEvt);

void TopFCNC_GenEvt::ReconstructEvt(const std::vector<TRootMCParticle*> &mcParticles){

  wLeptonicChannel_ = kNone;
  zLeptonicChannel_ = kNone;

  for(unsigned int i=0; i<mcParticles.size(); i++){
	if( mcParticles[i]->status() != 3) continue;
	if( fabs(mcParticles[i]->type()) == 6 ){
		if( fabs(mcParticles[i]->dauOneId()) == 24 ) {   smDecayTop_ = *mcParticles[i]; continue; }
		else if( fabs(mcParticles[i]->dauOneId()) == 23 ) { fcncDecayTop_ = *mcParticles[i]; continue; }
		else {
			cerr << "Cannot assign top quark according to its decay..." << endl;
			exit(1);
		}
	}

	if( fabs(mcParticles[i]->motherType()) == 6 )
	{
		if( fabs(mcParticles[i]->type()) == 24 ){
			W_ = *mcParticles[i];
			continue;
		}
		else if( fabs(mcParticles[i]->type()) ==  5 ){
			B_ = *mcParticles[i];
			continue;
		}
		else if( fabs(mcParticles[i]->type()) == 23 ){
			Z_ = *mcParticles[i];
			continue;
		}
		else if( fabs(mcParticles[i]->type()) >= 1  && fabs(mcParticles[i]->type()) <= 4 ){
			Q_ = *mcParticles[i];
			continue;
		}
/*
		else {
			cerr << "Cannot assign top quark decay product. Type =" << mcParticles[i]->type() << endl;
			exit(1);
		}
*/
	}
	else if( fabs(mcParticles[i]->motherType()) == 24 && fabs(mcParticles[i]->grannyType()) == 6 ){
		if( fabs(mcParticles[i]->type()) >= 1  && fabs(mcParticles[i]->type()) <= 4 ){
			if(quarkFromW_.P() == 0){
				quarkFromW_ = *mcParticles[i];
				continue;
			}
			else{
				quarkBarFromW_ = *mcParticles[i];
				continue;
			}
		}
		else if( fabs(mcParticles[i]->type()) == 11 || fabs(mcParticles[i]->type()) == 13 || fabs(mcParticles[i]->type()) == 15 ){
			if(      fabs(mcParticles[i]->type()) == 11 ) wLeptonicChannel_ = kElec;
			else if( fabs(mcParticles[i]->type()) == 13 ) wLeptonicChannel_ = kMuon;
			else if( fabs(mcParticles[i]->type()) == 15 ) wLeptonicChannel_ = kTau;

			leptonFromW_ = *mcParticles[i];
			continue;
		}
		else if( fabs(mcParticles[i]->type()) == 12 || fabs(mcParticles[i]->type()) == 14 || fabs(mcParticles[i]->type()) == 16 ){
			neutrino_ = *mcParticles[i];
			continue;
		}
/*
		else {
			cerr << "Cannot assign W boson decay. Type =" << mcParticles[i]->type() << endl;
			exit(1);
		}
*/
	}
	else if( fabs(mcParticles[i]->motherType()) == 23 && fabs(mcParticles[i]->grannyType()) == 6 ){
		if( fabs(mcParticles[i]->type()) == 11 || fabs(mcParticles[i]->type()) == 13 || fabs(mcParticles[i]->type()) == 15 ){
			if(      fabs(mcParticles[i]->type()) == 11 ) zLeptonicChannel_ = kElec;
			else if( fabs(mcParticles[i]->type()) == 13 ) zLeptonicChannel_ = kMuon;
			else if( fabs(mcParticles[i]->type()) == 15 ) zLeptonicChannel_ = kTau;

			if(lepton1FromZ_.E() == 0){
				lepton1FromZ_ = *mcParticles[i];
				continue;
			}
			else{
				lepton2FromZ_ = *mcParticles[i];
				continue;
			}
		}
/*
		else {
			cerr << "Cannot assign Z boson decay. Type =" << mcParticles[i]->type() << endl;
			exit(1);
		}
*/
	}
  }
}

void TopFCNC_GenEvt::MatchJetsToPartons(const std::vector<TRootJet*> &jets, const int algorithm, const bool useMaxDist, const bool useDeltaR, const double maxDist){
	if(zLeptonicChannel_ == kNone){
		cerr << "Cannot find the Z leptonic decay channel. Use the method ReconstructEvt first." << endl;
		exit(1);
	}

	JetPartonMatching *myJetPartonMatcher = 0;

	std::vector<TLorentzVector> partons;
	std::vector<TLorentzVector> tljets;

	for(unsigned int i=0;i<jets.size();i++){
		tljets.push_back((TLorentzVector)*jets[i]);
	}
/*
  quarkFromW_    = 0;
  quarkBarFromW_ = 0;
  B_ = 0;
  Q_ = 0;
*/
	if(wLeptonicChannel_ == kNone){ // dileptonic event tt->Wb+Zq->qqb+llq
		partons.push_back(quarkFromW_);
		partons.push_back(quarkBarFromW_);
		partons.push_back(B_);
		partons.push_back(Q_);
	}
	else{ // trileptonic events tt->Wb+Zq->lvb+llq
		partons.push_back(B_);
		partons.push_back(Q_);
	}
	myJetPartonMatcher = new JetPartonMatching(partons,tljets,algorithm,useMaxDist,useDeltaR,maxDist);
	//myJetPartonMatcher->print();
	if(wLeptonicChannel_ == kNone){ // dileptonic event tt->Wb+Zq->qqb+llq
		matchedQuarkFromW_    = (myJetPartonMatcher->getMatchForParton(0,0) >0 ? *jets[myJetPartonMatcher->getMatchForParton(0,0)] : TRootJet());
		matchedQuarkBarFromW_ = (myJetPartonMatcher->getMatchForParton(1,0) >0 ? *jets[myJetPartonMatcher->getMatchForParton(1,0)] : TRootJet());
		matchedB_             = (myJetPartonMatcher->getMatchForParton(2,0) >0 ? *jets[myJetPartonMatcher->getMatchForParton(2,0)] : TRootJet());
		matchedQ_             = (myJetPartonMatcher->getMatchForParton(3,0) >0 ? *jets[myJetPartonMatcher->getMatchForParton(3,0)] : TRootJet());
	}
	else{
		matchedB_ = (myJetPartonMatcher->getMatchForParton(0,0) >0 ? *jets[myJetPartonMatcher->getMatchForParton(0,0)] : TRootJet());
		matchedQ_ = (myJetPartonMatcher->getMatchForParton(1,0) >0 ? *jets[myJetPartonMatcher->getMatchForParton(1,0)] : TRootJet());
	}
	if(myJetPartonMatcher) delete myJetPartonMatcher;
}

void TopFCNC_GenEvt::MatchLeptonsToZ(const std::vector<TRootMuon*> &leptons, const int algorithm, const bool useMaxDist, const bool useDeltaR, const double maxDist){
 
	JetPartonMatching *myLeptonMatcher = 0;

	std::vector<TLorentzVector> genlept;
	std::vector<TLorentzVector> tleptons;

	for(unsigned int i=0;i<leptons.size();i++){
		tleptons.push_back((TLorentzVector)*leptons[i]);
	}

  genlept.push_back(lepton1FromZ_);
  genlept.push_back(lepton2FromZ_);

	myLeptonMatcher = new JetPartonMatching(genlept,tleptons,algorithm,useMaxDist,useDeltaR,maxDist);

  if(myLeptonMatcher->getMatchForParton(0,0) >0) matchedLepton1FromZ_ = *leptons[myLeptonMatcher->getMatchForParton(0,0)];
  if(myLeptonMatcher->getMatchForParton(1,0) >0) matchedLepton2FromZ_ = *leptons[myLeptonMatcher->getMatchForParton(1,0)];

	matchedZ_ = matchedLepton1FromZ_+matchedLepton2FromZ_;

}
void TopFCNC_GenEvt::MatchLeptonsToZ(const std::vector<TRootElectron*> &leptons, const int algorithm, const bool useMaxDist, const bool useDeltaR, const double maxDist){
 
	JetPartonMatching *myLeptonMatcher = 0;

	std::vector<TLorentzVector> genlept;
	std::vector<TLorentzVector> tleptons;

	for(unsigned int i=0;i<leptons.size();i++){
		tleptons.push_back((TLorentzVector)*leptons[i]);
	}

  genlept.push_back(lepton1FromZ_);
  genlept.push_back(lepton2FromZ_);

	myLeptonMatcher = new JetPartonMatching(genlept,tleptons,algorithm,useMaxDist,useDeltaR,maxDist);

  if(myLeptonMatcher->getMatchForParton(0,0) >0) matchedLepton1FromZ_ = *leptons[myLeptonMatcher->getMatchForParton(0,0)];
  if(myLeptonMatcher->getMatchForParton(1,0) >0) matchedLepton2FromZ_ = *leptons[myLeptonMatcher->getMatchForParton(1,0)];

	matchedZ_ = matchedLepton1FromZ_+matchedLepton2FromZ_;

}

void TopFCNC_GenEvt::FillResolutions(ResolutionFit* resFitLeptons, ResolutionFit* resFitLightJets, ResolutionFit* resFitBJets){

  if(foundZ_){
    resFitLeptons->Fill(&matchedLepton1FromZ_,&lepton1FromZ_);
    resFitLeptons->Fill(&matchedLepton2FromZ_,&lepton2FromZ_);
  }
  resFitBJets->Fill(&matchedB_,&B_);
  resFitLightJets->Fill(&matchedQ_,&Q_);
	if(wLeptonicChannel_ == kNone){ // dileptonic event tt->Wb+Zq->qqb+llq
		resFitLightJets->Fill(&matchedQuarkFromW_,&quarkFromW_);
		resFitLightJets->Fill(&matchedQuarkBarFromW_,&quarkBarFromW_);
  }
}

