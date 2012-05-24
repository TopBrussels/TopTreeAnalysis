#ifndef TopFCNC_Evt_h
#define TopFCNC_Evt_h

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include "TLorentzVector.h"
#include "TopTreeProducer/interface/TRootJet.h"
#include "TopTreeProducer/interface/TRootMET.h"
#include "TopTreeAnalysis/Reconstruction/interface/Combination.h"
#include "TopTreeAnalysis/Reconstruction/interface/MEzCalculator.h"
#include "TopTreeAnalysis/Selection/interface/Selection.h"

using namespace std;
using namespace stdcomb;

class TopFCNC_Evt : public TObject
{
	public: 
		// W/Z boson leptonic decay channel
		enum LeptonType {kNone, kMuon, kElec, kTau};

		TopFCNC_Evt() :
			TObject(),
      eventID_(0),
      runID_(0),
      lumiBlockID_(0),
      idParton1_(-1),
      idParton2_(-1),
      xParton1_(-1.),
      xParton2_(-1.),
      factorizationScale_(-1.),
      nPV_(0),
      nPUBXm1_(0),
      nPU_(0),
      nPUBXp1_(0),
      nTruePU_(-1.),
      eventWeight_(1.),
			zLeptonicChannel_(kNone),
			wLeptonicChannel_(kNone),
			isDiLeptonic_(false),
			isTriLeptonic_(false),
			wMass_(80.4),
			zMass_(91.2),
			topMass_(170.),
			hadWMassResol_(10.),
			hadTopMassResol_(17.4),
			lepWMassResol_(20.),
			lepZMassResol_(4.),
			lepTopMassResol_(17.2),
			fcncTopMassResol_(10.),
			massChi2_(-999.),
			massChi2_SM_(-999.),
			massChi2_FCNC_(-999.),
			verbose_(false)
			{;}

		TopFCNC_Evt(LeptonType type) :
			TObject(),
      eventID_(0),
      runID_(0),
      lumiBlockID_(0),
      idParton1_(-1),
      idParton2_(-1),
      xParton1_(-1.),
      xParton2_(-1.),
      factorizationScale_(-1.),
      nPV_(0),
      nPUBXm1_(0),
      nPU_(0),
      nPUBXp1_(0),
      nTruePU_(-1.),
      eventWeight_(1.),
			zLeptonicChannel_(type),
			wLeptonicChannel_(kNone),
			isDiLeptonic_(true),
			isTriLeptonic_(false),
			wMass_(80.4),
			zMass_(91.2),
			topMass_(170.),
			hadWMassResol_(10.),
			hadTopMassResol_(17.4),
			lepWMassResol_(20.),
			lepZMassResol_(4.),
			lepTopMassResol_(17.2),
			fcncTopMassResol_(10.),
			massChi2_(-999.),
			massChi2_SM_(-999.),
			massChi2_FCNC_(-999.),
			verbose_(false)
			{;}

		TopFCNC_Evt(LeptonType type1, LeptonType type2) :
			TObject(),
      eventID_(0),
      runID_(0),
      lumiBlockID_(0),
      idParton1_(-1),
      idParton2_(-1),
      xParton1_(-1.),
      xParton2_(-1.),
      factorizationScale_(-1.),
      nPV_(0),
      nPUBXm1_(0),
      nPU_(0),
      nPUBXp1_(0),
      nTruePU_(-1.),
      eventWeight_(1.),
			zLeptonicChannel_(type1),
			wLeptonicChannel_(type2),
			isDiLeptonic_(false),
			isTriLeptonic_(true),
			wMass_(80.4),
			zMass_(91.2),
			topMass_(170.),
			hadWMassResol_(10.),
			hadTopMassResol_(17.4),
			lepWMassResol_(20.),
			lepZMassResol_(4.),
			lepTopMassResol_(17.2),
			fcncTopMassResol_(10.),
			massChi2_(-999.),
			massChi2_SM_(-999.),
			massChi2_FCNC_(-999.),
			verbose_(false)
			{;}

		TopFCNC_Evt(const TopFCNC_Evt& evt) :
			TObject(evt),
      eventID_(evt.eventID_),
      runID_(evt.runID_),
      lumiBlockID_(evt.lumiBlockID_),
      idParton1_(evt.idParton1_),
      idParton2_(evt.idParton2_),
      xParton1_(evt.xParton1_),
      xParton2_(evt.xParton2_),
      factorizationScale_(evt.factorizationScale_),
      nPV_(evt.nPV_),
      nPUBXm1_(evt.nPUBXm1_),
      nPU_(evt.nPU_),
      nPUBXp1_(evt.nPUBXp1_),
      nTruePU_(evt.nTruePU_),
      eventWeight_(evt.eventWeight_),
			zLeptonicChannel_(evt.zLeptonicChannel_),
			wLeptonicChannel_(evt.wLeptonicChannel_),
			isDiLeptonic_(evt.isDiLeptonic_),
			isTriLeptonic_(evt.isTriLeptonic_),	
			smDecayTop_(evt.smDecayTop_),
			W_(evt.W_),
			B_(evt.Z_),
			fcncDecayTop_(evt.fcncDecayTop_),
			Z_(evt.Z_),
			Q_(evt.Q_),
			leptonFromW_(evt.leptonFromW_),
			neutrino_(evt.neutrino_),
			quarkFromW_(evt.quarkFromW_),
			quarkBarFromW_(evt.quarkBarFromW_),
			lepton1FromZ_(evt.lepton1FromZ_),
			lepton2FromZ_(evt.lepton2FromZ_),
			selectedJets_(evt.selectedJets_),
			ISR_(evt.ISR_),
			smDecayTopRadiation_(evt.smDecayTopRadiation_),
			fcncDecayTopRadiation_(evt.fcncDecayTopRadiation_),
			wMass_(evt.wMass_),
			zMass_(evt.zMass_),
			topMass_(evt.topMass_),
			hadWMassResol_(evt.hadWMassResol_),
			hadTopMassResol_(evt.hadTopMassResol_),
			lepWMassResol_(evt.lepWMassResol_),
			lepZMassResol_(evt.lepZMassResol_),
			lepTopMassResol_(evt.lepTopMassResol_),
			fcncTopMassResol_(evt.fcncTopMassResol_),
			massChi2_(evt.massChi2_),
			massChi2_SM_(evt.massChi2_SM_),
			massChi2_FCNC_(evt.massChi2_FCNC_),
			verbose_(evt.verbose_)
			{;}
	
		~TopFCNC_Evt(){;}

    const UInt_t eventID()     const { return eventID_;}
    const UInt_t runID()       const { return runID_;}
    const UInt_t lumiBlockID() const { return lumiBlockID_;}

    const Int_t idParton1() const { return idParton1_;}
    const Int_t idParton2() const { return idParton2_;}

    const Float_t xParton1()           const { return xParton1_;}
    const Float_t xParton2()           const { return xParton2_;}
    
    const Float_t factorizationScale() const { return factorizationScale_;}

    const UInt_t  nPV()     const { return nPV_;}
    const UInt_t  nPUBXm1() const { return nPUBXm1_;}
    const UInt_t  nPU()     const { return nPU_;}
    const UInt_t  nPUBXp1() const { return nPUBXp1_;}
    const Float_t nTruePU() const { return nTruePU_;}

    const Float_t eventWeight() const { return eventWeight_;}

		LeptonType zLeptonicChannel() const { return zLeptonicChannel_;}
		LeptonType wLeptonicChannel() const { return wLeptonicChannel_;}
	  
		Bool_t isDiLeptonic()  const {return isDiLeptonic_;}
		Bool_t isTriLeptonic() const {return isTriLeptonic_;}
		Bool_t isDiLeptonic(LeptonType t)                 const { return (zLeptonicChannel_==t ? true : false); }
		Bool_t isTriLeptonic(LeptonType t1,LeptonType t2) const { return((wLeptonicChannel_==t1 && zLeptonicChannel_==t2) ? true : false);}

		const TRootParticle leptonFromW()   const { return leptonFromW_;}
		const TRootParticle neutrino()      const { return neutrino_;}
		const TRootJet quarkFromW()         const { return quarkFromW_;}
		const TRootJet quarkBarFromW()      const { return quarkBarFromW_;}

		const TRootParticle lepton1FromZ()  const { return lepton1FromZ_;}
		const TRootParticle lepton2FromZ()  const { return lepton2FromZ_;}

		const TRootParticle smDecayTop()    const { return smDecayTop_;}
		const TRootParticle W()             const { return W_;}
		const TRootJet B()                  const { return B_;}
		const TRootParticle fcncDecayTop()  const { return fcncDecayTop_;}
		const TRootParticle Z()             const { return Z_;}
		const TRootJet Q()                  const { return Q_;}

    const std::vector<TRootJet> selectedJets()             const { return selectedJets_;}
		const std::vector<TRootJet> ISR()                      const { return ISR_;}
		const std::vector<TRootJet> smDecayTopRadiation()      const { return smDecayTopRadiation_;}
		const std::vector<TRootJet> fcncDecayTopRadiation()    const { return fcncDecayTopRadiation_;}

		const Float_t WMass()   const { return wMass_; }
		const Float_t ZMass()   const { return zMass_; }
		const Float_t TopMass() const { return topMass_; }

		const Float_t HadWMassResol()    const { return hadWMassResol_   ; }
		const Float_t HadTopMassResol()  const { return hadTopMassResol_ ; }
		const Float_t LepWMassResol()    const { return lepWMassResol_   ; }
		const Float_t LepZMassResol()    const { return lepZMassResol_   ; }
		const Float_t LepTopMassResol()  const { return lepTopMassResol_ ; }
		const Float_t FcncTopMassResol() const { return fcncTopMassResol_; }

		const Float_t MassChi2()      const { return massChi2_; }
		const Float_t MassChi2_SM()   const { return massChi2_SM_; }
		const Float_t MassChi2_FCNC() const { return massChi2_FCNC_; }

    void SetEventID(UInt_t eventID)         { eventID_ = eventID; }
    void SetRunID(UInt_t runID)             { runID_ = runID; }
    void SetLumiBlockID(UInt_t lumiBlockID) { lumiBlockID_ = lumiBlockID; }

    void SetIdParton1(Int_t idParton1) { idParton1_ = idParton1; }
    void SetIdParton2(Int_t idParton2) { idParton2_ = idParton2; }

    void SetxParton1(Float_t xParton1)           { xParton1_ = xParton1; }
    void SetxParton2(Float_t xParton2)           { xParton2_ = xParton2; }
    void SetFactorizationScale(Float_t factorizationScale) { factorizationScale_ = factorizationScale; }

    void SetnPV(UInt_t nPV)         { nPV_ = nPV; }
    void SetnPUBXm1(UInt_t nPUBXm1) { nPUBXm1_ = nPUBXm1; }
    void SetnPU(UInt_t nPU)         { nPU_ = nPU; }
    void SetnPUBXp1(UInt_t nPUBXp1) { nPUBXp1_ = nPUBXp1; }
    void SetnTruePU(Float_t nTruePU){ nTruePU_ = nTruePU; }

    void SetEventWeight(Float_t eventWeight) { eventWeight_ = eventWeight;}

		void SetDiLeptonicChannel(LeptonType type)
		{
			isDiLeptonic_ = true;
			isTriLeptonic_= false;
			zLeptonicChannel_ = type;
		}
		void SetTriLeptonicChannel(LeptonType type1, LeptonType type2)
		{
			isDiLeptonic_  = false;
			isTriLeptonic_ = true;
			zLeptonicChannel_ = type1;
			wLeptonicChannel_ = type2;
		}

		void SetTLorentzVector(TRootParticle &smDecayTop, TRootParticle &W, TRootJet &B, TRootParticle &fcncDecayTop, TRootParticle &Z, TRootJet &Q, TRootParticle &leptonFromW, TRootParticle &neutrino, TRootJet &quarkFromW, TRootJet &quarkBarFromW, TRootParticle &lepton1FromZ, TRootParticle &lepton2FromZ)
		{
			smDecayTop_    = smDecayTop;
			W_             = W;
			B_             = B;
			fcncDecayTop_  = fcncDecayTop;
			Z_             = Z;
			Q_             = Q;

			leptonFromW_   = leptonFromW;
			neutrino_      = neutrino;
			quarkFromW_    = quarkFromW;
			quarkBarFromW_ = quarkBarFromW;

			lepton1FromZ_  = lepton1FromZ;
			lepton2FromZ_  = lepton2FromZ;

		}

		void SetSmDecayTop(TRootParticle smDecayTop){ smDecayTop_ = smDecayTop; }
		void SetW(TRootParticle W){ W_ = W; }
		void SetB(TRootJet B){ B_ = B; }
		void SetFcncDecayTop(TRootParticle fcncDecayTop){ fcncDecayTop_ = fcncDecayTop; }
		void SetZ(TRootParticle Z){ Z_ = Z; }
		void SetQ(TRootJet Q){ Q_ = Q; }

		void SetLeptonFromW(TRootParticle leptonFromW){ leptonFromW_ = leptonFromW; }
		void SetNeutrino(TRootParticle neutrino){ neutrino_ = neutrino; }
		void SetQuarkFromW(TRootJet quarkFromW){ quarkFromW_ = quarkFromW; }
		void SetQuarkBarFromW(TRootJet quarkBarFromW){ quarkBarFromW_ = quarkBarFromW; }

		void SetLepton1FromZ(TRootParticle lepton1FromZ){ lepton1FromZ_ = lepton1FromZ; }
		void SetLepton2FromZ(TRootParticle lepton2FromZ){ lepton2FromZ_ = lepton2FromZ; }

		void SetRadiation(std::vector<TRootJet> smDecayTopRadiation, std::vector<TRootJet> fcncDecayTopRadiation, std::vector<TRootJet> ISR)
		{
			smDecayTopRadiation_   = smDecayTopRadiation; 
			fcncDecayTopRadiation_ = fcncDecayTopRadiation;
			ISR_ = ISR; 
		}
    void SetSelectedJets(std::vector<TRootJet> selectedJets){ selectedJets_ = selectedJets; }
		void SetISR(std::vector<TRootJet> ISR){ISR_ = ISR;}
		void SetSmDecayTopRadiation(std::vector<TRootJet> smDecayTopRadiation){ smDecayTopRadiation_ = smDecayTopRadiation; }
		void SetFcncDecayTopRadiation(std::vector<TRootJet> fcncDecayTopRadiation){ fcncDecayTopRadiation_ = fcncDecayTopRadiation; }

		void SetWMass(Float_t wMass){ wMass_ = wMass; }
		void SetZMass(Float_t zMass){ zMass_ = zMass; }
		void SetTopMass(Float_t topMass){ topMass_ = topMass; }

		void SetHadWMassResol(Float_t hadWMassResol)       { hadWMassResol_    = hadWMassResol; }
		void SetHadTopMassResol(Float_t hadTopMassResol)   { hadTopMassResol_  = hadTopMassResol; }
		void SetLepWMassResol(Float_t lepWMassResol)       { lepWMassResol_    = lepWMassResol; }
		void SetLepZMassResol(Float_t lepZMassResol)       { lepZMassResol_    = lepZMassResol; }
		void SetLepTopMassResol(Float_t lepTopMassResol)   { lepTopMassResol_  = lepTopMassResol; }
		void SetFcncTopMassResol(Float_t fcncTopMassResol) { fcncTopMassResol_ = fcncTopMassResol; }

		void SetVerbosity(Bool_t verbose)                  { verbose_ = verbose; }

		void ReconstructEvt()
		{
			if(isDiLeptonic_ == isTriLeptonic_)
			{
				cerr << "Specify if top FCNC event candidate is di- or tri-leptonic before event reconstruction" << endl;
				exit(1);
			}
			else if(isDiLeptonic_)  ReconstructDiLeptEvt();
			else if(isTriLeptonic_) ReconstructTriLeptEvt();
		}

		void ReconstructDiLeptEvt()
		{
			if(!isDiLeptonic_)
			{
				cerr << "Top FCNC event candidate is not di-leptonic. Cannot be reconstructed as such" << endl;
				exit(1);
			}
			W_ = quarkFromW_ + quarkBarFromW_ ;
			smDecayTop_ = B_ + W_ ;
			Z_ = lepton1FromZ_ + lepton2FromZ_ ;
			fcncDecayTop_ = Q_ + Z_ ;
			massChi2_SM_   = pow((W_.M()-wMass_)/hadWMassResol_,2)+pow((smDecayTop_.M()-topMass_)/hadTopMassResol_,2);
			massChi2_FCNC_ = pow((Z_.M()-zMass_)/lepZMassResol_,2)+pow((fcncDecayTop_.M()-topMass_)/fcncTopMassResol_,2);
			massChi2_      = pow((W_.M()-wMass_)/hadWMassResol_,2)+pow((smDecayTop_.M()-topMass_)/hadTopMassResol_,2)+pow((fcncDecayTop_.M()-topMass_)/fcncTopMassResol_,2);
		}

		void ReconstructTriLeptEvt()
		{
			if(!isTriLeptonic_)
			{
				cerr << "Top FCNC event candidate is not tri-leptonic. Cannot be reconstructed as such" << endl;
				exit(1);
			}
			W_ = leptonFromW_+ neutrino_ ;
			smDecayTop_ = B_ + W_ ;
			Z_ = lepton1FromZ_ + lepton2FromZ_ ;
			fcncDecayTop_ = Q_ + Z_ ;
			massChi2_SM_   = pow((W_.M()-wMass_)/lepWMassResol_,2)+pow((smDecayTop_.M()-topMass_)/lepTopMassResol_,2);
			massChi2_FCNC_ = pow((Z_.M()-zMass_)/lepZMassResol_,2)+pow((fcncDecayTop_.M()-topMass_)/fcncTopMassResol_,2);
			massChi2_      = pow((W_.M()-wMass_)/lepWMassResol_,2)+pow((smDecayTop_.M()-topMass_)/lepTopMassResol_,2)+pow((fcncDecayTop_.M()-topMass_)/fcncTopMassResol_,2);
		}

		void ReconstructDiLeptEvt(const TRootParticle* lept1, const TRootParticle* lept2, const std::vector<TopTree::TRootJet*> &init_jets, bool UseMinChi2 = true)
		{
			// Topology to reconstruct : tt-> bW + qZ -> bqq + qll
			lepton1FromZ_ = *lept1;
			lepton2FromZ_ = *lept2;
			selectedJets_.clear();
			for(UInt_t i=0;i<init_jets.size();++i) selectedJets_.push_back(*init_jets[i]);
			//std::vector<TRootJet*> jets = init_jets;

			if(UseMinChi2)
			{
				UInt_t NofJets = 4 ;
				if(selectedJets_.size()<NofJets){
          cout<<"*** Number of selected jets less than four ***"<<endl;
				  return;
				}
				UInt_t *numbers = new UInt_t[selectedJets_.size()];
				for(UInt_t i=0;i<selectedJets_.size();i++) numbers[i]=i;
				UInt_t *comb = new UInt_t[NofJets];
				for(UInt_t i=0;i <NofJets;i++) comb[i]=i;
				UInt_t *Permutation = new UInt_t[NofJets];
				for(UInt_t i=0;i<NofJets;i++) Permutation[i]=i;

				double massChi2_tmp_ = 99999999.;

				do
				{
					if(verbose_)
						cout<<"-- Jet combination considered : "<<comb[0]<<"/"<<comb[1]<<"/"<<comb[2]<<"/"<<comb[3]<<endl;

					massChi2_tmp_ = 99999999.;

					do
					{
						if(verbose_)
							cout<<"--- Permutations : "<<comb[0]<<"/"<<comb[1]<<"/"<<comb[2]<<"/"<<comb[3]<<endl;

						quarkFromW_    = selectedJets_[comb[0]];
						quarkBarFromW_ = selectedJets_[comb[1]];
						B_             = selectedJets_[comb[2]];
						Q_             = selectedJets_[comb[3]];
						ReconstructDiLeptEvt();

						if(massChi2_ < massChi2_tmp_)
						{
							for(UInt_t i=0;i<NofJets;i++) Permutation[i] = comb[i];
							massChi2_tmp_ = massChi2_ ;
						}
					}
					while(next_permutation(comb,comb+NofJets));
   			}
   			while(next_combination(numbers,numbers+selectedJets_.size(),comb,comb+NofJets));

				quarkFromW_    = selectedJets_[Permutation[0]];
				quarkBarFromW_ = selectedJets_[Permutation[1]];
				B_             = selectedJets_[Permutation[2]];
				Q_             = selectedJets_[Permutation[3]];
				ReconstructDiLeptEvt();

        delete Permutation;
        delete comb;
        delete numbers;
			}
			else
			{
				// Search for jet pairs with an invariant mass matching the W boson mass.
				UInt_t W_idx_1 = 0;
				UInt_t W_idx_2 = 1;
				Float_t W_massDiff = 99999.;
				for(UInt_t i=0;i<selectedJets_.size()-1;i++){
					for(UInt_t j=i+1;j<selectedJets_.size();j++){
						if(W_massDiff>fabs((selectedJets_[i]+selectedJets_[j]).M()-wMass_)){
							W_idx_1 = i;
							W_idx_2 = j;
							W_massDiff = fabs((selectedJets_[i]+selectedJets_[j]).M()-wMass_);
						}
					}
				}
				quarkFromW_    = selectedJets_[W_idx_1];
				quarkBarFromW_ = selectedJets_[W_idx_2];
				// Erase W boson jets candidates from the jet list
				selectedJets_.erase(selectedJets_.begin()+W_idx_2);
				selectedJets_.erase(selectedJets_.begin()+W_idx_1);
	
				Float_t Top_massDiff = 99999.;
				UInt_t B_idx = 0;
				UInt_t Q_idx = 0;

				for(UInt_t i=0;i<selectedJets_.size();i++){
					if(Top_massDiff>fabs((selectedJets_[i]+quarkFromW_+quarkBarFromW_).M()-topMass_)){
						B_idx = i;
						Top_massDiff = fabs((selectedJets_[i]+quarkFromW_+quarkBarFromW_).M()-topMass_);
					}
				}
				B_ = selectedJets_[B_idx];
				// Erase b-jet candidates from the jet list
				selectedJets_.erase(selectedJets_.begin()+B_idx);

				Top_massDiff = 99999.;
				// Search for the jet giving with the Z candidate an invariant mass close to the top mass
				for(UInt_t i=0;i<selectedJets_.size();i++){
					if(Top_massDiff>fabs((selectedJets_[i]+lepton1FromZ_+lepton2FromZ_).M()-topMass_)){
						Q_idx = i;
						Top_massDiff = fabs((selectedJets_[i]+lepton1FromZ_+lepton2FromZ_).M()-topMass_);
					}
				}
				Q_ = selectedJets_[Q_idx];
				ReconstructDiLeptEvt();
			}
		}

		void ReconstructTriLeptEvt(const TRootParticle* leptZ1, const TRootParticle* leptZ2, const TRootParticle* leptW, const std::vector<TopTree::TRootJet*> &init_jets, const TRootParticle *MET, bool UseMinChi2 = true)
		{
			// Topology to reconstruct : tt-> bW + qZ -> blv + qll
			lepton1FromZ_ = *leptZ1;
			lepton2FromZ_ = *leptZ2;
			leptonFromW_  = *leptW;
			neutrino_     = *MET;
			// Recover the kinematics of the neutrino from the MET
			MEzCalculator MyMEzCalc;
			MyMEzCalc.SetMET(*MET);
			MyMEzCalc.SetMuon(*leptW);
			Float_t MEz = MyMEzCalc.Calculate();
			neutrino_.SetPz(MEz);

			selectedJets_.clear();
			for(UInt_t i=0;i<init_jets.size();++i) selectedJets_.push_back(*init_jets[i]);
			//std::vector<TRootJet*> jets = init_jets;

			if(UseMinChi2)
			{
				UInt_t NofJets = 2 ;
				UInt_t *numbers = new UInt_t[selectedJets_.size()];
				for(UInt_t i=0;i<selectedJets_.size();i++) numbers[i]=i;
				UInt_t *comb = new UInt_t[NofJets];
				for(UInt_t i=0;i <NofJets;i++) comb[i]=i;
				UInt_t *Permutation = new UInt_t[NofJets];
				for(UInt_t i=0;i<NofJets;i++) Permutation[i]=i;

				double massChi2_tmp_ = 99999999.;

				do
				{
					if(verbose_)
						cout<<"-- Jet combination considered : "<<comb[0]<<"/"<<comb[1]<<endl;

					massChi2_tmp_ = 99999999.;

					do
					{
						if(verbose_)
							cout<<"--- Permutations : "<<comb[0]<<"/"<<comb[1]<<endl;

						B_             = selectedJets_[comb[0]];
						Q_             = selectedJets_[comb[1]];
						ReconstructTriLeptEvt();

						if(massChi2_ < massChi2_tmp_)
						{
							for(UInt_t i=0;i<NofJets;i++) Permutation[i] = comb[i];
							massChi2_tmp_ = massChi2_ ;
						}
					}
					while(next_permutation(comb,comb+NofJets));
   				}
   				while(next_combination(numbers,numbers+selectedJets_.size(),comb,comb+NofJets));

				B_             = selectedJets_[Permutation[0]];
				Q_             = selectedJets_[Permutation[1]];
				ReconstructTriLeptEvt();

			}
			else
			{
				// Search for the jet giving with the Z candidate an invariant mass close to the top mass
				Float_t Top_massDiff = 99999.;
				UInt_t B_idx = 0;
				UInt_t Q_idx = 0;

				for(UInt_t i=0;i<selectedJets_.size();i++){
					if(Top_massDiff>fabs((selectedJets_[i]+leptonFromW_+neutrino_).M()-topMass_)){
						B_idx = i;
						Top_massDiff = fabs((selectedJets_[i]+leptonFromW_+neutrino_).M()-topMass_);
					}
				}
				B_ = selectedJets_[B_idx];
				// Erase b-jet candidates from the jet list
				selectedJets_.erase(selectedJets_.begin()+B_idx);

				Top_massDiff = 99999.;
				for(UInt_t i=0;i<selectedJets_.size()-1;i++){
					if(Top_massDiff>fabs((selectedJets_[i]+lepton1FromZ_+lepton2FromZ_).M()-topMass_)){
						Q_idx = i;
						Top_massDiff = fabs((selectedJets_[i]+lepton1FromZ_+lepton2FromZ_).M()-topMass_);
					}
				}
				Q_ = selectedJets_[Q_idx];
				ReconstructTriLeptEvt();
			}
	}
	
		virtual TString typeName() const { return "TopFCNC_Evt"; }

		friend std::ostream& operator<< (std::ostream& stream, const TopFCNC_Evt& fcncEvt)
		{
			stream << "TopFCNC_Evt - is ";
			if(fcncEvt.isDiLeptonic()) {
				stream <<" DiLeptonic (";
				if(fcncEvt.zLeptonicChannel()==0) stream <<" no channel )";
				else if(fcncEvt.zLeptonicChannel()==1) stream <<" muonic )";
				else if(fcncEvt.zLeptonicChannel()==2) stream <<" electronic )";
				else if(fcncEvt.zLeptonicChannel()==3) stream <<" tauonic )";
			}
			else if(fcncEvt.isTriLeptonic()){
				stream <<" TriLeptonic (";
				if(fcncEvt.wLeptonicChannel()==0) stream << "W : no channel and ";
				else if(fcncEvt.wLeptonicChannel()==1) stream << "W : muonic and ";
				else if(fcncEvt.wLeptonicChannel()==2) stream << "W : electronic and ";
				else if(fcncEvt.wLeptonicChannel()==3) stream << "W : tauonic and ";

				if(fcncEvt.zLeptonicChannel()==0) stream << "Z : no channel )";
				else if(fcncEvt.zLeptonicChannel()==1) stream << "Z : muonic )";
				else if(fcncEvt.zLeptonicChannel()==2) stream << "Z : electronic )";
				else if(fcncEvt.zLeptonicChannel()==3) stream << "Z : tauonic )";
			}
			stream << std::endl;
			stream << "Nof ISR: "<< fcncEvt.ISR().size() << std::endl;
			stream << "Nof Top radiations: "<< fcncEvt.smDecayTopRadiation().size() + fcncEvt.fcncDecayTopRadiation().size() << std::endl;
			//stream << "TRootGenEvent - Charge=" << setw(2) << jet.charge() << " (Et,eta,phi)=("<< setw(10) << jet.Et() <<","<< setw(10) << jet.Eta() <<","<< setw(10) << jet.Phi() << ")"
			return stream;
		};


//	private:
	protected:

    UInt_t eventID_;
    UInt_t runID_;
    UInt_t lumiBlockID_;

    Int_t idParton1_;
    Int_t idParton2_;

    Float_t xParton1_;
    Float_t xParton2_;
    Float_t factorizationScale_;

    UInt_t  nPV_;
    UInt_t  nPUBXm1_;
    UInt_t  nPU_;
    UInt_t  nPUBXp1_;
    Float_t nTruePU_;

    Float_t eventWeight_;

		LeptonType zLeptonicChannel_;
		LeptonType wLeptonicChannel_;

		Bool_t isDiLeptonic_;
		Bool_t isTriLeptonic_;

		TRootParticle  smDecayTop_;
		TRootParticle  W_;
		TRootJet       B_;
		TRootParticle  fcncDecayTop_;
		TRootParticle  Z_;
		TRootJet       Q_;

		TRootParticle leptonFromW_;
		TRootParticle neutrino_;
		TRootJet      quarkFromW_;
		TRootJet      quarkBarFromW_;

		TRootParticle lepton1FromZ_;
		TRootParticle lepton2FromZ_;

		std::vector<TRootJet> selectedJets_;
		std::vector<TRootJet> ISR_;
		std::vector<TRootJet> smDecayTopRadiation_;
		std::vector<TRootJet> fcncDecayTopRadiation_;

		Float_t wMass_;
		Float_t zMass_;
		Float_t topMass_;

		Float_t hadWMassResol_;
		Float_t hadTopMassResol_;
		Float_t lepWMassResol_;
		Float_t lepZMassResol_;
		Float_t lepTopMassResol_;
		Float_t fcncTopMassResol_;

		Float_t massChi2_;
		Float_t massChi2_SM_;
		Float_t massChi2_FCNC_;

		Bool_t verbose_;

		ClassDef (TopFCNC_Evt,2);
	};

#endif
