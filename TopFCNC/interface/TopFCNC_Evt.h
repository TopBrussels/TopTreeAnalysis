#ifndef TopFCNC_Evt_h
#define TopFCNC_Evt_h

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include "TLorentzVector.h"

using namespace std;

class TopFCNC_Evt : public TObject
{
	public: 
		// W/Z boson leptonic decay channel
		enum LeptonType {kNone, kMuon, kElec, kTau};

		TopFCNC_Evt() :
			TObject(),
			WLeptonicChannel_(kNone),
			ZLeptonicChannel_(kNone),
			isDiLeptonic_(false),
			isTriLeptonic_(false)
			{;}

		TopFCNC_Evt(const TopFCNC_Evt& evt) :
			TObject(evt),
			WLeptonicChannel_(evt.WLeptonicChannel_),
			ZLeptonicChannel_(evt.ZLeptonicChannel_),
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
			ISR_(evt.ISR_),
			smDecayTopRadiation_(evt.smDecayTopRadiation_),
			fcncDecayTopRadiation_(evt.fcncDecayTopRadiation_)
			{;}
	
		~TopFCNC_Evt(){;}

		Bool_t isDiLeptonic() const {return isDiLeptonic_;}
		Bool_t isTriLeptonic() const {return isTriLeptonic_;}
		LeptonType WLeptonicChannel() const { return WLeptonicChannel_;}
		LeptonType ZLeptonicChannel() const { return ZLeptonicChannel_;}
		Bool_t isDiLeptonic(LeptonType type) const { return (ZLeptonicChannel()==type ? true : false); }
		Bool_t isTriLeptonic(LeptonType type1, LeptonType type2) const { return ((WLeptonicChannel()==type1 && ZLeptonicChannel()==type2) ? true : false); }
	  
		const TLorentzVector leptonFromW()   const { return leptonFromW_;}
		const TLorentzVector neutrino()      const { return neutrino_;}
		const TLorentzVector quarkFromW()    const { return quarkFromW_;}
		const TLorentzVector quarkBarFromW() const { return quarkBarFromW_;}

		const TLorentzVector lepton1FromZ()  const { return lepton1FromZ_;}
		const TLorentzVector lepton2FromZ()  const { return lepton2FromZ_;}

		const TLorentzVector smDecayTop()    const { return smDecayTop_;}
		const TLorentzVector W()             const {return W_;}
		const TLorentzVector B()             const {return B_;}
		const TLorentzVector fcncDecayTop()  const {return fcncDecayTop_;}
		const TLorentzVector Z()             const {return Z_;}
		const TLorentzVector Q()             const {return Q_;}

		const std::vector<TLorentzVector> smDecayTopRadiation() const { return smDecayTopRadiation_;}
		const std::vector<TLorentzVector> fcncDecayTopRadiation() const { return fcncDecayTopRadiation_;}
		const std::vector<TLorentzVector> ISR()const { return ISR_;}
 

		void SetDiLeptonicChannel(LeptonType type)
		{
			isDiLeptonic_ = true;
			isTriLeptonic_= false;
			ZLeptonicChannel_ = type;
		}
		void SetTriLeptonicChannel(LeptonType type1, LeptonType type2)
		{
			isDiLeptonic_  = false;
			isTriLeptonic_ = true;
			WLeptonicChannel_ = type1;
			ZLeptonicChannel_ = type2;
		}

		void SetTLorentzVector(TLorentzVector &smDecayTop, TLorentzVector &W, TLorentzVector &B, TLorentzVector &fcncDecayTop, TLorentzVector &Z, TLorentzVector &Q, TLorentzVector &leptonFromW, TLorentzVector &neutrino, TLorentzVector &quarkFromW, TLorentzVector &quarkBarFromW, TLorentzVector &lepton1FromZ, TLorentzVector &lepton2FromZ)
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

		void SetsmDecayTop(TLorentzVector smDecayTop){ smDecayTop_ = smDecayTop; }
		void SetW(TLorentzVector W){ W_ = W; }
		void SetB(TLorentzVector B){ B_ = B; }
		void SetfcncDecayTop(TLorentzVector fcncDecayTop){ fcncDecayTop_ = fcncDecayTop; }
		void SetZ(TLorentzVector Z){ Z_ = Z; }
		void SetQ(TLorentzVector Q){ Q_ = Q; }

		void SetLeptonFromW(TLorentzVector leptonFromW){ leptonFromW_ = leptonFromW; }
		void SetNeutrino(TLorentzVector neutrino){ neutrino_ = neutrino; }
		void SetQuarkFromW(TLorentzVector quarkFromW){ quarkFromW_ = quarkFromW; }
		void SetQuarkBarFromW(TLorentzVector quarkBarFromW){ quarkBarFromW_ = quarkBarFromW; }

		void SetLepton1FromZ(TLorentzVector lepton1FromZ){ lepton1FromZ_ = lepton1FromZ; }
		void SetLepton2FromZ(TLorentzVector lepton2FromZ){ lepton2FromZ_ = lepton2FromZ; }

		void SetRadiation(std::vector<TLorentzVector> smDecayTopRadiation, std::vector<TLorentzVector> fcncDecayTopRadiation, std::vector<TLorentzVector> ISR)
		{
			smDecayTopRadiation_   = smDecayTopRadiation; 
			fcncDecayTopRadiation_ = fcncDecayTopRadiation;
			ISR_ = ISR; 
		}
		void SetsmDecayTopRadiation(std::vector<TLorentzVector> smDecayTopRadiation){ smDecayTopRadiation_ = smDecayTopRadiation; }
		void SetfcncDecayTopRadiation(std::vector<TLorentzVector> fcncDecayTopRadiation){ fcncDecayTopRadiation_ = fcncDecayTopRadiation; }
		void SetISR(std::vector<TLorentzVector> ISR){ISR_ = ISR;}

		void ReconstructEvt()
		{
			if(isDiLeptonic_ == isTriLeptonic_)
			{
				cerr << "Specify if top FCNC event candidate is di- or tri-leptonic before event reconstruction" << endl;
				exit(1);
			}
			if(isDiLeptonic_) W_ = quarkFromW_ + quarkBarFromW_ ;
			else              W_ = leptonFromW_+ neutrino_ ;
			smDecayTop_ = B_ + W_ ;
			Z_ = lepton1FromZ_ + lepton2FromZ_ ;
			fcncDecayTop_ = Q_ + Z_ ;
		}

		virtual TString typeName() const { return "TopFCNC_Evt"; }

		friend std::ostream& operator<< (std::ostream& stream, const TopFCNC_Evt& fcncEvt)
		{
			stream << "TopFCNC_Evt - is ";
			if(fcncEvt.isDiLeptonic()) {
				stream <<" DiLeptonic (";
				if(fcncEvt.ZLeptonicChannel()==0) stream <<" no channel )";
				else if(fcncEvt.ZLeptonicChannel()==1) stream <<" muonic )";
				else if(fcncEvt.ZLeptonicChannel()==2) stream <<" electronic )";
				else if(fcncEvt.ZLeptonicChannel()==3) stream <<" tauonic )";
			}
			else if(fcncEvt.isTriLeptonic()){
				stream <<" TriLeptonic (";
				if(fcncEvt.WLeptonicChannel()==0) stream << "W : no channel and ";
				else if(fcncEvt.WLeptonicChannel()==1) stream << "W : muonic and ";
				else if(fcncEvt.WLeptonicChannel()==2) stream << "W : electronic and ";
				else if(fcncEvt.WLeptonicChannel()==3) stream << "W : tauonic and ";

				if(fcncEvt.ZLeptonicChannel()==0) stream << "Z : no channel )";
				else if(fcncEvt.ZLeptonicChannel()==1) stream << "Z : muonic )";
				else if(fcncEvt.ZLeptonicChannel()==2) stream << "Z : electronic )";
				else if(fcncEvt.ZLeptonicChannel()==3) stream << "Z : tauonic )";
			}
			stream << std::endl;
			stream << "Nof ISR: "<< fcncEvt.ISR().size() << std::endl;
			stream << "Nof Top radiations: "<< fcncEvt.smDecayTopRadiation().size() + fcncEvt.fcncDecayTopRadiation().size() << std::endl;
			//stream << "TRootGenEvent - Charge=" << setw(2) << jet.charge() << " (Et,eta,phi)=("<< setw(10) << jet.Et() <<","<< setw(10) << jet.Eta() <<","<< setw(10) << jet.Phi() << ")"
			return stream;
		};


//	private:
	protected:
		LeptonType WLeptonicChannel_;
		LeptonType ZLeptonicChannel_;

		Bool_t isDiLeptonic_;
		Bool_t isTriLeptonic_;

		TLorentzVector smDecayTop_;
		TLorentzVector W_;
		TLorentzVector B_;
		TLorentzVector fcncDecayTop_;
		TLorentzVector Z_;
		TLorentzVector Q_;

		TLorentzVector leptonFromW_;
		TLorentzVector neutrino_;
		TLorentzVector quarkFromW_;
		TLorentzVector quarkBarFromW_;

		TLorentzVector lepton1FromZ_;
		TLorentzVector lepton2FromZ_;

		std::vector<TLorentzVector> ISR_;
		std::vector<TLorentzVector> smDecayTopRadiation_;
		std::vector<TLorentzVector> fcncDecayTopRadiation_;

		ClassDef (TopFCNC_Evt,1);
	};

#endif
