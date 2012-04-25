#ifndef TopFCNC_GenEvt_h
#define TopFCNC_GenEvt_h

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include "TLorentzVector.h"
#include "TopTreeProducer/interface/TRootJet.h"
#include "TopTreeProducer/interface/TRootMET.h"
#include "TopTreeProducer/interface/TRootMCParticle.h"
#include "TopTreeProducer/interface/TRootParticle.h"
#include "TopTreeAnalysis/MCInformation/interface/JetPartonMatching.h"
#include "TopTreeAnalysis/MCInformation/interface/ResolutionFit.h"
#include "TopTreeAnalysis/Reconstruction/interface/MEzCalculator.h"
#include "TopTreeAnalysis/Selection/interface/Selection.h"

using namespace std;

class TopFCNC_GenEvt : public TObject
{
	public: 
		// W/Z boson leptonic decay channel
		enum LeptonType {kNone, kMuon, kElec, kTau};

		TopFCNC_GenEvt() :
			TObject(),
			zLeptonicChannel_(kNone),
			wLeptonicChannel_(kNone),
			foundZ_(false)
//			isDiLeptonic_(false),
//			isTriLeptonic_(false)
			{;}
		~TopFCNC_GenEvt(){;}

		Bool_t isDiLeptonic()  const {return (wLeptonicChannel_== kNone);}
		Bool_t isTriLeptonic() const {return (wLeptonicChannel_!= kNone);}
		Bool_t isDiLeptonic(LeptonType type)                     const { return  (zLeptonicChannel_==type ? true : false); }
		Bool_t isTriLeptonic(LeptonType type1, LeptonType type2) const { return ((wLeptonicChannel_==type1 && zLeptonicChannel_==type2) ? true : false); }
		LeptonType zLeptonicChannel() const { return zLeptonicChannel_;}
		LeptonType wLeptonicChannel() const { return wLeptonicChannel_;}
	  
		const TRootMCParticle smDecayTop()    const { return smDecayTop_;}
		const TRootMCParticle W()             const { return W_;}
		const TRootMCParticle B()             const { return B_;}
		const TRootMCParticle fcncDecayTop()  const { return fcncDecayTop_;}
		const TRootMCParticle Z()             const { return Z_;}
		const TRootMCParticle Q()             const { return Q_;}

		const TRootMCParticle leptonFromW()   const { return leptonFromW_;}
		const TRootMCParticle neutrino()      const { return neutrino_;}
		const TRootMCParticle quarkFromW()    const { return quarkFromW_;}
		const TRootMCParticle quarkBarFromW() const { return quarkBarFromW_;}

		const TRootMCParticle lepton1FromZ()  const { return lepton1FromZ_;}
		const TRootMCParticle lepton2FromZ()  const { return lepton2FromZ_;}

		const TRootJet matchedB()             const { return matchedB_;}
		const TRootJet matchedQ()             const { return matchedQ_;}
		const TRootJet matchedQuarkFromW()    const { return matchedQuarkFromW_;}
		const TRootJet matchedQuarkBarFromW() const { return matchedQuarkBarFromW_;}

		const TRootParticle matchedLepton1FromZ() const { return matchedLepton1FromZ_;}
		const TRootParticle matchedLepton2FromZ() const { return matchedLepton2FromZ_;}
		const TRootParticle matchedZ()            const { return matchedZ_;}

		void ReconstructEvt(const std::vector<TRootMCParticle*> &mcParticles);
		void MatchJetsToPartons(const std::vector<TRootJet*> &jets, const int algorithm = JetPartonMatching::totalMinDist, const bool useMaxDist = true, const bool useDeltaR = true, const double maxDist = 0.3);
		void MatchLeptonsToZ(const std::vector<TRootMuon*>     &leptons, const int algorithm = JetPartonMatching::totalMinDist, const bool useMaxDist = true, const bool useDeltaR = true, const double maxDist = 0.1);
		void MatchLeptonsToZ(const std::vector<TRootElectron*> &leptons, const int algorithm = JetPartonMatching::totalMinDist, const bool useMaxDist = true, const bool useDeltaR = true, const double maxDist = 0.1);

    void FillResolutions(ResolutionFit* resLeptons, ResolutionFit* resFitLightJets, ResolutionFit* resFitBJets);
	protected:
		LeptonType zLeptonicChannel_;
		LeptonType wLeptonicChannel_;

    Bool_t foundZ_;
//		Bool_t isDiLeptonic_;
//		Bool_t isTriLeptonic_;

		TRootMCParticle smDecayTop_;
		TRootMCParticle W_;
		TRootMCParticle B_;
		TRootMCParticle fcncDecayTop_;
		TRootMCParticle Z_;
		TRootMCParticle Q_;

		TRootMCParticle leptonFromW_;
		TRootMCParticle neutrino_;
		TRootMCParticle quarkFromW_;
		TRootMCParticle quarkBarFromW_;

		TRootMCParticle lepton1FromZ_;
		TRootMCParticle lepton2FromZ_;

		TRootJet matchedB_;
		TRootJet matchedQ_;
		TRootJet matchedQuarkFromW_;
		TRootJet matchedQuarkBarFromW_;

		TRootParticle matchedLepton1FromZ_;
		TRootParticle matchedLepton2FromZ_;
		TRootParticle matchedZ_;

//		ClassDef (TopFCNC_GenEvt,1);
	};

#endif
