#ifndef VJetEstPseudoExp_h
#define VJetEstPseudoExp_h

#include <iostream>
//#include <iomanip>
#include <vector>
#include <numeric>
#include <algorithm>

//#include "TRandom3.h"
#include "Math/Random.h"
#include "Math/GSLRndmEngines.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"

#include "RooMsgService.h"

#include "VJetEstimation.h"

using namespace std;
using namespace RooFit ;

class VJetEstPseudoExp {

  public:
	VJetEstPseudoExp(int, VJetEstimation *vjetEst);
	~VJetEstPseudoExp();

	void SetCondProb(double***** condProb) {condProb_ = condProb;};
	void RollTheDice(string &method, string &option, vector<int> SetFixedVarIdx, bool doMinos, bool UseUnBinMLE, bool UseMJLE, bool Print);
	/** Randomize the numbers of events with {0,1,2,3} b-jets and each time perform the VJetEstimation */

	void SetInitNbjets(double**** initNbjets);
	void ResetInitialValues();
	/** Reset the numbers of events with {0,1,2,3} b-jets to their initial values*/
	void GetPseudoExperiments(double****,double****,double****,double****);
	/** Fill the array pointers with the results of the pseudo-exp */
	
	double GetNvVar( int wp, bool doFit, double xlow, double xhigh) const;
	double GetNvVar( int wp, int njets, bool doFit, double xlow, double xhigh) const; 
	double GetNvVar( int wp, int njets, int ibjets, bool doFit, double xlow, double xhigh) const; 
	double GetNttVar(int wp, bool doFit, double xlow, double xhigh) const;
	double GetNttVar(int wp, int njets, bool doFit, double xlow, double xhigh) const; 
	double GetNttVar(int wp, int njets, int ibjets, bool doFit, double xlow, double xhigh) const;

	double GetNvBias( int wp, bool doFit, double xlow, double xhigh) const;
	double GetNvBias( int wp, int njets, bool doFit, double xlow, double xhigh) const;
	double GetNvBias( int wp, int njets, int ibjets, bool doFit, double xlow, double xhigh) const;
	double GetNttBias(int wp, bool doFit, double xlow, double xhigh) const;
	double GetNttBias(int wp, int njets, bool doFit, double xlow, double xhigh) const;
	double GetNttBias(int wp, int njets, int ibjets, bool doFit, double xlow, double xhigh) const;

	TH1F* Get_hNv( int wp, int nbjets, int njets);
	TH1F* Get_hNv( int wp, int njets);
	TH1F* Get_hNv( int wp);
	TH1F* Get_hDiffNv( int wp, int nbjets, int njets);
	TH1F* Get_hDiffNv( int wp, int njets);
	TH1F* Get_hDiffNv( int wp);
	TH1F* Get_hNtt(int wp, int nbjets, int njets);
	TH1F* Get_hNtt(int wp, int njets);
	TH1F* Get_hNtt(int wp);
	TH1F* Get_hDiffNtt(int wp, int nbjets, int njets);
	TH1F* Get_hDiffNtt(int wp, int njets);
	TH1F* Get_hDiffNtt(int wp);

	void Write(TFile*, string);
	void PrintResults();
	void PrintResults_latex(ofstream &ofile);

  private:
	int NbOfDataSets_;
	int NbOfJetsBins_;
	int NbOfBJetsBins_;
	int NbOfBtagWorkingPoint_;
	int NbOfPseudoExp_;
	int NbOfSuccessfulPE_;
	int NbOfRealPE_;
	
	VJetEstimation *vjEst_;
	//TRandom3* rand;
	ROOT::Math::Random<ROOT::Math::GSLRngMT> *rand_;
	double *****condProb_; /** Store the conditional probabilities needed to randomize the number of selected events with 0,1,2 and 3
	b-jets, in case of use of several b-tagging working points (per jet multiplicity and per dataset) */

	/// Initial numbers of events for each b-jet multiplicity
	/// -- per b-tagging working point, jet multiplicity and per dataset
	double****  initNbjets_;
	
	/** Randomized numbers of events with 0,1,2 and 3 b-jets
	-- per pseudo-exp, per b-tagging working points, per jet multiplicity and per dataset */
	double**** randN0bjets_; 
	double**** randN1bjets_;
	double**** randN2bjets_;
	double**** randN3bjets_;

	/** Histograms for the efficiency estimation
	-- per b-tagging working point and jet multiplicity */
	TH1F*** hEudsc_;
	TH1F*** hEuds_;
	TH1F*** hEb_;
	TH1F*** hPullEudsc_;
	TH1F*** hPullEuds_;
	TH1F*** hPullEb_;
	TH1F*** hDiffEudsc_;
	TH1F*** hDiffEuds_;
	TH1F*** hDiffEb_;
	/** Histograms for the parameter estimation
	-- per b-tagging working point, per jet and b-jet multiplicity */
	TH1F****  hNtt_;
	TH1F****  hNv_;
	TH1F****  hNvb_;
	TH1F****  hPullNtt_;
	TH1F****  hPullNv_;
	TH1F****  hPullNvb_;
	TH1F****  hDiffNtt_;
	TH1F****  hDiffNv_;
	TH1F****  hDiffNvb_;
	TH2F****  hDiffNtt_vs_MinLL_;
	TH2F****  hDiffNv_vs_MinLL_;
	TH2F****  hDiffNvb_vs_MinLL_;

	TH1F****  hRandN0bjets_;
	TH1F****  hRandN1bjets_;
	TH1F****  hRandN2bjets_;
	TH1F****  hRandN3bjets_;

	TH3F****  hRandNbjets_vs_MinLL_;
};

#endif
