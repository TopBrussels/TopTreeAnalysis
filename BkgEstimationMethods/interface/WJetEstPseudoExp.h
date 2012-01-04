#ifndef WJetEstPseudoExp_h
#define WJetEstPseudoExp_h

#include <iostream>
#include <iomanip>
#include <vector>

#include "TRandom3.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"

#include "WJetEstimation.h"

using namespace std;



class WJetEstPseudoExp {

  public:
	WJetEstPseudoExp(int, WJetEstimation *);
	//WJetEstPseudoExp(int, int, double**, double**, double**, double**, double**);
	~WJetEstPseudoExp();

	void RollTheDice(bool);
	void RollTheDice(double **, bool);
	
	void GetPseudoExperiments(double***,double***,double***,double***);
	
	double GetNWError(int ijets, bool doFit=false) const; //ijets = 0-> 4 jets - ijets [0-2]
	double GetNWError(bool doFit=false) const;
	double GetNWBias(int ijets, bool doFit=false) const;
	double GetNWBias(bool doFit=false) const;
	double GetNTTError(int ijets, bool doFit=false) const; //ijets = 0-> 4 jets - ijets [0-2]
	double GetNTTError(bool doFit=false) const;
	double GetNTTBias(int ijets, bool doFit=false) const;
	double GetNTTBias(bool doFit=false) const;

	TH1F* Get_hNw( int nbjets, int njets);
	TH1F* Get_hNw( int njets);
	TH1F* Get_hNw();
	TH1F* Get_hNtt(int nbjets, int njets);
	TH1F* Get_hNtt(int njets);
	TH1F* Get_hNtt();

	void Write(TFile*, string);
	void PrintResults();
	
  private:
	int NbOfDataSets_;
	int NbOfPseudoExp_;
	int NbOfSuccessfulPE_;
	
	WJetEstimation *wjEst_;
	TRandom3* rand;
	
	double**  initN0bjets_;
	double**  initN1bjets_;
	double**  initN2bjets_;
	double**  initN3bjets_;
	
	double*** randN0bjets_; 
	double*** randN1bjets_;
	double*** randN2bjets_;
	double*** randN3bjets_;
	
	double*  multijets_est_N0_;
	double*  multijets_est_N1_;
	double*  multijets_est_N2_;
	double*  multijets_est_N3_;
	
	TH1F**   hEudsc;
	TH1F**   hEb;
	TH1F**   hNtt_Njets;
	TH1F**   hNw_Njets;
	TH1F**   hNtt_Nbjets;
	TH1F**   hNw_Nbjets;
	TH1F**   hNtotal;

	TH1F*    hNtt_Incl;
	TH1F*    hNw_Incl;

	TH1F***  hRandN0bjets;
	TH1F***  hRandN1bjets;
	TH1F***  hRandN2bjets;
	TH1F***  hRandN3bjets;
};

#endif
