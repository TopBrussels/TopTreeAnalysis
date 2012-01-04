#ifndef JESTools_h
#define JESTools_h

#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <sys/stat.h>
#include <map>

#include "TROOT.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"
#include "TCanvas.h"

// RooFit librairies
#include "RooArgSet.h"
#include "RooAddition.h"
#include "RooCategory.h"
#include "RooConstVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHist.h"
#include "RooHistPdf.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooGenericPdf.h"
#include "RooMCStudy.h"
#include "RooMinuit.h"
#include "RooPlot.h"
#include "RooFitResult.h"

using namespace std;

pair<float,std::pair<float,float> > GetExpectation(TH1F* histo, string method = "gaus");
float* correctMTopErrors(float mTopError, float DElError); // for non unity pull distribution
float* correctJESErrors(float DElError, float DEbError); // for non unity pull distribution
float** calibrateMTop(float mTop, float mTopError, float DEl, float DElError, string calibCurveName = ""); // apply calibrationCurve
float** estimateGoodCombisMVA(vector<float> binning, TH1F* goodTemplate, TH1F* badTemplate, TH1F* totalMVA, string pathforPNG, TFile* f, string title = ""); // estimate the number of good and bad jet combi's as a function of the MVA value

#endif
