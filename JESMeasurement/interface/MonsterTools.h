#ifndef MonsterTools_h
#define MonsterTools_h

#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <cmath>

#include "TROOT.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"

#include "TopTreeAnalysis/JESMeasurement/interface/LightMonster.h"
#include "TopTreeAnalysis/JESMeasurement/interface/Monster.h"
#include "TopTreeAnalysis/JESMeasurement/interface/ExpCorrCalculator.h"

using namespace std;

float maxProb(Monster* monster, int JetCombi=0);
int nHoles(Monster* monster, bool topMass, int JetCombi=0);
float minDisHoleMax(TH2F* monster);
float mTopMax(TH2F* monster);
float probNoCorr(Monster* monster, TH2F* dummyMonster, bool topMass, int jetCombi=0);
bool hadrJetsMVAMatched(LightMonster* monster, int jetCombi=0);
bool hadrJetsMVAMatched(Monster* monster, int jetCombi=0);
bool wJetsMVAMatched(Monster* monster, int jetCombi=0);
int nWJetsMVAMatched(Monster* monster, int jetCombi=0);
int nHadrBJetsMVAMatched(Monster* monster, int jetCombi=0);
TH2F* histoMonster(Monster* monster, TH2F* dummyMonster, int jetCombi=0);

vector<double> CalculateParabola(TH1D* parabola, int size, bool fit);
vector< vector<float> > FitSuperMonster(TH2F* superMonster, string pathforPNG, TFile* f, string title, bool topMass);

void fillExpCorr(Monster *monster, ExpCorrCalculator *expCorr);

#endif
