#ifndef BTagJetSelection_CSV_h
#define BTagJetSelection_CSV_h

#include <string>
#include <vector>
#include <iostream>

using namespace std;

class BTagJetSelection_CSV{

  public:
   BTagJetSelection_CSV();
   ~BTagJetSelection_CSV();
   
   int HighestProbSelection(int bTagLoop, int ConsideredBTagger, vector<float> KinFitProb, vector<float> MlbProb, vector<float> btagTCHE, vector<float> btagTCHP, vector<float> btagSSVHE, vector<float> btagSSVHP, vector<float> btagCSV);
   
   private:

};

#endif
