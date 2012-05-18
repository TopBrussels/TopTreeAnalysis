#ifndef BTagJetSelection_h
#define BTagJetSelection_h

#include <string>
#include <vector>
#include <iostream>

using namespace std;

class BTagJetSelection{

  public:
   BTagJetSelection();
   ~BTagJetSelection();
   
   int HighestProbSelection(int bTagLoop, int ConsideredBTagger, vector<float> KinFitProb, vector<float> btagTCHE, vector<float> btagTCHP, vector<float> btagSSVHE, vector<float> btagSSVHP, vector<float> btagCSV);
   
   private:

};

#endif
