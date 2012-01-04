#ifndef Monster_h
#define Monster_h

#include <string>
#include <vector>
#include <map>
#include <iostream>

#include "Rtypes.h"
#include "TObject.h"
#include "TVector3.h"
#include "TRef.h"
#include "TH2F.h"

using namespace std;

class Monster : public TObject
{
 public:
	
  Monster():
    eventID_(0)
    ,runID_(0)
    ,lumiBlockID_(0)
    ,maxMVA_(0)
    ,eventWeight_(0)
    ,fitResults_()
   {;}
  
  Monster(unsigned int eventID, unsigned int runID, unsigned int lumiBlockID, float maxMVA, float eventWeight, TH2F fitResults):
    eventID_(eventID)
    ,runID_(runID)
    ,lumiBlockID_(lumiBlockID)
    ,maxMVA_(maxMVA)
    ,eventWeight_(eventWeight)
    ,fitResults_(fitResults)
   {;}
  
  unsigned int eventID() const { return eventID_; }
  unsigned int runID() const { return runID_; }
  unsigned int lumiBlockID() const { return lumiBlockID_; }
  float maxMVA() const { return maxMVA_; }
  float eventWeight() const { return eventWeight_; }
  TH2F fitResults() const { return fitResults_; }
  
  ~Monster() {;}
  
 protected:
  
  unsigned int eventID_;
  unsigned int runID_;
  unsigned int lumiBlockID_;
  float maxMVA_;
  float eventWeight_;
  TH2F fitResults_;
  
  ClassDef (Monster,2);
};

#endif


