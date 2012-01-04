#ifndef BTagName_h
#define BTagName_h

#include <string>

using namespace std;

class BTagName{

  public:
   BTagName();
   ~BTagName();
   
   std::string NameGiving(int TCHEbTagLoop, int NumberTCHEbTags,  int TCHPbTagLoop, int NumberTCHPbTags, int SSVHEbTagLoop, int NumberSSVHEbTags, int SSVHPbTagLoop, int NumberSSVHPbTags);
   std::string NameGivingPres(int TCHEbTagLoop, int NumberTCHEbTags,  int TCHPbTagLoop, int NumberTCHPbTags, int SSVHEbTagLoop, int NumberSSVHEbTags, int SSVHPbTagLoop, int NumberSSVHPbTags);
   
  private:
   std::string bTagFileOutput[14][14][14][14];
   std::string PresOutput[14][14][14][14];
};

#endif
