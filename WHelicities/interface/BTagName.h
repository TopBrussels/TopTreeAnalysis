#ifndef BTagName_h
#define BTagName_h

#include <string>

#include <iostream>

using namespace std;

class BTagName{

  public:
   BTagName();
   ~BTagName();
   
   std::string NameGiving(int TCHEbTagLoop, int NumberTCHEbTags,  int TCHPbTagLoop, int NumberTCHPbTags, int SSVHEbTagLoop, int NumberSSVHEbTags, int SSVHPbTagLoop, int NumberSSVHPbTags, int CSVbTagLoop, int NumberCSVbTags,int JPbTagLoop, int NumberJPbTags,int JBPbTagLoop, int NumberJBPbTags);
   std::string NameGivingPres(int TCHEbTagLoop, int NumberTCHEbTags,  int TCHPbTagLoop, int NumberTCHPbTags, int SSVHEbTagLoop, int NumberSSVHEbTags, int SSVHPbTagLoop, int NumberSSVHPbTags, int CSVbTagLoop, int NumberCSVbTags,int JPbTagLoop, int NumberJPbTags,int JBPbTagLoop, int NumberJBPbTags);
   
  private:
   std::string bTagFileOutput[14+14+14+14+14+14+14];
   std::string PresOutput[14+14+14+14+14+14+14];
};

#endif
