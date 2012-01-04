#ifndef BTagName_CSV_h
#define BTagName_CSV_h

#include <string>

using namespace std;

class BTagName_CSV{

  public:
   BTagName_CSV();
   ~BTagName_CSV();
   
   std::string NameGiving(int TCHEbTagLoop, int NumberTCHEbTags,  int TCHPbTagLoop, int NumberTCHPbTags, int SSVHEbTagLoop, int NumberSSVHEbTags, int SSVHPbTagLoop, int NumberSSVHPbTags, int CSVbTagLoop, int NumberCSVbTags);
   std::string NameGivingPres(int TCHEbTagLoop, int NumberTCHEbTags,  int TCHPbTagLoop, int NumberTCHPbTags, int SSVHEbTagLoop, int NumberSSVHEbTags, int SSVHPbTagLoop, int NumberSSVHPbTags, int CSVbTagLoop, int NumberCSVbTags);
   
  private:
   std::string bTagFileOutput[14+14+14+14+14];
   std::string PresOutput[14+14+14+14+14];
};

#endif
