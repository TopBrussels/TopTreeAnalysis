#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <map>
#include <vector>

using namespace std;

typedef struct {
	
	string type;
    string sel;
	float F;
	float uF;
	float diffF;
	float udiffF;
	float eb;
	float ueb;
	float diffeb;
    float nomdiffeb;
    float unomdiffeb;
	float udiffeb;
	float xs;
	float uxs;
	float diffxs;
	float udiffxs;
    float SF;
    float diffSF;

} systematic;

float max (float a, float b) {
	if (fabs(a)>fabs(b)) return fabs(a);
	return fabs(b);
}
	
float round(float depth, float n)
{
    float d;
    int    i;
	
    /* rescale 123.45678 to 12345.678 */ 
    d = n * depth;
    /* round off: 12345.678 + 0.5 = 12346.178 -> 12346 */ 
    i = d + 0.5;
    /* restore to its original scale: 12346 -> 123.46 */
    d = (float)i / depth;
	
    return d;
}

void printSystematics(std::string type, std::map<std::string, std::vector<float> > systematics, std::map<std::string,float> bias, int depth) {
    cout <<  "----------------- Systematics on "<< type <<" -----------------" << endl;
	float sumSyst=0;
	
	for (std::map<std::string, std::vector<float> >::iterator it = systematics.begin(); it != systematics.end(); ++it) {
		
		if (it->second.size() > 2)
			cout << "Error, " << it->first << " has more than two values. This is unexpected!" << endl;
		
		else {
            
            float sys=max(it->second[0],it->second[1]);
            
            //cout << max(it->second[0],it->second[1]) << endl;
            
            //if (bias.find(it->first) != bias.end()) {
                //cout << sys << " " << bias[it->first] << endl;
            //    sys = sqrt(pow(sys,2)+pow(bias[it->first],2));
            //}
            
            //if (sys > 1./(depth*10) ) {
                if (it->second.size() ==1) (it->second).push_back(0);
			
                cout << "+> Systematic " << it->first << " = " << round(depth,sys) << endl;
			
                sumSyst+=pow(sys,2);
                
			//}
		}	
		
	}
	
	cout << endl << "+> TOTAL: " << round(depth,sqrt(sumSyst)) << endl;

    ofstream f(("systematics/tmp/tmp_total_"+type).c_str(), ios::out | ios::trunc );
    
    f << round(1000,sqrt(sumSyst)/100);
    
    f.close();

}

int main (int argc, char* argv[]) {
	
	std::map<std::string, std::vector<float> > systematics_f;
	std::map<std::string, std::vector<float> > systematics_eb;
	std::map<std::string, std::vector<float> > systematics_xs;
	std::map<std::string, std::vector<float> > systematics_SF;
    std::map<std::string, float> bias_eb;
    
	if (argc < 2) exit(1);
	
	string fileName = (string)argv[1];
	
	fstream inFile(fileName.c_str(), ios::in);
	
	vector<systematic> systematics;
	
	cout << "Reading file " << fileName << endl;
    
    float nomdiffeb = -1;
    float unomdiffeb = -1;
	
	while (!inFile.eof()) {
		
		systematic sys;
		
		inFile >> sys.type >> sys.sel >> sys.F >> sys.uF >> sys.diffF >> sys.udiffF >> sys.eb >> sys.ueb >> sys.diffeb >> sys.udiffeb >> sys.xs >> sys.uxs >> sys.diffxs >> sys.udiffxs >> sys.SF >> sys.diffSF;
        
        if (sys.type.find("nominal") != string::npos) {
            nomdiffeb=sys.eb;
            unomdiffeb=sys.ueb;
            //cout << sys.type << " " << nomdiffeb << endl;
        } else if (sys.type.find("Data") == string::npos){
            
            //cout << sys.type << endl;
            
            sys.diffF = fabs(sys.diffF);
            sys.diffeb = fabs(sys.diffeb);
            sys.diffxs = fabs(sys.diffxs);
            
            sys.SF = fabs(sys.SF)*100;
            sys.diffSF = fabs(sys.diffSF)*100;
            
            sys.nomdiffeb=nomdiffeb;
            sys.unomdiffeb=unomdiffeb;
            
            //cout << sys.type << " " << sys.nomdiffeb << endl;

            systematics.push_back(sys);
        }
 		
	}
	
	////////////////////
	// calcing syst unc
	////////////////////
	
	//for PDF
	float sumplusF = 0;
	float summinF = 0;
	float sumplusE = 0;
	float summinE = 0;
	float sumplusXS = 0;
	float summinXS = 0;
    float summinSF = 0;
    float sumplusSF = 0;
	
	//for Backgrounds
    std::map<std::string, float> map_scaleXS_eb;
    std::map<std::string, float> map_scaleXS_SF;
    std::map<std::string, float> map_scaleXS_F;
    std::map<std::string, float> map_scaleXS_XS;

	//for matching
    std::map<std::string, float> map_match_eb;
    std::map<std::string, float> map_match_SF;
    std::map<std::string, float> map_match_F;
    std::map<std::string, float> map_match_XS;
	
	//for scale
    std::map<std::string, float> map_scale_eb;
    std::map<std::string, float> map_scale_SF;
    std::map<std::string, float> map_scale_F;
    std::map<std::string, float> map_scale_XS;
    
	int nPdf=0;
	for (unsigned int s=0;s<systematics.size();s++) {
	
		systematic sys = systematics[s];
        
        //sys.diffeb = (sys.diffeb/sys.eb)*100;
        
        string title = "Method bias";
        if (systematics_eb.find(title) == systematics_eb.end()) {
            if (sys.nomdiffeb >= sys.unomdiffeb) {
                systematics_eb[title].push_back(sys.nomdiffeb);
                systematics_SF[title].push_back(sys.nomdiffeb);
            }
            else {
                systematics_eb[title].push_back(sys.unomdiffeb);
                systematics_SF[title].push_back(sys.unomdiffeb);                
            }
            //cout << "lol" << endl;
        }
        
		//pdfunc
		if (sys.type.find("PDF-Uncert") != string::npos) {
		
			nPdf++;
			/*if (sys.diffF > 0) {
				sumplusF+=sys.diffF*sys.diffF;
			}	
			if (sys.diffeb < 0) {
				summinF+=sys.diffF*sys.diffF;
			}	

			if (sys.diffeb > 0) {
				sumplusE+=sys.diffeb*sys.diffeb;
			}	
			if (sys.diffeb < 0) {
				summinE+=sys.diffeb*sys.diffeb;
			}	
			
			if (sys.diffxs > 0) {
				sumplusXS+=sys.diffxs*sys.diffxs;
			}	
			if (sys.diffxs < 0) {
				summinXS+=sys.diffxs*sys.diffxs;
			}	*/
			if (nPdf%2 == 0) {
			
				//cout << nPdf << "  down" << endl;
				sumplusF+=sys.diffF*sys.diffF;
				sumplusE+=sys.diffeb*sys.diffeb;
				sumplusXS+=sys.diffxs*sys.diffxs;
				sumplusSF+=sys.diffSF*sys.diffSF;

			} else {
			
				//cout << nPdf << "  up" << endl;
				summinF+=sys.diffF*sys.diffF;
				summinE+=sys.diffeb*sys.diffeb;
				summinXS+=sys.diffxs*sys.diffxs;
				summinSF+=sys.diffSF*sys.diffSF;
			
			}
            
            bias_eb["PDF Uncertainties"]=sys.nomdiffeb;
		}
		// the simple systematics
		else if (sys.type.find("JES") != string::npos) {
			systematics_f["Jet Energy Scale"].push_back(sys.diffF);
			systematics_eb["Jet Energy Scale"].push_back(sys.diffeb);
			systematics_xs["Jet Energy Scale"].push_back(sys.diffxs);
			systematics_SF["Jet Energy Scale"].push_back(sys.diffSF);
            bias_eb["Jet Energy Scale"]=sys.nomdiffeb;
		}
		else if (sys.type.find("JER") != string::npos) {
			systematics_f["Jet Energy Resolution"].push_back(sys.diffF);
			systematics_eb["Jet Energy Resolution"].push_back(sys.diffeb);
			systematics_xs["Jet Energy Resolution"].push_back(sys.diffxs);
			systematics_SF["Jet Energy Resolution"].push_back(sys.diffSF);
            bias_eb["Jet Energy Resolution"]=sys.nomdiffeb;
		}
		else if (sys.type.find("PileUp") != string::npos) {
			systematics_f["PileUp"].push_back(sys.diffF);
			systematics_eb["PileUp"].push_back(sys.diffeb);
			systematics_xs["PileUp"].push_back(sys.diffxs);
			systematics_SF["PileUp"].push_back(sys.diffSF);
            bias_eb["PileUp"]=sys.nomdiffeb;
		}
		else if (sys.type.find("UE") != string::npos) {
			systematics_f["Underlying event"].push_back(sys.diffF);
			systematics_eb["Underlying event"].push_back(sys.diffeb);
			systematics_xs["Underlying event"].push_back(sys.diffxs);
			systematics_SF["Underlying event"].push_back(sys.diffSF);
            bias_eb["Underlying event"]=sys.nomdiffeb;
		}
		else if (sys.type.find("ISRFSR") != string::npos) {
			systematics_f["ISRFSR"].push_back(sys.diffF);
			systematics_eb["ISRFSR"].push_back(sys.diffeb);
			systematics_xs["ISRFSR"].push_back(sys.diffxs);
			systematics_SF["ISRFSR"].push_back(sys.diffSF);
            bias_eb["ISRFSR"]=sys.nomdiffeb;
		}
		
		else if (sys.type.find("RightReg") != string::npos) {
			systematics_f["Right Region definition"].push_back(sys.diffF);
			systematics_eb["Right Region definition"].push_back(sys.diffeb);
			systematics_xs["Right Region definition"].push_back(sys.diffxs);
			systematics_SF["Right Region definition"].push_back(sys.diffSF);
            bias_eb["Right Region definition"]=sys.nomdiffeb;
		}

        else if (sys.type.find("TopMass") != string::npos) { // divide by 9 to get the 1GeV/c2 effect
			systematics_f["Top Quark Mass"].push_back(sys.diffF/9.);
			systematics_eb["Top Quark Mass"].push_back(sys.diffeb/9.);
			systematics_xs["Top Quark Mass"].push_back(sys.diffxs/9.);
			systematics_SF["Top Quark Mass"].push_back(sys.diffSF/9.);
            bias_eb["Top Quark Mass"]=sys.nomdiffeb;
		}

		else if (sys.type.find("SF") != string::npos) {
			systematics_xs["Muon ID/Eff and Trigger eff SF"].push_back(sys.diffeb);
            //bias_eb["M"]=sys.nomdiffeb;
		}
		
		
		// background scaling
		else if (sys.type.find("XSDown") != string::npos || sys.type.find("XSUp") != string::npos) {
            
            if (map_scaleXS_eb.find(sys.sel) == map_scaleXS_eb.end()) {
                map_scaleXS_eb[sys.sel]=sys.diffeb;
                map_scaleXS_SF[sys.sel]=sys.diffSF;
                map_scaleXS_F[sys.sel]=sys.diffF;
                map_scaleXS_XS[sys.sel]=sys.diffxs;
            } else {
                map_scaleXS_eb[sys.sel]=max(map_scaleXS_eb[sys.sel],sys.diffeb);
                map_scaleXS_SF[sys.sel]=max(map_scaleXS_SF[sys.sel],sys.diffSF);
                map_scaleXS_F[sys.sel]=max(map_scaleXS_F[sys.sel],sys.diffF);
                map_scaleXS_XS[sys.sel]=max(map_scaleXS_XS[sys.sel],sys.diffxs);
                
            }
            bias_eb["Background Composition"]=sys.nomdiffeb;
            
		}
        
        // Q2 scaling
		else if (sys.type.find("ScaleDown") != string::npos || sys.type.find("ScaleUp") != string::npos) {
            
            if (map_scale_eb.find(sys.sel) == map_scale_eb.end()) {
                map_scale_eb[sys.sel]=sys.diffeb;
                map_scale_SF[sys.sel]=sys.diffSF;
                map_scale_F[sys.sel]=sys.diffF;
                map_scale_XS[sys.sel]=sys.diffxs;
            } else {
                map_scale_eb[sys.sel]=max(map_scale_eb[sys.sel],sys.diffeb);
                map_scale_SF[sys.sel]=max(map_scale_SF[sys.sel],sys.diffSF);
                map_scale_F[sys.sel]=max(map_scale_F[sys.sel],sys.diffF);
                map_scale_XS[sys.sel]=max(map_scale_XS[sys.sel],sys.diffxs);
                
            }
            bias_eb["Factorisation Scale"]=sys.nomdiffeb;
            
		}
        
        // ME-PS matching
		else if (sys.type.find("MatchingDown") != string::npos || sys.type.find("MatchingUp") != string::npos) {
            
            if (map_match_eb.find(sys.sel) == map_match_eb.end()) {
                map_match_eb[sys.sel]=sys.diffeb;
                map_match_SF[sys.sel]=sys.diffSF;
                map_match_F[sys.sel]=sys.diffF;
                map_match_XS[sys.sel]=sys.diffxs;
            } else {
                map_match_eb[sys.sel]=max(map_match_eb[sys.sel],sys.diffeb);
                map_match_SF[sys.sel]=max(map_match_SF[sys.sel],sys.diffSF);
                map_match_F[sys.sel]=max(map_match_F[sys.sel],sys.diffF);
                map_match_XS[sys.sel]=max(map_match_XS[sys.sel],sys.diffxs);
                
            }
            bias_eb["ME-PS Matching Threshold"]=sys.nomdiffeb;
            
		}
    }
    
    // BackGround composition
    
    float bgComp[4]={-1,-1,-1,-1};
    for (std::map<std::string, float>::iterator it = map_scaleXS_eb.begin(); it != map_scaleXS_eb.end(); ++it) {
        cout << "BGComp - EB - " << it->first << " " << it->second << "^2 * " << bias_eb["Background Composition"] << "^2 = " << sqrt(pow(it->second,2)+pow(bias_eb["Background Composition"],2)) << endl;
        cout << "BGComp - SF - " << it->first << " " << map_scaleXS_SF[it->first] << endl;
        cout << endl;
        if (bgComp[0] == -1) bgComp[0] = pow(map_scaleXS_F[it->first],2); else bgComp[0] += pow(map_scaleXS_F[it->first],2);
        if (bgComp[1] == -1) bgComp[1] = pow(map_scaleXS_eb[it->first],2); else bgComp[1] += pow(map_scaleXS_eb[it->first],2);
        if (bgComp[2] == -1) bgComp[2] = pow(map_scaleXS_SF[it->first],2); else bgComp[2] += pow(map_scaleXS_SF[it->first],2);
        if (bgComp[3] == -1) bgComp[3] = pow(map_scaleXS_XS[it->first],2); else bgComp[3] += pow(map_scaleXS_XS[it->first],2);
    }
	if (bgComp[0] != -1) systematics_f["Background Composition"].push_back(sqrt(bgComp[0]));
	if (bgComp[1] != -1) systematics_eb["Background Composition"].push_back(sqrt(bgComp[1]));
	if (bgComp[2] != -1) systematics_SF["Background Composition"].push_back(sqrt(bgComp[2]));
	if (bgComp[3] != -1) systematics_xs["Background Composition"].push_back(sqrt(bgComp[3]));

	// Q2 scaling composition
    
    float scale[4]={-1,-1,-1,-1};
    for (std::map<std::string, float>::iterator it = map_scale_eb.begin(); it != map_scale_eb.end(); ++it) {
        cout << "Q2 - EB - " << it->first << " " << it->second << "^2 * " << bias_eb["Factorisation Scale"] << "^2 = " << sqrt(pow(it->second,2)+pow(bias_eb["Factorisation Scale"],2)) << endl;
        cout << "Q2 - SF - " << it->first << " " << map_scale_SF[it->first] << endl;
        cout << endl;
        //cout << it->first << " " << it->second << endl;
        if (scale[0] == -1) scale[0] = pow(map_scale_F[it->first],2); else scale[0] += pow(map_scale_F[it->first],2);
        if (scale[1] == -1) scale[1] = pow(map_scale_eb[it->first],2); else scale[1] += pow(map_scale_eb[it->first],2);
        if (scale[2] == -1) scale[2] = pow(map_scale_SF[it->first],2); else scale[2] += pow(map_scale_SF[it->first],2);
        if (scale[3] == -1) scale[3] = pow(map_scale_XS[it->first],2); else scale[3] += pow(map_scale_XS[it->first],2);
    }
	if (scale[0] != -1) systematics_f["Factorisation Scale"].push_back(sqrt(scale[0]));
	if (scale[1] != -1) systematics_eb["Factorisation Scale"].push_back(sqrt(scale[1]));
	if (scale[2] != -1) systematics_SF["Factorisation Scale"].push_back(sqrt(scale[2]));
	if (scale[3] != -1) systematics_xs["Factorisation Scale"].push_back(sqrt(scale[3]));

    // ME-PS matching composition
    
    float match[4]={-1,-1,-1,-1};
    for (std::map<std::string, float>::iterator it = map_match_eb.begin(); it != map_match_eb.end(); ++it) {
        cout << "MEPS - EB - " << it->first << " " << it->second << "^2 * " << bias_eb["ME-PS Matching Threshold"] << "^2 = " << sqrt(pow(it->second,2)+pow(bias_eb["ME-PS Matching Threshold"],2)) << endl;
        cout << "MEPS - SF - " << it->first << " " << map_match_SF[it->first] << endl;
        cout << endl;
        //cout << it->first << " " << it->second << endl;
        if (match[0] == -1) match[0] = pow(map_match_F[it->first],2); else match[0] += pow(map_match_F[it->first],2);
        if (match[1] == -1) match[1] = pow(map_match_eb[it->first],2); else match[1] += pow(map_match_eb[it->first],2);
        if (match[2] == -1) match[2] = pow(map_match_SF[it->first],2); else match[2] += pow(map_match_SF[it->first],2);
        if (match[3] == -1) match[3] = pow(map_match_XS[it->first],2); else match[3] += pow(map_match_XS[it->first],2);
    }
	if (match[0] != -1) systematics_f["ME-PS Matching Threshold"].push_back(sqrt(match[0]));
	if (match[1] != -1) systematics_eb["ME-PS Matching Threshold"].push_back(sqrt(match[1]));
	if (match[2] != -1) systematics_SF["ME-PS Matching Threshold"].push_back(sqrt(match[2]));
	if (match[3] != -1) systematics_xs["ME-PS Matching Threshold"].push_back(sqrt(match[3]));


 	/*cout << "Pdf uncertainty Eb: " << max(sqrt(sumplusE),sqrt(summinE)) << endl;
 	cout << "Pdf uncertainty XS: " << max(sqrt(sumplusXS),sqrt(summinXS)) << endl;
	*/
	if (nPdf > 0) systematics_f["PDF Uncertainties"].push_back(max(sqrt(sumplusF),sqrt(summinF)));
	if (nPdf > 0) systematics_eb["PDF Uncertainties"].push_back(max(sqrt(sumplusE),sqrt(summinE)));
	if (nPdf > 0) systematics_xs["PDF Uncertainties"].push_back(max(sqrt(sumplusXS),sqrt(summinXS)));
	if (nPdf > 0) systematics_SF["PDF Uncertainties"].push_back(max(sqrt(sumplusSF),sqrt(summinSF)));
	
	std::map<std::string,float> empty;
    
	printSystematics("F_CS",systematics_f,empty,100);	
	printSystematics("E_b",systematics_eb,empty,10);	
	printSystematics("SF",systematics_SF,empty,10);	
	printSystematics("XS",systematics_xs,empty,10);	
	
	cout << endl;
	
	return 0;
	
}
