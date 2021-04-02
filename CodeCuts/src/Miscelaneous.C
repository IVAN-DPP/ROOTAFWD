/*******************************************/
// Author: Edwin Munevar Espitia
// Author: Ivan Dario Piernagorda Peña     
// Date:   15/02/2020                  
// Title:  Photon Polarization Tables                              
/*******************************************/

#ifndef MISCELANEOUS_H
#define MISCELANEOUS_H

#include "../include/Miscelaneous.h"
#include "../include/Libraries.h"


#define PI 3.141592654

using namespace std;

//Required by trip files functions
int dcustDebugFlag = 0;
vector<vector<int> > dcustTripIntervals; //1st & last num is INCLUSIVE (so if 1,4; covers 1,2,3,AND 4)


///////////////////////////////////
//  POLARIZATION TABLES          //
///////////////////////////////////

enum {
  E_ID    = 0,
  ENERGY  = 1,
  ENH     = 2,
  ENHERR  = 3,
  ENHFIT  = 4,
  PFIT    = 5,
  PCOR    = 6,
  PCORERR = 7,
  PSMOOTH = 8,
  EDGE    = 9,
};

enum {
  PARA1=0,
  PERP1=1,
  PARA2=2,
  PERP2=3,
  PARA3=4,
  PERP3=5
};


double polTable[73][500][385][10];  //where its [plane][edge][E_id][field]
int polTableN[73]={};               //No of entries for para and perp
char polFirstLines[73][500][250];   //to keep the 1st lines if needed (info from files)

int edgeEventLow[5000];            //hold the current table of edge positions for event ranges
int edgeEventHigh[5000];
double edgeEventEdge[5000];
int edgeEventPlane[5000];
int edgeEventN;
int edgeIndex=0;
int lastEdgeEvent=0;
double lastCohEdge=0.0;
int lastCohPlane=-1;


void ListFilesAtDir(string DirPathName,
		    vector<string> &ListFiles){
  DIR *d;
  struct dirent *dir;
  string File;
  if ( (d = opendir(DirPathName.c_str())) != NULL ){  
    // Do the push back in a vector ListFiles
    while ((dir = readdir(d)) != NULL){
      File=dir->d_name;
      if(File != "." && File!= "..")
	ListFiles.push_back(DirPathName+"/"+File);
    }
    closedir(d);
  } 
}

void ListFilesAtDir(string DirPathName,
		    string Extention,
		    vector<string> &ListFiles){
  
  DIR *d;
  char *p1,*p2;
  int ret;
  struct dirent *dir;
  string File;
  if ( (d = opendir(DirPathName.c_str())) != NULL ){  
    // Do the push back in a vector ListFiles
    while ((dir = readdir(d)) != NULL){
      File.clear();
      p1=strtok(dir->d_name,".");
      p2=strtok(NULL,".");
      
      if(p2!=NULL){
	ret=strcmp(p2,Extention.c_str());
	if(ret==0){
	  File.insert(0,DirPathName+"/");
	  File.insert(File.size(),p1);
	  File+="."+Extention;
	  ListFiles.push_back(File);
	}
      }
    }
    closedir(d);
  }
  
}

int LoadPolTable(int plane, char *PolTableList, map<vector<float>,int> &keysPlane){
  FILE *fplist,*fpfile;              //file pointers for filelist and each file therein
  char lline[250];                 //for reading in lines from filelist
  char fline[250];                 //for reading in lines from file
  char filename[250];              //file
  int  fcount=0;                   //counter for no of files read in
  int  chancount=0;                //counter for no of channels read in
  int  eid=0;
  //double edge=0.0;
  
  if((fplist=fopen(PolTableList,"r"))==NULL){ //open filelist
    cout << "Error Couldn't open file: " << PolTableList << endl;
    return -1;
  }
  cout << "List of Files: " << PolTableList << endl;
  fcount=0; 
  //for each file in the list
  while(fgets(lline,240,fplist) != NULL){
    if((lline[0] == '*')||(lline[0] == '#')) continue; //skip comments
    sscanf(lline,"%s",filename);                       //read in filename
    

    if((fpfile=fopen(filename,"r"))==NULL){              //open file
      cout << "Error Couldn't open file: " << filename << endl;
      return -1;
    }
    
    if(fgets(polFirstLines[plane][polTableN[plane]],240,fpfile) == NULL)
      perror("There isn't file"); //save the 1st line
    sscanf(strrchr(filename,'_')+1,"%lg",&polTable[plane][fcount][0][EDGE]);

    
    //File name:
    //Beam_4726_CohEdge_1.7_PARA_1687.0.fitted2
    //      ^            ^
    //      |            |
    //   Electron      Coh Edge
    
    string EBeam = filename;    float iEBeam = 0;
    string EPhot = filename;    float iEPhot = 0;
    string Pol   = filename;    float iPol = 0;
    
    EBeam = EBeam.substr(EBeam.find_first_of("_")+1,4);     //Electron energy from table name
    EPhot = EPhot.substr(EPhot.find_first_of(".")-1,3);     //Photon energy from table name
    Pol   = Pol.substr(Pol.find_first_of("P"),4);           //PARA or PERP string

    // Beam energy of file name, if you need other energy put here
    // to find the energy, see beam_en branch in the files .root
    switch(stoi(EBeam.c_str())){
    case 4199: iEBeam = 4.2; break;
    case 4072: iEBeam = 4.1; break;
    case 4482: iEBeam = 4.5; break;
    case 4726: iEBeam = 4.7; break;
    case 4756: iEBeam = 4.8; break;
    case 5052: iEBeam = 5.1; break;
    case 5166: iEBeam = 5.2; break;
    }

    iEPhot = stof(EPhot.c_str());
    
    if(Pol == "PARA"){iPol = 0;}
    else {iPol=1;}
    
    vector<float> Keys {iEBeam,iEPhot,iPol};
    
    keysPlane[Keys] = plane;

    chancount=0;                                             //starting array at 1 for consistency with E_ID
    
    while((fgets(fline,240,fpfile)) != NULL){
      if((fline[0] == '*')||(fline[0] == '#')) continue;     //skip comments    
      sscanf(fline,"%d",&eid);                               //first get the E_ID
      sscanf(fline,"%*d%lg%lg%lg%lg%lg%lg%lg%lg",
	     &polTable[plane][fcount][eid][ENERGY],
	     &polTable[plane][fcount][eid][ENH],
	     &polTable[plane][fcount][eid][ENHERR],
	     &polTable[plane][fcount][eid][ENHFIT],
	     &polTable[plane][fcount][eid][PFIT],
	     &polTable[plane][fcount][eid][PCOR],
	     &polTable[plane][fcount][eid][PCORERR],
	     &polTable[plane][fcount][eid][PSMOOTH]);
      chancount++;
    }
    
    fclose(fpfile); //close the file
    if(chancount!=384){
      cout << "Should be 384 lines in " << filename << " - only got " << chancount << endl;
      return -1;
    }

    polTableN[plane]++;
    
    fcount++;
  }
  
  
  fclose(fplist);
    
  return(0);
}


 
//double GetPol(int plane, double edge, int eid, int poltype = PSMOOTH, double lowThresh=0.2, double highThresh=0.3)
double GetPol(int plane, double edge, int eid, int poltype, double lowThresh, double highThresh)
{
  //get polarization based on eid and edge position
  
  int eIndex=0;
  double pol=-1.0;
  if((edge<polTable[plane][1][0][EDGE])||(edge>polTable[plane][polTableN[plane]-1][0][EDGE])) return -1.0;
  
  //find index
  for(eIndex=0;eIndex<polTableN[plane];eIndex++){
    if(polTable[plane][eIndex][0][EDGE]>=edge) break;
  }

  pol=polTable[plane][eIndex][eid][poltype];
  
  if((polTable[plane][eIndex][0][ENERGY]<edge)&&(pol<lowThresh)) pol = -1.0;
  if((polTable[plane][eIndex][0][ENERGY]>edge)&&(pol<highThresh)) pol = -1.0;
  
  return pol;
}

//double GetPol(int plane, double edge, double eg, int poltype = PSMOOTH, double lowThresh=0.2, double highThresh=0.3){
double GetPol(int plane, double edge, double eg, int poltype, double lowThresh, double highThresh){
  //get polarization based on ephoton energy and edge position
  
  int eIndex=0;
  double pol=-1.0;
  int eid=0;
  
  //Check edge in in range of tables
  if((edge<polTable[plane][1][0][EDGE])||(edge>polTable[plane][polTableN[plane]-1][0][EDGE])) return -1.0;
  
  //find index
  for(int eIndex=0;eIndex<polTableN[plane];eIndex++){
    if(polTable[plane][eIndex][0][EDGE]>=edge) break;
  }
  
  //find eid
  for(eid=1;eid<=384;eid++){
    if(polTable[plane][eIndex][eid][ENERGY]<=eg) break;
  }
  
  pol=polTable[plane][eIndex][eid][poltype];

  if((polTable[plane][eIndex][0][ENERGY]<edge)&&(pol<lowThresh)) pol = -1.0;
  if((polTable[plane][eIndex][0][ENERGY]>edge)&&(pol<highThresh)) pol = -1.0;
  
  return pol;
}


void GetPolAv(vector<float> Keys,vector<vector<int>> &ItP,vector<vector<double>> &AvP,double GetPol){


  if(Keys[1] == float(1.3) && Keys[2] == 0) {ItP[0][0]++; AvP[0][0]+=GetPol;}
  else if(Keys[0] == float(4.1) && Keys[1] == float(1.5) && Keys[2] == 0) {ItP[1][0]++; AvP[1][0]+=GetPol;}
  else if(Keys[0] == float(4.5) && Keys[1] == float(1.5) && Keys[2] == 0) {ItP[2][0]++; AvP[2][0]+=GetPol;}
  else if(Keys[0] == float(4.1) && Keys[1] == float(1.7) && Keys[2] == 0) {ItP[3][0]++; AvP[3][0]+=GetPol;}
  else if(Keys[0] == float(4.7) && Keys[1] == float(1.7) && Keys[2] == 0) {ItP[4][0]++; AvP[4][0]+=GetPol;}
  else if(Keys[0] == float(4.8) && Keys[1] == float(1.7) && Keys[2] == 0) {ItP[5][0]++; AvP[5][0]+=GetPol;}
  else if(Keys[1] == float(1.9) && Keys[2] == 0) {ItP[6][0]++; AvP[6][0]+=GetPol;}
  else if(Keys[0] == float(5.1) && Keys[1] == float(2.1) && Keys[2] == 0) {ItP[7][0]++; AvP[7][0]+=GetPol;}
  else if(Keys[0] == float(5.2) && Keys[1] == float(2.1) && Keys[2] == 0) {ItP[8][0]++; AvP[8][0]+=GetPol;}
  else if(Keys[0] == float(5.2) && Keys[1] == float(2.3) && Keys[2] == 0) {ItP[9][0]++; AvP[9][0]+=GetPol;}


  else if(Keys[1] == float(1.3) && Keys[2] == float(1)) {ItP[0][1]++; AvP[0][1]+=GetPol;}
  else if(Keys[0] == float(4.1) && Keys[1] == float(1.5) && Keys[2] == float(1)) {ItP[1][1]++; AvP[1][1]+=GetPol;}
  else if(Keys[0] == float(4.5) && Keys[1] == float(1.5) && Keys[2] == float(1)) {ItP[2][1]++; AvP[2][1]+=GetPol;}
  else if(Keys[0] == float(4.1) && Keys[1] == float(1.7) && Keys[2] == float(1)) {ItP[3][1]++; AvP[3][1]+=GetPol;}
  else if(Keys[0] == float(4.7) && Keys[1] == float(1.7) && Keys[2] == float(1)) {ItP[4][1]++; AvP[4][1]+=GetPol;}
  else if(Keys[0] == float(4.8) && Keys[1] == float(1.7) && Keys[2] == float(1)) {ItP[5][1]++; AvP[5][1]+=GetPol;}
  else if(Keys[1] == float(1.9) && Keys[2] == float(1)) {ItP[6][1]++; AvP[6][1]+=GetPol;}
  else if(Keys[0] == float(5.1) && Keys[1] == float(2.1) && Keys[2] == float(1)) {ItP[7][1]++; AvP[7][1]+=GetPol;}
  else if(Keys[0] == float(5.2) && Keys[1] == float(2.1) && Keys[2] == float(1)) {ItP[8][1]++; AvP[8][1]+=GetPol;}
  else if(Keys[0] == float(5.2) && Keys[1] == float(2.3) && Keys[2] == float(1)) {ItP[9][1]++; AvP[9][1]+=GetPol;}
  
}

void GetPolAvTable(vector<vector<int>> ItP,vector<vector<double>> AvP){
  cout << setprecision(2);
  cout << 1.3 << "\t" << 4199 << "\t" << AvP[0][0]/ItP[0][0] << "\t" << ItP[0][0] << "\t" << AvP[0][1]/ItP[0][1] << "\t" << ItP[0][1] << endl
       << 1.5 << "\t" << 4072 << "\t" << AvP[1][0]/ItP[1][0] << "\t" << ItP[1][0] << "\t" << AvP[1][1]/ItP[1][1] << "\t" << ItP[1][1] << endl
       << 1.5 << "\t" << 4482 << "\t" << AvP[2][0]/ItP[2][0] << "\t" << ItP[2][0] << "\t" << AvP[2][1]/ItP[2][1] << "\t" << ItP[2][1] << endl
       << 1.7 << "\t" << 4072 << "\t" << AvP[3][0]/ItP[3][0] << "\t" << ItP[3][0] << "\t" << AvP[3][1]/ItP[3][1] << "\t" << ItP[3][1] << endl
       << 1.7 << "\t" << 4726 << "\t" << AvP[4][0]/ItP[4][0] << "\t" << ItP[4][0] << "\t" << AvP[4][1]/ItP[4][1] << "\t" << ItP[4][1] << endl
       << 1.7 << "\t" << 4756 << "\t" << AvP[5][0]/ItP[5][0] << "\t" << ItP[5][0] << "\t" << AvP[5][1]/ItP[5][1] << "\t" << ItP[5][1] << endl
       << 1.9 << "\t" << 5052 << "\t" << AvP[6][0]/ItP[6][0] << "\t" << ItP[6][0] << "\t" << AvP[6][1]/ItP[6][1] << "\t" << ItP[6][1] << endl
       << 2.1 << "\t" << 5052 << "\t" << AvP[7][0]/ItP[7][0] << "\t" << ItP[7][0] << "\t" << AvP[7][1]/ItP[7][1] << "\t" << ItP[7][1] << endl
       << 2.1 << "\t" << 5166 << "\t" << AvP[8][0]/ItP[8][0] << "\t" << ItP[8][0] << "\t" << AvP[8][1]/ItP[8][1] << "\t" << ItP[8][1] << endl
       << 2.3 << "\t" << 5166 << "\t" << AvP[9][0]/ItP[9][0] << "\t" << ItP[9][0] << "\t" << AvP[9][1]/ItP[9][1] << "\t" << ItP[9][1] << endl;
        
}

void GetPolAvTableLatex(vector<vector<int>> ItP,vector<vector<double>> AvP,string Path,string Caption,string Label){
  cout << setprecision(2);
  fstream LTXT;
  LTXT.open(Path,ios::out);
  
  LTXT << "\\begin{table}[H]" << endl
       << "\t\\centering"       << endl
       << "\t\\begin{tabular}{|c|c|c|c|c|c|}" << endl
       << "\t\\hline" << endl
       << "\t\t$E_{\\gamma}$ Edge (GeV) & Energy Beam (GeV) & $\\bar{P}$ (Para) & N$^{\\circ}$ Events & $\\bar{P}$ (Perp) & N$^{\\circ}$ Events \\\\ \\hline\n"
       << "\t\t$" << 1.3 << "$\t&$" << 4.199 << "$\t&$" << AvP[0][0]/ItP[0][0] << "$\t&$" << ItP[0][0] << "$\t&$" << AvP[0][1]/ItP[0][1] << "$\t&$" << ItP[0][1] << "\\\\\n"
       << "\t\t$" << 1.5 << "$\t&$" << 4.072 << "$\t&$" << AvP[1][0]/ItP[1][0] << "$\t&$" << ItP[1][0] << "$\t&$" << AvP[1][1]/ItP[1][1] << "$\t&$" << ItP[1][1] << "\\\\\n"
       << "\t\t$" << 1.5 << "$\t&$" << 4.482 << "$\t&$" << AvP[2][0]/ItP[2][0] << "$\t&$" << ItP[2][0] << "$\t&$" << AvP[2][1]/ItP[2][1] << "$\t&$" << ItP[2][1] << "\\\\\n"
       << "\t\t$" << 1.7 << "$\t&$" << 4.072 << "$\t&$" << AvP[3][0]/ItP[3][0] << "$\t&$" << ItP[3][0] << "$\t&$" << AvP[3][1]/ItP[3][1] << "$\t&$" << ItP[3][1] << "\\\\\n"
       << "\t\t$" << 1.7 << "$\t&$" << 4.726 << "$\t&$" << AvP[4][0]/ItP[4][0] << "$\t&$" << ItP[4][0] << "$\t&$" << AvP[4][1]/ItP[4][1] << "$\t&$" << ItP[4][1] << "\\\\\n"
       << "\t\t$" << 1.7 << "$\t&$" << 4.756 << "$\t&$" << AvP[5][0]/ItP[5][0] << "$\t&$" << ItP[5][0] << "$\t&$" << AvP[5][1]/ItP[5][1] << "$\t&$" << ItP[5][1] << "\\\\\n"
       << "\t\t$" << 1.9 << "$\t&$" << 5.052 << "$\t&$" << AvP[6][0]/ItP[6][0] << "$\t&$" << ItP[6][0] << "$\t&$" << AvP[6][1]/ItP[6][1] << "$\t&$" << ItP[6][1] << "\\\\\n"
       << "\t\t$" << 2.1 << "$\t&$" << 5.052 << "$\t&$" << AvP[7][0]/ItP[7][0] << "$\t&$" << ItP[7][0] << "$\t&$" << AvP[7][1]/ItP[7][1] << "$\t&$" << ItP[7][1] << "\\\\\n"
       << "\t\t$" << 2.1 << "$\t&$" << 5.166 << "$\t&$" << AvP[8][0]/ItP[8][0] << "$\t&$" << ItP[8][0] << "$\t&$" << AvP[8][1]/ItP[8][1] << "$\t&$" << ItP[8][1] << "\\\\\n"
       << "\t\t$" << 2.3 << "$\t&$" << 5.166 << "$\t&$" << AvP[9][0]/ItP[9][0] << "$\t&$" << ItP[9][0] << "$\t&$" << AvP[9][1]/ItP[9][1] << "$\t&$" << ItP[9][1] << "\\\\\n"
       << "\t\t\\hline" << endl
       << "\t\\end{tabular}" << endl
       << "\t\\caption{"     << Caption << "}" << endl
       << "\t\\label{tab:"  << Label   << "}" << endl
       << "\\end{table}"     << endl;
  
  cout << "We Create: " << Path.substr(Path.find_last_of("/")+1,Path.size()-Path.find_last_of("/")) << endl;
  LTXT.close();
}

void GetEventPercent(vector<int> Events){
  cout << "\t"                      << "\t" << "Initial Events"               	  << "\t" << Events[0] << "\t" << Events[0]*100/Events[0] << endl
       << "Vertex Cut"              << "\t" << "Vertex \textit{z} of $K^{+}$"     << "\t" << Events[1] << "\t" << Events[1]*100/Events[0] << endl
       << ""                        << "\t" << "$\\Delta \\Beta$ to $\\pi^{-}$"   << "\t" << Events[2] << "\t" << Events[2]*100/Events[0] << endl
       << "Particle Identification" << "\t" << "$\\Delta \\Beta$ to $p$"          << "\t" << Events[3] << "\t" << Events[3]*100/Events[0] << endl
       << ""                        << "\t" << "$\\Delta \\Beta$ to $K^{+}$" 	  << "\t" << Events[4] << "\t" << Events[4]*100/Events[0] << endl
       << "$\\Delta t$ cut to most "
       << "likely photon"           << "\t" << "-1.0 < $\\Delta t_{K^{+}}$ < 1.0" << "\t" << Events[5]  << "\t" << Events[5]*100/Events[0]  << endl
       << "\t"			    << "\t" << "$\\phi_{p}$ cut" 		  << "\t" << Events[6]  << "\t" << Events[6]*100/Events[0]  << endl
       << "Fiduacial cuts"	    << "\t" << "$\\phi_{K^{+}}$ cut" 		  << "\t" << Events[7]  << "\t" << Events[7]*100/Events[0]  << endl
       << "\t"			    << "\t" << "$\\phi_{\\pi^{-}}$ cut"		  << "\t" << Events[8]  << "\t" << Events[8]*100/Events[0]  << endl
       << "\t" 			    << "\t" << "1" 				  << "\t" << Events[9]  << "\t" << Events[9]*100/Events[0]  << endl
       << "\t" 			    << "\t" << "2" 				  << "\t" << Events[10] << "\t" << Events[10]*100/Events[0] << endl
       << "Tager Ephp" 		    << "\t" << "3" 				  << "\t" << Events[11] << "\t" << Events[11]*100/Events[0] << endl
       << "\t" 			    << "\t" << "4" 				  << "\t" << Events[12] << "\t" << Events[12]*100/Events[0] << endl
       << "\t" 			    << "\t" << "5" 				  << "\t" << Events[13] << "\t" << Events[13]*100/Events[0] << endl
       << "Photo Pol"		    << "\t" << "Menor a 0.5"			  << "\t" << Events[14] << "\t" << Events[14]*100/Events[0] << endl
       << "Kaones no Piones"	    << "\t" << "" 				  << "\t" << Events[15] << "\t" << Events[15]*100/Events[0] << endl
       << "Background remove by "
       << "$K^{+}\\Lambda^{0}n $"   << "\t" << "Invariant Mass (8$\\sigma$)" 	  << "\t" << Events[16] << "\t" << Events[16]*100/Events[0] << endl
       << "Rescattering remove by "
       << "$K^{+}\\Lambda^{0}n $"   << "\t" << "Cut in the Missing Momentum" 	  << "\t" << Events[18] << "\t" << Events[18]*100/Events[0] << endl
       << "Background remove by "
       << "$K^{+}\\Lambda^{0}n $"   << "\t" << "it is in the code"	 	  << "\t" << Events[17] << "\t" << Events[17]*100/Events[0] << endl;
}

void GetEventPercentLatex(vector<int> Events,string Path,string Caption,string Label){
  fstream LTXT;
  LTXT.open(Path,ios::out);
  
  LTXT << "\\begin{table}[H]" 		     << endl
       << "\t\\centering"       	     << endl
       << "\t\\begin{tabular}{|c|c|c|c|}"    << endl
       << "\t\t\\hline"			     << endl
       << "\t\tSection" 		     << setw(43-7)  <<"\t&" << "Cut done"
       << setw(33-8)	<< "\t&" << "Events"   << "\t&" << "Events Percent $\\%$"  << "\\\\\\hline\n"
       << "\t\t"	                     << setw(43)    <<"\t&" << "Initial Events"
       << setw(33-14)   << "\t&" << Events[0]  << "\t&" << double(Events[0]*100)/double(Events[0]) << "\\\\\n"
       << "\t\tVertex Cut"                   << setw(43-10) <<"\t&" << "Vertex \\textit{z} of $K^{+}$"
       << setw(33-28)  	<< "\t&" << Events[1]  << "\t&" << double(Events[1]*100)/double(Events[0]) << "\\\\\n"
       << "\t\t"                             << setw(43)    <<"\t&" << "$\\Delta \\beta$ to $\\pi^{-}$"
       << setw(33-27)  	<< "\t&" << Events[2]  << "\t&" << double(Events[2]*100)/double(Events[0]) << "\\\\\n"
       << "\t\tParticle Identification"      << setw(43-23) <<"\t&" << "$\\Delta \\beta$ to $p$"
       << setw(33-21)  	<< "\t&" << Events[3]  << "\t&" << double(Events[3]*100)/double(Events[0]) << "\\\\\n"
       << "\t\t"                             << setw(43)    <<"\t&" << "$\\Delta \\beta$ to $K^{+}$"
       << setw(33-25)	<< "\t&" << Events[4]  << "\t&" << double(Events[4]*100)/double(Events[0]) << "\\\\\n"
       << "\t\t$\\Delta t$ cut to most "
       << "likely photon" 	             << setw(43-36) <<"\t&" << "$-1.0 < \\Delta t_{K^{+}} < 1.0$"
       << "\t&" << Events[5]  << "\t&" << double(Events[5]*100)/double(Events[0])   << "\\\\\n"
       << "\t\t"			     << setw(43)    <<"\t&" << "$\\phi_{p}$ cut"
       << setw(33-14) 		  	<< "\t&" << Events[6]  << "\t&" << double(Events[6]*100)/double(Events[0]) << "\\\\\n"
       << "\t\tFiduacial cuts"	    	     << setw(43-14) <<"\t&" << "$\\phi_{K^{+}}$ cut"
       << setw(33-18) 		<< "\t&" << Events[7]  << "\t&" << double(Events[7]*100)/double(Events[0]) << "\\\\\n"
       << "\t\t"			     << setw(43)    <<"\t&" << "$\\phi_{\\pi^{-}}$ cut"
       << setw(33-20)		<< "\t&" << Events[8]  << "\t&" << double(Events[8]*100)/double(Events[0]) << "\\\\\n"
       << "\t\t" 			     << setw(43)    <<"\t&" << "1"
       << setw(33-1)		<< "\t&" << Events[9]  << "\t&" << double(Events[9]*100)/double(Events[0]) << "\\\\\n"
       << "\t\t" 			     << setw(43)    <<"\t&" << "2"
       << setw(33-1)		<< "\t&" << Events[10] << "\t&" << double(Events[10]*100)/double(Events[0]) << "\\\\\n"
       << "\t\tTager Epho" 		     << setw(43-10) <<"\t&" << "3"
       << setw(33-1)   	        << "\t&" << Events[11] << "\t&" << double(Events[11]*100)/double(Events[0]) << "\\\\\n"
       << "\t\t" 			     << setw(43)    <<"\t&" << "4"
       << setw(33-1)		<< "\t&" << Events[12] << "\t&" << double(Events[12]*100)/double(Events[0]) << "\\\\\n"
       << "\t\t" 			     << setw(43)    <<"\t&" << "5"
       << setw(33-1)		<< "\t&" << Events[13] << "\t&" << double(Events[13]*100)/double(Events[0]) << "\\\\\n"
       << "\t\tPhoto Pol"		     << setw(43-9)  <<"\t&" << "Menor a 0.5"
       << setw(33-11)	        << "\t&" << Events[14] << "\t&" << double(Events[14]*100)/double(Events[0]) << "\\\\\n"
       << "\t\tKaones no Piones"	     << setw(43-16) <<"\t&" << ""
       << setw(33)	  	<< "\t&" << Events[15] << "\t&" << double(Events[15]*100)/double(Events[0]) << "\\\\\n"
       << "\t\tBackground remove by "
       << "$K^{+}\\Lambda^{0}$" 	     << setw(43-39) << "\t&" << "Invariant Mass (8$\\sigma$)"
       << setw(33-26) 	<< "\t&" << Events[16] << "\t&" << double(Events[16]*100)/double(Events[0]) << "\\\\\n"
       << "\t\tRescattering remove by "
       << "$K^{+}\\Lambda^{0}n$"   	     << "\t&"       << "Cut in the Missing Momentum"
       << setw(33-27)	<< "\t&" << Events[18] << "\t&" << double(Events[18]*100)/double(Events[0]) << "\\\\\n"
       << "\t\tBackground remove by "
       << "$K^{+}\\Lambda^{0}n$"   	     << "\t&" 	    << "it is in the code"
       << setw(33-17)	  	<< "\t&" << Events[17] << "\t&" << double(Events[17]*100)/double(Events[0]) << "\\\\\n"
       << "\t\t\\hline" 		     << endl
       << "\t\\end{tabular}" 		     << endl
       << "\t\\caption{"     		     << Caption << "}" << endl
       << "\t\\label{tab:"  		     << Label   << "}" << endl
       << "\\end{table}"     		     << endl;

  cout << "We Create: " << Path.substr(Path.find_last_of("/")+1,Path.size()-Path.find_last_of("/")) << endl;
  LTXT.close();
}

void GetEventBackgroundSigmaLaTeX(vector<vector<double>> Events,string Path,string Caption,string Label){
  fstream LTXT;
  LTXT.open(Path,ios::out);
  
  LTXT << "\\begin{table}[H]" 		     		<< endl
       << "\t\\centering"       	     		<< endl
       << "\t\\begin{tabular}{|c|c|c|}"		    	<< endl
       << "\t\t\\hline"			     		<< endl
       << "\t\t$\\Sigma^{-}$ Signal" 		     	<< setw(43-7)  		<< "\t&" << "Background"
       << "\t&"	<< "Events Background / Events $\\Sigma^{-}$ * 100 $\\%$"  	<< "\\\\\\hline\n"
       << "\t\t" << "$1.1-1.3$" << "\t&$" << Events[0][0] << "$\t&$" << Events[1][0] << "$\t&$" << double(Events[1][0])/double(Events[0][0]) << "$\\\\\n"
       << "\t\t" << "$1.3-1.5$" << "\t&$" << Events[0][1] << "$\t&$" << Events[1][1] << "$\t&$" << double(Events[1][1])/double(Events[0][1]) << "$\\\\\n"
       << "\t\t" << "$1.5-1.7$" << "\t&$" << Events[0][2] << "$\t&$" << Events[1][2] << "$\t&$" << double(Events[1][2])/double(Events[0][2]) << "$\\\\\n"
       << "\t\t" << "$1.7-1.9$" << "\t&$" << Events[0][3] << "$\t&$" << Events[1][3] << "$\t&$" << double(Events[1][3])/double(Events[0][3]) << "$\\\\\n"
       << "\t\t" << "$1.9-2.1$" << "\t&$" << Events[0][4] << "$\t&$" << Events[1][4] << "$\t&$" << double(Events[1][4])/double(Events[0][4]) << "$\\\\\n"
       << "\t\t" << "$2.1-2.3$" << "\t&$" << Events[0][5] << "$\t&$" << Events[1][5] << "$\t&$" << double(Events[1][5])/double(Events[0][5]) << "$\\\\\n"
       << "\t\t\\hline" 		     << endl
       << "\t\\end{tabular}" 		     << endl
       << "\t\\caption{"     		     << Caption << "}" << endl
       << "\t\\label{tab:"  		     << Label   << "}" << endl
       << "\\end{table}"     		     << endl;

  cout << "We Create: " << Path.substr(Path.find_last_of("/")+1,Path.size()-Path.find_last_of("/")) << endl;
  LTXT.close();
}
#endif
