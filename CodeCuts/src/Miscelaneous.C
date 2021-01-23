//Author: Edwin Munevar

/* Miscelaneous.C.  */
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

  
#endif
