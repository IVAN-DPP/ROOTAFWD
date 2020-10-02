#include "Histograms.h"
#include "CodeCuts.h"

#include <cstdlib>
#include <cstdio>

void Codigo(){

  Clase Cls;
  Cls.CodeCuts();

  remove("Codigo_C_ACLiC_dict_rdict.pcm");
  remove("Codigo_C.d");
  remove("Codigo_C.so");

}
