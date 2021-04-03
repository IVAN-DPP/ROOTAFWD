#include "Histograms.h"
#include "CodeCuts.h"
#include "include/Libraries.h"

void Codigo(){

  Codecuts Cls;
  // Only active the method if CodeCuts SIGMA particle was changed
  Cls.CodeCuts();
  // Cls.CodeCutsCosBin();
  // Cls.CodeCutsAsym();
  remove("Codigo_C_ACLiC_dict_rdict.pcm");
  remove("Codigo_C.d");
  remove("Codigo_C.so");
  remove("Codigo_C_ACLiC_dict.cxx_tmp_28880");
}
