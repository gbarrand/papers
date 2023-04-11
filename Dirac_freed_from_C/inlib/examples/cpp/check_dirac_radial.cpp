// Copyright (C) 2010, Guy Barrand. All rights reserved.
// See the file inlib.license for terms.

#include "../../tests/test/check_dirac_radial"

#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <iostream>

int main(int argc,char** argv) {
#ifdef INLIB_MEM
  inlib::mem::set_check_by_class(true);{
#endif //INLIB_MEM

  bool verbose = (argc==2 && !::strcmp(argv[1],"-verbose"))?true:false;

  if(!check_dirac_radial_x2(std::cout,verbose)) return EXIT_FAILURE;

  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  std::cout << "exit ..." << std::endl;

#ifdef INLIB_MEM
  }inlib::mem::balance(std::cout);
#endif //INLIB_MEM
  
  return EXIT_SUCCESS;
}
