/*
 * output.cc: 
 *
 *  provides a function to output a vector to a file
 *
 */

#include <fstream>

#include "output.hh"

// TODO: add error handling
void outputVector(std::vector<double> v, std::string filename) {
  std::ofstream out;
  out.open(filename.c_str());

  for (double d : v ) {
	out << d << std::endl;
  }
  out.close();
}
