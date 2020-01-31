/*
 * gammamain.cpp
 *
 *  Created on: 27-jun-2019
 *      Author: M. El-Kebir
 */

#include "utils.h"
#include "basetree.h"
#include <fstream>
#include <lemon/arg_parser.h>

int main(int argc, char** argv)
{
  typedef Digraph::NodeMap<int> IntNodeMap;
  typedef Digraph::NodeMap<double> DoubleNodeMap;
  typedef Digraph::ArcMap<int> IntArcMap;
  
  std::string host_filename, ptree_sol_filename;
  bool nonbinaryTree = false;
  int nrUnsampledHosts = 0;
  
  lemon::ArgParser ap(argc, argv);
  ap.refOption("b", "is the tree non binary (default: false)", nonbinaryTree);
  ap.refOption("u", "Number of unsampled hosts (default: 0)", nrUnsampledHosts);
  ap.other("<host>");
  ap.other("<ptree_sol>");
  ap.parse();
  
  if (ap.files().size() != 2)
  {
    std::cerr << "Error1: expected <host> <ptree_sol>" << std::endl;
    return 1;
  }
  
  host_filename = ap.files()[0];
  ptree_sol_filename = ap.files()[1];
  
  std::ifstream host_file(host_filename.c_str());
  if (!host_file.good())
  {
    std::cerr << "Error2: failed opening '" << host_filename << "' for reading" << std::endl;
    return 1;
  }
  
  std::ifstream ptree_sol_file(ptree_sol_filename.c_str());
  if (!ptree_sol_file.good())
  {
    std::cerr << "Error3: failed opening '" << ptree_sol_filename << "' for reading" << std::endl;
    return 1;
  }
  
  BaseTree B;
  try
  {
    BaseTree::readInputWithHosts(host_file, ptree_sol_file, nrUnsampledHosts, !nonbinaryTree, B);
  }
  catch (std::runtime_error& e)
  {
    std::cerr << "Error4: " << e.what() << std::endl;
    return 1;
  }
  host_file.close();
  ptree_sol_file.close();
  
//  Sankoff sankoff(B, B.getHostLabel(B.root()));
//  sankoff.setSolFromInput();
  
  int mu = B.getMu(B.getHostLabeling());
  std::cout << "mu = " << mu << std::endl;
  
  ArcMatrix N = B.getN(B.getHostLabeling());
  int gamma = N.size();
  std::cout << "gamma = " << gamma << std::endl;
  
  B.writeDOT(N, B.getHostLabeling(), std::cerr);
  
  return 0;
}
