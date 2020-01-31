/*
 * sankoffmain.cpp
 *
 *  Created on: 9-nov-2018
 *      Author: P. Sashittal
 */

#include "utils.h"
#include "naive.h"
#include "basetree.h"
#include <fstream>
#include <lemon/arg_parser.h>

int main(int argc, char** argv)
{
    
    typedef Digraph::NodeMap<int> IntNodeMap;
    typedef Digraph::NodeMap<double> DoubleNodeMap;
    typedef Digraph::ArcMap<int> IntArcMap;
    
    std::string host_filename, ptree_filename, ptree_sol_filename;
    int rootLabel = 0;
    int nrUnsampledHosts = 0;
    bool transTree = false;
    bool enumerate = false;
    bool consensus = false;
		bool parsimony = false;
		std::string contactmap_filename;
    int enumLimit = std::numeric_limits<int>::max();
    
    lemon::ArgParser ap(argc, argv);
    ap.refOption("r", "Root label (default: 0)", rootLabel);
		ap.refOption("p", "Parsimony flag (default: false)", parsimony);
		ap.refOption("m", "Contact map file: (default: empty)", contactmap_filename);
    ap.refOption("u", "Number of unsampled hosts (default: 0)", nrUnsampledHosts);
    ap.refOption("t", "Transmission tree instead of host file", transTree);
    ap.refOption("e", "Enumerate all the solutions (default: false)", enumerate);
    ap.refOption("c", "Find consensus Sankoff solution (deafault: false)", consensus);
    ap.refOption("l", "Enumeration solution number limit (default: intMax)", enumLimit);
    ap.other("<host> / <transmission_tree>");
    ap.other("<ptree>");
    ap.other("<output_ptree>");
    ap.parse();

    
    if (ap.files().size() != 3)
    {
        std::cerr << "Error1: expected <host>/<transmission_tree> <ptree> <output_file>" << std::endl;
        return 1;
    }

    host_filename = ap.files()[0];
    ptree_filename = ap.files()[1];
    ptree_sol_filename = ap.files()[2];
    
    std::ifstream host_file(host_filename.c_str());
    if (!host_file.good())
    {
        std::cerr << "Error2: failed opening '" << host_filename << "' for reading" << std::endl;
        return 1;
    }
    
    std::ifstream ptree_file(ptree_filename.c_str());
    if (!ptree_file.good())
    {
        std::cerr << "Error3: failed opening '" << ptree_filename << "' for reading" << std::endl;
        return 1;
    }
		
		std::ifstream contactmap_file;
		if (!contactmap_filename.empty())
		{
				contactmap_file.open(contactmap_filename.c_str());
				if (!contactmap_file.good())
				{
						std::cerr << "Error 4: failed openning '" << contactmap_filename << "' for reading" << std::endl;
						return 1;
				}
		}

    BaseTree B;
    try
    {
				if (!contactmap_filename.empty())
				{
						BaseTree::readInputWithHostsContactMap(host_file, ptree_file, contactmap_file, nrUnsampledHosts, true, B);
						
				}
				else
				{
						BaseTree::readInputWithHosts(host_file, ptree_file, nrUnsampledHosts, true, B);
				}
		}
    catch (std::runtime_error& e)
    {
        std::cerr << "Error4: " << e.what() << std::endl;
        return 1;
    }
    host_file.close();
    ptree_file.close();
		if (!contactmap_filename.empty())
		{
				contactmap_file.close();
		}
    
    std::ofstream ptree_sol(ptree_sol_filename.c_str());
    if (!ptree_sol.good())
    {
        std::cerr << "Error5: failed opening '" << ptree_sol_filename << "' for writing" << std::endl;
        return 1;
    }
    
    if (consensus && enumerate)
    {
        std::cerr << "Error6: enumerate and consensus does not work together. chose one" << std::endl;
        return 1;
    }
		
		if (consensus && !parsimony)
		{
				std::cerr << "Error7: can not find consensus solution for contact map constrained solutions" << std::endl;
				return 1;
		}
		
    Naive solver(B, rootLabel - 1);
		
		if (!parsimony)
		{
				solver.initContact();
				
				if (!enumerate)
				{
						solver.solveContact(ptree_sol, 1);
				}
				else
				{
						solver.solveContact(ptree_sol, enumLimit);
				}
		}
		else
		{
				solver.init();

				if (consensus)
				{
						solver.solve(ptree_sol);
				}
				else
				{
						solver.solve(ptree_sol, enumerate, enumLimit);
				}
		}
		
    ptree_sol.close();
    
    return 0;
}
