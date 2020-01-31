/*
 * fptmain.cpp
 *
 *  Created on: 9-dec-2019
 *      Author: P. Sashittal
 */

#include "utils.h"
#include "fptsolver.h"
#include "basetree.h"
#include <fstream>
#include <lemon/arg_parser.h>

int main(int argc, char** argv)
{

    std::string outputptree_filename, ptree_filename, host_filename;
    int rootLabel = 0;
    int nrUnsampledHosts = 0;
		std::string contactmap_filename;

    lemon::ArgParser ap(argc, argv);
		ap.refOption("r", "Root label (default: 0)", rootLabel);
		ap.refOption("c", "contact map file: (default: empty)", contactmap_filename);
    ap.other("<host>");
    ap.other("<ptree>");
    ap.other("<output_ptree_file>");
    ap.parse();


    if (ap.files().size() != 3)
    {
        std::cerr << "Error1: expected <host> <ptree> <output_ptree_file>" << std::endl;
        return 1;
    }

    host_filename = ap.files()[0];
    ptree_filename = ap.files()[1];
    outputptree_filename = ap.files()[2];

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
				std::cerr << "Error5: " << e.what() << std::endl;
				return 1;
		}
		host_file.close();
		ptree_file.close();
		if (!contactmap_filename.empty())
		{
				contactmap_file.close();
		}
		
		FptSolver solver(B, rootLabel - 1, outputptree_filename.c_str());
		
		//solver.solve();
		if (!solver.solveFull())
		{
				std::cout << "no solution found" << std::endl;
		}
		else
		{
				std::cout << "solution successfully found" << std::endl;
		}
		
    host_file.close();
    ptree_file.close();
		contactmap_file.close();		
		
    return 0;
}
