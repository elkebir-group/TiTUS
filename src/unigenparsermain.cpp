/*
 * samplesankoffmain.cpp
 *
 *  Created on: 9-nov-2018
 *      Author: P. Sashittal
 */

#include "utils.h"
#include "unigenparser.h"
#include "basetree.h"
#include <fstream>
#include <lemon/arg_parser.h>

int main(int argc, char** argv)
{
  
    typedef Digraph::NodeMap<int> IntNodeMap;
    typedef Digraph::NodeMap<double> DoubleNodeMap;
    typedef Digraph::ArcMap<int> IntArcMap;
		int muMax = std::numeric_limits<int>::max();
		
		std::string ptree_filename, host_filename, varlist_filename, unigen_filename, trans_prefix;
    std::string summary_filename;
		std::string contactmap_filename;
		std::string ttree_filename;
		
		lemon::ArgParser ap(argc, argv);
    ap.refOption("S", "summary file name", summary_filename);
		ap.refOption("muMax", "maximum value of mu : (default: max<int>)", muMax);
		ap.refOption("c", "contact map file : (default: empty)", contactmap_filename);
		ap.refOption("t", "transmission tree file: (default: empty)", ttree_filename);
		ap.other("<host>");
		ap.other("<ptree>");
		ap.other("<varlist>");
		ap.other("<unigen_file>");
		ap.other("<output_prefix>");
		ap.parse();

		if (ap.files().size() != 5)
		{
				std::cerr << "Error1: expected <host> <ptree> <varlist> <unigen_file> <output_prefix>" << std::endl;
				return 1;
		}
		
		if (!contactmap_filename.empty() && !ttree_filename.empty())
		{
				std::cerr << "cannot use both contact map and transmission tree input simulataneously" << std::endl;
		}
		
		host_filename = ap.files()[0];
		ptree_filename = ap.files()[1];
		varlist_filename = ap.files()[2];
		unigen_filename = ap.files()[3];
		trans_prefix = ap.files()[4];
		
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
		std::ifstream ttree_file;
		if (!contactmap_filename.empty())
		{
				contactmap_file.open(contactmap_filename.c_str());
				if (!contactmap_file.good())
				{
						std::cerr << "Error 4: failed openning '" << contactmap_filename << "' for reading" << std::endl;
						return 1;
				}
		}
		else if (!ttree_filename.empty())
		{
				ttree_file.open(ttree_filename.c_str());
				if (!ttree_file.good())
				{
						std::cerr << "Error 5: failed openning '" << contactmap_filename << "' for reading" << std::endl;
						return 1;
				}
		}
		
		BaseTree B;
		try
		{
				if (!contactmap_filename.empty())
				{
						BaseTree::readInputWithHostsContactMap(host_file, ptree_file, contactmap_file, 0, true, B);
				}
				else if (!ttree_filename.empty())
				{
						BaseTree::readInputWithHostsTransmissionTree(host_file, ptree_file, ttree_file, 0, true, B);
				}
				else
				{
						BaseTree::readInputWithHosts(host_file, ptree_file, 0, true, B);
				}
		}
		catch (std::runtime_error& e)
		{
				std::cerr << "Error6: " << e.what() << std::endl;
				return 1;
		}
		
		host_file.close();
		ptree_file.close();
		if (!contactmap_filename.empty())
		{
				contactmap_file.close();
		}
		if (!ttree_filename.empty())
		{
				ttree_file.close();
		}
		
		bool computeDist = false;
		if (!ttree_filename.empty())
		{
				computeDist = true;
		}

		UnigenParser parser(B, trans_prefix, muMax, computeDist);
		
		std::ifstream varlist_file(varlist_filename.c_str());
		if (!varlist_file.good())
		{
				std::cerr << "Error7: failed opening '" << varlist_filename << "' for reading" << std::endl;
				return 1;
		}
		
		parser.readVarFile(varlist_file);
		varlist_file.close();
		
		std::ifstream unigen_file(unigen_filename.c_str());
		if (!unigen_file.good())
		{
				std::cerr << "Error8: failed opening '" << unigen_filename << "' for reading" << std::endl;
				return 1;
		}
		
		parser.parseUnigen(unigen_file);
		unigen_file.close();
  
    std::ofstream outTT(trans_prefix);
    parser.writeTransTrees(outTT);
    outTT.close();
  
		if (!summary_filename.empty())
		{
				std::ofstream outS(summary_filename.c_str());
				if (computeDist)
				{
						parser.writeSummaryStatsWithDist(outS);
				}
				else
				{
						parser.writeSummaryStats(outS);	
				}
				outS.close();
		}
  
		return 0;
}
