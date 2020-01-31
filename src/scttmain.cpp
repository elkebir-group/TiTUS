/*
 * scttmain.cpp
 *
 *  Created on: 23-dec-2019
 *      Author: P. Sashittal
 */

#include "utils.h"
#include "scttsolver.h"
#include <fstream>
#include <lemon/arg_parser.h>

int main(int argc, char** argv)
{

    typedef Digraph::NodeMap<int> IntNodeMap;
    typedef Digraph::NodeMap<double> DoubleNodeMap;
    typedef Digraph::ArcMap<int> IntArcMap;

		std::string output_consensus_filename, host_filename, trans_filename_prefix;
		std::string contactmap_filename;
		
		lemon::ArgParser ap(argc, argv);
		ap.refOption("c", "contact map file: (default: empty)", contactmap_filename);
		ap.refOption("m", "Contact map file: (default: empty)", contactmap_filename);
		ap.other("<host>");
		ap.other("<transmission tree prefix>");
		ap.other("<output_ptree_file>");
		ap.parse();
		
		if (ap.files().size() != 3)
		{
				std::cerr << "Error1: expected <host>/<transmission_tree> <ptree> <output_file>" << std::endl;
				return 1;
		}
		
		host_filename = ap.files()[0];
		trans_filename_prefix = ap.files()[1];
		output_consensus_filename = ap.files()[2];

		ScttSolver solver(0);
		
		std::ifstream host_file(host_filename.c_str());
		if (!host_file.good())
		{
				std::cerr << "Error2: failed opening '" << host_filename << "' for reading" << std::endl;
				return 1;
		}
		
		solver.readHost(host_file);
		
		host_file.close();
		
		if (!contactmap_filename.empty())
		{
				std::ifstream contact_file(contactmap_filename.c_str());
				if (!contact_file.good())
				{
						std::cerr << "Error2: failed opening '" << contactmap_filename << "' for reading" << std::endl;
						return 1;
				}
				else
				{
						solver.readContactMap(contact_file);
						
						contact_file.close();
				}
		}
		else
		{
				solver.setFullContactMap();
		}

		{
        std::ifstream inTT(trans_filename_prefix.c_str());
        if (!inTT.good())
        {
            std::cerr << "Error4: could not open '" << trans_filename_prefix << "' for reading" << std::endl;
            return 1;
        }
        solver.readTransmissionTrees(inTT);
        inTT.close();
		}
		
		solver.buildParentChildGraph();
		solver.solve();
		
		std::ofstream output_file(output_consensus_filename.c_str());
		solver.writeConsensusTree(output_file);
		output_file.close();
		
		output_consensus_filename += ".dot";
		std::ofstream output_dot(output_consensus_filename.c_str());
		solver.writeDot(output_dot);
		output_dot.close();
		
    return 0;
}
