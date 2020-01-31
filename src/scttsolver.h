/*
 * scttsolver.h
 *
 *  Created on: 23-dec-2019
 *      Author: P. Sashittal
 */

#include "utils.h"
#include <fstream>

class ScttSolver
{
public:
		/// Default constructor
		ScttSolver();
		
		ScttSolver(int unhosts);
		
		/// read host file
		bool readHost(std::istream& in);
		
		/// read contact map
		bool readContactMap(std::istream& in);
		
		/// generate contact map
		void setFullContactMap();
		
		/// read transmission tree files
		bool readTransmissionTrees(std::istream& in);
		
		/// parentchild graph (edge label and edge cost on the contact map)
		void buildParentChildGraph();
		
		/// solve (find the minimum spanning tree on contact map)
		void solve();
		
		/// write solution file
		void writeConsensusTree(std::ostream& out) const;
		
		void writeConsensusTree(std::ostream& out, const std::string& msg) const;
		
		/// write DOT file
		void writeDot(std::ostream& out) const;
		
protected:
		
		// contact map variables
		Digraph _contactMap;
		StringNodeMap _nodeToId;
		StringToNodeMap _idToNode;
		StringVector _hostLabel;
		IntNodeMap _nodeToIndex;
		NodeVector _indexToNode;
		
		// host file variables
		int _nhosts;
		DoubleVector _enttime;
		DoubleVector _remtime;
		int _unsampledHosts;
		
		// transmission tree variables
		int _ktrees;
		IntVectorArcMap _treeEdgeWeights;
		
		// edge label and edge cost maps
		IntArcMap _edgeLabel;
		IntArcMap _edgeCost;
		
		// spanning tree solution
		BoolArcMap _mdst;
};
