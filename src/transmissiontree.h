/*
 * transmissiontree.h
 *
 *  Created on: 23-dec-2019
 *      Author: P. Sashittal
 */

#include "utils.h"
#include <fstream>

class TransmissionTree
{
public:
		/// Default constructor
		TransmissionTree();
		
		TransmissionTree(const Digraph& G);
		
		/// copy constructor
		TransmissionTree(const TransmissionTree& other);
		
		void writeTransmissionTree();
		
		/// DFS tree constructor
		TransmissionTree(const Digraph& G, const NodeVector indexToNodeG, const Node anchorNode);
		
		// set DFS tree
		bool setDFSTree(const Digraph& G, const NodeVector indexToNodeG, const Node anchorNode);
		
		bool setTreeFromMap(const Digraph& G, const IntNodeMap& nodeToIndexG, const BoolArcMap& selectedArcs);
		
		void changeTreeFromMap(const Digraph& G, const IntNodeMap& nodeToIndexG, const BoolArcMap& selectedArcs);
		
		void changeTreeFromPairs(const Digraph& G, const IntNodeMap& nodeToIndexG, const IntPairSet& selectedPairs);
		
		void changeTreeFromArcList(const Digraph& G, const IntNodeMap& nodeToIndexG, const ArcList& arcsL);
		
		/// read transmission trees file
		
		/// read host file
		
		/// get DFS tree
		//Digraph getDFSTree(const Digraph& G);
		
		// enumerate the spanning trees
		void enumerateTrees(const Digraph& G);
		
		const Digraph& tree() const
		{
				return _tTree;
		}
		
		const Node getContactNode(int s) const
		{
				return _indexToNode[s];
		}
		
		const int getContactIndex(Node u) const
		{
				return _nodeToIndex[u];
		}
				
protected:
		
		Digraph _tTree;
		int _nhosts;
		IntNodeMap _nodeToIndex;
		NodeVector _indexToNode;
};
