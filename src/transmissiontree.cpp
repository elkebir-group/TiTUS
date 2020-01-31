/*
 * transmissiontree.cpp
 *
 *  Created on: 23-dec-2019
 *      Author: P. Sashittal
 */

#include "transmissiontree.h"
#include <lemon/dfs.h>

TransmissionTree::TransmissionTree()
		: _tTree()
		, _nodeToIndex(_tTree)
		, _indexToNode()
{
}

TransmissionTree::TransmissionTree(const Digraph& G)
		: _tTree()
		, _nhosts()
		, _nodeToIndex(_tTree)
		, _indexToNode()
{
		_nhosts = lemon::countNodes(G);
		
		for (int i = 0; i < _nhosts; ++i)
		{
				Node u = _tTree.addNode();
				_nodeToIndex[u] = i;
				_indexToNode.push_back(u);
		}
}

/*
TransmissionTree::TransmissionTree(const TransmissionTree& other)
		: _tTree()
		, _nodeToIndex(_tTree)
		, _indexToNode()
{
		const Digraph& otherTree = other.tree();

}
*/

TransmissionTree::TransmissionTree(const Digraph& G,
																   const NodeVector indexToNodeG,
																	 const Node anchorNode)
		: _tTree()
		, _nodeToIndex(_tTree)
		, _indexToNode()
{
		// construct a DFS tree
		lemon::Dfs<Digraph> dfs(G);
		
		dfs.run(anchorNode);
		
		// make the DFS transmission tree
		// add nodes
		int nhosts = indexToNodeG.size();
		for (int i = 0; i < nhosts; ++i)
		{
				Node u = _tTree.addNode();
				_nodeToIndex[u] = i;
				_indexToNode.push_back(u);
		}
		
		// add arcs
		for (NodeIt v(G); v != lemon::INVALID; ++v)
		{
				Node u = dfs.predNode(v);
				if (u != lemon::INVALID)
				{
						_tTree.addArc(u, v);
				}
		}
		
		/*
		// write the transmission tree
		std::cout << "transmission tree is " << std::endl;
		for (ArcIt eij(_tTree); eij != lemon::INVALID; ++eij)
		{
				std::cout << _nodeToIndex[_tTree.source(eij)] << " -> " << _nodeToIndex[_tTree.target(eij)] << std::endl;
		}
		*/
}

bool TransmissionTree::setDFSTree(const Digraph& G, const NodeVector indexToNodeG, const Node anchorNode)
{
		// construct a DFS tree
		lemon::Dfs<Digraph> dfs(G);
		
		dfs.run(anchorNode);
		
		// make the DFS transmission tree
		// add nodes
		int nhosts = indexToNodeG.size();
		for (int i = 0; i < nhosts; ++i)
		{
				if (dfs.reached(indexToNodeG[i]))
				{
						Node u = _tTree.addNode();
						_nodeToIndex[u] = i;
						_indexToNode.push_back(u);
				}
				else
				{
						return false;
				}
		}
		
		_nhosts = nhosts;
		
		// add arcs
		for (NodeIt v(G); v != lemon::INVALID; ++v)
		{
				Node u = dfs.predNode(v);
				if (u != lemon::INVALID)
				{
						_tTree.addArc(u, v);
				}
		}
		
		
		writeTransmissionTree();
		/*
		// write the transmission tree
		std::cout << "transmission tree is " << std::endl;
		for (ArcIt eij(_tTree); eij != lemon::INVALID; ++eij)
		{
				std::cout << _nodeToIndex[_tTree.source(eij)] << " -> " << _nodeToIndex[_tTree.target(eij)] << std::endl;
		}
		*/
		
		return true;
}

bool TransmissionTree::setTreeFromMap(const Digraph& G, const IntNodeMap& nodeToIndexG, const BoolArcMap& selectedArcs)
{
		int nhosts = lemon::countNodes(G);
		for (int i = 0; i < nhosts; ++i)
		{
				Node u = _tTree.addNode();
				_nodeToIndex[u] = i;
				_indexToNode.push_back(u);
		}
		
		for (ArcIt eij(G); eij != lemon::INVALID; ++eij)
		{
				if (selectedArcs[eij])
				{
						Node u = G.source(eij);
						Node v = G.target(eij);
						
						int uIndex = nodeToIndexG[u];
						int vIndex = nodeToIndexG[v];
						
						_tTree.addArc(_indexToNode[uIndex], _indexToNode[vIndex]);
				}
		}
		
		//writeTransmissionTree();
		
		return true;
}

void TransmissionTree::changeTreeFromMap(const Digraph& G, const IntNodeMap& nodeToIndexG, const BoolArcMap& selectedArcs)
{
		for (ArcIt eij(_tTree); eij != lemon::INVALID; ++eij)
		{
				_tTree.erase(eij);
		}
		
		for (ArcIt eij(G); eij != lemon::INVALID; ++eij)
		{
				if (selectedArcs[eij])
				{
						Node u = G.source(eij);
						Node v = G.target(eij);
						
						int uIndex = nodeToIndexG[u];
						int vIndex = nodeToIndexG[v];
						
						_tTree.addArc(_indexToNode[uIndex], _indexToNode[vIndex]);
				}
		}
		
		//writeTransmissionTree();
}

void TransmissionTree::changeTreeFromPairs(const Digraph& G, const IntNodeMap& nodeToIndexG, const IntPairSet& selectedPairs)
{
		/*
		std::cout << "index pair set is: \n";
		for (IntPair nodeIndexPair : selectedPairs)
		{
				std::cout << nodeIndexPair.first << ", " << nodeIndexPair.second << std::endl;
		}
		*/
		
		lemon::DynArcLookUp<Digraph> transTreeLookUp(_tTree);
		
		for (ArcIt eij(G); eij != lemon::INVALID; ++eij)
		{
				Node u = G.source(eij);
				Node v = G.target(eij);
				
				int uIndex = nodeToIndexG[u];
				int vIndex = nodeToIndexG[v];
				
				Arc treeArc = transTreeLookUp(_indexToNode[uIndex], _indexToNode[vIndex]);
				
				if (selectedPairs.find(std::make_pair(uIndex, vIndex)) != selectedPairs.end())
				{
						if (treeArc == lemon::INVALID)
						{
								_tTree.addArc(_indexToNode[uIndex], _indexToNode[vIndex]);
						}
				}
				else
				{
						if (treeArc != lemon::INVALID)
						{
								_tTree.erase(treeArc);
						}
				}
		}
		
		//writeTransmissionTree();
}

void TransmissionTree::writeTransmissionTree()
{
		// write the transmission tree
		std::cout << "transmission tree is " << std::endl;
		for (ArcIt eij(_tTree); eij != lemon::INVALID; ++eij)
		{
				std::cout << _nodeToIndex[_tTree.source(eij)] << " -> " << _nodeToIndex[_tTree.target(eij)] << std::endl;
		}
}
