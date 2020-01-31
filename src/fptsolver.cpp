/*
 * fptsolver.cpp
 *
 *  Created on: 21-dec-2019
 *      Author: P. Sashittal
 */

#include "fptsolver.h"
#include <lemon/dfs.h>

FptSolver::FptSolver(const BaseTree& T,
										 const int rootLabel,
										 std::string psolFile)
		: _T(T)
		, _tTree(T.contactMap())
		, _rootLabel(rootLabel)
		, _psolFile(psolFile.c_str())
		, _solMap(_T.tree())
{
}

bool FptSolver::solveFull()
{
		const Digraph& G = tree();
		Node root = _T.root();
		int nrInfectedHosts = _T.getNrHost();
		
		// check if root label is out of bound
		// check if root host is possible
		if (_rootLabel >= nrInfectedHosts || _rootLabel < -1)
		{
				std::cerr << "root label is out of bounds" << std::endl;
				return false;
		}
		
		if (_rootLabel > -1)
		{
				if (entTime(_rootLabel) > time(root) || remTime(_rootLabel) < time(root))
				{
						std::cerr << "value of root label is " << _rootLabel + 1 << std::endl;
						std::cerr << "not a feasible root" << std::endl;
						return false;
				}
		}
		
		// check if leaf labels are feasible
		for (NodeIt vi(G); vi != lemon::INVALID; ++vi)
		{
				if (lemon::countOutArcs(G, vi) == 0)
				{
						int lhost = label(vi);
						if (entTime(lhost) > time(vi) || remTime(lhost) < time(vi))
						{
								std::cerr << "label of leaf " << _T.getName(vi) << " is wrong" << std::endl;
								return false;
						}
				}
		}

		/*
		for (int s = 0; s < nrInfectedHosts; ++s)
		{
				if (entTime(s) <= time(root) && remTime(s) >= time(root))
				{
						std::cout << "checking for root: " << s << std::endl;
						
						if (_tTree.setDFSTree(_T.contactMap(), _T.getContactNodeVector(), _T.getContactNode(s)))
						{
								if (enumeratedSolve())
								{
										return true;
								}
						}
						else
						{
								std::cout << "no DFS tree found" << std::endl;
						}
						
				}
		}
		*/
		
		for (int s = 0; s < nrInfectedHosts; ++s)
		{
				if (entTime(s) <= time(root) && remTime(s) >= time(root))
				{
						std::cout << "checking for root: " << s << std::endl;
						
						if (run(s))
						{
								std::cout << "solution found for root : " << s << std::endl;
								return true;
						}
						else
						{
								std::cout << "no solution for root : " << s << std::endl;
						}
				}
				else
				{
						std::cout << "no solution for root : " << s << std::endl;
				}
		}
		
		return false;
}


bool FptSolver::run(int s)
{
		const Digraph& contactMap = _T.contactMap();
		
		BoolNodeMap nodesG(contactMap, true);
		BoolNodeMap nodesT(contactMap, false);
		BoolNodeMap nodesL(contactMap, false);
		
		BoolArcMap arcsG(contactMap, true);
		BoolArcMap arcsT(contactMap, false);
		BoolArcMap arcsL(contactMap, false);
  
		SubDigraph G(contactMap, nodesG, arcsG);
		SubDigraph T(contactMap, nodesT, arcsT);
		SubDigraph L(contactMap, nodesL, arcsL);

		SubBfs bfsL(L);
		bfsL.init();
		
		Node r = _T.getContactNode(s);
		
		ArcList F;
		for (SubOutArcIt a(G, r); a != lemon::INVALID; ++a)
		{
				F.push_back(a);
		}
  
		T.enable(r);

		if (grow(r, G, T, L, bfsL, F))
		{
				return true;
		}
		else
		{
				return false;
		}
}


bool FptSolver::grow(Node root, SubDigraph& G, SubDigraph& T, SubDigraph& L, SubBfs& bfsL, ArcList& F)
{
		const Digraph& contactMap = _T.contactMap();
		int nNodes = lemon::countNodes(contactMap);
		
		// finds all spanning trees rooted at r containing T
		if (lemon::countArcs(T) ==  nNodes - 1)
		{
				// clear L
				for (NodeIt v(contactMap); v != lemon::INVALID; ++v)
				{
						L.disable(v);
				}
				for (ArcIt a(contactMap); a != lemon::INVALID; ++a)
				{
						L.disable(a);
				}
			
				// report T and copy to L
				IntPairSet selectedPairs;
				L.enable(root);
				for (SubArcIt a(T); a != lemon::INVALID; ++a)
				{
						selectedPairs.insert(std::make_pair(_T.getContactIndex(T.source(a)), _T.getContactIndex(T.target(a))));
						L.enable((Arc)a);
						L.enable((Node)T.target(a));
				}

				_tTree.changeTreeFromPairs(contactMap, _T.getContactIndexMap(), selectedPairs);
				
				if (solve())
				{
						return true;
				}
		}
		else
		{
				ArcList FF;
			
				bool done;
				do
				{
						//assert(!F.empty());
						if (F.empty())
						{
								return false;
						}
					
						Arc uv = F.back();
						F.pop_back();
						Node u = contactMap.source(uv);
						Node v = contactMap.target(uv);
					
						assert(T.status(u));
						assert(!T.status(v));
						assert(!T.status(uv));
					
						// add uv to T
						T.enable(uv);
						T.enable(v);
						//std::cout << _label[u] << " -> " << _label[v] << std::endl;
					
						ArcList newF = F;
					
						// push each arc vw where w not in V(T) onto F
						for (SubOutArcIt vw(G, v); vw != lemon::INVALID; ++vw)
						{
								Node w = G.target(vw);
								if (!T.status(w))
								{
										newF.push_back(vw);
								}
						}
					
						// remove each arc wv where w in T from F
						for (ArcListNonConstIt it = newF.begin(); it != newF.end();)
						{
								if (contactMap.target(*it) == v && T.status(contactMap.source(*it)))
								{
										it = newF.erase(it);
								}
								else
								{
										++it;
								}
						}
					
						if (grow(root, G, T, L, bfsL, newF))
						{
								return true;
						}
					
						G.disable(uv);
						T.disable(uv);
						T.disable(v);
					
						FF.push_back(uv);
					
						done = true;
						if (lemon::countNodes(L) != 0)
						{
								bfsL.run(v);
							
								for (SubInArcIt wv(G, v); wv != lemon::INVALID; ++wv)
								{
										Node w = G.source(wv);
										if (!bfsL.reached(w))
												done = false;
								}
						}
				} while (!done);
			
				for (ArcListRevIt it = FF.rbegin(); it != FF.rend(); ++it)
				{
						Arc a = *it;
						assert(!G.status(a));
					
						F.push_back(*it);
						G.enable(a);
				}
		}
		
		return false;
}

bool FptSolver::enumeratedSolve()
{
		if (solve())
		{
				return true;
		}
		
		const Digraph& contactMap = _T.contactMap();
		ArcSet nonBack, back;
		IntPairSet currentSelectedPairs;
		
		lemon::ArcLookUp<Digraph> transTreeLookUp(_tTree.tree());
		
		for (ArcIt eij(contactMap); eij != lemon::INVALID; ++eij)
		{
				int uIndex = _T.getContactIndex(contactMap.source(eij));
				int vIndex = _T.getContactIndex(contactMap.target(eij));
				
				Node u = _tTree.getContactNode(uIndex);
				Node v = _tTree.getContactNode(vIndex);
				
				if (transTreeLookUp(u, v) == lemon::INVALID)
				{
						if (isAncestor(_tTree.tree(), v, u))
						{
								back.insert(eij);
						}
						else
						{
								nonBack.insert(eij);
						}
				}
				else
				{
						currentSelectedPairs.insert(std::make_pair(uIndex, vIndex));
				}
		}
		
		_allSelectedPairs.insert(currentSelectedPairs);
		
		for (Arc nonbackArc : nonBack)
		{
				Node u = contactMap.source(nonbackArc);
				Node v = contactMap.target(nonbackArc);
				IntPairSet selectedPairs = currentSelectedPairs;

				for (auto it = selectedPairs.begin(); it != selectedPairs.end(); )
				{
						IntPair itPair = *it;
						if (itPair.second == _T.getContactIndex(v))
						{
								selectedPairs.erase(it);
								break;
						}
						else
						{
								++it;
						}
				}
				selectedPairs.insert(std::make_pair(_T.getContactIndex(u), _T.getContactIndex(v)));
				
				if (_allSelectedPairs.find(selectedPairs) == _allSelectedPairs.end())
				{
						_tTree.changeTreeFromPairs(contactMap, _T.getContactIndexMap(), selectedPairs);
						
						if (enumeratedSolve())
						{
								return true;
						}
				}
		}
		
		return false;
}

bool FptSolver::solve()
{
		const Digraph& G = tree();
		Node root = _T.root();
		const Digraph& C = _tTree.tree();
		
		// initialize all solMap to -1
		for (NodeIt vi(G); vi != lemon::INVALID; ++vi)
		{
				_solMap[vi] = -1;
		}
		
		// set all leaf labels
		for (NodeIt vi(G); vi != lemon::INVALID; ++vi)
		{
				if (OutArcIt(G, vi)  == lemon::INVALID)
				{
						_solMap[vi] = _T.getHostLabel(vi);
				}
		}
		
		// find the root label
		Node rootHost;
		
		for (NodeIt p(C); p != lemon::INVALID; ++p)
		{
				if (InArcIt(C, p) == lemon::INVALID)
				{
						rootHost = p;
				}
		}
		
		// label the root node with root host
		_solMap[root] = _tTree.getContactIndex(rootHost);
		//std::cout << "root host is labeled: " << _T.getContactIndex(rootHost) << std::endl;
		
		// solve for the vertex labeling
		if (!fixLCALabels(rootHost))
		{
				return false;
		}
		
		// fill the blanks
		for (NodeIt vi(G); vi != lemon::INVALID; ++vi)
		{
				if (OutArcIt(G, vi) == lemon::INVALID)
				{
						Node leafHost = _tTree.getContactNode(_solMap[vi]);
						Node leafHostParent;
						if (InArcIt(C,leafHost) != lemon::INVALID)
						{
								 leafHostParent = C.source(InArcIt(C, leafHost));
						}
						else
						{
								leafHostParent = leafHost;
						}
								
						
						Node climber = vi;
						
						while (InArcIt(G, climber) != lemon::INVALID)
						{
								climber = G.source(InArcIt(G,climber));
								if (_solMap[climber] == -1)
								{
										_solMap[climber] = _solMap[vi];
								}
								else if (_solMap[climber] == _tTree.getContactIndex(leafHost) || _solMap[climber] == _tTree.getContactIndex(leafHostParent))
								{
										break;
								}
								else
								{
										// no solution
										return false;
								}
										
						}
				}
		}
		
		// write the solution
		if (!_psolFile.good())
		{
				std::cerr << "error: failed to open output file" << std::endl;
				return false;
		}
		else
		{
				_T.writePtree(_psolFile, _solMap);
				return true;
		}
		
}

bool FptSolver::fixLCALabels(Node q)
{
		const Digraph& G = tree();
		const Digraph& C = _tTree.tree();

		for (NodeIt a(C); a != lemon::INVALID; ++a)
		{
				for (NodeIt b(C); b != lemon::INVALID; ++b)
				{
						if (_tTree.getContactIndex(a) < _tTree.getContactIndex(b))
						{
								NodeSet pairHosts {a, b};
								if ( q == getLCA(C, pairHosts))
								{
										NodeSet aLeafSet;
										NodeSet bLeafSet;
										
										for (NodeIt vi(G); vi != lemon::INVALID; ++vi)
										{
												if (OutArcIt(G, vi) == lemon::INVALID)
												{
														if (_T.getHostLabel(vi) == _tTree.getContactIndex(a))
														{
																aLeafSet.insert(vi);
														}
														else if (_T.getHostLabel(vi) == _tTree.getContactIndex(b))
														{
																bLeafSet.insert(vi);
														}
												}
										}
										
										for (Node aLeaf : aLeafSet)
										{
												NodeList aLeafPath = pathFromRoot(G, aLeaf);
												Node aAncestor;
												
												for (Node vine : aLeafPath)
												{
														if (_solMap[vine] == _tTree.getContactIndex(q))
														{
																aAncestor = vine;
																break;
														}
												}
												
												for (Node bLeaf : bLeafSet)
												{
														NodeList bLeafPath = pathFromRoot(G, bLeaf);
														Node bAncestor;
														
														for (Node vine : bLeafPath)
														{
																if (_solMap[vine] == _tTree.getContactIndex(q))
																{
																		bAncestor = vine;
																		break;
																}
														}
														
														if (aAncestor == bAncestor)
														{
																NodeSet leafNodes {aLeaf, bLeaf};
																Node anchor = getLCA(G, leafNodes);
																
																if (_T.getTime(anchor) > _T.getRemTime(_tTree.getContactIndex(q)) || _T.getTime(anchor) < _T.getEntTime(_tTree.getContactIndex(q)))
																{
																		//no solution
																		return false;
																}
																else
																{
																		NodeList anchorPath = pathFromRoot(G, anchor);
																		Node pathHost = _tTree.getContactNode(_solMap[root()]);
																		
																		for (Node vine : anchorPath)
																		{
																				if (_solMap[vine] == -1)
																				{
																						if (!isFeasibleLabel(vine, q))
																						{
																								return false;
																						}
																						else
																						{
																								_solMap[vine] = _T.getContactIndex(q);
																						}
																				}
																				else
																				{
																						Node currHost = _tTree.getContactNode(_solMap[vine]);
																						if (pathHost == currHost)
																						{
																								continue;
																						}
																						else if (isChild(C, currHost, pathHost))
																						{
																								pathHost = currHost;
																						}
																						else
																						{
																								// no solution
																								return false;
																						}
																				}
																		}
																}
														}
												}
										}
								}
						}
				}
		}
		
		for (OutArcIt qr(C, q); qr != lemon::INVALID; ++qr)
		{
				if (!seedChildTrees(C.target(qr)))
				{
						return false;
				}
		}

		for (OutArcIt qr(C, q); qr != lemon::INVALID; ++qr)
		{
				if (!fixLCALabels(C.target(qr)))
				{
						return false;
				}
		}
		
		return true;
}



bool FptSolver::seedChildTrees(Node r) const
{
		const Digraph& G = tree();
		const Digraph& C = _tTree.tree();

		Node rParent = C.source(InArcIt(C,r));
		
		NodeSet subTreeHosts;
		NodeListVector allPaths;
		double infectionTime = _T.getEntTime(_tTree.getContactIndex(r));
		
		for (NodeIt a(C); a != lemon::INVALID; ++a)
		{
				if (isAncestor(C, r, a))
				{
						subTreeHosts.insert(a);
				}
		}
		
		NodeSet anchorSet;
		
		for (NodeIt vi(G); vi != lemon::INVALID; ++vi)
		{
				if (isLeaf(G, vi))
				{
						Node leafHost = _tTree.getContactNode(_T.getHostLabel(vi));
						if (subTreeHosts.find(leafHost) != subTreeHosts.end())
						{
								Node climber = vi;
								Node anchor;
								while (InArcIt(G, climber) != lemon::INVALID)
								{
										anchor = climber;
										climber = G.source(InArcIt(G, climber));
										
										int currLabel = _solMap[climber];
										
										if (currLabel == -1)
										{
												continue;
										}
										else if (currLabel == _tTree.getContactIndex(rParent))
										{
												if (anchorSet.insert(anchor).second)
												{
														if (_T.getTime(climber) > infectionTime)
														{
																infectionTime = _T.getTime(climber);
														}
												}
												break;
										}
										else
										{
												//std::cout << "error for host " << _tTree.getContactIndex(r) << std::endl;
												return false;
										}
								}
						}
				}
		}
		
		if (infectionTime > _T.getRemTime(_tTree.getContactIndex(r)))
		{
				return false;
		}
		
		//std::cout << "Infection time for host " << _tTree.getContactIndex(r) << " is " << infectionTime << std::endl;
		
		for (Node anchor : anchorSet)
		{
				if (!isLeaf(G, anchor))
				{
						//auxSeedChildTrees(anchor, _tTree.getContactIndex(rParent), _tTree.getContactIndex(r), infectionTime);
						if (!auxSeedChildTrees(anchor, rParent, r, infectionTime))
						{
								return false;
						}
				}
		}
		
		return true;
}

bool FptSolver::auxSeedChildTrees(Node u, Node q, Node r, double infTime) const
{
		const Digraph& G = tree();
		int qIndex = _T.getContactIndex(q);
		int rIndex = _T.getContactIndex(r);
		
		if (_T.getTime(u) < infTime)
		{
				if (!isFeasibleLabel(u, q) || isLeaf(G, u))
				{
						return false;
				}
				else
				{
						_solMap[u] = qIndex;
				}
				
				
				for (OutArcIt uv(G, u); uv != lemon::INVALID; ++uv)
				{
						Node v = G.target(uv);
						
						if (!isLeaf(G, v))
						{
								if(!auxSeedChildTrees(v, q, r, infTime))
								{
										return false;
								}
						}
				}
		}
		else
		{
				if (!isFeasibleLabel(u, r))
				{
						return false;
				}
				else
				{
						if (!isLeaf(G, u))
						{
								_solMap[u] = rIndex;
						}
				}
		}
		
		return true;
}

/*
void FptSolver::auxSeedChildTrees(Node u, int q , int r, double infTime) const
{
		const Digraph& G = tree();
		
		if (_T.getTime(u) < infTime)
		{
				_solMap[u] = q;
				
				for (OutArcIt uv(G, u); uv != lemon::INVALID; ++uv)
				{
						Node v = G.target(uv);
						
						if (!isLeaf(G, v))
						{
								auxSeedChildTrees(v, q, r, infTime);
						}
				}
		}
		else
		{
				_solMap[u] = r;
		}
		
}
*/

bool FptSolver::isFeasibleLabel(const Node u, const Node q) const
{
		int qIndex = _T.getContactIndex(q);
		
		if (_T.getTime(u) < _T.getEntTime(qIndex) || _T.getTime(u) > _T.getRemTime(qIndex))
		{
				return false;
		}
		else
		{
				return true;
		}
}

Node FptSolver::getLCA(const Digraph& T, const NodeSet& nodes)
{
		if (nodes.size() == 1)
		{
				return *nodes.begin();
		}
		
		NodeListVector allPaths;
		NodeListItVector iteratorList;
		
		for (Node node : nodes)
		{
				allPaths.push_back(pathFromRoot(T, node));
				iteratorList.push_back(allPaths.back().begin());
		}
		
		Node lca;
		bool same = true;
		while (same)
		{
				lca = *iteratorList.front();
    
				for (int i = 0; i < iteratorList.size(); ++i)
				{
						if (++iteratorList[i] == allPaths[i].end())
						{
								same = false;
						}
				}
    
				if (same)
				{
						Node first = *iteratorList.front();
						for (const NodeListIt& it : iteratorList)
						{
								if (*it != first)
								{
										same = false;
								}
						}
				}
		}
		
		return lca;
}

NodeList FptSolver::pathFromRoot(const Digraph& T, Node u)
{
		NodeList path;
		
		while(InArcIt(T, u) != lemon::INVALID)
		{
				path.push_front(u);
				u = T.source(InArcIt(T, u));
		}
		
		path.push_front(u);
		
		return path;
}

bool FptSolver::isAncestor(const Digraph& T, Node u, Node v)
{
		while (v != lemon::INVALID)
		{
				if (u == v)
				{
						return true;
				}
				else
				{
						InArcIt a(T, v);
						if (a == lemon::INVALID)
						{
								v = lemon::INVALID;
						}
						else
						{
								v = T.source(a);
						}
				}
		}
		
		return u == v;
}

bool FptSolver::isLeaf(const Digraph& T, const Node u)
{
		if (OutArcIt(T, u) == lemon::INVALID)
		{
				return true;
		}
		else
		{
				return false;
		}
}

bool FptSolver::isChild(const Digraph& T, const Node v, const Node u)
{
		if (u == T.source(InArcIt(T, v)))
		{
				return true;
		}
		else
		{
				return false;
		}
}
