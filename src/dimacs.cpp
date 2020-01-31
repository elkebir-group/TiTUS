/*
 * dimacs.cpp
 *
 *  Created on: 9-nov-2018
 *      Author: P. Sashittal
 */

#include "dimacs.h"

Dimacs::Dimacs(const BaseTree& T,
               const int rootLabel,
               std::string psol,
               std::string varlist)
    : _T(T)
		, _rootLabel(rootLabel)
		, _bottleneck(false)
    , _Rv(_T.tree())
		, _Rst()
		, _solFile(psol.c_str())
    , _varFile(varlist.c_str())
{
}

bool Dimacs::solveDimacs()
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
		
    initVariables();
		
		for (NodeIt vi(G); vi != lemon::INVALID; ++vi)
		{
				int sum = 0;
				for (int s = 0; s < nrInfectedHosts; ++s)
				{
						sum += _Rv[vi][s];
				}
				
				if (sum == 0)
				{
						std::cerr << "no host labeling feasible for vertex " << name(vi) << " : no solution!" << std::endl;
						return false;
				}
		}
    
    writeVariableList();
		
    // write independent support
    writeIndependentSupport();
    std::cout << "independent support written" << std::endl;
		
    // vertex label time constraints
    writeVertexTimeConstraints();
    std::cout<< "vertex time constraints written" << std::endl;
  
		// vertex label contact map constraints
		writeVertexContactConstraints();
		std::cout << "vertex contact constraints written" << std::endl;

		// infection edge time constraint
		writeInfectionTimeConstraints();
		std::cout << "infection time constraint written" << std::endl;
		
		// infection window constraints
		writeInfectionWindowConstraints();
		std::cout << "infection window constraints written" << std::endl;
		
		if (_bottleneck)
		{
				writeStrongBottleneckCosntraints();
				std::cout << "strong bottleneck constraint written" << std::endl;
		}
		
		// direct infection constraints
		if (_rootLabel > -1)
		{
				writeRootedDirectInfectionConstraints(_rootLabel);
		}
		else
		{
				writeUnRootedDirectInfectionConstraint();
		}
		std::cout << "direct infection constraint written" << std::endl;
  
    _solFile << "p cnf " << _numVar << " " << _clauses.size() << std::endl;
  
    for (auto clause : _clauses)
    {
        for (auto var : clause)
        {
          _solFile << var << " ";
        }
        _solFile << "0" << std::endl;
    }
    _solFile.close();
    
    return true;
}

void Dimacs::initVariables()
{
    if (!_solFile.good())
    {
        std::cerr << "error: failed opening solution file" << std::endl;
    }
    if (!_varFile.good())
    {
        std::cerr << "error: failed opening solution file" << std::endl;
    }
    
    const Digraph& G = tree();
		const Digraph& C = contactMap();
    const int nrInfectedHosts = _T.getNrHost();
    
    // init vertex variables
    for (NodeIt vi(G); vi != lemon::INVALID; ++vi)
    {
        _Rv[vi] = IntVector(nrInfectedHosts, 0);
    }
    
    int rank = 1;
    
    for (NodeIt vi(G); vi != lemon::INVALID; ++vi)
    {
				if (vi != root() || _rootLabel == -1)
				{
						for (int s = 0; s < nrInfectedHosts; ++s)
						{
								if (time(vi) >= entTime(s) && time(vi) <= remTime(s))
								{
										_Rv[vi][s] = rank;
										++rank;
								}
						}
				}
				else if (vi == root())
				{
						_Rv[vi][_rootLabel] = rank;
						++rank;
				}
    }

		_contact = BoolMatrix(nrInfectedHosts);
		for (int s = 0; s < nrInfectedHosts; ++s)
		{
				_contact[s] = BoolVector(nrInfectedHosts, false);
		}
		
		for (ArcIt eij(C); eij != lemon::INVALID; ++eij)
		{
				int sourceIndex = contactIndex(C.source(eij));
				int targetIndex = contactIndex(C.target(eij));
				
				_contact[sourceIndex][targetIndex] = true;
		}
		
		_Rst = IntMatrix(nrInfectedHosts);
		for (int s = 0; s < nrInfectedHosts; ++s)
		{
				_Rst[s] = IntVector(nrInfectedHosts);
				for (int t = 0; t < nrInfectedHosts; ++t)
				{
						if (_contact[s][t])
						{
								_Rst[s][t] = rank;
								++rank;
						}
				}
		}
		
    _numVar = rank - 1;
}

void Dimacs::writeVariableList()
{
    const Digraph& G = tree();
    const int nrInfectedHosts = _T.getNrHost();
    
    int rank = 1;
    
    for (NodeIt vi(G); vi != lemon::INVALID; ++vi)
    {
        for (int s = 0; s < nrInfectedHosts; ++s)
        {
            if (_Rv[vi][s] > 0)
            {
                _varFile << _Rv[vi][s] << " " << _T.getName(vi) << " " << _T.getContactName(s) << " (v,s)"<< std::endl;
                ++rank;
            }
        }
    }
    
		for (int s = 0; s < nrInfectedHosts; ++s)
		{
				for (int t = 0; t < nrInfectedHosts; ++t)
				{
						if (_Rst[s][t] > 0)
						{
								_varFile << _Rst[s][t] << " " << _T.getContactName(s) << " " << _T.getContactName(t) << " (s,t)" << std::endl;
								++rank;
						}
				}
		}
    
    _varFile.close();
}

void Dimacs::writeIndependentSupport()
{
    const Digraph& G = tree();
    const int nrInfectedHosts = _T.getNrHost();
    
    int rank = 1;
    
    _solFile << "c ind ";
    
    for (NodeIt vi(G); vi != lemon::INVALID; ++vi)
    {
        for (int s = 0; s < nrInfectedHosts; ++s)
        {
            if (_Rv[vi][s] > 0)
            {
                _solFile << getVarIndex(vi, s) << " ";
                ++rank;
                
                if ((rank - 1)%10 == 0)
                {
                    _solFile << "0" << std::endl;
                    _solFile << "c ind ";
                }
            }
        }
    }
    
    _solFile << " 0" << std::endl;
}

void Dimacs::writeVertexTimeConstraints()
{
    const Digraph& G = tree();
    const int nrInfectedHosts = _T.getNrHost();
		
		for (NodeIt vi(G); vi != lemon::INVALID; ++vi)
		{
				if (lemon::countOutArcs(G, vi) == 0)
				{
						int lhost = label(vi);
						for (int s = 0; s < nrInfectedHosts; ++s)
						{
								if (lhost == s)
								{
                    _clauses.push_back(IntVector({getVarIndex(vi, s)}));
//                    _solFile << getVarIndex(vi, s) << " 0" << std::endl;
								}
								else if (_Rv[vi][s] > 0)
								{
                    _clauses.push_back(IntVector({-getVarIndex(vi, s)}));
//                    _solFile << -getVarIndex(vi, s) << " 0" << std::endl;
								}
						}
				}
				else
				{
						IntVector feasibleLabels;

						for (int s = 0; s < nrInfectedHosts; ++s)
						{
								if (_Rv[vi][s] > 0)
								{
										feasibleLabels.push_back(s);
								}
						}
				
						writeOneHotConstraints(vi, feasibleLabels);
				}
		}
}

void Dimacs::writeOneHotConstraints(Node u, IntVector feasibleLabels)
{
		int nfeasible = feasibleLabels.size();
		assert(nfeasible > 0);
		
		for (int i = 0; i < nfeasible; ++i)
		{
				for (int j = 0; j < i; ++j)
				{
            _clauses.push_back(IntVector({-getVarIndex(u, feasibleLabels[i]), -getVarIndex(u, feasibleLabels[j])}));
//            _solFile << -getVarIndex(u, feasibleLabels[i]) << " " << -getVarIndex(u, feasibleLabels[j]) << " 0" << std::endl;
				}
		}
		
  
    _clauses.push_back(IntVector({getVarIndex(u, feasibleLabels[0])}));
//    _solFile << getVarIndex(u, feasibleLabels[0]);
		for (int i = 1; i < nfeasible; ++i)
		{
        _clauses.back().push_back(getVarIndex(u, feasibleLabels[i]));
//        _solFile << " " << getVarIndex(u, feasibleLabels[i]);
		}
//    _solFile << " 0" << std::endl;
}

void Dimacs::writeVertexContactConstraints()
{
		const Digraph& G = tree();
		int nrInfectedHosts = _T.getNrHost();
		
		for (NodeIt vi(G); vi != lemon::INVALID; ++vi)
		{
				for (int s = 0; s < nrInfectedHosts; ++s)
				{
						if (_Rv[vi][s] > 0)
						{
								for (int t = 0; t < nrInfectedHosts; ++t)
								{
										if (s != t && !_contact[s][t])
										{
												for (OutArcIt eij(G,vi); eij != lemon::INVALID; ++eij)
												{
														Node vj = G.target(eij);
														
														if (_Rv[vj][t] > 0)
														{
                                _clauses.push_back(IntVector({-getVarIndex(vi, s), -getVarIndex(vj, t)}));
//                                _solFile << -getVarIndex(vi, s) << " " << -getVarIndex(vj, t) << " 0" << std::endl;
														}
												}
										}
										else if (s != t && _contact[s][t])
										{
												for (OutArcIt eij(G,vi); eij != lemon::INVALID; ++eij)
												{
														Node vj = G.target(eij);
														
														if (_Rv[vj][t] > 0)
														{
                                _clauses.push_back(IntVector({-getVarIndex(vi, s), -getVarIndex(vj, t), getVarIndex(s, t)}));
//                                _solFile << -getVarIndex(vi, s) << " " << -getVarIndex(vj, t) << " " << getVarIndex(s, t) << " 0" << std::endl;
														}
												}
										}
								}
						}
				}
		}
}

void Dimacs::writeInfectionTimeConstraints()
{
		const Digraph& G = tree();
		int nrInfectedHosts = _T.getNrHost();
		
		for (ArcIt eij(G); eij != lemon::INVALID; ++eij)
		{
				Node vi = G.source(eij);
				Node vj = G.target(eij);
				
				for (ArcIt ekl(G); ekl != lemon::INVALID; ++ekl)
				{
						Node vk = G.source(ekl);
						Node vl = G.target(ekl);
						
						if (arc_index(eij) < arc_index(ekl))
						{
								double ti = time(vi);
								double tj = time(vj);
								double tk = time(vk);
								double tl = time(vl);
								
								if ((ti > tl) || (tj < tk))
								{
										for (int s = 0; s < nrInfectedHosts; ++s)
										{
												for (int t = 0; t < nrInfectedHosts; ++t)
												{
														if (s != t && _contact[s][t])
														{
																if (_Rv[vi][s] > 0 && _Rv[vj][t] >0 && _Rv[vk][s] > 0 && _Rv[vl][t])
																{
                                    _clauses.push_back(IntVector({-getVarIndex(vi, s), -getVarIndex(vj, t), -getVarIndex(vk, s), -getVarIndex(vl, t)}));
//                                    _solFile << -getVarIndex(vi, s) << " " << -getVarIndex(vj, t) << " " << -getVarIndex(vk, s) << " " << -getVarIndex(vl, t) << " 0" << std::endl;
																}
														}
												}
										}
								}
						}
				}
		}
}

void Dimacs::writeStrongBottleneckCosntraints()
{
		const Digraph& G = tree();
		int nrInfectedHosts = _T.getNrHost();
		
		for (ArcIt eij(G); eij != lemon::INVALID; ++eij)
		{
				Node vi = G.source(eij);
				Node vj = G.target(eij);
				
				for (ArcIt ekl(G); ekl != lemon::INVALID; ++ekl)
				{
						Node vk = G.source(ekl);
						Node vl = G.target(ekl);
						
						if (arc_index(eij) < arc_index(ekl))
						{
								for (int t = 0; t < nrInfectedHosts; ++t)
								{
										if (t != _rootLabel)
										{
												for (int s = 0; s < nrInfectedHosts; ++s)
												{
														if (s != t && _contact[s][t])
														{
																if (_Rv[vi][s] > 0 && _Rv[vj][t] && _Rv[vk][s] > 0 && _Rv[vl][t])
																{
                                    _clauses.push_back(IntVector({-getVarIndex(vi, s), -getVarIndex(vj, t), -getVarIndex(vk, s), -getVarIndex(vl, t)}));
//                                    _solFile << -getVarIndex(vi, s) << " " << -getVarIndex(vj, t) << " " << -getVarIndex(vk, s) << " " << -getVarIndex(vl, t) << " 0" << std::endl;
																}
														}
												}
										}
								}
						}
				}
		}
}

void Dimacs::writeRootedDirectInfectionConstraints(int rootLabel)
{
		int nrInfectedHosts = _T.getNrHost();
		Node root = _T.root();
		
		// fix the root label
    _clauses.push_back(IntVector({getVarIndex(root, rootLabel)}));
//    _solFile << getVarIndex(root, rootLabel) << " 0" << std::endl;
  
		// no one can infect the root label
		for (int s = 0; s < nrInfectedHosts; ++s)
		{
				if (s != _rootLabel && _contact[s][_rootLabel])
				{
            _clauses.push_back(IntVector({-getVarIndex(s, _rootLabel)}));
//            _solFile << -getVarIndex(s, _rootLabel) << " 0" << std::endl;
				}
		}
		
		// rest of hosts can not be infected by more than one host
		for (int t = 0; t < nrInfectedHosts; ++t)
		{
				if (t != rootLabel)
				{

						IntVector feasibleSources;
						
						for (int s = 0; s < nrInfectedHosts; ++s)
						{
								if (s != t && _contact[s][t])
								{
										feasibleSources.push_back(s);
								}
						}
						
						int nfeasible = feasibleSources.size();
						for (int i = 0; i < nfeasible; ++i)
						{
								for (int j = 0; j < i; ++j)
								{
                    _clauses.push_back(IntVector({-getVarIndex(feasibleSources[i], t), -getVarIndex(feasibleSources[j], t)}));
//                    _solFile << -getVarIndex(feasibleSources[i], t) << " " << -getVarIndex(feasibleSources[j], t) << " 0" << std::endl;
								}
						}
				}
		}
}

void Dimacs::writeInfectionWindowConstraints()
{
		const Digraph& G = tree();
		int nrInfectedHosts = _T.getNrHost();
		
		for (int s = 0; s < nrInfectedHosts; ++s)
		{
				if (_T.infectionWindowExists(s))
				{
						std::cout << "infection window constraint for host " << _T.getContactName(s) << std::endl;
						
						double outTime = _T.getInfectionWindow(s);
						
						for ( ArcIt eij(G); eij != lemon::INVALID; ++eij)
						{
								Node vi = G.source(eij);
								Node vj = G.target(eij);
								
								double ti = time(vi);
								double tj = time(vj);
								
								if (ti > outTime && tj <= _T.getRemTime(s))
								{
										/*
										if( _Rv[vi][s] == 0 || _Rv[vj][s] == 0)
										{
												std::cout << _T.getName(vi) << "(" << ti << ")" << " -> " << _T.getName(vj) << "(" << tj << ")" << " for host " << _T.getContactName(s) << " with values " << _Rv[vi][s] << " and " << _Rv[vj][s] << " where outtime is " << outTime << " and removal time is " << remTime(s) << std::endl;
										}
										*/
										
                    _clauses.push_back(IntVector({-getVarIndex(vj, s), getVarIndex(vi, s)}));
//                    _solFile << -getVarIndex(vj, s) << " " << getVarIndex(vi, s) << " 0" << std::endl;
								}
						}
				}
		}
}

void Dimacs::writeOneHotConstraints(int t, IntVector feasibleSources)
{
		int nfeasible = feasibleSources.size();
		assert(nfeasible > 0);
		
		for (int i = 0; i < nfeasible; ++i)
		{
				for (int j = 0; j < i; ++j)
				{
            _clauses.push_back(IntVector({-getVarIndex(feasibleSources[i], t), -getVarIndex(feasibleSources[j], t)}));
//            _solFile << -getVarIndex(feasibleSources[i], t) << " " << -getVarIndex(feasibleSources[j], t) << " 0" << std::endl;
				}
		}
		
		_clauses.push_back(IntVector({feasibleSources[0], t}));
//    _solFile << getVarIndex(feasibleSources[0], t);
		for (int i = 1; i < nfeasible; ++i)
		{
        _clauses.back().push_back(getVarIndex(feasibleSources[i], t));
//        _solFile << " " << getVarIndex(feasibleSources[i], t);
		}
//    _solFile << " 0" << std::endl;
}


void Dimacs::writeOneHotConstraints(Node u, BoolVector feasbileLabels)
{
    const int nrInfectedHosts = _T.getNrHost();
    assert(feasbileLabels.size() == nrInfectedHosts);
    
    IntVector feasiblePaths;
    
    for (int s = 0; s < nrInfectedHosts; ++s)
    {
        if (feasbileLabels[s])
        {
            feasiblePaths.push_back(s);
        }
    }
    
    int nfeasible = feasiblePaths.size();
    
    for (int i = 0; i < nfeasible; ++i)
    {
        for (int j = 0; j < i; ++j)
        {
            _clauses.push_back(IntVector({-getVarIndex(u, feasiblePaths[i]), -getVarIndex(u, feasiblePaths[j])}));
//            _solFile << -getVarIndex(u, feasiblePaths[i]) << " " << -getVarIndex(u, feasiblePaths[j]) << " 0" << std::endl;
        }
    }
    
    _clauses.push_back(IntVector({getVarIndex(u, feasiblePaths[0])}));
//    _solFile << getVarIndex(u, feasiblePaths[0]);
    for (int i = 1; i < nfeasible; ++i)
    {
        _clauses.back().push_back(getVarIndex(u, feasiblePaths[i]));
//        _solFile << " " << getVarIndex(u, feasiblePaths[i]);
    }
//    _solFile << " 0" << std::endl;
}

void Dimacs::writeUnRootedDirectInfectionConstraint()
{
		int nrInfectedHosts = _T.getNrHost();
		Node root = _T.root();
		
		// no one can get infected by more than one host
		for (int t = 0; t < nrInfectedHosts; ++t)
		{
				IntVector feasibleSources;
						
				for (int s = 0; s < nrInfectedHosts; ++s)
				{
						if (s != t && _contact[s][t])
						{
								feasibleSources.push_back(s);
						}
				}

				int nfeasible = feasibleSources.size();

				for (int i = 0; i < nfeasible; ++i)
				{
						for (int j = 0; j < i; ++j)
						{
                _clauses.push_back(IntVector({-getVarIndex(feasibleSources[i], t), -getVarIndex(feasibleSources[j], t)}));
//                _solFile << -getVarIndex(feasibleSources[i], t) << " " << -getVarIndex(feasibleSources[j], t) << " 0" << std::endl;
						}
				}
				
				// root node can not get infected by anyone
				for (int i = 0; i < nfeasible; ++i)
				{
						if (_Rv[root][t] > 0 && _contact[feasibleSources[i]][t])
            {
                _clauses.push_back(IntVector({-getVarIndex(root, t), -getVarIndex(feasibleSources[i], t)}));
//                _solFile << -getVarIndex(root, t) << " " << -getVarIndex(feasibleSources[i], t) << " 0" << std::endl;
						}
				}
		}
}
