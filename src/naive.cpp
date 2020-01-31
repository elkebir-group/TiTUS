/*
 * naive.cpp
 *
 *  Created on: 9-nov-2018
 *      Author: P. Sashittal
 */

#include "naive.h"

Naive::Naive(const BaseTree& T, const int rootLabel = -1)
    :  _T(T)
    , _rootLabel(rootLabel)
    , _solMap(_T.tree())
    , _M(_T.tree())
    , _N(_T.tree())
    , _W(_T.tree())
    , _G(_T.tree())
		//, _contactLookUp(_T.contactMap())
{
    if (_rootLabel == -1)
    {
        std::cout << "no root label prescribed" << std::endl;
    }
    else
    {
        std::cout << "value of the root label is " << _rootLabel + 1 << std::endl;
    }
    std::cout << "Solving Sankoff Labeling" << std::endl;
}

void Naive::init()
{
    const Digraph& G = tree();
    const int nrInfectedHosts = _T.getNrHost();
    
    for (NodeIt vi(G); vi != lemon::INVALID; ++vi)
    {
        _M[vi] = UInt64Vector(nrInfectedHosts);
        _N[vi] = UInt64Vector(nrInfectedHosts);
        _W[vi] = UInt64Vector(nrInfectedHosts);
        _G[vi] = UInt64Vector(nrInfectedHosts);
    }
}

void Naive::initContact()
{
		const Digraph& G = tree();
		const int nrInfectedHosts = _T.getNrHost();
		
		for (NodeIt vi(G); vi != lemon::INVALID; ++vi)
		{
				_N[vi] = UInt64Vector(nrInfectedHosts);
		}
}

bool Naive::solveSankoff()
{
    Node root = _T.root();
    run(root);
    
    if (_rootLabel == -1)
    {
        _solMap[root] = std::min_element(_M[root].begin(), _M[root].end()) - _M[root].begin();
    }
    else
    {
        _solMap[root] = _rootLabel;
    }
    
    if (_M[root][_solMap[root]] == std::numeric_limits<uint64_t>::max())
    {
        std::cout << "no solution exists!" << std::endl;
        return false;
    }
    else
    {
        std::cout << "Root Label is: " << _solMap[root] + 1 << std::endl;
        std::cout << "Infection number solution is: " << _M[root][_solMap[root]] << std::endl;
        computeSol(root, _solMap[root]);
        
        int gamma = computeCoInf();
        std::cout << "Coinfection number solution is: " << gamma << std::endl;
        
        return true;
    }
}

uint64_t Naive::countSankoff()
{
    Node root = _T.root();
    int nrInfectedHosts = _T.getNrHost();
    runCount(root);
    
    uint64_t numSols = 0;
    
    if (_rootLabel == -1)
    {
        uint64_t min_cost = *std::min_element(_M[root].begin(), _M[root].end());
        
        if (min_cost == std::numeric_limits<uint64_t>::max())
        {
            std::cout << "no solutions exist!" << std::endl;
            return 0;
        }
        
        std::cout << "feasible roots are: ";
        for (int s = 0; s < nrInfectedHosts; ++s)
        {
            if (_M[root][s] == min_cost)
            {
                std::cout << s + 1 << ", ";
                numSols += _N[root][s];
            }
            else
            {
                _N[root][s] = 0;
            }
        }
        std::cout << std::endl;
        
        std::cout << "Infection number solution is: " << min_cost << std::endl;
    }
    else
    {
        if (_M[root][_rootLabel] == std::numeric_limits<uint64_t>::max())
        {
            std::cout << "no solutions exist!" << std::endl;
            return 0;
        }
        
        numSols = _N[root][_rootLabel];
        
        std::cout << "Infection number solution is: " << _M[root][_rootLabel] << std::endl;
        
        for (int s = 0; s < nrInfectedHosts; ++s)
        {
            if (s != _rootLabel)
            {
                _N[root][s] = 0;
            }
        }
    }
    
    _numSols = numSols;
    std::cout << "Number of Sankoff solutions is: " << _numSols << std::endl;
    
    computeWeight(root);
    
    /*
    std::cout << "node\tstate\tweight" << std::endl;
    for (NodeIt vi(_T.tree()); vi != lemon::INVALID; ++vi)
    {
        for (int s = 0; s < nrInfectedHosts; ++s)
        {
            if (_W[vi][s] > 0)
            {
                std::cout << name(vi) << "\t" << s + 1 << "\t" << _W[vi][s] << std::endl;
            }
        }
    }
    */
    
    return numSols;
}

uint64_t Naive::countContact()
{
		Node root = _T.root();
		int nrInfectedHosts = _T.getNrHost();
		
		runContactCount(root);
		
		if (_rootLabel == -1)
		{
				for (int s = 0; s < nrInfectedHosts; ++s)
				{
						if (_N[root][s] > 0)
						{
								_solMap[root] = s;
								break;
						}
				}
		}
		else
		{
				_solMap[root] = _rootLabel;
		}
		
		computeContactSol(root, _solMap[root]);
		
		uint64_t numSols = 0;
		
		if (_rootLabel == -1)
		{
				for (int s = 0; s < nrInfectedHosts; ++s)
				{
						numSols += _N[root][s];
				}
		}
		else
		{
				numSols = _N[root][_rootLabel];
				
				_solMap[root] = _rootLabel;
		}
		
		_numSols = numSols;
		std::cout << "Number of contact-constrained solutions is: " << _numSols << std::endl;
		
		return numSols;
}

bool Naive::solveConsensus()
{
    Node root = _T.root();
    int nrInfectedHosts = _T.getNrHost();
    runCount(root);
		
    uint64_t numSols = 0;
		
    if (_rootLabel == -1)
    {
        uint64_t min_cost = *std::min_element(_M[root].begin(), _M[root].end());
				
        if (min_cost == std::numeric_limits<uint64_t>::max())
        {
            std::cout << "no solutions exist!" << std::endl;
            return false;
        }
        
        std::cout << "feasible roots are: ";
        for (int s = 0; s < nrInfectedHosts; ++s)
        {
            if (_M[root][s] == min_cost)
            {
                std::cout << s + 1 << ", ";
                numSols += _N[root][s];
            }
            else
            {
                _N[root][s] = 0;
            }
        }
        std::cout << std::endl;
    }
    else
    {
        numSols = _N[root][_rootLabel];
        
        if (_M[root][_rootLabel] == std::numeric_limits<uint64_t>::max())
        {
            std::cout << "no solutions exist!" << std::endl;
            return false;
        }
        
        for (int s = 0; s < nrInfectedHosts; ++s)
        {
            if (s != _rootLabel)
            {
                _N[root][s] = 0;
            }
        }
    }
    
    _numSols = numSols;
    std::cout << "Number of Sankoff solutions is: " << _numSols << std::endl;
    
    computeWeight(root);
    
    runConsensus(root);
    
    int consensusRootLabel;
    
    if (_rootLabel == -1)
    {
        consensusRootLabel = std::min_element(_G[root].begin(), _G[root].end()) - _G[root].begin();
    }
    else
    {
        consensusRootLabel = _rootLabel;
    }
    
    _solMap[root] = consensusRootLabel;
    
    if (_M[root][_solMap[root]] == std::numeric_limits<uint64_t>::max())
    {
        std::cout << "no solution exists!" << std::endl;
        return false;
    }
    else
    {
        computeConsensus(root, consensusRootLabel);
        
        int gamma = computeCoInf();
        
        std::cout << "Coinfection number consensus solution is: " << gamma << std::endl;
        

        std::cout << "node\tstate\tweight" << std::endl;
        for (NodeIt vi(_T.tree()); vi != lemon::INVALID; ++vi)
        {
            for (int s = 0; s < nrInfectedHosts; ++s)
            {
                if (_W[vi][s] > 0)
                {
                    std::cout << name(vi) << "\t\t" << s + 1 << "\t\t" << _W[vi][s] << std::endl;
                }
            }
        }
        
        std::cout << "node\tstate\tpenalty" << std::endl;
        for (NodeIt vi(_T.tree()); vi != lemon::INVALID; ++vi)
        {
            for (int s = 0; s < nrInfectedHosts; ++s)
            {
                if (_W[vi][s] > 0)
                {
                    std::cout << name(vi) << "\t\t" << s + 1 << "\t\t" << _G[vi][s] << std::endl;
                }
            }
        }
        
        return true;
    }
    
}


void Naive::run(Node u)
{
    const Digraph& G = tree();
    const int nrInfectedHosts = _T.getNrHost();
    
    if (lemon::countOutArcs(G, u) == 0)
    {
        // is leaf
        for (int s = 0; s < nrInfectedHosts; ++s)
        {
            if (s == _T.getHostLabel(u))
            {
                _M[u][s] = 0;
            }
            else
            {
                _M[u][s] = std::numeric_limits<uint64_t>::max();
            }
        }
    }
    else
    {
        for (Digraph::OutArcIt a(G, u); a != lemon::INVALID; ++a)
        {
            Node v = G.target(a);
            
            run(v);
        }
        
        // compute _M[u][s]
        for (int s = 0; s < nrInfectedHosts; ++s)
        {
            if (_T.getEntTime(s) > _T.getTime(u) || _T.getRemTime(s) < _T.getTime(u))
            {
                _M[u][s] = std::numeric_limits<uint64_t>::max();
            }
            else
            {
                _M[u][s] = 0;
                
                for (Digraph::OutArcIt a(G, u); a != lemon::INVALID; ++a)
                {
                    Node v = G.target(a);
                    
                    uint64_t minimum_cost = std::numeric_limits<uint64_t>::max();
                    for (int t = 0; t < nrInfectedHosts; ++t)
                    {
                        uint64_t cost_vt;
                        
                        if (s != t)
                        {
                            if (_M[v][t] == std::numeric_limits<uint64_t>::max())
                            {
                                cost_vt = std::numeric_limits<uint64_t>::max();
                            }
                            else
                            {
                                cost_vt = _M[v][t] + 1;
                            }
                        }
                        else
                        {
                            cost_vt = _M[v][t];
                        }
                        
                        if (cost_vt < minimum_cost)
                        {
                            minimum_cost = cost_vt;
                        }
                    }
                    
                    _M[u][s] += minimum_cost;
                }
            }
        }
    }
}

void Naive::runCount(Node u)
{
    const Digraph& G = tree();
    const int nrInfectedHosts = _T.getNrHost();
    
    if (lemon::countOutArcs(G, u) == 0)
    {
        // is leaf
        for (int s = 0; s < nrInfectedHosts; ++s)
        {
            if (s == _T.getHostLabel(u))
            {
                _M[u][s] = 0;
                _N[u][s] = 1;
            }
            else
            {
                _M[u][s] = std::numeric_limits<uint64_t>::max();
                _N[u][s] = 0;
            }
        }
    }
    else
    {
        for (Digraph::OutArcIt a(G, u); a != lemon::INVALID; ++a)
        {
            Node v = G.target(a);
            
            runCount(v);
        }
        
        // compute _M[u][s] and _N[u][s]
        for (int s = 0; s < nrInfectedHosts; ++s)
        {
            if (_T.getEntTime(s) > _T.getTime(u) || _T.getRemTime(s) < _T.getTime(u))
            {
                _M[u][s] = std::numeric_limits<uint64_t>::max();
                _N[u][s] = 0;
            }
            else
            {
                _M[u][s] = 0;
                _N[u][s] = 1;
                
                for (Digraph::OutArcIt a(G, u); a != lemon::INVALID; ++a)
                {
                    Node v = G.target(a);
                    uint64_t childPaths = 0;
                    
                    uint64_t minimum_cost = computeMinCost(v, s);
                    
                    _M[u][s] += minimum_cost;
                    
                    for (int t = 0; t < nrInfectedHosts; ++t)
                    {
                        uint64_t cost_vt = computeCost(v, t, s);
                        
                        if (cost_vt == minimum_cost)
                        {
                            childPaths += _N[v][t];
                        }
                    }

                    _N[u][s] *= childPaths;
                }
            }
        }
    }
}

void Naive::runContactCount(Node u)
{
		const Digraph& G = tree();
		const Digraph& C = _T.contactMap();
		const int nrInfectedHosts = _T.getNrHost();
		
		if (lemon::countOutArcs(G, u) == 0)
		{
				// is leaf
				for (int s = 0; s < nrInfectedHosts; ++s)
				{
						if (s == _T.getHostLabel(u))
						{
								_N[u][s] = 1;
						}
						else
						{
								_N[u][s] = 0;
						}
				}
		}
		else
		{
				for (Digraph::OutArcIt a(G,u); a != lemon::INVALID; ++a)
				{
						Node v = G.target(a);
						
						runContactCount(v);
				}
				
				// compute _N[u][s] for all s
				for (int s = 0; s < nrInfectedHosts; ++s)
				{
						if (_T.getEntTime(s) > _T.getTime(u) || _T.getRemTime(s) < _T.getTime(u))
						{
								_N[u][s] = 0;
						}
						else
						{
								_N[u][s] = 1;
								
								for (Digraph::OutArcIt a(G,u); a != lemon::INVALID; ++a)
								{
										Node v = G.target(a);
										uint64_t childPaths = _N[v][s];
										
										Node sHost = _T.getContactNode(s);
										for (Digraph::OutArcIt e(C, sHost); e != lemon::INVALID; ++e)
										{
												Node tHost = C.target(e);
												
												childPaths += _N[v][_T.getContactIndex(tHost)];
										}
										
										_N[u][s] *= childPaths;
								}
						}
				}
		}
}

void Naive::runConsensus(Node u)
{
    const Digraph& G = tree();
    const int nrInfectedHosts = _T.getNrHost();
    
    if (lemon::countOutArcs(G, u) == 0)
    {
        // is leaf
        for (int s = 0; s < nrInfectedHosts; ++s)
        {
            if (s == _T.getHostLabel(u))
            {
                _G[u][s] = 0;
            }
            else
            {
                _G[u][s] = std::numeric_limits<uint64_t>::max();
            }
        }
    }
    else
    {
        for (Digraph::OutArcIt a(G,u); a != lemon::INVALID; ++a)
        {
            Node v = G.target(a);
            
            runConsensus(v);
        }
        
        // compute _G[u][s]
        for (int s = 0; s < nrInfectedHosts; ++s)
        {
            if (_T.getEntTime(s) > _T.getTime(u) || _T.getRemTime(s) < _T.getTime(u))
            {
                _G[u][s] = std::numeric_limits<uint64_t>::max();
            }
            else
            {
                _G[u][s] = _numSols - _W[u][s];
                
                for (Digraph::OutArcIt a(G,u); a !=lemon::INVALID; ++a)
                {
                    Node v = G.target(a);
                    
                    uint64_t minG_cost = std::numeric_limits<uint64_t>::max();
                    
                    uint64_t minimum_cost = computeMinCost(v, s);
                    
                    for (int t = 0; t < nrInfectedHosts; ++t)
                    {
                        uint64_t cost_vt = computeCost(v, t, s);
                                                
                        if (cost_vt == minimum_cost)
                        {
                            if (_G[v][t] < minG_cost)
                            {
                                minG_cost = _G[v][t];
                            }
                        }
                    }
                    _G[u][s] += minG_cost;
                }
            }
        }
    }
}


void Naive::computeSol(Node u, int label)
{
    const Digraph& G = tree();
    const int nrInfectedHosts = _T.getNrHost();
    
    for (Digraph::OutArcIt a(G, u); a != lemon::INVALID; ++a)
    {
        Node v = G.target(a);
        uint64_t minimum_cost = std::numeric_limits<uint64_t>::max();
        
        for (int t = 0; t < nrInfectedHosts; ++t)
        {
            uint64_t cost_vt = computeCost(v, t, label);
            
            if(cost_vt < minimum_cost)
            {
                minimum_cost = cost_vt;
                _solMap[v] = t;
            }
        }
        
        computeSol(v, _solMap[v]);
    }
    
}

void Naive::computeContactSol(Node u, int s)
{
		const Digraph& G = tree();
		const Digraph& C = _T.contactMap();
		
		Node sHost = _T.getContactNode(s);
		
		for (Digraph::OutArcIt a(G, u); a != lemon::INVALID; ++a)
		{
				Node v = G.target(a);

				if (_N[v][s] > 0)
				{
						_solMap[v] = s;
				}
				else
				{
						for (Digraph::OutArcIt e(C, sHost); e != lemon::INVALID; ++e)
						{
								Node tHost = C.target(e);
								int t = _T.getContactIndex(tHost);
								
								if (_N[v][t] > 0)
								{
										_solMap[v] = t;
										break;
								}
						}
				}
				
				computeContactSol(v, _solMap[v]);
		}
}

void Naive::computeConsensus(Node u, int label)
{
    const Digraph& G = tree();
    const int nrInfectedHosts = _T.getNrHost();
    
    for (Digraph::OutArcIt a(G, u); a != lemon::INVALID; ++a)
    {
        Node v = G.target(a);
        
        uint64_t minimum_cost = computeMinCost(v, label);
        uint64_t minG = std::numeric_limits<uint64_t>::max();

        for (int t = 0; t < nrInfectedHosts ; ++t)
        {
            uint64_t cost_vt = computeCost(v, t, label);
            
            if (cost_vt == minimum_cost)
            {
                if (_G[v][t] < minG)
                {
                    minG = _G[v][t];
                    _solMap[v] = t;
                }
            }
        }
        
        computeConsensus(v, _solMap[v]);
    }
}

/* TODO: uniform sampler (Palash)
*/

void Naive::getSample()
{
  const int nrInfectedHosts = _T.getNrHost();
  
  std::uniform_real_distribution<double> distribution(0.0,1.0);
  double p = distribution(g_rng);
  
  double cum = 0;
  for (int s = 0; s < nrInfectedHosts; ++s)
  {
    cum += (double(_N[_T.root()][s]) / double(_numSols));
    if (p <= cum)
    {
      _solMap[_T.root()] = s;
      getSample(_T.root(), s);
      break;
    }
  }
}

void Naive::getSample(Node u, int label)
{
    const Digraph& G = tree();
    const int nrInfectedHosts = _T.getNrHost();
    
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    
    for (Digraph::OutArcIt a(G,u); a != lemon::INVALID; ++a)
    {
        Node v = G.target(a);
        
        uint64_t minimum_cost = computeMinCost(v, label);
        
        uint64_t totalN = 0;
        UInt64IntPairVector feasibleLabels;
        
        for (int t = 0; t < nrInfectedHosts; ++t)
        {
            uint64_t cost_vt = computeCost(v, t, label);
            
            if(cost_vt == minimum_cost)
            {
                totalN += _N[v][t];
                feasibleLabels.push_back(std::make_pair(totalN, t));
            }
        }
        
        double random_number = totalN*distribution(g_rng);
        
        for (int i = 0; i < feasibleLabels.size(); ++i)
        {
            if (random_number < feasibleLabels[i].first)
            {
                _solMap[v] = feasibleLabels[i].second;
				break;
            }
        }
        
        getSample(v, _solMap[v]);
    }
}

void Naive::getContactSample()
{
  const int nrInfectedHosts = _T.getNrHost();
  
  std::uniform_real_distribution<double> distribution(0.0,1.0);
  double p = distribution(g_rng);
  
  double cum = 0;
  for (int s = 0; s < nrInfectedHosts; ++s)
	{
			cum += (double(_N[_T.root()][s]) / double(_numSols));
			if (p <= cum)
			{
					_solMap[_T.root()] = s;
					getContactSample(_T.root(), s);
					break;
			}
	}
}

void Naive::getContactSample(Node u, int s)
{
		const Digraph& G = tree();
		const Digraph& C = _T.contactMap();
		Node sHost = _T.getContactNode(s);
		
		std::uniform_real_distribution<double> distribution(0.0,1.0);
		
		for (Digraph::OutArcIt a(G,u); a != lemon::INVALID; ++a)
		{
				Node v = G.target(a);
				
				uint64_t totalN = 0;
				UInt64IntPairVector feasibleLabels;

				if (_N[v][s] > 0)
				{
						totalN += _N[v][s];
						feasibleLabels.push_back(std::make_pair(totalN, s));
				}
				
				for (Digraph::OutArcIt e(C, sHost); e != lemon::INVALID; ++e)
				{
						Node tHost = C.target(e);
						int t = _T.getContactIndex(tHost);
						
						if (_N[v][t] > 0)
						{
								totalN += _N[v][t];
								feasibleLabels.push_back(std::make_pair(totalN, t));
						}
				}
				
				double random_number = totalN*distribution(g_rng);
				
				for (int i = 0; i < feasibleLabels.size(); ++i)
				{
						if (random_number < feasibleLabels[i].first)
						{
								_solMap[v] = feasibleLabels[i].second;
								break;
						}
				}
				
				getContactSample(v, _solMap[v]);
		}
}

void Naive::computeWeight(Node u)
{
    const Digraph& G = tree();
    const int nrInfectedHosts = _T.getNrHost();
		
    if (lemon::countInArcs(G, u) == 0)
    {
        // is root
        if (_rootLabel == -1)
        {
            uint64_t min_cost = *std::min_element(_M[u].begin(), _M[u].end());
            
            for (int s = 0; s < nrInfectedHosts; ++s)
            {
                if (_M[u][s] == min_cost)
                {
                    _W[u][s] = _N[u][s];
                }
                else
                {
                    _W[u][s] = 0;
                }
            }
        }
        else
        {
            for (int s = 0; s < nrInfectedHosts; ++s)
            {
                if (s == _rootLabel)
                {
                    _W[u][s] = _N[u][s];
                }
                else
                {
                    _W[u][s] = 0;
                }
            }
        }
    }
    
    for (Digraph::OutArcIt a(G, u); a != lemon::INVALID; ++a)
    {
        Node v = G.target(a);
        
        for (int s = 0; s < nrInfectedHosts; ++s)
        {
            _W[v][s] = 0;
        }
        
        for (int s = 0; s < nrInfectedHosts; ++s)
        {
            uint64_t minimum_cost = computeMinCost(v, s);
            
            uint64_t totalN = 0;
            
            for (int t = 0; t < nrInfectedHosts; ++t)
            {
                uint64_t cost_vt = computeCost(v, t, s);
                
                if (cost_vt == minimum_cost)
                {
                    totalN += _N[v][t];
                }
            }
            
            for (int t = 0; t < nrInfectedHosts; ++t)
            {
                uint64_t cost_vt = computeCost(v, t, s);
                
                if (cost_vt == minimum_cost)
                {
                    _W[v][t] = _W[v][t] + _N[v][t]*_W[u][s]/totalN;
                }
            }
        }
        
        computeWeight(v);
    }
}

bool Naive::enumerateSankoff(std::ostream& psol)
{
    Node root = _T.root();
    const int nrInfectedHosts = _T.getNrHost();
    run(root);
    
    // enumeration chart
    // chart[sol_id][(vector of label and nodeID)]
    
    std::cout << "solId\troot\tCoInf" << std::endl;
    
    if (_rootLabel == -1)
    {
        uint64_t minInfNum = *std::min_element(_M[root].begin(), _M[root].end());
        
        for (int s = 0; s < nrInfectedHosts; ++s)
        {
            if (_M[root][s] == minInfNum)
            {
                IntStringPairMatrix solutionEnum = enumerateSol(root, s);
                writeSol(psol, solutionEnum);
            }
        }
    }
    else
    {
        IntStringPairMatrix solutionEnum = enumerateSol(root, _rootLabel);
        writeSol(psol, solutionEnum);
    }
    
    return true;
}

bool Naive::enumerateContact(std::ostream& psol)
{
		Node root = _T.root();
		const int nrInfectedHosts = _T.getNrHost();
		
		std::cout << "SolId\troot\tnInf\tCoInf" << std::endl;
		
		if (_rootLabel)
		{
				for (int s = 0; s < nrInfectedHosts; ++s)
				{
						if (_N[root][s] > 0)
						{
								IntStringPairMatrix solutionEnum = enumerateContactSol(root, s);
								writeContactSol(psol, solutionEnum);
						}
				}
		}
		else
		{
				if (_N[root][_rootLabel] > 0)
				{
						IntStringPairMatrix solutionEnum = enumerateContactSol(root, _rootLabel);
						writeContactSol(psol, solutionEnum);
				}
				else
				{
						return false;
				}
		}
		
		return true;
}

void Naive::writeSol(std::ostream& psol, const IntStringPairMatrix& sol)
{
    uint64_t numSols = sol.size();
    
    for (int i = 0; i < numSols; ++i)
    {
        int numVec = sol[i].size();
        
        for (int j = 0; j < numVec; ++j)
        {
            _solMap[_T.getNodeByLabel(sol[i][j].second)] = sol[i][j].first;
        }
        
        int gamma = computeCoInf();
        std::cout << i << "\t" << _solMap[_T.root()] + 1 << "\t" << gamma << std::endl;
        psol << "# idx = " << i << " -- root = " << _solMap[_T.root()] + 1 << " -- gamma = " << gamma << std::endl;
        
        _T.writePtree(psol, _solMap);
        
    }
}

void Naive::writeContactSol(std::ostream& psol, const IntStringPairMatrix& sol)
{
		uint64_t numSols = sol.size();
		
		for (int i = 0; i < numSols; ++i)
		{
				int numVec = sol[i].size();
				
				for (int j = 0; j < numVec; ++j)
				{
						_solMap[_T.getNodeByLabel(sol[i][j].second)] = sol[i][j].first;
				}
				
				int mu = computeInf();
				int gamma = computeCoInf();
				
				std::cout << i << "\t" << _solMap[_T.root()] + 1 << "\t" << mu << "\t" << gamma << std::endl;
				psol << "# idx = " << i << " -- root = " << _solMap[_T.root()] + 1 << " -- mu = " << mu << " -- gamma = " << gamma << std::endl;
				
				_T.writePtree(psol, _solMap);
				
		}
}

IntStringPairMatrix Naive::enumerateSol(Node u, int label)
{
    const Digraph& G = _T.tree();
    const int nrInfectedHosts = _T.getNrHost();
    
    const int nrChildren_u = lemon::countOutArcs(G, u);
    
    IntStringPairMatrix chart;
    
    // check if node is leaf and if it is then trivially return the
    // leaf label with leaf Id
    if (nrChildren_u == 0)
    {
        IntStringPairVector leafEntry{std::make_pair(_T.getHostLabel(u), _T.getName(u))};
        chart.push_back(leafEntry);
        
        return chart;
    }
    
    std::vector<IntStringPairMatrix> chartChildren(nrChildren_u);
//    IntStringPairMatrix chartChild1;
//    IntStringPairMatrix chartChild2;
    
    int childId = 0;
    
    for (Digraph::OutArcIt a(G,u); a != lemon::INVALID; ++a)
    {
        Node v = G.target(a);
        
        int minimum_cost = computeMinCost(v, label);
        
        for (int t = 0; t < nrInfectedHosts; ++t)
        {
            int cost_vt = computeCost(v, t, label);
            
            if (cost_vt == minimum_cost)
            {
                IntStringPairMatrix childTemp = enumerateSol(v, t);
                
                chartChildren[childId].insert(chartChildren[childId].end(),
                                              childTemp.begin(), childTemp.end());
						}
        }
        
        ++childId;
    }
    
    UInt64Vector sizes(nrChildren_u);
    for (int i = 0; i < nrChildren_u; ++i)
    {
        sizes[i] = chartChildren[i].size();
    }
    
     // merge chartChildren into chart
    UInt64Vector indices(nrChildren_u, 0);
    do
    {
        IntStringPairVector solInstance;
        
        for (int idx = 0; idx < nrChildren_u; ++idx)
        {
            solInstance.insert(solInstance.end(),
                               chartChildren[idx][indices[idx]].begin(),
                               chartChildren[idx][indices[idx]].end());
        }
        
        solInstance.push_back(std::make_pair(label, _T.getName(u)));
        
        chart.push_back(solInstance);
    } while (next(sizes, indices));
    
    // assert that size of chart is consistent
    
    return chart;
}

IntStringPairMatrix Naive::enumerateContactSol(Node u, int s)
{
		const Digraph& G = _T.tree();
		const Digraph& C = _T.contactMap();
		Node sHost = _T.getContactNode(s);
		
		IntStringPairMatrix chart;
		
		int nrChildren_u = lemon::countOutArcs(G, u);
		
		if (nrChildren_u == 0)
		{
				// is leaf
				IntStringPairVector leafEntry{std::make_pair(_T.getHostLabel(u), _T.getName(u))};
				chart.push_back(leafEntry);
				
				return chart;
		}
		
		std::vector<IntStringPairMatrix> chartChildren(nrChildren_u);
		
		int childId = 0;
		
		for (Digraph::OutArcIt a(G, u); a != lemon::INVALID; ++a)
		{
				Node v = G.target(a);
				
				if (_N[v][s] > 0)
				{
						IntStringPairMatrix childTemp = enumerateContactSol(v, s);
						
						chartChildren[childId].insert(chartChildren[childId].end(), childTemp.begin(), childTemp.end());
				}
				
				for (Digraph::OutArcIt e(C, sHost); e != lemon::INVALID; ++e)
				{
						Node tHost = C.target(e);
						int t = _T.getContactIndex(tHost);
						
						if (_N[v][t] > 0)
						{
								IntStringPairMatrix childTemp = enumerateContactSol(v, t);
								
								chartChildren[childId].insert(chartChildren[childId].end(), childTemp.begin(), childTemp.end());
						}
				}
				
				++childId;
		}
		
		UInt64Vector sizes(nrChildren_u);
		for (int i = 0; i < nrChildren_u; ++i)
		{
				sizes[i] = chartChildren[i].size();
		}
		
		// merge chartChildren into chart
		UInt64Vector indices(nrChildren_u, 0);
		do
		{
				IntStringPairVector solInstance;
				
				for (int idx = 0; idx < nrChildren_u; ++idx)
				{
						solInstance.insert(solInstance.end(),
															 chartChildren[idx][indices[idx]].begin(),
															 chartChildren[idx][indices[idx]].end());
				}
				
				solInstance.push_back(std::make_pair(s, _T.getName(u)));
				
				chart.push_back(solInstance);
		} while (next(sizes, indices));
		
		// assert that size of chart is consistent
		
		return chart;
}

bool Naive::next(const UInt64Vector& sizes, UInt64Vector& indices)
{
    assert(sizes.size() == indices.size());
    
    const int n = sizes.size();
    
    // identify first index that can be incremented
    int firstIndex = 0;
    for (; firstIndex < sizes.size(); ++firstIndex)
    {
        if (indices[firstIndex] < sizes[firstIndex] - 1)
            break;
    }
    
    if (firstIndex == n)
        return false;
    
    for (int i = 0; i < firstIndex; ++i)
    {
        indices[i] = 0;
    }
    
    ++indices[firstIndex];
    
    return true;
}

int Naive::computeCoInf()
{
    const Digraph& G = tree();
    const int nrInfectedHosts = _T.getNrHost();
    
    int gamma = 0;
    
    for (int s = 0; s < nrInfectedHosts; ++s)
    {
        for (int t = 0; t < nrInfectedHosts; ++t)
        {
            if (s != t)
            {
                // timeIntervals - (endtime, starttime)
                DoublePairVector timeIntervals;
                
                for (ArcIt eij(G); eij != lemon::INVALID; ++eij)
                {
                    Node vi = G.source(eij);
                    Node vj = G.target(eij);
                    
                    if (_solMap[vi] == s && _solMap[vj] == t) {
                        timeIntervals.push_back(std::make_pair(_T.getTime(vj), _T.getTime(vi)));
                    }
                }
                /*
                if(timeIntervals.size() > 0)
                {
                    std::cout << "infector -> infectee :: " << s << " -> " << t << std::endl;
                
                    std::cout << "number of infection edges: " << timeIntervals.size() << std::endl;
                }
                */
                while (timeIntervals.size() > 0)
                {
                    /*
                    for (int i = 0; i < timeIntervals.size(); ++i)
                    {
                        std::cout << timeIntervals[i].first << ", " << timeIntervals[i].second << std::endl;
                    }
                    */
                    
                    DoublePair minTime = *std::min_element(timeIntervals.cbegin(), timeIntervals.cend());
                    
                    timeIntervals.erase(std::remove_if(timeIntervals.begin(), timeIntervals.end(), [minTime](const DoublePair& x){ return x.second <= minTime.first;}), timeIntervals.end());
                    
                    ++gamma;
                    
                    //std::cout << "current size of batch: " << timeIntervals.size() << std::endl;
                }
            }
        }
    }
    
    return gamma;
}

int Naive::computeInf()
{
		const Digraph& G = tree();
		int mu = 0;
		
		for (ArcIt eij(G); eij != lemon::INVALID; ++eij)
		{
				Node vi = G.source(eij);
				Node vj = G.target(eij);
				
				if (_solMap[vi] != _solMap[vj])
				{
						++mu;
				}
		}
		
		return mu;
}

bool Naive::solve(std::ostream &psol, bool enumerate, uint64_t enumLimit)
{
    // check if root label is out of bound
    // check if root host is possible
    if (_rootLabel >= _T.getNrHost() || _rootLabel < -1)
    {
        std::cerr << "root label is out of bounds" << std::endl;
        return false;
    }
    
    if (_rootLabel > -1)
    {
        if (_T.getEntTime(_rootLabel) > _T.getTime(root()) || _T.getRemTime(_rootLabel) < _T.getTime(root()))
        {
            _T.writePtree(psol, "infeasible root");
            std::cerr << "value of root label is " << _rootLabel + 1 << std::endl;
            std::cerr << "not a feasible root" << std::endl;
            return false;
        }
    }
    
    // check if leaf labels are feasible
    for (NodeIt vi(_T.tree()); vi != lemon::INVALID; ++vi)
    {
        if (lemon::countOutArcs(_T.tree(), vi) == 0)
        {
            int lhost = _T.getHostLabel(vi);
            if (_T.getEntTime(lhost) > _T.getTime(vi) || _T.getRemTime(lhost) < _T.getTime(vi))
            {
                _T.writePtree(psol, "infeasible root labeling in ptree");
                std::cerr << "label of leaf " << _T.getName(vi) << " is wrong" << std::endl;
                return false;
            }
        }
    }
    
    // solve the problem
    if (!enumerate) {
        if (solveSankoff())
        {
            // write the psol
            _T.writePtree(psol, _solMap);
            return true;
        }
        else
        {
            // write empty file
            _T.writePtree(psol, "no solution");
            return false;
        }
    }
    else
    {
        uint64_t nsols = countSankoff();
        
        if (nsols <= enumLimit && nsols > 0)
        {
            return enumerateSankoff(psol);
        }
        else
        {
            return false;
        }
    }
}

bool Naive::solve(std::ostream &psol)
{
    // finds the sankoff consensus solution
    
    // check if root label is out of bound
    // check if root host is possible
    if (_rootLabel >= _T.getNrHost() || _rootLabel < -1)
    {
        std::cerr << "root label is out of bounds" << std::endl;
        return false;
    }
    
    if (_rootLabel > -1)
    {
        if (_T.getEntTime(_rootLabel) > _T.getTime(root()) || _T.getRemTime(_rootLabel) < _T.getTime(root()))
        {
            _T.writePtree(psol, "infeasible root");
            std::cerr << "value of root label is " << _rootLabel + 1 << std::endl;
            std::cerr << "not a feasible root" << std::endl;
            return false;
        }
    }
    
    // check if leaf labels are feasible
    for (NodeIt vi(_T.tree()); vi != lemon::INVALID; ++vi)
    {
        if (lemon::countOutArcs(_T.tree(), vi) == 0)
        {
            int lhost = _T.getHostLabel(vi);
            if (_T.getEntTime(lhost) > _T.getTime(vi) || _T.getRemTime(lhost) < _T.getTime(vi))
            {
                _T.writePtree(psol, "infeasible root labeling in ptree");
                std::cerr << "label of leaf " << _T.getName(vi) << " is wrong" << std::endl;
                return false;
            }
        }
    }
    
    // solve the problem
    if (solveConsensus())
    {
        // write the psol
        _T.writePtree(psol, _solMap);
        return true;
    }
    else
    {
        // write empty file
        _T.writePtree(psol, "no solution");
        return false;
    }
}

bool Naive::solveContact(std::ostream &psol, uint64_t enumLimit)
{
		// finds and enumerates the contact-constrained vertex labels
		
		// check if root label is out of bound
		// check if root host is possible
		if (_rootLabel >= _T.getNrHost() || _rootLabel < -1)
		{
				std::cerr << "root label is out of bounds" << std::endl;
				return false;
		}
		
		if (_rootLabel > -1)
		{
				if (_T.getEntTime(_rootLabel) > _T.getTime(root()) || _T.getRemTime(_rootLabel) < _T.getTime(root()))
				{
						_T.writePtree(psol, "infeasible root");
						std::cerr << "value of root label is " << _rootLabel + 1 << std::endl;
						std::cerr << "not a feasible root" << std::endl;
						return false;
				}
		}
		
		// check if leaf labels are feasible
		for (NodeIt vi(_T.tree()); vi != lemon::INVALID; ++vi)
		{
				if (lemon::countOutArcs(_T.tree(), vi) == 0)
				{
						int lhost = _T.getHostLabel(vi);
						if (_T.getEntTime(lhost) > _T.getTime(vi) || _T.getRemTime(lhost) < _T.getTime(vi))
						{
								_T.writePtree(psol, "infeasible root labeling in ptree");
								std::cerr << "label of leaf " << _T.getName(vi) << " is wrong" << std::endl;
								return false;
						}
				}
		}
		
		uint64_t nsols = countContact();
		
		if (nsols > 0)
		{
				if (enumLimit == 1)
				{
						// write the psol
						_T.writePtree(psol, _solMap);
				}
				else if (nsols <= enumLimit)
				{
						// enumerate all solutions
						enumerateContact(psol);
						return false;
				}
				else
				{
						std::cerr << "number of solutions exceeds enumLimit" << std::endl;
						return false;
				}
		}
		else
		{
				// write empty file
				_T.writePtree(psol, "no solution");
				return false;
		}
		
		return true;
}
