/*
 * unigenparser.cpp
 *
 *  Created on: 21-jan-2020
 *      Author: P. Sashittal
 */

#include "unigenparser.h"

UnigenParser::UnigenParser(const BaseTree& T,
													 std::string trans_prefix)
		: _T(T)
		, _solMap(_T.tree())
    , _sampledHostMap(_T.tree())
		, _nodeHostPair()
		, _transPrefix(trans_prefix)
		, _computeDist(false)
		, _numVar()
{
  initSampledHosts(_T.root());
}

UnigenParser::UnigenParser(const BaseTree& T,
													 std::string trans_prefix,
													 int muMax)
: _T(T)
, _solMap(_T.tree())
, _sampledHostMap(_T.tree())
, _nodeHostPair()
, _transPrefix(trans_prefix)
, _computeDist(false)
, _numVar()
, _muMax(muMax)
{
  initSampledHosts(_T.root());
}

UnigenParser::UnigenParser(const BaseTree& T,
													 std::string trans_prefix,
													 int muMax,
													 bool computeDist)
: _T(T)
, _solMap(_T.tree())
, _sampledHostMap(_T.tree())
, _nodeHostPair()
, _transPrefix(trans_prefix)
, _computeDist(computeDist)
, _numVar()
, _muMax(muMax)
{
  initSampledHosts(_T.root());
}

void UnigenParser::initSampledHosts(Node v)
{
  if (OutArcIt(_T.tree(), v) == lemon::INVALID)
  {
    // this is a leaf
    _sampledHostMap[v] = IntSet();
    _sampledHostMap[v].insert(_T.getHostLabel(v));
  }
  else
  {
    _sampledHostMap[v] = IntSet();
    for (OutArcIt a(_T.tree(), v); a != lemon::INVALID; ++a)
    {
      Node w = _T.tree().target(a);
      initSampledHosts(w);
      _sampledHostMap[v].insert(_sampledHostMap[w].begin(), _sampledHostMap[w].end());
    }
  }
}

bool UnigenParser::readVarFile(std::istream& in)
{
		int numvar = 0;
		
		while (in.good())
		{
				std::string line;
				getline(in, line);
				
				if (line.empty())
						break;
				
				StringVector s;
				boost::split(s, line, boost::is_any_of("\t "));
				
				if (s.size() != 4)
				{
						std::cerr << "Error: line '" << line << "' incorrect number of columns" << std::endl;
						return false;
				}
				
				if (s[3] == "(v,s)")
				{
						_nodeHostPair.push_back(std::make_pair(s[1], s[2]));
						++numvar;
				}
				else
				{
						break;
				}
		}
		
		_numVar = numvar;
  
		return true;
}

void UnigenParser::writeTransTrees(std::ostream& out) const
{
    out << _transmissionCountVector.size() << " # trans trees" << std::endl;
    for (int idx = 0; idx < _transmissionCountVector.size(); ++idx)
    {
        out << _transmissionCountVector[idx].size() << " # edges, tree " << idx << std::endl;
        for (const auto& kv : _transmissionCountVector[idx])
        {
            int s = kv.first.first;
            int t = kv.first.second;
          
            const std::string& sStr = _T.getContactName(s);
            const std::string& tStr = _T.getContactName(t);
            out << sStr << " " << tStr << " " << kv.second << std::endl;
        }
    }
}

void UnigenParser::parseUnigen(std::istream& in)
{
    _transmissionCountVector.clear();
		const Digraph& G = _T.tree();
		const Digraph& C = _T.contactMap();
		
    int solIdx = 0;
		while (in.good())
		{
				std::string line;
				getline(in, line);
				
				if (line.empty())
				{
						break;
				}
				
				for (NodeIt vi(G); vi != lemon::INVALID; ++vi)
				{
						if (lemon::countOutArcs(G, vi) == 0)
						{
								_solMap[vi] = _T.getHostLabel(vi);
						}
						else
						{
								_solMap[vi] = -1;
						}
				}
				
				StringVector s;
				boost::split(s, line, boost::is_any_of("\t "));
				
				for (int i = 0; i < _numVar; ++i)
				{
						std::string varId = s[i];
						
						if (varId.at(0) != '-')
						{
								int varIndex = std::stoi(varId) - 1;
								
								std::string nodeName = _nodeHostPair[varIndex].first;
								std::string hostName = _nodeHostPair[varIndex].second;
								
								Node u = _T.getNodeByLabel(nodeName);
								Node hostNode = _T.getContactNode(hostName);
								int hostIndex = _T.getContactIndex(hostNode);
								
								_solMap[u] = hostIndex;
						}
				}
      
        // we have a solution now, update transmission counts
				std::map<std::pair<int, int>, int> currTransmissionCount;
				int mu = 0;
        int lambda = 0;
        for (ArcIt a_ij(G); a_ij != lemon::INVALID; ++a_ij)
        {
          Node v_i = G.source(a_ij);
          Node v_j = G.target(a_ij);
						
          int s = _solMap[v_i];
          assert(s != -1);
          int t = _solMap[v_j];
          assert(t != -1);
          
          if (s == t) continue;
          
          IntPair st = std::make_pair(s, t);
						
          if (_sampledHostMap[v_j].count(t) == 0)
          {
            ++lambda;
          }
				
				  if (currTransmissionCount.count(st) == 0)
					{
							currTransmissionCount[st] = 0;
					}
				  ++currTransmissionCount[st];
				
				  ++mu;
        }
				
				int commonEdges = 0;
				// compute number of common edges
				for (ArcIt est(C); est != lemon::INVALID; ++est)
				{
						int sIndex = _T.getContactIndex(C.source(est));
						int tIndex = _T.getContactIndex(C.target(est));
						
						IntPair stContact = std::make_pair(sIndex, tIndex);
						
						if (currTransmissionCount.count(stContact) > 0)
						{
								++commonEdges;
						}
				}
				
				int weightedDist = getDistanceFromContact(currTransmissionCount);
				
				if (mu <= _muMax)
				{
						_transmissionCountVector.push_back(currTransmissionCount);
            _transmissionNumber.push_back(mu);
            _solIdx.push_back(solIdx);
            _maxBottleneckSize.push_back(0);
            for (const auto& kv : currTransmissionCount)
            {
                _maxBottleneckSize.back() = std::max(_maxBottleneckSize.back(),
                                                     kv.second);
            }
          
          _unsampledLineages.push_back(lambda);
						
				  _unsampledStrains.push_back(getNumUnsampledStrains());
						
				  _commonEdges.push_back(commonEdges);
						
				  _weightedDist.push_back(weightedDist);
				}
        ++solIdx;
		}
}

void UnigenParser::writeSummaryStats(std::ostream& out) const
{
  out << "solIdx" << "\t" << "transmissions" << "\t" << "maxBottleneckSize" << "\t" << "unsampledLineages" << "\t" << "unsampledStrains" << "\t" << "numCommonEdges" << std::endl;
  
  const int nrSols = _solIdx.size();
  for (int idx = 0; idx < nrSols; ++idx)
  {
    out << _solIdx[idx] << "\t" << _transmissionNumber[idx] << "\t" << _maxBottleneckSize[idx] << "\t" << _unsampledLineages[idx] << "\t" << _unsampledStrains[idx] << "\t" << _commonEdges[idx] << std::endl;
  }
}

void UnigenParser::writeSummaryStatsWithDist(std::ostream& out) const
{
  out << "solIdx" << "\t" << "transmissions" << "\t" << "maxBottleneckSize" << "\t" << "unsampledLineages" << "\t" << "unsampledStrains" << "\t" << "numCommonEdges" << "\t" << "weightedDist" << std::endl;
  
  const int nrSols = _solIdx.size();
  for (int idx = 0; idx < nrSols; ++idx)
	{
			out << _solIdx[idx] << "\t" << _transmissionNumber[idx] << "\t" << _maxBottleneckSize[idx] << "\t" << _unsampledLineages[idx] << "\t" << _unsampledStrains[idx] << "\t" << _commonEdges[idx] << "\t" << _weightedDist[idx] << std::endl;
	}
}

bool UnigenParser::isAncester(const Digraph& T, Node u, Node v)
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

int UnigenParser::getNumUnsampledStrains()
{
		const Digraph& G = _T.tree();
		std::vector<NodeSet> leafSetperHost(_T.getNrHost());
		
		for (NodeIt vi(G); vi != lemon::INVALID; ++vi)
		{
				if (lemon::countOutArcs(G, vi) == 0)
				{
						// is leaf
						leafSetperHost[_solMap[vi]].insert(vi);
				}
		}
		
		int gamma = 0;
		
		for (NodeIt vi(G); vi != lemon::INVALID; ++vi)
		{
				bool isUnsampled = true;
				
				for (Node leaf : leafSetperHost[_solMap[vi]])
				{
						if (isAncester(G, vi, leaf))
						{
								isUnsampled = false;
								break;
						}
				}
				
				if (isUnsampled)
				{
						++gamma;
				}
		}
		
		return gamma;
}

int UnigenParser::getDistanceFromContact(IntPairToIntMap currTransmissionCount)
{
		const Digraph& C = _T.contactMap();
		int countContact = 0, countTransmission = 0;
		int countIntersect = 0;
		
		for (ArcIt est(C); est != lemon::INVALID; ++est)
		{
				int sIndex = _T.getContactIndex(C.source(est));
				int tIndex = _T.getContactIndex(C.target(est));
				
				IntPair stContact = std::make_pair(sIndex, tIndex);
				
				if (currTransmissionCount.count(stContact) > 0)
				{
						countIntersect += std::min(currTransmissionCount[stContact], _T.getContactNumStrains(est));
				}
				
				countContact += _T.getContactNumStrains(est);
		}
		
		for (const auto& kv : currTransmissionCount)
		{
				countTransmission += kv.second;
		}
		
		return countContact + countTransmission - 2*countIntersect;
}
