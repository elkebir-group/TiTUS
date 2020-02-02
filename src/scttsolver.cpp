/*
 * scttsolver.cpp
 *
 *  Created on: 23-dec-2019
 *      Author: P. Sashittal
 */

#include "scttsolver.h"
#include <lemon/min_cost_arborescence.h>

ScttSolver::ScttSolver()
: _contactMap()
, _nodeToId(_contactMap)
, _idToNode()
, _nodeToIndex(_contactMap)
, _indexToNode()
, _nhosts()
, _enttime()
, _remtime()
, _unsampledHosts()
, _treeEdgeWeights(_contactMap)
, _edgeLabel(_contactMap)
, _edgeCost(_contactMap)
, _mdst(_contactMap)
{
		_unsampledHosts = 0;
}

ScttSolver::ScttSolver(int unhosts)
: _contactMap()
, _nodeToId(_contactMap)
, _idToNode()
, _nodeToIndex(_contactMap)
, _indexToNode()
, _nhosts()
, _enttime()
, _remtime()
, _unsampledHosts(unhosts)
, _treeEdgeWeights(_contactMap)
, _edgeLabel(_contactMap)
, _edgeCost(_contactMap)
, _mdst(_contactMap)
{
}


bool ScttSolver::readHost(std::istream& in)
{
		int idx = 0;
		_hostLabel.clear();
		while (in.good())
		{
				std::string line;
				getline(in, line);

				if (line.empty())
						break;

				StringVector s;
				boost::split(s, line, boost::is_any_of("\t "));

				if (s.size() != 3)
				{
						std::cerr << "Error: line '" << line << "' incorrect number of columns in host file" << std::endl;
						return false;
				}

				double etime = std::stod(s[1]);
				double rtime = std::stod(s[2]);

				_hostLabel.push_back(s[0]);
				_enttime.push_back(etime);
				_remtime.push_back(rtime);

				idx++;
		}

		int totalTime = *std::max_element(_remtime.begin(), _remtime.end());

		/// add unsampled hosts
		for (int s = 0; s < _unsampledHosts; ++s)
		{
				_hostLabel.push_back("Unsampled" + std::to_string(s+1));

				_enttime.push_back(0);
				_remtime.push_back(totalTime);
				idx++;
		}

		// initialize nhosts
		_nhosts = idx;

		std::cout << "number of hosts is: " << _nhosts << std::endl;
		/*
		 std::cout << "removal times are: " << std::endl;
		 for (int s = 0; s < _nhosts; ++s)
		 {
		 std::cout << "host " << s+1 << " has " << _enttime[s] << " and " << _remtime[s] << std::endl;
		 }
		 */
		std::cout<< "total time is " << totalTime << std::endl;

		return true;
}

bool ScttSolver::readContactMap(std::istream& in)
{
		for (int i = 0; i < _nhosts; ++i)
		{
				Node u = _contactMap.addNode();
				_nodeToIndex[u] = i;
				_indexToNode.push_back(u);
				_nodeToId[u] = _hostLabel[i];
				_idToNode[_hostLabel[i]] = u;
		}

		while (in.good())
		{
				std::string line;
				getline(in, line);

				if (line.empty())
						break;

				StringVector s;
				boost::split(s, line, boost::is_any_of("\t "));

				if (s.size() != 2)
				{
						std::cerr << "Error: line of contact map'" << line << "' incorrent number of columns" << std::endl;
						return false;
				}

				std::string source = s[0];
				std::string target = s[1];

				int sourceIndex, targetIndex;

				// find source Index
				StringVector::iterator itr = std::find(_hostLabel.begin(), _hostLabel.end(), source);
				if (itr != _hostLabel.cend())
				{
						sourceIndex = std::distance(_hostLabel.begin(), itr);
				}
				else
				{
						std::cout << "Error: line'" << line << "' hostlabel " << source << " not found in the host file" << std::endl;
						return false;
				}

				// find target Index
				itr = std::find(_hostLabel.begin(), _hostLabel.end(), target);
				if (itr != _hostLabel.cend())
				{
						targetIndex = std::distance(_hostLabel.begin(), itr);
				}
				else
				{
						std::cout << "Error: line'" << line << "' hostlabel " << target << " not found in the host file" << std::endl;
						return false;
				}

				if (sourceIndex == targetIndex)
				{
						std::cout << "Error: line'" << line << "' hostlabel " << source << " and hostLabel " << target << " return the same index" << std::endl;
						return false;
				}

				// if the time-intervals intersect, then add a directed arc from source-> target
				if (_remtime[sourceIndex] >= _enttime[targetIndex] && _enttime[sourceIndex] <= _remtime[targetIndex])
				{
						Node sourceHostNode = _idToNode.find(source)->second;
						Node targetHostNode = _idToNode.find(target)->second;

						_contactMap.addArc(sourceHostNode, targetHostNode);
				}
		}

		// we can now add the unsampled hosts
		for (int s = 0; s < _unsampledHosts; ++s)
		{
				Node unsampledHostNode = _idToNode.find("Unsampled" + std::to_string(s + 1))->second;

				for (int t = 0; t < _nhosts - _unsampledHosts; ++t)
				{
						Node sampledHostNode = _idToNode.find(_hostLabel[t])->second;

						_contactMap.addArc(unsampledHostNode, sampledHostNode);
						_contactMap.addArc(sampledHostNode, unsampledHostNode);
				}

				for (int t = 0; t < _unsampledHosts; ++t)
				{
						if (t != s)
						{
								Node anotherUnsampledHostNode = _idToNode.find("Unsampled" + std::to_string(t + 1))->second;

								_contactMap.addArc(unsampledHostNode, anotherUnsampledHostNode);
								_contactMap.addArc(anotherUnsampledHostNode, unsampledHostNode);
						}
				}
		}

		/*
		// print whole contact map
		std::cout << "contact map is " << std::endl;
		for (ArcIt eij(_contactMap); eij != lemon::INVALID; ++eij)
		{
				std::cout << _nodeToIndex[_contactMap.source(eij)] << " -> " << _nodeToIndex[_contactMap.target(eij)] << std::endl;
		}
		*/

		return true;
}

void ScttSolver::setFullContactMap()
{
		for (int i = 0; i < _nhosts; ++i)
		{
				Node u = _contactMap.addNode();
				_nodeToIndex[u] = i;
				_indexToNode.push_back(u);
				_nodeToId[u] = _hostLabel[i];
				_idToNode[_hostLabel[i]] = u;
		}

		for (int s = 0; s < _nhosts; ++s)
		{
				for (int t = 0; t < s; ++t)
				{
						if (_remtime[s] >= _enttime[t] && _enttime[s] <= _remtime[t])
						{
								Node host1 = _idToNode.find(_hostLabel[s])->second;
								Node host2 = _idToNode.find(_hostLabel[t])->second;

								_contactMap.addArc(host1, host2);
								_contactMap.addArc(host2, host1);
						}
				}
		}

		// print whole contact map
		std::cout << "contact map (interval graph) is " << std::endl;
		for (ArcIt eij(_contactMap); eij != lemon::INVALID; ++eij)
		{
				std::cout << _nodeToId[_contactMap.source(eij)] << " -> " << _nodeToId[_contactMap.target(eij)] << std::endl;
		}
}

bool ScttSolver::readTransmissionTrees(std::istream& in)
{
  // add to _treeEdgeWeights vector;
  lemon::ArcLookUp<Digraph> contactLookUp(_contactMap);

  std::string line;
  getline(in, line);

  std::stringstream ss(line);
  int nrTrees = -1;
  ss >> nrTrees;
  if (nrTrees < 0)
  {
    std::cerr << "error: incorrect number of transmission trees" << std::endl;
    return false;
  }

  for (int treeIdx = 0; treeIdx < nrTrees; ++treeIdx)
  {
    getline(in, line);

    std::stringstream ss(line);
    int nrEdges = -1;

    ss >> nrEdges;
    if (nrEdges < 0)
    {
      std::cerr << "error: incorrect number of edges for tree " << treeIdx << std::endl;
      return false;
    }

    for (int edgeIdx = 0; edgeIdx < nrEdges; ++edgeIdx)
    {
      getline(in, line);

      StringVector s;
      boost::split(s, line, boost::is_any_of("\t "));

      if (s.size() != 3)
      {
        std::cerr << "error: line '" << line << "' incorrect number of columns in ttree file" << std::endl;
        return false;
      }

      std::string source = s[0];
      std::string target = s[1];

      int weight = std::stoi(s[2]);

      int sourceIndex, targetIndex;

      StringVector::iterator itr = std::find(_hostLabel.begin(), _hostLabel.end(), source);
      if (itr != _hostLabel.cend())
      {
        sourceIndex = std::distance(_hostLabel.begin(), itr);
      }
      else
      {
        std::cout << "Error: line'" << line << "' hostlabel " << source << " not found in the host file" << std::endl;
        return false;
      }

      // find target Index
      itr = std::find(_hostLabel.begin(), _hostLabel.end(), target);
      if (itr != _hostLabel.cend())
      {
        targetIndex = std::distance(_hostLabel.begin(), itr);
      }
      else
      {
        std::cout << "Error: line'" << line << "' hostlabel " << target << " not found in the host file" << std::endl;
        return false;
      }

      if (sourceIndex == targetIndex)
      {
        std::cout << "Error: line'" << line << "' hostlabel " << source << " and hostLabel " << target << " return the same index" << std::endl;
        return false;
      }


      Node sourceNode = _indexToNode[sourceIndex];
      Node targetNode = _indexToNode[targetIndex];

      Arc contactArc = contactLookUp(sourceNode, targetNode);

      if (contactArc == lemon::INVALID)
      {
        std::cerr << "Warning: line '" << line << "' has arc that does not belong to contact map" << std::endl;
      }
      else
      {
        _treeEdgeWeights[contactArc].push_back(weight);
      }
    }
  }

		return true;
}

void ScttSolver::buildParentChildGraph()
{
		for (ArcIt a(_contactMap); a != lemon::INVALID; ++a)
		{
				IntVector treeWeights = _treeEdgeWeights[a];

				int nsize = treeWeights.size();

				if (nsize == 0)
				{
						_edgeLabel[a] = 1;
						// should I instead remove this edge?
				}
				else
				{
						int nhalf;

						if (nsize%2 == 0)
						{
								nhalf = nsize/2 - 1;
						}
						else
						{
								nhalf = nsize/2;
						}

						std::nth_element(treeWeights.begin(), treeWeights.begin() + nhalf, treeWeights.end());

						int solWeight = treeWeights[nhalf];

						/*
						Node u = _contactMap.source(a);
						Node v = _contactMap.target(a);
						std::cout << "optimal weight for " << _nodeToId[u] << " -> " << _nodeToId[v] << " is " << solWeight << " with nsize = " << nsize << " and nhalf = " << nhalf << std::endl;
						*/

						if (solWeight > 0)
						{
								_edgeLabel[a] = solWeight;
						}
						else
						{
								_edgeLabel[a] = 1;
						}

				}

				_edgeCost[a] = 0;
				for (int i = 0; i < nsize; i++)
				{
						_edgeCost[a] += std::abs(treeWeights[i] - _edgeLabel[a]) - treeWeights[i];
				}
		}

		/*
		for (ArcIt e(_contactMap); e != lemon::INVALID; ++e)
		{
				Node u = _contactMap.source(e);
				Node v = _contactMap.target(e);

				std::cout << _nodeToId[u] << " -> " << _nodeToId[v] << " " << _edgeLabel[e] << " " << _edgeCost[e] << std::endl;
		}
		*/
}

void ScttSolver::solve()
{
		int minCost = std::numeric_limits<int>::max();

		for (NodeIt v(_contactMap); v != lemon::INVALID; ++v)
		{
				BoolArcMap mdst(_contactMap, false);
				lemon::Bfs<Digraph> bfs(_contactMap);

				bfs.run(v);

				bool chosen = true;
				for (NodeIt u(_contactMap); u != lemon::INVALID; ++u)
				{
						if (!bfs.reached(u))
						{
								chosen = false;
								break;
						}
				}

				if (chosen)
				{
						int spanningCost = lemon::minCostArborescence(_contactMap, _edgeCost, v, mdst);

						if (spanningCost < minCost)
						{
								minCost = spanningCost;

								for (ArcIt e(_contactMap); e != lemon::INVALID; ++e)
								{
										_mdst[e] = mdst[e];
								}
						}
				}
		}
}

void ScttSolver::writeConsensusTree(std::ostream& out, const std::string& msg) const
{
		out << msg;
}

void ScttSolver::writeConsensusTree(std::ostream& out) const
{
		for (ArcIt e(_contactMap); e != lemon::INVALID; ++e)
		{
				if (_mdst[e])
				{
						Node u = _contactMap.source(e);
						Node v = _contactMap.target(e);

						out << _nodeToId[u] << "\t" << _nodeToId[v] << "\t" << _edgeLabel[e] << std::endl;
				}
		}
}

void ScttSolver::writeDot(std::ostream& out) const
{
		StringVector colorMap({
				"#3243BA",
				"#0363E1",
				"#0D75DC",
				"#1485D4",
				"#0998D1",
				"#06A7C6",
				"#15B1B4",
				"#38B99E",
				"#65BE86",
				"#92BF73",
				"#B7BD64",
				"#D9BA56",
				"#F8BB44",
				"#FCCE2E",
				"#F5E41D",
				"#F9FB0E"
		});

		out << "digraph N {" << std::endl;
		out << "\toverlap=\"false\"" << std::endl;
		out << "\trankdir=\"LR\"" << std::endl;
		for (int s = 0; s < _nhosts; ++s)
		{
				if (_nhosts < colorMap.size()) {
						out << "\t" << s << " [label=\"" << _hostLabel[s] << "\",width=1.2,height=1.2,style=\"\",penwidth=3,color=\"" << colorMap[int (s * colorMap.size() / _nhosts)] << "\"]" << std::endl;
				}
				else
				{
						out << "\t" << s << " [label=\"" << _hostLabel[s] << "\",width=1.2,height=1.2,style=\"\",penwidth=3,color=\"" << colorMap[s % colorMap.size()] << "\"]" << std::endl;
				}
		}

		for (ArcIt e(_contactMap); e != lemon::INVALID; ++e)
		{
				if (_mdst[e])
				{
						int s = _nodeToIndex[_contactMap.source(e)];
						int t = _nodeToIndex[_contactMap.target(e)];

						out << "\t" << s << " -> " << t << " [penwidth=1,color=black,label=\"" << _edgeLabel[e] << "\"]" << std::endl;
				}
		}

		out << "}" << std::endl;
}
