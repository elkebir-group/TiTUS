/*
 * unigenparser.h
 *
 *  Created on: 21-jan-2020
 *      Author: P. Sashittal
 */

#include "utils.h"
#include "basetree.h"
#include <fstream>

class UnigenParser
{
public:
		/// Default constructor
		UnigenParser();
		
		UnigenParser(const BaseTree& T, std::string transPrefix);
		
		UnigenParser(const BaseTree& T, std::string transPrefix, int muMax);

		UnigenParser(const BaseTree& T, std::string transPrefix, int muMax, bool computeDist);
		
		bool readVarFile(std::istream& in);
		
		void parseUnigen(std::istream& in);
  
    void writeTransTrees(std::ostream& out) const;
  
    void writeSummaryStats(std::ostream& out) const;

		void writeSummaryStatsWithDist(std::ostream& out) const;
		
		bool isAncester(const Digraph& T, Node u, Node v);
		
		int getNumUnsampledStrains();
		
		int getDistanceFromContact(IntPairToIntMap currTransmissionCount);
		  
private:
  void initSampledHosts(Node v);
  
protected:
    /// Contact map and ptree
		const BaseTree& _T;
		IntNodeMap _solMap;
  
    Digraph::NodeMap<IntSet> _sampledHostMap;
  
    /// Variable index to <v,s> pair
		StringPairVector _nodeHostPair;
  
    /// Number of transmitted strains between every pair of hosts
    std::vector<std::map<std::pair<int, int>, int>> _transmissionCountVector;
  
    IntVector _solIdx;

    IntVector _transmissionNumber;
    
    IntVector _maxBottleneckSize;
    
    IntVector _unsampledStrains;
		
		IntVector _unsampledLineages;
  
		IntVector _commonEdges;
		
		IntVector _weightedDist;
		
    std::string _transPrefix;
  
		bool _computeDist;
		
		int _numVar;
		
		int _muMax;
};
