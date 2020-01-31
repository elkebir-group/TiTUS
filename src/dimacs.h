/*
 * dimacs.h
 *
 *  Created on: 9-nov-2018
 *      Author: P. Sashittal
 */

#include "basetree.h"
#include "utils.h"
#include <fstream>

class Dimacs
{
public:
    Dimacs(const BaseTree& T,
           const int rootLabel,
           std::string psol,
           std::string varlist);
    
    /// initialize the exporter
    void initVariables();
    
    /// solver for dimacs
    virtual bool solveDimacs();
    
    /// variable list writing
    void writeVariableList();
		
    /// vertex labeling time constraints
    void writeVertexTimeConstraints();
		
		/// vertex labeling contact map constraints
		void writeVertexContactConstraints();
		
		/// rooted direct infection
		void writeRootedDirectInfectionConstraints(int rootLabel);
		
		/// unrooted direct infection
		void writeUnRootedDirectInfectionConstraint();
		
    /// write oneHot constraints (node)
    void writeOneHotConstraints(Node u, BoolVector feasibleLabels);
		
		/// write oneHot constraints (node)
		void writeOneHotConstraints(Node u, IntVector feasibleLabels);

		/// write oneHost constraints (contact)
		void writeOneHotConstraints(int t, IntVector feasibleSources);
		    
    /// write time constraints
    void writeInfectionTimeConstraints();
		
		/// write infection window constraints
		void writeInfectionWindowConstraints();
		
		/// write strong bottleneck constraints
		void writeStrongBottleneckCosntraints();
		
    /// write independent support
    void writeIndependentSupport();
		
		void setStrongBottleneck()
		{
				_bottleneck = true;
		}
		
    /// get variable index
    int getVarIndex(Node u, int label)
    {
        assert(_Rv[u][label] > 0);
        return _Rv[u][label];
    }
		
		/// get variable index
		int getVarIndex(int s, int t)
		{
				assert(_Rst[s][t] > 0);
				return _Rst[s][t];
		}
		
    /// get feasible labels
    IntVector getFeasibleLables(Node u)
    {
        IntVector feasibleLabels;
        
        for (int s = 0; s < _T.getNrHost(); ++s)
        {
            if (_Rv[u][s] > 0)
            {
                feasibleLabels.push_back(s);
            }
        }
        
        return feasibleLabels;
    }
		
		
		// BaseTree functions
		
		/// Return underlying LEMON tree
		virtual const Digraph& tree() const
		{
				return _T.tree();
		}
		
		/// Return the underlying contact map
		virtual const Digraph& contactMap() const
		{
				return _T.contactMap();
		}
		
		/// Return the index of a node
		virtual const int index(Node u) const
		{
				return _T.getIndex(u);
		}
		
		virtual const int contactIndex(Node u) const
		{
				return _T.getContactIndex(u);
		}
		
		/// Return the index of an arc
		virtual const int arc_index(Arc v) const
		{
				return _T.getArcIndex(v);
		}
		
		/// Return node identifier
		///
		/// @param v Node in T
		virtual const std::string& name(Node v) const
		{
				return _T.getName(v);
		}
		
		/// Return label of a node
		virtual const int label(Node u) const
		{
				return _T.getHostLabel(u);
		}
		
		/// Return timestamp of a node
		virtual const double time(Node u) const
		{
				return _T.getTime(u);
		}
		
		/// Return entry time of a host
		virtual const double entTime(int s) const
		{
				return _T.getEntTime(s);
		}
		
		/// Return removal time of a host
		virtual const double remTime(int s) const
		{
				return _T.getRemTime(s);
		}
		
		/// Return root node
		virtual Node root() const
		{
				return _T.root();
		}
		
protected:
		// the Tree and the epi data
		const BaseTree& _T;
		// root label
		const int _rootLabel;
		// bottleneck constraint
		bool _bottleneck;
		
    // vaiable rank in the list of
    // variables for SAT problem (vertices and edges)
    IntVectorNodeMap _Rv;
		IntMatrix _Rst;
		BoolMatrix _contact;
  
    // clauses
    IntMatrix _clauses;
		
		// number of boolean variables
    int _numVar;
		
		// output DIMACS file
    std::ofstream _solFile;
    // output variable file
    std::ofstream _varFile;
};
