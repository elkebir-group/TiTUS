/*
 * ilpsolver.h
 *
 *  Created on: 9-nov-2018
 *      Author: P. Sashittal
 */

#include "utils.h"
#include "basetree.h"
#include <fstream>

class Naive
{
public:
    Naive(const BaseTree& T, const int rootLabel);
    
    /// initialize the solver
    void init();
		
		/// initialize the solver for contact constrained solutions
		void initContact();
		
    /// call the dynamic program (solve or enumerate)
    /// and read the solution to solution map
    virtual bool solve(std::ostream& psol, bool enumerate, uint64_t enumLimit);
    
    /// call the dynamic program to get concensus solution
    virtual bool solve(std::ostream& psol);
  
		/// call the dynamic program to enumerate contact-constrained solution
		virtual bool solveContact(std::ostream& psol, uint64_t enumLimit);
		
    const IntNodeMap& getSolMap() const
    {
        return _solMap;
    }
 
    /// solve for co-infection number
    /// for given solMap
    int computeCoInf();
		
		// solve for infection number for given solMap
		int computeInf();
  
    /// get sample from uniform distribution
    /// over sankoff solution space
    void getSample();
		
		/// get sample from uniform distribution
		/// over contact-constrained solution space
		void getContactSample();
		
protected:
    bool next(const UInt64Vector& sizes, UInt64Vector& indices);

		/// contact-constrained functions
		uint64_t countContact();
		
		void computeContactSol(Node u, int label);
		
		bool enumerateContact(std::ostream& psol);
		
		IntStringPairMatrix enumerateContactSol(Node u, int label);
		
		void runContactCount(Node u);

		void writeContactSol(std::ostream& psol, const IntStringPairMatrix& sol);
		
		void getContactSample(Node u, int label);
		
    /// solve Sankoff
    bool solveSankoff();
		
    /// Solve for the parsimony score at given node
    /// (dynamic programming)
    /// @param u Node
    void run(Node u);
    
    /// solve for the optimal label of children
    /// for a given node and label
    /// (dynamic programming)
    /// $param u Node
    void computeSol(Node u, int label);
  
    /// count the number of Sankoff solutions
    uint64_t countSankoff();
		
    /// Solve for the number of Sankoff solutions
    /// (dynamic programming)
    /// @param u Node
    void runCount(Node u);
		
    /// compute the number of instances of Sankoff solutions
    /// with particular node label
    /// @param u Node
    void computeWeight(Node u);
    
    /// solve the sankoff consensus problem
    bool solveConsensus();

    /// solve for the consensus penalty chart
    /// (dynamic programming)
    /// @param u Node
    void runConsensus(Node u);
    
    /// solve for the optimal label for concensus
    //// @param u Node
    //// @param label Host
    void computeConsensus(Node u, int label);
    
    /// get sample from uniform distribution
    /// over sankoff solution space
    void getSample(Node u, int label);
    
    /// enumerate Sankoff
    bool enumerateSankoff(std::ostream& psol);
		
    /// write the enumerated solution
    void writeSol(std::ostream& psol, const IntStringPairMatrix& sol);
    
    /// compute subtree enumeration solution
    /// (dynamic programming)
    IntStringPairMatrix enumerateSol(Node u, int label);
    
    /// compute cost_vt
    uint64_t computeCost(Node v, int label, int plabel)
    {
        uint64_t cost_vt;
        
        if (plabel != label)
        {
            if (_M[v][label] == std::numeric_limits<uint64_t>::max())
            {
                cost_vt = std::numeric_limits<uint64_t>::max();
            }
            else
            {
                cost_vt = _M[v][label] + 1;
            }
        }
        else
        {
            cost_vt = _M[v][label];
        }
        
        return cost_vt;
    }
    
    /// compute minimum cost
    uint64_t computeMinCost(Node v, int plabel)
    {
        int nrInfectedHosts = _T.getNrHost();
        uint64_t minimum_cost = std::numeric_limits<uint64_t>::max();
        
        for (int t = 0; t < nrInfectedHosts; ++t)
        {
            uint64_t cost_vt = computeCost(v, t, plabel);
            
            if (cost_vt < minimum_cost)
            {
                minimum_cost = cost_vt;
            }
        }
        
        return minimum_cost;        
    }
    
    /// Return underlying LEMON tree
    virtual const Digraph& tree() const
    {
        return _T.tree();
    }
    
    /// Return the index of a node
    virtual const int index(Node u) const
    {
        return _T.getIndex(u);
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
    virtual const double enttime(int s) const
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
    const BaseTree& _T;
    /// root label
    const int _rootLabel;
    /// Solution vertex labeling
    IntNodeMap _solMap;
    /// dynamic programing chart (parsimony score)
    UInt64VectorNodeMap _M;
    /// number of paths chart
    UInt64VectorNodeMap _N;
    /// number of instances chart
    UInt64VectorNodeMap _W;
    /// penalty score chart
    UInt64VectorNodeMap _G;
    /// number of Sankoff solutions
    uint64_t _numSols;
		
		//lemon::ArcLookUp<Digraph> _contactLookUp;
};
