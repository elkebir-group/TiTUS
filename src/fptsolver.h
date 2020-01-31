/*
 * fptsolver.h
 *
 *  Created on: 21-dec-2019
 *      Author: P. Sashittal
 */

#include "basetree.h"
#include "transmissiontree.h"
#include "utils.h"
#include <fstream>

class FptSolver
{
public:
		
    FptSolver(const BaseTree& T, const int rootLabel, std::string psol);
		
		/// full solver using FPT
		bool solveFull();
		
		bool enumeratedSolve();
		
		bool run(int s);
		
		bool grow(Node r, SubDigraph& G, SubDigraph& T, SubDigraph& L, SubBfs& bfsL, ArcList& F);
		
    /// solver using FPT
    bool solve();
		
		bool solve(TransmissionTree S);
		
		static Node getLCA(const Digraph& G, const NodeSet& nodes);
		
		static NodeList pathFromRoot(const Digraph& G, Node node);
		
		static bool isAncestor(const Digraph& G, Node u, Node v);
		
		static bool isLeaf(const Digraph& G, const Node u);
		
		static bool isChild(const Digraph& G, const Node v, const Node u);
		
		/// fix LCA labels
		bool fixLCALabels(Node q);
		
		void fixLCALabels(Node q, TransmissionTree S) const;
		
		/// seed the child trees
		bool seedChildTrees(Node r) const;
		
		/// auxilary function for seed the child trees function
		void auxSeedChildTrees(Node u, int q, int r, double infTime) const;
		
		bool auxSeedChildTrees(Node u, Node q, Node r, double infTime) const;
		
		bool isFeasibleLabel(const Node u, const Node q) const;
		
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
		
		// the timed phylogeny and the epi data
		const BaseTree& _T;
		
		TransmissionTree _tTree;
		
		IntPairSetSet _allSelectedPairs;
		
		// root label
		const int _rootLabel;
		
    // output file
		std::ofstream _psolFile;
		
		// solution labeling
		mutable IntNodeMap _solMap;
};
