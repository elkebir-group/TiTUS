/*
 * basetree.h
 *
 *  Created on: 08-Jan-2019
 *      Author: P. Sashittal
 */

#ifndef BASETREE_H
#define BASETREE_H

#include "utils.h"

/// This class models a tree whose nodes are labeled by unique identifiers
class BaseTree
{
public:
    /// Default constructor
    BaseTree();

    /// constructor for initialization with unsampled hosts
    BaseTree(int unhosts);
	
    /// Destructor
    virtual ~BaseTree()
    {
    }
    
    /// read Ptree
    /// @param in Input stream
		bool readPtree(std::istream& in, bool binaryTree);

    /// read Ttree
    /// @param in Input stream
		bool readTtree(std::istream& in);
    
    /// read Host file
		bool readHost(std::istream& in);
		
		/// read contact map
		bool readContactMap(std::istream& in);
		
		/// set full contact map
		void setFullContactMap();
		
		/// read infection window file
		bool readInfectionWindow(std::istream& in);
    
    /// write Ptree with no strings attached
    void writePtree(std::ostream& out) const;

    /// write Ptree with a message
    void writePtree(std::ostream& out,
                    const std::string& msg) const;
    
    /// write Ptree with given vertex labeling
    void writePtree(std::ostream& out,
                    const IntNodeMap& label) const;

    /// write Ptree with given vertex labeling and arc labeling and time
    void writePtree(std::ostream& out,
                    const IntNodeMap& vlabel,
                    const IntArcMap& elabel,
                    const DoubleArcMap& etime) const;
    
    /// write Ttree
    void writeTtree(std::ostream& out) const;
        
    /// Decide whether two nodes occur in distinct branches
    ///
    /// @param u Node
    /// @param v Node
    bool areIncomparable(Node u, Node v) const
    {
        return !isAncestor(u, v) && !isAncestor(v, u);
    }
    
    /// Decide whether one node is an ancestor (or self) of another
    ///
    /// @param u Node
    /// @param v Node
    bool isAncestor(Node u, Node v) const
    {
        while (v != _root)
        {
            if (u == v)
                return true;
            else
                v = _tree.source(InArcIt(_tree, v));
        }
        
        return u == v;
    }
    
    static bool isAncestor(const Digraph& tree,
                           Node u,
                           Node v)
    {
        while (v != lemon::INVALID)
        {
            if (u == v)
            {
                return true;
            }
            else
            {
                InArcIt a(tree, v);
                if (a == lemon::INVALID)
                {
                    v = lemon::INVALID;
                }
                else
                {
                    v = tree.source(a);
                }
            }
        }
        
        return u == v;
    }
    
    /// Return parent node. If u is the root lemon::INVALID is returned.
    ///
    /// @param u Node
    Node parent(Node u) const
    {
        if (u == _root)
        {
            return lemon::INVALID;
        }
        else
        {
            return _tree.source(InArcIt(_tree, u));
        }
    }
    
    /// Return the root node
    Node root() const
    {
        return _root;
    }
    
    /// Return underlying LEMON tree
    const Digraph& tree() const
    {
        return _tree;
    }
		
		/// Return the underlying contact map
		const Digraph& contactMap() const
		{
				return _contactMap;
		}
		
    /// Return the node name map
    const StringNodeMap& getNameMap() const
    {
        return _nodeToId;
    }
    
    /// Return the node index map
    const IntNodeMap& getIndexMap() const
    {
        return _nodeToIndex;
    }
    
    /// Return the index of a node
    const int getIndex(Node u) const
    {
        return _nodeToIndex[u];
    }
		
		/// Return the index of the contact node
		const int getContactIndex(Node u) const
		{
				return _contactNodeToIndex[u];
		}
	
		const IntNodeMap& getContactIndexMap() const
		{
				return _contactNodeToIndex;
		}
		
		const NodeVector getContactNodeVector() const
		{
				return _contactIndexToNode;
		}
		
		const Node getNode(int idx) const
		{
				assert(0 <= idx && idx < _indexToNode.size());
				return _indexToNode[idx];
		}
		
		const Node getContactNode(int idx) const
		{
				assert(0 <= idx && idx < _contactIndexToNode.size());
				return _contactIndexToNode[idx];
		}
		
    /// Return the index of an arc
    const int getArcIndex(Arc v) const
    {
        return _arcToIndex[v];
    }
    
    /// Return the timestamp of a node
    const double getTime(Node u) const
    {
        return _timestamp[u];
    }
    
    /// Return entry time of a host
    const double getEntTime(int s) const
    {
        return _enttime[s];
    }
    
    /// Return removal time of a host
    const double getRemTime(int s) const
    {
        return _remtime[s];
    }

    /// Return total time of the tree
    const double getTotalTime() const
    {
        return _totalTime;
    }

    /// Return the number of hosts
    const int getNrHost() const
    {
        return _nhosts;
    }

    /// Return number of unsampled hosts
    const int getUnHosts() const
    {
        return _unshosts;
    }
  
    /// Return maximum number of infection events
    const int getMaxInf() const
    {
        return _maxInf;
    }
    
    /// Set the number of unsamples hosts
    void setUnHosts(int unshosts)
    {
        _unshosts = unshosts;
    }

    /// set maximum number of infection events
    void setMaxInf(int maxInf)
    {
        _maxInf = maxInf;
    }
    
    /// overloading set maximum infection events
    void setMaxInf()
    {
        _maxInf = lemon::countArcs(tree());
    }
    
    /// Return the identifier of a node
    ///
    /// @param u Node
    const std::string& getName(Node u) const
    {
        return _nodeToId[u];
    }

		/// Return identifier for the contact node
		const std::string& getContactName(int s) const
		{
				return _contactNodeToId[_contactIndexToNode[s]];
		}
		
    /// Return host label of a node
    ///
    /// @param u Node
    const int getHostLabel(Node u) const
    {
        return _label[u];
    }
  
    /// Return host labeling of a node
    ///
    /// @param u Node
    const IntNodeMap& getHostLabeling() const
    {
        return _label;
    }
		
    /// Return transmission number
    int getMu(const IntNodeMap& l) const;
  
    /// Return co-transmission number
    int getGamma(const IntNodeMap& l) const;
  
    /// Return minimum co-transmission network
    ArcMatrix getN(const IntNodeMap& l) const;
  
    /// Decides whether the given set of nodes is connected
    ///
    /// @param nodes Node set
    bool isConnected(const NodeSet& nodes) const;
    
    /// Return a node corresponding to the given identifier.
    /// If there is no node with the given identifier, lemon::INVALID is returned.
    ///
    /// @param lbl Identifier
    Node getNodeByLabel(const std::string& lbl) const
    {
        auto it = _idToNode.find(lbl);
        if (it == _idToNode.end())
        {
            return lemon::INVALID;
        }
        else
        {
            return it->second;
        }
    }
  
    const std::string& getHostLabel(int s) const
    {
      assert(0 <= s && s < _nhosts);
      return _hostLabel[s];
    }
  
		Node getContactNode(const std::string& name) const
		{
				auto it = _contactIdToNode.find(name);
				if (it == _contactIdToNode.end())
				{
						return lemon::INVALID;
				}
				else
				{
						return it->second;
				}
		}
		
		int getContactNumStrains(const Arc est) const
		{
				return _contactArcToNumStrains[est];
		}
		
		bool infectionWindowExists(const int s) const
		{
				if (s < 0 || s >= _nhosts)
				{
						return false;
				}
				else
				{
						std::string hostName = _hostLabel[s];
						
						if (_hostInfectionWindow.count(hostName) != 0)
						{
								return true;
						}
						else
						{
								return false;
						}
				}
		}
		
		double getInfectionWindow (const int s) const
		{
				assert(0 <= s && s < _nhosts);
				assert(infectionWindowExists(s));
				
				std::string hostName = _hostLabel[s];
				return _hostInfectionWindow.at(hostName);
		}
		
    void writeDOT(const ArcMatrix& N,
                  const IntNodeMap& l,
                  std::ostream& out) const;
  
    static void readInputWithHosts(std::istream& hfile,
                                   std::istream& pfile,
                                   int unsampled_hosts,
                                   bool binaryTree,
                                   BaseTree& B);
  
		static void readInputWithHostsContactMap(std::istream& hfile,
																						 std::istream& pfile,
																						 std::istream& cfile,
																						 int unsampled_hosts,
																						 bool binaryTree,
																						 BaseTree& B);
  
		static void readInputWithHostsTransmissionTree(std::istream& hfile,
																									 std::istream& pfile,
																									 std::istream& cfile,
																									 int unsampled_hosts,
																									 bool binaryTree,
																									 BaseTree& B);
		
protected:
    /// Tree
    Digraph _tree;
    /// Root
    Node _root;
    /// Node label (name)
    StringNodeMap _nodeToId;
    /// Get node by label (name)
    StringToNodeMap _idToNode;
    /// Node index
    IntNodeMap _nodeToIndex;
		/// Reverse node index
		NodeVector _indexToNode;
    /// Arc index
    IntArcMap _arcToIndex;
    /// timestamp of each node
    DoubleNodeMap _timestamp;
    /// vertex label (host) of nodes
    IntNodeMap _label;
    /// host label (string)
    StringVector _hostLabel;
    /// number of infected hosts
    int _nhosts;
    /// number of unsampled hosts
    int _unshosts;
    /// entry time for infected hosts
    DoubleVector _enttime;
    /// removal time of infected hosts
    DoubleVector _remtime;
    /// maximum number of infection events
    int _maxInf;
    /// total time of the tree
    double _totalTime;
		/// contact map
		Digraph _contactMap;
		/// contact map node id (name)
		StringNodeMap _contactNodeToId;
		/// reverse node id (name)
		StringToNodeMap _contactIdToNode;
		/// contact map node index
		IntNodeMap _contactNodeToIndex;
		/// contact reverse map index to node
		NodeVector _contactIndexToNode;
		/// number of strains transferred
		IntArcMap _contactArcToNumStrains;
		
		/// infection window data
		StringToDoubleMap _hostInfectionWindow;
};

#endif // BASETREE_H
