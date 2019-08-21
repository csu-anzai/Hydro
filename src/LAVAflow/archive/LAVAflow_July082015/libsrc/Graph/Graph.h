#ifndef LAVA_GRAPH_GRAPH_H
#define LAVA_GRAPH_GRAPH_H

#include <list>
#include <iostream>
#include <ostream>

#include "Node.h"
#include "Edge.h"


namespace Graph
{

template<class N, class E>
class Graph
{

private:
	// List of nodes in the graph holding type N
	std::list<Node<N,E>*> nodeList;
	// List of edges in the graph holding type E
	std::list<Edge<N,E>*> edgeList;
	// List of subgraph "root" Nodes
	std::list<Node<N,E>*> subgraphRoots;


public:
	Graph()
	{
		std::clog<<"[Graph::Graph] Creating generic graph"<<std::endl;
	}
	~Graph()
	{
		// TODO remove all nodes and edges properly to avoid memory leaks
	}


	int getNumberOfNodes(){ return this->nodeList.size(); }
	int getNumberOfEdges(){ return this->edgeList.size(); }
	int getNumberOfSubgraphs(){ return this->subgraphRoots.size(); }

	Node<N,E>* addNode(N nodeVal);
	Edge<N,E>* addEdge(Node<N,E>* first, Node<N,E>* second, E value);
	bool removeNode(Node<N,E>* nodeToRemove);
	bool removeEdge(Edge<N,E>* edgeToRemove);

	Node<N,E>* findNode(N nodeVal);

	typename std::list<Node<N,E>*>::iterator nodeBegin(){ return this->nodeList.begin();}
	typename std::list<Node<N,E>*>::iterator nodeEnd(){ return this->nodeList.end();}


	void computeSubgraphRoots();

};



/*==========================================

	Add a node to the graph.
	This function returns a pointer to the
	new node

==========================================*/
template<class N, class E>
Node<N,E>* Graph<N,E>::addNode(N nodeVal)
{
	Node<N,E>* newNode = new Node<N,E>(nodeVal);

	// Attempt to push the new node onto the nodeList
	this->nodeList.push_back(newNode);

	// Return the newNode
	return newNode;
}

/*==========================================

	Add an edge to the graph.
	This function returns a pointer to the
	new edge

==========================================*/
template<class N, class E>
Edge<N,E>* Graph<N,E>::addEdge(Node<N,E>* first, Node<N,E>* second, E value)
{
	Edge<N,E>* newEdge = new Edge<N,E>(first, second, value);

	// Add the new edge to the incident edge lists for the two nodes
	first->addEdge(newEdge);
	second->addEdge(newEdge);

	// Attempt to push the new edge onto the edgeList
	edgeList.push_back(newEdge);

	// Return the newEdge
	return newEdge;
}




/*==========================================

	Remove a node from the graph. This involves 
	also removing all incident edges from the graph

	(1)	Get all edges incident to the node
	(2)	For every incident edge, remove the edge from the OTHER node's incident list
	(3)	Delete all incident edges
	(4)	Delete the node

==========================================*/
template<class N, class E>
bool Graph<N,E>::removeNode(Node<N,E>* nodeToRemove)
{
	// Remove the edges from the opposing nodes' edge lists
	// as well as the Graph's edgeList
	typename std::list<Edge<N,E>*>::iterator edgeIt;
	for(edgeIt = nodeToRemove->begin(); edgeIt != nodeToRemove->end(); ++edgeIt)
	{
		// Get the opposite node of the edge
		Node<N,E>* otherNode = (*edgeIt)->getOtherNode(nodeToRemove);

		// Remove the edge from the otherNode's incident edge list
		otherNode->removeEdge(*edgeIt);

		// Remove from the Graph's edge list
		this->edgeList.remove(*edgeIt);
	}

	// Delete all incident edges
	while(!nodeToRemove->isDisconnected())
	{
		delete nodeToRemove->back();
		nodeToRemove->pop_back();
	}

	// Remove the node from the Graph's nodeList
	this->nodeList.remove(nodeToRemove);

	// Delete the node itself
	delete nodeToRemove;

}


/*==========================================

	Remove an edge from the graph. This involves
	removing it from the incident list of the 
	two nodes it connects

==========================================*/
template<class N, class E>
bool Graph<N,E>::removeEdge(Edge<N,E>* edgeToRemove)
{
	// Get the two affected nodes
	Node<N,E>* first = edgeToRemove->getFirst();
	Node<N,E>* second = edgeToRemove->getSecond();

	// Remove the edge from their incident lists
	first->removeEdge(edgeToRemove);
	second->removeEdge(edgeToRemove);

	// Remove the edge from the Graph's edge list
	this->edgeList.remove(edgeToRemove);

	// Delete the edge
	delete edgeToRemove;
}



/*==========================================

	Find the node in the graph that 
	has its value equal to the provided 
	search value

	This method will find ONE node in the graph.
	For graphs with unique nodes identifiers, this
	is fine. For other graphs, one needs to implement
	a new function to find ALL the nodes satisfying
	this criterion.

	This method traverses the list of nodes
	from the front and back simultaneously
	in an effort to get some speedup for large 
	graphs.

==========================================*/
template<class N, class E>
Node<N,E>* Graph<N,E>::findNode(N nodeVal)
{

	Node<N,E>* foundNode = NULL;

	typename std::list<Node<N,E>*>::iterator fIt;
	typename std::list<Node<N,E>*>::reverse_iterator rIt;

	// Traverse the list of nodes
	for(fIt = nodeList.begin(), rIt = nodeList.rbegin(); 
		fIt != nodeList.end(), rIt != nodeList.rend(); 
		++fIt, ++rIt)
	{
		if((*fIt)->getValue()==nodeVal)
		{
			foundNode = *fIt;
			break;
		}

		if((*rIt)->getValue()==nodeVal)
		{
			foundNode = *rIt;
			break;
		}
	}

	return foundNode;
}




/*==========================================

	Find the "root" nodes of all subgraphs.
	The word "root" implies the first node
	visited of a subgraph, such that traversal
	beginning at the "root" can find all nodes
	in the subgraph. This is not the same concept
	of "root" as used in trees. 

	There are two options when doing this:
	One is to actually store the subgraphs as
	Graph objects. This itself has two options.
	The first is to create new Node and Edge objects
	for each subgraph. The drawback here is memory usage.
	The other option is to use the pointers to the previously
	made objects (from the supergraph). This has the disadvantage
	that when the destructor is called on the subgraphs, it will
	deallocate the Nodes and Edges of the supergraph.

	The second option is to only store the individual "root" nodes
	for each subgraph in the supergraph. This Node object can then 
	be passed to other utility functions (e.g. search functions), 
	as the basis of traversal. This is the approach we use here.
	
==========================================*/

template<class N, class E>
void Graph<N,E>::computeSubgraphRoots()
{
	// Clear subgraph roots
	this->subgraphRoots.clear();

	// Set all flooded states to false
	typename std::list<Node<N,E>*>::iterator it;
	for(it=this->nodeList.begin(); it!=this->nodeList.end(); ++it)
	{
		(*it)->setVisitedState(false);
	}

	// Loop over all Nodes in the graph
	for(it=this->nodeList.begin(); it!=this->nodeList.end(); ++it)
	{
		// If the current Node has been visited, continue
		if((*it)->getVisitedState())
			continue;

		// Otherwise, add the current node to the Root list and
		// visit it's neighbors
		// The Node::visitNeighbors function is recursive and will
		// visit all Nodes in the subgraph
		subgraphRoots.push_back(*it);
		(*it)->visitNeighbors();
	}

}



}; // end namespace




#endif