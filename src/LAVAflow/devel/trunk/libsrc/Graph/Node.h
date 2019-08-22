#ifndef LAVA_GRAPH_NODE_H
#define LAVA_GRAPH_NODE_H

#include <list>
#include <iostream>


namespace Graph
{

template<class N, class E>
class Edge;


template<class N, class E>
class Node
{

private:

	N value;
	std::list<Edge<N,E>*> incidentEdges;

	// Auxilliary variables for various graph algorithms
	bool wasVisited;

public:
	Node()
	{
		std::clog<<"[Graph::Node] Creating empty node"<<std::endl;
	}
	Node(N valueIn)
	{
		this->value = valueIn;
	}



	N getValue(){ return this->value; }
	void setValue(N value){ this->value = value; }
	
	void setVisitedState(bool visit){ this->wasVisited = visit; }
	bool getVisitedState(){ return this->wasVisited; }
	void visitNeighbors();

	int getNumberOfEdges(){return this->incidentEdges.size();}
	bool isDisconnected(){return this->incidentEdges.empty();}

	void addEdge(Edge<N,E>* edge);
	void removeEdge(Edge<N,E>* edge);

	// Iterators over incident edges
	typename std::list<Edge<N,E>*>::iterator begin(){ return this->incidentEdges.begin();}
	typename std::list<Edge<N,E>*>::iterator end(){ return this->incidentEdges.end();}
	typename std::list<Edge<N,E>*>::reverse_iterator rbegin(){ return this->incidentEdges.rbegin();}
	typename std::list<Edge<N,E>*>::reverse_iterator rend(){ return this->incidentEdges.rend();}

	// List accessors
	Edge<N,E>* back(){ return this->incidentEdges.back(); }
	void pop_back(){ this->incidentEdges.pop_back(); }

	/* data */
};

/*==========================================

	Add an incident edge to the node

==========================================*/
template<class N, class E>
void Node<N,E>::addEdge(Edge<N,E>* edge)
{
	this->incidentEdges.push_back(edge);
}


/*==========================================

	Remove an incident edge from the node

==========================================*/
template<class N, class E>
void Node<N,E>::removeEdge(Edge<N,E>* edge)
{
	this->incidentEdges.remove(edge);
}


/*==========================================

	Visit all neighbors and set their 
	wasVisited flag

==========================================*/
template<class N, class E>
void Node<N,E>::visitNeighbors()
{
	// Visited the current node
	this->setVisitedState(true);

	typename std::list<Edge<N,E>*>::iterator edgeIt;
	for(edgeIt = this->incidentEdges.begin(); edgeIt != this->incidentEdges.end(); ++edgeIt)
	{
		// Get the opposite node of the edge
		Node<N,E>* W = (*edgeIt)->getOtherNode(this);

		// If the other node is already flooded, continue
		if(W->getVisitedState())
		{
			continue;
		}

		// Otherwise, visit all of W's neighbors
		W->visitNeighbors();

	}
}


}; // end namespace
#endif