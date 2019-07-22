#ifndef LAVA_GRAPH_EDGE_H
#define LAVA_GRAPH_EDGE_H

#include <cstdlib>

namespace Graph
{


template<class N, class E>
class Node;

template<class N, class E>
class Edge
{

private:
	Node<N,E>* first;
	Node<N,E>* second;

	E value;

public:
	Edge()
	{
		first = NULL;
		second = NULL;
		value = NULL;
	}


	Edge(Node<N,E>* firstNode, Node<N,E>* secondNode, E edgeValue)
	{
		first = firstNode;
		second = secondNode;
		value = edgeValue;
	}

	/* data */

	// Accessors
	E getValue(){ return this->value; }
	void setValue(E valueIn){ this->value = valueIn; }
	Node<N,E>* getFirst(){ return this->first; }
	Node<N,E>* getSecond(){ return this->second; }
	void setFirst(Node<N,E>* firstIn){ this->first = firstIn; }
	void setSecond(Node<N,E>* secondIn){ this->second = secondIn; }


	// Utilities
	Node<N,E>* getOtherNode(Node<N,E>* nodeStart);


};



/*==========================================

	Given a node, return the other node
	which comprises the edge

==========================================*/
template<class N, class E>
Node<N,E>* Edge<N,E>::getOtherNode(Node<N,E>* nodeStart)
{
	Node<N,E>* otherNode; 
	otherNode = NULL;

	if(nodeStart==this->first)
		otherNode = this->second;
	else if(nodeStart==this->second)
		otherNode = this->first;

	return otherNode;
}





}; // end namespace

#endif