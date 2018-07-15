// Coulombo Ⓒ 2018
// [Computer Physics Communications] Różański & Zieliński:
// Efficient computation of Coulomb and exchange integrals for very large systems
#ifndef COULOMBO_GRAPH_HPP
#define COULOMBO_GRAPH_HPP

#include <cassert>
#include <set>
#include <utility>
#include <vector>

//----------------------------------------------------------------------

using VertexIndex = int;
using Edge = std::pair<VertexIndex, VertexIndex>;
using AdjacencyList = std::vector<std::set<VertexIndex>>;

//----------------------------------------------------------------------

class GraphVertexComparator {
public:
	GraphVertexComparator(AdjacencyList& neighbors)
			:neighbors(neighbors) { }

	bool operator()(VertexIndex v, VertexIndex w) const
	{
		auto vs = neighbors[v].size(), ws = neighbors[w].size();
		return (vs<ws) || (vs==ws && v<w);
	}

private:
	AdjacencyList& neighbors;
};

//----------------------------------------------------------------------

class Graph {
public:
	Graph(VertexIndex vertexCount)
			:vertexCount(vertexCount), edges(edges_) { }

	void addEdge(VertexIndex start, VertexIndex end)
	{
		assert(start>=0 && start<vertexCount);
		assert(end>=0 && end<vertexCount);
		edges_.insert({start, end});
	}

	std::set<VertexIndex> computeVertexCover() const
	{
		std::set<VertexIndex> cover;
		for (const Edge& edge : edges_) {
			if (edge.second==edge.first) {
				cover.insert(edge.first);
			}
		}
		AdjacencyList neighbors(vertexCount);
		for (const Edge& edge : edges_) {
			if (edge.second!=edge.first && !cover.count(edge.first) && !cover.count(edge.second)) {
				neighbors[edge.first].insert(edge.second);
				neighbors[edge.second].insert(edge.first);
			}
		}
		GraphVertexComparator comparator(neighbors);
		std::set<VertexIndex, GraphVertexComparator> queue(comparator);
		for (VertexIndex v = 0; v<vertexCount; ++v) {
			if (!neighbors[v].empty()) {
				queue.insert(v);
			}
		}
		while (!queue.empty()) {
			VertexIndex v = *queue.begin();
			if (neighbors[v].size()==1) {
				v = *neighbors[v].begin();
			}
			else {
				v = *queue.rbegin();
			}
			queue.erase(v);
			cover.insert(v);
			for (VertexIndex w : neighbors[v]) {
				if (queue.erase(w)) {
					neighbors[w].erase(v);
					if (!neighbors[w].empty()) {
						queue.insert(w);
					}
				}
			}
		}
		return cover;
	}

private:
	VertexIndex vertexCount;
	std::set<Edge> edges_;

public:
	const std::set<Edge>& edges;
};

//----------------------------------------------------------------------

#endif
