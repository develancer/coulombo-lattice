// Coulombo Ⓒ 2018
// [Computer Physics Communications] Różański & Zieliński:
// Efficient computation of Coulomb and exchange integrals for multi-million atom nanostructures
#ifndef COULOMBO_PLANNER_HPP
#define COULOMBO_PLANNER_HPP

#include <algorithm>
#include <vector>
#include "Graph.hpp"

//----------------------------------------------------------------------

struct PlannerProduct {
	int product;
	bool conjugate;
};

struct PlannerIntegral {
	PlannerProduct left, right;
};

struct PlannerStep {
	int id; // ID of integral to be computed in this step
	PlannerProduct left, right;
};

//----------------------------------------------------------------------

class Planner {
public:
	Planner(int productCount)
			:graph(productCount) { }

	void addIntegral(int i1, int i2, int i3, int i4)
	{
		PlannerIntegral integral{getProduct(i1, i4), getProduct(i2, i3)};
		graph.addEdge(integral.left.product, integral.right.product);
		integrals.push_back(integral);
	}

	std::vector<PlannerStep> computePlan(void) const
	{
		auto cover = graph.computeVertexCover();
		std::vector<PlannerStep> steps;
		int id = 0;
		for (PlannerIntegral integral : integrals) {
			if (!cover.count(integral.left.product)) {
				std::swap(integral.left, integral.right);
			}
			steps.push_back({id++, integral.left, integral.right});
		}
		std::sort(steps.begin(), steps.end(), [](const PlannerStep& a, const PlannerStep& b) {
			if (a.left.product<b.left.product) return true;
			if (a.left.product>b.left.product) return false;
			if (a.right.product<b.right.product) return true;
			if (a.right.product>b.right.product) return false;
			return (a.left.conjugate==a.right.conjugate)
					&& (b.left.conjugate!=b.right.conjugate);
		});
		return steps;
	}

private:
	Graph graph;
	std::vector<PlannerIntegral> integrals;

	PlannerProduct getProduct(int iL, int iR)
	{
		bool conjugate = false;
		if (iL<iR) {
			conjugate = true;
			std::swap(iL, iR);
		}
		int iLR = iL*(iL-1)/2+iR-1;
		assert(iLR>=0);
		return PlannerProduct{iLR, conjugate};
	}
};

//----------------------------------------------------------------------

#endif
