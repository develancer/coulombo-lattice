// Coulombo Ⓒ 2018
// [Computer Physics Communications] Różański & Zieliński:
// Efficient computation of Coulomb and exchange integrals for multi-million atom nanostructures
#ifndef COULOMBO_PLANNER_HPP
#define COULOMBO_PLANNER_HPP

#include <algorithm>
#include <vector>
#include "mpi.hpp"
#include "Graph.hpp"

//----------------------------------------------------------------------

struct PlannerProduct {
	int product;
	unsigned char conjugate; // instead of bool since there is no MPI_BOOL
};

struct PlannerStep {
	int id; // ID of integral to be computed in this step
	PlannerProduct left, right;
};

//----------------------------------------------------------------------

class Planner {
public:
	Planner() : typePlannerProduct(createTypePlannerProduct()) { }

	virtual ~Planner() { }

	virtual void addIntegral(short, short, short, short)
	{
		// nothing here
	}

	virtual void computePlan()
	{
		// nothing here
	}

	virtual bool getNextStep(PlannerStep& step)
	{
		mpi::broadcast(&step.id, 1, MPI_INT);
		mpi::broadcast(&step.left, 1, typePlannerProduct);
		mpi::broadcast(&step.right, 1, typePlannerProduct);
		return step.id >= 0;
	}

protected:
	mpi::Type typePlannerProduct;

private:
	mpi::Type createTypePlannerProduct()
	{
		int dimensionBlockLengths[2] = {1, 1};
		MPI_Aint dimensionDisplacements[2] = {offsetof(PlannerProduct, product), offsetof(PlannerProduct, conjugate)};
		MPI_Datatype dimensionTypes[2] = {MPI_INT, MPI_UNSIGNED_CHAR};
		return mpi::structure(2, dimensionBlockLengths, dimensionDisplacements, dimensionTypes);
	}
};

//----------------------------------------------------------------------

class MasterPlanner : public Planner {
public:
	MasterPlanner(int productCount)
			:productCount(productCount), finalized(false)
	{ }

	void addIntegral(short i1, short i2, short i3, short i4)
	{
		if (finalized) {
			throwfl("plan is already computed");
		}
		auto id = steps.size();
		if (id > std::numeric_limits<int>::max()) {
			throwfr("too many integrals to compute");
		}
		steps.push_back({static_cast<int>(id), getProduct(i1, i4), getProduct(i2, i3)});
	}

	void computePlan()
	{
		if (finalized) {
			throwfl("plan is already computed");
		}
		auto cover = computeVertexCover();
		for (PlannerStep& step : steps) {
			if (!cover.count(step.left.product)) {
				std::swap(step.left, step.right);
			}
		}
		std::sort(steps.begin(), steps.end(), [](const PlannerStep& a, const PlannerStep& b) {
			if (a.left.product<b.left.product) return true;
			if (a.left.product>b.left.product) return false;
			if (a.right.product<b.right.product) return true;
			if (a.right.product>b.right.product) return false;
			return (a.left.conjugate==a.right.conjugate)
					&& (b.left.conjugate!=b.right.conjugate);
		});
		nextStep = steps.begin();
		finalized = true;
	}

	bool getNextStep(PlannerStep& step)
	{
		if (!finalized) {
			throwfl("plan is not yet computed");
		}
		if (nextStep == steps.end()) {
			step = {-1, {-1, false}, {-1, false}};
		} else {
			step = *nextStep++;
		}
		return Planner::getNextStep(step);
	}

private:
	const int productCount;
	std::vector<PlannerStep> steps;
	std::vector<PlannerStep>::iterator nextStep;
	bool finalized;

	std::set<VertexIndex> computeVertexCover() const
	{
		Graph graph(productCount);
		for (const PlannerStep& step : steps) {
			graph.addEdge(step.left.product, step.right.product);
		}
		return graph.computeVertexCover();
	}

	PlannerProduct getProduct(short iL, short iR) const
	{
		unsigned char conjugate = 0;
		if (iL<iR) {
			conjugate = 1;
			std::swap(iL, iR);
		}
		int iLR = static_cast<int>(iL)*static_cast<int>(iL-1)/2+iR-1;
		assert(iLR>=0);
		return PlannerProduct{iLR, conjugate};
	}
};

//----------------------------------------------------------------------

#endif
