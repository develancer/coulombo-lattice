// Coulombo Ⓒ 2018
// [Computer Physics Communications] Różański & Zieliński:
// Efficient computation of Coulomb and exchange integrals for multi-million atom nanostructures
#include <cppunit/TestCase.h>
#include <cppunit/TestCaller.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/ui/text/TestRunner.h>

//----------------------------------------------------------------------

#include "CoulombCalculator.hpp"

class CoulombCalculatorTest_Interaction: public Interaction {
	void map(Domain<real>& G) const
	{
		for (int iz = 0; iz<G.dimension.z; ++iz)
			for (int iy = 0; iy<G.dimension.y; ++iy)
				for (int ix = 0; ix<G.dimension.x; ++ix)
					G(ix, iy, iz) = exp(-ix-iy-iz);
	}
};

class CoulombCalculatorTest: public CppUnit::TestCase {
	CPPUNIT_TEST_SUITE( CoulombCalculatorTest );
	CPPUNIT_TEST( testTwoPointCharges );
	CPPUNIT_TEST_SUITE_END();

	void testTwoPointCharges()
	{
		Dimension dimension(5, 5, 10);
		CoulombCalculator cc(dimension);
		complex result;

		cc.initialize(CoulombCalculatorTest_Interaction());

		cc.input.zeros();
		cc.input(2, 2, 2) = 1.0;
		cc.prepare();
		result = cc.calculate();
		CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, result.real(), 1.0e-12);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, result.imag(), 1.0e-12);

		cc.input.zeros();
		cc.input(2, 2, 7) = 1.0;
		result = cc.calculate();
		CPPUNIT_ASSERT_DOUBLES_EQUAL(exp(-5), result.real(), 1.0e-12);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, result.imag(), 1.0e-12);
	}
};

//----------------------------------------------------------------------

#include "Dimension.hpp"

class DimensionTest: public CppUnit::TestCase {
	CPPUNIT_TEST_SUITE( DimensionTest );
	CPPUNIT_TEST( testSimple );
	CPPUNIT_TEST_EXCEPTION( testTooLarge, std::logic_error );
	CPPUNIT_TEST_EXCEPTION( testNegative, std::logic_error );
	CPPUNIT_TEST_SUITE_END();

	void testSimple()
	{
		Dimension dimension(20, 30, 40);
		uword cells = 20*30*40;
		CPPUNIT_ASSERT_EQUAL(20, dimension.x);
		CPPUNIT_ASSERT_EQUAL(30, dimension.y);
		CPPUNIT_ASSERT_EQUAL(40, dimension.z);
		CPPUNIT_ASSERT_EQUAL(cells, dimension.cells());

		dimension = dimension.twice();
		cells *= 8;
		CPPUNIT_ASSERT_EQUAL(40, dimension.x);
		CPPUNIT_ASSERT_EQUAL(60, dimension.y);
		CPPUNIT_ASSERT_EQUAL(80, dimension.z);
		CPPUNIT_ASSERT_EQUAL(cells, dimension.cells());
	}

	void testTooLarge()
	{
		int large = 1000000000;
		Dimension(large, large, large);
	}

	void testNegative()
	{
		Dimension(30, 10, -10);
	}
};

//----------------------------------------------------------------------

#include "Domain.hpp"
#include "DomainAllocator.hpp"

class DomainAllocatorTest: public CppUnit::TestCase {
	CPPUNIT_TEST_SUITE( DomainAllocatorTest );
	CPPUNIT_TEST( testAllocate );
	CPPUNIT_TEST_EXCEPTION( testTooLarge, std::bad_array_new_length );
	CPPUNIT_TEST( testWithVector );
	CPPUNIT_TEST_SUITE_END();

public:
	void testAllocate()
	{
		const size_t size = 1000;
		const DomainAllocator<unsigned> allocator;
		unsigned* p = allocator.allocate(size);
		CPPUNIT_ASSERT(p);
		for (unsigned i = 0; i<size; ++i) {
			p[i] = i;
		}
		for (unsigned i = 0; i<size; ++i) {
			CPPUNIT_ASSERT_EQUAL(i, p[i]);
		}
		allocator.deallocate(p, 1000);
	}

	void testTooLarge()
	{
		const size_t size = std::numeric_limits<size_t>::max()/2;
		const DomainAllocator<long> allocator;
		allocator.allocate(size);
	}

	void testWithVector()
	{
		const size_t size = 1000;
		DomainData<unsigned> p;
		p.resize(size/2);
		CPPUNIT_ASSERT(p.size()==size/2);
		p.resize(size*2);
		CPPUNIT_ASSERT(p.size()==size*2);
		p.resize(size);
		CPPUNIT_ASSERT(p.size()==size);
		for (unsigned i = 0; i<size; ++i) {
			p[i] = i;
		}
		for (unsigned i = 0; i<size; ++i) {
			CPPUNIT_ASSERT_EQUAL(i, p[i]);
		}
	}
};

//----------------------------------------------------------------------

#include "Domain.hpp"

class DomainTest: public CppUnit::TestCase {
	CPPUNIT_TEST_SUITE( DomainTest );
	CPPUNIT_TEST( testSimple );
	CPPUNIT_TEST_SUITE_END();

	void testSimple()
	{
		const int X = 2, Y = 3, Z = 4;
		const uword cells = X*Y*Z;
		DomainData<real> data(cells);
		Domain<real> domain(data.data(), Dimension(X, Y, Z));
		for (uword i = 0; i<cells; ++i) {
			data[i] = i;
		}
		CPPUNIT_ASSERT_EQUAL(X, static_cast<int>(domain.n_rows));
		CPPUNIT_ASSERT_EQUAL(Y, static_cast<int>(domain.n_cols));
		CPPUNIT_ASSERT_EQUAL(Z, static_cast<int>(domain.n_slices));
		CPPUNIT_ASSERT_EQUAL(X, domain.dimension.x);
		CPPUNIT_ASSERT_EQUAL(Y, domain.dimension.y);
		CPPUNIT_ASSERT_EQUAL(Z, domain.dimension.z);
		// testing column-major ordering
		for (int ix = 0; ix<X; ++ix)
			for (int iy = 0; iy<Y; ++iy)
				for (int iz = 0; iz<Z; ++iz) {
					real i = (iz*Y+iy)*X+ix;
					CPPUNIT_ASSERT_EQUAL(i, domain(ix, iy, iz));
				}
	}
};

//----------------------------------------------------------------------

#include "Domain.hpp"
#include "FourierPlan.hpp"

class FourierPlanTest: public CppUnit::TestCase {
	CPPUNIT_TEST_SUITE( FourierPlanTest );
	CPPUNIT_TEST( testSimpleTransformR2C );
	CPPUNIT_TEST_SUITE_END();

	void testSimpleTransformR2C()
	{
		const int inputSize = 1024;
		const int outputSize = inputSize/2+1;
		DomainData<real> input(inputSize);
		DomainData<complex> output(outputSize);
		FourierPlan plan(fftw_plan_dft_r2c_1d(
				inputSize,
				input.data(),
				reinterpret_cast<fftw_complex*>(output.data()),
				FFTW_ESTIMATE | FFTW_PRESERVE_INPUT
		));
		int k0 = 16;
		double phase0 = M_PI_4;
		for (int n = 0; n<inputSize; ++n) {
			input[n] = cos(2*M_PI*k0*n/inputSize+phase0);
		}
		plan.execute();
		for (int k = 0; k<outputSize; ++k) {
			if (k!=k0) {
				CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, std::abs(output[k]), 1.0e-12);
			}
		}
		CPPUNIT_ASSERT_DOUBLES_EQUAL(inputSize/2, std::abs(output[k0]), 1.0e-12);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(phase0, std::arg(output[k0]), 1.0e-12);
	}
};

//----------------------------------------------------------------------

#include "Graph.hpp"

class GraphTest: public CppUnit::TestCase {
	CPPUNIT_TEST_SUITE( GraphTest );
	CPPUNIT_TEST( testVertexCoverClique );
	CPPUNIT_TEST( testVertexCoverNoEdges );
	CPPUNIT_TEST( testVertexCoverLoops );
	CPPUNIT_TEST( testVertexCoverTwoEdges );
	CPPUNIT_TEST( testVertexCoverTriangleWithLoops );
	CPPUNIT_TEST_SUITE_END();

	void testVertexCoverClique()
	{
		Graph graph(4);
		for (int j = 0; j<4; ++j) {
			for (int i = 0; i<j; ++i) {
				graph.addEdge(i, j);
			}
		}
		auto cover = graph.computeVertexCover();
		CPPUNIT_ASSERT(cover.size()==3);
	}

	void testVertexCoverNoEdges()
	{
		Graph graph(3);
		auto cover = graph.computeVertexCover();
		CPPUNIT_ASSERT(cover.empty());
	}

	void testVertexCoverLoops()
	{
		Graph graph(3);
		graph.addEdge(0, 0);
		graph.addEdge(1, 1);
		graph.addEdge(2, 2);
		auto cover = graph.computeVertexCover();
		CPPUNIT_ASSERT(cover.size()==3);
	}

	void testVertexCoverTwoEdges()
	{
		Graph graph(3);
		graph.addEdge(0, 1);
		graph.addEdge(1, 2);
		auto cover = graph.computeVertexCover();
		CPPUNIT_ASSERT(cover.size()==1);
		CPPUNIT_ASSERT(cover.count(1));
	}

	void testVertexCoverTriangleWithLoops()
	{
		Graph graph(3);
		graph.addEdge(0, 0);
		graph.addEdge(0, 1);
		graph.addEdge(1, 1);
		graph.addEdge(1, 2);
		graph.addEdge(2, 0);
		auto cover = graph.computeVertexCover();
		CPPUNIT_ASSERT(cover.size()==2);
		CPPUNIT_ASSERT(cover.count(0));
		CPPUNIT_ASSERT(cover.count(1));
	}
};

//----------------------------------------------------------------------

#include "Round.hpp"

class RoundTest: public CppUnit::TestCase {
	CPPUNIT_TEST_SUITE( RoundTest );
	CPPUNIT_TEST( testPowerOf2 );
	CPPUNIT_TEST( testPowerOf235 );
	CPPUNIT_TEST( testPowerOf235div4 );
	CPPUNIT_TEST_SUITE_END();

	void testPowerOf2()
	{
		CPPUNIT_ASSERT_EQUAL(1, Round::up<2>(0, 1));
		CPPUNIT_ASSERT_EQUAL(16, Round::up<2>(13, 1));
		CPPUNIT_ASSERT_EQUAL(16, Round::up<2>(16, 1));
		CPPUNIT_ASSERT_EQUAL(32, Round::up<2>(17, 1));
		CPPUNIT_ASSERT_EQUAL(64, Round::up<2>(34, 1));
		CPPUNIT_ASSERT_EQUAL(128, Round::up<2>(101, 1));
	}

	void testPowerOf235()
	{
		CPPUNIT_ASSERT_EQUAL(1, (Round::up<2, 3, 5>(0, 1)));
		CPPUNIT_ASSERT_EQUAL(15, (Round::up<2, 3, 5>(13, 1)));
		CPPUNIT_ASSERT_EQUAL(16, (Round::up<2, 3, 5>(16, 1)));
		CPPUNIT_ASSERT_EQUAL(18, (Round::up<2, 3, 5>(17, 1)));
		CPPUNIT_ASSERT_EQUAL(36, (Round::up<2, 3, 5>(34, 1)));
		CPPUNIT_ASSERT_EQUAL(108, (Round::up<2, 3, 5>(101, 1)));
	}

	void testPowerOf235div4()
	{
		CPPUNIT_ASSERT_EQUAL(4, (Round::up<2, 3, 5>(0, 4)));
		CPPUNIT_ASSERT_EQUAL(16, (Round::up<2, 3, 5>(13, 4)));
		CPPUNIT_ASSERT_EQUAL(16, (Round::up<2, 3, 5>(16, 4)));
		CPPUNIT_ASSERT_EQUAL(20, (Round::up<2, 3, 5>(17, 4)));
		CPPUNIT_ASSERT_EQUAL(36, (Round::up<2, 3, 5>(34, 4)));
		CPPUNIT_ASSERT_EQUAL(108, (Round::up<2, 3, 5>(101, 4)));
	}
};

//----------------------------------------------------------------------

#include "mpi.hpp"

int main(int argc, char** argv)
{
	mpi::init(argc, argv);
	if (mpi::size()!=1) return 1;

	CppUnit::TextUi::TestRunner runner;
	runner.addTest(CoulombCalculatorTest::suite());
	runner.addTest(DimensionTest::suite());
	runner.addTest(DomainAllocatorTest::suite());
	runner.addTest(DomainTest::suite());
	runner.addTest(FourierPlanTest::suite());
	runner.addTest(GraphTest::suite());
	runner.addTest(RoundTest::suite());
	runner.run();

	mpi::finalize();
	return 0;
}
