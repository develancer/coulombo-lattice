// Coulombo Ⓒ 2018
// [Computer Physics Communications] Różański & Zieliński:
// Efficient computation of Coulomb and exchange integrals for multi-million atom nanostructures
#include "mpi.hpp"
#include "throwf.hpp"

#define THROW_IF_FAILS(COMMAND) if (COMMAND) { throwfr("mpi::%s failed", __FUNCTION__); }

namespace mpi {
//----------------------------------------------------------------------

Type::Type(MPI_Datatype datatype)
		:datatype(datatype) { }

Type::Type(Type&& source)
		:datatype(source.datatype)
{
	source.datatype = nullptr;
}

Type::~Type()
{
	if (datatype) {
		MPI_Type_free(&datatype);
		datatype = nullptr;
	}
}

void Type::operator=(Type&& source)
{
	if (datatype) {
		MPI_Type_free(&datatype);
		datatype = nullptr;
	}
	std::swap(datatype, source.datatype);
}

Type::operator MPI_Datatype() const
{
	return datatype;
}

//----------------------------------------------------------------------

void init(int& argc, char**& argv)
{
	int provided = 0;
	THROW_IF_FAILS(
			MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided)
	);
}

void finalize(void)
{
	MPI_Finalize();
}

int rank(void)
{
	int rank;
	THROW_IF_FAILS(MPI_Comm_rank(MPI_COMM_WORLD, &rank));
	return rank;
}

bool root(void)
{
	return !rank();
}

int size(void)
{
	int size;
	THROW_IF_FAILS(MPI_Comm_size(MPI_COMM_WORLD, &size));
	return size;
}

void allgather(void* sendbuf, int sendcount, MPI_Datatype datatype, void* recvbuf, int recvcount)
{
	THROW_IF_FAILS(
			MPI_Allgather(sendbuf, sendcount, datatype, recvbuf, recvcount, datatype, MPI_COMM_WORLD)
	);
}

void alltoallv(void* sendbuf, int* sendcounts, int* sdispls, MPI_Datatype datatype, void* recvbuf, int* recvcounts,
		int* rdispls)
{
	THROW_IF_FAILS(
			MPI_Alltoallv(sendbuf, sendcounts, sdispls, datatype, recvbuf, recvcounts, rdispls, datatype,
					MPI_COMM_WORLD)
	);
}

void barrier()
{
	THROW_IF_FAILS(MPI_Barrier(MPI_COMM_WORLD));
}

void broadcast(void* buffer, int count, MPI_Datatype datatype)
{
	THROW_IF_FAILS(
			MPI_Bcast(buffer, count, datatype, 0, MPI_COMM_WORLD)
	);
}

void gather(void* sendbuf, int sendcount, MPI_Datatype datatype, void* recvbuf, int recvcount)
{
	THROW_IF_FAILS(
			MPI_Gather(sendbuf, sendcount, datatype, recvbuf, recvcount, datatype, 0, MPI_COMM_WORLD)
	);
}

void reduce(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op)
{
	THROW_IF_FAILS(
			MPI_Reduce(sendbuf, recvbuf, count, datatype, op, 0, MPI_COMM_WORLD)
	);
}

Type contiguous(int count, MPI_Datatype datatype)
{
	MPI_Datatype result;
	THROW_IF_FAILS(
			MPI_Type_contiguous(count, datatype, &result) || MPI_Type_commit(&result)
	);
	return Type(result);
}

Type structure(int count, int* blockLengths, MPI_Aint* displacements, MPI_Datatype* types)
{
	MPI_Datatype result;
	THROW_IF_FAILS(
			MPI_Type_create_struct(count, blockLengths, displacements, types, &result) || MPI_Type_commit(&result)
	);
	return Type(result);
}

Type subarray(int dims, int* sizes, int* subsizes, int* starts, MPI_Datatype type)
{
	MPI_Datatype result;
	THROW_IF_FAILS(
			MPI_Type_create_subarray(dims, sizes, subsizes, starts, MPI_ORDER_FORTRAN, type, &result)
					|| MPI_Type_commit(&result)
	);
	return Type(result);
}

Type vector(int blockLength, int stride, MPI_Datatype type)
{
	MPI_Datatype result;
	THROW_IF_FAILS(
			MPI_Type_vector(1, blockLength, stride, type, &result) || MPI_Type_commit(&result)
	);
	return Type(result);
}

//----------------------------------------------------------------------
}

#undef THROW_IF_FAILS
