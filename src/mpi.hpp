// Coulombo Ⓒ 2018
// [Computer Physics Communications] Różański & Zieliński:
// Efficient computation of Coulomb and exchange integrals for very large systems
#ifndef COULOMBO_MPI_HPP
#define COULOMBO_MPI_HPP

#include <mpi.h>

namespace mpi {

class Type {
public:
	explicit Type(MPI_Datatype datatype = nullptr);

	Type(Type&& source);

	~Type();

	void operator=(Type&& source);

	operator MPI_Datatype();

private:
	MPI_Datatype datatype;

	void operator=(const Type&) = delete;
	Type(const Type&) = delete;
};

void init(int& argc, char**& argv);

void finalize(void);

int size(void);

int rank(void);

bool root(void);

void allgather(void* sendbuf, int sendcount, MPI_Datatype datatype, void* recvbuf, int recvcount);

void alltoallv(void* sendbuf, int* sendcounts, int* sdispls, MPI_Datatype datatype, void* recvbuf,
		int* recvcounts, int* rdispls);

void barrier();

void broadcast(void* buffer, int count, MPI_Datatype datatype);

void reduce(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op);

Type contiguous(int count, MPI_Datatype datatype);

Type structure(int count, int* blockLengths, MPI_Aint* displacements, MPI_Datatype* types);

Type subarray(int dims, int* sizes, int* subsizes, int* starts, MPI_Datatype type);

Type vector(int blockLength, int stride, MPI_Datatype type);

}

#endif
