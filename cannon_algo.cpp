# code
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <iostream>

using namespace std;
//#define N 1000


int main(int argc, char* argv[])
{
	MPI_Comm Snyder_comm;
	MPI_Status status;
	int rank, size;
	int shift;
	int i, j, k;
	int dims[2];
	int periods[2];
	int left, right, up, down;
	double* A, * B, * C;
	double* buf, * tmp;
	double start, end;
	unsigned int iseed = 0;
	int Nl, N;



	printf("\n");
	cout << "enter matrix size" << endl;
	cin >> N;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	srand(iseed);
	dims[0] = 0; dims[1] = 0;
	periods[0] = 1; periods[1] = 1;
	MPI_Dims_create(size, 2, dims);
	if (dims[0] != dims[1]) {
		if (rank == 0) printf("The number of processors must be a square.\n");
		MPI_Finalize();
		return 0;
	}
	Nl = N / dims[0];
	
	A = (double*)malloc(Nl * Nl *sizeof(double));
	B = (double*)malloc(Nl * Nl * sizeof(double));
	buf = (double*)malloc(Nl * Nl * sizeof(double));
	C = (double*)calloc(Nl * Nl, sizeof(double));
	for (i = 0; i < Nl; i++)
	{
		for (j = 0; j < Nl; j++) 
		{
			A[i * Nl + j] = rand() % 10;
			B[i * Nl + j] = rand() % 10;
			C[i * Nl + j] = (double)0;
		}
	}
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &Snyder_comm);
	MPI_Cart_shift(Snyder_comm, 0, 1, &left, &right);
	MPI_Cart_shift(Snyder_comm, 1, 1, &up, &down);
	start = MPI_Wtime();
	for (shift = 0; shift < dims[0]; shift++) 
	{
		// Matrix multiplication
		for (i = 0; i < Nl; i++)
		{
			for (k = 0; k < Nl; k++)
			{
				double val = 0;
				for (j = 0; j < Nl; j++)
				{
					val += (A[i * Nl + k] * B[k * Nl + j]);
				}
				C[i * Nl + j] = val;
				//cout << "result :" << C[i * Nl + j] << endl;
				if (shift == dims[0] - 1) break;
				// Communication
				MPI_Sendrecv(A, Nl * Nl, MPI_DOUBLE, left, 1, buf, Nl * Nl, MPI_DOUBLE, right, 1, Snyder_comm, &status);
				tmp = buf; buf = A; A = tmp;
				MPI_Sendrecv(B, Nl * Nl, MPI_DOUBLE, up, 2, buf, Nl * Nl, MPI_DOUBLE, down, 2, Snyder_comm, &status);
				tmp = buf; buf = B; B = tmp;
			}
		}
		
	}
	for (int i = 0; i < N - 1; i++)
	{
		for (int j = 0; j < N-1; j++) 
		{
			C[i * Nl + j];
			
		}
		cout << "Matrix_C[i] =" << C[i * Nl + j] << endl;
	}
	MPI_Barrier(Snyder_comm);
	end = MPI_Wtime();


	if (rank == 0)
		printf("Time: %.4fms\n", (end - start) / 1024.0);
	//cout << "Elapse_Time"<<(end - start) / 1024.0 << "ms" << endl;
	free(A); free(B); free(buf); //free(C);
	MPI_Finalize();
	return 0;
}
