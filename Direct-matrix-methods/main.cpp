#include <iostream>
#include "mpi.h"
#include "LU-decomposition.cpp"

int main(){
    int SIZE = 4;
    Matrix<double> Mat(SIZE,SIZE), L(SIZE,SIZE), U(SIZE,SIZE);
    Mat.fill_rand(1, 50);

    double start_time = MPI_Wtime();

    DirectLU(Mat, L, U);

    double end_time = MPI_Wtime();
    Mat.print();

    std::cout << end_time - start_time << std::endl;

    return 0;
}
