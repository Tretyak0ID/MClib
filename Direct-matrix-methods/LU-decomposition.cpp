#include <iostream>
#include <iomanip>
#include "../Matrix.cpp"

void GaussLU(Matrix<double>& A, Matrix<double>& L, Matrix<double>& U){
    if (A.get_m() > 0){
        double a = A.get_Mij(0, 0);
        Matrix<double> c(1, A.get_n() - 1);
        Matrix<double> b(A.get_m() - 1, 1);
        Matrix<double> A1(A.get_m() - 1, A.get_n() - 1), L1(A.get_m() - 1, A.get_n() - 1), U1(A.get_m() - 1, A.get_n() - 1);

        for(int i = 0; i < A.get_m() - 1; i++){
            b.set_Mij(i, 0, A.get_Mij(i + 1, 0));
            L.set_Mij(i + 1, 0, 1 / a * b.get_Mij(i,0));
            U.set_Mij(i + 1, 0, 0.0);
        }

        for(int i = 0; i < A.get_n() - 1; i++){
            c.set_Mij(0, i, A.get_Mij(0, i + 1));
            L.set_Mij(0, i + 1, 0.0);
            U.set_Mij(0, i + 1, c.get_Mij(0, i));
        }

        L.set_Mij(0, 0, 1.0);
        U.set_Mij(0, 0, a);
        A1 = A.getDownMinor() - (b * c) * (1 / a);

        GaussLU(A1, L1, U1);

        for (int i = 0; i <  A.get_m() - 1; i++){
            for (int j = 0; j <  A.get_n() - 1; j++){
                L.set_Mij(i + 1, j + 1, L1.get_Mij(i, j));
                U.set_Mij(i + 1, j + 1, U1.get_Mij(i, j));
            }
        }
    }
}

void DirectLU(Matrix<double>& A, Matrix<double>& L, Matrix<double>& U){

    for(int i = 0; i < A.get_m(); i++){
        for(int j = 0; j < A.get_n(); j++){
            if( i == j ) L.set_Mij(i, j, 1.0);
            else L.set_Mij(i, j, 0.0);
            U.set_Mij(i, j, 0.0);;
        }
    }

    for (int i = 0; i < A.get_m(); i++){
        for(int j = 0; j < A.get_n(); j++){
            double sumU = 0.0;
            double sumL = 0.0;

            if (i <= j){
              U.set_Mij(i, j, A.get_Mij(i, j);
              for (int k = 0; k < i-1; k++){
                  sumU += L.get_Mij(i, k)*U.get_Mij(k, j);
                  U.set_Mij(i, j, U.get_Mij(i, j) - L.get_Mij(i, k)*U.get_Mij(k, j));
              }
            }

            if (i > j){
              for (int k = 0; k < j-1; k++){
                  sumL += L.get_Mij(i, k)*U.get_Mij(k, j);
              }
              L.set_Mij(i, j, (A.get_Mij(i, j) - sumL) / U.get_Mij(j, j));
            }
        }
    }
}
