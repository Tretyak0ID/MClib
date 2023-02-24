#include <iostream>
#include <iomanip>
#include "Matrix.cpp"

Matrix<double> StrassenMult2x2double(const Matrix<double>& left, const Matrix<double>& right){
    Matrix<double> Comp(left.get_m(), left.get_n());

    double D  = (left.get_Mij(0, 0) + left.get_Mij(1, 1)) * (right.get_Mij(0, 0) + right.get_Mij(1, 1));
    double D1 = (left.get_Mij(0, 1) - left.get_Mij(1, 1)) * (right.get_Mij(1, 0) + right.get_Mij(1, 1));
    double D2 = (left.get_Mij(1, 0) - left.get_Mij(0, 0)) * (right.get_Mij(0, 0) + right.get_Mij(0, 1));
    double H1 = (left.get_Mij(0, 0) + left.get_Mij(0, 1)) * right.get_Mij(1, 1);
    double H2 = (left.get_Mij(1, 0) + left.get_Mij(1, 1)) * right.get_Mij(0, 0);
    double V1 = left.get_Mij(1, 1) * (right.get_Mij(1, 0) - right.get_Mij(0, 0));
    double V2 = left.get_Mij(0, 0) * (right.get_Mij(0, 1) - right.get_Mij(1, 1));

    Comp.set_Mij(0, 0, D + D1 + V1 - H1);
    Comp.set_Mij(0, 1, V2 + H1);
    Comp.set_Mij(1, 0, V1 + H2);
    Comp.set_Mij(1, 1, D + D2 + V2 - H2);

    return Comp;
}

Matrix<Matrix<double>> MatrixSplit(const Matrix<double>& A){
    int m = A.get_m() / 2;
    int n = A.get_n() / 2;
    Matrix<double> A11(m, n), A12(m, n), A21(m, n), A22(m, n);
    Matrix<Matrix<double>> Mat(2, 2);
    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            A11.set_Mij(i, j, A.get_Mij(i, j));
            A21.set_Mij(i, j, A.get_Mij(i + m, j));
            A12.set_Mij(i, j, A.get_Mij(i, j + n)); 
            A22.set_Mij(i, j, A.get_Mij(i + m, j + n));
        }
    }
    Mat.set_Mij(0, 0, A11);
    Mat.set_Mij(0, 1, A12);
    Mat.set_Mij(1, 0, A21);
    Mat.set_Mij(1, 1, A22);
    return Mat;
}

Matrix<double> MatrixConcat(const Matrix<double>& A11, const Matrix<double> A12, const Matrix<double> A21, const Matrix<double>& A22){
    int n = A11.get_n() + A12.get_n();
    int m = A11.get_m() + A21.get_m();
    Matrix<double> Mat(m, n);
    for (int i = 0; i < A11.get_m(); i++){
        for (int j = 0; j < A11.get_n(); j++){
            Mat.set_Mij(i, j, A11.get_Mij(i, j));
            Mat.set_Mij(i + A11.get_m(), j, A21.get_Mij(i, j));
            Mat.set_Mij(i, j + A11.get_n(), A12.get_Mij(i, j));
            Mat.set_Mij(i + A11.get_m(), j + A11.get_n(), A22.get_Mij(i, j));
        }
    }
    return Mat;
}

Matrix<double> StrassenMult(const Matrix<double>& left, const Matrix<double>& right){
    Matrix<double> Comp(left.get_m(), left.get_n());

    if (left.get_n() == 2){
        Comp = StrassenMult2x2double(left, right);
    }else{
        Matrix<Matrix<double>> block_mat_left = MatrixSplit(left);
        Matrix<Matrix<double>> block_mat_right = MatrixSplit(right);

        Matrix<double> D  = StrassenMult((block_mat_left.get_Mij(0, 0) + block_mat_left.get_Mij(1, 1)), (block_mat_right.get_Mij(0, 0) + block_mat_right.get_Mij(1, 1)));
        Matrix<double> D1 = StrassenMult((block_mat_left.get_Mij(0, 1) - block_mat_left.get_Mij(1, 1)), (block_mat_right.get_Mij(1, 0) + block_mat_right.get_Mij(1, 1)));
        Matrix<double> D2 = StrassenMult((block_mat_left.get_Mij(1, 0) - block_mat_left.get_Mij(0, 0)), (block_mat_right.get_Mij(0, 0) + block_mat_right.get_Mij(0, 1)));
        Matrix<double> H1 = StrassenMult((block_mat_left.get_Mij(0, 0) + block_mat_left.get_Mij(0, 1)), block_mat_right.get_Mij(1, 1));
        Matrix<double> H2 = StrassenMult((block_mat_left.get_Mij(1, 0) + block_mat_left.get_Mij(1, 1)), block_mat_right.get_Mij(0, 0));
        Matrix<double> V1 = StrassenMult(block_mat_left.get_Mij(1, 1), (block_mat_right.get_Mij(1, 0) - block_mat_right.get_Mij(0, 0)));
        Matrix<double> V2 = StrassenMult(block_mat_left.get_Mij(0, 0), (block_mat_right.get_Mij(0, 1) - block_mat_right.get_Mij(1, 1)));

        Matrix<double> C11 = D + D1 + V1 - H1;
        Matrix<double> C12 =  V2 + H1;
        Matrix<double> C21 =  V1 + H2;
        Matrix<double> C22 =  D + D2 + V2 - H2;

        Comp = MatrixConcat(C11, C12, C21, C22);
    }
    return Comp;
}