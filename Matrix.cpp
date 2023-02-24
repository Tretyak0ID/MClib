#include <iostream>

template <typename T>
class Matrix{
private:
    T* M;
    int m;
    int n;

public:
    // конструкторы
    Matrix(){
        n = m = 0;
        M = nullptr;
    }

    Matrix(int _m, int _n){
        m = _m;
        n = _n;

        M = (T*) new T[m * n];
    }

    Matrix(const Matrix& _M){
        m = _M.m;
        n = _M.n;

        M = (T*) new T[m * n];

    }

    ~Matrix(){
        
        if (m * n > 0)
            delete[] M;
    }

    // методы доступа
    T get_Mij (int i, int j) const {
            return M[i * n + j];
    }

    void set_Mij(int i, int j, T value){
        /*if ((i < 0) || (i >= m))
            return;
        if ((j < 0) || (j >= n))
            return;*/
        M[i * n + j] = value;
    }

    int get_m() const {
        return m;
    }

    int get_n() const {
        return n;
    }

    // Основные методы
    void fill_rand(int s, int e){
        for (int i = 0; i < m * n; i++)
            M[i] = (rand() % (e - s)) + s;
    }

    void fill_diag(double s){
        if(m >= n){
            for (int i = 0; i < n; i++)
                M[i * n + i] = s;
        }else{
            for (int i = 0; i < m; i++)
                M[i * n + i] = s;
        }
    }

    Matrix getDownMinor(){
        Matrix<double> Min(m - 1, n - 1);

        for (int i = 1; i < m; i++){
            for (int j = 1; j < n; j++){
                Min.M[(i - 1) * n + j - 1] = M[i * n + j];
            }
        }

        return Min;
    }

    void print(){
        std::cout << std::endl;
        for (int i = 0; i < m; i++){
            for (int j = 0; j < n; j++)
                std::cout << std::setw(10)<<M[i * n + j];
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    // Перегрузка операторов
    Matrix operator=(const Matrix& _M){
        

        if (m * n > 0){
            delete[] M;
        }

        m = _M.m;
        n = _M.n;

        M = (T*) new T[m * n];

        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                M[i * n + j] = _M.M[i * n + j];

        return *this;
    }

    Matrix operator*(const Matrix& Mr){
        Matrix Comp(m, Mr.n);
        
        for (int i = 0; i < m; i++){
            for (int j = 0; j < Mr.n; j++){
                for (int k = 0; k < n; k++){
                    Comp.M[i * Mr.n + j] = Comp.M[i * Mr.n + j] + M[i * n + k]*Mr.M[k * Mr.n + j];
                }
            }
        }
        return Comp;
    }

    Matrix operator*(double s){
        Matrix Comp(m, n);

        for (int i = 0; i < m; i++){
            for (int j = 0; j < n; j++){
                M[i * n + j] = M[i * n + j]*s;
            }
        }

        return Comp;
    }

    Matrix operator+(const Matrix& Mr){
        Matrix Sum(m, Mr.n);

        for (int i = 0; i < m; i++){
            for (int j = 0; j < Mr.n; j++){
                Sum.M[i * n + j] = M[i * n + j] + Mr.M[i * n + j];
            }
        }

        return Sum;
    }

    Matrix operator-(const Matrix& Mr){
        Matrix Sum(m, Mr.n);

        for (int i = 0; i < m; i++){
            for (int j = 0; j < Mr.n; j++){
                Sum.M[i * n + j] = M[i * n + j] - Mr.M[i * n + j];
            }
        }

        return Sum;
    }

};
