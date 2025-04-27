#include <iostream>
#include <vector>

class Matrix {
private:
    size_t rows, cols;
    std::vector<std::vector<double>> data;

public:
    // Конструктор для создания матрицы заданного размера
    Matrix(size_t rows, size_t cols, double initial_value = 0.0) 
        : rows(rows), cols(cols), data(rows, std::vector<double>(cols, initial_value)) {}

    // Конструктор для создания матрицы из двумерного вектора
    Matrix(const std::vector<std::vector<double>>& input_data) {
        if (input_data.empty() || input_data[0].empty()) {
            throw std::invalid_argument("Матрица не может быть пустой");
        }
        rows = input_data.size();
        cols = input_data[0].size();
        data = input_data;
    }

    // Оператор сложения матриц
    Matrix operator+(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols) {
            throw std::invalid_argument("Матрицы должны быть одного размера для сложения");
        }
        Matrix result(rows, cols);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                result.data[i][j] = data[i][j] + other.data[i][j];
            }
        }
        return result;
    }

    // Оператор вычитания матриц
    Matrix operator-(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols) {
            throw std::invalid_argument("Матрицы должны быть одного размера для вычитания");
        }
        Matrix result(rows, cols);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                result.data[i][j] = data[i][j] - other.data[i][j];
            }
        }
        return result;
    }

    // Оператор умножения матриц
    Matrix operator*(const Matrix& other) const {
        if (cols != other.rows) {
            throw std::invalid_argument("Количество столбцов первой матрицы должно совпадать с количеством строк второй");
        }
        Matrix result(rows, other.cols);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < other.cols; ++j) {
                for (size_t k = 0; k < cols; ++k) {
                    result.data[i][j] += data[i][k] * other.data[k][j];
                }
            }
        }
        return result;
    }

    // Умножение матрицы на скаляр
    Matrix operator*(double scalar) const {
        Matrix result(rows, cols);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                result.data[i][j] = data[i][j] * scalar;
            }
        }
        return result;
    }

    // Транспонирование матрицы
    Matrix transpose() const {
        Matrix result(cols, rows);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                result.data[j][i] = data[i][j];
            }
        }
        return result;
    }

    // Вывод матрицы на экран
    void print() const {
        for (const auto& row : data) {
            for (const auto& elem : row) {
                std::cout << elem << " ";
            }
            std::cout << std::endl;
        }
    }

    // Для того чтобы считать матрицы L и U. Нужно написать GETelem SETelem
    double GETelem(size_t get_i, size_t get_j) const {
        if (get_i >= rows || get_j >= cols) {
            throw std::out_of_range("Неверный индекс матрицы");
        }
        else{
            return data[get_i][get_j];
        }
    }
    
    double SETelem(size_t set_i, size_t set_j, double set_value) {
        if (set_i >= rows || set_j >= cols) {
            throw std::out_of_range("Неверный индекс матрицы");
        }
        else {
            data[set_i][set_j] = set_value;
        }
    }
    //======================================================================================================================//
    void LU_razl() {
        if (rows != cols) {
            throw std::invalid_argument("Матрица должна быть квадратной");
        }
        size_t n = rows;
    
        // Разложение
        for (size_t k = 0; k < n; k++) {
            // Разложение U
            for (size_t j = k; j < n; j++) {
                double sum_u = 0.0;
                for (size_t m = 0; m < k; m++) {
                    sum_u += GETelem(k, m) * GETelem(m, j);
                }
                SETelem(k, j, GETelem(k, j) - sum_u);
            }
    
            // Разложение L
            for (size_t i = k + 1; i < n; i++) {
                double sum_l = 0.0;
                for (size_t m = 0; m < k; m++) {
                    sum_l += GETelem(i, m) * GETelem(m, k);
                }
                SETelem(i, k, (GETelem(i, k) - sum_l) / GETelem(k, k));
            }
        }
    }
    //===============================================================================================================//
    // Заполняем L
    Matrix getL() const {
        Matrix L(rows, cols);
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                if (i > j) {
                    L.SETelem(i, j, data[i][j]);
                } else if (i == j) {
                    L.SETelem(i, j, 1.0);
                } else {
                    L.SETelem(i, j, 0.0);
                }
            }
        }
        return L;
    }
    // Заполняем U
    Matrix getU() const {
        Matrix U(rows, cols);
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                if (i <= j) {
                    U.SETelem(i, j, data[i][j]);
                } else {
                    U.SETelem(i, j, 0.0);
                }
            }
        }
        return U;
    }

    std::vector<double> LY(const std::vector<double>& b) const {
        size_t n = rows;
        std::vector<double> y(n, 0.0);

        for (size_t i = 0; i < n; i++) {
            double sum_y = 0.0;
            for (size_t j = 0; j < i; j++) {
                sum_y += GETelem(i, j) * y[j];
            }
            y[i] = b[i] - sum_y;
        }
        return y;
    }

    std::vector<double> UX(const std::vector<double>& y) const {
        size_t n = rows;
        std::vector<double> x(n, 0.0);

        for (int i = n - 1; i >= 0; i--) {
            double sum = 0.0;
            for (size_t j = i + 1; j < n; j++) {
                sum += GETelem(i, j) * x[j];
            }
            x[i] = (y[i] - sum) / GETelem(i, i);
        }
        return x;
    }

    std::vector<double> SYSTEM(const std::vector<double>& b) {
        LU_razl();
        std::vector<double> y = LY(b);
        std::vector<double> x = UX(y);
        return x;
    }
};

int main() {

    Matrix A({{1, 2, 1}, {4, 5, 6}, {7, 8, 9}});
    std::vector<double> b1 = {1, 2, 3};
    std::vector<double> b2 = {4, 5, 6};

    std::cout << "Матрица A:" << std::endl;
    A.print();

    std::cout << "Вектор b1:" << std::endl;
    for (const auto& elem : b1) {
        std::cout << elem << " ";
    }
    std::cout << std::endl;

    std::cout << "Вектор b2:" << std::endl;
    for (const auto& elem : b2) {
        std::cout << elem << " ";
    }
    std::cout << std::endl;

    std::vector<double> x = A.SYSTEM(b1);

    std::cout << "Вектор x1:" << std::endl;
    for (const auto& elem : x) {
        std::cout << elem << " ";
    }
    std::cout << std::endl;

    x = A.SYSTEM(b2);

    std::cout << "Вектор x2:" << std::endl;
    for (const auto& elem : x) {
        std::cout << elem << " ";
    }
    std::cout << std::endl;
    return 0;
}