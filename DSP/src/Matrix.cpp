// #include <cstdlib> // malloc, realloc, free
#include <string> // memcpy, memset, memmove
#include <sstream> // stringstream
#include <iomanip>
#include "DSP/Math/Matrix.hpp"

using namespace DSP::Math;

// Constructor
Matrix::Matrix():
rows(0),
columns(0),
elements(nullptr) 
{}

// Constructor
Matrix::Matrix(const size_t rows, const size_t columns):
rows(0),
columns(0),
elements(nullptr) 
{
    if (rows > 0 && columns > 0) {
        this->elements = new real_t[rows * columns];
        if (this->elements != nullptr) {
            this->rows = rows;
            this->columns = columns;
        }
    }
}

// Constructor
Matrix::Matrix(const size_t rows, const size_t columns, const real_t* elements):
rows(0),
columns(0),
elements(nullptr) 
{
    if (rows > 0 && columns > 0) {
        this->elements = new real_t[rows * columns];
        if (this->elements != nullptr) {
            this->rows = rows;
            this->columns = columns;

            if (elements != nullptr) { memcpy(this->elements, elements, rows * columns * sizeof(real_t)); }
            else { memset(this->elements, 0, rows * columns * sizeof(real_t)); }
        }
    }
}

Matrix::Matrix(const size_t rows, const size_t columns, const std::vector<real_t>& elements):
rows(0),
columns(0),
elements(nullptr) 
{
    if (rows > 0 && columns > 0 && (elements.size() == rows * columns)) {
        this->elements = new real_t[elements.size()];
        if (this->elements != nullptr) {
            this->rows = rows;
            this->columns = columns;
            memcpy(this->elements, elements.data(), elements.size() * sizeof(real_t));
        }
    }
}

// Constructor
Matrix::Matrix(const std::vector<real_t>& row_vec):
rows(0),
columns(0),
elements(nullptr) 
{
    if (row_vec.size() > 0) {
        this->elements = new real_t[row_vec.size()];
        if (this->elements != nullptr) {
            this->rows = row_vec.size();
            this->columns = 1;
            memcpy(this->elements, row_vec.data(), row_vec.size() * sizeof(real_t));
        }
    }
}

// Initilize
Matrix::Matrix(const std::initializer_list<real_t>& row_il):
rows(0),
columns(0),
elements(nullptr) 
{
    if (row_il.size() > 0) {
        this->elements = new real_t[row_il.size()];
        if (this->elements != nullptr) {
            this->rows = 1;
            this->columns = row_il.size();

            size_t k = 0;
            for (const auto& il_element : row_il) {
                this->elements[k] = il_element;
                k += 1;
            }
        }
    }
}

// Initilize
void Matrix::operator=(const std::initializer_list<real_t>& row_il) {

    // Release old array
    this->release_internal_array();

    if (row_il.size() > 0) {
        this->elements = new real_t[row_il.size()];
        if (this->elements != nullptr) {
            this->rows = 1;
            this->columns = row_il.size();

            size_t k = 0;
            for (const auto& il_element : row_il) {
                this->elements[k] = il_element;
                k += 1;
            }
        }
    }

}

// Copy
Matrix::Matrix(const Matrix& other):
rows(0),
columns(0),
elements(nullptr)
{
    if (other.rows > 0 && other.columns > 0 && other.elements != nullptr) {
        this->elements = new real_t[other.rows * other.columns];
        if (this->elements != nullptr) {
            this->rows = other.rows;
            this->columns = other.columns;
            memcpy(this->elements, other.elements, other.rows * other.columns * sizeof(real_t));
        }
    }
}

// Copy
void Matrix::operator=(const Matrix& other) {
    if (this != &other) {

        // Release old array
        this->release_internal_array();

        if (other.rows > 0 && other.columns > 0 && other.elements != nullptr) {
            this->elements = new real_t[other.rows * other.columns];
            if (this->elements != nullptr) {
                this->rows = other.rows;
                this->columns = other.columns;
                memcpy(this->elements, other.elements, other.rows * other.columns * sizeof(real_t));
            }
        }
    }
}

// Move
Matrix::Matrix(Matrix&& other):
rows(0),
columns(0),
elements(nullptr)
{
    // Move
    this->rows = other.rows;
    this->columns = other.columns;
    this->elements= other.elements;

    // Invalidate
    other.rows = 0;
    other.columns = 0;
    other.elements = nullptr;
}

void Matrix::operator=(Matrix&& other) {
    if (this != &other) {

        // Release old array
        this->release_internal_array();

        // Move
        this->rows = other.rows;
        this->columns = other.columns;
        this->elements= other.elements;

        // Invalidate
        other.rows = 0;
        other.columns = 0;
        other.elements = nullptr;
    }
}

// Destructor
void Matrix::release_internal_array() {
    // Release array
    if (this->elements != nullptr) { 
        delete[] this->elements; 
        this->elements = nullptr; 
    }

    // Reset 
    this->rows = 0;
    this->columns = 0;
}

// Destructor
Matrix::~Matrix() {
    this->release_internal_array();
}

// to bool
Matrix::operator bool() const {
    return (this->elements != nullptr);
}

// to std::vector
Matrix::operator std::vector<real_t>() const {
    return std::vector<real_t>(this->elements, this->elements + this->get_size());
}

// to string
std::string Matrix::to_string() const {
    std::stringstream buff;
    buff << "{" << std::endl;
    // buff << std::setw(12) << std::fixed;
    for (size_t i = 0; i < this->rows; ++i) {
        for (size_t j = 0; j < this->columns; ++j) {
            
            if (j < this->columns - 1) { buff << std::setw(12) << std::fixed << this->m(i, j) << ", "; }
            else { buff << std::setw(12) << std::fixed << this->m(i, j) << std::endl; }
        }
    }
    buff << "}" << std::endl;
    return buff.str();
}
Matrix::operator std::string() const {
    return this->to_string();
}

// reshape
Matrix Matrix::reshape(const size_t rows, const size_t columns) const {
    if (this->get_size() != (rows * columns)) { return Matrix(); }
    Matrix other(*this);
    other.rows = rows;
    other.columns = columns;
    return other;
}
Matrix& Matrix::reshape(const size_t rows, const size_t columns) {
    if (this->get_size() == (rows * columns)) {
        this->rows = rows;
        this->columns = columns;
    }
    return *this;
}

// Element access 1-based
real_t& Matrix::operator()(const size_t row, const size_t column) {
    return this->elements[(row - 1) * this->columns + (column - 1)];
}

// Element access 1-based
const real_t& Matrix::operator()(const size_t row, const size_t column) const {
    return this->elements[(row - 1) * this->columns + (column - 1)];
}

// Element access 0-based
real_t& Matrix::m(const size_t row_index, const size_t column_index) {
    return this->elements[row_index * this->columns + column_index];
}

// Element access 0-based
const real_t& Matrix::m(const size_t row_index, const size_t column_index) const {
    return this->elements[row_index * this->columns + column_index];
}

// Element access 1-based
real_t& Matrix::at(const size_t row, const size_t column) {
    return this->elements[(row - 1) * this->columns + (column - 1)];
}

const real_t& Matrix::at(const size_t row, const size_t column) const {
    return this->elements[(row - 1) * this->columns + (column - 1)];
}

// Element access
real_t& Matrix::operator[](const size_t index) {
    return this->elements[index];
}

// Element access
const real_t& Matrix::operator[](const size_t index) const {
    return this->elements[index];
}


// Getter
size_t Matrix::get_num_rows() const {
    return this->rows;
}
size_t Matrix::get_num_columns() const {
    return this->columns;
}
size_t Matrix::get_size() const {
    return this->rows * this->columns;
}
real_t* Matrix::get_elements() {
    return this->elements;
}
const real_t* Matrix::get_elements() const {
    return this->elements;
}


// Getter 0-based
Matrix Matrix::row_vector(const size_t row_index) const {
    Matrix vec(1, this->columns);
    for (size_t k = 0; k < this->columns; ++k) {
        vec[k] = this->m(row_index, k);
    }
    return vec;
}

// Getter 0-based
Matrix Matrix::column_vector(const size_t column_index) const {
    Matrix vec(this->rows, 1);
    for (size_t k = 0; k < this->rows; ++k) {
        vec[k] = this->m(k, column_index);
    }
    return vec;
}

// Getter 0-based
Matrix Matrix::submatrix(const size_t row_index, const size_t column_index) const {
    Matrix sub(this->rows - 1, this->columns - 1);

    size_t k = 0, l;
    for (size_t i = 0; i < this->rows; ++i) {
        if (i == row_index) { continue; }

        l = 0;
        for (size_t j = 0; j < this->columns; ++j) {
            if (j == column_index) { continue; }
                sub.m(k, l) = this->m(i, j);
            ++l;
        }
        ++k;
    }

    return sub;
}


// Arithmetic
Matrix Matrix::operator+() const {
    return *this;
}
Matrix Matrix::operator-() const {
    Matrix neg(this->rows, this->columns);
    const size_t size = this->rows * this->columns;
    for (size_t k = 0; k < size; ++k) {
        neg.elements[k] = -(this->elements[k]);
    }

    return neg;
}

// Arithmetic
Matrix Matrix::operator+(const Matrix& other) const {
    if (this->rows != other.rows || this->columns != other.columns) { return Matrix(); }

    Matrix res(this->rows, this->columns);
    for (size_t k = 0; k < (this->rows * this->columns); ++k) {
        res.elements[k] = this->elements[k] + other.elements[k];
    }

    return res;
}
Matrix Matrix::operator-(const Matrix& other) const {
    if (this->rows != other.rows || this->columns != other.columns) { return Matrix(); }

    Matrix res(this->rows, this->columns);
    for (size_t k = 0; k < (this->rows * this->columns); ++k) {
        res.elements[k] = this->elements[k] - other.elements[k];
    }

    return res;
}
Matrix Matrix::operator*(const Matrix& other) const {
    if (this->columns != other.rows) { return Matrix(); }

    // Create Matrix
    Matrix res(this->rows, other.columns);

    // Calculate product
    for (size_t i = 0; i < res.rows; ++i) {
        for (size_t j = 0; j < res.columns; ++j) {

            // dot product
            res.m(i, j) = this->m(i, 0) * other.m(0, j);
            for (size_t k = 1; k < this->columns; ++k) {
                res.m(i, j) += this->m(i, k) * other.m(k, j);
            }
        }
    }

    return res;
}
Matrix Matrix::operator/(const Matrix& other) const { // (A / B) == (inv(B) * A)
    if (this->rows != other.rows) { return Matrix(); }
    if (other.columns == 1) { return Matrix::solve(*this, other); }
    else { return other.pinv() * (*this); }
}

// Arithmetic
void Matrix::operator+=(const Matrix& other) {
    *this = ((*this) + other);
}
void Matrix::operator-=(const Matrix& other) {
    *this = ((*this) - other);
}
void Matrix::operator*=(const Matrix& other) {
    *this = ((*this) * other);
}
void Matrix::operator/=(const Matrix& other) {
    *this = ((*this) / other);
}

// Arithmetic
Matrix Matrix::operator*(const real_t scalar) const {
    Matrix res(this->rows, this->columns);

    const size_t size = this->rows * this->columns;
    for (size_t k = 0; k < size; ++k) {
        res.elements[k] = this->elements[k] * scalar;
    }

    return res;

}
Matrix Matrix::operator/(const real_t scalar) const {
    Matrix res(this->rows, this->columns);

    const size_t size = this->rows * this->columns;
    for (size_t k = 0; k < size; ++k) {
        res.elements[k] = this->elements[k] / scalar;
    }

    return res;

}
void Matrix::operator*=(const real_t scalar) {
    const size_t size = this->rows * this->columns;
    for (size_t k = 0; k < size; ++k) {
        this->elements[k] *= scalar;
    }
}
void Matrix::operator/=(const real_t scalar) {
    const size_t size = this->rows * this->columns;
    for (size_t k = 0; k < size; ++k) {
        this->elements[k] /= scalar;
    }
}


// Determinante
real_t Matrix::det() const {
    return Matrix::det(*this);
}

// Inverse 
// row == columns -> Gauß, 
// row > column -> ~M = ~(C^T * C) * C^T
// row < column -> ~M = C^T * ~(C * C^T)
Matrix Matrix::operator~() const {
    return Matrix::pinv(*this);
}

Matrix Matrix::inv() const {
    return Matrix::inv(*this);
}

Matrix Matrix::pinv() const {
    return Matrix::pinv(*this);
}

// Transpose
Matrix Matrix::transpose() const {
    Matrix res(this->columns, this->rows);

    for (size_t i = 0; i < this->rows; ++i) {
        for (size_t j = 0; j < this->columns; ++j) {
            res.m(j, i) = this->m(i, j);
        }
    }

    return res;
}

// Diagonal
std::vector<real_t> Matrix::diag() const {

    const size_t size = (this->rows <= this->columns ? this->rows : this->columns);

    std::vector<real_t> diagonal_elements;
    diagonal_elements.reserve(size);
    for (size_t k = 0; k < size; ++k) { diagonal_elements.push_back(this->m(k, k)); }

    return diagonal_elements;
}

// solve (A * x = b)
Matrix Matrix::solve(const Matrix& b) const {
    return Matrix::solve(*this, b);
}

// 0-based
void Matrix::swap_rows(const size_t row_a_index, const size_t row_b_index) {
    real_t temp;
    for (size_t k = 0; k < this->columns; ++k) {
        temp = this->m(row_a_index, k);
        this->m(row_a_index, k) = this->m(row_b_index, k);
        this->m(row_b_index, k) = temp;
    }
}

// 0-based
void Matrix::swap_columns(const size_t column_a_index, const size_t column_b_index) {
    real_t temp;
    for (size_t k = 0; k < this->rows; ++k) {
        temp = this->m(k, column_a_index);
        this->m(k, column_a_index) = this->m(k, column_b_index);
        this->m(k, column_b_index) = temp;
    }
}

// 0-based
void Matrix::scale_up_row(const size_t row_index, const size_t column_a_index, const size_t column_b_index, const real_t s) {
    for (size_t k = column_a_index; k <= column_b_index; ++k) {
        this->m(row_index, k) *= s;
    }
}

// 0-based
void Matrix::scale_up_column(const size_t column_index, const size_t row_a_index, const size_t row_b_index, const real_t s) {
    for (size_t k = row_a_index; k <= row_b_index; ++k) {
        this->m(k, column_index) *= s;
    }
}

// 0-based
void Matrix::scale_down_row(const size_t row_index, const size_t column_a_index, const size_t column_b_index, const real_t s) {
    for (size_t k = column_a_index; k <= column_b_index; ++k) {
        this->m(row_index, k) /= s;
    }
}

// 0-based
void Matrix::scale_down_column(const size_t column_index, const size_t row_a_index, const size_t row_b_index, const real_t s) {
    for (size_t k = row_a_index; k <= row_b_index; ++k) {
        this->m(k, column_index) /= s;
    }
}

// 0-based (row_a += s * row_b)
void Matrix::modify_row(const size_t row_a_index, const size_t row_b_index, const size_t column_a_index, const size_t column_b_index, const real_t s) {
    for (size_t k = column_a_index; k <= column_b_index; ++k) {
        this->m(row_a_index, k) += s * this->m(row_b_index, k);
    }
}

// 0-based (column_a += s * column_b)
void Matrix::modify_column(const size_t column_a_index, const size_t column_b_index, const size_t row_a_index, const size_t row_b_index, const real_t s) {
    for (size_t k = row_a_index; k <= row_b_index; ++k) {
        this->m(k, column_a_index) += s * this->m(k, column_b_index);
    }
}

real_t Matrix::gauss_jordan(Matrix* const A, Matrix* const B, const gauss_jordan_task_t task) {
    if (A == nullptr) { return 0; }
    if (A->rows == 0 || A->rows != A->columns) { return 0; } 
    if (task != calc_det && B == nullptr) { return 0; }
    if (task == solve_lse && B->columns != 1) { return 0; }

    // Vorbereitung
    const size_t size = A->rows;
    real_t alpha = 1, beta;
    if (task == invert_mat) { *B = Matrix::eye(size); }

    // Vorwärts-Elimination
    for (size_t k = 0; k < size - 1; ++k) {

        // Wenn das Element auf der Hauptdiagonalen Null ist
        if (A->m(k, k) == 0) {

            // Suche Zeile
            for (size_t j = k + 1; j < size; ++j) {

                // Wenn Element nicht Null
                if (A->m(j, k) != 0) {

                    // Tausche Zeilen
                    A->swap_rows(k, j);
                    if (task != calc_det) { B->swap_rows(k, j); }

                    // Beim Zeilen tauchen ändert sich das Vorzeichen der Determinante
                    alpha *= (-1); 
                    break;
                }
            }

            // Nicht lösbar
            if (A->m(k, k) == 0) { return 0; } 
        }

        // Nullen unterhalb der Haupdiagonalen erzeugen
        for (size_t j = k + 1; j < size; ++j) {
            beta = (-1) * (A->m(j, k) / A->m(k, k));

            // Null in Spalte k Zeile j erzeugen
            A->m(j, k) = 0;
            A->modify_row(j, k, k+1, size-1, beta);
            if (task != calc_det) { B->modify_row(j, k, 0, B->columns-1, beta); }
        }
    }

    // Determinante berechnen
    for (size_t k = 0; k < size; ++k) { alpha *= A->m(k, k); }

    // Falls nur die Determinante berechnet werden sollte 
    // oder die Determinante Null ist, stoppe hier.
    if (task == calc_det || alpha == 0) { return alpha;}

    // Rückwärts-Elimination
    for (size_t k = size; k > 1; --k) {

        // Nullen oberhalb der Haupdiagonalen erzeugen
        for (size_t j = k - 1; j > 0; --j) {
            beta = (-1) * (A->at(j, k) / A->at(k, k));

            // Null in Spalte k Zeile j erzeugen
            A->at(j, k) = 0;
            // A->modify_row(j-1, k-1, k, size-1, beta);
            B->modify_row(j-1, k-1, 0, B->columns-1, beta);
        }
    }

    // Elemente auf der Hauptdiagonalen zu 1 umformen
    for (size_t k = 0; k < size; ++k) {
        B->scale_down_row(k, 0, B->columns-1, A->m(k, k));
        A->m(k, k) = 1;
    }

    // Wir sind fertig.
    // A enthält jetzt die Einheitsmatrix
    // Und B die Lösung für das LGS bzw. die inverse Matrix
    return alpha; //    
}

// Static
real_t Matrix::det(Matrix M) {
    return gauss_jordan(&M, nullptr, calc_det);
}
Matrix Matrix::inv(Matrix M) {
    Matrix B;
    const real_t det_M = gauss_jordan(&M, &B, invert_mat);
    if (det_M == 0) { return Matrix(); }
    else { return B; }
}
Matrix Matrix::pinv(const Matrix& M) {
    if (M.rows > M.columns) {
        // Überbestimmt: row > column -> ~C = ~(C^T * C) * C^T
        const auto Mt = M.transpose();
        return Matrix::inv(Mt * M) * Mt;
    }
    else if (M.rows < M.columns) {
        // Unterbestimmt: row < column -> ~C = C^T * ~(C * C^T)
        const auto Mt = M.transpose();
        return Mt * Matrix::inv(M * Mt);
    }
    else {
        // row == columns -> Gauß-Jordan
        return Matrix::inv(M);
    }
}
Matrix Matrix::transpose(const Matrix& M) {
    return M.transpose();
}

// solve (A * x = b)
Matrix Matrix::solve(Matrix A, Matrix b) {
    if (b.columns != 1) { return Matrix(); }
    if (A.rows != b.rows) { return Matrix(); }
    if (A.rows == A.columns) {
        const real_t det_A = gauss_jordan(&A, &b, solve_lse);
        if (det_A != 0) { return b; }
        else { return Matrix(); }
    }
    else {
        return Matrix::pinv(A) * b;
    }
}

Matrix Matrix::zeros(const size_t dim) {
    return Matrix::zeros(dim, dim);
}
Matrix Matrix::zeros(const size_t rows, const size_t columns) {
    return Matrix(rows, columns, nullptr);
}
Matrix Matrix::eye(const size_t dim) {
    return Matrix::eye(dim, dim);
}
Matrix Matrix::eye(const size_t rows, const size_t columns) {

    Matrix M = Matrix::zeros(rows, columns);
    const size_t size = (rows <= columns ? rows : columns);
    for (size_t k = 0; k < size; ++k) { M.m(k, k) = 1; }
    return M;
}
Matrix Matrix::diag(const std::vector<real_t>& vec) {
    Matrix M = Matrix::zeros(vec.size());
    for (size_t k = 0; k < vec.size(); ++k) { M.m(k, k) = vec[k]; }
    return M;
}
std::vector<real_t> Matrix::diag(const Matrix& M) {
    return M.diag();
}
Matrix Matrix::vandermonde(const size_t order, const std::vector<real_t>& x_values) {

    Matrix V(x_values.size(), order+1);

    for (size_t i = 0; i < x_values.size(); ++i) {

        V.m(i, 0) = 1;
        V.m(i, 1) = x_values[i];
        for (size_t j = 2; j <= order; ++j) {
            V.m(i, j) = V.m(i, j-1) * x_values[i];
        }
    }

    return V;
}

Matrix Matrix::reshape(const Matrix& mat, const size_t rows, const size_t columns) {
    return mat.reshape(rows, columns);
}
Matrix Matrix::row_vector(const std::vector<real_t>& vec) {
    return Matrix(1, vec.size(), vec.data());
}
Matrix Matrix::column_vector(const std::vector<real_t>& vec) {
    return Matrix(vec.size(), 1, vec.data());
}

Matrix operator*(const real_t scalar, const Matrix& M) {
    return M * scalar;
}