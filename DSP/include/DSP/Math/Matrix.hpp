#ifndef SJ_DSP_MATRIX_HPP
#define SJ_DSP_MATRIX_HPP

#include <string> // 
#include <vector> // 
// #include <initializer_list>
// #include <utility> // pair, move
#include "DSP/dsp_types.h" // real_t


namespace DSP {

namespace Math {

DSP_CLASS Matrix {

    protected:
        size_t rows;
        size_t columns;
        real_t* elements;

        typedef enum  {
            calc_det,
            solve_lse,
            invert_mat
        } gauss_jordan_task_t;
        static real_t gauss_jordan(Matrix* const A, Matrix* const B, const gauss_jordan_task_t task);

        // 0-based
        void swap_rows(const size_t row_a_index, const size_t row_b_index);

        // 0-based
        void swap_columns(const size_t column_a_index, const size_t column_b_index);

        // 0-based
        void scale_up_row(const size_t row_index, const size_t column_a_index, const size_t column_b_index, const real_t s);

        // 0-based
        void scale_up_column(const size_t column_index, const size_t row_a_index, const size_t row_b_index, const real_t s);

        // 0-based
        void scale_down_row(const size_t row_index, const size_t column_a_index, const size_t column_b_index, const real_t s);

        // 0-based
        void scale_down_column(const size_t column_index, const size_t row_a_index, const size_t row_b_index, const real_t s);

        // 0-based (row_a += s * row_b)
        void modify_row(const size_t row_a_index, const size_t row_b_index, const size_t column_a_index, const size_t column_b_index, const real_t s);

        // 0-based (column_a += s * column_b)
        void modify_column(const size_t column_a_index, const size_t column_b_index, const size_t row_a_index, const size_t row_b_index, const real_t s);

    public:

        // Constructor
        Matrix();
        Matrix(const size_t rows, const size_t columns);
        Matrix(const size_t rows, const size_t columns, const real_t* elements);
        Matrix(const size_t rows, const size_t columns, const std::vector<real_t>& elements);
        Matrix(const std::vector<real_t>& row_vec);

        // Initilize
        Matrix(const std::initializer_list<real_t>& row_il);
        void operator=(const std::initializer_list<real_t>& row_il);

        // Copy
        Matrix(const Matrix& other);
        void operator=(const Matrix& other);

        // Move
        Matrix(Matrix&& other);
        void operator=(Matrix&& other);

        // Destructor
        void release_internal_array();
        virtual ~Matrix();

        // to bool
        operator bool() const;

        // to std::vector
        operator std::vector<real_t>() const;

        // to string
        std::string to_string() const;
        operator std::string() const;

        // reshape
        Matrix reshape(const size_t rows, const size_t columns) const;
        Matrix& reshape(const size_t rows, const size_t columns);


        // Element access 1-based
        real_t& operator()(const size_t row, const size_t column);
        const real_t& operator()(const size_t row, const size_t column) const;

        // Element access 0-based
        real_t& m(const size_t row_index, const size_t column_index);
        const real_t& m(const size_t row_index, const size_t column_index) const;

        // Element access 1-based
        real_t& at(const size_t row, const size_t column);
        const real_t& at(const size_t row, const size_t column) const;

        // Element access
        real_t& operator[](const size_t index);
        const real_t& operator[](const size_t index) const;


        // Getter
        size_t get_num_rows() const;
        size_t get_num_columns() const;
        size_t get_size() const;
        real_t* get_elements();
        const real_t* get_elements() const;

        // Getter 0-based
        Matrix row_vector(const size_t row_index) const;
        Matrix column_vector(const size_t column_index) const;
        Matrix submatrix(const size_t row_index, const size_t column_index) const;


        // Arithmetic
        Matrix operator+() const;
        Matrix operator-() const;

        // Arithmetic
        Matrix operator+(const Matrix& other) const;
        Matrix operator-(const Matrix& other) const;
        Matrix operator*(const Matrix& other) const;
        Matrix operator/(const Matrix& other) const; // (A / B) == (inv(B) * A)

        // Arithmetic
        void operator+=(const Matrix& other);
        void operator-=(const Matrix& other);
        void operator*=(const Matrix& other);
        void operator/=(const Matrix& other); // (A / B) == (inv(B) * A)

        // Arithmetic
        Matrix operator*(const real_t scalar) const;
        Matrix operator/(const real_t scalar) const;
        void operator*=(const real_t scalar);
        void operator/=(const real_t scalar);


        // Determinante
        real_t det() const;

        // Inverse 
        // row == columns -> Gauß, 
        // row > column -> ~M = ~(C^T * C) * C^T
        // row < column -> ~M = C^T * ~(C * C^T)
        Matrix operator~() const;
        Matrix inv() const; // 
        Matrix pinv() const; //

        // Transpose
        Matrix transpose() const;

        // Diagonal
        std::vector<real_t> diag() const;

        // solve (A * x = b)
        Matrix solve(const Matrix& b) const;

        // Static
        static real_t det(Matrix M);
        static Matrix inv(Matrix M);
        static Matrix pinv(const Matrix& M);
        static Matrix transpose(const Matrix& M);

        // solve (A * x = b)
        static Matrix solve(Matrix A, Matrix b);

        static Matrix zeros(const size_t dim);
        static Matrix zeros(const size_t rows, const size_t columns);
        static Matrix eye(const size_t dim);
        static Matrix eye(const size_t rows, const size_t columns);
        static Matrix diag(const std::vector<real_t>& vec);
        static std::vector<real_t> diag(const Matrix& M);
        static Matrix vandermonde(const size_t order, const std::vector<real_t>& x_values);


        static Matrix reshape(const Matrix& mat, const size_t rows, const size_t columns);
        static Matrix row_vector(const std::vector<real_t>& vec);
        static Matrix column_vector(const std::vector<real_t>& vec);

}; // class Matrix

}; // Namespace Math

}; // Namespace DSP

DSP_FUNCTION DSP::Math::Matrix operator*(const real_t scalar, const DSP::Math::Matrix& M);


#endif // SJ_DSP_MATRIX_HPP