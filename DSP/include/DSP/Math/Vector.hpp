#ifndef SJ_DSP_VECTOR_HPP
#define SJ_DSP_VECTOR_HPP

#include <string> // 
#include <vector> // 
#include <initializer_list>
#include <utility> // pair, move
#include "DSP/dsp_types.h" // real_t
#include "DSP/Math/Matrix.hpp"

namespace DSP {

namespace Math {

DSP_CLASS Vector : protected Matrix {

    public:

        // Constructor
        Vector();
        Vector(const size_t size);
        Vector(const size_t size, const real_t* vec);
        Vector(const std::vector<real_t> vec);

        // Initilize
        Vector(const std::initializer_list<real_t>& il);
        void operator=(const std::initializer_list<real_t>& il);

        // Copy
        Vector(const Matrix& colum_mat);
        Vector(const Vector& other);
        void operator=(const Vector& other);

        // Move
        Vector(Matrix&& colum_mat);
        Vector(Vector&& other);
        void operator=(Vector&& other);

        // Destructor
        using Matrix::release_internal_array;
        virtual ~Vector();

        // to bool
        using Matrix::operator bool;

        // to std::vector
        using Matrix::operator std::vector<real_t>;

        // to string
        using Matrix::to_string;
        using Matrix::operator std::string;

        // to Matrix
        Matrix to_Matrix() const;

        // to real_t
        operator real_t() const;
        Vector abs() const;
        real_t norm_squared() const;
        real_t norm() const;
        real_t length_squared() const;
        real_t length() const;
        Vector normalized() const;
        real_t cosphi(const Vector& other) const;
        real_t angle(const Vector& other) const;
        real_t project(const Vector& other) const;
        Vector projection(const Vector& other) const;


        // element access 1-based
        real_t& operator()(const size_t num_element);
        const real_t& operator()(const size_t num_element) const;

        // Element access 0-based
        real_t& v(const size_t index);
        const real_t& v(const size_t index) const;

        // Element access 1-based
        real_t& at(const size_t num_element);
        const real_t& at(const size_t num_element) const;

        // Element access 0-based
        using Matrix::operator[];

        // Getter
        using Matrix::get_size;
        using Matrix::get_elements;
        size_t size() const { return this->get_size(); }


        // Arithmetic
        Vector operator+() const;
        Vector operator-() const;

        // Arithmetic
        Vector operator+(const Vector& other) const;
        Vector operator-(const Vector& other) const;
        real_t operator*(const Vector& other) const;

        // Arithmetic
        void operator+=(const Vector& other);
        void operator-=(const Vector& other);
        void operator*=(const Vector& other) = delete;

        // Arithmetic
        Vector operator*(const real_t scalar) const;
        Vector operator/(const real_t scalar) const;
        void operator*=(const real_t scalar);
        void operator/=(const real_t scalar);

        // We can't divide Vectors
        Vector operator/(const Vector& other) const = delete; 
        void operator/=(const Vector& other) = delete; 

        // Dot product
        real_t dot(const Vector& other) const;

        // Cross product (for R^3 vector only)
        Vector cross(const Vector& other) const;

        // Convolution and polynomial multiplication
        // w = conv(u,v) returns the convolution of vectors u and v. 
        // If u and v are vectors of polynomial coefficients, 
        // convolving them is equivalent to multiplying the two polynomials.
        Vector conv(const Vector& other) const;

        // Deconvolution and polynomial division
        // [q,r] = deconv(u,v) deconvolves a vector v out of a vector u using long division, 
        // and returns the quotient q and remainder r such that u = conv(v,q) + r. 
        // If u and v are vectors of polynomial coefficients, then deconvolving them is equivalent to 
        // dividing the polynomial represented by u by the polynomial represented by v.
        std::pair<Vector, Vector> deconv(const Vector& other) const;


        // Dot product
        static real_t dot(const Vector& u, const Vector& v);
        
        // Cross product (for R^3 vector only)
        static Vector cross(const Vector& u, const Vector& v);

        // Convolution and polynomial multiplication
        // w = conv(u,v) returns the convolution of vectors u and v. 
        // If u and v are vectors of polynomial coefficients, 
        // convolving them is equivalent to multiplying the two polynomials.
        static Vector conv(const Vector& u, const Vector& v);

        // Deconvolution and polynomial division
        // [q,r] = deconv(u,v) deconvolves a vector v out of a vector u using long division, 
        // and returns the quotient q and remainder r such that u = conv(v,q) + r. 
        // If u and v are vectors of polynomial coefficients, then deconvolving them is equivalent to 
        // dividing the polynomial represented by u by the polynomial represented by v.
        static std::pair<Vector, Vector> deconv(const Vector& u, const Vector& v);


        static Vector abs(const Vector& vec);
        static real_t norm_squared(const Vector& vec);
        static real_t norm(const Vector& vec);
        static real_t length_squared(const Vector& vec);
        static real_t length(const Vector& vec);
        static Vector normalized(const Vector& vec);
        static real_t cosphi(const Vector& a, const Vector& b);
        static real_t angle(const Vector& a, const Vector& b);
        static real_t project(const Vector& a, const Vector& b);
        static Vector projection(const Vector& a, const Vector& b);


  public:

        // ---- iterators ----

        DSP_CLASS iterator {

            public:
                using iterator_category = std::random_access_iterator_tag;
                // using difference_type   = std::ptrdiff_t;
                using value_type        = real_t;
                using pointer           = real_t*;  // or also value_type*
                using reference         = real_t&;  // or also value_type&
                using const_value_type        = const real_t;
                using const_pointer           = const real_t*;  // or also value_type*
                using const_reference         = const real_t&;  // or also value_type&

                // Constructor
                iterator(Vector& vec, const int32_t index);
                iterator(Vector* const vec, const int32_t index);

                // Comparison
                bool operator==(const iterator& other) const;
                bool operator!=(const iterator& other) const;
                bool operator<(const iterator& other) const;
                bool operator<=(const iterator& other) const;
                bool operator>(const iterator& other) const;
                bool operator>=(const iterator& other) const;

                // Access
                reference operator*();
                const_reference operator*() const;
                reference operator->();
                const_reference operator->() const;
                reference operator[](const int32_t offset);
                const_reference operator[](const int32_t offset)const;

                // Prefix increment
                iterator& operator++();
                iterator& operator--();

                // Postfix increment
                iterator operator++(int);
                iterator operator--(int);

                // Arithmetic
                iterator operator+(const int32_t n) const;
                iterator operator-(const int32_t n) const;
                iterator& operator+=(const int32_t n);
                iterator& operator-=(const int32_t n);
                int32_t operator-(const iterator& other) const;

            private:
                Vector* vec;
                int32_t index;
        };
        
        DSP_CLASS const_iterator {

            public:
                using iterator_category = std::random_access_iterator_tag;
                // using difference_type   = std::ptrdiff_t;
                using const_value_type        = const real_t;
                using const_pointer           = const real_t*;  // or also value_type*
                using const_reference         = const real_t&;  // or also value_type&

                // Constructor
                const_iterator(const Vector& vec, const int32_t index);
                const_iterator(const Vector* const vec, const int32_t index);

                // Comparison
                bool operator==(const const_iterator& other) const;
                bool operator!=(const const_iterator& other) const;
                bool operator<(const const_iterator& other) const;
                bool operator<=(const const_iterator& other) const;
                bool operator>(const const_iterator& other) const;
                bool operator>=(const const_iterator& other) const;

                // Access
                const_reference operator*() const;
                const_reference operator->() const;
                const_reference operator[](const int32_t offset)const;

                // Prefix increment
                const_iterator& operator++();
                const_iterator& operator--();

                // Postfix increment
                const_iterator operator++(int);
                const_iterator operator--(int);

                // Arithmetic
                const_iterator operator+(const int32_t n) const;
                const_iterator operator-(const int32_t n) const;
                const_iterator& operator+=(const int32_t n);
                const_iterator& operator-=(const int32_t n);
                int32_t operator-(const const_iterator& other) const;

            private:
                const Vector* vec;
                int32_t index;
        };

        DSP_CLASS reverse_iterator {

            public:
                using iterator_category = std::random_access_iterator_tag;
                // using difference_type   = std::ptrdiff_t;
                using value_type        = real_t;
                using pointer           = real_t*;  // or also value_type*
                using reference         = real_t&;  // or also value_type&

                // Constructor
                reverse_iterator(Vector& vec, const int32_t index);
                reverse_iterator(Vector* const vec, const int32_t index);

                // Comparison
                bool operator==(const reverse_iterator& other) const;
                bool operator!=(const reverse_iterator& other) const;
                bool operator<(const reverse_iterator& other) const;
                bool operator<=(const reverse_iterator& other) const;
                bool operator>(const reverse_iterator& other) const;
                bool operator>=(const reverse_iterator& other) const;

                // Access
                reference operator*();
                const reference operator*() const;
                reference operator->();
                const reference operator->() const;
                reference operator[](const int32_t offset);
                const reference operator[](const int32_t offset)const;

                // Prefix increment
                reverse_iterator& operator++();
                reverse_iterator& operator--();

                // Postfix increment
                reverse_iterator operator++(int);
                reverse_iterator operator--(int);

                // Arithmetic
                reverse_iterator operator+(const int32_t n) const;
                reverse_iterator operator-(const int32_t n) const;
                reverse_iterator& operator+=(const int32_t n);
                reverse_iterator& operator-=(const int32_t n);
                int32_t operator-(const reverse_iterator& other) const;

            private:
                Vector* vec;
                int32_t index;
        };

        DSP_CLASS const_reverse_iterator {

            public:
                using iterator_category = std::random_access_iterator_tag;
                // using difference_type   = std::ptrdiff_t;
                using const_value_type        = const real_t;
                using const_pointer           = const real_t*;  // or also value_type*
                using const_reference         = const real_t&;  // or also value_type&

                // Constructor
                const_reverse_iterator(const Vector& vec, const int32_t index);
                const_reverse_iterator(const Vector* const vec, const int32_t index);

                // Comparison
                bool operator==(const const_reverse_iterator& other) const;
                bool operator!=(const const_reverse_iterator& other) const;
                bool operator<(const const_reverse_iterator& other) const;
                bool operator<=(const const_reverse_iterator& other) const;
                bool operator>(const const_reverse_iterator& other) const;
                bool operator>=(const const_reverse_iterator& other) const;

                // Access
                const_reference operator*() const;
                const_reference operator->() const;
                const_reference operator[](const int32_t offset)const;

                // Prefix increment
                const_reverse_iterator& operator++();
                const_reverse_iterator& operator--();

                // Postfix increment
                const_reverse_iterator operator++(int);
                const_reverse_iterator operator--(int);

                // Arithmetic
                const_reverse_iterator operator+(const int32_t n) const;
                const_reverse_iterator operator-(const int32_t n) const;
                const_reverse_iterator& operator+=(const int32_t n);
                const_reverse_iterator& operator-=(const int32_t n);
                int32_t operator-(const const_reverse_iterator& other) const;

            private:
                const Vector* vec;
                int32_t index;
        };


        iterator begin();
        const_iterator begin() const;
        iterator end();
        const_iterator end() const;
        reverse_iterator rbegin();
        const_reverse_iterator rbegin() const;
        reverse_iterator rend();
        const_reverse_iterator rend() const;
        const_iterator cbegin() const;
        const_iterator cend() const;
        const_reverse_iterator crbegin() const;
        const_reverse_iterator crend() const;

}; // class Vector

}; // Namespace Math

}; // Namespace DSP


DSP_FUNCTION DSP::Math::Vector operator*(const real_t scalar, const DSP::Math::Vector& v);
DSP_FUNCTION DSP::Math::Vector operator*(const DSP::Math::Matrix& M, const DSP::Math::Vector& v);

#endif // SJ_DSP_VECTOR_HPP