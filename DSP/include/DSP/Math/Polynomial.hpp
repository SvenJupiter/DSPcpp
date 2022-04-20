#ifndef SJ_DSP_POLYNOMIAL_HPP
#define SJ_DSP_POLYNOMIAL_HPP


#include <string> // 
#include <vector> // 
#include <initializer_list>
#include <utility> // pair, move
#include <iterator> // For std::forward_iterator_tag
#include <cstddef>  // For std::ptrdiff_t
#include "DSP/dsp_types.h" // real_t


namespace DSP {

namespace Math {

DSP_CLASS Polynomial {

    protected:
        mutable size_t order;
        mutable real_t* coeffs;


    public:

        // Constructor
        Polynomial();

        // Create a polynomial of size 'order' but don't initilize its coeffs
        Polynomial(const size_t order);

        // Create a polynomial for a givan array (if 'coeffs' is NULL, set all coeffs to 0)
        Polynomial(const size_t order, const real_t* const coeffs);

        // Create a polynomial from a give vector
        Polynomial(const std::vector<real_t>& coeffs);

        // Initilize
        Polynomial(const std::initializer_list<real_t>& il);
        void operator=(const std::initializer_list<real_t>& il);

        // Copy
        Polynomial(const Polynomial& other);
        void operator=(const Polynomial& other);

        // Move
        Polynomial(Polynomial&& other);
        void operator=(Polynomial&& other);

        // polyval:
        // Calculate value of polynomial for 'x'
        // y = p(x) = a0 + a1 * x + a2 * x^2 + ... an * x^n
        real_t operator()(const real_t x) const;

        // Getter
        size_t get_order() const;

        // access coeffs
        real_t& operator[](const size_t index);
        const real_t& operator[](const size_t index) const;

        // access coeffs (poly.a(0), poly->a(2), ...)
        real_t& a(const size_t index);
        const real_t& a(const size_t index) const;

        // convert to readable form
        std::string to_string() const;

        // Destructor
        virtual ~Polynomial();


        // operators
        Polynomial operator+() const;
        Polynomial operator-() const;
        Polynomial operator+(const Polynomial& other) const;
        Polynomial operator-(const Polynomial& other) const;
        Polynomial operator*(const Polynomial& other) const;
        void operator+=(const Polynomial& other);
        void operator-=(const Polynomial& other);
        void operator*=(const Polynomial& other);

        // scalar
        Polynomial operator+(const real_t scalar) const;
        Polynomial operator-(const real_t scalar) const;
        Polynomial operator*(const real_t scalar) const;
        Polynomial operator/(const real_t scalar) const;
        void operator+=(const real_t scalar);
        void operator-=(const real_t scalar);
        void operator*=(const real_t scalar);
        void operator/=(const real_t scalar);

        // polydiv
        std::pair<Polynomial, Polynomial> operator/(const Polynomial& other) const;

        // roots
        // std::vector<real_t> roots() const;

        // polyfit
        static Polynomial fit(const size_t order, const std::vector<real_t>& x, const std::vector<real_t>& y);
        static void fit(real_t* const p, const size_t order, const real_t* const x, const real_t* const y, const size_t size);

        // from roots
        // static Polynomial poly(const std::vector<real_t>& roots, const real_t an);

        void shrink_to_fit() const;
        void grow_to(const size_t new_order) const;

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
                iterator(Polynomial& poly, const int32_t index);
                iterator(Polynomial* const poly, const int32_t index);

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
                Polynomial* poly;
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
                const_iterator(const Polynomial& poly, const int32_t index);
                const_iterator(const Polynomial* const poly, const int32_t index);

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
                const Polynomial* poly;
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
                reverse_iterator(Polynomial& poly, const int32_t index);
                reverse_iterator(Polynomial* const poly, const int32_t index);

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
                Polynomial* poly;
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
                const_reverse_iterator(const Polynomial& poly, const int32_t index);
                const_reverse_iterator(const Polynomial* const poly, const int32_t index);

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
                const Polynomial* poly;
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
};

}; // Namespace Math

}; // Namespace DSP

DSP_FUNCTION DSP::Math::Polynomial operator+(const real_t scalar, const DSP::Math::Polynomial& poly);
DSP_FUNCTION DSP::Math::Polynomial operator-(const real_t scalar, const DSP::Math::Polynomial& poly);
DSP_FUNCTION DSP::Math::Polynomial operator*(const real_t scalar, const DSP::Math::Polynomial& poly);

#endif // SJ_DSP_POLYNOMIAL_HPP