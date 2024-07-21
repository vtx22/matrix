#ifndef _V_MATRIX_HPP
#define _V_MATRIX_HPP

#define _V_MATRIX_PRINT true

#include <cstdint>
#include <array>
#include <algorithm>

#if _V_MATRIX_PRINT
#include <iostream>
#include <iomanip>
#endif

enum MATRIX_TYPE
{
    IDENTITY,
    ZEROES,
    ONES,
} typedef MATRIX_TYPE;

template <uint8_t ROWS, uint8_t COLS>
class Matrix
{
public:
    Matrix(const std::array<float, ROWS * COLS> &data) : _data(data)
    {
        // static_assert(ROWS > 0 && COLS > 0, "Matrix dimensions must be greater than 0!");
    }

    Matrix()
    {
        // static_assert(ROWS > 0 && COLS > 0, "Matrix dimensions must be greater than 0!");
        _data.fill(0.0f);
    }

    Matrix(MATRIX_TYPE type)
    {
        // static_assert(ROWS > 0 && COLS > 0, "Matrix dimensions must be greater than 0!");

        switch (type)
        {
        case IDENTITY:
        {
            static_assert(ROWS == COLS, "Matrix must be square to be initialized as Identity Matrix!");

            _data.fill(0.0f);

            uint8_t col = 0;
            int16_t last_row = -1;
            for (uint16_t i = 0; i < ROWS * COLS; i++)
            {
                if (i % COLS == col && i / COLS > last_row)
                {
                    _data[i] = 1.0f;
                    col++;
                    last_row = i / COLS;
                }
            }
        }
        break;
        case ZEROES:
            _data.fill(0.0f);
            break;
        case ONES:
            _data.fill(1.0f);
            break;
        }
    }

    /*
    Calculate the trace of the matrix
    @return Trace value of matrix
    */
    float trace()
    {
        static_assert(COLS == ROWS, "Matrix must be square to calculate trace!");

        float trace = 0.0f;

        for (uint8_t i = 0; i < ROWS; i++)
        {
            trace += _data[i * COLS + i];
        }

        return trace;
    }

    /*
    Calculate the matrix determinant
    @return Value of the determinant
    */
    float det()
    {
        static_assert(ROWS == COLS, "Matrix must be square to calculate determinant!");

        return _determinant(*this);
    }

    /*
    Get the transposed matrix A^T
    @return Transposed matrix
    */
    Matrix<COLS, ROWS> transposed()
    {
        Matrix<COLS, ROWS> transposed;

        if (ROWS <= COLS)
        {
            for (uint8_t r = 0; r < ROWS; r++)
            {
                transposed.set_col(r, this->get_row(r));
            }
        }
        else
        {
            for (uint8_t c = 0; c < COLS; c++)
            {
                transposed.set_row(c, this->get_col(c));
            }
        }
        return transposed;
    }

    /*
    Get the inverse matrix A^(-1)
    @return Inverse matrix
    */
    Matrix inverse()
    {
        float d = det();

        // Inverse does not exist
        if (abs(d) < 1e-9)
        {
            return *this;
        }

        return adjugate() / d;
    }

    /*
    Get the adjugate matrix adj(A)
    @return Adjugate matrix
    */
    Matrix adjugate()
    {
        static_assert(COLS == ROWS, "Matrix must be square to calculate adjugate!");

        Matrix adj_t;

        for (uint16_t e = 0; e < ROWS * COLS; e++)
        {
            uint8_t row = e / COLS;
            uint8_t col = e % COLS;

            float d = submatrix(row, col).det();

            if ((row + col) % 2 != 0)
            {
                d *= -1.0f;
            }

            adj_t.set(e, d);
        }

        return adj_t.transposed();
    }

    /*
    Get the submatrix that does not contain the given row and column
    @param row Row to delete
    @param col Column to delete
    @return Reduced matrix
    */
    Matrix<ROWS - 1, COLS - 1> submatrix(uint8_t row, uint8_t col)
    {
        Matrix<ROWS - 1, COLS - 1> submatrix;

        uint16_t submatrix_index = 0;
        for (uint16_t e = 0; e < ROWS * COLS; e++)
        {
            // Current element is not part of submatrix
            if ((e % COLS) == col || (e / COLS) == row)
            {
                continue;
            }

            submatrix.set(submatrix_index++, _data[e]);
        }

        return submatrix;
    }

    /*
    Get the number of rows of the matrix
    @return The number of rows of the matrix
    */
    uint8_t rows() const
    {
        return ROWS;
    }

    /*
    Get the number of columns of the matrix
    @return The number of columns of the matrix
    */
    uint8_t cols() const
    {
        return COLS;
    }

    /*
    Set one element of the matrix
    @param index Index of the element
    @param value Value to set the element to
    */
    void set(uint16_t index, float value)
    {
        _data.at(index) = value;
    }

    /*
    Set one element of the matrix
    @param row Row of the element (starting at 0)
    @param col Column of the element (starting at 0)
    @param value Value to set the element to
    */
    void set(uint8_t row, uint8_t col, float value)
    {
        _data.at(row * COLS + col) = value;
    }

    /*
    Get the value of one element
    @param index Index of the element
    @return The value of the element
    */
    float get(uint16_t index) const
    {
        return _data.at(index);
    }

    /*
    Get the value of one element
    @param row Row of the element (starting at 0)
    @param col Column of the element (starting at 0)
    @return The value of the element
    */
    float get(uint8_t row, uint8_t col) const
    {
        return _data.at(row * COLS + col);
    }

    /*
    Get one row of the matrix as an array
    @param row Row number (starting at 0)
    @return Values of the row as std::array
    */
    std::array<float, COLS> get_row(uint8_t row) const
    {
        std::array<float, COLS> ra;
        uint16_t start = row * COLS;
        std::copy(_data.begin() + start, _data.begin() + start + COLS, ra.begin());
        return ra;
    }

    /*
    Get one column of the matrix as an array
    @param col Column number (starting at 0)
    @return Values of the column as std::array
    */
    std::array<float, ROWS> get_col(uint8_t col) const
    {
        std::array<float, ROWS> ca;

        for (uint8_t r = 0; r < ROWS; r++)
        {
            ca[r] = _data[r * COLS + col];
        }

        return ca;
    }

    /*
    Set the elements of one row
    @param row Row number (starting at 0)
    @param elements Row element array
    */
    void set_row(uint8_t row, std::array<float, COLS> elements)
    {
        uint16_t start = row * COLS;
        std::copy(elements.begin(), elements.end(), _data.begin() + start);
    }

    /*
    Set the elements of one column
    @param col Column number (starting at 0)
    @param elements Column element array
    */
    void set_col(uint8_t col, std::array<float, ROWS> elements)
    {
        for (uint8_t r = 0; r < ROWS; r++)
        {
            _data[r * COLS + col] = elements[r];
        }
    }

    template <uint8_t R, uint8_t C>
    Matrix<ROWS, C> operator*(const Matrix<R, C> &rhs) const
    {
        static_assert(COLS == R, "Matrix dimension mismatch at Matrix multiplication!");

        Matrix<ROWS, C> result;

        for (uint16_t res_i = 0; res_i < ROWS * C; res_i++)
        {
            float res_i_element = 0.0f;
            int row = res_i / C;
            int col = res_i % C;

            for (uint8_t rc_i = 0; rc_i < COLS; rc_i++)
            {
                res_i_element += _data[row * COLS + rc_i] * rhs[rc_i * C + col];
            }

            result.set(res_i, res_i_element);
        }

        return result;
    }

    Matrix operator*(float scalar) const
    {
        Matrix result = *this;
        return result *= scalar;
    }

    Matrix &operator*=(float scalar)
    {
        for (uint16_t i = 0; i < ROWS * COLS; i++)
        {
            _data[i] *= scalar;
        }

        return *this;
    }

    Matrix operator/(float scalar) const
    {
        Matrix result = *this;
        return result /= scalar;
    }

    Matrix &operator/=(float scalar)
    {
        for (uint16_t i = 0; i < ROWS * COLS; i++)
        {
            _data[i] /= scalar;
        }

        return *this;
    }

    template <uint8_t R, uint8_t C>
    Matrix operator-(const Matrix<R, C> &rhs) const
    {
        static_assert(ROWS == R && COLS == C, "Matrix dimension mismatch at Matrix addition!");

        Matrix<ROWS, COLS> result = *this;

        result -= rhs;

        return result;
    }

    template <uint8_t R, uint8_t C>
    Matrix &operator-=(const Matrix<R, C> &rhs)
    {
        static_assert(ROWS == R && COLS == C, "Matrix dimension mismatch at Matrix addition!");

        for (uint16_t i = 0; i < ROWS * COLS; i++)
        {
            _data[i] -= rhs[i];
        }

        return *this;
    }

    template <uint8_t R, uint8_t C>
    Matrix operator+(const Matrix<R, C> &rhs) const
    {
        static_assert(ROWS == R && COLS == C, "Matrix dimension mismatch at Matrix addition!");

        Matrix<ROWS, COLS> result = *this;

        result += rhs;

        return result;
    }

    template <uint8_t R, uint8_t C>
    Matrix &operator+=(const Matrix<R, C> &rhs)
    {
        static_assert(ROWS == R && COLS == C, "Matrix dimension mismatch at Matrix addition!");

        for (uint16_t i = 0; i < ROWS * COLS; i++)
        {
            _data[i] += rhs[i];
        }

        return *this;
    }

    float operator[](uint16_t i) const
    {
        return _data.at(i);
    }

#if _V_MATRIX_PRINT
    friend std::ostream &operator<<(std::ostream &out, const Matrix &m)
    {
        std::ios old_state(nullptr);
        old_state.copyfmt(out);

        out << std::setw(10) << std::fixed << std::showpos << std::setprecision(3);

        out << "\n";
        for (uint8_t i = 0; i < m.rows(); i++)
        {
            out << "( ";
            for (uint8_t j = 0; j < m.cols(); j++)
            {
                out << m[i * m.cols() + j] << " ";
            }
            out << ")\n";
        }

        out.copyfmt(old_state);
        return out;
    }

    /*
    Print matrix to std::cout
    */
    void print()
    {
        std::cout << *this;
    }

#endif

private:
    std::array<float, ROWS * COLS> _data;

    float _determinant(const Matrix<1, 1> &matrix)
    {
        return matrix[0];
    }

    float _determinant(const Matrix<2, 2> &matrix)
    {
        return matrix[0] * matrix[3] - matrix[1] * matrix[2];
    }

    template <uint8_t N>
    float _determinant(Matrix<N, N> matrix)
    {
        float det = 0.0f;

        int8_t sign = 1;

        // Develop submatrix for all rows of column 0
        for (uint8_t r = 0; r < N; r++)
        {
            Matrix<N - 1, N - 1> sm = matrix.submatrix(r, 0);

            det += sign * matrix[r * N] * _determinant(sm);
            sign *= -1;
        }

        return det;
    }
};

typedef Matrix<2, 2> Matrix2;
typedef Matrix<3, 3> Matrix3;
typedef Matrix<4, 4> Matrix4;

typedef Matrix<2, 1> Vector2;
typedef Matrix<3, 1> Vector3;
typedef Matrix<4, 1> Vector4;

#endif // _V_MATRIX_HPP