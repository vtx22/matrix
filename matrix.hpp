#ifndef _V_MATRIX_HPP
#define _V_MATRIX_HPP

#include <cstdint>
#include <array>
#include <iostream>
#include <iomanip>

template <uint8_t ROWS, uint8_t COLS>
class Matrix
{
public:
    Matrix(const std::array<float, ROWS * COLS> &data) : _data(data)
    {
        static_assert(ROWS > 0 && COLS > 0, "Matrix dimensions must be greater than 0!");
    }

    Matrix()
    {
        static_assert(ROWS > 0 && COLS > 0, "Matrix dimensions must be greater than 0!");
        _data.fill(0.0f);
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

    template <uint8_t R, uint8_t C>
    Matrix<ROWS, C> operator*(const Matrix<R, C> &rhs)
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

    template <uint8_t R, uint8_t C>
    Matrix operator-(const Matrix<R, C> &rhs)
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
    Matrix operator+(const Matrix<R, C> &rhs)
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

    friend std::ostream &operator<<(std::ostream &out, const Matrix &m)
    {
        std::streamsize ss = std::cout.precision();
        out << "\n";
        for (uint8_t i = 0; i < m.rows(); i++)
        {
            out << "( ";
            for (uint8_t j = 0; j < m.cols(); j++)
            {
                out << m[i * m.rows() + j] << " ";
            }
            out << ")\n";
        }

        out << std::setprecision(ss);
        return out;
    }

private:
    std::array<float, ROWS * COLS> _data;
};

typedef Matrix<2, 2> Matrix2;
typedef Matrix<3, 3> Matrix3;
typedef Matrix<4, 4> Matrix4;

typedef Matrix<2, 1> Vector2;
typedef Matrix<3, 1> Vector3;
typedef Matrix<4, 1> Vector4;

#endif // _V_MATRIX_HPP