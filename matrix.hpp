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
    Matrix(const std::array<float, ROWS * COLS> &data) : _data(data) {}
    Matrix()
    {
        _data.fill(0.0f);
    }

    uint8_t rows() const
    {
        return ROWS;
    }
    uint8_t cols() const
    {
        return COLS;
    }

    void set(uint16_t index, float value)
    {
        _data[index] = value;
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