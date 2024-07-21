# matrix
Minimal C++ Matrix library with fully static matrix dimensions known at compile time
# Features
 - Matrix initialization
 - Special Matricies (Identity, Ones, Zeroes)
 - Operators: Addition, Subtraction, Scalar Multiplication/Division, Matrix Multiplication
 - Element manipulation
 - Submatrix creation (row/column elimination)
 - Trace
 - Determinant
 - Transposed Matrix
 - Adjugate Matrix
 - Inverse Matrix
 - Matrix printing to console

# Usage
## Create Matrix
```C++
Matrix<3, 3> m;                               // Create a 3x3 matrix filled with zeroes per default
Matrix<3, 3> m({0, 1, 2, 3, 4, 5, 6, 7, 8});  // Create 3x3 matrix from array

// Sepcial types
Matrix<3, 3> m(ZEROES);   // 3x3 matrix explicetely filled with zeroes
Matrix<3, 3> m(ONES);     // 3x3 matrix explicetely filled with ones
Matrix<3, 3> m(IDENTITY); // 3x3 identity matrix (ones along diagonal)
```
Special predefined matricies:
```C++
typedef Matrix<2, 2> Matrix2;
typedef Matrix<3, 3> Matrix3;
typedef Matrix<4, 4> Matrix4;

const Matrix<2, 2> MatrixI2(IDENTITY);
const Matrix<3, 3> MatrixI3(IDENTITY);
const Matrix<4, 4> MatrixI4(IDENTITY);
```

## Set/Get elements
```C++
m.get(1, 2);   // Get the matrix element at row 1 column 2 (count starts at 0)
m.get(8);      // Get the matrix element at absolute index 8 (count starts at 0 and goes rowwise)
m[8];          // Get the matrix element at absolute index 8 (count starts at 0 and goes rowwise)    

```
```C++
m.set(1, 2, 10.0f);   // Set the matrix element at row 1 column 2 to 10.0
m.set(8, 10.0f);      // Set the matrix element at absolute index 8 to 10.0
```
```C++
m.get_row(1);       // Get all elements of row 1 as array
m.get_col(1);       // Get all elements of column 1 as array
```
```C++
m.set_row(1, {1, 2, 3});       // Set all elements of row 1 from array
m.set_col(1, {1, 2, 3});       // Set all elements of column 1 from array
```
```C++
m.rows();       // Get the number of matrix rows
m.cols();       // Get the number of matrix columns
```
## Submatrix
Get the submatrix by eliminating a row and a column
```C++
Matrix<2, 2> sm = m.submatrix(1, 2);  // Delete row 1 and column 2 from matrix to get submatrix
```
## Matrix Algebra
```C++
Matrix<2, 2> m1, m2;
```
```C++
Matrix<2, 2> sum = m1 + m2;      // Matrix elementwise addition
Matrix<2, 2> dif = m1 - m2;      // Matrix elementwise subtraction
Matrix<2, 2> pro = m1 * m2;      // Matrix multiplication
Matrix<2, 2> sca = m1 * 2.0f;    // (Elemetwise) multiplication by scalar
Matrix<2, 2> sca = m1 / 2.0f;    // (Elemetwise) division by scalar
```
## Matrix Tranformations
```C++
m.transposed(); // Get transposed matrix
```
```C++
m.adjugate();   // Get adjugate matrix
```
```C++
m.inverse();    // Get inverse matrix
```
## Matrix special values
```C++
m.trace(); // Get trace of the matrix 
```
```C++
m.det(); // Get determinant of the matrix 
```
## Printing
Matrix can be printed to console. `_V_MATRIX_PRINT` has to be set to `true` in `matrix.hpp`
```C++
std::cout << m << std::endl;  // Print using << operator
```
```C++
m.print();                    // Same as above
```
