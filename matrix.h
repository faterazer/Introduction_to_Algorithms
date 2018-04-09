// 此头文件包含了各种矩阵算法。
// This head file includes various algorithms of matrix.

#ifndef _MATRIX_H
#define _MATRIX_H

#include <initializer_list>
using std::initializer_list;

class matrix_int
{
public:
	matrix_int() = default;
	matrix_int(const matrix_int &rhs);
	matrix_int(int r, int c);
	matrix_int(initializer_list<initializer_list<int>> il);
	~matrix_int();

	matrix_int& operator=(const matrix_int &rhs);
	matrix_int &operator+=(const matrix_int &rhs);
	matrix_int &operator-=(const matrix_int &rhs);
	int* operator[](int r) { return *(room + r); }
	const int* operator[](int r) const { return *(room + r); }

	int Rows() const { return rows; }
	int Columns() const { return columns; }

	void debug_print() const;
private:
	int **room = nullptr;
	int rows = 0, columns = 0;
};

matrix_int operator+(const matrix_int &left, const matrix_int &right);
matrix_int operator-(const matrix_int &left, const matrix_int &right);

// 暴力求解，时间复杂度为 θ(n^3)。
// Matrix Multiply algorithm by brute force method, running time: θ(n^3).
matrix_int bf_matrix_multiply(const matrix_int &A, const matrix_int &B);

// 矩阵相乘的递归求法，时间复杂度为 Θ(n^3)。
// Matrix Multiply algorithm by recurive method. Running time : Θ(n^3).
matrix_int recursive_matrix_multiply(const matrix_int &A, const matrix_int &B);

// Strassen 算法，时间复杂度为 Θ(n^lg7）， lg7 的值介于 2.80 和 2.81 之间。
// Strassen method. Running time: Θ(n^lg7). Since lg7 lies between 2.80 and 2.81, 
// Strassen's algorithm runs in O(n^2.81) time.
matrix_int strassen_matrix_multiply(const matrix_int &A, const matrix_int &B);

#endif // !_MATRIX_H