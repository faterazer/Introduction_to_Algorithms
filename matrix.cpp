#include "matrix.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
using std::cout;
using std::endl;
using std::max;
using std::vector;

matrix_int::matrix_int(const matrix_int & rhs)
	: rows(rhs.rows), columns(rhs.columns)
{
	room = new int*[rows];
	for (int i = 0; i != rows; ++i) {
		*(room + i) = new int[columns];
		for (int j = 0; j != columns; ++j)
			room[i][j] = rhs.room[i][j];
	}
}

matrix_int::matrix_int(int r, int c)
{
	room = new int*[r];
	for (int i = 0; i != r; ++i)
		*(room + i) = new int[c]();

	rows = r;
	columns = c;
}

matrix_int::matrix_int(initializer_list<initializer_list<int>> il)
{
	rows = il.size();
	columns = il.begin()->size();

	room = new int*[rows];
	for (int i = 0; i != rows; ++i)
		*(room + i) = new int[columns];

	for (int i = 0; i != rows; ++i)
		for (int j = 0; j != columns; ++j)
			room[i][j] = *((il.begin() + i)->begin() + j);
}

matrix_int::~matrix_int()
{
	for (int i = 0; i != rows; ++i)
		delete[] *(room + i);
	delete[] room;
	rows = 0;
	columns = 0;
}

matrix_int& matrix_int::operator=(const matrix_int &rhs)
{
	int **newp = new int*[rhs.rows];
	for (int i = 0; i != rhs.rows; ++i) {
		*(newp + i) = new int[rhs.columns];
		for (int j = 0; j != rhs.columns; ++j)
			newp[i][j] = rhs.room[i][j];
	}

	for (int i = 0; i != rows; ++i)
		delete[] *(room + i);
	delete[] room;
	rows = 0;
	columns = 0;

	room = newp;
	rows = rhs.rows;
	columns = rhs.columns;

	return *this;
}

matrix_int & matrix_int::operator+=(const matrix_int & rhs)
{
	for (int i = 0; i != rows; ++i)
		for (int j = 0; j != columns; ++j)
			room[i][j] += rhs[i][j];
	return *this;
}

matrix_int & matrix_int::operator-=(const matrix_int & rhs)
{
	for (int i = 0; i != rows; ++i)
		for (int j = 0; j != columns; ++j)
			room[i][j] += rhs[i][j];
	return *this;
}

void matrix_int::debug_print() const
{
	cout << endl;
	for (int i = 0; i != rows; ++i) {
		for (int j = 0; j != columns; ++j)
			std::cout << room[i][j] << " ";
		cout << endl;
	}
	cout << endl;
}

matrix_int operator+(const matrix_int & left, const matrix_int & right)
{
	matrix_int result = left;
	result += right;
	return result;
}

matrix_int operator-(const matrix_int & left, const matrix_int & right)
{
	matrix_int result = left;
	result -= right;
	return result;
}

matrix_int bf_matrix_multiply(const matrix_int & A, const matrix_int & B)
{
	matrix_int result(A.Rows(), B.Columns());
	for (int i = 0; i != A.Rows(); ++i)
		for (int j = 0; j != B.Columns(); ++j)
			for (int k = 0; k != A.Columns(); ++k)
				result[i][j] += A[i][k] * B[k][j];
	return result;
}

/*
参数 *_row，*col 是子方阵的第一个元素的下标，参数 n 代表这个子方阵的阶数，方阵 A 与
方阵 B 的阶数应相同。

 Paraments *_row, *_col is a square submatrix's first entry's index, and parament 
 n means this square submatrix's rows(cols). A and B should has equal order.
 */
void recursive_multiply(const matrix_int &A, const matrix_int &B, matrix_int &result,
						int A_row, int A_col, int B_row, int B_col, int n)
{
	if (n == 1) {
		// 当子矩阵只有一个元素时，直接相乘并返回结果。
		// 如果元素的位置在矩阵 A 和 B 的范围外，则其值为 0。从而将非方矩阵虚拟成方阵。
		// While submaxtrix has only one entry, we multiply directly and return.
		// We can virtualize a non-square matrix as a square matrix through assume
		// all outrange entries' value is zero. 
		int A_vaule = A_row < A.Rows() && A_col < A.Columns() ? A[A_row][A_col] : 0;
		int B_vaule = B_row < B.Rows() && B_col < B.Columns() ? B[B_row][B_col] : 0;

		if (A_row < A.Rows() && B_col < B.Columns())
			result[A_row][B_col] += A_vaule * B_vaule;
		return;
	}
	else {
		// We will divide n * n matrices into four n/2 * n/2 matrices. 
		// result11 = A11 * B11 + A12 * B21 recursively.
		recursive_multiply(A, B, result, A_row, A_col, B_row, B_col, n / 2);
		recursive_multiply(A, B, result, A_row, A_col + n / 2, B_row + n / 2, B_col,
			n / 2);
		// result12 = A11 * B12 + A12 * B22 recursively.
		recursive_multiply(A, B, result, A_row, A_col, B_row, B_col + n / 2, n / 2);
		recursive_multiply(A, B, result, A_row, A_col + n / 2, B_row + n / 2,
			B_col + n / 2, n / 2);
		// result21 = A21 * B11 + A22 * B21 recursively.
		recursive_multiply(A, B, result, A_row + n / 2, A_col, B_row, B_col, n / 2);
		recursive_multiply(A, B, result, A_row + n / 2, A_col + n / 2,
			B_row + n / 2, B_col, n / 2);
		// result22 = A21 * B12 + A22 * B22 recursively.
		recursive_multiply(A, B, result, A_row + n / 2, A_col, B_row, B_col + n / 2,
			n / 2);
		recursive_multiply(A, B, result, A_row + n / 2, A_col + n / 2,
			B_row + n / 2, B_col + n / 2, n / 2);
	}
}

matrix_int recursive_matrix_multiply(const matrix_int & A, const matrix_int & B)
{
	// 如果运算的矩阵不是方阵，我们可以将其虚拟为方阵。我们可以比较矩阵 A 和 B 的行列
	// 数，得到最大值 m，然后计算得到 n = 2^m，n 应不小于 m，然后我们将矩阵 A 和 B 虚
	// 拟为 n 阶方阵，超出 A 或 B 范围的元素视为 0。
	// Some works for pretreatment:
	// We can compare matrix A's rows, columns and B's rows, columns, then we get
	// the max one.
	// So we can math a number: n = 2^m and numbwe n large or equal than the max
	// number we got in last step.
	// Now we can assume both of matrix A and B are n-dimensional matrices.
	// All entries out range of A or B assume value 0.

	matrix_int result(A.Rows(), B.Columns());
	int m1 = max(A.Rows(), A.Columns()), m2 = max(B.Rows(), B.Columns());
	int n = max(m1, m2);
	int i = 0;
	while (n > pow(2, i))
		++i;
	n = static_cast<int>(pow(2, i));

	recursive_multiply(A, B, result, 0, 0, 0, 0, n);

	return result;
}

// 子矩阵减法。
// Get two square submatrices' difference.
matrix_int submatrices_sub(const matrix_int &A, const matrix_int &B, int A_row,
						   int A_col, int B_row, int B_col, int n) 
{
	matrix_int result(n, n);
	for (int i = 0; i != n; ++i) {
		for (int j = 0; j != n; ++j) {
			int a = A_row + i < A.Rows() && A_col + j < A.Columns() ? A[A_row + i][A_col + j] : 0;
			int b = B_row + i < B.Rows() && B_col + j < B.Columns() ? B[B_row + i][B_col + j] : 0;
			result[i][j] += a - b;
		}
	}
	return result;
}

// 子矩阵加法。
// Get two square submatrices' sum.
matrix_int submatrices_add(const matrix_int &A, const matrix_int &B, int A_row,
						   int A_col, int B_row, int B_col, int n) {
	matrix_int result(n, n);
	for (int i = 0; i != n; ++i) {
		for (int j = 0; j != n; ++j) {
			int a = A_row + i < A.Rows() && A_col + j < A.Columns() ? A[A_row + i][A_col + j] : 0;
			int b = B_row + i < B.Rows() && B_col + j < B.Columns() ? B[B_row + i][B_col + j] : 0;
			result[i][j] += a + b;
		}
	}
	return result;
}

// 子矩阵递归乘法。
// Multiply recursively.
void submatrices_multiply(const matrix_int &A, const matrix_int &B, matrix_int &result, 
						  int A_row, int A_col, int B_row, int B_col, int base_row, int base_col, int n) 
{
	if (n == 1) {
		int A_value = A_row < A.Rows() && A_col < A.Columns() ? A[A_row][A_col] : 0;
		int B_value = B_row < B.Rows() && B_col < B.Columns() ? B[B_row][B_col] : 0;
		result[A_row - base_row][B_col - base_col] += A_value * B_value;
		return;
	}
	else {
		submatrices_multiply(A, B, result, A_row, A_col, B_row, B_col, base_row, base_col, n / 2);
		submatrices_multiply(A, B, result, A_row, A_col + n / 2, B_row + n / 2, B_col, base_row, base_col, n / 2);
		submatrices_multiply(A, B, result, A_row, A_col, B_row, B_col + n / 2, base_row, base_col, n / 2);
		submatrices_multiply(A, B, result, A_row, A_col + n / 2, B_row + n / 2, B_col + n / 2, base_row, base_col, n / 2);
		submatrices_multiply(A, B, result, A_row + n / 2, A_col, B_row, B_col, base_row, base_col, n / 2);
		submatrices_multiply(A, B, result, A_row + n / 2, A_col + n / 2, B_row + n / 2, B_col, base_row, base_col, n / 2);
		submatrices_multiply(A, B, result, A_row + n / 2, A_col, B_row, B_col + n / 2, base_row, base_col, n / 2);
		submatrices_multiply(A, B, result, A_row + n / 2, A_col + n / 2, B_row + n / 2, B_col + n / 2, base_row, base_col, n / 2);
	}
}

/*
Strassen 算法包含 4 个步骤：
1. 将矩阵 A，B 和输出矩阵 C，分解为 n/2 * n/2 的子矩阵。采用下标计算方法，此步骤
花费 Θ(1) 时间。
2. 创建 10 个 n/2 * n/2 的矩阵 S1, S2, ..., S10，每个矩阵保存步骤 1 中创建的两个
子矩阵的和或差。 花费时间为 Θ(n^2)。
3. 用步骤 1 中创建的子矩阵和步骤 2 中创建的 10 个矩阵，递归地计算 7 个矩阵积 P1，
P2，...，P7。每个矩阵 Pi 都是 n/2 * n/2 的。
4. 通过 Pi 矩阵的不同组合进行加减运算，计算出结果矩阵 C 的子矩阵 C11，C12，C21，C22。
花费时间 Θ(n^2)。

在步骤 2 中，创建如下 10 个矩阵：
S1 = B12 - B22
S2 = A11 + A12
S3 = A21 + A22
S4 = B21 - B11
S5 = A11 + A22
S6 = B11 + B22
S7 = A12 - A22
S8 = B21 + B22
S9 = A11 - A21
S10 = B11 + B12

在步骤 3 中，递归地计算 7 次 n/2 * n/2 矩阵的乘法，如下所示：
P1 = A11 * S1	(A11 * B12 - A11 * B22)
P2 = S2 * B22	(A11 * B22 + A12 * B22)
P3 = S3 * B11	(A21 * B11 + A22 * B11)
P4 = A22 * S4	(A22 * B21 - A22 * B11)
P5 = S5	* S6	(A11 * B11 + A11 * B22 + A22 * B11 + A22 * B22)
P6 = S7 * S8	(A12 * B21 + A12 * B22 - A22 * B21 - A22 * B22)
P7 = S9 * S10	(A11 * B11 + A11 * B12 - A21 * B11 - A21 * B12)

步骤 4 对步骤 3 创建的 Pi 矩阵进行加减法运算，计算出 C 的 4 个 n/2 * n/2 的子矩阵：
C11 = P5 + P4 - P2 + P6
C12 = P1 + P2
C21 = P3 + P4
C22 = P5 + P1 - P3 - P7

It has four steps:
1. Divide the input matrices A and B and output matrix C into n/2 * n/2
submatrices, This step takes Θ(1) time by index calculation.
2. Create 10 matrices S1, S2 ... S10, each of which is n/2 * n/2 and is the
sum or difference of two matrices created in step 1. We can create all 10
matrices in Θ(n^2) time.
3. Using the submatrices created in step 1 and the 10 matrices created in
step 2, recursively compute seven matrix products P1, P2 ... P7. Each
matrix Pi is n/2 * n/2.
4. Compute the desired submatrices C11, C12, C21, C22 of the result matrix C
by adding and subtracting various combinations of the Pi matrices. We can
compute all four submatrices in Θ(n^2) time.

In step2, we create the following 10 matrices:
S1 = B12 - B22
S2 = A11 + A12
S3 = A21 + A2
S4 = B21 - B11
S5 = A11 + A22
S6 = B11 + B22
S7 = A12 - A22
S8 = B21 + B22
S9 = A11 - A21
S10 = B11 + B12

In step3, we recursively multiply n/2 * n/2 matrices seven times to compute the
following n/2 * n/2 matrices, each of which is the sum or difference of products
of A and B submatrices:
P1 = A11 * S1	(A11 * B12 - A11 * B22)
P2 = S2 * B22	(A11 * B22 + A12 * B22)
P3 = S3 * B11	(A21 * B11 + A22 * B11)
P4 = A22 * S4	(A22 * B21 - A22 * B11)
P5 = S5	* S6	(A11 * B11 + A11 * B22 + A22 * B11 + A22 * B22)
P6 = S7 * S8	(A12 * B21 + A12 * B22 - A22 * B21 - A22 * B22)
P7 = S9 * S10	(A11 * B11 + A11 * B12 - A21 * B11 - A21 * B12)

Step4 adds and subtracts the Pi matrices created in step 3 to construct the four
n/2 * n/2 submatrices of the product C:
C11 = P5 + P4 - P2 + P6
C12 = P1 + P2
C21 = P3 + P4
C22 = P5 + P1 - P3 - P7
*/
matrix_int strassen_matrix_multiply(const matrix_int & A, const matrix_int & B)
{
	matrix_int result(A.Rows(), B.Columns());
	int m1 = max(A.Rows(), A.Columns()), m2 = max(B.Rows(), B.Columns());
	int n = max(m1, m2);
	int i = 1;
	while (n > pow(2, i)) ++i;
	n = static_cast<int>(pow(2, i));

	vector<matrix_int> S(10);
	S[0] = submatrices_sub(B, B, 0, n / 2, n / 2, n / 2, n / 2);
	S[1] = submatrices_add(A, A, 0, 0, 0, n / 2, n / 2);
	S[2] = submatrices_add(A, A, n / 2, 0, n / 2, n / 2, n / 2);
	S[3] = submatrices_sub(B, B, n / 2, 0, 0, 0, n / 2);
	S[4] = submatrices_add(A, A, 0, 0, n / 2, n / 2, n / 2);
	S[5] = submatrices_add(B, B, 0, 0, n / 2, n / 2, n / 2);
	S[6] = submatrices_sub(A, A, 0, n / 2, n / 2, n / 2, n / 2);
	S[7] = submatrices_add(B, B, n / 2, 0, n / 2, n / 2, n / 2);
	S[8] = submatrices_sub(A, A, 0, 0, n / 2, 0, n / 2);
	S[9] = submatrices_add(B, B, 0, 0, 0, n / 2, n / 2);

	vector<matrix_int> P(7, matrix_int(n / 2, n / 2));
	submatrices_multiply(A, S[0], P[0], 0, 0, 0, 0, 0, 0, n / 2);
	submatrices_multiply(S[1], B, P[1], 0, 0, n / 2, n / 2, 0, n / 2, n / 2);
	submatrices_multiply(S[2], B, P[2], 0, 0, 0, 0, 0, 0, n / 2);
	submatrices_multiply(A, S[3], P[3], n / 2, n / 2, 0, 0, n / 2, 0, n / 2);
	submatrices_multiply(S[4], S[5], P[4], 0, 0, 0, 0, 0, 0, n / 2);
	submatrices_multiply(S[6], S[7], P[5], 0, 0, 0, 0, 0, 0, n / 2);
	submatrices_multiply(S[8], S[9], P[6], 0, 0, 0, 0, 0, 0, n / 2);

	matrix_int result11 = P[4] + P[3] - P[1] + P[5];
	matrix_int result12 = P[0] + P[1];
	matrix_int result21 = P[2] + P[3];
	matrix_int result22 = P[4] + P[0] - P[2] - P[6];

	for (int i = 0; i != result.Rows(); ++i) {
		for (int j = 0; j != result.Columns(); ++j) {
			if (i < n / 2 && j < n / 2)
				result[i][j] = result11[i][j];
			else if (i < n / 2 && j >= n / 2)
				result[i][j] = result12[i][j - n / 2];
			else if (i >= n / 2 && j < n / 2)
				result[i][j] = result21[i - n / 2][j];
			else
				result[i][j] = result22[i - n / 2][j - n / 2];
		}
	}
	return result;
}
