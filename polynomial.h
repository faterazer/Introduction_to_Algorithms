// 此文件包含了多种多项式求值算法。
// This head file includes various algorithms for evaluating a polynomial.

#ifndef _POLYNOMIAL_H
#define _POLYNOMIAL_H

#include <vector>
using std::vector;

// 朴素多项式求值算法，时间复杂度为：Θ(n^2)。
// Native polynomial evaluation algorithm，runnning time is Θ(n^2).
int native_polynomial_evaluation(const vector<int> &coefficients, int x)
{
	int result = 0;
	for (int i = 0; i != coefficients.size(); ++i) {
		int x_i = 1;
		for (int j = 0; j != i; ++j) {
			x_i *= x;
		}
		result += coefficients[i] * x_i;
	}
	return result;
}

// 霍纳 - 秦九韶算法，时间复杂度为 Θ(n)。
// Horner-Qin Jiushao algorithm, running time is Θ(n).
int horner_polynomial_evaluation(const vector<int> &coefficients, int x)
{
	int result = 0;
	for (int i = coefficients.size() - 1; i >= 0; --i)
		result = coefficients[i] + result * x;
	return result;
}

#endif // !_POLYNOMIAL_H
