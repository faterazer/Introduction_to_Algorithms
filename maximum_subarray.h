/*
此文件包含了求解最大子数组的算法，包括暴力求解、分治法求解和动态规划求解。
时间复杂度：暴力求解为 Θ(n^2)，分治法为 Θ(nlgn)，动态规划为 Θ(n)。

This file includes a set of maximum_subarray algorithms.
Running times: brute-force method is Θ(n^2), divide-and-conquer 
method is Θ(nlgn) and dynamic method is Θ(n).
*/

#ifndef _MAXIMUM_SUBARRAY_H
#define	_MAXIMUM_SUBARRAY_H

#include <limits>

#define STANDARD_HOLD 2

class max_subarray
{
public:
	max_subarray() 
		: low(-1), high(-1), sum(std::numeric_limits<int>::min()) {}
	max_subarray(int b, int l, int s) 
		: low(b), high(l), sum(s) {}
	max_subarray(const max_subarray &right) 
		: low(right.low), high(right.high), sum(right.sum) {}

	int low, high, sum;
};

// 暴力求解法。
// Brute-force method.
max_subarray bf_find_maximum_subarray(const int arr[], int low, int high) 
{
	max_subarray result(low, low, arr[low]);
	for (int i = low; i <= high; ++i) {
		int sum = arr[i];
		for (int j = i + 1; j <= high; ++j) {
			sum += arr[j];
			if (sum > result.sum) {
				result.low = i;
				result.high = j;
				result.sum = sum;
			}
		}
	}
	return result;
}

// 分治法求解。
// divide-and-conquer method.
max_subarray find_max_crossing_subarray(int arr[], int low, int mid, int high)
{
	int left_sum = std::numeric_limits<int>::min();
	int max_left = mid;
	int sum = 0;
	for (int i = mid; i >= low; --i) {
		sum += arr[i];
		if (sum > left_sum) {
			left_sum = sum;
			max_left = i;
		}
	}

	int right_sum = std::numeric_limits<int>::min();
	int max_right = mid + 1;
	sum = 0;
	for (int i = mid + 1; i <= high; ++i) {
		sum += arr[i];
		if (sum > right_sum) {
			right_sum = sum;
			max_right = i;
		}
	}

	return max_subarray(max_left, max_right, left_sum + right_sum);
}

max_subarray recursive_find_maximum_subarray(int arr[], int low, int high)
{
	if (low == high)
		return max_subarray(low, high, arr[low]);
	else {
		int mid = (low + high) / 2;
		max_subarray left = recursive_find_maximum_subarray(arr, low, mid);
		max_subarray right = recursive_find_maximum_subarray(arr, mid + 1, high);
		max_subarray cross = find_max_crossing_subarray(arr, low, mid, high);
		if (left.sum >= right.sum && left.sum >= cross.sum)
			return left;
		else if (right.sum >= left.sum && right.sum >= cross.sum)
			return right;
		else
			return cross;
	}
}

max_subarray mixed_find_maximum_subarray(int arr[], int low, int high)
{
	if (high - low <= STANDARD_HOLD)
		return bf_find_maximum_subarray(arr, low, high);
	else {
		int mid = (low + high) / 2;
		max_subarray left = recursive_find_maximum_subarray(arr, low, mid);
		max_subarray right = recursive_find_maximum_subarray(arr, mid + 1, high);
		max_subarray cross = find_max_crossing_subarray(arr, low, mid, high);
		if (left.sum >= right.sum && left.sum >= cross.sum)
			return left;
		else if (right.sum >= left.sum && right.sum >= cross.sum)
			return right;
		else
			return cross;
	}
}

// 动态规划
// Dynamic method.
max_subarray dynamic_find_maximum_subarray(int arr[], int low, int high)
{
	max_subarray result(low, low, arr[low]);
	max_subarray temp(result);
	for (int i = low + 1; i <= high; ++i) {
		if (temp.sum < 0) {
			temp.low = i;
			temp.high = i;
			temp.sum = arr[i];
		}
		else {
			temp.high = i;
			temp.sum += arr[i];
		}

		if (temp.sum > result.sum) {
			result.low = temp.low;
			result.high = temp.high;
			result.sum = temp.sum;
		}
	}
	return result;
}

#endif // !_MAXIMUM_SUBARRAY_H