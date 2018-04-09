/*
此文件包含的算法为通过修改归并排序查找序列中的逆序对，具体描述见思考题2-4。
时间复杂度：Θ(nlgn)。

This file includes inversions count algorithm by modifing merge sort. Details
can refer to problem 2-4.
Running time is Θ(nlgn).
*/ 

#ifndef _INVERSION_H
#define _INVERSION_H

#include <limits>

int modified_merge(int arr[], int low, int mid, int high)
{
	int count = 0;

	if (arr[mid] < arr[mid + 1]) return count;

	int left_len = mid - low + 1, right_len = high - mid;
	int *left = new int[left_len + 1];
	int *right = new int[right_len + 1];
	for (int i = 0; i != left_len; ++i)
		left[i] = arr[low + i];
	for (int i = 0; i != right_len; ++i)
		right[i] = arr[mid + i + 1];
	// 加入哨兵。We place on the bottom of each pile a sentinel card.
	left[left_len] = std::numeric_limits<int>::max();
	right[right_len] = std::numeric_limits<int>::max();

	int i = 0, j = 0;
	for (int k = low; k <= high; ++k) {
		if (left[i] <= right[j])
			arr[k] = left[i++];
		else {
			arr[k] = right[j++];
			count += mid - low + 1 - i;
		}
	}

	delete[] left;
	delete[] right;
	return count;
}

int modified_merge_sort(int arr[], int low, int high)
{
	int count = 0;

	if (low < high) {
		int mid = (low + high) / 2;
		count += modified_merge_sort(arr, low, mid);
		count += modified_merge_sort(arr, mid + 1, high);
		count += modified_merge(arr, low, mid, high);
	}

	return count;
}

int inversions_count(const int arr[], int low, int high)
{
	int len = high - low + 1;
	int *temp_arr = new int[len];
	for (int i = low; i <= high; ++i)
		temp_arr[i] = arr[i];
	int count = modified_merge_sort(temp_arr, low, high);

	delete[] temp_arr;
	return count;
}

#endif // !_INVERSION_H