// 此头文件包含各种排序算法。
// This head file includes various sort algorithms.


#ifndef _SORT_H
#define _SORT_H

#include <limits>

#define THRESHOLD 16

// 插入排序。时间复杂度，在最好情况下：Θ(n)，最坏情况下：Θ(n^2)。
// Insertion sort. Running time: Best case: Θ(n), worst case: Θ(n^2).
void insertion_sort(int arr[], int low, int high)
{
	for (int i = low + 1; i <= high; ++i) {
		int key = arr[i];

		int j = i - 1;
		for (; j >= low && arr[j] > key; --j)
			arr[j + 1] = arr[j];
		arr[j + 1] = key;
	}
}


// 归并排序。时间复杂度为 Θ(nlgn)。
// Merge sort. Running time: Θ(nlgn), where lgn stands for log2^n.
void merge(int arr[], int low, int mid, int high)
{
	if (arr[mid] < arr[mid + 1]) return;  // 已经有序。The array is sorted.

	int left_len = mid - low + 1, right_len = high - mid;
	int *left = new int[left_len + 1];
	int *right = new int[right_len + 1];
	for (int i = 0; i != left_len; ++i)
		left[i] = arr[low + i];
	for (int i = 0; i != right_len; ++i)
		right[i] = arr[mid + i + 1];
	// 加入哨兵。We place on the bottom of each pile a sentinel card.
	left[left_len] = std::numeric_limits<int >::max();
	right[right_len] = std::numeric_limits<int>::max();

	int i = 0, j = 0;
	for (int k = low; k <= high; ++k) {
		if (left[i] <= right[j])
			arr[k] = left[i++];
		else
			arr[k] = right[j++];
	}

	delete[] left;
	delete[] right;
}

void merge_sort(int arr[], int low, int high)
{
	if (low < high) {
		int mid = (low + high) / 2;
		merge_sort(arr, low, mid);
		merge_sort(arr, mid + 1, high);
		merge(arr, low, mid, high);
	}
}

// 当待排序子序列足够小时，使用插入排序进行优化。
// Using insertion sort within merge sort when subproblems become sufficiently small.
void mixed_merge_sort(int arr[], int low, int high)
{
	if ((high - low) < THRESHOLD)
		insertion_sort(arr, low, high);
	else {
		int mid = (low + high) / 2;
		mixed_merge_sort(arr, low, mid);
		mixed_merge_sort(arr, mid + 1, high);
		merge(arr, low, mid, high);
	}
}


// 冒泡排序。时间复杂度在最好情况下为 Θ(n)，最坏情况下为 Θ(n^2)。
// Bubble sort algorithm. Running time: Best case: Θ(n), worst case: Θ(n^2).
void bubble_sort(int arr[], int len)
{
	bool isSorted = false;
	for (int i = 0; !isSorted && i != len - 1; ++i) {
		isSorted = true;
		for (int j = len - 1; j != i; --j) {
			if (arr[j] < arr[j - 1]) {
				isSorted = false;
				int temp = arr[j];
				arr[j] = arr[j - 1];
				arr[j - 1] = temp;
			}
		}
	}
}

#endif // !_SORT_H