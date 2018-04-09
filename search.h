// 此头文件包含各种排序算法。
// This head file includes various search algorithms.


#ifndef _SEARCH_H
#define _SEARCH_H

// 线性查找。时间复杂度：Θ(n)。
// Linear search algorithm. Running time is Θ(n).
int linear_search(const int arr[], int len, int key)
{
	for (int i = 0; i != len; ++i)
		if (arr[i] == key) return i;
	return -1;
}

// 二分查找。时间复杂度：Θ(lgn)。
// Binary search algorithm. Running time is Θ(lgn).
int binary_search(const int arr[], int len, int key)
{
	int begin = 0, mid = len / 2;

	while (begin != len) {
		mid = (begin + len) / 2;
		if (arr[mid] >= key)
			len = mid;
		else if (arr[mid] < key)
			begin = mid + 1;
	}

	if (arr[begin] == key)
		return begin;
	else
		return -1;
}

#endif // !_SEARCH_H