#include "Grid.h"
#include <algorithm>
#include <vector>

#pragma once
class Filter
{
public:
	Filter();
	~Filter();
	static int median(std::vector<int> values)
	{
		if (values.empty()) return 0;
		else {
			std::sort(values.begin(), values.end());
			if (values.size() % 2 == 0)
			{
				printf("even");
				return (values[values.size() / 2 - 1] + values[values.size() / 2]) / 2;
			}
			else
				return values[values.size() / 2];
		}
	}
	void m_Filter_zero(Grid &grid, int radius);
	//bubble sort the entire window
	static void m_Filter_fullsort(Grid &padded, Grid &image, int radius)
	{
			for (int y = radius; y < image.rows + radius; ++y)
			{
				int top = std::max(y - radius, 0);
				int bottom = std::min(y + radius, 4);
		
				for (int x = radius; x < image.cols + radius; ++x)
				{
					int left = std::max(x - radius, 0);
					int right = std::min(x + radius, 4);
					std::vector<int> values;
					for (int v = y - radius; v <= y + radius; ++v)
					{
						for (int u = x - radius; u <= x + radius; ++u)
						{
							values.push_back(padded.grid[v * padded.cols + u]);
						}
					}
					int med = median(values);
					image.grid[(y - radius) * image.cols + (x - radius)] = med;
				}
			}
	}
	//bubble sort for half the window
	static void m_Filter_half(Grid &padded, Grid &image, const int size)
	{
		int radius = (int)(size - 1) / 2;
		int windowSize = size*size;
		int halfWindow = ((windowSize + 1) / 2);
		for (int y = radius; y < image.rows + radius; ++y)
		{
			int top = std::max(y - radius, 0);
			int bottom = std::min(y + radius, 4);

			for (int x = radius; x < image.cols + radius; ++x)
			{
				int left = std::max(x - radius, 0);
				int right = std::min(x + radius, 4);
				int k = 0;
				int *window = new int[windowSize];
				for (int v = y - radius; v <= y + radius; ++v)
				{
					for (int u = x - radius; u <= x + radius; ++u)
					{
						window[k++] = padded.grid[v * padded.cols + u];
					}
				}
				for (int j = 0; j < halfWindow; ++j)
				{
					//   Find position of minimum element
					int min = j;
					for (int l = j + 1; l < windowSize; ++l)
						if (window[l] < window[min])
							min = l;
					//   Put found minimum element in its place
					const int temp = window[j];
					window[j] = window[min];
					window[min] = temp;
				}
				image.grid[(y - radius) * image.cols + (x - radius)] = window[halfWindow-1];
				delete[] window;
			}
		}
	}
	//quick select found on google 
	//http://blog.teamleadnet.com/2012/07/quick-select-algorithm-find-kth-element.html
	static void m_Filter_quickselect(Grid &padded, Grid &image, const int size)
	{
		int radius = (int)(size - 1) / 2;
		int windowSize = size*size;
		int halfWindow = ((windowSize + 1) / 2);
		for (int y = radius; y < image.rows + radius; ++y)
		{
			int top = std::max(y - radius, 0);
			int bottom = std::min(y + radius, 4);

			for (int x = radius; x < image.cols + radius; ++x)
			{
				int left = std::max(x - radius, 0);
				int right = std::min(x + radius, 4);
				int k = 0;
				int *window = new int[windowSize];
				for (int v = y - radius; v <= y + radius; ++v)
				{
					for (int u = x - radius; u <= x + radius; ++u)
					{
						window[k++] = padded.grid[v * padded.cols + u];
					}
				}
				int median = halfWindow - 1;
				int from = 0, to = windowSize - 1;
				while (from < to)
				{
					int r = from, w = to;
					int mid = window[(r + w) / 2];

					while (r < w)
					{
						if (window[r] >= mid)
						{
							int temp = window[w];
							window[w] = window[r];
							window[r] = temp;
							w--;
						}
						else
						{
							r++;
						}
					}
					if (window[r] > mid)
					{
						r--;
					}
					if (median <= r)
					{
						to = r;
					}
					else
					{
						from = r + 1;
					}
				}
				image.grid[(y - radius) * image.cols + (x - radius)] = window[median];
				delete[] window;
			}
		}
	}
	/*
	* This Quickselect routine is based on the algorithm described in
	* "Numerical recipes in C", Second Edition,
	* Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
	* This code by Nicolas Devillard - 1998. Public domain.
	*/
	static void m_Filter_quickselect2(Grid &padded, Grid &image, const int size)
	{
		int radius = (int)(size - 1) / 2;
		int windowSize = size*size;
		for (int y = radius; y < image.rows + radius; ++y)
		{
			int top = std::max(y - radius, 0);
			int bottom = std::min(y + radius, 4);

			for (int x = radius; x < image.cols + radius; ++x)
			{
				int left = std::max(x - radius, 0);
				int right = std::min(x + radius, 4);
				int k = 0;
				int *window = new int[windowSize];
				for (int v = y - radius; v <= y + radius; ++v)
				{
					for (int u = x - radius; u <= x + radius; ++u)
					{
						window[k++] = padded.grid[v * padded.cols + u];
					}
				}
				int low, high;
				int median;
				int middle, ll, hh;

				low = 0; high = windowSize - 1; median = (low + high) / 2;
				for (;;)
				{
					if (high <= low)
						break;
					if (high == low + 1)
					{
						if (window[low] > window[high])
						{
							int temp = window[low];
							window[low] = window[high];
							window[high] = temp;
						}
						break;
					}

					middle = (low + high) / 2;
					if (window[middle] > window[high])
					{
						int temp = window[middle];
						window[middle] = window[high];
						window[high] = temp;
					}
					if (window[low] > window[high])
					{
						int temp = window[low];
						window[low] = window[high];
						window[high] = temp;
					}
					if (window[middle] > window[low])
					{
						int temp = window[low];
						window[low] = window[middle];
						window[middle] = temp;
					}
					int temp = window[low + 1];
					window[low + 1] = window[middle];
					window[middle] = temp;

					ll = low + 1;
					hh = high;
					for (;;)
					{
						do{
							ll++;
						} while (window[low] > window[ll]);
						do
						{
							hh--;
						} while (window[hh] > window[low]);
						if (hh < ll)
							break;
						int temp = window[ll];
						window[ll] = window[hh];
						window[hh] = temp;
					}

					temp = window[low];
					window[low] = window[hh];
					window[hh] = temp;

					if (hh <= median)
						low = ll;
					if (hh >= median)
						high = hh - 1;
				}
				image.grid[(y - radius) * image.cols + (x - radius)] = window[median];
				delete[] window;
			}
		}
	}
	
};