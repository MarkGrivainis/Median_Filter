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

};