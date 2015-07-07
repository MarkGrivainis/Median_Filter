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
	static void m_Filter_extended(Grid &grid, int size)
	{
		int radius = (int)(size - 1) / 2;
			Grid bordered = Grid(grid.rows + (radius * 2), grid.cols + (radius * 2));
			bordered.Clear();
			for (int y = 0; y < grid.rows; ++y)
			{
				memcpy(bordered.grid + (y + radius) * bordered.cols + radius, grid.grid + y * grid.cols, grid.cols*sizeof(size_t));
				for (int x = 0; x < radius; ++x)
				{
					bordered.grid[(y + radius) * bordered.cols + x] = grid.grid[y * grid.cols];
					bordered.grid[(y + radius) * bordered.cols + grid.cols+radius+x] = grid.grid[y * grid.cols + grid.cols-1];
				}
			}
			for (int y = 0; y < radius; ++y)
			{
				memcpy(bordered.grid + y * bordered.cols, bordered.grid + radius * bordered.cols, bordered.cols*sizeof(size_t));
				memcpy(bordered.grid + (bordered.rows - 1 - y) * bordered.cols, bordered.grid + (bordered.rows - 1 - radius ) * bordered.cols, bordered.cols*sizeof(size_t));
			}
			//bordered.Print();
			for (int y = radius; y < grid.rows + radius; ++y)
			{
				int top = std::max(y - radius, 0);
				int bottom = std::min(y + radius, 4);
		
				for (int x = radius; x < grid.cols + radius; ++x)
				{
					int left = std::max(x - radius, 0);
					int right = std::min(x + radius, 4);
					std::vector<int> values;
					for (int v = y - radius; v <= y + radius; ++v)
					{
						for (int u = x - radius; u <= x + radius; ++u)
						{
							values.push_back(bordered.grid[v * bordered.cols + u]);
						}
					}
					int med = median(values);
					grid.grid[(y - radius) * grid.cols + (x - radius)] = med;
				}
			}
	}
};

