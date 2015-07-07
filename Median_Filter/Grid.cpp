#include "Grid.h"
#include <iostream>

Grid::Grid(int rows, int cols) :rows(rows), cols(cols)
{
	grid = new unsigned int[rows * cols];
	for (int i = 0; i < rows * cols; i++)
		grid[i] = 0;
}


Grid::~Grid()
{
	delete[] grid;
}

void Grid::Print()
{
	printf("\n");
	long counter = 0;
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; ++j)
		{
			printf("%d\t", grid[i * cols + j]);
			counter += grid[i * cols + j];
		}
		printf("\n");
	}
	printf("total = %d\n", counter);
}

int Grid::Count()
{
	long counter = 0;
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; ++j)
		{
			//printf("grid [%d][%d] = %d\n", i, j, grid[i][j]);
			counter += grid[i * cols + j];
		}
	return counter;
}

void Grid::Clear()
{
	for (int i = 0; i < rows * cols; i++)
		grid[i] = 0;
}

void Grid::Add(float x, float y )
{
	int x_i = (int)floor(x * cols);
	int y_i = (int)floor(y * rows);

	if (x_i == cols)
		x_i--;
	if (y_i == rows)
		y_i--;
	grid[y_i * cols + x_i]++;
}