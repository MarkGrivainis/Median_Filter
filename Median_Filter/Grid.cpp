#include "Grid.h"


Grid::Grid(int rows, int cols) :rows(rows), cols(cols)
{
	grid = new int[rows * cols];
	_size = rows * cols;
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
			printf("%-10d", grid[i * cols + j]);
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
Grid Grid::Pad(int radius)
{
	Grid bordered = Grid(rows + (radius * 2), cols + (radius * 2));
	bordered.Clear();
	for (int y = 0; y < rows; ++y)
	{
		memcpy(bordered.grid + (y + radius) * bordered.cols + radius, grid + y * cols, cols*sizeof(size_t));
		for (int x = 0; x < radius; ++x)
		{
			bordered.grid[(y + radius) * bordered.cols + x] = grid[y * cols];
			bordered.grid[(y + radius) * bordered.cols + cols + radius + x] = grid[y * cols + cols - 1];
		}
	}
	for (int y = 0; y < radius; ++y)
	{
		memcpy(bordered.grid + y * bordered.cols, bordered.grid + radius * bordered.cols, bordered.cols*sizeof(size_t));
		memcpy(bordered.grid + (bordered.rows - 1 - y) * bordered.cols, bordered.grid + (bordered.rows - 1 - radius) * bordered.cols, bordered.cols*sizeof(size_t));
	}
	return bordered;
	
}
void Grid::set(int number)
{
	for (int i = 0; i < rows * cols; i++)
		grid[i] = number;
}