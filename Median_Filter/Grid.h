#pragma once
#include <iostream>
class Grid
{
private:
	
public:
	int *grid;
	int _size;
	int cols, rows;
	Grid(int rows, int cols);
	Grid(const Grid &other) :rows(other.rows), cols(other.cols)
	{
		grid = new int[other._size];
		_size = other._size;
		memcpy(grid, other.grid, sizeof(int) * _size);
	}
	~Grid();
	void Print();
	void PrintRange(int x1, int x2, int y1, int y2);
	int Count();
	void Clear();
	void Add(float x, float y);
	Grid Pad(int radius);
	Grid ZeroPad(int radius);
	void set(int);
};

