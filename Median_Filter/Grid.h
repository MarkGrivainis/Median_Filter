#pragma once
class Grid
{
private:
	
public:
	unsigned int *grid;
	int cols, rows;
	Grid(int rows, int cols);
	~Grid();
	void Print();
	int Count();
	void Clear();
	void Add(float x, float y);
};

