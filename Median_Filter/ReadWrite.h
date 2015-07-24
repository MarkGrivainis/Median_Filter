#include "Grid.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <omp.h>

#pragma once
class ReadWrite
{
public:
	ReadWrite();
	~ReadWrite();
	static int LoadData_s(Grid &grid, char* name)
	{
		double r_t1, r_t2, b_t1, b_t2;
		double read = 0, bin = 0;
		grid.Clear();
		FILE *stream;
		const long num = 100000000;
		float *list = new float[num];
		__int64 filesize;
		unsigned int  numread = 1;
		int readcount = 0;
		if (fopen_s(&stream, name, "rb") == 0)
		{
			if (_fseeki64(stream, 0l, SEEK_END) == 0)
			{
				filesize = _ftelli64(stream);
				if (filesize == -1)
					std::cout << "ftell failed" << std::endl;
				else
					std::cout << "--- # bytes in file \t\t: " << filesize << std::endl;
				_fseeki64(stream, 0l, SEEK_SET); // to be tidy
			}
			else
			{
				std::cout << "fseek failed" << std::endl;
			}
			
			while (true)
			{
				r_t1 = omp_get_wtime();
				numread = fread(list, sizeof(float), num, stream);
				r_t2 = omp_get_wtime();
				read += r_t2 - r_t1;
				if (numread == 0)
					break;
				b_t1 = omp_get_wtime();
				for (int i = 0; i < numread; i += 2)
				{
					int x_i = (int)floor(list[i] * grid.cols);
					int y_i = (int)floor(list[i + 1] * grid.rows);

					if (x_i == grid.cols)
						x_i--;
					if (y_i == grid.rows)
						y_i--;
					++grid.grid[y_i * grid.cols + x_i];
				}
				b_t2 = omp_get_wtime();
				bin += b_t2 - b_t1;

				readcount++;

			}
			
			printf("--- # floats \t\t\t: %d\n", filesize / 8);
			printf("--- Time reading \t\t: %-5.5f sec.\n", read);
			printf("--- Time binning \t\t: %-5.5f sec.\n", bin);
			printf("--- # of floats binned \t\t: %d\n", grid.Count());

			fclose(stream);
		}
		else
		{
			printf("File could not be opened\n");
			return -1;
		}

		std::cout << "--- # of reads required \t: " << readcount << std::endl;
		return 0;
	}
	static int LoadData_omp(Grid &grid, char* name)
	{
		double r_t1, r_t2, b_t1, b_t2;
		double read = 0, bin = 0;
		grid.Clear();
		FILE *stream;
		const long num = 100000000;
		float *list = new float[num];
		__int64 filesize;
		unsigned int  numread = 1;
		int readcount = 0;
		if (fopen_s(&stream, name, "rb") == 0)
		{
			if (_fseeki64(stream, 0l, SEEK_END) == 0)
			{
				filesize = _ftelli64(stream);
				if (filesize == -1)
					std::cout << "ftell failed" << std::endl;
				else
					std::cout << "--- # bytes in file \t\t: " << filesize << std::endl;
				_fseeki64(stream, 0l, SEEK_SET); // to be tidy
			}
			else
			{
				std::cout << "fseek failed" << std::endl;
			}
			while (true)
			{
				r_t1 = omp_get_wtime();
				numread = fread(list, sizeof(float), num, stream);
				r_t2 = omp_get_wtime();
				read += r_t2 - r_t1;
				if (numread == 0)
					break;
				b_t1 = omp_get_wtime();
				#pragma omp parallel for schedule(static)
				for (int i = 0; i < numread; i += 2)
				{
					int x_i = (int)(list[i] * grid.cols);
					int y_i = (int)(list[i + 1] * grid.rows);
					
					if (x_i == grid.cols)
						x_i--;
					if (y_i == grid.rows)
						y_i--;
					int index = y_i * grid.cols + x_i;
					#pragma omp atomic
					++grid.grid[index];
				}
				b_t2 = omp_get_wtime();
				bin += b_t2 - b_t1;
				readcount++;
			}
			printf("--- # floats \t\t\t: %d\n", filesize / 8);
			printf("--- Time reading \t\t: %-5.5f sec.\n", read);
			printf("--- Time binning \t\t: %-5.5f sec.\n", bin);
			printf("--- # of floats binned \t\t: %d\n", grid.Count());

			fclose(stream);
		}
		else
		{
			printf("File could not be opened\n");
			return -1;
		}

		std::cout << "--- # of reads required \t: " << readcount << std::endl;
		return 0;
	}

	static void WriteData(Grid &grid, char* name)
	{
		float addition_x = (1 / (float)grid.cols) / 2;
		float addition_y = (1 / (float)grid.rows) / 2;
		std::ofstream myfile;
		myfile.open(name);
		std::string line;

		for (int y = 0; y < grid.rows + 1; ++y)
		{
			line = "";
			for (int x = 0; x < grid.cols + 1; ++x)
			{
				if (x == 0 && y == 0)
					line += "";
				else if (y == 0 && x > 0)
					line += ", " + std::to_string(((x - 1) / (float)grid.cols) + addition_x) + "";
				else if (x == 0 && y >0)
					line += std::to_string(((y - 1) / (float)grid.rows) + addition_y) + "";
				else
					line += ", " + std::to_string(grid.grid[(y - 1) * grid.cols + (x - 1)]) + "";

			}
			myfile << line + "\n";
		}
		myfile.close();
	}
};

