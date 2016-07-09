/*
* Defines configuration:
* _DEBUG: The application will display debug information for every iteration, success and failure messages
* _ITERATIONS=: How many iterations of the simulation shall be computed
* _DIM=: The amount of cells for each axis in the space
* _1D: The application will execute the simulation in one dimension
* _2D: The application will execute the simulation in two dimension
* _3D: The application will execute the simulation in three dimensions
* _STEP_ALTERNATIVE: Uses an alternative method for computing a step
* _INLINE: Forces the inline of the step function
* _OFFSET_BUFFER: One of the buffers will have an offset, with the objective of avoiding cache block sharing when cache associativity is too low
* _OFFSET_BUFFER_SIZE=: If _OFFSET_BUFFER is defined, this will be the offset size, in bytes
* _RESTRICT: Adds the restrict keyword, to warn the compiler about no memory aliasing
*/

typedef float cell;
#define PRINT_TYPE "%.1f"
#define DIFFUSE_FACTOR 0.1f

#ifdef _RESTRICT
#define _RESTRICT __restrict__
#else
#define _RESTRICT
#endif

#ifdef _INLINE
#define _INLINE __attribute__((always_inline))
#else
#define _INLINE __attribute__((noinline))
#endif

#if !(defined(_1D) ^ defined(_2D) ^ defined(_3D))
#error "Must define _1D xor _2D xor _3D"
#endif

#ifndef _DIM
#error "Must define _DIM"
#endif

#ifndef _ITERATIONS
#error "Must define _ITERATIONS"
#endif

#ifdef _OFFSET_BUFFER
#ifndef _OFFSET_BUFFER_SIZE
#error "Must define _OFFSET_BUFFER_SIZE"
#endif
#endif

#ifdef _1D
#define SIZE _DIM
int getOffset(int x);
#endif

#ifdef _2D
#define SIZE _DIM*_DIM
int getOffset(int x, int y);
#endif

#ifdef _3D
#define SIZE _DIM*_DIM*_DIM
int getOffset(int x, int y, int z);
#endif

#ifdef _DEBUG
void debugPrint(cell *data);
#include <stdio.h>
#endif

#include <stdlib.h>

#ifdef _1D
inline int getOffset(int x)
{
	return x;
}

#ifndef _STEP_ALTERNATIVE
_INLINE void step(cell *_RESTRICT read_buffer, cell *_RESTRICT write_buffer)
{
	int i;
	for (i = 1; i < _DIM - 1; i++)
	{
		write_buffer[getOffset(i)] =
			read_buffer[getOffset(i)] +
			(DIFFUSE_FACTOR *
			(read_buffer[getOffset(i + 1)] +
				read_buffer[getOffset(i - 1)] -
				2 * read_buffer[getOffset(i)]));
	}
}

#else
_INLINE void step(cell *_RESTRICT read_buffer, cell *_RESTRICT write_buffer)
{
	int i;
	for (i = 1; i < _DIM - 1; i++)
	{
		write_buffer[getOffset(i)] += read_buffer[getOffset(i)] - read_buffer[getOffset(i)] * 2 * DIFFUSE_FACTOR;
		cell temp = read_buffer[getOffset(i)] * DIFFUSE_FACTOR;
		write_buffer[getOffset(i + 1)] += temp;
		write_buffer[getOffset(i - 1)] += temp;
	}
}

#endif
#ifdef _DEBUG
void debugPrint(cell *data)
{
	int i;
	for (i = 0; i < _DIM; i++)
	{
		printf(PRINT_TYPE"\t", data[getOffset(i)]);
	}
	printf("\n");
}
#endif
#endif

#ifdef _2D
inline int getOffset(int x, int y)
{
	return (x*_DIM) + y;
}

#ifndef _STEP_ALTERNATIVE
_INLINE void step(cell *_RESTRICT read_buffer, cell *_RESTRICT write_buffer)
{
	int i, j;
	for (i = 1; i < _DIM - 1; i++)
	{
		for (j = 1; j < _DIM - 1; j++)
		{
			write_buffer[getOffset(i, j)] =
				read_buffer[getOffset(i, j)] +
				(DIFFUSE_FACTOR *
				(read_buffer[getOffset(i + 1, j)] +
					read_buffer[getOffset(i - 1, j)] +
					read_buffer[getOffset(i, j + 1)] +
					read_buffer[getOffset(i, j - 1)] -
					4 * read_buffer[getOffset(i, j)]));
		}
	}
}

#else
_INLINE void step(cell *_RESTRICT read_buffer, cell *_RESTRICT write_buffer)
{
	int i, j;
	for (i = 1; i < _DIM - 1; i++)
		for (j = 1; j < _DIM - 1; j++)
		{
			cell celula = read_buffer[getOffset(i, j)];
			write_buffer[getOffset(i, j)] += celula - (celula * 4 * DIFFUSE_FACTOR);
			cell temp = celula * DIFFUSE_FACTOR;
			write_buffer[getOffset(i + 1, j)] += temp;
			write_buffer[getOffset(i - 1, j)] += temp;
			write_buffer[getOffset(i, j + 1)] += temp;
			write_buffer[getOffset(i, j - 1)] += temp;
		}
}
#endif
#ifdef _DEBUG
void debugPrint(cell *data)
{
	int i, j;
	for (i = 0; i < _DIM; i++)
	{
		for (j = 0; j < _DIM; j++)
		{
			printf(PRINT_TYPE"\t", data[getOffset(i, j)]);
		}
		printf("\n");
	}
	printf("\n");
}
#endif
#endif

#ifdef _3D
inline int getOffset(int x, int y, int z)
{
	return (x*_DIM*_DIM) + y*_DIM + z;
}

#ifndef _STEP_ALTERNATIVE
_INLINE void step(cell *_RESTRICT read_buffer, cell *_RESTRICT write_buffer)
{
	int i, j, k;
	for (i = 1; i < _DIM - 1; i++)
		for (j = 1; j < _DIM - 1; j++)
			for (k = 1; k < _DIM - 1; k++)
			{
				write_buffer[getOffset(i, j, k)] =
					read_buffer[getOffset(i, j, k)] +
					(DIFFUSE_FACTOR *
					(read_buffer[getOffset(i + 1, j, k)] +
						read_buffer[getOffset(i - 1, j, k)] +
						read_buffer[getOffset(i, j + 1, k)] +
						read_buffer[getOffset(i, j - 1, k)] +
						read_buffer[getOffset(i, j, k - 1)] +
						read_buffer[getOffset(i, j, k + 1)] -
						6 * read_buffer[getOffset(i, j, k)]));
			}
}
#else
_INLINE void step(cell *_RESTRICT read_buffer, cell *_RESTRICT write_buffer)
{
	int i, j, k;
	for (i = 1; i < _DIM - 1; i++)
		for (j = 1; j < _DIM - 1; j++)
			for (k = 1; k < _DIM - 1; k++)
			{
				cell celula = read_buffer[getOffset(i, j, k)];
				write_buffer[getOffset(i, j, k)] += celula - (celula * 6 * DIFFUSE_FACTOR);
				cell temp = celula * DIFFUSE_FACTOR;
				write_buffer[getOffset(i + 1, j, k)] += temp;
				write_buffer[getOffset(i - 1, j, k)] += temp;
				write_buffer[getOffset(i, j + 1, k)] += temp;
				write_buffer[getOffset(i, j - 1, k)] += temp;
				write_buffer[getOffset(i, j, k + 1)] += temp;
				write_buffer[getOffset(i, j, k - 1)] += temp;
			}
}
#endif
#ifdef _DEBUG
void debugPrint(cell *data)
{
	int i, j;
	for (i = 0; i < _DIM; i++)
	{
		for (j = 0; j < _DIM; j++)
		{
			printf(PRINT_TYPE"\t", data[getOffset(i, j, _DIM / 2)]);
		}
		printf("\n");
	}
	printf("\n");
}
#endif
#endif

int main()
{
	cell *buffer[2];
	cell *aux;
	// We need a single contiguous memory block
#ifndef _OFFSET_BUFFER
	buffer[0] = (cell*)malloc(sizeof(cell) * SIZE);
#else // _OFFSET_BUFFER
	buffer[0] = (cell*)malloc(sizeof(cell) * SIZE + _OFFSET_BUFFER_SIZE);
	buffer[0] += _OFFSET_BUFFER_SIZE;
#endif
	buffer[1] = (cell*)malloc(sizeof(cell) * SIZE);

	if ((buffer[0] == NULL) || (buffer[1] == NULL))
	{
#ifdef _DEBUG
		printf("Error: Failed to allocate memory\n");
#endif
#ifdef _OFFSET_BUFFER
		if (buffer[0] != NULL)
			buffer[0] -= _OFFSET_BUFFER_SIZE;
#endif
		free(buffer[0]);
		free(buffer[1]);
		return -1;
	}

	// Clean the memory
	memset(buffer[0], 0, sizeof(cell)*SIZE);
	memset(buffer[1], 0, sizeof(cell)*SIZE);

	// Add a initial value to the simulation
#ifdef _1D
	buffer[0][getOffset(_DIM / 2)] = 1000;
#endif
#ifdef _2D
	buffer[0][getOffset(_DIM / 2, _DIM / 2)] = 10000;
#endif
#ifdef _3D
	buffer[0][getOffset(_DIM / 2, _DIM / 2, _DIM / 2)] = 10000;
#endif

	int iteration;
	for (iteration = 0; iteration < _ITERATIONS; iteration++)
	{
#ifdef _DEBUG
		debugPrint(buffer[0]);
#endif
		step(buffer[0], buffer[1]);
		aux = buffer[0];
		buffer[0] = buffer[1];
		buffer[1] = aux;
	}
#ifdef _OFFSET_BUFFER
	buffer[0] -= _OFFSET_BUFFER_SIZE;
#endif

	free(buffer[0]);
	free(buffer[1]);
#ifdef _DEBUG
	printf("Simulation completed successfully\n");
#endif
	return 0;
}
