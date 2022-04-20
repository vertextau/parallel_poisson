#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "params.h"
#include "function.h"

void set_boundaries(double *x, int x_size, int rank,  int rank_size, int *coord, double h)
{
    /*
     * map local coordinates to global
     */
    int global_j = coord[0];
    int global_i = abs(coord[1] - (rank_size - 1));

    /*
     * last row
     */
    if (coord[1] == 0) {
        int last_row = (x_size-1) * x_size;
        int y_i = (global_i*(x_size-2) + x_size-1);

        for (int j = 0; j < x_size; ++j) {
            int x_j = (global_j*(x_size-2)+j);
            x[last_row + j] = boundaries(x_j*h, y_i*h);
        }
    }

    /*
     * left column
     */
    if (coord[0] == 0) {
        int x_j = global_j*(x_size-2);

        for (int i = 0; i < x_size; ++i) {
            int y_i = global_i*(x_size-2)+i;
            x[i*x_size] = boundaries(x_j*h, y_i*h);
        }
    }

    /*
     * first row
     */
    if (coord[1] == rank_size-1) {
        int y_i = (global_i*(x_size-2));

        for (int j = 0; j < x_size; ++j) {
            int x_j = (global_j*(x_size-2)+j);
            x[j] = boundaries(x_j*h, y_i*h);
        }
    }

    /*
     * right column
     */
    if (coord[0] == rank_size-1) {
        int x_j = (global_j*(x_size-2)+x_size-1);

        for (int i = 0; i < x_size; ++i) {
            int y_i = (global_i*(x_size-2)+i);
            x[i*x_size + x_size-1] = boundaries(x_j*h, y_i*h);
        }
    }
}

void update_neighbors(double *x, int x_size, MPI_Comm grid_comm, int *neighbours, MPI_Datatype col_type)
{
    int tag = 0;

    MPI_Sendrecv(&x[(x_size-2)*x_size], x_size, MPI_DOUBLE, neighbours[SOUTH], tag,
            &x[0], x_size, MPI_DOUBLE, neighbours[NORTH], tag, grid_comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&x[x_size], x_size, MPI_DOUBLE, neighbours[NORTH], tag,
            &x[(x_size-1)*x_size], x_size, MPI_DOUBLE, neighbours[SOUTH], tag, grid_comm, MPI_STATUS_IGNORE);

    tag = 1;
    MPI_Sendrecv(&x[x_size-2], 1, col_type, neighbours[EAST], tag,
            &x[0], 1, col_type, neighbours[WEST], tag, grid_comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&x[1], 1, col_type, neighbours[WEST], tag,
            &x[x_size-1], 1, col_type, neighbours[EAST], tag, grid_comm, MPI_STATUS_IGNORE);

}

void red_nodes(double *x, int x_size, int *grid_coord, int q, double h, double omega)
{
    /*
     * map local coordinates to global
     */
    int global_j = grid_coord[0];
    int global_i = abs(grid_coord[1] - (q - 1));

    for (int i = 1; i < x_size-1; ++i) {
        for (int j = 1; j < x_size-1; ++j) {
            int x_i = (global_i*(x_size-2)+i);
            int y_j = (global_j*(x_size-2)+j);

            if ((x_i+y_j)%2 == 0) {
                double r = (1/4.0) * (x[(i-1)*x_size+j]
                        + x[(i+1)*x_size+j] + x[i*x_size+j+1]
                        + x[i*x_size+j-1] + h*h*f(x_i*h, y_j*h)) - x[i*x_size+j];
                x[i*x_size+j] += omega*r;
            }
        }
    }    
}

void black_nodes(double *x, int x_size, int *grid_coord, int q, double h, double omega)
{
    /*
     * map local coordinates to global
     */
    int global_j = grid_coord[0];
    int global_i = abs(grid_coord[1] - (q - 1));

    for (int i = 1; i < x_size-1; ++i) {
        for (int j = 1; j < x_size-1; ++j) {
            int x_i = (global_i*(x_size-2)+i);
            int y_j = (global_j*(x_size-2)+j);

            if ((x_i+y_j)%2 != 0) {
                double r = (1/4.0) * (x[(i-1)*x_size+j]
                        + x[(i+1)*x_size+j] + x[i*x_size+j+1]
                        + x[i*x_size+j-1] + h*h*f(x_i*h, y_j*h)) - x[i*x_size+j];
                x[i*x_size+j] += omega*r;
            }
        }
    }
}

void save_matrix(FILE *f, double *mat, int x_size, int rank, int q, MPI_Comm grid_comm)
{
    int grid_row, grid_col;
    int source;
    int coords[2];

    if (rank == 0) {
        double *temp = malloc(sizeof(double)*x_size*x_size);

        for (int mat_row = 0, zrow = 0; mat_row < x_size*q; ++mat_row) {

            grid_row = mat_row/x_size;
            coords[1] = abs(grid_row-(q-1));
            for (grid_col = 0; grid_col < q; grid_col++) {
                coords[0] = grid_col;
                MPI_Cart_rank(grid_comm, coords, &source);
                if (source == 0) {
                    for(int mat_col = 0; mat_col < x_size; ++mat_col) {
                        fprintf(f, "%f ", mat[zrow*x_size+mat_col]);
                    }
                    ++zrow;
                } else {
                    MPI_Recv(temp, x_size, MPI_DOUBLE, source, 0, grid_comm, MPI_STATUS_IGNORE);
                    for(int mat_col = 0; mat_col < x_size; ++mat_col)
                        fprintf(f, "%f ", temp[mat_col]);
                }
            }
            fprintf(f, "\n");
        }
        free(temp);
    } else {
        for (int mat_row = 0; mat_row < x_size; ++mat_row) {
            MPI_Send(&mat[mat_row*x_size], x_size, MPI_DOUBLE, 0, 0, grid_comm);
        }
    }
}

double calc_difference(double *x, double *temp, int x_size)
{
    double max = 0;

    for (int i = 1; i < x_size-1; ++i) {
        for (int j = 1; j < x_size-1; ++j) {
            double m = fabs(x[i*x_size+j] - temp[i*x_size+j]);
            if (m > max) {
                max = m;
            }
        }
    }

    return max;
}
