#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "function.h"
#include "utils.h"
#include "params.h"

#define MAX_ITER 10000

int main(int argc, char** argv)
{
    int rank_no, proc_num, grid_rank, q, n, N, row_num, col_num;

    MPI_Comm grid_comm, row_comm, column_comm;

    MPI_Datatype column_type;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_no);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_num);

    n = 1024;
    q = (int) sqrt((double) proc_num);
    N = n/q;
    double omega = 1.5;

    /*
     * column type
     */
    MPI_Type_vector(N+2, 1, N+2, MPI_DOUBLE, &column_type);
    MPI_Type_commit(&column_type);


    int dim[2], wrap[2], coord[2], free_coord[2];

    dim[0] = dim[1] = q;
    wrap[0] = wrap[1] = 0;

    MPI_Cart_create(MPI_COMM_WORLD, 2, dim, wrap, 1, &grid_comm);
    MPI_Comm_rank(grid_comm, &grid_rank);
    MPI_Cart_coords(grid_comm, grid_rank, 2, coord);

    /*
     * neighbours coord
     */
    int neighbours_ranks[4];

    double *result_matrix = malloc(sizeof(double)*(N+2)*(N+2));
    double *temp = malloc(sizeof(double)*(N+2)*(N+2));

    memset(result_matrix, 0, sizeof(double)*(N+2)*(N+2));
    memset(temp, 0, sizeof(double)*(N+2)*(N+2));

    MPI_Cart_shift(grid_comm, 0, 1, &neighbours_ranks[WEST], &neighbours_ranks[EAST]);
    MPI_Cart_shift(grid_comm, 1, 1, &neighbours_ranks[SOUTH], &neighbours_ranks[NORTH]);

    double h = 1.0/(n+1);
    
    set_boundaries(result_matrix, N+2, grid_rank, q, coord, h);

    int iter_num = 0;
    int convergence = 0;
    double result, local_diff;
    double eps = 0.0001;

    double io = 0, alg = 0;

    do {

#ifdef BENCHMARK
        MPI_Barrier(MPI_COMM_WORLD);
        double io0 = MPI_Wtime();
#endif
        memcpy(temp, result_matrix, sizeof(double)*(N+2)*(N+2));

#ifdef BENCHMARK
        MPI_Barrier(MPI_COMM_WORLD);
        double io0_e = MPI_Wtime() - io0;

        MPI_Barrier(MPI_COMM_WORLD);
        double alg1 = MPI_Wtime();
#endif
        red_nodes(result_matrix, N+2, coord, q, h, omega);

#ifdef BENCHMARK
        MPI_Barrier(MPI_COMM_WORLD);
        double alg1_e = MPI_Wtime() - alg1;

        MPI_Barrier(MPI_COMM_WORLD);
        double io1 = MPI_Wtime();
#endif
        update_neighbors(result_matrix, N+2, grid_comm, neighbours_ranks, column_type);

#ifdef BENCHMARK
        MPI_Barrier(MPI_COMM_WORLD);
        double io1_e = MPI_Wtime() - io1;

        MPI_Barrier(MPI_COMM_WORLD);
        double alg2 = MPI_Wtime();
#endif
        black_nodes(result_matrix, N+2, coord, q, h, omega);

#ifdef BENCHMARK
        MPI_Barrier(MPI_COMM_WORLD);
        double alg2_e = MPI_Wtime() - alg2;

        MPI_Barrier(MPI_COMM_WORLD);
        double io2 = MPI_Wtime();
#endif
        update_neighbors(result_matrix, N+2, grid_comm, neighbours_ranks, column_type);

#ifdef BENCHMARK
        MPI_Barrier(MPI_COMM_WORLD);
        double io2_e = MPI_Wtime() - io2;
#endif

        local_diff = calc_difference(result_matrix, temp, N+2);
        MPI_Allreduce(&local_diff, &result, 1, MPI_DOUBLE, MPI_MAX, grid_comm);

#ifdef BENCHMARK
        alg2_e = alg2_e + alg1_e;
        io2_e = io2_e + io1_e + io0_e;

        MPI_Reduce(&alg2_e, &alg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&io2_e, &io, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif

#ifdef DEBUG
        if (grid_rank == 0) {
            printf("Iteration: %d, Current difference: %f\n", iter_num, result);
        }
#endif

        if (result < eps) {
            convergence = 1;
#ifdef DEBUG
            if (grid_rank == 0) {
                printf("iter: %d\n", iter_num);
                printf("diff: %f\n", result);
            }
#endif
        }
    } while (!convergence && ++iter_num < MAX_ITER);

#ifdef BENCHMARK
    if (grid_rank == 0) {
        printf("alg: %lf\n", alg);
        printf("io: %lf\n", io);
    }
#endif

    /*
     * Prepare result for saving and
     * remove boundaries
     */
    double *wo_boundaries = malloc(sizeof(double)*N*N);
    for (int i = 1; i <= N; ++i) {
        for (int j = 1; j <= N; ++j) {
            wo_boundaries[(i-1)*N+(j-1)] = result_matrix[i*(N+2)+j];
        }
    }

    FILE *output_file = fopen("output.txt", "w");

    if (output_file  == NULL) {
        puts("ERROR: couldn't open a file and save the matrix");
        goto FINALIZE;
    } else {
        save_matrix(output_file, wo_boundaries, N+2, grid_rank, q, grid_comm);
    }

FINALIZE:

    free(result_matrix);
    free(temp);
    free(wo_boundaries);

    MPI_Finalize();
}

