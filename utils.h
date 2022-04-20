#ifndef UTILS_H
#define UTILS_H

void set_boundaries(double *x, int x_size, int rank,  int rank_size, int *coord, double h);
void update_neighbors(double *x, int x_size, MPI_Comm grid_comm, int *neighbours, MPI_Datatype col_type);
void red_nodes(double *x, int x_size, int *grid_coord, int q, double h, double omega);
void black_nodes(double *x, int x_size, int *grid_coord, int q, double h, double omega);
void save_matrix(FILE *f, double *mat, int x_size, int rank, int q, MPI_Comm grid_comm);
double calc_difference(double *x, double *temp, int x_size);

#endif /* UTILS_H */
