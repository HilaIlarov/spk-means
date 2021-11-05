#ifndef SPKMEANS_H_
#define SPKMEANS_H_
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <math.h>

typedef struct
{
    double size;
    double *sum;
} cluster;

typedef struct
{
    int index;
    double eigenvalue;
} eigenvalue_index;

double **callocation(int num_of_vectors, int num_of_cells);

/*k_means aid functions*/
double *add_vectors(double *v1, double *v2, int d);
double *sub_vectors(double *v1, double *v2, int d);
double pow_vec(double *v, int d);
double *mul_vectors(double *v, double s, int d);
double *build_v_zero(double *v_zero, int d);
cluster *allocate_s(cluster *s, double *v_zero, int k);
double **init_mu(double **mu, double **x, int k, int d);
int calc_min_index(double **x, double **mu, int i, int d, int K);
void update_mu(double **mu, double *v_zero, cluster *s, int d, int K);
double **build_mu(double **x, double **mu, cluster *s, double *v_zero, int N, int K, int d);
void print_mu(double **mu, int k, int d);
void free_array(double **arr, int n);
void free_s(cluster *s, int n);
void k_means(int k, double **x, int N, int d);

/*input aid functions*/
int *get_dim_and_N(char *file_name);
double **x_filling(FILE *file, double **x, int dim);

/*spk aid functions*/
double norm(double *v1, double *v2, int d);
double **w_filling(double **w, double **x, int dim, int N);

double sum_coordinates(double *v, int d);
double *d_filling(double *d, double **w, int N);
double *d_minus_filling(double *d_minus, double *d, int N);

double **mul_diagonal_matrix_left(double *d, double **w, int n);
double **mul_diagonal_matrix_right(double **w, double *d, int n);
double **sub_matrix_from_id(double **mat, int n);
double **form_l_norm(double **w, double *d_minus, int N);

double *max_off_diagonal_indexes(double **l_norm, int N);
double calc_t(double teta);
double *calc_max_i_j_and_c_and_s(double **l_norm, int N);
double **initial_V(double **V, double *ijcs, int N);
double sum_off(double **m, int N);
double **create_A_tag(double **A, double **A_tag, double *ijcs, int N);
double **v_mul_p(double **m, double *ijcs, int N);
double **get_V_and_A_tag(double **V, double **l_norm, double **A, double **A_tag, int N);

eigenvalue_index *get_eigenvalues_and_indexes(double **A_tag, eigenvalue_index *index_ev, int N);
int cmpfunc(const void *a, const void *b);
int get_k(eigenvalue_index *index_ev, int N);

double **build_U(double **U, double **V, eigenvalue_index *index_ev, int k, int N);
double **build_T(double **U, double **T, int k, int N, int from_c);

double **spk(int k, char *goal, double **vectors, int N, int dim, int from_c);

void print_arr_as_diagonal_matrix(double *arr, int n);
void print_matrix(double **m, int n, int d);
void print_tranpose(double **m, int n, int d);
void print_eigenvalues(eigenvalue_index *index_ev, int N);

#endif
