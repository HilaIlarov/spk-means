#include "spkmeans.h"
#define max_iter 300

/*calloc a 2-D array with "num_of_vectors" cells and each cell with size of "num_of_cells"*/
double **callocation(int num_of_vectors, int num_of_cells)
{
    int i;
    double **res = (double **)calloc(num_of_vectors, sizeof(double *));
    if (res == NULL)
    {
        printf("An Error Has Occured");
        assert(res != NULL);
    }
    for (i = 0; i < num_of_vectors; i++)
    {
        res[i] = (double *)calloc(num_of_cells, sizeof(double));
        if (res[i] == NULL)
        {
            printf("An Error Has Occured");
            assert(res[i] != NULL);
        }
    }
    return res;
}
/*takes v1 = <a1, a2,...,ad> and v2 = <b1, b2,..., bd> and return v = <a1 + b1, a2 + b2,...,ad + bd>*/
double *add_vectors(double *v1, double *v2, int d)
{
    int i = 0;
    double *res;
    res = (double *)calloc(d, sizeof(double));
    if (res == NULL)
    {
        printf("An Error Has Occured");
        assert(res != NULL);
    }
    for (; i < d; i++)
    {
        res[i] = v1[i] + v2[i];
    }
    return res;
}
/*takes v1 = <a1, a2,...,ad> and v2 = <b1, b2,..., bd> and return v = <a1 - b1, a2 - b2,...,ad - bd>*/
double *sub_vectors(double *v1, double *v2, int d)
{
    int i = 0;
    double *res;
    res = (double *)calloc(d, sizeof(double));
    if (res == NULL)
    {
        printf("An Error Has Occured");
        assert(res != NULL);
    }    
    for (; i < d; i++)
    {
        res[i] = v1[i] - v2[i];
    }
    return res;
}
/*takes v = <a1, a2,...,ad> and return v' = <a1 ^ 2, a2 ^ 2,..., ad ^ 2>*/
double pow_vec(double *v, int d)
{
    int i = 0;
    double sum = 0;
    for (; i < d; i++)
    {
        sum += v[i] * v[i];
    }
    return sum;
}
/*takes v = <a1, a2,...,ad> and return v' = <a1 * s, a2 * s,..., ad * s>*/
double *mul_vectors(double *v, double s, int d)
{
    int i = 0;
    double *res;
    res = (double *)calloc(d, sizeof(double));
    if (res == NULL)
    {
        printf("An Error Has Occured");
        assert(res != NULL);
    }    
    for (i = 0; i < d; i++)
    {
        res[i] = v[i] * s;
    }
    return res;
}
/*build vector v_zero = <0, 0,..., 0> : len(v) = d*/
double *build_v_zero(double *v_zero, int d)
{
    int i;
    v_zero = (double *)calloc(d, sizeof(double));
    if (v_zero == NULL)
    {
        printf("An Error Has Occured");
        assert(v_zero != NULL);
    }    
    for (i = 0; i < d; i++)
    {
        v_zero[i] = 0;
    }
    return v_zero;
}
/*build an array of k clusters and initialize for each cluster its fields with size = 0 and sum = v_zero*/
cluster *allocate_s(cluster *s, double *v_zero, int k)
{
    int i = 0;
    s = (cluster *)malloc(k * sizeof(cluster));
    if (s == NULL)
    {
        printf("An Error Has Occured");
        assert(s != NULL);
    }
    for (i = 0; i < k; i++)
    {
        s[i].size = 0;
        s[i].sum = v_zero;
    }
    return s;
}
/*build 2-D array (k * d) and initialize it with x's vectors by order : mu = <x_1, x_2,..., x_k>*/
double **init_mu(double **mu, double **x, int k, int d)
{
    int i;
    int j;
    for (i = 0; i < k; i++)
    {
        for (j = 0; j < d; j++)
        {
            mu[i][j] = x[i][j];
        }
    }
    return mu;
}
/*Calculate the index of the minimal mu that gives the minimal substraction value for x[i]*/
int calc_min_index(double **x, double **mu, int i, int d, int K)
{
    int j = 0;
    double *sub_x = NULL;
    double arg_min = 0;
    int min_index = 0;
    double val;
    sub_x = sub_vectors(x[i], mu[0], d);
    arg_min = pow_vec(sub_x, d);
    min_index = 0;
    for (j = 1; j < K; j++) /* go through all mus and find the minimal one */
    {
        sub_x = sub_vectors(x[i], mu[j], d);
        val = pow_vec(sub_x, d);
        if (val < arg_min)
        {
            arg_min = val;
            min_index = j;
        }
    }
    free(sub_x);
    return min_index;
}
/*update the mu's array*/
void update_mu(double **mu, double *v_zero, cluster *s, int d, int K)
{
    int j;
    for (j = 0; j < K; j++)
    {
        if (s[j].size != 0)
        {
            mu[j] = mul_vectors(s[j].sum, (1.0 / (s[j].size)), d);
        }
        else
        {
            mu[j] = v_zero;
        }
    }
}
/*calculate the final centroids and update mu array*/
double **build_mu(double **x, double **mu, cluster *s, double *v_zero, int N, int K, int d)
{
    int i = 1;
    int change = 1;
    int iter = 0;
    int count_x = 0;
    int min_index = 0;
    double curr_s_index;
    /*for each x in array of x find its match mu  and update S and x d index*/
    while (change && iter < max_iter)
    {
        count_x = 0;
        for (i = 0; i < N; i++) /* for each x in list of x's */
        {
            min_index = calc_min_index(x, mu, i, d, K);
            /* x vector should be in group s with the index of min_index */
            /* update group index*/
            curr_s_index = x[i][d];        /*x current s group */
            if (min_index != curr_s_index) /* if x need to move to another group */
            {
                count_x++; /*there is at least one vector that updated */
                x[i][d] = min_index;
                if (curr_s_index != -1) /* not first transition then remove from current */
                {
                    s[(int)curr_s_index].sum = sub_vectors(s[(int)curr_s_index].sum, x[i], d);
                    s[(int)curr_s_index].size--;
                }
                /* add vector x to new group */
                s[min_index].sum = add_vectors(s[min_index].sum, x[i], d);
                s[min_index].size++;
            }
        }
        if (count_x == 0)
        {
            change = 0;
        }
        update_mu(mu, v_zero, s, d, K);
        iter++;
    }
    return mu;
}
/*print the final centroids*/
void print_mu(double **mu, int k, int d)
{
    int i;
    int j;
    double val;
    for (i = 0; i < k; i++)
    {
        for (j = 0; j < d; j++)
        {
            val = mu[i][j];
            if (val < 0 && val > -0.00005){
                val = 0.0;
            }
            if (j == d - 1)
            {
                printf("%.4f", val);
            }
            else
            {
                printf("%.4f,", val);
            }
        }
        if (i != k - 1)
        {
            printf("\n");
        }
    }
}
/*free functions*/
void free_array(double **arr, int n)
{
    int i;
    for (i = 0; i < n; i++)
    {
        free(arr[i]);
    }
    free(arr);
}
void free_s(cluster *s, int n)
{
    int i;
    for (i = 0; i < n; i++)
    {
        free(s[i].sum);
    }
    free(s);
}
/*k means algorithm*/
void k_means(int k, double **x, int N, int d)
{
    /*init*/
    double *v_zero = NULL;
    double **mu = NULL;
    cluster *s = NULL;
    mu = callocation(k, d);
    /*Initializng mus_array*/
    mu = init_mu(mu, x, k, d);
    /*allocating S*/
    v_zero = build_v_zero(v_zero, d);
    /*allocating S*/
    s = allocate_s(s, v_zero, k);
    /*calc mu*/
    mu = build_mu(x, mu, s, v_zero, N, k, d);
    print_mu(mu, k, d);
    /*free them all*/
    /*free mu, s and v_zero */
    free_array(mu, k);
    free_s(s, k);
    free(v_zero);
    
}
/*get N (= num of vectors in input) and dim (= dimension of vector)*/
int *get_dim_and_N(char *file_name)
{
    char line[1000];
    char *result;
    int i = 1;
    int dim = 1;
    int N = 0;
    FILE *file = fopen(file_name, "r");
    int *d_N = (int *)calloc(2, sizeof(int));
    if (d_N == NULL)
    {
        printf("An Error Has Occured");
        assert(d_N != NULL);
    }
    while ((result = fgets(line, 1000, file)) != NULL)
    {
        if (i == 1)
        {
            while (*result)
            {
                if (*result == ',')
                {
                    dim++;
                }
                ++result;
            }
            i = 0;
        }
        N++;
    }
    d_N[0] = dim;
    d_N[1] = N;
    rewind(file); /*back to start of file*/
    return d_N;
}
/*build X (= vectors in input file)*/
double **x_filling(FILE *file, double **x, int dim)
{
    int i = 0;
    int j = 0;
    char line[1000];
    char *result;
    char *pt;
    while ((result = fgets(line, 1000, file)) != NULL)
    {
        pt = strtok(result, ",");
        while (pt != NULL && j < dim)
        {
            double val = atof(pt);
            x[i][j] = val;
            pt = strtok(NULL, ",");
            j++;
        }
        j = 0;
        i++;
    }
    return x;
}
/*calc the norm of v1 and v2 according to the norm formula*/
double norm(double *v1, double *v2, int d)
{
    int i = 0;
    double res = 0;
    for (; i < d; i++)
    {
        res += pow((v1[i] - v2[i]), 2);
    }
    res = sqrt(res);
    return res;
}
/*calc weights of W matrix*/
double **w_filling(double **w, double **x, int dim, int N)
{
    /*w matrix filling with weights*/
    int i = 0;
    int j = 0;
    double weight = 0;
    /*Initializing W - matrix of weights*/
    w = callocation(N, N);
    for (i = 0; i < N; i++)
    {
        for (j = i; j < N; j++)
        {
            if (i == j)
            {
                w[i][j] = 0;
            }
            else
            {
                weight = exp((-0.5) * norm(x[i], x[j], dim));
                w[i][j] = weight;
                w[j][i] = weight;
            }
        }
    }
    return w;
}
/*calc sum(v[i]) for all 1<=i<=d when v = <v1, v2,...,vd>*/
double sum_coordinates(double *v, int d)
{
    int i = 0;
    double res = 0;
    for (; i < d; i++)
    {
        res += v[i];
    }
    return res;
}
/*calc Diagonal matrix based on W*/
double *d_filling(double *d, double **w, int N)
{
    int i = 0;
    /*Initializing D - diagonal matrix*/
    d = (double *)calloc(N, sizeof(double));
    if (d == NULL)
    {
        printf("An Error Has Occured");
        assert(d != NULL);
    }
    /*D matrix filling with weights sum*/
    for (; i < N; i++){
        d[i] = sum_coordinates(w[i], N);
    }
    return d;
}
/*calc D^(-0.5) based on D*/
double *d_minus_filling(double *d_minus, double *d, int N)
{
    int i = 0;
    /*Initializing D sqrt -0.5 - diagonal matrix*/
    d_minus = (double *)calloc(N, sizeof(double));        
    if (d_minus == NULL)
    {
        printf("An Error Has Occured");
        assert(d_minus != NULL);
    }    
    /*D ^ -0.5 calculation and filling*/
    for (; i < N; i++){
        d_minus[i] = 1.0 / sqrt(d[i]);
    }
    return d_minus;
}
/*compute D * W*/
double **mul_diagonal_matrix_left(double *d, double **w, int n)
{
    int i;
    int j;
    double **res = (double **)calloc(n, sizeof(double *));
    if (res == NULL)
    {
        printf("An Error Has Occured");
        assert(res != NULL);
    }    
    /*init*/
    for (i = 0; i < n; i++)
    {
        res[i] = (double *)calloc(n, sizeof(double));
        if (res[i] == NULL)
        {
            printf("An Error Has Occured");
            assert(res[i] != NULL);
        }
    }
    /*mul*/
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            res[i][j] = d[i] * w[i][j];
        }
    }
    return res;
}
/*compute W * D*/
double **mul_diagonal_matrix_right(double **w, double *d, int n)
{
    int i;
    int j;
    double **res = (double **)calloc(n, sizeof(double *));
    if (res == NULL)
    {
        printf("An Error Has Occured");
        assert(res != NULL);
    }    
    /*init*/
    for (i = 0; i < n; i++)
    {
        res[i] = (double *)calloc(n, sizeof(double));
        if (res[i] == NULL)
        {
            printf("An Error Has Occured");
            assert(res[i] != NULL);
        }    
    }
    /*mul*/
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            res[i][j] = d[i] * w[j][i];
        }
    }
    return res;
}
/*build and return matrix m (= I - mat)*/
double **sub_matrix_from_id(double **mat, int n)
{
    int i = 0;
    int j = 0;
    double **m = (double **)calloc(n, sizeof(double *));
    if (m == NULL)
    {
        printf("An Error Has Occured");
        assert(m != NULL);
    }    
    for (; i < n; i++)
    {
        m[i] = (double *)calloc(n, sizeof(double));
        if (m[i] == NULL)
        {
            printf("An Error Has Occured");
            assert(m[i] != NULL);
        }
    }
    for (i = 0; i < n; i++)
    {  
        for (j = 0; j < n; j++)
        {
            if (i == j)
            {
                m[i][j] = 1.0 - mat[i][j];    
            }
            else
            {
                m[i][j] = 0.0 - mat[i][j];    
            }
        }
    }
    return m;
}
/*build lnorm based on W and D^(-0.5)*/
double **form_l_norm(double **w, double *d_minus, int N)
{
    double **l_left;
    double **l_right;
    double **l_norm;
    l_left = mul_diagonal_matrix_left(d_minus, w, N);
    l_right = mul_diagonal_matrix_right(l_left, d_minus, N);
    l_norm = sub_matrix_from_id(l_right, N);
    free_array(l_left, N);
    free_array(l_right, N);
    return l_norm;
}
/*find Aij, the off-diagonal element with the largest absolute value and return i,j in an array*/
double *max_off_diagonal_indexes(double **l_norm, int N)
{
    int i = 0;
    int j = i;
    double abs_val = 0;
    double max = fabs(l_norm[0][1]);
    double *max_index = (double *)calloc(2, sizeof(double));
    if (max_index == NULL)
    {
        printf("An Error Has Occured");
        assert(max_index != NULL);
    }    
    max_index[0] = 0;
    max_index[1] = 1;
    for (; i < N; i++)
    {
        for (j = i + 1; j < N; j++)
        {
            abs_val = fabs(l_norm[i][j]);
            if (max < abs_val)
            {
                max = abs_val;
                max_index[0] = i;
                max_index[1] = j;
            }
        }
    }
    return max_index;
}
/*calc t based on teta*/
double calc_t(double teta)
{
    double t = 1.0 / (fabs(teta) + sqrt(pow(teta, 2) + 1));
    if (teta < 0)
    {
        t *= -1;
    }
    return t;
}
/*build an array ijcs = [i, j, c ,s] with the function mentioned above*/
double *calc_max_i_j_and_c_and_s(double **l_norm, int N)
{
    double teta;
    double t = 0;
    double c = 0;
    double s = 0;
    int i = 0;
    int j = 0;
    double *tlav_indexes = max_off_diagonal_indexes(l_norm, N);
    double *ijcs = (double *)calloc(4, sizeof(double));
    if (ijcs == NULL)
    {
        printf("An Error Has Occured");
        assert(ijcs != NULL);
    }
    i = (int)tlav_indexes[0];
    j = (int)tlav_indexes[1];
    free(tlav_indexes);
    teta = (l_norm[j][j] - l_norm[i][i]) / (2.0 * l_norm[i][j]);
    t = calc_t(teta);
    c = 1.0 / (sqrt(pow(t, 2) + 1));
    s = t * c;
    ijcs[0] = i;
    ijcs[1] = j;
    ijcs[2] = c;
    ijcs[3] = s;
    return ijcs;
}
/*build the first Rotation matrix P and put V = P*/
double **initial_V(double **V, double *ijcs, int N)
{
    int i = (int)ijcs[0];
    int j = (int)ijcs[1];
    double c = ijcs[2];
    double s = ijcs[3];
    int k = 0;
    for (; k < N; k++)
    {
        V[k][k] = 1;
    }
    V[i][i] = c;
    V[j][j] = c;
    V[i][j] = s;
    V[j][i] = (-1) * s;
    return V;
}
/*calc sum of squares of all off-diagonal elements of m*/
double sum_off(double **m, int N)
{
    int i = 0;
    int j = 0;
    double res = 0;
    for (; i < N; i++)
    {
        for (j = i + 1; j < N; j++)
        {
            res += pow(m[i][j], 2);
        }
    }
    res *= 2;
    return res;
}
/*build A_tag matrix based on the relation between A and A_tag*/
double **create_A_tag(double **A, double **A_tag, double *ijcs, int N)
{
    int i = (int)(ijcs[0]);
    int j = (int)(ijcs[1]);
    double c = ijcs[2];
    double s = ijcs[3];
    int r;
    double a_val_r_i;
    double a_val_r_j;
    double a_val_i_i;
    a_val_i_i = A[i][i];

    for (r = 0; r < N; r++)
    {
        if (r != i && r != j)
        {
            a_val_r_i = A[r][i];
            a_val_r_j = A[r][j];
            A_tag[r][i] = c * a_val_r_i - s * a_val_r_j;
            A_tag[i][r] = A_tag[r][i];
            A_tag[r][j] = c * a_val_r_j + s * a_val_r_i;
            A_tag[j][r] = A_tag[r][j];
        }
    }
    A_tag[i][i] = c * c * a_val_i_i + s * s * A[j][j] - 2 * s * c * A[i][j];
    A_tag[j][j] = s * s * a_val_i_i + c * c * A[j][j] + 2 * s * c * A[i][j];
    A_tag[i][j] = 0;
    A_tag[j][i] = 0;
    return A_tag;
}
/*compute V * P when P is the rotation matrix*/
double **v_mul_p(double **m, double *ijcs, int N)
{
    int i = (int)(ijcs[0]);
    int j = (int)ijcs[1];
    double c = ijcs[2];
    double s = ijcs[3];
    int r;
    double m_r_i;
    double m_r_j;
    for (r = 0; r < N; r++)
    {
        m_r_i = m[r][i];
        m_r_j = m[r][j];
        m[r][i] = m_r_i * c + m_r_j * (-s);
        m[r][j] = m_r_i * s + m_r_j * c;
    }
    return m;
}
/*calculate A_tag matrix and V = P1* P2 * P3... with the functions above*/
double **get_V_and_A_tag(double **V, double **l_norm, double **A, double **A_tag, int N)
{
    double off_A;
    double off_A_tag;
    int convergence;
    int count = 1;
    const double epsilon = pow(10, -15);
    double *ijcs;
    V = callocation(N, N);
    A = l_norm;
    /*calc max i and j and get c and s*/
    ijcs = calc_max_i_j_and_c_and_s(l_norm, N);
    /*V = P1*/
    V = initial_V(V, ijcs, N);
    off_A = sum_off(A, N);
    A_tag = A;
    A_tag = create_A_tag(A, A_tag, ijcs, N);
    off_A_tag = sum_off(A_tag, N);
    convergence = off_A - off_A_tag <= epsilon;
    while (convergence != 1 && count < 100)
    {
        ijcs = calc_max_i_j_and_c_and_s(A_tag, N);
        V = v_mul_p(V, ijcs, N);
        A = A_tag;
        off_A = off_A_tag;
        A_tag = create_A_tag(A, A_tag, ijcs, N);
        off_A_tag = sum_off(A_tag, N);
        convergence = off_A - off_A_tag <= epsilon;
        count++;
    }
    free(ijcs);
    return V;
}
/*find eigenvalues (= the diagonal elements of A_tag) and match each one its own index by order and return as struct array*/
eigenvalue_index *get_eigenvalues_and_indexes(double **A_tag, eigenvalue_index *index_ev, int N)
{
    int i = 0;
    index_ev = (eigenvalue_index *)malloc(N * sizeof(eigenvalue_index));
    if (index_ev == NULL)
    {
        printf("An Error Has Occured");
        assert(index_ev != NULL);
    }
    for (; i < N; i++)
    {
        index_ev[i].index = i;
        index_ev[i].eigenvalue = A_tag[i][i];
    }
    return index_ev;
}
/*comparator function in order to sort the eigevalues*/
int cmpfunc(const void *a, const void *b)
{
    eigenvalue_index *A = (eigenvalue_index *)a;
    eigenvalue_index *B = (eigenvalue_index *)b;
    if (A->eigenvalue < B->eigenvalue)
    {
        return -1;
    }
    else
    {
        if (A->eigenvalue == B->eigenvalue)
        {
            if (A->index < B->index)
            {
                return -1;
            }
        }
        return 1;
    }
}
/*calc k = (argmax_i(delta_i), i = 1,...,floor(N/2))*/
int get_k(eigenvalue_index *index_ev, int N)
{
    int r;
    double delta;
    double max_delta = -1;
    int max_index = -1;
    for (r = 1; r <= floor(N / 2); r++)
    {
        delta = fabs(index_ev[r].eigenvalue - index_ev[r - 1].eigenvalue);
        if (delta > max_delta)
        {
            max_delta = delta;
            max_index = r;
        }
    }
    return max_index;
}
/*build U matrix from V and according to the sorted eigenvalues*/
double **build_U(double **U, double **V, eigenvalue_index *index_ev, int k, int N)
{
    int i = 0;
    int j = 0;
    int index;
    U = callocation(N, k);
    for (; i < k; i++)
    {
        index = index_ev[i].index;
        for (j = 0; j < N; j++)
        {
            U[j][i] = V[j][index];
        }
    }
    return U;
}
/*Build T matrix from U*/
double **build_T(double **U, double **T, int k, int N, int from_c)
{
    int i = 0;
    int j = 0;
    int l = 0;
    double sum = 0;
    if(from_c == 1)
    {
        T = callocation(N + 1, k + 1);
    }
    else
    {
        T = callocation(N + 1, k);
    }
    for (j = 0; j < N; j++)
    {
        for (i = 0; i < k; i++)
        {
            sum += pow(U[j][i], 2);
        }
        for (l = 0; l < k; l++)
        {
            T[j][l] = U[j][l] / (sqrt(sum));
        }
        if(from_c == 1)
        {
            T[j][k] = -1;
        }
        sum = 0;
    }
    T[N][0] = k;
    return T;
}
/*form spk process*/
double **spk(int k, char *goal, double **x, int N, int dim, int from_c)
{
    /*initializing*/
    double *d = NULL;
    double *d_minus_half = NULL;
    double **w = NULL;
    double **l_norm = NULL;
    double **A = NULL;
    double **A_tag = NULL;
    double **V = NULL;
    double **U = NULL;
    double **T = NULL;
    eigenvalue_index *index_ev = NULL;
    /*w matrix filling with weights*/
    w = w_filling(w, x, dim, N);
    if (strcmp(goal, "wam") == 0)
    {
        print_matrix(w, N, N);
        free_array(w, N);
        free_array(x, N);
        exit(0);
    }
    /*d matrix filling*/
    d = d_filling(d, w, N);
    if (strcmp(goal, "ddg") == 0)
    {
        print_arr_as_diagonal_matrix(d, N);
        free(d);
        free_array(w, N);
        free_array(x, N);
        exit(0);
    }
    /*d_minus_half matrix filling*/
    d_minus_half = d_minus_filling(d_minus_half, d, N);
    /*l_norm matrix filling*/
    l_norm = form_l_norm(w, d_minus_half, N);
    free_array(w, N);
    free(d);
    free(d_minus_half);
    if (strcmp(goal, "lnorm") == 0)
    {
        print_matrix(l_norm, N, N);
        free_array(l_norm, N);
        free_array(x, N);
        exit(0);
    }
    if (strcmp(goal, "jacobi") == 0)
    {
        l_norm = x;
    }
    A_tag = l_norm;
    /*Inistializing V */
    V = get_V_and_A_tag(V, l_norm, A, A_tag, N);
    index_ev = get_eigenvalues_and_indexes(A_tag, index_ev, N);
    free_array(l_norm, N);
    if (strcmp(goal, "jacobi") == 0)
    {
        print_eigenvalues(index_ev, N);
        print_tranpose(V, N, N);
        free_array(V, N);
        free(index_ev);
        exit(0);
    }
    /*Sort the eigenvalues*/
    qsort(index_ev, N, sizeof(eigenvalue_index), cmpfunc);
    if (k == 0)
    {
        k = get_k(index_ev, N);
    }
    /*Build U and T*/
    U = build_U(U, V, index_ev, k, N);
    T = build_T(U, T, k, N, from_c);
    free(index_ev);
    free_array(V, N);
    free_array(U, N);
    /*Run kmeans (HW1)*/
    if (from_c == 1)
    {
        k_means(k, T, N, k);
        free_array(T, N + 1);
        free_array(x, N);
        exit(0);
    }
    return T;
}
/*print a diagonal matrix, of which the diagonal values are given in the array*/
void print_arr_as_diagonal_matrix(double *arr, int n)
{
    int i = 0;
    int j = 0;
    double val = 0;
    for (; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            val = 0;
            if(i == j){
                val = arr[i];
            }
            if (val < 0 && val > -0.00005){
                val = 0.0;
            }
            if (j < n - 1)
            {
                printf("%.4f,", val);
            }
            else
            {
                printf("%.4f", val);
            }
        }
        if (i != n-1)
        {
            printf("\n");
        }
    }
}
/*print matrix function*/
void print_matrix(double **m, int n, int d)
{
    int i = 0;
    int j = 0;
    double val;
    for (; i < n; i++)
    {
        for (j = 0; j < d; j++)
        {
            val = m[i][j];
            if (val < 0 && val > -0.00005){
                val = 0.0;
            }
            if (j < d - 1)
            {
                printf("%.4f,", val);
            }
            else
            {
                printf("%.4f", val);
            }
        }
        if (i != n-1)
        {
            printf("\n");
        }
    }
}
/*print the transpose of a matrix*/
void print_tranpose(double **m, int num_of_vectors, int num_of_cells)
{
    int i = 0;
    int j = 0;
    double val;
    for (; i < num_of_cells; i++)
    {
        for (; j < num_of_vectors; j++)
        {
            val = m[j][i];
            if (val < 0 && val > -0.00005){
                val = 0.0;
            }            
            if (j < num_of_vectors - 1)
            {
                printf("%.4f,", val);
            }
            else
            {
                printf("%.4f", val);
            }
        }
        if (i != num_of_cells - 1)
        {
            printf("\n");
        }
        j = 0;
    }
}
/*print the eigenvalues*/
void print_eigenvalues(eigenvalue_index *index_ev, int N)
{
    int i = 0;
    double val;
    for (; i < N; i++)
    {
        val = index_ev[i].eigenvalue;
        if (val < 0 && val > -0.00005){
            val = 0.0;
        } 
        if (i != N - 1)
        {
            printf("%.4f,", val);
        }
        else
        {
            printf("%.4f", val);
        }
    }
    printf("\n");
}
/*The main function*/
int main(int argc, char *argv[])
{
    int k = 0;
    double **x;
    int *d_N;
    char *file_name; 
    char *goal;      
    int dim;
    int N;
    FILE *file;
    /*We assumed that the given k is a string that represents an integer*/
    k = atoi(argv[1]);
    goal = argv[2];
    file_name = argv[3];
    file = fopen(file_name, "r");
    if (file == NULL)
    {
        printf("An Error Has Occured");
        assert(file != NULL);
    }
    argc += 0;
    d_N = get_dim_and_N(file_name);
    dim = d_N[0];
    N = d_N[1];
    free(d_N);
    x = callocation(N, dim);
    /*x_array filling with vectors*/
    x = x_filling(file, x, dim);
    spk(k, goal, x, N, dim, 1);
    return 0;
}