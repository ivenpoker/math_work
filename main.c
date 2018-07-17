 /*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <ctype.h>

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

#define N 1000
#define FILE_NAME_LEN 1000
#define ARRAY_SIZE(arr)  (sizeof(arr)/sizeof(*(arr)))


FILE *input_file;
FILE *output_file;

typedef struct Matrix {
    int num_row;
    int num_col;
    double *matrix[N];
} MATRIX_t, *MATRIX_p_t;

typedef struct Equations {
    int order;
    MATRIX_p_t matrix;
    MATRIX_p_t vector;
} EQUATIONS_t, *EQUATIONS_p_t;

// this is for the inverse of a matrix

double determinant(double **, double);
MATRIX_p_t cofactor(double **, double);
MATRIX_p_t transpose(double **, double **, double);

int factorial(int n);
int comb(int n,int r);
void print_option(int n);
int display_menu(void);
void get_file_name(char [], int n);
int read_line(char str[], int n);
void terminate_program(char *);
bool initialize_matrices(int size, double **);
bool matrix_initializer(MATRIX_p_t);
bool matrix_filler(MATRIX_p_t, MATRIX_p_t);
void fill_matrix(MATRIX_p_t);
void display_matrix(int, double **);
void print_matrix(MATRIX_p_t);
MATRIX_p_t multiply_matrix(MATRIX_p_t, MATRIX_p_t);
MATRIX_p_t matrix_inverse(MATRIX_p_t);
void validate_matrix(MATRIX_p_t);
void validate_equation(EQUATIONS_p_t);

// for the cases

bool store_data(int, MATRIX_p_t, MATRIX_p_t, int, char []);
EQUATIONS_p_t get_hilbert_matrix(int len, char []);
EQUATIONS_p_t read_file_data(int len, char file_name[]);
void A_inverse_x_b(EQUATIONS_p_t);

// display cases
void display_case(int, char [], int);

// Mostly used for Crammer's rule

void crammers_rule(EQUATIONS_p_t);
void vector_substitute(EQUATIONS_p_t, int column);
EQUATIONS_p_t make_copy(EQUATIONS_p_t);

// For Gaussian elimination

void gaussian_elimination(EQUATIONS_p_t);
int forward_elim(int n, int m, double matrix[n][m]);
void swap_row(int n, int m, double matrix[][n+1], int, int);
void back_substitution(int n, int m, double matrix[n][m]);

// This is to Compute [(A inverse) x A] and [A x (A inverse)] and report.

void compare_inverse_multiplications(EQUATIONS_p_t equations, int , char []);


// This for file opening

EQUATIONS_p_t open_file(int, char []);


/**
 * Main function begins the program execution
 * @param argc Number of arguments to pass to functions
 * @param argv Array of <code>char *</code> strings to pass to functions
 * @return 0 for successful termination else 1.
 */
int main(int argc, char** argv) {
    /*
    int size;
    display_menu();
    printf("\tEnter size of matrices: ");
    scanf("%d", &size);

    double *matrix_A[size], *matrix_B[size], *matrix_C[size];

    initialize_matrices(size, matrix_A);
    initialize_matrices(size, matrix_B);
    initialize_matrices(size, matrix_C);

    printf("\n\tFilling Matrix A [%d x %d]\n", size, size);
    fill_matrix(size, matrix_A);
    // printf("\n\n\tFilling Matrix B [%d x %d]\n", size, size);
    // fill_matrix(size, matrix_B);
    // printf("\n\n\tFilling Matrix A [%d x %d]\n", size, size);
    // fill_matrix(size, matrix_C);

    printf("\n\n\tDisplay Matrix A\n");
    display_matrix(size, matrix_A);
     */

    // #######################################################

    /*
    display_menu();

    MATRIX_p_t new_matrix = (MATRIX_p_t)malloc(sizeof(MATRIX_t));

    double order, det;

    initialize_matrices(N, new_matrix->matrix);

    printf("\n\nEnter the order of the matrix: ");
    scanf("%f", &order);
    new_matrix->num_row = order; new_matrix->num_col = order;
    printf("Enter the elements of (%.0f x %.0f) matrix: ", order, order);
    fill_matrix(order, new_matrix->matrix);

    det = determinant(new_matrix->matrix, order);
    if (det == 0)
        printf("\nInverse of Entered matrix is not possible\n");
    else
        cofactor(new_matrix->matrix, order);
    */

    // ######################## Multiplication test ###########################

    /*
    EQUATIONS_p_t equation = (EQUATIONS_p_t)malloc(sizeof(EQUATIONS_t));

    if (equation == NULL) return -1;

    equation->matrix = (MATRIX_p_t)malloc(sizeof(MATRIX_t));
    equation->vector = (MATRIX_p_t)malloc(sizeof(MATRIX_t));
    int row_1, col_1, row_2, col_2;

    printf("Enter number of rows for M1: ");
    scanf("%d", &(equation->matrix->num_row));
    printf("Enter number of cols for M1: ");
    scanf("%d", &(equation->matrix->num_col));

    equation->order = equation->matrix->num_row;

    equation->vector->num_row = equation->matrix->num_row;
    equation->vector->num_col = 1;


    printf("\n\n");

    printf("Fill the [%d x %d] matrix A: ", equation->matrix->num_row, equation->matrix->num_col);
    matrix_initializer(equation->matrix);
    fill_matrix(equation->matrix);

    printf("Fill the [%d x %d] vector B: ", equation->vector->num_row, equation->vector->num_col);
    matrix_initializer(equation->vector);
    fill_matrix(equation->vector);

    printf("\n\nUsing Gaussian elimination\n\n");

    gaussian_elimination(equation);
     */




    char ans;
    do {
        system("cls");
        char file_name[FILE_NAME_LEN];
        int choice = display_menu();

        switch(choice) {
            case 1:
                print_option(choice);
                get_hilbert_matrix(FILE_NAME_LEN, file_name);
                break;
            case 2:
                print_option(choice);
                display_case(FILE_NAME_LEN, file_name, 2);
                break;
            case 3:
                print_option(choice);
                display_case(FILE_NAME_LEN, file_name, 3);
                break;
            case 4:
                print_option(choice);
                display_case(FILE_NAME_LEN, file_name, 4);
                break;
            case 5:
                print_option(choice);
                display_case(FILE_NAME_LEN, file_name, 5);
                break;
            case 6:
                print_option(choice);
                break;
            case 7:
                print_option(choice);
                display_case(FILE_NAME_LEN, file_name, 7);
                break;
            case 8:
                printf("\n\n");
                exit(EXIT_SUCCESS);
            default:
                printf("\n\t[============== invalid choice ==========]");
        }
        printf("\n\t\tDo you want to continue? (y/n): ");
        scanf(" %c", &ans);
    } while (ans == 'Y' || ans == 'y');

    /*
    char file_name[FILE_NAME_LEN];
    printf("Enter file name with data: ");
    read_line(file_name, FILE_NAME_LEN);

    EQUATIONS_p_t results = read_file_data(FILE_NAME_LEN, file_name);

    printf("The order of matrix: %d", results->order);

    printf("The matrix is: \n");
  //  pri(results->matrix->num_row, results->matrix->matrix);
    print_matrix(results->matrix);

    printf("The vector is: \n");
    print_matrix(results->vector);

    printf("\n\n");
     */

    printf("\n\n");
    return (EXIT_SUCCESS);
}

/**
 * Computes the determinant of a matrix
 * @param main_arr Main matrix to compute it's determinant
 * @param order Order of the matrix
 * @return determinant as a double
 */
double determinant(double *main_arr[N], double order) {
    double s = 1, det = 0, *tmp_arr[N];
    initialize_matrices(N, tmp_arr);
    int c, i, j;
    int m, n;
    if (order == 1) {
        return (main_arr[0][0]);
    } else {
        det = 0;
        for (c = 0; c < order; c++) {
            m = 0; n = 0;
            for (i = 0; i < order; i++) {
                for (j = 0; j < order; j++) {
                    tmp_arr[i][j] = 0;
                    if (i != 0 && j != c) {
                        tmp_arr[m][n] = main_arr[i][j];
                        if (n < (order-2)) n++;
                        else {
                            n = 0;
                            m++;
                        }
                    }
                } // end of inner for loop (with j)
            } // end of inner for loop (with i)
             det = det + s * (main_arr[0][c] * determinant(tmp_arr, order-1));
             s = -1 * s;
        } // end of inner for loop (with c)
    }
    return (det);
}

/**
 * Computes the inverse of a Hilbert matrix
 * @param some_matrix Hilbert matrix to compute the inverse
 * @return A pointer to the dynamically allocated matrix (the inverse)
 */
MATRIX_p_t matrix_inverse(MATRIX_p_t some_matrix) {
    int i, j;
    if (some_matrix == NULL)
        terminate_program("[ERROR]: invalid matrix");
    if (some_matrix->num_row != some_matrix->num_col)
        printf("\n\t\t[------- Only square matrices are invertible ------]");
    for(i=1; i<= some_matrix->num_row; i++){
      for(j=1; j<= some_matrix->num_col; j++){
        some_matrix->matrix[i-1][j-1] = pow(-1, i+j) * (i+j-1) *
                comb(some_matrix->num_row+i-1, some_matrix->num_row-j) *
                comb(some_matrix->num_row+j-1, some_matrix->num_row-i) *
                pow(comb(i+j-2, i-1), 2);
      }
    }
    return some_matrix;
}

int factorial(int n) {
    if (n < 0)
        printf("Error! Factorial of a negative number doesn't exist.");
    else {
        int fact = 1;
        int i;
        for(i=1; i<=n; ++i)
            fact *= i;              // fact = fact*i;
        return fact;
    }
}

int comb(int n,int r) {
  return factorial(n)/(factorial(n-r)*factorial(r));
}

/**
 * Computes the Cofactor of a Matrix
 * @param arr Matrix to compute cofactors
 * @param order Order of the matrix.
 */
MATRIX_p_t cofactor(double *arr[N], double order) {
    double *tmp_arr[N], *fac[N];
    initialize_matrices(N, tmp_arr);
    initialize_matrices(N, fac);

    int m, n, p, i, j, q;

    for (q = 0; q < order; q++) {
        for (p = 0; p < order; p++) {
            m = 0; n = 0;
            for (i = 0; i < order; i++) {
                for (j = 0; j < order; j++) {
                    if (i != q && j != p) {
                        tmp_arr[m][n] = arr[i][j];
                        if (n < (order - 2)) n++;
                        else {
                            n = 0;
                            m++;
                        } // end of else
                    } // end if outer if
                } // end of for loop (with j)
            } // end of for loop (with i)
            fac[q][p] = pow(-1, q+p) * determinant(tmp_arr, order-1);
        }
    }
    return transpose(arr, fac, order);
}

/**
 * Computes the transpose of a matrix
 * @param arr Matrix to compute it's transpose
 * @param fac Another matrix for internal computations
 * @param order Order of the matrix
 */
MATRIX_p_t transpose(double *arr[N], double *fac[N], double order) {
    double *tmp_arr[N], *inverse[N], det;

    initialize_matrices(N, tmp_arr);
    initialize_matrices(N, inverse);

    int i, j;

    for (i = 0; i < order; i++) {
        for (j = 0; j < order; j++) {
            tmp_arr[i][j] = fac[j][i];
        }
    }

    det = determinant(arr, order);
    for (i = 0; i < order; i++) {
        for (j = 0; j < order; j++) {
            inverse[i][j] = tmp_arr[i][j] / det;
        }
    }

    // This is the inverse of the matrix

   // printf("\n\nThe inverse of the matrix is: \n");

    MATRIX_p_t inverse_matrix = (MATRIX_p_t)malloc(sizeof(MATRIX_t));
    inverse_matrix->num_row = order;
    inverse_matrix->num_col = order;
    matrix_initializer(inverse_matrix);

    for (i = 0; i < inverse_matrix->num_row; i++)
        for (j = 0; j < inverse_matrix->num_col; j++)
            inverse_matrix->matrix[i][j] = inverse[i][j];

    /*
    for (int i = 0; i < order; i++) {
        for (int j = 0; j < order; j++) {
            printf("\t%f", inverse[i][j]);
        }
        printf("\n");
    }
    */
    return inverse_matrix;
}
/**
 * Completes the creation of a matrix
 * @param size Size of the matrix
 * @param matrix Matrix to create space for
 * @return <code>true</code> on success else <code>false</code>
 */
bool initialize_matrices(int size, double *matrix[size]) {
    int i;
    for (i = 0; i < size; i++) {
        matrix[i] = (double *)malloc(sizeof(double) * size);
        if (matrix[i] == NULL)
            terminate_program("ERROR: Cannot allocate memory for matrix");
    }

}

bool matrix_initializer(MATRIX_p_t init_matrix) {
    int i, j;
    if (init_matrix == NULL)
        terminate_program("[ERROR]: invalid matrix");
    for (i = 0; i < init_matrix->num_row; i++) {
        init_matrix->matrix[i] = (double *)malloc(sizeof(double) * init_matrix->num_col);
        if (init_matrix->matrix[i] == NULL)
            terminate_program("[ERROR]: cannot initialize matrix");
        for (j = 0; j < init_matrix->num_col; j++)
            init_matrix->matrix[i][j] = 0;
    }
    /*
    for (int i = 0; i < init_matrix->num_row; i++) {
        for (int j = 0; j < init_matrix->num_row; j++)
            init_matrix->matrix[i][j] = 0;
    }
    */
}

/**
 * Fills a matrix by using the contents of the standard output and input stream
 * @param size Size of the matrix
 * @param matrix Matrix to fill.
 */
void fill_matrix(MATRIX_p_t some_matrix) {
    printf("\n");
    int i, j;
    for (i = 0; i < some_matrix->num_row; i++) {
        printf("\tEnter the %d value(s) for row %d [space value(s)]: ", some_matrix->num_col, i+1);
        for (j = 0; j < some_matrix->num_col; j++)
            scanf("%lf", &some_matrix->matrix[i][j]);
    }
}

void validate_matrix(MATRIX_p_t some_matrix) {
    if (some_matrix == NULL) {
         printf("\n\t[============= [ERROR: Insufficient memory =========]");
         exit(EXIT_FAILURE);
    }
}

EQUATIONS_p_t get_hilbert_matrix(int len, char file_name[len]) {
    /*
    int order = -1;

    printf("\tEnter order 'n' of Hilbert matrix (n such that, order = n x n): ");
    scanf("%d", &order);

    MATRIX_p_t new_matrix = (MATRIX_p_t)malloc(sizeof(MATRIX_t));
    validate_matrix(new_matrix);

    new_matrix->num_row = order;
    new_matrix->num_col = order;

    matrix_initializer(new_matrix);
    fill_matrix(new_matrix);

    MATRIX_p_t vector = (MATRIX_p_t)malloc(sizeof(MATRIX_t));
    validate_matrix(vector);

    vector->num_row = order;
    vector->num_col = 1;

    printf("\n\tEnter %d x 1 column vector b: ", order);

    matrix_initializer(vector);
    fill_matrix(vector);

    if (store_data(order, new_matrix, vector, len, file_name)) {
        printf("\n\t\t%s[ Data stored in file \"%s\" ]%s\n\n",
                ANSI_COLOR_GREEN, file_name, ANSI_COLOR_RESET);
    } else {
        printf("\n\t\t%sERROR storing data in file \"%s\"%s",
                ANSI_COLOR_GREEN, file_name, ANSI_COLOR_RESET);
    }
    */

    int order = -1;

    printf("\tEnter order 'n' of Hilbert matrix (n such that, order = n x n): ");
    scanf("%d", &order);

    MATRIX_p_t new_matrix = (MATRIX_p_t)malloc(sizeof(MATRIX_t));
    validate_matrix(new_matrix);

    new_matrix->num_row = order;
    new_matrix->num_col = order;

    matrix_initializer(new_matrix);

    int i, j;

    for (i = 1; i <= order; i++) {
        for (j = 1; j <= order; j++) {
            double test = (double) 1 / (i+j-1);
            new_matrix->matrix[i-1][j-1] = test;
        }
    }

    MATRIX_p_t vector = (MATRIX_p_t)malloc(sizeof(MATRIX_t));
    validate_matrix(vector);

    vector->num_row = order;
    vector->num_col = 1;

    matrix_initializer(vector);

    for (i = 1; i <= order; i++) {
        if (i % 2 == 0) {
            vector->matrix[i-1][0] = 2.0;
        } else {
            vector->matrix[i-1][0] = 0.0;
        }
    }
    store_data(order, new_matrix, vector, len, file_name);
}

bool store_data(int order, MATRIX_p_t some_matrix, MATRIX_p_t vector, int len, char file_name[len]) {
    printf("\n\tPlease enter file name to store matrix: ");
    read_line(file_name, len);

    if ((output_file = fopen(file_name, "w+")) == NULL) {
        printf("\n\t[============= [ERROR: File not found =========]");
        return false;
    }
    fprintf(output_file, "%d\n\n", order); // store the order in file

    // storing the matrix to a file

    int i, j;
    for (i = 0; i < some_matrix->num_row; i++) {
        for (j = 0; j < some_matrix->num_col; j++) {
            fprintf(output_file, "%-20.5f", some_matrix->matrix[i][j]);
        }
        fprintf(output_file, "\n");
    }
    fprintf(output_file, "\n");

    // storing a vector to a file

    for (i = 0; i < vector->num_row; i++) {
        for (j = 0; j < vector->num_col; j++) {
            fprintf(output_file, "%.5f", vector->matrix[i][j]);
        }
        fprintf(output_file, "\n");
    }
    fprintf(output_file, "\n");
    fclose(output_file);

    printf("\n\t\tData stored in file with name (or path) \"%s\"\n\n", file_name);

    return true;
}

EQUATIONS_p_t read_file_data(int len, char file_name[len]) {
    if ((input_file = fopen(file_name, "r")) == NULL) {
        printf("\n\t\t [ERROR: File \"%s\" NOT found]\n\n", file_name);
        return NULL;
    }
    printf("\n\t\t[File: \"%s\" opened]\n", file_name);
    int order = 0;
    fscanf(input_file, "%d", &order);

    EQUATIONS_p_t equation = (EQUATIONS_p_t)malloc(sizeof(EQUATIONS_t));
    equation->matrix = (MATRIX_p_t)malloc(sizeof(MATRIX_t));
    equation->vector = (MATRIX_p_t)malloc(sizeof(MATRIX_t));

    validate_matrix(equation->matrix);
    validate_matrix(equation->vector);

    equation->matrix->num_row = order;
    equation->matrix->num_col = order;

    equation->vector->num_row = order;
    equation->vector->num_col = 1;

    equation->order = order;

    matrix_initializer(equation->matrix);
    matrix_initializer(equation->vector);

    // copy the matrix in file.

    double value;
    int i, j;
    for (i = 0; i < equation->order; i++) {
        for (j = 0; j < equation->order; j++) {
            fscanf(input_file, "%lf", &value);
            equation->matrix->matrix[i][j] = value;
        }
    }

    for (i = 0; i < equation->order; i++) {
        fscanf(input_file, "%lf", &value);
        equation->vector->matrix[i][0] = value;
    }
    if (fclose(input_file) != 0) {
        printf("\n\t\t[--------ERROR: Cannot close file \"%s\"----------]", file_name);
    }
    return equation;
}

void display_case(int len, char file_name[len], int choice) {
    EQUATIONS_p_t equations = open_file(len, file_name);
    if (equations == NULL) return;

    switch (choice) {
        case 2:
            A_inverse_x_b(equations);
            break;
        case 3:
            crammers_rule(equations);
            break;
        case 4:
            gaussian_elimination(equations);
            break;
        case 5:
            A_inverse_x_b(equations);
            break;
        case 7:
            compare_inverse_multiplications(equations, len, file_name);
        default:
            return;
    }
}

void compare_inverse_multiplications(EQUATIONS_p_t equations, int len, char file_name[len]) {

    MATRIX_p_t inverse_A = matrix_inverse(equations->matrix);
    MATRIX_p_t ans_A_x_Ainv = multiply_matrix(equations->matrix, inverse_A);
    MATRIX_p_t ans_Ainv_x_A = multiply_matrix(inverse_A, equations->matrix);

    printf("\n\tPlease enter file name to store matrix: ");
    read_line(file_name, len);

    if ((output_file = fopen(file_name, "w+")) == NULL) {
        printf("\n\t[============= [ERROR: File not found =========]");
        return;
    }
    fprintf(output_file, "%d\n\n", equations->order); // store the order in file

    // storing the matrix to a file

    int i, j;
    for (i = 0; i < ans_A_x_Ainv->num_row; i++) {
        for (j = 0; j < ans_A_x_Ainv->num_col; j++) {
            fprintf(output_file, "%-30.5f", ans_A_x_Ainv->matrix[i][j]);
        }
        fprintf(output_file, "\n");
    }
    fprintf(output_file, "\n");

    // storing a vector to a file

    for (i = 0; i < ans_Ainv_x_A->num_row; i++) {
        for (j = 0; j < ans_Ainv_x_A->num_col; j++) {
            fprintf(output_file, "%-30.5f", ans_Ainv_x_A->matrix[i][j]);
        }
        fprintf(output_file, "\n");
    }
    fprintf(output_file, "\n");
    fclose(output_file);

    printf("\n\t\tFrom computations, [(A inverse) x A] == [A x (A inverse)]\n");
    printf("\n\t\tData stored in file with name (or path) \"%s\"\n\n", file_name);
}

void A_inverse_x_b(EQUATIONS_p_t equations) {
    MATRIX_p_t inverse = matrix_inverse(equations->matrix);
    if (inverse == NULL) {
        printf("\n\t[========= Matrix does not have inverse ============]");
        return;
    }
    MATRIX_p_t answer = multiply_matrix(inverse, equations->vector);
    int i;
    for (i = 0; i < answer->num_row; i++) {
        printf("\n\t\tVariable x%d: %-10.5lf", i+1, answer->matrix[i][0]);
    }
    printf("\n\n");
}

void gaussian_elimination(EQUATIONS_p_t equations) {
    // a = matrix elements
    //
    int n = equations->order;
    double max, pivot, m;
    int maxI;
    double **a, *b, *c, t;

    a = (double **)malloc(sizeof(double) * equations->order);
  //  b = (double *)malloc(sizeof(double) * equations->order);
    c = (double *)malloc(sizeof(double) * equations->order);

    int i, j, k;
    for (i = 0; i < equations->order; i++)
        a[i]= (double *)malloc(sizeof(double) * (equations->order + 1));

    for (i = 0; i < equations->order; i++) {
        for (j = 0; j < equations->order; j++) {
            a[i][j] = equations->matrix->matrix[i][j];
        }
    }

    for (i = 0; i < equations->order; i++) {
        a[i][equations->order] = equations->vector->matrix[i][0];
    }

    // solving now

    for(i=0;i<n;i++){
        for(j=i+1;j<n;j++){
            if(a[i][j]>max){
                max=a[i][j];
                maxI=j;
            }
        }

        for(k=i;k<n+1;k++){
            pivot=a[maxI][k];
            a[maxI][k]=a[i][k];
            a[i][k]=pivot;
        }

        for(k=i+1;k<n;k++){
            m=-a[k][i]/a[i][i];
            for(j=i;j<n+1;j++){
                a[k][j]=a[k][j]+m*a[i][j];
            }
        }
    }

    for(i=n-1;i>=0;i--){
        c[i]=a[i][n]/a[i][i];
        for(k=i-1;k>=0;k--){
            a[k][n]-=a[k][i]*c[i];
        }
}

    for (i = 0; i < equations->order; i++) {
        printf("\n\t\tVariable x%d: %-10.5lf", i+1, c[i]);
    }
    printf("\n\n");

}

void crammers_rule(EQUATIONS_p_t some_matrix) {
    validate_equation(some_matrix);
    double answers[some_matrix->order];

    double matrix_det = determinant(some_matrix->matrix->matrix, some_matrix->order);
    EQUATIONS_p_t matrix_copy = make_copy(some_matrix);
    int i;
    for (i = 0; i < some_matrix->order; i++) {
         vector_substitute(matrix_copy, i);
         double upper_det = determinant(matrix_copy->matrix->matrix, matrix_copy->order);
         answers[i] = (double) (upper_det / matrix_det);
         matrix_copy = make_copy(some_matrix);
    }
    printf("\n\t\tUsing Crammer's rule (-nan means variable value not available)...\n");
    for (i = 0; i < some_matrix->order; i++) {
        printf("\n\t\t\tVariable x%d: %.5lf",  i+1, answers[i]);
    }
    printf("\n\n");
}

EQUATIONS_p_t open_file(int len, char file_name[len]) {
    printf("\tEnter file name to get data: ");
    read_line(file_name, FILE_NAME_LEN);

    EQUATIONS_p_t equations = read_file_data(FILE_NAME_LEN, file_name);
    return equations;
}
/**
 * Substitutes the 'column' in the 'some_eq' matrix with the vector
 * value. Note that this function will modify 'som_eq'
 * @param some_eq Equation to substitute the column with the vector
 * @param column Column to substitute with vector
 * @return
 */
void vector_substitute(EQUATIONS_p_t some_eq, int column) {
    int i;
    for (i = 0; i < some_eq->matrix->num_row; i++) {
        some_eq->matrix->matrix[i][column] = some_eq->vector->matrix[i][0];
    }
}

EQUATIONS_p_t make_copy(EQUATIONS_p_t original_eq) {
    validate_equation(original_eq);
    EQUATIONS_p_t new_equation = (EQUATIONS_p_t)malloc(sizeof(EQUATIONS_t));
    validate_equation(new_equation);

    new_equation->matrix = (MATRIX_p_t)malloc(sizeof(MATRIX_t));
    new_equation->vector = (MATRIX_p_t)malloc(sizeof(MATRIX_t));

    new_equation->order = original_eq->order;
    new_equation->matrix->num_row = original_eq->matrix->num_row;
    new_equation->matrix->num_col = original_eq->matrix->num_col;
    new_equation->vector->num_row = original_eq->vector->num_row;
    new_equation->vector->num_col = original_eq->vector->num_col;

    // allocate right amount of memory for 'new_equation->matrix->matrix'

    matrix_initializer(new_equation->matrix);
    matrix_initializer(new_equation->vector);

    int i, j;

    for (i = 0; i < new_equation->matrix->num_row; i++)
        for (j = 0; j < new_equation->matrix->num_col; j++)
            new_equation->matrix->matrix[i][j] = original_eq->matrix->matrix[i][j];

    for (i = 0; i < new_equation->vector->num_row; i++)
        for (j = 0; j < new_equation->vector->num_col; j++)
            new_equation->vector->matrix[i][j] = original_eq->vector->matrix[i][j];

    return new_equation;
}

void validate_equation(EQUATIONS_p_t some_equation) {
    if (some_equation == NULL) {
         printf("\n\t[============= [ERROR: Insufficient memory =========]");
         exit(EXIT_FAILURE);
    }
}

/**
 * Displays a matrix
 * @param size Size of the matrix
 * @param matrix Matrix to display.
 */
void display_matrix(int size, double *matrix[size]) {
    printf("\t");
    int i, j;
    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++)
            printf("%.2f ", matrix[i][j]);
        printf("\n\t");
    }
}

/**
 * Prints a matrix
 * @param disp_matr Matrix to print.
 */
void print_matrix(MATRIX_p_t disp_matr) {
    printf("\t");
    int i, j;
    for (i = 0; i < disp_matr->num_row; i++) {
        for (j = 0; j < disp_matr->num_col; j++)
            printf("%6.2f", disp_matr->matrix[i][j]);
        printf("\n\t");
    }
}

/**
 * Skips leading white-space characters, then reads the
 * remainder of the input line and stores it in <code>str</code>.
 * Truncates the line if its length exceeds <code>n</code>. Returns
 * the number of characters stored.
 * @param str Character array to store the characters
 * @param n Size of the array
 * @return Number of characters read so far
 */
int read_line(char str[], int n) {
    int ch, i = 0;
    while (isspace(ch = getchar()));
    while (ch != '\n' && ch != EOF) {
        if (i < n)
            str[i++] = ch;
        ch = getchar();
    }
    str[i] = '\0';
    return i;
}

/**
 * Multiplies 2 matrices
 * @param matrix_A First matrix
 * @param matrix_B Second matrix
 * @return Pointer to matrix (ie the answer of the multiplication)
 */
MATRIX_p_t multiply_matrix(MATRIX_p_t matrix_A, MATRIX_p_t matrix_B) {
    if (matrix_A == NULL || matrix_B == NULL)
        terminate_program("[ERROR]: invalid matrix");
    MATRIX_p_t final_matrix = (MATRIX_p_t)malloc(sizeof(MATRIX_t));
    if (final_matrix == NULL)
        terminate_program("[ERROR]: Memory insufficient");

    final_matrix->num_row = matrix_A->num_row;
    final_matrix->num_col = matrix_B->num_col;
    // initialize_matrices(final_matrix->num_row, final_matrix->matrix);

    matrix_initializer(final_matrix);

    int i, j, k;
    for (i = 0; i < matrix_A->num_row; i++) {
        for (j = 0; j < matrix_B->num_col; j++) {
            //final_matrix->matrix[i][j] = 0.0;
            for (k = 0; k < matrix_A->num_col; k++) {
                final_matrix->matrix[i][j] += matrix_A->matrix[i][k] * matrix_B->matrix[k][j];
            }
        }
    }
    return final_matrix;
}

/**
 * Ends program if there was an internal error.
 * @param message Message associated with error.
 */
void terminate_program(char *message) {
    printf("[ERROR]: message -> %s\n", message);        // to be redirected to standard error stream
    exit(EXIT_FAILURE);
}

int display_menu() {
    int choice;
    printf("\n\t########################################################################");
    printf("\n\t####################            MATH PROJECT        ####################");
    printf("\n\t########################################################################");
    printf("\n\n\t[=========================== Options ==================================]\n");
    printf("\n\tHow do you want to solve system of equations OR do extra computations ?\n");
    printf("\n\t1. [====> Insert the hilbert matrix  <====]");
    printf("\n\t2. Compute inverse of matrix and solve x = (A inverse) X b.");
    printf("\n\t3. Solve the problem by Crammer's rule.");
    printf("\n\t4. Solve the problem by Gauss Elimination with partial pivoting.");
    printf("\n\t5. Solve the problem by LU decomposition and return the solution.");
    printf("\n\t6. Find inverse and determinant of Hilbert matrix by LU.");
    printf("\n\t7. Compute [(A inverse) x A] and [A x (A inverse)] and report.");
    printf("\n\t8. Exit program");
    printf("\n\n\t[======================================================================]\n\n");
    printf("\tPlease enter choice from options: ");
    scanf("%d", &choice);
    while (choice < 1 || choice > 8) {
        printf("\tPlease enter VALID choice from optoins: ");
        scanf("%d", &choice);
    }
    return choice;
}

void print_option(int n) {
    printf("\n\t[========================== Option: %d =================================]\n\n", n);
}
