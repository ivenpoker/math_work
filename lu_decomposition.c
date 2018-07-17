#include<stdio.h>
#include<stdlib.h>
#include<math.h>

FILE *output_file;
FILE *results;

void print_menu(){
  printf("+----------------------+\n");
  printf("+   Welcome            +\n");
  printf("+  1) Lu decomposition +\n");
  printf("+  2) quit             +\n");
  printf("+----------------------+\n");
}
int main()
{
  output_file = fopen("lu_matrix.txt", "w");
  print_menu();
  int option;
  scanf("%d", &option);
  while(option != 2 ){
  int i,j,k,n, row;

  double **A, **IA, *b, *x; int *p;
  printf("Enter the order of the matrix: ");
  scanf("%d",&n);
  A = (double**) malloc(sizeof(double)*n);
  for(row = 0; row<n; row++) {
        A[row] = (double *) malloc(n * sizeof(double));
  }
  IA = (double**) malloc(sizeof(double)*n);
  for(row = 0; row<n; row++) {
        IA[row] = (double *) malloc(n * sizeof(double));
  }
  b = malloc(sizeof(double) * n);
  x = malloc(sizeof(double)*n);
  p = malloc(sizeof(int)*n);
  for(i=1; i<=n; i++)
  {
    for(j=1; j<=n; j++)
    {
      double val = (double) 1.0 / (i+j-1);
      A[i-1][j-1] = val ;
      printf("%lf ", A[i-1][j-1]);
    }
    printf("\n");
  }
  for ( i = 0; i < n; ++i ) {
    if (( i + 1 )%2 == 0) {
      b[i] = 2;
    } else {
      b[i] = 0;
    }
  }
  store_data(A, b, n, output_file);
  results = fopen("lu_result.txt", "w");
  int s = LUPDecompose(A, n, 0.000000001, p);
  if ( s != 1 ) {
    printf("Could not be decomposed");
  } else {
    LUPSolve(A, p, b, n, x);
    double det = A[0][0];
    int h;
    for (h = 1; h < n; h++)
        det *= A[h][h];

    if ((p[n] - n) % 2 != 0)
        det = -det;
    printf ("\nHere is the inverse of A:\n");
    fprintf (results, "Here is the inverse of A:\n");
    LUPInvert(A, p, n, IA);
    for(i=0; i<n; i++)
    {
      for (j = 0; j <n; ++j ) {
        printf("%.5e  ",IA[i][j]);
        fprintf(results, "%.5e  ",IA[i][j]);
      }
      printf("\n");
      fprintf(results, "\n");
    }
    printf("\nHere is the vector X(solution vector):\n[x]\n");
    fprintf(results, "\nHere is the vector X(solution vector):\n[x]\n");
    for ( i = 0; i < n; ++i ){
      printf("x[%d] = %.5e  \n",i, x[i]);
      fprintf(results, "x[%d] = %.5e  \n",i, x[i]);
    }
    printf("\nThe determinant of A is:\n\t %.5e\n", det);
    fprintf(results, "\nThe determinant of A is:\n\t %.5e\n", det);
  }
  fclose(results);
  print_menu();
  for(row = 0; row<n; row++) {
        free(A[row]);
        free(IA[row]);
  }
  free(A);
  free(IA);
  A = NULL;
  IA = NULL;
  scanf("%d", &option);
  }
  printf("Thanks and bye!!!!");
  return 0;
}

void LUPSolve(double **A, int *P, double *b, int N, double *x) {
    int i, j, k;
    for ( i = 0; i < N; i++) {
        x[i] = b[P[i]];

        for ( k = 0; k < i; k++)
            x[i] -= A[i][k] * x[k];
    }

    for (i = N - 1; i >= 0; i--) {
        for (k = i + 1; k < N; k++)
            x[i] -= A[i][k] * x[k];

        x[i] = x[i] / A[i][i];
    }
}


int LUPDecompose(double **A, int N, double Tol, int *P) {
    double **L = (double **) malloc(N*sizeof(double *));
    int i, j, k;
    for (i = 0; i < N; ++i ){
      L[i] = (double* ) malloc(N*sizeof(double));
    }
    double **U = (double **) malloc(N*sizeof(double *));
    for ( i = 0; i < N; ++i ){
      U[i] = (double* ) malloc(N*sizeof(double));
    }
    for (i = 0; i < N; i++){
      L[i][i] = 1;
      for ( j =0; j < N; ++j){
        U[i][j] = A[i][j];
      }
    }
    int imax;
    double maxA, *ptr, absA;

    for (i = 0; i <= N; i++)
        P[i] = i;

    for (i = 0; i < N; i++) {
        maxA = 0.0;
        imax = i;

        for (k = i; k < N; k++)
            if ((absA = fabs(A[k][i])) > maxA) {
                maxA = absA;
                imax = k;
            }

        if (maxA == 0) return 0;
        if (imax != i) {
            j = P[i];
            P[i] = P[imax];
            P[imax] = j;
            ptr = A[i];
            A[i] = A[imax];
            A[imax] = ptr;
            P[N]++;
        }

        for (j = i + 1; j < N; j++) {
            A[j][i] /= A[i][i];
            L[j][i] = A[j][i];
            for (k = i + 1; k < N; k++)
                A[j][k] -= A[j][i] * A[i][k];
        }
    }
    for (i = 0; i < N; i++){
      for (j =0; j < N; ++j){
        U[i][j] = A[i][j];
      }
    }
    for (i = 1; i < N; i++){
      for (j =0; j < i; ++j){
        U[i][j] =0;
      }
    }
    printf("\nHere is the matrix L\n[L]\n");
    for(i=0; i<N; i++)
    {
      for (j = 0; j <N; ++j )
        printf("%9.5lf  ",L[i][j]);
      printf("\n");
    }
    printf("\nHere is the matrix U\n[U]\n");
    for(i=0; i<N; i++)
    {
      for (j = 0; j <N; ++j )
        printf("%9.5lf  ", U[i][j]);
      printf("\n");
    }

    return 1;
}

void LUPInvert(double **A, int *P, int N, double **IA) {
    int i, j, k;
    for (j = 0; j < N; j++) {
        for (i = 0; i < N; i++) {
            if (P[i] == j)
                IA[i][j] = 1.0;
            else
                IA[i][j] = 0.0;

            for ( k = 0; k < i; k++)
                IA[i][j] -= A[i][k] * IA[k][j];
        }

        for (i = N - 1; i >= 0; i--) {
            for (k = i + 1; k < N; k++)
                IA[i][j] -= A[i][k] * IA[k][j];

            IA[i][j] = IA[i][j] / A[i][i];
        }
    }
}

void store_data(double **A, double*b, int n, FILE* output_file){
  int i, j;
  if (output_file != NULL){
        fprintf(output_file, "%d\n", n);
  for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            fprintf(output_file, "%.5f", A[i][j]);
        }
        fprintf(output_file, "\n");
    }
    fprintf(output_file, "\n");

    for (i = 0; i < n; i++) {
        fprintf(output_file, "%.5f", b[i]);
        fprintf(output_file, "\n");
    }
    fprintf(output_file, "\n");
    fclose(output_file);
  }
}
