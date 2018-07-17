// File name: main.c
// author: CSC403 Group II

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define INPUT_FILE "info.txt"    // define constant to hold the input file name
#define OUTPUT_FILE "output.txt"  // defines constant to hold the output file name

// prototypes
void regenerate_matrices(int n, double m[n][n], int b[n]);
int factorial(int);
int comb(int, int);
void inverse(int n, double m[n][n]);
void inverse_solution_of_matrix(int n, double m[n][n], int b[n]);
void matrix_mul(int n, double m1[n][n], double m2[n][n], double result[n][n]);
void matrix_mul_findings(int n, double m1[n][n], double m2[n][n]);
double determinant(int n, double m[n][n]);
void crammer(int n, double m[n][n], int b[n]);
void guassian_elimination(int n, double m[n][n], int b[n]);
void lu_decomposition(int n, double A[n][n], int b[n]);
int display_menu();

int main(void)
{
    FILE * fp;  // file pointer

    // open file and ensures that it actually opened
    if ((fp = fopen(INPUT_FILE, "r")) == NULL){
        puts("File could not be opened");
        exit(EXIT_FAILURE);
    }
    int n, i, j;

    fscanf(fp, "%d", &n);  // reads n from the file
    double A[n][n];
    int b[n];

    // read the matrix from the file
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            fscanf(fp, "%lf", &A[i][j]);
        }
    }

    // reads the vector from the file
    for(i = 0; i < n; i++){
        fscanf(fp, "%d", &b[i]);
    }

    fclose(fp); // fclose closes the file

  int option;
  char answer;
  /*
  This do while loop handles the menu display and calls the other methods to solve individual problems
  based on the options specification.
  It displays a menu and takes the choice from the user. It then calls the function responsible for
  solving that problem.
  */
  do{
    option = display_menu();

    switch(option){
    case 1: regenerate_matrices(n, A, b);
            inverse_solution_of_matrix(n, A, b);
            break;
    case 2: regenerate_matrices(n, A, b);
            crammer(n, A, b);
            break;
    case 3: regenerate_matrices(n, A, b);
            guassian_elimination(n, A, b);
            break;
    case 4: regenerate_matrices(n, A, b);
            lu_decomposition(n, A, b);
            break;
    case 5: exit(EXIT_SUCCESS);
            break;
    default: printf("Please enter a valuable option.");
            break;

    }

    printf("\nDo you want to continue (y or n): ");
    scanf(" %c", &answer);
  } while(answer == 'y' || answer == 'Y');

  return 0;
}


/*
    This function regenerates the matrix getting it back from the file.
    This is used to set the matrix and vector back to its defualt format.
*/
void regenerate_matrices(int n, double m[n][n], int b[n]){
    FILE * fp;  // file pointer

    // open file and ensures that it actually opened
    if ((fp = fopen(INPUT_FILE, "r")) == NULL){
        puts("File could not be opened");
        exit(EXIT_FAILURE);
    }
    int i, j;

    fscanf(fp, "%d", &n);  // reads n from the file

    // read the matrix from the file
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            fscanf(fp, "%lf", &m[i][j]);
        }
    }

    // reads the vector from the file
    for(i = 0; i < n; i++){
        fscanf(fp, "%d", &b[i]);
    }

    fclose(fp); // fclose closes the file

}


// This function computes the factorial of any number n
int factorial(int n)
{
    if (n < 0)
        printf("Error! Factorial of a negative number doesn't exist.");
    else
    {
        int fact = 1;
        int i;
        for(i=1; i<=n; ++i)
        {
            fact *= i;              // fact = fact*i;
        }
        return fact;
    }

}

// This function computes n combination r
int comb(int n,int r){
  return factorial(n)/(factorial(n-r)*factorial(r));
}

/* This function multiplies two matrices and store their result in a third matrix*/
void matrix_mul(int n, double m1[n][n], double m2[n][n], double result[n][n]){

    size_t i, j, k;

    // Initialise all entries of result to 0
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            result[i][j] = 0;
        }
    }

    // Calculate the multiplication of two matrices
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            for(k = 0; k < n; k++){
                result[i][j] += m1[i][k] * m2[k][j];
            }
        }
    }

}

/* The matrix takes two matrices[an inverse of a matrix and the matrix itself] and multiply two combinations
    inv(A) * A      and     A * inv(A).
    It then stores thier individual results in a file.
*/
void matrix_mul_findings(int n, double m1[n][n], double m2[n][n]){

    double result1[n][n], result2[n][n];  // two matrices to store respective results
    int i, j;

    FILE *fp;  // file pointer

    // open file and ensures that it actually opened
    if((fp = fopen(OUTPUT_FILE, "w")) == NULL){
        puts("File could not be opened");
        exit(EXIT_FAILURE);
    }

    // This computes the m * inv(m)
    matrix_mul(n, m1, m2, result1);
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            fprintf(fp, "%f ", result1[i][j]);
        }
        fprintf(fp, "\n");
    }

    fprintf(fp, "\n");

    // This computes the inv(m) * m
    matrix_mul(n, m2, m1, result2);
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            fprintf(fp, "%f ", result2[i][j]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);  // this closes the file

}


// This code fragment computes the inverse of a matrix by using combinations
void inverse(int n, double m[n][n]){
  size_t i, j;
  for(i=1;i<=n;i++){
    for(j=1;j<=n;j++){
      m[i-1][j-1]=pow(-1, i+j) * (i+j-1) * comb(n+i-1, n-j) * comb(n+j-1, n-i) * pow(comb(i+j-2, i-1), 2);
    }
  }
}

/*
    The function solves the equation x = inv(A) * b.
    It takes a matrix and its vector and uses the inverse function to get the solution.
*/
void inverse_solution_of_matrix(int n, double m[n][n], int b[n]){

    // This copies matrix inorder to have an unaltered matrix copy
    // be latered used by other functions
    double matrix[n][n];
    int i, j;
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            matrix[i][j] = m[i][j];
        }
    }

    inverse(n, m);   // computes the inverse of the matrix

    printf("\n\nThe inverse of the matrix is as follows: \n");
    /*display the inverse matrix*/
    for (i = 0; i < n; i++){
        for(j =0; j < n; j++){
            printf("%.5f ", m[i][j]);
        }
        printf("\n");
    }

    //solving the matrix i.e. x = inv(A)*b
    double x[n];
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            x[i] += m[i][j] * b[j];
        }
    }

    printf("\nThe solution of the matrix is as follows:\n");
    for(i = 0; i < n; i++){
        printf("x%d = %f\n", i+1, x[i]);
    }

    matrix_mul_findings(n, m, matrix); // solving question 2(b)

}


/*
    This function computes the determinant of a matrix
*/
double determinant(int n, double m[n][n]){
    double s = 1, det = 0, tmp_arr[n][n];
    int p, q, c, j, i;
    for(p = 0; p < n; p++){
        for(q = 0; q < n; q++){
            tmp_arr[p][q] = 0;
        }

    }
    if (n == 1) {
        return (m[0][0]);
    } else {
        det = 0;
        for ( c = 0; c < n; c++) {
            p = 0; q = 0;
            for (i = 0; i < n; i++) {
                for ( j = 0; j < n; j++) {
                    tmp_arr[i][j] = 0;
                    if (i != 0 && j != c) {
                        tmp_arr[p][q] = m[i][j];
                        if (q < (n-2)) q++;
                        else {
                            q = 0;
                            p++;
                        }
                    }
                } // end of inner for loop (with j)
            } // end of inner for loop (with i)
             det = det + s * (m[0][c] * determinant(n-1, tmp_arr));
             s = -1 * s;
        } // end of inner for loop (with c)
    }
    return (det);
}


// This code fragment solve the matrix by crammer's rule
void crammer(int n, double m[n][n], int b[n]){
    double det = determinant(n, m);
	printf("\nThe determinant of the given square matrix is \ndet=%f\n\n",det);
	int i, k, j, p;
	double result[n];
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
            m[j][i] = b[j];
        }
        result[i]= determinant(n, m) / det;
        regenerate_matrices(n, m, b);
    }
    for(i = 0; i < n; i++){
        printf("x%d = %f\n", i+1, result[i]);

    }

}


/*This method solves a linear system by guassian elimination*/
void guassian_elimination(int n, double m[n][n], int b[n]){

  int i, j, k;
  double c, x[n], sum, matrix[n][n+1];

  for(i = 0; i < n; i++){
    for(j = 0; j < n; j++){
        matrix[i][j] = m[i][j];
    }
    matrix[i][n] = b[n-1];
  }

  // loop for generating the upper triangular matrix
  for(j = 1; j <= n; j++){
    for(i = 1; i <= n; i++){
      if(i > j){
        c = matrix[i][j]/matrix[j][j];
        for(k = 1; k <= n+1; k++){
          matrix[i][k] = matrix[i][k] - c * matrix[j][k];
        }
      }
    }
  }

  x[n] = matrix[n][n+1] / matrix[n][n];

  /* This loop handles backward substitution*/
  for(i = n - 1; i >= 1; i--){
    sum = 0;
    for(j = i + 1; j <= n; j++){
      sum = sum + matrix[i][j] * x[j];
    }
    x[i] = (m[i][n+1] - sum) / matrix[i][i];
  }

  printf("\nThe solution is: \n");
  for(i = 1; i <= n; i++){
    printf("\nx%d = %f\t", i, x[i]); /* x1, x2, x3 are the required solutions */
  }

}

/*
    This function solves the linear system by LU decomposition
    It generates the uper and lower triangular matrices and used these matrices to
    get the solution of the linear system.
*/
void lu_decomposition(int n, double A[n][n], int b[n]){
     int i,j,k;
     int p=0;

     float L[n][n];//lower triangular matrix
     float U[n][n];//upper triangular matrix
     float sum=0.0;//holder
     float x[n];//
     float z[n];

     for(i = 0; i < n; i++){
        x[i] = 0;
        z[i] = 0;
     }

	//LU decomposition
	//upper triangular matrix
	for(i=1;i<=n;i++)
	{
	//setting main diagonal to 1
		L[i][i]=1;
		for(k=i;k<=n;k++)
		{
           sum=0;
           for (p=1;p<=i-1;p++)
           	 sum+=L[k][p]*U[p][k];
			L[i][j]=(A[i][j]-sum)/L[k][k];
		}

		for(j=i+1;j<=n;j++)
        {
            sum=0;
            for(p=1;p<=i-1;p++)
                sum+=L[i][p]*U[p][j];
            U[i][j]=(A[i][j]-sum)/L[i][i];
        }
	}

	//dispaly the LU matrix
	//input code later

	//LZ=b. find Z

	for (i=1;i<=n;i++)
	{
		sum+=L[i][p]*z[p];
	    z[i]=(b[i]-sum)/L[i][i];
    }

    //UX=z
    for (i=n;i>0;i--)
    {
    	sum=0;
    	for (p=n;p>i;p--)
    		sum+=U[i][p]*x[p];
    	x[i]=(z[i]-sum)/U[i][i];
    }

    //display solution x[i]

	for(i=0;i<n;i++)
	{

		printf("%f, ",x[i]);
	}
}


/*This function displays the menu containing the options to perform on the hilbert matrix
 It takes a choice from the user and returns it.
*/
int display_menu(){
    int choice;
    printf("------------------------------------------------------------------------------------------\n");
    printf("\nChoose one of the options below to solve the Hilbert matrix\n");
    printf("1 -- Solve the problem by calculating the inverse of A and then computing x = inv(A)*b\n");
    printf("2 --- Solve the problem by Cramers rule.\n");
    printf("3 --- Solve the problem by Guass elimination with partial pivoting.\n");
    printf("4 --- Solve the problem by LU decomposition.\n");
    printf("5 --- Exit");
    printf("\nEnter your option: ");
    scanf("%d", &choice);
    while(choice < 1 || choice > 5){
        printf("Please enter a VALID choice: ");
        scanf("%d", &choice);
    }
    return choice;
}
