//File name: hilbert.c
//Author: CSC403 Group II


/*              Purpose
                -------
  This programme generates the hilbert matrix and stores it in a file
  It also generates a vector and stores it in that same file
  It stores along n which belongs to the size of the matrix i.e. nxn matrix at the begining
*/

/*
 *  Hilbert matrix formula: A(I, J) = 1 / ( I + J - 1)
 *  Vector constrain: bi = -2 if i is even else 0 if i is odd for i = 0, 1, ... n
*/

#include<stdio.h>
#include<stdlib.h>

#define FILE_NAME "info.txt"          // This contains the name of our file

// prototypes
void hilbert(int n, float m[n][n]);
void vector(int n, int b[n]);
void save_matrices(int n, float m[n][n], int b[n]);

/* The main function where all the saving into our file processes take place
   It calls other functions to help with processing.
*/
int main(void)
{
   int size;

   printf("Enter the size of the hilbert matrix: ");
   scanf("%d", &size);

   float A[size][size];  // matrix of size nxn
   int b[size];
   size_t i, j;  // counter

   // Initialize elements of the array to zero
   for (i = 0; i < size; i++){
       for (j = 0; j < size; j++){
            A[i][j] = 0;
       }
   }

   hilbert(size, A);   // generates the hilbert matrix

   vector(size, b);    // generates the vector b

   save_matrices(size, A, b);


}

/*  This function takes a matrix and convert it to a hilbert matrix
    It uses the formula A(I, J) = 1 / (I + J -1)
*/
void hilbert(int n, float m[n][n])
{
   size_t i, j;
   for (i = 1; i <= n; i++){
      for (j = 1; j <= n; j++){
	    // Addition of a 1 instead of a subtracting, since array begins at (0,0) not (1,1) like a matrix
            m[i-1][j-1] = 1.0 / (float) (i + j - 1);
      }
   }

}

/*
    Function generates the vector b, by assigning -2 if position, i is even
    and assigning 0 if position, i is odd
*/
void vector(int n, int b[n])
{
   size_t i;
   for (i = 1; i <= n; i++){
      b[i-1] = i % 2 == 0 ? -2 : 0;
   }

}

/*This method create a file and saves the hilbert matrix, and its vector in that file along with n*/
void save_matrices(int n, float m[n][n], int b[n]){

  FILE *fp;        // file pointer

   // Try to open/create file, if file fails print message and exit
   // else write the size, matrix and vector to file
   if ( (fp = fopen(FILE_NAME, "w")) == NULL){
        puts("File could not be opened");
        exit(EXIT_FAILURE);
   }
   else {

        // Write the size of the matrix to the file
        fprintf(fp, "%d\n", n);

        size_t i, j;
        // Write each row of the matrix on a separate line
        for(i = 0; i < n; i++){
            for(j = 0; j < n; j++){
                fprintf(fp, "%.5f ", m[i][j]);

            }
            fprintf(fp, "\n");
        }

        // Write each row of the vector on a separate line
        for(i = 0; i < n; i++){
            fprintf(fp, "%d\n", b[i]);
        }
   }

   fclose(fp);    // fclose closes file
}
