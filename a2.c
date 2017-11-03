/**
  a2.c
  MPI implementation of the parallel merging algorithm.
  Contributors:
  -> Frank Khalil	  160226600
  -> Sarah Johnston	150139570
  -> Brad Katz		  130750210
  -> Gareth Sharpe	090361370
*/

// Includes
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdbool.h> // For bools in C! Works with C99+
#include <math.h>

// CONSTANTS
const int FIRST = 0;

const int RAND_STEP_SIZE  = 10;
const int INDEX_NOT_FOUND = -1;

// Static functions
/**
  Generates two array_size arrays of ints/
  
  @param arr_a Pointer to the first array
  @param arr_b Pointer to the second array
  @param The length of the arrays as given on the command line
*/
static void gen_arrays(int arr_a[], int arr_b[], const int array_size)
{  
    // Seed RNG
    srand(time(NULL));
  
    // Generate the first element of each array
    arr_a[0] = rand() % RAND_STEP_SIZE;
    arr_b[0] = rand() % RAND_STEP_SIZE;
  
    for (int i = 1; i < array_size; ++i)
    {
        arr_a[i] = (rand() % ((i * RAND_STEP_SIZE) - arr_a[i - 1])) + arr_a[i - 1];
        arr_b[i] = (rand() % ((i * RAND_STEP_SIZE) - arr_b[i - 1])) + arr_b[i - 1];
    }
}

/**
  Prints a given array out to the console.
  Ex: [1, 2, 3, 4] is printed as 1, 2, 3, 4 followed by a newline.
  
  @param array The array to be printed
  @param array_size The length of the array
*/
static void print_array(const int array[], const int array_size)
{
    for (int i = 0; i < array_size - 1; ++i)
        printf("%d, ", array[i]);
  
    printf("%d\n", array[array_size - 1]);
}

/**
  Runs a binary search on a given array and returns the index of the last occurrence of num_to_find if it exists in the array.
  
  @param array Pointer to the array to search
  @param num_to_find The integer to search for in the array
  @param array_size The length of the array to search
  
  @return index_found The index of the last occurrence of num_to_find, or -1 if num_to_find does not exist in the array
*/
static int binary_search(const int array[], const int num_to_find, const int array_size)
{
    int first_index = 0;
    int last_index = array_size - 1;
    int midpoint = floor((first_index + last_index) / 2);
  
    // Assume error case
    int index_found = -1;
  
    while (first_index <= last_index)
    {
        if (array[midpoint] < num_to_find)
        {
            first_index = midpoint + 1;
        }
        else if (array[midpoint] == num_to_find)
        {
            index_found = midpoint;
            first_index = midpoint + 1;
        }
        else
        {
            last_index = midpoint - 1;
        }
      
        midpoint = floor((first_index + last_index) / 2);
    }
  
    return index_found;
}

int main(int argc, char *argv[])
{
  int process_rank;
  int num_processors;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_processors);
  
  if (argc < 2)
  {
      if (process_rank == FIRST)
        printf("ERROR: Missing array length exponent. Usage: ./a2 [array length exponent]\n");
      
      MPI_Finalize();
      return 0;
  }
  
  // The array size is 2^arr_pow
  const int arr_pow = atoi(argv[1]);
  const int array_size = pow(2, arr_pow);
  
  int *arr_a = calloc(array_size, sizeof(int));
  int *arr_b = calloc(array_size, sizeof(int));
  
  MPI_Status status;
  MPI_Request request;
  
  // Only process 0 needs to generate the two arrays
  if (process_rank == FIRST)
      gen_arrays(arr_a, arr_b, array_size);
  
  const int arr_a_group_size = arr_pow;
  
  int *sub_arr_a = calloc(arr_a_group_size, sizeof(int));
  
  // Scatter a partion of arr_a to all other processors
  MPI_Scatter(arr_a, arr_a_group_size, MPI_INT, sub_arr_a, arr_a_group_size, MPI_INT, FIRST, MPI_COMM_WORLD);
  
  if (process_rank == FIRST)
  {
      int *sendcounts = calloc(num_processors, sizeof(int)); // How many elements to send to each of the n-1 processes
      int *displacements = calloc(num_processors, sizeof(int)); // The displacements where each segment begins
    
      // Partition arr_b and send the partitions to the n-1 processors.
      for (int current_processor = 1; current_processor < num_processors - 1; ++current_processor)
      {
          int last_index_of_sub_arr_b = binary_search(sub_arr_a, sub_arr_a[arr_a_group_size - 1], arr_a_group_size);
      }
      
          
      int *arr_c = calloc(array_size * 2, sizeof(int));
    
      // Receive merged partitions and shove them in arr_c
      
    
      print_array(arr_c, array_size * 2);
      free(arr_c);
  }
  else if (process_rank > FIRST)
  {
      // Receive sub_arr_b start and end indices from FIRST and fill sub_arr_b
      int sub_arr_b_start_index = 0;
      int sub_arr_b_end_index = 0;
    
      MPI_Recv(&sub_arr_b_start_index, 1, MPI_INT, FIRST, TAG_SUB_ARR_B_START_INDEX, MPI_COMM_WORLD, &status);
      MPI_Recv(&sub_arr_b_end_index, 1, MPI_INT, FIRST, TAG_SUB_ARR_B_END_INDEX, MPI_COMM_WORLD, &status);
      
      const int sub_arr_b_size = (sub_arr_b_end_index - sub_arr_b_start_index) + 1;
    
      int *sub_arr_b = calloc(sub_arr_b_size, sizeof(int));
      int sub_arr_index = 0;
    
      for (int i = sub_arr_b_start_index; i < sub_arr_b_end_index + 1; ++i)
      {
          sub_arr_b[sub_arr_index++] = arr_b[i]; 
      }
    
      // Merge local partitions of arr_a and arr_b
      const int sub_arr_c_size = sub_arr_b_size + arr_a_group_size;
      
      int *sub_arr_c = calloc(sub_arr_c_size, sizeof(int));
      int sub_arr_c_index = 0, i = 0, j = 0;
      
      while (i < arr_a_group_size && j < sub_arr_b_size)
      {
          if (sub_arr_a[i] < sub_arr_b[j])
              sub_arr_c[sub_arr_c_index++] = sub_arr_a[i++];
          else
              sub_arr_c[sub_arr_c_index++] = sub_arr_b[j++];
      }
    
      if (i >= arr_a_group_size)
          while (j < sub_arr_b_size)
              sub_arr_c[sub_arr_c_index++] = sub_arr_b[j++];
    
      if (j >= sub_arr_b_size)
          while (i < arr_a_group_size)
              sub_arr_c[sub_arr_c_index++] = sub_arr_a[i++];
    
      // Send merged partition back to the root process
      
    
    
      free(sub_arr_b);
      free(sub_arr_C);
  }
  
  MPI_Finalize();
  free(arr_a);
  free(arr_b);
  free(sub_arr_a);
  return 0;
}
