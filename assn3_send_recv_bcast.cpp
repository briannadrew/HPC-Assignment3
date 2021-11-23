// Name: assn3.cpp
// Author: Brianna Drew
// ID: #0622446
// Date Created: 2021-11-17
// Last Modified: 2021-11-19

// libraries
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <mpi.h>

int main (int argc, char** argv) {
  int procNum, procRank; // number of processes and process number
  const int arr_size = 100; // size of array (number of columns and number of rows)
  const int block_size = (arr_size / 2) + 2; // size of each sub-array with halos on all sides (4 sub-arrays for 4 threads)
  const float init_val = 10000000000000000000000000000000;
  MPI_Status status;

  // initialize MPI and get # of processes and process number
  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &procNum);
  MPI_Comm_rank(MPI_COMM_WORLD, &procRank);

  // print error and exit if number of processes is not 4
  if (procNum != 4) {
    fprintf(stderr,"ERROR: Please use 4 processes for this application.");
    MPI_Finalize();
    exit(-1);
  }

  /*
    GET NEIGHBOUR PROCESS ID's  
  */
  int leftProc, rightProc, topProc, bottomProc, topL_Proc, topR_Proc, bttmL_Proc, bttmR_Proc;
  rightProc = procRank + 1;
  if (rightProc >= procNum || procRank == 1) {rightProc = MPI_PROC_NULL;}
  leftProc = procRank - 1;
  if (leftProc < 0 || procRank == 2) {leftProc = MPI_PROC_NULL;}
  topProc = procRank - 2;
  if (topProc < 0) {topProc = MPI_PROC_NULL;}
  bottomProc = procRank + 2;
  if (bottomProc >= procNum) {bottomProc = MPI_PROC_NULL;}
  topL_Proc = procRank - 3;
  if (procRank != 3) {topL_Proc = MPI_PROC_NULL;}
  topR_Proc = procRank - 1;
  if (procRank != 2) {topR_Proc = MPI_PROC_NULL;}
  bttmL_Proc = procRank + 1;
  if (procRank != 1) {bttmL_Proc = MPI_PROC_NULL;}
  bttmR_Proc = procRank + 3;
  if (procRank != 0) {bttmR_Proc = MPI_PROC_NULL;}

  // create 2D array (with room for halos all around)
  float *arr = new float[(block_size) * (block_size)];

  // initialize all array values to 0
  for (int i = 0; i < block_size; i++) {
    for (int j = 0; j < block_size; j++) {
      arr[i*block_size+j] =  0;
    }
  }

  float max = 0; // largest number in array
  int count = 0, i, j = 1, max_i, max_j = 1; // count for each iteration of the simulation

  // update middle leftmost value (of total array) to largest number
  if (procRank == 0) {
    arr[(block_size - 2)*block_size+1] = 10000000000000000000000000000000;
    i = block_size - 2;
    max_i = block_size - 2;
  }
  else {
    i = 1;
    max_i = 1;
  }

  // create buffers to store all values of each process array in each process (easy access to neighbours)
  float *r1 = new float[(block_size) * (block_size)];
  float *r2 = new float[(block_size) * (block_size)];
  float *r3 = new float[(block_size) * (block_size)];
  float *r4 = new float[(block_size) * (block_size)];

  // initialize all neighbour array values to 0
  for (int i = 0; i < block_size; i++) {
    for (int j = 0; j < block_size; j++) {
      r1[i*block_size+j] =  0;
      r2[i*block_size+j] =  0;
      r3[i*block_size+j] =  0;
      r4[i*block_size+j] =  0;
    }
  }
  r1[(block_size - 2)*block_size+1] = init_val;

  // repeat simulation until no cell has a number larger than 12 starting at the cell with the largest value
  do {
    do {
      int neighbours = 0;
      float neighbour_sum = 0, neighbour_avg = 0, new_vals = 0;

      // top left corner
      if (i == 1 && j == 1) {
        // determine # of neighbours and sum their values
        neighbours = 3;
        neighbour_sum += arr[i*block_size+(j+1)];
        neighbour_sum += arr[(i+1)*block_size+j];
        neighbour_sum += arr[(i+1)*block_size+(j+1)];

        if (leftProc != MPI_PROC_NULL) {
          neighbours += 2;
          // SEND, RECEIVES, AND UPDATES HERE
          if (procRank == leftProc) {
            //MPI_Send(&arr[i*block_size+(block_size - 2)], 1, MPI_FLOAT, rightProc, 1, MPI_COMM_WORLD);
            //MPI_Send(&arr[(i + 1)*block_size+(block_size - 2)], 1, MPI_FLOAT, rightProc, 2, MPI_COMM_WORLD);
            switch (leftProc)
            {
            case 0:
              neighbour_sum += r1[i*block_size+0];
              neighbour_sum += r1[(i + 1)*block_size+0];
              break;
            case 1:
              neighbour_sum += r2[i*block_size+0];
              neighbour_sum += r2[(i + 1)*block_size+0];
              break;
            case 2:
              neighbour_sum += r3[i*block_size+0];
              neighbour_sum += r3[(i + 1)*block_size+0];
              break;
            default:
              break;
            }
          }
          //MPI_Recv(&arr[i*block_size+0], 1, MPI_FLOAT, leftProc, 1, MPI_COMM_WORLD, &status);
          //MPI_Recv(&arr[(i + 1)*block_size+0], 1, MPI_FLOAT, leftProc, 2, MPI_COMM_WORLD, &status);
          //neighbour_sum += arr[i*block_size+0];
          //neighbour_sum += arr[(i + 1)*block_size+0];
        }
        if (topProc != MPI_PROC_NULL) {
          neighbours += 2;
          // SEND, RECEIVES, AND UPDATES HERE
          if (procRank == topProc) {
            //MPI_Send(&arr[(block_size - 2)*block_size+j], 1, MPI_FLOAT, bottomProc, 3, MPI_COMM_WORLD);
            //MPI_Send(&arr[(block_size - 2)*block_size+(j + 1)], 1, MPI_FLOAT, bottomProc, 4, MPI_COMM_WORLD);
            switch (topProc)
            {
            case 0:
              neighbour_sum += r1[(i - 1)*block_size+j];
              neighbour_sum += r1[(i - 1)*block_size+(j + 1)];
              break;
            case 1:
              neighbour_sum += r2[(i - 1)*block_size+j];
              neighbour_sum += r2[(i - 1)*block_size+(j + 1)];
              break;
            default:
              break;
            }
          }
          //MPI_Recv(&arr[(i - 1)*block_size+j], 1, MPI_FLOAT, topProc, 3, MPI_COMM_WORLD, &status);
          //MPI_Recv(&arr[(i - 1)*block_size+(j + 1)], 1, MPI_FLOAT, topProc, 4, MPI_COMM_WORLD, &status);
          //neighbour_sum += arr[(i - 1)*block_size+j];
          //neighbour_sum += arr[(i - 1)*block_size+(j + 1)];
        }
        if (topL_Proc != MPI_PROC_NULL) {
          neighbours += 1;
          // SEND, RECEIVES, AND UPDATES HERE
          if (procRank == topL_Proc) {
            //MPI_Send(&arr[(block_size - 2)*block_size+(block_size - 2)], 1, MPI_FLOAT, bttmR_Proc, 5, MPI_COMM_WORLD);
            neighbour_sum += r1[(i - 1)*block_size+0];
          }
          //MPI_Recv(&arr[(i - 1)*block_size+0], 1, MPI_FLOAT, topL_Proc, 5, MPI_COMM_WORLD, &status);
          //neighbour_sum += arr[(i - 1)*block_size+0];
        }

        // calculate average of neighbours
        neighbour_avg = neighbour_sum / neighbours;
        // calculate new values for neighbours (5% of the difference between current cell and average of neighbours)
        new_vals = abs(0.05 * (arr[i*block_size+j] - neighbour_avg));

        // change values of neighbours to new value
        arr[i*block_size+(j+1)] = new_vals;
        arr[(i+1)*block_size+j] = new_vals;
        arr[(i+1)*block_size+(j+1)] = new_vals;

        if (leftProc != MPI_PROC_NULL) {
          // SEND, RECEIVES, AND UPDATES HERE
          MPI_Send(&new_vals, 1, MPI_FLOAT, leftProc, 1, MPI_COMM_WORLD);
          MPI_Send(&new_vals, 1, MPI_FLOAT, leftProc, 2, MPI_COMM_WORLD);
          if (procRank == leftProc) {
            MPI_Recv(&arr[i*block_size+(block_size - 2)], 1, MPI_FLOAT, rightProc, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(&arr[(i + 1)*block_size+(block_size - 2)], 1, MPI_FLOAT, rightProc, 2, MPI_COMM_WORLD, &status);
          }
        }
        if (topProc != MPI_PROC_NULL) {
          // SEND, RECEIVES, AND UPDATES HERE
          MPI_Send(&new_vals, 1, MPI_FLOAT, topProc, 3, MPI_COMM_WORLD);
          MPI_Send(&new_vals, 1, MPI_FLOAT, topProc, 4, MPI_COMM_WORLD);
          if (procRank == topProc) {
            MPI_Recv(&arr[(block_size - 2)*block_size+j], 1, MPI_FLOAT, bottomProc, 3, MPI_COMM_WORLD, &status);
            MPI_Recv(&arr[(block_size - 2)*block_size+(j + 1)], 1, MPI_FLOAT, bottomProc, 4, MPI_COMM_WORLD, &status);
          }
        }
        if (topL_Proc != MPI_PROC_NULL) {
          // SEND, RECEIVES, AND UPDATES HERE
          MPI_Send(&new_vals, 1, MPI_FLOAT, topL_Proc, 5, MPI_COMM_WORLD);
          if (procRank == topL_Proc) {
            MPI_Recv(&arr[(block_size - 2)*block_size+(block_size - 2)], 1, MPI_FLOAT, bttmR_Proc, 5, MPI_COMM_WORLD, &status);
          }
        }
      }

      // bottom left corner
      else if (i == block_size - 2 && j == 1) {
        // determine # of neighbours and sum their values
        neighbours = 3;
        neighbour_sum += arr[i*block_size+(j+1)];
        neighbour_sum += arr[(i-1)*block_size+(j+1)];
        neighbour_sum += arr[(i-1)*block_size+j];

        if (leftProc != MPI_PROC_NULL) {
          neighbours += 2;
          // SEND, RECEIVES, AND UPDATES HERE
          if (procRank == leftProc) {
            //MPI_Send(&arr[i*block_size+(block_size - 2)], 1, MPI_FLOAT, rightProc, 1, MPI_COMM_WORLD);
            //MPI_Send(&arr[(i - 1)*block_size+(block_size - 2)], 1, MPI_FLOAT, rightProc, 2, MPI_COMM_WORLD);
            switch (leftProc)
            {
            case 0:
              neighbour_sum += r1[i*block_size+0];
              neighbour_sum += r1[(i - 1)*block_size+0];
              break;
            case 1:
              neighbour_sum += r2[i*block_size+0];
              neighbour_sum += r2[(i - 1)*block_size+0];
              break;
            case 2:
              neighbour_sum += r3[i*block_size+0];
              neighbour_sum += r3[(i - 1)*block_size+0];
              break;
            default:
              break;
            }
          }
          //MPI_Recv(&arr[i*block_size+0], 1, MPI_FLOAT, leftProc, 1, MPI_COMM_WORLD, &status);
          //MPI_Recv(&arr[(i - 1)*block_size+0], 1, MPI_FLOAT, leftProc, 2, MPI_COMM_WORLD, &status);
          //neighbour_sum += arr[i*block_size+0];
          //neighbour_sum += arr[(i - 1)*block_size+0];
        }
        if (bottomProc != MPI_PROC_NULL) {
          neighbours += 2;
          // SEND, RECEIVES, AND UPDATES HERE
          if (procRank == bottomProc) {
            //MPI_Send(&arr[1*block_size+1], 1, MPI_FLOAT, topProc, 3, MPI_COMM_WORLD);
            //MPI_Send(&arr[1*block_size+(j + 1)], 1, MPI_FLOAT, topProc, 4, MPI_COMM_WORLD);
            switch (bottomProc)
            {
            case 2:
              neighbour_sum += r3[(block_size - 1)*block_size+j];
              neighbour_sum += r3[(block_size - 1)*block_size+(j + 1)];
              break;
            case 3:
              neighbour_sum += r4[(block_size - 1)*block_size+j];
              neighbour_sum += r4[(block_size - 1)*block_size+(j + 1)];
              break;
            default:
              break;
            }
          }
          //MPI_Recv(&arr[(block_size - 1)*block_size+j], 1, MPI_FLOAT, bottomProc, 3, MPI_COMM_WORLD, &status);
          //MPI_Recv(&arr[(block_size - 1)*block_size+(j + 1)], 1, MPI_FLOAT, bottomProc, 4, MPI_COMM_WORLD, &status);
          //neighbour_sum += arr[(block_size - 1)*block_size+j];
          //neighbour_sum += arr[(block_size - 1)*block_size+(j + 1)];
        }
        if (bttmL_Proc != MPI_PROC_NULL) {
          neighbours += 1;
          // SEND, RECEIVES, AND UPDATES HERE
          if (procRank == bttmL_Proc) {
            //MPI_Send(&arr[1*block_size+(block_size - 2)], 1, MPI_FLOAT, topR_Proc, 5, MPI_COMM_WORLD);
            neighbour_sum += r3[(block_size - 1)*block_size+0];
          }
          //MPI_Recv(&arr[(block_size - 1)*block_size+0], 1, MPI_FLOAT, bttmL_Proc, 5, MPI_COMM_WORLD, &status);
          //neighbour_sum += arr[(block_size - 1)*block_size+0];
        }

        // calculate average of neighbours
        neighbour_avg = neighbour_sum / neighbours;
        // calculate new values for neighbours (5% of the difference between current cell and average of neighbours)
        new_vals = abs(0.05 * (arr[i*block_size+j] - neighbour_avg));

        // change values of neighbours to new value
        arr[i*block_size+(j+1)] = new_vals;
        arr[(i-1)*block_size+(j+1)] = new_vals;
        arr[(i-1)*block_size+j] = new_vals;

        if (leftProc != MPI_PROC_NULL) {
          // SEND, RECEIVES, AND UPDATES HERE
          MPI_Send(&new_vals, 1, MPI_FLOAT, leftProc, 1, MPI_COMM_WORLD);
          MPI_Send(&new_vals, 1, MPI_FLOAT, leftProc, 2, MPI_COMM_WORLD);
          if (procRank == leftProc) {
            MPI_Recv(&arr[i*block_size+(block_size - 2)], 1, MPI_FLOAT, rightProc, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(&arr[(i - 1)*block_size+(block_size - 2)], 1, MPI_FLOAT, rightProc, 2, MPI_COMM_WORLD, &status);
          }
        }
        if (bottomProc != MPI_PROC_NULL) {
          // SEND, RECEIVES, AND UPDATES HERE
          MPI_Send(&new_vals, 1, MPI_FLOAT, bottomProc, 3, MPI_COMM_WORLD);
          MPI_Send(&new_vals, 1, MPI_FLOAT, bottomProc, 4, MPI_COMM_WORLD);
          if (procRank == bottomProc) {
            MPI_Recv(&arr[1*block_size+1], 1, MPI_FLOAT, topProc, 3, MPI_COMM_WORLD, &status);
            MPI_Recv(&arr[1*block_size+(j + 1)], 1, MPI_FLOAT, topProc, 4, MPI_COMM_WORLD, &status);
          }
        }
        if (bttmL_Proc != MPI_PROC_NULL) {
          // SEND, RECEIVES, AND UPDATES HERE
          MPI_Send(&new_vals, 1, MPI_FLOAT, bttmL_Proc, 5, MPI_COMM_WORLD);
          if (procRank == bttmL_Proc) {
            MPI_Recv(&arr[1*block_size+(block_size - 2)], 1, MPI_FLOAT, topR_Proc, 5, MPI_COMM_WORLD, &status);
          }
        }
      }

      // top right corner
      else if (i == 1 && j == block_size - 2) {
        // determine # of neighbours and sum their values
        neighbours = 3;
        neighbour_sum += arr[i*block_size+(j-1)];
        neighbour_sum += arr[(i+1)*block_size+(j-1)];
        neighbour_sum += arr[(i+1)*block_size+j];

        if (rightProc != MPI_PROC_NULL) {
          neighbours += 2;
          // SEND, RECEIVES, AND UPDATES HERE
          if (procRank == rightProc) {
            //MPI_Send(&arr[1*block_size+1], 1, MPI_FLOAT, leftProc, 1, MPI_COMM_WORLD);
            //MPI_Send(&arr[(i + 1)*block_size+1], 1, MPI_FLOAT, leftProc, 2, MPI_COMM_WORLD);
            switch (rightProc)
            {
            case 1:
              neighbour_sum += r2[1*block_size+(block_size - 1)];
              neighbour_sum += r2[(i + 1)*block_size+(block_size - 1)];
              break;
            case 2:
              neighbour_sum += r3[1*block_size+(block_size - 1)];
              neighbour_sum += r3[(i + 1)*block_size+(block_size - 1)];
              break;
            case 3:
              neighbour_sum += r4[1*block_size+(block_size - 1)];
              neighbour_sum += r4[(i + 1)*block_size+(block_size - 1)];
              break;
            default:
              break;
            }
          }
          //MPI_Recv(&arr[1*block_size+(block_size - 1)], 1, MPI_FLOAT, rightProc, 1, MPI_COMM_WORLD, &status);
          //MPI_Recv(&arr[(i + 1)*block_size+(block_size - 1)], 1, MPI_FLOAT, rightProc, 2, MPI_COMM_WORLD, &status);
          //neighbour_sum += arr[1*block_size+(block_size - 1)];
          //neighbour_sum += arr[(i + 1)*block_size+(block_size - 1)];
        }
        if (topProc != MPI_PROC_NULL) {
          neighbours += 2;
          // SEND, RECEIVES, AND UPDATES HERE
          if (procRank == topProc) {
            //MPI_Send(&arr[(block_size - 2)*block_size+j], 1, MPI_FLOAT, bottomProc, 3, MPI_COMM_WORLD);
            //MPI_Send(&arr[(block_size - 2)*block_size+(j - 1)], 1, MPI_FLOAT, bottomProc, 4, MPI_COMM_WORLD);
            switch (topProc)
            {
            case 0:
              neighbour_sum += r1[0*block_size+j];
              neighbour_sum += r1[0*block_size+(j - 1)];
              break;
            case 1:
              neighbour_sum += r2[0*block_size+j];
              neighbour_sum += r2[0*block_size+(j - 1)];
              break;
            default:
              break;
            }
          }
          //MPI_Recv(&arr[0*block_size+j], 1, MPI_FLOAT, topProc, 3, MPI_COMM_WORLD, &status);
          //MPI_Recv(&arr[0*block_size+(j - 1)], 1, MPI_FLOAT, topProc, 4, MPI_COMM_WORLD, &status);
          //neighbour_sum += arr[0*block_size+j];
          //neighbour_sum += arr[0*block_size+(j - 1)];
        }
        if (topR_Proc != MPI_PROC_NULL) {
          neighbours += 1;
          // SEND, RECEIVES, AND UPDATES HERE
          if (procRank == topR_Proc) {
            //MPI_Send(&arr[(block_size - 2)*block_size+1], 1, MPI_FLOAT, bttmL_Proc, 5, MPI_COMM_WORLD);
            neighbour_sum += r2[0*block_size+(block_size - 1)];
          }
          //MPI_Recv(&arr[0*block_size+(block_size - 1)], 1, MPI_FLOAT, topR_Proc, 5, MPI_COMM_WORLD, &status);
          //neighbour_sum += arr[0*block_size+(block_size - 1)];
        }

        // calculate average of neighbours
        neighbour_avg = neighbour_sum / neighbours;
        // calculate new values for neighbours (5% of the difference between current cell and average of neighbours)
        new_vals = abs(0.05 * (arr[i*block_size+j] - neighbour_avg));

        // change values of neighbours to new value
        arr[i*block_size+(j-1)] = new_vals;
        arr[(i+1)*block_size+(j-1)] = new_vals;
        arr[(i+1)*block_size+j] = new_vals;

        if (rightProc != MPI_PROC_NULL) {
          // SEND, RECEIVES, AND UPDATES HERE
          MPI_Send(&new_vals, 1, MPI_FLOAT, rightProc, 1, MPI_COMM_WORLD);
          MPI_Send(&new_vals, 1, MPI_FLOAT, rightProc, 2, MPI_COMM_WORLD);
          if (procRank == rightProc) {
            MPI_Recv(&arr[1*block_size+1], 1, MPI_FLOAT, leftProc, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(&arr[(i + 1)*block_size+1], 1, MPI_FLOAT, leftProc, 2, MPI_COMM_WORLD, &status);
          }
        }
        if (topProc != MPI_PROC_NULL) {
          // SEND, RECEIVES, AND UPDATES HERE
          MPI_Send(&new_vals, 1, MPI_FLOAT, topProc, 3, MPI_COMM_WORLD);
          MPI_Send(&new_vals, 1, MPI_FLOAT, topProc, 4, MPI_COMM_WORLD);
          if (procRank == topProc)
          {
            MPI_Recv(&arr[(block_size - 2)*block_size+j], 1, MPI_FLOAT, bottomProc, 3, MPI_COMM_WORLD, &status);
            MPI_Recv(&arr[(block_size - 2)*block_size+(j - 1)], 1, MPI_FLOAT, bottomProc, 4, MPI_COMM_WORLD, &status);
          }
        }
        if (topR_Proc != MPI_PROC_NULL) {
          // SEND, RECEIVES, AND UPDATES HERE
          MPI_Send(&new_vals, 1, MPI_FLOAT, topR_Proc, 5, MPI_COMM_WORLD);
          if (procRank == topR_Proc) {
            MPI_Recv(&arr[(block_size - 2)*block_size+1], 1, MPI_FLOAT, bttmL_Proc, 5, MPI_COMM_WORLD, &status);
          }
        }
      }

      // bottom right corner
      else if (i == block_size - 2 && j == block_size - 2) {
        // determine # of neighbours and sum their values
        neighbours = 3;
        neighbour_sum += arr[i*block_size+(j-1)];
        neighbour_sum += arr[(i-1)*block_size+(j-1)];
        neighbour_sum += arr[(i-1)*block_size+j];

        if (rightProc != MPI_PROC_NULL) {
          neighbours += 2;
          // SEND, RECEIVES, AND UPDATES HERE
          if (procRank == rightProc) {
            //MPI_Send(&arr[(block_size - 2)*block_size+1], 1, MPI_FLOAT, leftProc, 1, MPI_COMM_WORLD);
            //MPI_Send(&arr[(i - 1)*block_size+1], 1, MPI_FLOAT, leftProc, 2, MPI_COMM_WORLD);
            switch (rightProc)
            {
            case 1:
              neighbour_sum += r2[(block_size - 2)*block_size+(block_size - 1)];
              neighbour_sum += r2[(i - 1)*block_size+(block_size - 1)];
              break;
            case 2:
              neighbour_sum += r3[(block_size - 2)*block_size+(block_size - 1)];
              neighbour_sum += r3[(i - 1)*block_size+(block_size - 1)];
              break;
            case 3:
              neighbour_sum += r4[(block_size - 2)*block_size+(block_size - 1)];
              neighbour_sum += r4[(i - 1)*block_size+(block_size - 1)];
              break;
            default:
              break;
            }
          }
          //MPI_Recv(&arr[(block_size - 2)*block_size+(block_size - 1)], 1, MPI_FLOAT, rightProc, 1, MPI_COMM_WORLD, &status);
          //MPI_Recv(&arr[(i - 1)*block_size+(block_size - 1)], 1, MPI_FLOAT, rightProc, 2, MPI_COMM_WORLD, &status);
          //neighbour_sum += arr[(block_size - 2)*block_size+(block_size - 1)];
          //neighbour_sum += arr[(i - 1)*block_size+(block_size - 1)];
        }
        if (bottomProc != MPI_PROC_NULL) {
          neighbours += 2;
          // SEND, RECEIVES, AND UPDATES HERE
          if (procRank == bottomProc) {
            //MPI_Send(&arr[1*block_size+j], 1, MPI_FLOAT, topProc, 3, MPI_COMM_WORLD);
            //MPI_Send(&arr[i*block_size+(j - 1)], 1, MPI_FLOAT, topProc, 4, MPI_COMM_WORLD);
            switch (bottomProc)
            {
            case 2:
              neighbour_sum += r3[(block_size - 1)*block_size+j];
              neighbour_sum += r3[(block_size - 1)*block_size+(j - 1)];
              break;
            case 3:
              neighbour_sum += r4[(block_size - 1)*block_size+j];
              neighbour_sum += r4[(block_size - 1)*block_size+(j - 1)];
              break;
            default:
              break;
            }
          }
          //MPI_Recv(&arr[(block_size - 1)*block_size+j], 1, MPI_FLOAT, bottomProc, 3, MPI_COMM_WORLD, &status);
          //MPI_Recv(&arr[(block_size - 1)*block_size+(j-1)], 1, MPI_FLOAT, bottomProc, 4, MPI_COMM_WORLD, &status);
          //neighbour_sum += arr[(block_size - 1)*block_size+j];
          //neighbour_sum += arr[(block_size - 1)*block_size+(j-1)];
        }
        if (bttmR_Proc != MPI_PROC_NULL) {
          neighbours += 1;
          // SEND, RECEIVES, AND UPDATES HERE
          if (procRank == bttmR_Proc) {
            //MPI_Send(&arr[1*block_size+1], 1, MPI_FLOAT, topL_Proc, 5, MPI_COMM_WORLD);
            neighbour_sum += r4[(block_size - 1)*block_size+(block_size - 1)];
          }
          //MPI_Recv(&arr[(block_size - 1)*block_size+(block_size - 1)], 1, MPI_FLOAT, bttmR_Proc, 5, MPI_COMM_WORLD, &status);
          //neighbour_sum += arr[(block_size - 1)*block_size+(block_size - 1)];
        }

        // calculate average of neighbours
        neighbour_avg = neighbour_sum / neighbours;
        // calculate new values for neighbours (5% of the difference between current cell and average of neighbours)
        new_vals = abs(0.05 * (arr[i*block_size+j] - neighbour_avg));

        // change values of neighbours to new value
        arr[i*block_size+(j-1)] = new_vals;
        arr[(i-1)*block_size+(j-1)] = new_vals;
        arr[(i-1)*block_size+j] = new_vals;

        if (rightProc != MPI_PROC_NULL) {
          // SEND, RECEIVES, AND UPDATES HERE
          MPI_Send(&new_vals, 1, MPI_FLOAT, rightProc, 1, MPI_COMM_WORLD);
          MPI_Send(&new_vals, 1, MPI_FLOAT, rightProc, 2, MPI_COMM_WORLD);
          if (procRank == rightProc) {
            MPI_Recv(&arr[(block_size - 2)*block_size+1], 1, MPI_FLOAT, leftProc, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(&arr[(i - 1)*block_size+1], 1, MPI_FLOAT, leftProc, 2, MPI_COMM_WORLD, &status);
          }
        }
        if (bottomProc != MPI_PROC_NULL) {
          // SEND, RECEIVES, AND UPDATES HERE
          MPI_Send(&new_vals, 1, MPI_FLOAT, bottomProc, 3, MPI_COMM_WORLD);
          MPI_Send(&new_vals, 1, MPI_FLOAT, bottomProc, 4, MPI_COMM_WORLD);
          if (procRank == bottomProc) {
            MPI_Recv(&arr[1*block_size+j], 1, MPI_FLOAT, topProc, 3, MPI_COMM_WORLD, &status);
            MPI_Recv(&arr[i*block_size+(j - 1)], 1, MPI_FLOAT, topProc, 4, MPI_COMM_WORLD, &status);
          }
        }
        if (bttmR_Proc != MPI_PROC_NULL) {
          // SEND, RECEIVES, AND UPDATES HERE
          MPI_Send(&new_vals, 1, MPI_FLOAT, bttmR_Proc, 5, MPI_COMM_WORLD);
          if (procRank == bttmR_Proc) {
            MPI_Recv(&arr[1*block_size+1], 1, MPI_FLOAT, topL_Proc, 5, MPI_COMM_WORLD, &status);
          }
        }
      }

      // left side
      else if (i > 1 && i < block_size - 2 && j == 1) {
        // determine # of neighbours and sum their values
        neighbours = 5;
        neighbour_sum += arr[(i-1)*block_size+j];
        neighbour_sum += arr[(i-1)*block_size+(j+1)];
        neighbour_sum += arr[(i)*block_size+(j+1)];
        neighbour_sum += arr[(i+1)*block_size+j];
        neighbour_sum += arr[(i+1)*block_size+(j+1)];

        if (leftProc != MPI_PROC_NULL) {
          neighbours = 8;
          // GET VALUES SENT FROM NEIGHBOURS TO HALO
          if (procRank == leftProc) {
            //MPI_Send(&arr[i*block_size+(block_size - 2)], 1, MPI_FLOAT, rightProc, 1, MPI_COMM_WORLD);
            //MPI_Send(&arr[(i-1)*block_size+(block_size - 2)], 1, MPI_FLOAT, rightProc, 2, MPI_COMM_WORLD);
            //MPI_Send(&arr[(i+1)*block_size+(block_size - 2)], 1, MPI_FLOAT, rightProc, 3, MPI_COMM_WORLD);
            switch (leftProc)
            {
            case 0:
              neighbour_sum += r1[i*block_size+0];
              neighbour_sum += r1[(i - 1)*block_size+0];
              neighbour_sum += r1[(i + 1)*block_size+0];
              break;
            case 1:
              neighbour_sum += r2[i*block_size+0];
              neighbour_sum += r2[(i - 1)*block_size+0];
              neighbour_sum += r2[(i + 1)*block_size+0];
              break;
            case 2:
              neighbour_sum += r3[i*block_size+0];
              neighbour_sum += r3[(i - 1)*block_size+0];
              neighbour_sum += r2[(i + 1)*block_size+0];
              break;
            default:
              break;
            }
          }
          //MPI_Recv(&arr[i*block_size+0], 1, MPI_FLOAT, leftProc, 1, MPI_COMM_WORLD, &status);
          //MPI_Recv(&arr[(i-1)*block_size+0], 1, MPI_FLOAT, leftProc, 2, MPI_COMM_WORLD, &status);
          //MPI_Recv(&arr[(i+1)*block_size+0], 1, MPI_FLOAT, leftProc, 3, MPI_COMM_WORLD, &status);
          //neighbour_sum += arr[i*block_size+0];
          //neighbour_sum += arr[(i-1)*block_size+0];
          //neighbour_sum += arr[(i+1)*block_size+0];
        }

        // calculate average of neighbours
        neighbour_avg = neighbour_sum / neighbours;
        // calculate new values for neighbours (5% of the difference between current cell and average of neighbours)
        new_vals = abs(0.05 * (arr[i*block_size+j] - neighbour_avg));

        // change values of neighbours to new value
        arr[(i-1)*block_size+j] = new_vals;
        arr[(i-1)*block_size+(j+1)] = new_vals;
        arr[i*block_size+(j+1)] = new_vals;
        arr[(i+1)*block_size+j] = new_vals;
        arr[(i+1)*block_size+(j+1)] = new_vals;

        if (leftProc != MPI_PROC_NULL) {
          // UPDATE NEIGHBOUR VALUES IN TOP PROCESS
          MPI_Send(&new_vals, 1, MPI_FLOAT, leftProc, 1, MPI_COMM_WORLD);
          MPI_Send(&new_vals, 1, MPI_FLOAT, leftProc, 2, MPI_COMM_WORLD);
          MPI_Send(&new_vals, 1, MPI_FLOAT, leftProc, 3, MPI_COMM_WORLD);
          if (procRank == leftProc) {
            MPI_Recv(&arr[i*block_size+(block_size - 2)], 1, MPI_FLOAT, rightProc, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(&arr[(i-1)*block_size+(block_size - 2)], 1, MPI_FLOAT, rightProc, 2, MPI_COMM_WORLD, &status);
            MPI_Recv(&arr[(i+1)*block_size+(block_size - 2)], 1, MPI_FLOAT, rightProc, 3, MPI_COMM_WORLD, &status);
          }
        }
      }

      // right side
      else if (i < block_size - 2 && i > 1 && j == block_size - 2) {
        // determine # of neighbours and sum their values
        neighbours = 5;
        neighbour_sum += arr[(i-1)*block_size+j];
        neighbour_sum += arr[(i-1)*block_size+(j-1)];
        neighbour_sum += arr[i*block_size+(j-1)];
        neighbour_sum += arr[(i+1)*block_size+j];
        neighbour_sum += arr[(i+1)*block_size+(j-1)];

        if (rightProc != MPI_PROC_NULL) {
          neighbours = 8;
          // GET VALUES SENT FROM NEIGHBOURS TO HALO
          if (procRank == rightProc) {
            //MPI_Send(&arr[i*block_size+1], 1, MPI_FLOAT, leftProc, 1, MPI_COMM_WORLD);
            //MPI_Send(&arr[(i-1)*block_size+1], 1, MPI_FLOAT, leftProc, 2, MPI_COMM_WORLD);
            //MPI_Send(&arr[(i+1)*block_size+1], 1, MPI_FLOAT, leftProc, 3, MPI_COMM_WORLD);
            switch (rightProc)
            {
            case 1:
              neighbour_sum += r2[i*block_size+(block_size - 1)];
              neighbour_sum += r2[(i-1)*block_size+(block_size - 1)];
              neighbour_sum += r2[(i+1)*block_size+(block_size - 1)];
              break;
            case 2:
              neighbour_sum += r3[i*block_size+(block_size - 1)];
              neighbour_sum += r3[(i-1)*block_size+(block_size - 1)];
              neighbour_sum += r3[(i+1)*block_size+(block_size - 1)];
              break;
            case 3:
              neighbour_sum += r4[i*block_size+(block_size - 1)];
              neighbour_sum += r4[(i-1)*block_size+(block_size - 1)];
              neighbour_sum += r4[(i+1)*block_size+(block_size - 1)];
              break;
            default:
              break;
            }
          }
          //MPI_Recv(&arr[i*block_size+(block_size - 1)], 1, MPI_FLOAT, rightProc, 1, MPI_COMM_WORLD, &status);
          //MPI_Recv(&arr[(i-1)*block_size+(block_size - 1)], 1, MPI_FLOAT, rightProc, 2, MPI_COMM_WORLD, &status);
          //MPI_Recv(&arr[(i+1)*block_size+(block_size - 1)], 1, MPI_FLOAT, rightProc, 3, MPI_COMM_WORLD, &status);
          //neighbour_sum += arr[i*block_size+(block_size - 1)];
          //neighbour_sum += arr[(i-1)*block_size+(block_size - 1)];
          //neighbour_sum += arr[(i+1)*block_size+(block_size - 1)];
        }

        // calculate average of neighbours
        neighbour_avg = neighbour_sum / neighbours;
        // calculate new values for neighbours (5% of the difference between current cell and average of neighbours)
        new_vals = abs(0.05 * (arr[i*block_size+j] - neighbour_avg));

        // change values of neighbours to new value
        arr[(i-1)*block_size+j] = new_vals;
        arr[(i-1)*block_size+(j-1)] = new_vals;
        arr[i*block_size+(j-1)] = new_vals;
        arr[(i+1)*block_size+j] = new_vals;
        arr[(i+1)*block_size+(j-1)] = new_vals;

        if (rightProc != MPI_PROC_NULL) {
          // UPDATE NEIGHBOUR VALUES IN TOP PROCESS
          MPI_Send(&new_vals, 1, MPI_FLOAT, rightProc, 1, MPI_COMM_WORLD);
          MPI_Send(&new_vals, 1, MPI_FLOAT, rightProc, 2, MPI_COMM_WORLD);
          MPI_Send(&new_vals, 1, MPI_FLOAT, rightProc, 3, MPI_COMM_WORLD);
          if (procRank == rightProc) {
            MPI_Recv(&arr[i*block_size+1], 1, MPI_FLOAT, leftProc, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(&arr[(i-1)*block_size+1], 1, MPI_FLOAT, leftProc, 2, MPI_COMM_WORLD, &status);
            MPI_Recv(&arr[(i+1)*block_size+1], 1, MPI_FLOAT, leftProc, 3, MPI_COMM_WORLD, &status);
          }
        }
      }

      // top side
      else if (i == 1 && j < block_size - 2 && j > 1) {
        // determine # of neighbours and sum their values
        neighbours = 5;
        neighbour_sum += arr[i*block_size+(j-1)];
        neighbour_sum += arr[i*block_size+(j+1)];
        neighbour_sum += arr[(i+1)*block_size+(j-1)];
        neighbour_sum += arr[(i+1)*block_size+j];
        neighbour_sum += arr[(i+1)*block_size+(j+1)];

        if (topProc != MPI_PROC_NULL) {
          neighbours = 8;
          // GET VALUES SENT FROM NEIGHBOURS TO HALO
          if (procRank == topProc) {
            //MPI_Send(&arr[(block_size - 2)*block_size+j], 1, MPI_FLOAT, bottomProc, 1, MPI_COMM_WORLD);
            //MPI_Send(&arr[(block_size - 2)*block_size+(j-1)], 1, MPI_FLOAT, bottomProc, 2, MPI_COMM_WORLD);
            //MPI_Send(&arr[(block_size - 2)*block_size+(j+1)], 1, MPI_FLOAT, bottomProc, 3, MPI_COMM_WORLD);
            switch (topProc)
            {
            case 0:
              neighbour_sum += r1[(i-1)*block_size+j];
              neighbour_sum += r1[(i-1)*block_size+(j-1)];
              neighbour_sum += r1[(i-1)*block_size+(j+1)];
              break;
            case 1:
              neighbour_sum += r2[(i-1)*block_size+j];
              neighbour_sum += r2[(i-1)*block_size+(j-1)];
              neighbour_sum += r2[(i-1)*block_size+(j+1)];
              break;
            default:
              break;
            }
          }
          //MPI_Recv(&arr[(i-1)*block_size+j], 1, MPI_FLOAT, topProc, 1, MPI_COMM_WORLD, &status);
          //MPI_Recv(&arr[(i-1)*block_size+(j-1)], 1, MPI_FLOAT, topProc, 2, MPI_COMM_WORLD, &status);
          //MPI_Recv(&arr[(i-1)*block_size+(j+1)], 1, MPI_FLOAT, topProc, 3, MPI_COMM_WORLD, &status);
          //neighbour_sum += arr[(i-1)*block_size+j];
          //neighbour_sum += arr[(i-1)*block_size+(j-1)];
          //neighbour_sum += arr[(i-1)*block_size+(j+1)];
        }

        // calculate average of neighbours
        neighbour_avg = neighbour_sum / neighbours;
        // calculate new values for neighbours (5% of the difference between current cell and average of neighbours)
        new_vals = abs(0.05 * (arr[i*block_size+j] - neighbour_avg));

        // change values of neighbours to new value
        arr[i*block_size+(j-1)] = new_vals;
        arr[i*block_size+(j+1)] = new_vals;
        arr[(i+1)*block_size+(j-1)] = new_vals;
        arr[(i+1)*block_size+j] = new_vals;
        arr[(i+1)*block_size+(j+1)] = new_vals;

        if (topProc != MPI_PROC_NULL) {
          // UPDATE NEIGHBOUR VALUES IN TOP PROCESS
          MPI_Send(&new_vals, 1, MPI_FLOAT, topProc, 1, MPI_COMM_WORLD);
          MPI_Send(&new_vals, 1, MPI_FLOAT, topProc, 2, MPI_COMM_WORLD);
          MPI_Send(&new_vals, 1, MPI_FLOAT, topProc, 3, MPI_COMM_WORLD);
          if (procRank == topProc) {
            MPI_Recv(&arr[(block_size - 2)*block_size+j], 1, MPI_FLOAT, bottomProc, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(&arr[(block_size - 2)*block_size+(j-1)], 1, MPI_FLOAT, bottomProc, 2, MPI_COMM_WORLD, &status);
            MPI_Recv(&arr[(block_size - 2)*block_size+(j+1)], 1, MPI_FLOAT, bottomProc, 3, MPI_COMM_WORLD, &status);
          }
        }
      }

      // bottom side
      else if (i == block_size - 2 && j > 1 && j < block_size - 2) {
        // determine # of neighbours and sum their values
        neighbours = 5;
        neighbour_sum += arr[i*block_size+(j-1)];
        neighbour_sum += arr[i*block_size+(j+1)];
        neighbour_sum += arr[(i-1)*block_size+(j-1)];
        neighbour_sum += arr[(i-1)*block_size+j];
        neighbour_sum += arr[(i-1)*block_size+(j+1)];

        if (bottomProc != MPI_PROC_NULL) { 
          neighbours = 8;
          // GET VALUES SENT FROM NEIGHBOURS TO HALO
          if (procRank == bottomProc) {
            //MPI_Send(&arr[1*block_size+j], 1, MPI_FLOAT, topProc, 1, MPI_COMM_WORLD);
            //MPI_Send(&arr[1*block_size+(j-1)], 1, MPI_FLOAT, topProc, 2, MPI_COMM_WORLD);
            //MPI_Send(&arr[1*block_size+(j+1)], 1, MPI_FLOAT, topProc, 3, MPI_COMM_WORLD);
            switch (bottomProc)
            {
            case 2:
              neighbour_sum += r3[(i+1)*block_size+j];
              neighbour_sum += r3[(i+1)*block_size+(j-1)];
              neighbour_sum += r3[(i+1)*block_size+(j+1)];
              break;
            case 3:
              neighbour_sum += r4[(i+1)*block_size+j];
              neighbour_sum += r4[(i+1)*block_size+(j-1)];
              neighbour_sum += r4[(i+1)*block_size+(j+1)];
              break;
            default:
              break;
            }
          }
          //MPI_Recv(&arr[(i+1)*block_size+j], 1, MPI_FLOAT, bottomProc, 1, MPI_COMM_WORLD, &status);
          //MPI_Recv(&arr[(i+1)*block_size+(j-1)], 1, MPI_FLOAT, bottomProc, 2, MPI_COMM_WORLD, &status);
          //MPI_Recv(&arr[(i+1)*block_size+(j+1)], 1, MPI_FLOAT, bottomProc, 3, MPI_COMM_WORLD, &status);
          //neighbour_sum += arr[(i+1)*block_size+j];
          //neighbour_sum += arr[(i+1)*block_size+(j-1)];
          //neighbour_sum += arr[(i+1)*block_size+(j+1)];
        }

        // calculate average of neighbours
        neighbour_avg = neighbour_sum / neighbours;
        // calculate new values for neighbours (5% of the difference between current cell and average of neighbours)
        new_vals = abs(0.05 * (arr[i*block_size+j] - neighbour_avg));

        // change values of neighbours to new value
        arr[i*block_size+(j-1)] = new_vals;
        arr[i*block_size+(j+1)] = new_vals;
        arr[(i-1)*block_size+(j-1)] = new_vals;
        arr[(i-1)*block_size+j] = new_vals;
        arr[(i-1)*block_size+(j+1)] = new_vals;

        if (bottomProc != MPI_PROC_NULL) {
          // UPDATE NEIGHBOUR VALUES IN BOTTOM PROCESS
          MPI_Send(&new_vals, 1, MPI_FLOAT, bottomProc, 1, MPI_COMM_WORLD);
          MPI_Send(&new_vals, 1, MPI_FLOAT, bottomProc, 2, MPI_COMM_WORLD);
          MPI_Send(&new_vals, 1, MPI_FLOAT, bottomProc, 3, MPI_COMM_WORLD);
          if (procRank == bottomProc) {
            MPI_Recv(&arr[1*block_size+j], 1, MPI_FLOAT, topProc, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(&arr[1*block_size+(j-1)], 1, MPI_FLOAT, topProc, 2, MPI_COMM_WORLD, &status);
            MPI_Recv(&arr[1*block_size+(j+1)], 1, MPI_FLOAT, topProc, 3, MPI_COMM_WORLD, &status);
          }
          MPI_Sendrecv(&new_vals, 1, MPI_FLOAT, bottomProc, 0, &arr[1*block_size+j], 1, MPI_FLOAT, procRank, 0, MPI_COMM_WORLD, &status);
          MPI_Sendrecv(&new_vals, 1, MPI_FLOAT, bottomProc, 0, &arr[1*block_size+(j-1)], 1, MPI_FLOAT, procRank, 0, MPI_COMM_WORLD, &status);
          MPI_Sendrecv(&new_vals, 1, MPI_FLOAT, bottomProc, 0, &arr[1*block_size+(j+1)], 1, MPI_FLOAT, procRank, 0, MPI_COMM_WORLD, &status);
        }
      }

      // inner cells
      else {
        // determine # of neighbours and sum their values
        neighbours = 8;
        neighbour_sum += arr[(i-1)*block_size+(j-1)];
        neighbour_sum += arr[(i-1)*block_size+j];
        neighbour_sum += arr[(i-1)*block_size+(j+1)];
        neighbour_sum += arr[i*block_size+(j-1)];
        neighbour_sum += arr[i*block_size+(j+1)];
        neighbour_sum += arr[(i+1)*block_size+(j-1)];
        neighbour_sum += arr[(i+1)*block_size+j];
        neighbour_sum += arr[(i+1)*block_size+(j+1)];

        // calculate average of neighbours
        neighbour_avg = neighbour_sum / neighbours;
        // calculate new values for neighbours (5% of the difference between current cell and average of neighbours)
        new_vals = abs(0.05 * (arr[i*block_size+j] - neighbour_avg));

        // change values of neighbours to new value
        arr[(i-1)*block_size+(j-1)] = new_vals;
        arr[(i-1)*block_size+j] = new_vals;
        arr[(i-1)*block_size+(j+1)] = new_vals;
        arr[i*block_size+(j-1)] = new_vals;
        arr[i*block_size+(j+1)] = new_vals;
        arr[(i+1)*block_size+(j-1)] = new_vals;
        arr[(i+1)*block_size+j] = new_vals;
        arr[(i+1)*block_size+(j+1)] = new_vals;
      }

      // if at the end of a row, go to the next row
      if (i != block_size - 2 && j == block_size - 2) {
        i++;
        j = 1;
      }

      // if at the end of the matrix, go to the beginning
      else if (i == block_size - 2 && j == block_size - 2) {
        i = 1;
        j = 1;
      }

      // go to the next column
      else {
        j++;
      }
    } while (i != max_i || j != max_j);

    // determine largest number in array
    max = arr[1*block_size+1];
    for (int m = 1; m < block_size - 1; m++) {
      for (int n = 1; n < block_size - 1; n++) {
        if (arr[m*block_size+n] > max) {
          max = arr[m*block_size+n];
          i = m;
          max_i = m;
          j = n;
          max_j = n;
        }
      }
    }
    count += 1;
  } while (max > 12);

  printf("Process #%d had %d iterations.\n", procRank, count);
  // delete 2D array and exit
  delete [] arr;
  MPI_Finalize();
  return 0;
}