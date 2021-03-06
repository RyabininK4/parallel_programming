#include "mpi.h"
#include <stdlib.h> 
#include <stack>
#include <iostream>
#include <stdio.h>
#include <time.h>

using std::cout;
using std::endl;

int rank;
int size;

int MPI_Bcast_Tree(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm);
int parent(int your_num, int child);
int level(int your_num, int child);

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int const n = 30000000;
    int*  bufInt = new int[n];
    float* bufFloat = new float[n];
    double* bufDou = new double[n];
    for (int i = 0; i < n; i++) {
        bufInt[n] = rand() % 10;
        bufFloat[n] = ((float)rand() / RAND_MAX * (20 - 10) + -10);
        bufDou[n] = ((double)rand() / RAND_MAX * 10);
    }
    double time;
    if (rank == 0) {
        time = MPI_Wtime();

    }

    MPI_Bcast_Tree(&bufInt, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast_Tree(&bufFloat, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast_Tree(&bufDou, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (rank == 0) {
    time -= MPI_Wtime();
    printf("%lf15", time*(-1));
    cout << endl;
    time = MPI_Wtime();
    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Bcast(&bufInt, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&bufFloat, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&bufDou, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        time -= MPI_Wtime();
        printf("%lf15", time*(-1));
        cout << endl;
    }
    /* 
     if (rank == 0) {
        cout << "my rank : " << rank << endl;
        cout << "int : ";
        for (int i = 0; i < 1; i++) {
            cout << bufInt[i] << " ";
        }
        cout << endl << "float : ";
        for (int i = 0; i < 1; i++) {
            cout << bufFloat[i] << " ";
        }
        cout << endl << "double : ";
        for (int i = 0; i < 1; i++) {
            cout << bufDou[i] << " ";
        }
        cout << endl;
    }*/
    MPI_Finalize();
    return 0;
}

int parent(int your_num, int child) {
    return ((your_num - 1) / child);
}

int level(int your_num, int child) {
    int lvl = 2;
    int num = your_num;
    while (parent(num, child) != 0) {
        num = parent(num, child);
        lvl += 1;
    }
    return lvl;
}

int MPI_Bcast_Tree(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm) {
    MPI_Barrier(comm);
    int child = 1;
    int i;
    void* receiver = malloc(sizeof(datatype)*count);
    MPI_Status status;

    if (rank == root) {
        for (i = 1; ((i <= child) && (i < size)); i++) {
            if (i != root)
                MPI_Send(buffer, count, datatype, i, 0, comm);
            else
                MPI_Send(buffer, count, datatype, 0, 0, comm); 
        }
    } 
    else {
        if (rank != 0) {
            if (level(rank, child) == 2) {
                MPI_Recv(receiver, count, datatype, root, 0, comm, &status);
            }
            else {
                if (parent(rank, child) == root) 
                    MPI_Recv(receiver, count, datatype, 0, 0, comm, &status);
                else
                    MPI_Recv(receiver, count, datatype, parent(rank, child), 0, comm, &status); 
            }

            int dest = rank*child + 1;

            for (i = 1; ((i <= child) && (dest < size)); i++) {
                if (dest != root)
                    MPI_Send(receiver, count, datatype, dest, 0, comm);
                else
                    MPI_Send(receiver, count, datatype, 0, 0, comm);
                dest += 1;
            }
        }
        else {
            if (root <= child) {
                MPI_Recv(receiver, count, datatype, root, 0, comm, &status);
            }
            else {
                MPI_Recv(receiver, count, datatype, parent(root, child), 0, comm, &status);
            }
            int dest = root*child + 1;
            for (i = 1; ((i <= child) && (dest < size)); i++) {
                MPI_Send(receiver, count, datatype, dest, 0, comm);
                dest += 1;
            }
        }
        
       // cout << "my rank : " << rank << endl; 
    }
    return 0;
} 