#include "mpi.h"
#include <stdlib.h> 
#include <stack>
#include <iostream>
#include <stdio.h>

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
    int* bufInt = new int[5];
    float* bufFloat = new float[5];
    double* bufDou = new double[5];
    for (int i = 0; i < 5; i++) {
        bufInt[i] = i;
        bufFloat[i] = double(i) / 10;
        bufDou[i] = float(i) / 7;
    }
    
    MPI_Bcast_Tree(&bufInt, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast_Tree(&bufFloat, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast_Tree(&bufDou, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    //вывод из 1 процесса
    if (rank == 0) {
        cout << "int : ";
        for (int i = 0; i < 5; i++) {
            cout << bufInt[i] << " ";
        }
        cout << endl << "float : ";
        for (int i = 0; i < 5; i++) {
            cout << bufFloat[i] << " ";
        }
        cout << endl << "double : ";
        for (int i = 0; i < 5; i++) {
            cout << bufDou[i] << " ";
        }
        cout << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    //вывод из 2 процесса
    if (rank == 1) {
        cout << "int : ";
        for (int i = 0; i < 5; i++) {
            cout << bufInt[i] << " ";
        }
        cout << endl << "float : ";
        for (int i = 0; i < 5; i++) {
            cout << bufFloat[i] << " ";
        }
        cout << endl << "double : ";
        for (int i = 0; i < 5; i++) {
            cout << bufDou[i] << " ";
        }
        cout << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    //вывод из 3 процесса
    if (rank == 2) {
        cout << "int : ";
        for (int i = 0; i < 5; i++) {
            cout << bufInt[i] << " ";
        }
        cout << endl << "float : ";
        for (int i = 0; i < 5; i++) {
            cout << bufFloat[i] << " ";
        }
        cout << endl << "double : ";
        for (int i = 0; i < 5; i++) {
            cout << bufDou[i] << " ";
        }
        cout << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    //вывод из 4 процесса
    if (rank == 3) {
        cout << "int : ";
        for (int i = 0; i < 5; i++) {
            cout << bufInt[i] << " ";
        }
        cout << endl << "float : ";
        for (int i = 0; i < 5; i++) {
            cout << bufFloat[i] << " ";
        }
        cout << endl << "double : ";
        for (int i = 0; i < 5; i++) {
            cout << bufDou[i] << " ";
        }
        cout << endl;
    }
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
        cout << "my rank : " << rank << endl;
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
        cout << "my rank : " << rank << endl;
    }
    return 0;
} 