#include <mpi.h>
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

int procSize, procRank;

double *genMatrix(int sizeX, int sizeY) { // генерация случайной матрицы
    srand(time(0));
    double *matrix = new double[sizeX * sizeY];
    for (int i = 0; i < sizeX * sizeY; i++)
    {
        matrix[i] = (double)(rand() % 10);
    }
    for (int i = 0; i < sizeX*sizeY; i++)
    {
        if (matrix[i] == 0)
            matrix[i] = (double)(1);
    }
    return matrix;
}

void printMatrix(double* matrix, int sizeX, int sizeY) { //печатаем матрицу

    for (int i = 0; i < sizeY; i++) {
        for (int j = 0; j < sizeX; j++) {
            printf("%6.2f ", matrix[i*sizeX + j]);
        }
        std::cout << "\n";
    }
    std::cout << "\n";

}

// однопроцессорный вариант
double determOne(const double* _matrix, int size) {
    double* matrix = new double[size * size];
    memcpy(matrix, _matrix, sizeof(double)* size * size); 
    int swapCount = 0; 
    double *tempStr = new double[size]; 
    for (int i = 0; i < size - 1; i++) {
        int maxIndex = i; 
        double maxCoeff = matrix[i*size + i]; 
        for (int j = i; j < size; j++) { 
            if ((maxCoeff * maxCoeff) < (matrix[j*size + i] * matrix[j*size + i])) { 
                maxCoeff = matrix[j*size + i];
                maxIndex = j;
            }
        }
        if (maxIndex != i) { 
            swapCount++; 
            memcpy(tempStr, &matrix[i * size], sizeof(double)* size); 
            memcpy(&matrix[i * size], &matrix[maxIndex * size], sizeof(double)* size); 
            memcpy(&matrix[maxIndex * size], tempStr, sizeof(double)* size); 
        }
        double dEl = matrix[i*size + i];
        for (int j = i + 1; j < size; j++) {                                   
            double k = matrix[j*size + i] / dEl;
            for (int l = 0; l < size; l++) {
                matrix[j*size + l] -= matrix[i*size + l] * k;
            }
        }
    }
    double det = matrix[0];
    for (int i = 1; i < size; i++) {
        det *= matrix[i*size + i];
    }

    return det * (((int)swapCount % 2) ? -1 : 1);
}

// параллельный вариант расчета детерминанта
double determMPI(double* matrix, int size) {
    double det = 1;
    int *countMap = 0; 
    int *offsetMap = 0; 
    int swapCount = 0; 
    if (procRank == 0) { 
        countMap = new int[procSize]; 
        offsetMap = new int[procSize];
        for (int i = 0; i < procSize; i++) {
            countMap[i] = ((size / procSize) * size);                                               
            offsetMap[i] = (i * (((size / procSize) * size)));                                                      
        }
        countMap[procSize - 1] += ((size % procSize)* size); 
    }

    double *part, *tempStr;
    int pSize = (size / procSize); 
    int pSizeEnd = (size % procSize); 
    int pSizeThis = (procRank != (procSize - 1)) ? pSize : (pSize + pSizeEnd); 
    part = new double[pSizeThis * size];
    tempStr = new double[size];

    MPI_Scatterv(matrix, countMap, offsetMap, MPI_DOUBLE, part, (pSizeThis * size), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
    if (procRank == 0) {
        delete[] countMap;
        delete[] offsetMap;
    }

    for (int i = 0; i < size; i++) {                               
        int iProc = i / pSize;
        iProc = (iProc >(procSize - 1)) ? (procSize - 1) : iProc; 
        int *maxIndexTemp = 0; 
        double *maxElTemp = 0;
        if (procRank == iProc) 
        {
            maxIndexTemp = new int[procSize];
            maxElTemp = new double[procSize];
        }
        int startIndexP, maxIndexP;
        double maxElP;
        if (procRank >= iProc) 
        { 
            startIndexP = (iProc == procRank) ? (i - (iProc * pSize)) : 0;
            maxElP = part[startIndexP * size + i];
            maxIndexP = startIndexP; 
            for (int j = startIndexP; j < pSizeThis; j++) {
                if (maxElP * maxElP < part[j * size + i] * part[j * size + i]) {
                    maxElP = part[j * size + i];
                    maxIndexP = j;
                }
            }
        }
       
        MPI_Gather(&maxElP, 1, MPI_DOUBLE, maxElTemp, 1, MPI_DOUBLE, iProc, MPI_COMM_WORLD);
        MPI_Gather(&maxIndexP, 1, MPI_INT, maxIndexTemp, 1, MPI_INT, iProc, MPI_COMM_WORLD);

        
        int maxElIndex = i; 
        if (procRank == iProc) {
            double maxElVal = part[((i - (iProc * pSize)) * size) + i]; 
            int maxElProc = procRank; 
            for (int j = procRank; j < procSize; j++) {
                if (maxElVal * maxElVal < maxElTemp[j] * maxElTemp[j]) {
                    maxElVal = maxElTemp[j];
                    maxElIndex = j * pSize + maxIndexTemp[j];
                    maxElProc = j;
                }
            }
        }

        MPI_Bcast(&maxElIndex, 1, MPI_INT, iProc, MPI_COMM_WORLD);
        int iProcMax = maxElIndex / pSize;
        iProcMax = (iProcMax >(procSize - 1)) ? (procSize - 1) : iProcMax;
        if (iProcMax == iProc) {
            if (maxElIndex != i) {
                if (procRank == iProc) {
                    memcpy(tempStr, &part[(i - (iProc * pSize)) * size], sizeof(double)* size);
                    memcpy(&part[(i - (iProc * pSize)) * size], &part[(i - (iProcMax * pSize)) * size], sizeof(double)* size);
                    memcpy(&part[(i - (iProcMax * pSize)) * size], tempStr, sizeof(double)* size);
                }
            }
        }
        else {
            MPI_Status status;
            if (procRank == iProcMax) {
                MPI_Send(&part[(maxElIndex - (iProcMax * pSize)) * size], size, MPI_DOUBLE, iProc, i, MPI_COMM_WORLD);
            }
            if (procRank == iProc) {

                memcpy(tempStr, &part[(i - (iProc * pSize)) * size], sizeof(double)* size);
                MPI_Recv(&part[(i - (iProc * pSize)) * size], size, MPI_DOUBLE, iProcMax, i, MPI_COMM_WORLD, &status);
                MPI_Send(tempStr, size, MPI_DOUBLE, iProcMax, i, MPI_COMM_WORLD);
            }
            if (procRank == iProcMax) {
                swapCount++;
                MPI_Recv(&part[(maxElIndex - (iProcMax * pSize)) * size], size, MPI_DOUBLE, iProc, i, MPI_COMM_WORLD, &status);
            }  
        }
        
        countMap = new int[procSize];
        offsetMap = new int[procSize];
        
        for (int j = 0; j < procSize; j++) {
            countMap[j] = ((j < procRank) ? 1 : size);
            offsetMap[j] = ((i - (iProc * pSize)) * size); 
        }

        MPI_Scatterv(part, countMap, offsetMap, MPI_DOUBLE, tempStr, size, MPI_DOUBLE, iProc, MPI_COMM_WORLD);
        delete[] countMap;
        delete[] offsetMap;
        for (int j = 0; j < pSizeThis; j++) {
            if ((procRank * pSize + j) > i) { 
                double k = part[j*size + i] / tempStr[i];
                for (int l = 0; l < size; l++) {
                    part[j*size + l] -= tempStr[l] * k;
                }
            }
        }
    }
    double pDet = 1;
    for (int i = 0; i < size; i++) {

        int iProc = i / pSize;
        iProc = (iProc >(procSize - 1)) ? (procSize - 1) : iProc;

        if (procRank == iProc) {
            pDet *= part[(i - (pSize * procRank))*size + i];
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&pDet, &det, 1, MPI_DOUBLE, MPI_PROD, 0, MPI_COMM_WORLD);
    int allCount = 0;
    MPI_Reduce(&swapCount, &allCount, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    return det * ((allCount % 2) ? -1 : 1);
}

int main(int argc, char* argv[]) {
    int size = 5; 
    double *matrix = 0; 
    double timeOne, timeMPI; 
    double detOne, detMPI;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &procSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);

    MPI_Barrier(MPI_COMM_WORLD);
    if (procRank == 0) 
    { 
        matrix = genMatrix(size, size);

        fprintf(stdout, "Matrix: %d X %d\n", size, size);
        std::cout << std::endl;
        fflush(stdout);
        printMatrix(matrix, size, size);

        timeOne = MPI_Wtime();
        detOne = determOne(matrix, size);
        timeOne = MPI_Wtime() - timeOne;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    timeMPI = MPI_Wtime();
    detMPI = determMPI(matrix, size);
    timeMPI = MPI_Wtime() - timeMPI;

    if (procRank == 0) {

        fprintf(stdout, "The determinant of the matrix counted serial version: %f \n", detOne);
        fflush(stdout);

        fprintf(stdout, "Time serial version: %f \n", timeOne);
        fflush(stdout);

        fprintf(stdout, "The determinant of the matrix counted parallel version: %f \n", detMPI);
        fflush(stdout);

        fprintf(stdout, "Time parallel version: %f \n", timeMPI);
        fflush(stdout);
    }

    MPI_Finalize();
    return 0;
}