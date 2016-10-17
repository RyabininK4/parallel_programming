#include <mpi.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>

using std::ifstream;
using std::cout;
using std::string;

int main(int argc, char** argv)
{
    int procNum, procRank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);

    int lenStr = 0; // длина файла
    int count = 0; // количество символов в файле
    char* mas; //массив символов из файла

    double start, stop;

    if (procRank == 0) {
        ifstream file("book1.txt");

        while (!file.eof()) {
            lenStr++;
            file.get();
        }
        file.close();
    }

    MPI_Bcast(&lenStr, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    mas = new char[lenStr];

    if (procRank == 0) {
        ifstream file("book1.txt");

        for (int i = 0; i < lenStr; i++)
            mas[i] = file.get();

        file.close();
    }

    MPI_Bcast(mas, lenStr, MPI_CHAR, 0, MPI_COMM_WORLD);

    if (procRank == 0)
        start = MPI_Wtime();

    int i1, i2;
    int k = lenStr / procNum;

    i1 = (int)(k * procRank);
    i2 = (int)(k * (procRank + 1));

    if (procRank == procNum - 1)
        i2 = lenStr;

    for (int j = i1; j < i2; j++)
        if ((' ' != mas[j]) && ('\0' != mas[j]) && ('\n' != mas[j]) && ('\t' != mas[j]) && ((int)mas[j] >= 33) && ((int)mas[j] <= 255))
            count++;


    MPI_Reduce(&count, &mas, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (procRank == 0) {
        stop = MPI_Wtime();
        cout << "Time for search - " << std::fixed << std::setprecision(3) << stop - start << " sec\n";
    }

    MPI_Finalize();
    delete mas;
    return 0;
}