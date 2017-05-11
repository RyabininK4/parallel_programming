#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <string>
#include <math.h>
#include <stack>
#include <valarray>
#include <vector>
#include <time.h>
#include <algorithm>

using namespace std;

int findOverlappingPair(string str1, string str2, string& str)
{
    int max = INT_MIN;
    int len1 = str1.length();
    int len2 = str2.length();

    for (int i = 1; i <= min(len1, len2); i++)
    {
        if (str1.compare(len1 - i, i, str2, 0, i) == 0)
        {
            if (max < i)
            {
                max = i;
                str = str1 + str2.substr(i);
            }
        }
    }

    for (int i = 1; i <= min(len1, len2); i++)
    {
        if (str1.compare(0, i, str2, len2 - i, i) == 0)
        {
            if (max < i)
            {
                max = i;
                str = str2 + str1.substr(i);
            }
        }
    }

    return max;
}

int main(int argc, char* argv[])
{
    vector<string> arr = { "abc", "bcd"};
    for (int i = 0; i < 1000; i++) {
        arr.push_back("er");
    }

    int len = arr.size();
    double start, stop;

    int procNum, procRank;

    omp_set_num_threads(1);
    cout << omp_get_max_threads() << endl;
    MPI_Init(&argc, &argv);


    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);

    int border1;
    int border2;
    int k = arr.size() / procNum;

    border1 = (int)(k * procRank);
    border2 = (int)(k * (procRank + 1));

    int l, r;
    string resStr;
    int max;
    start = MPI_Wtime();
    while (len != 1)
    {
        max = INT_MIN;
#pragma omp parallel for 
        for (int i = border1; i < border2; i++)
        {
            for (int j = i + 1; j < len; j++)
            {
                string str;

                int res = findOverlappingPair(arr[i], arr[j], str);

                if (max < res)
#pragma omp critical 
                {
                    max = res;
                    resStr = str;
                    l = i, r = j;
                }
            }
        }
        len--;
        if (max == INT_MIN)
            arr[0] += arr[len];
        else
        {
            arr[l] = resStr;
            arr[r] = arr[len];
        }
    }
    stop = MPI_Wtime();

    if (procRank == 0) {
        cout << "Result = " << arr[0] << endl;
        cout << "Time =  " << stop - start << " sec\n";
    }

    MPI_Finalize();
    return 0;
}