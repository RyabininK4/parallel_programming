#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <ctime>
#include <stack>
#include <string>
#include <valarray>
#include <vector>
#include <time.h>
#include <algorithm>
#include "tbb/task_scheduler_init.h"
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/task_scheduler_init.h"
#include "tbb/spin_mutex.h"


using namespace std;
using namespace tbb;

static vector<string> Arr;
static int Max;
static string ResStr;
static int L, R;
static int result;

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

class CalcParallel
{
public:
    static spin_mutex myMutex;
    void operator() (const blocked_range<int>& range) const {
        Max = INT_MIN;
        for (int i = 0; i < range.end(); i++)
        {
            for (int j = i + 1; j < result; j++)
            {
                string str;
                int res = findOverlappingPair(Arr[i], Arr[j], str);
                
                if (Max < res)
                {
                    spin_mutex::scoped_lock lock(myMutex);
                    Max = res;
                    ResStr = str;
                    L = i, R = j;
                }
            }
        }
        
    }
};

void findShortestSuperstringParalell(vector<string> arr)
{
    result = Arr.size();

    while (result != 1)
    {
        parallel_for(blocked_range<int>(0, result), CalcParallel());
        result--;
        if (Max == INT_MIN)
            Arr[0] += Arr[result];
        else
        {
            Arr[L] = ResStr;
            Arr[R] = Arr[result];
        }
    }
    cout << Arr[0] << endl;
}

spin_mutex CalcParallel::myMutex;

int main()
{
    vector<string> arr = { "abc", "cde"};
    for (int i = 0; i < 150; i++) {
     arr.push_back("er");
    }
    Arr = arr;
    
    clock_t begin_time1 = clock();
    task_scheduler_init init(2);
    findShortestSuperstringParalell(Arr);
    clock_t stop_time1 = clock();
    cout << "Parallel time: " << double(stop_time1 - begin_time1) / CLOCKS_PER_SEC << endl;
    init.terminate();
    
    clock_t begin_time2 = clock();
    findShortestSuperstringParalell(Arr);
    clock_t stop_time2 = clock();
    cout << "Single time: " << double(stop_time2 - begin_time2) / CLOCKS_PER_SEC << endl;
    
    system("pause");
    return 0;
}
