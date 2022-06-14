// C++ program to print all combination
// of size r in an array of size n
#pragma once
#include <bits/stdc++.h>
#include <vector>
class Combination
{
private:
    // void combinationUtil(int arr[], int data[],
    //                      int start, int end,
    //                      int index, int r, std::vector<std::array<int, 2>> &result);
    /* arr[] ---> Input Array
data[] ---> Temporary array to
store current combination
start & end ---> Starting and
Ending indexes in arr[]
index ---> Current index in data[]
r ---> Size of a combination to be printed */
    void combinationUtil(int arr[], int data[],
                         int start, int end,
                         int index, int r, std::vector<std::array<int, 2>> &result)
    {
        // Current combination is ready
        // to be printed, print it
        if (index == r)
        {
            std::array<int, 2> tmp;
            for (int j = 0; j < r; j++)
            {
                // std::cout << data[j] << " ";
                tmp[j] = data[j];
            }
            result.push_back(tmp);
            std::cout << std::endl;
            return;
        }

        // replace index with all possible
        // elements. The condition "end-i+1 >= r-index"
        // makes sure that including one element
        // at index will make a combination with
        // remaining elements at remaining positions
        for (int i = start; i <= end &&
                            end - i + 1 >= r - index;
             i++)
        {
            data[index] = arr[i];
            combinationUtil(arr, data, i + 1,
                            end, index + 1, r, result);
        }
    }

public:
    // The main function that prints
    // all combinations of size r
    // in arr[] of size n. This function
    // mainly uses combinationUtil()
    void getCombination(int arr[], int n, int r, std::vector<std::array<int, 2>> &result)
    {
        // A temporary array to store
        // all combination one by one
        int data[r];

        // Print all combination using
        // temporary array 'data[]'
        combinationUtil(arr, data, 0, n - 1, 0, r, result);
    }
};