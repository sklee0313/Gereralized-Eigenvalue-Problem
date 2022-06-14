#include <vector>
#include <stack>
#include <iostream>
#pragma once
std::stack<clock_t> tictoc_stack;

void tic()
{
    tictoc_stack.push(clock());
}

void toc()
{
    std::cout << ((double)(clock() - tictoc_stack.top())) / CLOCKS_PER_SEC
              << " s" << std::endl;
    tictoc_stack.pop();
}