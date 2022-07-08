#include <vector>
#pragma once

class Element
{ // git branch merge https://goddaehee.tistory.com/275
public:
    std::vector<int> nodesIds;
    Element(std::vector<int> nodesIds) : nodesIds(nodesIds) {}
};