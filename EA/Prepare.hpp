#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <numeric>

std::vector<std::string> split(std::string, char);
double calculateResultDst(std::vector<int>);
double calculateResultDstVect(std::vector<int>*);

void createDistanceMatrix(int, bool); 
void Loader(std::string, bool);

std::vector<int> random(int);
//std::vector<int> greedy(int, int);