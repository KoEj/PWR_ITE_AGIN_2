#pragma once
#include <vector>
#include <set>
#include <algorithm>
#include <fstream>
#include "Prepare.hpp"

std::vector<int>* mutationSwap(std::vector<int>*);
std::vector<int>* mutationInversion(std::vector<int>*);
std::vector<int> orderedCrossover(std::vector<int>*, std::vector<int>*);
std::vector<int> pmxCrossover(std::vector<int>*, std::vector<int>*);
std::vector<double> evaluate(std::vector<std::vector<int>>*);
std::vector<std::vector<int>> initialise(int);
int tournament(std::vector<std::vector<int>>*);
int roulette(std::vector<std::vector<int>>*);
std::vector<int> evolution(int, bool, int, int);
