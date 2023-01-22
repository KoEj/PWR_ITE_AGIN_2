#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <numeric>

std::vector<std::string> split(std::string, char);
void Loader(std::string filename, bool loaderPrint);
void createDistanceMatrix(int N, bool printMatrix);
std::vector<std::vector<int>> random(int trucksNumber, int magasinCapacity, int citiesNumber);
std::vector<std::vector<int>> greedy(int trucksNumber, int magasinCapacity, int citiesNumber);
std::vector<double> calculateResultDstVect(int startCity, std::vector<std::vector<int>> vector);
double calculateResult(std::vector<double> vectorResult);
double calculateResult(std::string name, std::vector<double> vectorResult);
int antColonyOptimalization(int trucksNumber, int magasinCapacity, int citiesNumber);
void updatePheromoneMatrix(int citiesNumber);
double calculateProbability(int current_node, int next_node);
void calculateResultDstVect(int, std::vector<int>);
void getAntResult();

double cost(std::vector<std::vector<int>>);
std::vector<std::vector<int>> initializeRandomSolution(int, int, int);
bool is_valid_route(std::vector<std::vector<int>>);
std::vector < std::vector<int>> generateNeighbor(std::vector<std::vector<int>>);
std::vector<std::vector<int>> generateNeighbors(std::vector<std::vector<int>>);
void simulated_annealing();
