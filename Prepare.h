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
std::vector<std::vector<int>> antColonyOptimalization(int trucksNumber, int magasinCapacity, int citiesNumber);
double calculateProbability(double pheromone[][1000], double dimension[][1000], int i, int j, int alphaParam, int betaParam, int citiesNumber);
void updatePheromoneMatrix(double pheromone[][1000], std::vector<std::vector<int>>& ants, std::vector<double>& antDistance, int numAnts, double alpha, double beta, int citiesNumber, double evaporate);
double calculateTotalDistance(double dimension[][1000], std::vector<int>& ant);
std::vector<std::vector<double>> getResult(double dimension[][1000], std::vector<std::vector<int>>& ants, int startCity, int num_ants, int citiesNumber);