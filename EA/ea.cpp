#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <sstream> 
#include <map>
#include <algorithm>
#include <ctime>
#include <set>
#include "Evolution.hpp"
#include "Prepare.hpp"

double X[1000], Y[1000];
double dimensions[1000][1000];
double C[1000];
double capac[1000];
double visited[1000];
int number_of_trucks = 40;
int number_of_cities = 351;
int capacity = 436;

int pop_size;
int number_of_retesting = 10;

//std::string fileName = "berlin11_modified.tsp";
std::string fileName = "file6h.dat";


std::vector<double> calculateResultDstVect(int startCity, std::vector<std::vector<int>> vector) {
    double distance = 0;
    double result = 0;
    std::vector<double> distances;

    for (auto vec : vector) {
        // a - vector int
        result += dimensions[startCity][vec.at(0)];
        //std::cout << "result1: " << result << " ";
        for (int i = 1; i < vec.size(); i++) {
            distance = dimensions[vec.at(i)][vec.at(i - 1)];
            result += distance;
            //std::cout << "result2: " << result << " ";
        }

        result += dimensions[vec.at(vec.size() - 1)][startCity];
        //std::cout << "result3: " << result << " ";
        distances.push_back(result);
        result = 0;
    }

    //std::cout << std::endl << "calculateResultDstVect: " << std::endl;
    //for (auto a : distances) {
    //    std::cout << a << " ";
    //}

    return distances;
}

double cost(std::vector<std::vector<int>> routes) {
    double all_cost = 0;
    std::vector<double> truckDistances = calculateResultDstVect(0, routes);

    for (int i = 0; i < truckDistances.size(); i++)
    {
        all_cost += truckDistances[i];
    }
    return all_cost;
}

int main() {
    bool visted = false;
    std::vector<int> greedy_route;
    std::vector<std::vector<int>> greedy_routes;

    srand(time(NULL));
    double randomResultSum = 0, randomMinResult = DBL_MAX, randomMaxResult = 0;
    double greedyResultSum = 0, greedyMinResult = DBL_MAX, greedyMaxResult = 0;
    double geneticResultSum = 0, geneticMinResult = DBL_MAX, geneticMaxResult = 0;

    // For now hardcoded, can be changed later
    Loader(fileName, false);
    createDistanceMatrix(number_of_cities, false);
    
    for (int i = 1; i < number_of_retesting + 1; i++) {
        std::cout << "iteration: " << i << std::endl;

        do {
            std::vector<int> geneticPath = evolution(number_of_cities, false, 5, i);

            //std::cout << "Truck capacity: " << capacity << std::endl;

            for (int j = 0; j < number_of_trucks; j++) {
                capac[j] = capacity;
            }

            for (int j = 0; j < geneticPath.size(); j++) {
                if (C[geneticPath[j] - 1] != 0) {
                    visited[j] = 0;
                    //std::cout << geneticPath[j] << "\t[" << C[geneticPath[j] - 1] << "]\n";
                }
            }

            for (int j = 0; j < number_of_trucks; j++) {
                for (int k = 0; k < geneticPath.size(); k++) {
                    if (capac[j] >= C[geneticPath[k] - 1] && visited[k] == 0 && C[geneticPath[k] - 1] != 0) {
                        //std::cout << j << ' ' << capac[j] << ' ' << C[geneticPath[k] - 1] << "\t New cap: " << capac[j] - C[geneticPath[k] - 1] << std::endl;
                        capac[j] -= C[geneticPath[k] - 1];
                        visited[k] = 1;
                        greedy_route.push_back(geneticPath[k]);
                        //std::cout << geneticPath[k] << " ";
                    }
                    else if (C[geneticPath[k] - 1] == 0) {
                        visited[k] = 2;
                    }
                }
                //std::cout << "\n";
                greedy_routes.push_back(greedy_route);
                greedy_route.clear();
            }

            for (int j = 0; j < geneticPath.size(); j++) {
                if (visited[j] == 0) {
                    std::cout << "\nERROR, ONE OR MORE CITIES NOT VISITED!";
                    visted = false;
                }
                else
                {
                    visted = true;
                }
            }
        } while (visted == false);

        /*std::cout << std::endl << "TruckRoutesId: " << std::endl;
        for (auto a : greedy_routes) {
            for (auto b : a) {
                std::cout << b << " ";
            }
            std::cout << std::endl;
        }*/
        double calculateGreedy = cost(greedy_routes);
        
        std::cout << "\ncalculateGreedy: " << calculateGreedy << std::endl;

        //double calculatedGeneticResult = calculateResultDst(geneticPath);
        geneticResultSum += calculateGreedy;
        geneticMinResult = std::min(geneticMinResult, calculateGreedy);
        geneticMaxResult = std::max(geneticMaxResult, calculateGreedy);
        calculateGreedy = 0; 
        greedy_routes.clear();
    }

    std::cout << "\nGenetic result avg: " << geneticResultSum / number_of_retesting << std::endl;
    std::cout << "Genetic min result: " << geneticMinResult << std::endl;
    std::cout << "Genetic max result: " << geneticMaxResult << std::endl;

    std::cout;
	std::cout << "END" << std::endl;
    return 0;
}
