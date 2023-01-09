#include "Prepare.h"
#include "Solution.h"

extern double X[1000], Y[1000];
extern double dimensions[1000][1000];
extern double C[1000];
extern int number_of_trucks;
extern int number_of_cities;
extern int capacity;

void randomSolution(std::string filename, bool loaderPrint, bool distanceMatrixPrint) {
	Loader(filename, loaderPrint);
	createDistanceMatrix(number_of_cities, distanceMatrixPrint);
	std::vector<std::vector<int>> truckRoutesId = random(number_of_trucks, capacity, number_of_cities);
	std::vector<double> truckDistances = calculateResultDstVect(0, truckRoutesId);
	calculateResult(truckDistances);
}
