#include "Solution.h"

// Global parameters
double X[1000], Y[1000];                  // values to increse
double C[1000];                          // as above
double dimensions[1000][1000];            // as above

int number_of_trucks;
int number_of_cities;
int capacity;

std::string filename = "file.dat";  // path to change (depends on where the file/files are)

int main() {
    srand(time(NULL));
    randomSolution(filename, false, false);
    greedySolution(filename, false, false);
    return 0;
}