#include "Prepare.h"

extern double X[1000], Y[1000];
extern double dimensions[1000][1000];
extern double C[1000];
extern int number_of_trucks;
extern int number_of_cities;
extern int capacity;
double pheromone[1000][1000];

/// Function to split input information
std::vector<std::string> split(std::string str, char delimiter) {
    std::vector<std::string> internal;
    std::stringstream ss(str); // Turn the string into a stream. 
    std::string tok;

    while (std::getline(ss, tok, delimiter)) {
        internal.push_back(tok);
    }

    return internal;
}

/// Function Loader
//  filename - path to the file
void Loader(std::string filename, bool loaderPrint) {
    // To reading from a file   
    std::ifstream file(filename);
    std::string line;

    /// Values to be saved from file
    // Most important
    std::string file_name;
    double file_dimension = 0;
    double file_capacity = 0;

    // Less important
    std::string file_comment;
    std::string file_type;
    std::string file_edge_type;

    // Other parameters
    int i = 0;
    int k = 0;
    int no_line = 1;

    /// READING
    while (getline(file, line)) {
        // Reading items before data
        if (line.rfind("NAME", 0) == 0) {
            std::string line_backup = line;
            char line_temp = line_backup.back();
            line_backup.pop_back();
            if (int(line_backup.back()) >= 48 && int(line_backup.back()) <= 57)
                number_of_trucks = int(line_backup.back() - 48) * 10 + line_temp - 48;
            else
                number_of_trucks = int(line_temp) - 48;
        }
        if (line.rfind("COMMENT", 0) == 0) {
            file_comment = line;
        }
        if (line.rfind("TYPE", 0) == 0) {
            file_type = line;
        }
        if (line.rfind("DIMENSION", 0) == 0) {
            std::vector<std::string> stringVector = split(line, ' ');
            file_dimension = stod(stringVector[2]);
            number_of_cities = file_dimension;
        }
        if (line.rfind("EDGE_WEIGHT_TYPE", 0) == 0) {
            file_edge_type = line;
        }
        if (line.rfind("CAPACITY", 0) == 0) {
            std::vector<std::string> stringVector = split(line, ' ');
            file_capacity = stod(stringVector[2]);
            capacity = file_capacity;
        }
        // Reading data
        /* Reading x,y data of the cities */
        if (line == "NODE_COORD_SECTION") {
            // cout << line << endl;
        }
        if (no_line > 7 && no_line <= (7 + file_dimension)) {
            std::vector<std::string> stringVector = split(line, ' ');

            X[i] = stod(stringVector[1]);
            Y[i++] = stod(stringVector[2]);
        }
        /* Reading demand of cities */
        if (line == "DEPOT_SECTION") {
            // cout << line << endl;
        }
        if (no_line > (7 + file_dimension + 1) && (no_line <= (7 + file_dimension + 1) + file_dimension)) {
            std::vector<std::string> stringVector = split(line, ' ');

            C[k++] = stod(stringVector[1]);
        }
        // When reached the end of file -> break
        if (line == "EOF") {
            break;
        }
        no_line++;
    }
    file.close();

    if (loaderPrint) {
        // Printing information about file
        std::cout << std::endl;
        std::cout << "Number of trucks: " << number_of_trucks << std::endl;
        std::cout << "Number of cities: " << number_of_cities << std::endl;
        std::cout << "Capacity: " << capacity << std::endl;
        std::cout << "Dimension: " << file_dimension << std::endl;

        std::cout << filename << std::endl;
        std::cout << no_line << std::endl;

        // Validate 
        std::cout << "No" << ' ' << "X" << ' ' << " Y" << ' ' << "C" << std::endl;
        for (int i = 0; i < 16; i++)
        {
            std::cout << i << ' ' << X[i] << ' ' << Y[i] << ' ' << C[i] << std::endl;
        }
    }
}

void createDistanceMatrix(int N, bool printMatrix) {
    // start from index number 1
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double distance = 0;
            if (i != j) {
                distance = std::sqrt(std::pow(X[i] - X[j], 2) + std::pow(Y[i] - Y[j], 2));
            }
            dimensions[i][j] = distance;
        }
    }

    if (printMatrix) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                std::cout << dimensions[i][j] << "\t";
            }
            std::cout << std::endl;
        }
    }
}

std::vector<std::vector<int>> random(int trucksNumber, int magasinCapacity, int citiesNumber) {
    std::vector<std::vector<int>> trucks;
    std::vector<std::vector<int>> truckRoutesId;
    std::vector<int> idVector;
    std::vector<int> helperVector;
    std::vector<int> cValues;
    int i = 0;
    int valuesAdded = 0;
    bool alreadyIn = false;

    for (int j = 0; j < citiesNumber; j++) {
        idVector.push_back(j);
    }

    while (valuesAdded != citiesNumber - 1) {
        valuesAdded = 0;
        i = 0;

        cValues.clear();
        trucks.clear();
        helperVector.clear();
        truckRoutesId.clear();
        
        for (int j = 0; j < citiesNumber; j++) {
            cValues.push_back(C[j]);
        }

        for (int j = 0; j < trucksNumber; j++) {
            trucks.push_back((std::vector<int>)0);
        }

        std::random_shuffle(idVector.begin(), idVector.end());

        while (i < trucksNumber) {
            for (int j = 0; j < citiesNumber; j++) {
                int truckVectorSum = std::accumulate(trucks.at(i).begin(), trucks.at(i).end(), 0);

                if ((cValues.at(idVector.at(j)) + truckVectorSum) <= magasinCapacity) {
                    if ((cValues.at(idVector.at(j)) != 0)) {
                        trucks.at(i).push_back(cValues.at(idVector.at(j)));
                        helperVector.push_back(j);
                        cValues.at(idVector.at(j)) = 0;
                        valuesAdded++;
                    }
                }
            }

            truckRoutesId.push_back(helperVector);
            helperVector.clear();
            i++;
        }
       /* std::cout << "valuesAdded: " << valuesAdded << std::endl;*/
    }

    //std::cout << std::endl << "Wyniki: " << std::endl;
    //for (auto a : trucks) {
    //    for (auto b : a) {
    //        std::cout << b << " ";
    //    }
    //    std::cout << std::endl;
    //}

    //std::cout << std::endl << "TruckRoutesId: " << std::endl;
    //for (auto a : truckRoutesId) {
    //    for (auto b : a) {
    //        std::cout << b << " ";
    //    }
    //    std::cout << std::endl;
    //}

    return truckRoutesId;
}

std::vector<std::vector<int>> greedy(int trucksNumber, int magasinCapacity, int citiesNumber) {
    std::vector<std::vector<double>> dimensionsCopy;
    std::vector<std::vector<int>> trucks;
    std::vector<std::vector<int>> truckRoutesId;
    std::vector<int> cValues;
    std::vector<int> path;
    std::vector<int> helperVector;
    std::vector<double> tempVector;
    int n = 0;
    int valuesAdded = 0;
    int currentCity = 0;
    int id = currentCity;
    double distance = DBL_MAX;

    for (int j = 0; j < citiesNumber; j++) {
        cValues.push_back(C[j]);
    }

    for (int i = 0; i < citiesNumber; i++) {
        for (int j = 0; j < citiesNumber; j++) {
            tempVector.push_back(dimensions[i][j]);
        }
        dimensionsCopy.push_back(tempVector);
        tempVector.clear();
    }

    //std::cout << "============================" << std::endl;
    //for (auto a : dimensionsCopy) {
    //    for (auto b : a) {
    //        std::cout << b << "\t";
    //    }
    //    std::cout << std::endl;
    //}
    //std::cout << "============================" << std::endl;

    for (int j = 0; j < trucksNumber; j++) {
        trucks.push_back((std::vector<int>)0);
    }

    while (n < trucksNumber) {
        for (int i = 0; i < citiesNumber; i++) {
            for (int j = 0; j < citiesNumber; j++) {
                if (dimensionsCopy[currentCity][j] != 0 && dimensionsCopy[currentCity][j] < distance) {
                    if (std::find(path.begin(), path.end(), j) == path.end()) {
                        distance = dimensionsCopy[currentCity][j];
                        id = j;
                    }

                }
            }

            int truckVectorSum = std::accumulate(trucks.at(n).begin(), trucks.at(n).end(), 0);
            if (id != 0 && ((cValues.at(id) + truckVectorSum) <= magasinCapacity)) {
                if ((cValues.at(id) != 0)) {
                    trucks.at(n).push_back(cValues.at(id));
                    path.push_back(id);
                    helperVector.push_back(id);
                    //std::cout << "Distance: " << distance << "  ";
                    //std::cout << "cValue: " << cValues.at(id) << "  ";
                    //std::cout << "id: " << id << std::endl;
                }
            }

            currentCity = id;
            id = 0;
            distance = DBL_MAX;

            for (auto a : path) {
                for (int i = 0; i < citiesNumber; i++) {
                    dimensionsCopy[a][i] = 0;
                    dimensionsCopy[i][a] = 0;
                }
            }
        } 

        truckRoutesId.push_back(helperVector);

        //std::cout << "==============" << n <<"==============" << std::endl;
        //for (auto a : dimensionsCopy) {
        //    for (auto b : a) {
        //        std::cout << b << "\t";
        //    }
        //    std::cout << std::endl;
        //}
        //std::cout << "============================" << std::endl;

        currentCity = 0;
        helperVector.clear();
        n++;
    }

    //std::cout << std::endl;std::cout << std::endl;
    //for (auto a : truckRoutesId) {
    //    for (auto b : a) {
    //        std::cout << b << "\t";
    //    }
    //    std::cout << std::endl;
    //}

    //std::cout << std::endl;std::cout << std::endl;
    //for (auto a : trucks) {
    //    for (auto b : a) {
    //        std::cout << b << "\t";
    //    }
    //    std::cout << std::endl;
    //}
        

    return truckRoutesId;
}

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

double calculateResult(std::vector<double> vectorResult) {
    return calculateResult("", vectorResult);
}

double calculateResult(std::string name, std::vector<double> vectorResult) {
    double vectorSum = std::accumulate(vectorResult.begin(), vectorResult.end(), 0.0);
    std::string stringName = name.empty() ? "" : name;
    std::cout << std::endl << "VectorSum " << stringName  << ": " << vectorSum << std::endl;
    return vectorSum;
}

std::vector<std::vector<int>> antColonyOptimalization(int trucksNumber, int magasinCapacity, int citiesNumber) {
    int numberOfAnts = 20;
    int iterations = 100;
    double startingPheromoneValue = 0.1;
    double evaporationPheromone = 0.1;
    double Q = 100.0;
    double alpha = 2.0;
    double beta = 5.0;

    for (int i = 0; i < citiesNumber; i++) {
        for (int j = 0; j < citiesNumber; j++) {
            pheromone[i][j] = startingPheromoneValue;
        }
    }

    for (int it = 0; it < iterations; it++) {
        std::vector<std::vector<int>> ants(numberOfAnts, std::vector<int>(citiesNumber));
        std::vector<int> antLoadVec(numberOfAnts);
        std::vector<double> antDistanceVec(numberOfAnts);

        for (int i = 0; i < numberOfAnts; i++) {
            ants.at(i).at(0) = rand() % citiesNumber;
            antLoadVec.at(i) = C[ants.at(i).at(0)];
        }


        for (int i = 0; i < numberOfAnts; i++) {
            for (int j = 1; j < citiesNumber; j++) {
                std::vector<double> probabilities(citiesNumber);

                for (int k = 0; k < citiesNumber; k++) {
                    if (k != ants.at(i).at(j - 1)) {
                        if (antLoadVec[i] + C[k] <= magasinCapacity) {
                            probabilities[k] = calculateProbability(pheromone, dimensions, ants.at(i).at(j - 1), k, alpha, beta, citiesNumber);
                        }
                        else {
                            probabilities[k] = 0;
                        }
                    }
                }

                // Normalize the probabilities
                double sum_probabilities = 0.0;
                for (int k = 0; k < citiesNumber; k++) {
                    sum_probabilities += probabilities[k];
                }
                for (int k = 0; k < citiesNumber; k++) {
                    probabilities[k] /= sum_probabilities;
                }

                // Choose the next node based on the probabilities
                double random_number = (double)rand() / RAND_MAX;
                double cumulative_probability = 0.0;
                int nextCity = 0;
                for (int k = 0; k < citiesNumber; k++) {
                    cumulative_probability += probabilities[k];
                    if (random_number <= cumulative_probability) {
                        nextCity = k;
                        break;
                    }
                }

                // Move the ant to the next node
                ants[i][j] = nextCity;
                antLoadVec[i] += C[nextCity];
                antDistanceVec[i] += dimensions[ants.at(i).at(j - 1)][nextCity];

            }
        }
        if (it == (iterations - 1)) {
            return ants;
        } else {
            updatePheromoneMatrix(pheromone, ants, antDistanceVec, numberOfAnts, alpha, beta, citiesNumber, evaporationPheromone);
        }
        // next iteration
    }
}

double calculateProbability(double pheromone[][1000], double dimension[][1000], int i, int j, int alphaParam, int betaParam, int citiesNumber) {
    double numerator = pow(pheromone[i][j], alphaParam) * pow(1.0 / dimension[i][j], betaParam);
    double denominator = 0.0;
    for (int k = 0; k < citiesNumber; k++) {
        if (k != i) {
            denominator += pow(pheromone[i][k], alphaParam) * pow(1.0 / dimension[i][k], betaParam);
        }
    }
    return numerator / denominator;
}

void updatePheromoneMatrix(double pheromone[][1000], std::vector<std::vector<int>>& ants, std::vector<double>& antDistance, int numAnts, double alpha, double beta, int citiesNumber, double evaporate) {
    for (int i = 0; i < numAnts; i++) {
        double delta_pheromone = 1.0 / antDistance[i];
        for (int j = 0; j < citiesNumber - 1; j++) {
            pheromone[ants.at(i).at(j)][ants.at(i).at(j + 1)] += delta_pheromone;
        }
    }

    for (int i = 0; i < citiesNumber; i++) {
        for (int j = 0; j < citiesNumber; j++) {
            pheromone[i][j] *= (1 - evaporate);
        }
    }
}

double calculateTotalDistance(double dimension[][1000], std::vector<int>& ant) {
    double total_distance = 0.0;
    for (int i = 0; i < ant.size() - 1; i++) {
        total_distance += dimension[ant[i]][ant[i + 1]];
    }
    return total_distance;
}

std::vector<std::vector<int>> getResult(double dimension[][1000], std::vector<std::vector<int>>& ants, int startCity, int num_ants, int citiesNumber) {
    int bestAnt = 0;
    double bestDistance = std::numeric_limits<double>::max();
    std::vector<std::vector<int>> resultVec(number_of_trucks);

    // Find the best ant
    for (int i = 0; i < num_ants; i++) {
        std::cout << i << " ";
        double currentDistance = calculateTotalDistance(dimension, ants[i]);
        if (currentDistance < bestDistance) {
            bestDistance = currentDistance;
            bestAnt = i;
        }
    }

    // Get the result vector
    int currentCity = startCity;
    for (int i = 0; i < citiesNumber; i++) {
        resultVec[0].push_back(currentCity);
        currentCity = ants[bestAnt][currentCity];
    }

    return resultVec;
}