#include "Prepare.h"

extern double X[1000], Y[1000];
extern double dimensions[1000][1000];
extern double C[1000];
extern int number_of_trucks;
extern int number_of_cities;
extern int capacity;

//ant 
const int numberOfAnts = 100;
const int iterations = 100;
const double startingPheromoneValue = 0.1;
const double evaporationPheromone = 0.2;
const double Q = 50.0;
const double alpha = 1.0;
const double beta = 5.0;
std::string csvFileName = "test.csv";
bool toCSV = true;

double pheromone[1000][1000];
std::vector<std::vector<int>> ants(numberOfAnts, std::vector<int>(number_of_cities));
std::vector<double> ant_distance(numberOfAnts);
std::vector<int> ant_load(numberOfAnts);
std::vector<std::vector<int>> trucks;

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

//ant methods
double calculateProbability(int current_node, int next_node) {
    double pheromoneValue = pheromone[current_node][next_node];
    double distanceValue = dimensions[current_node][next_node];
    return pow(pheromoneValue, alpha) * pow(1.0 / distanceValue, beta);
}

void updatePheromoneMatrix(int citiesNumber) {
    for (int i = 0; i < numberOfAnts; i++) {
        double delta_pheromone = 1.0 / ant_distance[i];
        for (int j = 0; j < citiesNumber - 1; j++) {
            //std::cout << ants[i][j] << " " << ants[i][j + 1] << std::endl;
            pheromone[ants[i][j]][ants[i][j + 1]] += delta_pheromone;
        }
    }

    for (int i = 0; i < citiesNumber; i++) {
        for (int j = 0; j < citiesNumber; j++) {
            pheromone[i][j] *= (1 - evaporationPheromone);
        }
    }
}

void getAntResult() {
    int best_ant = 0;
    double best_distance = std::numeric_limits<double>::max();

    for (int i = 0; i < numberOfAnts; i++) {
        if (ant_distance[i] < best_distance) {
            best_distance = ant_distance[i];
            best_ant = i;
        }
    }
    trucks.clear();
    int start_city = 0;
    calculateResultDstVect(start_city, ants[best_ant]);
    std::cout << "Ant Colony Optimalization: " << best_distance << std::endl;
}

void calculateResultDstVect(int startCity, std::vector<int> route) {
    double distance = 0;
    double result = 0;
    int truckIndex = 0;
    trucks.push_back({ route[0] });

    for (int i = 1; i < route.size(); i++) {
        distance = dimensions[route.at(i)][route.at(i - 1)];
        result += distance;
        trucks[truckIndex].push_back(route[i]);
        if (result > capacity) {
            truckIndex++;
            trucks.push_back({});
            result = 0;
        }
    }
    if (trucks[truckIndex].size() != 0) {
        trucks[truckIndex].push_back(startCity);
    }

}

int antColonyOptimalization(int trucksNumber, int magasinCapacity, int citiesNumber) {
    std::ofstream csvFile;
    if (toCSV) {
        // std::cout << "toCSV init" << std::endl;
        csvFile.open(csvFileName, std::ios_base::app);
        csvFile << "numberOfAnts" + std::to_string(numberOfAnts) + "\n";
        csvFile << "\iterations: " + std::to_string(iterations);
        csvFile << "\startingPheromoneValue: " + std::to_string(startingPheromoneValue);
        csvFile << "\evaporationPheromone: " + std::to_string(evaporationPheromone);
        csvFile << "\Q: " + std::to_string(Q);
        csvFile << "\alpha: " + std::to_string(alpha);
        csvFile << "\beta: " + std::to_string(beta);
        csvFile << "\niteration,best,worst,avg\n";
    }


    std::vector<std::vector<std::vector<int>>> trucks;

    for (int i = 0; i < citiesNumber; i++) {
        for (int j = 0; j < citiesNumber; j++) {
            pheromone[i][j] = startingPheromoneValue;
        }
    }

    for (int it = 1; it <= iterations; it++) {
        for (int i = 0; i < numberOfAnts; i++) {
            int currentCity = 0;
            ants.at(i).push_back(currentCity);
            ant_load.at(i) = C[currentCity];

            for (int j = 1; j < citiesNumber; j++) {
                std::vector<double> probabilities;
                double probability_sum = 0.0;

                for (int k = 0; k < citiesNumber; k++) {
                    if (std::find(ants.at(i).begin(), ants.at(i).end(), k) != ants.at(i).end()) {
                        probabilities.push_back(0.0);
                        continue;
                    }

                    if (ant_load.at(i) + C[k] > magasinCapacity) {
                        probabilities.push_back(0.0);
                        continue;
                    }

                    double probability = calculateProbability(currentCity, k);
                    probability_sum += probability;
                    probabilities.push_back(probability);
                }

                std::vector<double> cumulative_probabilities;
                double cumulative_probability_sum = 0.0;

                for (int k = 0; k < citiesNumber; k++) {
                    probabilities[k] /= probability_sum;
                    cumulative_probability_sum += probabilities[k];
                    cumulative_probabilities.push_back(cumulative_probability_sum);
                }

                double random_value = double(std::rand()) / RAND_MAX;

                int nextCity = 0;
                for (int k = 0; k < citiesNumber; k++) {
                    if (random_value <= cumulative_probabilities[k]) {
                        nextCity = k;
                        break;
                    }
                }

                ants.at(i).push_back(nextCity);
                ant_load.at(i) += C[nextCity];
                ant_distance[i] += dimensions[currentCity][nextCity];
                currentCity = nextCity;
            }
        }
        updatePheromoneMatrix(citiesNumber);

        if (toCSV && it % 5 == 0) {
            int bestIndex = 0;
            int worstIndex = 0;
            double maxValue = DBL_MIN;
            double minValue = DBL_MAX;
            double sum = 0;
            
            for (int i = 0; i < numberOfAnts; i++) {
                sum += ant_distance[i];

                if (ant_distance[i] < minValue) {
                    minValue = ant_distance[i];
                    bestIndex = i;
                }
                if (ant_distance[i] > maxValue) {
                    maxValue = ant_distance[i];
                    worstIndex = i;
                }
            }

            //csvFile << "\ngeneration,best,worst,avg";
            csvFile << std::to_string(it) + "," +
                std::to_string(minValue) + "," +
                std::to_string(maxValue) + "," +
                std::to_string(sum / numberOfAnts) + "\n";
        }
    }

    return 0;
}
