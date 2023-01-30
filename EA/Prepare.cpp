#include "Prepare.hpp"

extern double X[1000], Y[1000];
extern double dimensions[1000][1000];
extern double C[1000];
extern int number_of_trucks;
extern int number_of_cities;
extern int capacity;


std::vector<std::string> split(std::string str, char delimiter) {
    std::vector<std::string> internal;
    std::stringstream ss(str); // Turn the string into a stream. 
    std::string tok;

    while (std::getline(ss, tok, delimiter)) {
        internal.push_back(tok);
    }

    return internal;
}

double calculateResultDst(std::vector<int> vector) {
    double result = 0;
    double dist;
    int startCity = vector.front();

    for (int i = 1; i < vector.size(); i++) {
        dist = dimensions[vector.at(i)][vector.at(i - 1)];
        result += dist;
    }
    result += dimensions[vector.at(vector.size() - 1)][startCity];

    return result;
}

double calculateResultDstVect(std::vector<int>* vector) {
    double result = 0;
    double dist;
    int startCity = vector -> front();

    for (int i = 1; i < vector -> size(); i++) {
        dist = dimensions[vector -> at(i)][vector -> at(i - 1)];
        result += dist;
    }
    result += dimensions[vector -> at(vector -> size() - 1)][startCity];

    return result;
}

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
		for (int i = 0; i < number_of_cities; i++)
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

std::vector<int> random(int N) {
    std::vector<int> randomPath;
    int arr[1000];
    bool duplicate = false;
    int i = 0;

    while (i < N) {
        int randomVal = std::rand() % N + 1;
        duplicate = false;

        for (int j = 0; j < i; j++) {
            if (randomVal == arr[j]) {
                duplicate = true;
            }
        }

        if (!duplicate) {
            arr[i] = randomVal;
            i++;
        }
    }

    for (int i = 1; i < N;i++) {
        randomPath.push_back(arr[i - 1]);

        if (i == N - 1) {
            randomPath.push_back(arr[i]);
        }
    }

    return randomPath;
}
