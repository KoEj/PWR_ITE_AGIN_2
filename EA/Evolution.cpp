#include "Evolution.hpp"

int populationSize = 150;
int generations = 5000;
int tournamentSize = 15;
double mutationProbability = 0.1;
double crossProbability = 0.8;

std::vector<double> evaluation;

std::vector<int>* mutationSwap(std::vector<int>* vector) {
    // zmienic na wieksza ilosc genow
    int firstVal = std::rand() % vector -> size();
    int secondVal;
    do {
        secondVal = std::rand() % vector -> size();
    } while (secondVal == firstVal);

    std::swap(vector -> at(firstVal), vector -> at(secondVal));
    return vector;
}

std::vector<int>* mutationInversion(std::vector<int>* vector) {
    int firstVal = std::rand() % vector -> size() + 1;
    int secondVal;
    do {
        secondVal = std::rand() % vector -> size() + 1;
    } while (secondVal == firstVal);

    // first value has to be lower than second one
    if (firstVal > secondVal) {
        std::swap(firstVal, secondVal);
    }

    std::vector<int> vectorTemp;;
    for (int i = firstVal - 1; i <= secondVal - 1; i++) {
        vectorTemp.push_back(vector -> at(i));
    }

    for (int i = firstVal - 1; i <= secondVal - 1;i++) {
        vector -> at(i) = vectorTemp.back();
        vectorTemp.pop_back();
    }

    return vector;
}

std::vector<int> orderedCrossover(std::vector<int>* firstParent, std::vector<int>* secondParent) {
    std::vector<int> child;
    int firstVal = std::rand() % firstParent -> size() + 1;
    int secondVal;
    do {
        secondVal = std::rand() % firstParent -> size() + 1;
    } while (secondVal == firstVal);

    if (firstVal > secondVal) {
        std::swap(firstVal, secondVal);
    }

    /*std::cout << firstVal << " " << secondVal;

    std::cout << std::endl << "FirstParent: ";
    for (auto x : firstParent) {
        std::cout << x << " ";
    }

    std::cout << std::endl <<"SecondParent: ";
    for (auto x : secondParent) {
        std::cout << x << " ";
    }
    std::cout << std::endl;*/

    for (int i = 0;i < firstParent -> size();i++) {
        child.push_back(0);
    }

    for (int i = firstVal - 1; i <= secondVal - 1; i++) {
        child.at(i) = firstParent -> at(i);
    }

    //std::cout << "child prepare: ";
    //for (auto x : child) {
    //    std::cout << x << " ";
    //}
    //std::cout << std::endl;

    int parent_index = 0;
    for (int i = 0; i < firstParent -> size(); i++) {
        if (child.at(i) == 0) {
            for (int j = 0; j < firstParent -> size(); j++) {
                if (secondParent -> at(parent_index) == child.at(j)) {
                    parent_index++;
                    j = -1;
                }
            }
            child.at(i) = secondParent -> at(parent_index);
        }
    }

    //std::cout << "child: ";
    //for (auto x : child) {
    //    std::cout << x << " ";
    //}
    //std::cout << std::endl;
    return child;
}

std::vector<int> pmxCrossover(std::vector<int>* firstParent, std::vector<int>* secondParent) {
    std::vector<int> child = *firstParent;
    bool duplicate;

    do {
        duplicate = false;
        int firstVal = std::rand() % firstParent -> size();
        int secondVal;
        do {
            secondVal = std::rand() % firstParent -> size();
        } while (secondVal == firstVal);

        if (firstVal > secondVal) {
            std::swap(firstVal, secondVal);
        }

        /*
        std::cout << firstVal + 1 << " " << secondVal + 1;

        std::cout << std::endl << "FirstParent:  ";
        for (auto x : firstParent) {
            std::cout << x << " ";
        }

        std::cout << std::endl << "SecondParent: ";
        for (auto x : secondParent) {
            std::cout << x << " ";
        }
        std::cout << std::endl;
        */

        for (int j = firstVal; j <= secondVal; j++) {
            for (int i = 0; i < firstParent -> size(); i++) {
                if (child.at(i) == secondParent -> at(j)) {
                    child.at(i) = firstParent -> at(j);
                }
            }
        }

        for (int i = firstVal; i <= secondVal; i++) {
            child.at(i) = secondParent -> at(i);
        }

        std::vector<int> tempVector = child;
        std::sort(tempVector.begin(), tempVector.end());

        for (int i = 1; i < tempVector.size();i++) {
            if (tempVector.at(i - 1) == tempVector.at(i)) {
                duplicate = true;
                break;
            }
        }
        //for (auto x : child) {
        //    std::cout << x << " ";
        //}
        //std::cout << std::endl;
    } while (duplicate == true);

    //for (auto x : child) {
    //    std::cout << x << " ";
    //}
    //std::cout << std::endl;

    return child;
}

int tournament(std::vector<std::vector<int>>* population) {
    std::set<int> numbers;
    double bestDistance = INT_MAX;
    int bestDistanceIndex = 0;

    do {
        numbers.insert(rand() % populationSize);
    } while (numbers.size() < tournamentSize);

    for (auto& number : numbers) {
        double currDistance = calculateResultDst(population -> at(number));
        if (bestDistance > currDistance) {
            bestDistance = currDistance;
            bestDistanceIndex = number;
        }
    }
    //std::cout << std::endl << bestDistanceIndex << " " << bestDistance;
    return bestDistanceIndex;
}

int roulette(std::vector<std::vector<int>>* population) {
    std::vector<double> evalua = evaluate(population);
    double distanceSum = 0;
    double probability = 0;
    double randomProbability = 0;
    double oldIndex = 0;
    double bestDistanceIndex = 0;

    for (auto dist : evalua) {
        distanceSum += dist;
    }

    for (int i = 0; i < populationSize; i++) {
        probability += (1 - (evalua.at(i) / distanceSum)) / populationSize;
    }

    randomProbability = probability * std::rand();
    probability = 0;
    std::cout << randomProbability;

    for (int i = 0; i < populationSize; i++) {
        std::cout << "evaluation.at(i) " << evalua.at(i);
        probability += (1 - (evalua.at(i) / distanceSum));
        std::cout << ", probability " << probability;

        std::cout << ", old Index" << oldIndex;
        if (probability > randomProbability && oldIndex <= randomProbability) {
            bestDistanceIndex = oldIndex;
            std::cout << ", bestIndex" << bestDistanceIndex;

            break;
        }
        std::cout << std::endl;
        oldIndex = probability;
    }

    std::cout << "***" << bestDistanceIndex;
    return bestDistanceIndex;
}

std::vector<std::vector<int>> initialise(int N) {
    std::vector<std::vector<int>> population;

    //for (int i = 0; i < populationSize; i++) {
    //    population.push_back(random(N));
    //}
    for (int i = 0; i < populationSize; i++) {
        population.push_back(random(N));
    }

    return population;
}

std::vector<double> evaluate(std::vector<std::vector<int>>* population) {
    for (int i = 0; i < population->size(); i++) {
        evaluation.push_back(calculateResultDstVect(&(population -> at(i))));
    }

    return evaluation;
}

std::vector<int> evolution(int N, bool toCSV, int csvIteration, int iteration) {
    std::ofstream csvFile;
    if (toCSV) {
       // std::cout << "toCSV init" << std::endl;
        csvFile.open("fl417_4_3.csv", std::ios_base::app);
        csvFile << "iteration" + std::to_string(iteration) + "\n";
        csvFile << "\npopulationSize: " + std::to_string(populationSize);
        csvFile << "\ngenerations: " + std::to_string(generations);
        csvFile << "\ntournamentSize: " + std::to_string(tournamentSize);
        csvFile << "\nmutationProbability: " + std::to_string(mutationProbability);
        csvFile << "\ncrossProbability: " + std::to_string(crossProbability);
        csvFile << "\ngeneration,best,worst,avg\n";
    }

    int currentGeneration = 0;
    std::vector<std::vector<int>> population = initialise(N);
    evaluate(&population);

    while (generations > currentGeneration) {
        std::vector<std::vector<int>> newPopulation = population;
        std::vector<int> firstParent;
        std::vector<int> secondParent;
        std::vector<int> child;

        for (int i = 0; i < populationSize; i++) {
            double crossRandom = ((std::rand() % 100) + 1) / 100.0;
            double mutationRandom = ((std::rand() % 100) + 1) / 100.0;

            int firstIndex = tournament(&population);
            firstParent = population.at(firstIndex);
            child = firstParent;

            if (crossRandom <= crossProbability) {
                int secondIndex = tournament(&population);
                secondParent = population.at(secondIndex);
                child = orderedCrossover(&firstParent, &secondParent);
            }

            if (mutationRandom <= mutationProbability) {
                mutationInversion(&child);
                //mutationSwap(&child);
            }

            newPopulation.at(i) = child;

            child.clear();
            firstParent.clear();
            secondParent.clear();
        }

        population = newPopulation;
        newPopulation.clear();

        evaluate(&population);
        currentGeneration++;

        if (toCSV && currentGeneration % csvIteration == 0) {
            int bestIndex = 0;
            int worstIndex = 0;
            double maxValue = DBL_MIN;
            double minValue = DBL_MAX;
            double sum = 0;
            for (int i = 0; i < population.size(); i++) {
                //best index
                double curValue = calculateResultDst(population.at(i));
                sum += curValue;
                if (curValue < minValue) {
                    bestIndex = i;
                    minValue = curValue;
                }
                if (curValue > maxValue) {
                    worstIndex = i;
                    maxValue = curValue;
                }

            }
            //csvFile << "\ngeneration,best,worst,avg";
            csvFile << std::to_string(currentGeneration) + "," +
                std::to_string(minValue) + "," + 
                std::to_string(maxValue) + "," + 
                std::to_string(sum/population.size()) + "\n";

            //std::cout << "toCSV saved" << std::endl;
        }
    }
    
    int bestIndex = 0;
    double minValue = DBL_MAX;
    for (int i = 0; i < population.size(); i++) {
        //best index
        double curValue = calculateResultDst(population.at(i));
        if (curValue < minValue) {
            bestIndex = i;
            minValue = curValue;
        }
    }

    if (toCSV) {
        csvFile.close();
        //std::cout << "************ toCSV closed" << std::endl;

    }

    return population.at(bestIndex);
}