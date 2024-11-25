// Thermoadaptive Ising model on d-dimensional lattices
// Creates global magnetization profiles and identifies its critical thresholds through bisections. 

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>

class Graph {
    std::vector<std::vector<int > > adjacencyList;

public:
    void addEdge(int src, int dest) {
        adjacencyList[src].push_back(dest);
        adjacencyList[dest].push_back(src);
    }

    void createGraphSLFullnD(int L, int D) {
        std::vector<int> coordinates;
        coordinates.resize(D);

        int totalNodes = pow(L, D);
        adjacencyList.resize(totalNodes);

        for (int i = 0; i < totalNodes; i++) {
            int temp = i;
            int target = 0;

            // get the coordinates
            for (int j = 0; j < D; j++) {
                coordinates[j] = temp % L;
                temp = int(temp / L);
            }

            // add edges to d neighbors
            for (int j = 0; j < D; j++) {
                // periodic
                if (coordinates[j] == L - 1)
                    addEdge(i, i - (L - 1) * pow(L, j));
                else
                    addEdge(i, i + pow(L, j));
            }
        }
    }

    const std::vector<std::vector<int > > & getAdjacencyList() const {
        return adjacencyList;
    }
};

// Function to calculate the magnetization of the Ising model
double calculateMagnetization(const std::vector<int>& spins) {
    int sum = 0;
    for (size_t i = 0; i < spins.size(); i++) {
        sum += spins[i];
    }
    return static_cast<double>(sum) / spins.size();
}

// Function to simulate the Ising model
double simulateIsingModel(const Graph& graph, int numSpins, double temperature, unsigned int seed) {

    std::mt19937 gen(seed); // Use the specified seed
    std::uniform_real_distribution<> dis(0.0, 1.0);
    std::vector<int> spins(numSpins, 1); // Initialize all spins to +1

    double magnetization = 1; 

    for (int iter = 0; iter < 100; iter ++) {
        
        double beta = magnetization / temperature;

        for (int step = 0; step < 1000; step++) {

            for (size_t i = 0; i < numSpins; i++) {
                // Select a random spin
                int spinIndex = dis(gen) * numSpins;

                // Calculate the local field
                double localField = 0.0;
                for (size_t j = 0; j < graph.getAdjacencyList()[spinIndex].size(); j++) {
                    int neighbor = graph.getAdjacencyList()[spinIndex][j];
                    localField += spins[neighbor];
                }

                // Calculate the energy difference
                double deltaE = 2 * spins[spinIndex] * localField;

                // Decide whether to flip the spin
                if (deltaE <= 0 || dis(gen) < exp(-beta * deltaE)) {
                    spins[spinIndex] *= -1;
                }
            }
            
        }

        magnetization = calculateMagnetization(spins);
        if (magnetization < 0.05) {
            break;
        }

    }
    double magnetization2 = calculateMagnetization(spins);

    return magnetization2;
}

/*
int main() {
    int L = 25; // Linear size of the lattice
    int D = 2;  // Number of dimensions
    Graph graph;
    graph.createGraphSLFullnD(L, D);

    double minTemperature = 1.5;
    double maxTemperature = 3.0;
    double temperatureStep = 0.1;

    for (double temperature = minTemperature; temperature <= maxTemperature; temperature += temperatureStep) {
        double averageMagnetization = simulateIsingModel(graph, pow(L, D), temperature);
        std::cout << "Temperature: " << temperature << ", Magnetization: " << averageMagnetization << std::endl;
    }

    return 0;
}
*/


// ...

int main() {

    int L = 100; // Linear size of the lattice
    int D = 2;  // Number of dimensions

    Graph graph;
    graph.createGraphSLFullnD(L, D);

    double minTemperature = 0.0;
    double maxTemperature = 4.0;
    double tolerance = 0.001;
    int numRuns = 10; // Number of runs
    std::vector<double> thresholds; // Vector to store the thresholds

    std::ofstream outFile("thresholds_ising_D" + std::to_string(D) + "_L" + std::to_string(L) + ".txt", std::ios_base::app);
    if (!outFile) {
        std::cerr << "Failed to create/open the file." << std::endl;
        return 1; // Return an error value
    }

    for (int run = 0; run < numRuns; ++run) {
        double left = minTemperature;
        double right = maxTemperature;

        // Generate a random seed for this run
        std::random_device rd;
        unsigned int seed = rd();

        while (right - left > tolerance) {
            double mid = (left + right) / 2.0;
            double averageMagnetization = simulateIsingModel(graph, pow(L, D), mid, seed);
        
            if (averageMagnetization < 0.05) {
                right = mid;
            } else {
                left = mid;
            }
        }

        double criticalTemperature = (left + right) / 2.0;
        std::cout << "Critical Temperature (Run " << run+1 << "): " << criticalTemperature << std::endl;
        thresholds.push_back(criticalTemperature); // Store the threshold in the vector

        outFile << thresholds[run] << std::endl;
    }

    outFile.close();

    return 0;
}

