// Find the pseudo-critical threshold Tc in thermoadaptive Ising models on d-dimensional lattices
// Sets the system below Tc and gathers clusters as consecutive variations in the total magnetization
// during the metastable relaxation to the paramagnetic phase

#include <iostream>
#include <vector>
#include <cmath>
#include <random>

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
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>

// ...

// Modified function to simulate the Ising model and save spin count differences
double simulateIsingModel(const Graph& graph, int numSpins, double temperature, unsigned int seed) {

    std::mt19937 gen(seed); // Use the specified seed
    std::uniform_real_distribution<> dis(0.0, 1.0);
    std::vector<int> spins(numSpins, 1); // Initialize all spins to +1

    double prevmag = 1;

    for (int iter = 0; iter < 100; iter++) {

        double beta = prevmag / temperature;
        //double beta = prevmag * prevmag / temperature;

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

        prevmag = calculateMagnetization(spins);
    }

    double magnetization2 = calculateMagnetization(spins);
    return magnetization2;
}


double simulateIsingModel_evoclusters(const Graph& graph, int L, int numSpins, double temperature, unsigned int seed, double delta) {

    std::mt19937 gen(seed); // Use the specified seed
    std::uniform_real_distribution<> dis(0.0, 1.0);
    std::vector<int> spins(numSpins, 1); // Initialize all spins to +1

    // Generate file with delta val for unique file saving
    std::string fileName = "cascades7D/spin_count_differences_L" + std::to_string(L) + "_delta" + std::to_string(delta) + ".txt";
    std::ofstream outFile(fileName);
    if (!outFile) {
        std::cerr << "Failed to create/open the file " << fileName << std::endl;
        return -1; // Return an error value
    }

    double prevmag = 1.0; // initial magnetizaiton
    int difference = 0; // initial difference in spin count
    int prevSpinCount = numSpins; // initial number of spins up

    for (int iter = 0; iter < 1000; iter++) {
        
        double beta = prevmag / temperature;
        //double beta = prevmag * prevmag / temperature;

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

        int currentSpinCount = std::count(spins.begin(), spins.end(), 1);
        difference = currentSpinCount - prevSpinCount;
        prevmag = calculateMagnetization(spins);
        prevSpinCount = currentSpinCount;

        //double currentmag = calculateMagnetization(spins);
        //difference = numSpins * (currentmag - prevmag); 
        // prevmag = currentmag;

        outFile << iter << " " << difference << ' ' << prevmag << std::endl;

        if (prevmag < 0.1) {
            break;
        }

    }

    outFile.close();

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


/*
int main() {

    int L = 200; // Linear size of the lattice
    int D = 2;  // Number of dimensions

    Graph graph;
    graph.createGraphSLFullnD(L, D);

    // Generate a random seed for this run
    std::random_device rd;
    unsigned int seed = rd();

    // Saves clusters' evolution to an external .txt file
    double temperature = 1.83; 
    double averageMagnetization = simulateIsingModel_evoclusters(graph, L, pow(L, D), temperature, seed);
    std::cout << "Average Magnetization: " << averageMagnetization << std::endl;

    return 0;
}
*/


int main() {

    int D = 7;  // Number of dimensions
    int numLvals = 5;
    int Lvals[] = {2, 3, 4, 5};

    for (int i = 0; i < numLvals; ++i) {
        int L = Lvals[i];

        Graph graph;
        graph.createGraphSLFullnD(L, D);

        double minTemperature = 0.0;
        double maxTemperature = 10.0;
        double tolerance = 0.001;
        double left = minTemperature;
        double right = maxTemperature;

        // Generate a random seed for this run
        std::random_device rd;
        unsigned int seed = rd();

        while (right - left > tolerance) {
            double mid = (left + right) / 2.0;
            double averageMagnetization = simulateIsingModel(graph, pow(L, D), mid, seed);
            
            if (averageMagnetization < 0.1) {
                right = mid;
            } else {
                left = mid;
            }
        }
        
        // Print the criticla temperature found via bisections
        double criticalTemperature = (left + right) / 2.0;
        std::cout << "Critical Temperature: " << criticalTemperature << std::endl;

        // Call the simulateIsingModel function with the critical temperature + delta
        // Saves clusters' evolution to an external .txt file

        double end = 0.00001;
        double factor = 10;
        double delta = 0.1;

        while (delta >= end) {

            double deltatemperature = criticalTemperature + delta; 
            double averageMagnetization = simulateIsingModel_evoclusters(graph, L, pow(L, D), deltatemperature, seed, delta);
            std::cout << "Average Magnetization: " << averageMagnetization << std::endl;

            delta *= 1 / factor; // Decrease the loop counter
        }

    }

    return 0;
}



// g++ -o isingmodelD_evoclusters isingmodelD_evoclusters.cpp
// ./isingmodelD_evoclusters 
// First create an appropriate folder