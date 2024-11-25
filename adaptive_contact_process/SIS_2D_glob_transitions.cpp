#include <iostream>
#include <vector>
#include <random>
#include <fstream>

const int tmax = 100;     // Maximum simulation time
const int k = 3;

// Function to count the number of active neighbors (4 nearest neighbors)
int countActiveNeighbors(const std::vector<std::vector<int > > & lattice, int x, int y) {

    int count = 0;
    int N = lattice.size();

    // Define the offsets for the 4 nearest neighbors (top, bottom, left, right)
    std::vector<std::pair<int, int > > offsets;
    offsets.push_back(std::make_pair(0, 1));
    offsets.push_back(std::make_pair(0, -1));
    offsets.push_back(std::make_pair(-1, 0));
    offsets.push_back(std::make_pair(1, 0));

    // Check each neighbor
    for (const auto& offset : offsets) {
        // Calculate the neighbor coordinates with periodic wrapping
        int nx = (x + offset.first + N) % N;
        int ny = (y + offset.second + N) % N;

        // Count the active neighbors
        if (lattice[nx][ny] == 1) {
            count++;
        }
    }

    return count;
}


// Function to count the number of active nodes (A) in the lattice
int countActiveNodes(int L, const std::vector<std::vector<int > > & lattice) {
    int numActiveNodes = 0;
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            if (lattice[i][j] == 1) {
                numActiveNodes++;
            }
        }
    }
    return numActiveNodes;
}

double simulate_SIS(double p_loc, int L, unsigned int seed) {

    // Create a 2D lattice of active nodes (A) and initialize all nodes to 1 (active)
    std::vector<std::vector<int > > lattice(L, std::vector<int>(L, 1));

    // Create a uniform_real_distribution for random numbers
    std::mt19937 gen(seed); // Use the specified seed
    std::uniform_real_distribution<> dis(0.0, 1.0);
    int numActiveNodes = countActiveNodes(L, lattice);

    // Simulation loop
    for (int t = 1; t <= tmax; t++) {
    
        // Loop through each lattice site
        for (int i = 0; i < L; i++) {

            for (int j = 0; j < L; j++) {

                // Infecting nodes
                if (lattice[i][j] == 0) {

                    double fractionA = static_cast<double>(numActiveNodes) / (L * L);
                    double prob = std::pow(fractionA, k);

                    int numActiveNeighbors = countActiveNeighbors(lattice, i, j);

                    // if k=0, then only lattice; 
                    // if k>0, then global influence + lattice
                    if (dis(gen) < prob * numActiveNeighbors / 4){
                        lattice[i][j] = 1; 
                        numActiveNodes += 1; 
                    }
                } 

                // Reviving nodes
                else if (lattice[i][j] == 1) {

                    if (dis(gen) < p_loc) {

                        lattice[i][j] = 0; // It is locally infected 
                        numActiveNodes += -1;
                    }
                }
            }
        }
    }

    int numNodes = countActiveNodes(L, lattice);
    double fraction = static_cast<double>(numNodes) / (L * L);

    return fraction; 

}


int main() {
    
    int numLvals = 0; 
    int L_vals[] = {100};

    for (int i = 0; i <= numLvals; ++i) {

        int L = L_vals[i];

        // Generate a random seed for this run
        std::random_device rd;
        unsigned int seed = rd();

        // Range of values to search the threshold
        double dp_steps = 100; 
        double p_inf = 0;
        double p_upper = 1;
        double dp = (p_upper - p_inf) / dp_steps; 

        // Open a file for writing the simulation data
        std::string fileName = "SIS_glob/SISlattice_K" + std::to_string(k) + "_L" + std::to_string(L) + ".txt";
        std::ofstream outputFile(fileName);

        while (p_inf <= p_upper) { 

            double fraction = simulate_SIS(p_inf, L, seed); 
            std::cout << "This is: " << L << ", for p=" << p_inf << "we have: S=" << fraction << std::endl;

            // Save the time step and the total number of active nodes to the file
            outputFile << p_inf << " " << fraction << std::endl;

            // Increase the loop counter
            p_inf += dp; 

        }

     // Close the file
    outputFile.close();
    
    }

    return 0;
}


// % g++ -std=c++11 SIS_glob_transitions.cpp -o SIS_transition
// ./SIS_transition