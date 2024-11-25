#include <iostream>
#include <vector>
#include <random>
#include <fstream>

const int tmax = 1000;     // Maximum simulation time
const int k = 5;

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
        if (lattice[nx][ny] == 0) {
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

bool simulate_SIS(double p_loc, int L, unsigned int seed) {

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

                // Infected node -- attempting recovery
                if (lattice[i][j] == 0) {

                    double fractionA = static_cast<double>(numActiveNodes) / (L * L);
                    double prob = std::pow(fractionA, k);

                    if (dis(gen) < prob){
                        lattice[i][j] = 1; // It is locally infected 
                        numActiveNodes += 1; 
                    }
                } 

                // Susceptible node -- attempting infection
                else if (lattice[i][j] == 1) {

                    if (dis(gen) < p_loc) {
                        lattice[i][j] = 0; // It is locally infected 
                        numActiveNodes += -1;
                    }
                }
            }
        }

        double fraction = static_cast<double>(numActiveNodes) / (L * L);
        if (fraction  < 1e-1) {
            return true; // Break the simulation loop if the condition is met
        }
    }

    return false; 

}


int main() {
    
    int numLvals = 32; 
    int L_vals[] = {10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 125, 150, 175, 200, 225, 250, 275, 300, 350, 400, 450, 500, 750, 1000};

    for (int i = 0; i <= numLvals; ++i) {

        int L = L_vals[i];

        // Open a file for writing the simulation data
        std::string fileName = "SIS_glob_thresholds/othersim_th_K" + std::to_string(k) + "_L" + std::to_string(L) + ".txt";
        std::ofstream outputFile(fileName);
        
        for (int n = 0; n < 100; n++){

            // Generate a random seed for this run
            std::random_device rd;
            unsigned int seed = rd();

            // Range of values to search the threshold
            double p_inf = 0;
            double p_upper = 1;
            
            // Bisection method to find the threshold
            while (p_upper - p_inf > 1e-7) { // Precision threshold for bisection
                double p_mid = (p_inf + p_upper) / 2.0;

                if (simulate_SIS(p_mid, L, seed)) {
                    // Condition met, reduce upper bound
                    p_upper = p_mid;
                } else {
                    // Condition not met, increase lower bound
                    p_inf = p_mid;
                }
            }

            // Printing the critical threshold
            double critical_p = (p_inf + p_upper) / 2.0;
            std::cout << "Critical threshold of p: " << critical_p << std::endl;

            // Save output to the external file
            outputFile << L << " " << critical_p << std::endl;

        }

         outputFile.close();

    }

    return 0;
}


// % g++ -std=c++11 SIS_glob_thresholds.cpp -o SIS_global_th     
// ./SIS_global