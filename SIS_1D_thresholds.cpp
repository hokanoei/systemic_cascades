#include <iostream>
#include <vector>
#include <random>
#include <fstream>

const int tmax = 10000; // MCsweeps
const int k = 1;

// Function to count the number of active neighbors (2 nearest neighbors in 1D)
int countActiveNeighbors(const std::vector<int>& lattice, int i) {
    int count = 0;
    int N = lattice.size();

    // Define the offsets for the 2 nearest neighbors (left, right)
    std::vector<int> offsets = { -1, 1 };

    // Check each neighbor
    for (const auto& offset : offsets) {
        // Calculate the neighbor index with periodic wrapping
        int ni = (i + offset + N) % N;

        // Count the active neighbors
        if (lattice[ni] == 1) {
            count++;
        }
    }

    return count;
}


// Function to count the number of active nodes (A) in the lattice
int countActiveNodes(const std::vector<int>& lattice) {
    int numActiveNodes = 0;
    for (int i = 0; i < lattice.size(); i++) {
        if (lattice[i] == 1) {
            numActiveNodes++;
        }
    }
    return numActiveNodes;
}



bool simulate_SIS_1D(double p_loc, int L, unsigned int seed) {

    // Create a 1D lattice of active nodes (A) and initialize all nodes to 1 (active)
    std::vector<int> lattice(L, 1);

    // Create a uniform_real_distribution for random numbers
    std::mt19937 gen(seed); // Use the specified seed
    std::uniform_real_distribution<> dis(0.0, 1.0);
    int numActiveNodes = countActiveNodes(lattice);

    // Setting sufficient MCsweeps to wash out self-correlations
    int tmax = L; 

    // Simulation loop
    for (int t = 1; t <= tmax; t++) {

        // Loop through each lattice site
        for (int i = 0; i < L; i++) {

            // Infecting nodes
            if (lattice[i] == 0) {

                double fractionA = static_cast<double>(numActiveNodes) / L;
                double prob = std::pow(fractionA, k);

                int numActiveNeighbors = countActiveNeighbors(lattice, i);

                // if k=0, then only lattice;
                // if k>0, then global influence + lattice

                if (dis(gen) < prob * numActiveNeighbors / 2) {
                    lattice[i] = 1;
                    numActiveNodes += 1;
                }
            }
            // Reviving nodes
            else if (lattice[i] == 1) {
                if (dis(gen) < p_loc) {
                    lattice[i] = 0; // It is locally infected
                    numActiveNodes += -1;
                }
            }
        }

        double fraction = static_cast<double>(numActiveNodes) / L;
        if (fraction  < 1e-1) {
            return true; // Break the simulation loop if the condition is met
        }
    }

    return false; 
}


int main() {
    
    int numLvals = 0; 
    int L_vals[] = {100000};

    for (int i = 0; i <= numLvals; ++i) {

        int L = L_vals[i];

        // Open a file for writing the simulation data
        std::string fileName = "SIS_1Dglob_thresholds/othersim_th_K" + std::to_string(k) + "_L" + std::to_string(L) + ".txt";
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

                if (simulate_SIS_1D(p_mid, L, seed)) {
                    // Condition met, reduce upper bound
                    p_upper = p_mid;
                } else {
                    // Condition not met, increase lower bound
                    p_inf = p_mid;
                }
            }

            // Printing the critical threshold
            double critical_p = (p_inf + p_upper) / 2.0;
            std::cout << "Critical threshold at L = " << L << " is p = " << critical_p << std::endl;

            // Save output to the external file
            outputFile << L << " " << critical_p << std::endl;

        }

         outputFile.close();

    }

    return 0;
}



// % g++ -std=c++11 SIS_global.cpp -o SIS_global      
// ./SIS_global