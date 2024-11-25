#include <iostream>
#include <vector>
#include <random>
#include <fstream>

const int tmax = 100; // Maximum simulation time
const int k = 1; // controlling the global revival coupling
int D = 3; // Dimensionality of the lattice

// Function to count the number of active neighbors in a d-dimensional lattice
int countActiveNeighbors(const std::vector<int>& lattice, const std::vector<int>& coordinates, int D, int L) {

    int count = 0;
    std::vector<std::vector<int>> offsets;

    // Generate offsets for neighbors in each dimension
    for (int i = 0; i < D; i++) {
        std::vector<int> offset(D, 0);
        offset[i] = 1;
        offsets.push_back(offset);
    }

    // Check each neighbor
    for (const auto& offset : offsets) {
        std::vector<int> neighborCoordinates(D);
        for (int i = 0; i < D; i++) {
            neighborCoordinates[i] = (coordinates[i] + offset[i] + L) % L;
        }

        // Convert d-dimensional coordinates to a 1D index
        int neighborIndex = 0;
        for (int i = 0; i < D; i++) {
            neighborIndex += neighborCoordinates[i] * pow(L, i);
        }

        // Count the active neighbors
        if (lattice[neighborIndex] == 0) {
            count++;
        }
    }

    return count;
}



int countActiveNodes(int L, const std::vector<int>& lattice) {
    int numActiveNodes = 0;
    for (int i = 0; i < lattice.size(); i++) {
        if (lattice[i] == 1) {
            numActiveNodes++;
        }
    }
    return numActiveNodes;
}


double simulate_SIS(double p_loc, int L, int D, unsigned int seed) {

    // Create a vector to represent the lattice
    std::vector<int> lattice(pow(L, D), 1);

    std::mt19937 gen(seed);
    std::uniform_real_distribution<> dis(0.0, 1.0);

    int numActiveNodes = countActiveNodes(L, lattice);

    // Simulation loop
    for (int t = 1; t <= tmax; t++) {

        // Loop through each lattice site
        for (int i = 0; i < pow(L, D); i++) {

            std::vector<int> coordinates(D);
            int temp = i;

            // Get the coordinates
            for (int j = 0; j < D; j++) {

                coordinates[j] = temp % L;
                temp = int(temp / L);
            }

            // Infected node -- attempting recovery
            if (lattice[i] == 0) {

                double fractionA = static_cast<double>(numActiveNodes) / pow(L, D);
                double prob = std::pow(fractionA, k);

                if (dis(gen) < prob) {

                    lattice[i] = 1; // It is locally infected
                    numActiveNodes += 1;
                }

            } else { // Susceptible node -- attempting infection

                int numActiveNeighbors = countActiveNeighbors(lattice, coordinates, D, L);

                if (dis(gen) < static_cast<double>(numActiveNeighbors) / pow(2, D)) {
                    lattice[i] = 0; // It is locally infected
                    numActiveNodes += -1;
                    
                } else if (dis(gen) < p_loc) {

                    if (lattice[i] == 1 ){
                        lattice[i] = 0; // It is locally infected
                        numActiveNodes += -1;
                    }
                }
            }
        }
    }

    // Output or analyze the results as needed
    int numNodes = countActiveNodes(L, lattice);
    double fraction = static_cast<double>(numNodes) / pow(L, D);

    return fraction; 

}




int main() {

    int numLvals = 1; // Number of lattice sizes in L_vals
    int L_vals[] = {25};

    for (int i = 0; i < numLvals; ++i) {

        int L = L_vals[i];

        // Assign a random seed
        std::random_device rd;
        unsigned int seed = rd();

        // Range of values to search the threshold
        double dp_steps = 100;
        double p_inf = 0;
        double p_upper = 1.0;
        double dp = (p_upper - p_inf) / dp_steps;

        // Open a file for writing the simulation data
        std::string fileName = "SIS_glob/othersimD" +std::to_string(D)+ "_K" + std::to_string(k) + "_L" + std::to_string(L) + "_rev1.txt";
        std::ofstream outputFile(fileName);

        while (p_inf <= p_upper) {

            double fraction = simulate_SIS(p_inf, L, D, seed); 
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