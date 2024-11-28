#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <string>
#include <string.h>
#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <bitset>
#include <stdexcept>
#include <exception>
#include <map>
#include <cassert>
#include <iomanip>
#include <limits>

#include <cmath>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <random>

#include <set>
#include <filesystem>

#include <unistd.h>
#include <sys/resource.h>
#include <chrono>
#include <thread>


#include "ising_multiplex_3_q.hpp"
#include "HyperGraph.hpp"
#include "Lattice_Multiplex.hpp"
#include "Graph.hpp"

using std::pair;
using std::make_pair;
using std::vector;
using std::logic_error;
using std::map;
using std::cout;
using std::cerr;
using std::setprecision;
using namespace std;




int main(int argc, char *argv[]){



    //### Spatial hypergraph
    int L; // linear size
    double d = 2; //dimension
    int kcore = 2;
    cin>>L;
    double n = pow(L,d); //graph size
    stringstream ssL;
    ssL<<L;
    int res_num = 1000;

    for(int i=0;i<res_num;i++){
        HyperGraph HG;
        HG.setHyperGraph(n);
        HG.createHyperGraphD(L,d);
        double pc = HG.scan_p_kcore_lion(kcore);
        cout<<"pc = "<<pc<<endl;
        //HG.scan_p_kcore_lion_avalanches_single(kcore, where_to); //get av size
}


    /*
//#########################################

    //### Percolation of Interdependent networks
    int L; // linear size
    double d = 2; //dimension
    cin>>L;
    double n = pow(L,d); //graph size
    stringstream ssL;
    ssL<<L;
    int res_num = 1000;
    for(int i=0;i<res_num;i++){
        string where_to = "/Users/bnayagross/Desktop/test_res/";
        Lattice_Multiplex LM(n, d);
        double pc = LM.scan_p_giant_lion(where_to, LM.G1.gen());
        cout<<"pc = "<<pc<<endl;
    }
*/
//#########################################
/*

//### Interdependent spin networks
double L;
double q = 1;
double d = 2;
cin>>L;
stringstream ssL;
ssL<<L;
double Tmin = 0.0001, Tmax = 12;
double NOI = 100, NOF = 10*pow(L,d);
for(int i=0;i<400;i++){
    vector<double> Tc_vec;
    string where_to = "/Users/bnayagross/Desktop/test_res/";

    ising_multiplex_3_q I(L,d,q);
    I.scan_beta_local_average_lion(Tmin,Tmax, where_to, I.G1.gen(), NOI,NOF);
    //I.scan_beta_mixed_layers_rmodel(1.655,2,0.001, where_to, I.G1.gen(), 100,100*pow(L,2));
    //Tc_vec.push_back(I.scan_beta_lion_not_mixed(Tmin,Tmax, where_to, I.G1.gen(), NOI,NOF));

    }
*/

}
