/*
 * ising_multiplex_3_q.hpp
 *
 *  Created on: Feb 20, 2018
 *      Author: benaya
 */

#ifndef SRC_ising_multiplex_3_q_HPP_
#define SRC_ising_multiplex_3_q_HPP_


#include "Graph.hpp"
#include <math.h>
#include <string>
#include <numeric>



class ising_multiplex_3_q{
public:

    double L,N;  // size, N=LXL
    double q;    // fraction of interdependent nodes
    double p;
    Graph G1, G2, G3; // The networks
    vector<bool> inter_connected12, inter_connected13, inter_connected23; //vector which state is the node is interdependent or not
    vector<int> spins1, spins2, spins3;  // The state of the spins in each network. [-1,1].
    vector<double> local_m1, local_m2, local_m3;

    vector<double> local_m1_thermal_avg, local_m2_thermal_avg, local_m3_thermal_avg;
    vector<double> local_m1_counter, local_m2_counter, local_m3_counter;

    double M1, M2, M3;
    double E1, E2, E3, E;

    double giant1, giant2, giant3;

    ising_multiplex_3_q(double L_, double d, double q_);
    void Create_interlinks();
    void initialize_spins_disordered();
    void initialize_spins_ordered_up();
    void initialize_spins_ordered_down();
    void initializing_local_m();
    void flip_spin(int spin, int net_idx, int step);
    void revive_ordered_up();
    void revive_ordered_down();
    void revive_disordered();

    double scan_beta_local_average_lion(double beta_i,double beta_f,string where_to, int TestNum, double NOI, double NOF);
    double scan_beta_lion_not_mixed(double beta_i,double beta_f,string where_to, int TestNum, double NOI, double NOF);

    void scan_beta_thermal(double Ti,double Tf,double dT, string where_to, int TestNum, double NOI, double NOF);
    void scan_beta_highorder(double Ti,double Tf,double dT,string where_to, int TestNum, double NOI, double NOF);

};

ising_multiplex_3_q::ising_multiplex_3_q(double L_, double d, double q_):L(L_),N(pow(L_,d)),q(q_){

    G1.setGraph(N);
    G1.createGraphSLFullnD(int(pow(N,(double)1.0/(double) d)),d);
    giant1 = G1.testConnectivity();
    for(double i=0;i<1000000000;i++) //pause
        i=i*1;

    G2.setGraph(N);
    G2.createGraphSLFullnD(int(pow(N,(double)1.0/(double)d)),d);
    giant2 = G2.testConnectivity();

    for(double i=0;i<1000000000;i++) //pause
        i=i*1;

    G3.setGraph(N);
    G3.createGraphSLFullnD(int(pow(N,(double)1.0/(double)d)),d);
    giant3 = G3.testConnectivity();


    Create_interlinks();
    initialize_spins_ordered_up();
    //initialize_spins_disordered();
    //initialize_spins_ordered_down();
    initializing_local_m();



    M1 = 0;
    M2 = 0;
    M3 = 0;
    for(int i=0;i<N;i++){
        if(G1.active[i]){
            M1+=spins1[i];
            M2+=spins2[i];
            M3+=spins3[i];
        }
    }


    cout<<"finish constructor"<<endl;

}

void ising_multiplex_3_q::Create_interlinks(){
    inter_connected12 = vector<bool>(N, 0); //create vector of size N full with zeroes
    inter_connected13 = vector<bool>(N, 0); //create vector of size N full with zeroes
    inter_connected23 = vector<bool>(N, 0); //create vector of size N full with zeroes

    double rand_num;
    for(int i=0;i<N;i++){ //for each node
        inter_connected12[i] = 1;

        rand_num = (double)G1.gen()/G1.gen.max();
        if(rand_num < q){
            inter_connected13[i] = 1;
            inter_connected23[i] = 1;
        }


    }
}

void ising_multiplex_3_q::initialize_spins_disordered(){
    double rand_num; //for each nodes intialized the spin +1 or -1 with equal probability

    for(int i=0;i<N;i++){ //network1
        rand_num = (double)G1.gen()/G1.gen.max();
        if(rand_num < 0.5){
            spins1.push_back(1);
        }
        else{
            spins1.push_back(-1);
        }
    }

    for(int i=0;i<N;i++){ //network2
        rand_num = (double)G1.gen()/G1.gen.max();
        if(rand_num < 0.5){
            spins2.push_back(1);
        }
        else{
            spins2.push_back(-1);
        }
    }

    for(int i=0;i<N;i++){ //network3
        rand_num = (double)G1.gen()/G1.gen.max();
        if(rand_num < 0.5){
            spins3.push_back(1);
        }
        else{
            spins3.push_back(-1);
        }
    }
}

void ising_multiplex_3_q::initialize_spins_ordered_up(){
    spins1 = vector<int>(N,1); //all the spins are up in both layers
    spins2 = vector<int>(N,1);
    spins3 = vector<int>(N,1);

}

void ising_multiplex_3_q::initialize_spins_ordered_down(){
    spins1 = vector<int>(N,-1); //all the spins are down in both layers
    spins2 = vector<int>(N,-1);
    spins3 = vector<int>(N,-1);

}


void ising_multiplex_3_q::initializing_local_m(){

    for(int i=0;i<N;i++){
        /*if(!G1.connected[i]){
            local_m1.push_back(0);
            local_m1_thermal_avg.push_back(0);
            local_m1_counter.push_back(0);
            continue;
        }*/

        double sum = 0, count = 0;
        for(int j=0;j<G1.v[i].size();j++){
            //if(G1.connected[G1.v[i][j]]){
            count++;
            sum+=spins1[G1.v[i][j]];
            //	}
        }

        if(count != 0)
            local_m1.push_back((double)sum/count);
        else
            local_m1.push_back(0);
        local_m1_thermal_avg.push_back(0);
        local_m1_counter.push_back(0);

    }

    for(int i=0;i<N;i++){
        /*if(!G2.connected[i]){
            local_m2.push_back(0);
            local_m2_thermal_avg.push_back(0);
            local_m2_counter.push_back(1);
            continue;
        }*/
        double sum = 0, count = 0;;
        for(int j=0;j<G2.v[i].size();j++){
            //	if(G2.connected[G2.v[i][j]]){
            sum+=spins2[G2.v[i][j]];
            count++;
            //}
        }

        if(count != 0){
            local_m2_thermal_avg.push_back((double)sum/count);
            local_m2.push_back((double)sum/count);
        }
        else{
            local_m2.push_back(0);
            local_m2_thermal_avg.push_back(0);
        }
        local_m2_counter.push_back(1);

    }

    for(int i=0;i<N;i++){
        /*	if(!G3.connected[i]){
                local_m3.push_back(0);
                local_m3_thermal_avg.push_back(0);
                local_m3_counter.push_back(1);
                continue;
            }*/
        double sum = 0, count = 0;;
        for(int j=0;j<G3.v[i].size();j++){
            //	if(G3.connected[G3.v[i][j]]){
            sum+=spins3[G3.v[i][j]];
            count++;
            //}
        }

        if(count != 0){
            local_m3_thermal_avg.push_back((double)sum/count);
            local_m3.push_back((double)sum/count);
        }
        else{
            local_m3.push_back(0);
            local_m3_thermal_avg.push_back(0);
        }
        local_m3_counter.push_back(1);

    }
}


void ising_multiplex_3_q::flip_spin(int spin, int net_idx, int step){

    //if(net_idx == 1 and G1.connected[spin]){
    if(net_idx == 1){
        spins1[spin]*=-1;

        for(int i=0;i<G1.v[spin].size();i++){
            int j = G1.v[spin][i];
            //if(!G1.connected[j])
            //	continue;

            double count = 0;
            for(int jj = 0; jj < G1.v[j].size(); jj++)
                //if(G1.connected[G1.v[j][jj]])
                count++;

            local_m1[j]+= (double)2*spins1[spin]/count;

            if(step > 10000000 and step%10000 == 0){
                local_m1_thermal_avg[j] += local_m1[j];
                local_m1_counter[j]++;
            }

            //E1 += -4*spins1[spin]*spins1[j];

        }
        M1+=2*spins1[spin];
    }
        //else if(net_idx == 2 and G2.connected[spin]){
    else if(net_idx == 2){
        spins2[spin]*=-1;

        for(int i=0;i<G2.v[spin].size();i++){
            int j = G2.v[spin][i];
            //if(!G2.connected[j])
            //	continue;

            double count = 0;
            for(int jj = 0; jj < G2.v[j].size(); jj++)
                //if(G2.connected[G2.v[j][jj]])
                count++;

            local_m2[j]+= (double)2*spins2[spin]/count;

            if(step > 10000000 and step%10000 == 0){
                local_m2_thermal_avg[j] += local_m2[j];
                local_m2_counter[j]++;
            }
            //E2 += -4*spins2[spin]*spins2[j];
        }

        M2+=2*spins2[spin];
    }
        //else if(net_idx == 3 and G3.connected[spin]){
    else if(net_idx == 3){
        spins3[spin]*=-1;

        for(int i=0;i<G3.v[spin].size();i++){
            int j = G3.v[spin][i];
            //if(!G3.connected[j])
            //	continue;

            double count = 0;
            for(int jj = 0; jj < G3.v[j].size(); jj++)
                //	if(G3.connected[G3.v[j][jj]])
                count++;

            local_m3[j]+= (double)2*spins3[spin]/count;

            if(step > 10000000 and step%10000 == 0){
                local_m3_thermal_avg[j] += local_m3[j];
                local_m3_counter[j]++;
            }
            //E2 += -4*spins2[spin]*spins2[j];
        }

        M3+=2*spins3[spin];
    }

}

double ising_multiplex_3_q::scan_beta_local_average_lion(double Tmin,double Tmax,string where_to, int TestNum, double NOI, double NOF){


    vector<double> T_vec, M_vec;
    double rand_num, avg_local_m;
    int u;
    double pi;
    double T = (Tmax + Tmin)/2;
    double epsilon = 0.2, epsilon_T = 0.001;
    while(Tmax - Tmin > epsilon_T){

        double beta = 1/T;
        for(int j=0;j<NOI;j++){
            //cout<<"NOI = "<<j<<"  ,M1="<<M1/giant1<<"  ,M2="<<M2/giant2<<",  beta="<<beta<<endl;

            //avg_local_m = calc_average_local_M(1);
            for(int i=0;i<3*NOF;i++){

                int ppi = G1.gen()%3;
                if(ppi == 0){
                    u = G1.gen()%int(N);

                    if(inter_connected12[u] and inter_connected13[u])
                        //pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*local_m2[u]*local_m3[u]*G1.v[u].size()));
                        pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*(M2/N)*(M3/N)*G1.v[u].size()));
                    else if(inter_connected12[u] and !inter_connected13[u])
                        pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*local_m2[u]*G1.v[u].size()));
                    else if(!inter_connected12[u] and inter_connected13[u])
                        pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*local_m3[u]*G1.v[u].size()));
                    else
                        pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*G1.v[u].size()));

                    rand_num = (double)G1.gen()/G1.gen.max();
                    if(rand_num < pi){
                        flip_spin(u,1,i);
                    }
                }
                else if (ppi == 1){
                    u = G2.gen()%int(N);

                    if(inter_connected12[u] and inter_connected23[u])
                        //pi = 1/(1+exp(+2*beta*spins2[u]*local_m1[u]*local_m2[u]*local_m3[u]*G2.v[u].size()));
                        pi = 1/(1+exp(+2*beta*spins2[u]*(M1/N)*local_m2[u]*(M3/N)*G2.v[u].size()));
                    else if(inter_connected12[u] and !inter_connected23[u])
                        pi = 1/(1+exp(+2*beta*spins2[u]*local_m1[u]*local_m2[u]*G2.v[u].size()));
                    else if(!inter_connected12[u] and inter_connected23[u])
                        pi = 1/(1+exp(+2*beta*spins2[u]*local_m2[u]*local_m3[u]*G2.v[u].size()));
                    else
                        pi = 1/(1+exp(+2*beta*spins2[u]*local_m2[u]*G2.v[u].size()));

                    rand_num = (double)G2.gen()/G2.gen.max();
                    if(rand_num < pi){
                        flip_spin(u,2,i);
                    }
                }
                else if (ppi == 2){
                    u = G3.gen()%int(N);

                    if(inter_connected13[u] and inter_connected23[u])
                        //pi = 1/(1+exp(+2*beta*spins3[u]*local_m1[u]*local_m2[u]*local_m3[u]*G3.v[u].size()));
                        pi = 1/(1+exp(+2*beta*spins3[u]*(M1/N)*(M2/N)*local_m3[u]*G3.v[u].size()));
                    else if(inter_connected13[u] and !inter_connected23[u])
                        pi = 1/(1+exp(+2*beta*spins3[u]*local_m1[u]*local_m3[u]*G3.v[u].size()));
                    else if(!inter_connected13[u] and inter_connected23[u])
                        pi = 1/(1+exp(+2*beta*spins3[u]*local_m2[u]*local_m3[u]*G3.v[u].size()));
                    else
                        pi = 1/(1+exp(+2*beta*spins3[u]*local_m3[u]*G3.v[u].size()));

                    rand_num = (double)G3.gen()/G3.gen.max();
                    if(rand_num < pi){
                        flip_spin(u,3,i);
                    }
                }


            }
        }

        T_vec.push_back(T);
        M_vec.push_back(M1/N);

        if(M1/giant1 > epsilon){
            Tmin = T;
            T = (Tmax + Tmin)/2;
        }
        else{
            Tmax = T;
            T = (Tmax + Tmin)/2;
        }

        revive_ordered_up();
        //revive_disordered();


    }

    ofstream dataInfo;
    stringstream ss;
    const char* cfile;
    string fileR = where_to + "/spins";
    ss<<TestNum;
    fileR = fileR + ss.str() + ".txt";
    cfile = fileR.c_str();
    dataInfo.open(cfile, ofstream::out);
    for(int i=0; i<T_vec.size();i++)
    {
        dataInfo<<T_vec[i]<<"	"<<(double) M_vec[i]<<endl;
    }

    dataInfo.close();
    cout<<"Tc="<<T<<",  q="<<q<<", L = "<<L<<endl;
    return T;
}


void ising_multiplex_3_q::scan_beta_highorder(double Ti,double Tf,double dT,string where_to, int TestNum, double NOI, double NOF){



    double rand_num, avg_local_m;
    int u;
    double pi;

    vector<double> T_vec, M_vec;
    double epsilon = 0.2, epsilon_T = 0.01;
    for(double T = Ti; T >=Tf;T-=dT){
        T_vec.push_back(T);
        double beta = 1/T;
        for(int j=0;j<NOI;j++){
            //cout<<"NOI = "<<j<<"  ,M1="<<M1/giant1<<"  ,M2="<<M2/giant2<<",  beta="<<beta<<endl;

            //avg_local_m = calc_average_local_M(1);
            for(int i=0;i<3*NOF;i++){

                int ppi = G1.gen()%3;
                if(ppi == 0){
                    u = G1.gen()%int(N);

                    if(inter_connected12[u] and inter_connected13[u])
                        pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*local_m2[u]*local_m3[u]*G1.v[u].size()));
                        //pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*(M2/N)*(M3/N)*G1.v[u].size()));
                    else if(inter_connected12[u] and !inter_connected13[u])
                        pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*local_m2[u]*G1.v[u].size()));
                    else if(!inter_connected12[u] and inter_connected13[u])
                        pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*local_m3[u]*G1.v[u].size()));
                    else
                        pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*G1.v[u].size()));

                    rand_num = (double)G1.gen()/G1.gen.max();
                    if(rand_num < pi){
                        flip_spin(u,1,i);
                    }
                }
                else if (ppi == 1){
                    u = G2.gen()%int(N);

                    if(inter_connected12[u] and inter_connected23[u])
                        pi = 1/(1+exp(+2*beta*spins2[u]*local_m1[u]*local_m2[u]*local_m3[u]*G2.v[u].size()));
                        //pi = 1/(1+exp(+2*beta*spins2[u]*(M1/N)*local_m2[u]*(M3/N)*G2.v[u].size()));
                    else if(inter_connected12[u] and !inter_connected23[u])
                        pi = 1/(1+exp(+2*beta*spins2[u]*local_m1[u]*local_m2[u]*G2.v[u].size()));
                    else if(!inter_connected12[u] and inter_connected23[u])
                        pi = 1/(1+exp(+2*beta*spins2[u]*local_m2[u]*local_m3[u]*G2.v[u].size()));
                    else
                        pi = 1/(1+exp(+2*beta*spins2[u]*local_m2[u]*G2.v[u].size()));

                    rand_num = (double)G2.gen()/G2.gen.max();
                    if(rand_num < pi){
                        flip_spin(u,2,i);
                    }
                }
                else if (ppi == 2){
                    u = G3.gen()%int(N);

                    if(inter_connected13[u] and inter_connected23[u])
                        pi = 1/(1+exp(+2*beta*spins3[u]*local_m1[u]*local_m2[u]*local_m3[u]*G3.v[u].size()));
                        //pi = 1/(1+exp(+2*beta*spins3[u]*(M1/N)*(M2/N)*local_m3[u]*G3.v[u].size()));
                    else if(inter_connected13[u] and !inter_connected23[u])
                        pi = 1/(1+exp(+2*beta*spins3[u]*local_m1[u]*local_m3[u]*G3.v[u].size()));
                    else if(!inter_connected13[u] and inter_connected23[u])
                        pi = 1/(1+exp(+2*beta*spins3[u]*local_m2[u]*local_m3[u]*G3.v[u].size()));
                    else
                        pi = 1/(1+exp(+2*beta*spins3[u]*local_m3[u]*G3.v[u].size()));

                    rand_num = (double)G3.gen()/G3.gen.max();
                    if(rand_num < pi){
                        flip_spin(u,3,i);
                    }
                }


            }
        }
        cout<<"M = "<<M1/giant1<<",    T="<<T<<",  q="<<q<<", L = "<<L<<endl;
        M_vec.push_back(M1/N);

        //revive_ordered_up();
        //revive_disordered();


    }

    ofstream dataInfo;
    stringstream ss;
    const char* cfile;
    string fileR = where_to + "/3-layers_highorder_DOIC_k6_test";
    ss<<TestNum;
    fileR = fileR + ss.str() + ".txt";
    cfile = fileR.c_str();
    dataInfo.open(cfile, ofstream::out);
    for(int i=0; i<T_vec.size();i++)
    {
        dataInfo<<T_vec[i]<<"	"<<(double) M_vec[i]<<endl;
    }

    dataInfo.close();
}

double ising_multiplex_3_q::scan_beta_lion_not_mixed(double Tmin,double Tmax,string where_to, int TestNum, double NOI, double NOF){



    double rand_num, avg_local_m;
    int u;
    double pi;
    double T = (Tmax + Tmin)/2;
    double epsilon = 0.2, epsilon_T = 0.01;
    while(Tmax - Tmin > epsilon_T){

        double beta = 1/T;
        for(int j=0;j<NOI;j++){
            //cout<<"NOI = "<<j<<"  ,M1="<<M1/giant1<<"  ,M2="<<M2/giant2<<",  beta="<<beta<<endl;

            //avg_local_m = calc_average_local_M(1);
            for(int i=0;i<NOF;i++){

                u = G1.gen()%int(N);

                if(inter_connected12[u] and inter_connected13[u]){
                    if(local_m2_counter[u] != 0 and local_m3_counter[u] != 0)
                        pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*(local_m2_thermal_avg[u]/local_m2_counter[u])*(local_m3_thermal_avg[u]/local_m3_counter[u])*G1.v[u].size()));
                    else if(local_m2_counter[u] != 0 and local_m3_counter[u] == 0)
                        pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*(local_m2_thermal_avg[u]/local_m2_counter[u])*local_m3[u]*G1.v[u].size()));
                    else if(local_m2_counter[u] == 0 and local_m3_counter[u] != 0)
                        pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*local_m2[u]*(local_m3_thermal_avg[u]/local_m3_counter[u])*G1.v[u].size()));
                    else
                        pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*local_m2[u]*local_m3[u]*G1.v[u].size()));
                }
                else if(inter_connected12[u] and !inter_connected13[u]){
                    if(local_m2_counter[u] != 0)
                        pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*(local_m2_thermal_avg[u]/local_m2_counter[u])*G1.v[u].size()));
                    else
                        pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*local_m2[u]*G1.v[u].size()));
                }
                else if(!inter_connected12[u] and inter_connected13[u]){
                    if(local_m3_counter[u] != 0)
                        pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*(local_m3_thermal_avg[u]/local_m3_counter[u])*G1.v[u].size()));
                    else
                        pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*local_m3[u]*G1.v[u].size()));
                }
                else
                    pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*G1.v[u].size()));

                rand_num = (double)G1.gen()/G1.gen.max();
                if(rand_num < pi){
                    flip_spin(u,1,i);
                }

            }

            for(int i = 0;i < N; i++){
                local_m2_thermal_avg[i] = 0;
                local_m2_counter[i] = 0;
            }

            for(int i=0;i<NOF;i++){

                u = G2.gen()%int(N);
                if(inter_connected12[u] and inter_connected23[u]){
                    if(local_m1_counter[u] != 0 and local_m3_counter[u] != 0)
                        pi = 1/(1+exp(+2*beta*spins2[u]*local_m2[u]*(local_m1_thermal_avg[u]/local_m1_counter[u])*(local_m3_thermal_avg[u]/local_m3_counter[u])*G2.v[u].size()));
                    else if(local_m1_counter[u] != 0 and local_m3_counter[u] == 0)
                        pi = 1/(1+exp(+2*beta*spins2[u]*local_m2[u]*(local_m1_thermal_avg[u]/local_m1_counter[u])*local_m3[u]*G2.v[u].size()));
                    else if(local_m1_counter[u] == 0 and local_m3_counter[u] != 0)
                        pi = 1/(1+exp(+2*beta*spins2[u]*local_m1[u]*local_m2[u]*(local_m3_thermal_avg[u]/local_m3_counter[u])*G2.v[u].size()));
                    else
                        pi = 1/(1+exp(+2*beta*spins2[u]*local_m1[u]*local_m2[u]*local_m3[u]*G2.v[u].size()));
                }
                else if(inter_connected12[u] and !inter_connected23[u]){
                    if(local_m1_counter[u] != 0)
                        pi = 1/(1+exp(+2*beta*spins2[u]*local_m2[u]*(local_m1_thermal_avg[u]/local_m1_counter[u])*G2.v[u].size()));
                    else
                        pi = 1/(1+exp(+2*beta*spins2[u]*local_m1[u]*local_m2[u]*G2.v[u].size()));
                }
                else if(!inter_connected12[u] and inter_connected23[u]){
                    if(local_m3_counter[u] != 0)
                        pi = 1/(1+exp(+2*beta*spins2[u]*local_m2[u]*(local_m3_thermal_avg[u]/local_m3_counter[u])*G2.v[u].size()));
                    else
                        pi = 1/(1+exp(+2*beta*spins2[u]*local_m2[u]*local_m3[u]*G2.v[u].size()));
                }
                else
                    pi = 1/(1+exp(+2*beta*spins2[u]*local_m2[u]*G2.v[u].size()));

                rand_num = (double)G2.gen()/G2.gen.max();
                if(rand_num < pi){
                    flip_spin(u,2,i);
                }
            }

            for(int i = 0;i < N; i++){
                local_m3_thermal_avg[i] = 0;
                local_m3_counter[i] = 0;
            }

            for(int i=0;i<NOF;i++){

                u = G3.gen()%int(N);

                if(inter_connected13[u] and inter_connected23[u]){
                    if(local_m1_counter[u] != 0 and local_m2_counter[u] != 0)
                        pi = 1/(1+exp(+2*beta*spins3[u]*local_m3[u]*(local_m1_thermal_avg[u]/local_m1_counter[u])*(local_m2_thermal_avg[u]/local_m2_counter[u])*G3.v[u].size()));
                    else if(local_m1_counter[u] != 0 and local_m2_counter[u] == 0)
                        pi = 1/(1+exp(+2*beta*spins3[u]*local_m3[u]*(local_m1_thermal_avg[u]/local_m1_counter[u])*local_m2[u]*G3.v[u].size()));
                    else if(local_m1_counter[u] == 0 and local_m2_counter[u] != 0)
                        pi = 1/(1+exp(+2*beta*spins3[u]*local_m1[u]*local_m3[u]*(local_m2_thermal_avg[u]/local_m2_counter[u])*G3.v[u].size()));
                    else
                        pi = 1/(1+exp(+2*beta*spins3[u]*local_m1[u]*local_m2[u]*local_m3[u]*G3.v[u].size()));
                }
                else if(inter_connected13[u] and !inter_connected23[u]){
                    if(local_m1_counter[u] != 0)
                        pi = 1/(1+exp(+2*beta*spins3[u]*local_m3[u]*(local_m1_thermal_avg[u]/local_m1_counter[u])*G3.v[u].size()));
                    else
                        pi = 1/(1+exp(+2*beta*spins3[u]*local_m1[u]*local_m3[u]*G3.v[u].size()));
                }
                else if(!inter_connected13[u] and inter_connected23[u]){
                    if(local_m2_counter[u] != 0)
                        pi = 1/(1+exp(+2*beta*spins3[u]*local_m3[u]*(local_m2_thermal_avg[u]/local_m2_counter[u])*G3.v[u].size()));
                    else
                        pi = 1/(1+exp(+2*beta*spins3[u]*local_m2[u]*local_m3[u]*G3.v[u].size()));
                }
                else
                    pi = 1/(1+exp(+2*beta*spins3[u]*local_m3[u]*G3.v[u].size()));

                rand_num = (double)G3.gen()/G3.gen.max();
                if(rand_num < pi){
                    flip_spin(u,3,i);
                }
            }

            for(int i = 0;i < N; i++){
                local_m1_thermal_avg[i] = 0;
                local_m1_counter[i] = 0;
            }
            /*
                for(int i=0;i<NOF;i++){

                    u = G1.gen()%int(N);

                    if(inter_connected12[u] and inter_connected13[u])
                        pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*(M2/N)*(M3/N)*G1.v[u].size()));
                    else if(inter_connected12[u] and !inter_connected13[u])
                        pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*(M2/N)*G1.v[u].size()));
                    else if(!inter_connected12[u] and inter_connected13[u])
                        pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*(M3/N)*G1.v[u].size()));
                    else
                        pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*G1.v[u].size()));

                    rand_num = (double)G1.gen()/G1.gen.max();
                    if(rand_num < pi){
                        flip_spin(u,1,i);
                    }

            }
            for(int i=0;i<NOF;i++){

                    u = G2.gen()%int(N);

                    if(inter_connected12[u] and inter_connected23[u])
                        pi = 1/(1+exp(+2*beta*spins2[u]*(M1/N)*local_m2[u]*(M3/N)*G2.v[u].size()));
                    else if(inter_connected12[u] and !inter_connected23[u])
                        pi = 1/(1+exp(+2*beta*spins2[u]*(M1/N)*local_m2[u]*G2.v[u].size()));
                    else if(!inter_connected12[u] and inter_connected23[u])
                        pi = 1/(1+exp(+2*beta*spins2[u]*local_m2[u]*(M3/N)*G2.v[u].size()));
                    else
                        pi = 1/(1+exp(+2*beta*spins2[u]*local_m2[u]*G2.v[u].size()));

                    rand_num = (double)G2.gen()/G2.gen.max();
                    if(rand_num < pi){
                        flip_spin(u,2,i);
                    }
                }
                for(int i=0;i<NOF;i++){

                    u = G3.gen()%int(N);

                    if(inter_connected13[u] and inter_connected23[u])
                        pi = 1/(1+exp(+2*beta*spins3[u]*(M1/N)*(M2/N)*local_m3[u]*G3.v[u].size()));
                    else if(inter_connected13[u] and !inter_connected23[u])
                        pi = 1/(1+exp(+2*beta*spins3[u]*(M1/N)*local_m3[u]*G3.v[u].size()));
                    else if(!inter_connected13[u] and inter_connected23[u])
                        pi = 1/(1+exp(+2*beta*spins3[u]*(M2/N)*local_m3[u]*G3.v[u].size()));
                    else
                        pi = 1/(1+exp(+2*beta*spins3[u]*local_m3[u]*G3.v[u].size()));

                    rand_num = (double)G3.gen()/G3.gen.max();
                    if(rand_num < pi){
                        flip_spin(u,3,i);
                    }
                }
                */

        }

        //cout<<"M = "<<M1/giant1<<",    T="<<T<<",  q="<<q<<",   k = "<<k<<", L = "<<L<<endl;
        if(M1/giant1 > epsilon){
            Tmin = T;
            T = (Tmax + Tmin)/2;
        }
        else{
            Tmax = T;
            T = (Tmax + Tmin)/2;
        }

        revive_ordered_up();
        //revive_disordered();

        for(int i = 0;i < N; i++){
            local_m1_thermal_avg[i] = 0;
            local_m1_counter[i] = 0;
            local_m2_thermal_avg[i] = local_m2[i];
            local_m2_counter[i] = 1;
            local_m3_thermal_avg[i] = local_m3[i];
            local_m3_counter[i] = 1;
        }


    }

    /*ofstream dataInfo;
    stringstream ss;
    const char* cfile;
    string fileR = where_to + "/spins";
    ss<<TestNum;
    fileR = fileR + ss.str() + ".txt";
    cfile = fileR.c_str();
    dataInfo.open(cfile, ofstream::out);
    for(int i=0; i<T_vec.size();i++)
    {
        dataInfo<<T_vec[i]<<"	"<<(double) M_vec[i]<<endl;
    }

    dataInfo.close();*/
    //cout<<"Tc="<<T<<",  q="<<q<<",   k = "<<k<<", L = "<<L<<endl;
    return T;
}


void ising_multiplex_3_q::scan_beta_thermal(double Ti,double Tf,double dT, string where_to, int TestNum, double NOI, double NOF){



    double rand_num, avg_local_m;
    int u;
    double pi;
    double epsilon = 0.2, epsilon_T = 0.01;
    vector<double> T_vec,  M_vec;
    for(double T = Ti; T>=Tf; T-=dT){
        T_vec.push_back(T);
        double beta = 1/T;
        for(int j=0;j<NOI;j++){
            //cout<<"NOI = "<<j<<"  ,M1="<<M1/giant1<<"  ,M2="<<M2/giant2<<",  beta="<<beta<<endl;

            //avg_local_m = calc_average_local_M(1);
            for(int i=0;i<NOF;i++){

                u = G1.gen()%int(N);

                if(inter_connected12[u] and inter_connected13[u]){
                    if(local_m2_counter[u] != 0 and local_m3_counter[u] != 0)
                        pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*(local_m2_thermal_avg[u]/local_m2_counter[u])*(local_m3_thermal_avg[u]/local_m3_counter[u])*G1.v[u].size()));
                    else if(local_m2_counter[u] != 0 and local_m3_counter[u] == 0)
                        pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*(local_m2_thermal_avg[u]/local_m2_counter[u])*local_m3[u]*G1.v[u].size()));
                    else if(local_m2_counter[u] == 0 and local_m3_counter[u] != 0)
                        pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*local_m2[u]*(local_m3_thermal_avg[u]/local_m3_counter[u])*G1.v[u].size()));
                    else
                        pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*local_m2[u]*local_m3[u]*G1.v[u].size()));
                }
                else if(inter_connected12[u] and !inter_connected13[u]){
                    if(local_m2_counter[u] != 0)
                        pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*(local_m2_thermal_avg[u]/local_m2_counter[u])*G1.v[u].size()));
                    else
                        pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*local_m2[u]*G1.v[u].size()));
                }
                else if(!inter_connected12[u] and inter_connected13[u]){
                    if(local_m3_counter[u] != 0)
                        pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*(local_m3_thermal_avg[u]/local_m3_counter[u])*G1.v[u].size()));
                    else
                        pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*local_m3[u]*G1.v[u].size()));
                }
                else
                    pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*G1.v[u].size()));

                rand_num = (double)G1.gen()/G1.gen.max();
                if(rand_num < pi){
                    flip_spin(u,1,i);
                }

            }

            for(int i = 0;i < N; i++){
                local_m2_thermal_avg[i] = 0;
                local_m2_counter[i] = 0;
            }

            for(int i=0;i<NOF;i++){

                u = G2.gen()%int(N);
                if(inter_connected12[u] and inter_connected23[u]){
                    if(local_m1_counter[u] != 0 and local_m3_counter[u] != 0)
                        pi = 1/(1+exp(+2*beta*spins2[u]*local_m2[u]*(local_m1_thermal_avg[u]/local_m1_counter[u])*(local_m3_thermal_avg[u]/local_m3_counter[u])*G2.v[u].size()));
                    else if(local_m1_counter[u] != 0 and local_m3_counter[u] == 0)
                        pi = 1/(1+exp(+2*beta*spins2[u]*local_m2[u]*(local_m1_thermal_avg[u]/local_m1_counter[u])*local_m3[u]*G2.v[u].size()));
                    else if(local_m1_counter[u] == 0 and local_m3_counter[u] != 0)
                        pi = 1/(1+exp(+2*beta*spins2[u]*local_m1[u]*local_m2[u]*(local_m3_thermal_avg[u]/local_m3_counter[u])*G2.v[u].size()));
                    else
                        pi = 1/(1+exp(+2*beta*spins2[u]*local_m1[u]*local_m2[u]*local_m3[u]*G2.v[u].size()));
                }
                else if(inter_connected12[u] and !inter_connected23[u]){
                    if(local_m1_counter[u] != 0)
                        pi = 1/(1+exp(+2*beta*spins2[u]*local_m2[u]*(local_m1_thermal_avg[u]/local_m1_counter[u])*G2.v[u].size()));
                    else
                        pi = 1/(1+exp(+2*beta*spins2[u]*local_m1[u]*local_m2[u]*G2.v[u].size()));
                }
                else if(!inter_connected12[u] and inter_connected23[u]){
                    if(local_m3_counter[u] != 0)
                        pi = 1/(1+exp(+2*beta*spins2[u]*local_m2[u]*(local_m3_thermal_avg[u]/local_m3_counter[u])*G2.v[u].size()));
                    else
                        pi = 1/(1+exp(+2*beta*spins2[u]*local_m2[u]*local_m3[u]*G2.v[u].size()));
                }
                else
                    pi = 1/(1+exp(+2*beta*spins2[u]*local_m2[u]*G2.v[u].size()));

                rand_num = (double)G2.gen()/G2.gen.max();
                if(rand_num < pi){
                    flip_spin(u,2,i);
                }
            }

            for(int i = 0;i < N; i++){
                local_m3_thermal_avg[i] = 0;
                local_m3_counter[i] = 0;
            }

            for(int i=0;i<NOF;i++){

                u = G3.gen()%int(N);

                if(inter_connected13[u] and inter_connected23[u]){
                    if(local_m1_counter[u] != 0 and local_m2_counter[u] != 0)
                        pi = 1/(1+exp(+2*beta*spins3[u]*local_m3[u]*(local_m1_thermal_avg[u]/local_m1_counter[u])*(local_m2_thermal_avg[u]/local_m2_counter[u])*G3.v[u].size()));
                    else if(local_m1_counter[u] != 0 and local_m2_counter[u] == 0)
                        pi = 1/(1+exp(+2*beta*spins3[u]*local_m3[u]*(local_m1_thermal_avg[u]/local_m1_counter[u])*local_m2[u]*G3.v[u].size()));
                    else if(local_m1_counter[u] == 0 and local_m2_counter[u] != 0)
                        pi = 1/(1+exp(+2*beta*spins3[u]*local_m1[u]*local_m3[u]*(local_m2_thermal_avg[u]/local_m2_counter[u])*G3.v[u].size()));
                    else
                        pi = 1/(1+exp(+2*beta*spins3[u]*local_m1[u]*local_m2[u]*local_m3[u]*G3.v[u].size()));
                }
                else if(inter_connected13[u] and !inter_connected23[u]){
                    if(local_m1_counter[u] != 0)
                        pi = 1/(1+exp(+2*beta*spins3[u]*local_m3[u]*(local_m1_thermal_avg[u]/local_m1_counter[u])*G3.v[u].size()));
                    else
                        pi = 1/(1+exp(+2*beta*spins3[u]*local_m1[u]*local_m3[u]*G3.v[u].size()));
                }
                else if(!inter_connected13[u] and inter_connected23[u]){
                    if(local_m2_counter[u] != 0)
                        pi = 1/(1+exp(+2*beta*spins3[u]*local_m3[u]*(local_m2_thermal_avg[u]/local_m2_counter[u])*G3.v[u].size()));
                    else
                        pi = 1/(1+exp(+2*beta*spins3[u]*local_m2[u]*local_m3[u]*G3.v[u].size()));
                }
                else
                    pi = 1/(1+exp(+2*beta*spins3[u]*local_m3[u]*G3.v[u].size()));

                rand_num = (double)G3.gen()/G3.gen.max();
                if(rand_num < pi){
                    flip_spin(u,3,i);
                }
            }

            for(int i = 0;i < N; i++){
                local_m1_thermal_avg[i] = 0;
                local_m1_counter[i] = 0;
            }
            /*
                for(int i=0;i<NOF;i++){

                    u = G1.gen()%int(N);

                    if(inter_connected12[u] and inter_connected13[u])
                        pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*(M2/N)*(M3/N)*G1.v[u].size()));
                    else if(inter_connected12[u] and !inter_connected13[u])
                        pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*(M2/N)*G1.v[u].size()));
                    else if(!inter_connected12[u] and inter_connected13[u])
                        pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*(M3/N)*G1.v[u].size()));
                    else
                        pi = 1/(1+exp(+2*beta*spins1[u]*local_m1[u]*G1.v[u].size()));

                    rand_num = (double)G1.gen()/G1.gen.max();
                    if(rand_num < pi){
                        flip_spin(u,1,i);
                    }

            }
            for(int i=0;i<NOF;i++){

                    u = G2.gen()%int(N);

                    if(inter_connected12[u] and inter_connected23[u])
                        pi = 1/(1+exp(+2*beta*spins2[u]*(M1/N)*local_m2[u]*(M3/N)*G2.v[u].size()));
                    else if(inter_connected12[u] and !inter_connected23[u])
                        pi = 1/(1+exp(+2*beta*spins2[u]*(M1/N)*local_m2[u]*G2.v[u].size()));
                    else if(!inter_connected12[u] and inter_connected23[u])
                        pi = 1/(1+exp(+2*beta*spins2[u]*local_m2[u]*(M3/N)*G2.v[u].size()));
                    else
                        pi = 1/(1+exp(+2*beta*spins2[u]*local_m2[u]*G2.v[u].size()));

                    rand_num = (double)G2.gen()/G2.gen.max();
                    if(rand_num < pi){
                        flip_spin(u,2,i);
                    }
                }
                for(int i=0;i<NOF;i++){

                    u = G3.gen()%int(N);

                    if(inter_connected13[u] and inter_connected23[u])
                        pi = 1/(1+exp(+2*beta*spins3[u]*(M1/N)*(M2/N)*local_m3[u]*G3.v[u].size()));
                    else if(inter_connected13[u] and !inter_connected23[u])
                        pi = 1/(1+exp(+2*beta*spins3[u]*(M1/N)*local_m3[u]*G3.v[u].size()));
                    else if(!inter_connected13[u] and inter_connected23[u])
                        pi = 1/(1+exp(+2*beta*spins3[u]*(M2/N)*local_m3[u]*G3.v[u].size()));
                    else
                        pi = 1/(1+exp(+2*beta*spins3[u]*local_m3[u]*G3.v[u].size()));

                    rand_num = (double)G3.gen()/G3.gen.max();
                    if(rand_num < pi){
                        flip_spin(u,3,i);
                    }
                }
                */

        }

        cout<<"M = "<<M1/giant1<<",    T="<<T<<",  q="<<q<<", L = "<<L<<endl;

        //revive_ordered_up();
        //revive_disordered();

        for(int i = 0;i < N; i++){
            local_m1_thermal_avg[i] = 0;
            local_m1_counter[i] = 0;
            local_m2_thermal_avg[i] = local_m2[i];
            local_m2_counter[i] = 1;
            local_m3_thermal_avg[i] = local_m3[i];
            local_m3_counter[i] = 1;
        }

        M_vec.push_back(M1/N);


    }

    ofstream dataInfo;
    stringstream ss;
    const char* cfile;
    string fileR = where_to + "/3-layers_thermal_DOIC_k6_test";
    ss<<TestNum;
    fileR = fileR + ss.str() + ".txt";
    cfile = fileR.c_str();
    dataInfo.open(cfile, ofstream::out);
    for(int i=0; i<T_vec.size();i++)
    {
        dataInfo<<T_vec[i]<<"	"<<(double) M_vec[i]<<endl;
    }

    dataInfo.close();
}

void ising_multiplex_3_q::revive_ordered_up(){
    for(int i=0;i<N;i++){    //flip all the spins back up
        if(spins1[i] == -1)
            flip_spin(i,1,i);
        if(spins2[i] == -1)
            flip_spin(i,2,i);
        if(spins3[i] == -1)
            flip_spin(i,3,i);
    }
}

void ising_multiplex_3_q::revive_ordered_down(){
    for(int i=0;i<N;i++){    //flip all the spins back up
        if(spins1[i] == 1)
            flip_spin(i,1,i);
        if(spins2[i] == 1)
            flip_spin(i,2,i);
        if(spins2[i] == 1)
            flip_spin(i,3,i);
    }
}

void ising_multiplex_3_q::revive_disordered(){
    double rand_num; //for each nodes intialized the spin +1 or -1 with equal probability

    for(int i=0;i<N;i++){ //network1
        rand_num = (double)G1.gen()/G1.gen.max();
        if(rand_num < 0.5){
            flip_spin(i,1,i);
        }
    }

    for(int i=0;i<N;i++){ //network2
        rand_num = (double)G1.gen()/G1.gen.max();
        if(rand_num < 0.5){
            flip_spin(i,2,i);
        }
    }

    for(int i=0;i<N;i++){ //network3
        rand_num = (double)G1.gen()/G1.gen.max();
        if(rand_num < 0.5){
            flip_spin(i,3,i);
        }
    }
}


#endif /* SRC_ising_multiplex_3_q_HPP_ */
