#include "Graph.hpp"
#include <algorithm>
using namespace std;
using std::pair;
using std::make_pair;
#define pi_val 3.14159265359


class Lattice_Multiplex{
public:
	Graph G1, G2;  //networks
	vector<int> dep1;  //dependency links from network1 to network2
	vector<int> dep2;  //vice versa
    Lattice_Multiplex(int n, int d);
	~Lattice_Multiplex();
	void Create_dependency_r_meanfield();
	void Random_attack(double p);
	int kill_survivors_giant(Graph &G, int net_idx);
    void deactivate(int node_id, int net_idx);
	bool cascade_failure(int r_h);
	void revive();
	void PBC();
	vector<double> cascade_failure_random_giant(double p);
	double scan_p_giant_lion(string where_to, int TestNum);


};

Lattice_Multiplex::Lattice_Multiplex(int n, int d):dep1(vector<int>(n,-1)), dep2(vector<int>(n,-1))
{
	G1.setGraph(n);

	for(double i=0;i<1000000000;i++) //pause
		i=i*1;

	G2.setGraph(n);
	PBC();

	G1.createGraphSLFullnD(int(pow(n,(double)1.0/(double) d)),d);
	G2.createGraphSLFullnD(int(pow(n,(double)1.0/(double) d)),d);

	Create_dependency_r_meanfield();


}
Lattice_Multiplex::~Lattice_Multiplex(){}

void Lattice_Multiplex::PBC(){
	G1.periodic = true;
	G2.periodic = true;
}

void Lattice_Multiplex::Create_dependency_r_meanfield(){
	vector<int> node_order;
	int L=sqrt(G1.n);
	int count_test, good_test=0;
	node_order.resize(G1.n);
	for(int i=0;i<G1.n;i++)
		node_order[i] = i;
	std::shuffle(node_order.begin(), node_order.end(), G1.gen);

	for(int i=0;i<G1.n;i++)
	{

		dep2[node_order[i]] = i;
		dep1[i] = node_order[i];
		good_test++;

	}
	//cout<<"dependency created!!"<<endl;
	//cout<<"how many nodes dont have dep link?"<<endl;
	cout<<"bad nodes: "<<G1.n-good_test<<"  percent: "<<(double)(G1.n-good_test)/G1.n<<endl;


}

void Lattice_Multiplex::deactivate(int node_id, int net_idx)
{


    int L = int(sqrt(G1.n));

    //deactivate_no_dep(node_id);

    if(net_idx == 1)
    {

        if(G1.active[node_id] == 1){
            for(int j = 0; j < G1.v[node_id].size(); j++)
                if(G1.active[G1.v[node_id][j]] == 0){
                    //D_corr_length.merge_components(node_id, G1.v[node_id][j]);
                }
        }

        G1.active[node_id] = 0;
        if(dep1[node_id]>=0)
            G2.active[dep1[node_id]] = 0;


    }
    else
    {


        G2.active[node_id] = 0;
        if(dep2[node_id]>=0){
            if(G1.active[dep2[node_id]] == 1){
                for(int j = 0; j < G1.v[dep2[node_id]].size(); j++)
                    if(G1.active[G1.v[dep2[node_id]][j]] == 0){
                        //D_corr_length.merge_components(dep2[node_id], G1.v[dep2[node_id]][j]);
                    }
            }
            G1.active[dep2[node_id]] = 0;
        }


    }

}

void Lattice_Multiplex::Random_attack(double p){
	vector<int> node_order;
	node_order.resize(G1.n);
	for(int i=0;i<G1.n;i++)
		node_order[i] = i;

	std::shuffle(node_order.begin(), node_order.end(), G1.gen);
	for(int i=0;i<int(G1.n*(1-p));i++){
		if(G1.active[node_order[i]] == 1)
			deactivate(node_order[i],1);
	}


}

int Lattice_Multiplex::kill_survivors_giant(Graph &G, int net_idx)
{
		//cout<<"in"<<endl;
		queue<int> Q;
		int j, size,expCount=0, cluster_id=0, giant_size=0, giant_id;
		int *explored = new int[G.n];//0-white, 1-grey, 2-black
		int c = 0;
		for (int i = 0; i < G.n; ++i) //initialization explored
		{
			    if (G.active[i]==0)
			    {
			    	explored[i] = 2;
			    	expCount++;
			    	G.components[i]=cluster_id++;
			    }
			    else if(G.getDegree(i)==0)
			    {
			    	explored[i] = 2;
			    	expCount++;
			    	G.components[i]=cluster_id++;
			    }
			    else
			    	explored[i] = 0;
		}
		j=0;
		while(expCount<G.n)
		{
			while(explored[j]!=0) //find a source for BFS algorithm
			{
				j++;
			}
			cluster_id++;
			Q.push(j);
			G.components[j]=cluster_id;
			size=1;
			explored[j] = 1;
			expCount++;

			while (!Q.empty()) {

				int u = Q.front();
				Q.pop();
				explored[u] = 2;

				for(int i=0;i<G.v[u].size();i++)
				{
					if(explored[G.v[u][i]]==0)
					{
						size++;
						Q.push(G.v[u][i]);
						G.components[G.v[u][i]] = cluster_id;
						explored[G.v[u][i]] = 1;
						expCount++;
					}
				}

			}
			if(size>giant_size)
			{
				giant_size=size;
				giant_id=cluster_id;
			}
		}
		int count = 0;
		for(int k=0;k<G.n;k++)
		{
			if(G.components[k] != giant_id){
				deactivate(k, net_idx);
				count++;
			}
		}
		//cout<<"out"<<endl;
		//cout<<"count = "<<count<<endl;
		delete [] explored;
		return giant_size;
}

void Lattice_Multiplex::revive()
{
	for(int i=0; i<G1.n ; i++)
	{
		G1.active[i] = 1;
		G2.active[i] = 1;

	}
}

vector<double> Lattice_Multiplex::cascade_failure_random_giant(double p){

	int giant_size_net1_old=G1.n, giant_size_net1_new;
	int steps = 0;
	int count = 0;
	double Test_num = G1.gen();
	vector<double> giant_vec;
	Random_attack(p);
	//Flat_interface();

	while(1)
	{
		giant_vec.push_back(giant_size_net1_old);

		//giant_size_net1_new = remove_kcore_from_single_layer(2,1);
		giant_size_net1_new = kill_survivors_giant(G1, 1);


		//holes_dist(steps);

		//cout<<"size: "<<giant_size_net1_new<<endl;
		if(giant_size_net1_new == giant_size_net1_old)
			break;

		giant_size_net1_old = giant_size_net1_new;
		//remove_kcore_from_single_layer(2,2);
		kill_survivors_giant(G2, 2);

		steps++;
	}
	vector<double> temp = {(double)giant_size_net1_new,(double)steps};

	return temp;
}

double Lattice_Multiplex::scan_p_giant_lion(string where_to, int TestNum){
	vector<double> p_vec;
	vector<double> giant_vec;
	vector<double> results;
	vector<double> NOI_vec;
	vector<double> critical_droplet_m;
	vector<int> av_temp, av_final;
	double av_final_mass;



	double pmax = 1, pmin = 0;
	double p = (pmax + pmin)/2;
	double epsilon = 0.005;
	double epsilon_p = 0.000001;

	double giant_final;
	while(pmax - pmin > epsilon_p){

		p_vec.push_back(p);
		//results = cascade_failure_random_giant_avalanches(p,av_temp);
		results = cascade_failure_random_giant(p);
		giant_vec.push_back((double)(results[0]/G1.n));
		NOI_vec.push_back(results[1]);
		cout<<"p="<<p<<", P="<<giant_vec.back()<<", NOI = "<<results[1]<<endl;
		if(giant_vec.back()<epsilon){
			pmin = p;
			p = (p + pmax)/2;
		}
		else{
			pmax = p;
			p = (p + pmin)/2;
			giant_final = giant_vec.back();

		}
		revive();

	}
	p = p + 0.05;
	results = cascade_failure_random_giant(p);
	giant_final = (double)(results[0]/G1.n);


	if(giant_final > epsilon){
		ofstream dataInfo;
		stringstream ss;
		const char* cfile;
		string fileR = where_to + "/g";
		ss<<TestNum;
		fileR = fileR + ss.str() + ".txt";
		cfile = fileR.c_str();
		dataInfo.open(cfile, ofstream::out);

		dataInfo<<p<<"	"<<giant_final<<endl;
		dataInfo.close();
	}


	return giant_final;

}