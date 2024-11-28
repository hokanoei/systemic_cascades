
using namespace std;
#include "my_structs.hpp"
#include <queue>
typedef std::mersenne_twister_engine< uint32_t, 32, 351, 175, 19, 0xccab8ee7, 11, 0xffffffff, 7, 0x31b6ab00, 15, 0xffe50000, 17, 1812433253 > mt11213b;

class HyperGraph {
public:
        int n, L;
        vector<int> components;
        vector<vector<pair<int,int> > > v;
        vector<int> active;
        vector <int> connected;
        bool periodic;
        mt11213b gen;

        HyperGraph();
        ~HyperGraph();
        void setHyperGraph(int size);
        void createHyperGraphD(int L, int D);
        bool isConnected(int u1, int u2, int u3);
        bool isConnected2(int u1, int u2);
        int testConnectivity();
        double k_core(int kcore);
        double k_core_avalanches(int kcore, vector<double>& av_vec);
        double k_core_avalanches2(int kcore, vector<double>& av_vec);
        void scan_p_kcore(int kcore, double dp, string where_to, int TestNum);
        double scan_p_kcore_lion(int kcore);
        void scan_p_kcore_lion_avalanches(int kcore, string where_to);
        void scan_p_kcore_lion_avalanches_single(int kcore, string where_to);


};
HyperGraph::HyperGraph() {n=0;}
void HyperGraph::setHyperGraph(int size) {
      if (size < 2)
          n = 2;
      else
          n = size;
      v.resize(n);
      connected.resize(n);
      components.resize(n);
      active = vector<int>(n,1);
      gen.seed(time(0));
}
HyperGraph::~HyperGraph() {
  for(int i=0;i<v.size();i++)
	  v[i].clear();
}

void HyperGraph::createHyperGraphD(int L_, int D) {

	L = L_;
	int s,t,u;
	vector<int> coordinates;
	coordinates.resize(D);
	if(pow(L,D) != n) {
        L = L + 1;
    }
	cout<<"realL = "<<L<<endl;

	//for each node connect it to its d neighbors in each axis
	for(int i=0;i<pow(L,D);i++)
	{
	    int temp = i;

		//find all the coordinates of i
	    for(int j = 0; j < D; j++){
			coordinates[j] = temp%L;
			temp = int(temp/L);
	    }

		// connect node i to each of his neighbors in each axis
	    for(int j = 0; j < D; j++){


			//edge of the network, coonect to the beginning
			if(coordinates[j] == L-1){
				t = i;
				u = i - (L-1)*pow(L,j);
				//addEdge(i,i - (L-1)*pow(L,j));
			}
			else{
				t = i;
				u = i + pow(L,j);
				//addEdge(i,i + pow(L,j));
			}

			do{
				s = gen()%n; //random node
			}while(s == t or s == u);

			pair <int, int> p = make_pair(s,t);
			v[u].push_back(p);
			p = make_pair(s,u);
			v[t].push_back(p);
			p = make_pair(u,t);
			v[s].push_back(p);
	    }
	}
}

bool HyperGraph::isConnected(int u1, int u2, int u3) {
    if(active[u1]==0||active[u2]==0||active[u3]==0)
    	return false;
    for(int i=0;i<v[u1].size();i++)
    	if((v[u1][i].first==u2 and v[u1][i].second == u3) or (v[u1][i].first==u3 and v[u1][i].second == u2))
    		return true;
    return false;
}
bool HyperGraph::isConnected2(int u1, int u2) {
    if(active[u1]==0||active[u2]==0)
    	return false;
    for(int i=0;i<v[u1].size();i++)
    	if(v[u1][i].first==u2 and v[u1][i].second == -1)
    		return true;
    return false;
}

void HyperGraph::scan_p_kcore(int kcore, double dp, string where_to, int TestNum){

	vector<int> node_order;
	vector<double> giant_vec, p_vec;

	for(int i=0;i<n;i++)
		node_order.push_back(i);

	std::shuffle(node_order.begin(), node_order.end(), gen);

	int node_count = int(dp*n);
	giant_vec.push_back(1);
	p_vec.push_back(1);
	int count = 0;
	while(count<n){
		//for(int j = 0;j < count; j++){
		//	active[node_order[j]] = 0;
		//}
		for(int j = 0;j < node_count; j++){
			active[node_order[count]] = 0;
			count++;
		}
		double p = 1 - (double)count/n;
		p_vec.push_back(p);
		double giant = k_core(kcore);

		//for(int j = 0;j < n; j++){
		//	active[j] = 1;
		//}
		giant_vec.push_back(giant);
		cout<<"p = "<<p<<", giant = "<<giant<<endl;

	}

	ofstream dataInfo;
	stringstream ss;
	const char* cfile;
	string fileR = where_to + "/2core";
	ss<<TestNum;
	fileR = fileR + ss.str() + ".txt";
	cfile = fileR.c_str();
	dataInfo.open(cfile, ofstream::out);
	for(int i=0; i<p_vec.size();i++)
	{
		dataInfo<<p_vec[i]<<"	"<<giant_vec[i]<<endl;
	}

	dataInfo.close();

}

double HyperGraph::k_core(int kcore){

	/*vector<int> TBR;
	vector<int> out_TBR;
	vector<int> degrees;
	double giant = n;

	for(int i = 0; i < n; i++){
		if(active[i] == 0){
			giant--;
			degrees.push_back(0);
			out_TBR.push_back(1);
			continue;
		}
		out_TBR.push_back(0);

		int degree = 0;
		for(int j = 0; j < v[i].size(); j++){
			if(active[v[i][j].first] == 1 and active[v[i][j].second] == 1)
				degree++;

		}

		degrees.push_back(degree);
		if(degree < kcore and active[i] == 1){
			TBR.push_back(i);
			active[i] = 0;
			giant--;
		}

	}

	while(TBR.size() != 0){

		int node = TBR.back();
		TBR.pop_back();
		out_TBR[node] = 1;
		for(int j = 0; j < v[node].size(); j++){
			if(active[v[node][j].first] == 1 and active[v[node][j].second] == 1){
				degrees[v[node][j].first]--;
				degrees[v[node][j].second]--;
			}
			if(active[v[node][j].first] == 0 and active[v[node][j].second] == 1){
				if(out_TBR[v[node][j].first] == 0){
					degrees[v[node][j].first]--;
					degrees[v[node][j].second]--;
				}
			}
			if(active[v[node][j].first] == 1 and active[v[node][j].second] == 0){
				if(out_TBR[v[node][j].second] == 0){
					degrees[v[node][j].first]--;
					degrees[v[node][j].second]--;
				}
			}

			if(degrees[v[node][j].first] < kcore and active[v[node][j].first] == 1){
				TBR.push_back(v[node][j].first);
				active[v[node][j].first] = 0;
				giant--;
			}

			if(degrees[v[node][j].second] < kcore and active[v[node][j].second] == 1){
				TBR.push_back(v[node][j].second);
				active[v[node][j].second] = 0;
				giant--;
			}
		}


	}

	//degree check
	for(int i = 0; i < n; i++){
		if(active[i] == 0)
			continue;
		int degree = 0;
		for(int j = 0; j < v[i].size(); j++){
			if(active[v[i][j].first] == 1 and active[v[i][j].second] == 1)
				degree++;
		}
		if(degree != degrees[i]){
			cout<<"Error"<<endl;
			cout<<"Real degree = "<<degree<<", curre degree = "<<degrees[i]<<endl;
		}
	}

	return giant/n;
*/

	vector<int> TBR;
	double giant = n;

	for(int i = 0; i < n; i++){
		if(active[i] == 0){
			giant--;
			TBR.push_back(i);
			continue;
		}

		int degree = 0;
		for(int j = 0; j < v[i].size(); j++){
			if(active[v[i][j].first] == 1 and active[v[i][j].second] == 1)
				degree++;

		}

		if(degree < kcore){
			TBR.push_back(i);
			active[i] = 0;
			giant--;
		}

	}

	while(TBR.size() != 0){

		int node = TBR.back();
		TBR.pop_back();
		for(int j = 0; j < v[node].size(); j++){


				int node1 = v[node][j].first;
				int node2 = v[node][j].second;

				if(active[node1] == 1){
					int degree = 0;
					for(int q = 0; q < v[node1].size(); q++){
						if(active[v[node1][q].first] == 1 and active[v[node1][q].second] == 1)
							degree++;
					}

					if(degree < kcore){
						TBR.push_back(node1);
						active[node1] = 0;
						giant--;
					}

				}

				if(active[node2] == 1){
					int degree = 0;
					for(int q = 0; q < v[node2].size(); q++){
						if(active[v[node2][q].first] == 1 and active[v[node2][q].second] == 1)
							degree++;
					}

					if(degree < kcore){
						TBR.push_back(node2);
						active[node2] = 0;
						giant--;
					}

				}
		}

	}

		for(int i = 0; i < n; i++){
			if(active[i] == 0){
				continue;
			}

			int degree = 0;
			for(int j = 0; j < v[i].size(); j++){
				if(active[v[i][j].first] == 1 and active[v[i][j].second] == 1)
					degree++;

			}
			if(degree < kcore)
				cout<<"Error"<<endl;
		}

	return giant/n;

}

double HyperGraph::k_core_avalanches2(int kcore, vector<double>& av_vec){

	vector<int> TBR;
	double giant = n;

	for(int i = 0; i < n; i++){
		if(active[i] == 0){
			giant--;
			continue;
		}

		int degree = 0;
		for(int j = 0; j < v[i].size(); j++){
			if(active[v[i][j].first] == 1 and active[v[i][j].second] == 1)
				degree++;

		}

		if(degree < kcore){
			TBR.push_back(i);
			active[i] = 0;
			giant--;
		}

	}

	av_vec.push_back(TBR.size());
	int counter_old = TBR.size(), counter_new = 0;
	while(TBR.size() != 0){

		int node = TBR.back();
		TBR.pop_back();
		counter_old--;
		for(int j = 0; j < v[node].size(); j++){


				int node1 = v[node][j].first;
				int node2 = v[node][j].second;

				if(active[node1] == 1){
					int degree = 0;
					for(int q = 0; q < v[node1].size(); q++){
						if(active[v[node1][q].first] == 1 and active[v[node1][q].second] == 1)
							degree++;
					}

					if(degree < kcore){
						TBR.push_back(node1);
						active[node1] = 0;
						counter_new++;
						giant--;
					}

				}

				if(active[node2] == 1){
					int degree = 0;
					for(int q = 0; q < v[node2].size(); q++){
						if(active[v[node2][q].first] == 1 and active[v[node2][q].second] == 1)
							degree++;
					}

					if(degree < kcore){
						TBR.push_back(node2);
						active[node2] = 0;
						counter_new++;
						giant--;
					}

				}
		}

		if(counter_old == 0){
			if(counter_new != 0)
				av_vec.push_back(counter_new);
			counter_old = counter_new;
			counter_new = 0;
		}

	}

		for(int i = 0; i < n; i++){
			if(active[i] == 0){
				continue;
			}

			int degree = 0;
			for(int j = 0; j < v[i].size(); j++){
				if(active[v[i][j].first] == 1 and active[v[i][j].second] == 1)
					degree++;

			}
			if(degree < kcore)
				cout<<"Error"<<endl;
		}

	return giant/n;

}
double HyperGraph::k_core_avalanches(int kcore, vector<double>& av_vec){

	vector<int> TBR;
	int degree;
	int source_idx = -1;
	while(source_idx < n){

		do{
			source_idx++;
			if(source_idx == n)
				break;

			if(active[source_idx] == 0){
				continue;
			}

			degree = 0;
			for(int j = 0; j < v[source_idx].size(); j++){
				if(active[v[source_idx][j].first] == 1 and active[v[source_idx][j].second] == 1)
					degree++;

			}

		}while(!(degree < kcore and active[source_idx]));

		if(source_idx == n)
				break;

		TBR.push_back(source_idx);
		int av_mass = 0;
		while(TBR.size() != 0){

			int node = TBR.back();
			TBR.pop_back();
			av_mass++;
			for(int j = 0; j < v[node].size(); j++){

				int node1 = v[node][j].first;
				int node2 = v[node][j].second;

				if(active[node1] == 1){
					degree = 0;
					for(int q = 0; q < v[node1].size(); q++){
						if(active[v[node1][q].first] == 1 and active[v[node1][q].second] == 1)
							degree++;
					}

					if(degree < kcore){
						TBR.push_back(node1);
						active[node1] = 0;
					}

				}

				if(active[node2] == 1){
					degree = 0;
					for(int q = 0; q < v[node2].size(); q++){
						if(active[v[node2][q].first] == 1 and active[v[node2][q].second] == 1)
							degree++;
					}

					if(degree < kcore){
						TBR.push_back(node2);
						active[node2] = 0;
					}

				}
			}
		}
		//if(av_mass > 1)
			av_vec.push_back(av_mass);
	}

	double giant = 0;
	for(int i = 0; i < n; i++){
			if(active[i]){
				giant++;
			}
	}

	return giant/n;
}
void HyperGraph::scan_p_kcore_lion_avalanches(int kcore, string where_to){


	double pmax = 1, pmin = 0;
	double p = (pmax + pmin)/2;
	double epsilon = 0.005;
	double epsilon_p = 0.000001;
	vector<double> av_vec_final, av_vec_temp;

	vector<int> node_order;
	for(int i=0;i<n;i++)
		node_order.push_back(i);
	std::shuffle(node_order.begin(), node_order.end(), gen);

	while(pmax - pmin > epsilon_p){

		for(int i = 0; i < (1-p)*n; i++)
			active[node_order[i]] = 0;

		av_vec_temp.clear();
		double giant = k_core_avalanches2(kcore, av_vec_temp);
		cout<<"p = "<<p<<", giant = "<<giant<<endl;
		for(int i = 0; i < n; i++)
			active[i] = 1;

		if(giant<epsilon){
			pmin = p;
			p = (p + pmax)/2;

		}
		else{
			pmax = p;
			p = (p + pmin)/2;
			av_vec_final = av_vec_temp;
		}
	}

	stringstream ss;
	ss<<gen();

	ofstream dataInfo;
	const char* cfile;
	string fileR = where_to + ss.str() + ".txt";
	cfile = fileR.c_str();
	dataInfo.open(cfile, ofstream::out);
	for(int i= 0; i < av_vec_final.size(); i++)
		dataInfo<<av_vec_final[i]<<endl;
	dataInfo.close();

}
void HyperGraph::scan_p_kcore_lion_avalanches_single(int kcore, string where_to){


  double pmax = 1, pmin = 0;
  double p = (pmax + pmin)/2;
  double epsilon = 0.005;
  double epsilon_p = 0.000001;
  vector<double> av_vec_final, av_vec_temp;

  vector<int> node_order;
  for(int i=0;i<n;i++)
    node_order.push_back(i);
  std::shuffle(node_order.begin(), node_order.end(), gen);

  while(pmax - pmin > epsilon_p){

        for(int i = 0; i < (1-p)*n; i++)
          active[node_order[i]] = 0;

        av_vec_temp.clear();
        double giant = k_core_avalanches2(kcore, av_vec_temp);
        cout<<"p = "<<p<<", giant = "<<giant<<endl;
        for(int i = 0; i < n; i++)
          active[i] = 1;

        if(giant<epsilon){
          pmin = p;
          p = (p + pmax)/2;

        }
        else{
              pmax = p;
              p = (p + pmin)/2;
              av_vec_final = av_vec_temp;
        }
  }

  double old_giant, new_giant;
  do{
    for(int i = 0; i < n; i++)
      active[i] = 1;
    std::shuffle(node_order.begin(), node_order.end(), gen);
    for(int i = 0; i < (1-p)*n; i++)
      active[node_order[i]] = 0;
    av_vec_temp.clear();
    old_giant = k_core_avalanches2(kcore, av_vec_temp);

    }while(old_giant<epsilon);


    int curr_idx = (1-p)*n;
    bool flag = true;
    do{
          if(active[node_order[curr_idx]] == 1){
            active[node_order[curr_idx]] = 0;
            new_giant = k_core_avalanches2(kcore, av_vec_temp);
            if(old_giant*n - new_giant*n > 1)
                flag = false;
            else
                old_giant = new_giant;
          }

        curr_idx++;
    }while(flag);

  cout<<"Av size: "<<old_giant*n - new_giant*n<<endl;
  stringstream ss;
  ss<<gen();

  ofstream dataInfo;
  const char* cfile;
  string fileR = where_to + ss.str() + ".txt";
  cfile = fileR.c_str();
  dataInfo.open(cfile, ofstream::out);
  dataInfo<<old_giant*n - new_giant*n<<endl;
  dataInfo.close();

}


double HyperGraph::scan_p_kcore_lion(int kcore){


	double pmax = 1, pmin = 0;
	double p = (pmax + pmin)/2;
	double epsilon = 0.005;
	double epsilon_p = 0.000001;

	vector<int> node_order;
	for(int i=0;i<n;i++)
		node_order.push_back(i);
	std::shuffle(node_order.begin(), node_order.end(), gen);

	while(pmax - pmin > epsilon_p){

		for(int i = 0; i < (1-p)*n; i++)
			active[node_order[i]] = 0;

		double giant = k_core(kcore);
		cout<<"p = "<<p<<", giant = "<<giant<<endl;
		for(int i = 0; i < n; i++)
			active[i] = 1;

		if(giant<epsilon){
			pmin = p;
			p = (p + pmax)/2;

		}
		else{
			pmax = p;
			p = (p + pmin)/2;
		}
	}
	return p;

}

int HyperGraph::testConnectivity()
{// we are going to change the BFS algorithem a little bit
	// we will remember the nodes of the giant and then we will make connected
	queue<int> Q;
	int j, size,expCount=0, cluster_id=0, giant_size=0, giant_id;
	int *explored = new int[n];//0-white, 1-grey, 2-black
	vector<int>* sources = new vector<int>;
	//vector<int>* giant_nodes;
	bool flag = true;

	for (int i = 0; i < n; ++i) //initialization explored
	{
		    if (active[i]==0)
		    {
		    	explored[i] = 2;
		    	expCount++;
		    	components[i]=cluster_id++;
		    }
		    else
		    	explored[i] = 0;
	}
	j=0;
	cluster_id--;
	while(expCount<n &&flag)
	{
	while(explored[j]!=0) //find a source for BFS algorithm
	{
	    j++;
	}
	cluster_id++;
	Q.push(j);
	components[j]=cluster_id;
	sources->push_back(j);
	size=1;
	explored[j] = 1;
	expCount++;

	while (!Q.empty()) {

	    int u = Q.front();
	    Q.pop();
	    explored[u] = 2;

	for(int i=0;i<v[u].size();i++)
	{
	    if(explored[v[u][i].first]==0)
	    {
	    	size++;
	    	Q.push(v[u][i].first);
	    	components[v[u][i].first] = cluster_id;
	    	explored[v[u][i].first] = 1;
	    	expCount++;
	    }

		if(v[u][i].second != -1){
			if(explored[v[u][i].second]==0)
			{
				size++;
				Q.push(v[u][i].second);
				components[v[u][i].second] = cluster_id;
				explored[v[u][i].second] = 1;
				expCount++;
			}
		}
	}

	}
	if(size>giant_size)
	{
		giant_size=size;
		giant_id=cluster_id;
	}
	}
	for(int k=0;k<n;k++)
	{
		if(components[k] == giant_id)
			connected[k] = 1;
		else
			connected[k] = 0;
	}
	//cout<<"unconnected: "<< n - giant_size<<endl;
	//cout<<giant_size<<endl;
	delete [] explored;
	return giant_size;
}


