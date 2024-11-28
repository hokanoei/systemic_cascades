

#ifndef SRC_MYGRAPH_HPP_
#define SRC_MYGRAPH_HPP_


using namespace std;
#include "my_structs.hpp"
#include <queue>
#include <pthread.h>
typedef std::mersenne_twister_engine< uint32_t, 32, 351, 175, 19, 0xccab8ee7, 11, 0xffffffff, 7, 0x31b6ab00, 15, 0xffe50000, 17, 1812433253 > mt11213b;

bool sortFunc(PitStruct Pit1, PitStruct Pit2)
{
	return (Pit1.length<Pit2.length);
}
class Graph {
public:
        int n;
        vector<int> components; //used for testconnectivity
        vector<vector<int> > v;
        vector<int> active;
        vector <int> connected;	//connected to the giant
        vector<PitStruct> pitVec;   //i ,j , r^2
        float p;  // 0<=p<=1
        bool periodic;
        mt11213b gen;

        //Forest Fire!!!!
    	//results vectors, updating at each time step
    	vector<double> r_t; //avarage distance from source at time t
    	vector<int> I_t;    // number of infected nodes at time t


    	vector<vector<int>> arti;
    	vector<int> cluster_id;
    	int count = 0;
    	int largest_cluster_id2 = 0;
    	int largest_cluster_size2 = 0;

	vector<bool> inter_nodes;
	vector<pair<int,int> > inter_vec;


	vector<pair<int,int> > edge_vec; // only if needed!!
        Graph();
        ~Graph();
        void showGraph();
        void setGraph(int size);
        bool isConnected(int, int);
        void addEdge(int , int );
        void killEdge(int, int);
        int getDegree(int);
        void showEdges(int);
        int getEdgesNum();
        int getNodesNum();
        void activate(int);
        void deactivate(int);
        int BFS();
	void add_ghost_field(double h);
        bool isAllExplored(bool *);
        void createGraphRR(double);
        void createGraphER_doubled(double);
        void createGraphER(double);
	void createGraphER_cliqueq(double k, int q_clique);
	void createGraphER_n_cliqueq(double k, int q_clique, int n_clique);
	void createGraphER_with_edge_vec(double);
        void createGraphSL(int);
        void createGraphSL_k8();
        void createGraphSL_k7();
        void createGraphSL_k6();
        void createGraphSL_k5();
        void createGraphSL_k3();
        void createGraphSLFull();
	void createGraphSLFull3D(int L);
	void createGraphSLFullnD(int L, int D);
        void createGraphSLZ_1D(float, float);
        void createGraphSLZ(float, float);
        void createGraphSLZ_uniform(float,float);
        void createGrpahSLZ_hetrogeny(float, float ,int);
        double createGraphSF_embedded_in_2D(double, int, double);
        double createGraphSF(float, int);
        double createGraphSF_test(float, int);
        double createGraphSF_test_doubled(float, int);
        void resetGraph();
        void activateNodes(float);
        void deactivateNodes(float);
        void changeToP(float);
        int binarySearch(int,int,double);
        bool isUp(int);
        bool isDown(int);
        bool isLeft(int);
        bool isRight(int);
        void SetPitVec();
        void resetPitVec();
        void showPitVec();
        void sortPitVec();
        int testConnectivity();
	int testConnectivity_deactivating();
        int SIR(double beta);
        vector<vector<double>> SIR_radius(double beta);
        int SIR_multiply(double beta, int source_num);
        void createGraphSLZcost(float, float,long double);
        void SIR_simulation(double beta);
        void Forest_fire(double beta);
        void SIR_node_perc(double beta);
        void BCCUtil(int u, int disc[], int low[], vector<vector<int> > *cluster_id, vector<pair<int,int> > *st, vector<bool> *arti, int parent[], int *count, int *largest_cluster_id, int *largest_cluster_size);
        void BCC();
        void BCCUtil_rec(int u, int disc[], int low[], vector<pair<int,int> > *st,int parent[]);
        void BCC_rec();
};
Graph::Graph() {n=0; p=0;}
void Graph::setGraph(int size) {
      if (size < 2) n = 2;
      else n = size;
      p=1;
      v.resize(n);
      connected.resize(n);
      components.resize(n);
      //v = new vector<vector<int> >(n);
      //components = new vector<int> (n);
      //connected = new vector<bool> (n);
      active = vector<int>(n,1);
      gen.seed(time(0));
}
Graph::~Graph() {
  for(int i=0;i<v.size();i++)
	  v[i].clear();
}
void Graph::showGraph()
{
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<v[i].size();j++)
			if(active[i]==1 && active[v[i][j]]==1)
			cout<<"("<<i<<","<<v[i][j]<<")"<<endl;
	}
}
bool Graph::isConnected(int u1, int u2) {
    if(active[u1]==0||active[u2]==0)
    	return false;
    for(int i=0;i<v[u1].size();i++)
    	if(v[u1][i]==u2)
    		return true;
    return false;
}
inline void Graph::addEdge(int u1, int u2) {
    v[u1].push_back(u2);
    v[u2].push_back(u1);
}
void Graph::killEdge(int u1, int u2){
	for(int i=0;i<v[u1].size();i++)
		if(v[u1][i] == u2){
			v[u1].erase(v[u1].begin() + i);
			break;
		}
	for(int i=0;i<v[u2].size();i++)
		if(v[u2][i] == u1){
			v[u2].erase(v[u2].begin() + i);
			break;
		}
}
int Graph::getDegree(int u)
{
	return v[u].size();

}
void Graph::showEdges(int u)
{
	for(int i=0;i<v[u].size();i++)
		cout<<"("<<u<<","<<v[u][i]<<")"<<endl;
}
int Graph::getEdgesNum()
{
	int num = 0;
	for(int i=0;i<n;i++)
	{
		num+=v[i].size();
	}
	return num/2;

}
int Graph::getNodesNum()
{
	int num = 0;
    for(int i=0;i<n;i++)
 	   if(active[i]==1)
 		   num++;
    return num;
}
void Graph::activate(int u)
{
	active[u] = 1;
}
void Graph::deactivate(int u)
{
	active[u] = 0;
}
int Graph::BFS() {
    queue<int> Q;
    int j, size, largestSize=0,expCount=0;
    int *explored = new int[n];//0-white, 1-grey, 2-black

    for (int i = 0; i < n; ++i) //initialization explored
    {

 	    if (active[i]==0)
 	    {
 		explored[i] = 2;
 		expCount++;
 	    }
 	   else if(getDegree(i)==0)
 	   {
 	   explored[i] = 2;
 	   expCount++;
 	   }
 	   else
 	         explored[i] = 0;
    }
	j=0;
    while(expCount<n)
    {
    while(explored[j]!=0) //find a source for BFS algorithm
    {
        j++;
    }
    Q.push(j);
    size=1;
    explored[j] = 1;
    expCount++;

   while (!Q.empty()) {

        int u = Q.front();
        Q.pop();
        explored[u] = 2;

    for(int i=0;i<v[u].size();i++)
    {
        if(explored[v[u][i]]==0)
        {
        	size++;
        	Q.push(v[u][i]);
        	explored[v[u][i]] = 1;
        	expCount++;
        }
    }

    }
    if(size>largestSize)
    {
    	largestSize = size;
    }
    }
    delete [] explored;
    return largestSize;
}
void Graph::add_ghost_field(double h){

	vector<int> temp;
	v.push_back(temp);
	double rand_num;
	for(int i = 0; i < n; i++){
		rand_num = (double)gen()/gen.max();
		if(rand_num < h)
			addEdge(i,n);
	}
	n++;
}
bool Graph::isAllExplored(bool * explored)
{
	for(int i=0;i<n;i++)
		if(!explored[i])
			return false;
	return true;
}
void Graph::createGraphRR(double k){
	vector<int> stub;

	for(int i=0; i<n;i++){
		for(int j=0;j<k;j++)
			stub.push_back(i);
	}

	std::shuffle(stub.begin(), stub.end(), gen);
	for(int i=0;i<stub.size()-1;i++){
		if(stub[i] != stub[i+1] and !isConnected(stub[i],stub[i+1]))
		{
			addEdge(stub[i],stub[i+1]);
			i++;
		}
	}


/*
	vector <int> distribution = vector<int> (n,0);
	for(int i=0;i<v.size();i++){
		distribution[v[i].size()]++;
	}

	ofstream dataInfo;
	stringstream ss;
	const char* cfile;
	string fileR = "/home/benaya/Documents/RR.txt";
	cfile = fileR.c_str();
	dataInfo.open(cfile, ofstream::out);
	for(int i=0; i<distribution.size();i++)
	{
		dataInfo<<i<<"	"<<distribution[i]<<endl;
	}
	dataInfo.close();
*/


}
void Graph::createGraphER(double k) {
	//k is average degree
	srand(time(0));
	int s,t,i=0;
	while(i<(n*k)/2)
	{
		s = gen()%n;
		t = gen()%n;
		if(s!=t && !isConnected(s,t))
		{
			addEdge(s,t);
			i++;
		}
	}
	//cout<<"start with "<<n<<" nodes and "<< (n*k)/2<<" edges."<<endl;
}

void Graph::createGraphER_cliqueq(double k, int q_clique) {
	//k is average degree
	srand(time(0));
	

	double Edges_num = (n*k)/2;
	cout<<"total edges = "<<Edges_num<<endl;
	vector<int> node_order;
	node_order.resize(n);
	for(int j=0;j<n;j++)
		node_order[j] = j;
	std::shuffle(node_order.begin(), node_order.end(), gen);

	for(int n1 = 0; n1 < q_clique-1; n1++){
		for(int n2 = n1 + 1; n2 < q_clique; n2++){
			addEdge(node_order[n1],node_order[n2]);
			Edges_num--;
		}
	}

	cout<<"Left edges = "<<Edges_num<<endl;
	int s,t,i=0;
	while(i<Edges_num)
	{
		s = gen()%n;
		t = gen()%n;
		if(s!=t && !isConnected(s,t))
		{
			addEdge(s,t);
			i++;
		}
	}
	//cout<<"start with "<<n<<" nodes and "<< (n*k)/2<<" edges."<<endl;
}


void Graph::createGraphER_n_cliqueq(double k, int q_clique, int n_clique) {
	//k is average degree
	srand(time(0));
	

	double Edges_num = (n*k)/2;
	cout<<"total edges = "<<Edges_num<<endl;

	for(int q_idx = 0; q_idx < n_clique; q_idx++){
		vector<int> node_order;
		node_order.resize(n);
		for(int j=0;j<n;j++)
			node_order[j] = j;
		std::shuffle(node_order.begin(), node_order.end(), gen);

		for(int n1 = 0; n1 < q_clique-1; n1++){
			for(int n2 = n1 + 1; n2 < q_clique; n2++){
				if(!isConnected(node_order[n1],node_order[n2])){
					addEdge(node_order[n1],node_order[n2]);
					Edges_num--;
				}
			}
		}
	}

	cout<<"Left edges = "<<Edges_num<<endl;
	int s,t,i=0;
	while(i<Edges_num)
	{
		s = gen()%n;
		t = gen()%n;
		if(s!=t && !isConnected(s,t))
		{
			addEdge(s,t);
			i++;
		}
	}
	//cout<<"start with "<<n<<" nodes and "<< (n*k)/2<<" edges."<<endl;
}

void Graph::createGraphER_with_edge_vec(double k) {
	//k is average degree
	srand(time(0));
	int s,t,i=0;
	while(i<(n*k)/2)
	{
		s = gen()%n;
		t = gen()%n;
		if(s!=t && !isConnected(s,t))
		{
			addEdge(s,t);
			edge_vec.push_back(make_pair(s,t));
			i++;
		}
	}
	//cout<<"start with "<<n<<" nodes and "<< (n*k)/2<<" edges."<<endl;
}
void Graph::createGraphER_doubled(double k) {
	int s,t,i=0;
	while(i<(n*k)/4)
	{
		s = gen()%int(n/2);
		t = gen()%int(n/2);
		if(s!=t && !isConnected(s,t))
		{
			addEdge(s,t);
			i++;
		}
	}
	i=0;
	while(i<(n*k)/4)
	{
		s = gen()%int(n/2) + int(n/2);
		t = gen()%int(n/2) + int(n/2);
		if(s!=t && !isConnected(s,t))
		{
			addEdge(s,t);
			i++;
		}
	}

}
void Graph::createGraphSL(int k) {
	//k is average degree
    int line = sqrt(n),rounds,nodeNum1,nodeNum2, num1, i=0;
    if(line*line!=n)
    	cout<<"not square lattice!!"<<endl;
    else{
    	if(k>4)
    	   rounds = n*4;
    	else
    	   rounds = n*k;
    	while(i!=rounds/2)
    	{
    	nodeNum1 = rand()%n;
    	num1 = (rand()%4) + 1;
    	if(num1 == 1)//Up
    	{
    		if(isUp(nodeNum1))
    			nodeNum2=line*(line-1) +nodeNum1;
    		else
    			nodeNum2=nodeNum1-line;
    	}
    	else if(num1 == 2)
    	{
    		if(isRight(nodeNum1))
    			nodeNum2=nodeNum1 - line +1;
    		else
    			nodeNum2 = nodeNum1 +1;
    	}
    	else if(num1 == 3)
    	{
    		if(isDown(nodeNum1))
    			nodeNum2=nodeNum1-(line*(line-1));
			else
				nodeNum2=nodeNum1+line;
    	}
    	else
    	{
    		if(isLeft(nodeNum1))
    			nodeNum2=nodeNum1+line-1;
    		else
    			nodeNum2=nodeNum1-1;
    	}
    	if(!isConnected(nodeNum1,nodeNum2))
    	{
    		addEdge(nodeNum1,nodeNum2);
    		i++;
    		//cout<<i<<" "<<nodeNum1<<" "<<nodeNum2<<" "<<num1<<endl;

    	}
    }
    }
    cout<<"start with "<<n<<" nodes and "<< (n*k)/2<<" edges."<<endl;
}
void Graph::createGraphSL_k5(){
	createGraphSLFull();
	int i=0, node,target, num,x,y, L = sqrt(n);
	while(i!=n/2){
		node = gen()%n;
		num = gen()%4;
		if(num == 0){
			x = ((node%L)+1) %L;
			y = ((int(node/L))-1+L)%L;
			target = x+L*y;
		}
		else if(num == 1){
			x = ((node%L)+1) %L;
			y = ((int(node/L))+1)%L;
			target = x+L*y;
		}
		else if(num == 2){
			x = ((node%L)-1+L) %L;
			y = ((int(node/L))+1)%L;
			target = x+L*y;
		}
		else{
			x = ((node%L)-1+L) %L;
			y = ((int(node/L))-1+L)%L;
			target = x+L*y;
		}
		if(!isConnected(node, target)){
			addEdge(node,target);
			i++;
		}
	}
}
void Graph::createGraphSL_k6(){
	createGraphSLFull();
	int i=0, node,target, num,x,y, L = sqrt(n);
	while(i!=n){
		node = gen()%n;
		num = gen()%4;
		if(num == 0){
			x = ((node%L)+1) %L;
			y = ((int(node/L))-1+L)%L;
			target = x+L*y;
		}
		else if(num == 1){
			x = ((node%L)+1) %L;
			y = ((int(node/L))+1)%L;
			target = x+L*y;
		}
		else if(num == 2){
			x = ((node%L)-1+L) %L;
			y = ((int(node/L))+1)%L;
			target = x+L*y;
		}
		else{
			x = ((node%L)-1+L) %L;
			y = ((int(node/L))-1+L)%L;
			target = x+L*y;
		}
		if(!isConnected(node, target)){
			addEdge(node,target);
			i++;
		}
	}
}
void Graph::createGraphSL_k7(){
	createGraphSLFull();
	int i=0, node,target, num,x,y, L = sqrt(n);
	while(i<(3*n/2)){
		node = gen()%n;
		num = gen()%4;
		if(num == 0){
			x = ((node%L)+1) %L;
			y = ((int(node/L))-1+L)%L;
			target = x+L*y;
		}
		else if(num == 1){
			x = ((node%L)+1) %L;
			y = ((int(node/L))+1)%L;
			target = x+L*y;
		}
		else if(num == 2){
			x = ((node%L)-1+L) %L;
			y = ((int(node/L))+1)%L;
			target = x+L*y;
		}
		else{
			x = ((node%L)-1+L) %L;
			y = ((int(node/L))-1+L)%L;
			target = x+L*y;
		}
		if(!isConnected(node, target)){
			addEdge(node,target);
			i++;
		}
	}
}
void Graph::createGraphSL_k8(){
	int L = sqrt(n), x ,y;
	for(int i=0;i<n;i++){
		x = i%L;
		y = int(i/L);
		addEdge(i,((x+1)%L) + ((y-1+L)%L)*L); // upper right
		addEdge(i,((x+1)%L) + y*L); // right
		addEdge(i,((x+1)%L) + ((y+1)%L)*L); // lower right
		addEdge(i,x + ((y+1)%L)*L); // down
	}
}
void Graph::createGraphSL_k3(){
	createGraphSLFull();
	int i=0, node,target, num,x,y, L = sqrt(n);
	while(i!=n/2){
		node = gen()%n;
		num = gen()%4;
		if(num == 0){
			x = ((node%L)+1) %L;
			y = int(node/L);
			target = x+L*y;
		}
		else if(num == 1){
			x = ((node%L)-1+L) %L;
			y = int(node/L);
			target = x+L*y;
		}
		else if(num == 2){
			x = node%L;
			y = ((int(node/L))+1)%L;
			target = x+L*y;
		}
		else{
			x = node%L;
			y = ((int(node/L))-1+L)%L;
			target = x+L*y;
		}
		if(isConnected(node, target)){
			killEdge(node,target);
			i++;
		}
	}
}
void Graph::createGraphSLFull(){
	int line = sqrt(n);
	for(int i=0;i<n;i++)
	{
		if(isRight(i))
			addEdge(i,i-line+1);
		else
			addEdge(i,i+1);
		if(isDown(i))
			addEdge(i,i-(line*(line-1)));
		else
            addEdge(i,i+line);
	}
}

void Graph::createGraphSLFull3D(int L){
	for(int i=0;i<L*L*L;i++)
	{
	    int x = i%L, y = (int(i/L))%L, z = int(i/(L*L));
            if(x == L-1)
		addEdge(i,z*L*L + y*L);
	    else
		addEdge(i,z*L*L + y*L + (x+1));
	    if(y == L-1)
		addEdge(i,z*L*L + x);
	    else
		addEdge(i,z*L*L + (y+1)*L + x);
            if(z == L-1)
		addEdge(i,y*L + x);
	    else
		addEdge(i,(z+1)*L*L + y*L + x);
	}
}
void Graph::createGraphSLFullnD(int L, int D){
	vector<int> coordinates;
	coordinates.resize(D);
	if(pow(L,D) != n)
		L = L+1;
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
			if(coordinates[j] == L-1)
				addEdge(i,i - (L-1)*pow(L,j));
			else
				addEdge(i,i + pow(L,j));
	    }
	}
}
void Graph::createGraphSLZ_1D(float k, float zeta){

	//periodic function
	double y,x;
	int i=0,source,target;
	std::uniform_real_distribution<double> randreal(0,1);
	randreal(gen);
	while(i<(n*k)/2)
	{
		source = gen()%n;

		y=randreal(gen);
		x = -zeta*log(1-y);

		if(gen()%2 == 0){
			target = (int(source + x))%n;
			if(source != target and !isConnected(source,target)){
				i++;
				addEdge(source,target);
			}
		}
		else{
			target = (int(source - x + n))%n;
			if(source != target and !isConnected(source,target)){
				i++;
				addEdge(source,target);
			}
		}
	}


}
void Graph::createGraphSLZ(float k, float zeta)
{

	bool out_of_boundaries;
	double y,x;
	int i=0,j,idx,source,target,c, line=sqrt(n),Px,Py;
	vector<PitStruct>::iterator low;
	std::uniform_real_distribution<double> randreal(0,1);
	randreal(gen);
	while(i<(n*k)/2)
	{


		out_of_boundaries=true;
		j=1;
		//cout<<"edge num: "<<i<<endl;

		source = gen()%n;
		//source = rand()%n;
		//y = ((double) rand() / (RAND_MAX));

		y=randreal(gen);
		x = -zeta*log(1-y);
		x*=x;
		if(x>=n){
			//cout<<"system too small!!!    x^2="<<x<<"n= "<<n<<endl;
			continue;
		}
		//idx = binarySearch(0,pitVec.size()-1,x);
		low =lower_bound(pitVec.begin(),pitVec.end(),x, MyCompareStruct());
		idx=low-pitVec.begin();

		if((pitVec[(idx+1)].length-x) < (x - (pitVec[(idx)].length)))
		{
			idx++;
			while(pitVec[(idx)].length==pitVec[(idx+j)].length)
				j++;
			j--;
			if(j!=0)
				{
					//j=rand()%(j+1);
					j=gen()%(j+1);
					idx+=j;
				}
		}
		else if(idx==0)
			j=0;
		else
		{
			while(pitVec[(idx)].length==pitVec[(idx-j)].length)
				j++;
			j--;
			if(j!=0)
			{
				//j=rand()%(j+1);
				j=gen()%(j+1);
				idx-=j;
			}
		}

		//cout<<"i= "<<pitVec->at(idx).at(0)<<"  j= "<<pitVec->at(idx).at(1)<<"  r^2= "<<pitVec->at(idx).at(2)<<" x^2= "<<x<<endl;
        //c = rand()%2;
		c = gen()%2;
        if(c==1)
        {
        	Px = pitVec[(idx)].i;
        	Py = pitVec[(idx)].j;
        }
        else
        {
        	Py = pitVec[(idx)].i;
        	Px = pitVec[(idx)].j;
        }
        //c = rand()%2;
        c = gen()%2;
        if(c==1)
        	Py*=-1;
        //c = rand()%2;
        c = gen()%2;
        if(c==1)
        	Px*=-1;
        /*cout<<i<<". Px= "<<Px<<" Py= "<<Py<<" source= "<<source;*/
        if((source%line)+Px>=line)
        {
        	out_of_boundaries = false;
        	target = ((source/line)*line) + (((source%line)+Px)%line);
        }
        else if(((source%line)+Px)<0)
        {
        	out_of_boundaries =false;
        	if(source<line)
        		target = line + source +Px;
        	else
        		target = ((source/line)*line) + ((source + Px)%line);
        }
        else
        	target = source+Px;
        /*cout<<" targetx= "<<target<<endl;*/
        if(((target/line) + Py )>=line)
        {
        	out_of_boundaries = false;
        	target = ((((target/line) + Py)%line)*line) + target%line;
        }
        else if(((target/line) + Py)<0)
        {
        	out_of_boundaries=false;
        	target = (((target/line) + Py + line)*line) + target%line;
        }
        else
        	target = target+ Py*line;
        if(!(periodic|| out_of_boundaries))
        {
        	continue;

        }
        /*cout<<" targety= "<<target<<endl;
        	cout<<"entered source="<<source<<"  target="<<target<<"degree= "<<getDegree(source)<<endl;*/
        if(source!=target && !isConnected(source,target))
        {
        	addEdge(source,target);
        	i++;
        }
	}
}
void Graph::createGraphSLZ_uniform(float k , float zeta)
{
	int i=0 ,source ,target ,dx ,dy;
	int L=sqrt(n);
	int t_x, t_y, s_x, s_y;

	while(i<(n*k)/2)
	{
		source = gen()%n;
		s_x = source%L;
		s_y = source/L;
		dx = round(fmod(gen(),2*zeta) - zeta); //[-zeta , +zeta] uniformly
		dy = round(fmod(gen(),2*zeta) - zeta);
		if((dx*dx + dy*dy)<zeta*zeta)
		{
			t_x = (s_x + dx + L)%L; //commet: working only for periodic!!!
			t_y = (s_y + dy + L)%L;
			target = t_y*L + t_x;
			if(source!=target && !isConnected(source,target))
			{
				addEdge(source,target);
				i++;
			}
		}

	}


}
void Graph::createGrpahSLZ_hetrogeny(float k_intra, float k_inter, int zeta){
	
	inter_nodes.resize(n,0);
	int L = sqrt(n);
	int num_of_ER = n/(zeta*zeta);
	int num_of_intra_edges = (zeta*zeta*k_intra)/2;
	//int num_of_inter_edges = num_of_ER*k_inter/2;
	int num_of_inter_edges = n*k_inter/2;
	int num_of_ER_L = L/zeta;

	vector<pair<int,int> > intra_vec;
	//cout<<"n = "<<num_of_inter_edges<<", zeta = "<<zeta<<endl;

	int i;

	int dx1, dx2, dy1, dy2;

	int n1, n2;

	for(int x=0;x<L;x+=zeta){
		for(int y=0;y<L;y+=zeta){

			i = 0;
			while(i<num_of_intra_edges){

				dx1 = gen()%zeta;
				dx2 = gen()%zeta;
				dy1 = gen()%zeta;
				dy2 = gen()%zeta;

				n1 = x+dx1 + (y+dy1)*L;
				n2 = x+dx2 + (y+dy2)*L;

				if(n1!=n2 and !isConnected(n1,n2)){
					addEdge(n1,n2);
					intra_vec.push_back(make_pair(n1,n2));
					i++;
				}
			}
		}
	}

	i = 0;
	int n;
	int op;
	while(i<num_of_inter_edges){

		n = gen()%num_of_ER; //choose ER
		op = gen()%4; // sides

		dx1 = gen()%zeta;
		dx2 = gen()%zeta;
		dy1 = gen()%zeta;
		dy2 = gen()%zeta;

		n1 = ((n%num_of_ER_L)*zeta + dx1) + (int(n/num_of_ER_L)*zeta + dy1)*L; //node location within n

		if(op == 0){ //up
			n= n%num_of_ER_L + ((int(n/num_of_ER_L)-1 + num_of_ER_L)%num_of_ER_L)*num_of_ER_L;
			n2 = ((n%num_of_ER_L)*zeta + dx2) + (int(n/num_of_ER_L)*zeta + dy2)*L;
		}
		else if(op == 1){ //right
			n= (n%num_of_ER_L+1)%num_of_ER_L + int(n/num_of_ER_L)*num_of_ER_L;
			n2 = ((n%num_of_ER_L)*zeta + dx2) + (int(n/num_of_ER_L)*zeta + dy2)*L;
		}
		else if(op == 2){ //down
			n= n%num_of_ER_L + ((int(n/num_of_ER_L)+1)%num_of_ER_L)*num_of_ER_L;
			n2 = ((n%num_of_ER_L)*zeta + dx2) + (int(n/num_of_ER_L)*zeta + dy2)*L;
		}
		else{// left
			n= (n%num_of_ER_L-1+num_of_ER_L)%num_of_ER_L + int(n/num_of_ER_L)*num_of_ER_L;
			n2 = ((n%num_of_ER_L)*zeta + dx2) + (int(n/num_of_ER_L)*zeta + dy2)*L;
		}

		if(!isConnected(n1, n2)){
			addEdge(n1,n2);
			inter_vec.push_back(make_pair(n1,n2));
			inter_nodes[n1] = 1;
			inter_nodes[n2] = 1;
			i++;
		}
	}


/*
	ofstream dataInfo;
	string where_to = "/home/benaya/Desktop/net_shlomo/inter.txt";
	const char* cfile = where_to.c_str();
	dataInfo.open(cfile, ofstream::out);
	for(int i=0;i<inter_vec.size();i++){
		dataInfo<<inter_vec[i].first<<"	"<<inter_vec[i].second<<endl;
	}
	dataInfo.close();

	ofstream dataInfo2;
	string where_to2 = "/home/benaya/Desktop/net_shlomo/intra.txt";
	const char* cfile2 = where_to2.c_str();
	dataInfo2.open(cfile2, ofstream::out);
	for(int i=0;i<intra_vec.size();i++){
		dataInfo2<<intra_vec[i].first<<"	"<<intra_vec[i].second<<endl;
	}
	dataInfo2.close();*/
}
double Graph::createGraphSF_embedded_in_2D(double gamma, int k_min,double zeta){
	std::uniform_real_distribution<double> randreal(0,1);
	vector<PitStruct>::iterator low;
	randreal(gen);
	vector<int> k_vec;
	double y, x, power = (double)1/(1-gamma);
	double E_expected = 0;

	for(int i=0; i<n;i++){
		y = randreal(gen);
		x = k_min*pow(y,power);

		if(int(x) < sqrt(n)){
			k_vec.push_back(int(x));
			E_expected+=int(x);
		}
		else
			i--;
	}
	E_expected/=2;

	vector<double> edge_vec(n,0);

	int count = 0, E_count = 0;
	int source;
	bool out_of_boundaries;
	int j, idx,target,c, line=sqrt(n), Px, Py;
	while((count < 100000) and (E_count < E_expected)){
		count++;
		source = gen()%n;
		if(k_vec[source] == 0){
			continue;
		}

		out_of_boundaries=true;
		j=1;
		y=randreal(gen);
		x = -zeta*log(1-y);

		x*=x;

		if(x>=n){
			continue;
		}

		low =lower_bound(pitVec.begin(),pitVec.end(),x, MyCompareStruct());
		idx=low-pitVec.begin();

		if((pitVec[(idx+1)].length-x) < (x - (pitVec[(idx)].length)))
		{
			idx++;
			while(pitVec[(idx)].length==pitVec[(idx+j)].length)
				j++;
			j--;
			if(j!=0)
				{
					j=gen()%(j+1);
					idx+=j;
				}
		}
		else if(idx==0)
			j=0;
		else
		{
			while(pitVec[(idx)].length==pitVec[(idx-j)].length)
				j++;
			j--;
			if(j!=0)
			{
				j=gen()%(j+1);
				idx-=j;
			}
		}

		c = gen()%2;
        if(c==1)
        {
        	Px = pitVec[idx].i;
        	Py = pitVec[idx].j;
        }
        else
        {
        	Py = pitVec[idx].i;
        	Px = pitVec[idx].j;
        }

        c = gen()%2;
        if(c==1)
        	Py*=-1;

        c = gen()%2;
        if(c==1)
        	Px*=-1;

        if((source%line)+Px>=line)
        {
        	out_of_boundaries = false;
        	target = ((source/line)*line) + (((source%line)+Px)%line);
        }
        else if(((source%line)+Px)<0)
        {
        	out_of_boundaries =false;
        	if(source<line)
        		target = line + source +Px;
        	else
        		target = ((source/line)*line) + ((source + Px)%line);
        }
        else
        	target = source+Px;

        if(((target/line) + Py )>=line)
        {
        	out_of_boundaries = false;
        	target = ((((target/line) + Py)%line)*line) + target%line;
        }
        else if(((target/line) + Py)<0)
        {
        	out_of_boundaries=false;
        	target = (((target/line) + Py + line)*line) + target%line;
        }
        else
        	target = target+ Py*line;
        if(!(periodic|| out_of_boundaries))
        {
        	continue;

        }

        if(source!=target and !isConnected(source,target) and k_vec[target] != 0)
        {
        	addEdge(source,target);
        	E_count++;
        	count = 0;
        	k_vec[source]--;
        	k_vec[target]--;
        	edge_vec[int(sqrt(Px*Px + Py*Py))]++;
        }
	}

	ofstream dataInfo;
	stringstream ss, ss2;
	string fileR = "/storage/home/benaya/SF_embedded_in_2D/link_length_dist/gamma";
	const char* cfile;
	ss<<gamma;
	ss2<<zeta;
	fileR = fileR + ss.str() +"/zeta" + ss2.str() + ".txt";
	cfile = fileR.c_str();
	dataInfo.open(cfile, ofstream::out);
	for(int j=0; j<edge_vec.size();j++)
	{
		if(edge_vec[j] != 0 )
			dataInfo<<j<<"	"<<edge_vec[j]<<endl;

	}
	dataInfo.close();


	cout<<"E_expected = "<<E_expected<<endl;
	cout<<"E_count = "<<E_count<<endl;
	cout<<"Connected percentage: "<<E_count/E_expected<<endl;
	return E_count/E_expected;
}
double Graph::createGraphSF(float gamma, int k_min){
	std::uniform_real_distribution<double> randreal(0,1);
	randreal(gen);
	vector<int> k_vec;
	vector<int> stub;
	double y, x, power = (double)1/(1-gamma);
	int counter = 0;
	double CS;

	int mink = 100, maxk = 0;
	for(int i=0; i<n;i++){
		y = randreal(gen);
		x = k_min*pow(y,power);

		if(int(x) < n){
			k_vec.push_back(int(x));
			vector <int> temp = vector<int>(int(x),i);
			stub.insert( stub.end(), temp.begin(), temp.end() );
		}
		else
			i--;
	}



	CS = stub.size();
	int idx, t, s;
	while(counter < 1){
		idx = 0;
		std::shuffle(stub.begin(), stub.end(), gen);

		while(idx<stub.size()-1 and stub.size() !=0){
			//cout<<"idx= "<<idx<<",  size=" << stub.size()<<endl;
			t = stub[idx];
			s = stub[idx+1];

	        if(s!=t && !isConnected(s,t))
	        {
	        	addEdge(s,t);
	        	stub.erase(stub.begin()+idx,stub.begin() + idx+2);
	        }
	        else
	        	idx++;
		}
		if(stub.size() < CS){
			CS = stub.size();
			counter = 0;
		}
		else
			counter++;
	}

	vector <int> distribution = vector<int> (n,0);
	for(int i=0;i<v.size();i++){
		distribution[v[i].size()]++;
	}
	float sum = 0, count = 0;
	for(int i=0;i<distribution.size();i++){
		sum+=i*distribution[i];
		count+=distribution[i];
	}

	ofstream dataInfo;
	stringstream ss;
	const char* cfile;
	string fileR = "/home/benaya/Documents/SF.txt";
	cfile = fileR.c_str();
	dataInfo.open(cfile, ofstream::out);
	for(int i=0; i<distribution.size();i++)
	{
		dataInfo<<i<<"	"<<distribution[i]<<endl;
	}
	dataInfo.close();
	return (double)(sum/count);

}
double Graph::createGraphSF_test(float gamma, int k_min){
	std::uniform_real_distribution<double> randreal(0,1);
	randreal(gen);
	vector<int> k_vec;
	vector<int> stub;
	double y, x, power = (double)1/(1-gamma);

	for(int i=0; i<n;i++){
		y = randreal(gen);
		x = k_min*pow(y,power);

		if(int(x) < sqrt(n)){
			k_vec.push_back(int(x));
			vector <int> temp = vector<int>(int(x),i);
			stub.insert( stub.end(), temp.begin(), temp.end() );
		}
		else
			i--;
	}
	std::shuffle(stub.begin(), stub.end(), gen);

	int s,t,counter = 0;
	for(int i=0;i<stub.size()-1;i++){
			t = stub[i];
			s = stub[i+1];
		if(s!=t && !isConnected(s,t)){
        	addEdge(s,t);
        	i++;
		}
	}

/*	vector <int> distribution = vector<int> (n,0);
	double sum = 0,count = 0;;
	for(int i=0;i<v.size();i++){
		distribution[v[i].size()]++;
		if(v[i].size() > 0){
			sum+=(1/(double)(1+v[i].size()));
			count++;
		}
	}



	ofstream dataInfo;
	stringstream ss;
	const char* cfile;
	string fileR = "/home/benaya/Documents/SF.txt";
	cfile = fileR.c_str();
	dataInfo.open(cfile, ofstream::out);
	for(int i=0; i<distribution.size();i++)
	{
		dataInfo<<i<<"	"<<distribution[i]<<endl;
	}
	dataInfo.close();*/
	return (double)4;

}
double Graph::createGraphSF_test_doubled(float gamma, int k_min){
	std::uniform_real_distribution<double> randreal(0,1);
	randreal(gen);
	vector<int> k_vec;
	vector<int> stub;
	double y, x, power = (double)1/(1-gamma);

	for(int i=0; i<int(n/2);i++){
		y = randreal(gen);
		x = k_min*pow(y,power);

		if(int(x) < sqrt(int(n/2))){
			k_vec.push_back(int(x));
			vector <int> temp = vector<int>(int(x),i);
			stub.insert( stub.end(), temp.begin(), temp.end() );
		}
		else
			i--;
	}
	std::shuffle(stub.begin(), stub.end(), gen);

	int s,t,count = 0;
	for(int i=0;i<stub.size()-1;i++){
			t = stub[i];
			s = stub[i+1];
		if(s!=t && !isConnected(s,t)){
        	addEdge(s,t);
        	i++;
        	count++;
		}
	}
	cout<<count<<endl;
	k_vec.clear();
	stub.clear();

	for(int i=0; i<int(n/2);i++){
			y = randreal(gen);
			x = k_min*pow(y,power);

			if(int(x) < sqrt(int(n/2))){
				k_vec.push_back(int(x));
				vector <int> temp = vector<int>(int(x),i);
				stub.insert( stub.end(), temp.begin(), temp.end() );
			}
			else
				i--;
		}
	std::shuffle(stub.begin(), stub.end(), gen);
	count = 0;
	for(int i=0;i<stub.size()-1;i++){
			t = stub[i] + int(n/2);
			s = stub[i+1] + int(n/2);
		if(s!=t && !isConnected(s,t)){
	        addEdge(s,t);
	        i++;
	        count++;
		}
	}
	cout<<count<<endl;
	return (double)4;

}
void Graph::resetGraph()
{

	active.clear();
	components.clear();
	connected.clear();
	n=0;
	p=0;
}
void Graph::activateNodes(float p1)
{
	//cout<<"rise to p="<<p+p1<<endl;
     int num, number = n*p1;
     while(number>0)
     {
    	 num = rand()%n;
    	 if(active[num]==0)
    	 {
    		 active[num]=1;
    		 number--;
    	 }
     }
     //cout<<"nodes number:"<<getNodesNum()<<endl;
    // cout<<"edges number:"<<getEdgesNum()<<endl;
}
void Graph::deactivateNodes(float p1)
{
	//cout<<"reduce to p=" <<p-p1<<endl;
     int num, number = n*p1;
     while(number>0)
     {
    	 num = rand()%n;
    	 if(active[num]==1)
    	 {
    		 active[num]=0;
    		 number--;
    	 }
     }
    // cout<<"nodes number:"<<getNodesNum()<<endl;
    // cout<<"edges number:"<<getEdgesNum()<<endl;
}
void Graph::changeToP(float p1)
{
	if(p1>p)
		activateNodes(p1-p);
	else
		deactivateNodes(p-p1);
	p=p1;
}
int Graph::binarySearch(int min, int max, double l)
{

    if(min==max || max-min==1)
    	return min;
	if(pitVec[((min+max)/2)].length>l)
		return binarySearch(min,(max+min)/2,l);
	else if(pitVec[((min+max)/2)].length<l)
		return binarySearch((max+min)/2,max,l);
	else
		return (max+min)/2;
}
bool Graph::isUp(int num)
{
	int line = sqrt(n);
	if(num<line)
		return true;
	return false;
}
bool Graph::isDown(int num)
{
	int line = sqrt(n);
	if(num>=n-line)
		return true;
	return false;
}
bool Graph::isLeft(int num)
{
	int line = sqrt(n);
	if(num%line==0)
		return true;
	return false;
}
bool Graph::isRight(int num)
{
	int line = sqrt(n);
	if((num+1)%line==0)
		return true;
	return false;
}
void Graph::SetPitVec()
{
	/*
	bool flag;
	int j;
	vector<int> temp;
	pitVec = new vector<vector<int> >;
	for(int i=0;i<sqrt(n);i++)
	{
		j=i;
		flag = true;
		while(flag)
		{
			//cout<<"i= "<<i<<"  j="<<j<<" l^2="<<(i*i + j*j)<<endl;
			if((i*i + j*j) > n)
				flag = false;
			else
			{
				temp.push_back(i);
				temp.push_back(j);
				temp.push_back(i*i + j*j);
				pitVec->push_back(temp);
				temp.clear();
				j=j+1;
			}
		}
		//cout<<GetTime()<<" i= "<<i<<endl;
	}*/

	for(int i=0;i<=sqrt(n);i++)
	{
		for(int j=i; (i*i + j*j) <= n; j++){
				pitVec.push_back(PitStruct(i,j));


		}
		//cout<<GetTime()<<" i= "<<i<<endl;
	}	
}
void Graph::resetPitVec()
{
	pitVec.clear();
}
void Graph::showPitVec()
{
	for(int i=0;i<pitVec.size();i++)
		cout<<"i is: "<<pitVec[i].i<<"  j is: "<<pitVec[i].j<<" r^2 is: "<<pitVec[i].length<<endl;
}
void Graph::sortPitVec()
{
	sort(pitVec.begin(),pitVec.end(),sortFunc);
}
int Graph::testConnectivity_deactivating()
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
		    else if(getDegree(i)==0)
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
	    if(explored[v[u][i]]==0)
	    {
	    	size++;
	    	Q.push(v[u][i]);
	    	components[v[u][i]] = cluster_id;
	    	explored[v[u][i]] = 1;
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
	for(int k=0;k<n;k++)
	{
		if(components[k] == giant_id)
			active[k] = 1;
		else
			active[k] = 0;
	}
	//cout<<"unconnected: "<< n - giant_size<<endl;
	//cout<<giant_size<<endl;
	delete [] explored;
	return giant_size;
}
int Graph::testConnectivity()
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
		    else if(getDegree(i)==0)
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
	    if(explored[v[u][i]]==0)
	    {
	    	size++;
	    	Q.push(v[u][i]);
	    	components[v[u][i]] = cluster_id;
	    	explored[v[u][i]] = 1;
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
int Graph::SIR(double beta){
	cout<<"Enter SIR function"<<endl;
	int r_counter = 0;
	vector<int> infected;
	int current_size;
	int t=0;
	int n1, n2;

	vector<bool> recovered = vector<bool>(n,0);
	vector<bool> infected_bool = vector<bool>(n,0);

	//omp_set_num_threads(threads);
	//#pragma omp parallel for
	vector<vector<int> > nodes = v;
		do{
			n1 = gen()%n;
		}while(!connected[n1]);
		n2 = nodes[n1][gen()%nodes[n1].size()];


	//infect single bond
	infected.push_back(n1);
	infected.push_back(n2);
	infected_bool[n1] = 1;
	infected_bool[n2] = 1;


	while(!infected.empty())
	{
		current_size = infected.size();
		for(int inf_idx=0;inf_idx<current_size;inf_idx++)
		{
			n1 = infected[inf_idx];
			for(int i=0;i<nodes[n1].size();i++)
			{
				if(!recovered[nodes[n1][i]] and !infected_bool[nodes[n1][i]])
				{
					if((double)gen()/gen.max()<beta)
					{
						infected.push_back(nodes[n1][i]);
						infected_bool[nodes[n1][i]] = 1;

					}
				}
			}
			recovered[n1] = 1;
			infected_bool[n1] = 0;

		}
		infected.erase (infected.begin(),infected.begin()+current_size);
		t++;

	}
	cout<<"Exit SIR function"<<endl;
	return t;

}
vector<vector<double>> Graph::SIR_radius(double beta){
	vector<double> radius_vec, full_radius_vec;
	vector<int> infected;
	vector<double> branching_factor;
	int current_size;
	int t=0;
	int n1, n2;
	long double radius_sum, radius_sum_full = 0;
	double max_r=0;
	int radius_sum_full_counter = 0;
	int L = int(sqrt(n));
	int infected_counter_currently, infected_counter_past = 1;
	//omp_set_num_threads(threads);
	//#pragma omp parallel for
	//vector<vector<int> > nodes = v;
	do{
		n1 = gen()%n;
		}while(!connected[n1]);
	n2 = v[n1][gen()%v[n1].size()];

	int sx = n1%L, sy = int(n1/L), tx, ty;

	//infect single bond
	infected.push_back(n1);
	infected.push_back(n2);
	for(int i=0;i<v[n1].size();i++)
		if(v[n1][i] == n2)
		{
			v[n1].erase(v[n1].begin() + i);
			break;
		}
	for(int i=0;i<v[n2].size();i++)
		if(v[n2][i] == n1)
		{
			v[n2].erase(v[n2].begin() + i);
			break;
		}
	while(!infected.empty())
	{
		infected_counter_currently = 0;
		radius_sum = 0;
		current_size = infected.size();
		for(int inf_idx=0;inf_idx<current_size;inf_idx++)
		{
			n1 = infected[inf_idx];

			tx = n1%L;
			ty = int(n1/L);
			int dx = sx - tx;
			int dy = sy - ty;
			if (dx > L/2)
				dx=dx-L;
			else if (dx < -L/2)
				dx=L+dx;
			if (dy > L/2)
			    dy=dy-L;
			else if (dy < -L/2)
				dy = L +dy;
			if(sqrt(dx*dx + dy*dy) > max_r)
				max_r = sqrt(dx*dx + dy*dy);
			radius_sum+= sqrt(dx*dx + dy*dy);
			radius_sum_full+=sqrt(dx*dx + dy*dy);
			radius_sum_full_counter++;
			for(int i=0;i<v[n1].size();i++)
			{
				if(v[v[n1][i]].size()!=0)
				{
					if((double)gen()/gen.max()<beta)
					{
						infected_counter_currently++;
						if(!(find(infected.begin(), infected.end(), v[n1][i]) != infected.end())) {
							infected.push_back(v[n1][i]);
						}
					}
				}
			}
			v[n1].clear();
		}
		radius_vec.push_back((double)radius_sum/current_size);
		full_radius_vec.push_back((double)radius_sum_full/radius_sum_full_counter);
		infected.erase (infected.begin(),infected.begin()+current_size);

		//branching_factor.push_back((double)infected_counter_currently/infected_counter_past);
		branching_factor.push_back((double)infected_counter_currently);
		//cout<<infected_counter_currently<<endl;
		infected_counter_past = infected_counter_currently;
		t++;
	}
	vector<vector<double> > results;
	vector<double> time;
	time.push_back((double)t);
	time.push_back(max_r);
	results.push_back(radius_vec);
	results.push_back(full_radius_vec);
	results.push_back(time);
	results.push_back(branching_factor);
	return results;
}
int Graph::SIR_multiply(double beta, int source_num){
	vector<int> infected;
	int current_size;
	int t=0;
	int n1, n2;
	int counter=0;
	//omp_set_num_threads(threads);
	//#pragma omp parallel for
	vector<vector<int> > nodes = v;

	//infect cource_num bonds
	while(counter<source_num){
		n1 = gen()%n;
		if(connected[n1]){
			n2 = nodes[n1][gen()%nodes[n1].size()];
			bool flag = find(infected.begin(), infected.end(), n1) == infected.end();
			if(flag){
				infected.push_back(n1);
				if(find(infected.begin(), infected.end(), n2) == infected.end()){
					infected.push_back(n2);
				}
				counter++;
			}
			else{
				if(find(infected.begin(), infected.end(), n2) == infected.end()){
					infected.push_back(n2);
					counter++;
				}
			}
		}
	}

	while(!infected.empty())
	{
		current_size = infected.size();
		for(int inf_idx=0;inf_idx<current_size;inf_idx++)
		{
			n1 = infected[inf_idx];
			for(int i=0;i<nodes[n1].size();i++)
			{
				if(nodes[nodes[n1][i]].size()!=0)
				{
					if((double)gen()/gen.max()<beta)
					{
							if((find(infected.begin(), infected.end(), nodes[n1][i]) == infected.end())) {
								infected.push_back(nodes[n1][i]);
							}
					}
				}
			}
			nodes[n1].clear();
		}
		infected.erase (infected.begin(),infected.begin()+current_size);
		t++;
	}


	return t;

}
void Graph::createGraphSLZcost(float k, float zeta, long double length_sum){

	bool out_of_boundaries;
	double y,x;
	int i=0,j,idx,source,target,c, line=sqrt(n),Px,Py;
	vector<PitStruct>::iterator low;
	std::uniform_real_distribution<double> randreal(0,1);
	randreal(gen);
	long double sum = 0;
	while(sum < length_sum)
	{
		out_of_boundaries=true;
		j=1;
		//cout<<"edge num: "<<i<<endl;

		source = gen()%n;
		//source = rand()%n;
		//y = ((double) rand() / (RAND_MAX));

		y=randreal(gen);
		x = -zeta*log(1-y);
		x*=x;
		if(x>=n){
			//cout<<"system too small!!!    x^2="<<x<<"n= "<<n<<endl;
			continue;
		}
		//idx = binarySearch(0,pitVec.size()-1,x);
		low =lower_bound(pitVec.begin(),pitVec.end(),x, MyCompareStruct());
		idx=low-pitVec.begin();
		if((pitVec[(idx+1)].length-x) < (x - (pitVec[(idx)].length)))
		{
			idx++;
			while(pitVec[(idx)].length==pitVec[(idx+j)].length)
				j++;
			j--;
			if(j!=0)
				{
					//j=rand()%(j+1);
					j=gen()%(j+1);
					idx+=j;
				}
		}
		else if(idx==0)
			j=0;
		else
		{
			while(pitVec[(idx)].length==pitVec[(idx-j)].length)
				j++;
			j--;
			if(j!=0)
			{
				//j=rand()%(j+1);
				j=gen()%(j+1);
				idx-=j;
			}
		}

		//cout<<"i= "<<pitVec->at(idx).at(0)<<"  j= "<<pitVec->at(idx).at(1)<<"  r^2= "<<pitVec->at(idx).at(2)<<" x^2= "<<x<<endl;
        //c = rand()%2;
		c = gen()%2;
        if(c==1)
        {
        	Px = pitVec[(idx)].i;
        	Py = pitVec[(idx)].j;
        }
        else
        {
        	Py = pitVec[(idx)].i;
        	Px = pitVec[(idx)].j;
        }
        //c = rand()%2;
        c = gen()%2;
        if(c==1)
        	Py*=-1;
        //c = rand()%2;
        c = gen()%2;
        if(c==1)
        	Px*=-1;
        /*cout<<i<<". Px= "<<Px<<" Py= "<<Py<<" source= "<<source;*/
        if((source%line)+Px>=line)
        {
        	out_of_boundaries = false;
        	target = ((source/line)*line) + (((source%line)+Px)%line);
        }
        else if(((source%line)+Px)<0)
        {
        	out_of_boundaries =false;
        	if(source<line)
        		target = line + source +Px;
        	else
        		target = ((source/line)*line) + ((source + Px)%line);
        }
        else
        	target = source+Px;
        /*cout<<" targetx= "<<target<<endl;*/
        if(((target/line) + Py )>=line)
        {
        	out_of_boundaries = false;
        	target = ((((target/line) + Py)%line)*line) + target%line;
        }
        else if(((target/line) + Py)<0)
        {
        	out_of_boundaries=false;
        	target = (((target/line) + Py + line)*line) + target%line;
        }
        else
        	target = target+ Py*line;
        if(!(periodic|| out_of_boundaries))
        {
        	continue;

        }
        /*cout<<" targety= "<<target<<endl;
        	cout<<"entered source="<<source<<"  target="<<target<<"degree= "<<getDegree(source)<<endl;*/
        if(source!=target && !isConnected(source,target))
        {
        	sum +=sqrt(pitVec[(idx)].length);
        	addEdge(source,target);
        	i++;
        }
	}
}
void Graph::SIR_simulation(double beta){
	vector<double> radius_vec, full_radius_vec;
	vector<int> infected;
	vector<double> branching_factor;
	int current_size;
	int t=0;
	int n1, n2;
	long double radius_sum, radius_sum_full = 0;
	double max_r=0;
	int radius_sum_full_counter = 0;
	int L = int(sqrt(n));
	int infected_counter_currently, infected_counter_past = 1;

	string fileR = "/home/benaya/Desktop/pics_epidemic/SIR.txt";
	const char *cfile = fileR.c_str();
	ofstream dataInfo;
	dataInfo.open(cfile, ofstream::out);


	//omp_set_num_threads(threads);
	//#pragma omp parallel for
	//vector<vector<int> > nodes = v;
	do{
		n1 = gen()%n;
		}while(!connected[n1]);
	n2 = v[n1][gen()%v[n1].size()];
	int sx = n1%L, sy = int(n1/L), tx, ty;

	//infect single bond
	infected.push_back(n1);
	infected.push_back(n2);
	dataInfo<<-1<<endl;

	for(int i=0;i<v[n1].size();i++)
		if(v[n1][i] == n2)
		{
			v[n1].erase(v[n1].begin() + i);
			break;
		}
	for(int i=0;i<v[n2].size();i++)
		if(v[n2][i] == n1)
		{
			v[n2].erase(v[n2].begin() + i);
			break;
		}
	t++;
	while(!infected.empty())
	{

		infected_counter_currently = 0;
		radius_sum = 0;
		current_size = infected.size();
		for(int ee = 0;ee<current_size;ee++){
			dataInfo<<infected[ee]<<endl;
		}
		dataInfo<<-1<<endl;

		for(int inf_idx=0;inf_idx<current_size;inf_idx++)
		{
			n1 = infected[inf_idx];

			tx = n1%L;
			ty = int(n1/L);
			int dx = sx - tx;
			int dy = sy - ty;
			if (dx > L/2)
				dx=dx-L;
			else if (dx < -L/2)
				dx=L+dx;
			if (dy > L/2)
			    dy=dy-L;
			else if (dy < -L/2)
				dy = L +dy;
			if(sqrt(dx*dx + dy*dy) > max_r)
				max_r = sqrt(dx*dx + dy*dy);
			radius_sum+= sqrt(dx*dx + dy*dy);
			radius_sum_full+=sqrt(dx*dx + dy*dy);
			radius_sum_full_counter++;
			for(int i=0;i<v[n1].size();i++)
			{
				if(v[v[n1][i]].size()!=0)
				{
					if((double)gen()/gen.max()<beta)
					{
						infected_counter_currently++;
						if(!(find(infected.begin(), infected.end(), v[n1][i]) != infected.end())) {
							infected.push_back(v[n1][i]);
						}
					}
				}
			}
			v[n1].clear();
		}
		radius_vec.push_back((double)radius_sum/current_size);
		full_radius_vec.push_back((double)radius_sum_full/radius_sum_full_counter);
		infected.erase (infected.begin(),infected.begin()+current_size);

		branching_factor.push_back((double)infected_counter_currently/infected_counter_past);
		infected_counter_past = infected_counter_currently;

		t++;
		cout<<"t="<<t<<endl;
	}
	dataInfo.close();

}
void Graph::Forest_fire(double beta){

	//there is no need to return the variable t since it is simply the size of r_t and I_t
	int t = 0;
	int n1, n2;

	int L = sqrt(n); // system length

	int sx,sy; //source coordinates
	int tx,ty; //target coordinates
	int dx, dy;

	vector<bool> infected = vector<bool>(n,0);
	vector<bool> recovered = vector<bool>(n,0);

	queue<int> Q; // queue of currently infected nodes

	// immunate each node with probability beta

	for(int i=0;i<n;i++){
		if((double)gen()/gen.max()<beta){
			recovered[i] = 1;
		}
	}

	//infect single node

	do{
		n1 = gen()%n;
	}while(!connected[n1] or recovered[n1] == 1);

	Q.push(n1);
	sx = n1%L;
	sy = int(n1/L);


	int currently_infected = 1; // how many nodes are infected currently,
								//this variable updates itself each time step

	I_t.push_back(1); //single node is infected at t=0
	r_t.push_back(0); // average <r> at t=0 is 0

	double R_sum = 0; // sum of all distances from source at time t

	while(!Q.empty()){ //main loop


		if(currently_infected == 0){ // advance to t+1

			currently_infected = Q.size();
			I_t.push_back(currently_infected);
			r_t.push_back((double)(R_sum/currently_infected));
			R_sum = 0;
			t++;
		}

		n1 = Q.front();
		Q.pop();
		infected[n1] = 0;
		recovered[n1] = 1;
		currently_infected--; //one less...

		for(int i=0;i<v[n1].size();i++){

			n2 = v[n1][i];
			if(recovered[n2] == 0 and infected[n2] == 0){
				Q.push(n2);
				infected[n2] = 1;

				tx = n2%L;
				ty = int(n2/L);

				dx = sx - tx;
				dy = sy - ty;
				R_sum+= sqrt(dx*dx + dy*dy);
			}

		}

	}
}
void Graph::SIR_node_perc(double beta){

	//there is no need to return the variable t since it is simply the size of r_t and I_t
	int t = 0;
	int n1, n2;

	int L = sqrt(n); // system length

	int sx,sy; //source coordinates
	int tx,ty; //target coordinates
	int dx, dy;

	vector<bool> infected = vector<bool>(n,0);
	vector<bool> recovered = vector<bool>(n,0);

	queue<int> Q; // queue of currently infected nodes

	//infect single node

	do{
		n1 = gen()%n;
	}while(!connected[n1]);

	Q.push(n1);
	sx = n1%L;
	sy = int(n1/L);


	int currently_infected = 1; // how many nodes are infected currently,
								//this variable updates itself each time step

	I_t.push_back(1); //single node is infected at t=0
	r_t.push_back(0); // average <r> at t=0 is 0

	double R_sum = 0; // sum of all distances from source at time t

	while(!Q.empty()){ //main loop


		if(currently_infected == 0){ // advance to t+1

			currently_infected = Q.size();
			I_t.push_back(currently_infected);
			r_t.push_back((double)(R_sum/currently_infected));
			R_sum = 0;
			t++;
		}

		n1 = Q.front();
		Q.pop();
		infected[n1] = 0;
		recovered[n1] = 1;
		currently_infected--; //one less...

		for(int i=0;i<v[n1].size();i++){

			n2 = v[n1][i];
			if(recovered[n2] == 0 and infected[n2] == 0){

				if((double)gen()/gen.max()<beta){
					Q.push(n2);
					infected[n2] = 1;

					tx = n2%L;
					ty = int(n2/L);

					dx = sx - tx;
					dy = sy - ty;
					R_sum+= sqrt(dx*dx + dy*dy);
				}
			}

		}

	}
}

void Graph::BCCUtil(int node, int disc[], int low[], vector<vector<int> > *cluster_id, vector<pair<int,int> > *st,vector<bool> *arti,
					int parent[], int *count, int *largest_cluster_id, int *largest_cluster_size)
{
	// A static variable is used for simplicity, we can avoid use
	// of static variable by passing a pointer.
	int time = 0;

	// Initialize discovery time and low value
	disc[node] = low[node] = ++time;

	// Go through all vertices adjacent to this
	vector<int> nodes;
	vector<int> indexs;
	nodes.push_back(node);
	indexs.push_back(0);

	while(!nodes.empty()){
		bool new_node_flag = false;
		int u = nodes.back();
		int children = 0;

		for (int i = indexs.back(); i < v[u].size(); i++)
		{
			int w = v[u][i]; // w is current adjacent of 'u'

			// If v is not visited yet, then recur for it
			if (disc[w] == -1)
			{
				children++;
				parent[w] = u;
				//store the edge in stack
				st->push_back(make_pair(u,w));
				nodes.push_back(w);
				indexs[indexs.size()-1]++;
				indexs.push_back(0);
				disc[w] = low[w] = ++time;
				new_node_flag = true;
				break;
			}

			// Update low value of 'u' only of 'v' is still in stack
			// (i.e. it's a back edge, not cross edge).
			// Case 2 -- per Strongly Connected Components Article
			else if(w != parent[u] && disc[w] < low[u])
			{
				low[u] = min(low[u], disc[w]);
				st->push_back(make_pair(u,w));
			}
		}

		if(new_node_flag)
			continue;

		nodes.pop_back();
		indexs.pop_back();
		low[parent[u]] = min(low[parent[u]], low[u]);


		if( (disc[parent[u]] == 1 && children > 1) ||
			(disc[parent[u]] > 1 && low[u] >= disc[parent[u]]) )
		{
			(*arti)[parent[u]] = 1;
			//(*articulations).push_back(parent[u]);
			//(*articulations)[parent[u]].push_back(*count);
			//cout<<"node "<<parent[u]<<" is arti"<<endl;
			int size = 0;
			while(st->back().first != parent[u] || st->back().second != u)
			{
				//cout << st->back().first << "--" << st->back().second << " ";
				st->pop_back();
				if((*arti)[st->back().first] or (!(*arti)[st->back().first] and (*cluster_id)[st->back().first].size() == 0)){
					(*cluster_id)[st->back().first].push_back(*count);
					size++;
				}
				if((*arti)[st->back().second] or (!(*arti)[st->back().second] and (*cluster_id)[st->back().second].size() == 0)){
					(*cluster_id)[st->back().second].push_back(*count);
					size++;
				}
			}
			//cout << st->back().first << "--" << st->back().second;
			if((*arti)[st->back().first] or (!(*arti)[st->back().first] and (*cluster_id)[st->back().first].size() == 0)){
				(*cluster_id)[st->back().first].push_back(*count);
				size++;
			}
			if((*arti)[st->back().second] or (!(*arti)[st->back().second] and (*cluster_id)[st->back().second].size() == 0)){
				(*cluster_id)[st->back().second].push_back(*count);
				size++;
			}
			st->pop_back();
			//cout << endl;
			if(size > *largest_cluster_size){
				*largest_cluster_size = size;
				*largest_cluster_id = *count;
			}

			(*count)++;
		}
	}
}
void Graph::BCC()
{
	int count = 0;
	vector<bool> *arti = new vector<bool>(n,0);
	int *disc = new int[n];
	int *low = new int[n];
	int *parent = new int[n];
	vector<pair<int,int> >  *st = new vector<pair<int,int> >;
	vector<vector<int> > *cluster_id = new vector<vector<int>>;
	for(int i=0;i<n;i++){
		vector<int> temp;
		(*cluster_id).push_back(temp);
	}
	int largest_cluster_id = -1;
	int largest_cluster_size = 0;
	// Initialize disc and low, and parent arrays
	for (int i = 0; i < n; i++)
	{
		disc[i] = -1;
		low[i] = -1;
		parent[i] = -1;
	}

	for (int i = 0; i < n; i++)
	{
		if (disc[i] == -1)
			BCCUtil(i, disc, low, cluster_id,st, arti,parent,&count,&largest_cluster_id,&largest_cluster_size);

		int j = 0;
		//If stack is not empty, pop all edges from stack
		int size = 0;
		while(st->size() > 0)
		{
			j = 1;
			//cout << st->back().first << "--" << st->back().second << " ";
			if((*arti)[st->back().first] or (!(*arti)[st->back().first] and (*cluster_id)[st->back().first].size() == 0)){
				(*cluster_id)[st->back().first].push_back(count);
				size++;
			}
			if((*arti)[st->back().second] or (!(*arti)[st->back().second] and (*cluster_id)[st->back().second].size() == 0)){
				(*cluster_id)[st->back().second].push_back(count);
				size++;
			}

			st->pop_back();
		}
		if(j == 1)
		{
			//cout << endl;
			if(size > largest_cluster_size){
				largest_cluster_size = size;
				largest_cluster_id = count;
			}
			count++;
		}
	}

/*	cout<<"largest_cluster_id="<<largest_cluster_id<<endl;
	cout<<"largest_cluster_size="<<largest_cluster_size<<endl;
	for(int i=0;i<count;i++){
		int size = 0;
		for(int k=0;k<n;k++){
			if(cluster_id[k] == i){
				cout<<k<<",";
				size++;
			}

			else{
				if((*articulations)[k] == 1)
					for(int w=0;w<v[k].size();w++){
						if(cluster_id[v[k][w]] == i){
							if((*articulations)[v[k][w]] == 0){
								size++;
								cout<<k<<",";
							}
							break;
						}
					}
				}
			}

		cout<<endl;
		cout<<"cluser id = "<<i<<", cluster size="<<size<<endl;

	}
	cout<<"The articulations points are: ";
	for(int q=0;q<n;q++){
		if((*articulations)[q] == 1)
			cout<<q<<",";
	}
	cout<<endl;
	cout << "Above are " << count << " biconnected components in graph";*/


/*	for(int i=0;i<n;i++){
		cout<<"node "<<i<<endl;
		for(int j=0;j<(*cluster_id)[i].size();j++)
			cout<<(*cluster_id)[i][j]<<" ";
		cout<<endl;
	}*/

	cout<<"largest_cluster_id="<<largest_cluster_id<<endl;
	cout<<"largest_cluster_size="<<largest_cluster_size<<endl;
	int new_size = 0;
	for(int i=0;i<n;i++){
		if((*cluster_id)[i].size() == 1 and (*cluster_id)[i][0] == largest_cluster_id){
			new_size++;
			continue;
		}
		if((*cluster_id)[i].size() > 1){
			//cout<<"node: "<<i<<" is arti"<<endl;
			bool flag = true;
			for(int j=0;j<(*cluster_id)[i].size();j++){
				//cout<<(*cluster_id)[i][j]<<endl;
				if((*cluster_id)[i][j] == largest_cluster_id)
						flag = false;
			}
			if(flag){
				active[i] = 0;
				//cout<<"kill node "<<i<<endl;
			}
			else
				new_size++;
			continue;
		}


		active[i] = 0;
		//cout<<"kill node "<<i<<endl;
	}
	for(int i=0;i<n;i++){
		if(active[i] == 0)
			continue;
		int temp = 0;
		for(int j=0;j<v[i].size();j++){
			if(active[v[i][j]] == 1)
				temp++;
		}
		if(temp<2){
			cout<<"kill one"<<endl;
			active[i] = 0;
		}
	}
	cout<<"new_size="<<new_size<<endl;




/*		else{
			if((*articulations)[i] == 0)
				active[i] = 0;
			else{
				bool flag = true;
				for(int w=0;w<v[i].size();w++){
					if(cluster_id[v[i][w]] == largest_cluster_id){
						if((*articulations)[v[i][w]] == 0){
							flag = false;
						}
						break;
					}
				}

				if(flag)
					active[i] = 0;
			}
		}*/


	ofstream dataInfo;
	stringstream ss;
	const char* cfile;
	string fileR = "/storage/home/benaya/Resistors/zeta_interdependent/ER_example/k";
	ss<<3;
	fileR = fileR + ss.str() + ".txt";
	cfile = fileR.c_str();
	dataInfo.open(cfile, ofstream::out);
	for(int i=0; i<n;i++)
	{
		dataInfo<<i<<" ";
		if(active[i] and (*cluster_id)[i].size()>0)
			dataInfo<<-2<<" ";
		else if(active[i] and (*cluster_id)[i].size() == 0)
			dataInfo<<-1<<" ";
		for(int j=0;j<v[i].size();j++)
			dataInfo<<v[i][j]<<" ";
		dataInfo<<endl;
	}

	dataInfo.close();

/*	ofstream dataInfo;
	stringstream ss;
	const char* cfile;
	string fileR = "/storage/home/benaya/Resistors/zeta_interdependent/ER_example/k";
	ss<<4;
	fileR = fileR + ss.str() + ".txt";
	cfile = fileR.c_str();
	dataInfo.open(cfile, ofstream::out);
	for(int i=0; i<p_vec.size();i++)
	{
		dataInfo<<p_vec[i]<<"	"<<backbone_vec[i]<<"	"<<conductivity_vec[i]<<endl;
	}

	dataInfo.close();*/

	for(int i=0;i<n;i++){
		if(active[i] == 0)
			continue;

		int temp = 0;
		for(int j=0;j<v[i].size();j++){
			if(active[v[i][j]] == 1)
				temp++;
		}
		if(temp < 2){
			cout<<"bad node: "<<i<<endl;
			if((*cluster_id)[i].size() > 1){
				cout<<"arti point"<<endl;
				cout<<"cluster size"<<(*cluster_id)[i].size()<<endl;
			}
			else
				cout<<"not arti point"<<endl;
			cout<<"cluster id="<<(*cluster_id)[i][0]<<endl;
			for(int j=0;j<v[i].size();j++){
				if(active[v[i][j]] == 1){
					cout<<"negibor "<<v[i][j]<<endl;
					cout<<"neigbor of neigbor:"<<endl;
					for(int q=0;q<v[v[i][j]].size();q++){
						if(active[v[v[i][j]][q]] == 1)
							cout<<"negibor "<<v[v[i][j]][q]<<endl;
					}
				}

			}
		}
	}
}
void Graph::BCCUtil_rec(int u, int disc[], int low[], vector<pair<int,int> > *st,int parent[])
{
    // A static variable is used for simplicity, we can avoid use
    // of static variable by passing a pointer.
    static int time = 0;

    // Initialize discovery time and low value
    disc[u] = low[u] = ++time;
    int children = 0;

    // Go through all vertices adjacent to this

    for (int i=0; i<v[u].size(); ++i)
    {
        int q = v[u][i];  // v is current adjacent of 'u'

        // If v is not visited yet, then recur for it
        if (disc[q] == -1)
        {
            children++;
            parent[q] = u;
            //store the edge in stack
            st->push_back(make_pair(u,q));
            BCCUtil_rec(q, disc, low, st, parent);

            // Check if the subtree rooted with 'v' has a
            // connection to one of the ancestors of 'u'
            // Case 1 -- per Strongly Connected Components Article
            low[u]  = min(low[u], low[q]);

            //If u is an articulation point,
            //pop all edges from stack till u -- v
            arti[u].push_back(count);
            if( (disc[u] == 1 && children > 1) ||
                (disc[u] > 1 && low[q] >= disc[u]) )
            {
            	int size = 0;
                while(st->back().first != u || st->back().second != q)
                {
                    cout << st->back().first << "--" << st->back().second << " ";
                    if(cluster_id[st->back().first] == -1){
                    	cluster_id[st->back().first] = count;
                    	size++;
                    }
                    if(cluster_id[st->back().second] == -1){
                    	cluster_id[st->back().second] = count;
                    	size++;
                    }
                    st->pop_back();
                }
                cout << st->back().first << "--" << st->back().second;
                if(cluster_id[st->back().first] == -1){
                	cluster_id[st->back().first] = count;
                	size++;
                }
                if(cluster_id[st->back().second] == -1){
                	cluster_id[st->back().second] = count;
                	size++;
                }
                st->pop_back();
                cout << endl;
    			if(size > largest_cluster_size2){
    				largest_cluster_size2 = size;
    				largest_cluster_id2 = count;
    			}
                count++;
            }
        }

        // Update low value of 'u' only of 'v' is still in stack
        // (i.e. it's a back edge, not cross edge).
        // Case 2 -- per Strongly Connected Components Article
        else if(q != parent[u] && disc[q] < low[u])
        {
            low[u]  = min(low[u], disc[q]);
            st->push_back(make_pair(u,q));
        }
    }
}

// The function to do DFS traversal. It uses BCCUtil()
void Graph::BCC_rec()
{
	for(int i=0;i<n;i++){
		vector<int> temp;
		arti.push_back(temp);
		cluster_id.push_back(-1);
	}

    int *disc = new int[n];
    int *low = new int[n];
    int *parent = new int[n];
    vector<pair<int,int> > *st = new vector<pair<int,int> >[n];

    // Initialize disc and low, and parent arrays
    for (int i = 0; i < n; i++)
    {
        disc[i] = -1;
        low[i] = -1;
        parent[i] = -1;
    }

    for (int i = 0; i < n; i++)
    {
        if (disc[i] == -1)
            BCCUtil_rec(i, disc, low, st, parent);

        int j = 0;
        int size = 0;
        //If stack is not empty, pop all edges from stack
        while(st->size() > 0)
        {
            j = 1;
            cout << st->back().first << "--" << st->back().second << " ";
            if(cluster_id[st->back().first] == -1){
            	cluster_id[st->back().first] = count;
            	size++;
            }
            if(cluster_id[st->back().second] == -1){
            	cluster_id[st->back().second] = count;
            	size++;
            }
            st->pop_back();
        }
        if(j == 1)
        {
            cout << endl;
			if(size > largest_cluster_size2){
				largest_cluster_size2 = size;
				largest_cluster_id2 = count;
			}
            count++;
        }
    }

	cout<<"largest_cluster_id="<<largest_cluster_id2<<endl;
	cout<<"largest_cluster_size="<<largest_cluster_size2<<endl;
	for(int i=0;i<n;i++){
		if(cluster_id[i] == largest_cluster_id2)
			continue;
		if(find(arti[i].begin(),arti[i].end(),largest_cluster_id2) != arti[i].end()){
			cout<<"node: "<<i<<endl;
			cout<<"clusters: "<<endl;
			for(int j=0;j<arti[i].size();j++)
				cout<<arti[i][j]<<endl;
			continue;
		}
		active[i] = 0;
	}

	ofstream dataInfo;
	stringstream ss;
	const char* cfile;
	string fileR = "/storage/home/benaya/Resistors/zeta_interdependent/ER_example/k";
	ss<<3;
	fileR = fileR + ss.str() + ".txt";
	cfile = fileR.c_str();
	dataInfo.open(cfile, ofstream::out);
	for(int i=0; i<n;i++)
	{
		dataInfo<<i<<" ";
		if(active[i] and arti[i].size()>0)
			dataInfo<<-2<<" ";
		else if(active[i] and arti[i].size() == 0)
			dataInfo<<-1<<" ";
		for(int j=0;j<v[i].size();j++)
			dataInfo<<v[i][j]<<" ";
		dataInfo<<endl;
	}

	dataInfo.close();
}
#endif
