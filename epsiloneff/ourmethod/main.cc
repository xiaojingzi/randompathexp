#include<stdio.h>
#include<stdlib.h>
#include<sys/time.h>
#include<time.h>
#include<iostream>
#include<set>
#include<map>
#include<vector>
#include<algorithm>
#include <cmath>
#include <cstring>
//#include <cstdlib>
using namespace std;

//#define N 1023170
//#define M 3172790
//#define N 112411
//#define M 417869
//#define T 10
//#define episilon 0.005
//int N = 0;
//int M = 0;
//int R = 0;
//#define R (int)(((double)0.5 / (episilon * episilon)) * (log2(T * (T-1)) + 1 + log((double)1/0.1)))

static int N = 112411;
static int M = 417869;
static int T = 10;
static double Epsilon = 0.005;

static int R = -1 ;

inline int calR() {
    return (int)(((double)0.5 / (Epsilon * Epsilon )) * (   (log2( T * (T-1)/2 ) + 1) + log((double)1/0.1)));
}



double timestamp (){
  struct timeval tv;
  gettimeofday (&tv, 0);
  return tv.tv_sec + 1e-6*tv.tv_usec;
}

typedef struct 
{
  int id;
  char name[1024];
  double weight;
  /* data */
}Node,*_Node;


int D;
char dataset[100];
char network_file[100];
char map_file[100];

char similarPathfile[100];
char similarStructurefile[100];


char ** Gnode;
//char Gnode[N][1024];

vector<pair<int,double> > * G ;
vector<pair<int,double> > GWeightedDegree;
vector<pair<int,double> > * GAccumulatedNeighborVector;
//vector<pair<int,double>> *G = new vector<pair<int,double>>[N];

//int Gedges[M*2];
//int Gedgee[M*2];

//int start[M];
//int length[M];

int * Gedges;
int * Gedgee;

int * start;
int * length;


//vector<int> node2path[N];
vector<int> * node2path;
//vector<int> path2node[R];
//vector<int> * path2node = new vector<int>[R];
vector<int> * path2node;


vector<pair<int, int> > * similar_path_list;
vector<int> * similar_structure_list;
// double t1 = 0.0;
// double t2 = 0.0;


void init(char* data, int t, int d, double epsilon){
  T=t;//10
  D=d;//50
  Epsilon = epsilon; 

  int displayEpsilon = round(epsilon*10000);
  printf("epsilon=%lf, displayed epsilonn = %d\n", epsilon, displayEpsilon);
  cout<< epsilon << " "<< displayEpsilon << endl;
  sprintf(dataset, "%s", data);
  sprintf(network_file, "data/%s.graph", dataset);
  sprintf(map_file, "data/%s.dict", dataset);
  sprintf(similarPathfile,"result/%s_%d_%d_%d.pathsim",dataset, T, D, displayEpsilon);
  sprintf(similarStructurefile,"result/%s_%d_%d_%d.pathvec",dataset,T, D, displayEpsilon);
  printf("%s, %s, %s\n", dataset, network_file, map_file);
}

void loadNode(){
  //  FILE *fp = fopen(network_file, "r");
  FILE *fp = fopen(map_file, "r");
  char line[1000];
  int d=0;
  printf("map file : %s \n", map_file);
  if(fp == NULL ) {
	fprintf(stderr, "open file : %s failed !!!! \n", map_file);
	fflush(NULL);
  }

  while(!feof(fp)){
    d++;
   fgets(line,1000, fp);
    //printf("%s\n", line);
    int a;
    char tem[1024];
    char * k = strtok(line,"\t");
    memcpy(tem,k,1023);
    if( k == NULL) {
       fprintf(stderr,"unknown val : %s\n",line);
       continue;
    }
    k =  strtok(NULL,"\t");
    if( k == NULL) {
       fprintf(stderr,"unknown val : %s\n",line);
       continue;
    }
    a = atoi(k);
    //printf("%s -- %d \n" , tem,a);
    if(a > N ) {
    	fprintf(stderr,"a val : %d , N val : %d \n", a , N);
    } else {
    	memcpy(Gnode[a],tem,1023);
    }
    //printf("val : %s %d \n",Gnode[a],a); 
   //Gnode[a]=tem;
    if(d % 10000 == 0 ) {
 	printf("load : %d ... \n", d);
	fflush(NULL);
    }
//       printf("%d %s\n",a ,tem);
//       if(d>20){ printf("wan\n"); exit(0); }
  }
  printf("Load nodes done.\n");
  fclose(fp);
}

void readMandN(char * line , int &n, int &m) {
    char * k = strtok(line,"\t");
    if( k == NULL) {
       fprintf(stderr,"unknown val : %s\n",line);
       return;
    }
    n = atoi(k); 
    k =  strtok(NULL,"\t");
    if( k == NULL) {
       fprintf(stderr,"unknown val : %s\n",line);
        return ;
    }
    m = atoi(k);
    printf("read N: %d M: %d \n",n,m);
}


// N , M first line ,Node & Edge number
void loadLink(){
  FILE *fp = fopen(network_file, "r");
  char line[100];
  fgets(line,1000, fp);
  readMandN(line,N,M);
  R = calR();
  printf("R=%d\n", R);
  Gnode = new char*[N];
  for(int i = 0 ; i < N ; i++ )
      Gnode[i] = new char[1024];
 
  G = new vector<pair<int,double> >[N];
  GAccumulatedNeighborVector =  new vector<pair<int,double> >[N];

  Gedges = new int[M * 2];
  Gedgee = new int[M * 2];

  start = new int[M];
  length = new int[M];

  node2path = new vector<int>[N];
  path2node = new vector<int>[R];

  similar_path_list = new vector<pair<int,int> >[N];
  similar_structure_list = new vector<int>[N];


  for (int i = 0 ;i < N ; i++){
    GWeightedDegree.push_back(make_pair(i,0));
  }
  while(!feof(fp)){
    int a, b;
    double c;
    fgets(line,100, fp);
    sscanf(line, "%d\t%d\t%lf", &a, &b,&c);
    G[a].push_back(make_pair(b,c));
    
    G[b].push_back(make_pair(a,c));
    //if(a==b)printf("a=%d, b=%d\n",a,b);
  }
  printf("load edges done.\n");
}

//binary search
int BinarySearch(vector<pair<int, double> > array, int value)  
{  

    int low = 0;  
    int high = array.size() - 1;  
    while (low <= high)  
    {  
        int mid = low + (high - low) / 2;  
        if (array[mid].second == value)  
            return array[mid].first;  
        else if (array[mid].second > value)  
            high = mid - 1;  
        else  
            low = mid + 1;  
    } 
    if (high < 0 )
      high = 0;
    if (high >  array.size()-1 ) 
      high = array.size() -1;
    return array[high].first;  
}  

struct CmpByValue{
  bool operator()(const pair<int, int>&lhs, const pair<int, int>& rhs){
    return lhs.second < rhs.second;
  }
};

void generateRandomPath(){

  printf("sample size = %d\n", R);
  // double t1 = timestamp();

 
  sort(GWeightedDegree.begin(), GWeightedDegree.end(), CmpByValue());
  
  vector<pair<int, double> > accumated_degee;
  accumated_degee.push_back(make_pair(GWeightedDegree[0].first, GWeightedDegree[0].second));
  for(unsigned int i = 1 ; i < N ; i++) {
    accumated_degee.push_back( make_pair(GWeightedDegree[i].first, GWeightedDegree[i].second + accumated_degee[i-1].second) ) ;
  } 

  for(unsigned int i = 0 ; i < N ; i++) {
    for(vector<pair<int,double> >::const_iterator it = G[i].begin(); it != G[i].end() ; it ++) {
      GWeightedDegree[i].second += it->second;
    }
  }

  for(unsigned int i = 0 ; i < N ; i++) {
    int neighborlen = G[i].size();
    if (neighborlen > 0){
      // sort(G[i].begin(), G[i].end(), CmpByValue());
      GAccumulatedNeighborVector[i].push_back(make_pair(G[i][0].first, G[i][0].second/ GWeightedDegree[i].second));

      if(neighborlen > 1){
        for(unsigned int j = 1 ; j < neighborlen ; j++) {
          GAccumulatedNeighborVector[i].push_back( make_pair(G[i][j].first, G[i][j].second / GWeightedDegree[i].second+ GAccumulatedNeighborVector[i][j-1].second));
        }
      }
    }
  }

  // double t2 = timestamp();
  // printf("Time of preparing accumulated vectors = %lf s\n", t2-t1);


  // R = 1;
  // int randid = 0;
  
  // t1 = timestamp();
  for(int r=0; r<R; r++){
    //for(int j=0; j<R; j++){
      int current_node = rand() % N ;// 0 - N-1
      // int current_node = randid;
      // select by degree begin ----
      // double rand_degree = (double)rand() / RAND_MAX * accumated_degee[N-1].second;
      // int current_node = BinarySearch(accumated_degee, rand_degree);
      // select by degree end ---
      int t=0;
      //int pathid=i*R+j;
      int pathid=r;
      path2node[pathid].push_back(current_node);
      node2path[current_node].push_back(pathid);
      while(t<T){
        // printf("current node = %d\n", current_node);
        int neighborlen = G[current_node].size();
        if(neighborlen == 0)
          break;
      
        // printf("neighbors\n");
        // for(unsigned int i = 0 ; i < neighborlen ; i++) {
        //     printf("%d ", G[current_node][i].first);
        // }
        int father = current_node;
        
        // double total_weight = 0;
        // for(vector<pair<int,double> >::const_iterator it = G[current_node].begin(); it != G[current_node].end() ; it ++) {
        //     total_weight += it->second;
        // }
        // printf("total weight = %lf\n", total_weight);
        // printf("neighbor number = %d\n", neighborlen);
        // printf("original weight\n");
        // for(unsigned int i = 0 ; i < neighborlen ; i++) {
        //     printf("%lf ", G[current_node][i].second);
        // }
        // printf("\n");
        // vector<double> ws;
        // ws.push_back(G[current_node][0].second / total_weight);
        // for(unsigned int i = 1 ; i < neighborlen ; i++) {
        //     ws.push_back( G[current_node][i].second / total_weight + ws[i-1]);
        // }

        // printf("accumulated weight\n");
        // for(int i = 0 ; i < ws.size(); i++){
        //   printf("%lf ", ws[i]);
        // }
        // printf("\n");

        double rand_val = (double)rand() / RAND_MAX;
        // printf("random value = %lf\n", rand_val);

        // int neighbor_index = 0;
        // for(unsigned int i = 0 ; i < neighborlen ; i++) {
        //     if(ws[i] >= rand_val ) {
        //         neighbor_index = i;
        //         break;
        //     }
        // }
        // current_node = G[father][neighbor_index].first;
        
        int neighbor_index = 0;
        for(unsigned int i = 0 ; i < neighborlen ; i++) {
            if(GAccumulatedNeighborVector[current_node][i].second >= rand_val ) {
                neighbor_index = i;
                break;
            }
        }
        current_node = GAccumulatedNeighborVector[current_node][neighbor_index].first;
        
        // printf("neighbor length = %d\n", neighborlen);
        // printf("accumulated neighbor length = %d\n", GAccumulatedNeighborVector[current_node].size());
        
        // double rand_weight = (double) rand() / RAND_MAX  ;//* GAccumulatedNeighborVector[current_node][neighborlen-1].second;
        // int neighbor_index = BinarySearch(GAccumulatedNeighborVector[current_node], rand_weight);
        // current_node = neighbor_index;

        // printf("neighbor node = %d, pathid = %d\n", current_node, pathid);
      
        path2node[pathid].push_back(current_node);
        node2path[current_node].push_back(pathid);//can not parallel
        t++;
      }
      // if((r % 1000) == 0){
      //  t2 = timestamp();
      //  printf("path id = %d, time = %lf s\n", r, (t2 - t1));
      //  t1 = t2;
      // }
    
  }
}





struct CmpByValueReversed{
  bool operator()(const pair<int, int>&lhs, const pair<int, int>& rhs){
    return lhs.second > rhs.second;
  }
};

void calculatePathSim(){
  for(int i=0; i<N; i++){
    map<int, int> nodesInSamePath;
   
    // t1 = timestamp() ;
    vector<int>::iterator iter;
    for(iter = node2path[i].begin(); iter != node2path[i].end(); iter++){
      set<int> nodeset;
      nodeset.insert(path2node[*iter].begin(), path2node[*iter].end());
      set<int>::iterator it;
      for(it = nodeset.begin(); it != nodeset.end(); it++){
        if(!nodesInSamePath.count(*it))
          nodesInSamePath[*it]=0;
        nodesInSamePath[*it]++;
      }
    }
    vector<pair<int, int> > sortedNodes(nodesInSamePath.begin(), nodesInSamePath.end());
    sort(sortedNodes.begin(), sortedNodes.end(), CmpByValueReversed());

    // t2 = timestamp();
    // printf("search top K nodes with path similarity = %lf s\n", (t2 - t1));

    int tempt =sortedNodes.size(); 
    int dim = (D<tempt)? D: tempt;
    
    similar_path_list[i].insert(similar_path_list[i].begin(), sortedNodes.begin(), sortedNodes.begin()+dim);
    for(int j = 0; j<dim; j++){
      int tem = sortedNodes[j].second;
      similar_structure_list[i].push_back(tem);
    }
    if(similar_structure_list[i].size()<D){
      for(int j=dim; j<D; j++)
        similar_structure_list[i].push_back(0);

    }

    //if((i % 10000)==0)
    //  printf("calculate %d\n", i);
  }

}

void savePathSim(){
  //sprintf(similarPathfile, "result/%s/similarPath", dataset);
  FILE *fp = fopen(similarPathfile, "w");
  for(int i=0; i<N; i++){
    vector<pair<int, int> >::iterator iter;
    for(iter = similar_path_list[i].begin(); iter != similar_path_list[i].end(); iter++){
      fprintf(fp, "%d ", iter->first);
    }
    fprintf(fp, "\n");
    //if((i%10000) == 0)
    //  printf("%d\n", i);
  }
  fclose(fp);

}

void saveStructureSim(){
  FILE *fp = fopen(similarStructurefile, "w");
  for(int i=0; i<N; i++){
    vector<int>::iterator iter;
    for(iter = similar_structure_list[i].begin(); iter != similar_structure_list[i].end(); iter++){
      fprintf(fp, "%d ", *iter);
    }
    fprintf(fp, "\n");
    //if((i%10000) == 0)
    //  printf("%d\n", i);
  }
  fclose(fp);
}


void precompute(char* dataset, int T, int D, double Epsilon){
  init(dataset, T, D, Epsilon);
  double t1 = timestamp();
  //loadNode();
  // double t2 = timestamp();
  // printf("Time of loading nodes = %lf s\n", t2-t1);
  loadLink();
  double t3 = timestamp();
  printf("Time of loading edges = %lf s\n", t3-t1);
  generateRandomPath();
  double t4= timestamp();
  printf("Time of generating random paths = %lf s\n", t4-t3);
  calculatePathSim();
  double t5= timestamp();
  printf("Time of calcualating pathsim = %lf s\n", t5-t4);
  savePathSim();
  saveStructureSim();
  double t6= timestamp();
  printf("Time of saving  = %lf s\n", t6-t5);

}

int main(int argc, char* argv[]){
  srand((unsigned long)time(NULL));
  int T=10; // path length
  int D=50;  // vector dimension
  double Epsilon = 0.01;
  char *data = argv[1];
  if(argc > 2)
    T = atoi(argv[2]);
  if(argc > 3)
    D = atoi(argv[3]);
  if(argc > 4)
    Epsilon = atof(argv[4]);
  printf("data: %s, T: %d, D:  %d,  epsilon: %lf\n", data, T, D, Epsilon);
  precompute(data, T, D, Epsilon);

  return 0;
}
