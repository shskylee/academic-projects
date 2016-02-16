// parallel monte carlo simultion for SIS model on random networks


#include<iostream>
#include<vector>
#include<string>
#include<cmath>
#include<cstdlib>
#include<numeric>
#include<time.h>
#include<fstream>
#include<ctime>
#include<queue>
#include<iomanip>
#include<algorithm>
#include<mpi.h>

using namespace std;

double lambda;        //spreading rate
int    N;
int    E;
double time_max=1.0e+6;
double t_delta=1.0;

template<typename T>
struct Node{
    T data;
    Node<T> *next;
    Node& operator=(const Node* a)
    {
        data=a->data;
        next=a->next;
        return this;
    }
};

double randomG(){
    return (double)rand()/(double)RAND_MAX;
}


vector<Node<int>*> adjacency;

vector<Node<int>*> ER(){
    
    srand(time(NULL));
    
    vector<pair<int, int> > edges(E);
    for (int j=0; j<E; j++) {
        double rnd1=randomG();
        double rnd2=randomG();
        int e1=floor(rnd1*N);
        int e2=floor(rnd2*N);
        while(e1==e2)
        {
            // avoiding the self loop
            rnd2=randomG();
            e2=floor(rnd2*N);
        }
        edges[j].first=e1;
        edges[j].second=e2;
    }
    
    vector<Node<int>* > adj(N);
    for (int i=0; i<N; i++) {
        adj[i]=new Node<int>;
        adj[i]->next=NULL;
    }
    
    for(int j=0;j<edges.size();j++){
        int row1=edges[j].first;
        int row2=edges[j].second;
        Node<int>* ptr=adj[row1];
        if (ptr==NULL) {
            ptr->data=edges[j].second;
        }else{
            Node<int>* ptr1;
            ptr1=new Node<int>;
            ptr1->next=NULL;
            ptr1->data=edges[j].second;
            while (ptr->next!=NULL) ptr=ptr->next;
            ptr->next=ptr1;
        }
        
        ptr=adj[row2];
        if (ptr==NULL) {
            ptr->data=edges[j].first;
        }else{
            Node<int>* ptr2;
            ptr2=new Node<int>;
            ptr2->next=NULL;
            ptr2->data=edges[j].first;
            while (ptr->next!=NULL) ptr=ptr->next;
            ptr->next=ptr2;
        }
    }
    
    return adj;
}

int main(int argc, char* argv[]){

 
    N=atoi(argv[1]);
    E=atof(argv[2]);
    lambda=atof(argv[3]);
    
    
    
    //MPI initializes
    MPI_Init(&argc, &argv);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    double start, end;
    

   // MPI_Barrier(MPI_COMM_WORLD); /* IMPORTANT */
   // start = MPI_Wtime();
    
    int length=floor(time_max/t_delta);
    vector<double> p_mean(length,0.0);
    vector<double> times(length,0.0);
     
   
    adjacency=ER();
    vector<int>    states;
    vector<int>    active_nodes;
    vector<int>    links(N,0);
    int active_links=0;
    int actives=N;
    double t=0.0;
    srand(time(NULL)*world_rank);

    //continuous time MC
    
    for (int i=0; i<N; i++) {
            active_nodes.push_back(i);
            states.push_back(1);
    }
        
    for (int i=0; i<actives; i++) {
        Node<int>* ptr=adjacency[i]->next;
        while (ptr!=NULL) {
            links[i]++;
            ptr=ptr->next;}
            active_links += links[i];
            
    }
    
    
    vector<int> cnt(length,0);
    while (actives>1) {
            
        int pos=floor(actives*randomG());
        int nd=active_nodes[pos];
        
        double rnd=randomG();
        t += 1.0/(actives+lambda*active_links);
        
        if (rnd<1.0*actives/(actives+lambda*active_links)) {
            states[nd]=0;
            active_nodes.erase(active_nodes.begin()+pos);
            actives--;
            active_links -= links[nd];
        }
        else{
            vector<int> inactive_neighbors;
            Node<int>*  ptr=adjacency[nd]->next;
            while (ptr!=NULL) {
                if (states[ptr->data]==0){
                    inactive_neighbors.push_back(ptr->data);
                    
                }
            
                ptr=ptr->next;
            }
            int s1=inactive_neighbors.size();
            if (s1>1) {
                int pos1=floor(s1*randomG());
                int nd1=inactive_neighbors[pos1];
                active_nodes.push_back(nd1);
                states[nd1]=1;
                actives++;
                active_links +=links[nd1];
            }
                
            }
     
      if(t>time_max){
             break;}
        
      int count=floor(t/t_delta);
      times[count] +=t;
      p_mean[count]+= 1.0*actives/N;
      cnt[count]++;
        
    }
    
    
    vector<double> global_times(length);
    vector<double> global_p_mean(length);
    vector<int> global_cnt(length);
    MPI_Reduce(&times[0],&global_times[0],length,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&p_mean[0],&global_p_mean[0],length,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&cnt[0],&global_cnt[0],length,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    
    if (world_rank==0) {
        
        for (int i=0; i<length; i++) {
            if (global_cnt[i]>0) {
                global_p_mean[i]/= global_cnt[i];
                global_times[i] /= global_cnt[i];
            }
        }
         
        string filename1;
        string filename2;
        filename1=(string)"sis_hmn_N"+argv[1]+"E"+argv[2] +"l"+argv[3]+".txt";
        filename2=(string)"time_cost_N"+argv[1]+"E"+argv[2] +"l"+argv[3]+".txt";
        ofstream ofs1(filename1.c_str());
        ofstream ofs2(filename2.c_str());
        ofs1<<fixed;
        ofs2<<fixed;
        for (int i=0; i<length; i++){
          if (global_cnt[i]>0) {
            ofs1<<setprecision(10)<<global_times[i]<<"\t"<<global_p_mean[i]<<"\n";
          }}
        ofs2<<"elapsed time: "<<setprecision(10)<<end-start<<"seconds"<<"\n";
        
    }
    
    
    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();
   
    MPI_Finalize();

    return 0;
}

