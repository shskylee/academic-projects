// parallel monte carlo simulation for SIS model on 1D modular networks

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

int       g;
int       cn;                      // number of nodes in one supernode
int       N;
int       num;
double    lambda;                  // spreading rate
double    time_max=1.0e+5;
double    t_b=0.001;     // logrithmic bin parameter x_(i+1)=x_i*exp(t_b)
double    time_min=1.0e-3;
double    t_delta=1.0;


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

vector<Node<int>*> OneD(){
    
    srand(time(NULL));
    vector<pair<int, int> > Gedges;
    
    for(int i=0;i<num;i++){
        for(int j=0;j<cn;j++){
            for(int k=j+1;k<cn;k++){
                Gedges.push_back(make_pair(i*cn+j,i*cn+k));
            }}}
    
    
    for (int i=0; i<num-1; i++) {
        double rnd1=randomG();
        double rnd2=randomG();
        int e1=i*cn+floor(rnd1*cn);
        int e2=(i+1)*cn+floor(rnd2*cn);
        Gedges.push_back(make_pair(e1,e2));
    }
    double rnd1=randomG();
    double rnd2=randomG();
    int    e1=floor(rnd1*cn);
    int    e2=(num-1)*cn+floor(rnd2*cn);
    Gedges.push_back(make_pair(e1,e2));
    
    vector<Node<int>* > adj(N);
    for (int i=0; i<N; i++) {
        adj[i]=new Node<int>;
        adj[i]->next=NULL;
    }
    
    for(int j=0;j<Gedges.size();j++){
        int row1=Gedges[j].first;
        int row2=Gedges[j].second;
        Node<int>* ptr=adj[row1];
        if (ptr==NULL) {
            ptr->data=Gedges[j].second;
        }else{
            Node<int>* ptr1;
            ptr1=new Node<int>;
            ptr1->next=NULL;
            ptr1->data=Gedges[j].second;
            while (ptr->next!=NULL) ptr=ptr->next;
            ptr->next=ptr1;
        }
        
        ptr=adj[row2];
        if (ptr==NULL) {
            ptr->data=Gedges[j].first;
        }else{
            Node<int>* ptr2;
            ptr2=new Node<int>;
            ptr2->next=NULL;
            ptr2->data=Gedges[j].first;
            while (ptr->next!=NULL) ptr=ptr->next;
            ptr->next=ptr2;
        }
    }
    
    return adj;
}


int main(int argc, char* argv[]){
    //MPI initializes
    MPI_Init(&argc, &argv);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    
    g=atoi(argv[1]);
    cn=atoi(argv[2]);
    lambda=atof(argv[3]);
    num=pow(2,g);
    N=num*cn;
    double start, end;
    

    MPI_Barrier(MPI_COMM_WORLD); /* IMPORTANT */
    start = MPI_Wtime();
    
    int length=ceil(log(time_max/time_min)/t_b)+1;
    vector<double> p_mean(length,0.0);
    vector<double> times(length,0.0);
    vector<int> cnt(length,0);
     vector<double> time_axis(length);
    time_axis[0]=time_min;
    for(int i=1;i<length;i++){
        time_axis[i]=time_min*exp(i*t_b);
    }
    
   // cout<<"length is: "<<length<<endl;
    adjacency=OneD();
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
    
    
    
    int count=0;
    
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
        
      if(t>time_axis[count]){
            count++;
            times[count] += t;
            p_mean[count]+= 1.0*actives/N;
            cnt[count]++;
            
        }else{
            times[count] += t;
            p_mean[count]+= 1.0*actives/N;
            cnt[count]++;
        }

    }

   
    vector<double> global_times(length);
    vector<double> global_p_mean(length);
    vector<int> global_cnt(length);
    MPI_Reduce(&times[0],&global_times[0],length,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&p_mean[0],&global_p_mean[0],length,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&cnt[0],&global_cnt[0],length,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
  
    
  

    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();
    
    if (world_rank==0) {
        
        string filename1;
        string filename2;
        filename1=(string)"sis_1D_s"+argv[1]+"cn"+argv[2]+"l"+argv[3]+"run"+argv[4]+".txt";
        filename2=(string)"time_cost_1D_s"+argv[1]+"cn"+argv[2]+"l"+argv[3]+"run"+argv[4]+".txt";
        ofstream ofs1(filename1.c_str());
        ofstream ofs2(filename2.c_str());
        ofs1<<fixed;
        ofs2<<fixed;
        for (int i=1; i<length; i++){
            if (global_cnt[i]>0) {
                ofs1<<setprecision(10)<<global_times[i]/global_cnt[i]<<"\t"<<global_p_mean[i]/global_cnt[i]<<"\n";
            }}
        ofs2<<"elapsed time: "<<setprecision(10)<<end-start<<"seconds"<<"\n";
        
    }
     
     
    
    

    
    MPI_Finalize();

    return 0;
}

