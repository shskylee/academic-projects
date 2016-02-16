// parallel monte carlo simulation for sis model on hierarhical modular network in Munoz paper

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

int    s;             // the largest hierarhical level
int    M0;            // the number of nodes in the lowest level module
double p;             // the probability for a link existing between the the lowest level modules
int    alpha;         // expected number of links between the lowest level modules
double lambda;        //spreading rate
int    N;
double time_max=1.0e+5;
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

vector<Node<int>*>  hmn(){
    
    
    vector<pair<int, int> > Gedges;
    
    srand(time(NULL));
    
    for (int l=0;l<=s; ++l){
        if (l==0){
            int m=M0*pow(2,l);
            int n=N/m;
            for (int i=0; i<n;++i) {
                for (int j=0; j<m;++j) {
                    for (int k=j+1;k<m;++k) {
                        Gedges.push_back(make_pair(i*m+j,i*m+k));
                    }}}}
        else{
            int m=M0*pow(2,l-1);
            int n=N/(2*m);
            for (int i=0; i<n; i++) {
                int j=2*i+1;
                int count=0;
                while (count<1) {
                    for (int h=0; h<m; h++) {
                        for (int k=0; k<m; k++) {
                            double rnd=randomG();
                            if (rnd<=alpha*pow(p,l)) {
                                Gedges.push_back(make_pair((2*i)*m+h,j*m+k));
                                ++count;}}}}}}
    }
    
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

 
    s=atoi(argv[1]);
    lambda=atof(argv[2]);
    M0=2;
    p=1.0/4.0;
    alpha=1;
    N=M0*pow(2,s);
    
    
    //MPI initializes
    MPI_Init(&argc, &argv);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    double start, end;
    

    MPI_Barrier(MPI_COMM_WORLD); /* IMPORTANT */
    start = MPI_Wtime();
    
    int length=floor(time_max/t_delta);
    double* p_mean= new double[length];
    double* times=new double[length];
    
  
    adjacency=hmn();
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
    
    int count=0,cnt=0;
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
        }else{
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
        
        if (floor(t/t_delta)==count) {
            double average=(1.0*actives)/N;
            times[count]  += t;
            p_mean[count] += average;
            cnt += 1;
        }else{
            times[count] /= cnt;
            p_mean[count] /=cnt;
            cnt=1;
            count += 1;
            double average=(1.0*actives)/N;
            times[count]  += t;
            p_mean[count] += average;
            
        }
        
        if (t>time_max) {
            break;
        }
    }
    
    
    int* elements_procs=new int[world_size];
    
    MPI_Allgather(&count, 1, MPI_INT, elements_procs, 1, MPI_INT, MPI_COMM_WORLD);
    int elements_min=*min_element(elements_procs,elements_procs+world_size);
    int elements_max=*max_element(elements_procs,elements_procs+world_size);

    
    int     extra_num=elements_max-elements_min;
    int*    extra_per_proc=new int[extra_num]();
    int*    extra_procs=new int[extra_num];
    double* extra_local_times= new double[extra_num]();
    double* extra_local_p_mean=new double[extra_num]();
    double* extra_global_times= new double[extra_num];
    double* extra_global_p_mean=new double[extra_num];
    
    for (int i=elements_min+1; i<=count; i++) {
        int i0=i-elements_min-1;
        extra_per_proc[i0]=1;
        extra_local_p_mean[i0]=p_mean[i];
        extra_local_times[i0]=times[i];
    }
    
    double* global_times= new double[elements_min+1];
    double* global_p_mean= new double[elements_min+1];
    

    
    
    MPI_Reduce(times,global_times,elements_min+1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(p_mean,global_p_mean,elements_min+1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

    MPI_Allreduce(extra_per_proc,extra_procs,extra_num,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(extra_local_times,extra_global_times,extra_num,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(extra_local_p_mean,extra_global_p_mean,extra_num,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    
    
    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();
    
    if (world_rank==0) {
        
        for (int i=0; i<elements_min+1; i++) {
            global_p_mean[i]/= world_size;
            global_times[i] /= world_size;
        }
        for (int i=0; i<extra_num; i++) {
            
            extra_global_p_mean[i]/= world_size;
            extra_global_times[i] /= extra_procs[i];
        }
        
        
        
        string filename1;
        string filename2;
        filename1=(string)"sis_hmn_s"+argv[1] +"l"+argv[2]+".txt";
        filename2=(string)"time_cost_s"+argv[1] +"l"+argv[2]+".txt";
        ofstream ofs1(filename1.c_str());
        ofstream ofs2(filename2.c_str());
        ofs1<<fixed;
        ofs2<<fixed;
        for (int i=0; i<elements_min+1; i++){
         ofs1<<setprecision(10)<<global_times[i]<<"\t"<<global_p_mean[i]<<"\n";
        }
        for (int i=0; i<extra_num; i++){
         ofs1<<setprecision(10)<<extra_global_times[i]<<"\t"<<extra_global_p_mean[i]<<"\n";
        }
        

        
        
        ofs2<<"elapsed time: "<<setprecision(10)<<end-start<<"seconds"<<"\n";
        
    }
    


    

    

    
    delete[] global_times;
    delete[] global_p_mean;
    delete[] extra_local_times;
    delete[] extra_local_p_mean;
    delete[] extra_global_p_mean;
    delete[] extra_global_times;
    delete[] extra_per_proc;
    delete[] extra_procs;


    MPI_Finalize();

    return 0;
}

