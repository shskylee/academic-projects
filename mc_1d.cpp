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


using namespace std;

double    time_max=1.0e+7;
double    t_b=0.001;     // logrithmic bin parameter x_(i+1)=x_i*exp(t_b)
double    time_min=1.0;

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

int main(int argc, char* argv[]){
  
  
    int N=atoi(argv[1]);
    int E=atoi(argv[2]);
    double lambda=atof(argv[3]);
    

    clock_t start, end;
    double used_time;
    
    start=clock();
    
    srand(time(NULL));
    
    // contruct adjacency matrix
    
    int max_edges=100;
    int** adjacency=new int*[N];
    for (int i=0; i<N; i++) {
        adjacency[i]=new int[max_edges];
    }
    
    int* edges_count=new int[N]();
    for(int j=0;j<E;j++){
        int e1=rand()%N;
        int e2=rand()%N;
        while(e1==e2){
           e2=rand()%N;
        }
        adjacency[e1][edges_count[e1]]=e2;
        adjacency[e2][edges_count[e2]]=e1;
        edges_count[e1]++;
        edges_count[e2]++;
     }
    
    //continuous step MC
    
     int length=ceil(log(time_max/time_min)/t_b)+1;
     double* p_mean=new double[length]();
     double* times=new double[length]();
     int*    cnt=new int[length]();
     double* time_axis=new double[length]();
     time_axis[0]=time_min;
     for(int i=1;i<length;i++){
     time_axis[i]=time_min*exp(i*t_b);
     }
     
     int* states=new int[N];
     int* active_nodes=new int[N];
     int active_links=0;
     
     for (int i=0; i<N; i++) {
     active_nodes[i]=i;
     states[i]=1;
     }
     
     for (int i=0; i<N; i++) {
     active_links += edges_count[i];
     }
     
     
     int count=0;
     int actives=N;
     double t=0.0;
     
     while (actives>1) {
     
     int pos=rand()%actives;
     int nd=active_nodes[pos];
     
     
     double rnd=randomG();
         
     t += 1.0/(double)(actives+lambda*active_links);
     
     
     
     if (rnd<1.0*actives/(double)(actives+lambda*active_links)) {
     states[nd]=0;
     for (int i=pos; i<actives-1; i++) {
     active_nodes[i]=active_nodes[i+1];
     }
     active_nodes[actives-1]=-1;
     actives--;
     active_links -= edges_count[nd];
     }else{
     
     int* inactive_neighbors=new int[edges_count[nd]];
     int ct=0;
     for (int i=0; i<edges_count[nd]; i++) {
     if (states[adjacency[nd][i]]==0) {
     inactive_neighbors[ct]=adjacency[nd][i];
     ct++;
     }
     }
     
     
     
     int pos1;
     int nd1;
     
     if (ct>=1) {
     pos1=rand()%ct;
     nd1=inactive_neighbors[pos1];
     active_nodes[actives]=nd1;
     states[nd1]=1;
     actives++;
     active_links +=edges_count[nd1];
     }
     
     delete[] inactive_neighbors;
     }
     
     if(t>time_max){
     break;
     }
     
     if(t>time_axis[count]){
     count++;
     if(count<length){
     times[count] += t;
     p_mean[count]+= 1.0*actives/N;
     cnt[count]++;
     }else{
     break;
     }
     }else{
     times[count] += t;
     p_mean[count]+= 1.0*actives/N;
     cnt[count]++;
     }
     
     end=clock();
     used_time=(end-start)/CLOCKS_PER_SEC;
     if(used_time>3600.0){
     break;}
     }
     
    string filename1;
    string filename2;
    filename1=(string)"sis_ER_N"+argv[1]+"E"+argv[2]+"l"+argv[3]+"run"+argv[4]+".txt";
    filename2=(string)"time_cost_ER_N"+argv[1]+"E"+argv[2]+"l"+argv[3]+"run"+argv[4]+".txt";
    ofstream ofs1(filename1.c_str());
    ofstream ofs2(filename2.c_str());
    ofs1<<fixed;
    ofs2<<fixed;
    
    for (int i=1; i<length; i++){
        if (cnt[i]>0) {
            ofs1<<setprecision(10)<<times[i]/cnt[i]<<"\t"<<p_mean[i]/cnt[i]<<"\n";
        }}
    
    end=clock();
    ofs2<<"elapsed time: "<<setprecision(10)<<(end-start)/CLOCKS_PER_SEC<<"seconds"<<"\n";

    
    
    // discrete step MC
    /*
    int     length=(int)time_max;
    double* p_mean=new double[length]();
    double* times=new double[length]();
    int*    cnt=new int[length]();
    double* time_axis=new double[length]();
    
    int* states=new int[N];
    int* active_nodes=new int[N];
    
    for (int i=0; i<N; i++) {
        active_nodes[i]=i;
        states[i]=1;
    }
    
    int actives=N;
    double t=0.0;
    
    // initially
    
    times[0]=0.0;
    p_mean[0]=1.0;
    
    int t0=0;
    while (actives>1) {
        
        int actives_old=actives;
        
        for (int count=0; count<actives_old; count++) {
            
            int pos=rand()%actives;
            int nd=active_nodes[pos];
            if (states[nd]==1) {
                double rnd=randomG();
                if (rnd<1.0/(1.0+lambda)) {
                    states[nd]=0;
                    for (int i=pos; i<actives-1; i++) {
                        active_nodes[i]=active_nodes[i+1];
                    }
                    active_nodes[actives-1]=-1;
                    actives--;
                }else{
                    
                    int* inactive_neighbors=new int[edges_count[nd]];
                    int  ct=0;
                    for (int i=0; i<edges_count[nd]; i++) {
                        if (states[adjacency[nd][i]]==0) {
                            inactive_neighbors[ct]=adjacency[nd][i];
                            ct++;
                        }
                    }
                    
                    int nd1;
                    
                    if (ct>=1) {
                        for (int i=0; i<ct; i++) {
                            double rnd1=randomG();
                            if (rnd1<lambda/(1.0+lambda)/ct) {
                                nd1=inactive_neighbors[i];
                                active_nodes[actives]=nd1;
                                states[nd1]=1;
                                actives++;
                            }
                        }
                        
                    }
                    
                    delete[] inactive_neighbors;
                }
            }
        }
        
        t += 1.0;
        t0++;
        
        
        if(t>time_max){
            break;
        }
        
        times[t0] += t;
        p_mean[t0]+= 1.0*actives/N;
        
        
        
        end=clock();
        used_time=(end-start)/CLOCKS_PER_SEC;
        if(used_time>1000.0){
            break;}
    }
     

    string filename1;
    string filename2;
    filename1=(string)"sis_ER_N"+argv[1]+"E"+argv[2]+"l"+argv[3]+"run"+argv[4]+".txt";
    filename2=(string)"time_cost_ER_N"+argv[1]+"E"+argv[2]+"l"+argv[3]+"run"+argv[4]+".txt";
    ofstream ofs1(filename1.c_str());
    ofstream ofs2(filename2.c_str());
    ofs1<<fixed;
    ofs2<<fixed;
    
    for (int i=1; i<length; i++){
        if (p_mean[i]>0) {
            ofs1<<setprecision(10)<<times[i]<<"\t"<<p_mean[i]<<"\n";
        }
        
    }
    
    end=clock();
    ofs2<<"elapsed time: "<<setprecision(10)<<(end-start)/CLOCKS_PER_SEC<<"seconds"<<"\n";
    */
    
    for(int i = 0; i<N; ++i) {
        delete [] adjacency[i];
    }
    delete[] adjacency;
    delete[] p_mean;
    delete[] times;
    delete[] cnt;
    delete[] time_axis;
    delete[] states;
    delete[] active_nodes;
 
    return 0;
}

