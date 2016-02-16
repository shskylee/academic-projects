function DF_eigen(l,g)
tic
b=2;
NewNodes=l;
NewBranches=b*NewNodes;
Generation=g;
Nodes=zeros(1,Generation);
for k=1:Generation
    Nodes(k)=2+NewNodes*sum(NewBranches.^[0:1:k-1]);
end

N=Nodes(Generation);

%edges=[[1,3];[1,4];[2,3];[2,4]]; %l=2
%edges=[[1,3];[1,4];[1,5];[2,3];[2,4];[2,5]];   %l=3
%edges=[[1,3];[1,4];[1,5];[1,6];[2,3];[2,4];[2,5];[2,6]];   %l=4
edges=[[1,3];[1,4];[1,5];[1,6];[1,7];[2,3];[2,4];[2,5];[2,6];[2,7]];
%%l=5


for i=2:Generation
    edges_new=[];
    for k=1:NewBranches^(i-1)
        for j=1:NewNodes
            edges_new=[edges_new;[edges(k,1),Nodes(i-1)+(k-1)*NewNodes+j]];
            edges_new=[edges_new;[edges(k,2),Nodes(i-1)+(k-1)*NewNodes+j]];
        end
    end
    edges=edges_new;
end


% graph the network && reduce the network

%{
r=size(edges,1);
row=[];
col=[];
v=[];
for i=1:r
    row=[new_row,edges(i,1),edges(i,2)];
    col=[col,edges(i,2),edges(i,1)];
    v=[v,1,1];
end
% adjacency matrix is
adj=sparse(new_row,col,v,N,N);
S=1; 
% S: starting node for graph travese
order=graphtraverse(adj, S,'Method','BFS');
neworder=1:N;
orderMap = containers.Map(order,neworder);

new_row=[];
new_col=[];
for i=1:r
    m=edges(i,1);
    n=edges(i,2);
    new_row=[new_row,orderMap(m)];
    new_col=[new_col,orderMap(n)];
 end
G=graph(new_row,new_col);
plot(G)

adj=adjacency(G);
adj=full(adj);
Dmat=zeros(N,N);
for i=1:N
    Dmat(i,i)=sum(adj(i,:));
end
lap=Dmat-adj;

%construct reduction vector reduce_v
repeat=[1];
cnt=1;
for i=2:N
    if lap(i,i)~=lap(i-1,i-1)
        repeat=[repeat,1];
        cnt=cnt+1;
    else
        repeat(cnt)=repeat(cnt)+1;
    end
end

rs=size(repeat,2);
reduce_v=zeros(rs,N);
csum=cumsum(repeat);
reduce_v(1,1)=1;

for i=2:rs
    for j=csum(i-1)+1:csum(i)
        reduce_v(i,j)=1/sqrt(repeat(i));
    end
end

% reduced laplacian is
oracle=zeros(N,N);
oracle(1,1)=1;


reduce_lap=reduce_v*lap*pinv(reduce_v);
% prove the reduced laplacian is idential to a one dimensional hierarchical
% line
reduce_N=2^g+1;
hierarchical_1D=zeros(reduce_N,reduce_N);
for i=1:g-1
    m=2^(g-1-i);
    for j=1:m
        node=2^i*(2*j-1)+1;
        hierarchical_1D(node,node)=2*l^i;
        hierarchical_1D(node,node-1)=-sqrt(l^i);
        hierarchical_1D(node-1,node)=-sqrt(l^i);
        hierarchical_1D(node,node+1)=-sqrt(l^i);
        hierarchical_1D(node+1,node)=-sqrt(l^i);
    end
end
for j=1:2^(g-1)
    node=2*j-1+1;
    hierarchical_1D(node,node)=2;
end
node=0+1;
hierarchical_1D(node,node)=l^g;
hierarchical_1D(node,node+1)=-sqrt(l^g);
hierarchical_1D(node+1,node)=-sqrt(l^g);
node=2^g+1;
hierarchical_1D(node,node)=l^g;
hierarchical_1D(node,node-1)=-sqrt(l^g);
hierarchical_1D(node-1,node)=-sqrt(l^g);

isequal(hierarchical_1D,reduce_lap)

(hierarchical_1D-reduce_lap)*(1e+12)
%}    
    
% construct laplacian matrix & find eigenvectors

r=size(edges,1);
v_d=zeros(1,N);
row=[];
col=[];
v=[];
for i=1:r
    row=[row,edges(i,1),edges(i,2)];
    col=[col,edges(i,2),edges(i,1)];
    v_d(edges(i,1)) = v_d(edges(i,1))+1;
    v_d(edges(i,2)) = v_d(edges(i,2))+1;
    v=[v,-1,-1];
end
row=[row,[1:1:N]];
col=[col,[1:1:N]];
v=[v,v_d];

lap=sparse(row,col,v,N,N);
%[evec,eval]=eigs(lap,4,'sm');

%[evec,eval]=eig(full(lap));

eval=eig(full(lap));

filename = sprintf('%s%d%s%d%s','MK_eval_b',l,'g',g,'.txt');
fileID = fopen(filename,'w');
fprintf(fileID,'%s\t%s\t%s\n','N','i','eigenvalues[i]');
for i=1:N
    fprintf(fileID,'%d\t%d\t%12.8f\n',N,i,abs(eval(i)));
end

fclose(fileID);



%laplacian=full(lap);

%[evec,eval]=eig(full(lap));

%{
average=mean(sqrt(evec.^2),2);
a=sqrt(evec.^2);
average=a(:,3*N/4);
filename = sprintf('%s%d%s%d%s','DF_evec_mean_l',l,'g',g,'.txt')
fileID = fopen(filename,'w');
fprintf(fileID,'%14.12f\n',average);
fclose(fileID);
%}

toc

end




    
    
