function [Lb,marked_state,initial_state] = Reduced_L(b,g)

% this script is to reduce the Laplacian matrix of MK lattices to hierarchical one dimensional line
% input:
% g --> generation for MK to be reduced
% b --> branching factor for MK to be reduced
%
% output:
% Lb --> the graph Laplacian on the comb (--> "L bar")
% marked_state --> vector |w> in the reduced system
% initial_state --> vector |s> in the reduced system

l=2;
NewNodes=9;
NewBranches=l*NewNodes;
N=2+NewNodes*sum(NewBranches.^[0:1:g-1]);
reduce_N=2^g+1;
Lb=zeros(reduce_N,reduce_N);
marked_state = zeros(reduce_N,1);
marked_state(1,1) = 1;

initial_state =ones(reduce_N,1);

for i=1:g-1
    m=2^(g-1-i);
    for j=1:m
        node=2^i*(2*j-1)+1;
        Lb(node,node)=2*b^i;
        initial_state(node)=sqrt(b^(g-i));
        Lb(node,node-1)=-sqrt(b^i);
        Lb(node-1,node)=-sqrt(b^i);
        Lb(node,node+1)=-sqrt(b^i);
        Lb(node+1,node)=-sqrt(b^i);
    end
end
for j=1:2^(g-1)
    node=2*j-1+1;
    initial_state(node)=sqrt(b^(g));
    Lb(node,node)=2;
end
node=0+1;
Lb(node,node)=b^g;
Lb(node,node+1)=-sqrt(b^g);
Lb(node+1,node)=-sqrt(b^g);
node=2^g+1;
Lb(node,node)=b^g;
Lb(node,node-1)=-sqrt(b^g);
Lb(node-1,node)=-sqrt(b^g);
initial_state = 1/sqrt(N).* initial_state;

Lb = sparse(Lb);
marked_state = sparse(marked_state);
initial_state = sparse(initial_state);

end




    
    
