function HN3_eigen(g) 
 tic,
 N=2^g;
 edges=[[1,N/2+1]];

 for i=1:g-1
     if i==1
        for j=1:2^(g-1-i)
            k=2^(i-1)*(4*j-3)+1;
            l=2^(i-1)*(4*j-1)+1;
            edges=[edges;[k,l];[k,k-1];[k,k+1];[l,l-1]];
            if l+1==N+1
               edges=[edges;[l,1]];
            else
               edges=[edges;[l,l+1]];
            end
        end
     else
         for j=1:2^(g-1-i)
            k=2^(i-1)*(4*j-3)+1;
            l=2^(i-1)*(4*j-1)+1;
            edges=[edges;[k,l]];
         end
     end
 end
         
 
 row=[];
 col=[];
 v=[];
 r=size(edges,1);
 for i=1:r
     row=[row,edges(i,1),edges(i,2)];
     col=[col,edges(i,2),edges(i,1)];
     v=[v,-1,-1];
 end
 
 row=[row,1:N];
 col=[col,1:N];
 v=[v,3*ones(1,N)];
 lap=sparse(row,col,v,N,N);
 
 % [evec,eval]=eigs(lap,4,'sm');
 % eigenval=diag(eval);
 
 eval=eig(full(lap));
 
 filename = sprintf('%s%d%s%d%s%d%s','HN3_eval_g',g,'.txt');
 fileID = fopen(filename,'w');
 fprintf(fileID,'%s\t%s\t%s\n','N','i','eigenvalues[i]');
 for i=1:N
     fprintf(fileID,'%d\t%d\t%12.8f\n',N,i,abs(eval(i)));
 end
 
 fclose(fileID);
    
 
 toc
 
end
    
    
