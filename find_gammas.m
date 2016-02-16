%clear;
% find quantum walk parameters gamma and corresponding success probability on Mk lattices
function [gammas,max_p_values]=find_gammas(gs)



g_list=gs; 	            % list of generations for MKs
b =9;                   % branching factor
n_list= [];
gammas=cell(length(g_list),1);

for cnt=1:length(g_list)
    n_list=[n_list,2+2*sum((2*b).^[0:1:g_list(cnt)-1])];  % N list
    gammas{cnt}=[0.000000001:0.000000001:0.00000001];
end



% max_p_values saves the data
max_p_values = cell(length(g_list),1);


for cnt =1:length(g_list)	
	
	tic
    
	
	% use equally spaced time points or random ones
	max_time = n_list(cnt);	
    num_times = n_list(cnt)/10000000;	% number of time samples
	times = [ 1:max_time/num_times:max_time ]; 
	
	

    % generate reduced Hamiltonian and vectors
	[ Lb, marked_state, initial_state ] = Reduced_L(b,g_list(cnt));
	

	for count = 1:length(gammas{cnt}())

		Hb=gammas{cnt}(count).*Lb-diag(marked_state);
		p_values = zeros(1,length(times));		
		t_count = 1;

		for t = times
			
			p_values(1,t_count) = abs( dot( expm(-i*Hb*t) * initial_state, marked_state ) )^2;
            % careful when using this 'definition' of the optimal gamma (which is actually the
            % correct one) since one has to exclude t = 0 on the one hand but not start to late on
            % the other hand in order to not miss the minimum); also swap comments for the command
            % for max_p_values below
            %p_values(1,t_count) = t/abs( dot( expm(-i*Hb*t) * initial_state, marked_state ) )^2;
			
            t_count = t_count + 1;

		end

		max_p_values{cnt}=[max_p_values{cnt},max(p_values)];
        

	end

	
	toc
	

end


%loglog(gammas{1}(:),max_p_values{1}(:),'r',gammas{2}(:),max_p_values{2}(:),'m',gammas{3}(:),max_p_values{3}(:),'b')
%semilogy(gammas{1}(:),max_p_values{1}(:),'r')

end
