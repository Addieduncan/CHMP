function[CHMP_out] = CHMP(hypergraph,parameters)
    %% get parameters

    % CHMP parameters
    beta_init = parameters.beta_init;
    beta_max = parameters.beta_max;
    rate = parameters.rate;
    sampling = parameters.sampling; % percent of cycles to sample


    % hypergraph parameters
    Ind = hypergraph.Ind;
    RijkMat = hypergraph.RijkMat;
    ErrVec = hypergraph.ErrVec;

    n = parameters.n; 
    num_hyperedges = hypergraph.num_hyperedges; 

 
    %% start CHMP initialization 
    disp('Initializing CHMP')
    tic;

    %% building the graph   
    Ind_i = Ind(:,1);
    Ind_j = Ind(:,2);
    Ind_k = Ind(:,3);

    %% find hypergraph cycles

    RijkMat5d = zeros(3,3,n,n,n,2);
    IndTen = zeros(n,n,n); 

    for l = 1:num_hyperedges
        i=Ind_i(l);j=Ind_j(l);k=Ind_k(l);
        RijkMat5d(:,:,i,j,k,:)=RijkMat(:,:,l,:); 
        RijkMat5d(:,:,i,k,j,1)=RijkMat(:,:,l,1)*RijkMat(:,:,l,2);
        RijkMat5d(:,:,i,k,j,2)=RijkMat(:,:,l,2)';
        RijkMat5d(:,:,j,i,k,1)=RijkMat(:,:,l,1)';
        RijkMat5d(:,:,j,i,k,2)=RijkMat(:,:,l,1)*RijkMat(:,:,l,2);
        RijkMat5d(:,:,j,k,i,1)=RijkMat(:,:,l,2);
        RijkMat5d(:,:,j,k,i,2)=(RijkMat(:,:,l,1)*RijkMat(:,:,l,2))';
        RijkMat5d(:,:,k,i,j,1)=(RijkMat(:,:,l,1)*RijkMat(:,:,l,2))';
        RijkMat5d(:,:,k,i,j,2)=RijkMat(:,:,l,1);
        RijkMat5d(:,:,k,j,i,1)=RijkMat(:,:,l,2)';
        RijkMat5d(:,:,k,j,i,2)=RijkMat(:,:,l,1)';
        IndTen(i,j,k)=l; 
        IndTen(i,k,j)=l;
        IndTen(j,i,k)=l;
        IndTen(j,k,i)=l;
        IndTen(k,i,j)=l;
        IndTen(k,j,i)=l;
    end

    % create cycle tracking vectors 
    edge_cycle_count = zeros(num_hyperedges,1); 
    cycles = []; 
    Ind_ijk = [];
    Ind_jkl = [];
    Ind_kli = [];
    Ind_lij = [];
    
   
    %% sample cycles

    % list all possible cycles
    possible_cycles = nchoosek(1:n,4);
    max_cycle = nchoosek(n,4);

    k_par = round(sampling * max_cycle); % number of cycles to sample, if sampling=1 then all cycles will be considered

    % sample cycles
    random_cycles = randperm(max_cycle,k_par); 
    sampled_cycles = sortrows(possible_cycles(random_cycles, :));

    % find all cycles and track indices
    for s = 1:k_par
        edges_found = all(ismember(Ind, sampled_cycles(s,:)), 2);
        edge_indices = find(edges_found);
        if size(edge_indices,1) == 4
            edge_cycle_count(edge_indices(1)) = edge_cycle_count(edge_indices(1)) + 1;
            edge_cycle_count(edge_indices(2)) = edge_cycle_count(edge_indices(2)) + 1;
            edge_cycle_count(edge_indices(3)) = edge_cycle_count(edge_indices(3)) + 1;
            edge_cycle_count(edge_indices(4)) = edge_cycle_count(edge_indices(4)) + 1;
            cycles = [cycles; edge_indices'];
            i = sampled_cycles(s,1);
            j = sampled_cycles(s,2);
            k = sampled_cycles(s,3);
            l = sampled_cycles(s,4);
            Ind_ijk = [Ind_ijk, [IndTen(i,j,k);l],[IndTen(j,k,l);i],[IndTen(k,l,i);j],[IndTen(l,i,j);k]];
            Ind_jkl = [Ind_jkl, IndTen(j,k,l),IndTen(k,l,i),IndTen(l,i,j),IndTen(i,j,k)];
            Ind_kli = [Ind_kli, IndTen(k,l,i),IndTen(l,i,j),IndTen(i,j,k),IndTen(j,k,l)];
            Ind_lij = [Ind_lij, IndTen(l,i,j),IndTen(i,j,k),IndTen(j,k,l),IndTen(k,l,i)];
        else
            continue;
        end
    end


    %% check for cycles
    if isempty(cycles)
        disp('No cycles in hypergraph. Terminating CHMP.')
        success = 0;

        SVec = ones(1,num_hyperedges);
        maxErr = ones(1,num_hyperedges);
        meanErr = ones(1,num_hyperedges);
        medianErr = ones(1,num_hyperedges);
        

        initialization_time = toc;

        iteration_time = 0;

    else   
        success = 1;

        %% initialize cycle inconsitency and weights

        % clean cycle data
        edge_cycle = unique(cycles); 
        edge_cycle_count = edge_cycle_count(edge_cycle); 
        
        cum_ind = [0;cumsum(edge_cycle_count)]; 
        m_pos = length(edge_cycle); 
        m_cycle = cum_ind(end); 
        
        
        % group according to edge ijk
        [~,ijksort]=sort(Ind_ijk(1,:));
        Ind_ijk = Ind_ijk(:,ijksort);
        Ind_jkl = Ind_jkl(ijksort);
        Ind_kli = Ind_kli(ijksort);
        Ind_lij = Ind_lij(ijksort);
        
        
        % initialize cycle potentials
        Rijk0Mat = zeros(3,3,m_cycle,2);
        Rjkl0Mat = zeros(3,3,m_cycle,2);
        Rkli0Mat = zeros(3,3,m_cycle,2);
        Rlij0Mat = zeros(3,3,m_cycle,2);

        % generate cycle inconsistency measures
        for s = 1:m_cycle
            t = Ind_ijk(1,s);
            i = Ind(t,1);
            j = Ind(t,2);
            k = Ind(t,3);
            l = Ind_ijk(2,s);
            Rijk0Mat(:,:,s,:) = RijkMat5d(:,:,i,j,k,:);
            Rjkl0Mat(:,:,s,:) = RijkMat5d(:,:,j,k,l,:);
            Rkli0Mat(:,:,s,:) = RijkMat5d(:,:,k,l,i,:);
            Rlij0Mat(:,:,s,:) = RijkMat5d(:,:,l,i,j,:);
        end
        
        R_cycle1 = zeros(3,3,m_cycle);
        R_cycle2 = zeros(3,3,m_cycle);
        R_cycle = zeros(3,3,m_cycle); % stores cycle product
        
        for j = 1:3
          R_cycle1 = R_cycle1 + bsxfun(@times,Rijk0Mat(:,j,:,:),Rjkl0Mat(j,:,:,:));
        end
        
        for j = 1:3
          R_cycle2 = R_cycle2 + bsxfun(@times,R_cycle1(:,j,:,:),Rkli0Mat(j,:,:,:));
        end
        
        for j = 1:3
          R_cycle = R_cycle + bsxfun(@times,R_cycle2(:,j,:,:),Rlij0Mat(j,:,:,:));
        end
        
        R_trace_one = reshape(R_cycle(1,1,:,1)+R_cycle(2,2,:,1)+R_cycle(3,3,:,1),[1,m_cycle]);
        R_trace_two = reshape(R_cycle(1,1,:,2)+R_cycle(2,2,:,2)+R_cycle(3,3,:,2),[1,m_cycle]);
        S0_long = sqrt(((1/sqrt(2))*abs((acos((R_trace_one-1)./2)/pi))).^2 + ((1/sqrt(2))*abs((acos((R_trace_two-1)./2)/pi))).^2);
        S0_vec = ones(1,num_hyperedges); % stores cycle inconsistency measure
        
        % generate initial weights
        Weight_vec = ones(1,m_cycle);
        S0_weight = S0_long.*Weight_vec;
        
        for l=1:m_pos
            IJK = edge_cycle(l);
            S0_vec(IJK) = sum(S0_weight((cum_ind(l)+1):cum_ind(l+1)))/sum(Weight_vec((cum_ind(l)+1):cum_ind(l+1)));
        end
        
        disp('Initialization completed');

        initialization_time = toc;

        %% iterate CHMP
        
        disp('Starting reweighting procedure');
        
        tic;
        
        iter = 0;
        
        SVec = S0_vec;
        beta = beta_init;
    
        maxErr = [];
        meanErr = [];
        medianErr = [];
    
        while beta <= beta_max
            Sjkl = SVec(Ind_jkl);
            Skli = SVec(Ind_kli);
            Slij = SVec(Ind_lij);
            S_sum = Sjkl+Skli+Slij;
            
            Weight_vec = exp(-beta*S_sum);
            S0_weight = S0_long.*Weight_vec;
        
            for l=1:m_pos
                IJK = edge_cycle(l);
                SVec(IJK) = sum(S0_weight((cum_ind(l)+1):cum_ind(l+1)))/sum(Weight_vec((cum_ind(l)+1):cum_ind(l+1)));
            end
   
    
            % track errors
            t_maxErr = log(max(abs(ErrVec-SVec)));
            maxErr = [maxErr;iter,t_maxErr];

            t_meanErr = log(sum(abs(ErrVec-SVec))/size(ErrVec,2));
            meanErr = [meanErr;iter,t_meanErr];

            t_medianErr = log(median(abs(ErrVec-SVec)));
            medianErr = [medianErr;iter,t_medianErr];
            
 
            iter = iter+1;

            % parameter controling the decay rate of reweighting function
            beta = beta*rate;
        end
        iteration_time = toc;
        disp('CHMP completed')
    end


    
    CHMP_out.SVec = SVec;
    CHMP_out.maxErr = maxErr;
    CHMP_out.meanErr = meanErr;
    CHMP_out.medianErr = medianErr;
    CHMP_out.success = success;
    CHMP_out.initialization_time = initialization_time;
    CHMP_out.iteration_time = iteration_time;
    


end