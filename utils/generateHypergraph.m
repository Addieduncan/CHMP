function[hypergraph]=generateHypergraph(parameters,rotations,noise_model)

    %% get parameters
    n = parameters.n;
    p = parameters.p;
    q = parameters.q;
    sigma = parameters.sigma;

    R_orig = rotations.R_orig;


    %% generate hypergraph
    G = rand(n,n,n) < p;
    AdjTen = zeros(n,n,n);
    Ind = []; 
    
    for i = 1:n
        for j = i+1:n
            for k = j+1:n
                if G(i,j,k) == 1
                    Ind = [Ind;i,j,k];
                    AdjTen(i,j,k) = 1;
                    AdjTen(i,k,j) = 1;
                    AdjTen(j,i,k) = 1;
                    AdjTen(j,k,i) = 1;
                    AdjTen(k,i,j) = 1;
                    AdjTen(k,j,i) = 1;
                end
            end
        end
    end
  
    num_hyperedges = size(Ind,1);
    
   %% generate hyperedge potentials from rotations
    Rijk_orig = zeros(3,3,num_hyperedges,2);
    for t = 1:num_hyperedges
        i=Ind(t,1); j=Ind(t,2); k=Ind(t,3);
        Rijk_orig(:,:,t,1)=R_orig(:,:,i)*(R_orig(:,:,j)');
        Rijk_orig(:,:,t,2)=R_orig(:,:,j)*(R_orig(:,:,k)');
    end
    
    %% add noise
    RijkMat = Rijk_orig;

    noiseIndLog = rand(1,num_hyperedges)>=q;
    noiseInd=find(noiseIndLog); 
    
    switch noise_model
        case 'uniform'
            RijkMat(:,:,noiseInd,:)= RijkMat(:,:,noiseInd,:)+sigma*(rand(3,3,length(noiseInd),2)-0.5); % Uniform noise in range [-sigma/2, sigma/2]
        case 'gaussian'
            RijkMat(:,:,noiseInd,:)= RijkMat(:,:,noiseInd,:)+sigma*randn(3,3,length(noiseInd),2); % add Gaussian noise
        otherwise
            error('Invalid noise type. Choose ''uniform'' or ''gaussian''.');
    end

    % project back to SO(3)
    for k = noiseInd
        for i =1:2
            [U, ~, V]= svd(RijkMat(:,:,k,i));
            S0 = diag([1,1,det(U*V')]);
            RijkMat(:,:,k,i) = U*S0*V';
        end
    end 
    
    %% add corruption

    corrIndLog = logical(1-noiseIndLog); 
    corrInd=find(corrIndLog); 

    for k = corrInd
        for i =1:2
            Q=randn(3);
            [U, ~, V]= svd(Q);
            S0 = diag([1,1,det(U*V')]);
            RijkMat(:,:,k,i) = U*S0*V';  
        end
    end
    
    %% hyperedge corruption error (using l2 product metric with geodesic distance on SO(3)
    
    % error for first element
    R_err_one = zeros(3,3,num_hyperedges);
    for j = 1:3
      R_err_one = R_err_one + bsxfun(@times,Rijk_orig(:,j,:,1),RijkMat(:,j,:,1));
    end
    
    R_err_trace_one = (reshape(R_err_one(1,1,:)+R_err_one(2,2,:)+R_err_one(3,3,:), [num_hyperedges,1]))';
    
    % error for second element
    R_err_two = zeros(3,3,num_hyperedges);
    for j = 1:3
      R_err_two = R_err_two + bsxfun(@times,Rijk_orig(:,j,:,2),RijkMat(:,j,:,2));
    end
    
    R_err_trace_two = (reshape(R_err_two(1,1,:)+R_err_two(2,2,:)+R_err_two(3,3,:), [num_hyperedges,1]))';

    % total error
    
    ErrVec = sqrt(((1/sqrt(2))*abs((acos((R_err_trace_one-1)./2)/pi))).^2 + ((1/sqrt(2))*abs((acos((R_err_trace_two-1)./2)/pi))).^2);

    %% save output
    hypergraph.Ind = Ind;
    hypergraph.RijkMat = RijkMat;
    hypergraph.ErrVec = ErrVec;
    hypergraph.num_hyperedges = num_hyperedges;
    hypergraph.AdjTen = AdjTen;
    
    
end