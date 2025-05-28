function[R_est,recovery_time] = GCW_hypergraph(hypergraph,parameters,CHMP_out)
    
    n = parameters.n;
    RijkMat = hypergraph.RijkMat;
    Ind = hypergraph.Ind;
    SVec = CHMP_out.SVec;
    beta_T = parameters.beta_max;

    tic;
    MijMat = zeros(3*n,3*n);
    AdjMat = zeros(n,n);
    SVecMat = zeros(n,n);
    mat_size = ones(1,n)*3;
    cum_ind = [0,cumsum(mat_size)];

    for i =1:n
        for j = i+1:n
            containing_edges = find(any(Ind == i,2) & any(Ind == j,2));
            if ~isempty(containing_edges)
                [min_weight,idx] = min(SVec(containing_edges));
                SVecMat(i,j) = min_weight;
                edge = containing_edges(idx);
                I = Ind(edge,1);J = Ind(edge,2);K = Ind(edge,3);
                if (i == I) & (j == J)
                    MijMat((cum_ind(i)+1):cum_ind(i+1), (cum_ind(j)+1):cum_ind(j+1)) = reshape(RijkMat(:,:,edge,1),[3,3]);
                elseif (i == I) & (j == K)
                    MijMat((cum_ind(i)+1):cum_ind(i+1), (cum_ind(j)+1):cum_ind(j+1)) = (reshape(RijkMat(:,:,edge,1),[3,3]))*(reshape(RijkMat(:,:,edge,2),[3,3]));
                elseif (i == J) & (j == K)
                    MijMat((cum_ind(i)+1):cum_ind(i+1), (cum_ind(j)+1):cum_ind(j+1)) = reshape(RijkMat(:,:,edge,2),[3,3]);
                elseif (i == J) & (j == I)
                    MijMat((cum_ind(i)+1):cum_ind(i+1), (cum_ind(j)+1):cum_ind(j+1)) = (reshape(RijkMat(:,:,edge,1),[3,3]))';
                elseif (i == K) & (j == I)
                    MijMat((cum_ind(i)+1):cum_ind(i+1), (cum_ind(j)+1):cum_ind(j+1)) = ((reshape(RijkMat(:,:,edge,1),[3,3]))*(reshape(RijkMat(:,:,edge,2),[3,3])))';
                elseif (i == K) & (j == J)
                    MijMat((cum_ind(i)+1):cum_ind(i+1), (cum_ind(j)+1):cum_ind(j+1)) = (reshape(RijkMat(:,:,edge,2),[3,3]))';
                end
                AdjMat(i,j) = 1;
            else
                AdjMat(i,j) = 0;
            end
        end
    end

    

    MijMat = MijMat + MijMat';

    AdjMat = AdjMat + AdjMat';

    SVecMat = SVecMat + SVecMat';

    Weights = exp(-beta_T.*SVecMat).*AdjMat;
    Weights = diag(1./sum(Weights,2))*Weights;
    Weights = kron(Weights, ones(3));  


    MMat = MijMat.*Weights;

    for i = 1:n
        MMat((cum_ind(i)+1):cum_ind(i+1), (cum_ind(i)+1):cum_ind(i+1)) = eye(3);
    end


    [V,~] = eigs(MMat,3,'la');
    V(:,1) = V(:,1)*sign(det(V(1:3,:))); 
    R_est = zeros(3,3,n);
    for i=1:n
       Ri = V((cum_ind(i)+1):cum_ind(i+1), :); 
       [Ur,~,Vr] = svd(Ri);
       S0 = diag([ones(1,3-1),det(Ur*Vr')]);
       R_est(:,:,i) = Ur*S0*Vr';
       
    end
    recovery_time = toc;

end