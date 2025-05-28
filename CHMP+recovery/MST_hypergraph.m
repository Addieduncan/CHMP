function [R_est,recovery_time] = MST_hypergraph(hypergraph,parameters,CHMP_out)

    % get data
    adjTensor = hypergraph.AdjTen;
    n = parameters.n;
    SVec = CHMP_out.SVec;
    Ind = hypergraph.Ind;
    RijkMat = hypergraph.RijkMat;
   
    
    tic;
    adjMat = zeros(n,n);
    
    edge_adjMat = zeros(n,n);
    
    for t = 1:size(Ind,1)
        I = Ind(t,1);
        J = Ind(t,2);
        K = Ind(t,3);
        if adjMat(I,J) == 0
            adjMat(I,J) = SVec(t);
            edge_adjMat(I,J) = t;
        elseif SVec(t) < adjMat(I,J)
            adjMat(I,J) = SVec(t);
            edge_adjMat(I,J) = t;
        end
    
        if adjMat(J,K) == 0
            adjMat(J,K) = SVec(t);
            edge_adjMat(J,K) = t;
        elseif SVec(t) < adjMat(J,K)
            adjMat(J,K) = SVec(t);
            edge_adjMat(J,K) = t;
        end
    
        if adjMat(I,K) == 0
            adjMat(I,K) = SVec(t);
            edge_adjMat(I,K) = t;
        elseif SVec(t) < adjMat(I,K)
            adjMat(I,K) = SVec(t);
            edge_adjMat(I,K) = t;
        end
    end
    
    
    [Ind_i,Ind_j,W_ij] = find(triu(adjMat));
    
    edge_adjMat = edge_adjMat + edge_adjMat';
    
    
    G = graph(Ind_i,Ind_j,W_ij);
    MST = minspantree(G);
    
    AdjTree = adjacency(MST);
    AdjTree = full(AdjTree);
    
   
    % compute Ri by multiplying Rij along the spanning tree
    rootnodes = 1;
    added=zeros(1,n);
    R_est = zeros(3,3,n);
    R_est(:,:,rootnodes)=eye(3);
    added(rootnodes)=1;
    newroots = [];
    
    %% MST
    
    while sum(added)<n
        for node_root = rootnodes
            leaves = find((AdjTree(node_root,:).*(1-added))==1);
            newroots = [newroots, leaves];
            for node_leaf=leaves
                edge_leaf = edge_adjMat(node_leaf,node_root);
                I = Ind(edge_leaf,1);J = Ind(edge_leaf,2);K = Ind(edge_leaf,3);
                if (node_root == I) & (node_leaf == J)
                    R_est(:,:,node_leaf)=(RijkMat(:,:,edge_leaf,1))'*R_est(:,:,node_root);
                elseif (node_root == I) & (node_leaf == K)
                    R_est(:,:,node_leaf)=(RijkMat(:,:,edge_leaf,1)*RijkMat(:,:,edge_leaf,2))'*R_est(:,:,node_root);
                elseif (node_root == J) & (node_leaf == K)
                    R_est(:,:,node_leaf)=(RijkMat(:,:,edge_leaf,2))'*R_est(:,:,node_root);
                elseif (node_root == J) & (node_leaf == I)
                    R_est(:,:,node_leaf)=(RijkMat(:,:,edge_leaf,1))*R_est(:,:,node_root);
                elseif (node_root == K) & (node_leaf == I)
                    R_est(:,:,node_leaf)=(RijkMat(:,:,edge_leaf,1)*RijkMat(:,:,edge_leaf,2))*R_est(:,:,node_root);
                elseif (node_root == K) & (node_leaf == J)
                    R_est(:,:,node_leaf)=(RijkMat(:,:,edge_leaf,2))*R_est(:,:,node_root);
                end
                added(node_leaf)=1;
            end
        end
        rootnodes = newroots;
    end
    recovery_time = toc;
end