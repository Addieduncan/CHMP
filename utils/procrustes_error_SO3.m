function [mean_err] = procrustes_error_SO3(R_est, R_orig, parameters)

    n = parameters.n;

    est = [];
    gt = [];
    for i = 1:n
        est = [est; R_est(:,:,i)]; 
        gt = [gt; R_orig(:,:,i)];    
    end

    % Procrustes alignment to find optimal orthogonal transformation
    [U, ~, V] = svd(gt' * est);
    Q = U * V';

    
    % Compute Frobenius norm error
    err = norm(est - gt * Q, 'fro');
    mean_err = err / n;

end