function[rotations]=generateRotations(parameters)

% get parameters
n = parameters.n;


% generate 3D rotations
R_orig = zeros(3,3,n);

for i = 1:n
    Q=randn(3);
    [U, ~, V]= svd(Q);
    S0 = diag([1,1,det(U*V')]);  
    R_orig(:,:,i)=U*S0*V';
end

% save output
rotations.R_orig = R_orig;