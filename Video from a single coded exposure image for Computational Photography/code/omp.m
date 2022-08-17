function [theta] = omp(A,y,noise_variance)
%omp Summary of this function goes here
%   Detailed explanation goes here

[~,n] = size(A);
[m,~] = size(y);

%initialize index vector to keep track of max index
index_vector = [];

%initiaize matrix to keep track of basis in every iteration
basis_matrix = [];

%initialize the residual vector to y
residue = y;

%initialize theta vector to 0
theta = zeros(n,1);

iteration_val = 0;
while(norm(residue)^2 > 3*m*noise_variance && iteration_val <= m)
   
    %Find the index of maximum projection of y along columns of A
    [~, current_max_index] = max(abs(residue'*normc(A)));

    %Keep track of this index and add the basis vector array with the column
    %having this index in A
    index_vector = [index_vector current_max_index];
    
    %Find the columns in A with respect to the index vector
    basis_matrix = A(:,index_vector);
    
    %Find the estimate of theta_temp 
    theta_temp = pinv(basis_matrix)*y;
    
    %Find the new residue vector
    residue = y - basis_matrix*theta_temp;
    
    %Increase the iteration counter
    iteration_val = iteration_val + 1;
end

%Create the sparse prediction vector
theta(index_vector,:) = theta_temp;
end



