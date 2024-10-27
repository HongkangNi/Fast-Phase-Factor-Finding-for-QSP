function y = half_chol(G,v)
% Input:
% G: an matrix of size n*r, where n>>r
%
% Output:
% y = L\v, where L is the LDL decomposition of A. The matrix A is described as follows.
% 
% Let Z = diag(ones(n-1,1),-1) be the lower shifting matrix.
% There is an unique n*n matrix A satisfying the displacement equation
% A - Z*A*Z'=G*G'
% This function is an O(rn^2) fast solver for calculating the LDL
% decomposition of this A, i.e. A = L*D*L', and then output the vector 
% y = L\v
%
% Remarks:
% The matrix A need not be given as an input, since we will not use it
% explicitly.
% The matrix L is also not constructed explicitly, since only the vector y = L\v is needed.
%
% Time complexity: O(rn^2)
% Space complexity: O(rn)

[n,~] = size(G);
Gi = G';
y = v';
for i = 1:n-1
    [~,Gi] = qr(Gi);
    u1 = Gi(1,2:end);
    y(i+1:end) = y(i+1:end) - u1/Gi(1,1)*y(i);
    Gi = [Gi(1,1:end-1);Gi(2:end,2:end)];    
end
y = y';
    
    
