function ddx = nbodyVect(x,mus)
%NBODYVECT - Fully vectorized n-body equation
% 
% INPUTS
% x     3n x 1 stacked initial position vectors:
%       [r1(1);r1(2);r1(3);r2(1);r2(2)r2(3);...;rn(1);rn(2);rn(3)]
% mus   gravitational parameters (G*m_i) where G is the
%       gravitational constant and m_i is the mass of the ith body.
% OUTPUTS
% ddx   second derivative of x

mus = mus(:); %make mus a column
n = length(mus); %this should be the number of planets

%input check
if (size(x,1) ~= 3*n)
    error('nbodyVect:inputError',['x0 must be 3n x 1 vectors',...
        'and mus must have length n for n bodies']);
end

% precompute index arrays
% k index must be circularly permuted to fill 1 x 3n
% j index must be repeated to fill 1 x 3n
% inds is 1 - 3*n in 3 x n
% inds2 is 1 - 3*n in n x n-1
inds = reshape(1:3*n,3,n);
cols = repmat(1:n,1,n);
ks = cols(logical(ones(n) - eye(n)));
js = reshape(repmat(1:n,n-1,1),1,n*(n-1));
inds2 = reshape(1:n*(n-1),n-1,n).';

% r_k/j = r_k/o - r_j/o
rkj = x(inds(:,ks)) - x(inds(:,js));

%d^2/dt^2(r_k/j) = G*sum_{k ~= j}(mu_k*r_k/j)/(r_k/j^3)
rkj = rkj * diag(1./sqrt(sum(rkj.^2)).^3)*diag(mus(ks));
ddx = sum(reshape(rkj(:,inds2(:)),3*n,n-1),2);
  