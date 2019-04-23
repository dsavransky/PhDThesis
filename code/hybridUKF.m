function [xhist,Phist] = hybridUKF(t,z,x0,P0,R,f,h)
%hybridEKF - Hybrid Unscented Kalman Filter
%
% INPUTS
% t     (N x 1) Array of times for N observations
% z     (N x m) matrix of observations for m element h function
% x0    (n x 1) vector of nitial conditions for n element state
% P0    (n x n) initial covariance matrix
% R     (m x m) observation noise covariance
% f     Function handle to (possible nonlinear) dynamics function. 
%       Must accept 2 arguments: scalar time and (n x 2n+1) 
%       matrix of sigma points so that f(t,x) = dx (n x 2n+1)
% h     Function handle to (possibly nonlinear) observation function. 
%       Input argument is an (n x 2n+1) matrix of sigma points 
%       so that h(x) = y (m x 2n+1)
%
% OUTPUTS
% xhist     (n x N) time history of filter state estimates
% Phist     (n x n x N) time history of filter covariance estimates


%basic dimensions
n = length(x0);
d = 2*n+1;

%prepare arrays
xminus = x0;
Pminus = P0;
xhist = zeros(n,length(t));
Phist = zeros(n,n,length(t));

%form weight matrices
%alph = 1e-3;
alph = 0.5;
bet = 2;
%kappa = 1;
kappa = 3-n;
lambda = alph^2*(n+kappa) - n;
c = alph^2*(n+kappa);
w_m = ones(d,1)*1/(2*(n+lambda));
w_c = w_m;
w_m(1) = lambda/(n+lambda);
w_c(1) = lambda/(n+lambda) + (1 - alph^2 + bet);
temp = eye(d) - repmat(w_m,1,d);
W = temp*diag(w_c)*temp.';

%do one update to prime the filter
A = chol(Pminus,'lower');
Xminus = repmat(xminus,1,d) + sqrt(c)*[zeros(n,1),A,-A];
Yminus = h(Xminus);
mu = Yminus*w_m;
S = Yminus*W*Yminus.' + R;
C = Xminus*W*Yminus.';
K = C/S;
x = xminus + K*(z(1,:).' - mu);
P = Pminus - K*S*K.';
xhist(:,1) = x;
Phist(:,:,1) = P;

%iterate
for k = 2:length(t)
    disp(k);
    %propagate
    [~,Zo] = ode45(@hybridUKF_deq,[t(k-1),t(k)],[x;P(:)]);
    xminus = Zo(end,1:n).';
    Pminus = reshape(Zo(end,n+1:end),n,n);
    
    %update
    A = chol(Pminus,'lower');
    Xminus = repmat(xminus,1,d) + sqrt(c)*[zeros(n,1),A,-A];
    Yminus = h(Xminus);
    mu = Yminus*w_m;
    S = Yminus*W*Yminus.' + R;
    C = Xminus*W*Yminus.';
    K = C/S;
    x = xminus + K*(z(k,:).' - mu);
    P = Pminus - K*S*K.';
end

    function dZ = hybridUKF_deq(T,Z)
        %unpackage input
        xi = Z(1:n);
        Pi = reshape(Z(n+1:end),n,n);
        
        %form sigma-point matrix
        A = chol(Pi);
        Xi = repmat(xi,1,d) + sqrt(c)*[zeros(n,1),A,-A];
        
        %integrate
        fXi = f(T,Xi);
        dxi = fXi*w_m;
        dPi = Xi*W*fXi.' + fXi*W*Xi.';
        
        %repackage
        dZ = [dxi;dPi(:)];
    end

end