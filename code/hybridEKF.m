function [xhist,Phist] = hybridEKF(t,z,x0,P0,R,Q,f,F,h,H)
%hybridEKF - Hybrid Extended Kalman Filter
%
% INPUTS
% t     (N x 1) Array of times for N observations
% z     (N x m) matrix of observations for m element h function
% x0    (n x 1) vector of nitial conditions for n element state
% P0    (n x n) initial covariance matrix
% R     (m x m) observation noise covariance
% Q     (n x n) process noise covariance
% f     Function handle to (possible nonlinear) dynamics function. 
%       Must accept 2 arguments: scalar time and (n x 1) state vector 
%       so that f(t,x) = dx (n x 1)
% F     Function handle to Jacobian of f. Same inputs as f and 
%       returns F(t,x) = df/dx (n x n)
% h     Function handle to (possibly nonlinear) observation function. 
%       Input argument is an (n x 1) state vector so that
%       h(x) = y (m x 1)
% H     Function handle to Jacobian of h.  Same inputs as h and 
%       returns
%       H(x) = dh/dx (m x n)
%
% OUTPUTS
% xhist     (n x N) time history of filter state estimates
% Phist     (n x n x N) time history of filter covariance estimates

n = length(x0); %basic dimension
xminus = x0;
Pminus = P0;
xhist = zeros(n,length(t));
Phist = zeros(n,n,length(t));

%do one update to prime the filter
Hk = H(xminus);
K = (Pminus*Hk.')/(Hk*Pminus*Hk.' + R);
x = xminus + K*(z(1,:).' - h(xminus));
P = (eye(length(Pminus)) - K*Hk)*Pminus;
xhist(:,1) = x;
Phist(:,:,1) = P;

for k = 2:length(t)

    disp(k);
    %propogate
    [~,Zo] = ode45(@hybridEKF_deq,[t(k-1),t(k)],[x;P(:)]);
    xminus = Zo(end,1:n).';
    Pminus = reshape(Zo(end,n+1:end),n,n);
    
    %update
    Hk = H(xminus);
    K = (Pminus*Hk.')/(Hk*Pminus*Hk.' + R);
    x = xminus + K*(z(k,:).' - h(xminus));
    P = (eye(length(Pminus)) - K*Hk)*Pminus;  
    
    %add to history
    xhist(:,k) = x;
    Phist(:,:,k) = P;
end

    %combined x and P differential equations
    function dZ = hybridEKF_deq(T,Z)
        xi = Z(1:n);
        Pi = reshape(Z(n+1:end),n,n);
        dxi = f(T,xi);
        Fi = F(T,xi);
        dPi = Fi*Pi + Pi*Fi.' + Q;
        dZ = [dxi;dPi(:)];
    end
end