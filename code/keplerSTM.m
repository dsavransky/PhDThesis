function Phi = keplerSTM(x0,dt,mu)
% keplerSTM - State transition matrix for Keplerian orbits.
%
% INPUTS
% x0    6n x 1 vector of stacked positions and velocities:
%       [r1_1,r1_2,r1_3,v1_1,v1_2,v1_3,r2_1,r2_2,r2_3,...].'
%       positions and velocities are assumed heliocentric
% dt    Scalar representing propagation time
% mu    n x 1 vector of gravitational parameters (mu_S + mu_Pi)
% 
% OUTPUT
% Phi   State transition matrix such that x = Phi*x0
% 
% Example:
%   %propogate a planet 10 days along its orbit:
%   x0 = [-0.8684, -0.6637, 0.7207, 0.0039, 0.0056, 0.0120].';
%   x1 = keplerSTM(x0,10,2.6935e-04)*x0;

%determine number of planets and validate input
nplanets = length(x0)/6;
if nplanets - floor(nplanets) > 0
    error('keplerSTM:inputError',['The length of x0 ',...
        'must be a multiple of 6.']);
end
if length(mu) ~= nplanets
    error('keplerSTM:inputError',['The length of mu vector ',...
        'must be the length of x0 divided by 6']);
end

%orient mu properly
mu = mu(:).';

%create position and velocity matrices
x0 = reshape(x0,6,nplanets);
r0 = x0(1:3,:);
v0 = x0(4:6,:);

%constants
r0norm = sqrt(sum(r0.^2));
nu0 = dot(r0,v0);
bet = 2*mu./r0norm - dot(v0,v0);

%initialization
u = zeros(1,nplanets);

%For elliptic orbits, calculate period effects
deltaU = zeros(1,length(bet));
eorbs = bet > 0;
if any(eorbs)
    P = 2*pi*mu(eorbs).*bet(eorbs).^(-3/2);
    n = floor((dt + P/2 - 2*nu0(eorbs)./bet(eorbs))./P);
    deltaU(eorbs) = 2*pi*n.*bet(eorbs).^(-5/2);
end

%continued fraction constants
a = 5;
b = 0;
c = 5/2;

%kepler iteration loop
t = zeros(1,nplanets);
counter = 0;
%loop until convergence of the time array to the time step
while max(abs(t-dt)) > max(eps(t))*4 && counter < 1000
    q = bet.*u.^2./(1+bet.*u.^2);
    U0w2 = 1 - 2*q;
    U1w2 = 2*(1-q).*u;
    U = 16/15*U1w2.^5 .* contFrac(a,b,c,q) + deltaU;
    U0 = 2*U0w2.^2-1;
    U1 = 2*U0w2.*U1w2;
    U2 = 2*U1w2.^2;
    U3 = bet.*U + U1.*U2/3;
    r = r0norm.*U0 + nu0.*U1 + mu.*U2;
    t = r0norm.*U1 + nu0.*U2 + mu.*U3;
    u = u - (t-dt)./(4*(1-q).*r);
    counter = counter+1;
end
if counter == 1000
     error('keplerSTM:convError','Failed to converge on t.');
end

%kepler solution
f = 1 - mu./r0norm.*U2;
g = r0norm.*U1 + nu0.*U2;
F = -mu.*U1./r./r0norm;
G = 1 - mu./r.*U2;

%construct the state transition matrix
Phi = zeros(6*nplanets);
for j=1:nplanets
    st = (j-1)*6+1;
    Phi(st:st+5,st:st+5) = [eye(3)*f(j)  eye(3)*g(j);...
                            eye(3)*F(j) eye(3)*G(j)];
end

    %calculate continued fraction
    function G = contFrac(a,b,c,x)
        
        %initialize
        k = 1 - 2*(a-b);
        l = 2*(c-1);
        d = 4*c*(c-1);
        n = 4*b*(c-a);
        A = ones(1,length(x));
        B = ones(1,length(x));
        G = ones(1,length(x));

        Gprev = zeros(1,length(x))+2;
        counter2 = 0;
        %loop until convergence of continued fraction
        while max(abs(G-Gprev)) > max(eps(G)) && counter2 < 1000
            k = -k;
            l = l+2;
            d = d+4*l;
            n = n+(1+k)*l;
            A = d./(d - n*A.*x);
            B = (A-1).*B;
            Gprev = G;
            G = G + B;
            counter2 = counter2+1;
        end
        if counter2 == 1000
            error('keplerSTM:convError','Failed to converge on G.');
        end
    end
end