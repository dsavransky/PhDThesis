function [a,e,E,I,omega,Omega,P,tau,A,B] = vec2orbElem(rs,vs,mus)
% vec2orbElem(rs,vs,mus) converts positions (rs) and velocities (vs)
% of bodies with gravitational parameters (mus) to Keplerian orbital
% elements.
% 
% INPUTS:
% rs    3n x 1 stacked initial position vectors:
%       [r1(1);r1(2);r1(3);r2(1);r2(2)r2(3);...;rn(1);rn(2);rn(3)]
%       or 3 x n matrix of position vecotrs.
% vs    3n x 1 stacked initial velocity vectors or 3 x n matrix
% mus   gravitational parameters (G*m_i) where G is the
%       gravitational constant and m_i is the mass of the ith body.
%       if all vectors represent the same body, mus may be a
%       scalar.
%
% OUTPUTS:
% a        semi-major axes
% e        eccentricities
% E        eccentric anomalies
% I        inclinations
% omega    arguments of pericenter
% Omega    longitudes of ascending nodes
% P        orbital periods
% tau      time of periapsis crossing
% A, B     orientation matrices
%
%Example:
%  %solar system oribtal elements
%  ssdat = load('solarSystemData.mat');
%  rs = ssdat.p0(1:end-3) - repmat(ssdat.p0(end-2:end),9,1);
%  vs = ssdat.v0(1:end-3) - repmat(ssdat.v0(end-2:end),9,1);
%  mus = ssdat.mus(1:9) + ssdat.mus(10);
%  [a,e,E,I,omega,Omega,P,tau,A,B] = vec2orbElem(rs,vs,mus);
%  %convert back:
%  r = A*diag(cos(E) - e) + B*diag(sin(E));
%  rdot = (-A*diag(sin(E))+B*diag(cos(E)))*...
%          diag(sqrt(mus(:).'.*a.^-3)./(1 - e.*cos(E)));

%condition inputs
mus = mus(:).';
nplanets = numel(rs)/3;
if mod(nplanets,1) ~= 0 || numel(vs) ~= nplanets*3 ||... 
    (length(mus) ~= nplanets && length(mus) ~= 1)
    error('vec2orbElem:inputError',['rs and vs must contain 3n ',...
        'elements and mus must have length n or 1 for n bodies']);
end
if length(rs) == numel(rs)
    rs = reshape(rs,3,nplanets);
end
if length(vs) == numel(vs)
    vs = reshape(vs,3,nplanets);
end

v2s = sum(vs.^2); 
r = sqrt(sum(rs.^2)); %orbital separation
Ws = 0.5*v2s - mus./r;
a = -mus/2./Ws; %semi-major axis

L = [rs(2,:).*vs(3,:) - rs(3,:).*vs(2,:);...
     rs(3,:).*vs(1,:) - rs(1,:).*vs(3,:);...
     rs(1,:).*vs(2,:) - rs(2,:).*vs(1,:)]; %angular momentum
L2s = sum(L.^2);
p = L2s./mus; %semi-latus rectum
e = sqrt(1 - p./a); %eccentricity

%ecentric anomaly
cosE = (1 - r./a)./e;
sinE = sum(rs.*vs)./(e.*sqrt(mus.*a));
E = atan2(sinE,cosE);

if nargout < 4, return; end

%inclination
sinI = sqrt(L(1,:).^2 + L(2,:).^2)./sqrt(L2s);
cosI = L(3,:)./sqrt(L2s);
I = atan2(sinI,cosI);

%argument of pericenter
sinw = ((vs(1,:).*L(2,:) - vs(2,:).*L(1,:))./mus - ...
    rs(3,:)./r)./(e.*sinI);
cosw = ((sqrt(L2s).*vs(3,:))./mus - (L(1,:).*rs(2,:) - ...
    L(2,:).*rs(1,:))./(sqrt(L2s).*r))./(e.*sinI);
omega = atan2(sinw,cosw);

%longitude of ascending node
cosO = -L(2,:)./(sqrt(L2s).*sinI);
sinO = L(1,:)./(sqrt(L2s).*sinI);
Omega = atan2(sinO,cosO);

if nargout < 7, return; end

%orbital periods
P = 2*pi*sqrt(a.^3./mus);

%time of periapsis crossing
tau = -(E - e.*sin(E))./sqrt(mus.*a.^-3);

if nargout < 9, return; end

%A and B vectors
A = [a.*(cos(Omega).*cos(omega) - sin(Omega).*cos(I).*sin(omega));...
     a.*(sin(Omega).*cos(omega) + cos(Omega).*cos(I).*sin(omega));...
     a.*sin(I).*sin(omega)];

B = [-a.*sqrt(1-e.^2).*(cos(Omega).*sin(omega) + ...
                sin(Omega).*cos(I).*cos(omega));...
      a.*sqrt(1-e.^2).*(-sin(Omega).*sin(omega) + ...
                cos(Omega).*cos(I).*cos(omega));...
      a.*sqrt(1-e.^2).*sin(I).*cos(omega)];