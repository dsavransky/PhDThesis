function [E,nu] = invKepler(M,e)
%invKepler - Newton-Raphson inversion of the Kepler Equation
%
%  INPUTS
%  M    n x 1 or 1 x n vector of mean anomalies.
%  e    vector of n eccentricity values, or scalar eccentricity
% 
%  OUTPUTS
%  E    n x 1 vector of eccentric anomalies
%  nu   n x 1 vector of true anomalies
% 
%  Example:
%  %one year of highly eccentric orbit
%  [E, nu] = invKepler(0:pi/100:2*pi,0.8);

%condition inputs
if (numel(M) ~= numel(e)) && numel(e) ~= 1
    error('invKepler:inputError',['For M of length n, e must',...
        'have 1 or n elements']);
end
M = M(:); e = e(:);

%initialize
counter = 0;
del = 1;
E = M./(1-e);
inds = E > sqrt(6*(1-e)./e);
if length(e) == 1
    einds = 1;
else
    einds = inds;
end
E(inds) = (6*M(inds)./e(einds)).^(1/3);

%iterate
while ((del > eps(2*pi)) && (counter <1000))
    E = E - (M - E + e.*sin(E))./(e.*cos(E)-1);
    del = max(abs(M - (E - e.*sin(E))));
    counter = counter+1;
end
if (counter == 1000) %convergence failure
    error('invKepler:overIteration',...
        'Maximum number of iterations exceeded');
end

%calculate true anomaly if needed
if (nargout > 1)
    nu = 2*atan(sqrt((1+e)./(1-e)).*tan(E/2));
    inds = nu < 0;
    nu(inds) = nu(inds)+2*pi;
end