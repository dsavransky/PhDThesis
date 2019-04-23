function [tfound,nfound] = compSim(nweeks,nvisits,projIWA,dMagLim)

Nplanets = 1e6;

arange = [0.7,1.5];
erange = [0,0.35];
p = 0.26;
R = 4.26349283e-5;

%sample from uniform distribution of semi-major axes,eccentricities and
%initial mean anomolies
a = arange(1) + (arange(2)-arange(1))*rand(Nplanets,1);
M = 2*pi*rand(Nplanets,1);
e = erange(1) + (erange(2)-erange(1))*rand(Nplanets,1);

%random psi, theta, phi: rotation about z,x,z axes
psi = rand(Nplanets,1)*2*pi;
theta = rand(Nplanets,1)*pi;
phi = rand(Nplanets,1)*2*pi;

if ~exist('nvisits','var'), nvisits = 52; end
nfound = zeros(nvisits,1);
tfound = zeros(nvisits,1);
found = zeros(Nplanets,1);

for visit = 1:nvisits
    visit
    %newton-raphson to figure out E:
    counter = 0;
    del = 1;
    E = M./(1-e);
    inds = E > sqrt(6*(1-e)./e);
    E(inds) = (6*M(inds)./e(inds)).^(1/3);
    while ((del > eps(2*pi)) && (counter <1000))
        E = E - (M - E + e.*sin(E))./(e.*cos(E)-1);
        del = max(abs(M - (E - e.*sin(E))));
        counter = counter+1;
    end
    
    %calculate true anomaly
    nu = 2*atan(sqrt((1+e)./(1-e)).*tan(E/2));
    inds = find(nu <0);
    nu(inds) = nu(inds)+2*pi;

    %calculate planetary distance
    r = a.*(1-e.*cos(E));
    rps_proj = [r.*(cos(nu).*sin(theta).*sin(phi)+sin(nu).*(cos(theta).*cos(psi).*sin(phi)+cos(phi).*sin(psi))),...
        r.*(cos(nu).*cos(phi).*sin(theta)+sin(nu).*(cos(theta).*cos(phi).*cos(psi)-sin(phi).*sin(psi))),...
        r.*(cos(theta).*cos(nu)-cos(psi).*sin(theta).*sin(nu))];

    %compute Lambert phase and delta mag
    beta = acos(rps_proj(:,3)./r);
    Phi = (sin(beta)+(pi - beta).*cos(beta))/pi;

    term1 = -2.5*log10(R.^2.*p);
    dMag = term1 - 2.5*log10(Phi) + 5*log10(r);

    %apparent seperation
    s = sqrt(rps_proj(:,1).^2 + rps_proj(:,2).^2);
    f = find(dMag < dMagLim & s >= projIWA);
    found(f) = 1;
    nfound(visit) = length(f);
    tfound(visit) = sum(found);
    
    %update M
    if nweeks == 0
    M = M + rand(1)*2*pi;
    else
    M = M + 2*pi*nweeks/52;
    end
end
