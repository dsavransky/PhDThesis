function [alpha,cexp] = planContrCurve(a,e,p,R,theta,d)

E = linspace(0,2*pi,1000).';
nu = 2*atan(sqrt((1+e)./(1-e)).*tan(E/2));
r = a.*(1-e.*cos(E));
rps_proj = ([r.*cos(nu),r.*cos(theta).*sin(nu),r.*sin(theta).*sin(nu)]);
term1 = -2.5*log10(R.^2.*p);
beta = acos(rps_proj(:,3)./r);
Phi = (sin(beta)+(pi - beta).*cos(beta))/pi;
dMag = term1 - 2.5*log10(Phi) + 5*log10(r);
s = sqrt(rps_proj(:,1).^2 + rps_proj(:,2).^2);
alpha = s./d;
cexp = -dMag/2.5;