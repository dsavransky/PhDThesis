function [H,xedges,yedges] = planContrBlob(a,e,p,R)

nbins = 1000;
xedges = logspace(-5,log10(a),nbins);
yedges = linspace(floor(log10(p*(R/(a*(1+e)))^2*1e-11)), ceil(log10(p*(R/(a*(1-e)))^2)),nbins);

H = zeros(nbins,nbins);
N = 1e6;
ddist = @(d) -13.8136 + 5.8901*d;

for j=1:100
    j
    E = linspace(0,2*pi,N).';
    theta = acos(rand(N,1)*2 - 1);
    d = sampleDist(ddist,N,[1,30]);
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
    
    h = hist2(alpha,cexp,xedges,yedges);
        
    H = H + h;
end