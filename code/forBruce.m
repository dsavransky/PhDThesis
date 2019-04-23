as = logspace(-1,2,150);
Ms = logspace(-1,2.9,200);

rocks = Ms < 10;
rmf=0.67;
rads = zeros(size(Ms));
rads(rocks) = (0.0592*rmf+0.0975)*(log10(Ms(rocks)).^2)+(0.2337*rmf+0.4938)*log10(Ms(rocks))+0.3102*rmf+0.7932;
mgiant = [10,17,28,46,77,129,215,318,464,774,1292,2154,3594];
rgiant = [0.29,0.727,0.921,1.004,1.063,1.116,1.141,1.146,1.152,1.158,1.173,1.165,1.134]*11.2;
%from Fortney et al 2007 with the 10Me values being for rocky planets. 
rads(~rocks) = interp1(mgiant,rgiant,Ms(~rocks));
ps = zeros(size(Ms));
ps(rocks) = 0.21;
ps(~rocks) = 0.3;
rads = rads*4.26349283e-5;
term1 = -2.5*log10(rads.^2.*ps);

N = 1e6;
edist = @(e) (e/0.0281).*exp((-e.^2.)/2/0.0281);
e = sampleDist(edist,N,[0,1]);
M = rand(N,1)*2*pi; %mean anomaly
[E,nu] = invKepler(M,e(:)); %eccentric & true anomaly
beta = acos(rand(N,1)*2 - 1);
Phi = (sin(beta)+(pi - beta).*cos(beta))/pi;
term2 = 2.5*log10(Phi);

res = zeros(length(as),length(Ms));
dists = feh.data.stars.DIST;
IWA = 75/1000;

for j = 1:length(as)
    j
    r = as(j).*(1-e.*cos(E));
    s = r.*sin(beta);
    for k = 1:length(Ms)
        dMag = term1(k) - term2 + 5*log10(r);
        tmp = 0;
        for i = 1:length(dists), tmp = tmp + length(find((s >= dists(i)*IWA) & (dMag <= 26))); end
        %tmp = tmp + length(find((s >= max(dists)*IWA) & (dMag <= 26)));
        res(j,k) = tmp;
    end
end





pcolor(as,Ms,log10(res.'/N/length(dists)))
shading flat
%imagesc(Ms,as,log10(res/N/length(dists)))
set(gca,'YScale','log','Xscale','log','YDir','normal','FontSize',14,'FontName','Times')
ylabel('Mass (M_\oplus)')
xlabel('Semi-major axis (AU)')
l = colorbar()
set(get(l,'Title'),'String','log p','FontSize',14,'FontName','Times')