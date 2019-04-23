afun0 = @(a) a.^(-0.62).*exp(-(a./60).^2);

% gravitational parameter:
G = 6.67428e-11;%m^3/kg/s^2
mAU = 1.495978707e11;%AU in meters
mEarth = 5.9736e24; %kg
G = G/mAU^3*86400^2*mEarth; %AU^3/mEarth/day^2
msun = 333060.402*G;
a85 = ((85/2/pi)^2*msun)^(1/3);

k = 1/quad(afun0,0,1e6);

afun = @(a) k*afun0(a);

f1 = quad(afun,0,a85);

fRa = @(R,a) 1/f1*afun(a).*fR_Kepler(R);

%%
f2 = @(R,a) k.*afun(a).*fR_Kepler(R);

q = quad2d(fRa,1,25,0,1e6); 

Rin = logspace(0,log10(22.6),1000);
ain = logspace(-1,2,1000);
[Rs,as] = meshgrid(Rin,ain);

fs = fRa(Rs,as);

d = fs;
good = find(d ~= 0);
d(good) = log10(d(good));

minexp = floor(min(d(good)));
maxexp = ceil(max(d(good)));
if maxexp-minexp >= 8
    L = minexp:2:maxexp;
else
    L = minexp:maxexp;
end

mn = min(d(good));
rng = max(d(good)) - mn;
d(good) = 1+63*(d(good)-mn)/rng;
l = 1+63*(L-mn)/rng;

figure(1)
clf
imagesc(Rin.',ain.',d);
set(gca,'YDir','normal')
set(gca,'Layer','top')
shading flat
cmap = colormap;
cmap(1,:) = 1;
colormap(cmap);

%draw colorbar
C = colorbar();
set(C,'Ytick',l,'YTickLabel',L,'FontSize',12,'FontName','Times');

set(gca,'FontName','Times','FontSize',14,'Xscale','log','Yscale','log')

xlabel('R (R_\oplus)')
ylabel('a (AU)')