fid = fopen('~/Downloads/SNRstarshade.csv');
fstring = '%f %f %f %s %f %f %f %f %f %f %f %f';
data = textscan(fid,fstring,'delimiter', ',','headerlines',1);
fclose(fid);

%%
h = 6.62606896e-34; %J*s
c = 299792458; %m/s
k = 1.3806504e-23; %J/K
Temp = 5770; %K
%J/s/m^2/sr/Hz
Idnu = @(nu) (2*h*nu.^3./c.^2)./(exp(h*nu./k/Temp) - 1);

l = data{5}(1)*1e-6;
dl = data{6}(1)*1e-6;

l = 0.55e-6;
dl = 0.089e-6

F = quadgk(Idnu,c/(l+dl/2),c/(l-dl/2))*4*pi

d = 3.08568025e17; %m
F/(4*pi*d.^2)


%%
Pdl = @(lambda) (2*h*c^2./lambda.^5)./(exp(h*c./lambda/k/Temp) - 1)./(h*c./lambda);
PV = quadgk(Pdl,l-dl/2,l+dl/2); %photons/s/m^2/sr

Rsun = 6.955e8; %m
PV*Rsun.^2/d.^2/100^2/(dl/1e-9)/10^(-4.505/2.5)
