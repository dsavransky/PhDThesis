%Rayleigh eccentricity distribution
s = sqrt(2/pi)*0.21; %sigma as func of mean eccentricity
f = @(x) (x/s^2).*exp(-x.^2/2/s^2);
erange = [0,0.5];

%power law semi-major axis
b = 0.61;
arange = [0.7,1.5];
%adist = @(a) (b-1)/arange(1)*(a/arange(1)).^-b;
adist = @(a) b*exp(-b*a);

%%
[H1,xedges1,yedges1] = genComp(1e6,'erange',erange);
[H2,xedges2,yedges2] = genComp(1e6,'erange',erange,'edist',f);
[H3,xedges3,yedges3] = genComp(1e6,'erange',erange,'adist',adist);

for j = 2:100
    j
    [h1] = genComp(1e6,'erange',erange);
    [h2] = genComp(1e6,'erange',erange,'edist',f);
    [h3] = genComp(1e6,'erange',erange,'adist',adist);
    H1 = H1+h1;
    H2 = H2+h2;
    H3 = H3+h3;
end

%%
plotSepvMag(H1/1e8,xedges1,yedges1,false,[],10)
plotSepvMag(H2/1e8,xedges2,yedges2,false,[],11)
plotSepvMag(H3/1e8,xedges3,yedges3,false,[],12)


