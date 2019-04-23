% fname = '~/Documents/NRO/kepler_rdist.csv';
% dat = csvread(fname,2);
% 
% %%
% figure(1)
% clf
% hold on
% for j=1:length(Rs)-1
%     semilogx([Rs(j),Rs(j+1)],[Rvals85(j),Rvals85(j)])
% end
% 
% %% sanity checks
% %giants - 6 to 22
%  0.0071*(8.0 - 6.0)+0.0102*(11.3-8.0)+ 0.0049*(16.0-11.3)+ 0.0014*(22-16)
%  
%  %large neps - 4 - 6
%  0.0187*(5.7-4) + 0.0071*(6 - 5.7)
%  
%  %small neps - 2 to 4
%  0.1739*(2.8-2) + 0.0609*(4-2.8)
%  
%  %% sample  < 85 day periods
% logfun = @(x,a,b) 1./(x*(log(b)-log(a)));
% 
% X = [];
% N = 1e6;
% for j=1:length(Rs)-1
%     j
%     numsamp = poissrnd(Rvals85(j)*N);
%     sfun = @(x) logfun(x,Rs(j),Rs(j+1));
%     X = [X;sampleDist(sfun,numsamp,[Rs(j),Rs(j+1)])];
% end
% 
% %%
% edges = Rs;
% pstr = 'b';
% [n,x] = histc(X,edges); 
% 
% figure(1)
% clf
% hold on;
% for j=1:length(n)-1
%     patch([edges(j) edges(j) edges(j+1) edges(j+1)], ...
%           [0 n(j) n(j) 0]/N,pstr,'FaceAlpha',0.0);
% end 
% hold off;
% 
% set(gca,'XScale','log','FontName','Times','FontSize',14)
% xlim([min(edges),max(edges)])
% %ylim([0,max(n/sum(n))])

%%
Rs = [1,1.4,2.0,2.8,4.0,5.7,8.0,11.3,16,22.6];
Rvals85 = [0.1555, 0.1671, 0.1739, 0.0609, 0.0187, 0.0071, 0.0102, 0.0049, 0.0014] ;

%%
afun = @(a) a.^(-0.62).*exp(-(a./70).^2);
%afun =  @(a)0.0517.^0.809/gamma(0.809).*(a.^(0.809-1)).*exp(-0.0517*a);

G = 6.67428e-11;%m^3/kg/s^2
mAU = 1.495978707e11;%AU in meters
mEarth = 5.9736e24; %kg
G = G/mAU^3*86400^2*mEarth; %AU^3/mEarth/day^2
msun = 333060.402*G;
a85 = ((85/2/pi)^2*msun)^(1/3);
ap8 = ((0.8/2/pi)^2*msun)^(1/3);
a145 = ((145/2/pi)^2*msun)^(1/3);
a245 = ((245/2/pi)^2*msun)^(1/3);
a418 = ((418/2/pi)^2*msun)^(1/3);

fgiants = quad(afun,0,100)/quad(afun,ap8,a418)*5.24/100;
fgiants85 = quad(afun,0,100)/quad(afun,ap8,a85)*2.0/100;
flneps = quad(afun,0,100)/quad(afun,ap8,a418)*3.18/100; 
fsneps = quad(afun,0,100)/quad(afun,ap8,a245)*30.9/100;
fses = quad(afun,0,100)/quad(afun,ap8,a145)*29.6/100;
inds = Rs >= 4;

fac1 = quad(afun,0,a85);
Rvals = quad(afun,0,100)*(Rvals85./fac1);

%account for longer orbital baseline data
Rvals(6:end) = Rvals(6:end)*2.75%(sum(diff(Rs(6:end)).*Rvals(6:end))/(Rs(end)-Rs(6)));
%% sample  all periods
logfun = @(x,a,b) 1./(x*(log(b)-log(a)));

X = [];
N = 1e6;
for j=1:length(Rs)-1
    j
    numsamp = poissrnd(Rvals(j)*N);
    sfun = @(x) logfun(x,Rs(j),Rs(j+1));
    X = [X;sampleDist(sfun,numsamp,[Rs(j),Rs(j+1)])];
end

%%
%ice/rock
R_ir = @(frac,M) (0.0912*frac + 0.1603).*log10(M).^2 + (0.333*frac + 0.7387).*log10(M) + (0.4639*frac + 1.1193);
M_ir = @(frac,R) 10.^((-(0.333*frac + 0.7387) + sqrt((0.333*frac + 0.7387).^2 - 4*(0.0912*frac + 0.1603).*...
    (0.4639*frac + 1.1193 - R)))./(2*(0.0912*frac + 0.1603)));

%rock/iron
R_ri = @(frac,M) (0.0592*frac + 0.0975).*log10(M).^2 + (0.2337*frac + 0.4938).*log10(M) + (0.3102*frac + 0.7932);
M_ri = @(frac,R) 10.^((-(0.2337*frac + 0.4938) + sqrt((0.2337*frac + 0.4938).^2 - 4*(0.0592*frac + 0.0975).*...
    (0.3102*frac + 0.7932 - R)))./(2*(0.0592*frac + 0.0975)));

ggdat = load('fortney_table4');

rockironfracs = R_ri(linspace(0,1,100),10);
icerockfracs = R_ir(linspace(0,1,100),100);
%%
R = X;
M = zeros(length(R),1);

%first, the things that can be rock/iron or ice/rock (up to 10 M_\oplus)
inds = R <= R_ri(0,10);
Rtmp = R(inds);
Mtmp = zeros(length(Rtmp),1);

% mnfracs = interp1(rockironfracs,linspace(0,1,100),Rtmp);
% mnfracs(mnfracs ~= mnfracs) = 0;
% %make probability of being rock/iron proportional to 1 - mnfrac
% rockiron = rand(length(Rtmp),1) < (1 - mnfracs);
% rifracs = rand(length(find(rockiron)),1).*(1 - mnfracs(rockiron))+mnfracs(rockiron);
% Mtmp(rockiron) = M_ri(rifracs,Rtmp(rockiron));
% icerock = ~rockiron;
% irfracs = rand(length(find(icerock)),1).*(1 - mnfracs(icerock))+mnfracs(icerock);
% Mtmp(icerock) = M_ir(irfracs,Rtmp(icerock));
% 

fracs = rand(length(Rtmp),1)*2-1;
icerock = fracs < 0;
Mtmp(icerock) = M_ir(abs(fracs(icerock)),Rtmp(icerock));
rockiron = fracs >= 0;

mnfracs = interp1(rockironfracs,linspace(0,1,100),Rtmp(rockiron));
bad = fracs(rockiron) < mnfracs;
tmp = fracs(rockiron);
tmp(bad) = tmp(bad).*(1 - mnfracs(bad))+mnfracs(bad);
fracs(rockiron) = tmp;
Mtmp(rockiron) = M_ri(fracs(rockiron),Rtmp(rockiron));

M(inds) = Mtmp;

%next, things that can be icy/watery
inds = (R > R_ri(0,10)) & (R <= min(ggdat.radii(:))*4.77894089e-4/4.26349283e-5);
Rtmp = R(inds);

mnfracs = interp1(icerockfracs,linspace(0,1,100),Rtmp);
mnfracs(mnfracs ~= mnfracs) = 0;
fracs = rand(length(Rtmp),1).*(1 - mnfracs)+mnfracs;
Mtmp = M_ir(fracs,Rtmp);

M(inds) = Mtmp;

% and now things that are gassy
a = sampleDist(afun,length(Rtmp),[0.05,100]);
a(a > 9) = 9;

inds = R > min(ggdat.radii(:))*4.77894089e-4/4.26349283e-5;
Rtmp = R(inds);
Mtmp = zeros(length(Rtmp),1);

radii = ggdat.radii*4.77894089e-4/4.26349283e-5;
radii(radii ~= radii) = 0;
Rtmp(Rtmp > max(radii(:))) = max(radii(:))*0.95;

closein = Rtmp >= max(max(radii(:,:,2)));
a(closein) = sampleDist(afun,length(find(closein)),[0.02,0.05]);

F = TriScatteredInterp(ggdat.x1(:),ggdat.x3(:),radii(:),ggdat.x2(:));
bad = 1:length(Rtmp);   
for j = 1:10
    disp(length(bad));
    coremass = rand(length(bad),1)*100;
    %Mtmp(bad) = griddata(ggdat.x1,ggdat.x3,radii,ggdat.x2,coremass,a(bad),Rtmp(bad));
    Mtmp(bad) = F(coremass,a(bad),Rtmp(bad));
    bad = find(Mtmp ~= Mtmp);
end
Mtmp(bad) = F(rand(length(bad),1)*15+10,a(bad),Rtmp(bad));
bad = find(Mtmp ~= Mtmp);
Mtmp(bad) = 300;

M(inds) = Mtmp;

