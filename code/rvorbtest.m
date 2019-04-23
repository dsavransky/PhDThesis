%% generate orbits
etas = [0.92,0.7,0.32];
Mranges = {[1,10],[10,100],[100,3000]};
Mfun_giants = @(x) (x/317.8).^(-1.3);
Mfun_neps =  @(x) (x/317.8).^(-1.3);
Mfun_ses = @(x) 1./(x*(log(10)-log(1)));
Mfuns = {Mfun_ses,Mfun_neps,Mfun_giants};

nPlanSys = 1e6;
types = cell(3,1);
lens = [0,cumsum(round(etas/sum(etas)*nPlanSys))];
for j=1:length(types)
    types{j} = lens(j)+1:lens(j+1);
end

% Rayleigh eccen dist
sedist = sqrt(2/pi)*0.21; %sigma as func of mean eccentricity
edist = @(x) (x/sedist^2).*exp(-x.^2/2/sedist^2);

% gravitational parameter:
G = 6.67428e-11;%m^3/kg/s^2
mAU = 1.495978707e11;%AU in meters
mEarth = 5.9736e24; %kg
G = G/mAU^3*86400^2*mEarth; %AU^3/mEarth/day^2

data = struct(...
    'adist',@(a) a.^(-0.62).*exp(-(a./70).^2),...
    'a_min',0.1,...
    'a_max',100,...
    'p_min',0.2,...
    'p_max',0.5,...
    'edist',edist,...
    'e_min',0,...
    'e_max',1);

%generate parameter distributions
a = sampleDist(data.adist,nPlanSys,[data.a_min,data.a_max]);
p = data.p_min + (data.p_max - data.p_min)*rand(nPlanSys,1);
e = sampleDist(data.edist,nPlanSys,[data.e_min,data.e_max]);

I = acos(rand(nPlanSys,1)*2 - 1);
Omega = rand(nPlanSys,1)*2*pi;
omega = rand(nPlanSys,1)*2*pi;
MeanAnom = 2*pi*rand(nPlanSys,1);

[EccAnom,nu] = invKepler(MeanAnom,e);

%figure out radii and Masses
%ice/rock
R_ir = @(frac,M) (0.0912*frac + 0.1603).*log10(M).^2 + (0.333*frac + 0.7387).*log10(M) + (0.4639*frac + 1.1193);
%rock/iron
R_ri = @(frac,M) (0.0592*frac + 0.0975).*log10(M).^2 + (0.2337*frac + 0.4938).*log10(M) + (0.3102*frac + 0.7932);

%load Fortney et al. table 4
ggdat = load('fortney_table4');

M = zeros(nPlanSys,1);
R = zeros(nPlanSys,1);
rtype = zeros(nPlanSys,1);

for j=1:3
    M(types{j}) = sampleDist(Mfuns{j},length(types{j}),Mranges{j});
end
  
%SuperEarth radii
Mtmp = M(types{1});
Rtmp = zeros(length(Mtmp),1);
fracs = rand(length(Rtmp),1)*2 - 1;
icerock = fracs < 0;
Rtmp(icerock) = R_ir(abs(fracs(icerock)),Mtmp(icerock));
rockiron = fracs >= 0;
Rtmp(rockiron) = R_ri(fracs(rockiron),Mtmp(rockiron));
R(types{1}) = Rtmp*4.26349283e-5; %in AU
rtype(types{1}) = fracs;
        
%Neptune radii
atmp = a;
atmp(atmp > 9) = 9;
Mtmp = M;
Mtmp(Mtmp < 17) = 17;
corerat = ggdat.x1./ggdat.x2;

coremass = rand(length(Mtmp(types{2})),1)*max(corerat((ggdat.x2 < 100)& (ggdat.radii == ggdat.radii)));
coremass = Mtmp(types{2}).*coremass;
coremass(Mtmp(types{2})< 28) = rand(length(find(Mtmp(types{2}) < 28)),1)*10;
coremass(Mtmp(types{2})< 46 & Mtmp(types{2}) > 28) = rand(length(find(Mtmp(types{2}) < 46 & Mtmp(types{2}) > 28)),1)*25;
coremass(Mtmp(types{2})< 77 & Mtmp(types{2}) > 46) = rand(length(find(Mtmp(types{2}) < 77 & Mtmp(types{2}) > 46)),1)*25;
coremass(Mtmp(types{2})< 129 & Mtmp(types{2}) > 77) = rand(length(find(Mtmp(types{2}) < 129 & Mtmp(types{2}) > 77)),1)*50;
R(types{2}) = interpn(ggdat.x1,ggdat.x2,ggdat.x3,ggdat.radii,coremass,Mtmp(types{2}),atmp(types{2}),'linear',1)*4.77894089e-4; %in AU

%giant radii
coremass = rand(length(Mtmp(types{3})),1)*100;
coremass(Mtmp(types{3})< 129) = rand(length(find(Mtmp(types{3}) < 129)),1)*50;
R(types{3}) = interpn(ggdat.x1,ggdat.x2,ggdat.x3,ggdat.radii,coremass,Mtmp(types{3}),atmp(types{3}),'linear',1)*4.77894089e-4; %in AU

%gravitational parameters & period
mus = G*(M + 333060.402);
P = 2*pi*sqrt(a.^3./mus);
term1 = -2.5*log10(R.^2.*p);

A = [a.*(cos(Omega).*cos(omega) - sin(Omega).*cos(I).*sin(omega)),...
     a.*(sin(Omega).*cos(omega) + cos(Omega).*cos(I).*sin(omega)),...
     a.*sin(I).*sin(omega)];

B = [-a.*sqrt(1-e.^2).*(cos(Omega).*sin(omega) + sin(Omega).*cos(I).*cos(omega)),...
      a.*sqrt(1-e.^2).*(-sin(Omega).*sin(omega) + cos(Omega).*cos(I).*cos(omega)),...
      a.*sqrt(1-e.^2).*sin(I).*cos(omega)];

dMagLim = 26;
projSepLim = 100/1000*10;
%%
%figure(1)
%clf
%hold on
%set(1,'Visible','off','Renderer','zbuffer')

MeanAnom = linspace(0,2*pi,1000);
%visorbs = [];
%maxacvis = [];
%maxacvis = 0;
%zeroacvis = [];
%maxvelvis = [];
%minvelvis = [];
%zerovelvis = [];

%1 percent period error
di = 0.05*length(MeanAnom);
maxvelpmdivis = [];
minvelpmdivis = [];
zerovelpmdivis = [];


visorbsrv = unique([onlymax;onlymin;onlyzero]).';

for j = visorbsrv%1:length(e) 

if mod(j,100) == 0, disp(j); end

[EccAnom,nu] = invKepler(MeanAnom,e(j));

n = sqrt(mus(j)/a(j)^3);
Eccdot = n./(1 - e(j)*cos(EccAnom));
Eccddot = -e(j)*n^2*sin(EccAnom)./(1 - e(j)*cos(EccAnom)).^3;

r = A(j,:).'*(cos(EccAnom) - e(j)).' + B(j,:).'*(sin(EccAnom)).';
v = A(j,:).'*(-sin(EccAnom).*Eccdot).' + B(j,:).'*(cos(EccAnom).*Eccdot).';
%ac = A(j,:).'*(-sin(EccAnom).*Eccddot - cos(EccAnom).*Eccdot.^2).' + B(j,:).'*(cos(EccAnom).*Eccddot - sin(EccAnom).*Eccdot.^2).';

rmag = sqrt(sum(r.^2,1));
%vmag = sqrt(sum(v.^2,1));

beta = acos(r(3,:)./rmag);
Phi = (sin(beta)+(pi - beta).*cos(beta))/pi;
dMag = term1(j) - 2.5*log10(Phi) + 5*log10(rmag);
s = sqrt(r(1,:).^2 + r(2,:).^2);

vis = (s > projSepLim) & (dMag < dMagLim);

if ~any(vis)
    disp('bink')
    continue
end

%which of the max/mins are visible at +/- 1%
if any(onlymax == j)
    [~,ind] = max(v(3,:));
    indp = ind+di;
    if indp > length(vis), indp = indp - length(vis); end
    indm = ind - di;
    if indm < 1, indm = indm+length(vis); end
    if vis(indp) && vis(indm),maxvelpmdivis = [maxvelpmdivis,j];end
end

if any(onlymin == j)
    [~,ind] = min(v(3,:));
    indp = ind+di;
    if indp > length(vis), indp = indp - length(vis); end
    indm = ind - di;
    if indm < 1, indm = indm+length(vis); end
    if vis(indp) && vis(indm),minvelpmdivis = [minvelpmdivis,j];end
end

if any(onlyzero == j)
try
    zvind1 = interp1(v(3,1:length(MeanAnom)/2),MeanAnom(1:length(MeanAnom)/2),0);
    [~,ind] = min(abs(MeanAnom - zvind1));
    if vis(ind)
        indp = ind+di;
        if indp > length(vis), indp = indp - length(vis); end
        indm = ind - di;
        if vis(indp) && vis(indm)
            zerovelpmdivis = [zerovelpmdivis,j];
            continue
        end
    end
catch
end
try
    zvind2 = interp1(v(3,length(MeanAnom)/2:end),MeanAnom(length(MeanAnom)/2:end),0);
    [~,ind] = min(abs(MeanAnom - zvind2));
    if vis(ind)
        indp = ind+di;
        if indp > length(vis), indp = indp - length(vis); end
        indm = ind - di;
        if vis(indp) && vis(indm)
            zerovelpmdivis = [zerovelpmdivis,j];
            continue
        end
    end
catch
end
end

% % calculate visibility at max,min, zero velocity
% [~,ind] = min(v(3,:));
% if vis(ind),minvelvis = [minvelvis,j];end
% [~,ind] = max(v(3,:));
% if vis(ind),maxvelvis = [maxvelvis,j];end
% 
% [~,mvind] = min(v(3,:));
% try
%     zvind1 = interp1(v(3,1:length(MeanAnom)/2),MeanAnom(1:length(MeanAnom)/2),0);
%     [~,ind] = min(abs(MeanAnom - zvind1));
%     if vis(ind)
%         zerovelvis = [zerovelvis,j];
%         continue
%     end
% catch
% end
% try
%     zvind2 = interp1(v(3,length(MeanAnom)/2:end),MeanAnom(length(MeanAnom)/2:end),0);
%     [~,ind] = min(abs(MeanAnom - zvind2));
%     if vis(ind) zerovelvis = [zerovelvis,j];end
% catch
% end

% % which orbits are ever visible
%visorbs = [visorbs,j];

% %location of minimum acceleration
%[~,ind] = min(ac(3,:));
%if vis(ind),maxacvis = [maxacvis,j];end

% %location of zero acceleration
% [~,ind] = min(abs(ac(3,:)));
% if vis(ind),zeroacvis = [zeroacvis,j];end


% plot3(r(1,:),r(2,:),r(3,:))
% r2 = r;
% r2(:,~vis) = nan;
% plot3(r2(1,:),r2(2,:),r2(3,:),'r')

% if ~vis(ind)
% 
% figure(1)
% clf
% subplot(2,1,1)
% set(gca,'FontName','Times','FontSize',12)
% plot(MeanAnom,r(1,:),'b',MeanAnom,r(2,:),'b')
% ylabel('Position (AU)')
% xlim([0,2*pi])
% r2 = r;
% r2(:,~vis) = nan;
% hold on
% plot(MeanAnom,r2(1,:),'r',MeanAnom,r2(2,:),'r','linewidth',2)
% hold off
% subplot(2,1,2)
% set(gca,'FontName','Times','FontSize',12)
% plot(MeanAnom,ac(3,:),'b')
% ylabel('Acceleration (AU/day^2)')
% xlabel('Mean Anomaly (rad)')
% xlim([0,2*pi])
% 
% pause
% end
% pause 

%[~,zind] = min(abs(ac(3,:)));
%[~,mind] = min(ac(3,:));


% figure(33)
% clf
% subplot(3,1,1)
% set(gca,'FontName','Times','FontSize',12)
% plot(MeanAnom,r(1,:),'r',MeanAnom,r(2,:),'b', MeanAnom,r(3,:),'g')
% ylabel('Position (AU)')
% xlim([0,2*pi])
% r2 = r;
% r2(:,~vis) = nan;
% hold on
% plot(MeanAnom,r2(1,:),'r',MeanAnom,r2(2,:),'b',MeanAnom,r2(3,:),'g','linewidth',2)
% plot([MeanAnom(mind),MeanAnom(mind)],get(gca,'ylim'),'k--')
% plot([MeanAnom(zind),MeanAnom(zind)],get(gca,'ylim'),'k--')
% hold off
% 
% subplot(3,1,2)
% set(gca,'FontName','Times','FontSize',12)
% plot(MeanAnom,v(1,:),'r',MeanAnom,v(2,:),'b', MeanAnom,v(3,:),'g')
% hold on
% plot([MeanAnom(mind),MeanAnom(mind)],get(gca,'ylim'),'k--')
% plot([MeanAnom(zind),MeanAnom(zind)],get(gca,'ylim'),'k--')
% hold off
% ylabel('Velocity (AU/day^2)')
% xlabel('Mean Anomaly (rad)')
% xlim([0,2*pi])
% 
% subplot(3,1,3)
% set(gca,'FontName','Times','FontSize',12)
% plot(MeanAnom,ac(1,:),'r',MeanAnom,ac(2,:),'b', MeanAnom,ac(3,:),'g')
% ylabel('Acceleration (AU/day^2)')
% xlabel('Mean Anomaly (rad)')
% xlim([0,2*pi])
% hold on
% plot([MeanAnom(mind),MeanAnom(mind)],get(gca,'ylim'),'k--')
% plot([MeanAnom(zind),MeanAnom(zind)],get(gca,'ylim'),'k--')
% hold off
% 
% pause


end


%%
tmp = e;
[n1,x1] = hist(tmp(maxvelvis));
[n2,x2] = hist(tmp(minvelvis),x1);
[n3,x3] = hist(tmp(zerovelvis),x1);
bar(x1,[n1/sum(n1);n2/sum(n2);n3/sum(n3)].')

%%
tmp = e;
[n1,x1] = hist(tmp(onlymax));
[n2,x2] = hist(tmp(onlymin),x1);
[n3,x3] = hist(tmp(onlyzero),x1);
bar(x1,[n1/sum(n1);n2/sum(n2);n3/sum(n3)].')

set(gca,'FontName','Times','FontSize',14)
legend({'Max','Min','Zero'})
xlim([0,max(tmp)])
ylabel('Normalized Frequency')

xlabel('Argument of Periapsis (rad)')
%%
tmp = zeros(nPlanSys,1);
tmp(maxvelvis) = 1;
tmp(minvelvis) = 0;
tmp(zerovelvis) = 0;
onlymax = find(tmp);

tmp = zeros(nPlanSys,1);
tmp(minvelvis) = 1;
tmp(maxvelvis) = 0;
tmp(zerovelvis) = 0;
onlymin = find(tmp);

tmp = zeros(nPlanSys,1);
tmp(zerovelvis) = 1;
tmp(minvelvis) = 0;
tmp(maxvelvis) = 0;
onlyzero = find(tmp);

%%
tmp1 = zeros(nPlanSys,1);
tmp1(maxvelvis) = 1;
tmp2 = zeros(nPlanSys,1);
tmp2(minvelvis) = 1;
tmp3 = zeros(nPlanSys,1);
tmp3(zerovelvis) = 1;