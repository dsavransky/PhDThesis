%instrument comparison
dr = '~/Documents/NRO/';
outdr = '~/Documents/NRO/keplerdist/';

%%
case1 = analyzeDRM_v3([dr,'nro_small_kepler_case1.mat']);
case2 = analyzeDRM_v3([dr,'nro_small_kepler_case2.mat']);
case2b = analyzeDRM_v3([dr,'nro_small_kepler_case2b.mat']);
case3 = analyzeDRM_v3([dr,'nro_small_kepler_case3.mat']);
case4 = analyzeDRM_v3([dr,'nro_small_kepler_case4.mat']);

%%
[case1m,mm1,datam1] = analyzeDRM_v3([dr,'nro_small_kepler_mod_case1.mat']);
[case2m,mm2,datam2] = analyzeDRM_v3([dr,'nro_small_kepler_mod_case2b.mat']);
[case3m,mm3,datam3] = analyzeDRM_v3([dr,'nro_small_kepler_mod_case3.mat']);
[case4m,mm4,datam4] = analyzeDRM_v3([dr,'nro_small_kepler_mod_case4.mat']);
%%
group1 = [case1,case2b,case3,case4];
group1_names = {'Case 1','Case 2','Case 3','Case 4'};
outname = [outdr,'keplerdist_hists'];
%outname = [];

plotResDistsMult(group1,1,'AuDETs',group1_names,outname,1,true,false)
plotResDistsMult(group1,1,'ADETs',group1_names,outname,2,true,false)
plotResDistsMult(group1,1,'fullspectra',group1_names,outname,3,true,false)


%%
group1 = [case1,case2m,case3m,case4m];
group1_names = {'Case 1','Case 2','Case 3','Case 4'};
outname = [outdr,'keplerdist_hists2'];
%outname = [];

plotResDistsMult(group1,1,'AuDETs',group1_names,outname,1,true,false)
plotResDistsMult(group1,1,'ADETs',group1_names,outname,2,true,false)
plotResDistsMult(group1,1,'fullspectra',group1_names,outname,3,true,false)

%%
S = load('../NRO/SIMBAD300.mat');
S = S.S;
haveL = (S.L == S.L);
haved = (S.DIST == S.DIST);
pc30 = S.DIST <= 30;

S30 = filtstruct(S,haveL & haved & pc30);

stars = S30;
mainseq =  ((stars.BVNEW < 0.74 & stars.MV < 6*stars.BVNEW + 1.8)...
            | (stars.BVNEW >= 0.74 & stars.BVNEW < 1.37 & stars.MV < 4.3*stars.BVNEW + 3.05)...
            | (stars.BVNEW >= 1.37 & stars.MV < 18*stars.BVNEW - 15.7))...
            & ((stars.BVNEW < 0.87 & stars.MV > -8*(stars.BVNEW - 1.35).^2 + 7.01)...
            | (stars.BVNEW >= 0.87 & stars.BVNEW < 1.45 & stars.MV > 5*stars.BVNEW + 0.81)...
            | (stars.BVNEW >= 1.45 & stars.MV > 18*stars.BVNEW - 18.04));

mults = ~strcmpi(stars.BINARY_CUT,'cut');

S30msnm = filtstruct(S30,mainseq & mults);

%%
%m2 = load([dr,'nro_small_kepler_case1.mat']);
%data = m2.data;
data = datam2;
data.stars = S30msnm;
G = 6.67428e-11;%m^3/kg/s^2
mAU = 1.495978707e11;%AU in meters
mSun = 1.9891e30; %kg
G = G/mAU^3*86400^2*mSun; %AU^3/Msun/day^2
data.Ms_est = ...
    10.^(0.002456*data.stars.MV.^2-0.09711*...
    data.stars.MV+0.4365)*G;
err = (rand(length(data.stars.MV),1)*2 - 1)*0.07;
data.Ms_true = data.Ms_est + err.*data.Ms_est;

universe = genUniverse_Kepler_mod(data,1);

rs = universe.rv(:,1:3).';
vs = universe.rv(:,4:6).';
[a,e,E,I,omega,Omega,P,tau,A,B] = vec2orbElem(rs(:),vs(:),universe.Mp+...
        data.Ms_true(universe.planInds));

universe.a = a.';
universe.r = sqrt(sum(rs.^2)).';
universe.s = sqrt(rs(1,:).^2 + rs(2,:).^2).';
        
beta = acos(rs(3,:).'./universe.r);
Phi = (sin(beta)+(pi - beta).*cos(beta))/pi;
term1 = -2.5*log10(universe.Rp.^2.*universe.p);
universe.dMag = term1 - 2.5*log10(Phi) + 5*log10(universe.r);

G = 6.67428e-11;%m^3/kg/s^2
mAU = 1.495978707e11;%AU in meters
mEarth = 5.9736e24; %kg
G = G/mAU^3*86400^2*mEarth; %AU^3/mEarth/day^2
universe.Mp = universe.Mp/G;
universe.Rp = universe.Rp/4.26349283e-5;

universe.s = universe.s./data.stars.DIST(universe.planInds);

%%
cs = colormap();

pntcs = log10(universe.Mp);
mn = min(pntcs);
mx = max(pntcs);

minexp = floor(mn);
maxexp = ceil(mx);
if maxexp-minexp >= 8
    L = minexp:2:maxexp;
else
    L = minexp:maxexp;
end 
rng = mx - mn;
pntcs = 1+63*(pntcs-mn)/rng;
l = 1+63*(L-mn)/rng;


[~,inds] = sort(universe.Mp);
figure(1)
clf
loglog(universe.s(1),10.^(-universe.dMag(1)/2.5),'o','MarkerSize',log10(universe.Rp(1))*8+5,'MarkerFaceColor',interp1(1:64,cs,pntcs(1)),'MarkerEdgeColor','k')
hold on
for j = inds(end:-1:1).'
    loglog(universe.s(j),10.^(-universe.dMag(j)/2.5),'o','MarkerSize',log10(universe.Rp(j))*8+5,'MarkerFaceColor',interp1(1:64,cs,pntcs(j)),'MarkerEdgeColor','k')
end
plot([0.103,0.103],[1e-9,1],'k--','LineWidth',2)
plot([0.103,100],[1e-9,1e-9],'k--','LineWidth',2)
plot([0.2063,0.2063],[1e-9,1],'k','LineWidth',2)
plot([0.2063,100],[1e-9,1e-9],'k','LineWidth',2)
plot([0.2063,0.2063],[1e-8,1],'k-.','LineWidth',2)
plot([0.2063,100],[1e-8,1e-8],'k-.','LineWidth',2)
plot([0.2063,0.2063],[1e-10,1],'k:','LineWidth',2)
plot([0.2063,100],[1e-10,1e-10],'k:','LineWidth',2)
hold off

C = colorbar();
set(C,'Ytick',l,'YTickLabel',L,'FontSize',12,'FontName','Times');
set(get(C,'Title'),'String','Mass (log10 M_\oplus)','FontName','Times','FontSize',12)
set(gca,'FontName','Times','FontSize',16)
xlabel('Separation (arcsec)')
ylabel('Contrast')




%%
tags = {'planInds','Mp','Rp','p','I','a','r','s','dMag','rtype'};

fid = fopen([outdr,'sample_planets_keplerdist2.txt'],'w');
len = 0;
fprintf(fid,'%s','# ');
for k = 1:length(tags)
    len = len + length(tags{k});
    if k < length(tags)
        fprintf(fid,'%s\t',tags{k});
    else
        fprintf(fid,'%s\n',tags{k});
    end
end
brk = '#';
for j = 1:len+5*length(tags),brk = [brk,'-'];end
fprintf(fid,'%s\n',brk);

for j = 1:length(universe.(tags{1}))
    for k = 1:length(tags)
        if strcmp(class(universe.(tags{k})(1)),'double')
            tmp = num2str(universe.(tags{k})(j));
        else
            tmp = char(universe.(tags{k})(j));
        end
        if k < length(tags) 
            fprintf(fid,'%s,\t',tmp);
        else
             fprintf(fid,'%s\n',tmp);
        end
    end
end

fclose(fid);


%%

[Msin1,Msout1,Rsin1,Rsout1,asin1,asout1] = getInOutDists(mm1,datam1);
[Msin2,Msout2,Rsin2,Rsout2,asin2,asout2] = getInOutDists(mm2,datam2);
[Msin3,Msout3,Rsin3,Rsout3,asin3,asout3] = getInOutDists(mm3,datam3);
[Msin4,Msout4,Rsin4,Rsout4,asin4,asout4] = getInOutDists(mm4,datam4);

%%
cs = [   0         0    1.0000
    1.0000         0         0
         0         0    0.2069
    1.0000         0    0.6897
         0    0.3103         0
    0.3828    0.7241    0.1034
         0    0.6207    1.0000
    0.6207    0.3103    0.2759
    0.4483    0.2069    0.7241
    0.0690    0.4828    0.5172
    1.0000    0.6897    1.0000
    1.0000    0.8276         0
    0.8276    0.4483         0
    0.7586    0.6897    0.4138
    0.8966         0    1.0000
    0.1379    0.1034         0
    0.9310    0.0345    0.3448
    0.4828         0    0.3448
    0.4828    0.4138    0.5862
    0.3448    0.6552    0.4483
    0.4138    0.3103         0
    0.6207         0         0
    1.0000    0.6897    0.6897
         0    0.2759    0.6207
    0.9310    0.4483    1.0000];  

%%
G = 6.67428e-11;%m^3/kg/s^2
mAU = 1.495978707e11;%AU in meters
mEarth = 5.9736e24; %kg
G = G/mAU^3*86400^2*mEarth; %AU^3/mEarth/day^2

nin = loghist([Msin1;Msin2;Msin3;Msin4]/G,[],1,'k--');
hold on
n1 = loghist(Msout1/G,[],[],'o-',cs(1,:));
n2 = loghist(Msout2/G,[],[],'s-',cs(2,:));
n3 = loghist(Msout3/G,[],[],'d-',cs(3,:));
n4 = loghist(Msout4/G,[],[],'^-',cs(4,:));
hold off
set(gca,'FontName','Times','FontSize',16)
xlabel('Mass (M$_\oplus$)','Interpreter','Latex')
ylabel('Normalized Frequency','Interpreter','Latex')
legend({'Input','Case 1','Case 2','Case 3','Case 4'},'Interpreter','Latex','FontSize',14)
set(gca,'XTick',[10,100,1000])
ylim([0,0.09])
%print('-depsc',[outdr,'massdists.eps'])
%%
nin = loghist([Rsin1;Rsin2;Rsin3;Rsin4]/4.26349283e-5,[],2,'k--');
hold on
n1 = loghist(Rsout1/4.26349283e-5,[],[],'o-',cs(1,:));
n2 = loghist(Rsout2/4.26349283e-5,[],[],'s-',cs(2,:));
n3 = loghist(Rsout3/4.26349283e-5,[],[],'d-',cs(3,:));
n4 = loghist(Rsout4/4.26349283e-5,[],[],'^-',cs(4,:));
hold off
set(gca,'FontName','Times','FontSize',16)
xlabel('Radius(R$_\oplus$)','Interpreter','Latex')
ylabel('Normalized Frequency','Interpreter','Latex')
ylim([0,0.1])
legend({'Input','Case 1','Case 2','Case 3','Case 4'},'Interpreter','Latex','FontSize',14)

%print('-depsc',[outdr,'radiusdists.eps'])

%%
nin = loghist([asin1;asin2;asin3;asin4],[],3,'k--');
hold on
n1 = loghist(asout1,[],[],'o-',cs(1,:));
n2 = loghist(asout2,[],[],'s-',cs(2,:));
n3 = loghist(asout3,[],[],'d-',cs(3,:));
n4 = loghist(asout4,[],[],'^-',cs(4,:));
hold off
set(gca,'FontName','Times','FontSize',16)
xlabel('Semi-major Axis (AU)','Interpreter','Latex')
ylabel('Normalized Frequency','Interpreter','Latex')
legend({'Input','Case 1','Case 2','Case 3','Case 4'},'Interpreter','Latex','FontSize',14)
ylim([0,0.07])
xlim([0.2,30])

%print('-depsc',[outdr,'smadists.eps'])