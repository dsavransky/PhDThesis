plan_fname = 'allPlans.csv';
[headers,units,data] = readIPAC(plan_fname);

star_fname = 'allStars.csv';
[Sheaders,Sunits,Sdata] = readIPAC(star_fname);

%resort by discovery year
year = strmatch('plndisc',headers);
[~,ord] = sort(str2double(data(:,year)));
data = data(ord,:);

%get breakdown by years
[years,inds] = unique(str2double(data(:,year)),'last');

%fill in missing years
tmp = [];
for j = 2:length(years)
    if years(j) - years(j - 1) > 1
        tmp = [tmp,[years(j-1)+1:years(j)-1;ones(1,(years(j) - years(j-1)-1))*inds(j-1)]];
    end
end
years = [years;tmp(1,:).'];
inds = [inds;tmp(2,:).'];
[years,tmp] = sort(years);
inds = inds(tmp);

%find coords and identifiers
ra = strmatch('ra',Sheaders);
dec = strmatch('dec',Sheaders);
pid = strmatch('uniqueid',headers);
sid = strmatch('uniqueid',Sheaders);
method = strmatch('plnmethod',headers);

pprops = zeros(3,1);
mass = strmatch('plnmassj',headers, 'exact');
smax = strmatch('plnorbsmax',headers, 'exact');
eccen = strmatch('plnorbeccen',headers, 'exact');
%%
figure(1)
clf
colormap('gray')
pos = get(1,'Position');
set(1,'Position',pos+[0,0,0,pos(4)*0.25]);
masses = str2double(data(:,mass));
[nm1,xm1] = hist(masses(masses<5),25);
[nm2,xm2] = hist(masses(masses>=5 & masses <=25),25);
subplot(2,1,1)
plot(0,0,'w.','Visible','off')
hold on
bar(xm1,nm1,'EdgeColor','k','FaceColor','w','LineWidth',1)
hold off
set(gca,'FontName','Times','FontSize',16)
ylim([0,max(nm1)+1])
grid on
legend({'$$m_P < 5$$ M$$_J$$'},'Interpreter','LaTeX')
subplot(2,1,2)
plot(5,0,'w.','Visible','off')
hold on
bar(xm2,nm2,'EdgeColor','k','FaceColor','w','LineWidth',1)
hold off
set(gca,'FontName','Times','FontSize',16)
ylim([0,max(nm2)+1])
grid on
legend({'$$m_P \ge 5$$ M$$_J$$'},'Interpreter','LaTeX')
xlabel('$$m_P$$ (M$$_J$$)','Interpreter','LaTeX')
%%
print('-depsc', '-painters','../figures/currMassHist.eps');

%%
figure(2)
clf
colormap('gray')
pos = get(2,'Position');
set(2,'Position',pos+[0,0,0,pos(4)*0.25]);
as = str2double(data(:,smax));
as = as(as == as);
[na1,xa1] = hist(as(as<0.1),25);
[na2,xa2] = hist(as(as>=0.1 & as < 6),25);
subplot(2,1,1)
plot(pi/100,pi/100,'w.')
hold on
bar(xa1,na1,'EdgeColor','k','FaceColor','w','LineWidth',1)
hold off
set(gca,'FontName','Times','FontSize',16)
ylim([0,max(na1)+1])
grid on
legend({'$$a < 0.1$$ AU'},'Interpreter','LaTeX')
subplot(2,1,2)
plot(pi/10,pi/10,'w.')
hold on
bar(xa2,na2,'EdgeColor','k','FaceColor','w','LineWidth',1)
hold off
set(gca,'FontName','Times','FontSize',16)
ylim([0,max(na2)+1])
grid on
legend({'$$a \ge 0.1$$ AU'},'Interpreter','LaTeX')
xlabel('$$a$$ (AU)','Interpreter','LaTeX')

%%
print('-depsc', '-painters','../figures/currSMaxHist.eps');

%%
figure(3)
clf
colormap('gray')
es = str2double(data(:,eccen));
es = es(es == es);
[ne1,xe1] = hist(es,25);
bar(xe1,ne1,'EdgeColor','k','FaceColor','w','LineWidth',1)
set(gca,'FontName','Times','FontSize',14)
xlabel('$$e$$','Interpreter','LaTeX')
grid on
%%
print('-depsc', '-painters','../figures/currEccenHist.eps');

%%
figure(4)
clf
colormap('gray')
bar(years,inds,'EdgeColor','k','FaceColor',[211,211,211]/255,'LineWidth',1)
set(gca,'FontName','Times','FontSize',14)
xlabel('Year','Interpreter','LaTeX')
ylabel('Known Exoplanets','Interpreter','LaTeX')
xlim([min(years),max(years)+1])
ylim([0,550])
grid on
%%
print('-depsc', '-painters','../figures/planDiscHist.eps');

