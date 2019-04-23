plan_fname = 'all_planets.csv';
star_fname = 'all_stars_w_planets.csv';

[headers,units,data] = readIPAC(plan_fname);
[Sheaders,Sunits,Sdata] = readIPAC(star_fname);

%resort by discovery year
year = find(strcmp('plndisc',headers));
[~,ord] = sort(str2double(data(:,year)));
data = data(ord,:);

%get breakdown by years
[years,inds] = unique(str2double(data(:,year)),'first');

%find coords and identifiers
ra = find(strcmp('ra',Sheaders));
dec = find(strcmp('dec',Sheaders));
pid = find(strcmp('uniqueid',headers));
sid = find(strcmp('uniqueid',Sheaders));
method = find(strcmp('plnmethod',headers));

M = find(strcmp('plnmassj',headers));
a = find(strcmp('plnorbsmax',headers));
e = find(strcmp('plnorbeccen',headers));

%find unique methods
mthds = cellstr(unique(char(data{:,method}),'rows'));
blanks = strcmp('',mthds);
mthds = mthds(~blanks);
lgnd = mthds;
if any(blanks)
    mthds{end+1} = '';
    lgnd{end+1} = 'Other'; 

end

for j = 1:length(lgnd)
    switch lgnd{j}
        case 'ima'
            lgnd{j} = 'Imaging';
        case 'micro'
            lgnd{j} = 'Microlensing';
        case 'rv'
            lgnd{j} = 'Radial Velocity';
        case 'tran'
            lgnd{j} = 'Transits';
        case 'pul'
            lgnd{j} = 'Pulsar Timing';
        case 'ast'
            lgnd{j} = 'Astrometry';
    end
end

%define colors and symbols
cmap = colormap('lines');
%cmap = cmap(1:length(mthds),:);
cmap = unique(cmap,'rows');
%cmap = cmap([1:3,7,5,8:9],:);
syms = 'p^dvos<';

%generate relative sizes
szs = log10(str2double(data(:,M)));
szs = round((szs - min(szs)).^1.8+6);


%%
%create figure
h = figure(337);
clf
set(h,'Position',[70,200,800,400]);

ax1 = subplot(1,1,1);
ax1pos = get(ax1,'Position');
ax1pos = ax1pos + [-0.12,0,0,0];
set(ax1,'Position',ax1pos);

%create legend
hold on
for i = 1:length(mthds)
    plot(0,0,syms(i),'Color',cmap(i,:),'MarkerSize',8,'Visible','off');
end
set(gca,'FontSize',14,'FontName','Times','Box','off','XTick',[],'YTick',[])
l = legend(lgnd);
lpos = get(l,'Position');
set(l,'Position',[ax1pos(1) + ax1pos(3)+0.005,ax1pos(2) + ax1pos(4)/2 - lpos(4)/2,lpos(3:4)])

%create skymap
axesm hammer;
framem; gridm; mlabel;
setm(gca,'MLabelParallel','equator', 'MLabelLocation',-120:60:120,'FontName','Times','FontSize',14)
axlim = axis();
tx = text(axlim(1)+0.25,axlim(3)+0.25,num2str(years(3)),'FontSize',16);
tx2 = text(axlim(2)-0.5,axlim(3)+0.25,num2str(4),'FontSize',16);

msizes = zeros(length(Sdata),1)+50;
%plot initial pionts up to 1994
for j = 1:4
    sind = find(strcmp(data(j,pid),Sdata(:,sid)));
    mthd = find(strcmp(data{j,method},mthds));
    scatterm(str2double(Sdata(sind,dec)),str2double(Sdata(sind,ra)),...
        msizes(sind),cmap(mthd,:),syms(mthd));
    msizes(sind) = msizes(sind)*1.5^2;
end
j = j+1;

%%
doMov = false;

%loop through years
if doMov
    mov = avifile('exoplanetHistory.avi');
    set(gcf,'Renderer','zbuffer')
    mov = addframe(mov,getframe(gcf));
end

fpyFull = 20;
currMonth = clock;
currMonth = currMonth(2);
xmxs = [20,10,1];
for k = 3:length(inds)
    %if normal year, use total fpy, otherwise, truncate
    if k ~= length(inds)
        fpy = fpyFull;
    else
        fpy = floor(fpyFull*currMonth/12);
    end
    
    %figure out how many planets will be in this year
    if k < length(inds)
        nplans = inds(k+1) - inds(k);
    else
        nplans = length(data) - inds(k) + 1;
    end
    
    if nplans == 1
        frames = round(fpyFull/2);
    elseif nplans == 2;
        frames = round([fpyFull/3,2*fpyFull/3]);
    else
        frames = round(linspace(1,fpy,nplans));
    end

    for nframe = 1:fpy
        if k ~= length(inds)
            set(tx,'String',num2str(max([years(k),years(k)+(nframe/fpy)-0.1]),'%4.1f'));
        else
            set(tx,'String',num2str(years(k)+(nframe*(currMonth - 1)/12/fpy),'%4.1f'));
        end
        while numel(frames) > 0 && frames(1) == nframe
            sind = find(strcmp(data(j,pid),Sdata(:,sid)));
            mthd = find(strcmp(data{j,method},mthds));
            scatterm(str2double(Sdata(sind,dec)),str2double(Sdata(sind,ra)),...
                msizes(sind),cmap(mthd,:),syms(mthd));
            msizes(sind) = msizes(sind)+50;
            frames = frames(2:end);
            set(tx2,'String',num2str(j));

            j = j+1;
        end
        if doMov,mov = addframe(mov,getframe(gcf));else pause(1/15); end
    end
end

if doMov
for i = 1:3,mov = addframe(mov,getframe(gcf));end
mov = close(mov);
end


%%
for k = 1:2
figure(k)
clf
for j = 1:length(mthds)
if k == 1,loglog(0,0,'.','Color',cmap(j,:));else semilogx(0,0,'.','Color',cmap(j,:));end
if j == 1, hold on, end
end
set(gca,'FontName','Times','FontSize',16)

for j=1:length(data)
    mthd = find(strcmp(data{j,method},mthds));
    if szs(j) == szs(j)
        if k==1
        loglog(str2double(data(j,a)),str2double(data(j,M)),'.','Color',cmap(mthd,:),'MarkerSize',20);
        else
            semilogx(str2double(data(j,a)),str2double(data(j,e)),'.','Color',cmap(mthd,:),'MarkerSize',szs(j));
        end
    end
end
hold off
if k==1
    legend(lgnd,'Location','SouthEast','FontSize',14);
    ylim([1e-3,1e2])
    ylabel('Mass (M_J)')
else
    legend(lgnd,'Location','NorthEast','FontSize',14);
    ylabel('Eccentricity')
end
xlabel('Semi-major axis (AU)')
xlim([1e-3,1e3])
hold off

if k==1
    print('-depsc','M_v_a.eps')
else
    print('-depsc','e_v_a.eps')
end
end

%%
mass = M;
smax = a;
eccen = e;

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
print('-depsc','currMassHist.eps');

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
print('-depsc','currSMaxHist.eps');

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
print('-depsc','currEccenHist.eps');

