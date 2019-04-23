%instrument comparison
dr = '~/Documents/InstrumentComparison/';
odr = '~/Documents/InstrumentComparison/ExopagFigs/';
%%
hybrid43 = analyzeDRM([dr,'hybrid_4m_3lD_SDO_Earths_HZ_noBinaries_v2b_1kruns.mat']);
hybrid44 = analyzeDRM([dr,'hybrid_4m_4lD_SDO_Earths_HZ_noBinaries_v2b_1kruns.mat']);
hybrid83 = analyzeDRM([dr,'hybrid_8m_3lD_SDO_Earths_HZ_noBinaries_v2b_1kruns.mat']);
hybrid83_4 = analyzeDRM([dr,'hybrid_8m_3lD_4mSDO_Earths_HZ_noBinaries_v2b_1kruns.mat']);
hybrid84 = analyzeDRM([dr,'hybrid_8m_4lD_SDO_Earths_HZ_noBinaries_v2b_1kruns.mat']);
hybrid84_4 = analyzeDRM([dr,'hybrid_8m_4lD_4mSDO_Earths_HZ_noBinaries_v2b_1kruns.mat']);
hybrid43_giants = analyzeDRM([dr,'hybrid_4m_3lD_SDO_Giants_noBinaries_v2b_1kruns.mat']);

mdo16 = analyzeDRM([dr,'mdo_16m_Earths_HZ_noBinaries_v2a_1kruns.mat']);
mdo8 = analyzeDRM([dr,'mdo_8m_Earths_HZ_noBinaries_v2a_1kruns.mat']);
mdo4 = analyzeDRM([dr,'mdo_4m_Earths_HZ_noBinaries_v2a_1kruns.mat']);

sdo8 = analyzeDRM([dr,'sdo_8m_Earths_HZ_noBinaries_v2a_1kruns.mat']);
sdo4 = analyzeDRM([dr,'sdo_4m_Earths_HZ_noBinaries_v2a_1kruns.mat']);

%%
coron16_3 = analyzeDRM([dr,'piaa_16m_3lD_Earths_HZ_noBinaries_v2a_1kruns.mat']);
coron16_4 = analyzeDRM([dr,'piaa_16m_4lD_Earths_HZ_noBinaries_v2a_1kruns.mat']);

coron4_2 = analyzeDRM([dr,'piaa_4m_2lD_Earths_HZ_noBinaries_v2a_1kruns.mat']);
coron4_3 = analyzeDRM([dr,'piaa_4m_3lD_Earths_HZ_noBinaries_v2a_1kruns.mat']);
coron4_4 = analyzeDRM([dr,'piaa_4m_4lD_Earths_HZ_noBinaries_v2a_1kruns.mat']);

coron8_2 = analyzeDRM([dr,'piaa_8m_2lD_Earths_HZ_noBinaries_v2a_1kruns.mat']);
coron8_3 = analyzeDRM([dr,'piaa_8m_3lD_Earths_HZ_noBinaries_v2a_1kruns.mat']);
coron8_4 = analyzeDRM([dr,'piaa_8m_4lD_Earths_HZ_noBinaries_v2a_1kruns.mat']);

%%
outdr = [odr,'Instruments/'];
names = {'hybrid_4m_3lD','hybrid_4m_4lD','hybrid_4m_3lD_Giants',...
    'mdo_4m','sdo_4m',...
    'coron_4m_2lD','coron_4m_3lD','coron_4m_4lD',...
    'coron_8m_2lD','coron_8m_3lD','coron_8m_4lD',...
    'hybrid_8m_3lD','hybrid_8m_4lD',...
    'hybrid_8m4m_3lD','hybrid_8m4m_4lD'};

obj = [hybrid43,hybrid44,hybrid43_giants,...
    mdo4,sdo4,...
    coron4_2,coron4_3,coron4_4,...
    coron8_2,coron8_3,coron8_4,...
    hybrid83,hybrid84,...
    hybrid83_4,hybrid84_4];

ld = '$^\lambda\hspace{-0.5ex}/\hspace{-0.5ex}_D$';
titles = {['4m Hybrid, 3',ld],['4m Hybrid, 4',ld],['4m Hybrid, 3',ld,', Giants'],...
    '4m MDO','4m SDO',...
    ['4m Coronagraph, 2',ld],['4m Coronagraph, 3',ld],['4m Coronagraph, 4',ld],...
    ['8m Coronagraph, 2',ld],['8m Coronagraph, 3',ld],['8m Coronagraph, 4',ld],...
    ['8m Hybrid, 3',ld],['8m Hybrid, 4',ld],...
    ['8m/4m Hybrid, 3',ld],['8m/4m Hybrid, 4',ld],};

for j=1:length(names)
    plotResDists(obj(j),'AuDETs',[outdr,names{j}],1,true,titles{j})
    plotResDists(obj(j),'ADETs',[outdr,names{j}],3,true,titles{j})
    plotResDists(obj(j),'fullspectra',[outdr,names{j}],5,true,titles{j})
    plotResDists(obj(j),'norbs',[outdr,names{j}],7,true,titles{j})
end

%%
outdr = [odr,'Comparisons/'];
%%
group = [hybrid43,hybrid44,hybrid83,hybrid83_4,hybrid84,hybrid84_4];
group_names = {['4m SDO, 3',ld],['4m SDO, 4',ld],['8m SDO, 3',ld],['8m w/4m SDO, 3',ld],['8m SDO, 4',ld],['8m w/4m SDO, 4',ld]};
outname = 'compare_hybrids';

plotResDistsMult(group,3,'AuDETs',group_names,[outdr,outname],1,true,false,'Hybrids')
plotResDistsMult(group,3,'ADETs',group_names,[outdr,outname],2,true,false,'Hybrids')
plotResDistsMult(group,3,'fullspectra',group_names,[outdr,outname],3,true,false,'Hybrids')
plotResDistsMult(group,3,'norbs',group_names,[outdr,outname],4,true,false,'Hybrids')

%%
group = [coron4_2,coron4_3,coron4_4,mdo4,sdo4,hybrid43,hybrid44];
group_names = {['2',ld],['3',ld],['4',ld],'MDO','SDO',['Hybrid, 3',ld],['Hybrid, 4',ld]};
outname = 'compare_4m';

plotResDistsMult(group,3,'AuDETs',group_names,[outdr,outname],1,true,false,'4m Designs')
plotResDistsMult(group,3,'ADETs',group_names,[outdr,outname],2,true,false,'4m Designs')
plotResDistsMult(group,3,'fullspectra',group_names,[outdr,outname],3,true,false,'4m Designs')
plotResDistsMult(group,3,'norbs',group_names,[outdr,outname],4,true,false,'4m Designs')
plotResDistsMult(group(1:4),3,'partspectra',group_names(1:4),[outdr,outname],5,true,false,'4m Designs')

%%
group = [coron8_2,coron8_3,coron8_4,mdo8,sdo8,hybrid83,hybrid84,hybrid83_4,hybrid84_4];
group_names = {['2',ld],['3',ld],['4 ',ld],'MDO','SDO',['Hybrid 3 ',ld],['Hybrid 4',ld],['8m/4m 3',ld],['8m/4m 4',ld]};
outname = 'compare_8m';

plotResDistsMult(group,3,'AuDETs',group_names,[outdr,outname],1,true,false,'8m Designs')
plotResDistsMult(group,3,'ADETs',group_names,[outdr,outname],2,true,false,'8m Designs')
plotResDistsMult(group,3,'fullspectra',group_names,[outdr,outname],3,true,false,'8m Designs')
plotResDistsMult(group,3,'norbs',group_names,[outdr,outname],4,true,false,'8m Designs')
plotResDistsMult(group(1:4),3,'partspectra',group_names(1:4),[outdr,outname],5,true,false,'8m Designs')

%%
group = [coron16_3,coron16_4,mdo16];
group_names = {['3',ld],['4 ',ld],'MDO'};
outname = 'compare_16m';

plotResDistsMult(group,3,'AuDETs',group_names,[outdr,outname],1,true,false,'16m Designs')
plotResDistsMult(group,3,'ADETs',group_names,[outdr,outname],2,true,false,'16m Designs')
plotResDistsMult(group,3,'fullspectra',group_names,[outdr,outname],3,true,false,'16m Designs')
plotResDistsMult(group,3,'norbs',group_names,[outdr,outname],4,true,false,'16m Designs')
plotResDistsMult(group,3,'partspectra',group_names,[outdr,outname],5,true,false,'16m Designs')

%%
outdr = [odr,'varMissionPortions/'];
dr = '~/Documents/InstrumentComparison/varMissionPortion/';
%%
basen = 'piaa_4m_3lD_Earths_HZ_noBinaries_v2a_1kruns_eta0.3_mp';
res = [];
for j=1:10
    tmp = analyzeDRM([dr,basen,num2str(j/10,'%1.1f'),'.mat']);
    res = [res,tmp];
end

group = res;
outname = 'coron_4m_3lD';

tt = ['4m Coronagraph, 3',ld,', $\eta_{\oplus} = 0.3$'];

plotResDists2(group,'AuDETs',[outdr,outname],1,true,tt)
plotResDists2(group,'ADETs',[outdr,outname],3,true,tt)
plotResDists2(group,'fullspectra',[outdr,outname],5,true,tt)
plotResDists2(group,'norbs',[outdr,outname],7,true,tt)
plotResDists2(group,'Auvisits',[outdr,outname],9,true,tt)

%%
tmp = load([dr,'piaa_4m_2lD_Earths_HZ_noBinaries_v2a_1kruns_eta0.3_varmp.mat']);
data = tmp.data;
all = tmp.res.universes;
tmp2.eta = tmp.res.eta;
out = [];
for j=1:10
    tmp2.universes = all{j};
    tmp = analyzeDRM(data,tmp2);
    out = [out,tmp];
end

group = out;
outname = 'coron_4m_2lD';
tt = ['4m Coronagraph, 2',ld,', $\eta_{\oplus} = 0.3$'];

plotResDists2(group,'AuDETs',[outdr,outname],1,true,tt)
plotResDists2(group,'ADETs',[outdr,outname],3,true,tt)
plotResDists2(group,'fullspectra',[outdr,outname],5,true,tt)
plotResDists2(group,'norbs',[outdr,outname],7,true,tt)
plotResDists2(group,'Auvisits',[outdr,outname],9,true,tt)

%%
for k = 2:4
    tmp = load([dr,'piaa_8m_',num2str(k),'lD_Earths_HZ_noBinaries_v2a_1kruns_eta0.3_varmp.mat']);
    data = tmp.data;
    all = tmp.res.universes;
    tmp2.eta = tmp.res.eta;
    out = [];
    for j=1:10
        tmp2.universes = all{j};
        tmp = analyzeDRM(data,tmp2);
        out = [out,tmp];
    end
    
    group = out;
    outname = ['coron_8m_',num2str(k),'lD'];
    tt = ['8m Coronagraph, ',num2str(k),ld,', $\eta_{\oplus} = 0.3$'];
    
    plotResDists2(group,'AuDETs',[outdr,outname],1,true,tt)
    plotResDists2(group,'ADETs',[outdr,outname],3,true,tt)
    plotResDists2(group,'fullspectra',[outdr,outname],5,true,tt)
    plotResDists2(group,'norbs',[outdr,outname],7,true,tt)
    plotResDists2(group,'Auvisits',[outdr,outname],9,true,tt)
end