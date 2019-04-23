%instrument comparison
dr = '~/Documents/InstrumentComparison/';
%d = dir([dr,'*.mat']]);
%d = {d.name};

outdr = '~/Documents/Conferences and Talks/AAS219/';

hybrid4 = analyzeDRM([dr,'hybrid_4m_3lD_SDO_Earths_HZ_noBinaries_v2a_1kruns.mat']);
hybrid83 = analyzeDRM([dr,'hybrid_8m_3lD_SDO_Earths_HZ_noBinaries_v2a_1kruns.mat']);
hybrid83_4 = analyzeDRM([dr,'hybrid_8m_3lD_4mSDO_Earths_HZ_noBinaries_v2a_1kruns.mat']);
hybrid84 = analyzeDRM([dr,'hybrid_8m_4lD_SDO_Earths_HZ_noBinaries_v2a_1kruns.mat']);
hybrid84_4 = analyzeDRM([dr,'hybrid_8m_4lD_4mSDO_Earths_HZ_noBinaries_v2a_1kruns.mat']);

hybrid4_giants = analyzeDRM([dr,'hybrid_4m_3lD_SDO_Giants_noBinaries_v2a_1kruns.mat']);
%%
plotResDists(hybrid4,'AuDETs',[outdr,'hybrid4m_3lD'],1,false,[0,65,0,0.25])
plotResDists(hybrid4,'ADETs',[outdr,'hybrid4m_3lD'],2,false)
plotResDists(hybrid4,'fullspectra',[outdr,'hybrid4m_3lD'],3,true,[0,60,0,0.25])

plotResDists(hybrid4_giants,'AuDETs',[outdr,'hybrid4m_3lD_Giants'],4,false,[0,65,0,0.2])
plotResDists(hybrid4_giants,'ADETs',[outdr,'hybrid4m_3lD_Giants'],5,false)
plotResDists(hybrid4_giants,'fullspectra',[outdr,'hybrid4m_3lD_Giants'],6,false,[0,60,0,0.2])

%%
plotResDists(hybrid83,'AuDETs',[outdr,'hybrid8m_3lD'],7)
plotResDists(hybrid83,'ADETs',[outdr,'hybrid8m_3lD'],8,false)
plotResDists(hybrid83,'fullspectra',[outdr,'hybrid8m_3lD'],9,false)


%%
group = [hybrid4,hybrid83,hybrid83_4,hybrid84,hybrid84_4];
group_names = {'4m SDO, 3$^\lambda\hspace{-0.5ex}/\hspace{-0.5ex}_D$','8m SDO, 3$^\lambda\hspace{-0.5ex}/\hspace{-0.5ex}_D$','8m w/4m SDO, 3$^\lambda\hspace{-0.5ex}/\hspace{-0.5ex}_D$','8m SDO, 4$^\lambda\hspace{-0.5ex}/\hspace{-0.5ex}_D$','8m w/4m SDO, 4$^\lambda\hspace{-0.5ex}/\hspace{-0.5ex}_D$'};

plotResDistsMult(group,3,'AuDETs',group_names,[outdr,'compare_hybrids'],1)
plotResDistsMult(group,3,'ADETs',group_names,[outdr,'compare_hybrids'],2,false)
plotResDistsMult(group,3,'fullspectra',group_names,[outdr,'compare_hybrids'],3,false)


%%
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
group1 = [coron4_2,coron4_3,coron4_4,mdo4,sdo4,hybrid4];
group1_names = {'2$^\lambda\hspace{-0.5ex}/\hspace{-0.5ex}_D$','3$^\lambda\hspace{-0.5ex}/\hspace{-0.5ex}_D$','4$^\lambda\hspace{-0.5ex}/\hspace{-0.5ex}_D$','MDO','SDO','Hybrid'};
outname = [outdr,'compare_4m'];

plotResDistsMult(group1,3,'AuDETs',group1_names,outname,1,true,false,'4 m')
plotResDistsMult(group1,3,'ADETs',group1_names,outname,2,false,false,'4 m')
plotResDistsMult(group1,3,'fullspectra',group1_names,outname,3,false,false,'4 m')
plotResDistsMult(group1(1:4),3,'partspectra',group1_names(1:4),outname,4,true,true,'4 m',[0.7098    0.5429    0.1893    0.3179])

%%
group2 = [coron8_2,coron8_3,coron8_4,mdo8,sdo8,hybrid83,hybrid84,hybrid83_4,hybrid84_4];
group2_names = {'2$^\lambda\hspace{-0.5ex}/\hspace{-0.5ex}_D$','3$^\lambda\hspace{-0.5ex}/\hspace{-0.5ex}_D$','4 $^\lambda\hspace{-0.5ex}/\hspace{-0.5ex}_D$','MDO','SDO','Hybrid 3 $^\lambda\hspace{-0.5ex}/\hspace{-0.5ex}_D$','Hybrid 4$^\lambda\hspace{-0.5ex}/\hspace{-0.5ex}_D$','8m/4m 3$^\lambda\hspace{-0.5ex}/\hspace{-0.5ex}_D$','8m/4m 4$^\lambda\hspace{-0.5ex}/\hspace{-0.5ex}_D$'};
outname = [outdr,'compare_8m'];

plotResDistsMult(group2,3,'AuDETs',group2_names,outname,1,true,false,'8 m',[0.7123    0.2905    0.2839    0.7095])
plotResDistsMult(group2,3,'ADETs',group2_names,outname,2,false,false,'8 m')
plotResDistsMult(group2,3,'fullspectra',group2_names,outname,3,false,false,'8 m')
plotResDistsMult(group2(1:4),3,'partspectra',group2_names(1:4),outname,4,false,true,'8 m')


%%
plotResDists(mdo4,'AuDETs','~/InstrumentComparison/mdo_4m')
plotResDists(mdo4,'fullspectra','~/InstrumentComparison/mdo_4m')
plotResDists(coron4_3,'AuDETs','~/InstrumentComparison/coron_4m_3lD')
plotResDists(coron4_3,'fullspectra','~/InstrumentComparison/coron_4m_3lD')
plotResDists(hybrid4,'fullspectra','~/InstrumentComparison/hybrid_4m_3lD')
plotResDists(hybrid4,'AuDETs','~/InstrumentComparison/hybrid_4m_3lD')

res = [mdo4,coron4_3,hybrid4];
plotResDistsMult(res,3,'fullspectra',{'MDO','3$^\lambda\hspace{-0.5ex}/\hspace{-0.5ex}_D$','Hybrid'},'~/InstrumentComparison/compare_4m_eta0.3')
plotResDistsMult(res,3,'AuDETs',{'MDO','3$^\lambda\hspace{-0.5ex}/\hspace{-0.5ex}_D$','Hybrid'},'~/InstrumentComparison/compare_4m_eta0.3')
plotResDistsMult(res,3,'ADETs',{'MDO','3$^\lambda\hspace{-0.5ex}/\hspace{-0.5ex}_D$','Hybrid'},'~/InstrumentComparison/compare_4m_eta0.3')

%%
plotResDistsMult([coron4_3,coron4,mdo4,sdo4,hybrid4],1,'AuDETs',{'3 $^\lambda\hspace{-0.5ex}/\hspace{-0.5ex}_D$','2 $^\lambda\hspace{-0.5ex}/\hspace{-0.5ex}_D$','MDO','SDO','Hybrid'},'~/InstrumentComparison/compare_4m_eta0.3',1)
plotResDistsMult([coron4_3,coron4,mdo4,sdo4,hybrid4],1,'ADETs',{'3 $^\lambda\hspace{-0.5ex}/\hspace{-0.5ex}_D$','2 $^\lambda\hspace{-0.5ex}/\hspace{-0.5ex}_D$','MDO','SDO','Hybrid'},'~/InstrumentComparison/compare_4m_eta0.3',3)
plotResDistsMult([coron4_3,coron4,mdo4,sdo4,hybrid4],1,'fullspectra',{'3 $^\lambda\hspace{-0.5ex}/\hspace{-0.5ex}_D$','2 $^\lambda\hspace{-0.5ex}/\hspace{-0.5ex}_D$','MDO','SDO','Hybrid'},'~/InstrumentComparison/compare_4m_eta0.3',2)
plotResDistsMult([coron4_3,coron4,mdo4,sdo4,hybrid4],1,'partspectra',{'3 $^\lambda\hspace{-0.5ex}/\hspace{-0.5ex}_D$','2 $^\lambda\hspace{-0.5ex}/\hspace{-0.5ex}_D$','MDO','SDO','Hybrid'},'~/InstrumentComparison/compare_4m_eta0.3',4)


%%
res = [coron8_4,coron8_3,coron8,sdo8,mdo8,hybrid8_4,hybrid8];
names = {'4 $^\lambda\hspace{-0.5ex}/\hspace{-0.5ex}_D$','3 $^\lambda\hspace{-0.5ex}/\hspace{-0.5ex}_D$','2 $^\lambda\hspace{-0.5ex}/\hspace{-0.5ex}_D$','SDO','MDO','Hybrid 4 $^\lambda\hspace{-0.5ex}/\hspace{-0.5ex}_D$','Hybrid 3 $^\lambda\hspace{-0.5ex}/\hspace{-0.5ex}_D$'};
fname = '~/InstrumentComparison/compare_8m_eta0.3';
%fname = [];
plotResDistsMult(res,1,'AuDETs',names,fname,1)
plotResDistsMult(res,1,'ADETs',names,fname,3)
plotResDistsMult(res,1,'fullspectra',names,fname,2)
plotResDistsMult(res,1,'partspectra',names,fname,4)

%%
res = [coron8_4,coron8_3,coron8,sdo8,mdo8,hybrid84_4,hybrid84];
names = {'4 $^\lambda\hspace{-0.5ex}/\hspace{-0.5ex}_D$','3 $^\lambda\hspace{-0.5ex}/\hspace{-0.5ex}_D$','2 $^\lambda\hspace{-0.5ex}/\hspace{-0.5ex}_D$','SDO','MDO','Hybrid 4 $^\lambda\hspace{-0.5ex}/\hspace{-0.5ex}_D$','Hybrid 3 $^\lambda\hspace{-0.5ex}/\hspace{-0.5ex}_D$'};
fname = '~/InstrumentComparison/compare_8m_4mHybrid_eta0.3';
%fname = [];
plotResDistsMult(res,1,'AuDETs',names,fname,1)
plotResDistsMult(res,1,'ADETs',names,fname,3)
plotResDistsMult(res,1,'fullspectra',names,fname,2)
plotResDistsMult(res,1,'partspectra',names,fname,4)

%%
res = coron8_4;
fname = '~/InstrumentComparison/coron_8m_4lD';
plotResDists(res,'AuDETs',fname)
plotResDists(res,'fullspectra',fname)
plotResDists(res,'ADETs',fname)