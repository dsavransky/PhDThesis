%instrument comparison

hybrid4b = analyzeDRM('/home/ds/InstrumentComparison/hybrid_4m_3lD_SDO_Earths_HZ_noBinaries_v2a_1kruns_100daysInt.mat');
hybrid4 = analyzeDRM('/home/ds/InstrumentComparison/hybrid_4m_3lD_SDO_Earths_HZ_noBinaries_v2a_1kruns.mat');
hybrid8 = analyzeDRM('/home/ds/InstrumentComparison/hybrid_8m_3lD_SDO_Earths_HZ_noBinaries_v2a_1kruns.mat');

%%
mdo16 = analyzeDRM('/home/ds/InstrumentComparison/mdo_16m_Earths_HZ_noBinaries_v2a_1kruns.mat');
mdo8 = analyzeDRM('/home/ds/InstrumentComparison/mdo_8m_Earths_HZ_noBinaries_v2a_1kruns.mat');
mdo4 = analyzeDRM('~/InstrumentComparison/theia_Earths_HZ_noBinaries_v2a_1kruns.mat');

%%
coron16_3 = analyzeDRM('/home/ds/InstrumentComparison/piaa_16m_3ld_Earths_HZ_noBinaries_v2a_1kruns.mat');
coron4_3 = analyzeDRM('/home/ds/InstrumentComparison/piaa_3ld_Earths_HZ_noBinaries_v2a_1kruns.mat');
coron8_3 = analyzeDRM('/home/ds/InstrumentComparison/piaa_8m_3ld_Earths_HZ_noBinaries_v2a_1kruns.mat');
coron8 = analyzeDRM('/home/ds/InstrumentComparison/piaa_8m_Earths_HZ_noBinaries_v2a_1kruns.mat');
coron4 = analyzeDRM('/home/ds/InstrumentComparison/piaa_Earths_HZ_noBinaries_v2a_1kruns.mat');


%%
plotResDists(mdo4,'AuDETs','~/InstrumentComparison/mdo_4m')
plotResDists(mdo4,'fullspectra','~/InstrumentComparison/mdo_4m')
plotResDists(coron4_3,'AuDETs','~/InstrumentComparison/coron_4m_3lD')
plotResDists(coron4_3,'fullspectra','~/InstrumentComparison/coron_4m_3lD')
plotResDists(hybrid4,'fullspectra','~/InstrumentComparison/hybrid_4m_3lD')
plotResDists(hybrid4,'AuDETs','~/InstrumentComparison/hybrid_4m_3lD')

res = [mdo4,coron4_3,hybrid4];
plotResDistsMult(res,3,'fullspectra',{'MDO','3$\lambda/D$','Hybrid'},'~/InstrumentComparison/compare_4m_eta0.3')
plotResDistsMult(res,3,'AuDETs',{'MDO','3$\lambda/D$','Hybrid'},'~/InstrumentComparison/compare_4m_eta0.3')
plotResDistsMult(res,3,'ADETs',{'MDO','3$\lambda/D$','Hybrid'},'~/InstrumentComparison/compare_4m_eta0.3')


%%
mdo4 = analyzeDRM('/home/ds/InstrumentComparison/theia_Earths_HZ_noBinaries_v2a_1kruns.mat','etas',3);
sdo4 = analyzeDRM('/home/ds/InstrumentComparison/sdo_Earths_HZ_noBinaries_v2a_1kruns_0.3only.mat');
hybrid4 = analyzeDRM('/home/ds/InstrumentComparison/hybrid_4m_3lD_SDO_Earths_HZ_noBinaries_v2a_1kruns.mat','etas',3);
coron4 = analyzeDRM('/home/ds/InstrumentComparison/piaa_Earths_HZ_noBinaries_v2a_1kruns_0.3only.mat');
coron4_3 = analyzeDRM('/home/ds/InstrumentComparison/piaa_3ld_Earths_HZ_noBinaries_v2a_1kruns_0.3only.mat');

%%
coron4_3B = analyzeDRM('/home/ds/InstrumentComparison/piaa_3ld_0.25t_Earths_HZ_noBinaries_v2a_1kruns_0.3only.mat');
mdo4B = analyzeDRM('/home/ds/InstrumentComparison/theia_0.5bp_Earths_HZ_noBinaries_v2a_1kruns_0.3only.mat');
%%
plotResDistsMult([coron4_3,coron4,mdo4,sdo4,hybrid4],1,'AuDETs',{'3 $\lambda/D$','2 $\lambda/D$','MDO','SDO','Hybrid'},'~/InstrumentComparison/compare_4m_eta0.3',1)
plotResDistsMult([coron4_3,coron4,mdo4,sdo4,hybrid4],1,'ADETs',{'3 $\lambda/D$','2 $\lambda/D$','MDO','SDO','Hybrid'},'~/InstrumentComparison/compare_4m_eta0.3',3)
plotResDistsMult([coron4_3,coron4,mdo4,sdo4,hybrid4],1,'fullspectra',{'3 $\lambda/D$','2 $\lambda/D$','MDO','SDO','Hybrid'},'~/InstrumentComparison/compare_4m_eta0.3',2)
plotResDistsMult([coron4_3,coron4,mdo4,sdo4,hybrid4],1,'partspectra',{'3 $\lambda/D$','2 $\lambda/D$','MDO','SDO','Hybrid'},'~/InstrumentComparison/compare_4m_eta0.3',4)


%%
coron8_4 = analyzeDRM('/home/ds/InstrumentComparison/piaa_8m_4ld_Earths_HZ_noBinaries_v2a_1kruns_0.3only.mat');
coron8_3 = analyzeDRM('/home/ds/InstrumentComparison/piaa_8m_3ld_Earths_HZ_noBinaries_v2a_1kruns_0.3only.mat');
coron8 = analyzeDRM('/home/ds/InstrumentComparison/piaa_8m_Earths_HZ_noBinaries_v2a_1kruns_0.3only.mat');
hybrid8 = analyzeDRM('/home/ds/InstrumentComparison/hybrid_8m_3lD_SDO_Earths_HZ_noBinaries_v2a_1kruns_0.3only.mat');
hybrid8_4 = analyzeDRM('/home/ds/InstrumentComparison/hybrid_8m_4lD_SDO_Earths_HZ_noBinaries_v2a_1kruns_0.3only.mat');
mdo8 = analyzeDRM('/home/ds/InstrumentComparison/mdo_8m_Earths_HZ_noBinaries_v2a_1kruns.mat','etas',3);
sdo8 = analyzeDRM('/home/ds/InstrumentComparison/sdo_8m_Earths_HZ_noBinaries_v2a_1kruns_0.3only.mat');
%%
hybrid84 = analyzeDRM('/home/ds/InstrumentComparison/hybrid_8m_3lD_4mSDO_Earths_HZ_noBinaries_v2a_1kruns_0.3only.mat');
hybrid84_4 = analyzeDRM('/home/ds/InstrumentComparison/hybrid_8m_4lD_4mSDO_Earths_HZ_noBinaries_v2a_1kruns_0.3only.mat');
%%
res = [coron8_4,coron8_3,coron8,sdo8,mdo8,hybrid8_4,hybrid8];
names = {'4 $\lambda/D$','3 $\lambda/D$','2 $\lambda/D$','SDO','MDO','Hybrid 4 $\lambda/D$','Hybrid 3 $\lambda/D$'};
fname = '~/InstrumentComparison/compare_8m_eta0.3';
%fname = [];
plotResDistsMult(res,1,'AuDETs',names,fname,1)
plotResDistsMult(res,1,'ADETs',names,fname,3)
plotResDistsMult(res,1,'fullspectra',names,fname,2)
plotResDistsMult(res,1,'partspectra',names,fname,4)

%%
res = [coron8_4,coron8_3,coron8,sdo8,mdo8,hybrid84_4,hybrid84];
names = {'4 $\lambda/D$','3 $\lambda/D$','2 $\lambda/D$','SDO','MDO','Hybrid 4 $\lambda/D$','Hybrid 3 $\lambda/D$'};
fname = '~/InstrumentComparison/compare_8m_4mHybrid_eta0.3';
%fname = [];
plotResDistsMult(res,1,'AuDETs',names,fname,1)
plotResDistsMult(res,1,'ADETs',names,fname,3)
plotResDistsMult(res,1,'fullspectra',names,fname,2)
plotResDistsMult(res,1,'partspectra',names,fname,4)
