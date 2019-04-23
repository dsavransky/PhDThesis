dr = '~/Documents/InstrumentComparison/O3_Jan2012/';
DI22 = analyzeDRM([dr,'DI22_1p5m_Earths_HZ_noBinaries_v2a_1kruns_5000initMass.mat']);
probe75_SuperEarths = analyzeDRM([dr,'probe_p75m_SuperEarths_HZ_noBinaries_v2a_1kruns.mat']);
probe75_Jupiters = analyzeDRM([dr,'probe_p75m_Jupiters_noBinaries_v2a_1kruns.mat']);


outdr = [dr,'SmallScopeRes/'];
names = {'O3_Earths','probe_p75m_SuperEarths','probe_p75m_Jupiters'};
obj = [DI22,probe75_SuperEarths,probe75_Jupiters];
titles = {'O$_3$ - Earth Twins','75 cm Probe - SuperEarths', '75 cm Probe - Jupiters'};

for j=1:length(names)
    plotResDists(obj(j),'AuDETs',[outdr,names{j}],1,true,titles{j})
    plotResDists(obj(j),'ADETs',[outdr,names{j}],3,true,titles{j})
    plotResDists(obj(j),'fullspectra',[outdr,names{j}],5,true,titles{j})
    plotResDists(obj(j),'norbs',[outdr,names{j}],7,true,titles{j})
    plotResDists(obj(j),'slewfuel',[outdr,names{j}],7,true,titles{j})
    plotResDists(obj(j),'skfuel',[outdr,names{j}],7,true,titles{j})
end