[res,raw] = analyzeDRM3('theia_Earths_HZ_noBinaries_v2a_baseline_500runs.mat');

%%
x = (1:500).';

cummean = @(y) cumsum(y)./x./mean(y);

figure(1)
clf
subplot(2,1,1)
plot(x,cummean(res.Avisits),x,cummean(res.Auvisits),x,cummean(res.ADETs),x,cummean(res.AuDETs))
grid on
set(gca,'FontName','Times','FontSize',12)
ylabel('Mean','Interpreter','LaTex')
ylim([0.9,1.1])

subplot(2,1,2)
plot(x,cumstd(res.Avisits),x,cumstd(res.Auvisits),x,cumstd(res.ADETs),x,cumstd(res.AuDETs))
grid on
ylim([0.5,1.5])
set(gca,'FontName','Times','FontSize',12)
ylabel('Standard Deviation','Interpreter','LaTex')
xlabel('Ensemble Size','Interpreter','LaTex')
l = legend({'All Visits','Unique Visits','All Detection','Unique Detections'},'Interpreter','LaTex');
pos = get(l,'Position');
set(l,'Position',pos+[0 0 0.05 0])
%%

print('-depsc','../Thesis/figures/burnin.eps');