[tfound1,nfound1] = compSim(1,100,0.7,26);
[tfound2,nfound2] = compSim(2,100,0.7,26);
[tfound4,nfound4] = compSim(4,100,0.7,26);
[tfound1d,nfound1d] = compSim(1/7,1000,0.7,26);
[tfoundr,nfoundr] = compSim(0,100,0.7,26);

%%
figure(253)
clf
n = 1:1:100;
pc =  mean([nfound1;nfound2;nfound4])/1e6;
Pn = 1 - (1 - pc).^(1:1000);
meh = semilogx([n,1000],[tfound1;tfound1(end)]/1e6,[n,1000],[tfound2;tfound2(end)]/1e6,[n,1000],[tfound4;tfound4(end)]/1e6,1:1000,tfound1d/1e6,[n,1000],[tfoundr;tfoundr(end)]/1e6,1:1000,Pn,'LineWidth',2);
set(meh(end-2),'Color',[0.4,0.4,0.4])
legend({'1 week','2 weeks','4 weeks','1 day','Random','$$P_n(k>0)$$'},'Interpreter','LaTeX','Location','SouthEast')
grid on
set(gca,'FontName','Times','FontSize',14)
xlabel('$$n$$','Interpreter','LaTeX')
ylabel('Fraction of Planets Found','Interpreter','LaTex')
%%
print('-depsc', '-r0','-painters', '../figures/compSim.eps')
%%
pf0 = get(gcf,'Position');
ylim([0.99,1])
set(gcf,'Position',pf0+[0 0 0 -pf0(end)*0.75])
set(gca,'XTickLabel',[])
xlabel('')
legend off
ylabel('')
set(gcf,'PaperPositionMode','auto')
%%
print('-depsc', '-r0', '-painters', '../figures/compSimZoom.eps')
