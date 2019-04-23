%ice/rock
R_ir = @(frac,M) (0.0912*frac + 0.1603).*log10(M).^2 + (0.333*frac + 0.7387).*log10(M) + (0.4639*frac + 1.1193);
%rock/iron
R_ri = @(frac,M) (0.0592*frac + 0.0975).*log10(M).^2 + (0.2337*frac + 0.4938).*log10(M) + (0.3102*frac + 0.7932);

fracs = [0,0.5,1];
M = 1:0.01:10;

[f1,m1] = meshgrid(fracs,M);

icerock = R_ir(f1,m1);
rockiron = R_ri(f1,m1);

hold off
plot(M,icerock(:,[3,2,1]),M,rockiron(:,[2,1]),'LineWidth',2)
axis([1,10.01,0.75,3.01])
set(gca,'FontSize',14,'Box','on','FontName','Times','FontWeight','normal')

% ggdat = load('fortney_table4');
% Mn = 10:1:100;
% coremass1 = zeros(length(Mn),1)+0.25;
% coremass2 = zeros(length(Mn),1)+0.5;
% 
% Rn1 = interpn(ggdat.x1,ggdat.x2,ggdat.x3,ggdat.radii,coremass1,Mn.',ones(length(Mn),1)*0.5,'linear',1)*11.209;
% Rn2 = interpn(ggdat.x1,ggdat.x2,ggdat.x3,ggdat.radii,coremass2,Mn.',ones(length(Mn),1)*0.5,'linear',1)*11.209;
% 
% hold on
% semilogx(Mn,[Rn1,Rn2])

legend({'Pure Ice','50\% Rock','Pure Rock', '50\% Iron', 'Pure Iron'},'Location','NorthWest','Interpreter','LaTex')
ylabel('Radius (R$$_\oplus$$)','Interpreter','LaTex')
xlabel('Mass (M$$_\oplus$$)','Interpreter','LaTex')
grid on
