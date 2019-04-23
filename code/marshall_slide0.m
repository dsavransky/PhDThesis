rng('shuffle')

nbins = 1000;
Nplanets = 1e6;

%earth
xedges = logspace(-2,log10(0.2025),nbins);
yedges = logspace(-11.5,log10(5e-9),nbins);
H = zeros(nbins,nbins);

arange = [0.7,1.5];
erange = [0,0.35];

p = 0.3;
R = 4.26349283e-5;
term1 = -2.5*log10(R.^2.*p);

%neptune
xedges = logspace(-2,log10(0.2025),nbins);
yedges = logspace(-10.5,log10(5e-8),nbins);
H = zeros(nbins,nbins);

p = 0.3;
R = 0.000165537115;
term1 = -2.5*log10(R.^2.*p);


for j=1:10
    j
    M = 2*pi*rand(Nplanets,1);
    a = arange(1) + (arange(2)-arange(1))*rand(Nplanets,1);
    e = erange(1) + (erange(2)-erange(1))*rand(Nplanets,1);
    d = 5 + (30-5)*rand(Nplanets,1);
    
    %newton-raphson to figure out E:
    counter = 0;
    del = 1;
    E = M;
    while ((del > 1e-8) && (counter <1000))
        E = E - (M - E + e.*sin(E))./(e.*cos(E)-1);
        del = sum(abs(M - (E - e.*sin(E))));
        counter = counter+1;
    end
    %calculate true anomaly
    nu = 2*atan(sqrt((1+e)./(1-e)).*tan(E/2));
    inds = find(nu <0);
    nu(inds) = nu(inds)+2*pi;
    
    %calculate planetary distance
    r = a.*(1-e.*cos(E));
    
    %random psi, theta, phi: rotation about z,x,z axes
    psi = rand(Nplanets,1)*2*pi;
    theta =  acos(rand(Nplanets,1)*2 - 1);
    %theta = rand(Nplanets,1)*2*pi;
    phi = rand(Nplanets,1)*2*pi;
    
    rps_proj = [r.*(cos(nu).*sin(theta).*sin(phi)+sin(nu).*(cos(theta).*cos(psi).*sin(phi)+cos(phi).*sin(psi))),...
        r.*(cos(nu).*cos(phi).*sin(theta)+sin(nu).*(cos(theta).*cos(phi).*cos(psi)-sin(phi).*sin(psi))),...
        r.*(cos(theta).*cos(nu)-cos(psi).*sin(theta).*sin(nu))];
    
    %compute Lambert phase and delta mag
    beta = acos(rps_proj(:,3)./r);
    Phi = (sin(beta)+(pi - beta).*cos(beta))/pi;
    
    %smaple albedos and planetary radii
    dMag = term1 - 2.5*log10(Phi) + 5*log10(r);
    %apparent seperation
    s = sqrt(rps_proj(:,1).^2 + rps_proj(:,2).^2);
    
    C = 10.^(-dMag/2.5);
    alpha = s./d;
    
    h = hist2(alpha,C,xedges,yedges);
    
    H = H + h;
end

buh = H;

%buh = h;
%buh(buh ~= 0) = 1;

buh(buh ~= 0) = log(buh(buh ~= 0));
imagesc(xedges',yedges',buh);
set(gca,'YDir','normal')
ylim([1e-12,1])
xlim([1e-2,1.25])
set(gca,'XScale','log','YScale','log')




%% lines
%Earth
p = 0.3;
R = 4.26349283e-5;
a = 1.;
e = 0.01671123;
d = 10.;
theta = 90*pi/180.;
[alpha1,cexp1] = planContrCurve(a,e,p,R,theta,d);

%rock/iron
R_ri = @(frac,M) (0.0592*frac + 0.0975).*log10(M).^2 + (0.2337*frac + 0.4938).*log10(M) + (0.3102*frac + 0.7932);
R_2E = R_ri(2/3,2)*R;
R_10E = R_ri(2/3,10)*R;
[alpha2,cexp2] = planContrCurve(a,e,p,R_2E,theta,d);
[alpha3,cexp3] = planContrCurve(a,e,p,R_10E,theta,d);

%neptune
pN = 0.4;
RN = 0.000165537115;
[alpha4,cexp4] = planContrCurve(1.,0.016,pN,RN,theta,d); %earth orbit
[alpha5,cexp5] = planContrCurve(30.44125206,0.0011214269,pN,RN,theta,d); %neptune orbit

%Jupiter
aJ = 5.204267;
eJ = 0.048775;
pJ = 0.52;
RJ = 0.000477894503;
[alpha6,cexp6] = planContrCurve(aJ,eJ,pJ,RJ,theta,d);

figure(1)
clf
semilogx(alpha1,cexp1,alpha2,cexp2,alpha3,cexp3,alpha4,cexp4,alpha5,cexp5,alpha6,cexp6,'linewidth',2)
set(gca,'FontName','Times','FontSize',12)
l = legend({'Earth-Twin','2$M_\oplus$, 2/3 Rock, Earth Orbit',...
    '10$M_\oplus$, 2/3 Rock, Earth Orbit','Neptune, Earth Orbit', 'Neptune-Twin','Jupiter-Twin'},...
    'Interpreter','LaTex','Location','best');
lpos = get(l,'Position');
set(l,'Position',lpos+[0,0,0.05,0])
xlabel('Apparent Separation (arcsec)','Interpreter','LaTex')
ylabel('Contrast ($\log_{10}(C)$)','Interpreter','LaTex')
title('Reflected light in V band, $i = 90^\circ$, $d = 10$pc','Interpreter','LaTex')

%print('-dpng','reflected_light_curves.png')

%% blobs

[H,xedges,yedges] = planContrBlob(a,e,p,R);

[HN1,xedgesN1,yedgesN1] = planContrBlob(a,e,pN,RN);

[HN,xedgesN,yedgesN] = planContrBlob(30.44125206,0.0011214269,pN,RN);

[HJ,xedgesJ,yedgesJ] = planContrBlob(aJ,eJ,pJ,RJ);

save 'blobs.mat' H HN1 HN HJ xedges xedgesN xedgesJ xedgesN1 yedges yedgesJ yedgesN yedgesN1

cleanH = @(H) (H.*(log10(H/1e8) > -7))/1e8;

figure(2)
clf
surface(xedges,yedges,log10(cleanH(H)))
shading flat
set(gca,'XScale','log','FontName','Times','FontSize',12)
hold on
surface(xedgesN1,yedgesN1,log10(cleanH(HN1)))
shading flat
surface(xedgesN,yedgesN,log10(cleanH(HN)))
shading flat
surface(xedgesJ,yedgesJ,log10(cleanH(HJ)))
shading flat
cmap = colormap;
cmap(1,:) = 1;
colormap(cmap);
xlabel('Apparent Separation (arcsec)')
ylabel('Contrast (log(C))')

[i,j] = ind2sub(size(H),find(cleanH(H) ~= 0));
text(xedges(max(i)),yedges(max(j)),'Earth-twin','FontName','Times','FontSize',12,'HorizontalAlignment','right')
[i,j] = ind2sub(size(H),find(cleanH(HN1) ~= 0));
text(xedgesN1(max(i)),yedgesN1(max(j)),'Neptune, Earth orbit','FontName','Times','FontSize',12,'HorizontalAlignment','right')
[i,j] = ind2sub(size(H),find(cleanH(HN) ~= 0));
text(xedgesN(max(i)),yedgesN(max(j)),'Neptune-Twin','FontName','Times','FontSize',12,'HorizontalAlignment','right')
[i,j] = ind2sub(size(H),find(cleanH(HJ) ~= 0));
text(xedgesJ(max(i)),yedgesJ(max(j)),'Jupiter-Twin','FontName','Times','FontSize',12,'HorizontalAlignment','right')
title('Reflected light in V band, d \in (1,30)')

print('-dpng','reflected_light_blobs.png')