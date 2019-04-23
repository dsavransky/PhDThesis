npoints = 1000;
%create orbit with 1 AU semi-major axis
a = 0.5; 
%generate mean anomoly and eccentricity
e = 0.35;
M = linspace(0,2*pi,npoints)';

%newton-raphson to figure out E:
counter = 0;
del = 1;
E = M;
while ((del > 1e-8) & (counter <1000))
    E = E - (M - E + e.*sin(E))./(e.*cos(E)-1);
    del = sum(abs(M - (E - e.*sin(E))));
    counter = counter+1;
end
%calculate planetary distance, semi-minor axis and true anomoly
r = a.*(1-e.*cos(E));
b = a.*sqrt(1-e);
nu = 2*atan(sqrt((1+e)./(1-e)).*tan(E/2));
nup = nu;
buh = find(nup < 0);
nup(buh) = nup(buh)+2*pi;
rps = [zeros(npoints,1), r.*sin(nu), r.*cos(nu)];

%create sphere and set rotations
[cx,cy,cz] = sphere(100);
psi = pi/4;
theta = 70*pi/180;
phi = pi/6;

poi = 900;

Ei = [0,pi,pi/2,-pi/2]';
ri = a.*(1-e.*cos(Ei));
nui = 2*atan(sqrt((1+e)./(1-e)).*tan(Ei/2));
rpsi = [zeros(4,1),ri.*sin(nui),ri.*cos(nui)];

%psi rotation about z axis:
rps_z = [r.*sin(nu).*sin(psi), r.*sin(nu).*cos(psi),r.*cos(nu)];
rpsi_z = [ri.*sin(nui).*sin(psi), ri.*sin(nui).*cos(psi),ri.*cos(nui)];

%theta rotation about x axis:
rps_zx = [r.*sin(nu).*sin(psi), r.*(cos(nu).*sin(theta)+cos(theta).*cos(psi).*sin(nu)),...
    r.*(cos(theta).*cos(nu) - cos(psi).*sin(theta).*sin(nu))];
rpsi_zx = [ri.*sin(nui).*sin(psi), ri.*(cos(nui).*sin(theta)+cos(theta).*cos(psi).*sin(nui)),...
    ri.*(cos(theta).*cos(nui) - cos(psi).*sin(theta).*sin(nui))];

%phi rotation about z axis
rps_zxz = [r.*(cos(nu).*sin(theta).*sin(phi)+sin(nu).*(cos(theta).*cos(psi).*sin(phi)+cos(phi).*sin(psi))),...
            r.*(cos(nu).*cos(phi).*sin(theta)+sin(nu).*(cos(theta).*cos(phi).*cos(psi)-sin(phi).*sin(psi))),...
            r.*(cos(theta).*cos(nu)-cos(psi).*sin(theta).*sin(nu))];
rpsi_zxz = [ri.*(cos(nui).*sin(theta).*sin(phi)+sin(nui).*(cos(theta).*cos(psi).*sin(phi)+cos(phi).*sin(psi))),...
            ri.*(cos(nui).*cos(phi).*sin(theta)+sin(nui).*(cos(theta).*cos(phi).*cos(psi)-sin(phi).*sin(psi))),...
            ri.*(cos(theta).*cos(nui)-cos(psi).*sin(theta).*sin(nui))];
                   

d = 5;

f1 = figure(232);
clf
set(f1,'Position',[15,50,800,700]);
hold on

%plot orbit
plot3(rps_zxz(:,1),rps_zxz(:,2),rps_zxz(:,3),'LineWidth',2)
set(gca,'XTick',[],'YTick',[],'ZTick',[])

%star and planet
l1 = light('Position',[0 0 0],'Style','local');
sun = surface(cx*0.05,cy*0.05,cz*0.05,'FaceLighting','none','FaceColor',[1,1,0],'LineStyle','none');
planet = surface(cx*0.03+rps_zxz(poi,1),cy*0.03+rps_zxz(poi,2),cz*0.03+rps_zxz(poi,3),'FaceColor',[0,0,1],'FaceLighting','phong','LineStyle','none','SpecularExponent',2);

plot3([0 rps_zxz(poi,1)],[0 rps_zxz(poi,2)],[0 rps_zxz(poi,3)],'k')

%orig orbit
ell(1) = fill3(rps(:,1),rps(:,2),rps(:,3),'r');

%plane of sky
ell(2) = fill3(rps_zxz(:,1),rps_zxz(:,2),zeros(npoints,1),'g');
set(ell,'EdgeColor','none','FaceAlpha',0.5,'FaceLighting','none')

%star-planet lines
plot3([0 rps_zxz(poi,1)],[0 rps_zxz(poi,2)],[0 0],'k')
plot3([rps_zxz(poi,1) rps_zxz(poi,1)],[rps_zxz(poi,2) rps_zxz(poi,2)],[0 rps_zxz(poi,3)],'k')

%plot3([0 0], [0 0], [0 -d],'b--','LineWidth',2)
plot3([0 rps_zxz(poi,1)],[0 rps_zxz(poi,2)],[-d rps_zxz(poi,3)],'k')

plot3([0 0],[0 0],[0 0.5],'k')
plot3([0 0],[0 0.5],[0 0],'k')
plot3([0 0.5],[0 0],[0 0],'k')

view(-25,15)
axis equal
xlim([min(rps_zxz(:,1)) max(rps_zxz(:,1))]*1.2)
ylim([min(rps(:,2)) max(rps(:,2))]*1.2)
zlim([min(rps(:,3)) max(rps(:,3))]*1.1)
hold off