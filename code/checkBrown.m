SNR = 5;
dMag = 26;
mag_z = 23;
Rspec = 5;
eta2 = 0.52;
T = 1;
Psi = 0.04;
D = 4;
F0 = 9500*100^2*1e9;
lambda = 550e-9;
mu = 1.5;
q = 2;
T_r = 3600;
DR = 0.001;
sigma_r = 2;
C = 1e-10;


A = pi*D^2/4;
S = T*eta2*F0*10.^(-(V+dMag)/2.5)*lambda/Rspec*A;
N1 = 2*mu*T*eta2*F0*10^(-mag_z/2.5)*lambda^3*(180*3600)^2/(pi*Psi*Rspec*16);
N2 = C*T*eta2*F0*10.^(-V/2.5)*lambda*D^2*pi^2/(Psi*Rspec*64);
N3 = DR/Psi;
N4 = sigma_r^2/Psi/T_r;

N = N1+N2+N3+N4;

t = (SNR./S).^2.*(S+q*N)


Ip = F0*10.^(-(V+dMag)/2.5);
cp = eta2*(lambda/Rspec)*Ip*T*A;
cs = eta2*(lambda/Rspec)*A*F0*10.^(-V/2.5)*(C/2);
Z = N1;
t2 = (SNR./cp.*((1/Psi)*((sigma_r/T_r + DR)...
    *(1+ 1/10) + Z) + cp + cs)).^2