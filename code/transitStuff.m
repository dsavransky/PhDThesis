
stream = RandStream('mt19937ar','Seed',sum(101*clock));
trans = [];
Rs = 0.00464912633; %AU
N = 1e6;
fa = @(a) a.^-0.62;
amin = 0.1;
amax = 100;
%fa = @(a) 1/(amax-amin);

s = [];
for j=1:10
    j
a = sampleDist(fa,N,[amin,amax]);
e = rand(stream,N,1);
M = 2*pi*rand(stream,N,1);
[E,nu] = invKepler(M,e);
r = a.*(1-e.*cos(E));
psi = rand(stream,N,1)*2*pi;
theta =  acos(rand(stream,N,1)*2 - 1);
phi = rand(stream,N,1)*2*pi;

rps_proj = [r.*(cos(nu).*sin(theta).*sin(phi)+sin(nu).*(cos(theta).*cos(psi).*sin(phi)+cos(phi).*sin(psi))),...
            r.*(cos(nu).*cos(phi).*sin(theta)+sin(nu).*(cos(theta).*cos(phi).*cos(psi)-sin(phi).*sin(psi))),...
            r.*(cos(theta).*cos(nu)-cos(psi).*sin(theta).*sin(nu))];
        
s = [s;sqrt(rps_proj(:,1).^2 + rps_proj(:,2).^2)];
end

trans = []
while length(trans) < 1e6
    length(trans)
    R = rand(stream,N,1)*0.5;
    trans = [trans; R(s(floor(rand(N,1)*1e7+1))/Rs < 1 + R)];
end

 [n,x] = hist(trans,1000);
 n = n/length(trans)/mean(diff(x));


rs = linspace(0,max(r),100);
ns = hist(s,rs);    
ps = zeros(length(rs),1);
for j=1:length(rs),j,ps(j) = sint_euniform(rs(j),fa,amin,amax);end


