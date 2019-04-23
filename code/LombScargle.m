function [f,P,prob] = LombScargle(t,h,ofac,hifac)
%LombScargle - Compute normalized Lomb-Scargle Periodogram
%
% INPUTS
% t     (N x 1) Array of times for N observations
% h     (N x 1) Array of observations
% ofac  Frequency oversampling rate
% hifac Highest frequency to consider (multiple of Nyquist frequency)
%
% OUTPUTS
% f     Array of M frequencies sampled
% P     Array of M periodogram values
% prob  Array of M probabilities of null hypothesis false alarm

%condition input
t = t(:);
h = h(:);

%sample length and time span
N = length(h);
T = max(t) - min(t);

%mean and variance 
mu = mean(h);
s2 = var(h);

%calculate sampling frequencies
f = (1/(T*ofac):1/(T*ofac):hifac*N/(2*T)).';

%angular frequencies and constant offsets
w = 2*pi*f;
tau = atan2(sum(sin(2*w*t.'),2),sum(cos(2*w*t.'),2))./(2*w);

%spectral power
cterm = cos(w*t.' - repmat(w.*tau,1,length(t)));
sterm = sin(w*t.' - repmat(w.*tau,1,length(t)));
P = (sum(cterm*diag(h-mu),2).^2./sum(cterm.^2,2) + ...
     sum(sterm*diag(h-mu),2).^2./sum(sterm.^2,2))/(2*s2);

%estimate of the number of independent frequencies
M=2*length(f)/ofac;

%statistical significane of power
prob = M*exp(-P);
inds = prob > 0.01;
prob(inds) = 1-(1-exp(-P(inds))).^M;
