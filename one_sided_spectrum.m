function [f,Xhat1] = one_sided_spectrum(t,X)
% This function provides the one-sided spectrum Xhat1 of a real signal X.

L = length(X);
dt = t(2) - t(1);

fs = 1/dt;
f = fs*(0:(L/2))/L;

Xhat = fft(X);
Xhat2 = abs(Xhat)/L; % two-sided spectrum
Xhat1 = Xhat2(1:L/2+1); Xhat1(2:end-1) = 2*Xhat1(2:end-1); % one-sided spectrum

end