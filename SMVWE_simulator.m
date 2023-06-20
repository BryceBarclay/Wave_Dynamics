function [X,XS] = SMVWE_simulator(dt,T,omegaS,omegaN,phiS,phiN,sigma,f)

% variables
t = 0:dt:T;
X = zeros(length(t),1);
XS = zeros(length(t),1);
dB = zeros(length(t),length(omegaN));

%% simulate sensor time series
[SS,TT] = meshgrid(t,t);
integrand = zeros([length(t),length(t),length(omegaN)]);
for i = 1:length(omegaN)
    integrand(:,:,i) = phiN(i)*(sigma(i)/omegaN(i)) * sin(omegaN(i)*(TT-SS));
end

% IC
for j = 1:length(omegaS)
    X(1) = X(1) + phiS(j)*f(j)*cos(omegaS(j)*t(1));
    XS(1) = XS(1) + phiS(j)*f(j)*cos(omegaS(j)*t(1));
end

for i = 2:length(t)
    z = normrnd(0,1,[1,length(omegaN)]);
    dB(i,:) = sqrt(dt)*z;

    for j = 1:length(omegaN)
        X(i) = X(i) + sum(integrand(i,1:i-1,j)'.*dB(2:i,j));
    end
    for j = 1:length(omegaS)
        X(i) = X(i) + phiS(j)*f(j)*cos(omegaS(j)*t(i));
        XS(i) = XS(i) + phiS(j)*f(j)*cos(omegaS(j)*t(i));
    end
end

end
