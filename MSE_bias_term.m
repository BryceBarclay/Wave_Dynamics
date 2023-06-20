function [BIAS] = MSE_bias_term(dt,T,omegaS,omega,phiS)
% computes phi*D_CC
%% variables
t = 0:dt:T;

integrand = zeros(length(t),length(omegaS));
for k = 1:length(omegaS)
    integrand(:,k) = phiS(k) * cos(omegaS(k)*t) .* cos(omega*t);
end

BIAS = zeros(length(t),length(omegaS));
for i = 2:length(t)
    for k = 1:length(omegaS)
        BIAS(i,k) = (1/t(i))*trapz(t(1:i),integrand(1:i,k)); 
    end
end

end