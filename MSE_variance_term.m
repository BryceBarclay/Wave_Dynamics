function [VAR] = MSE_variance_term(dt,T,MC,omegaN,omega,phiN,sigma)
% computes sum E[(sigma*phi*S_C)^2]
%% variables
t = 0:dt:T;
stoch = zeros(length(t),MC,length(omegaN)); % sigma*phi*S_C

%% simulate solution

% integrand
[SS,TT] = meshgrid(t,t);
integrand = zeros([length(t),length(t),length(omegaN)]);
for k = 1:length(omegaN)
    integrand(:,:,k) = phiN(k)*(sigma(k)/omegaN(k)) * sin(omegaN(k)*(TT-SS)) .* cos(omega*TT); 
end

dB = sqrt(dt)*normrnd(0,1,[length(t),MC,length(omegaN)]);
dB(1,:,:) = zeros(size(dB(1,:,:)));

int = zeros(length(t),length(omegaN));
for mc = 1:MC
    for i = 2:length(t)
        for j = 2:i
            for k = 1:length(omegaN)
                int(j,k) = sum(integrand(j,1:j-1,k)'.*dB(2:j,mc,k));
            end
        end
        for k = 1:length(omegaN)
            stoch(i,mc,k) = (t(i-1)/t(i))*stoch(i-1,mc,k) + (1/t(i))*trapz(t(i-1:i),int(i-1:i,k));
        end
    end
end

%% analysis
VAR = zeros(length(t),1);
for k = 1:length(omegaN)
    VAR = VAR + mean(stoch(:,:,k).^2,2);
end

end