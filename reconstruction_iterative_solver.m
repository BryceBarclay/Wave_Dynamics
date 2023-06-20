%% Iterative algorithm to solve for Fourier Coefficients

% set seed of random number generator
rng('default');
rng(1);

% numerical parameters
dt = 0.05;
T = 40; % T = 5 takes 100s and T = 20 takes 1600s (MC=200, omegaN = 1:20)
t = 0:dt:T;
MC = 500;

% 3D simulation parameters
numeig = 8; % in each dimension. Total is numeig^3
L = 1; W = pi; H = exp(1);
slidx = [10,20,30]; % index of sensor location
x = linspace(0,L,100); y = linspace(0,W,100); z = linspace(0,H,100);
[PHI,lam] = eigen_analytical(numeig,x,y,z,L,W,H);
omegaSidx = [2,10,15];
omegaNidx = 1:20;
omegaS = sqrt(lam(omegaSidx));
omegaN = sqrt(lam(omegaNidx));
f = [-0.5,1,1.5];
sigma = 0.2*ones(length(omegaN),1);
phiS = PHI(slidx(1),slidx(2),slidx(3),omegaSidx);
phiN = PHI(slidx(1),slidx(2),slidx(3),omegaNidx);


%% simulate noisy signal

[X,XS] = SMVWE_simulator(dt,T,omegaS,omegaN,phiS,phiN,sigma,f);

%% frequency analysis
% omega = freq*2pi
[freq,Xhat1] = one_sided_spectrum(t,X);
[freq,XShat1] = one_sided_spectrum(t,XS);

% find signal omegas (highest peaks)
[ampsort,indsort] = sort(Xhat1);
peakindex = [];
while size(peakindex) < 3
    peakindex = [peakindex,indsort(end)];
    indsort = indsort(abs(indsort-peakindex(end))>1);
end
omegaS_approx1 = 2*pi*freq(peakindex);
[M,I] = min(abs(omegaS_approx1 - sqrt(lam)),[],1); 
omegaS_approx2 = sqrt(lam(I));

figure;
plot(freq,Xhat1,'LineWidth',2.0);
title("Single-Sided Amplitude Spectrum of X(t)",'FontSize',18);
xlabel("f (Hz)",'FontSize',18);
ylabel("|Xhat1(f)|",'FontSize',18);

figure;
plot(freq,XShat1,'LineWidth',2.0);
title("Single-Sided Amplitude Spectrum of X(t)",'FontSize',18);
xlabel("f (Hz)",'FontSize',18);
ylabel("|XShat1(f)|",'FontSize',18);

%% simulation for error terms

BIAS = zeros(length(t),length(omegaS),length(omegaS));
VAR = zeros(length(t),length(omegaS));

tic;
for i = 1:length(omegaS)
    VAR(:,i) = MSE_variance_term(dt,T,MC,omegaN,omegaS(i),phiN,sigma);
    BIAS(:,:,i) = MSE_bias_term(dt,T,omegaS,omegaS(i),phiS);
end
toc;


%% iterative solver

iter = 10;

% initialize coefficients
fhat = zeros(iter,length(omegaS));
fhat(1,:) = ones(1,length(omegaS));
taumin = zeros(iter,length(omegaS));

for i = 2:iter
    for k = 1:length(omegaS)
        idx = [1:k-1,k+1:length(omegaS)]; % ~k
        error = BIAS(:,k,k).^(-2) .* ((BIAS(:,idx,k)*fhat(i-1,idx)').^2 + VAR(:,k));
        [errormin,ind] = min(error(2:end));
        taumin(i,k) = t(ind+1);
        fhat(i,k) = mean(X(1:ind+1)'.*cos(omegaS(k)*t(1:ind+1))) / (BIAS(ind+1,k,k));
    end
end

%% analytical error terms (for comparison)
% compare simple simulations to analytical error terms
% 1 stochastic term and 2 deterministic terms 
omega1 = omegaS(1);
omega2 = omegaS(2);
f1 = f(1);
f2 = f(2);
sigma1 = sigma(1);
phi1 = phiS(1);
phi2 = phiS(2);

p = t*omega1;

normalize = zeros(1,length(t)) + 1;
for i = 2:length(t)
    normalize(i) = trapz(t(1:i),cos(omega1*t(1:i)).^2)/t(i);
end

% Deterministic error
Dmn = phi2*omega1./(2*(omega1^2-omega2^2)*p).*((omega1-omega2)*sin(p*(1+omega2/omega1)) + (omega1+omega2)*sin(p*(1-omega2/omega1)));
det_error = 1./(abs(phi1*f1)*normalize) .* sqrt(f2^2*Dmn.^2);

% Stochastic error
ESnn2_cosine = phi1^2/(192*omega1^3)*(8*p + 12*sin(2*p) - 12./p + (3./(p.^2)).*sin(4*p));
ESnn2_cosine(1) = 0;% remove NaN
stoch_error = 1./(abs(phi1*f1)*normalize) .* sqrt(sigma1^2*ESnn2_cosine);

% MSE vs tau
Dcc12 = 0.5 * (sin((omega1 + omega2)*t)./((omega1 + omega2)*t) + sin((omega1 - omega2)*t)./((omega1 - omega2)*t));
Dcc11 = 0.5 * (sin((omega1 + omega1)*t)./((omega1 + omega1)*t) + 1);

BIAS_analytic = f2 * (phi2/phi1) * Dcc12./Dcc11;
BIAS2_analytic = BIAS_analytic.^2;

ESc2 = 1/(192*omega1^3) * (8*omega1*t + 12*sin(2*omega1*t) - 12./(omega1*t) + 3./(omega1^2*t.^2).*sin(4*omega1*t));

VAR_analytic = (sigma1./Dcc11).^2 .* ESc2;

MSE_analytic = BIAS2_analytic + VAR_analytic;


%% Analysis of simulation
% theoretical error
error_analytic = zeros(length(t),length(omegaS));
for k = 1:length(omegaS)
    idx = [1:k-1,k+1:length(omegaS)];
    error_analytic(:,k) = BIAS(:,k,k).^(-2) .* ((BIAS(:,idx,k)*f(idx)').^2 + VAR(:,k));
end

% Fourier coefficients vs observation time
fhatsim = zeros(length(t)-1,length(omegaS));
for i = 2:length(t)
    for k = 1:length(omegaS)
        fhatsim(i-1,k) = mean(X(1:i)'.*cos(omegaS(k)*t(1:i))) / BIAS(i,k,k);
    end
end

% solutions
exact_soln = zeros(1,length(t));
approx_soln = zeros(1,length(t));
tf_soln = zeros(1,length(t));
for k = 1:length(omegaS)
    exact_soln = exact_soln + phiS(k)*f(k)*cos(omegaS(k)*t);
    approx_soln = approx_soln + phiS(k)*fhat(end,k)*cos(omegaS(k)*t);
    tf_soln = tf_soln + phiS(k)*fhatsim(end,k)*cos(omegaS(k)*t); % approximate solution using final time
end

% L2 error
approx_error = trapz(t,(exact_soln - approx_soln).^2);
tf_error = trapz(t,(exact_soln - tf_soln).^2);
soln_norm = trapz(t,(exact_soln).^2);
approx_rel_error = sqrt(approx_error/soln_norm)
tf_rel_error = sqrt(tf_error/soln_norm)
tf_rel_error/approx_rel_error

norm(exact_soln - tf_soln,2)/norm(exact_soln - approx_soln,2)

%% plot

% comparison with theoretical expressions
figure;
plot(t(2:end),error_analytic(2:end,1),'LineWidth',1.5);
hold on;
plot(t(2:end),MSE_analytic(2:end),'LineWidth',1.5);
xlabel('t','FontSize',18);
ylabel('error','FontSize',18);

figure;
for k = 1:length(omegaS)
    plot(t(2:end),error_analytic(2:end,k),'LineWidth',1.5);
    hold on;
end
xlabel('t','FontSize',18);
ylabel('error','FontSize',18);

figure;
for k = 1:length(omegaS)
    plot(2:iter,taumin(2:end,k),'LineWidth',1.5);
    hold on;
end
xlabel('iteration','FontSize',18);
ylabel('optimal \tau','FontSize',18);

figure;
plot(t,X'-exact_soln,'LineWidth',2.0)
xlabel('t','FontSize',18);
ylabel('Noisy - noiseless signal','FontSize',18);

% coefficient estimate vs observation time
figure;
hold on;
for i = 1:length(omegaS)
    plot(t(2:end),abs(fhatsim(:,i) - f(i)),'LineWidth',2.0)
    plot(taumin(end,i),abs(fhat(end,i) - f(i)),'o')
end
xlabel('\tau','FontSize',18);
ylabel('fhat error','FontSize',18);

% all plotted
figure;
hold on;
plot(t,exact_soln,'-','LineWidth',2.0)
plot(t,approx_soln(end,:),'--','LineWidth',2.0)
plot(t,tf_soln,'--','LineWidth',2.0)
p = plot(t,X,'-','LineWidth',2.0);
p.Color = [p.Color 0.5];
xlabel('t','FontSize',18);
ylabel('E_1','FontSize',18);
legend('exact noiseless','reconstructed','final time reconstructed','exact noisy');

% partial plots
figure;
hold on;
plot(t,exact_soln,'-','LineWidth',2.0)
plot(t,approx_soln(end,:),'-','LineWidth',2.0)
p = plot(t,tf_soln,'-','LineWidth',2.0);
p.Color = [p.Color 0.5];
xlabel('t','FontSize',18);
ylabel('E_1','FontSize',18);
legend('\langle E_1\rangle','Alg 5.1 approximation','final time approximation','FontSize',16);
ylim([-0.6,0.4]);

figure;
hold on;
plot(t,exact_soln,'-','LineWidth',2.0)
p = plot(t,X,'-','LineWidth',2.0);
p.Color = [p.Color 0.5];
xlabel('t','FontSize',18);
ylabel('E_1','FontSize',18);
legend('\langle E_1\rangle','E_1','FontSize',16);
ylim([-0.6,0.4]);


