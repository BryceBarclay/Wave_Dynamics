%% Sparse Reconstruction of Solutions to Stratified Wave equation
% Medium is isotropic and stratified in z
%% parameters

% medium parameters
cm20 = 0.005;              % c^(-2)(0) 
gamma = 0.1;               % c^(-2)(x3) = cm20*exp(gamma*x3)
c0 = 1/sqrt(cm20);         % c(0)
dcm2dx0 = cm20*gamma;      % d(c^(-2))/dx at x3=0

% wave parameters
A = [1+0*1i,-0.0-0.2*1i];  % coefficients
omega = [40,47];           % frequencies
khat = [3;4]; khat = khat/norm(khat);
kabs = omega/c0;
k = khat*kabs;             % wavenumber at x3=0

% eigenvalues
lambda = omega.^2; 
xi = k(1,:).^2; 
xi3 = k(2,:).^2;           % cm2*lambda - xi

% set seed of random number generator
rng('default');
rng(1);

%% variables
x1dom = [0 10];
x3dom = [0 10];
tdom  = [0 30];

x1step = 0.1;
x3step = 0.1;
tstep = 0.02;
 
x1 = x1dom(1):x1step:x1dom(2);
x3 = x3dom(1):x3step:x3dom(2);
t = tdom(1):tstep:tdom(2);

[X1,X3] = meshgrid(x1,x3);

%% Eigenfunctions in z
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
beta1 = zeros(length(lambda),length(x3));
beta2 = zeros(length(lambda),length(x3));
BETA1 = zeros(length(x3),length(x1),length(lambda));
BETA2 = zeros(length(x3),length(x1),length(lambda));
for li = 1:length(lambda)
    k30 = sqrt(cm20*lambda(li) - xi(li));
    betabc = [1;-lambda(li)*dcm2dx0/(4*k30^2) + 1i*k30];

    [x31,beta1_vec] = ode45(@(x3,beta) beta_ev_exp_wave(x3,beta,lambda(li),xi(li),cm20,gamma),x3dom,real(betabc),opts);
    [x32,beta2_vec] = ode45(@(x3,beta) beta_ev_exp_wave(x3,beta,lambda(li),xi(li),cm20,gamma),x3dom,imag(betabc),opts);
    
    beta1(li,:) = interp1(x31,beta1_vec(:,1),x3);
    beta2(li,:) = interp1(x32,beta2_vec(:,1),x3);
    
    % grid needed for computation of E3
    [X1,BETA1(:,:,li)] = meshgrid(x1,beta1(li,:));
    [X1,BETA2(:,:,li)] = meshgrid(x1,beta2(li,:));
end

%% Sensor trajectory 

% circular trajectory
r = 3;
htraj = 4; ktraj = 4;
l = 1; % iterations around same loop
traj = [r*cos(l*2*pi*t/t(end)).' + htraj,r*sin(l*2*pi*t/(t(end))).' + ktraj];

% straight line acceleration
% dir = [1,2]; dir = dir/norm(dir,2);
% accel_max = min([x1(end),x3(end)]./dir)/t(end)^2;
% accel = accel_max;
% traj = accel*t'.^2*dir;

%% solution

tic;
% 2D eigenfunctions
phi = zeros(length(x1),length(x3),length(lambda));
for li = 1:length(lambda)
    phi(:,:,li) = exp(1i*(sqrt(xi(li))*X1))  .*  (reshape(BETA1(:,:,li),size(X1)) + 1i*reshape(BETA2(:,:,li),size(X1)));
end

Rs = zeros(1,length(t));
phi_s = zeros(length(t),length(lambda));
Rs0 = zeros(1,length(t)); % stationary sensor
%figure;
for ti = 1:length(t)
    % wave field
    E3 = zeros(size(X1));
    for li = 1:length(lambda)
        E3 = E3 + A(li) * exp(1i*(sqrt(xi(li))*X1 - sqrt(lambda(li))*t(ti)))  .*  (reshape(BETA1(:,:,li),size(X1)) + 1i*reshape(BETA2(:,:,li),size(X1)));
    end

    % sensor data
    Rs(ti) = interp2(X1,X3,(E3),traj(ti,1),traj(ti,2));
    Rs0(ti) = interp2(X1,X3,(E3),traj(1,1),traj(1,2));

    % eigenfunction data along sensor path
    for li = 1:length(lambda)
        phi_s(ti,li) = interp2(X1,X3,reshape(phi(:,:,li),size(X1)),traj(ti,1),traj(ti,2));
    end

    % plot
%     surf(x1,x3,real(E3));
%     title('E_3','FontSize',20);
%     xlabel('x_1','FontSize',20);
%     ylabel('x_3','FontSize',20);
%     zlim([-10,10]);
%     drawnow;
%     shg;
end
toc;

%% plots
% sensor trajectory
figure;
plot(traj(:,1),traj(:,2),'LineWidth',2.0);
%title('Sensor Trajectory','FontSize',18);
xlabel('x_1','FontSize',18);
ylabel('x_3','FontSize',18);
xlim([x1(1),x1(end)]);
ylim([x3(1),x3(end)]);

% eigenfunctions
% figure;
% hold on;
% plot(x3,beta1(1,:),'LineWidth',2.0);
% %plot(x3,real(exp(1i*(sqrt(cm20*exp(gamma*x3)*lambda(1) - xi(1))).*x3)),'LineWidth',2.0);
% xlabel('x_3','FontSize',18);
% ylabel('\beta_1','FontSize',18);
% 
% figure;
% hold on;
% plot(x3,beta2(1,:),'LineWidth',2.0);
% %plot(x3,imag(exp(1i*(sqrt(lambda - xi))*x3)),'LineWidth',2.0);
% xlabel('x_3','FontSize',18);
% ylabel('\beta_2','FontSize',18);

% sensor time series
% figure;
% plot(t,real(Rs),'LineWidth',2.0);
% title('Sensor Time Series','FontSize',20);
% xlabel('t','FontSize',20);
% ylabel('E_3(x(t),t)','FontSize',20);

%% ------------------------------------------------------------------------
%% -------------------------Sparse Reconstruction--------------------------
%% ------------------------------------------------------------------------
%% Fourier analysis
% sqrt(lambda) = frequency * 2pi
[f,Rshat] = one_sided_spectrum(t,real(Rs));
[f,Rs0hat] = one_sided_spectrum(t,real(Rs0));

% Max Doppler shift: omega -> (1 +/- max(v/c))omega
c = (cm20*exp(gamma*traj(1:end-1,2))).^(-1/2);
v = sqrt((traj(2:end,1) - traj(1:end-1,1)).^2 + ((traj(2:end,2) - traj(1:end-1,2)).^2))/tstep;
velratio_max = max(v)./min(c);

% exact signal frequencies (unknown for reconstruction)
exact_signal_freq = omega/(2*pi);
exact_freq_range = [exact_signal_freq*(1 - velratio_max);exact_signal_freq*(1 + velratio_max)];

% estimated signal frequencies
amplitude_tol = 5e-2;
spike_freq_bool = (Rshat > amplitude_tol);
spike_freq = f(spike_freq_bool); % spikes of FT

est_freq_range = [spike_freq*(1 + velratio_max).^(-1); spike_freq*(1 - velratio_max).^(-1)]';
freq_domain = unionofintervals(est_freq_range); % disjoint frequency intervals

%% Reconstruction

% random frequencies
N = 2000; % number of random guesses
omega_rand = 2*pi*unifrndintervals(freq_domain,N*length(lambda));
omega_rand = reshape(omega_rand,[length(lambda),N]);

% random eigenvalues (one direction of propagation)
lambda_rand = omega_rand.^2;
xi_rand = (omega_rand*khat(1,1)/c0).^2;

% length of time series used for least squares
ts_length = 500;

beta1_rand = zeros(length(lambda),length(x3));
beta2_rand = zeros(length(lambda),length(x3));
BETA1_rand = zeros(length(x3),length(x1),length(lambda));
BETA2_rand = zeros(length(x3),length(x1),length(lambda));

comp = zeros(length(x1),length(x3));
Rs_rand_comp = zeros(length(t),length(lambda));
A_rand = zeros(length(lambda),N);
error = zeros(N,1);

for i = 1:N
    %----create eigenfunctions----
    for li = 1:length(lambda)
        k30 = sqrt(cm20*lambda_rand(li,i) - xi_rand(li,i));
        betabc = [1;-lambda_rand(li,i)*dcm2dx0/(4*k30^2) + 1i*k30];
        [x31,beta1_vec] = ode45(@(x3,beta) beta_ev_exp_wave(x3,beta,lambda_rand(li,i),xi_rand(li,i),cm20,gamma),x3dom,real(betabc),opts);
        [x32,beta2_vec] = ode45(@(x3,beta) beta_ev_exp_wave(x3,beta,lambda_rand(li,i),xi_rand(li,i),cm20,gamma),x3dom,imag(betabc),opts);
    
        beta1_rand(li,:) = interp1(x31,beta1_vec(:,1),x3);
        beta2_rand(li,:) = interp1(x32,beta2_vec(:,1),x3);
    
        % grid needed for computation of E3
        [X1,BETA1_rand(:,:,li)] = meshgrid(x1,beta1_rand(li,:));
        [X1,BETA2_rand(:,:,li)] = meshgrid(x1,beta2_rand(li,:));
    end

    %----evaluate eigenfunctions on sensor trajectory----
    for ti = 1:length(t)
        for li = 1:length(lambda)
            comp = exp(1i*(sqrt(xi_rand(li,i))*X1 - sqrt(lambda_rand(li,i))*t(ti)))  .*  (reshape(BETA1_rand(:,:,li),size(X1)) + 1i*reshape(BETA2_rand(:,:,li),size(X1)));
            % sensor data
            Rs_rand_comp(ti,li) = interp2(X1,X3,comp,traj(ti,1),traj(ti,2));
        end
    end

    %----find coefficients minimizing L^2 error----
    A_rand(:,i) = Rs_rand_comp(1:ts_length,:)\Rs(1:ts_length).';
    Rs_rand = Rs_rand_comp*A_rand(:,i);
    error(i) = norm(Rs_rand(1:ts_length) - Rs(1:ts_length).',2)/norm(Rs(1:ts_length),2);

    if(error(i) == min(error(error > 0)))
        Rs_rand_best = Rs_rand; % retain lowest MSE time series
    end

    % display progress
    disp(i);
    disp(error(i));

end
toc;

%% Analyze results
% closest frequencies
[min_freq_diff,min_freq_diff_index] = min([vecnorm(omega_rand-[omega(2);omega(1)],2,1),vecnorm(omega_rand-omega',2,1)]);
if(min_freq_diff_index <= N)
    disp(omega_rand(:,min_freq_diff_index));
else
    disp(omega_rand(:,min_freq_diff_index-N));
end

% minimum error
[min_error,min_index] = min(error);
disp(omega_rand(:,min_index));
disp(min_error);

%% Plots

% original signal spectrum
figure;
hold on;
p = plot(2*pi*f,Rs0hat,'LineWidth',2.0);
plot(2*pi*f,Rshat,'LineWidth',2.0);
p.Color = [p.Color 0.5];
%title('One-sided spectrum','FontSize',18);
xlabel('\omega','FontSize',18);
ylabel('amplitude spectrum','FontSize',18);
legend('R_0','R_s','FontSize',18);

% comparison of result to original time series
figure;
plot(t(1:ts_length),real(Rs(1:ts_length)),t(1:ts_length),real(Rs_rand_best(1:ts_length)),'LineWidth',2.0);
%title('Time Series Comparison','FontSize',18);
xlabel('t','FontSize',18);

