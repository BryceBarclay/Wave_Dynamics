function dbetadt = beta_ev_exp_wave(x3,beta,lambda,xi,cm20,gamma)
% ODE function for the eigenvalue problem in x_3.

% c^(-2)(x_3):
cm2 = cm20*exp(gamma*x3);
k = sqrt(cm2*lambda - xi);

dbetadt = [beta(2); -(k.^2).*beta(1)];

end