function optParams = optimize_SC_MCMC(m_known, lb, ub, xopt, xtilt, yopt, ytilt, zopt, u1d, ...
    invStdPWRL, tiltstd, tiltreduced, nanstatbeginning)

GPS_std = 1./invStdPWRL;
priormeans = (ub + lb)/2;
priorvariances = (ub - lb)/3;
paramNames = ["HMM volume", "dpHMM", "vert semi-diameter", "horiz semi-diameter", "dip", "dpSC"];

bnds = [lb(1), ub(1); lb(2), ub(2); lb(3), ub(3); lb(4), ub(4); lb(5), ub(5)];
sigma = ones(size(u1d));
for i = 1:3
    sigma(:, i) = sigma(:, i) * GPS_std(i);
end

xstep = 0.0013717;
ntrials = 10000;
[x_keep, L_keep, count] = mcmc('create_MCMC_data',u1d,priormeans,xstep, bnds, sigma, ntrials, ...
    m_known, [xopt, xtilt], [yopt, ytilt], [zopt, 0], u1d, ...
    [invStdPWRL(1), invStdPWRL(2), invStdPWRL(3),1/((tiltstd)), 1/((tiltstd))], tiltreduced(1:2), nanstatbeginning, priormeans, priorvariances);
% Perform the optimization
disp("Success rate: " + num2str(count) + ", xstep: " + xstep);


% disp(size(x_keep))
burn = 1;

% Plot histograms of paramter distribution
figure(1);
for i = 1:length(paramNames)
    subplot(320 + i); 
    hist(x_keep(i,burn:end));
    title(paramNames(i));
end

figure(2)
for i = 1:length(paramNames)
    subplot(320 + i); 
    plot(x_keep(i,burn:end));
    title(paramNames(i));
end

optParams = x_keep(:, end);
