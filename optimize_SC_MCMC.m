function optParams = optimize_SC_MCMC(m_known, lb, ub, xopt, xtilt, yopt, ytilt, zopt, u1d, ...
    invStdPWRL, tiltstd, tiltreduced, nanstatbeginning, ntrials, saveFigs)

GPS_std = 1./invStdPWRL;
vol = (4/3 * pi * m_known(1) * m_known(2)^2);
priormeans = [vol, m_known(8), m_known(9), m_known(10), m_known(11), m_known(end)];
priorvariances = (ub - lb)/3;
paramNames = ["HMM volume", "dpHMM", "vert semi-diameter", "horiz semi-diameter", "dip", "dpSC"];

bnds = [lb; ub]';
sigma = ones(size(u1d));
for i = 1:3
    sigma(:, i) = sigma(:, i) * GPS_std(i);
end

% for i = 1:30
% xstep = 5^(-i)
xstep = 0.007;
[x_keep, L_keep, count] = mcmc('create_MCMC_data',u1d,priormeans,xstep, bnds, sigma, ntrials, ...
    m_known, [xopt, xtilt], [yopt, ytilt], [zopt, 0], u1d, ...
    [invStdPWRL(1), invStdPWRL(2), invStdPWRL(3),1/((tiltstd)), 1/((tiltstd))], tiltreduced(1:2), nanstatbeginning);
% Perform the optimizationdisp("Success rate: " + num2str(100*count/ntrials) + "%, xstep: " + xstep + " L_keep: " + max(L_keep));
% end

% disp(size(x_keep))
burn = 1000;

% Plot histograms of paramter distribution
figure(1);
for i = 1:length(paramNames)
    subplot(320 + i); 
    hist(x_keep(i,burn:end));
    title(paramNames(i));
end
if(saveFigs); saveas(1, "Figures/MCMC_hist_" + num2str(ntrials, "%.1e") + "trials.fig"); end

figure(2);
clf;
hold on
sgtitle("Step size = " + num2str(xstep));
for i = 1:length(paramNames)
    hold on;
    subplot(320 + i); 
    plot(x_keep(i,burn:end));
    yline(lb(i), 'LineWidth', 3);
    yline(ub(i), 'LineWidth', 3);
    ylim([lb(i), ub(i)]);
    title(paramNames(i));
    hold off;
end
if(saveFigs); saveas(2, "Figures/param_evolution_" + num2str(ntrials, "%.1e") + "trials.fig"); end

figure(3);
clf;
plot(L_keep);
xlabel("Iterations")
ylabel("Log Likelihood")
if(saveFigs); saveas(3, "Figures/likelihood_" + num2str(ntrials, "%.1e") + "trials.fig"); end

[~, idx] = max(L_keep);  % Index of the highest log-likelihood
optParams = x_keep(:, idx);
