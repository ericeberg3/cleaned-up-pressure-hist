function [optParams, posterior, gps_l2, insar_l2] = optimize_SC_MCMC(m_known, lb, ub, xopt, yopt, zopt, u1d, ...
     insarx, insary, insaru, look, insarweight, invStdPWRL, tiltstd, tiltreduced, nanstatbeginning, ntrials, gps_weight, saveFigs)

GPS_std = 1./invStdPWRL;
vol = (4/3 * pi * m_known(1) * m_known(2)^2);
priormeans = [vol, m_known(8), m_known(9), m_known(10), m_known(11), m_known(end)];
priormeans = priormeans + 0.1*randn(1,6) .* priormeans;
paramNames = ["HMM volume", "dpHMM", "vert semi-diameter", "horiz semi-diameter", "dip", "dpSC"];

bnds = [lb; ub]';
sigma = ones(size([u1d(:);insaru])) .* insarweight; % Insar stddev ~0.2m
% sigma = ones(size(insaru));
for i = 1:length(u1d(:))
    sigma(i) = GPS_std(ceil(i / size(u1d, 1)));
end

% [u1d(:);insaru]
% xstep = 0.02*ones(1,6); % 0.02
% xstep(4) = 0.007; % horiz semi-diam 

xstep = 7e-2*ones(1,6); % 0.02

% gps_weight = 4;
% xstep(4) = 0.02; % horiz semi-diam 
[x_keep, L_keep, count, gps_l2, insar_l2] = mcmc('create_MCMC_data',[u1d(:);insaru],priormeans,xstep, bnds, sigma, ntrials, gps_weight, ...
    m_known, xopt, yopt, zopt, u1d, insarx, insary, insaru, look, insarweight, ...
    [invStdPWRL(1), invStdPWRL(2), invStdPWRL(3),1/((tiltstd)), 1/((tiltstd))], tiltreduced(1:2), nanstatbeginning);

% [x_keep, L_keep, count] = mcmc_adaptive('create_MCMC_data',[u1d(:);insaru],priormeans,xstep, bnds, sigma, ntrials, ...
%     adapt_start, adapt_interval, scaling, epsilon, ...
%     m_known, xopt, yopt, zopt, u1d, insarx, insary, insaru, look, insarweight, ...
%     [invStdPWRL(1), invStdPWRL(2), invStdPWRL(3),1/((tiltstd)), 1/((tiltstd))], tiltreduced(1:2), nanstatbeginning);

gps_l2 = 0;
insar_l2 = 0;

burn = 5000;

% Store results:
[~, idx] = max(L_keep);  % Index of the highest log-likelihood
optParams = x_keep(:, idx);
posterior = x_keep(:, burn:end);

%% Plotting

% Plot histograms of paramter distribution
figure(5);
clf;
for i = 1:length(paramNames)
    subplot(320 + i); 
    numBins = 100;
    [counts, edges] = histcounts(x_keep(i, burn:end), numBins, 'Normalization', 'pdf');
    binCenters = edges(1:end-1) + diff(edges)/2;

    plot(binCenters, counts, 'LineWidth', 1.5);
    hold on;

    xline(optParams(i), '--r', 'LineWidth', 1.5);
    xline(priormeans(i), '-b', 'LineWidth', 1.5);

    title(paramNames(i));
    hold off;
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


%% Get top 10% of estimates and plot the volume of HMM vs. volume of SC
% Find the number of elements corresponding to the top 10%
num_top_elements = ceil(0.1 * length(L_keep));

% Sort the list in descending order and get the indices
[~, sorted_indices] = sort(L_keep, 'descend');

% Get the indices of the top 10% values
top_indices = sorted_indices(1:num_top_elements);

% figure(4);
% clf;
% plot(x_keep(1, top_indices) .* x_keep(2, top_indices), (4/3) * pi * x_keep(3, top_indices) .* (x_keep(4, top_indices)).^2 .* x_keep(6, top_indices), '.')
% xlabel("dpHMM * vHMM");
% ylabel("dpSC * vSC");
% if(saveFigs); saveas(4, "Figures/volume_correl_" + num2str(ntrials, "%.1e") + "trials.fig"); end

end