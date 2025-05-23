function [dp_low, dp_high] = GetErrors(N_draws, N_noise, posterior, ntime, ux, uy, uz, tiltx, tilty, ...
    dispstd, GPSNameList, rw_stddev, dp_weight, taiyi_parameters, npitloc, gps_sigma, nanstatbeginning, optParams)

%% First construct the covariance matrices to generate added nosie
%%% Delete dp, optimizedM when it works
N = length(tiltx);          % number of measurements
C_rw = zeros(N,N);
parfor i = 1:N
    for j = 1:N
        C_rw(i,j) = rw_stddev^2 * min(i,j);
    end
end
gps_sigma = 1./gps_sigma;

%% Use posterior distribution to generate a list of many candidate geometries
dtheta = 0;
randIdx = randi(size(posterior, 1), [1, N_draws]);
geo_samples = posterior(:, randIdx);
geo_samples(:, 1) = optParams;
dp_dist = zeros(N_draws, N_noise, ntime*2);
% h = waitbar(0,'Analyzing errors...');

if isempty(gcp('nocreate'))
    parpool('local');  % or a specific number of workers: parpool('local',4)
end

%% TODO - Incorporate uncertainties in shear modulus (and maybe poisson's ratio) into error analysis
for i = 1:N_draws
    % Format new geometry into optimizedM array
    m_samp = get_full_m(taiyi_parameters, geo_samples(:,i), true);

    % Create new green's functions
    [gHMM_samp, gSC_samp] = creategreens(m_samp(1:8), m_samp(9:end));
    [gTiltHMM_samp, gTiltSC_samp] = createtiltgreens(m_samp(1:8), m_samp(9:end), dtheta, false);
    gHMMflat_samp = gHMM_samp';
    gHMMflat_samp = gHMMflat_samp(:);
    gSCflat_samp = gSC_samp';
    gSCflat_samp = gSCflat_samp(:);

    temp = TimeDependentLSQtilt(gHMMflat_samp, gSCflat_samp, gTiltHMM_samp, gTiltSC_samp, ux, uy, uz, tiltx, tilty, ...
        dispstd, GPSNameList, rw_stddev, dp_weight, true);

    dp_ideal = temp(1:2*length(tiltx)); %[dp(:, 1), dp(:, 2)];
    offsets = temp((2*length(tiltx) + 1):end);
    offsets = reshape(offsets, [], 3);
    
    %% Now take the inverted pressure history and generate synthetic data
    % Use that synthetic data and add noise. Then re-invert to get the
    % variance of those individual estimates
    for j = 1:N_noise
        ux_pred = gHMM_samp(1,:)' * dp_ideal(1:length(tiltx)) + gSC_samp(1,:)' * dp_ideal(length(tiltx)+1:2*length(tiltx));
        ux_pred(nanstatbeginning, :) = ux_pred(nanstatbeginning, :) + offsets(:, 1);

        uy_pred = gHMM_samp(2,:)' * dp_ideal(1:length(tiltx)) + gSC_samp(2,:)' * dp_ideal(length(tiltx)+1:2*length(tiltx));
        uy_pred(nanstatbeginning, :) = uy_pred(nanstatbeginning, :) + offsets(:, 2);

        uz_pred = gHMM_samp(3,:)' * dp_ideal(1:length(tiltx)) + gSC_samp(3,:)' * dp_ideal(length(tiltx)+1:2*length(tiltx));
        uz_pred(nanstatbeginning, :) = uz_pred(nanstatbeginning, :) + offsets(:, 3);

        tilte_pred = gTiltHMM_samp(1) .* dp_ideal(1:length(tiltx)) + gTiltSC_samp(1) .* dp_ideal(length(tiltx)+1:2*length(tiltx));
        tiltn_pred = gTiltHMM_samp(2) .* dp_ideal(1:length(tiltx)) + gTiltSC_samp(2) .* dp_ideal(length(tiltx)+1:2*length(tiltx));
    
        % Generate noise for each data
        rw_east_noise = mvnrnd(zeros(length(C_rw), 1), C_rw);
        rw_north_noise = mvnrnd(zeros(length(C_rw), 1), C_rw);
        gps_e_noise = normrnd(0,gps_sigma(1),[1,length(ux)]);
        gps_n_noise = normrnd(0,gps_sigma(2),[1,length(ux)]);
        gps_u_noise = normrnd(0,gps_sigma(3),[1,length(ux)]);
    
        % Add generated noise to each data
        tilte_pred = tilte_pred + rw_east_noise;
        tiltn_pred = tiltn_pred + rw_north_noise;
        ux_pred = ux_pred + gps_e_noise;
        uy_pred = uy_pred + gps_n_noise;
        uz_pred = uz_pred + gps_u_noise;

        % Make NaN values for ux, uy, uz to match original data
        nanx = isnan(ux); nanx(:, 1) = 0;
        nany = isnan(ux); nany(:, 1) = 0;
        nanz = isnan(ux); nanz(:, 1) = 0;

        ux_pred(nanx) = "NaN"; 
        uy_pred(nany) = "NaN";
        uz_pred(nanz) = "NaN";
    
        %% Re-invert for pressure of now re-noised synthetic data
        temp = TimeDependentLSQtilt(gHMMflat_samp, gSCflat_samp, gTiltHMM_samp, gTiltSC_samp, ux_pred, uy_pred, uz_pred, ...
             tilte_pred, tiltn_pred, dispstd, GPSNameList, rw_stddev, dp_weight, true);
        dp_dist(i, j, :) = [temp(1:length(tiltx)) .* m_samp(8)/(1e6), temp((length(tiltx) + 1):2*length(tiltx)) .* m_samp(16)/(1e6)];
    end
end
dp_dist = reshape(dp_dist, N_draws * N_noise, ntime*2);

clear gHMM_samp gSC_samp gTiltHMM_samp gTiltSC_samp temp

%% Now pick the 10th and 90th percentile points for each sample

% pLow = zeros(size(dp_dist));
% pHigh = zeros(size(dp_dist));
% for i = 1:size(dp_dist, 2)
%     pLow(i) = prctile(dp_dist(:, i), 10);
%     pHigh(i) = prctile(dp_dist(:, i), 90);
% end

pLow  = prctile(dp_dist, 10, 1)';  % 10th percentile across columns
pHigh = prctile(dp_dist, 90, 1)';  % 90th percentile
dp_low = real([pLow(1:length(tiltx)), pLow((length(tiltx) + 1):(2*length(tiltx)))]);
dp_high = real([pHigh(1:length(tiltx)), pHigh((length(tiltx) + 1):(2*length(tiltx)))]);

%% Make a plot of all the diff pressure histories
figure(1);
clf;
for i = 2:(N_draws*N_noise)
    plot(dp_dist(i, length(tiltx) + 1:end), 'HandleVisibility','off');
    hold on;
end

plot(dp_dist(1, length(tiltx) + 1:end), "LineWidth", 4, "DisplayName", "MLE");
plot(dp_low(:, 2), "LineWidth", 6, "DisplayName", "10th percentile");
plot(dp_high(:, 2), "LineWidth", 6, "DisplayName", "90th percentile");
legend();
end