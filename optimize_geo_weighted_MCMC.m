%% This script uses green's functions created from an optimized HMM and SC to invert for the pressure history

% Loading in data and setting the downsample rate (from 5s to hourly)
downsamplerate = 720;
load(['data/u_mm_full_duration_' int2str(downsamplerate) '_interp.mat']);
load('data/gps_time_series.mat');
load('data/GPS_seism_locations.mat', 'GPSNameList', 'GPS_llh');
load('data/hawaii_line_new.mat', 'coast_new', 'new_pit')
load('data/sdh_clean.mat')
daily_GPS = readtable('data/daily_gps_mult_stations.csv');

%% Creating GPS names and setting up x and y coordinates
GPSNameList = ["69FL","92YN","AHUP","BDPK","BYRL","CNPK","CRIM","DEVL","OUTL","PUHI","PWRL","UWEV","V120","VSAS", "CALS"];
tiltNameList = ['SDH'];

t = downsample(universalt, downsamplerate);
ntime = length(t);

% convert u from mm to m
u = u .* 1e-3;

% Set outliers to NaN
for i = 1:size(u, 1)
    u(i, :, :) = filloutliers(squeeze(u(i, :, :)), NaN, "movmedian", 300, 2);
end 

% Center the displacements at the beginning of each time series at 0. Done 
% by taking the median of the first 20 non-nan data points for each station
u = CenterDisps(1, 20, u, GPSNameList, ntime);

clearvars interpu unfilledu nanu meansize rtnet_ka windowsize startcind new_pit endcind disps;

nstat = length(GPSNameList); % Number of GPS stations

% station x and y positions
xy = llh2local(GPS_llh, [-155.2784, 19.4073]);
x = xy(1,:) * 1000;
y = xy(2,:) * 1000;
z = zeros(1, nstat - 1);
%% Cleaning up stations with too many nan values

% Removing CALS
u(end, :, :) = [];
GPSNameList(end) = [];
nstat = nstat - 1;

%% Set up daily data to create u1d from June 12th - July 30th
nanstatend = false(1, length(GPSNameList));
startind = 43;
finalind = 91;

PWRLdata = [daily_GPS{startind:finalind, "PWRL_east"}, daily_GPS{startind:finalind, "PWRL_north"}, daily_GPS{startind:finalind, "PWRL_up"}];
daily_inv_std = 1./std(PWRLdata, "omitnan");

u1d = zeros(length(GPSNameList), 3);
for i = 1:length(GPSNameList)
    station = char(GPSNameList(i));
    check = all(ismember(station(1), '0123456789'));

    if(GPSNameList(i) == "CRIM" || GPSNameList(i) == "UWEV" || GPSNameList(i) == "BYRL") nanstatend(i) = true; end
    try
        if(check == 0)
            u1d(i, 1) = daily_GPS{finalind, GPSNameList(i) + "_east"} - daily_GPS{startind, GPSNameList(i) + "_east"};
            u1d(i, 2) = daily_GPS{finalind, GPSNameList(i) + "_north"} - daily_GPS{startind, GPSNameList(i) + "_north"};
            u1d(i, 3) = daily_GPS{finalind, GPSNameList(i) + "_up"} - daily_GPS{startind, GPSNameList(i) + "_up"};
        else
            u1d(i, 1) = daily_GPS{finalind, "x" + GPSNameList(i) + "_east"} - daily_GPS{startind, "x" + GPSNameList(i) + "_east"};
            u1d(i, 2) = daily_GPS{finalind, "x" + GPSNameList(i) + "_north"} - daily_GPS{startind, "x" + GPSNameList(i) + "_north"};
            u1d(i, 3) = daily_GPS{finalind, "x" + GPSNameList(i) + "_up"} - daily_GPS{startind, "x" + GPSNameList(i) + "_up"};
        end
    catch
        nanstatend(i) = true;
    end
    if(isnan(u1d(i, 1))); nanstatend(i) = true; end
end
u1d = u1d(~nanstatend, :);

%% COMMENT LATER
% u1d = u(:, :, end-finalind);
% u1d = u1d(~nanstatend, :);

clear PWRLdata daily_GPS

% Plotting this stuff
% Get original u1d to compare the two in a vector plot;

% finalindex = 114;
% u1d_5s = u(:, :, end-finalindex);
% u1d_5s = u1d_5s(~nanstatend, :);
% 
% hold off;
% figure(1);
% radscale = 1e4;
% clf;
% daily_quiver = quiver3(x(~nanstatend)', y(~nanstatend)', zeros(size(y(~nanstatend)))', u1d(:, 1) * radscale, u1d(:, 2) * radscale, u1d(:, 3) * radscale, 'AutoScale', 'off', 'LineWidth',2.75, 'MaxHeadSize', 0.3, 'Color', '#4DBEEE', 'DisplayName', 'Daily');
% hold on;
% secondly_quiver = quiver3(x(~nanstatend)', y(~nanstatend)', zeros(size(y(~nanstatend)))', u1d_5s(:, 1) * radscale, u1d_5s(:, 2) * radscale, u1d_5s(:, 3) * radscale, 'AutoScale', 'off', ...
%     'LineWidth',2.75, 'MaxHeadSize', 0.3, 'Color', '#A2142F', 'DisplayName', 'High frequency');
% cxy = llh2local(coast_new', [-155.2784, 19.4073]);
% cxy = cxy * 1000;
% plot(cxy(1, :)', cxy(2, :)', 'k.', 'HandleVisibility','off');
% xlabel('x (m)');
% ylabel('y (m)');
% zlabel('Scaled Displacement (m)');
% title("Displacements at various stations (y oriented N, x oriented E)")
% xlim([-8000, 8000]);
% ylim([-8000, 8000]);
% zlim([-8000, 8000]);
% legend;
% 
% text(x', y', z', GPSNameList(1:end-1));
% hold off;

%% Get random walk noise value from data
rw_stddev = GetRandomWalk(sdh_clean);

%% Align tilt data to displacment data
tilts = interp1(decyear(sdh_clean.t), sdh_clean.d, universalt);
tilts = interp1(universalt, tilts, t, 'linear');

% Tilt standard dev. is calculated by taking the mean of the E and N
% components from April 11th - April 18th 2018
tilts = (tilts - tilts(1, :));
tilts(1, :) = 0;
tiltstd = std(sdh_clean.d(1e5:1.1e5, 1:2), 1);
tiltstd = mean(tiltstd);

clear universalt

%% Manually set theta offset
% dtheta = -deg2rad(0);
% % dphi = tiltoffsets(2);
% 
% disp("Theta offset: " + dtheta)
% disp("phi offset: " + dphi)

%% Setting station locations and eliminating faulty stations 

% Only get the displacement data from the first and last timestamp for
% use in chamber optimization. u1d & tiltreduced are the long term
% displacements / tilts
finalindex = 149; % 114

% u1d = u(:, :, end-finalindex);
% 
% nanstatend = isnan(u1d(:, 2));
% 
% %%% Uncomment if we don't want to include CRIM, UWEV, BYRL %%% 
% nanstatend(find(contains(GPSNameList,'CRIM')), :, :) = 1;
% nanstatend(find(contains(GPSNameList,'UWEV')), :, :) = 1;
% nanstatend(find(contains(GPSNameList,'BYRL')), :, :) = 1;

% Delete nan stations from the data
% u1d = u1d(~nanstatend, :);
tiltreduced = tilts(end-finalindex, :) - tilts(1, :);
xopt = x(~nanstatend);
yopt = y(~nanstatend);
zopt = zeros(1, length(xopt));
% nanstatbeginning = nanstatbeginning(~nanstatend);
% offsets = 0.5 .* ones(1, sum(nanstatbeginning == 1));

% Splitting overall u into more readable ux, uy, uz
ux = squeeze(u(:, 1, :));
uy = squeeze(u(:, 2, :));
uz = squeeze(u(:, 3, :));



%% Setting the tilt location and NPIT location for optimization

% Tilt location
xyll = [-155.2937; 19.3905];
xy = llh2local(xyll, [-155.2784, 19.4073]) * 1000;
xtilt = xy(1);
ytilt = xy(2);

% NPIT location
npitloc = coord('NPIT', 'llh');
npitloc = llh2local(npitloc(1:2), [-155.2784, 19.4073]) * 1000;
vertscale = 5e3;
radscale = 5e3;

% Setting up weights for optimization
invStdPWRL = 1./std(squeeze(u(11, :, :)), 0, 2, "omitmissing");
invStdPWRL = invStdPWRL(:);

%% Optimizing SC geometry

% Use taiyi's priors except for strike angle and delta p
taiyi_parameters = [1600.79, 914.47, 90, 0, 70.3, 183, -2.18e3, -3e6, ... 
    277.01, 1621.47, 63, 136, npitloc(1) + 1890, npitloc(2) - 3030, -3630, -10e6];
taiyi_parameters_centered = [1600.79, 914.47, 90, 0, 0, 0, -2.18e3, -3e6, ... 
    277.01, 1621.47, 63, 136, npitloc(1) + 1890, npitloc(2) - 3030, -3630, -10e6];

aspect_ratio_HMM = 1;
aspect_ratio_SC = 1;
roman_HMM_vol = 9.1 * 1e9; 
roman_SC_vol = 10.5 * 1e9;
vert_sd_HMM = (3/(4*pi) * roman_HMM_vol * (aspect_ratio_HMM^2))^(1/3);
vert_sd_SC = (3/(4*pi) * roman_SC_vol * (aspect_ratio_SC^2))^(1/3);
roman_parameters = [vert_sd_HMM, vert_sd_HMM/(aspect_ratio_HMM), 90, 0, 0, 0, -900 - vert_sd_HMM, -3e6, ...
    vert_sd_SC, vert_sd_SC/(aspect_ratio_SC), 90, 136, -200, -1.9e3,-3.2e3 - vert_sd_SC, -1e6];

% ["HMM volume", "dpHMM", "vert semi-diameter", "horiz semi-diameter", "dip", "dpSC"];
lb = [4e9, -8e6, 150, 700, 40, -3e7];
ub = [9.5e9, 0, 500, 2500, 90, -1e6];
saveFigs = false;
ntrials = 1e5;
nwalkers = 100;


%%
% Make synthetic insar data
% insarx = (rand(1, 1e1) - 0.5) * 2 * 8e3;
% insary = (rand(1, 1e1) - 0.5) * 2 * 8e3;
% insarz = zeros(size(insary));
% insaru = spheroid(taiyi_parameters(1:8), [insarx; insary; insarz], 0.25, 3.08*10^9) + ...
%     spheroid(taiyi_parameters(9:end), [insarx; insary; insarz], 0.25, 3.08*10^9);
% insaru = insaru(3, :)';
insarweight = 1e10;
look = [0;0;1];

insar_data124 = readmatrix('data/track124_20180508_20180806.txt', 'Delimiter', '\t');
insarx = insar_data124(:,1);
insary = insar_data124(:,2);
insaru = insar_data124(:,3);

%% UNCOMMENT FOR MCMC
% delete(gcp('nocreate'));
[optParams, posterior] = optimize_SC_MCMC_noinsar(taiyi_parameters, lb, ub, xopt, yopt, zopt, u1d, ...
    daily_inv_std, tiltstd, nanstatend, ntrials, saveFigs);
% [optParams, posterior] = optimize_SC_MCMC(taiyi_parameters, lb, ub, xopt, ...
%     yopt, zopt, u1d, insarx, insary, insaru, look, insarweight, daily_inv_std, tiltstd, tiltreduced, nanstatend, ntrials, saveFigs);
% [optParams, posterior] = optimize_SC_MCMC_par(taiyi_parameters, lb, ub, xopt, xtilt, yopt, ytilt, ...
%     zopt, u1d, daily_inv_std, tiltstd, tiltreduced, nanstatend, ntrials, nwalkers, saveFigs);
disp(table(optParams(1),optParams(2),optParams(3),optParams(4),optParams(5),optParams(6),'VariableNames',{'HMM Vol', 'dpHMM', 'vert semi-diam', 'horiz semi-diam', 'dip', 'dpSC'}))

posterior = posterior';
aspect_ratio = 1.7496;
opt_vert_sd = (3/(4*pi) * optParams(1) * (aspect_ratio^2))^(1/3);
opt_horiz_sd = opt_vert_sd/(aspect_ratio);
optimizedM = [opt_vert_sd, opt_horiz_sd, taiyi_parameters(3:6), taiyi_parameters(7) - abs(opt_vert_sd - taiyi_parameters(1)), optParams(2:5)', 136, npitloc(1) + 1890, npitloc(2) - 3030, -3630, optParams(6)];
% optimizedM = taiyi_parameters;
disp("MCMC done");

save("figures/optimizedM.mat");

% optimizedM = taiyi_parameters;

% clear opt_vert_sd opt_horiz_sd

%% Now optimize for the tilt orientation:
% f = @(t)green_residuals_tilt(t, [optimizedM, offsets], [xopt, xtilt], [yopt, ytilt], [zopt, 0], u1d, ...
%     [invStdPWRL(1), invStdPWRL(2), invStdPWRL(3),1/((tiltstd)), 1/((tiltstd))], tiltreduced(1:2), nanstatbeginning);
% problem = createOptimProblem('fmincon','objective', f,'x0', 0,'lb',-pi/4, ...
%     'ub',pi/4,'options',options);
% ms = MultiStart;
% [dtheta, res] = run(ms,problem,20);
dtheta = deg2rad(0);
disp("Optimized theta = " + rad2deg(dtheta))

% tilts(:, 1) = 1e6 * atan(tan(1e-6 * tilts(:, 1)) * cos(dtheta) - tan(1e-6 * tilts(:, 2)) * sin(dtheta));
% tilts(:, 2) = 1e6 * atan(tan(1e-6 * tilts(:, 2)) * sin(dtheta) + tan(1e-6 * tilts(:, 2)) * cos(dtheta));

% Adding the optimized offsets to the GPS data

clearvars i j xopt yopt zopt options lb ub resnorm resdpos dsize dangle A b res;

%% Creating greens functions using the new parameters for M
% dtheta = -0.9 * pi/4; % In radians the offset angle of the tiltmeter
% dphi = 3 * pi/4;
[gHMM, gSC] = creategreens(optimizedM(1:8), optimizedM(9:end));
[gTiltHMM, gTiltSC] = createtiltgreens(optimizedM(1:8), optimizedM(9:end), dtheta, false);

%% Setting up offsets
nanstatbeginning = [isnan(ux(:, 1)); isnan(uy(:, 1)); isnan(uz(:, 1))];

%% Calculate the pressure at each time step
% Setting up green's functions and formatting displacements/tilts properly
gHMMflat = gHMM';
gHMMflat = gHMMflat(:);
gSCflat = gSC';
gSCflat = gSCflat(:);

ntime = max(size(u(1,1,:)));

tiltx = tilts(:, 1); % * (cos(dtheta) - sin(dtheta))/cos(dphi); 
tilty = tilts(:, 2); % * (cos(dtheta) + sin(dtheta))/cos(dphi);

i = 1;
dispstd = [invStdPWRL(1) .* ones(1, length(ux(:, i))), invStdPWRL(2) .* ones(1, length(ux(:, i))), invStdPWRL(3) .* ones(1, length(ux(:, i))), 1/((tiltstd)), 1/((tiltstd))];

% dispstd([1, length(ux(:, i)) + 1, length(ux(:, i))*2 + 1]) = [0, 0, 0];
dp_weight = 1e6; %% TODO - WHY DOES THIS NUMBER NOT IMPACT AS MUCH AS G
temp = TimeDependentLSQtilt(gHMMflat, gSCflat, gTiltHMM, gTiltSC, ux, uy, uz, tiltx, tilty, ...
    dispstd, GPSNameList, rw_stddev, dp_weight, true);
dp = real([temp(1:length(tiltx)); temp((length(tiltx) + 1):(2*length(tiltx)))])';
offsets = temp((2*length(tiltx) + 1):end);
ones_temp = zeros(1, length(nanstatbeginning));
ones_temp(nanstatbeginning) = offsets;
offsets = ones_temp;
nanstatbeginning = nanstatbeginning(1:nstat);
clear temp i ones_temp;

%% Modify ux, uy, uz to reflect solved offsets
offsets = reshape(offsets, [], 3);
% ux(nanstatbeginning, :) = ux(nanstatbeginning, :) + offsets(:, 1);
% uy(nanstatbeginning, :) = uy(nanstatbeginning, :) + offsets(:, 2);
% uz(nanstatbeginning, :) = uz(nanstatbeginning, :) + offsets(:, 3);

%%% NOT SURE IF I SHOULD KEEP THIS
% Fill outliers, which occur due to gaps in data throwing off LSQ suddenly
[dp(:, 2), TF] = filloutliers(dp(:, 2), "makima", "movmedian", 100, 1);
dp(TF, 1) = NaN;
dp(:, 1) = fillmissing(dp(:, 1), "makima", 1);

%% Error analysis

N_draws = 10;
N_noise = 5;
% dp_weight = 0; %2e4;
[dp_low, dp_high] = GetErrors(N_draws, N_noise, posterior, ntime, ux, uy, uz, tiltx, tilty, dispstd, ...
    GPSNameList, rw_stddev, dp_weight, taiyi_parameters, npitloc, invStdPWRL, nanstatbeginning, optParams);
 
%% Creating a matrix to store all predicted displacements + tilts in var called usim
%all errors are stored in u_error_low/high
usim = zeros(max(size(dp)), size(gHMM, 1), size(gHMM, 2));
u_error_low = zeros(max(size(dp)), size(gHMM, 1), size(gHMM, 2));
u_error_high = zeros(max(size(dp)), size(gHMM, 1), size(gHMM, 2));
uobs = zeros(max(size(dp)), size(gHMM, 1), size(gHMM, 2));

for i = 1:max(size(dp))
    % **** SETTING THE SIMULATED DISPLACEMENTS ****
    usim(i, :, 1:nstat) =  (gHMM .* dp(i, 1)) + (gSC .* dp(i, 2));
    usim(i, :, nstat + 1) = [gTiltHMM .* dp(i, 1) + gTiltSC .* dp(i, 2), 0];
end


%% Get chi^2 statistic

% Making convenient matrix for comparing observed and predicted
% displacements
uobs(:, :, 1:nstat) = permute(u(1:nstat, :, :), [3, 2, 1]);
uobs(:, :, nstat + 1) = tilts;

residuals = usim - uobs;

chi2 = 0;

for i = 1:size(usim, 2) % Component of GPS station
    for j = 1:nstat % GPS station number (tilt is the 15th element)
        if(j < nstat + 1); stddev = 1/invStdPWRL(i);
        else; stddev = tiltstd; end

        for k = 1:size(usim, 1) % Number of timesteps for each station
            % Compute contribution to chi2
            s = ((uobs(k, i, j) - usim(k, i, j))^2)/(stddev^2);
            if(~isnan(s))
                chi2 = chi2 + ((uobs(k, i, j) - usim(k, i, j))^2)/(stddev^2);
            end
        end
    end
end

%%% FIX FREE PARAM NUMBER %%%
DOF = nnz(~isnan(uobs)) - length(t)*2; % Number of non-nan observations - free params

disp("Reduced chi^2 (no tilt) = " + chi2/DOF);

%% Extract collapse amplitudes
collapset = D.gps_info.t_events(26:end);
collapset = decyear(datetime(collapset, 'ConvertFrom', 'datenum', 'Format', 'dd-MM-yy'));
% [ampHMM, ampSC] = ExtractCollapseAmplitude([dp(:, 1)'; dp(:, 2)'], t, collapset, t(3) - t(2));


%% Make synthetic insar data

insar_data124 = readmatrix('data/track124_20180508_20180806.txt', 'Delimiter', '\t');

insarbnd = [-1e4, -1e4; 1e4, 1e4];
insarxy=llh2local(insar_data124(:,1:2)', [-155.2784, 19.4073])'.*1e3;
in_bounds = (insarxy(:,1) >= insarbnd(1,1)) & (insarxy(:,1) <= insarbnd(2,1)) & ...
            (insarxy(:,2) >= insarbnd(1,2)) & (insarxy(:,2) <= insarbnd(2,2));

insardata_cropped = [insarxy(in_bounds,:), insar_data124(in_bounds, 3)];
% Filter the data
% data_filtered = data_out(in_bounds, :);

insarx = insardata_cropped(:,1);
insary = insardata_cropped(:,2);
insaru = insardata_cropped(:,3);
inc_angle = deg2rad(39);
look = [sin(inc_angle); 0; cos(inc_angle)];
%%
makeplots(x, y, z, u, u1d, ux, uy, uz, insarx, insary, insaru, look, tiltx, tilty, usim, t, nanstatend, nanstatbeginning, finalindex, collapset, dp, dp_low, dp_high, ...
    optimizedM, GPSNameList, gTiltHMM, gTiltSC, xtilt, ytilt, tiltreduced, radscale, coast_new, dtheta, 3, ntrials, offsets, saveFigs);
 
%%% TODO - Compare insar with tilt data to see if the tilt has an offsetf