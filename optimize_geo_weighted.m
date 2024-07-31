%% This script uses green's functions created from an optimized HMM and SC to invert for the pressure history

% Loading in data and setting the downsample rate (from 5s to hourly)
downsamplerate = 720;
load(['data/u_mm_full_duration_' int2str(downsamplerate) '_interp.mat']);
load('data/gps_time_series.mat');
load('data/GPS_seism_locations.mat', 'GPSNameList', 'GPS_llh');
load('data/hawaii_line_new.mat', 'coast_new', 'new_pit')
load('data/sdh_clean.mat')

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
u = CenterDisps(20, u, GPSNameList, ntime);

clearvars interpu unfilledu nanu meansize rtnet_ka windowsize startcind new_pit endcind disps;

nstat = length(GPSNameList); % Number of GPS stations

% station x and y positions
xy = llh2local(GPS_llh, [-155.2784, 19.4073]);
x = xy(1,:) * 1000;
y = xy(2,:) * 1000;
z = zeros(1, nstat);

%% Align tilt data to displacment data
tilts = interp1(decyear(sdh_clean.t), sdh_clean.d, universalt);
tilts = medfilt1(tilts, downsamplerate/40, [], 1, 'omitnan', 'truncate');
tilts = interp1(universalt, tilts, t, 'linear');

% Tilt standard dev. is calculated by taking the mean of the E and N
% components from April 11th - April 18th 2018
tilts = tilts - tilts(1, :);
tilts(1, :) = 0;
tiltstd = std(sdh_clean.d(1e5:1.1e5, 1:2), 1);
tiltstd = mean(tiltstd);

%% Manually set theta offset
% dtheta = -deg2rad(0);
% % dphi = tiltoffsets(2);
% 
% disp("Theta offset: " + dtheta)
% disp("phi offset: " + dphi)


%% Cleaning up stations with too many nan values

% Removing CALS
u(end, :, :) = [];
GPSNameList(end) = [];
nstat = nstat - 1;

%% Setting station locations and eliminating faulty stations 

% Only get the displacement data from the first and last timestamp for
% use in chamber optimization. u1d & tiltreduced are the long term
% displacements / tilts
finalindex = 100;
u1d = u(:, :, end-finalindex);
tiltreduced = tilts(end-finalindex, :) - tilts(1, :);

nanstatend = isnan(u1d(:, 2));
nanstatbeginning = isnan(u(: , 1, 1));

%%% Uncomment if we don't want to include CRIM, UWEV, BYRL %%% 
nanstatend(find(contains(GPSNameList,'CRIM')), :, :) = 1;
nanstatend(find(contains(GPSNameList,'UWEV')), :, :) = 1;
nanstatend(find(contains(GPSNameList,'BYRL')), :, :) = 1;

% Delete nan stations from the data
u1d = u1d(~nanstatend, :);
xopt = x(~nanstatend);
yopt = y(~nanstatend);
zopt = zeros(1, length(xopt));
nanstatbeginning = nanstatbeginning(~nanstatend);
offsets = 0.5 .* ones(1, sum(nanstatbeginning == 1));

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
vertscale = 4e3;
radscale = 4e3;

%% Making a guess and optimizing the parameters to create an estimate of the HMM and SC geometries 

mHMMguess = [1600.79, 914.47, 90, 0, 70.3, 183, -1940, 1e7]; % [1124.79, 914.47, 90, 0, 70.3, 183, -1940, dp];
mSCguess = [277.01, 1621.47, 63, 136, npitloc(1) + 1890, npitloc(2) - 3030, -3630, 1e7];

dposen = 2e3; % meters
dposup = 1000; % meters
dsize = 4e2; % meters
dangle = 63; % degrees
deltap = 1e8; % Pa

offsetslb = -1 * ones(size(offsets));
offsetsub = 1 * ones(size(offsets));

tiltlb = [-pi/2, 0]; % east, north
tiltub = [pi/2, pi/4];

% Setting lower and upper bounds
lb = [mHMMguess(1:7), 0, mSCguess(1:2), mSCguess(3), mSCguess(4), mSCguess(5:6) - dposen, mSCguess(7) - dposup, 0, offsetslb];
ub = [mHMMguess(1:7), 3e6, mSCguess(1:2), 90, mSCguess(4), mSCguess(5:6) + dposen, mSCguess(7), mSCguess(8) + deltap, offsetsub];

% Setting up weighting scheme
invStdPWRL = 1./std(squeeze(u(11, :, :)), 0, 2, "omitmissing");
invStdPWRL = invStdPWRL(:);

% Setting function for optimization
f = @(m)green_residuals(m, [xopt, xtilt], [yopt, ytilt], [zopt, 0], u1d, ...
    [invStdPWRL(1), invStdPWRL(2), invStdPWRL(3),1/((tiltstd)), 1/((tiltstd))], tiltreduced(1:2), nanstatbeginning);

%% Optimization of chamber geometry
% This is commented out so that it doesn't run every time the code is run.
% If you'd like to re-run the optimization, uncomment lines 140-142 and
% comment 144 and 145
options = optimoptions(@fmincon, 'PlotFcn', 'optimplotfval', 'OptimalityTolerance', 1e-20, ...
    'ConstraintTolerance', 1e-20, 'MaxIter', 1000, 'MaxFunctionEvaluations', 1e4, 'StepTolerance', 1e-20); % 'display', 'iter'
% [optimizedM, res] = fmincon(f, [mHMMguess, mSCguess, offsets], [], [], [], [], lb, ub, [], options);

% Use multistart
% problem = createOptimProblem('fmincon','objective', f,'x0',[mHMMguess, mSCguess, rand * pi/4, rand * pi/4, offsets],'lb',lb, ...
%     'ub',ub,'options',options);
% ms = MultiStart;
% [optimizedM, res] = run(ms,problem,10);

optimizedM = [1600.79000000000	914.470000000000	90	0	70.3000000000000	183	-1940	2999999.99992412	277.010000000000	1621.47000000000	89.9999774942710	136	162.745313651300	-479.303502915132	-4620.30457577829	12317057.4035003];
res = 79.436;

disp("Final Residual = " + res);
offsets = optimizedM(17:end);
optimizedM = optimizedM(1:16);

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

tiltreduced = tilts(end-finalindex, :) - tilts(1, :);


% Adding the optimized offsets to the GPS data
j = 1;
for i = 1:length(nanstatbeginning)
    if(nanstatbeginning(i) == 1)
        ux(i, :) = ux(i, :) + offsets(j);
        uy(i, :) = uy(i, :) + offsets(j);
        uz(i, :) = uz(i, :) + offsets(j);
        j = j + 1;
    end
end

clearvars i j xopt yopt zopt options lb ub resnorm resdpos dsize dangle A b res;

%% Creating greens functions using the new parameters for M
% dtheta = -0.9 * pi/4; % In radians the offset angle of the tiltmeter
% dphi = 3 * pi/4;
[gHMM, gSC] = creategreens(optimizedM(1:8), optimizedM(9:end));
[gTiltHMM, gTiltSC] = createtiltgreens(optimizedM(1:8), optimizedM(9:end), dtheta, false);

%% Calculate the pressure at each time step
% Setting up green's functions and formatting displacements/tilts properly
gHMMflat = gHMM';
gHMMflat = gHMMflat(:);
gSCflat = gSC';
gSCflat = gSCflat(:);
gTiltHMMflat = gTiltHMM';

ntime = max(size(u(1,1,:)));
dpHMM = zeros(1, ntime);
dpSC = zeros(1, ntime);
dp = [dpHMM', dpSC'];
tiltx = tilts(:, 1); % * (cos(dtheta) - sin(dtheta))/cos(dphi); 
tilty = tilts(:, 2); % * (cos(dtheta) + sin(dtheta))/cos(dphi);

i = 1;
dispstd = [invStdPWRL(1) .* ones(1, length(ux(:, i))), invStdPWRL(2) .* ones(1, length(ux(:, i))), invStdPWRL(3) .* ones(1, length(ux(:, i))), 1/((tiltstd)), 1/((tiltstd))];
temp = TimeDependentLSQtilt(gHMMflat, gSCflat, gTiltHMM, gTiltSC, ux, uy, uz, tiltx, tilty, dispstd, GPSNameList);
dp = real([temp(1:length(tiltx)); temp(length(tiltx) + 1:end)])';
clear temp i;

% Fill outliers, which occur due to gaps in data throwing off LSQ suddenly
[dp(:, 2), TF] = filloutliers(dp(:, 2), "makima", "movmedian", 100, 1);
dp(TF, 1) = NaN;
dp(:, 1) = fillmissing(dp(:, 1), "makima", 1);

%% Creating a matrix to store all predicted displacements + tilts in var called usim
usim = zeros(max(size(dp)), size(gHMM, 1), size(gHMM, 2));

for i = 1:max(size(dp))
    % **** SETTING THE SIMULATED DISPLACEMENTS ****
    usim(i, :, 1:nstat) =  (gHMM .* dp(i, 1)) + (gSC .* dp(i, 2)) ;
    usim(i, :, nstat + 1) = [gTiltSC .* dp(i, 1) + gTiltHMM .* dp(i, 2), 0];
end

%% Extract collapse amplitudes
collapset = D.gps_info.t_events(26:end);
collapset = decyear(datetime(collapset, 'ConvertFrom', 'datenum', 'Format', 'dd-MM-yy'));
% [ampHMM, ampSC] = ExtractCollapseAmplitude([dp(:, 1)'; dp(:, 2)'], t, collapset, t(3) - t(2));

%% Plots
% Convert time into matlab dateyear
year = floor(t);
partialYear = mod(t,1);
date0 = datenum(num2str(year),'yyyy');
date1 = datenum(num2str(year+1),'yyyy');
daysInYear = date1 - date0;
t = date0 + partialYear .* daysInYear;
t = datetime(t, 'ConvertFrom', 'datenum', 'Format', 'dd-MM-yy');
collapset = datetime(collapset, 'ConvertFrom', 'datenum', 'Format', 'dd-MM-yy');

%% 
% Pressure vs. time figure
figure(10);
plot(t(1:end-finalindex), dp(1:end-finalindex, 1) * optimizedM(8)/(1e6), 'DisplayName', 'HMM Pressure');
hold on
plot(t(1:end-finalindex), dp(1:end-finalindex, 2) * optimizedM(16)/(1e6), 'DisplayName', 'SC Pressure');
l2 = xline(t(end - finalindex), 'red', 'Label', 'Time', 'HandleVisibility', 'off');
title("Pressure vs. Time");
xlabel('Time');
ylabel('Pressure (MPa)');
% xlim([collapset(3) collapset(5)]);
legend();
hold off

%% Making grid of displacements and tilt
figure(7);
disptype = 1; % 1 = x, 2 = y, 3 = z

tlo = tiledlayout(4,4);
if(disptype == 1)
    title(tlo, "East Displacement vs. Time", 'FontSize', 24);
elseif(disptype == 2)
    title(tlo, "North Displacement vs. Time", 'FontSize', 24);
else
    title(tlo, "Vertical Displacement vs. Time", 'FontSize', 24);
end

for i = 1:max(size(GPSNameList) + 1)
    nexttile
    if(i < length(GPSNameList) + 1)
        if(disptype == 1)
            plot(t(1:end-finalindex), ux(i, 1:end-finalindex), '-', 'DisplayName', 'GPS', 'LineWidth', 1.2);
            hold on;
            plot(t(1:end-finalindex), usim(1:end-finalindex, 1, i), 'DisplayName', 'LSQ', 'LineWidth', 1.6);
            plot(t(end - finalindex), ux(i, end - finalindex), 'o-', 'MarkerFaceColor','red');
        elseif(disptype == 2)
            plot(t(1:end-finalindex), uy(i, 1:end-finalindex), '-', 'DisplayName', 'GPS', 'LineWidth', 1.2);
            hold on;
            plot(t(1:end-finalindex), usim(1:end-finalindex, 2, i), 'DisplayName', 'LSQ', 'LineWidth', 1.6);
            plot(t(end - finalindex), uy(i, end - finalindex), 'o-', 'MarkerFaceColor','red');
        else
            plot(t(1:end-finalindex), uz(i, 1:end-finalindex), '-', 'DisplayName', 'GPS', 'LineWidth', 1.2);
            hold on;
            plot(t(1:end-finalindex), squeeze(usim(1:end-finalindex, 3, i)), 'DisplayName', 'LSQ', 'LineWidth', 1.6);
            plot(t(end - finalindex), uz(i, end - finalindex), 'o-', 'MarkerFaceColor','red');
            % xline(collapset);
        end
        if(GPSNameList(i) == "UWEV" || GPSNameList(i) == "BYRL" || GPSNameList(i) == "CRIM")
                set(gca,'Color','k');
        end
        ylabel("Displacement (m)");
        title(GPSNameList(i))
    else
        simtiltx = (gTiltHMM(1) .* dp(:, 1)) + (gTiltSC(1) .* dp(:, 2));
        plot(t(1:end-finalindex),tiltx(1:end-finalindex), '-', 'DisplayName', 'Data', 'LineWidth', 1.2);
        hold on;
        plot(t(1:end-finalindex), simtiltx(1:end-finalindex), '-', 'DisplayName', 'LSQ', 'LineWidth', 1.6);
        title("Tilt e");
        ylim([-30, 250]);
        ylabel("Tilt (µrad)")
        hold off;
        nexttile;
        simtilty = (gTiltHMM(2) .* dp(:, 1)) + (gTiltSC(2) .* dp(:, 2));
        plot(t(1:end-finalindex), tilty(1:end-finalindex), '-', 'DisplayName', 'Data', 'LineWidth', 1.2);
        hold on;
        plot(t(1:end-finalindex), simtilty(1:end-finalindex), 'DisplayName', 'LSQ', 'LineWidth', 1.6);
        ylim([-90, 130]);
        ylabel("Tilt (µrad)")
        title("Tilt n");
    end
    hold off;
end
leg = legend('Orientation', 'Horizontal');
leg.Layout.Tile = 'north';
leg.FontSize = 14;

%% Print out statistics

disp("Net HMM dp: " + (dp(1, 1) - dp(end - finalindex, 1)))
disp("Net SC dp: " + (dp(1, 2) - dp(end - finalindex, 2)))

GPSrms = u - permute(usim(:, :, 1:14), [3, 2, 1]);
GPSrms = GPSrms(:);
GPSrms = rms(GPSrms, 'omitnan');
disp("GPS RMS Misfit: " + GPSrms)

tiltrms = rms(simtiltx - tiltx) + rms(simtilty - tilty);
disp("Tilt Unweighted RMS Misfit: " + tiltrms);

clear GPSrms tiltrms
%% Quiver Plot

mHMM = optimizedM(1:8);
mSC = optimizedM(9:end);

figure(4);
spheroid(mHMM);
hold on;
spheroid(mSC);
hold off;

[gHMM, ~, ~, ~] = spheroid(mHMM, [x; y; z(1:length(x))], 0.25, 3.08*10^9);
[gSC, ~, ~, ~] = spheroid(mSC, [x; y; z(1:length(x))], 0.25, 3.08*10^9);

[~, dHMM, ~, ~] = spheroid(mHMM, [xtilt; ytilt; 0], 0.25, 3.08*10^9);
[~, dSC, ~, ~] = spheroid(mSC, [xtilt; ytilt; 0], 0.25, 3.08*10^9);

tiltscale = 1e-3;
u1d = u(:, :, end-100);
u1d(end + 1, 1) = tiltreduced(1) .* tiltscale;
u1d(end, 2) = tiltreduced(2) .* tiltscale;
u1d(end, 3) = zeros(length(tiltreduced(1)), 1);
x(end+1) = xtilt;
y(end+1) = ytilt;

gTilt = createtiltgreens(mHMM, mSC, dtheta, false);
gtot = gHMM + gSC;
gtot(:, end + 1) = [gTilt, 0] .* tiltscale;


hold off;
figure(6);
realquiver = quiver3(x', y', zeros(size(x))', u1d(:, 1) * radscale, u1d(:, 2) * radscale, u1d(:, 3) * radscale, 'AutoScale', 'off', 'LineWidth',2.75, 'MaxHeadSize', 0.3, 'Color', '#4DBEEE', 'DisplayName', 'Data');
hold on;
plot3(mSCguess(5),mSCguess(6), mSCguess(7), '.', 'MarkerSize', 20, 'Color', '#f77036', 'HandleVisibility','off');
quiver3(mSCguess(5),mSCguess(6), mSCguess(7), optimizedM(13) - mSCguess(5), optimizedM(14) - mSCguess(6), optimizedM(15) - mSCguess(7), 'AutoScale', 'off', 'LineWidth', 2.75, 'MaxHeadSize', 0.5, 'Color', '#f77036', 'DisplayName','SC Center Shift');
simquiver = quiver3(x', y', zeros(size(x))', gtot(1, :)' * radscale, gtot(2, :)' * radscale, gtot(3, :)' * radscale, 'AutoScale', 'off', 'LineWidth',2.75, 'MaxHeadSize', 0.3, 'Color', '#A2142F', 'DisplayName','Optimization Result');
xlabel('x (m)');
ylabel('y (m)');
zlabel('Scaled Displacement (m)');
title("Displacements at various stations (y oriented N, x oriented E)")
xlim([-8000, 8000]);
ylim([-8000, 8000]);
zlim([-8000, 8000]);
GPSNameList(end + 1) = "Tilt";

cxy = llh2local(coast_new', [-155.2784, 19.4073]);
cxy = cxy * 1000;
plot(cxy(1, :)', cxy(2, :)', 'k.', 'HandleVisibility','off');

legend;

text(x', y', z', GPSNameList);


hold off;
