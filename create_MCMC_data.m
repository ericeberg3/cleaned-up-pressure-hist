function data = create_MCMC_data(m, m_guess, x, y, z, u, insarx, insary, ...
    insaru, look, insarweight, recstds, tilt, nanstatsbeg, priormeans, priorvariances)

npitloc = coord('NPIT', 'llh');
npitloc = llh2local(npitloc(1:2), [-155.2784, 19.4073]) * 1000;

% aspect_ratio = 1.7496;
% opt_vert_sd = (3/(4*pi) * m(1) * (aspect_ratio^2))^(1/3);
% opt_horiz_sd = opt_vert_sd/(aspect_ratio);

m_full = get_full_m(m_guess, m, true);
mHMM = m_full(1:8);
mSC = m_full(9:end);

% tilt(1) = tilt(1) * cos(tiltoffset(1)) - sin(tiltoffset(1));
% tilt(2) = tilt(2) * cos(tiltoffset(1)) + sin(tiltoffset(1));

j = 1;
for i = 1:length(nanstatsbeg)
    if(nanstatsbeg(i) == 1)
        % constsfull(i, :) = consts(j);
        j = j + 1;
    end
end
clear j i consts

% Forward model computations
% GPS data generation
[gHMM, ~, ~, ~] = spheroid(mHMM, [x(1:end); y(1:end); z(1:end)], 0.25, 3.08*10^9);
[gSC, ~, ~, ~] = spheroid(mSC, [x(1:end); y(1:end); z(1:end)], 0.25, 3.08*10^9);

gtot = gHMM + gSC;
GPS_data = gtot(:);

% InSAR data generation
[gHMM, ~, ~, ~] = spheroid(mHMM, [insarx; insary; zeros(size(insarx))], 0.25, 3.08*10^9);
[gSC, ~, ~, ~] = spheroid(mSC, [insarx; insary; zeros(size(insarx))], 0.25, 3.08*10^9);

gtot = gHMM + gSC;
insar_data = gtot' * look;

% data = [insar_data];
data = [GPS_data; insar_data];
end