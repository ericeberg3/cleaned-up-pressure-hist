function data = create_MCMC_data(m, m_guess, x, y, z, u, recstds, tilt, nanstatsbeg, priormeans, priorvariances)

npitloc = coord('NPIT', 'llh');
npitloc = llh2local(npitloc(1:2), [-155.2784, 19.4073]) * 1000;
mHMM = [m_guess(1:7), m(1)];
% mHMM = [1600.79, 914.47, 90, 0, 70.3, 183, -1940, m(1)];
mSCguess = m_guess(9:end); % [277.01, 1621.47, 63, 136, npitloc(1) + 1890, npitloc(2) - 3030, -3630, 1e7];
mSC = [m(2:4), mSCguess(4:end-1), m(end)];
% consts = m(17:end);

% tilt(1) = tilt(1) * cos(tiltoffset(1)) - sin(tiltoffset(1));
% tilt(2) = tilt(2) * cos(tiltoffset(1)) + sin(tiltoffset(1));

constsfull = zeros(length(nanstatsbeg), 3);
j = 1;
for i = 1:length(nanstatsbeg)
    if(nanstatsbeg(i) == 1)
        % constsfull(i, :) = consts(j);
        j = j + 1;
    end
end
clear j i consts

% Forward model computations
[gHMM, ~, ~, ~] = spheroid(mHMM, [x(1:end-1); y(1:end-1); z(1:end-1)], 0.25, 3.08*10^9);
[gSC, ~, ~, ~] = spheroid(mSC, [x(1:end-1); y(1:end-1); z(1:end-1)], 0.25, 3.08*10^9);

gtot = gHMM + gSC;
data = gtot';
end