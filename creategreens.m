function [gHMM, gSC] = creategreens(mHMM, mSC)
    load('data/GPS_seism_locations.mat', 'GPSNameList', 'GPS_llh');
    
    xy = llh2local(GPS_llh, [-155.2784, 19.4073]);
    x = xy(1,:) * 1000;
    y = xy(2,:) * 1000;
    z = zeros(size(x));
    
    % Sign convention negative dp is positive pressure increase
    
    [gHMM, ~, ~, ~] = spheroid(mHMM, [x; y; z], 0.25, 3.08*10^9);
    [gSC, ~, ~, ~] = spheroid(mSC, [x; y; z], 0.25, 3.08*10^9);
    
    gHMM = real(gHMM);
    gSC = real(gSC);
end
