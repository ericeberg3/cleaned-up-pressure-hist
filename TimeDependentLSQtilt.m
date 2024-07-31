%% This function calculates the pressure at all timesteps using a single large least squares step

function dp = TimeDependentLSQtilt(gHMMflat, gSCflat, gTiltHMM, gTiltSC, ux, uy, uz, tiltx, tilty, recstds, ...
    GPSNameList)
    gtotHMM = cat(1, gHMMflat);
    gtotSC = cat(1, gSCflat);
    
    % Getting rid of nan stations
    nanstations = [isnan(ux(:, 1200)); isnan(ux(:, 1200)); isnan(ux(:, 1200)); isnan(tiltx(1200)); isnan(tilty(1200))];
    nanstations(find(contains(GPSNameList,'UWEV')):length(GPSNameList):end) = 1;
    nanstations(find(contains(GPSNameList,'BYRL')):length(GPSNameList):end) = 1;
    nanstations(find(contains(GPSNameList,'CRIM')):length(GPSNameList):end) = 1;
    nanstations = logical(nanstations);
    
    gtotHMM(nanstations) = [];
    gtotSC(nanstations) = [];
    gtotTilt = cat(1, gTiltHMM', gTiltSC');
    nstat = size(ux, 1);
    ux(nanstations(1:nstat), :) = [];
    uy(nanstations(nstat+1:nstat*2), :) = [];
    uz(nanstations(nstat*2+1:nstat*3), :) = [];

    % Set up weight matrix  & green's functions matrix G
    % Wdiag = the diagonal of the GPS components of the Sigma^-1 matrix
    % tiltW = the tilt block of the Sigma^-1 matrix
    ntime = length(ux);
    greensize = length(gtotHMM);
    tiltsize = length(gTiltHMM);
    G = zeros((greensize + 2) * ntime, ntime * 2);
    ureal = zeros(1, (greensize + 2) * ntime);
    Wdiag = ones(size(greensize*ntime, 2), 1);
    tiltW = zeros(tiltsize*ntime, tiltsize*ntime);

    for i = 1:ntime
        % Fill in GPS stations greens functions into G and weights into
        % Wdiag
        G((((i-1) * greensize) + 1):(i * greensize), i) = gtotHMM';
        G(((i-1) * greensize + 1):(i * greensize), i + ntime) = gtotSC';
        ureal((((i -1) * greensize) + 1):((i - 1) * greensize + greensize)) = cat(1, ux(:, i), uy(:, i), uz(:, i));
        sizegps = greensize;
        Wdiag((greensize*(i - 1) + 1) : (greensize*(i - 1) + sizegps/3)) = recstds(end - 3*(length(recstds) - 2)/3); % weight x component
        Wdiag((greensize*(i - 1) + 1 + sizegps/3) : (greensize*(i - 1) + 2*sizegps/3)) = recstds(end - 2*(length(recstds) - 2)/3); % weight y component
        Wdiag((greensize*(i - 1) + 1 + 2*sizegps/3) : (greensize*(i - 1) + 3*sizegps/3)) = recstds(end - 1*(length(recstds) - 2)/3); % weight z component

        % Fill in tilt stations greens functions into G and weights into
        % tiltW
        gpsendind = greensize*ntime + 1;
        % G((gpsendind + tiltsize*(i-1)):(gpsendind + tiltsize*(i-1) + 1), i) = gTiltHMM';
        % G((gpsendind + tiltsize*(i-1)):(gpsendind + tiltsize*(i-1) + 1), i + ntime) = gTiltSC';
        G((gpsendind + i-1), i) = gTiltHMM(1);
        G((gpsendind + length(tiltx) + i - 1), i) = gTiltHMM(2);
        G((gpsendind + i-1), i + ntime) = gTiltSC(1);
        G((gpsendind + length(tiltx) + i-1), i + ntime) = gTiltSC(2);

        ureal(gpsendind + i-1) = tiltx(i);
        ureal(gpsendind + length(tiltx) + i-1) = tilty(i);
        % ureal((gpsendind + tiltsize*(i-1)):(gpsendind + tiltsize*(i-1) + 1)) = [tiltx(i), tilty(i)];
        deltat = 1; % in hours
        % The random walk standard deviation sigma that I am using is the
        % standard deviation from the no-signal period (incorrect, change
        % later)
        sigma = 0.8;
        % Sigma_{i, i}
        tiltW(i, i) = 1/(2*sigma^2) * (2/deltat);
        tiltW(i + ntime, i+ntime) = 1/(2*sigma^2) * (2/deltat);

        % % Sigma_{i - 1, i}
        tiltW(i -1 + ntime, i + ntime) = -1/(2*sigma^2) * (1/deltat);
        if(i > 1)
            tiltW(i-1, i) = -1/(2*sigma^2) * (1/deltat);
        end

        % Sigma_{i + 1, i}
        tiltW(i + 1, i) = -1/(2*sigma^2) * (1/deltat);
        if(i < ntime)
            tiltW(i + ntime + 1, i + ntime) = -1/(2*sigma^2) * (1/deltat);
        end
    end
    disp("RW Sigma = " + sigma)
    % make the random walk "start over"
    tiltW(ntime, ntime + 1) = 0;
    tiltW(ntime + 1, ntime) = 0;
    % Changing the tilt weighting 
    % tiltW = tiltW;

    % Filling in the Sigma^-1 matrix, called W
    W = zeros(size(ureal, 2), size(ureal, 2));
    W(1:length(Wdiag), 1:length(Wdiag)) = diag(Wdiag);
    W((length(Wdiag) + 1):end, (length(Wdiag) + 1):end) = tiltW;
    
    W = sparse(W);
    G = sparse(G);

    ureal = ureal';

    % Clean up NaN values in preparation of least squares inversion
    G(isnan(ureal), :) = [];
    W(isnan(ureal), :) = [];
    W(:, isnan(ureal)) = [];
    ureal(isnan(ureal)) = [];

    % Doing the weighted least squares
    A = (G' * W * G);
    dp = (A \ (G' * W * ureal))';
end