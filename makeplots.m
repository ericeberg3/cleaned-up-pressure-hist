function makeplots(x, y, z, u, ux, uy, uz, tiltx, tilty, usim, t, finalindex, collapset, ...
    dp, optimizedM, GPSNameList, gTiltHMM, gTiltSC, xtilt, ytilt, tiltreduced, radscale, mSCguess, coast_new, dtheta)
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
    disptype = 2; % 1 = x, 2 = y, 3 = z
    
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
    
    % [gHMM, ~, ~, ~] = spheroid(mHMM, [x; y; z(1:length(x))], 0.25, 3.08*10^9);
    % [gSC, ~, ~, ~] = spheroid(mSC, [x; y; z(1:length(x))], 0.25, 3.08*10^9);

    [gHMM, gSC] = creategreens(mHMM, mSC);
    [gTiltHMM, gTiltSC] = createtiltgreens(mHMM, mSC, 0, false);
    
    tiltscale = 1e-3;
    u1d = u(:, :, end-finalindex);
    u1d(end + 1, 1) = tiltreduced(1) .* tiltscale;
    u1d(end, 2) = tiltreduced(2) .* tiltscale;
    u1d(end, 3) = zeros(length(tiltreduced(1)), 1);
    x(end+1) = xtilt;
    y(end+1) = ytilt;
   
    gtot = gHMM + gSC;
    gtot(:, end + 1) = [(gTiltHMM + gTiltSC), 0] .* tiltscale;
    
    
    hold off;
    figure(6);
    clf;
    realquiver = quiver3(x', y', zeros(size(x))', u1d(:, 1) * radscale, u1d(:, 2) * radscale, u1d(:, 3) * radscale, 'AutoScale', 'off', 'LineWidth',2.75, 'MaxHeadSize', 0.3, 'Color', '#4DBEEE', 'DisplayName', 'Data');
    hold on;
    plot3(mSCguess(5),mSCguess(6), mSCguess(7), '.', 'MarkerSize', 20, 'Color', '#f77036', 'HandleVisibility','off');
    quiver3(mSCguess(5),mSCguess(6), mSCguess(7), optimizedM(13) - mSCguess(5), optimizedM(14) - mSCguess(6), optimizedM(15) - mSCguess(7), 'AutoScale', 'off', 'LineWidth', 2.75, 'MaxHeadSize', 0.5, 'Color', '#f77036', 'DisplayName','SC Center Shift');
    simquiver = quiver3(x', y', zeros(size(x))', gtot(1, :)' * radscale, gtot(2, :)' * radscale, gtot(3, :)' * radscale, 'AutoScale', 'off', 'LineWidth',2.75, 'MaxHeadSize', 0.3, 'Color', '#A2142F', 'DisplayName','Optimization Result');
    quiver3(x(end), y(end), 0, gTiltHMM(1) * 3e1, gTiltHMM(2) * 3e1, 0, 'AutoScale', 'off', 'LineWidth',2.75, 'MaxHeadSize', 0.3, 'Color', '#023799', 'DisplayName','HMM Greens fcn');
    quiver3(x(end), y(end), 0, gTiltSC(1) * 1e1, gTiltSC(2) * 1e1, 0, 'AutoScale', 'off', 'LineWidth',2.75, 'MaxHeadSize', 0.3, 'Color', '#910299', 'DisplayName','SC Greens fcn');
    % ellipsoid(mHMM(5),mHMM(6),mHMM(7),mHMM(2) * sind(mHMM(3)) * cosd(mHMM(4)),mHMM(2) * sind(mHMM(3)) * cosd(mHMM(4)),mHMM(1) * sind(mHMM(3)))
    % ellipsoid(0, 0, -1940, 1000, 1000, 1000)
    
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


end