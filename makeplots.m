function makeplots(x, y, z, u, u1d, ux, uy, uz, insarx, insary, insaru, look, tiltx, tilty, usim, t, nanstat, nanstatbeginning, finalindex, collapset, ...
    dp, dp_low, dp_high, optimizedM, GPSNameList, gTiltHMM, gTiltSC, xtilt, ytilt, tiltreduced, radscale, ...
    coast_new, dtheta, disptype, ntrials, offsets, saveFigs)

    u1d = squeeze(u(:, :, end-finalindex));
    nanstat = isnan(u1d(:, 1));
    u1d = u1d(~nanstat, :);
    %% Plotting with solved offsets
    ux = ux(:, :) + offsets(:, 1);
    uy = uy(:, :) + offsets(:, 2);
    uz = uz(:, :) + offsets(:, 3);

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
    
    %% Pressure vs. time figure with dual y-axes
    figure(10); 
    clf;
    
    % Define scale factor and offset for SC pressure (adjust these as needed)
    sc_scale = 2;   % Scaling factor for SC data
    sc_offset = 5;  % Downward offset for SC data
    
    % Plot HMM pressure on the right y-axis
    ax1 = axes;  % This axis will host the HMM data and the transformed SC data
    hold(ax1, 'on');
    set(ax1, 'YColor', '#0072BD')

    h1 = plot(ax1, t(1:end-finalindex), dp(1:end-finalindex, 1) * optimizedM(8)/(1e6), ...
        'Color', '#0072BD', 'DisplayName', 'HMM Pressure', 'LineWidth', 4);
    ylabel('HMM Pressure (MPa)');
    hold on;
    
    % Plot HMM confidence intervals on the right axis
    xPatch = [t(1:end-finalindex); flipud(t(1:end-finalindex))];
    yPatch_HMM = [dp_low(1:end-finalindex, 1); flipud(dp_high(1:end-finalindex, 1))];

    patch(ax1, xPatch, yPatch_HMM, 'blue', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % HMM conf interval
    
    % Adjust the SC pressure data with scaling and offset
    sc_pressure = dp(1:end-finalindex, 2)*optimizedM(16)/(1e6);
    sc_pressure_adj = sc_scale * sc_pressure - sc_offset;
    h2 = plot(ax1, t(1:end-finalindex), sc_pressure_adj, 'DisplayName', 'SC Pressure', 'LineWidth', 4);
    
    % Plot SC confidence intervals on the left axis
    yPatch_SC = [sc_scale*(dp_low(1:end-finalindex, 2)) - sc_offset; ...
                 flipud(sc_scale*(dp_high(1:end-finalindex, 2)) - sc_offset)];
    patch(ax1, xPatch, yPatch_SC, 'red', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    ylim([-20, 5]);
    
    % Create axes label for SC
    ax1_pos = get(ax1, 'Position');
    ax1_pos(1) = ax1_pos(1) - 0.05;
    ax2 = axes('Position', ax1_pos, ...
           'Color', 'none', ...
           'YAxisLocation', 'left', ...
           'XAxisLocation', 'top', ... 
           'XTick', [], ...       
           'Box', 'off');

    set(ax2, 'YLim', [-20, 5], 'YColor','#D95319');
    ax2.XAxis.Visible = 'off';
    
    ticks = get(ax2, 'YTick');
    tickLabels = (ticks + sc_offset) / sc_scale;
    set(ax2, 'YTickLabel', tickLabels);
    ylabel(ax2, 'SC Pressure (MPa)');

    % Set common figure properties
    set(ax2, 'FontSize', 14);
    set(ax1, 'FontSize', 14);
    title("Pressure vs. Time");
    xlabel('Time');
    legend([h1, h2], 'Location', 'northeast', "FontSize", 18);
    hold off
    
    if(saveFigs)
        saveas(10, "Figures/pvt_" + num2str(ntrials, "%.1e") + "trials.fig");
    end

    %% Making grid of displacements and tilt
    for j = 3
        % disptype = j;
        figure(7);
        % disptype = 2; % 1 = x, 2 = y, 3 = z
        
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

                    plot([0, 0], [0, offsets(i, j)], "LineWidth", 8);
                elseif(disptype == 2)
                    plot(t(1:end-finalindex), uy(i, 1:end-finalindex), '-', 'DisplayName', 'GPS', 'LineWidth', 1.2);
                    hold on;
                    plot(t(1:end-finalindex), usim(1:end-finalindex, 2, i), 'DisplayName', 'LSQ', 'LineWidth', 1.6);
                    plot(t(end - finalindex), uy(i, end - finalindex), 'o-', 'MarkerFaceColor','red');

                    plot([0, 0], [0, offsets(i, j)], "LineWidth", 8);
                else
                    plot(t(1:end-finalindex), uz(i, 1:end-finalindex), '-', 'DisplayName', 'GPS', 'LineWidth', 1.2);
                    hold on;
                    plot(t(1:end-finalindex), squeeze(usim(1:end-finalindex, 3, i)), 'DisplayName', 'LSQ', 'LineWidth', 1.6);
                    plot(t(end - finalindex), uz(i, end - finalindex), 'o-', 'MarkerFaceColor','red');

                    plot([0, 0], [0, offsets(i, j)], "LineWidth", 8);
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
        if(saveFigs); saveas(7, "Figures/displacements_" + disptype + "_" + num2str(ntrials, "%.1e") + "trials.fig"); end;
    end
    
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
    
    %% Geometry plot 
    mHMM = optimizedM(1:8);
    mSC = optimizedM(9:end);
    % 
    % figure(4);
    % clf;
    % spheroid(mHMM);
    % hold on;
    % spheroid(mSC);
    % hold off;
    % if(saveFigs); saveas(4, "Figures/geo_" + num2str(ntrials, "%.1e") + "trials.fig"); end

    %% Quiver Plot
    % [gHMM, ~, ~, ~] = spheroid(mHMM, [x; y; z(1:length(x))], 0.25, 3.08*10^9);
    % [gSC, ~, ~, ~] = spheroid(mSC, [x; y; z(1:length(x))], 0.25, 3.08*10^9);

    [gHMM, gSC] = creategreens(mHMM, mSC);
    [gTiltHMM, gTiltSC] = createtiltgreens(mHMM, mSC, 0, false);
    
    tiltscale = 1e-3;
    u1d(end + 1, 1) = tiltreduced(1) .* tiltscale;
    u1d(end, 2) = tiltreduced(2) .* tiltscale;
    u1d(end, 3) = zeros(length(tiltreduced(1)), 1);
    x(end+1) = xtilt;
    y(end+1) = ytilt;
    nanstat(end + 1) = false;
   
    gtot = gHMM + gSC;
    gtot(:, end + 1) = [(gTiltHMM + gTiltSC), 0] .* tiltscale;

    % Set up optimization result vectors from LSQ solution
    u1d_LSQ = usim(end-finalindex, :, :);
    u1d_LSQ = squeeze(u1d_LSQ(1, :, :))';
    u1d_LSQ(end, 1) = u1d_LSQ(end, 1) .* tiltscale;
    u1d_LSQ(end, 2) = u1d_LSQ(end, 2) .* tiltscale;
    
    mSCguess = optimizedM(9:end);
    hold off;
    figure(6);
    clf;
    hold on;
    %%% IF WE WANT 3D Plot change the scaling of Z back to radscale
    realquiver = quiver3(x(~nanstat)', y(~nanstat)', zeros(size(x(~nanstat)))', u1d(:, 1) * radscale, u1d(:, 2) * radscale, ...
        u1d(:, 3) * 0, 'AutoScale', 'off', 'LineWidth',2.75, 'MaxHeadSize', 0.3, 'Color', '#4DBEEE', 'DisplayName', 'Data');
    % quiver3(mSCguess(5),mSCguess(6), mSCguess(7), optimizedM(13) - mSCguess(5), optimizedM(14) - mSCguess(6), optimizedM(15) - mSCguess(7), 'AutoScale', 'off', 'LineWidth', 2.75, 'MaxHeadSize', 0.5, 'Color', '#f77036', 'DisplayName','SC Center Shift');
    simquiver = quiver3(x', y', zeros(size(x))', u1d_LSQ(:, 1) * radscale, u1d_LSQ(:, 2) * radscale, ...
        u1d_LSQ(:, 3) * 0, 'AutoScale', 'off', 'LineWidth',2.75, 'MaxHeadSize', 0.3, 'Color', '#A2142F', 'DisplayName','Optimization Result');

    % Plot vert. displacement circles:
    theta = linspace(0, 2*pi, 50);  
    verticalScale = 0.3e4; 
    u1d_ind = 1;
    for k = 1:length(u1d_LSQ)
        if(~nanstat(k))
            line_style = '--';
            if(u1d(u1d_ind,3) < 0); line_style = '-'; end

            r = abs(u1d(u1d_ind,3)) * verticalScale;
            x_data = x(k) + r*cos(theta);
            y_data = y(k) + r*sin(theta);

            plot(x_data, y_data, line_style, 'Color', '#4DBEEE', 'LineWidth', 3, 'HandleVisibility','off');
            u1d_ind = u1d_ind + 1;
        end

        line_style = '--';
        if(u1d_LSQ(k,3) < 0); line_style = '-'; end

        r_lsq = abs(u1d_LSQ(k,3)) * verticalScale;

        x_sim = x(k) + r_lsq*cos(theta);
        y_sim = y(k) + r_lsq*sin(theta);

        plot(x_sim, y_sim, line_style, 'Color', '#A2142F', 'LineWidth', 2.5, 'HandleVisibility','off');
    end

    % Add reference circle for scale
    % Define reference location and scale
    x_ref = 6000;
    y_ref = 6000;
    r_ref = 0.1 * verticalScale;   % Reference vertical displacement magnitude (e.g., 1 unit)
    x_ref_circle = x_ref + r_ref*cos(theta);
    y_ref_circle = y_ref + r_ref*sin(theta);
    
    % Plot the reference circle
    plot(x_ref_circle, y_ref_circle, 'r', 'LineWidth', 2, 'HandleVisibility','off');
    text(x_ref - 270, y_ref - r_ref - 200, "10 cm", 'Color', 'red', 'FontSize', 12, 'HandleVisibility','off');
    quiver(x_ref, y_ref, -radscale * 0.5, 0, 'LineWidth',2.75, 'MaxHeadSize', 0.3, 'Color', 'red', 'HandleVisibility','off'); % Vector scale reference
    text(x_ref - radscale*0.4, y_ref -radscale*0.1, "50 cm", 'Color', 'red', 'FontSize', 12);

    % Plot ellipsoids
    for i = 1:2
        if(i == 1); m = mHMM; name = "HMM"; col = "#ffa07d"; end
        if(i == 2); m = mSC; name = "SC"; col = "#ff7d7d"; end

        % Define ellipsoid parameters
        x_c = m(5); % x-coordinate of the ellipsoid center
        y_c = m(6); % y-coordinate of the ellipsoid center
        z_c = m(7); % z-coordinate of the ellipsoid center
        a = m(2); % Semi-axis along the x-direction
        b = m(1); % Semi-axis along the y-direction
        c = m(2); % Semi-axis along the z-direction
        
        % Generate ellipsoid points
        [x_ell, y_ell, z_ell] = ellipsoid(0, 0, 0, a, b, c, 50);


        % Flatten the matrices into vectors for transformation
        points = [x_ell(:), y_ell(:), z_ell(:)]';
        
        % Define rotation angles (degrees)
        strike = m(4); % Rotation around z-axis (strike angle)
        dip = m(3);    % Rotation around x-axis (dip angle)
        
        % Convert angles to radians
        strike_rad = deg2rad(strike);
        dip_rad = deg2rad(dip);
        
        % Rotation matrices
        Rz = [cos(strike_rad), sin(strike_rad), 0; % Rotation around z-axis
              -sin(strike_rad), cos(strike_rad),  0;
              0,              0,               1];
        Rx = [1, 0,               0;              % Corrected rotation around x-axis
              0, cos(dip_rad), -sin(dip_rad);
              0, sin(dip_rad),  cos(dip_rad)];
        
        % Combined rotation matrix
        R = Rx * Rz;
        
        % Rotate points
        rotated_points = R * points;
        
        % Reshape rotated points back to grid format
        x_rot = reshape(rotated_points(1, :), size(x_ell)) + x_c;
        y_rot = reshape(rotated_points(2, :), size(y_ell)) + y_c;
        z_rot = reshape(rotated_points(3, :), size(z_ell)) + z_c;
        
        % Plot the rotated ellipsoid
        surf(x_rot, y_rot, z_rot, 'EdgeColor', 'none', 'FaceColor', col, 'FaceAlpha', 0.8, 'DisplayName', name);

        hold on;
    end

    % Add lighting effects
    light('Position', [0, 0, 200], 'Style', 'infinite'); % Add light source
    lighting gouraud; % Smooth shading
    material shiny; % Enhance reflectivity for a polished appearance

    xlabel('Easting (m)', "FontSize", 18);
    ylabel('Northing (m)', "FontSize", 18);
    zlabel('Vertical (m)',  "FontSize", 18);
    title("Displacements at various stations (y oriented N, x oriented E)", "FontSize", 18)
    xlim([-8000, 8000]);
    ylim([-8000, 8000]);
    zlim([-8000, 8000]);
    GPSNameList(end + 1) = "SDH";
    
    cxy = llh2local(coast_new', [-155.2784, 19.4073]);
    cxy = cxy * 1000;
    plot(cxy(1, :)', cxy(2, :)', 'k.', 'HandleVisibility','off');
    
    text((x-300)', (y+300)', (zeros(size(x)) + 400)', GPSNameList, 'FontSize', 16, 'HandleVisibility','off');
    plot3(x', y', zeros(size(x))', '.', "Color", "#f76a05", "MarkerSize", 20, 'HandleVisibility','off')
    hold off;
    legend("FontSize", 18, "Location", "southwest");

    if(saveFigs); saveas(6, "Figures/quiver_" + num2str(ntrials, "%.1e") + "trials.fig"); end


    %% InSAR Plot
    
    figure(8);
    clf;
    hold on;
    % Plot insar data
    % Create a regular grid from your x and y data (assuming insarx and insary are your original arrays)
    xi = linspace(min(insarx), max(insarx), 200);
    yi = linspace(min(insary), max(insary), 200);
    [X, Y] = meshgrid(xi, yi);
    
    % Interpolate vertical displacement data (assumed to be in insaru(3,:))
    insaru_grid = griddata(insarx, insary, insaru(3, :), X, Y, 'cubic');
    
    % Plot using pcolor for a 2D color map
    pcolor(X, Y, insaru_grid);
    shading interp;           % Smooth out color transitions
    colormap jet;             % Set the colormap (try other maps if needed)
    colorbar;                 % Add a colorbar for reference
    
    % Optionally, adjust the color limits to reveal gradients
    clim([min(insaru_grid(:)), max(insaru_grid(:))]);

    % Plot GPS on top
    realquiver = quiver(x(~nanstat)', y(~nanstat)', u1d(:, 1) * radscale, u1d(:, 2) * radscale,  ...
        'AutoScale', 'off', 'LineWidth',2.75, 'MaxHeadSize', 0.3, 'Color', '#4DBEEE', 'DisplayName', 'Data');
    % quiver3(mSCguess(5),mSCguess(6), mSCguess(7), optimizedM(13) - mSCguess(5), optimizedM(14) - mSCguess(6), optimizedM(15) - mSCguess(7), 'AutoScale', 'off', 'LineWidth', 2.75, 'MaxHeadSize', 0.5, 'Color', '#f77036', 'DisplayName','SC Center Shift');
    simquiver = quiver(x', y', u1d_LSQ(:, 1) * radscale, u1d_LSQ(:, 2) * radscale, ...
        'AutoScale', 'off', 'LineWidth',2.75, 'MaxHeadSize', 0.3, 'Color', '#A2142F', 'DisplayName','Optimization Result');

    % Plot vert. displacement circles:
    theta = linspace(0, 2*pi, 50);  
    verticalScale = 0.3e4; 
    u1d_ind = 1;
    for k = 1:length(u1d_LSQ)
        if(~nanstat(k))
            line_style = '--';
            if(u1d(u1d_ind,3) < 0); line_style = '-'; end

            r = abs(u1d(u1d_ind,3)) * verticalScale;
            x_data = x(k) + r*cos(theta);
            y_data = y(k) + r*sin(theta);

            plot(x_data, y_data, line_style, 'Color', '#4DBEEE', 'LineWidth', 3, 'HandleVisibility','off');
            u1d_ind = u1d_ind + 1;
        end

        line_style = '--';
        if(u1d_LSQ(k,3) < 0); line_style = '-'; end

        r_lsq = abs(u1d_LSQ(k,3)) * verticalScale;

        x_sim = x(k) + r_lsq*cos(theta);
        y_sim = y(k) + r_lsq*sin(theta);

        plot(x_sim, y_sim, line_style, 'Color', '#A2142F', 'LineWidth', 2.5, 'HandleVisibility','off');
    end

    % Add reference circle for scale
    % Define reference location and scale
    x_ref = 6000;
    y_ref = 6000;
    r_ref = 0.1 * verticalScale;   % Reference vertical displacement magnitude (e.g., 1 unit)
    x_ref_circle = x_ref + r_ref*cos(theta);
    y_ref_circle = y_ref + r_ref*sin(theta);
    
    % Plot the reference circle
    plot(x_ref_circle, y_ref_circle, 'r', 'LineWidth', 2, 'HandleVisibility','off');
    text(x_ref - 270, y_ref - r_ref - 200, "10 cm", 'Color', 'red', 'FontSize', 12, 'HandleVisibility','off');
    quiver(x_ref, y_ref, -radscale * 0.5, 0, 'LineWidth',2.75, 'MaxHeadSize', 0.3, 'Color', 'red', 'HandleVisibility','off'); % Vector scale reference
    text(x_ref - radscale*0.4, y_ref -radscale*0.1, "50 cm", 'Color', 'red', 'FontSize', 12);
    xlabel('Easting (m)', "FontSize", 18);
    ylabel('Northing (m)', "FontSize", 18);
    zlabel('Vertical (m)',  "FontSize", 18);
    title("Displacements at various stations (y oriented N, x oriented E)", "FontSize", 18)
    xlim([-8000, 8000]);
    ylim([-8000, 8000]);
    zlim([-8000, 8000]);
    GPSNameList(end + 1) = "SDH";
    
    cxy = llh2local(coast_new', [-155.2784, 19.4073]);
    cxy = cxy * 1000;
    plot(cxy(1, :)', cxy(2, :)', 'k.', 'HandleVisibility','off');
    
    % text((x-300)', (y+300)', (zeros(size(x)) + 400)', GPSNameList, 'FontSize', 16, 'HandleVisibility','off');
    % plot3(x', y', zeros(size(x))', '.', "Color", "#f76a05", "MarkerSize", 20, 'HandleVisibility','off')
    hold off;
    legend("FontSize", 18, "Location", "southwest");

    if(saveFigs); saveas(6, "Figures/quiver_" + num2str(ntrials, "%.1e") + "trials.fig"); end
end