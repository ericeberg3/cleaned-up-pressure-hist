function m = get_full_m(prior, opt_params, forward)
    % Forward = true means that we want to get a full 1x16 m state vector
    % based on our incomplete guess with 6 parameters
    if(forward)
        taiyi_parameters = prior;
        % ["HMM volume", "dpHMM", "SC Volume", "SC aspect ratio", "dip", "dpSC"];
        aspect_ratio_HMM = 1.7496;
        opt_vert_sd = (3/(4*pi) * opt_params(1) * (aspect_ratio_HMM^2))^(1/3);
        opt_horiz_sd = opt_vert_sd/(aspect_ratio_HMM);
        
        mHMM = [opt_vert_sd, opt_horiz_sd, prior(3:6), prior(7) - abs(opt_vert_sd - prior(1)), opt_params(2)];
        
        aspect_ratio_SC = opt_params(4);
        opt_vert_sd = (3/(4*pi) * opt_params(3) * (aspect_ratio_SC^2))^(1/3);
        opt_horiz_sd = opt_vert_sd/(aspect_ratio_SC);
        prior_SC = prior(9:end);
        mSC = [opt_vert_sd, opt_horiz_sd, opt_params(5), prior_SC(4:7), opt_params(6)];
        
        m = [mHMM, mSC];
    % Otherwise we want to get the smaller m guess vector from a complete m
    % vector
    else
        m = [(4/3 * pi * prior(1) * prior(2)^2), prior(8), ...
            (4/3 * pi * prior(9) * prior(10)^2), prior(9)/prior(10), prior(11), prior(end)];
    end

end