function check_inputs(fnam, pass_band, stop_band, t_bounds, gamma, v1, v2, thresh)
    %     Step 0.1: Check class validity of all parameters
    if ~isa(fnam, "string") || ~isa(t_bounds, "double") || ~isa(stop_band, "double") ...
            || ~isa(pass_band, "double") || ~isa(gamma, "double") || ...
            ~isa(v1, "double") || ~isa(v2, "double") || ~isa(thresh, "double")
        error("Incorrect class of a parameter. Revise and re-run")
    end
    %     Step 0.2: Read in .wav file
    [data, samp_rate] = audioread(fnam);
    sound.d = data;
    sound.s = samp_rate;

    %     Step 0.3: Check validity of t_bounds
    if t_bounds(1) > t_bounds(2)
        t_bounds = [t_bounds(2), t_bounds(1)];
    end
    if t_bounds(1) == t_bounds(2) || ~isa(t_bounds, "double")
        t_bounds = [0, length(data) / samp_rate];
    end
    if t_bounds(1) < 0
        t_bounds(1) = 0;
    end
    if t_bounds(2) > length(data) / samp_rate
        t_bounds(2) = length(data) / samp_rate;
    end

    %     Step 0.4: Check validity of stop_band
    if stop_band(1) > stop_band(2)
        stop_band = [stop_band(2), stop_band(1)];
    end
    if stop_band(1) == stop_band(2) || ~isa(stop_band, "double")
        stop_band = [0, samp_rate / 2];
    end
    if stop_band(1) < 0
        stop_band(1) = 0;
    end
    if stop_band(2) > samp_rate / 2
        stop_band(2) = samp_rate / 2;
    end

    %     Step 0.5: Check validity of pass_band
    if pass_band(1) > pass_band(2)
        pass_band = [pass_band(2), pass_band(1)];
    end
    if pass_band(1) == pass_band(2) || ~isa(pass_band, "double")
        pass_band = [0, samp_rate / 2];
    end
    if pass_band(1) < 0
        pass_band(1) = 0;
    end
    if pass_band(2) > samp_rate / 2
        pass_band(2) = samp_rate / 2;
    end

    %     Step 0.6: Check validity of pass_band compared to stop_band
    if stop_band(1) > pass_band(1)
        stop_band(1) = pass_band(1);
    end
    if stop_band(2) < pass_band(2)
        stop_band(2) = pass_band(2);
    end
