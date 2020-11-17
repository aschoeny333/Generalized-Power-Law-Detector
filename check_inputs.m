function [t_bounds, stop_band, pass_band, eta_thresh, eta_noise, t_min, ...
    noise_thresh, max_noise_dur, min_noise_dur] = check_inputs(fnam, wav_dir, programs_dir, ... 
    pass_band, stop_band, t_bounds, gamma, v1, v2, eta_thresh, eta_noise, t_min, ...
    noise_thresh, max_noise_dur, min_noise_dur)
    %     Step 0.1: Check class validity of all parameters
    
    if ~(isa(fnam, 'string') || isa(fnam, 'char'))
        error('Incorrect class of fnam. Revise and re-run')
    end
    
    if ~(isa(wav_dir, 'string') || isa(wav_dir, 'char'))
        error('Incorrect class of wav_dir. Revise and re-run')
    end
    
    if ~(isa(programs_dir, 'string') || isa(programs_dir, 'char'))
        error('Incorrect class of programs_dir. Revise and re-run')
    end
    
    if ~isa(t_bounds, 'double')
        error('Incorrect class of t_bounds. Revise and re-run')
    end
    
    if ~isa(stop_band, 'double')
        error('Incorrect class of stop_band. Revise and re-run')
    end
    
    if ~isa(pass_band, 'double')
        error('Incorrect class of pass_band. Revise and re-run')
    end
    
    if ~isa(gamma, 'double')
        error('Incorrect class of gamma. Revise and re-run')
    end
    
    if ~isa(v1, 'double')
        error('Incorrect class of v1. Revise and re-run')
    end
    
    if ~isa(v2, 'double')
        error('Incorrect class of v2. Revise and re-run')
    end
    
    if ~isa(eta_thresh, 'double')
        error('Incorrect class of eta_thresh. Revise and re-run')
    end
    
    if ~isa(eta_noise, 'double')
        error('Incorrect class of eta_noise. Revise and re-run')
    end
    
    if ~isa(noise_thresh, 'double')
        error('Incorrect class of noise_thresh. Revise and re-run')
    end
    
    if ~isa(t_min, 'double')
        error('Incorrect class of t_min. Revise and re-run')
    end
    
    if ~isa(min_noise_dur, 'double')
        error('Incorrect class of min_noise_dur. Revise and re-run')
    end
    
    if ~isa(max_noise_dur, 'double')
        error('Incorrect class of max_noise_dur. Revise and re-run')
    end    
    
    %     Step 0.2: Read in .wav file
    cd(wav_dir);
    [data, samp_rate] = audioread(fnam);

    %     Step 0.3: Check validity of t_boundsd
    if t_bounds(1) > t_bounds(2)
        t_bounds = [t_bounds(2), t_bounds(1)];
    end
    if t_bounds(1) == t_bounds(2) || ~isa(t_bounds, 'double')
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
    if stop_band(1) == stop_band(2) || ~isa(stop_band, 'double')
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
    if pass_band(1) == pass_band(2) || ~isa(pass_band, 'double')
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
    
    %     Step 0.7: Check that several paramaters are non-negative. If so,
    %     change value to that used in Helble et al. (2012, 2015) or
    %     standard value suggested by Dr. Spiesberger
    
    if eta_thresh < 0
        eta_thresh = 2.62 * 10^-4;
    end
    
    if eta_noise < 0
        eta_noise = 2.07 * 10^-5;
    end
    
    if t_min < 0
        t_min = 0.35;
    end
    
    if noise_thresh < 0
        noise_thresh = eta_thresh / 2;
    end
    
    if max_noise_dur < 0
        max_noise_dur = 5;
    end
    
    if min_noise_dur < 0
        min_noise_dur = 1;
    end
    
    
    
