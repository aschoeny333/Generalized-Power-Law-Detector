% Check Classes and Values of Input Parameters
% 
% Author: Alex Schoeny, Dr. John Spiesberger
% 
% Goal: Ensure that the classes of input parameters are all valid. Where
% possible, ensure that values of the parameters are also valid
%
% Inputs:
% 
%     fnam - String, .wav file name (e.g. 'test.wav')
% 
%     pass_band - 1 x 2 Double array, Min/max passband frequencies (Hz).
%     Minimum in pass_band(1). For a lowpass filter, set pass_band(1)=0. If
%     pass_band is empty, output data are not filtered.
% 
%     stop_band - 1 x 2 Double array, Min/max stop band (Hz) for bandpass
%     filter. Minimum in stop_band(1). stop_band(1) only used when
%     pass_band(1)>0. then stop_band(1) <pass_band(1). Only used when
%     pass_band is not empty.
% 
%     t_bounds - 1 x 2 Double array, time interval of interest of form [start
%     time, end time] in sec (e.g. [360, 450]). If start time > end time,
%     values are swapped. If any other improper entry (e.g. start time =
%     end time, start time < 0, end time > time series length), t_bounds is
%     set to [0, time series length]
% 
%     gamma - 1 x 1 Double, paramater used to determine A and B (e.g. 1) as
%     used in Helble et al. (2012) equations 7 and 8, p. 2685
% 
%     v1 - 1 x 1 Double, exponential parameter applied to A to determine N
%     (e.g. 1) as defined in Helble et al. (2012) equation 6, p. 2685
% 
%     v2 - 1 x 1 Double, exponential parameter applied to B to determine N
%     (e.g. 2)as defined in Helble et al. (2012) equation 6, p. 2685
% 
%     eta_thresh - 1 x 1 Double, minimum test statistic threshold a signal
%     must reach in at least one time bin in order to be identified by
%     detector.m. Derived largely in sections IIIB-D in Helble et al.
%     (2012) Exact value given in Section IV (2.62 * 10^-4), Helble et al.
%     (2012) p. 2690
% 
%     eta_noise - 1 x 1 Double, minimum test statistic threshold a signal
%     must reach in all time bins in order to have its time interval
%     determined by detector.m; that is, test statistics below this value
%     imply a signal does not exist in that time bin. Derived largely in
%     sections IIIB-D in Helble et al. Exact value given in Section IV
%     (2.07 * 10^-5), Helble et al. (2012) p. 2690
%
%     noise_thresh - 1 x 1 double, value of test statistic above which
%     noise interval iteration procedure stops. Not referenced in either
%     Helble paper, but consider using eta_noise
%
%     max_noise_dur - 1 x 1 double, value in seconds that the length of the
%     noise bounds cannot exceed
%  
%     min_noise_dur - 1 x 1 double, value in seconds that the length of the
%     noise bounds must exceed
%
%     filter_order - 1 x 1 double, input argument for Matlab function
%     designfilt, suggested value of 10
%
%     detector_type - 1 x 1 double. If 1, runs detector.m, which replaces
%     each individual detection with columns containing noise conditions
%     and then recalculates test statistic (hence is better suited for data
%     with many signals, has longer runtime). If any other value, runs
%     detector_simple.m, which does not replace columns of detections
%
% Outputs: t_bounds, stop_band, pass_band, eta_thresh, eta_noise, t_min,
% noise_thresh, max_noise_dur, min_noise_dur, filter_order all have same
% description as above, but may have had value adjusted if original value
% was invalid

function [t_bounds, stop_band, pass_band, eta_thresh, eta_noise, t_min, ...
    noise_thresh, max_noise_dur, min_noise_dur, filter_order] = check_inputs(fnam, ...
    wav_dir, pass_band, stop_band, t_bounds, gamma, v1, v2, eta_thresh, ...
    eta_noise, t_min, noise_thresh, max_noise_dur, min_noise_dur, filter_order, detector_type)
    
% Step 0.1: Check class validity of all parameters
    if ~(isa(fnam, 'string') || isa(fnam, 'char'))
        error('Incorrect class of fnam. Revise and re-run')
    end
    
    if ~(isa(wav_dir, 'string') || isa(wav_dir, 'char'))
        error('Incorrect class of wav_dir. Revise and re-run')
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
    
    if ~isa(filter_order, 'double')
        error('Incorrect class of filter_order. Revise and re-run')
    end
    
    if ~isa(detector_type, 'double')
        error('Incorrect class of detector_type. Revise and re-run')
    end 
    
% Step 0.2: Read in .wav file
    [data, samp_rate] = audioread(strcat(wav_dir, fnam));

% Step 0.3: Check validity of t_boundsd
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

% Step 0.4: Check validity of stop_band
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

% Step 0.5: Check validity of pass_band
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

% Step 0.6: Check validity of pass_band compared to stop_band
    if stop_band(1) > pass_band(1)
        stop_band(1) = pass_band(1);
    end
    if stop_band(2) < pass_band(2)
        stop_band(2) = pass_band(2);
    end
    
% Step 0.7: Check that several paramaters are non-negative. If so, change
% value to that used in Helble et al. (2012, 2015) or standard value
% suggested by Dr. Spiesberger
    
    if eta_thresh < 0
        eta_thresh = 2.62 * 10^-4;
        disp('eta_thresh is negative, value changed by check_inputs.m to 2.62 * 10^-4');
    end
    
    if eta_noise < 0
        eta_noise = 2.07 * 10^-5;
        disp('eta_noise is negative, value changed by check_inputs.m to 2.07 * 10^-5');
    end
    
    if t_min < 0
        t_min = 0.35;
        disp('t_min is negative, value changed by check_inputs.m to 0.35');
    end
    
    if noise_thresh < 0
        noise_thresh = eta_thresh / 2;
        disp('noise_thresh is negative, value changed by check_inputs.m to eta_thresh / 2');
    end
    
    if max_noise_dur < 0
        max_noise_dur = 5;
        disp('max_noise_dur is negative, value changed by check_inputs.m to 5');
    end
    
    if min_noise_dur < 0
        min_noise_dur = 1;
        disp('min_noise_dur is negative, value changed by check_inputs.m to 1');
    end
    
    if filter_order < 0
        filter_order = 10;
        disp('filter_order is negative, value changed by check_inputs.m to 10');
    end
    
% Step 0.8: Ensure that filter_order is an integer

filter_order = round(filter_order);
    
    
    
