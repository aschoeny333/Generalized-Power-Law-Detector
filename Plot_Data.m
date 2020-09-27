% Generalized Power Law Algorithm plotting with test statistics 
%  
% Author: Alex Schoeny
% 
% Goal: Produce spectrograms of each major step of GPL implementation
% alongside corresponding test statistic plots 
% 
% References:
% Helble, Tyler A et al. ?A generalized power-law detection algorithm for
% humpback whale vocalizations.? The Journal of the Acoustical Society of
% America vol. 131,4 (2012): 2682-99. doi:10.1121/1.3685790
% 
% Inputs: X, mu, N, N_thresh from GPL.m
%
%     A note on size definitions: Some matrices and arrays sizes are
%     defined using m, n, and m'. The values of these will change based on
%     the values of the inputs of GPL.m, especially pass_band, stop_band,
%     and t_bounds. Note also that m' <= m. Matrices all have rows with
%     constant frequency and columns with constant time.
% 
%     X - Variable name: X. Double matrix of size m x n as described in the
%     note above; absolute value of fourier_trimmed. Not identical to the X
%     defined in Helble et al (2012) p.2684, but the X in Helble et al
%     (2012) is only ever used in a context where its absolute value is
%     taken, hence the shorthand used here
% 
%     mu - Variable name: mu. Double column vector of size m as described
%     in note above; contains each value of mu_k (eq 9 in Helble et al.
%     (2012)) p.2685
%
%     N - Double matrix of size m x n as described in note above; defined
%     by eq 6 in Helble et al (2012) p.2685
%     
%     N_thresh - Double matrix of size m' x n as described in note above;
%     cells are either identical to N or equal to thresh if the
%     corresponding value in N was less than thresh
%
%     eta_thresh - 1 x 1 Double, minimum test statistic threshold a signal
%     must reach in at least one time bin in order to be identified by
%     detector.m. Derived largely in sections IIIB-D in Helble et al. Exact
%     value given in Section IV (2.62 * 10^-4) p.2690. Used here to draw
%     threshold lines on test statistic plots
% 
%     eta_noise - 1 x 1 Double, minimum test statistic threshold a signal
%     must reach in all time bins in order to have its time interval
%     determined by detector.m; that is, test statistics below this value
%     imply a signal does not exist in that time bin. Derived largely in
%     sections IIIB-D in Helble et al. Exact value given in Section IV
%     (2.07 * 10^-5) p.2690. Used here to draw threshold lines on test
%     statistic plots
% 
% Outputs
% 
%     X_sum - Double row vector of size n as described in note above; test
%     statistic for original spectrogram defined by eq 2 in Helble et al
%     (2012) p.2684, where v = 1
%     
%     Xw_sum - Double row vector of size n as described in X note above;
%     test statistic for whitened spectrogram defined by eq 2 in Helble et
%     al (2012) p.2684, where v = 1 and X is replaced with X - mu
%     
%     N_sum - Double row vector of size n as described in X note above;
%     test statistic for spectrogram after GPL defined by eq 6 in Helble et
%     al (2012) p.2685
% 
%     Nt_sum - Double row vector of size n as described in X note above;
%     test statistic for spectrogram after GPL and thresholding defined by
%     eq 6 in Helble et al (2012) p.2685

% Plotting original spectrogram

function [X_sum, Xw_sum, N_sum, Nt_sum] = Plot_Data(original, mu, N, N_thresh, ...
    t_bounds, pass_band, stop_band, gamma, v1, v2, eta_thresh, eta_noise, ...
    time_intervals, X_masked, freq_intervals)

    % Unpack structured arrays
    X = original.X;
    freqs_trimmed = original.f;
    times = original.t;

    % First plot: Original vs Modified (N) spectrograms
    %     First subplot: Original spectrogram
    subplot(2,1,1);
    surf(times + t_bounds(1), freqs_trimmed, log10(X), 'EdgeColor', 'none'); % If power argument is changed, change title
    axis xy; 
    axis tight; 
    colormap(jet); view(0,90);
    xlabel('Time (secs)');
    c = colorbar;
    ylabel('Frequency(HZ)');
    ylabel(c, "Power");
    title("Original Spectrogram: pass_band = [" + num2str(pass_band(1)) + ... 
        ", " + num2str(pass_band(2)) + "], stop_band = [" + num2str(stop_band(1)) + ...
        ", " + num2str(stop_band(2)) + "], (\gamma, v_1, v_2) = (" + ...
        num2str(gamma) + ", " + num2str(v1) + ", " + num2str(v2) + ...
        "), power as log10(X)");
    
    %    Second subplot: Modified (N) Spectrogram
        subplot(2,1,2);
    surf(times + t_bounds(1), freqs_trimmed, log10(N), 'EdgeColor', 'none'); % If power argument is changed, change title 
    axis xy; 
    axis tight; 
    colormap(jet); view(0,90);
    xlabel('Time (secs)');
    c = colorbar;
    ylabel('Frequency(HZ)');
    ylabel(c, "Power");
    title("Spectrogram After GPL Applied: pass_band = [" + num2str(pass_band(1)) + ... 
        ", " + num2str(pass_band(2)) + "], stop_band = [" + num2str(stop_band(1)) + ...
        ", " + num2str(stop_band(2)) + "], (\gamma, v_1, v_2) = (" + ...
        num2str(gamma) + ", " + num2str(v1) + ", " + num2str(v2) + ...
        "), power as log10(X)");
    
    % Plot developments of spectrograms in the style of Helble et al.
    % (2012) figure 5
    
    figure;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subplot(4,2,1);
    surf(times + t_bounds(1), freqs_trimmed, log10(X), 'EdgeColor', 'none'); 
    axis xy; 
    axis tight; 
    colormap(jet); view(0,90);
    xlabel('Time (secs)');
    c = colorbar;
    ylabel('Frequency(HZ)');
    ylabel(c, "Power");
    title('Original Spectrogram');

    % Test statistic for original spectrogram
    subplot(4,2,2);
    X_sum = sum(X .^ 2);
    plot(times + t_bounds(1), X_sum);
    yline(eta_thresh, "b");
    yline(eta_noise, "k");
    xlabel("(Time (secs)");
    ylabel("\Sigma |X|^2");

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Plotting conditionally whitened spectrogram (Original minus mu)
    subplot(4,2,3);
    surf(times + t_bounds(1), freqs_trimmed, log10(abs(X - mu)), 'EdgeColor', 'none'); 
    axis xy; 
    axis tight; 
    colormap(jet); view(0,90);
    xlabel('Time (secs)');
    c = colorbar;
    ylabel('Frequency(HZ)');
    ylabel(c, "Power");
    title('Spectrogram After Conditional Whitening');

    % Test statistic for original spectrogram
    subplot(4,2,4);
    Xw_sum = sum((X-mu) .^ 2);
    plot(times + t_bounds(1), Xw_sum);
    yline(eta_thresh, "b");
    yline(eta_noise, "k");
    xlabel("(Time (secs)");
    ylabel("\Sigma |X_w|^2");

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Plotting spectrogram after GPL
    subplot(4,2,5);
    surf(times + t_bounds(1), freqs_trimmed, log10(N), 'EdgeColor', 'none');
    axis xy; 
    axis tight; 
    colormap(jet); view(0,90);
    xlabel('Time (secs)');
    c = colorbar;
    ylabel('Frequency(HZ)');
    ylabel(c, "Power");
    title('Spectrogram After GPL Applied');

    % Plotting T^g, the test statistic defined in Helble et al
    subplot(4,2,6);
    N_sum = sum(N);
    plot(times + t_bounds(1), N_sum * 1000);
    yline(eta_thresh * 1000, "b");
    yline(eta_noise * 1000, "k");
    xlabel("(Time (secs)");
    ylabel("T^g(X) \times 10^3");

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Plotting spectrogram after GPL and thresholding
    subplot(4,2,7);
    surf(times + t_bounds(1), freqs_trimmed, log10(N_thresh), 'EdgeColor', 'none'); 
    axis xy; 
    axis tight; 
    colormap(jet); view(0,90);
    xlabel('Time (secs)');
    c = colorbar;
    ylabel('Frequency(HZ)');
    ylabel(c, "Power");
    title("Spectrogram After GPL Applied and Thresholding");

    % Plotting T^g, the test statistic defined in Helble et al
    subplot(4,2,8);
    Nt_sum = sum(N_thresh);
    plot(times + t_bounds(1), N_sum * 1000);
    yline(eta_thresh * 1000, "b");
    yline(eta_noise * 1000, "k");
    xlabel("(Time (secs)");
    ylabel("T^g(X)' \times 10^3");
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Plotting spectrogram and T^g (row 3 above) with time bounds displayed
    
    figure;
    
    subplot(2,1,1);
    surf(times + t_bounds(1), freqs_trimmed, log10(N), 'EdgeColor', 'none');
    if numel(time_intervals) == 0
        % Don't plot any time bound lines
    elseif numel(time_intervals) == 2
        xline(time_intervals(1, 1), "color", "w", "LineWidth", 5);
        xline(time_intervals(2, 1), "color", "w", "LineWidth", 5);
    else
        for i = 1 : length(time_intervals(1, :))
            xline(time_intervals(1, i), "color", "w", "LineWidth", 5);
            xline(time_intervals(2, i), "color", "w", "LineWidth", 5);
        end
    end
    axis xy; 
    axis tight; 
    colormap(jet); view(0,90);
    xlabel('Time (secs)');
    c = colorbar;
    ylabel('Frequency(HZ)');
    ylabel(c, "Power");
    title('Spectrogram After GPL Applied With Bounds Determined By Detector');

    % Plotting T^g, the test statistic defined in Helble et al
    
    subplot(2,1,2);
    N_sum = sum(N);
    plot(times + t_bounds(1), N_sum * 1000);
    if numel(time_intervals) == 0
        % Don't plot any time bound lines
    elseif numel(time_intervals) == 2
        xline(time_intervals(1, 1), "color", "g", "LineWidth", 5);
        xline(time_intervals(2, 1), "color", "r", "LineWidth", 5);
    else
        for i = 1 : length(time_intervals(1, :))
            xline(time_intervals(1, i), "color", "g", "LineWidth", 5);
            xline(time_intervals(2, i), "color", "r", "LineWidth", 5);
        end
    end
    yline(eta_thresh * 1000, "b");
    yline(eta_noise * 1000, "k");
    xlabel("(Time (secs)");
    ylabel("T^g(X) \times 10^3");
    title("Test Statistic vs. Time With Bounds Determined By Detector");
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Plotting original vs. masked spectrogram

    figure;
    
    subplot(2,1,1)
    surf(times + t_bounds(1), freqs_trimmed, X, 'EdgeColor', 'none');
    axis xy; 
    axis tight; 
    colormap(jet); view(0,90);
    if numel(time_intervals) == 0
        % Don't plot any time bound lines
    elseif numel(time_intervals) == 2
        xline(time_intervals(1, 1), "color", "g", "LineWidth", 2);
        xline(time_intervals(2, 1), "color", "r", "LineWidth", 2);
        yline(freq_intervals(1, 1), "color", "g", "LineWidth", 2);
        yline(freq_intervals(2, 1), "color", "r", "LineWidth", 2);
%         rectangle("Position", [time_intervals(1,1), freq_intervals(1,1), ...
%             time_intervals(2,1) - time_intervals(1,1), freq_intervals(2,1) - freq_intervals(1,1)], ...
%             "EdgeColor", "w", "LineWidth", 3); 
    else
        for i = 1 : length(time_intervals(1, :))
            xline(time_intervals(1, i), "color", "g", "LineWidth", 2);
            xline(time_intervals(2, i), "color", "r", "LineWidth", 2);
            yline(freq_intervals(1, i), "color", "g", "LineWidth", 2);
            yline(freq_intervals(2, i), "color", "r", "LineWidth", 2);
%             rectangle("Position", [time_intervals(1,i), freq_intervals(1,i), ...
%             time_intervals(2,i) - time_intervals(1,i), freq_intervals(2,i) - freq_intervals(1,i)], ...
%             "EdgeColor", "k", "LineWidth", 3); 
        end
    end
    xlabel('Time (secs)');
    c = colorbar;
    ylabel('Frequency(HZ)');
    ylabel(c, "Power");
    title('Spectrogram Before any Manipulation, with Detection Windows');

    subplot(2,1,2)
    surf(times + t_bounds(1), freqs_trimmed, X_masked, 'EdgeColor', 'none');
    axis xy; 
    axis tight; 
    colormap(jet); view(0,90);
    if numel(time_intervals) == 0
        % Don't plot any time bound lines
    elseif numel(time_intervals) == 2
        rectangle("Position", [time_intervals(1,1), freq_intervals(1,1), ...
            time_intervals(2,1) - time_intervals(1,1), freq_intervals(2,1) - freq_intervals(1,1)], ...
            "EdgeColor", "w", "LineWidth", 3); 
    else
        for i = 1 : length(time_intervals(1, :))
            rectangle("Position", [time_intervals(1,i), freq_intervals(1,i), ...
            time_intervals(2,i) - time_intervals(1,i), freq_intervals(2,i) - freq_intervals(1,i)], ...
            "EdgeColor", "w", "LineWidth", 3); 
        end
    end
    xlabel('Time (secs)');
    c = colorbar;
    ylabel('Frequency(HZ)');
    ylabel(c, "Power");
    title('Spectrogram After GPL, Detector, and Masking Procedure');

    
    
    
    
    
    
    
    
    
    
    
    