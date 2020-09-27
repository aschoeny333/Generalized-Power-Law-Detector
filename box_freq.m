function [bin_intervals] = box_freq(X_masked, sig_intervals)
    rows = length(X_masked(:, 1));
    bin_intervals = zeros(size(sig_intervals));
    if ~isempty(sig_intervals)
        for i = 1:length(sig_intervals(1, :))
            min_freq = rows + 1;
            max_freq = 0;
            X_sig = X_masked(:, sig_intervals(1, i):sig_intervals(2, i));
            for r = 1:rows
                for c = 1:length(X_sig(1, :))
                    if X_sig(r,c) > 0
                        if r < min_freq
                            min_freq = r;
                        end
                        if r > max_freq
                            max_freq = r;
                        end
                    end
                end
            end
            bin_intervals(:, i) = [min_freq; max_freq];
        end
    end
end