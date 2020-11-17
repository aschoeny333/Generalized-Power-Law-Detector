% Testing the whitener function:
% The objective of whitener.m is to determine a value mu_k for each 
% k^th row of a spectrogram matrix
% 
% To test the efficacy of whitener.m, a row vector will be assigned
% a known distribution, and its distribution will be compared before and after
% subtracting the output of whitener.m. 
% 
% Expectation is to see a 
% 
% (Typing this and referring to Helble et al I am beginning to realize that 
% some understanding of the statistics section is necessary to properly
% test the performance of GPL.m and related programs)

X_test = zeros(rows, cols);
for r = 1:rows
    real = randn(cols, 1);
    imag = randn(cols, 1);
    Xk = complex(real, imag) / sqrt(2);
end



ray_x = linspace(0, 10, 250);
rayleigh = 4*ray_x/pi .* exp(-2*ray_x .^ 2 / pi);

figure; 
histogram(abs(Xk), 200, 'FaceColor', 'b', 'Normalization', 'pdf'); hold;
plot(ray_x, rayleigh); 
title('Comparing pdf of random variable X_k to predicted Rayleigh distribution with mean \surd\pi/2');

[mu, j_star, rows, cols] = whitener(abs(X_test));
Z = (abs(X_test) - mu) .^ 6;


