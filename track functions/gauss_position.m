function [x2, range] = gauss_position(x,L)
% smooths particle trajectories using Gaussian kernel following Mordant et al (2004)
% x - position values (px, mm, etc.)
% L - kernel size (frames). must be odd. enter L=0 for no smoothing.

% check if track is long enough to smooth
if L > 0
    
    % halfwidth
    T = L/2 - 0.5;

    % filter width parameter, value from Mordant et al (2004)
    w = T/1.5; 

    % Construct the Gaussian kernel
    t = -T:T;
    G = 1/(w*sqrt(2*pi))*exp(-t.^2/(2*w^2));
    G = G / sum(G); % Normalize area to 1.

    % Convolve kernel with position values
    x2 = conv(x, G, 'valid');
    range = (T+1:length(x2)+T); 

else
    % track not long enough to smooth
    x2 = x;
    range = 1:length(x2);

end
end