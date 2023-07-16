function [u, range] = gauss_velocity(x,L,dt)
% calculates smoothed particle velocities from position trajectories using
% differentiated Gaussian kernel (Mordant et al, 2004)
% x - position values (px, mm, etc.)
% L - kernel size (frames)
% fs - sampling frequency (1/s)

% check if track is long enough to smooth
if L > 0
    
    % halfwidth
    T = L/2 - 0.5;
    
    % filter width parameter, value from Mordant et al (2004)
    w = T/1.5; 
    
    % first order derivative of Gaussian:
    t = -T:T;
    dG = -t/(w^3*sqrt(2*pi)) .* exp(-t.^2/(2*w^2));
    % normalize as A*dG+B (see Springer Handbook of Fluid Mech - Tropea p. 795)
    % using the conditions sum(A*dG+B) = 0, sum(t*(A*dG+B)) = 1
    B = 0;
    A = -1/sum((-T:T).*dG);
    dG = A*dG+B;
    
    % Convolve kernel with position values
    u = conv(x, dG, 'valid')/dt;
    range = (T:length(u)+T-1);

else
    % track is not long enough to smooth
    u = gradient(x)/dt;
    range = 1:length(u);
    
end
end