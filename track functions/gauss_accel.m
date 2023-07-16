function [a, range] = gauss_accel(x,L,dt)
% calculates smoothed particle accelerations from position trajectories using
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
    
    % second order derivative of Gaussian:
    t = -T:T;
    ddG = - (w^2-t.^2)/(w^5*sqrt(2*pi)) .* exp(-t.^2/(2*w^2));
    % normalize as A*ddG+B (see Springer Handbook of Fluid Mech - Tropea p. 795)
    % using the conditions sum(A*ddG+B) = 0, sum((-T:T)^2/2*(A*ddg+B)) = 1
    A = 2/(sum(t.^2.*ddG) - 1/numel(t) * sum(ddG) * sum(t.^2));
    B = -A/numel(t) * sum(ddG);
    ddG = A*ddG+B;
    
    % Convolve kernel with position values
    a = conv(x, ddG, 'valid')/dt^2;
    range = (T:length(a)+T-1);

else
    % track not long enough to smooth
    a = gradient(gradient(x))/dt^2;
    range = 1:length(a);

end    
end