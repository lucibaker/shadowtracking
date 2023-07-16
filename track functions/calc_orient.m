function pxyz = calc_orient(th, d, Dp, particle_type)
% calculate 3d particle orientation from 2d orientation and projected axis length
% th: 2d orientation
% d: projected axis length
% Dp: nominal axis length
% particle type: 'r' or 'd' for rod or disk
% pxyz: particle orientation vector

if strncmp(particle_type,'r',1)
    pxyz = [d/Dp.*cos(th), ... % p_x (streamwise)
        d/Dp.*sin(th), ... % p_y (vertical)
        sqrt(1 - (d/Dp).^2)]; % p_z (spanwise)

elseif strncmp(particle_type,'d',1)
    pxyz = [sin(th).*sqrt(1 - (d/Dp).^2), ... % p_x (streamwise)
        cos(th).*sqrt(1 - (d/Dp).^2).*-sign(th), ... % p_y (vertical)
        d/Dp]; % p_z (spanwise)
else
    error('invalid particle type')
end