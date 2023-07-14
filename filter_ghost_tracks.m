function [tracks, tracklength] = filter_ghost_tracks(tracks, tracklength, edges, tol)
% remove tracks that start or end at a frame edge ("ghost tracks")
%
% inputs
% smtracks
% smtracklength
% edges: vector containing x-coord of frame edges [m]
% tol: tolerance on particle position wrt frame edge position
%
% outputs
% smtracks, smtracklength: with tracks that start or end at a frame edge removed

Ntracks = length(tracklength);
is_ghost = false(Ntracks,1);

for i = 1:Ntracks
    idx = find(tracks(:,5) == i);
    i_start = idx(1);
    i_end = idx(end);
    % check if the track start or end x-coord is within tol of any frame edge
    is_ghost(i) = any( abs(tracks(i_start,1) - edges) < tol | abs(tracks(i_end,1) - edges) < tol );
end

% remove ghost tracks
[tracks,tracklength] = remove_track(tracks,tracklength,find(is_ghost));

end
        
