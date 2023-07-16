function [tracks,tracklength] = remove_track(tracks, tracklength, track_id)
% remove the track corresponding to track_id, i.e. tracks(:,5)==track_id

for i = length(track_id):-1:1
    % indices of track to remove
    idx = tracks(:,5)==track_id(i);
    
    % remove track
    tracks(idx,:) = [];
    
    % update subsequent track_ids so they are all consecutive numbers
    tracks(tracks(:,5) > track_id(i), 5) = tracks(tracks(:,5) > track_id(i), 5) - 1;
    
    % remove track length entry
    tracklength(track_id(i)) = [];
end