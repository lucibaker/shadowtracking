function uid = find_track(tracks, cols, vals)
% find track that contains the values vals in the columns cols.
% for example, to find the track containing [xp, yp] = [0.535, -.219]:
% uid = find_track(tracks, [1, 2], [0.535, -.219])
% 
% inputs 
% tracks: array containing track info
% cols: columns of tracks array to search in 
% vals: values to find
%
% outputs
% uid: track id containing the value

% tolerance 
tol = 1e-6; 

% find entries matching values within tolerance
idx = true(size(tracks,1),1);
for i = 1:length(cols)
    idx = idx & abs( tracks(:,cols(i)) - vals(i) ) < tol;
end

% track id
uid = unique(tracks(idx, 5));

end