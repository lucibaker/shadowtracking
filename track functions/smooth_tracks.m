function [smtracks, smtracklength, avar_k] = smooth_tracks(tracks0,kernel,dt)
% apply gaussian smoothing to tracks and get position, velocity,
% acceleration. if length(kernel)>1, returns acceleration variance vs
% kernel size.


avar_k = zeros(size(kernel));
Nsmtracks = zeros(size(kernel));

for k = 1:length(kernel)
    fprintf(['kernel = ' num2str(kernel(k)) '...'])
    mintracks = kernel(k)*2+3;  % minimum length of tracks for smoothing (at least three points)
    j = 0;
    smtracks = zeros(0,size(tracks0,2));
    smtracklength = [];

    for i = 1:max(tracks0(:,5))
        idx_i = find(tracks0(:,5)==i);

        if any(idx_i) && length(tracks0(idx_i,1))>=mintracks 
            % sort track by lifetime counter (in case it is out of order after fixing a broken track)
            [~,idx_sort] = sort(tracks0(idx_i,6),'ascend'); 
            idx_i = idx_i(idx_sort);
            
            xp = tracks0(idx_i,1);
            yp = tracks0(idx_i,2);
            [x2, range] = gauss_position(xp',kernel(k));
            [y2, ~] = gauss_position(yp',kernel(k));
            [u2, ~] = gauss_velocity(xp',kernel(k),dt);  
            [v2, ~] = gauss_velocity(yp',kernel(k),dt);

            j = j + 1; 
            UID2 = ones(1,length(x2))*j;
            if ~isempty(x2)
                [Udot_ms2, ~] = gauss_accel(xp',kernel(k),dt);
                [Vdot_ms2, ~] = gauss_accel(yp',kernel(k),dt);
                smtracks = [smtracks; x2',y2',u2',v2',UID2',range', ...
                range'+tracks0(idx_i(1),7)',Udot_ms2',Vdot_ms2',tracks0(idx_i(range),10:end)]; 
                smtracklength = [smtracklength; length(x2)];
            end
        end
    end
    
    avar_k(k) = std(sqrt(smtracks(:,8).^2 + smtracks(:,9).^2),'omitnan')^2;
    Nsmtracks(k) = length(unique(smtracks(:,5)));
    fprintf([num2str(Nsmtracks(k)) ' tracks\n'])
end

end