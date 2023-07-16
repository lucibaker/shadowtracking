function [tracks,tracklength] = get_particle_tracks(centers,fs,searchrad,varargin)
% Track particles using kNN (k nearest neighbors) search. Previous particle
% displacement used as predictor for next timestep.
%
% Inputs:
% centers: cell array containing the coordinates of the particle centroids,
%   one cell per frame [px]
% kernel: smoothing kernel size
% fs: sampling frequency of particle centroids
% searchrad: search radius [m]
% (optional) angles, errchk:  cell array containing the orientation and
% error check info of the particles, one cell per frame
%
% Outputs:
% smtracks: smoothed particle tracks in the format [X_m Y_m U_ms V_ms UID lifetime frameno Ax_ms2 Ay_ms2 (th_p) (d_p) (errchk)]
% smtracklength: length of each track (# frames)

img_nt = length(centers);

snapshots = cell(img_nt,1); 
for i = 1:img_nt 
    snapshots{i} = [centers{i} zeros(size(centers{i},1),4)];
    for j = 1:nargin - 3 
        snapshots{i} = [snapshots{i} varargin{j}{i}];
    end
end

% get tracks
tracks = trackp(snapshots,searchrad,fs);

% convert to array
filled = tracks(~cellfun('isempty',tracks))';
if isempty(filled)
    error('no tracks found')
end
for i=1:length(filled)
    UID(i,1) = mode(filled{i}(:,5));
    UID(i,2) = length(filled{i}(:,5));
end
UID = sortrows(UID,2);
UID = flipud(UID);

tracklength = zeros(length(tracks),1);
for i = 1:length(tracks)
    tracklength(i) = size(tracks{i},1);
end   
filled_tracklength = tracklength(~cellfun('isempty',tracks));

filled_array = zeros(0,size(filled{1},2));
for i = 1:length(filled)
    filled_array = [filled_array; filled{i}];
end

tracks = filled_array;
tracklength = filled_tracklength;

end


function tracks = trackp(snapshots,searchrad,fs)
% Creates Lagrangian tracks from particle position and velocity data.
% Author : Kee Onn Fong (fongx065@umn.edu)
% edited by Luci Baker (bake0616@umn.edu)

%%%%%%%%%%%%%%%%%%%%%%%
fprintf('finding particle tracks...')
% tstart = tic;

nframes = numel(snapshots);
dt = 1/fs; % T is in s
%%%%%%%%%%%%%%%%%%%%%%%
%% Particle UID Assigning using knnsearch
% [X_m Y_m U_ms V_ms UID lifetime];
% tic
fA = snapshots{1};
fA(:,5) = (1:numel(fA(:,1)))'; %Unique ID assigned for first frame
fA(:,6) = 0; % lifetime of particle
snapshots{1} = fA;
maxmaxUID = 0; % overall max UID

for i = 1:nframes-1
    if isempty(snapshots{i}) || isempty(snapshots{i+1})
        continue;
    end
    fA = snapshots{i};
    fB = snapshots{i+1};
    flag = 0;
    if isempty(fB)
        disp(['frame #',num2str(i),' is empty'])
        continue;
    end
    if isempty(fA)
        count = i;
        flag = 1;
        while isempty(fA)
            fA = snapshots{count};
            count = count-1;
        end
    end
    fA = sortrows(fA,1);
    fB = sortrows(fB,1);
    [np,~]=size(fB);
    maxUIDfA = max([max(fA(:,5)),maxmaxUID]);
    maxmaxUID = maxUIDfA;
    minUIDfA = min(fA(:,5));
    fB(:,5) = ones(numel(fB(:,1)),1)*-1; %-1 is no ID assigned yet
    fB(:,6) = -1; % lifetime of particle 
    
    % predicted displacement
    pred_disp = fA(:,3:4);  

    % use k-nearest neighbor search to associate particles in fA and fB
    [idx,d] = knnsearch(fA(:,1:2) + pred_disp, fB(:,1:2));
    
    idx = [idx,(1:np)',d];
    idx = sortrows(idx,3); idxuncut = numel(idx(:,1));
    idx = idx((idx(:,3)<searchrad),:);
    idx = sortrows(idx,1);

    % if a particle in frame A matches with more than 1 particle in frame B
    rowdel = [];
    for k = 1:numel(idx(:,1))
        rowno = idx(k,1);
        if sum(idx(:,1)==rowno)>1
%             disp(['duplicate found at ',num2str(rowno),' in frame ',num2str(i)]);
            rownos = find(idx(:,1)==rowno);
            rowdel = vertcat(rowdel,rownos(2:end)); % keep only the closest match
        end
    end

    rowdel = unique(rowdel);cutted = 0;
    for k=1:numel(rowdel)
        idx(rowdel(k)-cutted,:)=[];
        cutted = cutted+1;
    end
    
    % check if particle fields are time-resolved
    if flag == 1
%         assign new UID for particles in fB and new lifetime count
        maxUIDfB = maxUIDfA;
        for j=1:numel(fB(:,1))
            if fB(j,5)==-1
                maxUIDfB = maxUIDfB+1;
                fB(j,5) = maxUIDfB;%+1;
                fB(j,6) = 0;
            end
        end
    else
        %associate particles in fA with itself in fB and add lifetime count
        for j=1:numel(idx(:,1))          
            UID = fA(idx(j,1),5);
            fB(idx(j,2),5) = UID;
            fB(idx(j,2),6) = fA(idx(j,1),6)+1;
            % update displacement
            fB(idx(j,2),3) = fB(idx(j,2),1) - fA(idx(j,1),1);
            fB(idx(j,2),4) = fB(idx(j,2),2) - fA(idx(j,1),2);
        end
        %assign new UID for particles in fB and new lifetime count
        maxUIDfB = maxUIDfA;
        for j=1:numel(fB(:,1))
            if fB(j,5)==-1
                maxUIDfB = maxUIDfB +1;
                fB(j,5) = maxUIDfB;%+1;
                fB(j,6) = 0;
            end
        end
    end        
    
    snapshots{i+1} = fB;
end


% Particle Velocity
% calculate velocity for connected particles

snapminUID = zeros(nframes,1);
snapmaxUID = zeros(nframes,1);
for i=1:nframes-1
    if isempty(snapshots{i}) || isempty(snapshots{i+1})
        snapminUID(i) = max(snapminUID);
        snapmaxUID(i) = max(snapmaxUID);
        continue;
    end
    fA = snapshots{i};
    fB = snapshots{i+1};
    for fbx=1:numel(fB(:,1))
       if fB(fbx,6)>0
          UID = fB(fbx,5);
          fax = find(fA(:,5)==UID);
          U_ms = (fB(fbx,1)-fA(fax,1))/dt;
          V_ms = (fB(fbx,2)-fA(fax,2))/dt;
          fA(fax,3) = U_ms;
          fA(fax,4) = V_ms;
       end
    end  
    vels = fA(:,3:4);
    vels(vels==0) = nan; % remove zero velocities on track ends
    fA(:,3:4) = vels;
    snapshots{i} = fA;
    snapminUID(i) = min(snapshots{i}(:,5));
    snapmaxUID(i) = max(snapshots{i}(:,5));
end
if isempty(snapshots{nframes})
    snapminUID(nframes) = max(snapminUID);
    snapmaxUID(nframes) = max(snapmaxUID);
else
    snapminUID(nframes) = min(snapshots{nframes}(:,5));
    snapmaxUID(nframes) = max(snapshots{nframes}(:,5));
end

% Forming Particle Tracks
% Previously all particle locations are sorted by frame number, 
% now they are going to be sorted by UID - particle tracks
% [X_m Y_m U_ms V_ms UID lifetime frame#];

ntracks = maxUIDfB;
tracks = cell(1,ntracks);
notempty = ones(1,ntracks);

for i=1:ntracks
    snapmin = find(snapmaxUID>=i,1,'first');
    snapmax = find(snapminUID<=i,1,'last');
    if isempty(snapmax < snapmin) || snapmax < snapmin     % track is empty (only one particle)
        notempty(i) = 0;
    else
        for j=snapmin:snapmax % frame number
            if ~isempty(snapshots{j})
                index = find(snapshots{j}(:,5)==i); % locate UID in frame
                if ~isempty(index)  
                    if size(snapshots{j},2) > 6
                        tracks{i} = vertcat(tracks{i},[snapshots{j}(index,1:6),j,snapshots{j}(index,7:end)]);
                    else
                        tracks{i} = vertcat(tracks{i},[snapshots{j}(index,:),j]);
                    end
                end
            end
        end
        if ~isempty(tracks{i})
            Udot_ms2 = gradient(tracks{i}(:,3))/dt;
            Vdot_ms2 = gradient(tracks{i}(:,4))/dt;
        else
            Udot_ms2 = [];
            Vdot_ms2 = [];
        end
        if size(tracks{i},2) > 7
            tracks{i} = [tracks{i}(:,1:7) Udot_ms2 Vdot_ms2 tracks{i}(:,8:end)];
        else
            tracks{i} = [tracks{i} Udot_ms2 Vdot_ms2];
        end
        tracks{i}(:,7) = tracks{i}(:,7) - 2;
    end
end
tracks = tracks(logical(notempty));  % remove empty tracks
ntracks = length(tracks);

% toc(tstart)


end
