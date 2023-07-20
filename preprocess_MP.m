function [avar_k] = preprocess_MP(n)
% preprocess particle tracks:
% -clean up data
% -compute orientation (p vector)
% -apply temporal smoothing to tracks and orientations

%% SETUP - set these variables first
expt_path = 'H:\My Drive\MATLAB\shadowtracking\';  % path to experiment folder
expt_name = 'demo';  % experiment or dataset name 

plot_on = false;  % make plots?

kernel = 3; % width of gaussian kernel for temporal smoothing


%% get experiment parameters
% add path to shadowtracking functions
addpath([expt_path 'track functions'])

% load experiment params
warning off
run_params = readtable(sprintf('%s\\data_%s\\run_parameters_%s.xlsx', expt_path, expt_name, expt_name));
warning on

fprintf('\nwindspeed = %2.f m/s, particle type = %s\n', run_params.WindSpeed_m_s(n), ...
        run_params.ParticleType{n});

nonsphere = strncmp(run_params.ParticleType{n},'d',1) || strncmp(run_params.ParticleType{n},'r',1);

load(sprintf('%s\\outputs_%s\\tracks_run%02d.mat', expt_path, expt_name, n))

%% remove center bright spot
spot_x = [0 0.015]; 
spot_z = [-.3 -.275]; 

i = 1;
while i <= max(tracks0(:,5))
    idx = tracks0(:,5) == i;
    if sum(idx)
        if nanmean(tracks0(idx,1)) > spot_x(1) && nanmean(tracks0(idx,1)) < spot_x(2) && ...
                nanmean(tracks0(idx,2)) > spot_z(1) && nanmean(tracks0(idx,2)) < spot_z(2)
            % remove the track corresponding to the center bright spot
            [tracks0, tracklength0] = remove_track(tracks0, tracklength0, i);  
        else
            % don't remove track, move to next track id
            i = i+1;            
        end
    end
end

%% repair broken tracks  
searchrad = 10e-3; 
skipframes = 5;
[tracks0,tracklength0] = fix_tracks(tracks0,tracklength0,searchrad,1/run_params.imagingFreq_Hz(n),skipframes);

%% remove non-shadow particles whose tracks break at frame edges
frame_edges = [-0.22, -0.195, -0.017, 0.008, 0.206, 0.238]; % frame edges  
tol = 0.005; % tolerance (max dist from frame edge)
[tracks0, tracklength0] = filter_ghost_tracks(tracks0, tracklength0, frame_edges, tol);

%% flip streamwise coord and orientation (WASIRF-specific)
tracks0(:,1) = -tracks0(:,1);  % x
tracks0(:,3) = -tracks0(:,3);  % u
tracks0(:,10) = -tracks0(:,10);  % theta

%% smooth tracks
sm_fn = sprintf('%s\\outputs_%s\\smtracks_run%02d.mat', expt_path, expt_name, n);
[smtracks, smtracklength, avar_k] = smooth_tracks(tracks0,kernel,1/run_params.imagingFreq_Hz(n));
if length(kernel) > 1
    figure; semilogy(kernel,avar_k,'k.'); 
    xlabel('kernel [frames]'), ylabel('var(a) [m^2/s^4]')
    k_opt = 5; k_opt_idx = logical(kernel >= k_opt);
    P = polyfit(kernel(k_opt_idx),log(avar_k(k_opt_idx)),1);
    hold on;
    plot(kernel,exp(kernel*P(1)+P(2)),'r-');
    plot(k_opt,avar_k(find(k_opt_idx,1,'first')),'ro','markersize',6);
    ylim([0 1.5])
end

save(sm_fn,'smtracks','smtracklength','kernel');

ntracks = length(smtracklength);

% get smoothed angles
if nonsphere
    [smangles, smangles_cont] = get_smangles(tracks0,kernel,1/run_params.imagingFreq_Hz(n),run_params.ParticleType{n},run_params.Dp_m(n), ...
        run_params.d(n),run_params.K(n));
    save(sm_fn,'smangles','smangles_cont','-append');
end

%% preview tracks
if plot_on

    % track lengths
    figure; histogram(smtracklength,100)
    xlabel('track length [frames]'); ylabel('count')
    
    figure;
    track_ids = 1:200;% find(smtracklength>300); % round(linspace(1,ntracks,100));  % 1:30; %
    c = jet(length(track_ids));
    for i = 1:length(track_ids)
        idx = smtracks(:,5)==track_ids(i);
        c_idx = ceil(rand*length(track_ids)); % round(smtracklength(track_ids(i))/max(smtracklength(track_ids))*length(track_ids));
        plot(smtracks(idx,1),smtracks(idx,2),'.','color',c(c_idx,:));
        hold on
    end
    axis equal; axis([-.5 .5 -.45 .05]);
    xlabel('x [m]'); ylabel('y [m]')
    
    fidx = smtracks(:,7)==1;
    pts = plot(smtracks(fidx,1),smtracks(fidx,2),'ko');
    pause(1/10)
    for i = 2:1000
        fidx = smtracks(:,7)==i;
        delete(pts)
        pts = plot(smtracks(fidx,1),smtracks(fidx,2),'ko');
        pause(1/10)
    end
end
