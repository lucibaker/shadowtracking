% preprocess_MP sensitivity
function mean_dist_error_mm = preprocess_MP_sens(n)

% clear
% close all

gdrive_path = 'C:\Users\ljbak\My Drive\';  'H:\My Drive\'; % %  'G:\My Drive\';  % 
addpath([gdrive_path 'MATLAB\fm-toolbox'])
expt_string = '220613';  % expt set
% n = 8;  % run number

% load experiment params
warning off
run_params = readtable(sprintf('%sMP in OSBL\\imaging expts\\run%s\\sens_parameters_%s.xlsx',gdrive_path,expt_string,expt_string));
warning on

fprintf('particle type = %s\n', run_params.ParticleType{n});

nonsphere = strncmp(run_params.ParticleType{n},'d',1) || strncmp(run_params.ParticleType{n},'r',1);

load(sprintf('centers_sens%02d.mat',n))

if nonsphere
    %% compute angles
    
    pxyz = cell(size(centers));
    Ntrials = 1; %10;
    % del_lr_array = linspace(0, 1.5, Ntrials); % linspace(0, 5e-3, Ntrials); %
    % del_dr_array = linspace(0, 1.5, Ntrials); % linspace(0, 2e-3, Ntrials);
    % figure;
    
    for k = 1:Ntrials
        Dp = run_params.Dp_m(n);
        if strncmp(run_params.ParticleType{n},'r',1)
            for i = 1:length(centers)
                if isempty(angles{i}); continue; end
                % Rod angles (from raw tracks)
                th0r = angles{i}(:,1); %tracks0(:,10);
                lr = angles{i}(:,2); %tracks0(:,11);
            
                % subtract 1st-percentile l_r, scaled by rod angle
            %     [N_lr,edges_lr] = histcounts(lr,200,'Normalization','cdf'); 
            %     del_lr0 = mean(tracks0(:,12)); %5e-3; %edges_lr(find(N_lr>.99,1)+1) - Dp;  %min(lr); % [m]
            %     del_lr = 1e-3;%mean(tracks0(:,12));
            %     % del_lr = edges_lr(find(N_lr>.01,1)+1);  %min(lr); % [m]
            % %     lr_shifted = max([(lr-del_lr.*(Dp-lr)/(Dp)), zeros(size(lr))],[],2);
            %     lr_shifted = max([(lr - del_lr0), zeros(size(lr))],[],2);
            %     lr_shifted = max([(lr_shifted - del_lr.*(Dp-lr_shifted)/(Dp)), zeros(size(lr_shifted))],[],2);
    %             del_lr0 = del_lr_array(k)*mean(errchk{i}); %tracks0(:,12)); % 
                del_lr0 = run_params.K(n)*mean(errchk{i});
                del_lr = run_params.d(n); %0; %del_lr_array(k); % 2.5e-3;
                lr_shifted = lr - del_lr0;
                lr_shifted = lr_shifted - del_lr.*(Dp-lr_shifted)/Dp;
            
                % compute angle cosines
                pxyz{i} = [lr_shifted/(Dp).*cos(th0r), ... % p_x
                lr_shifted/(Dp).*sin(th0r), ... % p_z
                sqrt(1 - (lr_shifted/(Dp)).^2)]; % p_y
            end
            
        elseif strncmp(run_params.ParticleType{n},'d',1)
            for i = 1:length(centers)
                if isempty(angles{i}); continue; end
                % Disk angles (from raw tracks)
                th0primer = angles{i}(:,1); %tracks0(:,10);
                dr = angles{i}(:,2); %tracks0(:,11);
            
                % subtract 1st-percentile d_r, scaled by disk angle
            %     [N_dr,edges_dr] = histcounts(dr,200,'Normalization','cdf');
            %     del_dr0 = mean(tracks0(:,12)) - Dp; 
            %     del_dr = 1e-3; %edges_dr(find(N_dr>.01,1)+1);  %min(dr);
            % %     del_dr = edges_dr(find(N_dr>.01,1)+1);  %min(dr);
            %     dr_shifted = max([(dr - del_dr0), zeros(size(dr))],[],2);
            %     dr_shifted = max([(dr_shifted - del_dr.*(Dp-dr_shifted)/(Dp)), zeros(size(dr_shifted))],[],2);
    %             del_dr0 = del_dr_array(k)*mean(errchk{i}) - Dp; %mean(tracks0(:,12)) - Dp; 
                del_dr0 = run_params.K(n)*mean(errchk{i}) - Dp;
                del_dr = run_params.d(n); %0; %del_dr_array(k); % 0.5e-3; % 
                dr_shifted = dr - del_dr0;
                dr_shifted = dr_shifted - del_dr.*(Dp-dr_shifted)/Dp;
                
                % compute angle cosines
                pxyz{i} = [sin(th0primer).*sqrt(1 - (dr_shifted/(Dp)).^2), ... % p_x
                    cos(th0primer).*sqrt(1 - (dr_shifted/(Dp)).^2).*-sign(th0primer), ... % p_z
                    dr_shifted/(Dp)]; % p_y
            end
            
        end
        
        %% calculate error
        % read known orientations from spreadsheet
        pxyz_true = table2array( readtable(sprintf('%sMP in OSBL\\imaging expts\\run%s\\sens_parameters_%s.xlsx',gdrive_path,expt_string,expt_string), ...
            'Sheet','orientations','Range',sprintf('run_%d',n)) );
        
        % error between measured and known orientations
        pxyz_error = zeros(0,3);
        pxyz_imag = 0;
        for i = 1:length(centers)
            if ~isempty(pxyz{i})
                for j = 1:size(pxyz{i},1)
                    pxyz_error = [pxyz_error; abs(pxyz{i}(j,:)) - pxyz_true(i,:)];
                end
                pxyz_imag = pxyz_imag + sum(logical(imag(pxyz{i}(:,3))));
            end
        end
        
        % average error and % imaginary in each image
        pxyz_error_avg = mean(abs(pxyz_error),1);
        pxyz_imag = pxyz_imag/size(pxyz_error,1);
    
        disp('mean error  px  pz  py:'); disp(pxyz_error_avg)
    
    %     subplot(121)
    %     if strncmp(run_params.ParticleType{n},'r',1)
    %         plot(del_lr_array(k), pxyz_error_avg(1), 'r+', ...
    %             del_lr_array(k), pxyz_error_avg(2), 'b+', ...
    %             del_lr_array(k), pxyz_error_avg(3), 'k+')
    %         xlabel('K multiplied to del lr0'); hold on
    %     else
    %         plot(del_dr_array(k), pxyz_error_avg(1), 'r+', ...
    %             del_dr_array(k), pxyz_error_avg(2), 'b+', ...
    %             del_dr_array(k), pxyz_error_avg(3), 'k+')
    %         xlabel('K multiplied to del dr0'); hold on
    %     end
    %     legend('px','pz', 'py');
    %     ylabel('avg error'); ylim([0 1])
    % 
    %     subplot(122)
    %     if strncmp(run_params.ParticleType{n},'r',1)
    %         plot(del_lr_array(k), pxyz_imag, 'k*')
    %         xlabel('del lr'); hold on
    %     else
    %         plot(del_dr_array(k), pxyz_imag, 'k*')
    %         xlabel('del dr'); hold on
    %     end
    %     ylabel('fraction imaginary'); ylim([0 1])
    %     sgtitle(run_params.ParticleType(n))
    
    %     keyboard
    end
end

%% error on distances
% plot to see the order of detected particles (for entering the
% interparticle distances in spreadsheet)
% figure; plot(centers{1}(:,1),centers{1}(:,2),'*-'); axis equal 

Np_true = 10; % true number of particles
Np = size(centers{1}(:,1),1); % detected number of particles

% measured distances
xp = repmat(centers{1}(:,1),[1,Np]);
dxp = xp - xp';
yp = repmat(centers{1}(:,2),[1,Np]);
dyp = yp - yp';
dists = sqrt(dxp.^2 + dyp.^2); % [m]
dists = padarray(dists,Np_true-Np*ones(2,1),'post');

% true distances
dists_true = table2array( readtable(sprintf('%sMP in OSBL\\imaging expts\\run%s\\sens_parameters_%s.xlsx',gdrive_path,expt_string,expt_string), ...
    'Sheet','distances','Range',sprintf('run_%d_dist',n)) )/1000; % [m]

dist_error = dists - dists_true;
mean_dist_error_mm = mean(abs(dist_error(:)),'omitnan')*1000;
disp('mean error distance mm')
disp(mean_dist_error_mm)

% keyboard