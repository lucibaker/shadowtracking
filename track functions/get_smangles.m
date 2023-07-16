function [smangles, smangles_cont] = get_smangles(tracks0, kernel, dt, particle_type, Dp, varargin)
% compute particle orientations and time derivatives. if length(kernel)>1,
% returns acceleration variance vs kernel size.
% 
% inputs:
% tracks0: unsmoothed tracks (array)
% kernel: smoothing kernel width (scalar or vector)
% dt: time interval between frames
% particle_type: particle type id code (eg r5, d10)
% Dp: particle major axis length [m]
% (d: scaled offset value del_dr or del_lr [m])
% (K: multiplicative factor for constant offset value ie. K*del_dr0 or K*del_lr0)

%% compute angles

% load rescaling parameters
if ~isempty(varargin)
    d = varargin{1};
    K = varargin{2};
else
    d = nan;
    K = 1;
end

if strncmp(particle_type,'r',1)
    % Rod angles (from raw tracks)
    th0r = tracks0(:,10);
    lr = tracks0(:,11);

    % subtract offset del_lr0, scaled by rod angle
    del_lr0 = K*mean(tracks0(:,12));
    del_lr = 2.5e-3; 
    if ~isnan(d); del_lr = d; end %  choose to use d0 if passed in
    lr_shifted = lr - del_lr0;
    lr_shifted = lr_shifted - del_lr.*(Dp-lr_shifted)/Dp;
    
    % plot lr histogram
    [lr_pdf,lr_range] = pdf_var(lr*1000,100,0);
    [lr_pdf2,lr_range2] = pdf_var(lr_shifted*1000,100,0);
    figure; plot(lr_range,lr_pdf,'k.-',lr_range2,lr_pdf2,'r.-')
    ylabel('PDF'); xlabel('$d$ [mm]'); grid on; %ylim([0 1]); %legend('$d$','$d_{corr}$','location','nw'); xlim([0 3.5])
    title('Rods')
    goodplot([4 3.5])

    % compute angle cosines
    pxyz = calc_orient(th0r,lr_shifted,Dp,particle_type); % px, pz, py
        
    rir = logical(~imag(pxyz(:,3)) & ~imag(pxyz(:,1))); % real idx, raw
    fprintf('%2.1f%% of raw p-hats real\n',sum(rir)/numel(rir)*100);   
    
else
    % Disk angles (from raw tracks)
    th0primer = tracks0(:,10);
    dr = tracks0(:,11);

    % subtract offset del_dr, scaled by disk angle
    del_dr0 = K*mean(tracks0(:,12)) - Dp; 
    del_dr = 0.5e-3; % 0;
    if ~isnan(d); del_dr = d; end %  choose to use d0 if passed in
    dr_shifted = dr - del_dr0;
    dr_shifted = dr_shifted - del_dr.*(Dp-dr_shifted)/Dp;
    
    % plot dr histogram
    [dr_pdf,dr_range] = pdf_var(dr*1000,100,0);
    [dr_pdf2,dr_range2] = pdf_var(dr_shifted*1000,100,0);
    figure; plot(dr_range,dr_pdf,'k.-',dr_range2,dr_pdf2,'r.-')
    xlabel('$d$ [mm]'); legend('$d$','$d_{corr}$'); grid on; 
    ylabel('PDF'); 
    title('Disks')
    goodplot([4 3.5])
    
    % compute angle cosines
    pxyz = calc_orient(th0primer,dr_shifted,Dp,particle_type);  % px, pz, py
    
    rir = logical(~imag(pxyz(:,2)) & ~imag(pxyz(:,1))); % real idx, raw
    fprintf('%2.1f%% of raw p-hats real\n',sum(rir)/numel(rir)*100);
    
end

%% resolve sign ambiguities
if strcmp(particle_type,'rods')
    sigma_Lp = 0.19e-3; % stdev of rod lengths [m]
    p1 = 0.15;
    p2 = 1 - sigma_Lp/(Dp);
else
    p1 = 0.15;
    p2 = 0.97;
end

plot_on = 0;
[pxyz, amb_list] = resolve_orientations(tracks0,pxyz,p1,p2,plot_on);

                        
%% temporal smoothing
% smangles: [px py pz px_dot py_dot pz_dot px_dd py_dd pz_dd] 
avar_k_ang = zeros(size(kernel));

for k = 1:length(kernel)
    fprintf(['kernel = ' num2str(kernel(k)) '...\n'])
    if kernel(k) > 1
        mintracks = kernel(k)*2+3;  % minimum length of tracks for smoothing (at least three points)
        smangles = zeros(0,9);
        for i = 1:max(tracks0(:,5))
            idx_i = find(tracks0(:,5)==i);
            if any(idx_i) && length(tracks0(idx_i,1))>=mintracks

                % sort track by lifetime counter (in case it is out of order after fixing a broken track)
                [~,idx_sort] = sort(tracks0(idx_i,6),'ascend'); 
                idx_i = idx_i(idx_sort);

                px_i = gauss_position(pxyz(idx_i,1),kernel(k));
                py_i = gauss_position(pxyz(idx_i,2),kernel(k));
                pz_i = gauss_position(pxyz(idx_i,3),kernel(k));
                px_d_i = gauss_velocity(pxyz(idx_i,1),kernel(k),dt);
                py_d_i = gauss_velocity(pxyz(idx_i,2),kernel(k),dt);
                pz_d_i = gauss_velocity(pxyz(idx_i,3),kernel(k),dt);
                px_dd_i = gauss_accel(pxyz(idx_i,1),kernel(k),dt);
                py_dd_i = gauss_accel(pxyz(idx_i,2),kernel(k),dt);
                pz_dd_i = gauss_accel(pxyz(idx_i,3),kernel(k),dt);

                smangles = [smangles; px_i, py_i, pz_i, px_d_i, py_d_i, pz_d_i, px_dd_i, py_dd_i, pz_dd_i];
            end
        end
        avar_k_ang(k) = std(sqrt(smangles(:,7).^2 + smangles(:,8).^2 + smangles(:,9).^2),'omitnan')^2;
    
    else
        smangles(:,1) = pxyz(:,1);
        smangles(:,2) = pxyz(:,2);
        smangles(:,3) = pxyz(:,3);
        avar_k_ang(k) = nan;
        warning('Kernel <= 1: no angular velocities calculated')
    end
    
end

if length(kernel) > 1
    n_end = 15;
    figure; semilogy(kernel*dt/t_plus,avar_k_ang*t_plus^2,'k+','linewidth',1,'markersize',4); hold on
    xlabel('$t_k^+$'); ylabel('var($\ddot{p}^+$)'); 
    n0 = find(kernel==17); % index of minimum kernel size for good exponential fit
    P = polyfit(kernel(n0:n_end)*dt/t_plus,log(avar_k_ang(n0:n_end)*t_plus^2),1); 
    semilogy(kernel*dt/t_plus, exp(P(1)*kernel*dt/t_plus + P(2)), 'k-','linewidth',1);
    semilogy(kernel(n0)*dt/t_plus,avar_k_ang(n0)*t_plus^2,'ro','linewidth',1,'markersize',8);
    xlim([0 (kernel(n_end)+1)*dt/t_plus]); %set(gca,'XTick',2:2:22); 
    goodplot([5 3.5])

    fprintf(['Best kernel size = ' num2str(kernel(n0)) '\n'])
    savefig('avar_k_ang.fig')
    keyboard
end

if length(kernel) == 1
    % remove particles whose p_hat vector deviates from unit length (by more than a threshold) 
    dev_thres = 0.1; % deviation allowed from unit length
    p_length = sqrt(sum(smangles(:,1:3).^2, 2));
    angle_check_idx1 = logical(abs(1 - p_length) > dev_thres);
    smangles(angle_check_idx1,:) = nan; 
 
    if kernel > 1
        % discard spurious angular accelerations
        pdd_thres = inf; % max p_ddot threshold
        smangles(smangles(:,7) > pdd_thres | smangles(:,8) > pdd_thres | smangles(:,9) > pdd_thres,7:9) = nan;
        
        % smangles: [ [p] [p_dot] [p_dd] [om_t] [al_t] str sta]
        % tumbling components of angular velocity and angular acceleration
        if size(smangles,2) == 9
            smangles = [smangles, zeros(size(smangles,1),6)];
        end
        smangles(:,10:12) = cross(smangles(:,1:3),smangles(:,4:6)); % tumbling ang. vel.
        smangles(:,13:15) = cross(smangles(:,1:3),smangles(:,7:9)); % tumbling ang. accel. 
        smangles(:,16) = sum(smangles(:,4:6).^2,2); % sq_tumb_rate
        smangles(:,17) = sum(smangles(:,7:9).^2,2); % sq_tumb_accel
    end
  
    % flip angles back to original signs
    smangles_cont = smangles; % continuous version for autocorrs
    smangles(:,2) = smangles(:,2).*sign(smangles(:,1));
    smangles(:,1) = abs(smangles(:,1));

    % remove imaginary and negative p_hats
    imag_remove = logical(imag(smangles(:,1)) | imag(smangles(:,2)) | imag(smangles(:,3)));
    neg_remove = any(smangles(:,[1 3]) < 0, 2);
    fprintf('%2.1f%% of smpx and smpy neg; %2.1f%% of smphats imag\n', sum(neg_remove)/numel(neg_remove)*100, sum(imag_remove)/numel(imag_remove)*100)
    smangles(imag_remove | neg_remove,:) = nan;
    smangles_cont(imag_remove | neg_remove,:) = nan;
end

end