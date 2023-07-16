function [pxyz, amb_list] = resolve_orientations(tracks0,pxyz,p1,p2,plot_on)
% resolve ambiguities in particle orientations when orientations hit the
% limits of their range 
% tracks0: array containing unsmoothed particle track data 
% pxyz: array containing unsmoothed particle orientation data [px py pz]
% amb_list: vector with same length as tracks0 containing the ambiguity 
%   flag for each observation (0 = no change, 1 = flip px & py, 2 = flip 
%   pz, 3 = flip px py & pz)
% p1: lower orientation range limit
% p2: upper orientation range limit
% plot_on: set to true to plot cases where particle orientation was flipped

amb_list = zeros(size(tracks0,1),1);

if plot_on
    figure; set(gcf,'Position',[680 92 1205 886]); subplot(311); 
end

for i = 1:max(tracks0(:,5))
    idx_i = find(tracks0(:,5)==i);
    if any(idx_i)
        % sort track by lifetime counter (in case it is out of order after fixing a broken track)
        [~,idx_sort] = sort(tracks0(idx_i,6),'ascend'); 
        idx_i = idx_i(idx_sort);
        
        if all(isreal(pxyz(idx_i,:))) 
            % original track angular acceleration
            pdd0 = [nan(1,3); diff(tracks0(idx_i,10:12),2,1); nan(1,3)]; 
            pdd0_mag = sum(pdd0(:,1:3).^2,2);
    
            % check for angle ambiguities
            amb_flag = find( (tracks0(idx_i,10) < p1 & islocalmin(tracks0(idx_i,10))) | ...
                (pxyz(idx_i,1) > p2 & islocalmax(pxyz(idx_i,1))) | ...
                (pxyz(idx_i,3) < p1 & islocalmin(pxyz(idx_i,3))) | ...
                (pxyz(idx_i,3) > p2 & islocalmax(pxyz(idx_i,3))) | ...
                (pxyz(idx_i,2) < -p2 & islocalmin(pxyz(idx_i,2))) | ...
                (pxyz(idx_i,2) > p2 & islocalmax(pxyz(idx_i,2))) | ...
                diff(sign([pxyz(idx_i,2);1])) );
    
            % if there are ambiguities, resolve by choosing the option with minimum angular acceleration
            if ~isempty(amb_flag)
                pxyz1 = pxyz(idx_i,:);
                pxyz2 = pxyz(idx_i,:);
                pxyz3 = pxyz(idx_i,:);
    
                for j = 1:length(amb_flag)
                    if amb_flag(j) < length(idx_i)
                        pxyz1(amb_flag(j)+1:end,1:2) = -pxyz1(amb_flag(j)+1:end,1:2); % change sign of px, py
                        pxyz2(amb_flag(j)+1:end,3) = -pxyz2(amb_flag(j)+1:end,3); % change sign of pz
                        pxyz3(amb_flag(j)+1:end,1:3) = -pxyz3(amb_flag(j)+1:end,1:3); % change sign of px, py, pz
    
                        % angular acceleration of altered tracks
                        pdd1 = [nan(1,3); diff(pxyz1,2,1); nan(1,3)];
                        pdd2 = [nan(1,3); diff(pxyz2,2,1); nan(1,3)];
                        pdd3 = [nan(1,3); diff(pxyz3,2,1); nan(1,3)];
                        pdd1_mag = sum(pdd1.^2,2);
                        pdd2_mag = sum(pdd2.^2,2);
                        pdd3_mag = sum(pdd3.^2,2);
    
                        % choose track with minimum pdd (angular acceleration)
                        [~,I] = min([pdd0_mag(amb_flag(j)+1), pdd1_mag(amb_flag(j)+1), pdd2_mag(amb_flag(j)+1), pdd3_mag(amb_flag(j)+1)]);
                        
                        % plot pdd and symmetry axis
                        if plot_on && I-1 > 0
                            subplot(311); plot(idx_i,tracks0(idx_i,10:12),idx_i(amb_flag(j)),0,'mo')
                            legend('p_x','p_y','p_z','Current frame','location','northeastoutside')
                            ylabel('p_i');  title(['Track ' num2str(i) ': choose option ' num2str(I-1)]); grid on; ylim([-1 1])
                            
                            subplot(312); plot(idx_i,pdd0_mag,'r',idx_i,pdd1_mag,'b',idx_i,pdd2_mag,'k',idx_i,pdd3_mag,'c',idx_i(amb_flag(j)),0,'mo','linewidth',1);
                            legend('0: Original','1: p_x, p_y reversed','2: p_z reversed','3: p_x, p_y, p_z reversed','location','northeastoutside')
                            xlabel('t [frames]'); ylabel('|d^2{\bfp}/dt^2| [frames^{-2}]'); %ylim([0 1e-3])
                            
                            subplot(313); 
                            if strcmp(particle_type,'rods')    
                                line([tracks0(idx_i,1) - Rp*tracks0(idx_i,10), tracks0(idx_i,1) + Rp*tracks0(idx_i,10)]'*1000, ...
                                    [tracks0(idx_i,2) - Rp*tracks0(idx_i,11), tracks0(idx_i,2) + Rp*tracks0(idx_i,11)]'*1000, ...
                                    [0 - Rp*tracks0(idx_i,12), 0 + Rp*tracks0(idx_i,12)]'*1000, ...
                                    'color','k','linewidth',1); hold on
                                line([tracks0(idx_i(amb_flag(j)),1) - Rp*tracks0(idx_i(amb_flag(j)),10), tracks0(idx_i(amb_flag(j)),1) + Rp*tracks0(idx_i(amb_flag(j)),10)]'*1000, ...
                                    [tracks0(idx_i(amb_flag(j)),2) - Rp*tracks0(idx_i(amb_flag(j)),11), tracks0(idx_i(amb_flag(j)),2) + Rp*tracks0(idx_i(amb_flag(j)),11)]'*1000, ...
                                    [0 - Rp*tracks0(idx_i(amb_flag(j)),12), 0 + Rp*tracks0(idx_i(amb_flag(j)),12)]'*1000, ...
                                    'color','m','linewidth',2); hold off
                                pause
                                cla(gca);
                            else
                                plotCircle3D([tracks0(idx_i,1:2),zeros(size(idx_i))]*1000,tracks0(idx_i,10:12),Rp*1000,[.5 .5 .5],1); hold on
                                plotCircle3D([tracks0(idx_i(amb_flag(j)),1:2),0]*1000,tracks0(idx_i(amb_flag(j)),10:12),Rp*1000,'m',2); hold on
                                line([tracks0(idx_i,1) - Rp*tracks0(idx_i,10), tracks0(idx_i,1) + Rp*tracks0(idx_i,10)]'*1000, ...
                                    [tracks0(idx_i,2) - Rp*tracks0(idx_i,11), tracks0(idx_i,2) + Rp*tracks0(idx_i,11)]'*1000, ...
                                    [0 - Rp*tracks0(idx_i,12), 0 + Rp*tracks0(idx_i,12)]'*1000, ...
                                    'color','k','linewidth',1); hold on
                                line([tracks0(idx_i(amb_flag(j)),1) - Rp*tracks0(idx_i(amb_flag(j)),10), tracks0(idx_i(amb_flag(j)),1) + Rp*tracks0(idx_i(amb_flag(j)),10)]'*1000, ...
                                    [tracks0(idx_i(amb_flag(j)),2) - Rp*tracks0(idx_i(amb_flag(j)),11), tracks0(idx_i(amb_flag(j)),2) + Rp*tracks0(idx_i(amb_flag(j)),11)]'*1000, ...
                                    [0 - Rp*tracks0(idx_i(amb_flag(j)),12), 0 + Rp*tracks0(idx_i(amb_flag(j)),12)]'*1000, ...
                                    'color','m','linewidth',2); hold off
                                view(2); axis equal
                                pause
                                cla(gca);
                            end
                        end
    
                        % apply changes to tracks
                        switch(I-1)
                            case 0  
                                % revert to original track
                                pxyz1 = pxyz(idx_i,:);
                                pxyz2 = pxyz(idx_i,:);
                                pxyz3 = pxyz(idx_i,:);
    
                            case 1 
                                % change sign of px, py only
                                pxyz(idx_i,:) = pxyz1;
                                pxyz2 = pxyz1;
                                pxyz3 = pxyz1;
                                pdd0 = [nan(1,3); diff(pxyz(idx_i,:),2,1); nan(1,3)]; 
                                pdd0_mag = sum(pdd0.^2,2);
    
                            case 2 
                                % change sign of pz only
                                pxyz(idx_i,:) = pxyz2;
                                pxyz1 = pxyz2;
                                pxyz3 = pxyz2;
                                pdd0 = [nan(1,3); diff(pxyz(idx_i,:),2,1); nan(1,3)];
                                pdd0_mag = sum(pdd0.^2,2);
                                
                            case 3 
                                % change sign of px, py, pz
                                pxyz(idx_i,:) = pxyz3;
                                pxyz1 = pxyz3;
                                pxyz2 = pxyz3;
                                pdd0 = [nan(1,3); diff(pxyz(idx_i,:),2,1); nan(1,3)];
                                pdd0_mag = sum(pdd0.^2,2);
                        end
                        amb_list(idx_i(amb_flag(j))) = I-1;
                        
                    end
                end
            end
        end 
    end
end

end