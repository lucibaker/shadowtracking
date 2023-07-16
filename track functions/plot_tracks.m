% plot tracks
figure;
track_ids = 1:500;
c = jet(length(track_ids));
for i = 1:length(track_ids)
    idx = tracks1(:,5)==track_ids(i);
    c_idx = ceil(rand*length(track_ids)); 
    plot(tracks1(idx,1),tracks1(idx,2),'.','color',c(c_idx,:));
    hold on
end
axis equal; axis([-.5 .5 -.45 .05]);
xlabel('x [m]'); ylabel('y [m]')