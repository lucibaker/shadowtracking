% surface elevation
n = 1;
load(['C:\Users\ljbak\Documents\BACKUPS\Surface Elevation\5 separated dewarped runs\run_' num2str(n) '_dewarped.mat']);
run_dewarp = cat(1,run_dewarp(:,1:end/4), run_dewarp(:,end/2+1:3*end/4), run_dewarp(:,end/4+1:end/2), run_dewarp(:,3*end/4+1:end));

