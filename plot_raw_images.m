A = cam_imread('CamB2-1000.tif',0);
figure; imagesc((A),[0 80]); axis equal tight; colormap gray
% set(gca,'YTick',[],'XTick',[]); 

A_sub = A(1640:1940, 230:530); %A(1175:1375, 935:1135); %
figure; imagesc(flipud(A_sub),[0 80]); axis equal tight; colormap gray
set(gca,'YTick',[],'XTick',[]); 

goodplot([3 3])

fn_str = 'C:\Users\ljbak\My Drive\MP in OSBL\results\ExiF paper\figs\CamB2-1000-r20';
% savefig([fn_str '.fig'])
% print([fn_str '.png'],'-dpng')