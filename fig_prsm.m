figure(1)
plot(TM_sto(1:length(TM_sto)), sto_PS,'- k','LineWidth',1);
hold on;
plot(TM_ASP(1:length(TM_ASP)), ASP_PS,'Color',[0.1,0.4,0.8],'LineWidth',1);
hold on;
plot(TM_SARAH(1:length(TM_SARAH)), SARAH_PS,'- m','LineWidth',1);
hold on;
plot(TM_PDHG(1:length(TM_PDHG)), PDHG_PS,'- g','LineWidth',1);
hold on;
plot(TM_AS(1:length(TM_AS)), AS_PS,'- r','LineWidth',1);
legend('sto-ADMM','ASP-ADMM','SARAH-ADMM','S-PDHG','AS-PRSM','Location','southoutside','Orientation','horizontal');
legend('boxoff')
xlabel('time (second)');
ylabel('PSNR');

zoom_x = [1700, 2100]; 
zoom_y = [24.2, 25.3];  

ax_main = gca;
main_pos = get(ax_main, 'Position'); 

new_main_width = 0.65; 
set(ax_main, 'Position', [main_pos(1), main_pos(2), new_main_width, main_pos(4)]);

main_pos = get(ax_main, 'Position');

zoom_height = main_pos(4); 
zoom_width = 0.08; 
zoom_left = main_pos(1) + main_pos(3) + 0.02; 
zoom_bottom = main_pos(2); 


ax_zoom = axes('Position', [zoom_left, zoom_bottom, zoom_width, zoom_height]);
hold on;

plot(TM_sto(1:length(TM_sto)), sto_PS,'- k','LineWidth',1.5);
plot(TM_ASP(1:length(TM_ASP)), ASP_PS,'Color',[0.1,0.4,0.8],'LineWidth',1.5);

plot(TM_SARAH(1:length(TM_SARAH)), SARAH_PS,'- m','LineWidth',1.5);

plot(TM_PDHG(1:length(TM_PDHG)), PDHG_PS,'- g','LineWidth',1.5);
plot(TM_AS(1:length(TM_AS)), AS_PS,'- r','LineWidth',1);
xlim(zoom_x);
ylim(zoom_y);

grid off;
box on;
set(gca, 'FontSize', 9,'YAxisLocation', 'right');

axes(ax_main);
rectangle('Position', [zoom_x(1), zoom_y(1), diff(zoom_x), diff(zoom_y)], ...
          'EdgeColor', 'black', 'LineWidth', 2, 'LineStyle', '-');
x_lim_main = xlim(ax_main);
y_lim_main = ylim(ax_main);

main_rect_top_x = (zoom_x(2) - x_lim_main(1)) / diff(x_lim_main);
main_rect_top_y = (zoom_y(2) - y_lim_main(1)) / diff(y_lim_main);
main_rect_bottom_x = (zoom_x(2) - x_lim_main(1)) / diff(x_lim_main);
main_rect_bottom_y = (zoom_y(1) - y_lim_main(1)) / diff(y_lim_main);

main_point_top_x = main_pos(1) + main_rect_top_x * main_pos(3);
main_point_top_y = main_pos(2) + main_rect_top_y * main_pos(4);

zoom_point_top_x = zoom_left;
zoom_point_top_y = zoom_bottom + zoom_height; 

annotation('line', [main_point_top_x, zoom_point_top_x], [main_point_top_y, zoom_point_top_y], ...
           'Color', 'black', 'LineWidth', 1.5, 'LineStyle', '-');
%%
rxs = reshape(xs,N);
idd = 10;

xt = rxs(:,:,idd);
x1 = x_ASP(:,:,idd);
x3 = x_sto(:,:,idd);
x5 = x_SARAH(:,:,idd);
x7 = x_PDHG(:,:,idd);
x9 = x_AS(:,:,idd);


psnr_ASP = 10*log10(256^2./(norm(x1(:) - xt(:))).^2);
psnr_sto = 10*log10(256^2./(norm(x3(:) - xt(:))).^2);
psnr_SARAH = 10*log10(256^2./(norm(x5(:) - xt(:))).^2);
psnr_PDHG = 10*log10(256^2./(norm(x7(:) - xt(:))).^2);
psnr_AS = 10*log10(256^2./(norm(x9(:) - xt(:))).^2);

%%
figure(2);
ha = tight_subplot(2, 3, [0.05 0.01], [0.05 0.05], [0.05 0.05]);

axes(ha(1));
imshow(rxs(:,:,idd), []);
title('ground truth', 'FontSize', 10);

axes(ha(2));
imshow(x_ASP(:,:,idd), []);
title(['ASP-ADMM=', num2str(psnr_ASP)], 'FontSize', 9);

axes(ha(3));
imshow(x_AS(:,:,idd), []);
title(['AS-PRSM=', num2str(psnr_AS)], 'FontSize', 9);

axes(ha(4));
imshow(x_sto(:,:,idd), []);
title(['sto-ADMM=', num2str(psnr_sto)], 'FontSize', 9);

axes(ha(5));
imshow(x_SARAH(:,:,idd), []);
title(['SARAH-ADMM=', num2str(psnr_SARAH)], 'FontSize', 9);

axes(ha(6));
imshow(x_PDHG(:,:,idd), []);
title(['S-PDHG=', num2str(psnr_PDHG)], 'FontSize', 9);