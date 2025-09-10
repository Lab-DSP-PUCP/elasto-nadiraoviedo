%% pre-load
clear, clc
close all
%%
freq = 200:20:360;
load stats2.mat
colors = [1 0 0; 0 0.6 0; 0 0 0.9; 0.5 0.5 0.5; 0.3 0.3 0.3];

%% SWS mean values comparison

close all

% --------------------- Inclusion -----------------------
figure('Position',[200 200 400 400]),%1
% t_fig = tiledlayout(1,2);
% nexttile;
TOF= yline(5.1,'--', 'LineWidth',1, 'Color',uint8([17 17 17]));
ylim([2 13]), xlim([180 380])
grid on, hold on
p1 = errorbar(freq,Info_inc_AMFM.SWS,Info_inc_AMFM.std,'o-', 'MarkerSize',5, 'LineWidth',2, 'MarkerFaceColor','auto');
p2 = errorbar(freq,Info_inc_FSST.SWS,Info_inc_FSST.std,'o-', 'MarkerSize',5, 'LineWidth',2, 'MarkerFaceColor','auto');
hold off
xlabel('Frequency [Hz]'), ylabel('SWS [m/s]')
legend([p1,p2,TOF],{'AMFM','FSST','TOF'}, ...
    'Location','northoutside','Orientation','horizontal');
%title('STFT')
title('Inclusion')
ax = gca; ax.FontSize = 25;
ax.Position = [0.2 0.2 0.7 0.7];

figure('Position',[200 200 400 400]),%2
% nexttile;
ylim([0 35]), xlim([180 380]);
grid on, hold on
p2 = plot(freq,Info_inc_AMFM.CV,'o-','MarkerSize',5,'LineWidth',2, "MarkerFaceColor","auto");
p1 = plot(freq,Info_inc_FSST.CV,'o-','MarkerSize',5,'LineWidth',2, "MarkerFaceColor","auto");
hold off
xlabel('Frequency [Hz]'), ylabel('CV [%]')
legend([p2,p1],{'AMFM','FSST'}, ...
    'Location','northoutside','Orientation','horizontal');
ax = gca; ax.FontSize = 25;
ax.Position = [0.2 0.2 0.7 0.7];
title('Inclusion')
% saveas(gcf,'CV_inc.jpg')



%---------------- Background -----------
figure('Position',[200 200 400 400]),%3
% t_fig = tiledlayout(1,2);
% nexttile;
TOF= yline(3.45,'--', 'LineWidth',1, 'Color',uint8([17 17 17]));
ylim([2 6]), xlim([180 380])
grid on, hold on
p1 = errorbar(freq,Info_back_AMFM.SWS,Info_back_AMFM.std,'o-', 'MarkerSize',5, 'LineWidth',2, 'MarkerFaceColor','auto');
p2 = errorbar(freq,Info_back_FSST.SWS,Info_back_FSST.std,'o-', 'MarkerSize',5, 'LineWidth',2, 'MarkerFaceColor','auto');
hold off
xlabel('Frequency [Hz]'), ylabel('SWS [m/s]')
legend([p1,p2,TOF],{'AMFM','FSST','TOF'}, ...
    'Location','northoutside','Orientation','horizontal');
%title('STFT')
title('Background')
ax = gca; ax.FontSize = 25;
ax.Position = [0.2 0.2 0.7 0.7];


figure('Position',[200 200 400 400]),%4
% nexttile;
ylim([0 35]), xlim([180 380]);
grid on, hold on
p1 = plot(freq,Info_back_AMFM.CV,'o-','MarkerSize',5,'LineWidth',2, "MarkerFaceColor","auto");
p2 = plot(freq,Info_back_FSST.CV,'o-','MarkerSize',5,'LineWidth',2, "MarkerFaceColor","auto");
hold off
xlabel('Frequency [Hz]'), ylabel('CV [%]')
legend([p1,p2],{'AMFM','FSST'}, ...
    'Location','northoutside','Orientation','horizontal');
title('Background')
ax = gca; ax.FontSize = 25;
ax.Position = [0.2 0.2 0.7 0.7];


%% Bias (valor absoluto)
% Info_inc_PD.bias = abs(Info_inc_PD.bias);
% Info_inc_RW.bias = abs(Info_inc_RW.bias);
% Info_inc_CWT.bias = abs(Info_inc_CWT.bias);
% Info_back_PD.bias = abs(Info_back_PD.bias);
% Info_back_RW.bias = abs(Info_back_RW.bias);
% Info_back_CWT.bias = abs(Info_back_CWT.bias);
% 
% Info_inc_STFT.bias = abs(Info_inc_STFT.bias);
% Info_back_STFT.bias = abs(Info_back_STFT.bias);


% --------------------- Inclusion -----------------------
figure('Position',[200 200 400 400]),%5
% t_fig = tiledlayout(1,2);
% nexttile;
ylim([-10 70]), xlim([180 380]);
grid on, hold on
yline(0,'k--', 'LineWidth',1.5)
% p_PD = plot(freq,Info_inc_PD.bias,'o-','MarkerSize',5,'LineWidth',2,'Color',"#EDB120", "MarkerFaceColor","auto");
p1 = plot(freq,Info_inc_AMFM.bias,'o-','MarkerSize',5,'LineWidth',2, "MarkerFaceColor","auto");
p2 = plot(freq,Info_inc_FSST.bias,'o-','MarkerSize',5,'LineWidth',2, "MarkerFaceColor","auto");
hold off
xlabel('Frequency [Hz]'), ylabel('Bias [%]')
legend([p1,p2],{'AMFM','FSST'}, ...
    'Location','northoutside','Orientation','horizontal');
ax = gca; ax.FontSize = 25;
ax.Position = [0.2 0.2 0.7 0.7];
title('Inclusion')
%saveas(gcf,'bias_inc.jpg')
% nexttile;


%---------------- Background -----------
figure('Position',[200 200 400 400])%6
ylim([-5 30]), xlim([180 380]);
grid on, hold on
yline(0,'k--', 'LineWidth',1.5)
p1 = plot(freq,Info_back_AMFM.bias,'o-','MarkerSize',5,'LineWidth',2, "MarkerFaceColor","auto");
p2 = plot(freq,Info_back_FSST.bias,'o-','MarkerSize',5,'LineWidth',2, "MarkerFaceColor","auto");
hold off
xlabel('Frequency [Hz]'), ylabel('Bias [%]')
legend([p1,p2],{'AMFM','FSST'}, ...
    'Location','northoutside','Orientation','horizontal');
ax = gca; ax.FontSize = 25;
ax.Position = [0.2 0.2 0.7 0.7];
title('Background')
% saveas(gcf,'bias_back.jpg')

%% CNR

figure('Position',[200 200 400 400]),%8
% t_fig = tiledlayout(1,2);
% nexttile;
ylim([7 80]), xlim([180 380]);
grid on, hold on

CNR_AMFM = 20*log10(2*(Info_inc_AMFM.SWS - Info_back_AMFM.SWS).^2 ./ ...
    (Info_inc_AMFM.std.^2 + Info_back_AMFM.std.^2) );
CNR_FSST = 20*log10(2*(Info_inc_FSST.SWS - Info_back_FSST.SWS).^2 ./ ...
    (Info_inc_FSST.std.^2 + Info_back_FSST.std.^2) );
% 
p1 = plot(freq,CNR_AMFM,'o-','MarkerSize',5,'LineWidth',2, "MarkerFaceColor","auto");
p2 = plot(freq,CNR_FSST,'o-','MarkerSize',5,'LineWidth',2, "MarkerFaceColor","auto");
hold off
xlabel('Frequency [Hz]'), ylabel('CNR [dB]')
title('Contrast-to-Noise Ratio')
legend([p1,p2],{'AMFM','FSST'}, ...
    'Location','northoutside','Orientation','horizontal');
ax = gca; ax.FontSize = 25;
ax.Position = [0.2 0.2 0.7 0.7];

% figure('Position',[200 200 400 400]),%9
% % nexttile;
% ylim([0 8]),
% xlim([180 380]);
% grid on, hold on
% % 
% p1 = plot(freq,resAMFM,'o-','MarkerSize',5,'LineWidth',2, "MarkerFaceColor","auto");
% p2 = plot(freq,resFSST,'o-','MarkerSize',5,'LineWidth',2, "MarkerFaceColor","auto");
% hold off
% xlabel('Frequency [Hz]'), ylabel('R_{2080} [mm]')
% %title('Resolution')
% legend([p1,p2],{'AMFM','FSST'}, ...
%     'Location','northoutside','Orientation','horizontal');
% ax = gca; ax.FontSize = 25;
% ax.Position = [0.2 0.2 0.7 0.7];


%% Resolution
figure('Position',[200 200 400 400])%
ylim([0 15]),
xlim([180 380]);
grid on, hold on

% 
p1 = plot(freq,resAMFM,'o-','MarkerSize',5,'LineWidth',2, "MarkerFaceColor","auto");
p2 = plot(freq,resFSST,'o-','MarkerSize',5,'LineWidth',2, "MarkerFaceColor","auto");
hold off
xlabel('Frequency [Hz]'), ylabel('R_{2080} [mm]')
%title('Resolution')
legend([p1,p2,TOF],{'AMFM','FSST','TOF'}, ...
    'Location','northoutside','Orientation','horizontal');
ax = gca; ax.FontSize = 25;
ax.Position = [0.2 0.2 0.7 0.7];


%% one estimator per figure
%close all

%nexttile
figure('Position',[200 200 400 400])%11
yline(5.1,'--', 'LineWidth',1, 'Color',"#D95319");
ylim([3 6]), xlim([180 380])
grid on, hold on
yline(3.45,'--', 'LineWidth',1, 'Color',"#0072BD");
p1 = errorbar(freq,Info_back_AMFM.SWS,Info_back_AMFM.std,'o-', 'MarkerSize',5, 'LineWidth',2, 'MarkerFaceColor','auto');
p2 = errorbar(freq,Info_inc_AMFM.SWS,Info_inc_AMFM.std,'o-', 'MarkerSize',5, 'LineWidth',2, 'MarkerFaceColor','auto');
hold off
xlabel('Frequency [Hz]'), ylabel('Shear Wave Speed [m/s]')
legend([p2,p1],{'Inc.','Back.'},'Location','northeast');
title('AMFM')
ax = gca; ax.FontSize = 12;
ax.Position = [0.2 0.2 0.7 0.7];

%nexttile
figure('Position',[200 200 400 400])%12
yline(5.1,'--', 'LineWidth',1.5, 'Color',"#D95319");
ylim([3 6]), xlim([180 380])
grid on, hold on
yline(3.45,'--', 'LineWidth',1.5, 'Color',"#0072BD");
p1 = errorbar(freq,Info_back_FSST.SWS,Info_back_FSST.std,'o-', 'MarkerSize',5, 'LineWidth',2, 'MarkerFaceColor','auto');
p2 = errorbar(freq,Info_inc_FSST.SWS,Info_inc_FSST.std,'o-', 'MarkerSize',5, 'LineWidth',2, 'MarkerFaceColor','auto');
hold off
xlabel('Frequency [Hz]'), ylabel('Shear Wave Speed [m/s]')
legend([p2,p1],{'Inc.','Back.'},'Location','northeast');
title('FSST')
ax = gca; ax.FontSize = 12;
ax.Position = [0.2 0.2 0.7 0.7];
