close all;
clear, clc;
%% 
freq = 400:100:800;
colors = [1 0 0; 0 0.6 0; 0 0 0.9; 0.5 0.5 0.5];
%%
load ('stats 8-4.mat')
%%
load ('stats 8-2.mat')
%%
load ('stats 4-2.mat')
%%
load ('stats 8-8.mat')
%%
load ('stats 4-4.mat')
%%
load ('stats 2-2.mat')
%% comparison between estimators in mean values
close all
figure('Position',[200 200 400 400])
ylim([3 5]), xlim([350 850])
grid on, hold on
TOF=yline(8,'--', 'LineWidth',5, 'Color',uint8([17 17 17]));

p2 = errorbar(freq,Info_left_AMFM.SWS,Info_left_AMFM.std,'o-', 'MarkerSize',5, 'LineWidth',3, 'MarkerFaceColor','auto');
p4 = errorbar(freq,Info_left_FSST.SWS,Info_left_FSST.std,'o-', 'MarkerSize',5, 'LineWidth',3, 'MarkerFaceColor','auto');
hold off
xlabel('Frequency [Hz]'), ylabel('SWS [m/s]')
%legend([p1,p2,p3,p4,p5,TOF],{'STFT','CWT','ASLT','FSST','WSST','TOF'},'Location','northeast','NumColumns',2);
legend([p2,p4,TOF],{'AMFM','FSST','TOF'}, ...
    'Location','northoutside','Orientation','horizontal');
% title('WSST')
ax = gca; ax.FontSize = 28;
ax.LineWidth = 1.5;
ax.Position = [0.2 0.2 0.7 0.7];
%saveas(gcf, 'C:\Users\Joaquin\OneDrive\Escritorio\Joaquin\Journal\Journal Standing Waves final 2.0\miFigura.jpg')

%%
figure('Position',[200 200 400 400])
ylim([0.6 3.4]), xlim([350 850])
grid on, hold on
TOF=yline(2,'--', 'LineWidth',5, 'Color',uint8([17 17 17]));
p2 = errorbar(freq,Info_right_AMFM.SWS,Info_right_AMFM.std,'o-', 'MarkerSize',5, 'LineWidth',3,'MarkerFaceColor','auto');
p4 = errorbar(freq,Info_right_FSST.SWS,Info_right_FSST.std,'o-', 'MarkerSize',5, 'LineWidth',3, 'MarkerFaceColor','auto');
hold off
xlabel('Frequency [Hz]'), ylabel('SWS [m/s]')
legend([p2,p4,TOF],{'AMFM','FSST','TOF'}, ...
    'Location','northoutside','Orientation','horizontal');
% title('WSST')
ax = gca; ax.FontSize = 28;
ax.LineWidth = 1.5;
ax.Position = [0.2 0.2 0.7 0.7];
% saveas(gcf, 'C:\Users\Joaquin\OneDrive\Escritorio\Joaquin\Journal\Journal Standing Waves final 2.0\miFigura.jpg')
% 

%% Bias (valor absoluto)
% --------------------- Left -----------------------
figure('Position',[200 200 400 400])
ylim([-16 16]), xlim([350 850]);
grid on, hold on
yline(0,'k--', 'LineWidth',3)
p_AMFM = plot(freq,Info_left_AMFM.bias,'o-','MarkerSize',5,'LineWidth',3, "MarkerFaceColor","auto");
p_FSST = plot(freq,Info_left_FSST.bias,'o-','MarkerSize',5,'LineWidth',3, "MarkerFaceColor","auto");

hold off
xlabel('Frequency [Hz]'), ylabel('Bias [%]')
legend([p_AMFM ,p_FSST],{'AMFM','FSST'},...
    'Location','northoutside','Orientation','horizontal');
ax = gca; ax.FontSize = 28;
ax.LineWidth = 1.5;
ax.Position = [0.2 0.2 0.7 0.7];
%title('Left')
%saveas(gcf,'bias_inc.jpg')

%---------------- Right -----------
figure('Position',[200 200 400 400])
ylim([-16 16]), xlim([350 850]);
grid on, hold on
yline(0,'k--', 'LineWidth',3)
p_AMFM = plot(freq,Info_right_AMFM.bias,'o-','MarkerSize',5,'LineWidth',3,"MarkerFaceColor","auto");
p_FSST = plot(freq,Info_right_FSST.bias,'o-','MarkerSize',5,'LineWidth',3, "MarkerFaceColor","auto");
hold off
xlabel('Frequency [Hz]'), ylabel('Bias [%]')
legend([p_AMFM ,p_FSST],{'AMFM','FSST'},...
    'Location','northoutside','Orientation','horizontal');
ax = gca; ax.FontSize = 28;
ax.LineWidth = 1.5;
ax.Position = [0.2 0.2 0.7 0.7];
%title('Right')
%saveas(gcf,'bias_back.jpg')

%% CV
% --------------------- left -----------------------
figure('Position',[200 200 400 400])
ylim([0 30]), xlim([350 850]);
grid on, hold on
p_AMFM = plot(freq,Info_left_AMFM.CV,'o-','MarkerSize',5,'LineWidth',3, "MarkerFaceColor","auto");
p_FSST = plot(freq,Info_left_FSST.CV,'o-','MarkerSize',5,'LineWidth',3, "MarkerFaceColor","auto");
hold off
xlabel('Frequency [Hz]'), ylabel('CV [%]')
legend([p_AMFM,p_FSST],{'AMFM','FSST'},...
    'Location','northoutside','Orientation','horizontal');
ax = gca; ax.FontSize = 28;
ax.LineWidth = 1.5;
ax.Position = [0.2 0.2 0.7 0.7];
%title('Left')
%saveas(gcf,'CV_inc.jpg')

%---------------- right -----------
figure('Position',[200 200 400 400])
ylim([0 30]), xlim([350 850]);
grid on, hold on
p_AMFM = plot(freq,Info_right_AMFM.CV,'o-','MarkerSize',5,'LineWidth',3, "MarkerFaceColor","auto");
p_FSST = plot(freq,Info_right_FSST.CV,'o-','MarkerSize',5,'LineWidth',3, "MarkerFaceColor","auto");
hold off
xlabel('Frequency [Hz]'), ylabel('CV [%]')
legend([p_AMFM,p_FSST],{'AMFM','FSST'},...
    'Location','northoutside','Orientation','horizontal');
ax = gca; ax.FontSize = 28;
ax.LineWidth = 1.5;
ax.Position = [0.2 0.2 0.7 0.7];
% title('Right')
%saveas(gcf,'CV_back.jpg')

%% CNR
figure('Position',[200 200 400 400])
ylim([10 100]), xlim([350 850]);
grid on, hold on

CNR_AMFM = 20*log10(2*(Info_left_AMFM.SWS - Info_right_AMFM.SWS).^2 ./ ...
    (Info_left_AMFM.std.^2 + Info_right_AMFM.std.^2) );
CNR_FSST = 20*log10(2*(Info_left_FSST.SWS - Info_right_FSST.SWS).^2 ./ ...
    (Info_left_FSST.std.^2 + Info_right_FSST.std.^2) );
% 
p_AMFM = plot(freq,CNR_AMFM,'o-','MarkerSize',5,'LineWidth',3,"MarkerFaceColor","auto");
p_FSST = plot(freq,CNR_FSST,'o-','MarkerSize',5,'LineWidth',3, "MarkerFaceColor","auto");
hold off
xlabel('Frequency [Hz]'), ylabel('CNR [dB]')
%title('Contrast-to-Noise Ratio')
legend([p_AMFM,p_FSST],{'AMFM','FSST'},...
    'Location','northoutside','Orientation','horizontal');
ax = gca; ax.FontSize = 28;
ax.LineWidth = 1.5;
ax.Position = [0.2 0.2 0.7 0.7];

%% Resolution
figure('Position',[200 200 400 400])
ylim([0 5]),
xlim([350 850]);
grid on, hold on

% 
p_AMFM = plot(freq,resAMFM,'o-','MarkerSize',5,'LineWidth',3,"MarkerFaceColor","auto");
p_FSST = plot(freq,resFSST,'o-','MarkerSize',5,'LineWidth',3, "MarkerFaceColor","auto");
hold off
xlabel('Frequency [Hz]'), ylabel('R_{2080} [mm]')
%title('Resolution')
legend([p_AMFM,p_FSST],{'AMFM','FSST'},...
    'Location','northoutside','Orientation','horizontal');
ax = gca; ax.FontSize = 28;
ax.LineWidth = 1.5;
ax.Position = [0.2 0.2 0.7 0.7];



%-----------------------------------------------------------------------------------------------------------------------------------------------------------
% %% HOMOGENEUS
% %-----------------------------------------------------------------------------------------------------------------------------------------------------------
% %% one estimator per figure
% % close all
% figure('Position',[200 200 400 400])
% ylim([1 3]), xlim([350 850])
% grid on, hold on
% TOF=yline(2,'--', 'LineWidth',5, 'Color',uint8([17 17 17]));
% p2 = errorbar(freq,Info_AMFM.SWS,Info_AMFM.std,'o-', 'MarkerSize',5, 'LineWidth',3, 'MarkerFaceColor','auto');
% p4 = errorbar(freq,Info_FSST.SWS,Info_FSST.std,'o-', 'MarkerSize',5, 'LineWidth',3, 'MarkerFaceColor','auto');
% hold off
% xlabel('Frequency [Hz]'), ylabel('SWS [m/s]')
% legend([p2,p4,TOF],{'AMFM','FSST','TOF'}, ...
%     'Location','northoutside','Orientation','horizontal');
% % title('Mean SWS')
% ax = gca; ax.FontSize = 28;
% ax.LineWidth = 1.5;
% ax.Position = [0.2 0.2 0.7 0.7];
% 
% %% Bias (valor absoluto)
% figure('Position',[200 200 400 400])
% ylim([-16 16]), xlim([350 850]);
% grid on, hold on
% yline(0,'k--', 'LineWidth',3)
% p_AMFM = plot(freq,Info_AMFM.bias,'o-','MarkerSize',5,'LineWidth',3,"MarkerFaceColor","auto");
% p_FSST = plot(freq,Info_FSST.bias,'o-','MarkerSize',5,'LineWidth',3,"MarkerFaceColor","auto");
% hold off
% xlabel('Frequency [Hz]'), ylabel('Bias [%]')
% legend([p_AMFM,p_FSST],{'AMFM','FSST'},...
%     'Location','northoutside','Orientation','horizontal');
% ax = gca; ax.FontSize = 28;
% ax.LineWidth = 1.5;
% ax.Position = [0.2 0.2 0.7 0.7];
% %title('Bias')
% 
% %% CV
% figure('Position',[200 200 400 400])
% ylim([0 30]), xlim([350 850]);
% grid on, hold on
% p_AMFM = plot(freq,Info_AMFM.CV,'o-','MarkerSize',5,'LineWidth',3,"MarkerFaceColor","auto");
% p_FSST = plot(freq,Info_FSST.CV,'o-','MarkerSize',5,'LineWidth',3, "MarkerFaceColor","auto");
% 
% hold off
% xlabel('Frequency [Hz]'), ylabel('CV [%]')
% legend([p_AMFM,p_FSST],{'AMFM','FSST'},...
%     'Location','northoutside','Orientation','horizontal');
% ax = gca; ax.FontSize = 28;
% ax.LineWidth = 1.5;
% ax.Position = [0.2 0.2 0.7 0.7];
% %title('Coefficient Variation')