clear;
clc;

close all;

load('MyColormaps.mat');
%%
move='right';
setup='parallel';
rd=[1.5,9];
BaseDir2 = ['D:\Universidah\AM-FM_demod\SWS_amfm_stimator\AM_FM\comparacion8\Paralelo\MatlabProcessed'];
BaseDir = 'D:\Universidah\AM-FM_demod\SWS_amfm_stimator\AM_FM';
swsDir = [BaseDir,'\sws'];
rd1 = [1.5 9];  

%%

for id = 1:9
%% Normalizado x vs tiempo
    directory = [BaseDir2,'\Image',num2str(id),'\sono.mat'];
    load(directory);
    % load('sono.mat')
    Properties.dx=3.08e-4;
    Properties.pitch=3.08e-04;

    [sono_filt_mov,sono_filt,mask]=process_sono_data(sono,Properties,move,rd);
    
    % Restar el componente DC de cada traza x(t)
    sono_filt_mov = sono_filt_mov - mean(sono_filt_mov, 3);

%% AM-FM
    tic
    [vshearsin,vshearcos] = AMFM_demod(sono_filt_mov,Properties,0);
    toc

    figure;
    % vshearsin=medfilt2(vshearsin,[18 6]);
    % vshearsin_im = mean (vshearsin,3);
    
    vshearcos=medfilt2(vshearcos,[18 6]);
    vshearsin_im = mean (vshearcos,3);

    imagesc( 10^3*Properties.Width_S,10^3*Properties.Depth_S,vshearsin_im);
    h = colorbar;
    ylabel(h, 'SWS m/s','FontSize',14);
    xlabel('Width [mm]','fontsize',14);ylabel('Depth [mm]','fontsize',14)
    title(['SWS AM-FM Vib Freq = ' num2str(Properties.VibFreq) ' Hz'],'fontsize',14)
    colormap turbo;
    set (gca,'clim',rd1);


    hold on
    [X,Z] = meshgrid(1000*Properties.Width_S,1000*Properties.Depth_S);
    L = 10; C = [20.5,15.6]; sep = 4;

    % ROI (inc)
    x_inc = [C(1)-L/2 C(1)+L/2];
    z_inc = [C(2)-L/2 C(2)+L/2];

    % ROI (back)   3.5 mm from inc
    x_back = [x_inc(1)-sep-L/2 x_inc(1)-sep x_inc(2)+sep x_inc(2)+sep+L/2];
    z_back = z_inc;


    plot([x_inc(1),x_inc(2),x_inc(2),x_inc(1),x_inc(1)],...
        [z_inc(1),z_inc(1),z_inc(2),z_inc(2),z_inc(1)],'w--','LineWidth',1.5)
    plot([x_back(1),x_back(2),x_back(2),x_back(1),x_back(1)],...
        [z_back(1),z_back(1),z_back(2),z_back(2),z_back(1)],'k--','LineWidth',1.5)
    plot([x_back(3),x_back(4),x_back(4),x_back(3),x_back(3)],...
        [z_back(1),z_back(1),z_back(2),z_back(2),z_back(1)],'k--','LineWidth',1.5)
    hold ofF

    % subplot(211);
    %% fSST
    tic
    SWS_FSST = process_FSST(sono_filt_mov,Properties,rd,12);
    % SWS_FSST = process_FSST_normal(sono_filt_mov,Properties,rd,1,a1,a2);
    toc
    % SWS_STFT_filt = zeros(size(SWS_STFT));
    % for t = 1:Properties.nframes
    %     SWS_STFT_filt(:,:,t) = medfilt2(SWS_STFT(:,:,t),[20,5],'symmetric');
    % end
    % SWS_FSST_im = mean(SWS_STFT_filt,3);
    SWS_FSST_im = mean(SWS_FSST,3);

    figure,imagesc(10^3*Properties.Width_S,10^3*Properties.Depth_S,SWS_FSST_im)
    set (gca,'clim',rd1),
    h = colorbar;
    ylabel(h, 'SWS m/s','FontSize',14);
    xlabel('Width [mm]','fontsize',14);ylabel('Depth [mm]','fontsize',14)
    title(['SWS FSST Vib Freq = ' num2str(Properties.VibFreq) ' Hz'],'fontsize',14)
    colormap turbo

    hold on
    [X,Z] = meshgrid(1000*Properties.Width_S,1000*Properties.Depth_S);
    L = 10; C = [20.5,15.6]; sep = 4;

    % ROI (inc)
    x_inc = [C(1)-L/2 C(1)+L/2];
    z_inc = [C(2)-L/2 C(2)+L/2];

    % ROI (back)   3.5 mm from inc
    x_back = [x_inc(1)-sep-L/2 x_inc(1)-sep x_inc(2)+sep x_inc(2)+sep+L/2];
    z_back = z_inc;


    plot([x_inc(1),x_inc(2),x_inc(2),x_inc(1),x_inc(1)],...
        [z_inc(1),z_inc(1),z_inc(2),z_inc(2),z_inc(1)],'w--','LineWidth',1.5)
    plot([x_back(1),x_back(2),x_back(2),x_back(1),x_back(1)],...
        [z_back(1),z_back(1),z_back(2),z_back(2),z_back(1)],'k--','LineWidth',1.5)
    plot([x_back(3),x_back(4),x_back(4),x_back(3),x_back(3)],...
        [z_back(1),z_back(1),z_back(2),z_back(2),z_back(1)],'k--','LineWidth',1.5)
    hold off
    
    % subplot(212);imagesc(vshearsin);

%%
    save([swsDir,'\',num2str(id),'.mat'],'SWS_FSST','SWS_FSST_im',...
        'vshearsin','vshearsin_im','Properties')
end

%%
[X,Z] = meshgrid(1000*Properties.Width_S,1000*Properties.Depth_S);
L = 10; C = [20.5,15.6]; sep = 4;

% ROI (inc)
% x_inc = [C(1)-L/2 C(1)+L/2];
% z_inc = [C(2)-L/2 C(2)+L/2];

x_inc = [15.5 25.5];
z_inc = [9.5 20.5];
ROI_inc = x_inc(1)<X & X<x_inc(2) & z_inc(1)<Z & Z<z_inc(2);

% ROI (back)   3.5 mm from inc
% x_back = [x_inc(1)-sep-L/2 x_inc(1)-sep x_inc(2)+sep x_inc(2)+sep+L/2];
% z_back = z_inc;

x_back = [6 10.2 31 35.2];
z_back = z_inc;
ROI_back = (( x_back(1)<X & X<x_back(2) ) | ( x_back(3)<X & X<x_back(4) ))...
    & z_back(1)<Z & Z<z_back(2);
x = 10^3*Properties.Width_S;
z = 10^3*Properties.Depth_S;

SWS_im_range = [2,6];
figure('Position', [200 200 320 320]),
imoverlay2(Properties.Bmode,vshearsin_im,[-70 0],SWS_im_range,0.5,x,z,ROI_inc);
hold on;
plot([x_inc(1),x_inc(2),x_inc(2),x_inc(1),x_inc(1)],...
    [z_inc(1),z_inc(1),z_inc(2),z_inc(2),z_inc(1)],'w--','LineWidth',4)
plot([x_back(1),x_back(2),x_back(2),x_back(1),x_back(1)],...
    [z_back(1),z_back(1),z_back(2),z_back(2),z_back(1)],'w--','LineWidth',4)
plot([x_back(3),x_back(4),x_back(4),x_back(3),x_back(3)],...
    [z_back(1),z_back(1),z_back(2),z_back(2),z_back(1)],'w--','LineWidth',4)
hold off;
xlabel('Width [mm]'), ylabel('Depth [mm]')
ax = gca; ax.FontSize = 30;