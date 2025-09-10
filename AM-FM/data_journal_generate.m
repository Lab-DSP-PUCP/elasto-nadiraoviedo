% Generates images for the three methods
close all;
clear, clc;
load('MyColormaps.mat');
%% Initialization
SWS_range = [1.5 8.5];
SWS_range_normative = [1 9.5];
SWS_im_range = [1.5 8.5];
% SWS_window = [45 15];
SWS_window = [65 16]; %[60 16] originalmente
wind_size =0;
freq = 400:100:800;
id = [2 4 8];

% BaseDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...nmaa
%     'Elastography\heterogeneo_sono'];
BaseDir2 = ['D:\Universidah\AM-FM_demod\SWS_amfm_stimator\AM_FM\comparacion8\Paralelo\MatlabProcessed'];
BaseDir = 'D:\Universidah\AM-FM_demod\SWS_amfm_stimator\AM_FM';

%%
for iq = 8 
    %% Data and image
    for id = 2:2:4
        %FOR HOMOGENEOUS
        %--------------------------------------------------------------------------------------------------------------
        iq =4;
        id =2;
        %--------------------------------------------------------------------------------------------------------------
        % BaseDir2 = ['/home/dsplab/Desktop/Joaquin/2024-2/Tesis/Tesis/data_journal/',num2str(iq),'-',num2str(id)];
        % swsDir = [BaseDir,'/sws_two_layers/',num2str(iq),'-',num2str(id);
        swsDir = [BaseDir,'/sws_two_layers/',num2str(iq),'-',num2str(id)];
        %% all frequencies
        for is = 1:5
            %is=3;
            BaseDir2 = ['D:\Universidah\AM-FM_demod\SWS_amfm_stimator\AM_FM\dataset_rf_iq_data\data_f_',num2str(freq(is)),'_sws_',num2str(iq),'.00_',num2str(id),'.00.mat'];
            % directory = [BaseDir2,'/data_f_',int2str(freq(is)),'_sws_',num2str(iq),'.00_',num2str(id),'.00.mat'];
            
            load(BaseDir2);
            %%
            IQ = matrix_IQ;
            
            IQ_dims = size(IQ); % IQ-data ---> 3D matrix of size [M,N,P]
             % if  IQ_dims(2) == 128
             %    IQ = IQ(:,33:96,:); 
             % end
            N_pv = 48; % # PV frames 48!!
            PW_ens = 2; % ensemble length (# PW/ensemble), PW: Plane Waves
            % % % % N_pv = 6; % # PV frames 48!!
            % % % % PW_ens = 2; % ensemble length (# PW/ensemble), PW: Plane Waves
            % PW_ens = 16; % ensemble length (# PW/ensemble), PW: Plane Waves
            N_frames = N_pv + PW_ens; % # IQ frames required
            %57->47
            IQ = IQ(:,:,1:N_frames);
            % IQ = IQ(:,:,26:38);
            N_angles = 1;
            [v,dinf] = pv_cal(IQ,dinf,N_angles,PW_ens);
            v_dim = size(v); % pv-video ---> 3D matrix of size [M,N,P]
            % lateral-dimension vector
            % x = ((0:v_dim(2) - 1) - v_dim(2)/2)*dinf.dx;
            x = x_index;
            % x_size = size(x);
            % if x_size(2) == 128
            %     x = x(:,33:96);
            % end
            dinf.x = x;
            % axial-dimension vector
            % z = (0:v_dim(1)- 1 )*dinf.dz;
            z = z_index_new;
            dinf.z=z;
            PRF = dinf.PRF;
            frame = 10;
            B_mode = db(IQ);
            B_mode = B_mode - max(B_mode(:));
            dinf.B_mode=B_mode(:,:,frame);
            v_abs = abs(v); % absolute value is taken since we're interested in the peak value
            N_avg = 4;
            sono_frames = 10; %10!!!!!
            
            for i = 1:sono_frames
                sono_video = sum(v_abs(:,:,i:i+N_avg-1),3)*(1/N_avg);
            end
            
            sono = sono_video(:,:,1);
            % sono_window = [7,2];
            sono_window = [4 2];%[6 2][9,3]; 
            [sono_norm,~,~] = normalize(sono,2); % depth normalization (for each row)
            sono_filt = medfilt2(sono_norm, sono_window); % median filtering 
            sono_filt_mov = sono_filt;
           % sono_filt_mov = BM3D(sono_filt,0.1);
            % sono_filt_mov = filter_bandpass(sono_filt,SWS_range_normative,dinf);
            % [sono_filt_mov,~,~] = normalize(sono_filt_mov,2); % depth normalization (for each row)
            % sono_filt_mov = imgaussfilt(sono_filt);
            % [gradThresh,numIter]=imdiffuseest(sono_filt_mov,'ConductionMethod','exponential','Connectivity','maximal');
            % sono_filt_diff = imdiffusefilt(sono_filt_mov,'GradientThreshold', ...
            %    gradThresh,'NumberOfIterations',numIter);
            % sono_filt_mov = sono_filt_diff;

            % figure; set(gcf,'Position',[200 200 600 250]);  
            % t_fig = tiledlayout(1,3);
            % nexttile; imagesc(x*1e3,z*1e3,sono); set(gcf,'Colormap',sonomap); title('Raw Sonoelastogram');
            % % nexttile; imagesc(x*1e3,z*1e3,sono_norm); set(gcf,'Colormap',sonomap); title(['Normalized']);
            % nexttile; imagesc(x*1e3,z*1e3,sono_filt); set(gcf,'Colormap',sonomap); title(['Medfilt']);
            % nexttile; imagesc(x*1e3,z*1e3,sono_filt_mov); set(gcf,'Colormap',sonomap); title('Filtered');

            figure, imagesc(x*1e3,z*1e3,sono_filt_mov); set(gcf,'Colormap',sonomap);
            xlabel('Width [mm]','fontsize',28);ylabel('Depth [mm]','fontsize',28)
            ax = gca; ax.FontSize = 28;
            axis image

            title({'Sonoelastogram',['f_{vib} = ',num2str(dinf.f_vib), 'Hz']}); xlabel('Lateral Distance [mm]'); ylabel('Depth [mm]');

        
         

           %% AM-FM
            tic
            [vshearsin,vshearcos] = AMFM_demod_StWS(sono_filt_mov,dinf);
            toc
       
            vshearsin_filt = medfilt2(vshearsin,SWS_window,'symmetric');
            vshearsin_im = mean(vshearsin_filt,3);
        
            SWS_im_range = [1.5 8.5];
            figure,imagesc(x*1e3,z*1e3,vshearsin_im,SWS_im_range)
            set (gca,'clim',SWS_im_range),
            h = colorbar;
            ylabel(h, 'SWS m/s','FontSize',45);
            xlabel('Width [mm]','fontsize',45);ylabel('Depth [mm]','fontsize',45)
            ax = gca; ax.FontSize = 45;
            title(['SWS AMFM Vib Freq = ' num2str(dinf.f_vib) ' Hz a ' num2str(iq),'.0-',num2str(id),'.0 m/s'],'fontsize',14)
            % title(['SWS FSST Vib Freq = ' num2str(dinf.f_vib) ' Hz a ' num2str(iq),'.0 m/s'],'fontsize',14)
            colormap turbo
             % 
            hold on
            % x_left = [-8 -2];
            x_left = [-14 -5];
            z_left = [18 42];
            plot([x_left(1),x_left(2),x_left(2),x_left(1),x_left(1)],...
                [z_left(1),z_left(1),z_left(2),z_left(2),z_left(1)],'w--','LineWidth',4)

            % x_right = [2 8];
            x_right = [5 14];
            z_right = [18 42];
            plot([x_right(1),x_right(2),x_right(2),x_right(1),x_right(1)],...
                [z_right(1),z_right(1),z_right(2),z_right(2),z_right(1)],'w--','LineWidth',4)

            % %x_ROI = [-8 8];
            % x_ROI = [-12 12];
            % % z_ROI = [37 63];
            % z_ROI = [18 42];
            % plot([x_ROI(1),x_ROI(2),x_ROI(2),x_ROI(1),x_ROI(1)],...
            %     [z_ROI(1),z_ROI(1),z_ROI(2),z_ROI(2),z_ROI(1)],'w--','LineWidth',4)
            hold off

      
            %% FSST
            wind_size = 0;
            tic
            SWS_FSST = process_FSST_loupas(sono_filt_mov,dinf,SWS_range,wind_size);
            toc
            SWS_FSST_filt = medfilt2(SWS_FSST,SWS_window,'symmetric');
            SWS_FSST_im = mean(SWS_FSST_filt,3);
    
            SWS_im_range = [1.5 8.5];
            figure,imagesc(x*1e3,z*1e3,SWS_FSST_im,SWS_im_range)
            set (gca,'clim',SWS_im_range),
            h = colorbar;
            ylabel(h, 'SWS m/s','FontSize',45);
            xlabel('Width [mm]','fontsize',45);ylabel('Depth [mm]','fontsize',45)
            ax = gca; ax.FontSize = 45;
            title(['SWS FSST Vib Freq = ' num2str(dinf.f_vib) ' Hz a ' num2str(iq),'.0-',num2str(id),'.0 m/s'],'fontsize',14)
            % title(['SWS FSST Vib Freq = ' num2str(dinf.f_vib) ' Hz a ' num2str(iq),'.0 m/s'],'fontsize',14)
            colormap turbo
             % 
            hold on
            % x_left = [-8 -2];
            x_left = [-14 -5];
            z_left = [18 42];
            plot([x_left(1),x_left(2),x_left(2),x_left(1),x_left(1)],...
                [z_left(1),z_left(1),z_left(2),z_left(2),z_left(1)],'w--','LineWidth',4)

            % x_right = [2 8];
            x_right = [5 14];
            z_right = [18 42];
            plot([x_right(1),x_right(2),x_right(2),x_right(1),x_right(1)],...
                [z_right(1),z_right(1),z_right(2),z_right(2),z_right(1)],'w--','LineWidth',4)

            % %x_ROI = [-8 8];
            % x_ROI = [-12 12];
            % % z_ROI = [37 63];
            % z_ROI = [18 42];
            % plot([x_ROI(1),x_ROI(2),x_ROI(2),x_ROI(1),x_ROI(1)],...
            %     [z_ROI(1),z_ROI(1),z_ROI(2),z_ROI(2),z_ROI(1)],'w--','LineWidth',4)
            hold off
    
    
            %%
            save([swsDir,'/',int2str(freq(is)),'.mat'],'vshearsin','vshearsin_im',...
                'SWS_FSST','SWS_FSST_im','dinf','wind_size');
        end
    end
end

%% B-MODE HETEROGENEOUS

x_0 = 1000*dinf.x;
z_0 = 1000*dinf.z;
[X,Z] = meshgrid(x_0,z_0);
% ROI (left)
x_left = [-14 -5];
z_left = [19 41];
ROI_left = x_left(1)<X & X<x_left(2) & z_left(1)<Z & Z<z_left(2);

% ROI (right)   3.5 mm from inc
x_right = [5 14];
z_right = [19  41];
ROI_right = x_right(1)<X & X<x_right(2) & z_right(1)<Z & Z<z_right(2);

image = dinf.B_mode;
figure('Position', [200 200 350 350]),
imagesc(x_0,z_0,image, [-60 0])
axis equal
xlim([x_0(1), x_0(end)]); xlabel('Width [mm]')
ylim([z_0(1), z_0(end)]); ylabel('Depth [mm]')
%ylabel(h, 'SWS m/s','FontSize',14);|
%xlabel('Width [mm]','fontsize',14);ylabel('Depth [mm]','fontsize',14)
%title(['SWS Vib Freq = ' num2str(Properties.VibFreq) ' Hz'],'fontsize',14)
colormap gray, 
c = colorbar; c.Label.String = 'dB';
ax = gca; ax.FontSize = 28;
hold on
plot([x_left(1),x_left(2),x_left(2),x_left(1),x_left(1)],...
    [z_left(1),z_left(1),z_left(2),z_left(2),z_left(1)],'r--','LineWidth',7)
plot([x_right(1),x_right(2),x_right(2),x_right(1),x_right(1)],...
    [z_right(1),z_right(1),z_right(2),z_right(2),z_right(1)],'r--','LineWidth',7)
hold off

%% B.MODE HOMOGENEOUS

x_0 = 1000*dinf.x;
z_0 = 1000*dinf.z;
[X,Z] = meshgrid(x_0,z_0);
% ROI 
x_ROI = [-14 14];
z_ROI = [19 41];
ROI = x(1)<X & X<x(2) & z(1)<Z & Z<z(2);

image = dinf.B_mode;
figure('Position', [200 200 350 350]),
imagesc(x_0,z_0,image, [-60 0])
axis equal
xlim([x_0(1), x_0(end)]); xlabel('Width [mm]')
ylim([z_0(1), z_0(end)]); ylabel('Depth [mm]')
%ylabel(h, 'SWS m/s','FontSize',14);
%xlabel('Width [mm]','fontsize',14);ylabel('Depth [mm]','fontsize',14)
%title(['SWS Vib Freq = ' num2str(Properties.VibFreq) ' Hz'],'fontsize',14)
colormap gray, 
c = colorbar; c.Label.String = 'dB';
ax = gca; ax.FontSize = 28;
hold on
plot([x_ROI(1),x_ROI(2),x_ROI(2),x_ROI(1),x_ROI(1)],...
                [z_ROI(1),z_ROI(1),z_ROI(2),z_ROI(2),z_ROI(1)],'r--','LineWidth',7)
hold off

%% SWS OVERLAY BMODE HETEROGENEOUS

figure('Position', [200 200 320 320]),
imoverlay2(dinf.B_mode,vshearsin_im,[-60 0],SWS_im_range,0.5,x_0,z_0,1);
hold on;
plot([x_left(1),x_left(2),x_left(2),x_left(1),x_left(1)],...
    [z_left(1),z_left(1),z_left(2),z_left(2),z_left(1)],'w--','LineWidth',5)
plot([x_right(1),x_right(2),x_right(2),x_right(1),x_right(1)],...
    [z_right(1),z_right(1),z_right(2),z_right(2),z_right(1)],'w--','LineWidth',5)
hold off;
xlabel('Width [mm]'), ylabel('Depth [mm]')
c = colorbar; c.Label.String = 'm/s';c.Label.FontSize = 28;
ax = gca; ax.FontSize = 28;

%% SWS OVERLAY BMODE HOMOGENEOUS

figure('Position', [200 200 320 320]),
imoverlay2(dinf.B_mode,vshearsin_im,[-60 0],SWS_im_range,0.5,x_0,z_0,1);
hold on;
plot([x_ROI(1),x_ROI(2),x_ROI(2),x_ROI(1),x_ROI(1)],...
    [z_ROI(1),z_ROI(1),z_ROI(2),z_ROI(2),z_ROI(1)],'w--','LineWidth',5)
hold off;
xlabel('Width [mm]'), ylabel('Depth [mm]')
c = colorbar; c.Label.String = 'm/s';c.Label.FontSize = 28;
ax = gca; ax.FontSize = 28;