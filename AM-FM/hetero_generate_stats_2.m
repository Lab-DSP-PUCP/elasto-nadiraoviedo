%% Initialization
Info_inc_AMFM = array2table(zeros(9,5),'VariableNames',{'SWS','std','CV','bias','bias_er'});
Info_back_AMFM = array2table(zeros(9,5),'VariableNames',{'SWS','std','CV','bias','bias_er'});
Info_inc_FSST = array2table(zeros(9,5),'VariableNames',{'SWS','std','CV','bias','bias_er'});
Info_back_FSST = array2table(zeros(9,5),'VariableNames',{'SWS','std','CV','bias','bias_er'});
resAMFM = zeros(9,1);
resFSST = zeros(9,1);

rd1 = [2,5];
% BaseDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
%     'Elastography\heterogeneo_sono'];
BaseDir = 'D:\Universidah\AM-FM_demod\SWS_amfm_stimator\AM_FM';
BaseDir2 = 'D:\Universidah\AM-FM_demod\SWS_amfm_stimator\AM_FM\comparacion8\Paralelo\MatlabProcessed';
swsDir = [BaseDir,'/sws'];

%% Loop for every frequency
%%
for id = 1:9
    %% loading relevant info
    load([BaseDir2,'/Image',num2str(id),'/sono.mat']);
    load([swsDir,'/',num2str(id),'.mat'])

    %% Generate ROI mask
    % Here we can plot the ROI and calculate parameters
    x = 1000*Properties.Width_S;
    z = 1000*Properties.Depth_S;
    [X,Z] = meshgrid(x,z);
    % lado = 11    C = (20.7,15.6)

    % ROI (inc)
    % x_inc = [15.2 26.2];
    % z_inc = [10.1 20.1];

    x_inc = [15.5 25.5];
    z_inc = [9.5 20.5];
    ROI_inc = x_inc(1)<X & X<x_inc(2) & z_inc(1)<Z & Z<z_inc(2);

    % ROI (back)   3.5 mm from inc
    % x_back = [6.2 9.7 31.7 35.2];

    x_back = [6 10.2 31 35.2];
    z_back = z_inc;
    ROI_back = (( x_back(1)<X & X<x_back(2) ) | ( x_back(3)<X & X<x_back(4) ))...
        & z_back(1)<Z & Z<z_back(2);

    %% Calculate bias and CV
    % SWS_CWT_im = mean(SWS_CWT,3);
    %SWS_CWT_im = SWS_CWT_unfilt(:,:,1);
    % SWS_FSST_im = mean(SWS_FSST,3);
    %SWS_STFT_im = SWS_STFT_unfilt(:,:,1);
    Info_inc_AMFM(id,:) = calc_param(vshearsin_im,ROI_inc,5.1);
    Info_back_AMFM(id,:) = calc_param(vshearsin_im,ROI_back,3.45);
    Info_inc_FSST(id,:) = calc_param(SWS_FSST_im,ROI_inc,5.1);
    Info_back_FSST(id,:) = calc_param(SWS_FSST_im,ROI_back,3.45);


    % resCWT(id) = getFWHM(SWS_CWT_im,x);
    % resSTFT(id) = getFWHM(SWS_STFT_im,x); 
    % resPD(id) = getFWHM(SWS_PD,x);
    % resRW(id) = getFWHM(SWS_RW,x);
    resAMFM(id) = getFWHM(vshearsin_im,x,...
        Info_back_AMFM(id,:).SWS, Info_inc_AMFM(id,:).SWS);
    resFSST(id) = getFWHM(SWS_FSST_im,x,...
        Info_back_FSST(id,:).SWS, Info_inc_FSST(id,:).SWS);


end
save('stats2.mat','Info_inc_FSST','Info_inc_AMFM',...
    'Info_back_FSST','Info_back_AMFM',...
    "resFSST","resAMFM");
%save('stats2.mat','Info_inc_FSST','Info_inc_CWT',...
 %   'Info_back_FSST','Info_back_CWT',...
  %  "resFSST","resCWT");

%% Plotting ROI
% Here we can plot the ROI and calculate parameters
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
image = Properties.Bmode;
figure('Position', [200 200 350 350]),
imagesc(x,z,image, [-70 0])
axis equal
xlim([x(1), x(end)]); xlabel('Width [mm]')
ylim([z(1), z(end)]); %ylabel('Depth [mm]')
%ylabel(h, 'SWS m/s','FontSize',14);
%xlabel('Width [mm]','fontsize',14);ylabel('Depth [mm]','fontsize',14)
%title(['SWS Vib Freq = ' num2str(Properties.VibFreq) ' Hz'],'fontsize',14)
colormap gray, 
c = colorbar; c.Label.String = 'dB';
hold on
plot([x_inc(1),x_inc(2),x_inc(2),x_inc(1),x_inc(1)],...
    [z_inc(1),z_inc(1),z_inc(2),z_inc(2),z_inc(1)],'w--','LineWidth',1.5)
plot([x_back(1),x_back(2),x_back(2),x_back(1),x_back(1)],...
    [z_back(1),z_back(1),z_back(2),z_back(2),z_back(1)],'w--','LineWidth',1.5)
plot([x_back(3),x_back(4),x_back(4),x_back(3),x_back(3)],...
    [z_back(1),z_back(1),z_back(2),z_back(2),z_back(1)],'w--','LineWidth',1.5)
hold off

set(gca,'FontSize',25);

%%

SWS_im_range = [2,8];
figure('Position', [200 200 320 320]),
% t_fig = tiledlayout(1,2);
% nexttile;
imoverlay2(Properties.Bmode,SWS_FSST_im,[-70 0],SWS_im_range,0.5,x,z,ROI_inc);
hold on;
plot([x_inc(1),x_inc(2),x_inc(2),x_inc(1),x_inc(1)],...
    [z_inc(1),z_inc(1),z_inc(2),z_inc(2),z_inc(1)],'w--','LineWidth',4)
plot([x_back(1),x_back(2),x_back(2),x_back(1),x_back(1)],...
    [z_back(1),z_back(1),z_back(2),z_back(2),z_back(1)],'w--','LineWidth',4)
plot([x_back(3),x_back(4),x_back(4),x_back(3),x_back(3)],...
    [z_back(1),z_back(1),z_back(2),z_back(2),z_back(1)],'w--','LineWidth',4)
hold off;
xlabel('Width [mm]'), ylabel('Depth [mm]')
ax = gca; ax.FontSize = 25;

%%
nexttile;
imoverlay2(Properties.Bmode,vshearsin_im,[-70 0],SWS_im_range,0.5,x,z,ROI_inc);
hold on;
plot([x_inc(1),x_inc(2),x_inc(2),x_inc(1),x_inc(1)],...
    [z_inc(1),z_inc(1),z_inc(2),z_inc(2),z_inc(1)],'w--','LineWidth',4)
plot([x_back(1),x_back(2),x_back(2),x_back(1),x_back(1)],...
    [z_back(1),z_back(1),z_back(2),z_back(2),z_back(1)],'w--','LineWidth',4)
plot([x_back(3),x_back(4),x_back(4),x_back(3),x_back(3)],...
    [z_back(1),z_back(1),z_back(2),z_back(2),z_back(1)],'w--','LineWidth',4)
hold off;
xlabel('Width [mm]'), %ylabel('Depth [mm]')
colormap turbo
c = colorbar; c.Label.String = 'm/s';c.Label.FontSize = 30;
ax = gca; ax.FontSize = 30;


%%
function Info = calc_param(image,mask,real_SWS)
    SWS = mean(image(mask));
    SD = std(image(mask));
    CV = SD/SWS * 100;
    bias = (SWS-real_SWS)/real_SWS *100;
    Info = array2table([SWS,SD,CV,bias,0]);
end

function out = getFWHM(SWS,x, c_back, c_inc)
% SWS = medfilt2(SWS,[20,5],'symmetric');
iz0 = 149; izf = 169;
swsRegion = SWS(iz0:izf,:);
meanProfile = mean(swsRegion);

% ix0Up = 1; ixfUp = 80;
% ix0Down = 48; ixfDown = 128; 
% % figure,
% % plot(x,meanProfile),
% % xline(x(ix0Up),'b')
% % xline(x(ixfUp),'b')
% % xline(x(ixfDown),'r')
% % xline(x(ix0Down),'r')
% % axis tight
% swsDiff = diff(meanProfile(ix0Up:ixfUp));
% xDiff = x(ix0Up:ixfUp-1);
% fwhmUp = fwhm(xDiff,swsDiff);
% 
% swsDiff = -diff(meanProfile(ix0Down:ixfDown));
% xDiff = x(ix0Down:ixfDown-1);
% fwhmDown = fwhm(xDiff,swsDiff);
% 
% out = (fwhmUp + fwhmDown)/2;

[fitresult, ~] = createFitSws(x, meanProfile, c_back, c_inc);
out = (fitresult.lambda1 + fitresult.lambda2) *log(4);
end