close all;
clear, clc;
%% Initialization
rd=[1.5,8.5];
freq = 400:100:800;
BaseDir = 'D:\Universidah\AM-FM_demod\SWS_amfm_stimator\AM_FM';
%%
for iq = 2
    %% Data and image
        for id = 2%:2:4
            iq = 8;
            id = 2;
            Info_left_AMFM = array2table(zeros(5,5),'VariableNames',{'SWS','std','CV','bias','bias_er'});
            Info_right_AMFM = array2table(zeros(5,5),'VariableNames',{'SWS','std','CV','bias','bias_er'});
            Info_left_FSST = array2table(zeros(5,5),'VariableNames',{'SWS','std','CV','bias','bias_er'});
            Info_right_FSST = array2table(zeros(5,5),'VariableNames',{'SWS','std','CV','bias','bias_er'});

            resAMFM = zeros(5,1);
            resFSST = zeros(5,1);

            for is = 1:5
                load([BaseDir,'/sws_two_layers/',num2str(iq),'-',num2str(id),'/',int2str(freq(is)),'.mat']);
                freq = 400:100:800;
                 %% Generate ROI mask
                % Here we can plot the ROI and calculate parameters
                x_0 = 1000*dinf.x;
                z_0 = 1000*dinf.z;
                [X,Z] = meshgrid(x_0,z_0);
                % ROI (inc)
                % x_left = [-8 -2];
                % x_right = [2 8];
                x_left = [-14 -5];
                x_right = [5 14];
                % z_left = [35 65];
                z_left = [19 41];%[18 42];
                ROI_left = x_left(1)<X & X<x_left(2) & z_left(1)<Z & Z<z_left(2);
        
                % ROI (back)  
                z_right = z_left;
                ROI_right   =x_right(1)<X & X<x_right(2) & z_right(1)<Z & Z<z_right(2);
                %% Calculate bias and CV
                
                Info_left_AMFM(is,:) = calc_param(vshearsin_im,ROI_left,iq);
                Info_right_AMFM(is,:) = calc_param(vshearsin_im,ROI_right,id);

                Info_left_FSST(is,:) = calc_param(SWS_FSST_im,ROI_left,iq);
                Info_right_FSST(is,:) = calc_param(SWS_FSST_im,ROI_right,id);

                resAMFM(is) = getFWHM(vshearsin_im,x_0,...
                    Info_left_AMFM(is,:).SWS, Info_right_AMFM(is,:).SWS);
                resFSST(is) = getFWHM(SWS_FSST_im,x_0,...
                    Info_left_FSST(is,:).SWS, Info_right_FSST(is,:).SWS);

            end
            SaveDir = ['stats ',num2str(iq),'-',num2str(id),'.mat'];
            save(SaveDir,'Info_left_AMFM','Info_left_FSST',...
                'Info_right_AMFM','Info_right_FSST',...
                'resFSST','resAMFM');
        end
end

%%
for iq = 4
    %% Data and image
        for id = 4
            Info_AMFM = array2table(zeros(5,5),'VariableNames',{'SWS','std','CV','bias','bias_er'});
            Info_FSST = array2table(zeros(5,5),'VariableNames',{'SWS','std','CV','bias','bias_er'});
            for is = 1:5
                load([BaseDir,'/sws_two_layers/',num2str(iq),'-',num2str(id),'/',int2str(freq(is)),'.mat']);
                freq = 400:100:800;
                 %% Generate ROI mask
                % Here we can plot the ROI and calculate parameters
                x_0 = 1000*dinf.x;
                z_0 = 1000*dinf.z;
                [X,Z] = meshgrid(x_0,z_0);
                % lado = 11    C = (20.7,15.6)
                % ROI (inc)
                % x = [-8 8];
                x = [-14 14];
                %x_inc = [-7.8 6.6];
                %7-5=2,21-5=16
                % z = [35 65];
                z = [19 41];%[18 42];
                ROI = x(1)<X & X<x(2) & z(1)<Z & Z<z(2);
        
                %% Calculate bias and CV
                
                Info_AMFM(is,:) = calc_param(vshearsin_im,ROI,iq);
                Info_FSST(is,:) = calc_param(SWS_FSST_im,ROI,iq);
            
            end
            SaveDir = ['stats ',num2str(iq),'-',num2str(id),'.mat'];
            save(SaveDir,'Info_AMFM','Info_FSST');
        end
end

%%
function Info = calc_param(image,mask,real_SWS)
    SWS = mean(image(mask));
    SD = std(image(mask));
    CV = SD/SWS * 100;
    bias = (SWS-real_SWS)/real_SWS *100;
    bias_er = SD/real_SWS * 100;
    Info = array2table([SWS,SD,CV,bias,bias_er]);
end

function out = getFWHM(SWS,x, c_right, c_left)
% SWS = medfilt2(SWS,[20,5],'symmetric');
iz0 = 216; izf = 228;
swsRegion = SWS(iz0:izf,:);
meanProfile = mean(swsRegion);

[fitresult, ~] = createFit1(x, meanProfile,c_right,c_left);
out = (fitresult.lambda1) *log(4);
end