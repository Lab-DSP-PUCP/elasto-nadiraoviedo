clear;
clc;
close all;

load('MyColormaps.mat');
%%
move='right';
setup='parallel';
rd=[1.5,9];
BaseDir2 = 'D:\Universidah\AM-FM_demod\SWS_amfm_stimator\AM_FM\comparacion8\Paralelo\MatlabProcessed';
BaseDir = 'D:\Universidah\AM-FM_demod\SWS_amfm_stimator\AM_FM';
swsDir = [BaseDir,'\sws'];
    

%%
%tic
for id = 8
    
%% Normalizado x vs tiempo
    
    directory = [BaseDir2,'\Image',num2str(id),'\sono.mat'];
    load(directory);
    % load('sono.mat')
    Properties.dx=3.08e-4;
    Properties.pitch=3.08e-04;

    [sono_filt_mov,sono_filt,mask]=process_sono_data(sono,Properties,move,rd);
    
    % Restar el componente DC de cada traza x(t)
    sono_filt_mov = sono_filt_mov - mean(sono_filt_mov, 3);
    %slicen=180;
    
    figure;
    %subplot(211);imagesc(squeeze(sono_filt_mov(slicen,:,:)));title("x vs t pre")
    %subplot(212);
    imagesc( 10^3*Properties.Width_S,10^3*Properties.Depth_S,sono_filt_mov(:,:,10));title("frame pre");set(gcf,'colormap',sonomap);
    title('Interference pattern','FontSize',14);
    xlabel('Width [mm]','fontsize',14);ylabel('Depth [mm]','fontsize',14)
    %% Seleccion de senial 1D para analizar
    
    % vshearsin_pre = zeros(size(sono_filt_mov,1),size(sono_filt_mov,2),size(sono_filt_mov,3));
    % vshearcos_pre = zeros(size(sono_filt_mov,1),size(sono_filt_mov,2),size(sono_filt_mov,3));
    %tic
    vshearsin_pre = zeros(size(sono_filt_mov,1),size(sono_filt_mov,2));
    vshearcos_pre = zeros(size(sono_filt_mov,1),size(sono_filt_mov,2));


            %% Parámetros de los filtros Gabor
            numChannels = 4;                   % Número de canales
            r0 = 0.1;                            % Frecuencia base
            commonRatio = 1.5;                 % Razón geométrica
            octaveBandwidth = 1.5;             % Ancho de banda de 1.5 octavas
            filterSize = size(sono_filt_mov,2);         % Coincidir con la longitud de la señal
            filters = cell(1, numChannels);
            
            % Generar banco de filtros de Gabor
            xFilt = linspace(-1, 1, filterSize);
            for k = 1:numChannels
                f = r0 * commonRatio^(k - 1); % Frecuencia escalada geométricamente
                sigma = f / (2 * pi * sqrt(log(2) / 2) * (2^octaveBandwidth - 1));
                % Filtro Gabor
                %gabor = exp(-xFilt.^2 / (2 * sigma^2)) .* cos(2 * pi * f * xFilt);
                % Filtro Gabor COMPLEJO (correcto para análisis AM-FM)
                gabor = exp(-xFilt.^2 / (2 * sigma^2)) .* exp(1j * 2 * pi * f * xFilt);
            
                % Normalizar el filtro a norma L2
                gabor = gabor / norm(gabor);
                filters{k} = gabor;
            end

    tic
    for framen = 10%:size(sono_filt_mov,3)

        % figure;
        % subplot(211);imagesc(squeeze(sono_filt_mov(150,:,:)));title("x vs t pre")
        % subplot(212);imagesc(sono_filt_mov(:,:,10));title("frame pre");set(gcf,'colormap',sonomap);

        %% Seleccion de senial 1D para analizar

        for slice = 1:size(sono_filt_mov,1)
        
            vect = sono_filt_mov(slice,:,framen);
            vect_A = hilbert(vect); %Crear senial analitica
            
            x_esp = (0:size(sono_filt_mov,2)-1)*Properties.dx;
            % fs = 1/Properties.dx;
            % N = length(vect); % Número de puntos
            % f = (0:N-1)*(fs/N); % Vector de frecuencias
            % VECT= fftshift(abs(fft(vect))/N);
            % f_shifted= (-N/2:N/2-1)*(fs/N); % Reordenar el eje de frecuencias
            % f_norm =(f_shifted/fs);

            % 
            % figure;
            % title("Senial elegida y descompuesta");
            % plot(x_esp*10e2,vect); xlabel("Width [mm]"); ylabel("Displacement");
            % subplot(211);plot(x_esp*10e2,vect); xlabel("Time [s]"); ylabel("Amplitude");
            % subplot(212);stem(f_shifted,VECT);  xlabel("Frecuency [Hz]"); ylabel("Magnitude");
           
            % 

            
            
            %% Aplicar filtros Gabor a la señal
            filteredSignals = zeros(numChannels, length(vect_A));
            for k = 1:numChannels
                % Convolución con el filtro Gabor
                filteredSignals(k, :) = conv(vect_A, filters{k}, 'same');
            end
            
            
            %% Normalizacion (Segun paper) 
            %Relaciona la amplitud con la magnitud
            psi= zeros(numChannels, length(vect_A));
            valmax=0;
            for i = 1:size(filters,2)
                G = abs(fft(filters{i}, 10 * length(filters{i})));
                mx=max(abs(G));
                psi(i, :) = abs(filteredSignals(i, :))/mx;
                E = sum(abs(psi(i, :)).^2);
                %Hallar el canal con mayor energia
                if E > valmax
                    domchan=psi(i,:);
                    valmax = E;
                    index = i;
                end
            end
            
            
            %% QEA
            [grad_phi_sin, grad_phi_cos] = compute_phase_gradient(filteredSignals(index, :), 1);
            % figure;
            % subplot(211);plot(x_esp * 1e3, (pi*Properties.VibFreqOffset)./grad_phi_cos); xlabel("Lateral distance [mm]");ylabel("Vshear");
            % title(sprintf('Altura: %d Frame: %d', slicen, framen));
            % subplot(212);plot(x_esp * 1e3,grad_phi_sin); xlabel("Lateral distance [mm]");ylabel("QEA [mm-1]");
            
            % vshearsin_pre(slice,:,framen)=(pi*Properties.VibFreqOffset)./grad_phi_sin;
            % vshearcos_pre(slice,:,framen)=(pi*Properties.VibFreqOffset)./grad_phi_cos;

            vshearsin_pre(slice,:)=(pi*Properties.VibFreq)./(grad_phi_sin*10^3)*1.2;
            %vshearcos_pre(slice,:)=(pi*Properties.VibFreq)./(grad_phi_cos*10^3);

            vshearcos_pre(slice,:)=(pi*Properties.VibFreqOffset)./grad_phi_cos;
        end
    end
    
    
    toc
    
    % vshears=mean(vshearsin_pre,3);
    % vshearc=mean(vshearcos_pre,3);

    vshears=vshearsin_pre;
    vshearc=vshearcos_pre;
    %toc
    
    figure;
    
    vshears=medfilt2(vshears,[9 3]);
    vshearsin_im = vshears;
    
    imagesc( 10^3*Properties.Width_S,10^3*Properties.Depth_S,vshears);
    
    h = colorbar;
    ylabel(h, 'SWS m/s','FontSize',14);
    xlabel('Width [mm]','fontsize',14);ylabel('Depth [mm]','fontsize',14)
    title(['SWS AM-FM Vib Freq = ' num2str(Properties.VibFreq) ' Hz sin'],'fontsize',14)
    colormap turbo;
    set (gca,'clim',[2 6]);
   
    % 
    % figure;
    % 
    % vshearc=medfilt2(vshearc,[15 5]);
    % vshearcos_im = vshearc;
    % 
    % imagesc( 10^3*Properties.Width_S,10^3*Properties.Depth_S,vshearc);
    % 
    % h = colorbar;
    % ylabel(h, 'SWS m/s','FontSize',14);
    % xlabel('Width [mm]','fontsize',14);ylabel('Depth [mm]','fontsize',14)
    % title(['SWS AM-FM Vib Freq = ' num2str(Properties.VibFreq) ' Hz cos'],'fontsize',14)
    % colormap turbo;
    % set (gca,'clim',[2 6]);
    
    % subplot(212);imagesc(vshearsin);

%%
end
%toc

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