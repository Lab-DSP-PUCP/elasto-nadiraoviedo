
function [vshears,vshearc]= AMFM_demod (sono_filt_mov,Properties,num_frame)
%   Function that return a SWS frame stimator
%   Author: Nadira Oviedo
%   
%   Inputs:
%       sono_filt_mov       Filtered sonoelasticity video
%       Properties          Sonoelasticity video properties
%       num_frame           Frame to evaluate, if 0, it calculates the mean
%       of all frame

%   Outputs:
%        vshearsin           vshearcos

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
    

    if num_frame == 0
        vshearsin_pre = zeros(size(sono_filt_mov,1),size(sono_filt_mov,2),size(sono_filt_mov,3));
        vshearcos_pre = zeros(size(sono_filt_mov,1),size(sono_filt_mov,2),size(sono_filt_mov,3));
        for framen = 4:size(sono_filt_mov,3)
    
            %% Seleccion de senial 1D para analizar
            for slice = 1:size(sono_filt_mov,1)
                vect = sono_filt_mov(slice,:,framen);
                vect_A = hilbert(vect); %Crear senial analitica 
                x_esp = (0:size(sono_filt_mov,2)-1)*Properties.dx;

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
                

                vshearsin_pre(slice,:,framen)=(pi*Properties.VibFreq)./(grad_phi_sin*10^3)*1.2;
                vshearcos_pre(slice,:,framen)=(pi*Properties.VibFreq)./(grad_phi_cos*10^3)*1.2;

                % vshearsin_pre(slice,:,framen)=(pi*Properties.VibFreqOffset)./grad_phi_sin;
                % vshearcos_pre(slice,:,framen)=(pi*Properties.VibFreqOffset)./grad_phi_cos;
            end
        end
    
        vshears=mean(vshearsin_pre,3);
        vshearc=mean(vshearcos_pre,3);


    elseif num_frame <= size(sono_filt_mov,3) && num_frame>=5
        vshearsin_pre = zeros(size(sono_filt_mov,1),size(sono_filt_mov,2));
        vshearcos_pre = zeros(size(sono_filt_mov,1),size(sono_filt_mov,2));
        for framen = num_frame
            %% Seleccion de senial 1D para analizar
            for slice = 1:size(sono_filt_mov,1)
                vect = sono_filt_mov(slice,:,framen);
                vect_A = hilbert(vect); %Crear senial analitica
                
                x_esp = (0:size(sono_filt_mov,2)-1)*Properties.dx;
                
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

                vshearsin_pre(slice,:)=(pi*Properties.VibFreq)./(grad_phi_sin*10^3)*1.2;
                vshearcos_pre(slice,:)=(pi*Properties.VibFreq)./(grad_phi_cos*10^3)*1.2;

                % vshearsin_pre(slice,:)=(pi*Properties.VibFreqOffset)./grad_phi_sin;
                % vshearcos_pre(slice,:)=(pi*Properties.VibFreqOffset)./grad_phi_cos;
            end
        end
        vshears=vshearsin_pre;
        vshearc=vshearcos_pre;

    else
        disp('seleccione un frame dentro del rango')
    end
end