function SWS_FSST = process_FSST(sono_filt_mov,Properties,rd,window_length)
%   Function that return a SWS video using the short time Fourier Transform
%   Author: Joaquin Sanchez
%   
%   Inputs:
%       sono_filt_mov       Filtered sonoelasticity video
%       Properties          Sonoelasticity video properties
%       rd                  SWS range (recommendation: no higher than 8 m/s)
%
%   Outputs:
%       SWS_FSST            SWS video

% Window size estimation

max_lambda = rd(2)./ (2*Properties.VibFreq); % using cmax
win_len = floor(max_lambda/Properties.pitch *0.5)*2+window_length;
window = hann(win_len);
Fo = 1/Properties.pitch;
%extend_len = floor(win_len/2); % for zero padding
%overlap = length(window)-1;

%c = linspace(rd(1),rd(2),128); % range for c
%fk = flip(2*Properties.VibFreq ./ c); % Posible values for spatial frequency fk

[Nz,Nx,Nt] = size(sono_filt_mov);
%Nx_1=Nx+win_len;
SWS_FSST = zeros([Nz,Nx,Nt]);
vibf = Properties.VibFreq;
parfor t_id = 1:Nt
    %U = wextend('addcol','ppd',sono_filt_mov(:,:,t_id),extend_len); % zero padding
    U = sono_filt_mov(:,:,t_id);
    fk_mat = zeros([Nz,Nx]);%Nx
    for i=1:Nz
        [s,w] = fsst(U(i,:),Fo, window, 'yaxis');
        %threshold = 1; % Define un umbral de energía
        %s(abs(s) < threshold) = 0; % Zeros para valores muy pequeños de s
        [~,id_f] = max(s);
        fk_mat(i,:) = w(id_f)';
        for j = 2:Nx  % Comenzamos desde la segunda columna
            if fk_mat(i, j) == 0
                fk_mat(i, j) = fk_mat(i, j-1);  % Reemplazar 0 con el valor anterior
            end
        end
    end

    % Definir el índice de las columnas que se van a seleccionar
    %indices = round(linspace(1, Nx_1, Nx));

    % Crear la nueva matriz seleccionando las columnas correspondientes
    %fk_mat_1 = fk_mat(:, indices);

    % c = 2f/fk
    % SWS_frame = 2*vibf ./ fk_mat;
    % SWS_FSST(:,:,t_id) = SWS_frame;

    modified_fk_mat = change_zeros(fk_mat);
    %modified_fk_mat = new_matrix(fk_mat,rd(1),rd(2));
    SWS_frame = 2*vibf ./ modified_fk_mat;
    SWS_FSST(:,:,t_id) = SWS_frame;
end

end