function [sono_filt_mov,filt_mov]=moving_filterv2(sono_filt,FR,pitch,f,df,sd,mov)

% Desired frequency response
[f1,f2] = freqspace([500 500],'meshgrid');
ft = df/FR*2;
fs1 = f/sd(2)*pitch*2*2;
fs2 = f/sd(1)*pitch*2*2;
Hd = ( abs(f1)<ft + 0.05 )&( abs(f1)>ft - 0.05 )&( abs(f2) > fs1 )&(abs(f2) < fs2);

% Filter direction
L=bwlabel(Hd);
if strcmp(mov,'right')
    L(L==1)=0;
    L(L==4)=0;
elseif strcmp(mov,'left')
    L(L==2)=0;
    L(L==3)=0;
else
    disp('Incorrect direction');
end
% Hd = Hd & L;

% 2D Hamming window
im_size = [size(sono_filt,2) size(sono_filt,3)];
win_dim = [floor(im_size(1)*0.1)*2+1,floor(im_size(2)*0.25)*2+1];
win = hamming(win_dim(1)).*hamming(win_dim(2))';

% Windowing filter
hd = ifft2(ifftshift(Hd));
hd_cent = fftshift(hd);
N1 = floor(size(hd,1)/2) + 1;
N2 = floor(size(hd,2)/2) + 1;
Nwin1 = floor(size(win,1)/2);
Nwin2 = floor(size(win,2)/2);
hd_tr = hd_cent(N1-Nwin1:N1+Nwin1, N2-Nwin2:N2+Nwin2);
h_final = (hd_tr.*win);

% Normalization
H_final = fftshift(fft2(h_final,500,500));
norm_factor = abs(max(H_final(:)));
filt_mov = h_final/norm_factor;

% surf(f1(1,:)*FR/2,f2(:,1)/pitch/2,abs(fftshift(fft2(filt_mov,500,500))),'EdgeColor','none')
% xlabel("Temporal freq [Hz]")
% ylabel("Spatial freq [m^{-1}]")
% colorbar


sono_filt_mov = zeros(size(sono_filt));
for z_id = 1:size(sono_filt,1)
    sono_slice = squeeze(sono_filt(z_id,:,:));
    sono_filt_mov(z_id,:,:) = filter2(filt_mov,sono_slice,'same');
    % if z_id == 50
    %     subplot(131),imagesc(abs(fftshift(fft2(sono_slice,500,500)))),axis off;
    %     subplot(132),imagesc(abs(fftshift(fft2(filt_mov,500,500)))),axis off;
    %     subplot(133),imagesc(abs(fftshift(fft2(squeeze(sono_filt_mov(z_id,:,:)),500,500)))),axis off;
    % 
    % end
end

filt_mov = fft2(filt_mov);

