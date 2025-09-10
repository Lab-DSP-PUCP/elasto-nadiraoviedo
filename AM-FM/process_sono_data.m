 function [sono_filt_mov,sono_filt,filt_mov]=process_sono_data(sono,Properties,move,rd)

% Median filter
sono_med = zeros(size(sono));
for t = 1:size(sono,3)
    sono_med(:,:,t) = medfilt2(sono(:,:,t), [9,3], 'symmetric');%[21,7]
end

% Depth normalization
[sono_filt,C,S] = normalize(sono_med,2);

% Normalized sono
sono_norm = normalize(sono,"center",C,"scale",S);

% 2D directional filter
[sono_filt_mov,filt_mov]=moving_filterv2(sono_filt,Properties.FrameRate, ...
    Properties.pitch,Properties.VibFreq,Properties.VibFreqOffset,rd,move);

end