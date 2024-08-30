function [PSD,f_sel,h_new] = lap_psd(s,h,lap,f_min,f_max)

% i,j per debug

%{
    Function to load and concatenate EEGdata.gdf

    Inputs:    - s:      matrix [samples x channels] (see also: 'sload')
               - h:      struct containing the main information related to
                         s (see also: 'sload')
               - lap:    matrix containing the laplacian kernel
               - f_min:  minimum frequecy
               - f_max:  maximum frequency

    Outputs:   - PSD:    matrix [windows x frequencies x channels]
               - f_sel:  vector contaning the selected frequencies
               - h_new:  h matrix computed on PSD

    SEE ALSO: 'proc_spectrogram', 'proc_pos2win'
%}

    s_lap = s*lap;
    
    wlength = 0.5;              % [s]
    pshift = 0.25;              % [s]
    wshift = 0.0625;            % [s]
    samplerate = h.SampleRate;
    mlength = 1;                % [s]

    [PSD, f] = proc_spectrogram(s_lap, wlength, wshift, pshift, samplerate, mlength);

    idx_min = find(f == f_min);
    idx_max = find(f == f_max);
    
    f_sel = f(idx_min:idx_max);     % [Hz]

    PSD = PSD(:,idx_min:idx_max,:);

    winconv = 'backward';
    h_new = h;
    h_new.EVENT.POS = proc_pos2win(h.EVENT.POS, wshift*h.SampleRate, winconv, mlength*h.SampleRate);
    
    end_POS_samp = h.EVENT.POS + h.EVENT.DUR -1;    
    end_POS_wind = proc_pos2win(end_POS_samp, wshift*h.SampleRate, winconv, mlength*h.SampleRate);    
    
    h_new.EVENT.DUR = end_POS_wind - h_new.EVENT.POS +1;
    
end