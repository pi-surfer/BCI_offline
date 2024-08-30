function [S,H,s_conc,h_conc,fns] = conc_sload(fns,w)

%{
    Function to load and concatenate EEGdata.gdf

    Inputs:    - fns:   'char' type variable containing the file.gdf names
                         separated one from the other by ' ' or \n
                         E.g.: 'filename1.gdf filename2.gdf filename3.gdf'
               - w:      parameter which allows to choose which EEG data 
                         must be loaded and concatenated:

                                 w = 0       ALL EEG data provided
                                 w = 1       OFFLINE EEG data
                                 w = 2       ONLINE EEG data

    Outputs:   - S:      structcontaining the s-matrices from 'sload'
                         of the selected files (according to w)
               - H:      struct containing the h-structures from 'sload'
                         of the selected files (according to w)
               - s_conc: matrix resulting from the concatenation of the
                         s-matrices of the selected files (according to w)
               - h_conc: struct containing the main informations
                         associated to s_conc
%}

fns = sort(split(fns));

% fns is ordered so that offline files occupy the first n_off positions 
% and online files occupy the others:

c = strfind(fns,'online');

for i = 1:length(c)
    c{i} = num2str(c{i});
end

n_off = length(find(strcmp(c,'') == 1));
n_files = length(fns);

j = 0;
k = 0;
fns_ordered = cell(n_files,1);

for i = 1:n_files
    if strcmp(c{i},'') == 1
        j = j+1;
        fns_ordered{j} = fns{i};
    else
        k = k+1;
        fns_ordered{n_off+k} = fns{i};
    end
end

fns = fns_ordered;

% Selection of which files must be concatenated:

if w == 0
    n_start = 1;
    n_stop = length(fns);

elseif w == 1
    n_start = 1;
    n_stop = n_off;

elseif w == 2
    n_start = n_off+1;
    n_stop = length(fns);

end    

% Loading and concatenation of the selected files one by one:

EVENT_TYP_conc = [];
EVENT_POS_conc = [];
EVENT_DUR_conc = [];
POS_offset = 0;
s = [];

idx = 0;
S_fields = cell(1,n_files);
H_fields = cell(1,n_files);

for f = 1:n_off
    
    idx = idx+1;
    Si_field = ['s_off_', num2str(f)];
    S_fields(idx) = {Si_field};
    
    Hi_field = ['h_off_', num2str(f)];
    H_fields(idx) = {Hi_field};
end
    
for f = 1:(n_files - n_off)
    
    idx = idx+1;
    Si_field = ['s_on_', num2str(f)];
    S_fields(idx) = {Si_field};
    
    Hi_field = ['h_on_', num2str(f)];
    H_fields(idx) = {Hi_field};
end

idx = n_start-1;
s_conc = [];

for n = n_start:n_stop
    
    idx = idx+1;
    
    fn = char(fns(n));
    [sn,hn] = sload(fn);
    sn = sn(:,1:16);
    
    [S.(S_fields{idx})] = sn;
    [H.(H_fields{idx})] = hn;
    
    s_conc = cat(1,s_conc,sn);
    
    EVENT_TYP_conc = cat(1,EVENT_TYP_conc,hn.EVENT.TYP);
    EVENT_POS_conc = cat(1,EVENT_POS_conc,(hn.EVENT.POS) + POS_offset);
    EVENT_DUR_conc = cat(1,EVENT_DUR_conc,hn.EVENT.DUR);
    POS_offset = POS_offset + size(sn,1);
    
end

EVENT_conc = struct('TYP',EVENT_TYP_conc,'POS',EVENT_POS_conc,'DUR',EVENT_DUR_conc);

ch_labels = ["Fz" num2str(1);"FC3" num2str(2);"FC1" num2str(3);"FCz" num2str(4);"FC2" num2str(5);"FC4" num2str(6);...
            "C3" num2str(7);"C1" num2str(8);"Cz" num2str(9);"C2" num2str(10);"C4" num2str(11);"CP3" num2str(12);...
            "CP1" num2str(13);"CPz" num2str(14);"CP2" num2str(15);"CP4" num2str(16);"GND" num2str(17)];

h_conc = struct('EVENT',EVENT_conc,'SampleRate',hn.SampleRate,'ChannelLabels',ch_labels);

h_conc.OfflineSessions = n_off;
h_conc.OnlineSessions = length(fns) - n_off;

H.OfflineSessions = n_off;
H.OnlineSessions = length(fns) - n_off;
H.Sessions = length(fns);

end