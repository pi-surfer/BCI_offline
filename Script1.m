%% PROJECT 3 SCRIPT 1 (offline runs)

%{
    Data recorded with 16-channel EEG amplifier (g.USBamp, g.Tec) at 512 Hz. 
    Electrodes layout: 10-20 international.

    Each subject: 2 or 3 "offline" runs (calibration, no real feedback)
                  6 "online" runs (with real feedback)
%}

close all;clear;clc  

 % NOTE: It's necessary to add the path according to the PC where eeglab is saved.
%addpath('C:\Users\andre\Desktop\Neuro\laboratorio neuro\eeglab_current\eeglab2022.1');
eeglab

% NOTE: It's necessary to change this path according to the PC and
% where the "conc_sload()" function is saved.
addpath 'C:\Users\franc\OneDrive\Documenti\Università\(1-1) Neurorobotics and Neurorihabilitation'

%% Importing and intra-subject concatention of EEG data, only OFFLINE runs

% NOTE: the EEG data of each subject are allocated in a folder named "a.._micontinuos".
% The 8 folders "a.._micontinuos" are allocated in another folder named "micontinus",
% located in turn inside the current folder.

% The code searches which files are allocated inside the data directory and
% creates a variabile with the filenames used in the function "conc_sload":

if ispc                     % True for the PC (Windows) version of MATLAB
    d = dir('micontinuous');   %dir lists the files in the folder micontinuous
    dir_content = {};
    %save within dir_content the names of the subfolders of micontinuous referring to each subject
    for i = 3:10
        subject = d(i).name;
        dir_content{i-2} = subject;
    end
    cd('micontinuous');
    n_sbjs = length(dir_content);
    
else                        % For the PC (macOS or Linux) version of MATLAB
    
    cd('micontinuous');

    dir_content = ls;
    dir_content(end) = [];
    dir_content = sort(split(dir_content));

    n_sbjs = length(dir_content);
end


% The outputs of this section will be 2 structures S and H which contain
% the s array [samples x channels] and the struct h of each run (only offline,
% only online or both according to the last input of "conc_sload") of each subject.

%defining labels
fields_S = {'s1','s2','s3','s4','s5','s6','s7','s8'};
fields_H = {'h1','h2','h3','h4','h5','h6','h7','h8'};

fields_S_conc = {'s1_conc','s2_conc','s3_conc','s4_conc','s5_conc','s6_conc',...
                 's7_conc','s8_conc'};
fields_H_conc = {'h1_conc','h2_conc','h3_conc','h4_conc','h5_conc','h6_conc',...
                'h7_conc','h8_conc'};
            

for d = 1:n_sbjs %For each subject, we analyze the contents of the folder provided
    
    if ispc
    
        c = dir(dir_content{d});
        fin = size(c,1);
        for a = 3:fin    %save within eeg_sbj_i the names of the files contained in the reference folder
             file_name = c(a).name;
             eeg_sbj_i{a-2} = file_name;
        end        
        cd(dir_content{d})
        
    else
        cd(dir_content{d});
        eeg_sbj_i = ls;
        eeg_sbj_i(end) = [];
    end
    
 
    
    [si,hi,si_conc,hi_conc] = conc_sload(eeg_sbj_i,1);
    
    %Saving within different structures of the data extracted by conc_sload
    [S.(fields_S{d})] = si;
    [H.(fields_H{d})] = hi;
    [S_conc.(fields_S_conc{d})] = si_conc;
    [H_conc.(fields_H_conc{d})] = hi_conc;
    
    %Computing the ERD/ERS for each subject
    Computing_ERD(si_conc,hi_conc,fields_S(d))
    
    cd ..      
end
    
cd ..