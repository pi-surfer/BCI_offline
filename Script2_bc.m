%% PROJECT 3 SCRIPT 2b-c

%  EVALUATION PHASE (ONLINE runs)

%{
    Data recorded with 16-channel EEG amplifier (g.USBamp, g.Tec) at 512 Hz. 
    Electrodes layout: 10-20 international.

    Each subject: 2 or 3 "offline" runs (calibration, no real feedback)
                  6 "online" runs (with real feedback)
%}

close all
clear
clc

addpath 'C:\Users\franc\OneDrive\Documenti\Università\(1-1) Neurorobotics and Neurorihabilitation'

% NOTE: it's necessary to add the path to the conc_sload function

%profile on          % measurement and analysis of the code's performance


%% Importing and intra-subject concatention of EEG data, only ONLINE runs

% NOTE: the EEG data of each subject are allocated in a folder named "a.._micontinuos".
% The 8 folders "a.._micontinuos" are allocated in another folder named "micontinus", 
% located in turn inside the current folder.

if ispc                     % True for the PC (Windows) version of MATLAB
    
    d = dir('micontinuous'); 
    dir_content = {};
    
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

% Defining labels:
fields_S = {'s1','s2','s3','s4','s5','s6','s7','s8'};
fields_H = {'h1','h2','h3','h4','h5','h6','h7','h8'};

for d = 1:n_sbjs 
    
    if ispc
    
        c = dir(dir_content{d});
        fin = size(c,1);
        for a = 3:fin    
             file_name = c(a).name;
             eeg_sbj_i{a-2} = file_name;
        end        
        cd(dir_content{d})
        
    else
        cd(dir_content{d});
        eeg_sbj_i = ls;
        eeg_sbj_i(end) = [];
    end
       
    [si,hi,~,~] = conc_sload(eeg_sbj_i,2);
        
    [S.(fields_S{d})] = si;
    [H.(fields_H{d})] = hi;
    
    cd ..      
    clear eeg_sbj_i
end
    
cd ..

clear d eeg_sbj_i i si hi si_conc hi_conc dir_content
%clear fields_H fields_S n_sbjs

% NOTE: The choice of using struct S and H was necessary in order to obtain 
% one code for all subjects that can automate the procedure while mantaing
% a good code usability, code evaluation and debug by sections.


%% Spatial filtering and PSD computation

% Spatial kernel loading:
load('laplacian16.mat');

% Frequency band selection:
f_min = 4;          % [Hz]
f_max = 48;         % [Hz]

fields_PSD = {'PSD1','PSD2','PSD3','PSD4','PSD5','PSD6','PSD7','PSD8'};
fields_H_psd = {'H_psd1','H_psd2','H_psd3','H_psd4','H_psd5','H_psd6',...
                'H_psd7','H_psd8'};

for i = 1:n_sbjs                % For each i-th subject
    
    fields_S = fieldnames(S);
    fields_H = fieldnames(H);
    
    si = S.(fields_S{i});
    hi = H.(fields_H{i});
    
    fields_sj = fieldnames(si);
    fields_hj = fieldnames(hi);
    
    n_files = hi.Sessions;
    n_on = hi.OnlineSessions;
    
    % Struct fields creation:   
    idx = 0;
    PSDi_fields = cell(1,n_files);
    H_psdi_fields = cell(1,n_files);
    
    for file = 1:n_on
        
        idx = idx+1;
        PSDi_field = ['PSD_on_', num2str(file)];
        PSDi_fields(idx) = {PSDi_field};

        H_psdi_field = ['hpsd_on_', num2str(file)];
        H_psdi_fields(idx) = {H_psdi_field};
    end

    % Spatial filtering and PSD computation for each j-th run:    
    for j = 1:hi.OnlineSessions
        
        sj = si.(fields_sj{j});
        sj = sj(:,1:16);            % Reference channel exclusion
        hj = hi.(fields_hj{j});
                
        [PSDj,f_sel,h_PSDj] = lap_psd(sj,hj,lap,f_min,f_max);
        
        [PSDi.(PSDi_fields{j})] = PSDj;
        [H_psdi.(H_psdi_fields{j})] = h_PSDj;
    end
    
    H_psdi.OfflineSessions = hi.OfflineSessions;
    H_psdi.OnlineSessions = hi.OnlineSessions;
    H_psdi.Sessions = hi.Sessions;
    
    [H_psd.(fields_H_psd{i})] = H_psdi;
    [PSD.(fields_PSD{i})] = PSDi;

end

H_psd.SelectedFrequencies = f_sel;

clear file f_min f_max i j lap si sj hi hj PSDj h_PSDj PSDi n_off n_files idx H_psdi
%clear fields_PSD fields_H_psd fields_S fields_H fields_sj fields_hj 
clear PSDi_field PSDi_fields H_psdi_field H_psdi_fields


%% Feature extraction

% Trial extraction and creation of the array
% [windows x frequencies x channels x trials]

% Decimal codes:
Fix_cross = 786;
both_hands = 773;
both_feet = 771;
Rest = 783;
Continuous_feedback = 781;
Target_hit = 897;
Target_miss = 898;

% Creation of figures for the results:
fig1 = figure('Name',['Test Accuracy'],...
              'NumberTitle','off','Position',[120 60 1200 700]);
    
fig2 = figure();
          
fig4 = figure('Name',['Scatter plots subjects 5 and 7'],...
              'NumberTitle','off','Position',[300 200 800 600]);
          
fields_H_psd_conc = {'H_psd1','H_psd2','H_psd3','H_psd4','H_psd5','H_psd6',...
                     'H_psd7','H_psd8'};
                 
times_m = [];
    
for i = 1:n_sbjs
   
    hi_psd = H_psd.(fields_H_psd{i});
    PSDi = PSD.(fields_PSD{i});
    PSDi_fields = fieldnames(PSDi);
    hi_psd_fields = fieldnames(hi_psd);
    
    %% PSD_conc computation
    
    EVENT_TYP_conc = [];
    EVENT_POS_conc = [];
    EVENT_DUR_conc = [];
    POS_offset = 0;
    PSDi_conc = [];

    for j = 1:hi_psd.OnlineSessions
        
        PSDi_j = PSDi.(PSDi_fields{j});
        hi_psd_j = hi_psd.(hi_psd_fields{j});       
        
        PSDi_conc = vertcat(PSDi_conc,PSDi_j);
        
        EVENT_TYP_conc = cat(1,EVENT_TYP_conc,hi_psd_j.EVENT.TYP);
        EVENT_POS_conc = cat(1,EVENT_POS_conc,(hi_psd_j.EVENT.POS) + POS_offset);
        EVENT_DUR_conc = cat(1,EVENT_DUR_conc,hi_psd_j.EVENT.DUR);
        
        POS_offset = POS_offset + size(PSDi_j,1);
    end

    EVENT_conc = struct('TYP',EVENT_TYP_conc,'POS',EVENT_POS_conc,'DUR',EVENT_DUR_conc);   
    hi_psd_conc = struct('EVENT',EVENT_conc);
    
    [H_psd_conc.(fields_H_psd_conc{i})] = hi_psd_conc;
    
    
    %% Trial extraction
    
    % Each trial starts with the fixation cross period and it ends with the
    % last window of the continuous feedback period.
    
    POS = hi_psd_conc.EVENT.POS;
    TYP = hi_psd_conc.EVENT.TYP;
    DUR = hi_psd_conc.EVENT.DUR;
    
    start_pos = POS(TYP == Fix_cross); 
    ck_pos = POS(TYP == Continuous_feedback);
    ck_dur = DUR(TYP == Continuous_feedback);

    end_pos = ck_pos + ck_dur -1;
    durata = end_pos - start_pos +1;

    length_min = min(durata);
    channels = size(PSDi_conc,3);
    ntrials = length(start_pos);
    freq = size(PSDi_conc,2);

    % Creation of the array [windows x frequencies x channels x trials]:  
    TrialPSDi = zeros(length_min,freq,channels,ntrials);

    for t = 1:ntrials
        TrialPSDi(:,:,:,t) = PSDi_conc(start_pos(t):(start_pos(t)+length_min-1),:,:);
    end
    
    
    %% Extraction of Ck := trials cue information vector (labels)

    val_cue1 = both_hands;
    val_cue2 = both_feet;

    idx_cue1 = find(TYP == val_cue1);
    idx_cue2 = find(TYP == val_cue2);

    a = 1;      b = 1;
    curr_idx1 = idx_cue1(a);
    curr_idx2 = idx_cue2(b);
    Cki = zeros(ntrials,1);

    for t = 1:ntrials

        if a+1 > length(idx_cue1)
            Cki(t:end) = val_cue2;

        elseif b+1 > length(idx_cue2)
            Cki(t:end) = val_cue1;

        elseif curr_idx1 < curr_idx2
            Cki(t) = val_cue1;
            a = a+1;
            curr_idx1 = idx_cue1(a);

        else
            Cki(t) = val_cue2;
            b = b+1;
            curr_idx2 = idx_cue2(b);
        end
    end

    clear a b curr_idx1 curr_idx2 t
    
    
    %% Continuos feedback extraction
    
    start_pos = POS(TYP == Continuous_feedback); 
    cf_dur = DUR(TYP == Continuous_feedback);
    
    typ_cf_idxs = find(TYP == Continuous_feedback);

    end_pos = start_pos + cf_dur -1;
    durata = end_pos - start_pos +1;

    length_min = min(durata);
    channels = size(PSDi_conc,3);
    ntrials = length(start_pos);
    freq = size(PSDi_conc,2);

    Continuos_feedback_PSDi = zeros(length_min,freq,channels,ntrials);

    for t = 1:ntrials
        Continuos_feedback_PSDi(:,:,:,t) = PSDi_conc(start_pos(t):(start_pos(t)+length_min-1),:,:);
    end

    
    %% Classifier
    
    ch = [9 9 9 9 1;...         %subject 1
        13 13 15 13 15;...
        11 7 9 9 3;...
        11 7 7 11 2;...
        7 12 2 7 7;...
        8 8 7 7 11;...
        2 4 14 13 16;...
        12 16 16 8 10];         %subject 8
    
    fr = [10 9 20 8 9;...
        5 6 5 4 6;...
        6 6 12 11 6;...
        5 4 5 4 5;...
        6 6 6 5 12;...
        5 6 6 5 5;...
        12 12 6 6 4;...
        5 5 6 5 5];         %subject 8

    % we take again all the data and extract the 5 features
    holder = [];
    labels = Cki;
    %We take the informations about the PSD
    holder_big = Continuos_feedback_PSDi;
    holder_size = size(holder_big,1);
    di_max = size(holder_big,4);
    %We repeat the labels for each window in the PSD
    labels = repelem(labels,holder_size);


     for fi=1:di_max
        sqzd = squeeze(holder_big(:,:,:,fi));
        % and we take a sample for each element in that trial
        for di=1:holder_size
            holder(end+1,:) = [sqzd(di,fr(i,1),ch(i,1)), sqzd(di,fr(i,2),ch(i,2)),...
                sqzd(di,fr(i,3),ch(i,3)), sqzd(di,fr(i,4),ch(i,4)), sqzd(di,fr(i,5),ch(i,5))];
        end
     end


    % We provide a scatterplot for some specific subjects
    if i == 5 || i == 7
        figure(fig4)
            if i == 5
                subplot(211)
            else
                subplot(212)
            end
        idx_hand = holder(find(labels==773),:);
        idx_feet = find(labels == 771);
        scatter3(idx_hand(:,1), idx_hand(:,2),...
                 idx_hand(:,3) ,'.', 'MarkerFaceColor','#0072BD')
        hold on
        scatter3(holder(idx_feet,1), holder(idx_feet,2),...
            holder(idx_feet,3), '.', 'MarkerFaceColor','#D95319')
        title(['Scatter for subject ', num2str(i)])
    end
    
    name = ['Classifiers\Classifier_', num2str(i), '.mat'];
    Classifier = load(name);
    Classifier = Classifier.Classifier;

    %obtain results
    [Gk, pp] = predict(Classifier,holder);

    logi_val = Gk == labels;
    hands_indices = find(labels == both_hands); 
    feet_indices = find(labels == both_feet);
    
    test_acc = sum(logi_val)/length(logi_val);
    test_acc_hands = sum(logi_val(hands_indices))/length(hands_indices);
    test_acc_feet = sum(logi_val(feet_indices))/length(feet_indices);


    %plot the results
    figure(fig1);
    subplot(4,4,i)
    sgtitle('Results on online data','FontSize',16,'FontWeight','bold')
    X = categorical({'Overall','Hands','Feet'});
    X = reordercats(X,{'Overall','Hands','Feet'});
    Y = [test_acc test_acc_hands test_acc_feet];
    bar(X,Y, 0.5,'grouped','g')
    ylim([0 1])
    title(['Single sample accuracy on subject ' num2str(i)])
    ylabel('Accuracy [%]')
   
    subplot(4,4,i+8)
    cm = confusionchart(Gk, labels);
    title(['Confusion matrix for subject ' num2str(i)])


    % We initialize the variables for the posterior probability
    %figure(fig2);
    alpha = 0.5;
    ex = 0;
    D = 0.5 * ones(size(labels,1), 2);

    % Compute the pp
    for wId = 2:size(labels,1)-1
        if ex == holder_size    % at the start of each trial
            ex = 0; 
        else 
            % We use the formula to compute it
            D(wId) = D(wId - 1) * alpha + pp(wId) * (1 - alpha);
            ex = ex + 1;
        end
    end

    
    % We plot a part of the Posterior Probability just for subject i as a
    % control
    
    if i == 5
        figure(fig2)
        fig2.Name= sprintf('Posterior probability for subject %d', i)
        plot(0:1:200-1,D(1:200,:))
        hold on
        plot(0:6:200-1,pp((1:6:200),1),'o')
        yline(0.2)
        yline(0.8)
        hold off
    
        yticks([0.2 0.5 0.8])
        yticklabels({'both hands = 0.2 ','0.5','both feet = 0.8'})
        title(['A posteriori probability for subject' num2str(i)])
    end
    
    
    % We now computer the average time to deliver a command
    nwindows = 0;
    time_w = [];
    flag = 0;
    
    for d = 1:size(D,1)
        % When the pp reaches the threshold we stop counting until we have
        % a new window
        if D(d,1) >= 0.8 || D(d,1) <= 0.2
            if flag == 0              
                nwindows = nwindows + 1;
                time_w(end+1) = nwindows;
                nwindows = 0;
                flag = 1;
            else
                continue
            end
            
        elseif D(d,1) == 0.5
            flag = 0;
            
        else
            nwindows = nwindows + 1;
        end
    end

    % These informations come from lap_psd
    wlength = 0.5;              % [s]
    wshift = 0.0625;            % [s]
    
    % We remove the time due to the overlap
    time = time_w*wlength - wshift*(time_w-1);
    time_m = mean(time);
    times_m(end+1) = time_m;
    
end

%% We plot the average time for each subject
fig3 = figure('NumberTitle','off','Position',[300 200 800 600]);
fig3.Name = sprintf('Avg time to deliver a command')
sgtitle('Average time to deliver a command','FontSize',16,'FontWeight','bold')
stem(times_m,'linewidth',2)
grid on
xlabel('Subject','FontSize',14)
ylabel('Time [s]','FontSize',14)


%% Code performance

%profile off
%profile viewer

