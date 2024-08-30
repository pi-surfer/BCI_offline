%% PROJECT 3 SCRIPT 2a

%  CALIBRATION PHASE (OFFLINE runs)

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

%profile on          % performance analysis of the code


%% EEG data loading, only OFFLINE runs

% NOTE: the EEG data of each subject are allocated in a folder named "a.._micontinuous".
% The 8 folders "a.._micontinuous" are allocated in another folder named "micontinuous",
% located in turn inside the current folder.

% The code searches which files are allocated inside the data directory and
% creates a variabile with the filenames, then used in the function "conc_sload":

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
       
    [si,hi,~,~] = conc_sload(eeg_sbj_i,1);
        
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
    n_off = hi.OfflineSessions;
    
    % Struct fields creation:    
    idx = 0;
    PSDi_fields = cell(1,n_files);
    H_psdi_fields = cell(1,n_files);
    
    for file = 1:n_off
        
        idx = idx+1;
        PSDi_field = ['PSD_on_', num2str(file)];
        PSDi_fields(idx) = {PSDi_field};

        H_psdi_field = ['hpsd_on_', num2str(file)];
        H_psdi_fields(idx) = {H_psdi_field};
    end
    
    % Spatial filtering and PSD computation for each j-th run:    
    for j = 1:hi.OfflineSessions
        
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


%% Feature extraction and selection

% Trial extraction and creation of the array
% [windows x frequencies x channels x trials].

% Decimal codes:
Fix_cross = 786;
both_hands = 773;
both_feet = 771;
Rest = 783;
Continuous_feedback = 781;
Target_hit = 897;
Target_miss = 898;

% Creation of figures for the results:
fig1 = figure('Name',['Classifier with 3 features'],...
              'NumberTitle','off','Position',[120 60 1200 700]);
    
fig2 = figure('Name',['Scatter plots subjects 4 and 7'],...
              'NumberTitle','off','Position',[300 200 800 600]);
          
fig4 = figure('Name',['Classifier with 5 features'],...
              'NumberTitle','off','Position',[120 40 1200 700]);

feature_definition
load("feature_selected.mat", 'F_5')

for i = 1:n_sbjs
   

    hi_psd = H_psd.(fields_H_psd{i});
    PSDi = PSD.(fields_PSD{i});
    PSDi_fields = fieldnames(PSDi);
    hi_psd_fields = fieldnames(hi_psd);
    
    %% Intra-subject concatenation of the PSD of the different runs
    
    EVENT_TYP_conc = [];
    EVENT_POS_conc = [];
    EVENT_DUR_conc = [];
    POS_offset = 0;
    PSDi_conc = [];

    for j = 1:hi_psd.OfflineSessions
        
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
    
    clear EVENT_TYP_conc EVENT_POS_conc EVENT_DUR_conc POS_offset 
    clear PSDi_j hi_psd_j j
    
    
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
    
    
    %% Feature selection using Fisher Score for visualization
    
    % Corresponding label extraction:    
    PSDi_fields = fieldnames(PSDi);
    hi_psd_fields = fieldnames(hi_psd);
    Fi_3D = zeros(channels,freq,hi_psd.OfflineSessions);
    
    
    for j = 1:hi_psd.OfflineSessions   % for each calibration run 
        
        % Corresponding label extraction:        
        PSDi_j = PSDi.(PSDi_fields{j});
        hi_psd_j = hi_psd.(hi_psd_fields{j});
        
        PSDj = PSDi_j;
        hj = hi_psd_j;

        % Extracting only the windows from the Cue to the Continuous
        % feedback periods for each class:
        
        start_pos = hj.EVENT.POS(hj.EVENT.TYP == both_hands | hj.EVENT.TYP == both_feet); 
        cf_pos = hj.EVENT.POS(hj.EVENT.TYP == Continuous_feedback);
        cf_dur = hj.EVENT.DUR(hj.EVENT.TYP == Continuous_feedback); 

        end_pos = cf_pos + cf_dur -1;
        duration = end_pos - start_pos +1;

        length_min = min(duration);
        ntrials = size(duration,1);
        flag = zeros(size(PSDj,1),1);
        classes = zeros(size(PSDj,1),1);

        
        for t = 1:ntrials
            flag(start_pos(t):end_pos(t),1) = 1;
            
            % Value assignment 1 to the flag vector to determine the windows to
            % be extracted for each trial and assignment .
            
            if hj.EVENT.TYP(hj.EVENT.POS == start_pos(t)) == both_hands
                classes(start_pos(t):end_pos(t),1) = both_hands;
                
            else
                classes(start_pos(t):end_pos(t),1) = both_feet;
            end
        end

        P_feet = PSDj(classes == both_feet,:,:);
        P_hands = PSDj(classes == both_hands,:,:);
        
        % Reshaping data in the format [windows x features]:       
        num_features = 23*16;
        num_windows_feet = size(P_feet,1);
        num_windows_hands = size(P_hands,1);
        class_feet = nan(num_windows_feet,num_features);
        class_hands = nan(num_windows_hands,num_features);

        val=1;
        for c = 1:channels
            for f = 1:freq
                class_feet(:,val)=P_feet(:,f,c);
                class_hands(:,val)=P_hands(:,f,c);
                val=val+1;
            end
        end
        
        % Computing Fisher score for each feature:        
        mean_feet = mean(class_feet);
        mean_hands = mean(class_hands);

        var_feet = var(class_feet);
        var_hands = var(class_hands);
        FS = nan(num_features,1);
        
        for f = 1:num_features
            FS(f) = (abs(mean_feet(f)-mean_hands(f)))/(sqrt(var_feet(f)+var_hands(f)));
        end

        % Converting back features index to channel/frequency:
        Fj = nan(channels,freq);
        index_start = 1:23:num_features;
        index_end = 23:23:num_features;
        for c = 1:channels
            Fj(c,:) = FS(index_start(c):index_end(c));
        end
        
        Fi_3D(:,:,j) = Fj;
    end

    n_runs = size(Fi_3D,3);
    
    % Visualizing the features maps for each calibration run:
    f1 = figure('Name',['Fisher Score - Subject ',num2str(i)],...
                'NumberTitle','off','Position',[120 280 1200 500]);
    sgtitle(['Fisher Score - Subject ',num2str(i)],'FontSize',16,'FontWeight','bold')
            
    xlabels = '';
    for l = 4:2:48      
        xlabels = append(xlabels,[num2str(l),' ']);
    end 
    xlabels(end) = [];
    xlabels = split(xlabels);
    
    for f = 1:n_runs
        subplot(1,n_runs+1,f)
        imagesc(Fi_3D(:,:,f))
        
        title(['Calibration run ',num2str(f)],'FontSize',14)
        ylabel('Channel','FontSize',12)
        xlabel('Frequency [Hz]','FontSize',12)
        xticks([1:1:23])
        xticklabels(xlabels)
        xtickangle(90)
        
    end
    
    % Visualizing the features selected:
    subplot(1,n_runs+1,n_runs+1)
    sel_feature = F_5.(fields_S{i});
    imagesc(sel_feature)
    
    title('Selected features')
    ylabel('Channel','FontSize',12)
    xlabel('Frequency [Hz]','FontSize',12)
    xticks([1:1:23])
    xticklabels(xlabels)
    xtickangle(90)


    %% Classifier
    % Vectors of chosen features for each subject
    ch = [9 9 9;...         %subject 1
        13 13 15;...
        11 7 9;...
        11 7 7;...
        7 12 2;...
        8 8 7;...
        2 4 14;...
        12 16 8 ];          %subject 8
    fr = [10 9 20;...
        5 6 5;...
        6 6 12;...
        5 4 5;...
        6 6 6;...
        5 6 6;...
        12 12 6;...
        5 5 5];             %subject 8
    
    % We initialize two vectors that will keep our samples and our labels
    holder = [];
    labels = Cki;
    
    % We take the informations about the PSD
    holder_big = Continuos_feedback_PSDi;
    holder_size = size(holder_big,1);
    di_max = size(holder_big,4);
    
    % We repeat the labels for each window in the PSD
    labels = repelem(labels,holder_size);

    % We unfold the matrix trial-by-trial
    for fi = 1:di_max
        sqzd = squeeze(holder_big(:,:,:,fi));
        % and we take a sample for each element in that trial
        
        for di=1:holder_size
            holder(end+1,:) = [sqzd(di,fr(i,1),ch(i,1)), sqzd(di,fr(i,2),ch(i,2)), sqzd(di,fr(i,3),ch(i,3))];
        end
    end


    % We provide a scatterplot for some specific subjects
    if i == 4 || i == 7
        figure(fig2)
            if i == 4
                subplot(211)
            else
                subplot(212)
            end
        idx_hand = holder(labels==773,:);
        idx_feet = find(labels == 771);
        scatter3(idx_hand(:,1), idx_hand(:,2),...
                 idx_hand(:,3) ,'.', 'MarkerFaceColor','#0072BD')
        hold on
        scatter3(holder(idx_feet,1), holder(idx_feet,2),...
            holder(idx_feet,3), '.', 'MarkerFaceColor','#D95319')
        title(['Scatter plot for subject ', num2str(i)],'FontSize',16,'FontWeight','bold')
    end

    Classifier = fitcdiscr(holder, labels,'DiscrimType','quadratic');

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
    sgtitle('Results with 3 features','FontSize',16,'FontWeight','bold')
    X = categorical({'Overall','Hands','Feet'});
    X = reordercats(X,{'Overall','Hands','Feet'});
    Y = [test_acc test_acc_hands test_acc_feet];
    bar(X,Y, 0.5,'grouped','g')
    title(['Single sample accuracy on subject ' num2str(i)])
    ylabel('Accuracy [%]')
    ylim([0,1])
   
    subplot(4,4,i+8)
    cm = confusionchart(Gk, labels);
    title(['Confusion matrix for subject ' num2str(i)])

    
    %% Classification using 5 features
    % The choice of the number of features is made manually because doing
    % it via forward selection is computationally infasible
    
       ch = [9 9 9 9 1;...         %subject 1
        13 13 15 13 15;...
        11 7 9 9 3;...
        11 7 7 11 2;...
        7 12 2 7 3;...
        8 8 7 7 11;...
        2 4 14 13 11;...
        12 16 16 8 10];         %subject 8
    
    fr = [10 9 20 8 9;...
        5 6 5 4 6;...
        6 6 12 11 6;...
        5 4 5 4 5;...
        6 6 6 5 6;...
        5 6 6 5 5;...
        12 12 6 6 5;...
        5 5 6 5 5];         %subject 8

    % we take again all the data and extract the 5 features
    holder = [];
    holder_big = Continuos_feedback_PSDi;
    holder_size = size(holder_big,1);
    di_max = size(holder_big,4);

     for fi=1:di_max
        sqzd = squeeze(holder_big(:,:,:,fi));
        % and we take a sample for each element in that trial
        for di=1:holder_size
            holder(end+1,:) = [sqzd(di,fr(i,1),ch(i,1)), sqzd(di,fr(i,2),ch(i,2)),...
                sqzd(di,fr(i,3),ch(i,3)), sqzd(di,fr(i,4),ch(i,4)), sqzd(di,fr(i,5),ch(i,5))];
        end
     end

    % we build new classifiers and plot the results in order to compare 
    % them with the ones obtained before
    Classifier = fitcdiscr(holder, labels,'DiscrimType','quadratic');
    [Gk, ~] = predict(Classifier, holder);

    logi_val = Gk == labels;
    hands_indices = find(labels == both_hands); 
    feet_indices = find(labels == both_feet);
    
    test_acc = sum(logi_val)/length(logi_val);
    test_acc_hands = sum(logi_val(hands_indices))/length(hands_indices);
    test_acc_feet = sum(logi_val(feet_indices))/length(feet_indices);


    %plot the results
    figure(fig4); 
    subplot(4,4,i)
    sgtitle('Results with 5 features','FontSize',16,'FontWeight','bold')
    X = categorical({'Overall','Hands','Feet'});
    X = reordercats(X,{'Overall','Hands','Feet'});
    Y = [test_acc test_acc_hands test_acc_feet];
    bar(X,Y, 0.5,'grouped','g')
    title(['Single sample accuracy on subject ' num2str(i)])
    ylabel('Accuracy [%]')
    ylim([0,1])
   
    subplot(4,4,i+8)
    cm = confusionchart(Gk, labels);
    title(['Confusion matrix for subject ' num2str(i)])

    % We save the classifiers in order to use them for the online sessions
    if not(isfolder("Classifiers"))
    mkdir("Classifiers")   
    end
    
    name = ['Classifiers\Classifier_', num2str(i), '.mat'];
    save(name,'Classifier');
    
end


%% Code performance

%profile off
%profile viewer

