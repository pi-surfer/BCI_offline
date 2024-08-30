function Computing_ERD(s,h,subject)
%{
    Inputs:    - s:         matrix [samples x channels] (see also: 'sload')
               - h:         struct containing the main information related to
                            s (see also: 'sload')
               - subject    label of the subject analysed

    Outputs:   - figure:    ERD/ERS over Ce,Cz,C4 in mu e beta band
               - figure:    topoplot concerning the subject analysed of
               ERD/ERS during fixation and during activty, in mu and beta band, for
               both hands and feet
%}

    subject=char(subject);
   
    load('laplacian16.mat');   %loading the Laplacian mask provided
  
   
    %Application of the Laplacian filter
    s_lap=s*lap;

    %Filtering the signal in mu and beta bands
    mu_low=8;
    mu_high=12;
    beta_low=18;
    beta_high=22;
    sample_rate=512;
    N=5;
    [b a]=butter(N,2*[mu_low mu_high]/sample_rate);
    filter_mu=filtfilt(b,a,s_lap);
    [b a]=butter(N,2*[beta_low beta_high]/sample_rate);
    filter_beta=filtfilt(b,a,s_lap);
    
    %signal rectification
    srect_mu=filter_mu.^2;
    srect_beta=filter_beta.^2;
   
    %Application of the moving average (1-second window)
    sample_rate=512;
    sec=1;     
    window=sec*sample_rate; 
    kernel=ones(window,1)/window;
    mavg_mu=filter(kernel,1,srect_mu);
    mavg_beta=filter(kernel,1,srect_beta);
    
    %Logarithm transform
    log_mu=log(mavg_mu);
    log_beta=log(mavg_beta);
 
    %%  ERD/ERS
    
    %%Creation of the matrix trialdata_mu and trialdata_beta [samples x channels x trials]
    %Each trial starts from the event related to the fixation cross to the end of the event related to the continuous feedback
   
    start_pos=h.EVENT.POS(h.EVENT.TYP == 786); 
    ck_pos=h.EVENT.POS(h.EVENT.TYP == 781);
    ck_dur=h.EVENT.DUR(h.EVENT.TYP == 781);

    end_pos=ck_pos+ck_dur -1;
    durata=end_pos-start_pos+1;

    length_min=min(durata);
    channels=size(s,2);
    ntrials=length(start_pos);

    trialdata_mu=zeros(length_min,channels,ntrials);
    trialdata_beta=zeros(length_min,channels,ntrials);

    for t=1:ntrials
        trialdata_mu(:,:,t) = log_mu(start_pos(t):(start_pos(t)+length_min-1),:);
        trialdata_beta(:,:,t) = log_beta(start_pos(t):(start_pos(t)+length_min-1),:);
     
    end

    %Creation of the matrix fixdata_mu and fixdata_beta [samples x channels x trials] related only to the fixation period
    event=786;
    startPOS=h.EVENT.POS(h.EVENT.TYP==event); 
    finishDUR=h.EVENT.DUR(h.EVENT.TYP==event);
    L=min(finishDUR);
    ntrials=length(startPOS);

    for i=1:ntrials
        sample_start=startPOS(i);
        sample_stop=sample_start+L-1;
        fixdata_mu(:,:,i)=log_mu(sample_start : sample_stop,:);
        fixdata_beta(:,:,i)=log_beta(sample_start : sample_stop,:);
    end

    %Computing the ERD/ERS for mu and beta band
    Baseline_mu=repmat(mean(fixdata_mu), [size(trialdata_mu,1) 1 1]);
    ERD_mu=100*(trialdata_mu-Baseline_mu)./Baseline_mu;
    Baseline_beta=repmat(mean(fixdata_beta), [size(trialdata_beta,1) 1 1]);
    ERD_beta=100*(trialdata_beta-Baseline_beta)./Baseline_beta;
%%  Temporal visualization
    c3=7; cz=9; c4=11; %meaningful channels
    
    Cue_data(:)=h.EVENT.TYP(find(h.EVENT.TYP==773 | h.EVENT.TYP==771));
    Cue_data=Cue_data';

    %averaging across the two MI classes in mu band
    ERD_mu_771=ERD_mu(:,:,Cue_data==771);
    ERD_mu_773=ERD_mu(:,:,Cue_data==773);
    mean_ERD_mu_771=mean(ERD_mu_771,3);
    mean_ERD_mu_773=mean(ERD_mu_773,3);

    %averaging across the two MI classes in beta band
    ERD_beta_771=ERD_beta(:,:,Cue_data==771);
    ERD_beta_773=ERD_beta(:,:,Cue_data==773);
    mean_ERD_beta_771=mean(ERD_beta_771,3);
    mean_ERD_beta_773=mean(ERD_beta_773,3);
    
    %definition of the time vector
    sample_rate=h.SampleRate;
    n_sample=size(mean_ERD_beta_773,1);
    T=n_sample/sample_rate;  
    t=0:1/sample_rate:T-1/sample_rate;
    
%% plotting the first figure

    fig1 = figure();
    fig1.Name=sprintf('ERD/ERS over Ce,Cz,C4, subject %s',subject);
    sgtitle(['ERD/ERS over Ce,Cz,C4, subject: ',subject])
    
    subplot(321)
    plot(t,mean_ERD_mu_771(:,7))
    hold on
    plot(t,mean_ERD_mu_773(:,7))
    legend('both feet','both hands','Location','northwest')
    xlabel('Time [s]'); ylabel('[ERD/ERS]');
    title('channel C3 : \mu')

    subplot(322)
    plot(t,mean_ERD_beta_771(:,7))
    hold on
    plot(t,mean_ERD_beta_773(:,7))
    legend('both feet','both hands','Location','northwest')
    xlabel('Time [s]'); ylabel('[ERD/ERS]');
    title('channel C3 : \beta')

    subplot(323)
    plot(t,mean_ERD_mu_771(:,9))
    hold on
    plot(t,mean_ERD_mu_773(:,9))
    legend('both feet','both hands','Location','northwest')
    xlabel('Time [s]'); ylabel('[ERD/ERS]');
    title('channel Cz : \mu')

    subplot(324)
    plot(t,mean_ERD_beta_771(:,9))
    hold on
    plot(t,mean_ERD_beta_773(:,9))
    legend('both feet','both hands','Location','northwest')
    xlabel('Time [s]'); ylabel('[ERD/ERS]');
    title('channel Cz : \beta')

    subplot(325)
    plot(t,mean_ERD_mu_771(:,11))
    hold on
    plot(t,mean_ERD_mu_773(:,11))
    legend('both feet','both hands','Location','northwest')
    xlabel('Time [s]'); ylabel('[ERD/ERS]');
    title('channel C4 : \mu')

    subplot(326)
    plot(t,mean_ERD_beta_771(:,11))
    hold on
    plot(t,mean_ERD_beta_773(:,11))
    legend('both feet','both hands','Location','northwest')
    xlabel('Time [s]'); ylabel('[ERD/ERS]');
    title('channel C4 : \beta')

%% 
    %For each band and the two MI class, averaging the signal in the fixation period and in the continuous feedback period
    %fixation period
    event=786;
    fixcross_dur=h.EVENT.DUR(h.EVENT.TYP == event);
    L=min(fixcross_dur);
    startpos=1;
    endpos=startpos+L-1;
    fixcross_mu_771=mean_ERD_mu_771(startpos:endpos,:);
    fixcross_mu_773=mean_ERD_mu_773(startpos:endpos,:);
    fixcross_beta_771=mean_ERD_beta_771(startpos:endpos,:);
    fixcross_beta_773=mean_ERD_beta_773(startpos:endpos,:);
    
    %continuous feedback period
    event=781;
    contfeed_dur=h.EVENT.DUR(h.EVENT.TYP == event);
    L=min(contfeed_dur);
    start_pos=size(mean_ERD_mu_771,1)-L+1;

    confeed_mu_771=mean_ERD_mu_771(start_pos:end,:);
    confeed_mu_773=mean_ERD_mu_773(start_pos:end,:);
    confeed_beta_771=mean_ERD_beta_771(start_pos:end,:);
    confeed_beta_773=mean_ERD_beta_773(start_pos:end,:);

    load('chanlocs16.mat') %Loading the chanlocs provided for this channel layout

    p1=mean(fixcross_mu_771,1);
    p2=mean(fixcross_mu_773,1);
    p3=mean(confeed_mu_771,1);
    p4=mean(confeed_mu_773,1);

    p5=mean(fixcross_beta_771,1);
    p6=mean(fixcross_beta_773,1);
    p7=mean(confeed_beta_771,1);
    p8=mean(confeed_beta_773,1);

    val_min=-200;
    val_max=200;

%% Topoplot visualization

    fig2 = figure();
    fig2.Name=sprintf('Topoplot subject %s',subject);
    sgtitle(['ERD/ERS subject: ',subject])
    
    subplot(241)
    topoplot(p2,chanlocs16);
    caxis([val_min val_max])
    colorbar
    title('ERD/ERS during fixation (\mu band) - both hands')

    subplot(242)
    topoplot(p4,chanlocs16);
    caxis([val_min val_max])
    colorbar
    title('ERD/ERS during activity (\mu band) - both hands')

    subplot(245)
    topoplot(p1,chanlocs16);
    caxis([val_min val_max])
    colorbar
    title('ERD/ERS during fixation (\mu band) - both feet')

    subplot(246)
    topoplot(p3,chanlocs16);
    caxis([val_min val_max])
    colorbar
    title('ERD/ERS during activity (\mu band) - both feet')

    subplot(243)
    topoplot(p6,chanlocs16);
    caxis([val_min val_max])
    colorbar
    title('ERD/ERS during fixation (\beta band) - both hands')

    subplot(244)
    topoplot(p8,chanlocs16);
    caxis([val_min val_max])
    colorbar
    title('ERD/ERS during activity (\beta band) - both hands')

    subplot(247)
    topoplot(p5,chanlocs16);
    caxis([val_min val_max])
    colorbar
    title('ERD/ERS during fixation (\beta band) - both feet')

    subplot(248)
    topoplot(p7,chanlocs16);
    caxis([val_min val_max])
    colorbar
    title('ERD/ERS during activity (\beta band) - both feet')
end