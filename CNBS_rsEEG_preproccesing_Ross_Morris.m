% Berensson-Allen Center for Non-invasive Brain Stimulaiton 
% CNBS Resting-State Semi-Automated Pipeline (R-SAP)

% Authors: Jessica M. Ross, PhD, jross4@bidmc.harvard.edu, Timothy P. Morris, PhD, tpmorris@bidmc.harvard.edu

% This pipeline uses EEGlab plugins (including clean_raw
% (ASR)), MARA ICA component selection, and TESA component selection to
% preproccess resting state EEG data:

% The pipeline initially cleans channels and removes artefacts (using ASR),
% ICA is then performed and MARA is used with a custom threshold to automatically delete only
% components with a very high probability of being artefact (reducing the
% possibilty for false positives = deleting components that are infact
% neural). Components that are kept as "likely neural' are sorted by their explained
% variance and then visualized in TESA (Rogasch et, al 2017) to be manually
% confirmed as either neural or artefact (those missed by MARA). 

% Two .txt files are generated with output varibles from the cleaning of
% each dataset. If cleaning in batch, results are appended to these .txt
% files. 

% Sections of the script need to be changed based on indivudal aspects or requirements of
% each individual dataset(or data batch):

% Lines 41-42, 309-311, 396-398, 449-451: homefolder, datafolder and savefolder should be identified in each step, and in the MATLAB list of paths (ex: addpath(genpath(''))
% Line 69-70: Bandpass filters (if low pass cutoff is over 50Hz, a notch filter for line noise will also need to be applied)
% Line 73-75: Epoch length 
% Line 79: Re-reference
% Line 81 : Use PCA compression or do not
% Line 85: PCA dimensions (if using PCA)
% Line 65, 100: File type is .vhdr, but if you have a different raw data file format, these lines will need to change
% Line 113: EMG channel deletions 

% After step 1 check for % of data kept after ASR (see command window after
% running step 1)

%% 

clear all; close all; clc;
homefolder = 'D:\SAGESII_P034_GB\SAGESII_P034_GB_V2_1_14_2021\Analysis\Andrei_training\Resting';
datafolder = 'D:\SAGESII_P034_GB\SAGESII_P034_GB_V2_1_14_2021\Analysis\Andrei_training\Resting';
cd(datafolder);

% Set up output file
tic
% Make a file for output
mkdir('auto_output');
addpath(genpath('auto_output/'));

% Open summary file for output
summ_path = strcat('auto_output/','auto_pipeline_summary');
file_name = [summ_path sprintf('%d') '.txt'];
fopen(file_name,'a+');

% Make header
fid = fopen(file_name,'a+');
fprintf(fid,'\n \n %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s' ,'file name',',','high pass (Hz)',',','low pass (Hz)',',','epoch length',',','PCA (0=no 1=yes)',',','downsampled to (Hz)',',','number of channels removed',',','reference',',','number of ICs before MARA',',','number of MARA comps deleted',',','[MARA comps deleted]');
  
% ========================================================================
% Step 1: Load Data, Upload Channels, Delete non-EEG Channels, Save -
% BRAINVISION
% ========================================================================

% Files = dir('*.vhdr');
Files = dir('*S01.set')
gate_channel = [];
channelsremoved = [];

bandlow=1; % VAR - BANDPASS LOWER EDGE
bandhigh=49; % VAR - BANDPASS UPPER EDGE
filtord = 4; %VAR - FILTER ORDER

EpochLength = 3; %VAR - EPOCH LENGTH TO DIVIDE CONTINUOUS DATA INTO
epoch_begin = 0; %Duration of time before marker to begin epoch, in seconds
epoch_end = 3; %Duration of time after marker to end epoch, in seconds

eventmarker = 'epochmark';

Rereference = 1; % VAR - Set to 0 for no rereferencing, 1 for average reference

DoPCA = 1; %VAR: Whether to do dimensionality reduction w PCA prior to ICA: 0 for no, 1 for yes
CalculateDimensions = 0; %VAR: Whether to calculate number of dimensions to reduce data to prior to ICA
PercentVar = 99; %VAR: If calculating dimensions, percent of variance to explain in PCA reduction
MinComp = 20; %Minimum number of components to include if calculating variance
PCAdimensions = 30; %VAR: Number of dimensions to reduce to if doing PCA but not calculating dimensions
ICAmethod = 'fast ICA'; %VAR: ICA type
% ICAmethod = 'infomax';

cd(homefolder);

 %-------------------------------------------------------------------------
for i = 1 : length(Files)
    %Load the file
    [ALLEEG, EEG, ~, ALLCOM] = eeglab;
    FileName = Files(i).name;
    cd(datafolder);
    thisfilename = FileName;
    Ind1 = find(FileName == '.');
    basefilename = thisfilename(1:Ind1-1);
%     EEG = pop_loadbv([datafolder '/'], thisfilename, [],[]);
    EEG = pop_loadset('filename', thisfilename, 'filepath', datafolder);
    EEG.setname= [basefilename '.set'];
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'setname',EEG.setname,'gui','off');
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);
    eeglab redraw
    counter = 1;
    
%     %Now define channel locations
%       EEG=pop_chanedit(EEG, 'lookup', ...
%                 '/Users/timothymorris/MATLAB-Drive/eeglab2019_0/plugins/dipfit3.0/standard_BESA/standard-10-5-cap385.elp');

    %Now define channel locations -->IF DEFINED BEFORE ASR, CHANNEL
%     CORRELATION WILL BE COMPARED TO RANSAC RECONSTRUCTION. IF DEFINED
%     AFTER, A COMMON CORRELATION BETWEEN ALL CHANNELS WILL BE USED
    EEG=pop_chanedit(EEG, 'lookup', ...
        'C:\Users\jross4\Desktop\eeglab-develop\plugins\dipfit3.3\standard_BESA\standard-10-5-cap385.elp');

    
    % now delete EMG channels IF THEY EXIST
    EEG = pop_select(EEG, 'nochannel', [64:EEG.nbchan]);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    eeglab redraw
     
%     DS = 5000;
    % Now downsample to 1000Hz
    DS = 1000;
    EEG = pop_resample(EEG,DS);

    %Save resulting dataset
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG.setname = [basefilename '_A01.set'];
    pop_saveset( EEG, 'filename',EEG.setname,'filepath', datafolder);
    
    %Bandpass filter
    EEG = tesa_filtbutter( EEG, bandlow, bandhigh, filtord, 'bandpass' );
    
    % Run clean_rawdata (ASR) 
    originalEEG = EEG; % copy EEG before clean_raw?
%     EEG = clean_rawdata(EEG, 5, [0.25 0.75], 0.8, 4, 5, 0.5 );  % Run clean_raw
    EEG = clean_rawdata(EEG, 5, 'off', 0.8, 'off', 5, 0.5 );  % Run clean_raw NEW!!!!!!!!!!!!!!!!
    Channels_Removed = setdiff({originalEEG.chanlocs.labels},{EEG.chanlocs.labels}, 'stable');  % Make a variable to show what's different between the original and new channel labels
    save([basefilename '_Channels_Removed.mat'], 'Channels_Removed')
 
% %  ------------   
% % MANUAL OVERRIDE IN CASE TOO MANY CHANNELS ARE DELETED BY ASR:
%     EEG = originalEEG;
%     eeglab redraw
% Channels_Removed = [{'AF7'}]
% save([basefilename '_Channels_Removed.mat'], 'Channels_Removed')
% %  ------------
    
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG )
   
    eeglab redraw 
    
    
    %interpolate missing channels 
    EEG = pop_interp(EEG, originalEEG.chanlocs, 'spherical');% Interpolate channels removed. 
    
    
    %And average reference (if desired)%     
    if Rereference == 1
        reference_type = 'average';
        EEG = pop_reref(EEG, []);
    end

    % Introduce event markers
    EEG.event=[];
    for j=1:round(EEG.pnts/(EpochLength*EEG.srate))
        EEG.event(j).latency=EpochLength*EEG.srate*j;
        EEG.event(j).duration=0;
        EEG.event(j).type='epochmark';
        EEG.event(j).source = basefilename;
        EEG.event(j).eventnumber = j;
    end

    % Epoch data
    newfilename = [basefilename 'A02.set']
    EEG = pop_epoch( EEG, { eventmarker }, [epoch_begin  epoch_end], 'newname', newfilename, 'epochinfo', 'yes');
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );
  
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG.setname = [basefilename '_A02.set'];
    pop_saveset( EEG, 'filename',EEG.setname,'filepath', datafolder);
   
    eeglab redraw
  
    
    % -------------------------------------------------------------
    %               6a)Perform ICA - fastica method
    % -------------------------------------------------------------
    if strcmp(ICAmethod,'fast ICA')==1
        if DoPCA == 1
            if CalculateDimensions == 1
                PCAdimensions = [];
                PCAdimensions = fcn_EstimateNrICAComp(EEG, PercentVar, MinComp);
                str = ['For file ' basefilename ', ' num2str(PCAdimensions) 'components are necessary to keep 99%% of the variance in the data'];
                disp(str);
                EEG = pop_runica(EEG,'icatype','fastica', 'approach', 'symm', 'g', ...
                'tanh', 'firsteig', 1, 'lasteig', PCAdimensions); %Does fastica with fixed PCA decomposition first
            else
                EEG = pop_runica(EEG,'icatype','fastica', 'approach', 'symm', 'g', ...
                'tanh', 'firsteig', 1, 'lasteig', PCAdimensions); %Does fastica with fixed PCA decomposition first
            end
            EEG.ica2_dimensions = PCAdimensions;
        else
            EEG = pop_runica(EEG,'icatype','fastica', 'approach', 'symm', 'g', 'tanh');
        end

        EEG = eeg_checkset( EEG );
        EEG.BadCmp=[];
        if DoPCA==1
            EEG.setname = [basefilename '_A03_fICA' num2str(PCAdimensions) '.set'];
        else
            EEG.setname = [basefilename '_A03_fICA.set'];
        end
        
        EEG = pop_saveset( EEG, 'filename',EEG.setname,'filepath',datafolder);
        EEG = eeg_checkset( EEG );
    end
    
    % -------------------------------------------------------------
    %               6b)Perform ICA - runica method
    % -------------------------------------------------------------
    if strcmp(ICAmethod,'infomax')==1
        if DoPCA == 1
            if CalculateDimensions == 1
                PCAdimensions = [];
                PCAdimensions = fcn_EstimateNrICAComp_RunPCA_Var(EEG, PercentVar);
                EEG = pop_runica(EEG, 'extended', 1, 'pca', PCAdimensions); %Does fastica with PCA decomposition first
            else
                EEG = pop_runica(EEG, 'extended', 1, 'pca', PCAdimensions); %Does fastica with PCA decomposition first
            end
        else
            EEG = pop_runica(pop_runica(EEG, 'extended',1));
        end

        EEG = eeg_checkset( EEG );
        EEG.BadCmp=[];
        if DoPCA==1
            EEG.setname = [basefilename '_A03_rICA' num2str(PCAdimensions) '.set'];
        else
            EEG.setname = [basefilename '_A03_rICA.set'];
        end
        EEG = pop_saveset( EEG, 'filename',EEG.setname,'filepath',datafolder);
        EEG = eeg_checkset( EEG );
    end
    
    string = ['Finished running ICA on ' FileName];
    disp(string);
   
    eeglab redraw
   
    
%     varRem = 99;
%     for k  = 1:length(EEG.varsPerc)
%         if varRem-EEG.varsPerc(k) > 0
%             sigComps(k) = k;
%         end
%         varRem = varRem-EEG.varsPerc(k);
%     end
%     numSigComps = length(sigComps);
    
    % Run MARA
    % Note: Select "plot and select components for removal" and select OK
    % on all plots, then type dbcont in command window
%     EEG = pop_processMARA ( ALLEEG,EEG,CURRENTSET );
%     
%     [EEG, EEG.varsPerc] = tesa_sortcomps(EEG);
%     allcomps = length(EEG.varsPerc);    
    
    [ALLEEG,EEG] = processMARA ( ALLEEG,EEG,CURRENTSET );
    EEG = eeg_checkset( EEG );
    
    pause(2)
    
    MARAmarked = EEG.reject.gcompreject;
    MARAmarked = find(MARAmarked); 
    numMARAmarked = length(MARAmarked);
    
    % Set threshold to 30%
    EEG.reject.gcompreject = zeros(size(EEG.reject.gcompreject)); 
    EEG.reject.gcompreject(EEG.reject.MARAinfo.posterior_artefactprob > 0.30) = 1;
    
    AUTOdelete = EEG.reject.gcompreject;
    AUTOdelete = find(AUTOdelete);

    % Remove selected components
    EEG = pop_subcomp(EEG,[],0);
    eeglab redraw
    
    % Convert AUTOdelete to a string
    AUTOdelete_str = num2str(AUTOdelete);
    
    % Save everything
    fid = fopen(file_name,'a+');
    fprintf(fid,'\n %s %s %d %s %d %s %d %s %d %s %d %s %d %s %s %s %d %s %d %s %s %s %s' ,basefilename,',',bandlow,',',bandhigh,',',EpochLength,',',DoPCA,',',DS,',',length(Channels_Removed),',',reference_type,',',PCAdimensions,',',length(AUTOdelete),',','[',AUTOdelete_str,']');

    % Save
    EEG = eeg_checkset( EEG );
%     [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG.setname = [basefilename '_A04_MARA.set'];
    pop_saveset( EEG, 'filename',EEG.setname,'filepath', datafolder);
%     clear('EEG')
    
end 

beep
toc

% Save run time
fid = fopen(file_name,'a+');
fprintf(fid,'\n %s %s %f' ,'Runtime (seconds)',',',toc);
eeglab redraw
fclose('all');  

%%
% % =======================================================================
% % Step 2: Component selection using TESA 
% % ======================================================================== 
clear all; close all; clc; 
tic
homefolder = 'D:\SAGESII_P034_GB\SAGESII_P034_GB_V2_1_14_2021\Analysis\Andrei_training\Resting';
datafolder = 'D:\SAGESII_P034_GB\SAGESII_P034_GB_V2_1_14_2021\Analysis\Andrei_training\Resting';
cd(datafolder);

% Open summary file for output
summ_path = strcat('auto_output/','auto_pipeline_summary');
file_name = [summ_path sprintf('%d') '.txt'];
fopen(file_name,'a+');

% Make header
fid = fopen(file_name,'a+');
fprintf(fid,'\n \n %s %s %s %s %s' ,'file name',',','number of manually deleted ICs using TESA',',','number of ICs remaining');
  

Files = dir('*_A04_MARA.set');

for i = 1 : length(Files)
    FileName = Files(i).name;
    
    clear EEG ALLEEG CURRENTSET ALLCOM;
    [ALLEEG, EEG, CURRENTSET ALLCOM] = eeglab;
    EEG = pop_loadset('filename', FileName, 'filepath', datafolder);
    Ind1 = strfind(FileName, '_A04_');
    basefilename = FileName(1:Ind1-1);
    eeglab redraw
    
    %First sort components by percent of variance explained
    [EEG, EEG.varsPerc] = tesa_sortcomps(EEG);
    num_comps = length(EEG.varsPerc);
    
%     varRem = 99;
%     for k  = 1:length(EEG.varsPerc)
%         if varRem-EEG.varsPerc(k) > 0
%             sigComps(k) = k;
%         end
%         varRem = varRem-EEG.varsPerc(k);
%     end
%     numSigComps = length(sigComps);

    %Open Component activations in eegplot window
	  pop_eegplot(EEG, 0, 1, 0, [], 'winlength', 5, 'dispchans', 5);
    
    %Create filename for after component marking
    thisfile = [basefilename '_A05b_MARAmarked.set'];
    EEG.setname= thisfile; EEG.filepath = datafolder;

    %Now call TESA for automatic component selection
    %Does TMS pulse, blink, lateral eye movements, muscle and electrode
    EEG = tesa_compselect( EEG,'compCheck','on','figSize','medium','plotTimeX',[0 2998],'plotFreqX',[1 100], ...
        'tmsMuscle','off', ...
        'blink','on','blinkThresh',2.5, 'blinkElecs',{'Fp1','Fp2'}, 'blinkFeedback','off',...
        'move','on','moveThresh',2, 'moveElecs',{'F7','F8'},'moveFeedback','off', ...
        'muscle','on','muscleThresh',0.6,'muscleFreqWin',[30 100],'muscleFreqIn',[7,75],'muscleFreqEx',[48,52],'muscleFeedback','off',...
        'elecNoise','on','elecNoiseThresh',4,'elecNoiseFeedback','off' );
    %'muscleFreqIn',[7,75]
    %'muscleFreqEx',[48,52],
    
    comps_final = length(EEG.icaweights(:,1));
    comps_removed = num_comps - comps_final;
    
    % Save everything
    fid = fopen(file_name,'a+');
    fprintf(fid,'\n %s %s %d %s %d' ,basefilename,',',comps_removed,',',comps_final);  
    EEG.filepath = []; %Removes filepath so it is not stored for later   
    thisfile = [basefilename '_A05_TESA.set'];
    EEG.setname= thisfile;
    EEG = pop_saveset( EEG, 'filename', thisfile, 'filepath',datafolder);

end

beep
toc

% Save run time
fid = fopen(file_name,'a+');
fprintf(fid,'\n %s %s %f' ,'Runtime (seconds)',',',toc);

fclose('all');

%%
%=======================================
%Step 3 Plot data before and after 
%=======================================

clear all; close all;

homefolder = 'D:\SAGESII_P034_GB\SAGESII_P034_GB_V2_1_14_2021\Analysis\Andrei_training\Resting';
datafolder = 'D:\SAGESII_P034_GB\SAGESII_P034_GB_V2_1_14_2021\Analysis\Andrei_training\Resting';

FilesPre = dir('*_A02.set');
FilesPost = dir('*_A05_TESA.set');

for p = 1:length(FilesPost)
    cd(datafolder);
    [ALLEEG, EEG, ~, ALLCOM] = eeglab;
    
    FileNamePre = FilesPre(p).name;
    thisfilename = FileNamePre;
    Ind1Pre = find(FileNamePre == '.');
    basefilenamePre = thisfilename(1:Ind1Pre-1);
    EEG = pop_loadset('filename', thisfilename, 'filepath', datafolder);
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);
    eeglab redraw
    
    FileNamePost = FilesPost(p).name;
    thisfilename = FileNamePost;
    Ind1Post = find(FileNamePost == '.');
    basefilenamePost = thisfilename(1:Ind1Post-1);
    EEG = pop_loadset('filename', thisfilename, 'filepath', datafolder);
    EEG = eeg_checkset( EEG );
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);
    eeglab redraw
    
    plot_title = strcat(basefilenamePost,': Original (black) and cleaned (red)');
    eegplot([ALLEEG(1,1).data], 'data2', [ALLEEG(1,2).data], 'title', plot_title,'eloc_file', EEG.chanlocs,'color','off')
%     n = get(gcf,'Number');
%     file_fig = strcat(basefilenamePost,'_pre_and_post.fig');
end

% NOTE: These do not automatically save!!! 
%%
%=======================================
%Step 4 Frequency calculations and plots
%=======================================

% Calculates sum of frequency band absolute and relative power within each frequency range, and peak alpha frequency and power from posterior electrodes.
% Output includes 3 plots and summary appended to frequency_band_power.txt.
% 
% window size = sampling rate
% overlap = 50%
%   TRT: FS = 5000 Hz with epochs 3 seconds long
%
% NOTE: Make sure your posterior electrode names match the electrode
% numbers listed below (lines 565-573).

clear all; close all; clc; 

homefolder = 'D:\SAGESII_P034_GB\SAGESII_P034_GB_V2_1_14_2021\Analysis\Andrei_training\Resting';
datafolder = 'D:\SAGESII_P034_GB\SAGESII_P034_GB_V2_1_14_2021\Analysis\Andrei_training\Resting';
data_loc =  datafolder;

% addpath(genpath(homefolder));

% Get files
cd(data_loc);
Files = dir('*_A05_TESA.set');

% % Remove hidden files
% Files = Files(arrayfun(@(x) ~strcmp(x.name(1),'.'),Files));
cd(homefolder)
% Make a file for spectral analysis
mkdir('spect_analysis');
addpath(genpath('spect_analysis/'));

% Open summary file for output
summ_path = strcat('spect_analysis/','frequency_band_power');
file_name = [summ_path sprintf('%d') '.txt'];
fopen(file_name,'a+');

cd(data_loc);
for i = 1:length(Files)
    [ALLEEG, EEG, ~, ALLCOM] = eeglab;
    FileName = Files(i).name; 
    EEG = pop_loadset('filename', FileName, 'filepath', data_loc);
    Ind1 = strfind(FileName, '_A05_TESA.set');
    basefilename = FileName(1:Ind1-1);
    eeglab redraw
    
    % Frequency domain parameters
    nwinsize = EEG.srate;
    noverlap = (ceil(nwinsize * 0.50)); % 50% overlap in calculating power
    
    figure(2)
    [spectra2,freqs2,speccomp2,contrib2,specstd2]=spectopo(EEG.data, size(EEG.data,2), ...
        EEG.srate, 'winsize', nwinsize, 'overlap', noverlap, 'chanlocs', EEG.chanlocs, 'freqrange', [1 50]);
    hold on
%     axis([0 30 -20 30]);
    title(basefilename,'interpreter','none');
    hold off
    
    file_fig1_m = strcat('spect_analysis/',basefilename,'_spect1');
    file_fig1_jpg = strcat('spect_analysis/',basefilename,'_spect1.jpg');
%     saveas(2,file_fig1_m);
    saveas(2,file_fig1_jpg);

    % Power and frequencies  
    Results.powerspectra = db2pow(spectra2);
    Results.powerfreqs = freqs2;
    whole48 = Results.powerspectra(:,1:48)';
    sum48 = sum(whole48);
    avgPSD = mean(sum48);
    
    figure(3)
    plot(Results.powerfreqs,Results.powerspectra)
    hold on
    xlim([1 50]);
    xlabel('Frequency');
    ylabel('Power');
    title(basefilename,'interpreter','none');
    hold off
    
    file_fig2_m = strcat('spect_analysis/',basefilename,'_spect2.fig');
    file_fig2_jpg = strcat('spect_analysis/',basefilename,'_spect2.jpg');
%     saveas(3,file_fig2_m);
    saveas(3,file_fig2_jpg);

    % ABSOLUTE, all electrodes: Find delta=1-4, theta=4-8, alpha=8-13,
    % beta=13-30, gamma=30-48
    deltaA = Results.powerspectra(:,1:4)';
    sumdeltaA = sum(deltaA);
    deltaPowerA = mean(sumdeltaA);
    
    thetaA = Results.powerspectra(:,4:8)';
    sumthetaA = sum(thetaA);
    thetaPowerA = mean(sumthetaA);
    
    alphaA = Results.powerspectra(:,8:13)';
    sumalphaA = sum(alphaA);
    alphaPowerA = mean(sumalphaA);
    
    betaA = Results.powerspectra(:,13:30)';
    sumbetaA = sum(betaA);
    betaPowerA = mean(sumbetaA);
    
    gammaA = Results.powerspectra(:,30:48)';
    sumgammaA = sum(gammaA);
    gammaPowerA = mean(sumgammaA);
    
    % RELATIVE: Find delta=1-4, theta=4-8, alpha=8-13, beta=13-30,
    % gamma=30-48
    deltaPowerR = mean(deltaPowerA)/avgPSD;
    thetaPowerR = mean(thetaPowerA)/avgPSD;
    alphaPowerR = mean(alphaPowerA)/avgPSD;
    betaPowerR = mean(betaPowerA)/avgPSD;
    gammaPowerR = mean(gammaPowerA)/avgPSD;
    
    % Find the frequency of the alpha peak from posterior electrodes only
    % CHECK TO MAKE SURE THESE ELECTRODE NUMBERS ARE CORRECT! 
    % SAGES: (Pz=13, P3=14, O1=16, Oz=17, O2=18, P4=19, P1=44, PO3=47, POz=48,
    % PO4=49, P2=52)
    alphaI = Results.powerspectra(:,8:13);
    alphaIpost = alphaI([13,14,16,17,18,19,44,47,48,49,52],:);
    postalpha = mean(alphaIpost);
    alphaFreqAll = Results.powerfreqs(8+1:13+1);
    [M,I] = max(postalpha);
    alphaPeakPower = M;
    alphaPeakFreq = alphaFreqAll(I);
    
    figure(4) 
    plot(Results.powerfreqs,Results.powerspectra([13,14,16,17,18,19,44,47,48,49,52],:))
    hold on
    xlim([1 50]);
    xlabel('Frequency');
    ylabel('Power');
    title_post = strcat(basefilename,': Posterior electrodes only');
    title(title_post,'interpreter','none');
    hold off
    
    file_fig4_m = strcat('spect_analysis/',basefilename,'_spect2_posterior.fig');
    file_fig4_jpg = strcat('spect_analysis/',basefilename,'_spect2_posterior.jpg');
%     saveas(3,file_fig2_m);
    saveas(4,file_fig4_jpg);
    
    % Save everything
    fid = fopen([homefolder '\' file_name],'a+');
    fprintf(fid,'\n \n %s %s %s %s %f',basefilename,',','delta range (1-4 Hz) absolute power',',',deltaPowerA);
    fprintf(fid,'\n %s %s %s %s %f',basefilename,',','theta range (4-8 Hz) absolute power',',',thetaPowerA);
    fprintf(fid,'\n %s %s %s %s %f',basefilename,',','alpha range (8-13) absolute power',',',alphaPowerA);
    fprintf(fid,'\n %s %s %s %s %f',basefilename,',','beta range (13-30) absolute power',',',betaPowerA);
    fprintf(fid,'\n %s %s %s %s %f',basefilename,',','gamma range (30-48 Hz) absolute power',',',gammaPowerA);
    
    fprintf(fid,'\n %s %s %s %s %f',basefilename,',','delta range (1-4 Hz) relative power',',',deltaPowerR);
    fprintf(fid,'\n %s %s %s %s %f',basefilename,',','theta range (4-8 Hz) relative power',',',thetaPowerR);
    fprintf(fid,'\n %s %s %s %s %f',basefilename,',','alpha range (8-13) relative power',',',alphaPowerR);
    fprintf(fid,'\n %s %s %s %s %f',basefilename,',','beta range (13-30) relative power',',',betaPowerR);
    fprintf(fid,'\n %s %s %s %s %f',basefilename,',','gamma range (30-48 Hz) relative power',',',gammaPowerR);
    
    fprintf(fid,'\n %s %s %s %s %d',basefilename,',','alpha peak frequency (from average of posterior electrodes)', ',',alphaPeakFreq);
    fprintf(fid,'\n %s %s %s %s %f',basefilename,',','alpha peak absolute power (from average of posterior electrodes)',',', alphaPeakPower);
%     close(2);
%     close(3);
%     
%     clear ALLEEG EEG;

eeglab redraw

end
  
fclose('all');  
















%% EXTRA OPTIONAL STEPS===============================

%%====================================================
%STEP 1.25: Manual reject channels if needed
%=====================================================

%If during the inspection of the frequency plots bad channels still exsit,
%edit this step to remove those channels. Enter the channel numbers in the bad_chan variable below. This step will then overwrite the
%previous file and you will need to go back and re run step 3. 

bad_chan = [50];

EEG = pop_select(EEG,'nochannel',bad_chan);
EEG = pop_interp(EEG, ALLEEG(1).chanlocs, 'spherical');

EEG = eeg_checkset(EEG);
eeglab redraw
   
thisfile = [basefilename '_A05_TESA.set'];
EEG.setname= thisfile;
EEG = pop_saveset( EEG, 'filename', thisfile, 'filepath',datafolder);

clear all; close all;

%%
%====================================================
%STEP 1.5: Manual epoch rejection and ICA with MARA
%====================================================

%IF MANUAL EPOCH REJECTION REQUIRED

%Run this step if after running intial step 1, there is epochs that need to
%be deleted manually (i.e ASR did not clean the data sufficiently). Remove
%.set files that were created in first step (A03_fICA30, A04_MARA and
%A05_TESA and save elsewhere for comparison). After running manually
%component selction, this step will re-run ICA with MARA saving the file
%that can be loaded into step 2. 


clear all; close all; clc; 

homefolder = '/Users/timothymorris/MATLAB-Drive';
datafolder = '/Users/timothymorris/Desktop/BBHI_TMSEEG/Restingstate';
data_loc =  datafolder;

% addpath(genpath(homefolder));

% % Open summary file for output
% summ_path = strcat('auto_output/','auto_pipeline_summary');
% file_name = [summ_path sprintf('%d') '.txt'];
% fopen(file_name,'a+');
% 
% % Make header
% fid = fopen(file_name,'a+');
% fprintf(fid,'\n \n %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s' ,'file name',',','high pass (Hz)',',','low pass (Hz)',',','epoch length',',','PCA (0=no, 1=yes)',',','downsampled to (Hz)',',','number of channels removed',',','reference',',','number of ICs before MARA',',','number of MARA comps deleted');
%   

DoPCA = 1; %VAR: Whether to do dimensionality reduction w PCA prior to ICA: 0 for no, 1 for yes
CalculateDimensions = 0; %VAR: Whether to calculate number of dimensions to reduce data to prior to ICA
PercentVar = 99; %VAR: If calculating dimensions, percent of variance to explain in PCA reduction
MinComp = 20; %Minimum number of components to include if calculating variance
PCAdimensions = 30; %VAR: Number of dimensions to reduce to if doing PCA but not calculating dimensions
ICAmethod = 'fast ICA'; %VAR: ICA type
% ICAmethod = 'infomax';

% Get files
cd(data_loc);
Files = dir('*_A02.set');


sinchan_prob_thresh = 3.5; %VAR: Single-channel threshold for rejection based
                            % on probability; default use 3.5
allchan_prob_thresh = 3; %VAR: All-channel threshold for rejection based
                            % on probability; default use 3
sinchan_kurt_thresh = 5; %VAR: Single-channel threshold for rejection based
                            % on kurtosis; default use 5
allchan_kurt_thresh = 3; %VAR: All-channel threshold for rejection based
                            % on kurtosis; default use 3
voltage_reject = 1; %VAR: Whether to reject based on a voltage threshold;
                        % Use 1 for yes, 0 otherwise
voltage_thresh =  100; %VAR: Voltage threshold for epoch rejection
ChansExclude = {'EOG1', 'EOG2' 'Fp1', 'Fpz', 'Fp2', 'AF7', 'AF8'}; %VAR: %Channels to exclude from voltage thresholding - Brainvision
% ChansExclude = {'EOG', 'FP1', 'FPZ', 'FP2', 'AF1', 'AFZ', 'AF2'}; %VAR: Channels to exclude from voltage thresholding - Nexstim
BeginTimeInclude = -0.5; %VAR: Beginning of time period to do voltage thresholding on
EndTimeInclude = 1; %VAR: End of time period to do voltage thresholding on
BeginTimeExclude = 0; %VAR: Beginning of time period to exclude from voltage thresholding
EndTimeExclude = 0.05; %VAR: End of time period to exclude from voltage thresholding (in seconds)


% %remove hidden files
% Files = Files(arrayfun(@(x) ~strcmp(x.name(1),'.'),Files));

for i = 1 : length(Files)
    FileName = Files(i).name;
    clear EEG ALLEEG CURRENTSET ALLCOM count2 chanind;
    [ALLEEG, EEG, CURRENTSET ALLCOM] = eeglab;
    EEG = pop_loadset('filename', FileName, 'filepath', datafolder);
    Ind1 = strfind(FileName, '_A02');
    basefilename = FileName(1:Ind1-1);
    eeglab redraw
    
% ------------------------------------------------------------------------
%           5a) Copy original data into new variable, then bandpass and
%           average reference data. Note that bandpassing and rereferecing 
%           is strictly for visualization to help determine which epochs to 
%           delete, and the original unfiltered data is replaced at the end
%           of this step
% ------------------------------------------------------------------------ 
    EEG.origdata = EEG.data;

    filtord = 4;
    stoplow = 57; %58?
    stophigh = 63;
    bandlow = 1;    %Lower edge of bandpass filter
    bandhigh = 50;   %Upper edge of bandpass filter
    % -------------- Filtering
    % Uses fourth-order Butterworth filter
    % Backwards and forwards filtering
    %Bandstop (notch) filter
    EEG = tesa_filtbutter( EEG, stoplow, stophigh, filtord, 'bandstop' );
    %Bandpass filter
    EEG = tesa_filtbutter( EEG, bandlow, bandhigh, filtord, 'bandpass' );
    % -------------- Average Referencing
    EEG = pop_reref(EEG, []);
    
% ------------------------------------------------------------------------
%           5b) Tag trials based on amplitude, probability, and kurtosis  
% ------------------------------------------------------------------------ 
    %Use EEGLAB function to calcuate Kurtosis and tag bad channels
    search_array = [];
    
    % Now identify abnormal trials based on probability and kurtosis
    EEG = pop_jointprob(EEG,1,[1:size(EEG.chanlocs,2)],sinchan_prob_thresh,allchan_prob_thresh,1,0);
    EEG = pop_rejkurt(EEG,1,[1:size(EEG.chanlocs,2)],sinchan_kurt_thresh,allchan_kurt_thresh,1,0);
    
    % Below will also identify trials with activity greater than the
    % voltage threshold, but will exclude eye channels
    if (voltage_reject==1)
        for count = 1 : EEG.nbchan
            electrodes{count} = EEG.chanlocs(count).labels;
            search_array = [search_array count];
        end
        
        for count2 = 1 : length(ChansExclude);
            tempchan = find(strcmp(electrodes, ChansExclude{count2}));
            if ~isempty(tempchan) 
                chanind(count2) = tempchan;
                search_array = search_array(search_array~=chanind(count2));
            end
        end
        

        EEG = pop_eegthresh(EEG,1,search_array,-voltage_thresh,voltage_thresh,...
            [BeginTimeInclude EndTimeExclude],[BeginTimeExclude EndTimeInclude],1,0);
    end
            
% ------------------------------------------------------------------------
%  5c)  Scroll, visualize tagged trials, manually tag trials with lots of 
%       muscle and other artifacts, untag what is not noisy, finalize
%       trials to delete, and save resulting dataset
% ------------------------------------------------------------------------ 
    eeglab redraw;
    disp('Visually select trials with lots of muscle and other artifacts, using pop_rejmenu');
    disp('Update marked trials, but do NOT delete. Then type dbcont when done');
    pop_rejmenu(EEG,1);
    keyboard;

%             Run pop_rejmenu from the GUI here, update marked trials, but
%             do NOT delete yet (needs to be saved)

    EEG.badtr = union(union(union(find(EEG.reject.rejmanual>0), ...
        find(EEG.reject.rejjp>0)), union(find(EEG.reject.rejkurt>0), ...
        find(EEG.reject.rejthresh>0))),find(EEG.reject.rejconst>0));
    EEG.setname = [basefilename '_marked.set'];
    EEG = pop_saveset( EEG, 'filename',[basefilename '_A03b_marked.set'],...
            'filepath',datafolder);

% ------------------------------------------------------------------------
%   5d) Delete Bad Trials
% ------------------------------------------------------------------------
            
            EEG.data = EEG.origdata;  %Copy original unfiltered data back
            EEG.origdata = [];
            EEG = pop_rejepoch( EEG, EEG.badtr ,0); % EEGLAB function to delete all tagged trials
            EEG = eeg_checkset( EEG );
            EEG.data=double(EEG.data);

    thisfile = [basefilename '_A03_ClEp.set'];
    EEG.setname= thisfile;
    EEG = pop_saveset( EEG, 'filename', thisfile, 'filepath',datafolder);

    eeglab redraw 

%

 % -------------------------------------------------------------
    %               6a)Perform ICA - fastica method
    % -------------------------------------------------------------
    if strcmp(ICAmethod,'fast ICA')==1
        if DoPCA == 1
            if CalculateDimensions == 1
                PCAdimensions = [];
                PCAdimensions = fcn_EstimateNrICAComp(EEG, PercentVar, MinComp);
                str = ['For file ' basefilename ', ' num2str(PCAdimensions) 'components are necessary to keep 99%% of the variance in the data'];
                disp(str);
                EEG = pop_runica(EEG,'icatype','fastica', 'approach', 'symm', 'g', ...
                'tanh', 'firsteig', 1, 'lasteig', PCAdimensions); %Does fastica with fixed PCA decomposition first
            else
                EEG = pop_runica(EEG,'icatype','fastica', 'approach', 'symm', 'g', ...
                'tanh', 'firsteig', 1, 'lasteig', PCAdimensions); %Does fastica with fixed PCA decomposition first
            end
            EEG.ica2_dimensions = PCAdimensions;
        else
            EEG = pop_runica(EEG,'icatype','fastica', 'approach', 'symm', 'g', 'tanh');
        end

        EEG = eeg_checkset( EEG );
        EEG.BadCmp=[];
        if DoPCA==1
            EEG.setname = [basefilename '_A03_fICA' num2str(PCAdimensions) '.set'];
        else
            EEG.setname = [basefilename '_A03_fICA.set'];
        end
        
        EEG = pop_saveset( EEG, 'filename',EEG.setname,'filepath',datafolder);
        EEG = eeg_checkset( EEG );
    end
    
    % -------------------------------------------------------------
    %               6b)Perform ICA - runica method
    % -------------------------------------------------------------
    if strcmp(ICAmethod,'infomax')==1
        if DoPCA == 1
            if CalculateDimensions == 1
                PCAdimensions = [];
                PCAdimensions = fcn_EstimateNrICAComp_RunPCA_Var(EEG, PercentVar);
                EEG = pop_runica(EEG, 'extended', 1, 'pca', PCAdimensions); %Does fastica with PCA decomposition first
            else
                EEG = pop_runica(EEG, 'extended', 1, 'pca', PCAdimensions); %Does fastica with PCA decomposition first
            end
        else
            EEG = pop_runica(pop_runica(EEG, 'extended',1));
        end

        EEG = eeg_checkset( EEG );
        EEG.BadCmp=[];
        if DoPCA==1
            EEG.setname = [basefilename '_A03_rICA' num2str(PCAdimensions) '.set'];
        else
            EEG.setname = [basefilename '_A03_rICA.set'];
        end
        EEG = pop_saveset( EEG, 'filename',EEG.setname,'filepath',datafolder);
        EEG = eeg_checkset( EEG );
    end
    
    string = ['Finished running ICA on ' FileName];
    disp(string);
   
    eeglab redraw
   
    
%     varRem = 99;
%     for k  = 1:length(EEG.varsPerc)
%         if varRem-EEG.varsPerc(k) > 0
%             sigComps(k) = k;
%         end
%         varRem = varRem-EEG.varsPerc(k);
%     end
%     numSigComps = length(sigComps);
    
    % Run MARA
    % Note: Select "plot and select components for removal" and select OK
    % on all plots, then type dbcont in command window
%     EEG = pop_processMARA ( ALLEEG,EEG,CURRENTSET );
%     
%     [EEG, EEG.varsPerc] = tesa_sortcomps(EEG);
%     allcomps = length(EEG.varsPerc);    
    
    [ALLEEG,EEG] = processMARA ( ALLEEG,EEG,CURRENTSET );
    EEG = eeg_checkset( EEG );
    
    pause(2)
    
    MARAmarked = EEG.reject.gcompreject;
    MARAmarked = find(MARAmarked); 
    numMARAmarked = length(MARAmarked);
    
    % Set threshold to 30%
    EEG.reject.gcompreject = zeros(size(EEG.reject.gcompreject)); 
    EEG.reject.gcompreject(EEG.reject.MARAinfo.posterior_artefactprob > 0.30) = 1;
    
    AUTOdelete = EEG.reject.gcompreject;
    AUTOdelete = find(AUTOdelete);

    % Remove selected components
    EEG = pop_subcomp(EEG,[],0);
    eeglab redraw
    
%     Save everything
%     fid = fopen(file_name,'a+');
%     fprintf(fid,'\n %s %s %d %s %d %s %d %s %d %s %d %s %d %s %s %s %d %s %d' ,basefilename,',',bandlow,',',bandhigh,',',EpochLength,',',DoPCA,',',DS,',',length(Channels_Removed),',',reference_type,',',PCAdimensions,',',length(AUTOdelete));

    % Save
    EEG = eeg_checkset( EEG );
%     [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG.setname = [basefilename '_manepoch_A04_MARA.set'];
    pop_saveset( EEG, 'filename',EEG.setname,'filepath', datafolder);
    clear('EEG')
    
end 


% beep
% toc
% 
% % Save run time
% fid = fopen(file_name,'a+');
% fprintf(fid,'\n %s %s %f' ,'Runtime (seconds)',',',toc);
% eeglab redraw
% fclose('all');  

%%

% %Do manual component seleciton here if needed 
% 
%    %First sort components by percent of variance explained
%     [EEG, EEG.varsPerc] = tesa_sortcomps(EEG);
% 
%     %Open Component activations in eegplot window
% 	  pop_eegplot(EEG, 0, 1, 0, [], 'winlength', 5, 'dispchans', 5);
%     
%     %Create filename for after component marking
%     thisfile = [basefilename '_S04b_manmarked.set'];
%     EEG.setname= thisfile; EEG.filepath = datafolder;
% 
%       %Now call TESA for automatic component selection
%     %Does TMS pulse, blink, lateral eye movements, muscle and electrode
%     EEG = tesa_compselect( EEG,'compCheck','on','figSize','medium','plotTimeX',[0 2998],'plotFreqX',[1 100], ...
%         'tmsMuscle','off', ...
%         'blink','on','blinkThresh',2.5, 'blinkElecs',{'Fp1','Fp2'}, 'blinkFeedback','off',...
%         'move','on','moveThresh',2, 'moveElecs',{'F7','F8'},'moveFeedback','off', ...
%         'muscle','on','muscleThresh',0.6,'muscleFreqWin',[30 100],'muscleFeedback','off',...
%         'elecNoise','on','elecNoiseThresh',4,'elecNoiseFeedback','off' );
%     %'muscleFreqIn',[7,75]
%     %'muscleFreqEx',[48,52],
%     EEG.filepath = []; %Removes filepath so it is not stored for later   
%     thisfile = [basefilename '_S04_ICA_MAN.set'];
%     EEG.setname= thisfile;
%     EEG = pop_saveset( EEG, 'filename', thisfile, 'filepath',datafolder);

%% 
% % =======================================================================
% % Step 6.5 Artifact rejection using MARA for ICs with greater than 30%
% % artifact probabilities
% %     05/2019, jross4@bidmc.harvard.edu, tpmorris@bidmc.harvard.edu
% % =======================================================================
% 
% clear all; close all; clc; 
% homefolder = '/Users/timothymorris/MATLAB-Drive';
% datafolder = '/Users/timothymorris/Desktop/BBHI_TMSEEG/EEGConcuss';
% savefolder = '/Users/timothymorris/Desktop/BBHI_TMSEEG/EEGConcuss/results/';
% cd(datafolder);
% 
% % Open summary file for output
% summ_path = strcat(savefolder,'MAN_AUTO_comparison');
% file_name = [summ_path sprintf('%d') '.txt'];
% fopen(file_name,'a+');
% 
% Files = dir('*S04_fICA*'); 
% for i = 1 : length(Files)
%     FileName = Files(i).name;
%     
%     clear EEG ALLEEG CURRENTSET ALLCOM;
%     [ALLEEG, EEG, CURRENTSET ALLCOM] = eeglab;
%     EEG = pop_loadset('filename', FileName, 'filepath', datafolder);
%     Ind1 = strfind(FileName, '_S04_');
%     basefilename = FileName(1:Ind1-1);
%     eeglab redraw
%     
%     [EEG, EEG.varsPerc] = tesa_sortcomps(EEG);
%     allcomps = length(EEG.varsPerc);
%     
% %     varRem = 99;
% %     for k  = 1:length(EEG.varsPerc)
% %         if varRem-EEG.varsPerc(k) > 0
% %             sigComps(k) = k;
% %         end
% %         varRem = varRem-EEG.varsPerc(k);
% %     end
% %     numSigComps = length(sigComps);
%     
%     % Run MARA
%     % Note: Select "plot and select components for removal" and select OK
%     % on all plots, then type dbcont in command window
% %     EEG = pop_processMARA ( ALLEEG,EEG,CURRENTSET );
%     [ALLEEG,EEG] = processMARA ( ALLEEG,EEG,CURRENTSET );
% 
%     EEG = eeg_checkset( EEG );
%     
%     pause(2)
%     
%     MARAmarked = EEG.reject.gcompreject;
%     MARAmarked = find(MARAmarked); 
%     numMARAmarked = length(MARAmarked);
%     
%     % Set threshold to 30%
%     EEG.reject.gcompreject = zeros(size(EEG.reject.gcompreject)); 
%     EEG.reject.gcompreject(EEG.reject.MARAinfo.posterior_artefactprob > 0.30) = 1;
%     
%     AUTOdelete = EEG.reject.gcompreject;
%     AUTOdelete = find(AUTOdelete);
% 
%     % Remove selected components
%     EEG = pop_subcomp(EEG,[],0);
%     eeglab redraw
%     
%     % Save
%     EEG = eeg_checkset( EEG );
% %     [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
%     EEG.setname = [basefilename '_S04_MARA.set'];
%     pop_saveset( EEG, 'filename',EEG.setname,'filepath', datafolder);
%     clear('EEG')
% end
    
%     
%     [ALLEEG, EEG] = eeglab;
%     FileName2 = [basefilename '_S04_ICA_MAN.set'];
%     EEG = pop_loadset('filename', FileName2, 'filepath', datafolder);
%     eeglab redraw
%     
%     MANmarked = [EEG.icaCompClass.TESA1.reject EEG.icaCompClass.TESA1.rejectBlink EEG.icaCompClass.TESA1.rejectEyeMove EEG.icaCompClass.TESA1.rejectMuscle EEG.icaCompClass.TESA1.rejectElectrodeNoise EEG.icaCompClass.TESA1.rejectSensory ];
%     MANmarked = sort(MANmarked);
%     
%     AUTOsigComps = find(AUTOdelete < numSigComps+1);
%     AUTOsigComps = AUTOdelete(AUTOsigComps);
%     numAUTOsigComps = length(AUTOsigComps);
%     
%     for j = 1:length(AUTOsigComps)
%         x(j) = ismember(AUTOsigComps(j),MANmarked);
%     end
%     
%     match = find(x);
%     for m = 1:length(match)
%         matchComp(m) = AUTOsigComps(match(m));
%     end
%     
%     numCorrect = length(matchComp);
%     numFalsePos = length(x)-numCorrect;
%     FalsePosComp = find(x == 0);
%     FalsePosComp = AUTOsigComps(FalsePosComp);
%     strFalsePosComp = num2str(FalsePosComp);
%     
%     for n = 1:length(MANmarked)
%         y(n) = ismember(MANmarked(n),AUTOdelete);
%     end
%     
%     FalseNeg = find(y==0);
%     numFalseNeg = length(FalseNeg);
%     FalseNeg = MANmarked(FalseNeg);
%     strFalseNegComp = num2str(FalseNeg);
%     
%     % Save everything
%     fid = fopen(file_name,'a+');
%     fprintf(fid,'\n %s','      ');
%     fprintf(fid,'\n %s',['File: ' basefilename]);
%     fprintf(fid,'\n %-100s %.0f','Total number of ICs:',allcomps);
%     fprintf(fid,'\n %-100s %.0f','Total number of significant ICs (contribute to 99% of the variance):',numSigComps);
%     fprintf(fid,'\n %-100s %.0f','Total number of ICs removed by auto:',numMARAmarked);
%     fprintf(fid,'\n %-100s %.0f','Number of significant auto-removed ICs:',numAUTOsigComps);
%     fprintf(fid,'\n %-100s %.0f','Number of significant ICs labeled as artifact by both manual and auto (true positive):',numCorrect);
%     fprintf(fid,'\n %-100s %.0f','Number of significant ICs labeled as neural by manual but artifact by auto (false positives):',numFalsePos);
%     fprintf(fid,'\n %-100s %s','False positive ICs in auto:',strFalsePosComp);
%     fprintf(fid,'\n %-100s %.0f','Number of significant ICs labeled as artifact by manual but neural by auto (false negatives):',numFalseNeg);
%     fprintf(fid,'\n %-100s %s','False negative ICs in auto:',strFalseNegComp);
%     fprintf(fid,'\n %-100s %.0f','Number of significant ICs labeled as neural by both manual and auto (true negative):',(numSigComps-numCorrect-numFalsePos-numFalseNeg));
%     
%     clear FalsNeg
%     
%     clear ALLEEG EEG;
% 
% end
% fclose('all');

%%
% % =======================================================================
% % Step 7: Component selection using TESA 
% % MODIFIED JR and TM 05/2019
% % ======================================================================== 
clear all; close all; clc; 
homefolder = '/Users/timothymorris/MATLAB-Drive';
datafolder = '/Users/timothymorris/Desktop/BBHI_TMSEEG/EEGConcuss';
savefolder = '/Users/timothymorris/Desktop/BBHI_TMSEEG/EEGConcuss/results/';
cd(datafolder);

Files = dir('*_S04_MARA.set');
% Files = dir('*S06_*'); 
for i = 1 : length(Files)
    FileName = Files(i).name;
    
    clear EEG ALLEEG CURRENTSET ALLCOM;
    [ALLEEG, EEG, CURRENTSET ALLCOM] = eeglab;
    EEG = pop_loadset('filename', FileName, 'filepath', datafolder);
    Ind1 = strfind(FileName, '_S04_');
    basefilename = FileName(1:Ind1-1);
    eeglab redraw
    
    %First sort components by percent of variance explained
    [EEG, EEG.varsPerc] = tesa_sortcomps(EEG);
    
%     varRem = 99;
%     for k  = 1:length(EEG.varsPerc)
%         if varRem-EEG.varsPerc(k) > 0
%             sigComps(k) = k;
%         end
%         varRem = varRem-EEG.varsPerc(k);
%     end
%     numSigComps = length(sigComps);

    %Open Component activations in eegplot window
	  pop_eegplot(EEG, 0, 1, 0, [], 'winlength', 5, 'dispchans', 5);
    
    %Create filename for after component marking
    thisfile = [basefilename '_S04b_MARAmarked.set'];
    EEG.setname= thisfile; EEG.filepath = datafolder;

    %Now call TESA for automatic component selection
    %Does TMS pulse, blink, lateral eye movements, muscle and electrode
    EEG = tesa_compselect( EEG,'compCheck','on','figSize','medium','plotTimeX',[0 2998],'plotFreqX',[1 100], ...
        'tmsMuscle','off', ...
        'blink','on','blinkThresh',2.5, 'blinkElecs',{'Fp1','Fp2'}, 'blinkFeedback','off',...
        'move','on','moveThresh',2, 'moveElecs',{'F7','F8'},'moveFeedback','off', ...
        'muscle','on','muscleThresh',0.6,'muscleFreqWin',[30 100],'muscleFeedback','off',...
        'elecNoise','on','elecNoiseThresh',4,'elecNoiseFeedback','off' );
    %'muscleFreqIn',[7,75]
    %'muscleFreqEx',[48,52],
    EEG.filepath = []; %Removes filepath so it is not stored for later   
    thisfile = [basefilename '_S04_ICA_MARA.set'];
    EEG.setname= thisfile;
    EEG = pop_saveset( EEG, 'filename', thisfile, 'filepath',datafolder);
end
% 

%%
% Calculates sum of frequency band absolute and relative power within each frequency range, and
% peak alpha frequency and power from posterior electrodes.
% Output includes 2 plots and summary appended to frequency_band_power.txt.
% 
% window size = sampling rate
% overlap = 50%
%   TRT: FS = 5000 Hz with epochs 3 seconds long

clear all; close all; clc; 

homefolder = '/Users/timothymorris/MATLAB-Drive';
datafolder = '/Users/timothymorris/Desktop/BBHI_TMSEEG/EEGConcuss';
data_loc =  datafolder;

% addpath(genpath(homefolder));

% Get files
cd(data_loc);
Files = dir('*_S04_ICA_MARA.set');

% % Remove hidden files
% Files = Files(arrayfun(@(x) ~strcmp(x.name(1),'.'),Files));

% Make a file for spectral analysis
mkdir('spect_analysis');
addpath(genpath('spect_analysis/'));

% Open summary file for output
summ_path = strcat('spect_analysis/','frequency_band_power');
file_name = [summ_path sprintf('%d') '.txt'];
fopen(file_name,'a+');

for i = 1:length(Files)
    [ALLEEG, EEG, ~, ALLCOM] = eeglab;
    FileName = Files(i).name; 
    EEG = pop_loadset('filename', FileName, 'filepath', data_loc);
    Ind1 = strfind(FileName, '_S04_ICA_MARA.set');
    basefilename = FileName(1:Ind1-1);
    eeglab redraw
    
    % Frequency domain parameters
    nwinsize = EEG.srate;
    noverlap = (ceil(nwinsize * 0.50)); % 50% overlap in calculating power
    
    figure(2)
    [spectra2,freqs2,speccomp2,contrib2,specstd2]=spectopo(EEG.data, size(EEG.data,2), ...
        EEG.srate, 'winsize', nwinsize, 'overlap', noverlap, 'chanlocs', EEG.chanlocs, 'freqrange', [1 50]);
    hold on
%     axis([0 30 -20 30]);
    title(basefilename,'interpreter','none');
    hold off
    
    file_fig1_m = strcat('spect_analysis/',basefilename,'_spect1');
    file_fig1_jpg = strcat('spect_analysis/',basefilename,'_spect1.jpg');
%     saveas(2,file_fig1_m);
    saveas(2,file_fig1_jpg);
    
    % Power and frequencies
    EEG.powerspectra=db2pow(spectra2);
    EEG.Powerfreqs=freqs2;
    
    Results.powerspectra = db2pow(spectra2);
    Results.powerfreqs = freqs2;
    avgPSD = mean(EEG.powerspectra(:));
    sumPSD = sum(EEG.powerspectra(:));
    
    figure(3)
    plot(Results.powerfreqs,Results.powerspectra)
    hold on
    xlim([1 50]);
    xlabel('Frequency');
    ylabel('Power');
    title(basefilename,'interpreter','none');
    hold off
    
    file_fig2_m = strcat('spect_analysis/',basefilename,'_spect2.fig');
    file_fig2_jpg = strcat('spect_analysis/',basefilename,'_spect2.jpg');
%     saveas(3,file_fig2_m);
    saveas(3,file_fig2_jpg);

    % ABSOLUTE, all electrodes: Find delta=1-4, theta=4-8, alpha=8-13,
    % beta=13-30, gamma=30-48
    deltaPowerA = sum(Results.powerspectra(1+1:4+1));
    thetaPowerA = sum(Results.powerspectra(4+1:8+1));
    alphaPowerA = sum(Results.powerspectra(8+1:13+1));
    betaPowerA = sum(Results.powerspectra(13+1:30+1));
    gammaPowerA = sum(Results.powerspectra(30+1:48+1));

    % RELATIVE: Find delta=1-4, theta=4-8, alpha=8-13, beta=13-30,
    % gamma=30-48
    deltaPowerR = sum(Results.powerspectra(1+1:4+1))/sumPSD;
    thetaPowerR = sum(Results.powerspectra(4+1:8+1))/sumPSD;
    alphaPowerR = sum(Results.powerspectra(8+1:13+1))/sumPSD;
    betaPowerR = sum(Results.powerspectra(13+1:30+1))/sumPSD;
    gammaPowerR = sum(Results.powerspectra(30+1:48+1))/sumPSD;
    
    % Find the frequency of the alpha peak (8-13 Hz)
%     alphaPowerAll = Results.powerspectra(8+1:13+1);
%     alphaFreqAll = Results.powerfreqs(8+1:13+1);
%     [M,I] = max(alphaPowerAll);
%     alphaPeakPower = M;
%     alphaPeakFreq = alphaFreqAll(I);
    
    % Find the frequency of the alpha peak from posterior electrodes only
    % CHECK TO MAKE SURE THESE ELECTRODE NUMBERS ARE CORRECT! 
    % (Pz=14, P3=15, O1=17, Oz=18, O2=19, P4=20, P1=44, PO3=47, POz=48,
    % PO4=49, P2=52)
    alphaPowerAll = mean(Results.powerspectra(8+1:13+1,[14,15,17,18,19,20,44,47,48,49,52]));
    alphaFreqAll = Results.powerfreqs(8+1:13+1);
    [M,I] = max(alphaPowerAll);
    alphaPeakPower = M;
    alphaPeakFreq = alphaFreqAll(I);
    
    % Save everything
    fid = fopen(file_name,'a+');
    fprintf(fid,'\n \n %s %s %s %s %f',basefilename,',','delta range (1-4 Hz) absolute power',',',deltaPowerA);
    fprintf(fid,'\n %s %s %s %s %f',basefilename,',','theta range (4-8 Hz) absolute power',',',thetaPowerA);
    fprintf(fid,'\n %s %s %s %s %f',basefilename,',','alpha range (8-13) absolute power',',',alphaPowerA);
    fprintf(fid,'\n %s %s %s %s %f',basefilename,',','beta range (13-30) absolute power',',',betaPowerA);
    fprintf(fid,'\n %s %s %s %s %f',basefilename,',','gamma range (30-48 Hz) absolute power',',',gammaPowerA);
    
    fprintf(fid,'\n %s %s %s %s %f',basefilename,',','delta range (1-4 Hz) relative power',',',deltaPowerR);
    fprintf(fid,'\n %s %s %s %s %f',basefilename,',','theta range (4-8 Hz) relative power',',',thetaPowerR);
    fprintf(fid,'\n %s %s %s %s %f',basefilename,',','alpha range (8-13) relative power',',',alphaPowerR);
    fprintf(fid,'\n %s %s %s %s %f',basefilename,',','beta range (13-30) relative power',',',betaPowerR);
    fprintf(fid,'\n %s %s %s %s %f',basefilename,',','gamma range (30-48 Hz) relative power',',',gammaPowerR);
    
    fprintf(fid,'\n %s %s %s %s %f',basefilename,',','alpha peak frequency (from average of posterior electrodes)', ',',alphaPeakFreq);
    fprintf(fid,'\n %s %s %s %s %f',basefilename,',','alpha peak absolute power (from average of posterior electrodes)',',', alphaPeakPower);
    close(2);
    close(3);
    
    clear ALLEEG EEG;
end
  


%%
% %%
% % =======================================================================
% %  Step 8: Low-Pass Filter at 50Hz, then interpolate missing channels
% % =======================================================================
% 
% clear all; close all; clc; 
% homefolder = '/Users/timothymorris/MATLAB-Drive';
% datafolder = '/Users/timothymorris/Desktop/BBHI_TMSEEG/EEGConcuss';
% savefolder = '/Users/timothymorris/Desktop/BBHI_TMSEEG/EEGConcuss/results/';
% UpperEdge = 50; % VAR - LOWPASS EDGE
% 
% cd(datafolder);
% 
% EEG = pop_loadset('filename', 'BV_ElectrodeTemplate_New.set', 'filepath', datafolder);
% eeglab redraw
% [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
% Files = dir('*S04_ICA_MARA.set');
% for ii = 1 : length(Files)
%     FileName = Files(ii).name;
%     Ind1 = strfind(FileName, 'S04_ICA_MARA');
%     basefilename = FileName(1:Ind1-2);
%     EEG = pop_loadset('filename', FileName, 'filepath', datafolder);
%     [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
%     
%       %Low-pass filter the data
%     Fs=EEG.srate;ord=4;
%     [xall,yall]=butter(ord,UpperEdge/(Fs/2),'low');
%     
%     for trial = 1:size(EEG.data,3)
%         for ch=1:size(EEG.data,1)
%             EEG.data(ch,:,trial) = double(filtfilt(xall,yall, double(EEG.data(ch,:,trial))));
%         end
%     end
%     
%     %And interpolate missing channels
%     EEG = pop_interp(EEG, ALLEEG(1).chanlocs, 'spherical');
% %     
% %     %And delete EOG leads
% %     EEG = pop_select( EEG,'nochannel',{'EOG1' 'EOG2'});
%     
%     EEG.setname = [basefilename '_S05_Final_MARA'];
%     EEG = pop_saveset( EEG, 'filename', EEG.setname, 'filepath', datafolder);
%     ALLEEG = pop_delset( ALLEEG, [2] );
%     [EEG ALLEEG CURRENTSET] = eeg_retrieve(ALLEEG,1);
% end
% 
% %% 
% % =======================================================================
% %  Step 9: GMFA Analysis: Calculate & Plot GMFA and peaks
% % =======================================================================
% % 
% clear all; close all; clc; 
% 
% homefolder = '/Users/timothymorris/MATLAB-Drive';
% datafolder = '/Users/timothymorris/Desktop/BBHI_TMSEEG/EEGConcuss';
% savefolder = '/Users/timothymorris/Desktop/BBHI_TMSEEG/EEGConcuss/results/';
% 
% cd(datafolder);
% Linetime = 0; %Time of pulse
% Graph1_start = -200; %VAR Time at which to begin topoplot
% Graph1_end = 400; %VAR Time at which to end topoplot
% TOI_start = 20; %VAR Time of interest start for peakfinding
% TOI_end = 400; %VAR Time of interest end for peakfinding
% threshval = 2; %VAR Number of standard deviations above the mean to be a 
%                 % minimum peak after stimulation; note that this does not
%                 % need to be defined at all, but probably should be, at
%                 % least as the max baseline voltage
% selval = 3; %VAR Criteria for labeling something a peak; does not need to 
%             % be set manually
% first_basetime = -450; %VAR Beginning of official baseline period
% last_basetime = -50; %VAR End of official baseline period
% FStitle = 16; %VAR Figure title font size
% FS = 14; %VAR Figure axes font size
% 
% Files = dir('*_S05_*.set'); 
% for i = 1 : length(Files)
%     FileName = Files(i).name;
%     clear EEG ALLEEG CURRENTSET ALLCOM;
%     [ALLEEG, EEG, CURRENTSET ALLCOM] = eeglab;
%     EEG = pop_loadset('filename', FileName, 'filepath', datafolder);
%     Ind1 = strfind(FileName, '.set');
%     basefilename = FileName(1:Ind1-1);
%     eeglab redraw
%     
%     %Define variables used in script
%     graphstart_bin = find((EEG.times>=Graph1_start), 1, 'first');
%     graphend_bin = find((EEG.times>=Graph1_end),1, 'first');
%     firstbasebin = find((EEG.times>=first_basetime),1,'first');
%     lastbasebin = find((EEG.times<=last_basetime),1,'last');
%     pulse1_latency = EEG.event(1,1).latency;
%     TOIstart_bin = find((EEG.times>=TOI_start), 1, 'first');
%     TOIend_bin = find((EEG.times>=TOI_end),1, 'first');
% %     First_peak_time = ceil(EEG.tmscut.cutTimesTMS(2));
%     
%     % GMFA calculations
%     %----------------------------------------------------------------------
%     %Do Standard GMFA analysis here using TESA
%     EEG = pop_tesa_tepextract( EEG, 'GMFA'); %Standard GMFA analysis
%     
%     % Calculate the normalized GMFA, which is GMFA divided by the mean GMFA
%     % in the baseline period
%     AvgbaseGMFA=mean(EEG.GMFA.R1.tseries(firstbasebin:lastbasebin));
%     EEG.GMFA.norm.tseries = EEG.GMFA.R1.tseries/AvgbaseGMFA;
%     EEG.GMFA.norm.time = EEG.GMFA.R1.time;
%     
%     %Now detect peaks automatically (note that this is NOT a TESA function,
%     %but rather uses the peakfinder script, and is my own analysis)
%     %----------------------------------------------------------------------
%     %First define the threshold - threshval SD above mean baseline value
%     tempdata = [];
%     tempdata = EEG.GMFA.R1.tseries;
%     peakLoc = []; peakMag = []; peaks=[]; 
%     mean_baseval = mean(tempdata(firstbasebin:lastbasebin));
%     std_baseval = std(tempdata(firstbasebin:lastbasebin));
%     maxbaseval = max(tempdata(firstbasebin:lastbasebin));
%     thresh = max(((std_baseval * threshval) + mean_baseval),maxbaseval); 
%    
%     %And selection criteria for something to be defined a peak
%     sel = (max(tempdata(firstbasebin:lastbasebin))-min(tempdata(firstbasebin:lastbasebin)))/selval;
%     
%     %Now find and save the peaks
%     [peakLoc, peakMag] = peakfinder(tempdata, sel, thresh);
%     peakLoc=EEG.GMFA.R1.time(peakLoc);
%     peaks(:,1)=peakLoc';
%     peaks(:,2)=peakMag';
%     
%     %Delete peaks before time region of interest
%     numpeaks = size(peaks,1);
%     for count=numpeaks:-1:1
%         if ((peaks(count,1)<TOI_start))
%             peaks(count,:) = [];
%         end
%     end
%     
%     %And Delete peaks after time region of interest
%     numpeaks = size(peaks,1);
%     for count=numpeaks:-1:1
%         if ((peaks(count,1))>TOI_end)
%             peaks(count,:) = [];
%         end
%     end
%     
%     EEG.GMFA.R1.peaks(:,1)= peaks(:,1);
%     EEG.GMFA.R1.peaks(:,2)= peaks(:,2);
%     
%     %Plot nGMFP
% %     figure;plot(EEG.times(graphstart_bin:graphend_bin), EEG.GMFA.norm.tseries(1,graphstart_bin:graphend_bin)); hold on;
% %     titlename = [basefilename ' Normalized GMFP'];
% %     title(titlename);
% %     ymax = max(EEG.GMFA.norm.tseries);
% %     axis([Graph1_start Graph1_end 0 (ymax + (ymax/10))]);
% %     line([Linetime, Linetime], [0 ymax], 'Color', 'r'); hold off; 
% %     hgexport(gcf, titlename, hgexport('factorystyle'), 'Format', 'png'); 
% % %     saveas(gcf,titlename);
% %     saveas(gcf,titlename, 'png');
%     
% %     %Save the resulting figure
% %     exportfig(gcf, titlename,'format', 'png', 'Color', 'cmyk', ...
% %         'Resolution', 600, 'FontMode', 'scaled', 'Bounds', 'tight', ...
% %         'LineMode','scaled');
% %     clf;
%     
%     %Plot ERPtopo
% %     topotitle = [basefilename ' ERPTopo'];
% %     figure; pop_timtopo(EEG, [Graph1_start Graph1_end], peaks(:,1));
% %     gtext(topotitle, 'fontsize', 12);
% %     hgexport(gcf, topotitle, hgexport('factorystyle'), 'Format', 'png');
%     
% %     %Save the resulting figure
% %     exportfig(gcf, topotitle,'format', 'png', 'Color', 'cmyk', ...
% %         'Resolution', 600, 'FontMode', 'scaled', 'FontSize', 10, ...
% %         'FontSizeMin', 8, 'FontSizeMax', 12, 'Bounds', 'tight', ...
% %         'LineMode','scaled', 'Width', 7);
% %     clf;
% 
%     
%     %And save the resulting file
%     EEG.setname = [basefilename '_S06_GMFA.set'];
%     EEG = pop_saveset( EEG, 'filename', EEG.setname, 'filepath', datafolder);
%     close all;
%     eeglab redraw;   
% end
