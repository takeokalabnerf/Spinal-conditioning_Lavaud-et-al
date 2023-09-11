% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
% Creation of table with electrophysiology information related to:
% - Controls : frY
% - Learners : frM
%
% Last updated: 11/09/23
% Simon Lavaud, Modified by Charlotte Bichara
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


% - - - Help on how to use it if you already have table:
%
% If you have frM frY table from preivous recording
% The best is to load your previous table, rename it as frM_old frY_old.
% Then you can run for a specific recording the frM frY algorithm.
% It will return a frM frY only for one recording.
% You can copy and past this new frM frY in the old version.
% Then delete the new frM and frY, and rename the frM_old and frY_old
% as frM and frY. You can now save it in the right folder.

% Load the database of interest with the directory of the different
% recording we want to use for creating the informative table
recording_database_updated_20220718;

% Initialize the variable
nch = 385; % Number of channels
Fs = 30000; % Freqency of the acquisition
analysis_type = 1; 

% Column of the table for analysis
Recording = [];
Id_spont = [];
Fr_spont = [];
Id_horr = [];
Fr_early = [];
Fr_middle = [];
Fr_late = [];
class = [];


%%
% % % % % % % % % % % % % % % % % % % % % % % %
% - - -
% - - - - - frM
% - - -
% % % % % % % % % % % % % % % % % % % % % % % %

clearvars frM

% Initialize inner variable
cont = 0;

% Iteration on the recording
% If only one recording to add : i = 1
% If all the recording: [7 8 10 11 12 13 14 15 16 17 18 19 20 21 22]

% i = number of the recording
for i = 8
    QQ=i;
    
    % Display the recording that is analysed
    disp(sprintf('Learners recording number %d',i))
    
    % Load the recording in a specified directory
    spk = loadKSdir(sort_masters_horridge{i});
    % From the recording, extract only the single units
    clusters = spk.cids(spk.cgs == 2);
    
    % Extract the event for each recording
    load(fullfile(sort_masters_horridge{i},'events.mat'),'events');
    evnames = {events.type};
    evtid = find(contains(evnames,'horridge')); % Recording with only conditionning experiment
    if isempty(evtid)==1
        evtid = find(contains(evnames,'opto')); % Recording with conditionning and optotagging
    end
    
    % Event description are saved in an event file created with the
    % recording

    if i==8 
        %End of master horridge number 8 at 32min
        lim_start=1; % Beginning of the conditioning task start with the first event
        lim_event=156; % Last event of the conditioning task    
    else % Initialization parameter in case of problem
        lim_start=2;
        lim_event=length(events(evtid).onsets);
    end
    
    % Computation of the timeline Early, mid, late
    % Pre criterion is: Early + Mid
    % Post criterion is: late
    % Mid not use anymore
    
    clearvars dstim tdist
    a=1;
    
    %Frequency
    for j = (lim_start+1):lim_event
        tdist(j) = events(evtid).onsets(j) - events(evtid).onsets(j-1);
    end
    
    dstim = smoothdata(tdist);
    idmax = find(dstim>1,1);
    if isempty(idmax)
        idmax=lim_start+10;
    end
    filter_non_onset=1;
    
    % Tbin
    % 0 to 600 - spont
    % 600 to time early (early phase)
    % time early to time mid (mid phase)
    % time mid to end horridge (late phase)
    % end horridge - end+600 (after phase)
    
    if ~isempty(idmax)
        tbin = [events(evtid).onsets(lim_start)-600, events(evtid).onsets(lim_start)-1;... %Spont
            events(evtid).onsets(lim_start), events(evtid).onsets(idmax);... %Early
            events(evtid).onsets(idmax), events(evtid).onsets(lim_event);... %Mid
            events(evtid).onsets(lim_event), events(evtid).onsets(lim_start)+600]; %Late
    else
        tbin = [events(evtid).onsets(1)-600, events(evtid).onsets(1)-1;...
            events(evtid).onsets(1), 60;...
            60, events(evtid).onsets(end);...
            events(evtid).onsets(end), events(evtid).onsets(1)+600];
    end
    
    % Recording number
    % Id of the single unit cluster for a specific recording 
    L = 1;
    step = 0.1;
    for j = 1:length(clusters)
        cont = cont+1;
        disp(sprintf('Learner %d, cluster %d',i,j))
        Recording(cont) = i;
        Id_horr(cont) = clusters(j);
        
        % Create the FR column for each phases depending on Tbin
        for k = 1:size(tbin,1)
            win = [tbin(k,1) tbin(k,1)+L];
            frmov = [];
            edges = tbin(k,1):L:tbin(k,2);
            frmov = histcounts(spk.st(spk.clu == clusters(j)),edges)/L;
            switch k
                case 1
                    Fr_spont{cont} = frmov;
                case 2
                    Fr_early{cont} = frmov;
                case 3
                    Fr_middle{cont} = frmov;
                case 4
                    Fr_late{cont} = frmov;
            end
        end
        
        k=4;
        win = [1500 1500+L];
        frmov = [];
        while win(2)<=2100
            stj = spk.st(spk.clu == clusters(j) & spk.st>win(1) & spk.st<=win(2));
            frmov = [frmov numel(stj)/L];
            win = win+L;
        end
        Fr_after{cont} = frmov;
        
        Fr_average_after(cont)=mean(frmov);
        
        Fr_average_spont(cont)=mean(Fr_spont{cont});
        Tbin_tot{cont}=tbin;
    end

end


% Create the frM table
frM = table();
frM.Recording = Recording';
frM.Tbin = Tbin_tot';
frM.Fr_spont = Fr_spont';
frM.Fr_average_spont = Fr_average_spont';
frM.Id_horr = Id_horr';
frM.Fr_early = Fr_early';
frM.Fr_middle = Fr_middle';
frM.Fr_late = Fr_late';
frM.Fr_after = Fr_after';
frM.Fr_average_after = Fr_average_after';




%%

% % % % % % % % % % % % % % % % % % % % % % % %
% - - -
% - - - - - frY
% - - -
% % % % % % % % % % % % % % % % % % % % % % % %

% Work exactly as the Leaner code

RecordingY = [];
Fr_spontY = [];
Id_horrY = [];
Fr_earlyY = [];
Fr_middleY = [];
Fr_lateY = [];
classY = [];

cont = 0;
for i = 7 
    QQ=i;
    
    disp(sprintf('Controls %d',i))
    
    spk = loadKSdir(sort_yokes_horridge{i});
    clusters = spk.cids(spk.cgs == 2);
    
    %crete time binning
    load(fullfile(sort_yokes_horridge{i},'events.mat'),'events');
    evnames = {events.type};
    evtid = find(contains(evnames,'horridge'));
    if isempty(evtid)==1
        evtid = find(contains(evnames,'opto'));
    end
    
    if i==7   
        lim_start=1; 
        lim_event=1603; 
    end
    
    clearvars dstim tdist a
    a=1;
    for j = (lim_start+1):lim_event
        tdist(a) = events(evtid).onsets(j) - events(evtid).onsets(j-1);
        a=a+1;
    end
    dstim = smoothdata(tdist);
    
    idmax = find(dstim>1,1);
    if isempty(idmax)
        idmax=lim_start+10;
    end
    
    
    if ~isempty(idmax)
        tbin = [events(evtid).onsets(lim_start)-600, events(evtid).onsets(lim_start)-1;... %Spont
            events(evtid).onsets(lim_start), events(evtid).onsets(idmax);... %Early
            events(evtid).onsets(idmax), events(evtid).onsets(lim_event);... %Mid
            events(evtid).onsets(lim_event), events(evtid).onsets(lim_start)+600]; %Late
    else
        tbin = [events(evtid).onsets(1)-600, events(evtid).onsets(1)-1;...
            events(evtid).onsets(1), 60;...
            60, events(evtid).onsets(end);...
            events(evtid).onsets(end), events(evtid).onsets(1)+600];
    end
    
    L = 1;
    step = 0.1;
    for j = 1:length(clusters)
        disp(sprintf('Controls %d, cluster %d',i,j))
        cont = cont+1;
        RecordingY(cont) = i;
        Id_horrY(cont) = clusters(j);
        Id_spontY(cont) = clusters(j);
        Id_afterY(cont) = clusters(j);
        
        for k = 1:size(tbin,1)
            win = [tbin(k,1) tbin(k,1)+L];
            frmov = [];
            edges = tbin(k,1):L:tbin(k,2);
            frmov = histcounts(spk.st(spk.clu == clusters(j)),edges)/L;
            switch k
                case 1
                    Fr_spontY{cont} = frmov;
                case 2
                    Fr_earlyY{cont} = frmov;
                case 3
                    Fr_middleY{cont} = frmov;
                case 4
                    Fr_lateY{cont} = frmov;
            end
        end
        k=4;
        win = [1500 1500+L];
        frmov = [];
        while win(2)<=2100
            stj = spk.st(spk.clu == clusters(j) & spk.st>win(1) & spk.st<=win(2));
            frmov = [frmov numel(stj)/L];
            win = win+L;
        end
        Fr_afterY{cont} = frmov;
        
        Fr_average_afterY(cont)=mean(frmov);
        Fr_average_spontY(cont)=mean(Fr_spontY{cont});
        Tbin_totY{cont}=tbin;
        
    end
end


frY = table();
frY.Recording = RecordingY';
frY.Tbin = Tbin_totY';
frY.Fr_average_spont = Fr_average_spontY';
frY.Fr_spont = Fr_spontY';
frY.Id_horr = Id_horrY';
frY.Fr_early = Fr_earlyY';
frY.Fr_middle = Fr_middleY';
frY.Fr_late = Fr_lateY';
frY.Fr_after = Fr_afterY';
frY.Fr_average_after = Fr_average_afterY';

%%
% % % % % % % % % % % % % % % % % % % % % % % %
% - - -
% - - - - - Up-regulation / Down-regulation frM frY
% - - -
% - - -
% - - -
% - - - - - 0=non participatiig; 1:increased; 2=decreased; 3=disappeared
% - - -
% - - -
% % % % % % % % % % % % % % % % % % % % % % % %

% - - - - - - - - - -
% - - - - - Initialize all the variable
% - - - - - - - - - -
clearvars compa_after zdata_early zdata_mid zdata_late zdata
clearvars evolution class2 evolution_cdf_pre_post evolution_fi evolution_cdf_pre_post_fi
clearvars compa_after zdata_early zdata_early_up zdata_early_down
clearvars zdata_mid zdata_mid_up zdata_mid_down zdata_late zdata
clearvars zdata_late zdata_late_up zdata_late_down
clearvars zdata zdata2 zdata3 zgroup_early zgroup_mid zgroup_late
clearvars mod_early_up_zscore mod_early_down_zscore
clearvars mod_mid_up_zscore mod_mid_down_zscore
clearvars mod_late_up_zscore mod_late_down_zscore
clearvars zdata_EandMmerge_up zdata_EandMmerge_down zdata_EandMmerge
clearvars mod_EandMmerge_up_zscore mod_EandMmerge_down_zscore zgroup_EandMmerge
clearvars EandMmerge


% - - - Parameter for up/down regulation
% Threshold of z-score detection
threshold_zscore=2;
% Threshold for p-value comparison 
% (not used anymore)
pvalue_threshold=0.001;
% Smoothing factor of the Fr (0.001 is almost no smoothing)
smoothing_factor=0.001;
% Threshold for persistent/intermitent modulation (0 = 0%, 1=100% of Fr modulated)
% (not used anymore)
threshold_per_zscore=0.5;
threshold_weak=0.2;

%%
% - - - - - - - -
% - - - Start calculating Firing rate
% - - - - - - -
clc
clearvars zgroup_EandMmerge zdata_EandMmerge zdata_EandMmerge_up zdata_EandMmerge_down
clearvars zgroup_late zdata_late zdata_late_up zdata_late_down

% Iterating on the size of previously calculated frM
for i = 1:size(frM,1)
    
    %Initialization
    evolution=0;
    evolution_cdf_pre_post=0;
    
    % Compute all the Fr of all the phase for futur calculation
    FrM_spont_s{i}=smoothdata(frM.Fr_spont{i},'gaussian', smoothing_factor);
    FrM_early_s{i}=smoothdata(frM.Fr_early{i},'gaussian', smoothing_factor);
    FrM_middle_s{i}=smoothdata(frM.Fr_middle{i},'gaussian', smoothing_factor);
    FrM_late_s{i}=smoothdata(frM.Fr_late{i},'gaussian', smoothing_factor);
    FrM_after_s{i}=smoothdata(frM.Fr_after{i},'gaussian', smoothing_factor);
    
    % Mean of the firing rate spont
    base = mean(frM.Fr_spont{i}(300:599));
    % Std of firing rate spont
    basestd = std(frM.Fr_spont{i}(300:599));
    
    % If std spont = 0 then assign std=0.1
    % Used for computing zscore
    if basestd==0
        basestd=0.1;
    end
    
    % - - -
    % Z-score calculation for each phase : EARLY / MID / LATE
    % - - -
    
    % intialization of the variable
    b=0;
    c=0;
    d=0;
    clearvars zdata
    mod_EandMmerge_up=[];
    mod_EandMmerge_down=[];
    
    clearvars zdata
    zdata=[];
    for j=1:size(FrM_early_s{i},2)
        zdata(j)= (FrM_early_s{i}(j) - base) / basestd;
    end
    zdata_early_full=zdata;
    
	clearvars zdata
    zdata=[];
    for j=1:size(FrM_middle_s{i},2)
        zdata(j)= (FrM_middle_s{i}(j) - base) / basestd;
    end
    zdata_mid_full=zdata;
    
    % Merge EARLY and MID = PRE-CRITERION PHASE
    EandMmerge=[FrM_early_s{i} FrM_middle_s{i}];
    
    for j=1:size(EandMmerge,2)
        
        % Z-score calculation
        zdata(j)= (EandMmerge(j) - base) / basestd;
        
        % If crossing the z-score threshold (either + or -) then b value increase
        if zdata(j)>threshold_zscore | zdata(j)<-threshold_zscore
            b=b+1;
        end
        
        % If crossing the z-score threshold (+) then c value increase
        if zdata(j)>threshold_zscore
            c=c+1;
            % If extraction of the modulation of the unit during the phase
            % Here is the calculation (here juste Fr phase - Fr spont)
            mod_EandMmerge_up=[(EandMmerge(j))-(base) mod_EandMmerge_up];
        end
        
        % If crossing the z-score threshold (-) then c value increase
        if zdata(j)<-threshold_zscore
            d=d+1;
            mod_EandMmerge_down=[(EandMmerge(j))-(base) mod_EandMmerge_down];
        end
        
    end
    
    % Percentage of zscore significant for + and - the threshold
    zdata_EandMmerge(i)=b/size(EandMmerge,2);
    % Percentage of zscore significant for + the threshold
    zdata_EandMmerge_up(i)=c/size(EandMmerge,2);
    % Percentage of zscore significant for - the threshold
    zdata_EandMmerge_down(i)=d/size(EandMmerge,2);
    % Raw zscore
    zdata_EandMmerge_full=zdata;
    
    % - - - - - - -
    % - - - Classify up/down regulation
    % - - - - - - -
    
    % - - Normal increase /// CAT = 1
    % example : z-score up > percentage threshold persistent % AND z-score down < percentage threshold intermitent
    if zdata_EandMmerge_up(i)>threshold_per_zscore & zdata_EandMmerge_down(i)<threshold_weak % ==0 %& mean(frM.Fr_early{i})>base
        zgroup_EandMmerge(i)=1;
        
    % - - Weak increase /// CAT = 10
    elseif zdata_EandMmerge_up(i)>threshold_weak & zdata_EandMmerge_down(i)<threshold_weak %==0 %& mean(frM.Fr_early{i})>base
        zgroup_EandMmerge(i)=10;
        
    % - - Not participating /// CAT = 0
    elseif zdata_EandMmerge_up(i)<threshold_weak & zdata_EandMmerge_down(i)<threshold_weak
        zgroup_EandMmerge(i)=0;
        
    % - - Normal decrease /// CAT = 2
    elseif  zdata_EandMmerge_up(i)<threshold_weak & zdata_EandMmerge_down(i)>threshold_per_zscore %& mean(frM.Fr_early{i})<base
        zgroup_EandMmerge(i)=2;
        
    % - - Weak decrease /// CAT = 20
    elseif zdata_EandMmerge_up(i)<threshold_weak & zdata_EandMmerge_down(i)>threshold_weak
        zgroup_EandMmerge(i)=20;
        
    % - - Increase AND Decrease /// CAT = 3 or 4, depending on the order of increase/decrease
    elseif zdata_EandMmerge_up(i)>threshold_per_zscore & zdata_EandMmerge_down(i)>threshold_per_zscore
        if zdata_EandMmerge_up(i)>zdata_EandMmerge_down(i)
            zgroup_EandMmerge(i)=3;
        else
            zgroup_EandMmerge(i)=4;
        end
        
    % - - Not able to categorise /// CAT = 5, neurons that behave wierdly
    else
        zgroup_EandMmerge(i)=5;
    end
    
    
    % LATE
    % /!\ same methods as for the merge EARLY MID
    a=1;
    b=0;
    c=0;
    d=0;
    zdata=[];
    mod_late_up=[];
    mod_late_down=[];
    for j=1:size(FrM_late_s{i},2)
        zdata(j)= (FrM_late_s{i}(j) - base) / basestd;
        
        if zdata(j)>threshold_zscore | zdata(j)<-threshold_zscore
            b=b+1;
        end
        if zdata(j)>threshold_zscore
            c=c+1;
            mod_late_up=[(frM.Fr_late{i}(j))-(base) mod_late_up];
        end
        if zdata(j)<-threshold_zscore
            d=d+1;
            mod_late_down=[(frM.Fr_late{i}(j))-(base) mod_late_down];
        end
        
    end
    
    zdata_late(i)=b/size(FrM_late_s{i},2);
    zdata_late_up(i)=c/size(FrM_late_s{i},2);
    zdata_late_down(i)=d/size(FrM_late_s{i},2);
    zdata_late_full=zdata;
    
    % Normal increase
    if zdata_late_up(i)>threshold_per_zscore & zdata_late_down(i)<threshold_weak % ==0 %& mean(frM.Fr_early{i})>base
        zgroup_late(i)=1;
        
        % Weak increase
    elseif zdata_late_up(i)>threshold_weak & zdata_late_down(i)<threshold_weak %==0 %& mean(frM.Fr_early{i})>base
        zgroup_late(i)=10;
        
        % Not participating
    elseif zdata_late_up(i)<threshold_weak & zdata_late_down(i)<threshold_weak
        zgroup_late(i)=0;
        
        % Normal decrease
    elseif  zdata_late_up(i)<threshold_weak & zdata_late_down(i)>threshold_per_zscore %& mean(frM.Fr_early{i})<base
        zgroup_late(i)=2;
        
        % Weak decrease
    elseif zdata_late_up(i)<threshold_weak & zdata_late_down(i)>threshold_weak
        zgroup_late(i)=20;
        
        % Increase AND Decrease
    elseif zdata_late_up(i)>threshold_per_zscore & zdata_late_down(i)>threshold_per_zscore
        if zdata_late_up(i)>zdata_late_down(i)
            zgroup_late(i)=3;
        else
            zgroup_late(i)=4;
        end
        
        % Not able to categorise
    else
        zgroup_late(i)=5;
    end
    
    % Raw zscore data
  	frM.zscore_early_full{i}=zdata_early_full;
    frM.zscore_mid_full{i}=zdata_mid_full;
    frM.zscore_late_full{i}=zdata_late_full;

    
    clearvars zdata_EandMmerge_full zdata_late_full
end

frM.zscore_group_EandMmerge=zgroup_EandMmerge';
frM.zscore_EandMmerge_up=zdata_EandMmerge_up';
frM.zscore_EandMmerge_down=zdata_EandMmerge_down';

frM.zscore_group_late=zgroup_late';
frM.zscore_late=zdata_late';
frM.zscore_late_up=zdata_late_up';
frM.zscore_late_down=zdata_late_down';

%%
% Same methods as for frM


clearvars zgroup_EandMmerge zdata_EandMmerge zdata_EandMmerge_up zdata_EandMmerge_down
clearvars zgroup_late zdata_late zdata_late_up zdata_late_down
% - - - - - - - -
% - - - Start calculating FrY
% - - - - - - - -

% Iterating on the size of previously calculated frM
for i = 1:size(frY,1)
    
    % Compute all the Fr of all the phase for futur calculation
    FrY_spont_s{i}=smoothdata(frY.Fr_spont{i},'gaussian', smoothing_factor);
    FrY_early_s{i}=smoothdata(frY.Fr_early{i},'gaussian', smoothing_factor);
    FrY_middle_s{i}=smoothdata(frY.Fr_middle{i},'gaussian', smoothing_factor);
    FrY_late_s{i}=smoothdata(frY.Fr_late{i},'gaussian', smoothing_factor);
    FrY_after_s{i}=smoothdata(frY.Fr_after{i},'gaussian', smoothing_factor);
    
    % Mean of the FR spont
    base = mean(frY.Fr_spont{i}(300:599));
    % Std of Fr spont
    basestd = std(frY.Fr_spont{i}(300:599));
    
    % If std spont = 0 then assign std=0.1
    if basestd==0
        basestd=0.1;
    end
    
    % - - -
    % Z-score calculation for each phase : EARLY / LATE
    % - - -
    
    % Early and Mid MERGE = EARLY
    
    clearvars zdata
    zdata=[];
    for j=1:size(FrY_early_s{i},2)
        zdata(j)= (FrY_early_s{i}(j) - base) / basestd;
    end
    zdata_early_full=zdata;
    
	clearvars zdata
    zdata=[];
    for j=1:size(FrY_middle_s{i},2)
        zdata(j)= (FrY_middle_s{i}(j) - base) / basestd;
    end
    zdata_mid_full=zdata;
    
    
    % intialization of the variable
    b=0;
    c=0;
    d=0;
    clearvars zdata
    mod_EandMmerge_up=[];
    mod_EandMmerge_down=[];
    
    % Merge EARLY and MID
    EandMmerge=[FrY_early_s{i} FrY_middle_s{i}];
    
    
    for j=1:size(EandMmerge,2)
        
        % Z-score calculation
        zdata(j)= (EandMmerge(j) - base) / basestd;
        
        % If crossing the z-score threshold (either + or -) then b value increase
        if zdata(j)>threshold_zscore | zdata(j)<-threshold_zscore
            b=b+1;
        end
        
        % If crossing the z-score threshold (+) then c value increase
        if zdata(j)>threshold_zscore
            c=c+1;
            % If extraction of the modulation of the unit during the phase
            % Here is the calculation (here juste Fr phase - Fr spont)
            mod_EandMmerge_up=[(EandMmerge(j))-(base) mod_EandMmerge_up];
        end
        
        % If crossing the z-score threshold (-) then c value increase
        if zdata(j)<-threshold_zscore
            d=d+1;
            %mod_early_down=[(frY.Fr_early{i}(a)+1)/(base+1) mod_early_down];
            mod_EandMmerge_down=[(EandMmerge(j))-(base) mod_EandMmerge_down];
        end
        
    end
    
    % Percentage of zscore significant for + and - the threshold
    zdata_EandMmerge(i)=b/size(EandMmerge,2);
    % Percentage of zscore significant for + the threshold
    zdata_EandMmerge_up(i)=c/size(EandMmerge,2);
    % Percentage of zscore significant for - the threshold
    zdata_EandMmerge_down(i)=d/size(EandMmerge,2);
    % Raw zscore
    zdata_EandMmerge_full=zdata;
    
    % - - - - - - -
    % - - - Classify up/down regulation
    % - - - - - - -
    
    
    % Normal increase /// CAT = 1
    % example : z-score up > percentage threshold persistent % AND z-score down < percentage threshold intermitent
    if zdata_EandMmerge_up(i)>threshold_per_zscore & zdata_EandMmerge_down(i)<threshold_weak % ==0 %& mean(frM.Fr_early{i})>base
        zgroup_EandMmerge(i)=1;
        
        % Weak increase /// CAT = 10
    elseif zdata_EandMmerge_up(i)>threshold_weak & zdata_EandMmerge_down(i)<threshold_weak %==0 %& mean(frM.Fr_early{i})>base
        zgroup_EandMmerge(i)=10;
        
        % Not participating /// CAT = 0
    elseif zdata_EandMmerge_up(i)<threshold_weak & zdata_EandMmerge_down(i)<threshold_weak
        zgroup_EandMmerge(i)=0;
        
        % Normal decrease /// CAT = 2
    elseif  zdata_EandMmerge_up(i)<threshold_weak & zdata_EandMmerge_down(i)>threshold_per_zscore %& mean(frM.Fr_early{i})<base
        zgroup_EandMmerge(i)=2;
        
        % Weak decrease /// CAT = 20
    elseif zdata_EandMmerge_up(i)<threshold_weak & zdata_EandMmerge_down(i)>threshold_weak
        zgroup_EandMmerge(i)=20;
        
        % Increase AND Decrease /// CAT = 3 or 4, depending on the order of increase/decrease
    elseif zdata_EandMmerge_up(i)>threshold_per_zscore & zdata_EandMmerge_down(i)>threshold_per_zscore
        if zdata_EandMmerge_up(i)>zdata_EandMmerge_down(i)
            zgroup_EandMmerge(i)=3;
        else
            zgroup_EandMmerge(i)=4;
        end
        
        % Not able to categorise /// CAT = 5, neurons that behave wierdly
    else
        zgroup_EandMmerge(i)=5;
    end
    
    
    % LATE
    % /!\ same methods as for the merge EARLY MID
    a=1;
    b=0;
    c=0;
    d=0;
    zdata=[];
    mod_late_up=[];
    mod_late_down=[];
    for j=1:size(FrY_late_s{i},2)
        zdata(j)= (FrY_late_s{i}(j) - base) / basestd;
        
        if zdata(j)>threshold_zscore | zdata(j)<-threshold_zscore
            b=b+1;
        end
        if zdata(j)>threshold_zscore
            c=c+1;
            mod_late_up=[(frY.Fr_late{i}(j))-(base) mod_late_up];
        end
        if zdata(j)<-threshold_zscore
            d=d+1;
            mod_late_down=[(frY.Fr_late{i}(j))-(base) mod_late_down];
        end
        
    end
    
    zdata_late(i)=b/size(FrY_late_s{i},2);
    zdata_late_up(i)=c/size(FrY_late_s{i},2);
    zdata_late_down(i)=d/size(FrY_late_s{i},2);
    zdata_late_full=zdata;
    
    % Normal increase
    if zdata_late_up(i)>threshold_per_zscore & zdata_late_down(i)<threshold_weak % ==0 %& mean(frM.Fr_early{i})>base
        zgroup_late(i)=1;
        
        % Weak increase
    elseif zdata_late_up(i)>threshold_weak & zdata_late_down(i)<threshold_weak %==0 %& mean(frM.Fr_early{i})>base
        zgroup_late(i)=10;
        
        % Not participating
    elseif zdata_late_up(i)<threshold_weak & zdata_late_down(i)<threshold_weak
        zgroup_late(i)=0;
        
        % Normal decrease
    elseif  zdata_late_up(i)<threshold_weak & zdata_late_down(i)>threshold_per_zscore %& mean(frM.Fr_early{i})<base
        zgroup_late(i)=2;
        
        % Weak decrease
    elseif zdata_late_up(i)<threshold_weak & zdata_late_down(i)>threshold_weak
        zgroup_late(i)=20;
        
        % Increase AND Decrease
    elseif zdata_late_up(i)>threshold_per_zscore & zdata_late_down(i)>threshold_per_zscore
        if zdata_late_up(i)>zdata_late_down(i)
            zgroup_late(i)=3;
        else
            zgroup_late(i)=4;
        end
        
        % Not able to categorise
    else
        zgroup_late(i)=5;
    end
    
    % Raw zscore data
  	frY.zscore_early_full{i}=zdata_early_full;
    frY.zscore_mid_full{i}=zdata_mid_full;
    frY.zscore_late_full{i}=zdata_late_full;
    
    clearvars zdata_EandMmerge_full zdata_late_full
end

frY.zscore_group_EandMmerge=zgroup_EandMmerge';
frY.zscore_EandMmerge_up=zdata_EandMmerge_up';
frY.zscore_EandMmerge_down=zdata_EandMmerge_down';

frY.zscore_group_late=zgroup_late';
frY.zscore_late=zdata_late';
frY.zscore_late_up=zdata_late_up';
frY.zscore_late_down=zdata_late_down';

%%
% - - - - - - - -
% - - - Computation of the Depth for Learners
% - - - - - - -

% beforehand you need to have a channel map file which is coming with your
% recording setup. It tells you the position of the different active site
% on your recording probe and allow you to determine the depth.

for i = 1:size(frM,1)
    if isnan(frM.Id_horr(i))
        load(fullfile(sort_masters_horridge{frM.Recording(i)},'chanMap.mat'));
        clu = frM.Id_horr(i);
        idepth = get_cluster_depth(sort_masters_horridge{frM.Recording(i)}, clu);
    else
        load(fullfile(sort_masters_before{frM.Recording(i)},'chanMap.mat'));
        clu = frM.Id_horr(i);
        idepth = get_cluster_depth(sort_masters_before{frM.Recording(i)}, clu);
    end
    depths_M(i) = abs(idepth.depth' - ycoords(masters_depth(frM.Recording(i))));
end
frM.depth = depths_M';

% - - - - - - - -
% - - - Computation of the Depth for Controls
% - - - - - - -
for i = 1:size(frY,1)
    if isnan(frY.Id_horr(i))
        load(fullfile(sort_yokes_horridge{frY.Recording(i)},'chanMap.mat'));
        clu = frY.Id_horr(i);
        idepth = get_cluster_depth(sort_yokes_horridge{frY.Recording(i)}, clu);
    else
        clu = frY.Id_horr(i);
        load(fullfile(sort_yokes_before{frY.Recording(i)},'chanMap.mat'));
        idepth = get_cluster_depth(sort_yokes_before{frY.Recording(i)}, clu);
    end
        depths_Y(i) = abs(idepth.depth' - ycoords(yokes_depth(frY.Recording(i))));

    
end
frY.depth = depths_Y';

%%
% - - - - - - - - 
% Identify Dead neurons for Learners
% - - - - - - - - - 

% Dead neurons display an activity that continously decreasing during the
% time of the recording until it reaches 0 and disappears for the following
% time of the recording

% Number of the recording studied
recordingsM=find(frM.Recording==8);

% Vizualisation display to be able to manually curate alive neurons from
% dead neurons
for i=1:length(recordingsM)
    j=recordingsM(i);
    
    if frM.Fr_average_after(j) < 1
        
            figure(j)
            plot([frM.Fr_spont{j} frM.Fr_early{j} frM.Fr_middle{j} frM.Fr_late{j} frM.Fr_after{j}])
            hold on
            line([600 600],[0 30],'Color','red')
            iter=600+length(frM.Fr_early{j});
            line([iter iter],[0 30],'Color','black')
            iter=iter+length(frM.Fr_middle{j});
            line([iter iter],[0 30],'Color','black')
            iter=iter+length(frM.Fr_late{j});
            line([iter iter],[0 30],'Color','red')
            iter=iter+600;
            line([iter iter],[0 30],'Color','red')
            hold off
        
    end
    
end


% - - - - - - - - 
% Identify Dead neurons for Controls
% - - - - - - - - - 
% Number of the recording studied
recordingsY=find(frY.Recording==7);

for i=1:length(recordingsY)
    j=recordingsY(i);
    
    if frY.Fr_average_after(j) < 1
        
            figure(j)
            plot([frY.Fr_spont{j} frY.Fr_early{j} frY.Fr_middle{j} frY.Fr_late{j} frY.Fr_after{j}])
            hold on
            line([600 600],[0 30],'Color','red')
            iter=600+length(frY.Fr_early{j});
            line([iter iter],[0 30],'Color','black')
            iter=iter+length(frY.Fr_middle{j});
            line([iter iter],[0 30],'Color','black')
            iter=iter+length(frY.Fr_late{j});
            line([iter iter],[0 30],'Color','red')
            iter=iter+600;
            line([iter iter],[0 30],'Color','red')
            hold off
        
    end
    
end