% Takeoka Lab - NERF empowered by imec, KU Leuven and VIB
% Author: Charlotte Bichara & Simon Lavaud
% 2023

% Load the data base or recording : ---------------------------------------
% 'sort_masters_horridge' (WT or CNO) are the directory of the sorted data.
% use 'loadKSdir' to extract the 'spk' structure containing information on
% spike time clusters etc... Available upon request. 
% -------------------------------------------------------------------------

% load the frM or frY tables of interest. if table  _En1, _vGat or _vGlut,
% rename the table to frM or frY, but load the appropriate database

%Understanding the subplot(1,2,2)
% red line = zscore, black line = firing rate
% if square = green : response passes jitter and reliability criterium
% (direct)
% if square = blue : response passes only reliability criterium (indirect)
% if square = red : response passes only zscore increase criteria (upregulation)

% PARAMETER FOR DIRECT RESPONSE
threshold_zscore=5; % threshold zscore to be considered as a potential response
threshold_jitter= 0.001; %1ms
threshold_relia= 0.05; % Reliability minimal
fenetre=5; % +/- window for Jitter (total 1ms)

% - - - - - - - - - - - -
% ANALYSIS MASTER/Yoked
% - - - - - - - - - - - -

reliabilitylat=[];
meanlat=[];
jitterlat=[];

spike_count_all={};
spike_latency_all={};
spike_count=[];
spike_latency=[];

first_direct=0;
peak_reponse_all=[];
modulation_all=[];
modulation_value_all=[];
modulation_p_value_all=[];

Con_switch = 0; % Not Used in the paper
Condition = 1; % Select if units from : 1 Master, 0 Yoked, 3 CNO

% Choose unit(s) of interest : Line number for frM or frY tables, Value in
% column 'Line' in frM_UOI or frY_UOI-----------------------------------
UOI= 984; % (ex: 984 = illustration for the paper, Figure 4B Proprioception)
% ----------------------------------------------------------------------

Rec_Court_ONLY_conditioning = 0; % Length of recording different for recording 3, = 3 for units of rec 3 of Ptf1a_CNO

% If analysis of the electro=tagging
tagging =0; % Not used for the paper

  
close all
for k=1:numel(UOI)

    window_test=[-0.01 0.06]; % Window of search peri stimulus (0)
    gw = gausswin(round(20),3);
    nb_direct_response=0;
    latt_focus_sum=[];

    i= UOI(k);
    QQ=i;
    clc

    % - - - - -  - - - - -
    % - - - - -  - - - - -

    %£££%
    if Condition == 1
        spk = loadKSdir(sort_masters_horridge{frM.Recording(i)});
    elseif Condition == 0
        spk = loadKSdir(sort_yokes_horridge{frY.Recording(i)});
    elseif Condition == 3
        spk = loadKSdir(sort_masters_horridge_CNO{frM_Ptf1a_CNO.Recording(i)});
    end

    %£££%
    if Condition == 1
        NOI=frM.Id_horr(i);
    elseif Condition == 0
        NOI=frY.Id_horr(i);
    elseif Condition == 3
        NOI=frM_Ptf1a_CNO.Id_horr(i);
    end

    [a b]=find(spk.clu==NOI);

    %£££%
    if Condition == 1
        load(fullfile(sort_masters_horridge{frM.Recording(i)},'events.mat'),'events');
    elseif Condition == 0
        load(fullfile(sort_yokes_horridge{frY.Recording(i)},'events.mat'),'events');
    elseif Condition == 3
        load(fullfile(sort_masters_horridge_CNO{frM_Ptf1a_CNO.Recording(i)},'events.mat'),'events');
    end

    evnames = {events.type};
    evtid = find(contains(evnames,'horridge'));
    if isempty(evtid)==1
        evtid = find(contains(evnames,'opto'));
    end

    E=events.onsets;
    if Rec_Court_ONLY_conditioning == 1
        pos_cut=length(events.onsets);
    else
        [cut pos_cut]=max(diff(E));
    end

    peaks_diffE = find(diff(E)>500);

    %£££%
    if Condition == 1
        if long_long_masters(frM.Recording(i))==1 %long_long_yokes(frY.Recording(i))==1 %
            if tagging == 0
                if Con_switch == 0
                    Ev=E(1:peaks_diffE(1));
                else
                    Ev=E(peaks_diffE(2):peaks_diffE(3));
                end
            else
                Ev=E(StartTAG:StopTAG);
            end
        else
            Ev=E(1:pos_cut);
        end
    elseif Condition == 0
        if long_long_yokes(frY.Recording(i))==1
            if tagging == 0
                if Con_switch == 0
                    Ev=E(1:peaks_diffE(1));
                else
                    Ev=E(peaks_diffE(2):peaks_diffE(3));
                end
            else
                Ev=E(StartTAG:StopTAG);
            end
        else
            Ev=E(1:pos_cut);
        end
    elseif Condition == 3
        if long_long_masters_CNO(frM_Ptf1a_CNO.Recording(i))==1
            if tagging == 0
                if Con_switch == 0
                    Ev=E(1:peaks_diffE(1));
                else
                    Ev=E(peaks_diffE(2):peaks_diffE(3));
                end
            else
                Ev=E(StartTAG:StopTAG);
            end
        else
            Ev=E(1:pos_cut);
        end
    end

    spikeTimes = spk.st(a);
    eventTimes = Ev;

    % CALCULATION OF THE PSTH
    [psth, bins, rasterX, rasterY, spikeCounts, ba] = psthAndBA(spikeTimes, eventTimes, window_test, 0.0001);

    opt = statset('Maxiter',1000,'Display','final');

    [W H]=nnmf(ba,2,'options',opt,'algorithm','als');

    % plot(H(1,:)+H(2,:))

    [tr,b] = find(ba);
    [rasterX,yy] = rasterize(bins(b));
    rasterY = yy+reshape(repmat(tr',3,1),1,length(tr)*3); % yy is of the form [0 1 NaN 0 1 NaN...] so just need to add trial number to everything

    % scale the raster ticks
    rasterScale = floor(numel(eventTimes)/100);
    rasterY(2:3:end) = rasterY(2:3:end)+rasterScale;

    figure(i)


    sgtitle(['Unit number : ' num2str(i)]);



    % -------------

    subplot(2,1,1)
    plot(rasterX,rasterY, 'k');
    ylim([0 numel(spikeCounts)])
    xlim(window_test)
    %ylabel('event number');
    %xlabel('time (sec)');
    makepretty;
    h = rectangle('Position', [0, 0, 0.001,  numel(spikeCounts)], ...
        'Curvature', 0.2, ...
        'FaceColor', [1, 0, 0, 0.5], ...
        'EdgeColor', [1, 0, 0, 0.5]);
    [t,s] = title('Raster plot') ;
    t.FontSize = 10;
    grid on
    grid minor

    % -------------
    subplot(2,1,2)
    [t,s] =title('PSTH') ;
    t.FontSize = 10;
    smWin = gw./sum(gw);
    hold on
    % smooth ba
    baSm = conv2(smWin,1,ba', 'same')'./0.001;
    psthSm = mean(baSm);
    plot(bins, psthSm,'k','LineWidth',1)

    %Value Psth before shock excluding 2ms before shock
    val_bef_shock_mean=mean(psthSm(1,1:80));
    val_bef_shock_std=std(psthSm(1,1:80));

    if val_bef_shock_std==0
        val_bef_shock_std=0.1;
    end

    % ZSCORE OF THE PSTH
    psthSm_zscore=(psthSm(110:700)-val_bef_shock_mean) / val_bef_shock_std;
    zscore_mean(k) = mean(psthSm_zscore);
    zscore_max(k) = max(psthSm_zscore);
    zscore_AUC(k) = trapz(psthSm_zscore);
    hold on
    plot(bins(110:700), psthSm_zscore,'r','LineWidth',1)
    hold off
    occur=[];
    for i=110:679 %To avoid index to get out of range
        occur(i)=0;
        if psthSm_zscore(i-109)>threshold_zscore            
            occur(i)=1;
        end
    end
    
 
    X=occur;
    d = [true, diff(X) ~= 0, true];  % TRUE if values change
    n = diff(find(d));
    
    jitterlat=[];
    meanlat=[];
    reliabilitylat=[];
    meanlat2=[];
    jitterlat2=[];
    reliabilitylat2=[];
    reliabilitylat_focus=[];
    jitterlat_focus=[];
    peak_reponse=[];
    density_response=[];
    
    %COMPUTE HOW MANY CHANGE IN ZSCORE
    if numel(n)>1 & numel(n)<699
        
        n2(1)=n(1);
        for i=2:numel(n)
            n2(i)=n(i)+n2(i-1);
        end
        
        idx = find(bins>0.002);
        for i=1:floor(numel(n)/2)
            
            %Start of the psth zscorepeak FOR EACH WINDOW
            inter_focus=n2(i*2-1):n2(i*2);
            
            % If window too small
            if numel(inter_focus)<3
                
                idx2 = inter_focus;
                idx1 = inter_focus;
                a=[];
                
            else
                
                % FIND PSTH PEAK WITHIN THE WINDOW
                [a b]=findpeaks(psthSm(1,inter_focus));
                
            end
            
            % IF NO PEAK FOUND
            if isempty(a)
                
                idx2 = inter_focus;
                idx1 = inter_focus;
            else
                
                % MAX COMPUTED IF MULTIPLE PEAK
                [aprime bprime]=max(a);
                
                if isempty(bprime)
                    bprime=1;
                end
                
                %LATENCE OF THE PEAK
                peak_reponse(i)=bins(inter_focus(1)+round(b(bprime)));
                
                % WINDOW FOR RELIABILITY
                idx1 = [inter_focus(1)+round(b(bprime))-10:inter_focus(1)+round(b(bprime))+10];
                % WINDOW FOR JITTER
                idx2 = [inter_focus(1)+round(b(bprime))-fenetre:inter_focus(1)+round(b(bprime))+fenetre];
                
            end
            
            h = rectangle('Position', [bins(idx2(1)), 0, 0.002,  3], ...
                'Curvature', 0.2, ...
                'FaceColor', [1, 0, 0, 0.9], ...
                'EdgeColor', [1, 0, 0, 0.9]);
            
            
            latt_focus=[];
            o=0;
            for j = 1:size(ba,1) %first on all trials
                ll = find(ba(j,idx1)==1);
                if ~isempty(ll)
                    o=o+1;
                    latt_focus(j) = bins(idx1(ll(1)));
                    density_response(j)=j;
                else
                    latt_focus(j) = NaN;
                    density_response(j)=NaN;
                end
                
            end
            
            
            latt_focus_jitter=[];
            o=0;
            for j = 1:size(ba,1) %first on all trials
                ll = find(ba(j,idx2)==1);
                if ~isempty(ll)
                    o=o+1;
                    latt_focus_jitter(j) = bins(idx2(ll(1)));
                    density_response_jitter(j)=j;
                else
                    latt_focus_jitter(j) = NaN;
                    density_response_jitter(j)=NaN;
                end
                
            end
            
            
            if isempty(a)
                
                reliabilitylat(i)=0;
                jitterlat(i)=0;
                reliabilitylat_focus(i)=0;
                jitterlat_focus(i)=0;
                
            else
                
                [c d]=max(ksdensity(density_response,1:size(ba,1)));
                
                if d>size(ba,1)-50
                    d=size(ba,1)-51;
                    
                elseif  d<51
                    
                    d=51;
                end
                
                if size(ba,1)<100
                    
                     reliabilitylat_focus(i)=sum(~isnan(latt_focus(:)))/numel(latt_focus);
                else

                    reliabilitylat_focus(i)=sum(~isnan(latt_focus(d-50:d+50)))/numel(d-50:d+50);
                
                end

                jitterlat_focus(i)=std(latt_focus_jitter(:),'omitnan');
                
                latt_focus_sum(i,:)=latt_focus(:);
                
                reliabilitylat(i)=sum(~isnan(latt_focus))/ numel(latt_focus);
                jitterlat(i)=std(latt_focus,'omitnan');
                meanlat(i)=mean(latt_focus,'omitnan');
                
                
            end
            
            latt_focus_sum(i,:)=latt_focus(:);
            
            reliabilitylat(i)=sum(~isnan(latt_focus))/ numel(latt_focus);
            jitterlat(i)=std(latt_focus,'omitnan');
            meanlat(i)=mean(latt_focus,'omitnan');
            
            % - - - - - - - - - - - - - - - - - -
            % - - - - - - - - -
            % Indirect response
            % - - - - - - - - -
            % - - - - - - - - - - - - - - - - - -
            
            if reliabilitylat_focus(i)>threshold_relia & jitterlat_focus(i)>threshold_jitter
                cue_response_indirect(k)=1;
                
                h = rectangle('Position', [bins(idx2(1)), 0, 0.002,  3], ...
                    'Curvature', 0.2, ...
                    'FaceColor', [0, 0, 1, 0.9], ...
                    'EdgeColor', [0, 0, 1, 0.9]);
                
            end
            
            % - - - - - - - - - - - - - - - - - -
            % - - - - - - - - -
            % Direct response
            % - - - - - - - - -
            % - - - - - - - - - - - - - - - - - -
            
            if reliabilitylat_focus(i)>threshold_relia & jitterlat_focus(i)<threshold_jitter
                
                % Iterate the number of direct response found
                nb_direct_response=nb_direct_response+1;
                
                peak_reponse_all(nb_direct_response)=peak_reponse(i);
                
                
                if first_direct==0
                    lat_first_direct_response=peak_reponse(i);
                    jitter_first_direct_response=jitterlat_focus(i);
                    rel_first_direct_response=reliabilitylat_focus(i);
                    first_direct=1;
                end
                
                h = rectangle('Position', [bins(idx2(1)), 0, 0.002,  3], ...
                    'Curvature', 0.2, ...
                    'FaceColor', [0, 1, 0, 0.9], ...
                    'EdgeColor', [0, 1, 0, 0.9]);
                
                clearvars spikeCounts_mdl psth_mdl bins_mdl ba_mdl
                % Bin of intest
                window_mdl=bins(idx1);
                
                % Compute the number of spikes elicited on the specific
                % windows of interest
                [psth_mdl, bins_mdl, rasterX_mdl, rasterY_mdl, spikeCounts_mdl, ba_mdl] = psthAndBA(spikeTimes, eventTimes, [window_mdl(1) window_mdl(end)], 0.0001);

            end
            
        end
        
    end
    
    grid on
    grid minor
    
    
end
