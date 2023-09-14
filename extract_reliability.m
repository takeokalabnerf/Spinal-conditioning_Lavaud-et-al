% Takeoka Lab - NERF empowered by imec, KU Leuven and VIB
% Authors: Simon Lavaud (main) and Charlotte Bichara
% 2023

% Creates additional columns in frM_UOI and frY_UOI tables containing:
% data used for the paper :
% SetOfSpikeStart : % of spikes/stimulation for the first X% (X=perc_shock)
% SetOfSpikeEnd : % of spikes/stimulation for the last X% (X=perc_shock)
% ReliabilityDelta : difference SetOfSpikeEnd-SetOfSpikeStart

clc


% Load :

% Load database of interest, contain the paths to the sorted and raw data files 
recording_database_updated_Engrailed
recording_database_updated_vGlut
recording_database_updated_vGat
recording_database_updated_20220718 % WT and Ptf1a mice

% Load frM_UOI and frY_UOI, frM and frY for the genotype of interest

% Percentage of shock to take into accout at the beginning/end 
% for reliability calculation, % of shocks analysed (firsts and lasts)
perc_shock=10;


% Parameters to set :

% Velocities and distance of conduction parameters: -----------------------
% -------- > Used for paper
sABmin = 12; sABmax = 20;
sAdmin = 0.1; sAdmax = 12; 
sCmin = 0.01; sCmax = 0.09;
sPromin = 20; sPromax = 120; 
Lenght_leg=5; % average length from recording site to the furthest stimulation electrode (data not shown)
% -------------------------------------------------------------------------


% Choose the population to analyse. if all = 0, analyse is done on full ---
% second order set:
Ptf1aON_analysis = 0;
vGlutON_analysis = 0;
vGatON_analysis = 0;
% -------------------------------------------------------------------------





% Window of responses to sensory stimulation
ABstart = (Lenght_leg/(sABmax*100))+0.001; ABstop = (Lenght_leg/(sABmin*100))+0.001;
Adstart = (Lenght_leg/(sAdmax*100))+0.001; Adstop = (Lenght_leg/(sAdmin*100))+0.001;
Cstart =  (Lenght_leg/(sCmax*100))+0.001; Cstop = (Lenght_leg/(sCmin*100))+0.001;
Prostart = (Lenght_leg/(sPromax*100))+0.001; Prostop = (Lenght_leg/(sPromin*100))+0.001;

% - - - - - - - - - - - 
% - - -  EXTRACT UOI MASTER - - - - 
% - - - - - - - - - - - 

clearvars unit_responsive
unit_responsive=find(frM_UOI.depth);

iter_on_lat=1;
iter_on_C=1;
iter_on_P=1;
iter_on_AB=1;
iter_on_Ad=1;

spikeCounts_AB=NaN(10000, 100);
spikeCounts_C=NaN(10000, 100);
spikeCounts_P=NaN(10000, 100);
spikeCounts_Ad=NaN(10000, 100);

spikeCounts_AB_Y=NaN(10000, 100);
spikeCounts_C_Y=NaN(10000, 100);
spikeCounts_P_Y=NaN(10000, 100);
spikeCounts_Ad_Y=NaN(10000, 100);

Lat_hist=[];
Mod_hist=[];
LabelFiber_hist=[];
Value_hist=[];
Cat_hist=[];
Slope_hist=[];
Pvalue_hist=[];
ValueSU_hist=[];
Cross_hist=[];
barPTslope_M=[];
Fr_aver_hist=[];
Depth_hist=[];
Ptf1a_hist=[];
vGlutON_hist=[];
vGatON_hist=[];

for i=1:numel(unit_responsive)
    
    Lat_iter=cell2mat(frM_UOI.latency_all(unit_responsive(i)));
    Mod_iter=sign(cell2mat(frM_UOI.modulation_all(unit_responsive(i))));
    Value_iter=cell2mat(frM_UOI.modulation_all(unit_responsive(i)));

    Cat_iter=ones(1,numel(Lat_iter))*frM_UOI.Category(unit_responsive(i));
    Slope_iter=cell2mat(frM_UOI.modulation_value_all(unit_responsive(i)));
    Pvalue_iter=cell2mat(frM_UOI.modulation_p_value_all(unit_responsive(i)));
    ValueSU_iter=ones(1,numel(Mod_iter))*unit_responsive(i);
    Ptf1aON_iter=ones(1,numel(Lat_iter))*frM_UOI.Ptf1aON(unit_responsive(i));
    vGlutON_iter=ones(1,numel(Lat_iter))*frM_UOI.vGlutON(unit_responsive(i));
    vGatON_iter=ones(1,numel(Lat_iter))*frM_UOI.vGatON(unit_responsive(i));

    Depth_iter=ones(1,numel(Mod_iter))*frM_UOI.depth(unit_responsive(i));

    Lat_hist=[Lat_hist Lat_iter];
    Mod_hist=[Mod_hist Mod_iter];
    Value_hist=[Value_hist Value_iter];
    Ptf1a_hist=[Ptf1a_hist Ptf1aON_iter];
    vGlutON_hist=[vGlutON_hist vGlutON_iter];
    vGatON_hist=[vGatON_hist vGatON_iter];
    Cat_hist=[Cat_hist Cat_iter];
    Slope_hist=[Slope_hist Slope_iter];
    Pvalue_hist=[Pvalue_hist Pvalue_iter];
    ValueSU_hist=[ValueSU_hist ValueSU_iter];

    Depth_hist=[Depth_hist Depth_iter];
    
end 



% - - - - - - - - - - - 
% - - -  EXTRACT UOI YOKED - - - - 
% - - - - - - - - - - -

unit_responsive=find( frY_UOI.depth);

Lat_histY=[];
Mod_histY=[];
LabelFiber_histY=[];
Value_histY=[];
Cat_histY=[];
Slope_histY=[];
Pvalue_histY=[];
ValueSU_histY=[];
Depth_histY=[];
Cross_histY=[];
Ptf1a_histY=[];
vGlutON_histY=[];
vGatON_histY=[];
clc
for i=1:numel(unit_responsive)
    
    Lat_iter=cell2mat(frY_UOI.latency_all(unit_responsive(i)));
    Mod_iter=sign(cell2mat(frY_UOI.modulation_all(unit_responsive(i))));
    Value_iter=cell2mat(frY_UOI.modulation_all(unit_responsive(i)));

    Cat_iter=ones(1,numel(Lat_iter))*frY_UOI.Category(unit_responsive(i));
    Slope_iter=cell2mat(frY_UOI.modulation_value_all(unit_responsive(i)));
    Pvalue_iter=cell2mat(frY_UOI.modulation_p_value_all(unit_responsive(i)));
    ValueSU_iter=ones(1,numel(Mod_iter))*unit_responsive(i);
    Depth_iter=ones(1,numel(Mod_iter))*frY_UOI.depth(unit_responsive(i));
    Ptf1aON_iter=ones(1,numel(Lat_iter))*frY_UOI.Ptf1aON(unit_responsive(i));
    vGlutON_iter=ones(1,numel(Lat_iter))*frY_UOI.vGlutON(unit_responsive(i));
    vGatON_iter=ones(1,numel(Lat_iter))*frY_UOI.vGatON(unit_responsive(i));

    Lat_histY=[Lat_histY Lat_iter];
    Mod_histY=[Mod_histY Mod_iter];
    Value_histY=[Value_histY Mod_iter];
    
    Cat_histY=[Cat_histY Cat_iter];
    Slope_histY=[Slope_histY Slope_iter];
    Pvalue_histY=[Pvalue_histY Pvalue_iter];
    ValueSU_histY=[ValueSU_histY ValueSU_iter];
    Depth_histY=[Depth_histY Depth_iter];
    Ptf1a_histY=[Ptf1a_histY Ptf1aON_iter];
    vGlutON_histY=[vGlutON_histY vGlutON_iter];
    vGatON_histY=[vGatON_histY vGatON_iter];
    
end


%%

% - - - - - - - - - - - 
% - - -  EXTRACTION RELIABILITY MASTER - - - - 
% - - - - - - - - - - -

Lat_AB=[];
Rec_AB=[];
Id_AB=[];
Lat_AB_Y=[];
Rec_AB_Y=[];
Id_AB_Y=[];

Lat_C=[];
Rec_C=[];
Id_C=[];
Lat_C_Y=[];
Rec_C_Y=[];
Id_C_Y=[];

Lat_P=[];
Rec_P=[];
Id_P=[];
Lat_P_Y=[];
Rec_P_Y=[];
Id_P_Y=[];

Lat_Ad=[];
Rec_Ad=[];
Id_Ad=[];
Lat_Ad_Y=[];
Rec_Ad_Y=[];
Id_Ad_Y=[];

SoP_C_Y=[];
SoP_C=[];
SoP_Ad_Y=[];
SoP_Ad=[];
SoP_AB_Y=[];
SoP_AB=[];
SoP_P_Y=[];
SoP_P=[];
CheckPtf1aAd=[];
CheckPtf1aAB=[];
CheckPtf1aAdY=[];
CheckPtf1aABY=[];

CheckvGlutONAd=[];
CheckGlutONAB=[];
CheckGlutONAdY=[];
CheckGlutONABY=[];

CheckGatONAd=[];
CheckGatONAB=[];
CheckGatONAdY=[];
CheckGatONABY=[];



for i=1:ValueSU_hist(end) %Number SU responsive iterate on i
    
    %Is that unit responsive more that once ? 
    % If yes number of response = HowManyResponse
    HowManyResponse = sum(ismember(ValueSU_hist,i));
    
    ReliabilityBinary_array=[];
    PerOfChangeM_indiv_array=[];
    Delta_RelM_indiv_array=[];
    ReliabilityBinary_arrayCha=[];
    SetOfSpikeStart=[];
    SetOfSpikeEnd=[];

        
    % Iterate on HowManyResponse on j
    for j=0:(HowManyResponse-1)
        
        % Get the line of the UOI from the Line list
        Line=frM_UOI.Line(i);
        
        % - - - - Check which recording:
        
        if frM_UOI.Recording(i)<99 % Ptf1a recording
            
            % Extract PSTH
            % /!\ name is sort_masters_horridge_classic !!
            spk = loadKSdir(sort_masters_horridge_classic{frM.Recording(Line)});
            
            NOI=frM.Id_horr(Line);
            [a , ~]=find(spk.clu==NOI);
            
            load(fullfile(sort_masters_horridge_classic{frM.Recording(Line)},'events.mat'),'events');
            evnames = {events.type};
            evtid = find(contains(evnames,'horridge'));
            if isempty(evtid)==1
                evtid = find(contains(evnames,'opto'));
            end
            
        elseif frM_UOI.Recording(i)<999 % En1 recording
            
            % Extract PSTH
            spk = loadKSdir(sort_masters_horridge_en1{frM_En1.Recording(Line)});
            
            NOI=frM_En1.Id_horr(Line);
            [a b]=find(spk.clu==NOI);
            
            load(fullfile(sort_masters_horridge_en1{frM_En1.Recording(Line)},'events.mat'),'events');
            evnames = {events.type};
            evtid = find(contains(evnames,'horridge'));
            if isempty(evtid)==1
                evtid = find(contains(evnames,'opto'));
            end
            
        elseif frM_UOI.Recording(i)<9999 % Vglut recording
            
            % Extract PSTH
            spk = loadKSdir(sort_masters_horridge_vglut{frM_vGlut.Recording(Line)});
            
            NOI=frM_vGlut.Id_horr(Line);
            [a b]=find(spk.clu==NOI);
            
            load(fullfile(sort_masters_horridge_vglut{frM_vGlut.Recording(Line)},'events.mat'),'events');
            evnames = {events.type};
            evtid = find(contains(evnames,'horridge'));
            if isempty(evtid)==1
                evtid = find(contains(evnames,'opto'));
            end
         
        else % Vgat recording
            
            % Extract PSTH
            spk = loadKSdir(sort_masters_horridge_vgat{frM_vGat.Recording(Line)/1});
            
            NOI=frM_vGat.Id_horr(Line);
            [a b]=find(spk.clu==NOI);
            
            load(fullfile(sort_masters_horridge_vgat{frM_vGat.Recording(Line)/1},'events.mat'),'events');
            evnames = {events.type};
            evtid = find(contains(evnames,'horridge'));
            if isempty(evtid)==1
                evtid = find(contains(evnames,'opto'));
            end
            
        end
        
        % Extract event
        E=events.onsets;
        
        if E(end)<6000 % Recording shorter than 1hour40 - only horridge without switch
            
            remove_outlier=find(diff(E)<50);
            E2=E(remove_outlier(1):end);
            
            [cut pos_cut]=max(diff(E2));
            %Start_of_E=find(E>850);
            
            eventTimes=E2(1:pos_cut);
            spikeTimes = spk.st(a);
            
        else % Recording with switch
            
            [cut pos_cut]=max(diff(E));
            
            peaks_diffE = find(diff(E)>500);
            
            Ev=E(1:peaks_diffE(1));

            spikeTimes = spk.st(a);
            eventTimes = Ev;
            
        end

        window_mdl=Lat_hist(iter_on_lat)-0.001:0.0001:Lat_hist(iter_on_lat)+0.001;
        
        [psth_mdl, bins_mdl, rasterX_mdl, rasterY_mdl, spikeCounts_mdl, ba_mdl] = psthAndBA(spikeTimes, eventTimes, [window_mdl(1) window_mdl(end)], 0.0001);
        
        % Just 0/1 (like here) or real number of spike
        spikeCounts_mdl=double(spikeCounts_mdl>0);
        
        %Here separated the shock in the first X and the last X%
        howlong=floor(numel(spikeCounts_mdl)/(100/perc_shock));
        
        % Store the spike elicited mean of the first X% of the shock
        % First set
        set_of_spike1=mean(spikeCounts_mdl(1:howlong,1));
        % Last set
        set_of_spike2=mean(spikeCounts_mdl(end-howlong:end,1));
    	
        SetOfSpikeStart=[SetOfSpikeStart set_of_spike1];
        SetOfSpikeEnd=[SetOfSpikeEnd set_of_spike2];
        
        frM_UOI.SetOfSpikeStart{i}=SetOfSpikeStart;
        frM_UOI.SetOfSpikeEnd{i}=SetOfSpikeEnd;
        
        PerOfChangeM_indiv=set_of_spike2*100./set_of_spike1-100;
        
        PerOfChangeM_indiv_array=[PerOfChangeM_indiv_array PerOfChangeM_indiv];
        
        frM_UOI.ReliabilityChange{i}=PerOfChangeM_indiv_array;
        
        if PerOfChangeM_indiv>10
            ReliabilityBinary=1;
        %frM_UOI.ReliabilityBinary(iter_on_lat)=1;
        elseif PerOfChangeM_indiv<-10
            ReliabilityBinary=-1;
        %frM_UOI.ReliabilityBinary(iter_on_lat)=-1;
        else
            ReliabilityBinary=0;
            %frM_UOI.ReliabilityBinary(iter_on_lat)=0;
        end
        ReliabilityBinary_array=[ReliabilityBinary_array ReliabilityBinary];

        frM_UOI.ReliabilityBinary{i}=ReliabilityBinary_array;


        % Calculation for Delta (Cha)
        PerOfChangeM_indivDelta=set_of_spike2-set_of_spike1;
        
        Delta_RelM_indiv_array=[Delta_RelM_indiv_array PerOfChangeM_indivDelta];
        
        frM_UOI.ReliabilityDelta{i}=Delta_RelM_indiv_array;

        if PerOfChangeM_indiv>0
            ReliabilityBinaryCha=1;
        %frM_UOI.ReliabilityBinary(iter_on_lat)=1;
        elseif PerOfChangeM_indiv<0
            ReliabilityBinaryCha=-1;
        %frM_UOI.ReliabilityBinary(iter_on_lat)=-1;
        else
            ReliabilityBinaryCha=0;
            %frM_UOI.ReliabilityBinary(iter_on_lat)=0;
        end
        ReliabilityBinary_arrayCha=[ReliabilityBinary_arrayCha ReliabilityBinaryCha];

        frM_UOI.ReliabilityBinaryCha{i}=ReliabilityBinary_arrayCha;

        % Store based on the latency for further analysis
        % Fiber type = 1 if C / 2 if Ad / 3 if AB / 4 if P
        if Lat_hist(iter_on_lat)>=Cstart
            spikeCounts_C(1:numel(spikeCounts_mdl),iter_on_C)=spikeCounts_mdl;
            Lat_C(iter_on_C)=Lat_hist(iter_on_lat);
            Rec_C(iter_on_C)=frM_UOI.Recording(i);
            Id_C(iter_on_C)=frM.Id_horr(i);
            Id_respC(iter_on_C)=frM_UOI.Line(i);
            SoP_C(1,iter_on_C)=set_of_spike1;
            SoP_C(2,iter_on_C)=set_of_spike2;
            CheckPtf1aC(iter_on_C)=Ptf1a_hist(iter_on_lat);
            CheckvGatONC(iter_on_C)=vGatON_hist(iter_on_lat);
            CheckvGlutONC(iter_on_C)=vGlutON_hist(iter_on_lat);
            iter_on_C=iter_on_C+1;
            CheckFiber(iter_on_lat)=1;

            
        elseif Lat_hist(iter_on_lat)>= Adstart & Lat_hist(iter_on_lat) < Adstop
            spikeCounts_Ad(1:numel(spikeCounts_mdl),iter_on_Ad)=spikeCounts_mdl;
            Lat_Ad(iter_on_Ad)=Lat_hist(iter_on_lat);
            Rec_Ad(iter_on_Ad)=frM_UOI.Recording(i);
            Id_Ad(iter_on_Ad)=frM.Id_horr(i);
            Id_respAd(iter_on_Ad)=frM_UOI.Line(i);
            SoP_Ad(1,iter_on_Ad)=set_of_spike1;
            SoP_Ad(2,iter_on_Ad)=set_of_spike2;
            CheckPtf1aAd(iter_on_Ad)=Ptf1a_hist(iter_on_lat);
            CheckvGatONAd(iter_on_Ad)=vGatON_hist(iter_on_lat);
            CheckvGlutONAd(iter_on_Ad)=vGlutON_hist(iter_on_lat);
            iter_on_Ad=iter_on_Ad+1;
            CheckFiber(iter_on_lat)=2;
            

        elseif Lat_hist(iter_on_lat)>=ABstart & Lat_hist(iter_on_lat) < ABstop
            spikeCounts_AB(1:numel(spikeCounts_mdl),iter_on_AB)=spikeCounts_mdl;
            Lat_AB(iter_on_AB)=Lat_hist(iter_on_lat);
            Rec_AB(iter_on_AB)=frM_UOI.Recording(i);
            Id_AB(iter_on_AB)=frM.Id_horr(i);
            Id_respAB(iter_on_AB)=frM_UOI.Line(i);
            SoP_AB(1,iter_on_AB)=set_of_spike1;
            SoP_AB(2,iter_on_AB)=set_of_spike2;
            CheckPtf1aAB(iter_on_AB)=Ptf1a_hist(iter_on_lat);
            CheckvGatONAB(iter_on_AB)=vGatON_hist(iter_on_lat);
            CheckvGlutONAB(iter_on_AB)=vGlutON_hist(iter_on_lat);
            iter_on_AB=iter_on_AB+1;
            CheckFiber(iter_on_lat)=3;
            
        elseif Lat_hist(iter_on_lat)>=Prostart & Lat_hist(iter_on_lat) < Prostop
            spikeCounts_P(1:numel(spikeCounts_mdl),iter_on_P)=spikeCounts_mdl;
            Lat_P(iter_on_P)=Lat_hist(iter_on_lat);
            Rec_P(iter_on_P)=frM_UOI.Recording(i);
            Id_P(iter_on_P)=frM.Id_horr(i);
            Id_respP(iter_on_P)=frM_UOI.Line(i);
            SoP_P(1,iter_on_P)=set_of_spike1;
            SoP_P(2,iter_on_P)=set_of_spike2;
            CheckPtf1aP(iter_on_P)=Ptf1a_hist(iter_on_lat);
            CheckvGatONP(iter_on_P)=vGatON_hist(iter_on_lat);
            CheckvGlutONP(iter_on_P)=vGlutON_hist(iter_on_lat);
            iter_on_P=iter_on_P+1;
            CheckFiber(iter_on_lat)=4;
        end
        
        iter_on_lat=iter_on_lat+1
        
    end
    
end



%%



% - - - - - - - - - - - 
% - - -  EXTRACTION RELIABILITY YOKED - - - - 
% - - - - - - - - - - -

iter_on_lat=1;
iter_on_C=1;
iter_on_P=1;
iter_on_AB=1;
iter_on_Ad=1;




for i=1:ValueSU_histY(end) %Number SU responsive
    
    %Is that unit responsive more that once ?
    HowManyResponse = sum(ismember(ValueSU_histY,i));
    
    ReliabilityBinary_array=[];
    PerOfChangeM_indiv_array=[];
    Delta_RelY_indiv_array=[];
    ReliabilityBinary_arrayCha=[];
    SetOfSpikeStart=[];
    SetOfSpikeEnd=[];
        
    for j=0:(HowManyResponse-1)
        

        Line=frY_UOI.Line(i);
        
        if frY_UOI.Recording(i)<99 %Engrailed
            
            % Extract PSTH
            spk = loadKSdir(sort_yokes_horridge_classic{frY.Recording(Line)/1});
            
            NOI=frY.Id_horr(Line);
            [a b]=find(spk.clu==NOI);
            
            load(fullfile(sort_yokes_horridge_classic{frY.Recording(Line)/1},'events.mat'),'events');
            evnames = {events.type};
            evtid = find(contains(evnames,'horridge'));
            if isempty(evtid)==1
                evtid = find(contains(evnames,'opto'));
            end
            
        elseif frY_UOI.Recording(i)<999
            
            % Extract PSTH
            spk = loadKSdir(sort_yokes_horridge_en1{frY_En1.Recording(Line)/1});
            
            NOI=frY_En1.Id_horr(Line);
            [a b]=find(spk.clu==NOI);
            
            load(fullfile(sort_yokes_horridge_en1{frY_En1.Recording(Line)/1},'events.mat'),'events');
            evnames = {events.type};
            evtid = find(contains(evnames,'horridge'));
            if isempty(evtid)==1
                evtid = find(contains(evnames,'opto'));
            end
            
        elseif frY_UOI.Recording(i)<9999
            
            % Extract PSTH
            spk = loadKSdir(sort_yokes_horridge_vglut{frY_vGlut.Recording(Line)/1});
            
            NOI=frY_vGlut.Id_horr(Line);
            [a b]=find(spk.clu==NOI);
            
            load(fullfile(sort_yokes_horridge_vglut{frY_vGlut.Recording(Line)/1},'events.mat'),'events');
            evnames = {events.type};
            evtid = find(contains(evnames,'horridge'));
            if isempty(evtid)==1
                evtid = find(contains(evnames,'opto'));
            end
            
        else
            
            % Extract PSTH
            spk = loadKSdir(sort_yokes_horridge_vgat{frY_vGat.Recording(Line)/1});
            
            NOI=frY_vGat.Id_horr(Line);
            [a b]=find(spk.clu==NOI);
            
            load(fullfile(sort_yokes_horridge_vgat{frY_vGat.Recording(Line)/1},'events.mat'),'events');
            evnames = {events.type};
            evtid = find(contains(evnames,'horridge'));
            if isempty(evtid)==1
                evtid = find(contains(evnames,'opto'));
            end
            
        end
        
        % Extract event
        E=events.onsets;
        
        if E(end)<6000 % Recording shorter than 1hour - only horridge without switch
            
            remove_outlier=find(diff(E)<50);
            E2=E(remove_outlier(1):end);
            
            [cut pos_cut]=max(diff(E2));
            %Start_of_E=find(E>850);
            
            eventTimes=E2(1:pos_cut);
            spikeTimes = spk.st(a);
            
        else
            
            [cut pos_cut]=max(diff(E));
            
            peaks_diffE = find(diff(E)>500);
            
            %if long_long_masters(frM.Recording(Line)/1)==1 %long_long_yokes(frY.Recording(i))==1 %
            Ev=E(1:peaks_diffE(1));
            %else
            %    Ev=E(1:pos_cut);
            %end
            
            spikeTimes = spk.st(a);
            eventTimes = Ev;
            
        end
        
        % CALCULATION OF THE PSTH BASED ON THE LATENCY
        %[psth, bins, rasterX, rasterY, spikeCounts, ba] = psthAndBA(spikeTimes, eventTimes, window_test, 0.0001);
        
        window_mdl=Lat_histY(iter_on_lat)-0.001:0.0001:Lat_histY(iter_on_lat)+0.001;
        
        [psth_mdl, bins_mdl, rasterX_mdl, rasterY_mdl, spikeCounts_mdl, ba_mdl] = psthAndBA(spikeTimes, eventTimes, [window_mdl(1) window_mdl(end)], 0.0001);
        
        % Just 0/1 or real number of spike
        spikeCounts_mdl=double(spikeCounts_mdl>0);
        
        %Here separated the shock in the first 50 and the last 50%
        howlong=floor(numel(spikeCounts_mdl)/(100/perc_shock));
        
        % Store the spike elicited mean of the first X% of the shock
        % First set
        set_of_spike1=mean(spikeCounts_mdl(1:howlong,1));
        % Last set
        set_of_spike2=mean(spikeCounts_mdl(end-howlong:end,1));
        
        SetOfSpikeStart=[SetOfSpikeStart set_of_spike1];
        SetOfSpikeEnd=[SetOfSpikeEnd set_of_spike2];
        
        frY_UOI.SetOfSpikeStart{i}=SetOfSpikeStart;
        frY_UOI.SetOfSpikeEnd{i}=SetOfSpikeEnd;
        
        PerOfChangeM_indiv=set_of_spike2*100./set_of_spike1-100;
        
        PerOfChangeM_indiv_array=[ PerOfChangeM_indiv_array PerOfChangeM_indiv];
        
        frY_UOI.ReliabilityChange{i}=PerOfChangeM_indiv_array;
        
        if PerOfChangeM_indiv>10
            ReliabilityBinary=1;
        %frM_UOI.ReliabilityBinary(iter_on_lat)=1;
        elseif PerOfChangeM_indiv<-10
            ReliabilityBinary=-1;
        %frM_UOI.ReliabilityBinary(iter_on_lat)=-1;
        else
            ReliabilityBinary=0;
            %frM_UOI.ReliabilityBinary(iter_on_lat)=0;
        end
        ReliabilityBinary_array=[ReliabilityBinary_array ReliabilityBinary];
        
        frY_UOI.ReliabilityBinary{i}=ReliabilityBinary_array;

        % Calculation for Delta (Cha) =============
        PerOfChangeY_indivDelta=set_of_spike2-set_of_spike1;
        
        Delta_RelY_indiv_array=[Delta_RelY_indiv_array PerOfChangeY_indivDelta];
        
        frY_UOI.ReliabilityDelta{i}=Delta_RelY_indiv_array;

        if PerOfChangeY_indivDelta>0
            ReliabilityBinaryCha=1;
        %frM_UOI.ReliabilityBinary(iter_on_lat)=1;
        elseif PerOfChangeY_indivDelta<0
            ReliabilityBinaryCha=-1;
        %frM_UOI.ReliabilityBinary(iter_on_lat)=-1;
        else
            PerOfChangeY_indivDelta=0;
            %frM_UOI.ReliabilityBinary(iter_on_lat)=0;
        end
        ReliabilityBinary_arrayCha=[ReliabilityBinary_arrayCha ReliabilityBinaryCha];

        frY_UOI.ReliabilityBinaryCha{i}=ReliabilityBinary_arrayCha;
        %=============================================

        % Store based on the latency
        if Lat_histY(iter_on_lat)>=Cstart
            spikeCounts_C_Y(1:numel(spikeCounts_mdl),iter_on_C)=spikeCounts_mdl;
            Lat_C_Y(iter_on_C)=Lat_histY(iter_on_lat);
            Rec_C_Y(iter_on_C)=frY_UOI.Recording(i);
            Id_C_Y(iter_on_C)=frY.Id_horr(i);
            Id_respC_Y(iter_on_C)=frY_UOI.Line(i);
            SoP_C_Y(1,iter_on_C)=set_of_spike1;
            SoP_C_Y(2,iter_on_C)=set_of_spike2;
            CheckPtf1aCY(iter_on_C)=Ptf1a_histY(iter_on_lat);
            CheckvGatONCY(iter_on_C)=vGatON_histY(iter_on_lat);
            CheckvGlutONCY(iter_on_C)=vGlutON_histY(iter_on_lat);
            iter_on_C=iter_on_C+1;
            CheckFiberY(iter_on_lat)=1;
            

        elseif Lat_histY(iter_on_lat)>=Adstart & Lat_histY(iter_on_lat)< Adstop
            spikeCounts_Ad_Y(1:numel(spikeCounts_mdl),iter_on_Ad)=spikeCounts_mdl;
            Lat_Ad_Y(iter_on_Ad)=Lat_histY(iter_on_lat);
            Rec_Ad_Y(iter_on_Ad)=frY_UOI.Recording(i);
            Id_Ad_Y(iter_on_Ad)=frY.Id_horr(i);
            Id_respAd_Y(iter_on_Ad)=frY_UOI.Line(i);
            SoP_Ad_Y(1,iter_on_Ad)=set_of_spike1;
            SoP_Ad_Y(2,iter_on_Ad)=set_of_spike2;
%             Ptf1aON_SoP_Ad_Y(1,iter_on_Ad)=set_of_spike1;
%             Ptf1aON_SoP_Ad_Y(2,iter_on_Ad)=set_of_spike2;
%             vGlutON_SoP_Ad_Y(1,iter_on_Ad)=set_of_spike1;
%             vGlutON_SoP_Ad_Y(2,iter_on_Ad)=set_of_spike2;
%             vGatON_SoP_Ad_Y(1,iter_on_Ad)=set_of_spike1;
%             vGatON_SoP_Ad_Y(2,iter_on_Ad)=set_of_spike2;
            CheckPtf1aAdY(iter_on_Ad)=Ptf1a_histY(iter_on_lat);
            CheckvGatONAdY(iter_on_Ad)=vGatON_histY(iter_on_lat);
            CheckvGlutONAdY(iter_on_Ad)=vGlutON_histY(iter_on_lat);
            iter_on_Ad=iter_on_Ad+1;
            CheckFiberY(iter_on_lat)=2;
            

        elseif Lat_histY(iter_on_lat)>=ABstart & Lat_histY(iter_on_lat)< ABstop
            spikeCounts_AB_Y(1:numel(spikeCounts_mdl),iter_on_AB)=spikeCounts_mdl;
            Lat_AB_Y(iter_on_AB)=Lat_histY(iter_on_lat);
            Rec_AB_Y(iter_on_AB)=frY_UOI.Recording(i);
            Id_AB_Y(iter_on_AB)=frY.Id_horr(i);
            Id_respAB_Y(iter_on_AB)=frY_UOI.Line(i);
            SoP_AB_Y(1,iter_on_AB)=set_of_spike1;
            SoP_AB_Y(2,iter_on_AB)=set_of_spike2;
            CheckPtf1aABY(iter_on_AB)=Ptf1a_histY(iter_on_lat);
            CheckvGatONABY(iter_on_AB)=vGatON_histY(iter_on_lat);
            CheckvGlutONABY(iter_on_AB)=vGlutON_histY(iter_on_lat);
            iter_on_AB=iter_on_AB+1;
            CheckFiberY(iter_on_lat)=3;
            

        elseif Lat_histY(iter_on_lat)>=Prostart & Lat_histY(iter_on_lat)< Prostop
            spikeCounts_P_Y(1:numel(spikeCounts_mdl),iter_on_P)=spikeCounts_mdl;
            Lat_P_Y(iter_on_P)=Lat_histY(iter_on_lat);
            Rec_P_Y(iter_on_P)=frY_UOI.Recording(i);
            Id_P_Y(iter_on_P)=frY.Id_horr(i);
            Id_respP_Y(iter_on_P)=frY_UOI.Line(i);
            SoP_P_Y(1,iter_on_P)=set_of_spike1;
            SoP_P_Y(2,iter_on_P)=set_of_spike2;
            CheckPtf1aPY(iter_on_P)=Ptf1a_histY(iter_on_lat);
            CheckvGatONPY(iter_on_P)=vGatON_histY(iter_on_lat);
            CheckvGlutONPY(iter_on_P)=vGlutON_histY(iter_on_lat);
            iter_on_P=iter_on_P+1;
            CheckFiberY(iter_on_lat)=4;
        end
        


        iter_on_lat=iter_on_lat+1
        
    end
    
end


% - - - - - - - - - -
%% ---- C FIBER ANALYSIS 
% - - - - - - - - - -

% - - - - - - -  - - - 
% - - - - - - -  - - - 
threshold=0.05;
IfInfWhichValue=1500;
% - - - - - - -  - - - 
% - - - - - - -  - - - 

DecM=[];
DecY=[];

incrM=[];
incrY=[];

PerOfChangeM=[];
PerOfChangeY=[];
ValueI_Y=[];
ValueI=[];
ValueD_Y=[];
ValueD=[];

figure
sgtitle('C fibers')
hold on

if vGlutON_analysis == 1
    PerOfChangeM=(SoP_C(2,find(CheckvGlutONC==1))-SoP_C(1,find(CheckvGlutONC==1)));
    PerOfChangeY=(SoP_C_Y(2,find(CheckvGlutONCY==1))-SoP_C_Y(1,find(CheckvGlutONCY==1)));
elseif vGatON_analysis == 1
    PerOfChangeM=(SoP_C(2,find(CheckvGatONC==1))-SoP_C(1,find(CheckvGatONC==1)));
    PerOfChangeY=(SoP_C_Y(2,find(CheckvGatONCY==1))-SoP_C_Y(1,find(CheckvGatONCY==1)));
elseif Ptf1aON_analysis == 1
    PerOfChangeM=(SoP_C(2,find(CheckPtf1aC==1))-SoP_C(1,find(CheckPtf1aC==1)));
    PerOfChangeY=(SoP_C_Y(2,find(CheckPtf1aCY==1))-SoP_C_Y(1,find(CheckPtf1aCY==1)));
else
    PerOfChangeM=(SoP_C(2,1:end)-SoP_C(1,1:end));
    PerOfChangeY=(SoP_C_Y(2,1:end)-SoP_C_Y(1,1:end));
end

a=1;
b=1;


for i=1:numel(PerOfChangeM)
    
    if isinf(PerOfChangeM(i))
    PerOfChangeM(i)=IfInfWhichValue;
    end

    if PerOfChangeM(i)>threshold % Increase
        subplot(2,3,1)
        line([0 1],[0 (PerOfChangeM(i))],'Color','r')
        axis([0 1 0 400])
        incrM(a)=PerOfChangeM(i);
        ValueI(a)=PerOfChangeM(i);
        a=a+1;
    elseif PerOfChangeM(i)<-threshold % Decrease
        subplot(2,3,2)
        line([0 1],[0 -(-PerOfChangeM(i))],'Color','r')
        axis([0 1 -100 0])
        DecM(b)=PerOfChangeM(i);
        ValueD(b)=PerOfChangeM(i);
        b=b+1;
    else % No Change
        subplot(2,3,3)
        line([0 1],[0 PerOfChangeM(i)],'Color','r')
        axis([0 1 -10 10])
    end
    
end

a=1;
b=1;
for i=1:numel(PerOfChangeY)
    if isinf(PerOfChangeY(i))
    PerOfChangeY(i)=IfInfWhichValue;
    end
    if PerOfChangeY(i)>threshold % Increase
        subplot(2,3,4)
        line([0 1],[0 (PerOfChangeY(i))],'Color','b')
        axis([0 1 0 400])
        ValueI_Y(a)=PerOfChangeY(i);
        incrY(a)=PerOfChangeY(i);
        a=a+1;
    elseif PerOfChangeY(i)<-threshold % Decrease
        subplot(2,3,5)
        line([0 1],[0 -(-PerOfChangeY(i))],'Color','b')
        axis([0 1 -100 0])
        ValueD_Y(b)=PerOfChangeY(i);
        DecY(b)=PerOfChangeY(i);
        b=b+1;
    else
        subplot(2,3,6)
        line([0 1],[0 PerOfChangeY(i)],'Color','b')
        axis([0 1 -10 10])
    end
    
end

PerOfChangeY_C=PerOfChangeY;
PerOfChangeM_C=PerOfChangeM;

if length(incrM)>1 & length(incrY)>1
    [a b]=ranksum(incrM,incrY);
end

figure
title(' reliability')
ax1 = nexttile;
pie(ax1,[size(SoP_C,2)-numel(DecM)-numel(incrM),numel(DecM),numel(incrM)]);
title('Master C')

ax2 = nexttile;
pie(ax2,[size(SoP_C_Y,2)-numel(DecY)-numel(incrY),numel(DecY),numel(incrY)]);
title('Yoked C')
legend('No change','Decrease','Increase')

chisquare(size(SoP_C,2)-numel(DecM)-numel(incrM), ...
    size(SoP_C,2), ...
    size(SoP_C_Y,2)-numel(DecY)-numel(incrY), ...
    size(SoP_C_Y,2))

chisquare(numel(DecM), ...
    size(SoP_C,2), ...
    numel(DecY), ...
    size(SoP_C_Y,2))

chisquare(numel(incrM), ...
    size(SoP_C,2), ...
    numel(incrY), ...
    size(SoP_C_Y,2))


figure
hold on
scatter(ones(numel(PerOfChangeM_C),1),PerOfChangeM_C,'MarkerEdgeColor','k','MarkerFaceColor','m','jitter','on','jitterAmount', .1);
scatter(2*ones(numel(PerOfChangeY_C),1),PerOfChangeY_C,'MarkerEdgeColor','k','MarkerFaceColor','b','jitter','on','jitterAmount', .1);
xlim([0 3])
alpha(0.5)


% Violin Plot
figure
subplot(1,2,1)
x=ValueD;
y=ValueD_Y;
scatter(ones(numel(x),1),x,'MarkerEdgeColor','k','MarkerFaceColor','m','jitter','on','jitterAmount', .1);
hold on
scatter(2*ones(numel(y),1),y,'MarkerEdgeColor','k','MarkerFaceColor','b','jitter','on','jitterAmount', .1);
xlim([0 3])
ylim([-1 0])
alpha(0.5)
title('C Decrease')
hold off

if length(x) >1 & length(y)> 1
    [h1 b1]=kstest( (x-mean(x))/std(x) );
    [h2 b2]=kstest( (y-mean(y))/std(y) );
    
    
    
    subplot(1,2,1)
    if h1==0 & h2==0 % data normal
        [c d]=ttest2(x,y);
        if d<0.001
            text(1,min(x)/2,3,'***','FontSize',20,'Color','red')
        elseif d<0.01
            text(1,min(x)/2,3,'**','FontSize',20,'Color','red')
            elseif d<0.05
            text(1,min(x)/2,3,'*','FontSize',20,'Color','red')
        else
            text(1,min(x)/2,3,'NS')
        end
    else
        [c d]=ranksum(x,y);
        if c<0.001
            text(1,min(x)/2,3,'***','FontSize',20,'Color','red')
        elseif c<0.01
            text(1,min(x)/2,3,'**','FontSize',20,'Color','red')
            elseif c<0.05
            text(1,min(x)/2,3,'*','FontSize',20,'Color','red')
            else
            text(1,max(x)/2,3,'NS')
        end
    end
end



subplot(1,2,2)
x=ValueI;
y=ValueI_Y;
scatter(ones(numel(x),1),x,'MarkerEdgeColor','k','MarkerFaceColor','m','jitter','on','jitterAmount', .1);
hold on
scatter(2*ones(numel(y),1),y,'MarkerEdgeColor','k','MarkerFaceColor','b','jitter','on','jitterAmount', .1);
xlim([0 3])
ylim([0 1])
alpha(0.5)
title('C Increase')

if length(x) > 1 & length(y) > 1
    [h1 b1]=kstest( (x-mean(x))/std(x) );
    [h2 b2]=kstest( (y-mean(y))/std(y) );
    

    
    
    if h1==0 & h2==0 % data normal
        [c d]=ttest2(x,y);
        if d<0.001
            text(1,max(x)/2,3,'***','FontSize',20,'Color','red')
        elseif d<0.01
            text(1,max(x)/2,3,'**','FontSize',20,'Color','red')
            elseif d<0.05
            text(1,max(x)/2,3,'*','FontSize',20,'Color','red')
        else
            text(1,max(x)/2,3,'NS')
        end
    else
        [c d]=ranksum(x,y);
        if c<0.001
            text(1,max(x)/2,3,'***','FontSize',20,'Color','red')
        elseif c<0.01
            text(1,max(x)/2,3,'**','FontSize',20,'Color','red')
            elseif c<0.05
            text(1,max(x)/2,3,'*','FontSize',20,'Color','red')
            else
            text(1,max(x)/2,3,'NS')
        end
    end
end

sgtitle('Reliability change')
hold off

figure
hold on
for i = 1:length(ValueD)
    plot([0, ValueD(i)],'Color','m')
    alpha(0.5)
end
for i = 1:length(ValueD_Y)
    plot([0, ValueD_Y(i)],'Color','b')
    alpha(0.5)
end
hold on
LINE=[];
LINEY=[];
LINE = [zeros(1,length(ValueD))' ValueD'];
LINEY = [zeros(1,length(ValueD_Y))' ValueD_Y'];
stdshade(LINE,0.5,'m');
hold on
stdshade(LINEY,0.5,'b');

hold on
for i = 1:length(ValueI)
    plot([0, ValueI(i)],'Color','m')
    alpha(0.5)
end
for i = 1:length(ValueI_Y)
    plot([0, ValueI_Y(i)],'Color','b')
    alpha(0.5)
end
hold on
LINE=[];
LINEY=[];
LINE = [zeros(1,length(ValueI))' ValueI'];
LINEY = [zeros(1,length(ValueI_Y))' ValueI_Y'];
stdshade(LINE,0.5,'m');
hold on
stdshade(LINEY,0.5,'b');
title('C')
ylim([-1 1])
% - - - - - - - - - 
%% ---- Ad FIBER
% - - - - - - - - - 


DecM=[];
DecY=[];

incrM=[];
incrY=[];

IfInfWhichValue=4860;

PerOfChangeM=[];
PerOfChangeY=[];
ValueI_Y=[];
ValueI=[];
ValueD_Y=[];
ValueD=[];

figure
sgtitle('Ad')
hold on


if vGlutON_analysis == 1
    PerOfChangeM=(SoP_Ad(2,find(CheckvGlutONAd==1))-SoP_Ad(1,find(CheckvGlutONAd==1)));
    PerOfChangeY=(SoP_Ad_Y(2,find(CheckvGlutONAdY==1))-SoP_Ad_Y(1,find(CheckvGlutONAdY==1)));
elseif vGatON_analysis == 1
    PerOfChangeM=(SoP_Ad(2,find(CheckvGatONAd==1))-SoP_Ad(1,find(CheckvGatONAd==1)));
    PerOfChangeY=(SoP_Ad_Y(2,find(CheckvGatONAdY==1))-SoP_Ad_Y(1,find(CheckvGatONAdY==1)));
elseif Ptf1aON_analysis == 1
    PerOfChangeM=(SoP_Ad(2,find(CheckPtf1aAd==1))-SoP_Ad(1,find(CheckPtf1aAd==1)));
    PerOfChangeY=(SoP_Ad_Y(2,find(CheckPtf1aAdY==1))-SoP_Ad_Y(1,find(CheckPtf1aAdY==1)));
else
    PerOfChangeM=(SoP_Ad(2,1:end)-SoP_Ad(1,1:end));
    PerOfChangeY=(SoP_Ad_Y(2,1:end)-SoP_Ad_Y(1,1:end));
end


a=1;
b=1;
for i=1:numel(PerOfChangeM)
    if isinf(PerOfChangeM(i))
    PerOfChangeM(i)=IfInfWhichValue;
    end
    if PerOfChangeM(i)>threshold %Increase
        subplot(2,3,1)
        line([0 1],[0 (PerOfChangeM(i))],'Color','r')
        ValueI(a)=PerOfChangeM(i);
        axis([0 1 0 400])
        incrM(a)=PerOfChangeM(i);
        a=a+1;
    elseif PerOfChangeM(i)<-threshold %Increase
        subplot(2,3,2)
        line([0 1],[0 -(-PerOfChangeM(i))],'Color','r')
        axis([0 1 -100 0])
        DecM(b)=PerOfChangeM(i);
        ValueD(b)=PerOfChangeM(i);
        b=b+1;
    else
        subplot(2,3,3)
        line([0 1],[0 PerOfChangeM(i)],'Color','r')
        axis([0 1 -10 10])
    end    
end

a=1;
b=1;
for i=1:numel(PerOfChangeY)
    if isinf(PerOfChangeY(i))
    PerOfChangeY(i)=IfInfWhichValue;
    end
    if PerOfChangeY(i)>threshold %Increase
        subplot(2,3,4)
        line([0 1],[0 (PerOfChangeY(i))],'Color','b')
        axis([0 1 0 400])
        incrY(a)=PerOfChangeY(i);
        ValueI_Y(a)=PerOfChangeY(i);
        a=a+1;
    elseif PerOfChangeY(i)<-threshold %Increase
        subplot(2,3,5)
        line([0 1],[0 -(-PerOfChangeY(i))],'Color','b')
        axis([0 1 -100 0])
        DecY(b)=PerOfChangeY(i);
        ValueD_Y(b)=PerOfChangeY(i);
        b=b+1;
    else
        subplot(2,3,6)
        line([0 1],[0 PerOfChangeY(i)],'Color','b')
        axis([0 1 -10 10])
    end
    
end

PerOfChangeY_Ad=PerOfChangeY;
PerOfChangeM_Ad=PerOfChangeM;

%[a b]=ranksum(incrM,incrY);
%[a b]=ranksum(DecM,DecY);


figure
sgtitle(' reliability')
ax1 = nexttile;
pie(ax1,[length(PerOfChangeM_Ad)-numel(DecM)-numel(incrM),numel(DecM),numel(incrM)]);
title('Master Ad')

ax2 = nexttile;
pie(ax2,[length(PerOfChangeY_Ad)-numel(DecY)-numel(incrY),numel(DecY),numel(incrY)]);
title('Yoked Ad')
legend('No change','Decrease','Increase')

% No change
chisquare(size(SoP_Ad,2)-numel(DecM)-numel(incrM), ...
    size(SoP_Ad,2), ...
    size(SoP_Ad_Y,2)-numel(DecY)-numel(incrY), ...
    size(SoP_Ad_Y,2))
% Decrease
chisquare(numel(DecM), ...
    size(SoP_Ad,2), ...
    numel(DecY), ...
    size(SoP_Ad_Y,2))
%increase
chisquare(numel(incrM), ...
    size(SoP_Ad,2), ...
    numel(incrY), ...
    size(SoP_Ad_Y,2))


figure
subplot(1,2,1)
hold on
for i = 1:length(PerOfChangeM_Ad)
    plot([0, PerOfChangeM_Ad(i)],'Color','m')
end
ylim([-1 1])
subplot(1,2,2)
hold on 
for i = 1:length(PerOfChangeY_Ad)
    plot([0, PerOfChangeY_Ad(i)],'Color','b')
end
ylim([-1 1])


figure; hold on
LINE = [zeros(1,length(PerOfChangeM_Ad))' PerOfChangeM_Ad'];
LINEY = [zeros(1,length(PerOfChangeY_Ad))' PerOfChangeY_Ad'];
stdshade(LINE,0.5,'m');
hold on
stdshade(LINEY,0.5,'b');
ylim([-1 1])
title('Ad')

x=ValueD;
y=ValueD_Y;

figure
scatter(ones(numel(PerOfChangeM_Ad),1),PerOfChangeM_Ad,'MarkerEdgeColor','k','MarkerFaceColor','m','jitter','on','jitterAmount', .1);
hold on
scatter(2*ones(numel(PerOfChangeY_Ad),1),PerOfChangeY_Ad,'MarkerEdgeColor','k','MarkerFaceColor','b','jitter','on','jitterAmount', .1);
xlim([0 3])

figure
subplot(1,2,1)
scatter(ones(numel(x),1),x,'MarkerEdgeColor','k','MarkerFaceColor','m','jitter','on','jitterAmount', .1);
hold on
scatter(2*ones(numel(y),1),y,'MarkerEdgeColor','k','MarkerFaceColor','b','jitter','on','jitterAmount', .1);
xlim([0 3])
ylim([-1 0])
% alpha(0.5)
title('Ad Decrease')
if length(x) > 1 & length(y) > 1
    [h1 b1]=kstest( (x-mean(x))/std(x) );
    [h2 b2]=kstest( (y-mean(y))/std(y) );
    

    
    
    if h1==0 & h2==0 % data normal
        [c d]=ttest2(x,y);
        if d<0.001
            text(1.5,min(x)/2,3,'***','FontSize',20,'Color','red')
        elseif d<0.01
            text(1.5,min(x)/2,3,'**','FontSize',20,'Color','red')
            elseif d<0.05
            text(1.5,min(x)/2,3,'*','FontSize',20,'Color','red')
        else
            text(1.5,min(x)/2,3,'NS')
        end
    else
        [c d]=ranksum(x,y);
        if c<0.001
            text(1.5,min(x)/2,3,'***','FontSize',20,'Color','red')
        elseif c<0.01
            text(1.5,min(x)/2,3,'**','FontSize',20,'Color','red')
            elseif c<0.05
            text(1.5,min(x)/2,3,'*','FontSize',20,'Color','red')
            else
            text(1.5,max(x)/2,3,'NS')
        end
    end
end

% violinplot(x,ones(numel(x),1))
% %scatter(ones(numel(x),1),x)
% ylim([-100 0])
% subplot(1,2,2)
% violinplot(y,ones(numel(y),1))
% %scatter(2*ones(numel(y),1),y)
% ylim([-100 0])

x=ValueI;
y=ValueI_Y;

subplot(1,2,2)

scatter(ones(numel(x),1),x,'MarkerEdgeColor','k','MarkerFaceColor','m','jitter','on','jitterAmount', .1);
hold on
scatter(2*ones(numel(y),1),y,'MarkerEdgeColor','k','MarkerFaceColor','b','jitter','on','jitterAmount', .1);
xlim([0 3])
ylim([0 1])
title('Ad Increase')
% alpha(0.5)
if length(x) > 1 & length(y) > 1
    [h1 b1]=kstest( (x-mean(x))/std(x) );
    [h2 b2]=kstest( (y-mean(y))/std(y) );
    
    if h1==0 & h2==0 % data normal
        [c d]=ttest2(x,y);
        if d<0.001
            text(1.5,max(x)/2,3,'***','FontSize',20,'Color','red')
        elseif d<0.01
            text(1.5,max(x)/2,3,'**','FontSize',20,'Color','red')
            elseif d<0.05
            text(1.5,max(x)/2,3,'*','FontSize',20,'Color','red')
        else
            text(1.5,max(x)/2,3,'NS')
        end
    else
        [c d]=ranksum(x,y);
        if c<0.001
            text(1.5,max(x)/2,3,'***','FontSize',20,'Color','red')
        elseif c<0.01
            text(1.5,max(x)/2,3,'**','FontSize',20,'Color','red')
            elseif c<0.05
            text(1.5,max(x)/2,3,'*','FontSize',20,'Color','red')
            else
            text(1.5,max(x)/2,3,'NS')
        end
    end
end

sgtitle('Reliability change')


figure
hold on
for i = 1:length(ValueD)
    plot([0, ValueD(i)],'Color','m')
    alpha(0.5)
end
for i = 1:length(ValueD_Y)
    plot([0, ValueD_Y(i)],'Color','b')
    alpha(0.5)
end
hold on
LINE=[];
LINEY=[];
LINE = [zeros(1,length(ValueD))' ValueD'];
LINEY = [zeros(1,length(ValueD_Y))' ValueD_Y'];
stdshade(LINE,0.5,'m');
hold on
stdshade(LINEY,0.5,'b');

hold on
for i = 1:length(ValueI)
    plot([0, ValueI(i)],'Color','m')
    alpha(0.5)
end
for i = 1:length(ValueI_Y)
    plot([0, ValueI_Y(i)],'Color','b')
    alpha(0.5)
end
hold on
LINE=[];
LINEY=[];
LINE = [zeros(1,length(ValueI))' ValueI'];
LINEY = [zeros(1,length(ValueI_Y))' ValueI_Y'];
stdshade(LINE,0.5,'m');
hold on
stdshade(LINEY,0.5,'b');
title('Ad')
ylim([-1 1])
% - - - - - - - - - 
%% ---- AB FIBER
% - - - - - - - - - 
clc

DecM=[];
DecY=[];

incrM=[];
incrY=[];

PerOfChangeM=[];
PerOfChangeY=[];

ValueI_Y=[];
ValueI=[];
ValueD_Y=[];
ValueD=[];

figure
hold on

if vGlutON_analysis == 1
    PerOfChangeM=(SoP_AB(2,find(CheckvGlutONAB==1))-SoP_AB(1,find(CheckvGlutONAB==1)));
    PerOfChangeY=(SoP_AB_Y(2,find(CheckvGlutONABY==1))-SoP_AB_Y(1,find(CheckvGlutONABY==1)));
elseif vGatON_analysis == 1
    PerOfChangeM=(SoP_AB(2,find(CheckvGatONAB==1))-SoP_AB(1,find(CheckvGatONAB==1)));
    PerOfChangeY=(SoP_AB_Y(2,find(CheckvGatONABY==1))-SoP_AB_Y(1,find(CheckvGatONABY==1)));
elseif Ptf1aON_analysis == 1
    PerOfChangeM=(SoP_AB(2,find(CheckPtf1aAB==1))-SoP_AB(1,find(CheckPtf1aAB==1)));
    PerOfChangeY=(SoP_AB_Y(2,find(CheckPtf1aABY==1))-SoP_AB_Y(1,find(CheckPtf1aABY==1)));
else
    PerOfChangeM=(SoP_AB(2,1:end)-SoP_AB(1,1:end));
    PerOfChangeY=(SoP_AB_Y(2,1:end)-SoP_AB_Y(1,1:end));
end

a=1;
b=1;
for i=1:numel(PerOfChangeM)
    if isinf(PerOfChangeM(i))
    PerOfChangeM(i)=IfInfWhichValue;
    end
    if PerOfChangeM(i)>threshold %Increase
        subplot(2,3,1)
        line([0 1],[0 (PerOfChangeM(i))],'Color','r')
        axis([0 1 0 400])
        incrM(a)=PerOfChangeM(i);
        ValueI(a)=PerOfChangeM(i);
        IncreaseAB(a)=Id_respAB(i);
        a=a+1;
    elseif PerOfChangeM(i)<-threshold %Decrease
        subplot(2,3,2)
        line([0 1],[0 -(-PerOfChangeM(i))],'Color','r')
        axis([0 1 -100 0])
        ValueD(b)=PerOfChangeM(i);
        DecM(b)=PerOfChangeM(i);
        DecreaseAB(b)=Id_respAB(i);
        b=b+1;
    else
        subplot(2,3,3)
        line([0 1],[0 PerOfChangeM(i)],'Color','r')
        axis([0 1 -10 10])
    end
    
end

a=1;
b=1;
for i=1:numel(PerOfChangeY)
    if isinf(PerOfChangeY(i))
    PerOfChangeY(i)=IfInfWhichValue;
    end
    if PerOfChangeY(i)>threshold %Increase
        subplot(2,3,4)
        line([0 1],[0 (PerOfChangeY(i))],'Color','b')
        axis([0 1 0 400])
        ValueI_Y(a)=PerOfChangeY(i);
        incrY(a)=PerOfChangeY(i);
        IncreaseAB_Y(a)=Id_respAB_Y(i);
        a=a+1;
    elseif PerOfChangeY(i)<-threshold %Increase
        subplot(2,3,5)
        line([0 1],[0 -(-PerOfChangeY(i))],'Color','b')
        axis([0 1 -100 0])
        ValueD_Y(b)=PerOfChangeY(i);
        DecY(b)=PerOfChangeY(i);
        DecreaseAB_Y(b)=Id_respAB_Y(i);
        b=b+1;
    else
        subplot(2,3,6)
        line([0 1],[0 PerOfChangeY(i)],'Color','b')
        axis([0 1 -10 10])
    end
    
end
PerOfChangeY_AB=PerOfChangeY;
PerOfChangeM_AB=PerOfChangeM;


% [a b]=ranksum(incrM,incrY);

figure
sgtitle('AB reliability')
ax1 = nexttile;
pie(ax1,[length(PerOfChangeM_AB)-numel(DecM)-numel(incrM),numel(DecM),numel(incrM)]);
title('Master AB')

ax2 = nexttile;
pie(ax2,[length(PerOfChangeY_AB)-numel(DecY)-numel(incrY),numel(DecY),numel(incrY)]);
title('Yoked AB')
legend('No change','Decrease','Increase')

chisquare(size(SoP_AB,2)-numel(DecM)-numel(incrM), ...
    size(SoP_AB,2), ...
    size(SoP_AB_Y,2)-numel(DecY)-numel(incrY), ...
    size(SoP_AB_Y,2))

chisquare(numel(DecM), ...
    size(SoP_AB,2), ...
    numel(DecY), ...
    size(SoP_AB_Y,2))

chisquare(numel(incrM), ...
    size(SoP_AB,2), ...
    numel(incrY), ...
    size(SoP_AB_Y,2))


figure

subplot(1,2,1)

x=ValueD;
y=ValueD_Y;

scatter(ones(numel(x),1),x,'MarkerEdgeColor','k','MarkerFaceColor','m','jitter','on','jitterAmount', .1);
hold on
scatter(2*ones(numel(y),1),y,'MarkerEdgeColor','k','MarkerFaceColor','b','jitter','on','jitterAmount', .1);
xlim([0 3])
ylim([-1 0])
alpha(0.5)
title('Ab Decrease')

if length(x) > 1 & length(y) > 1
    [h1 b1]=kstest( (x-mean(x))/std(x) );
    [h2 b2]=kstest( (y-mean(y))/std(y) );

    if h1==0 & h2==0 % data normal
        [c d]=ttest2(x,y);
        if d<0.001
            text(1,min(x)/2,3,'***','FontSize',20,'Color','red')
        elseif d<0.01
            text(1,min(x)/2,3,'**','FontSize',20,'Color','red')
            elseif d<0.05
            text(1,min(x)/2,3,'*','FontSize',20,'Color','red')
        else
            text(1,min(x)/2,3,'NS')
        end
    else
        [c d]=ranksum(x,y);
        if c<0.001
            text(1,min(x)/2,3,'***','FontSize',20,'Color','red')
        elseif c<0.01
            text(1,min(x)/2,3,'**','FontSize',20,'Color','red')
            elseif c<0.05
            text(1,min(x)/2,3,'*','FontSize',20,'Color','red')
            else
            text(1,max(x)/2,3,'NS')
        end
    end
end


x=ValueI;
y=ValueI_Y;

subplot(1,2,2)


scatter(ones(numel(x),1),x,'MarkerEdgeColor','k','MarkerFaceColor','m','jitter','on','jitterAmount', .1);
hold on
scatter(2*ones(numel(y),1),y,'MarkerEdgeColor','k','MarkerFaceColor','b','jitter','on','jitterAmount', .1);
xlim([0 3])
ylim([0 1])
alpha(0.5)
title('Ab Increase')

if length(x) > 1 & length(y) > 1
    
    [h1 b1]=kstest( (x-mean(x))/std(x) );
    [h2 b2]=kstest( (y-mean(y))/std(y) );
   
    if h1==0 & h2==0 % data normal
        [c d]=ttest2(x,y);
        if d<0.001
            text(1,max(x)/2,3,'***','FontSize',20,'Color','red')
        elseif d<0.01
            text(1,max(x)/2,3,'**','FontSize',20,'Color','red')
            elseif d<0.05
            text(1,max(x)/2,3,'*','FontSize',20,'Color','red')
        else
            text(1,max(x)/2,3,'NS')
        end
    else
        [c d]=ranksum(x,y);
        if c<0.001
            text(1,max(x)/2,3,'***','FontSize',20,'Color','red')
        elseif c<0.01
            text(1,max(x)/2,3,'**','FontSize',20,'Color','red')
            elseif c<0.05
            text(1,max(x)/2,3,'*','FontSize',20,'Color','red')
            else
            text(1,max(x)/2,3,'NS')
        end
    end
end

sgtitle('Reliability change')


figure
hold on
for i = 1:length(ValueD)
    plot([0, ValueD(i)],'Color','m')
    alpha(0.5)
end
for i = 1:length(ValueD_Y)
    plot([0, ValueD_Y(i)],'Color','b')
    alpha(0.5)
end
hold on
LINE=[];
LINEY=[];
LINE = [zeros(1,length(ValueD))' ValueD'];
LINEY = [zeros(1,length(ValueD_Y))' ValueD_Y'];
stdshade(LINE,0.5,'m');
hold on
stdshade(LINEY,0.5,'b');

hold on
for i = 1:length(ValueI)
    plot([0, ValueI(i)],'Color','m')
    alpha(0.5)
end
for i = 1:length(ValueI_Y)
    plot([0, ValueI_Y(i)],'Color','b')
    alpha(0.5)
end
hold on
LINE=[];
LINEY=[];
LINE = [zeros(1,length(ValueI))' ValueI'];
LINEY = [zeros(1,length(ValueI_Y))' ValueI_Y'];
stdshade(LINE,0.5,'m');
hold on
stdshade(LINEY,0.5,'b');
title('Ab')
ylim([-1 1])

% - - - - - - - - -
%% ---- P FIBER ANALYSIS 
% - - - - - - - - - -

% - - - - - - -  - - - 
% - - - - - - -  - - - 

IfInfWhichValue=1500;
% - - - - - - -  - - - 
% - - - - - - -  - - - 

DecM=[];
DecY=[];

incrM=[];
incrY=[];

PerOfChangeM=[];
PerOfChangeY=[];
ValueI_Y=[];
ValueI=[];
ValueD_Y=[];
ValueD=[];

figure
sgtitle('P fibers')
hold on

if vGlutON_analysis == 1
    PerOfChangeM=(SoP_P(2,find(CheckvGlutONP==1))-SoP_P(1,find(CheckvGlutONP==1)));
    PerOfChangeY=(SoP_P_Y(2,find(CheckvGlutONPY==1))-SoP_P_Y(1,find(CheckvGlutONPY==1)));
elseif vGatON_analysis == 1
    PerOfChangeM=(SoP_P(2,find(CheckvGatONP==1))-SoP_P(1,find(CheckvGatONP==1)));
    PerOfChangeY=(SoP_P_Y(2,find(CheckvGatONPY==1))-SoP_P_Y(1,find(CheckvGatONPY==1)));
elseif Ptf1aON_analysis == 1
    PerOfChangeM=(SoP_P(2,find(CheckPtf1aP==1))-SoP_P(1,find(CheckPtf1aP==1)));
    PerOfChangeY=(SoP_P_Y(2,find(CheckPtf1aPY==1))-SoP_P_Y(1,find(CheckPtf1aPY==1)));
else
    PerOfChangeM=(SoP_P(2,1:end)-SoP_P(1,1:end));
    PerOfChangeY=(SoP_P_Y(2,1:end)-SoP_P_Y(1,1:end));
end
a=1;
b=1;


for i=1:numel(PerOfChangeM)
    
    if isinf(PerOfChangeM(i))
    PerOfChangeM(i)=IfInfWhichValue;
    end

    if PerOfChangeM(i)>threshold % Increase
        subplot(2,3,1)
        line([0 1],[0 (PerOfChangeM(i))],'Color','r')
        axis([0 1 0 400])
        incrM(a)=PerOfChangeM(i);
        ValueI(a)=PerOfChangeM(i);
        a=a+1;
    elseif PerOfChangeM(i)<-threshold % Decrease
        subplot(2,3,2)
        line([0 1],[0 -(-PerOfChangeM(i))],'Color','r')
        axis([0 1 -100 0])
        DecM(b)=PerOfChangeM(i);
        ValueD(b)=PerOfChangeM(i);
        b=b+1;
    else % No Change
        subplot(2,3,3)
        line([0 1],[0 PerOfChangeM(i)],'Color','r')
        axis([0 1 -10 10])
    end
    
end

a=1;
b=1;
for i=1:numel(PerOfChangeY)
    if isinf(PerOfChangeY(i))
    PerOfChangeY(i)=IfInfWhichValue;
    end
    if PerOfChangeY(i)>threshold % Increase
        subplot(2,3,4)
        line([0 1],[0 (PerOfChangeY(i))],'Color','b')
        axis([0 1 0 400])
        ValueI_Y(a)=PerOfChangeY(i);
        incrY(a)=PerOfChangeY(i);
        a=a+1;
    elseif PerOfChangeY(i)<-threshold % Decrease
        subplot(2,3,5)
        line([0 1],[0 -(-PerOfChangeY(i))],'Color','b')
        axis([0 1 -100 0])
        ValueD_Y(b)=PerOfChangeY(i);
        DecY(b)=PerOfChangeY(i);
        b=b+1;
    else
        subplot(2,3,6)
        line([0 1],[0 PerOfChangeY(i)],'Color','b')
        axis([0 1 -10 10])
    end
    
end

PerOfChangeY_P=PerOfChangeY;
PerOfChangeM_P=PerOfChangeM;


% [a b]=ranksum(incrM,incrY);

figure
title(' reliability')
ax1 = nexttile;
pie(ax1,[size(SoP_P,2)-numel(DecM)-numel(incrM),numel(DecM),numel(incrM)]);
title('Master P')

ax2 = nexttile;
pie(ax2,[size(SoP_P_Y,2)-numel(DecY)-numel(incrY),numel(DecY),numel(incrY)]);
title('Yoked P')
legend('No change','Decrease','Increase')

chisquare(size(SoP_P,2)-numel(DecM)-numel(incrM), ...
    size(SoP_P,2), ...
    size(SoP_P_Y,2)-numel(DecY)-numel(incrY), ...
    size(SoP_P_Y,2))

chisquare(numel(DecM), ...
    size(SoP_P,2), ...
    numel(DecY), ...
    size(SoP_P_Y,2))

chisquare(numel(incrM), ...
    size(SoP_P,2), ...
    numel(incrY), ...
    size(SoP_P_Y,2))


% Violin Plot
x=ValueD;
y=ValueD_Y;

figure
subplot(1,2,1)

scatter(ones(numel(x),1),x,'MarkerEdgeColor','k','MarkerFaceColor','m','jitter','on','jitterAmount', .1);
hold on
scatter(2*ones(numel(y),1),y,'MarkerEdgeColor','k','MarkerFaceColor','b','jitter','on','jitterAmount', .1);
xlim([0 3])
ylim([-1 0])
alpha(0.5)
title('P Decrease')

if length(x) > 1 & length(y) > 1
    [h1 b1]=kstest( (x-mean(x))/std(x) );
    [h2 b2]=kstest( (y-mean(y))/std(y) );
    
    subplot(1,2,1)
    if h1==0 & h2==0 % data normal
        [c d]=ttest2(x,y);
        if d<0.001
            text(1,min(x)/2,3,'***','FontSize',20,'Color','red')
        elseif d<0.01
            text(1,min(x)/2,3,'**','FontSize',20,'Color','red')
            elseif d<0.05
            text(1,min(x)/2,3,'*','FontSize',20,'Color','red')
        else
            text(1,min(x)/2,3,'NS')
        end
    else
        [c d]=ranksum(x,y);
        if c<0.001
            text(1,min(x)/2,3,'***','FontSize',20,'Color','red')
        elseif c<0.01
            text(1,min(x)/2,3,'**','FontSize',20,'Color','red')
            elseif c<0.05
            text(1,min(x)/2,3,'*','FontSize',20,'Color','red')
            else
            text(1,max(x)/2,3,'NS')
        end
    end
end

x=ValueI;
y=ValueI_Y;

subplot(1,2,2)

scatter(ones(numel(x),1),x,'MarkerEdgeColor','k','MarkerFaceColor','m','jitter','on','jitterAmount', .1);
hold on
scatter(2*ones(numel(y),1),y,'MarkerEdgeColor','k','MarkerFaceColor','b','jitter','on','jitterAmount', .1);
xlim([0 3])
ylim([0 1])
alpha(0.5)
title('P Increase')

if length(x) > 1 & length(y) > 1
    [h1 b1]=kstest( (x-mean(x))/std(x) );
    [h2 b2]=kstest( (y-mean(y))/std(y) );
    
    if h1==0 & h2==0 % data normal
        [c d]=ttest2(x,y);
        if d<0.001
            text(1,max(x)/2,3,'***','FontSize',20,'Color','red')
        elseif d<0.01
            text(1,max(x)/2,3,'**','FontSize',20,'Color','red')
            elseif d<0.05
            text(1,max(x)/2,3,'*','FontSize',20,'Color','red')
        else
            text(1,max(x)/2,3,'NS')
        end
    else
        [c d]=ranksum(x,y);
        if c<0.001
            text(1,max(x)/2,3,'***','FontSize',20,'Color','red')
        elseif c<0.01
            text(1,max(x)/2,3,'**','FontSize',20,'Color','red')
            elseif c<0.05
            text(1,max(x)/2,3,'*','FontSize',20,'Color','red')
            else
            text(1,max(x)/2,3,'NS')
        end
    end
end

sgtitle('Reliability change')



figure
hold on
for i = 1:length(ValueD)
    plot([0, ValueD(i)],'Color','m')
    alpha(0.5)
end
for i = 1:length(ValueD_Y)
    plot([0, ValueD_Y(i)],'Color','b')
    alpha(0.5)
end
hold on
LINE=[];
LINEY=[];
LINE = [zeros(1,length(ValueD))' ValueD'];
LINEY = [zeros(1,length(ValueD_Y))' ValueD_Y'];
stdshade(LINE,0.5,'m');
hold on
stdshade(LINEY,0.5,'b');

hold on
for i = 1:length(ValueI)
    plot([0, ValueI(i)],'Color','m')
    alpha(0.5)
end
for i = 1:length(ValueI_Y)
    plot([0, ValueI_Y(i)],'Color','b')
    alpha(0.5)
end
hold on
LINE=[];
LINEY=[];
LINE = [zeros(1,length(ValueI))' ValueI'];
LINEY = [zeros(1,length(ValueI_Y))' ValueI_Y'];
stdshade(LINE,0.5,'m');
hold on
stdshade(LINEY,0.5,'b');
title('P')
ylim([-1 1])
% - - - - - - - - - 
%% Phase tuning analysis
% - - - - - - - - -


a=1;
for i=1:size(frM_UOI,1)
    

    if numel(frM_UOI.latency_all{i})==1
        
        frM_UOI.FiberType{i}=[CheckFiber(a)];
        a=a+1;
        
    else
        frM_UOI.FiberType{i}=[CheckFiber(a:a+numel(frM_UOI.latency_all{i})-1)];
        a=a+numel(frM_UOI.latency_all{i});

    end

end


a=1;
for i=1:size(frY_UOI,1)
    

    if numel(frY_UOI.latency_all{i})==1
        
        frY_UOI.FiberType{i}=[CheckFiberY(a)];
        a=a+1;
        
    else
        frY_UOI.FiberType{i}=[CheckFiberY(a:a+numel(frY_UOI.latency_all{i})-1)];
        a=a+numel(frY_UOI.latency_all{i});

    end

end



%
a=1;

c=1;
ad=1;
ab=1;

for i=1:size(frM_UOI,1)
    
    tabR=[];
    
    for j=1:numel(frM_UOI.latency_all{i})
        
        if CheckFiber(a)==1
            
            if PerOfChangeM_C(c)>threshold
                Rel=1;
            elseif PerOfChangeM_C(c)<-threshold
                Rel=-1;
            else
                Rel=0;
            end
            tabR=[tabR Rel]
            c=c+1;
            
        elseif CheckFiber(a)==2
            
            if PerOfChangeM_Ad(ad)>threshold
                Rel=1;
            elseif PerOfChangeM_Ad(ad)<-threshold
                Rel=-1;
            else
                Rel=0;
            end

            tabR=[tabR Rel]
            ad=ad+1;
            
        elseif CheckFiber(a)==3
            
            if PerOfChangeM_AB(ab)>threshold
                Rel=1;
            elseif PerOfChangeM_AB(ab)<-threshold
                Rel=-1;
            else
                Rel=0;
            end

            tabR=[tabR Rel]
            ab=ab+1;

        elseif CheckFiber(a)==4
                Rel=4;
                tabR=[tabR Rel]
        end
    
        a=a+1; 
    end
    
    %numel(frM_UOI.latency_all{i})
    frM_UOI.ReliabilityType{i}=tabR;
    
end






a=1;

c=1;
ad=1;
ab=1;

for i=1:size(frY_UOI,1)
    
    tabR=[];
    
    for j=1:numel(frY_UOI.latency_all{i})
        
        if CheckFiberY(a)==1
            
            if PerOfChangeY_C(c)>threshold
                Rel=1;
            elseif PerOfChangeY_C(c)<-threshold
                Rel=-1;
            else
                Rel=0;
            end
            tabR=[tabR Rel]
            c=c+1;
            
        elseif CheckFiberY(a)==2
            
            if PerOfChangeY_Ad(ad)>threshold
                Rel=1;
            elseif PerOfChangeY_Ad(ad)<-threshold
                Rel=-1;
            else
                Rel=0;
            end

            tabR=[tabR Rel]
            ad=ad+1;
            
        elseif CheckFiberY(a)==3
            
            if PerOfChangeY_AB(ab)>threshold
                Rel=1;
            elseif PerOfChangeY_AB(ab)<-threshold
                Rel=-1;
            else
                Rel=0;
            end

            tabR=[tabR Rel]
            ab=ab+1;

        elseif CheckFiberY(a)==4
                Rel=4;
                tabR=[tabR Rel]
        end
    
        a=a+1; 
    end
    
    %numel(frM_UOI.latency_all{i})
    frY_UOI.ReliabilityType{i}=tabR;
    
end



% - - - - - - - - - - -

tabReliabCurated_M=[];
tabPhaseTuning_M=[];

for i=1:size(frM_UOI,1)

    tabReliabCurated_M=[tabReliabCurated_M frM_UOI.ReliabilityType{i}];
    tabPhaseTuning_M=[tabPhaseTuning_M transpose(ones(numel(frM_UOI.ReliabilityType{i}),1)*frM_UOI.Category(i))];
end


tabReliabCurated_Y=[];
tabPhaseTuning_Y=[];

for i=1:size(frY_UOI,1)

    tabReliabCurated_Y=[tabReliabCurated_Y frY_UOI.ReliabilityType{i}];
    tabPhaseTuning_Y=[tabPhaseTuning_Y transpose(ones(numel(frY_UOI.ReliabilityType{i}),1)*frY_UOI.Category(i))];
end


% %  - - - - - - - - - -MASTER

IncreaseM_C=find(tabReliabCurated_M==1 & CheckFiber==1);

IncreaseM_C_E=numel(find(tabPhaseTuning_M(IncreaseM_C)==3 | tabPhaseTuning_M(IncreaseM_C)==4));
IncreaseM_C_S=numel(find(tabPhaseTuning_M(IncreaseM_C)==1 | tabPhaseTuning_M(IncreaseM_C)==2));
IncreaseM_C_L=numel(find(tabPhaseTuning_M(IncreaseM_C)==5 | tabPhaseTuning_M(IncreaseM_C)==6));
IncreaseM_C_D=numel(find(tabPhaseTuning_M(IncreaseM_C)==7 | tabPhaseTuning_M(IncreaseM_C)==8 | tabPhaseTuning_M(IncreaseM_C)==9));
IncreaseM_C_NC=numel(find(tabPhaseTuning_M(IncreaseM_C)==10 | isnan(tabPhaseTuning_M(IncreaseM_C)) ));


DecreaseM_C=find(tabReliabCurated_M==-1 & CheckFiber==1);

DecreaseM_C_E=numel(find(tabPhaseTuning_M(DecreaseM_C)==3 | tabPhaseTuning_M(DecreaseM_C)==4));
DecreaseM_C_S=numel(find(tabPhaseTuning_M(DecreaseM_C)==1 | tabPhaseTuning_M(DecreaseM_C)==2));
DecreaseM_C_L=numel(find(tabPhaseTuning_M(DecreaseM_C)==5 | tabPhaseTuning_M(DecreaseM_C)==6));
DecreaseM_C_D=numel(find(tabPhaseTuning_M(DecreaseM_C)==7 | tabPhaseTuning_M(DecreaseM_C)==8 | tabPhaseTuning_M(DecreaseM_C)==9));
DecreaseM_C_NC=numel(find(tabPhaseTuning_M(DecreaseM_C)==10 | isnan(tabPhaseTuning_M(DecreaseM_C)) ));


NoChangeM_C=find(tabReliabCurated_M==0 & CheckFiber==1);

NoChangeM_C_E=numel(find(tabPhaseTuning_M(NoChangeM_C)==3 | tabPhaseTuning_M(NoChangeM_C)==4));
NoChangeM_C_S=numel(find(tabPhaseTuning_M(NoChangeM_C)==1 | tabPhaseTuning_M(NoChangeM_C)==2));
NoChangeM_C_L=numel(find(tabPhaseTuning_M(NoChangeM_C)==5 | tabPhaseTuning_M(NoChangeM_C)==6));
NoChangeM_C_D=numel(find(tabPhaseTuning_M(NoChangeM_C)==7 | tabPhaseTuning_M(NoChangeM_C)==8 | tabPhaseTuning_M(NoChangeM_C)==9));
NoChangeM_C_NC=numel(find(tabPhaseTuning_M(NoChangeM_C)==10 | isnan(tabPhaseTuning_M(NoChangeM_C)) ));


% % % % % % % % % % Ad

IncreaseM_Ad=find(tabReliabCurated_M==1 & CheckFiber==2);

IncreaseM_Ad_E=numel(find(tabPhaseTuning_M(IncreaseM_Ad)==3 | tabPhaseTuning_M(IncreaseM_Ad)==4));
IncreaseM_Ad_S=numel(find(tabPhaseTuning_M(IncreaseM_Ad)==1 | tabPhaseTuning_M(IncreaseM_Ad)==2));
IncreaseM_Ad_L=numel(find(tabPhaseTuning_M(IncreaseM_Ad)==5 | tabPhaseTuning_M(IncreaseM_Ad)==6));
IncreaseM_Ad_D=numel(find(tabPhaseTuning_M(IncreaseM_Ad)==7 | tabPhaseTuning_M(IncreaseM_Ad)==8 | tabPhaseTuning_M(IncreaseM_Ad)==9));
IncreaseM_Ad_NC=numel(find(tabPhaseTuning_M(IncreaseM_Ad)==10 | isnan(tabPhaseTuning_M(IncreaseM_Ad)) ));


DecreaseM_Ad=find(tabReliabCurated_M==-1 & CheckFiber==2);

DecreaseM_Ad_E=numel(find(tabPhaseTuning_M(DecreaseM_Ad)==3 | tabPhaseTuning_M(DecreaseM_Ad)==4));
DecreaseM_Ad_S=numel(find(tabPhaseTuning_M(DecreaseM_Ad)==1 | tabPhaseTuning_M(DecreaseM_Ad)==2));
DecreaseM_Ad_L=numel(find(tabPhaseTuning_M(DecreaseM_Ad)==5 | tabPhaseTuning_M(DecreaseM_Ad)==6));
DecreaseM_Ad_D=numel(find(tabPhaseTuning_M(DecreaseM_Ad)==7 | tabPhaseTuning_M(DecreaseM_Ad)==8 | tabPhaseTuning_M(DecreaseM_Ad)==9));
DecreaseM_Ad_NC=numel(find(tabPhaseTuning_M(DecreaseM_Ad)==10 | isnan(tabPhaseTuning_M(DecreaseM_Ad)) ));


NoChangeM_Ad=find(tabReliabCurated_M==0 & CheckFiber==2);

NoChangeM_Ad_E=numel(find(tabPhaseTuning_M(NoChangeM_Ad)==3 | tabPhaseTuning_M(NoChangeM_Ad)==4));
NoChangeM_Ad_S=numel(find(tabPhaseTuning_M(NoChangeM_Ad)==1 | tabPhaseTuning_M(NoChangeM_Ad)==2));
NoChangeM_Ad_L=numel(find(tabPhaseTuning_M(NoChangeM_Ad)==5 | tabPhaseTuning_M(NoChangeM_Ad)==6));
NoChangeM_Ad_D=numel(find(tabPhaseTuning_M(NoChangeM_Ad)==7 | tabPhaseTuning_M(NoChangeM_Ad)==8 | tabPhaseTuning_M(NoChangeM_Ad)==9));
NoChangeM_Ad_NC=numel(find(tabPhaseTuning_M(NoChangeM_Ad)==10 | isnan(tabPhaseTuning_M(NoChangeM_Ad)) ));


% % % % % % % % % % AB


IncreaseM_AB=find(tabReliabCurated_M==1 & CheckFiber==3);

IncreaseM_AB_E=numel(find(tabPhaseTuning_M(IncreaseM_AB)==3 | tabPhaseTuning_M(IncreaseM_AB)==4));
IncreaseM_AB_S=numel(find(tabPhaseTuning_M(IncreaseM_AB)==1 | tabPhaseTuning_M(IncreaseM_AB)==2));
IncreaseM_AB_L=numel(find(tabPhaseTuning_M(IncreaseM_AB)==5 | tabPhaseTuning_M(IncreaseM_AB)==6));
IncreaseM_AB_D=numel(find(tabPhaseTuning_M(IncreaseM_AB)==7 | tabPhaseTuning_M(IncreaseM_AB)==8 | tabPhaseTuning_M(IncreaseM_AB)==9));
IncreaseM_AB_NC=numel(find(tabPhaseTuning_M(IncreaseM_AB)==10 | isnan(tabPhaseTuning_M(IncreaseM_AB)) ));


DecreaseM_AB=find(tabReliabCurated_M==-1 & CheckFiber==3);

DecreaseM_AB_E=numel(find(tabPhaseTuning_M(DecreaseM_AB)==3 | tabPhaseTuning_M(DecreaseM_AB)==4));
DecreaseM_AB_S=numel(find(tabPhaseTuning_M(DecreaseM_AB)==1 | tabPhaseTuning_M(DecreaseM_AB)==2));
DecreaseM_AB_L=numel(find(tabPhaseTuning_M(DecreaseM_AB)==5 | tabPhaseTuning_M(DecreaseM_AB)==6));
DecreaseM_AB_D=numel(find(tabPhaseTuning_M(DecreaseM_AB)==7 | tabPhaseTuning_M(DecreaseM_AB)==8 | tabPhaseTuning_M(DecreaseM_AB)==9));
DecreaseM_AB_NC=numel(find(tabPhaseTuning_M(DecreaseM_AB)==10 | isnan(tabPhaseTuning_M(DecreaseM_AB)) ));


NoChangeM_AB=find(tabReliabCurated_M==0 & CheckFiber==3);

NoChangeM_AB_E=numel(find(tabPhaseTuning_M(NoChangeM_AB)==3 | tabPhaseTuning_M(NoChangeM_AB)==4));
NoChangeM_AB_S=numel(find(tabPhaseTuning_M(NoChangeM_AB)==1 | tabPhaseTuning_M(NoChangeM_AB)==2));
NoChangeM_AB_L=numel(find(tabPhaseTuning_M(NoChangeM_AB)==5 | tabPhaseTuning_M(NoChangeM_AB)==6));
NoChangeM_AB_D=numel(find(tabPhaseTuning_M(NoChangeM_AB)==7 | tabPhaseTuning_M(NoChangeM_AB)==8 | tabPhaseTuning_M(NoChangeM_AB)==9));
NoChangeM_AB_NC=numel(find(tabPhaseTuning_M(NoChangeM_AB)==10 | isnan(tabPhaseTuning_M(NoChangeM_AB)) ));



% - - - - - -  YOKED

IncreaseY_C=find(tabReliabCurated_Y==1 & CheckFiberY==1);

IncreaseY_C_E=numel(find(tabPhaseTuning_Y(IncreaseY_C)==3 | tabPhaseTuning_Y(IncreaseY_C)==4));
IncreaseY_C_S=numel(find(tabPhaseTuning_Y(IncreaseY_C)==1 | tabPhaseTuning_Y(IncreaseY_C)==2));
IncreaseY_C_L=numel(find(tabPhaseTuning_Y(IncreaseY_C)==5 | tabPhaseTuning_Y(IncreaseY_C)==6));
IncreaseY_C_D=numel(find(tabPhaseTuning_Y(IncreaseY_C)==7 | tabPhaseTuning_Y(IncreaseY_C)==8 | tabPhaseTuning_Y(IncreaseY_C)==9));
IncreaseY_C_NC=numel(find(tabPhaseTuning_Y(IncreaseY_C)==10 | isnan(tabPhaseTuning_Y(IncreaseY_C)) ));


DecreaseY_C=find(tabReliabCurated_Y==-1 & CheckFiberY==1);

DecreaseY_C_E=numel(find(tabPhaseTuning_Y(DecreaseY_C)==3 | tabPhaseTuning_Y(DecreaseY_C)==4));
DecreaseY_C_S=numel(find(tabPhaseTuning_Y(DecreaseY_C)==1 | tabPhaseTuning_Y(DecreaseY_C)==2));
DecreaseY_C_L=numel(find(tabPhaseTuning_Y(DecreaseY_C)==5 | tabPhaseTuning_Y(DecreaseY_C)==6));
DecreaseY_C_D=numel(find(tabPhaseTuning_Y(DecreaseY_C)==7 | tabPhaseTuning_Y(DecreaseY_C)==8 | tabPhaseTuning_Y(DecreaseY_C)==9));
DecreaseY_C_NC=numel(find(tabPhaseTuning_Y(DecreaseY_C)==10 | isnan(tabPhaseTuning_Y(DecreaseY_C)) ));


NoChangeY_C=find(tabReliabCurated_Y==0 & CheckFiberY==1);

NoChangeY_C_E=numel(find(tabPhaseTuning_Y(NoChangeY_C)==3 | tabPhaseTuning_Y(NoChangeY_C)==4));
NoChangeY_C_S=numel(find(tabPhaseTuning_Y(NoChangeY_C)==1 | tabPhaseTuning_Y(NoChangeY_C)==2));
NoChangeY_C_L=numel(find(tabPhaseTuning_Y(NoChangeY_C)==5 | tabPhaseTuning_Y(NoChangeY_C)==6));
NoChangeY_C_D=numel(find(tabPhaseTuning_Y(NoChangeY_C)==7 | tabPhaseTuning_Y(NoChangeY_C)==8 | tabPhaseTuning_Y(NoChangeY_C)==9));
NoChangeY_C_NC=numel(find(tabPhaseTuning_Y(NoChangeY_C)==10 | isnan(tabPhaseTuning_Y(NoChangeY_C)) ));


% % % % % % % % % % Ad

IncreaseY_Ad=find(tabReliabCurated_Y==1 & CheckFiberY==2);

IncreaseY_Ad_E=numel(find(tabPhaseTuning_Y(IncreaseY_Ad)==3 | tabPhaseTuning_Y(IncreaseY_Ad)==4));
IncreaseY_Ad_S=numel(find(tabPhaseTuning_Y(IncreaseY_Ad)==1 | tabPhaseTuning_Y(IncreaseY_Ad)==2));
IncreaseY_Ad_L=numel(find(tabPhaseTuning_Y(IncreaseY_Ad)==5 | tabPhaseTuning_Y(IncreaseY_Ad)==6));
IncreaseY_Ad_D=numel(find(tabPhaseTuning_Y(IncreaseY_Ad)==7 | tabPhaseTuning_Y(IncreaseY_Ad)==8 | tabPhaseTuning_Y(IncreaseY_Ad)==9));
IncreaseY_Ad_NC=numel(find(tabPhaseTuning_Y(IncreaseY_Ad)==10 | isnan(tabPhaseTuning_Y(IncreaseY_Ad)) ));


DecreaseY_Ad=find(tabReliabCurated_Y==-1 & CheckFiberY==2);

DecreaseY_Ad_E=numel(find(tabPhaseTuning_Y(DecreaseY_Ad)==3 | tabPhaseTuning_Y(DecreaseY_Ad)==4));
DecreaseY_Ad_S=numel(find(tabPhaseTuning_Y(DecreaseY_Ad)==1 | tabPhaseTuning_Y(DecreaseY_Ad)==2));
DecreaseY_Ad_L=numel(find(tabPhaseTuning_Y(DecreaseY_Ad)==5 | tabPhaseTuning_Y(DecreaseY_Ad)==6));
DecreaseY_Ad_D=numel(find(tabPhaseTuning_Y(DecreaseY_Ad)==7 | tabPhaseTuning_Y(DecreaseY_Ad)==8 | tabPhaseTuning_Y(DecreaseY_Ad)==9));
DecreaseY_Ad_NC=numel(find(tabPhaseTuning_Y(DecreaseY_Ad)==10 | isnan(tabPhaseTuning_Y(DecreaseY_Ad)) ));


NoChangeY_Ad=find(tabReliabCurated_Y==0 & CheckFiberY==2);

NoChangeY_Ad_E=numel(find(tabPhaseTuning_Y(NoChangeY_Ad)==3 | tabPhaseTuning_Y(NoChangeY_Ad)==4));
NoChangeY_Ad_S=numel(find(tabPhaseTuning_Y(NoChangeY_Ad)==1 | tabPhaseTuning_Y(NoChangeY_Ad)==2));
NoChangeY_Ad_L=numel(find(tabPhaseTuning_Y(NoChangeY_Ad)==5 | tabPhaseTuning_Y(NoChangeY_Ad)==6));
NoChangeY_Ad_D=numel(find(tabPhaseTuning_Y(NoChangeY_Ad)==7 | tabPhaseTuning_Y(NoChangeY_Ad)==8 | tabPhaseTuning_Y(NoChangeY_Ad)==9));
NoChangeY_Ad_NC=numel(find(tabPhaseTuning_Y(NoChangeY_Ad)==10 | isnan(tabPhaseTuning_Y(NoChangeY_Ad)) ));


% % % % % % % % % % AB


IncreaseY_AB=find(tabReliabCurated_Y==1 & CheckFiberY==3);

IncreaseY_AB_E=numel(find(tabPhaseTuning_Y(IncreaseY_AB)==3 | tabPhaseTuning_Y(IncreaseY_AB)==4));
IncreaseY_AB_S=numel(find(tabPhaseTuning_Y(IncreaseY_AB)==1 | tabPhaseTuning_Y(IncreaseY_AB)==2));
IncreaseY_AB_L=numel(find(tabPhaseTuning_Y(IncreaseY_AB)==5 | tabPhaseTuning_Y(IncreaseY_AB)==6));
IncreaseY_AB_D=numel(find(tabPhaseTuning_Y(IncreaseY_AB)==7 | tabPhaseTuning_Y(IncreaseY_AB)==8 | tabPhaseTuning_Y(IncreaseY_AB)==9));
IncreaseY_AB_NC=numel(find(tabPhaseTuning_Y(IncreaseY_AB)==10 | isnan(tabPhaseTuning_Y(IncreaseY_AB)) ));


DecreaseY_AB=find(tabReliabCurated_Y==-1 & CheckFiberY==3);

DecreaseY_AB_E=numel(find(tabPhaseTuning_Y(DecreaseY_AB)==3 | tabPhaseTuning_Y(DecreaseY_AB)==4));
DecreaseY_AB_S=numel(find(tabPhaseTuning_Y(DecreaseY_AB)==1 | tabPhaseTuning_Y(DecreaseY_AB)==2));
DecreaseY_AB_L=numel(find(tabPhaseTuning_Y(DecreaseY_AB)==5 | tabPhaseTuning_Y(DecreaseY_AB)==6));
DecreaseY_AB_D=numel(find(tabPhaseTuning_Y(DecreaseY_AB)==7 | tabPhaseTuning_Y(DecreaseY_AB)==8 | tabPhaseTuning_Y(DecreaseY_AB)==9));
DecreaseY_AB_NC=numel(find(tabPhaseTuning_Y(DecreaseY_AB)==10 | isnan(tabPhaseTuning_Y(DecreaseY_AB)) ));


NoChangeY_AB=find(tabReliabCurated_Y==0 & CheckFiberY==3);

NoChangeY_AB_E=numel(find(tabPhaseTuning_Y(NoChangeY_AB)==3 | tabPhaseTuning_Y(NoChangeY_AB)==4));
NoChangeY_AB_S=numel(find(tabPhaseTuning_Y(NoChangeY_AB)==1 | tabPhaseTuning_Y(NoChangeY_AB)==2));
NoChangeY_AB_L=numel(find(tabPhaseTuning_Y(NoChangeY_AB)==5 | tabPhaseTuning_Y(NoChangeY_AB)==6));
NoChangeY_AB_D=numel(find(tabPhaseTuning_Y(NoChangeY_AB)==7 | tabPhaseTuning_Y(NoChangeY_AB)==8 | tabPhaseTuning_Y(NoChangeY_AB)==9));
NoChangeY_AB_NC=numel(find(tabPhaseTuning_Y(NoChangeY_AB)==10 | isnan(tabPhaseTuning_Y(NoChangeY_AB)) ));


%  - - - - 
close all

% Cfibers
figure
subplot(3,2,1)
pie([IncreaseM_C_E IncreaseM_C_S IncreaseM_C_L IncreaseM_C_D IncreaseM_C_NC])
title(['Master Increase:'  ...
    num2str(IncreaseM_C_E+IncreaseM_C_S+IncreaseM_C_L+IncreaseM_C_D+IncreaseM_C_NC) ' SUs'])
subplot(3,2,2)
pie([IncreaseY_C_E IncreaseY_C_S IncreaseY_C_L IncreaseY_C_D IncreaseY_C_NC])
title(['Yoked Increase:'  ...
    num2str(IncreaseY_C_E+IncreaseY_C_S+IncreaseY_C_L+IncreaseY_C_D+IncreaseY_C_NC) ' SUs'])

TotM=IncreaseM_C_E+IncreaseM_C_S+IncreaseM_C_L+IncreaseM_C_D+IncreaseM_C_NC;
TotY=IncreaseY_C_E+IncreaseY_C_S+IncreaseY_C_L+IncreaseY_C_D+IncreaseY_C_NC;
stat_In_C(1)=chisquare(IncreaseM_C_E,TotM,IncreaseY_C_E,TotY);
stat_In_C(2)=chisquare(IncreaseM_C_S,TotM,IncreaseY_C_S,TotY);
stat_In_C(3)=chisquare(IncreaseM_C_L,TotM,IncreaseY_C_L,TotY);
stat_In_C(4)=chisquare(IncreaseM_C_D,TotM,IncreaseY_C_D,TotY);
stat_In_C(5)=chisquare(IncreaseM_C_NC,TotM,IncreaseY_C_NC,TotY);


subplot(3,2,3)
pie([DecreaseM_C_E DecreaseM_C_S DecreaseM_C_L DecreaseM_C_D DecreaseM_C_NC])
title(['Master Decrease:'  ...
    num2str(DecreaseM_C_E+DecreaseM_C_S+DecreaseM_C_L+DecreaseM_C_D+DecreaseM_C_NC) ' SUs'])

subplot(3,2,4)
pie([DecreaseY_C_E DecreaseY_C_S DecreaseY_C_L DecreaseY_C_D DecreaseY_C_NC])
title(['Yoked Decrease:'  ...
    num2str(DecreaseY_C_E+DecreaseY_C_S+DecreaseY_C_L+DecreaseY_C_D+DecreaseY_C_NC) ' SUs'])

TotM=DecreaseM_C_E+DecreaseM_C_S+DecreaseM_C_L+DecreaseM_C_D+DecreaseM_C_NC;
TotY=DecreaseY_C_E+DecreaseY_C_S+DecreaseY_C_L+DecreaseY_C_D+DecreaseY_C_NC;
stat_Dec_C(1)=chisquare(DecreaseM_C_E,TotM,DecreaseY_C_E,TotY);
stat_Dec_C(2)=chisquare(DecreaseM_C_S,TotM,DecreaseY_C_S,TotY);
stat_Dec_C(3)=chisquare(DecreaseM_C_L,TotM,DecreaseY_C_L,TotY);
stat_Dec_C(4)=chisquare(DecreaseM_C_D,TotM,DecreaseY_C_D,TotY);
stat_Dec_C(5)=chisquare(DecreaseM_C_NC,TotM,DecreaseY_C_NC,TotY);

subplot(3,2,5)
pie([NoChangeM_C_E NoChangeM_C_S NoChangeM_C_L NoChangeM_C_D NoChangeM_C_NC])
title(['Master NC:'  ...
    num2str(NoChangeM_C_E+NoChangeM_C_S+NoChangeM_C_L+NoChangeM_C_D+NoChangeM_C_NC) ' SUs'])

subplot(3,2,6)
pie([NoChangeY_C_E NoChangeY_C_S NoChangeY_C_L NoChangeY_C_D NoChangeY_C_NC])
title(['Yoked NC:'  ...
    num2str(NoChangeY_C_E+NoChangeY_C_S+NoChangeY_C_L+NoChangeY_C_D+NoChangeY_C_NC) ' SUs'])

TotM=NoChangeM_C_E+NoChangeM_C_S+NoChangeM_C_L+NoChangeM_C_D+NoChangeM_C_NC;
TotY=NoChangeY_C_E+NoChangeY_C_S+NoChangeY_C_L+NoChangeY_C_D+NoChangeY_C_NC;
stat_NC_C(1)=chisquare(NoChangeM_C_E,TotM,NoChangeY_C_E,TotY);
stat_NC_C(2)=chisquare(NoChangeM_C_S,TotM,NoChangeY_C_S,TotY);
stat_NC_C(3)=chisquare(NoChangeM_C_L,TotM,NoChangeY_C_L,TotY);
stat_NC_C(4)=chisquare(NoChangeM_C_D,TotM,NoChangeY_C_D,TotY);
stat_NC_C(5)=chisquare(NoChangeM_C_NC,TotM,NoChangeY_C_NC,TotY);


legend('Early','Steady','Late','Decrease','NC')

% Ad fibers
figure
subplot(3,2,1)
pie([IncreaseM_Ad_E IncreaseM_Ad_S IncreaseM_Ad_L IncreaseM_Ad_D IncreaseM_Ad_NC])
title(['Master Increase:'  ...
    num2str(IncreaseM_Ad_E+IncreaseM_Ad_S+IncreaseM_Ad_L+IncreaseM_Ad_D+IncreaseM_Ad_NC) ' SUs'])

subplot(3,2,2)
pie([IncreaseY_Ad_E IncreaseY_Ad_S IncreaseY_Ad_L IncreaseY_Ad_D IncreaseY_Ad_NC])
title(['Yoked Increase:'  ...
    num2str(IncreaseY_Ad_E+IncreaseY_Ad_S+IncreaseY_Ad_L+IncreaseY_Ad_D+IncreaseY_Ad_NC) ' SUs'])

TotM=IncreaseM_Ad_E+IncreaseM_Ad_S+IncreaseM_Ad_L+IncreaseM_Ad_D+IncreaseM_Ad_NC;
TotY=IncreaseY_Ad_E+IncreaseY_Ad_S+IncreaseY_Ad_L+IncreaseY_Ad_D+IncreaseY_Ad_NC;
stat_In_Ad(1)=chisquare(IncreaseM_Ad_E,TotM,IncreaseY_Ad_E,TotY);
stat_In_Ad(2)=chisquare(IncreaseM_Ad_S,TotM,IncreaseY_Ad_S,TotY);
stat_In_Ad(3)=chisquare(IncreaseM_Ad_L,TotM,IncreaseY_Ad_L,TotY);
stat_In_Ad(4)=chisquare(IncreaseM_Ad_D,TotM,IncreaseY_Ad_D,TotY);
stat_In_Ad(5)=chisquare(IncreaseM_Ad_NC,TotM,IncreaseY_Ad_NC,TotY);

subplot(3,2,3)
pie([DecreaseM_Ad_E DecreaseM_Ad_S DecreaseM_Ad_L DecreaseM_Ad_D DecreaseM_Ad_NC])
title(['Master Decrease:'  ...
    num2str(DecreaseM_Ad_E+DecreaseM_Ad_S+DecreaseM_Ad_L+DecreaseM_Ad_D+DecreaseM_Ad_NC) ' SUs'])

subplot(3,2,4)
pie([DecreaseY_Ad_E DecreaseY_Ad_S DecreaseY_Ad_L DecreaseY_Ad_D DecreaseY_Ad_NC])
title(['Yoked Decrease:'  ...
    num2str(DecreaseY_Ad_E+DecreaseY_Ad_S+DecreaseY_Ad_L+DecreaseY_Ad_D+DecreaseY_Ad_NC) ' SUs'])

TotM=DecreaseM_Ad_E+DecreaseM_Ad_S+DecreaseM_Ad_L+DecreaseM_Ad_D+DecreaseM_Ad_NC;
TotY=DecreaseY_Ad_E+DecreaseY_Ad_S+DecreaseY_Ad_L+DecreaseY_Ad_D+DecreaseY_Ad_NC;
stat_Dec_Ad(1)=chisquare(DecreaseM_Ad_E,TotM,DecreaseY_Ad_E,TotY);
stat_Dec_Ad(2)=chisquare(DecreaseM_Ad_S,TotM,DecreaseY_Ad_S,TotY);
stat_Dec_Ad(3)=chisquare(DecreaseM_Ad_L,TotM,DecreaseY_Ad_L,TotY);
stat_Dec_Ad(4)=chisquare(DecreaseM_Ad_D,TotM,DecreaseY_Ad_D,TotY);
stat_Dec_Ad(5)=chisquare(DecreaseM_Ad_NC,TotM,DecreaseY_Ad_NC,TotY);

subplot(3,2,5)
pie([NoChangeM_Ad_E NoChangeM_Ad_S NoChangeM_Ad_L NoChangeM_Ad_D NoChangeM_Ad_NC])
title(['Master NC:'  ...
    num2str(NoChangeM_Ad_E+NoChangeM_Ad_S+NoChangeM_Ad_L+NoChangeM_Ad_D+NoChangeM_Ad_NC) ' SUs'])

subplot(3,2,6)
pie([NoChangeY_Ad_E NoChangeY_Ad_S NoChangeY_Ad_L NoChangeY_Ad_D NoChangeY_Ad_NC])
title(['Yoked NC:'  ...
    num2str(NoChangeY_Ad_E+NoChangeY_Ad_S+NoChangeY_Ad_L+NoChangeY_Ad_D+NoChangeY_Ad_NC) ' SUs'])

TotM=NoChangeM_Ad_E+NoChangeM_Ad_S+NoChangeM_Ad_L+NoChangeM_Ad_D+NoChangeM_Ad_NC;
TotY=NoChangeY_Ad_E+NoChangeY_Ad_S+NoChangeY_Ad_L+NoChangeY_Ad_D+NoChangeY_Ad_NC;
stat_NC_Ad(1)=chisquare(NoChangeM_Ad_E,TotM,NoChangeY_Ad_E,TotY);
stat_NC_Ad(2)=chisquare(NoChangeM_Ad_S,TotM,NoChangeY_Ad_S,TotY);
stat_NC_Ad(3)=chisquare(NoChangeM_Ad_L,TotM,NoChangeY_Ad_L,TotY);
stat_NC_Ad(4)=chisquare(NoChangeM_Ad_D,TotM,NoChangeY_Ad_D,TotY);
stat_NC_Ad(5)=chisquare(NoChangeM_Ad_NC,TotM,NoChangeY_Ad_NC,TotY);

legend('Early','Steady','Late','Decrease','NC')


% AB fibers
figure
subplot(3,2,1)
pie([IncreaseM_AB_E IncreaseM_AB_S IncreaseM_AB_L IncreaseM_AB_D IncreaseM_AB_NC])
title(['Master Increase:'  ...
    num2str(IncreaseM_AB_E+IncreaseM_AB_S+IncreaseM_AB_L+IncreaseM_AB_D+IncreaseM_AB_NC) ' SUs'])

subplot(3,2,2)
pie([IncreaseY_AB_E IncreaseY_AB_S IncreaseY_AB_L IncreaseY_AB_D IncreaseY_AB_NC])
title(['Yoked Increase:'  ...
    num2str(IncreaseY_AB_E+IncreaseY_AB_S+IncreaseY_AB_L+IncreaseY_AB_D+IncreaseY_AB_NC) ' SUs'])

TotM=IncreaseM_AB_E+IncreaseM_AB_S+IncreaseM_AB_L+IncreaseM_AB_D+IncreaseM_AB_NC;
TotY=IncreaseY_AB_E+IncreaseY_AB_S+IncreaseY_AB_L+IncreaseY_AB_D+IncreaseY_AB_NC;
stat_In_AB(1)=chisquare(IncreaseM_AB_E,TotM,IncreaseY_AB_S,TotY);
stat_In_AB(2)=chisquare(IncreaseM_AB_S,TotM,NoChangeY_Ad_S,TotY);
stat_In_AB(3)=chisquare(IncreaseM_AB_L,TotM,IncreaseY_AB_L,TotY);
stat_In_AB(4)=chisquare(IncreaseM_AB_D,TotM,IncreaseY_AB_D,TotY);
stat_In_AB(5)=chisquare(IncreaseM_AB_NC,TotM,IncreaseY_AB_NC,TotY);


subplot(3,2,3)
pie([DecreaseM_AB_E DecreaseM_AB_S DecreaseM_AB_L DecreaseM_AB_D DecreaseM_AB_NC])
title(['Master Decrease:'  ...
    num2str(DecreaseM_AB_E+DecreaseM_AB_S+DecreaseM_AB_L+DecreaseM_AB_D+DecreaseM_AB_NC) ' SUs'])

subplot(3,2,4)
pie([DecreaseY_AB_E DecreaseY_AB_S DecreaseY_AB_L DecreaseY_AB_D DecreaseY_AB_NC])
title(['Yoked Decrease:'  ...
    num2str(DecreaseY_AB_E+DecreaseY_AB_S+DecreaseY_AB_L+DecreaseY_AB_D+DecreaseY_AB_NC) ' SUs'])

TotM=DecreaseM_AB_E+DecreaseM_AB_S+DecreaseM_AB_L+DecreaseM_AB_D+DecreaseM_AB_NC;
TotY=DecreaseY_AB_E+DecreaseY_AB_S+DecreaseY_AB_L+DecreaseY_AB_D+DecreaseY_AB_NC;
stat_Dec_AB(1)=chisquare(DecreaseM_AB_E,TotM,DecreaseY_AB_E,TotY);
stat_Dec_AB(2)=chisquare(DecreaseM_AB_S,TotM,DecreaseY_AB_S,TotY);
stat_Dec_AB(3)=chisquare(DecreaseM_AB_L,TotM,DecreaseY_AB_L,TotY);
stat_Dec_AB(4)=chisquare(DecreaseM_AB_D,TotM,DecreaseY_AB_D,TotY);
stat_Dec_AB(5)=chisquare(DecreaseM_AB_NC,TotM,DecreaseY_AB_NC,TotY);

subplot(3,2,5)
pie([NoChangeM_AB_E NoChangeM_AB_S NoChangeM_AB_L NoChangeM_AB_D NoChangeM_AB_NC])
title(['Master NC:'  ...
    num2str(NoChangeM_AB_E+NoChangeM_AB_S+NoChangeM_AB_L+NoChangeM_AB_D+NoChangeM_AB_NC) ' SUs'])

subplot(3,2,6)
pie([NoChangeY_AB_E NoChangeY_AB_S NoChangeY_AB_L NoChangeY_AB_D NoChangeY_AB_NC])
title(['Yoked NC:'  ...
    num2str(NoChangeY_AB_E+NoChangeY_AB_S+NoChangeY_AB_L+NoChangeY_AB_D+NoChangeY_AB_NC) ' SUs'])

TotM=NoChangeM_AB_E+NoChangeM_AB_S+NoChangeM_AB_L+NoChangeM_AB_D+NoChangeM_AB_NC;
TotY=NoChangeY_AB_E+NoChangeY_AB_S+NoChangeY_AB_L+NoChangeY_AB_D+NoChangeY_AB_NC;
stat_NC_AB(1)=chisquare(NoChangeM_AB_E,TotM,NoChangeY_AB_E,TotY);
stat_NC_AB(2)=chisquare(NoChangeM_AB_S,TotM,NoChangeY_AB_S,TotY);
stat_NC_AB(3)=chisquare(NoChangeM_AB_L,TotM,NoChangeY_AB_L,TotY);
stat_NC_AB(4)=chisquare(NoChangeM_AB_D,TotM,NoChangeY_AB_D,TotY);
stat_NC_AB(5)=chisquare(NoChangeM_AB_NC,TotM,NoChangeY_AB_NC,TotY);

legend('Early','Steady','Late','Decrease','NC')



numel(find(CheckFiberY==3))

%%

% plot reliability with bars  DOESNT WORK WELL

PerOfChangeM=(SoP_Ad(2,1:end)-SoP_Ad(1,1:end))*100;
PerOfChangeY=(SoP_Ad_Y(2,1:end)-SoP_Ad_Y(1,1:end))*100;

figure
subplot(1,2,1)
hold on
[Ad_P,idx_Ad_P] = intersect(Id_respAd, Id_respP);
Ad_P_rec = Rec_Ad(idx_Ad_P);
[Ad_Ab,idx_Ad_Ab] = intersect(Id_respAd, Id_respAB);
Ad_Ab_rec = Rec_Ad(idx_Ad_Ab);
[Ad_C,idx_Ad_C] = intersect(Id_respAd, Id_respC);
Ad_C_rec = Rec_Ad(idx_Ad_C);

for i = 1:length(PerOfChangeM)  
    if find(intersect(Id_respAd, Id_respP) == Id_respAd(i))~=0 & Ad_P_rec(find(intersect(Id_respAd, Id_respP) == Id_respAd(i))) == Rec_Ad(i)
        plot([1, PerOfChangeM(i)],'m','LineStyle','--', 'LineWidth',2)
    elseif find(intersect(Id_respAd, Id_respAB) == Id_respAd(i))~=0 & Ad_Ab_rec(find(intersect(Id_respAd, Id_respAB) == Id_respAd(i))) == Rec_Ad(i)
        plot([1, PerOfChangeM(i)],'r','LineStyle','-.', 'LineWidth',2)
    elseif find(intersect(Id_respAd, Id_respC) == Id_respAd(i))~=0 & Ad_C_rec(find(intersect(Id_respAd, Id_respC) == Id_respAd(i))) == Rec_Ad(i)
        plot([1, PerOfChangeM(i)],'r','LineStyle',':', 'LineWidth',2)
    else
        plot([1, PerOfChangeM(i)],'m')
    end
end
title('Learner')

ylim([-100 100])
subplot(1,2,2)
hold on
[Ad_PY,idx_Ad_PY] = intersect(Id_respAd_Y, Id_respP_Y);
Ad_P_recY = Rec_Ad_Y(idx_Ad_PY);
for i = 1:length(PerOfChangeY)  
    if find(intersect(Id_respAd_Y, Id_respP_Y) == Id_respAd_Y(i))~=0 & Ad_P_recY(find(intersect(Id_respAd_Y, Id_respP_Y) == Id_respAd_Y(i))) == Rec_Ad_Y(i)
        plot([1, PerOfChangeY(i)],'b','LineStyle','--', 'LineWidth',2)
    else
        plot([1, PerOfChangeY(i)],'b')
    end
end
ylim([-100 100])
title('Control')


sgtitle('Ad Reliability change DOESNT WORK WELL')