% Takeoka Lab - NERF empowered by imec, KU Leuven and VIB
% Author: Charlotte Bichara 
% 2023

%% Tables required, in : GitHub\takeokalab\Updated fr table
% load frM (every genotype tables) 
% frM_UOI 
% frY (every genotype tables) 
% frY_UOI 
% frM_Ptf1a_CNO 
% frM_UOI_CNO 


% Measure the change of firing rate or z-score during the conditioning phase
% or full 10 min of test, and check the proportion of each amongst Up, Dow
% regulated and  non changing units (see methods)

% Choose the % of phase for measurments : ex : 10 : compare the first and last 10% of the selected phase
P = 10; % Used for the paper
% Threshold to group increased/decreased slopes (Slopes not used in the paper)
Tresh = 0.1; 

% choose to work on Firing rates or z-scores
firing_rate = 0; % 0 = Z-scores, 1 = FR
% Choose to work on pre criterion phase or the whole 10min
Earlyphase = 0; % 0 = All 10 min trial, 1 = Only pre-criterion

% Data source
masterO = find(frM.optotagged==1 & frM.Recording>10);% & frM.depth < 500); % Ptf1aON Learner
masterI = find(frM.inhibition==1 & frM.Recording>10 & frM.Category ~= 0); % Targets of Ptf1aON Learner
controlO = find(frY.optotagged==1 & frY.Recording>9);% & frY.depth < 500); % Ptf1aON control
controlI = find(frY.inhibition==1 & frY.Recording>9 & frY.Category ~= 0); % Targets of Ptf1aON control

% create the list of second order neurons for each mouse lines - to only
% focus on the non-second-order units

a=0; b=0;c=0;d=0;
for i = 1:length(frM_UOI.Recording)
    if frM_UOI.Recording(i) < 100
        a=a+1;
        R_WT_Ptf1a_tlx3_M(a)= frM_UOI.Line(i);
    elseif frM_UOI.Recording(i)  >= 100 & frM_UOI.Recording(i)  < 1000
        b=b+1;
        R_En1_M(b) = frM_UOI.Line(i);
    elseif frM_UOI.Recording(i)  >= 1000 & frM_UOI.Recording(i)  < 10000
        c=c+1;
        R_vGlut_M(c) = frM_UOI.Line(i);
    elseif frM_UOI.Recording(i)  >= 10000 
        d=d+1;
        R_vGat_M(d) = frM_UOI.Line(i);
    end
end
a=0; b=0;c=0;d=0;
for i = 1:length(frY_UOI.Recording)
    if frY_UOI.Recording(i)  < 100
        a=a+1;
        R_WT_Ptf1a_tlx3_Y(a)= frY_UOI.Line(i);
    elseif frY_UOI.Recording(i)  >= 100 & frY_UOI.Recording(i)  < 1000
        b=b+1;
        R_En1_Y(b) = frY_UOI.Line(i);
    elseif frY_UOI.Recording(i)  >= 1000 & frY_UOI.Recording(i)  < 10000
        c=c+1;
        R_vGlut_Y(c) = frY_UOI.Line(i);
    elseif frY_UOI.Recording(i)  >= 10000 
        d=d+1;
        R_vGat_Y(d) = frY_UOI.Line(i);
    end
end

% Creates list of the non-second-order units
%Learner
NR_WT_Ptf1a_tlx3_M=setdiff(find(frM.Recording), R_WT_Ptf1a_tlx3_M);
NR_En1_M=setdiff(find(frM_En1.Recording), R_En1_M);
NR_vGlut_M=setdiff(find(frM_vGlut.Recording), R_vGlut_M);
NR_vGat_M=setdiff(find(frM_vGat.Recording), R_vGat_M);
%Control
NR_WT_Ptf1a_tlx3_Y=setdiff(find(frY.Recording), R_WT_Ptf1a_tlx3_Y);
NR_En1_Y=setdiff(find(frY_En1.Recording), R_En1_Y);
NR_vGlut_Y=setdiff(find(frY_vGlut.Recording), R_vGlut_Y);
NR_vGat_Y=setdiff(find(frY_vGat.Recording), R_vGat_Y);

% Remove Ptf1aON as well ---Can Remove other cell type if required ----
%Learner
NR_WT_Ptf1a_tlx3_M=NR_WT_Ptf1a_tlx3_M(find(frM.optotagged(NR_WT_Ptf1a_tlx3_M) == 0));
% NR_En1_M=NR_En1_M(find(frM_En1.optotagged(NR_En1_M) == 0));
% NR_vGlut_M=NR_vGlut_M(find(frM_vGlut.optotagged(NR_vGlut_M) == 0));
% NR_vGat_M=NR_vGat_M(find(frM_vGat.optotagged(NR_vGat_M) == 0));
%Control
NR_WT_Ptf1a_tlx3_Y=NR_WT_Ptf1a_tlx3_Y(find(frY.optotagged(NR_WT_Ptf1a_tlx3_Y) == 0));
% NR_En1_Y=NR_En1_Y(find(frY_En1.optotagged(NR_En1_Y) == 0));
% NR_vGlut_Y=NR_vGlut_Y(find(frY_vGlut.optotagged(NR_vGlut_Y) == 0));
% NR_vGat_Y=NR_vGat_Y(find(frY_vGat.optotagged(NR_vGat_Y) == 0));
%---------------------------------------------------------

% Concatenated list
masterAll = [NR_WT_Ptf1a_tlx3_M', NR_En1_M(find(frM_En1.Recording(NR_En1_M) <3 | frM_En1.Recording(NR_En1_M) == 5 | frM_En1.Recording(NR_En1_M) >6))', NR_vGat_M', NR_vGlut_M'];
controlAll = [NR_WT_Ptf1a_tlx3_Y', NR_En1_Y(find(frY_En1.Recording(NR_En1_Y) <3 | frY_En1.Recording(NR_En1_Y) == 5 | frY_En1.Recording(NR_En1_Y) >6))', NR_vGat_Y', NR_vGlut_Y'];

% All = all units but the ones removed above. O = Opto, I = Inhibited :
% ONLY 1: All, used for the paper

listM = {masterAll,  masterO, masterI};
titlesM ={'All', 'Ptf1aON', 'Ptf1a targets'};

listY = {controlAll,  controlO, controlI};
titlesY ={'All', 'Ptf1aON', 'Ptf1a targets'};

count = 0;

for l = 1%:length(listM)
    %Learner

    RatioD_M = [];
    RatioD_Y  = [];
    increasedUnits  = [];
    increasedUnitsY  = [];
    decreasedUnits  = [];
    decreasedUnitsY  = [];
    DeltaD_M  = [];
    DeltaD_Y  = [];
    DeltaI_M  = [];
    DeltaI_Y = [];
    depthD_Y  = [];
    depthI_Y  = [];
    depthNP_M  = [];
    depthNP_Y  = [];
    depthD_M  = [];
    depthI_M = [];
    DeltaNP_Y = [];
    DeltaNP_M = [];
    NPUnitsY= [];
    NPUnitsM= [];
    Cat_M =[];
    Cat_Y =[];
    CatD_M=[]; CatI_M=[];CatNP_M=[];
    CatD_Y=[]; CatI_Y=[];CatNP_Y=[];
    Z_IY= [];
    Z_I= [];
    a=0;
    b=0;
    c=0;
    depthI_Up_Y=[];
    depthI_Up_M=[];
    idxDUp=0;
    idxIUp=0;
    idxINC=0;
    idxDUpY=0;
    idxIUpY=0;
    idxINCY=0;

    for j = 1:length(listM{l})
        i=listM{l}(j);
        clear Fr_EarlyMid Fr_Late Z_EarlyMid Z_EarlyMidY
 % 10% early vs 10% late =========    

        if Earlyphase == 0
            if firing_rate == 0
                if j<= length(NR_WT_Ptf1a_tlx3_M)  % Tlx3 wt Ptf1a
                    Z_EarlyMid = [frM.zscore_early_full{i,1}, frM.zscore_mid_full{i,1}, frM.zscore_late_full{i,1}];
                elseif j> length(NR_WT_Ptf1a_tlx3_M) && j <= length(NR_WT_Ptf1a_tlx3_M)+length(NR_En1_M)  % Engrailed
                    Z_EarlyMid = [frM_En1.zscore_early_full{i,1}, frM_En1.zscore_mid_full{i,1}, frM_En1.zscore_late_full{i,1}];
                elseif j> length(NR_WT_Ptf1a_tlx3_M)+length(NR_En1_M) && j <= length(NR_WT_Ptf1a_tlx3_M)+length(NR_En1_M)+length(NR_vGat_M) % vGat
                    Z_EarlyMid = [frM_vGat.zscore_early_full{i,1}, frM_vGat.zscore_mid_full{i,1}, frM_vGat.zscore_late_full{i,1}];
                elseif j> length(NR_WT_Ptf1a_tlx3_M)+length(NR_En1_M) && j <= length(NR_WT_Ptf1a_tlx3_M)+length(NR_En1_M)+length(NR_vGat_M)+length(NR_vGlut_M) % vGlut2
                    Z_EarlyMid = [frM_vGlut.zscore_early_full{i,1}, frM_vGlut.zscore_mid_full{i,1}, frM_vGlut.zscore_late_full{i,1}];
                end
            else
                if j<= length(NR_WT_Ptf1a_tlx3_M)  % Tlx3 wt Ptf1a
                    Z_EarlyMid = [frM.Fr_early{i,1}, frM.Fr_middle{i,1}, frM.Fr_late{i,1}];
                elseif j> length(NR_WT_Ptf1a_tlx3_M) && j <= length(NR_WT_Ptf1a_tlx3_M)+length(NR_En1_M)  % Engrailed
                    Z_EarlyMid = [frM_En1.Fr_early{i,1}, frM_En1.Fr_middle{i,1}, frM_En1.Fr_late{i,1}];
                elseif j> length(NR_WT_Ptf1a_tlx3_M)+length(NR_En1_M) && j <= length(NR_WT_Ptf1a_tlx3_M)+length(NR_En1_M)+length(NR_vGat_M) % vGat
                    Z_EarlyMid = [frM_vGat.Fr_early{i,1}, frM_vGat.Fr_middle{i,1}, frM_vGat.Fr_late{i,1}];
                elseif j> length(NR_WT_Ptf1a_tlx3_M)+length(NR_En1_M) && j <= length(NR_WT_Ptf1a_tlx3_M)+length(NR_En1_M)+length(NR_vGat_M)+length(NR_vGlut_M) % vGlut2
                    Z_EarlyMid = [frM_vGlut.Fr_early{i,1}, frM_vGlut.Fr_middle{i,1}, frM_vGlut.Fr_late{i,1}];
                end
            end
        else
            if firing_rate == 0
                if j<= length(NR_WT_Ptf1a_tlx3_M)  % Tlx3 wt Ptf1a
                    Z_EarlyMid = [frM.zscore_early_full{i,1}, frM.zscore_mid_full{i,1}];
                elseif j> length(NR_WT_Ptf1a_tlx3_M) && j <= length(NR_WT_Ptf1a_tlx3_M)+length(NR_En1_M)  % Engrailed
                    Z_EarlyMid = [frM_En1.zscore_early_full{i,1}, frM_En1.zscore_mid_full{i,1}];
                elseif j> length(NR_WT_Ptf1a_tlx3_M)+length(NR_En1_M) && j <= length(NR_WT_Ptf1a_tlx3_M)+length(NR_En1_M)+length(NR_vGat_M) % vGat
                    Z_EarlyMid = [frM_vGat.zscore_early_full{i,1}, frM_vGat.zscore_mid_full{i,1}];
                elseif j> length(NR_WT_Ptf1a_tlx3_M)+length(NR_En1_M) && j <= length(NR_WT_Ptf1a_tlx3_M)+length(NR_En1_M)+length(NR_vGat_M)+length(NR_vGlut_M) % vGlut2
                    Z_EarlyMid = [frM_vGlut.zscore_early_full{i,1}, frM_vGlut.zscore_mid_full{i,1}];
                end
            else
                if j<= length(NR_WT_Ptf1a_tlx3_M)  % Tlx3 wt Ptf1a
                    Z_EarlyMid = [frM.Fr_early{i,1}, frM.Fr_middle{i,1}, frM.Fr_late{i,1}];
                elseif j> length(NR_WT_Ptf1a_tlx3_M) && j <= length(NR_WT_Ptf1a_tlx3_M)+length(NR_En1_M)  % Engrailed
                    Z_EarlyMid = [frM_En1.Fr_early{i,1}, frM_En1.Fr_middle{i,1}];
                elseif j> length(NR_WT_Ptf1a_tlx3_M)+length(NR_En1_M) && j <= length(NR_WT_Ptf1a_tlx3_M)+length(NR_En1_M)+length(NR_vGat_M) % vGat
                    Z_EarlyMid = [frM_vGat.Fr_early{i,1}, frM_vGat.Fr_middle{i,1}];
                elseif j> length(NR_WT_Ptf1a_tlx3_M)+length(NR_En1_M) && j <= length(NR_WT_Ptf1a_tlx3_M)+length(NR_En1_M)+length(NR_vGat_M)+length(NR_vGlut_M) % vGlut2
                    Z_EarlyMid = [frM_vGlut.Fr_early{i,1}, frM_vGlut.Fr_middle{i,1}];
                end
            end        
        end
        if j<= length(NR_WT_Ptf1a_tlx3_M)  % Tlx3 wt Ptf1a
            Cat_M(j) = frM.Category(i);
        elseif j> length(NR_WT_Ptf1a_tlx3_M) && j <= length(NR_WT_Ptf1a_tlx3_M)+length(NR_En1_M)  % Engrailed
            Cat_M(j) = frM_En1.Category(i);
        elseif j> length(NR_WT_Ptf1a_tlx3_M)+length(NR_En1_M) && j <= length(NR_WT_Ptf1a_tlx3_M)+length(NR_En1_M)+length(NR_vGat_M) % vGat
            Cat_M(j) = frM_vGat.Category(i);
        elseif j> length(NR_WT_Ptf1a_tlx3_M)+length(NR_En1_M) && j <= length(NR_WT_Ptf1a_tlx3_M)+length(NR_En1_M)+length(NR_vGat_M)+length(NR_vGlut_M) % vGlut2
            Cat_M(j) = frM_vGlut.Category(i);
        end

        perc = round(P*(length(Z_EarlyMid)/100));
        Z_initial = mean(Z_EarlyMid(1:perc));
        Z_terminal = mean(Z_EarlyMid((length(Z_EarlyMid)-perc):end));
                
%       =====================        

        if Z_terminal+std(Z_EarlyMid((length(Z_EarlyMid)-perc):end))<Z_initial%-std(Z_EarlyMid(1:perc))
            a=a+1;
            decreasedUnits(a)=i;
            RatioD_M(a) = Z_terminal/Z_initial;
            DeltaD_M(a) = Z_terminal - Z_initial;
            depthD_M(a) = frM.depth(i);

            CatD_M(1,a) = Cat_M(j);
            CatD_M(2,a) = 2;

            if Cat_M(j) > 0 & Cat_M(j) < 7
                idxDUp = idxDUp+1;
                depthD_Up_M(idxDUp) = frM.depth(i);
            end

        elseif Z_terminal>Z_initial+std(Z_EarlyMid(1:perc)) %Z_terminal-std(Z_EarlyMid((length(Z_EarlyMid)-perc):end))>Z_initial+std(Z_EarlyMid(1:perc))
            b=b+1;
            increasedUnits(b)=i;
            DeltaI_M(b) = Z_terminal - Z_initial;
            depthI_M(b) = frM.depth(i);
            CatI_M(1,b) = Cat_M(j);
            CatI_M(2,b) = 1;
            if Cat_M(j) > 0 & Cat_M(j) < 7
                idxIUp = idxIUp+1;
                depthI_Up_M(idxIUp) = frM.depth(i);
            end
        else
            c=c+1;
            DeltaNP_M(c) = Z_terminal - Z_initial;
            NPUnitsM(c)=i;
            depthNP_M(c) = frM.depth(i);
            CatNP_M(1,c) = Cat_M(j);
            CatNP_M(2,c) = 0;

        end
    end
    

    % yoked
    
    a=0;
    b=0;
    c=0;

    for j = 1:length(listY{l}); 

        i=listY{l}(j); 
        clear Fr_EarlyMidY Fr_LateY
        if Earlyphase == 0
            if firing_rate == 0
                if j<= length(NR_WT_Ptf1a_tlx3_Y)  % Tlx3 wt Ptf1a
                    Z_EarlyMidY = [frY.zscore_early_full{i,1}, frY.zscore_mid_full{i,1}, frY.zscore_late_full{i,1}];
                elseif j> length(NR_WT_Ptf1a_tlx3_Y) && j <= length(NR_WT_Ptf1a_tlx3_Y)+length(NR_En1_Y)  % Engrailed
                    Z_EarlyMidY = [frY_En1.zscore_early_full{i,1}, frY_En1.zscore_mid_full{i,1}, frY_En1.zscore_late_full{i,1}];
                elseif j> length(NR_WT_Ptf1a_tlx3_Y)+length(NR_En1_Y) && j <= length(NR_WT_Ptf1a_tlx3_Y)+length(NR_En1_Y)+length(NR_vGat_Y) % vGat
                    Z_EarlyMidY = [frY_vGat.zscore_early_full{i,1}, frY_vGat.zscore_mid_full{i,1}, frY_vGat.zscore_late_full{i,1}];
                elseif j> length(NR_WT_Ptf1a_tlx3_Y)+length(NR_En1_Y)+length(NR_vGat_Y) && j <= length(NR_WT_Ptf1a_tlx3_Y)+length(NR_En1_Y)+length(NR_vGat_Y)+length(NR_vGlut_Y) % vGlut2
                    Z_EarlyMidY = [frY_vGlut.zscore_early_full{i,1}, frY_vGlut.zscore_mid_full{i,1}, frY_vGlut.zscore_late_full{i,1}];
                end
            else
                if j<= length(NR_WT_Ptf1a_tlx3_Y)   % Tlx3 wt Ptf1a
                    Z_EarlyMidY = [frY.Fr_early{i,1}, frY.Fr_middle{i,1}, frY.Fr_late{i,1}];
                elseif j> length(NR_WT_Ptf1a_tlx3_Y) && j <= length(NR_WT_Ptf1a_tlx3_Y)+length(NR_En1_Y)   % Engrailed
                    Z_EarlyMidY = [frY_En1.Fr_early{i,1}, frY_En1.Fr_middle{i,1}, frY_En1.Fr_late{i,1}];
                elseif j> length(NR_WT_Ptf1a_tlx3_Y)+length(NR_En1_Y) && j <= length(NR_WT_Ptf1a_tlx3_Y)+length(NR_En1_Y)+length(NR_vGat_Y) % vGat
                    Z_EarlyMidY = [frY_vGat.Fr_early{i,1}, frY_vGat.Fr_middle{i,1}, frY_vGat.Fr_late{i,1}];
                elseif j> length(NR_WT_Ptf1a_tlx3_Y)+length(NR_En1_Y)+length(NR_vGat_Y) && j <= length(NR_WT_Ptf1a_tlx3_Y)+length(NR_En1_Y)+length(NR_vGat_Y)+length(NR_vGlut_Y) % vGlut2
                    Z_EarlyMidY = [frY_vGlut.Fr_early{i,1}, frY_vGlut.Fr_middle{i,1}, frY_vGlut.Fr_late{i,1}];
                end
            end
        else
        
            if firing_rate == 0
                if j<= length(NR_WT_Ptf1a_tlx3_Y)  % Tlx3 wt Ptf1a
                    Z_EarlyMidY = [frY.zscore_early_full{i,1}, frY.zscore_mid_full{i,1}];
                elseif j> length(NR_WT_Ptf1a_tlx3_Y) && j <= length(NR_WT_Ptf1a_tlx3_Y)+length(NR_En1_Y)  % Engrailed
                    Z_EarlyMidY = [frY_En1.zscore_early_full{i,1}, frY_En1.zscore_mid_full{i,1}];
                elseif j> length(NR_WT_Ptf1a_tlx3_Y)+length(NR_En1_Y) && j <= length(NR_WT_Ptf1a_tlx3_Y)+length(NR_En1_Y)+length(NR_vGat_Y) % vGat
                    Z_EarlyMidY = [frY_vGat.zscore_early_full{i,1}, frY_vGat.zscore_mid_full{i,1}];
                elseif j> length(NR_WT_Ptf1a_tlx3_Y)+length(NR_En1_Y)+length(NR_vGat_Y) && j <= length(NR_WT_Ptf1a_tlx3_Y)+length(NR_En1_Y)+length(NR_vGat_Y)+length(NR_vGlut_Y) % vGlut2
                    Z_EarlyMidY = [frY_vGlut.zscore_early_full{i,1}, frY_vGlut.zscore_mid_full{i,1}];
                end
            else
                if j<= length(NR_WT_Ptf1a_tlx3_Y)   % Tlx3 wt Ptf1a
                    Z_EarlyMidY = [frY.Fr_early{i,1}, frY.Fr_middle{i,1}];
                elseif j> length(NR_WT_Ptf1a_tlx3_Y) && j <= length(NR_WT_Ptf1a_tlx3_Y)+length(NR_En1_Y)   % Engrailed
                    Z_EarlyMidY = [frY_En1.Fr_early{i,1}, frY_En1.Fr_middle{i,1}];
                elseif j> length(NR_WT_Ptf1a_tlx3_Y)+length(NR_En1_Y) && j <= length(NR_WT_Ptf1a_tlx3_Y)+length(NR_En1_Y)+length(NR_vGat_Y) % vGat
                    Z_EarlyMidY = [frY_vGat.Fr_early{i,1}, frY_vGat.Fr_middle{i,1}];
                elseif j> length(NR_WT_Ptf1a_tlx3_Y)+length(NR_En1_Y)+length(NR_vGat_Y) && j <= length(NR_WT_Ptf1a_tlx3_Y)+length(NR_En1_Y)+length(NR_vGat_Y)+length(NR_vGlut_Y) % vGlut2
                    Z_EarlyMidY = [frY_vGlut.Fr_early{i,1}, frY_vGlut.Fr_middle{i,1}];
                end
            end
        end

        if j<= length(NR_WT_Ptf1a_tlx3_Y)  % Tlx3 wt Ptf1a
            Cat_Y(j) = frY.Category(i);
        elseif j> length(NR_WT_Ptf1a_tlx3_Y) && j <= length(NR_WT_Ptf1a_tlx3_Y)+length(NR_En1_Y)  % Engrailed
            Cat_Y(j) = frY_En1.Category(i);
        elseif j> length(NR_WT_Ptf1a_tlx3_Y)+length(NR_En1_Y) && j <= length(NR_WT_Ptf1a_tlx3_Y)+length(NR_En1_Y)+length(NR_vGat_Y) % vGat
            Cat_Y(j) = frY_vGat.Category(i);
        elseif j> length(NR_WT_Ptf1a_tlx3_Y)+length(NR_En1_Y) && j <= length(NR_WT_Ptf1a_tlx3_Y)+length(NR_En1_Y)+length(NR_vGat_Y)+length(NR_vGlut_Y) % vGlut2
            Cat_Y(j) = frY_vGlut.Category(i);
        end
        
        percY = round(P*(length(Z_EarlyMidY)/100));
        Z_initialY = mean(Z_EarlyMidY(1:percY));
        Z_terminalY = mean(Z_EarlyMidY((length(Z_EarlyMidY)-percY):end));
        
        if Z_terminalY+std(Z_EarlyMidY((length(Z_EarlyMidY)-percY):end))<Z_initialY%-std(Z_EarlyMidY(1:percY))
            a=a+1;
            decreasedUnitsY(a)=i;
            RatioD_Y(a) = Z_terminalY/Z_initialY;
            DeltaD_Y(a) = Z_terminalY - Z_initialY;
            depthD_Y(a) = frY.depth(i);
            CatD_Y(1,a) = Cat_Y(j);
            CatD_Y(2,a) = 2;
            if Cat_Y(j) > 0 & Cat_Y(j) < 7
                idxDUpY = idxDUpY+1;
                depthD_Up_Y(idxDUpY) = frY.depth(i);
             end
        elseif Z_terminalY>Z_initialY+std(Z_EarlyMidY(1:percY)) %Z_terminalY-std(Z_EarlyMidY((length(Z_EarlyMidY)-percY):end))>Z_initialY+std(Z_EarlyMidY(1:percY))
            b=b+1;
            increasedUnitsY(b)=i;
            DeltaI_Y(b) = Z_terminalY - Z_initialY;
            depthI_Y(b) = frY.depth(i);
            Z_IY(b) = mean([frY.zscore_late_full{i,1},  frY.zscore_after_full{i,1}(1:600-length(frY.zscore_late_full{i,1}))]);
            CatI_Y(1,b) = Cat_Y(j);
            CatI_Y(2,b) = 1;
           
             if Cat_Y(j) > 0 & Cat_Y(j) < 7
                idxIUpY = idxIUpY+1;
                depthI_Up_Y(idxIUpY) = frY.depth(i);
             end

        else
            c=c+1;
            DeltaNP_Y(c) = Z_terminalY - Z_initialY;
            NPUnitsY(c)=i;
            depthNP_Y(c) = frY.depth(i);
            CatNP_Y(1,c) = Cat_Y(j);
            CatNP_Y(2,c) = 0;
        end
    end

    % Figures giving % of each cat in Increase/Decrease grps.
    catmod_M = [CatI_M CatD_M CatNP_M];
    catmod_Y = [CatI_Y CatD_Y CatNP_Y];
         Up = find(catmod_M(1,:) > 0 & catmod_M(1,:) <=6);
         Down = find(catmod_M(1,:) > 6 & catmod_M(1,:) <=9);
         NC = find(catmod_M(1,:) == 10 | isnan(catmod_M(1,:)));
     figure
     CatDelta = {Up, Down, NC};
     Cattitle = {'UP', 'Down', 'NC'};
     for i = 1:length(CatDelta)
         j=CatDelta{i};
         if i==1
            UP_I_M = numel(find(catmod_M(2,j) == 1));
            UP_D_M = numel(find(catmod_M(2,j) == 2));
            UP_NC_M = numel(find(catmod_M(2,j) == 0));
         elseif i==2
            D_I_M = numel(find(catmod_M(2,j) == 1));
            D_D_M = numel(find(catmod_M(2,j) == 2));
            D_NC_M = numel(find(catmod_M(2,j) == 0));
        elseif i==3
            NC_I_M = numel(find(catmod_M(2,j) == 1));
            NC_D_M = numel(find(catmod_M(2,j) == 2));
            NC_NC_M = numel(find(catmod_M(2,j) == 0));
         end
     end
         subplot(1,2,1)
         pie([NC_NC_M, NC_D_M, NC_I_M, D_NC_M, D_D_M, D_I_M, UP_NC_M, UP_D_M, UP_I_M])
         title('Learner')
         Up_Y = find(catmod_Y(1,:) > 0 & catmod_Y(1,:) <=6);
         Down_Y = find(catmod_Y(1,:) > 6 & catmod_Y(1,:) <=9);
         NC_Y = find(catmod_Y(1,:) == 10 | isnan(catmod_Y(1,:)));

     CatDeltaY = {Up_Y, Down_Y, NC_Y};
     for i = 1:length(CatDeltaY)
         j=CatDeltaY{i};
         if i==1
            UP_I_Y = numel(find(catmod_Y(2,j) == 1));
            UP_D_Y = numel(find(catmod_Y(2,j) == 2));
            UP_NC_Y = numel(find(catmod_Y(2,j) == 0));
         elseif i==2
            D_I_Y = numel(find(catmod_Y(2,j) == 1));
            D_D_Y = numel(find(catmod_Y(2,j) == 2));
            D_NC_Y = numel(find(catmod_Y(2,j) == 0));
        elseif i==3
            NC_I_Y = numel(find(catmod_Y(2,j) == 1));
            NC_D_Y = numel(find(catmod_Y(2,j) == 2));
            NC_NC_Y = numel(find(catmod_Y(2,j) == 0));
         end
     end
     subplot(1,2,2)
     pie([NC_NC_Y, NC_D_Y, NC_I_Y, D_NC_Y, D_D_Y, D_I_Y, UP_NC_Y, UP_D_Y, UP_I_Y])
     title('Control')
     sgtitle('Phase tuning for each groups changing FR during phase')
     legend('Non Changing NC', 'Decreasing NC', 'Increasing NC','Non Changing Down', 'Decreasing Down', 'Increasing Down','Non Changing Up', 'Decreasing Up', 'Increasing Up')
end

