% Takeoka Lab - NERF empowered by imec, KU Leuven and VIB
% Authors: Simon Lavaud 
% 2023

% Isolate the Ad/C and Ad/C+P responses for comparison
% data used for the paper :
% SetOfSpikeStart : % of spikes/stimulation for the first X% (X=perc_shock)
% SetOfSpikeEnd : % of spikes/stimulation for the last X% (X=perc_shock)
% ReliabilityDelta : difference SetOfSpikeEnd-SetOfSpikeStart

% Load frM_UOI and frY_UOI, frM and frY for the genotype of interest

% Parameters to set :

% Velocities and distance of conduction parameters: -----------------------
% -------- > Used for paper
sABmin = 12; sABmax = 20;
sAdmin = 0.1; sAdmax = 12; 
sCmin = 0.01; sCmax = 0.09;
sPromin = 20; sPromax = 120; 
Lenght_leg=5; % average length from recording site to the furthest stimulation electrode (data not shown)
% -------------------------------------------------------------------------

% Window of responses to sensory stimulation
ABstart = (Lenght_leg/(sABmax*100))+0.001; ABstop = (Lenght_leg/(sABmin*100))+0.001;
Adstart = (Lenght_leg/(sAdmax*100))+0.001; Adstop = (Lenght_leg/(sAdmin*100))+0.001;
Cstart = (Lenght_leg/(sCmax*100))+0.001; Cstop = (Lenght_leg/(sCmin*100))+0.001;
Prostart = (Lenght_leg/(sPromax*100))+0.001; Prostop = (Lenght_leg/(sPromin*100))+0.001;

% Units that show Ad/C and P responses for the fr UOI tables
respM_all_ADCP = [20 24 103 124 139]
respY_all_ADCP = [50 69]


% Initialization
Relia_Delta_Y_ADCP=[];
Relia_Delta_M_ADCP=[];

Type_of_lat_M_ADCP=[];
Type_of_lat_Y_ADCP=[];

Type_of_cat_M_ADCP=[];
Type_of_cat_Y_ADCP=[];

Lat_M_ADCP=[];
Lat_Y_ADCP=[];

cat_lat_check_M_ADCP=[];
cat_lat_check_Y_ADCP=[];


% - - - - - - - - - - - 
% - - -  EXTRACT PARAM AD/C+P LEARNER - - - - 
% - - - - - - - - - - -

for i=1:size(respM_all_ADCP,2)
    
    lineUOI=respM_all_ADCP(i);
    
    for k=1:numel(frM_UOI.latency_all{lineUOI})
        
        if frM_UOI.latency_all{lineUOI}(k)>=0.0052 | frM_UOI.latency_all{lineUOI}(k)<=0.0035
        Relia_Delta_M_ADCP=[Relia_Delta_M_ADCP frM_UOI.ReliabilityDelta{lineUOI}(k)];
        Type_of_cat_M_ADCP=[Type_of_cat_M_ADCP frM_UOI.Category(lineUOI)];
        Lat_M_ADCP=[Lat_M_ADCP frM_UOI.latency_all{lineUOI}(k)];
            if frM_UOI.latency_all{lineUOI}(k)>=0.0052
                cat_lat_check_M_ADCP=[cat_lat_check_M_ADCP 4];
            elseif frM_UOI.latency_all{lineUOI}(k)<=0.0035
                cat_lat_check_M_ADCP=[cat_lat_check_M_ADCP 1];
            end
            
        end
        
    end
    
end

% - - - - - - - - - - - 
% - - -  EXTRACT PARAM AD/C+P CONTROL  - - - - 
% - - - - - - - - - - -

for i=1:size(respY_all_ADCP,2)
    
    lineUOI=respY_all_ADCP(i);
    
    for k=1:numel(frY_UOI.latency_all{lineUOI})
        
        if frY_UOI.latency_all{lineUOI}(k)>=0.0052 | frY_UOI.latency_all{lineUOI}(k)<=0.0035
        Relia_Delta_Y_ADCP=[Relia_Delta_Y_ADCP frY_UOI.ReliabilityDelta{lineUOI}(k)];
        Type_of_cat_Y_ADCP=[Type_of_cat_Y_ADCP frY_UOI.Category(lineUOI)];
        Lat_Y_ADCP=[Lat_Y_ADCP frY_UOI.latency_all{lineUOI}(k)];
            if frY_UOI.latency_all{lineUOI}(k)>=0.0052
                cat_lat_check_Y_ADCP=[cat_lat_check_Y_ADCP 4];
            elseif frY_UOI.latency_all{lineUOI}(k)<=0.0035
                cat_lat_check_Y_ADCP=[cat_lat_check_Y_ADCP 1];
                
            end
        end
        
    end
    
end


%%
% - - - - - - - - - - - 
% - - -  FIND THE SET OF AD/C without P  - - - - 
% - - - - - - - - - - -

% Initialization
Relia_Delta_Y_ALL=[];
Relia_Delta_M_ALL=[];
Relia_Delta_M_CNO_ALL=[];

Type_of_lat_M_ALL=[];
Type_of_lat_Y_ALL=[];
Type_of_lat_M_CNO_ALL=[];

Type_of_cat_M_ALL=[];
Type_of_cat_Y_ALL=[];
Type_of_cat_M_CNO_ALL=[];

Lat_M_ALL=[];
Lat_Y_ALL=[];
Lat_M_CNO_ALL=[];

cat_lat_check_M_ALL=[];
cat_lat_check_Y_ALL=[];
cat_lat_check_M_CNO_ALL=[];

index_ALL_UOI_M=1:size(frM_UOI,1);
index_ALL_UOI_Y=1:size(frY_UOI,1);
index_ALL_UOI_M_CNO=1:size(frM_UOI_CNO,1);

% Remove the units with Ad/C and P responses
for i=size(frM_UOI_CNO,1):-1:1
    
    if ismember(i,[45 43 26 20])
       index_ALL_UOI_M_CNO(i)=[];
    end
    
end

for i=size(frM_UOI,1):-1:1
    
    if ismember(i,[139 124 103 24 20])
       index_ALL_UOI_M(i)=[];
    end
    
end

for i=1:size(frY_UOI,1):-1:1
    
    if ismember(i,[69 50])
       index_ALL_UOI_Y(i)=[];
    end
    
end

Relia_Delta_M_ALL=[];
Relia_Delta_Y_ALL=[];
Type_of_lat_M_ALL=[];
Type_of_lat_Y_ALL=[];
Type_of_cat_M_ALL=[];
Type_of_cat_Y_ALL=[];
cat_lat_check_M_ALL=[];
cat_lat_check_Y_ALL=[];


% Run the extraction code on the set free of P responses
% For Master CNO units
for i=1:numel(index_ALL_UOI_M_CNO)
    
    lineUOI=index_ALL_UOI_M_CNO(i);
    
    for k=1:numel(frM_UOI_CNO.latency_all{lineUOI})
        
        if frM_UOI_CNO.latency_all{lineUOI}(k)>=0.0052
        Relia_Delta_M_CNO_ALL=[Relia_Delta_M_CNO_ALL frM_UOI_CNO.ReliabilityDelta{lineUOI}(k)];
        Lat_M_CNO_ALL=[Lat_M_CNO_ALL frM_UOI_CNO.latency_all{lineUOI}(k)];
            if frM_UOI_CNO.latency_all{lineUOI}(k)>=0.0052
                cat_lat_check_M_CNO_ALL=[cat_lat_check_M_CNO_ALL 4];
            end
            
        end
        
    end
    
end

% For Master without CNO manipulation units
for i=1:numel(index_ALL_UOI_M)
    
    lineUOI=index_ALL_UOI_M(i);
    
    for k=1:numel(frM_UOI.latency_all{lineUOI})
        
        if frM_UOI.latency_all{lineUOI}(k)>=0.0052
        Relia_Delta_M_ALL=[Relia_Delta_M_ALL frM_UOI.ReliabilityDelta{lineUOI}(k)];
        Lat_M_ALL=[Lat_M_ALL frM_UOI.latency_all{lineUOI}(k)];
            if frM_UOI.latency_all{lineUOI}(k)>=0.0052
                cat_lat_check_M_ALL=[cat_lat_check_M_ALL 4];
            end
            
        end
        
    end
    
end

% For Controls without CNO manipulation units
for i=1:size(index_ALL_UOI_Y,2)
    
    lineUOI=index_ALL_UOI_Y(i);
    
    for k=1:numel(frY_UOI.latency_all{lineUOI})
        
        if frY_UOI.latency_all{lineUOI}(k)>=0.0052
        Relia_Delta_Y_ALL=[Relia_Delta_Y_ALL frY_UOI.ReliabilityDelta{lineUOI}(k)];
        Lat_Y_ALL=[Lat_Y_ALL frY_UOI.latency_all{lineUOI}(k)];
            if frY_UOI.latency_all{lineUOI}(k)>=0.0052
                cat_lat_check_Y_ALL=[cat_lat_check_Y_ALL 4];
            end
            
        end
        
    end
    
end


%% All reliability together

%AdC+P
Relia_Delta_M_Inc_ADCP=Relia_Delta_M_ADCP(find(Relia_Delta_M_ADCP>=0.05));
Relia_Delta_M_Dec_ADCP=Relia_Delta_M_ADCP(find(Relia_Delta_M_ADCP<=-0.05));
Relia_Delta_Y_Inc_ADCP=Relia_Delta_Y_ADCP(find(Relia_Delta_Y_ADCP>=0.05));
Relia_Delta_Y_Dec_ADCP=Relia_Delta_Y_ADCP(find(Relia_Delta_Y_ADCP<=-0.05));

%All AdC without AdC+P
Relia_Delta_M_Inc_ALL=Relia_Delta_M_ALL(find(Relia_Delta_M_ALL>=0.05));
Relia_Delta_M_Dec_ALL=Relia_Delta_M_ALL(find(Relia_Delta_M_ALL<=-0.05));
Relia_Delta_Y_Inc_ALL=Relia_Delta_Y_ALL(find(Relia_Delta_Y_ALL>=0.05));
Relia_Delta_Y_Dec_ALL=Relia_Delta_Y_ALL(find(Relia_Delta_Y_ALL<=-0.05));

% All AdC from CNO
Relia_Delta_M_Inc_CNO_ALL=Relia_Delta_M_CNO_ALL(find(Relia_Delta_M_CNO_ALL>=0.05));
Relia_Delta_M_Dec_CNO_ALL=Relia_Delta_M_CNO_ALL(find(Relia_Delta_M_CNO_ALL<=-0.05));


% Ad/C+P
figure;
pie([numel(Relia_Delta_M_Dec_ADCP), numel(Relia_Delta_M_Inc_ADCP), numel(Relia_Delta_M_ADCP)-numel(Relia_Delta_M_Inc_ADCP)-numel(Relia_Delta_M_Dec_ADCP)]);
title('Ad/C + P - Learner')
legend('decreasing_Reliability','increasing_Reliability', 'NonChanging_Reliability')

% All Ad/C without P
figure;
pie([numel(Relia_Delta_M_Dec_ALL), numel(Relia_Delta_M_Inc_ALL), numel(Relia_Delta_M_ALL)-numel(Relia_Delta_M_Inc_ALL)-numel(Relia_Delta_M_Dec_ALL)]);
title('All Ad/C without Ad/C+P - Learner')
legend('decreasing_Reliability','increasing_Reliability', 'NonChanging_Reliability')

% All Ad/C without P with CNO manipulation
figure;
pie([numel(Relia_Delta_M_Dec_CNO_ALL), numel(Relia_Delta_M_Inc_CNO_ALL), numel(Relia_Delta_M_CNO_ALL)-numel(Relia_Delta_M_Dec_CNO_ALL)-numel(Relia_Delta_M_Inc_CNO_ALL)]);
title('All Ad/C without Ad/C+P - Learner')
legend('decreasing_Reliability','increasing_Reliability', 'NonChanging_Reliability')



%% Extent of Reliability change

% Increase, Ad/C v.s Ad/C+P
Relia_Delta_M_Inc_ADCP=Relia_Delta_M_ADCP(find(Relia_Delta_M_ADCP>=0.05 & cat_lat_check_M_ADCP==4 )) *100;
Relia_Delta_M_Inc_ALL=Relia_Delta_M_ALL(find(Relia_Delta_M_ALL>=0.05 & cat_lat_check_M_ALL==4 )) *100;

% Decrease, Ad/C v.s Ad/C+P
Relia_Delta_M_Ded_ADCP=Relia_Delta_M_ADCP(find(Relia_Delta_M_ADCP<=-0.05 & cat_lat_check_M_ADCP==4 )) *100;
Relia_Delta_M_Dec_ALL=Relia_Delta_M_ALL(find(Relia_Delta_M_ALL<=-0.05 & cat_lat_check_M_ALL==4 )) *100;



