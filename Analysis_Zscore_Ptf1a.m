% Takeoka Lab - NERF empowered by imec, KU Leuven and VIB
% Author: Charlotte Bichara
% 2023


% Ptf1a Analysis

% Parameters setting
thresh =80; % % of zscore>z
z = 2;
SS=20;


idxoptotaggedPtfaY = find(frY.optotagged==1 & frY.Recording>10); % Ptf1aON Controls
idxoptotaggedPtfa = find(frM.optotagged==1 & frM.Recording>9); % Ptf1aON Learners

thresh =20; % % of zscore>z corresponds to units defind as upregulated
z = 2; % zscore threshold
SS=10; % Smoothing factor


% Learner

a=0;
Perc_above_T=[];
meanZ_sel=[];
Zscores_allH=[];
Smth_Zscores_allH=[];
meanZ=[];
for j=1:length(idxoptotaggedPtfa)
    i=idxoptotaggedPtfa(j);
    
    Zscores_allH(j,:) = [frM.zscore_early_full{i} frM.zscore_mid_full{i} frM.zscore_late_full{i} frM.zscore_after_full{i}(1:600-length([frM.zscore_early_full{i} frM.zscore_mid_full{i} frM.zscore_late_full{i}]))];
    Smth_Zscores_allH(j,:)=smooth(Zscores_allH(j,:),SS);
    meanZ(j) = mean(Zscores_allH(j,:));
    Early_ZscoresM = [frM.zscore_early_full{i} frM.zscore_mid_full{i}];  
     
    Perc_above_T(j) = (numel(find(Zscores_allH(j,:) > z))*100)/numel(find(Zscores_allH(j,:)));
    if Perc_above_T(j)>thresh
        a=a+1;
        meanZ_sel(a) = mean(Zscores_allH(j,:));
        Zscores_allH_sel(a,:) = smooth(Zscores_allH(j,:),SS);
        hold on
    end
end

% Control

b=0;
Perc_above_TY=[];
meanZY_sel=[];
Zscores_allHY=[];
Smth_Zscores_allHY=[];
meanZY=[];
for j=1:length(idxoptotaggedPtfaY)
    i=idxoptotaggedPtfaY(j);
    
    Zscores_allHY(j,:) = [frY.zscore_early_full{i} frY.zscore_mid_full{i} frY.zscore_late_full{i} frY.zscore_after_full{i}(1:600-length([frY.zscore_early_full{i} frY.zscore_mid_full{i} frY.zscore_late_full{i}]))];
    Smth_Zscores_allHY(j,:)=smooth(Zscores_allHY(j,:),SS);
    meanZY(j) = mean(Zscores_allHY(j,:));
    Perc_above_TY(j) = (numel(find(Zscores_allHY(j,:) > z))*100)/numel(find(Zscores_allHY(j,:)));

    if Perc_above_TY(j)>thresh
        b=b+1;
        meanZY_sel(b) = mean(Zscores_allHY(j,:));
        Zscores_allHY_sel(b,:) = smooth(Zscores_allHY(j,:),SS);
        hold on
    end
end


figure;
stdshade(Zscores_allH_sel(:,:),0.5,'m');
hold on; 
stdshade(Zscores_allHY_sel(:,:),0.5,'b');