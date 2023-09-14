% Takeoka Lab - NERF empowered by imec, KU Leuven and VIB
% Author: Charlotte Bichara 
% 2023

% Run depth_latency_tables
% Warning Fig S7C separates Ad and C . To see it separated, run
% depth_latency_tables with the variable C_Ad_separated = 1. Otherwise, run
% it with C_Ad_separated = 0.

%%

multiples=[];
multiplesY=[];
multiplesID_M=[];
multiplesID_Y=[];
multiplesIDX_M=[];
multiplesIDX_Y=[];
AB_Single_unit=[];
AB_SingleY_unit=[];
Ad_Single_unit=[];
Ad_SingleY_unit=[];
C_Single_unit=[];
C_SingleY_unit=[];
multiplesRC_M =[];
multiplesRC_Y=[];
Ad_Single_unit_RC=[];
Ad_Single_unit_RC_Y=[];
P_Single_unit=[];
P_SingleY_unit=[];
multiples_ad_RC ={};
multiples_ad_RC_Y={};
a=0;
AB_Single = 0;
Ad_Single = 0;
C_Single = 0;
P_Single = 0;

% Learner
for i = 1:length(frM_UOI.Line)   
    if length(frM_UOI.latency_all{i}) > 1
        a = a + 1;
        multiples{a} = frM_UOI.latency_all{i};
        multiplesID_M(a) = frM_UOI.Line(i);
        multiplesIDX_M(a) = i;
        multiplesRC_M{a} = frM_UOI.ReliabilityDelta(i);
    else
        if frM_UOI.latency_all{i}<Prostop
            P_Single = P_Single + 1;
            P_Single_unit(P_Single) = i;
        end
        if frM_UOI.latency_all{i}>=ABstart & frM_UOI.latency_all{i}<ABstop
            AB_Single = AB_Single + 1;
            AB_Single_unit(AB_Single) = i;
        end
        if frM_UOI.latency_all{i}>=Adstart & frM_UOI.latency_all{i}<Adstop
            Ad_Single = Ad_Single + 1;
            Ad_Single_unit(Ad_Single) = i;
            Ad_Single_unit_RC{Ad_Single} = frM_UOI.ReliabilityDelta(i);
        end
        if frM_UOI.latency_all{i}>=Cstart
            C_Single = C_Single + 1;
            C_Single_unit(C_Single) = i;
        end
    end
end


multiples_ab=[];
multiples_abY=[];
multiples_ab_idx=[];
multiples_abY_idx=[];
multiples_proY=[];
multiples_pro=[];
multiples_pro_idx=[];
multiples_proY_idx=[];
multiples_adY=[];
multiples_ad=[];
multiples_adY_idx=[];
multiples_ad_idx=[];
multiples_cY=[];
multiples_c=[];
multiples_cY_idx=[];
multiples_c_idx=[];
multiples_adY_RC=[];
multiples_ab_RC=[];


ab=0;ad=0;c=0;pro=0;
for j = 1:length(multiples)
    for i = 1:length(multiples{j})
        if multiples{j}(i)>=ABstart & multiples{j}(i)<ABstop
            ab = ab +1;
            multiples_ab(ab)=multiplesID_M(j);
            multiples_ab_idx(ab)=multiplesIDX_M(j); % line nb in frM_UOI of multiple ab

        elseif multiples{j}(i)<Prostop
            pro = pro + 1;
            multiples_pro(pro)=multiplesID_M(j);
            multiples_pro_idx(pro)=multiplesIDX_M(j); % line nb in frM_UOI of multiple P

        elseif multiples{j}(i)>=Adstart & multiples{j}(i)<Adstop
            ad=ad+1;
            multiples_ad(ad)=multiplesID_M(j);
            multiples_ad_idx(ad)=multiplesIDX_M(j);  % line nb in frM_UOI of multiple ad
            multiples_ad_RC(1, ad) = multiplesRC_M{j}; % 1. Reliability of responses / units (multi Ad)
            multiples_ad_RC(2, ad) = {i}; % 2. Position of ad response(s) in 1.
        elseif multiples{j}(i)>=Cstart
            c=c+1;
            multiples_c(c)=multiplesID_M(j);
            multiples_c_idx(c)=multiplesIDX_M(j);
        end
    end
end
AB_all = AB_Single + ab;
Ad_all = Ad_Single + ad;
C_all = C_Single + c;
P_all = P_Single + pro;


% Control

a=0;
AB_SingleY = 0;
Ad_SingleY = 0;
C_SingleY = 0;
P_SingleY = 0;
for i = 1:length(frY_UOI.Line)   
    if length(frY_UOI.latency_all{i}) > 1
        a = a + 1;
        multiplesY{a} = frY_UOI.latency_all{i};
        multiplesID_Y(a) = frY_UOI.Line(i);
        multiplesIDX_Y(a) = i;
        multiplesRC_Y{a} = frY_UOI.ReliabilityDelta(i);
    else
        if frY_UOI.latency_all{i}<Prostop
            P_SingleY = P_SingleY + 1;
            P_SingleY_unit(P_SingleY) = i;
        end
        if frY_UOI.latency_all{i}>=ABstart & frY_UOI.latency_all{i}<ABstop
            AB_SingleY = AB_SingleY + 1;
            AB_SingleY_unit(AB_SingleY) = i;
        end
        if frY_UOI.latency_all{i}>=Adstart & frY_UOI.latency_all{i}<Adstop
            Ad_SingleY = Ad_SingleY + 1;
            Ad_SingleY_unit(Ad_SingleY) = i;
            Ad_Single_unit_RC_Y{Ad_SingleY} = frY_UOI.ReliabilityDelta(i);
        end
        if frY_UOI.latency_all{i}>=Cstart
            C_SingleY = C_SingleY + 1;
            C_SingleY_unit(C_SingleY) = i;
        end
    end
end

abY=0;adY=0;cY=0;proY=0;
for j = 1:length(multiplesY)
    for i = 1:length(multiplesY{j})
        if multiplesY{j}(i)<Prostop
            proY = proY +1;
            multiples_proY(proY)=multiplesID_Y(j);
            multiples_proY_idx(proY)=multiplesIDX_Y(j);
        elseif multiplesY{j}(i)>=ABstart & multiplesY{j}(i)<ABstop
            abY = abY +1;
            multiples_abY(abY)=multiplesID_Y(j);
            multiples_abY_idx(abY)=multiplesIDX_Y(j);
        elseif multiplesY{j}(i)>=Adstart & multiplesY{j}(i)<Adstop
            adY=adY+1;
            multiples_adY(adY)=multiplesID_Y(j);
            multiples_adY_idx(adY)=multiplesIDX_Y(j);
            multiples_ad_RC_Y(1, adY) = multiplesRC_Y{j};
            multiples_ad_RC_Y(2, adY) = {i};
        elseif multiplesY{j}(i)>=Cstart
            cY=cY+1;
            multiples_cY(cY)=multiplesID_Y(j);
            multiples_cY_idx(cY)=multiplesIDX_Y(j);
        end
    end
end

AB_allY = AB_SingleY + abY;
Ad_allY = Ad_SingleY + adY;
C_allY = C_SingleY + cY;
P_allY = P_SingleY + proY;


% Stats L vs C
Stat_resp(1,1) = chisquare(AB_all, AB_all+ Ad_all+C_all+P_all, AB_allY, AB_allY+ Ad_allY+C_allY+P_allY);
Stat_resp(1,2) = chisquare(Ad_all, AB_all+ Ad_all+C_all+P_all, Ad_allY, AB_allY+ Ad_allY+C_allY+P_allY);
Stat_resp(1,3) = chisquare(C_all, AB_all+ Ad_all+C_all+P_all, C_allY, AB_allY+ Ad_allY+C_allY+P_allY);
Stat_resp(1,4) = chisquare(P_all, AB_all+ Ad_all+C_all+P_all, P_allY, AB_allY+ Ad_allY+C_allY+P_allY);


figure,
hold on
subplot(1,2,1)
pie([AB_Single, ab, Ad_Single,ad, C_Single, c, P_Single, pro]);
title('Master')
subplot(1,2,2)
pie([AB_SingleY, abY, Ad_SingleY,adY, C_SingleY, cY, P_SingleY, proY]);
title('Yoked')
legend('AB Single', 'AB multi', 'Ad Single', 'Ad multi', 'C Single', 'C multi','P Single', 'P multi');
sgtitle('Separated single and multi responses')

figure,
hold on
subplot(1,2,1)
pie([AB_all, Ad_all, C_all, P_all]);
title('Master', AB_all+ Ad_all+C_all+P_all)
subplot(1,2,2)
pie([AB_allY, Ad_allY, C_allY, P_allY]);
title('Yoked', AB_allY+ Ad_allY+C_allY+P_allY)
legend('AB', 'Ad', 'C','P');
sgtitle('Merged single and multi responses')
