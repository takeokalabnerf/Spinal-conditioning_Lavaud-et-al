
% Takeoka Lab - NERF empowered by imec, KU Leuven and VIB
% Author: Charlotte Bichara 
% 2023

% This codes Computes the probability of response for each second order
% type, for learner controls and CNO condition.

% Load frM_UOI, frY_UOI frM_Ptf1a_CNO, frM_UOI_CNO

% Run code depth_latency_tables (Create Tables with info for all responses and sets of usefull variables)


% Parameters to set :

% Velocities and distance of conduction parameters: -----------------------
% -------- > Used for paper
sABmin = 12; sABmax = 20;
sAdmin = 0.1; sAdmax = 12; 
sCmin = 0.01; sCmax = 0.09;
sPromin = 20; sPromax = 120; 
Lenght_leg=5; % average length from recording site to the furthest stimulation electrode (data not shown)
% -------------------------------------------------------------------------

%Threshold for Increased, decreased or no change in probability of response
threshold = 0.05;

% Standart deviation for outliers 
SD =3;

%% Create Tables with info for all responses and sets of usefull variables FOR CNO

Resp_M =find(frM_UOI_CNO.Recording);
All_latencies = cat(2, frM_UOI_CNO.latency_all{Resp_M,:});
  
tableauCNO=table;

Depths=[];
Units=[];
Rec=[];
Category=[];
d=0;
WF = [];

for i=1:length(Resp_M)
    j=Resp_M(i);
    k=frM_UOI_CNO.Line(j);
    
    
    depths=frM_UOI_CNO.depth(j)*ones(1,length(frM_UOI_CNO.latency_all{j,:}));  
    Depths(d+1:d+length(depths))=depths;
     
    category=frM_UOI_CNO.Category(j)*ones(1,length(frM_UOI_CNO.latency_all{j,:}));  
    Category(d+1:d+length(category))=category;

    units=k*ones(1,length(frM_UOI_CNO.latency_all{j,:})); 
    Units(d+1:d+length(units))=units;
    
    rec=frM_UOI_CNO.Recording(j)*ones(1,length(frM_UOI_CNO.latency_all{j,:})); 
    Rec(d+1:d+length(rec))=rec;
    
    d=d+length(depths);
end

tableauCNO.Latence=All_latencies';
tableauCNO.Unit=Units';
tableauCNO.Depth=-Depths';
tableauCNO.Recording=Rec';
tableauCNO.Category=Category';

%% Analysis second order neurons


cno_all_Ad = find(tableauCNO.Latence >= Adstart & tableauCNO.Latence < Adstop);
cno_all_Ab = find(tableauCNO.Latence >= ABstart & tableauCNO.Latence < ABstop);
cno_all_P = find(tableauCNO.Latence >= Prostart & tableauCNO.Latence < Prostop);
cno_all_C = find(tableauCNO.Latence >= Cstart);
cno_all_Ad_C = find(tableauCNO.Latence >= Adstart & tableauCNO.Latence < Adstop & tableauCNO.Latence >= Cstart);

rel_delta_cno_All = cat(2, frM_UOI_CNO.ReliabilityDelta{:});

rel_delta_cno_Ad = rel_delta_cno_All(cno_all_Ad);
rel_delta_cno_Ab = rel_delta_cno_All(cno_all_Ab);
rel_delta_cno_P = rel_delta_cno_All(cno_all_P);
rel_delta_cno_C = rel_delta_cno_All(cno_all_C);
rel_delta_cno_Ad_C = [rel_delta_cno_All(cno_all_Ad) rel_delta_cno_All(cno_all_Ad_C)];

rel_delta_All = cat(2, frM_UOI.ReliabilityDelta{:});
rel_delta_Ad = rel_delta_All(Ad);
rel_delta_Ab = rel_delta_All(Ab);
rel_delta_P = rel_delta_All(Pro);
rel_delta_C = rel_delta_All(C);
rel_delta_Ad_C = [rel_delta_All(Ad) rel_delta_All(C)];

depth_delta_Ad = tableau.Depth(Ad);
depth_delta_Ab = tableau.Depth(Ab);

rel_delta_All_Y = cat(2, frY_UOI.ReliabilityDelta{:});

rel_delta_Ad_Y = rel_delta_All_Y(Ady);
rel_delta_Ab_Y = rel_delta_All_Y(Aby);
rel_delta_P_Y = rel_delta_All_Y(Proy);
rel_delta_C_Y = rel_delta_All_Y(Cy);
rel_delta_Ad_C_Y = [rel_delta_All_Y(Ady) rel_delta_All_Y(Cy)];

depth_delta_Ad_Y = tableauY.Depth(Ady);
depth_delta_Ab_Y = tableauY.Depth(Aby);


%% remove outliers (3sd), plot decreased only, stats


for i = 1:4
    REL = []; RELY = []; REL_CNO = []; REL_D=[]; REL_CNO_D=[]; RELY_D=[];REL_I=[]; REL_CNO_I=[]; RELY_I=[]; RELY_NC=[];REL_NC=[]; REL_CNO_NC=[]; 
    if i == 1
       REL = rel_delta_Ad;RELY = rel_delta_Ad_Y;REL_CNO = rel_delta_cno_Ad;
       titre = 'Ad';
    elseif i == 2
        REL = rel_delta_Ab;RELY = rel_delta_Ab_Y;REL_CNO = rel_delta_cno_Ab;
        titre = 'Ab';
    elseif i == 3
        REL = rel_delta_P;RELY = rel_delta_P_Y;REL_CNO = rel_delta_cno_P;
        titre = 'P';
%     elseif i == 4
%         REL = rel_delta_C;RELY = rel_delta_C_Y;REL_CNO = rel_delta_cno_C;
%         titre = 'C';
    elseif i == 4 % Same as i == 1 if Ad and C already merged in the velocities window
        REL = rel_delta_Ad_C;RELY = rel_delta_Ad_C_Y;REL_CNO = rel_delta_cno_Ad_C;
        titre = 'merged Ad C';
    end
    
    % Removes outlier X std (X = variable 'SD')
    % decrease 
    if kstest(([RELY(find(RELY<-threshold))]-mean([RELY(find(RELY<-threshold))]))/std([RELY(find(RELY<-threshold))])) == 0
        RELY_D = RELY(find(RELY<-threshold  & RELY>mean(RELY(find(RELY<-threshold)))-SD*std(RELY(find(RELY<-threshold))) & RELY<mean(RELY(find(RELY<-threshold)))+SD*std(RELY(find(RELY<-threshold)))));
    else
        sprintf('%.0f non param Y', i)
        RELY_D = RELY(find(RELY<-threshold));
    end
    if kstest(([REL(find(REL<-threshold))]-mean([REL(find(REL<-threshold))]))/std([REL(find(REL<-threshold))])) == 0
        REL_D = REL(find(REL<-threshold  & REL>mean(REL(find(REL<-threshold)))-SD*std(REL(find(REL<-threshold))) & REL<mean(REL(find(REL<-threshold)))+SD*std(REL(find(REL<-threshold)))));
    else
        sprintf('%.0f non param M', i)
        REL_D = REL(find(REL<-threshold));
    end
    if length(REL_CNO(find(REL_CNO<-threshold))) > 1
        if kstest(([REL_CNO(find(REL_CNO<-threshold))]-mean([REL_CNO(find(REL_CNO<-threshold))]))/std([REL_CNO(find(REL_CNO<-threshold))])) == 0
            REL_CNO_D = REL_CNO(find(REL_CNO<-threshold  & REL_CNO>mean(REL_CNO(find(REL_CNO<-threshold)))-SD*std(REL_CNO(find(REL_CNO<-threshold))) & REL_CNO<mean(REL_CNO(find(REL_CNO<-threshold)))+SD*std(REL_CNO(find(REL_CNO<-threshold)))));
        else
            sprintf('%.0f non param CNO', i)
            REL_CNO_D = REL_CNO(find(REL_CNO<-threshold));
        end
    end
    
    %Increase

    if kstest(([RELY(find(RELY>threshold))]-mean([RELY(find(RELY>threshold))]))/std([RELY(find(RELY>threshold))])) == 0
        RELY_I = RELY(find(RELY>threshold  & RELY>mean(RELY(find(RELY>threshold)))-SD*std(RELY(find(RELY>threshold))) & RELY<mean(RELY(find(RELY>threshold)))+SD*std(RELY(find(RELY>threshold)))));
    else
        sprintf('%.0f non param Y', i)
        RELY_I = RELY(find(RELY<-threshold));
    end
    if kstest(([REL(find(REL>threshold))]-mean([REL(find(REL>threshold))]))/std([REL(find(REL>threshold))])) == 0
        REL_I = REL(find(REL>threshold  & REL>mean(REL(find(REL>threshold)))-SD*std(REL(find(REL>threshold))) & REL<mean(REL(find(REL>threshold)))+SD*std(REL(find(REL>threshold)))));
    else
        sprintf('%.0f non param M', i)
        REL_I = REL(find(REL>threshold));
    end
    if length(REL_CNO(find(REL_CNO>threshold))) > 1
        if kstest(([REL_CNO(find(REL_CNO>threshold))]-mean([REL_CNO(find(REL_CNO>threshold))]))/std([REL_CNO(find(REL_CNO>threshold))])) == 0
            REL_CNO_I = REL_CNO(find(REL_CNO>threshold  & REL_CNO>mean(REL_CNO(find(REL_CNO>threshold)))-SD*std(REL_CNO(find(REL_CNO>threshold))) & REL_CNO<mean(REL_CNO(find(REL_CNO>threshold)))+SD*std(REL_CNO(find(REL_CNO>threshold)))));
        else
            sprintf('%.0f non param CNO', i)
            REL_CNO_I = REL_CNO(find(REL_CNO>threshold));
        end
    end

        %No Change

    if kstest(([RELY(find(RELY<=threshold & RELY>=-threshold))]-mean([RELY(find(RELY<=threshold & RELY>=-threshold))]))/std([RELY(find(RELY<=threshold & RELY>=-threshold))])) == 0
        RELY_NC = RELY(find(RELY<=threshold & RELY>=-threshold & RELY>mean(RELY(find(RELY<=threshold & RELY>=-threshold)))-SD*std(RELY(find(RELY<=threshold & RELY>=-threshold))) & RELY<mean(RELY(find(RELY<=threshold & RELY>=-threshold)))+SD*std(RELY(find(RELY<=threshold & RELY>=-threshold)))));
    else
        sprintf('%.0f non param Y', i)
        RELY_NC = RELY(find(RELY<=threshold & RELY>=-threshold));
    end
    if kstest(([REL(find(REL<=threshold & REL>=-threshold))]-mean([REL(find(REL<=threshold & REL>=-threshold))]))/std([REL(find(REL<=threshold & REL>=-threshold))])) == 0
        REL_NC = REL(find(REL<=threshold & REL>=-threshold & REL>mean(REL(find(REL<=threshold & REL>=-threshold)))-SD*std(REL(find(REL<=threshold & REL>=-threshold))) & REL<mean(REL(find(REL<=threshold & REL>=-threshold)))+SD*std(REL(find(REL<=threshold & REL>=-threshold)))));
    else
        sprintf('%.0f non param Y', i)
        REL_NC = REL(find(REL<=threshold & REL>=-threshold));
    end

    if length(REL_CNO(find(REL_CNO<=threshold & REL_CNO>=-threshold))) > 1
        if kstest(([REL_CNO(find(REL_CNO<=threshold & REL_CNO>=-threshold))]-mean([REL_CNO(find(REL_CNO<=threshold & REL_CNO>=-threshold))]))/std([REL_CNO(find(REL_CNO<=threshold & REL_CNO>=-threshold))])) == 0
            REL_CNO_NC = REL_CNO(find(REL_CNO<=threshold & REL_CNO>=-threshold & REL_CNO>mean(REL_CNO(find(REL_CNO<=threshold & REL_CNO>=-threshold)))-SD*std(REL_CNO(find(REL_CNO<=threshold & REL_CNO>=-threshold))) & REL_CNO<mean(REL_CNO(find(REL_CNO<=threshold & REL_CNO>=-threshold)))+SD*std(REL_CNO(find(REL_CNO<=threshold & REL_CNO>=-threshold)))));
        else
            sprintf('%.0f non param Y', i)
            REL_CNO_NC = REL_CNO(find(REL_CNO<=threshold & REL_CNO>=-threshold));
        end
    end


    %Decrease
    figure
    boxchart(ones(numel(REL_D),1), REL_D);  
    if length(REL_CNO(find(REL_CNO<-threshold))) > 1
        hold on
        boxchart(2*ones(numel(REL_CNO_D),1), REL_CNO_D);
    end
    hold on
    boxchart(3*ones(numel(RELY_D),1), RELY_D);
    
    hold on
    scatter(ones(numel(REL_D),1),REL_D,'MarkerEdgeColor','k','MarkerFaceColor','m','jitter','on','jitterAmount', .1);
    if length(REL_CNO(find(REL_CNO<-threshold))) > 1
        hold on
        scatter(2*ones(numel(REL_CNO_D),1),REL_CNO_D,'MarkerEdgeColor','k','MarkerFaceColor','r','jitter','on','jitterAmount', .1);
    end
    hold on
    scatter(3*ones(numel(RELY_D),1),RELY_D,'MarkerEdgeColor','k','MarkerFaceColor','b','jitter','on','jitterAmount', .1);
    legend('Learner', 'Control');
    title(sprintf('Reliability %s',titre))
    xlim([0 4])
    ylim([-1 0])
    
    %Increase
    figure
    boxchart(ones(numel(REL_I),1), REL_I);  
    if length(REL_CNO(find(REL_CNO>threshold))) > 1
        hold on
        boxchart(2*ones(numel(REL_CNO_I),1), REL_CNO_I);
    end
    hold on
    boxchart(3*ones(numel(RELY_I),1), RELY_I);
    
    hold on
    scatter(ones(numel(REL_I),1),REL_I,'MarkerEdgeColor','k','MarkerFaceColor','m','jitter','on','jitterAmount', .1);
    if length(REL_CNO(find(REL_CNO>threshold))) > 1
        hold on
        scatter(2*ones(numel(REL_CNO_I),1),REL_CNO_I,'MarkerEdgeColor','k','MarkerFaceColor','r','jitter','on','jitterAmount', .1);
    end
    hold on
    scatter(3*ones(numel(RELY_I),1),RELY_I,'MarkerEdgeColor','k','MarkerFaceColor','b','jitter','on','jitterAmount', .1);
    legend('Learner', 'Control');
    title(sprintf('Reliability %s',titre))
    xlim([0 4])
    ylim([0 1])
      
    %Decrease
    if length(REL_D) > 0 & length(RELY_D) > 0
        if kstest(([REL_D]-mean([REL_D]))/std([REL_D])) == 0 & kstest(([RELY_D]-mean([RELY_D]))/std([RELY_D])) == 0 
            [a,Stat_RelValueD(i,1)] = ttest2(REL_D,RELY_D,'Tail','both');
        else
            Stat_RelValueD(i,1) = ranksum(REL_D,RELY_D,'Tail','both');
        end
        EffectSizeD(i,1)=(mean(REL_D)-mean(RELY_D))/std(RELY_D);
    end
    if length(REL_D) > 0 & length(REL_CNO_D) > 0
        if kstest(([REL_CNO_D]-mean([REL_CNO_D]))/std([REL_CNO_D])) == 0 & kstest(([REL_D]-mean([REL_D]))/std([REL_D])) == 0 
            [a,Stat_RelValueD(i,2)] = ttest2(REL_CNO_D,REL_D,'Tail','both');
        else
            Stat_RelValueD(i,2) = ranksum(REL_CNO_D,REL_D,'Tail','both');
        end
        EffectSizeD(i,2)=(mean(REL_D)-mean(REL_CNO_D))/std(REL_D);
    end
    if length(RELY_D) > 0 & length(REL_CNO_D) > 0
        if kstest(([REL_CNO_D]-mean([REL_CNO_D]))/std([REL_CNO_D])) == 0 & kstest(([RELY_D]-mean([RELY_D]))/std([RELY_D])) == 0 
            [a,Stat_RelValueD(i,3)] = ttest2(REL_CNO_D,RELY_D,'Tail','both');
        else
            Stat_RelValueD(i,3) = ranksum(REL_CNO_D,RELY_D,'Tail','both');
        end
        EffectSizeD(i,3)=(mean(RELY_D)-mean(REL_CNO_D))/std(RELY_D);
    end


    figure
    hold on
    subplot(1,3,1)
    RelGPM = [numel(REL_I), numel(REL_D), numel(REL_NC)]/sum([numel(REL_D),numel(REL_I), numel(REL_NC)]);
    spider_plot(RelGPM, 'AxesLabels', {'Increase', 'Decrease', 'Non Changing'}, 'AxesLimits', [0, 0, 0; 0.5, 0.5, 0.5],'AxesOffset',0);
    title(sprintf('Learner %.0f', numel(REL)))
    subplot(1,3,2)
    RelGPCNO = [numel(REL_CNO_I), numel(REL_CNO_D),numel(REL_CNO_NC)]/sum([numel(REL_CNO_D),numel(REL_CNO_I), numel(REL_CNO_NC)]);
    spider_plot(RelGPCNO, 'AxesLabels', {'Increase', 'Decrease', 'Non Changing'}, 'AxesLimits', [0, 0, 0; 0.7, 0.7, 0.7],'AxesOffset',0);
    title(sprintf('CNO %.0f', numel(REL_CNO)))
    subplot(1,3,3)
    RelGPY = [numel(RELY_I), numel(RELY_D),numel(RELY_NC)]/sum([numel(RELY_D),numel(RELY_I), numel(RELY_NC)]);
    spider_plot(RelGPY, 'AxesLabels', {'Increase', 'Decrease', 'Non Changing'}, 'AxesLimits', [0, 0, 0; 0.5, 0.5, 0.5],'AxesOffset',0);
    title(sprintf('Control %.0f', numel(RELY)))
    sgtitle(sprintf('No Outliers ... Reliability %s',titre))
    legend('decrease', 'increase', 'no change');

end

