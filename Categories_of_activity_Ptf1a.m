% Takeoka Lab - NERF empowered by imec, KU Leuven and VIB
% Author: Charlotte Bichara
% 2023

% Plot depths and modulation of the units during conditioning

Anesthesia = 0; % 0 corresponds to Anesthesia condition as described in the methods

if Anesthesia == 1
    IdxoptotaggedPtfa = find(frM.optotagged==1 & frM.Recording>30);
    IdxoptotaggedPtfaY = find(frY.optotagged==1 & frY.Recording>28);
elseif Anesthesia == 0
    IdxoptotaggedPtfa = find(frM.optotagged==1  & frM.Recording>10);
    IdxoptotaggedPtfaY = find(frY.optotagged==1  & frY.Recording>9);
elseif Anesthesia == 2
    IdxoptotaggedPtfa = find(frM.optotagged==1 & frM.Recording>10 & frM.Recording<=30);
    IdxoptotaggedPtfaY = find(frY.optotagged==1 & frY.Recording>9 & frY.Recording<=28);
end

%% Depths of ptf1a neurons

figure;
hold on
scatter(1*ones(1,length(frM.depth(find(frM.optotagged==1  & frM.Recording>10)))), -frM.depth(find(frM.optotagged==1  & frM.Recording>10)),'MarkerEdgeColor','k','MarkerFaceColor','m','jitter','on','jitterAmount', .2);
scatter(2*ones(1,length(frY.depth(find(frY.optotagged==1 & frY.Recording>9)))), -frY.depth(find(frY.optotagged==1 & frY.Recording>9)),'MarkerEdgeColor','k','MarkerFaceColor','b','jitter','on','jitterAmount', .2);
alpha(0.5)
xlim([0 3])
ylim([-1800 0])
title('depths of Ptf1aON units')

%% Learner
%Photoactivated

% Categories :
% Up-regulation (zscore >2 for >20% of phase)
% E = Early, pre-criterion
% L = Late, post-criterion
% EL = Early + Late
% Down-regulation (zscore <-2 for >20% of phase)
% Same ith index "d"
% Other = mix of Up and down, minor, excluded from analysis

% Divided in sub category of activity not used in the paper.
% In the paper, only Up (E+L+EL), Down (Ed, ELd, Ld) and no change (NP) are used

CategoryM = [0 0 0 0 0 0 0 0]; % 1: EL, 2 = E, 3= L, 4= ELd, 5= Ed, 6=Ld, 7= NP, 8= Other

% Not used in paper --------
weak_strong_both_1=[0 0 0 0];
weak_strong_both_2=[0 0];
weak_strong_both_3=[0 0];
% ---------------------------

for j = 1:length(IdxoptotaggedPtfa)
    i=IdxoptotaggedPtfa(j);
    
    if frM.Strong_firing(i) == 1

        if (frM.zscore_group_EandMmerge(i)==1 | frM.zscore_group_EandMmerge(i)==10 ) & ( frM.zscore_group_late(i)==1 | frM.zscore_group_late(i)==10)
            CategoryM(1)=CategoryM(1)+1;
            a=1;
            if frM.zscore_group_EandMmerge(i)==1 & frM.zscore_group_late(i)==1
                weak_strong_both_1(1)=weak_strong_both_1(1)+1;
            end
            if frM.zscore_group_EandMmerge(i)==1 & frM.zscore_group_late(i)==10
                weak_strong_both_1(2)=weak_strong_both_1(2)+1;
            end
            if frM.zscore_group_EandMmerge(i)==10 & frM.zscore_group_late(i)==1
                weak_strong_both_1(3)=weak_strong_both_1(3)+1;
            end
            if frM.zscore_group_EandMmerge(i)==10 & frM.zscore_group_late(i)==10
                weak_strong_both_1(4)=weak_strong_both_1(4)+1;
            end
        end

        if (frM.zscore_group_EandMmerge(i)==1 | frM.zscore_group_EandMmerge(i)==10 ) & ( frM.zscore_group_late(i)==0)
            CategoryM(2)=CategoryM(2)+1;
            a=1;
            if frM.zscore_group_EandMmerge(i)==1
                weak_strong_both_2(1)=weak_strong_both_2(1)+1;
            end
            if frM.zscore_group_EandMmerge(i)==10
                weak_strong_both_2(2)=weak_strong_both_2(2)+1;
            end
        end

        if (frM.zscore_group_EandMmerge(i)== 0) & (frM.zscore_group_late(i)==1 | frM.zscore_group_late(i)==10)
            CategoryM(3)=CategoryM(3)+1;
            a=1;
            if frM.zscore_group_late(i)==1
                weak_strong_both_3(1)=weak_strong_both_3(1)+1;
            end
            if frM.zscore_group_late(i)==10
                weak_strong_both_3(2)=weak_strong_both_3(2)+1;
            end
        end


        if (frM.zscore_group_EandMmerge(i)==2 | frM.zscore_group_EandMmerge(i)==20 ) & ( frM.zscore_group_late(i)==2 | frM.zscore_group_late(i)==20)
            CategoryM(4)=CategoryM(4)+1;
            a=1;
        end

        if (frM.zscore_group_EandMmerge(i)==2 | frM.zscore_group_EandMmerge(i)==20 ) & ( frM.zscore_group_late(i)==0)
            CategoryM(5)=CategoryM(5)+1;
            a=1;
        end

        if (frM.zscore_group_EandMmerge(i)==0 ) & ( frM.zscore_group_late(i)==2 | frM.zscore_group_late(i)==20)
            CategoryM(6)=CategoryM(6)+1;
            a=1;
        end
            if (frM.zscore_group_EandMmerge(i)==0 ) & ( frM.zscore_group_late(i)==0)
            CategoryM(7)=CategoryM(7)+1;
            a=1;
            end
    else
        CategoryM(7)=CategoryM(7)+1;
    end
end

CategoryM(8)= (numel(IdxoptotaggedPtfa) - sum(CategoryM));

%% Control
%Photoactivated

CategoryY = [0 0 0 0 0 0 0 0];
weak_strong_bothY_1=[0 0 0 0];
weak_strong_bothY_2=[0 0];
weak_strong_bothY_3=[0 0];

for j = 1:length(IdxoptotaggedPtfaY)
    i=IdxoptotaggedPtfaY(j);
    
    if frY.Strong_firing(i) == 1

        if (frY.zscore_group_EandMmerge(i)==1 | frY.zscore_group_EandMmerge(i)==10 ) & ( frY.zscore_group_late(i)==1 | frY.zscore_group_late(i)==10)
            CategoryY(1)=CategoryY(1)+1;
            a=1;
            if frY.zscore_group_EandMmerge(i)==1 & frY.zscore_group_late(i)==1
                weak_strong_bothY_1(1)=weak_strong_bothY_1(1)+1;
            end
            if frY.zscore_group_EandMmerge(i)==1 & frY.zscore_group_late(i)==10
                weak_strong_bothY_1(2)=weak_strong_bothY_1(2)+1;
            end
            if frY.zscore_group_EandMmerge(i)==10 & frY.zscore_group_late(i)==1
                weak_strong_bothY_1(3)=weak_strong_bothY_1(3)+1;
            end
            if frY.zscore_group_EandMmerge(i)==10 & frY.zscore_group_late(i)==10
                weak_strong_bothY_1(4)=weak_strong_bothY_1(4)+1;
            end
        end

        if (frY.zscore_group_EandMmerge(i)==1 | frY.zscore_group_EandMmerge(i)==10 ) & ( frY.zscore_group_late(i)==0)
            CategoryY(2)=CategoryY(2)+1;
            a=1;
            if frY.zscore_group_EandMmerge(i)==1
                weak_strong_bothY_2(1)=weak_strong_bothY_2(1)+1;
            end
            if frY.zscore_group_EandMmerge(i)==10
                weak_strong_bothY_2(2)=weak_strong_bothY_2(2)+1;
            end
        end

        if (frY.zscore_group_EandMmerge(i)==0 ) & ( frY.zscore_group_late(i)==1 | frY.zscore_group_late(i)==10)
            CategoryY(3)=CategoryY(3)+1;
            a=1;
            if frY.zscore_group_late(i)==1
                weak_strong_bothY_3(1)=weak_strong_bothY_3(1)+1;
            end
            if frY.zscore_group_late(i)==10
                weak_strong_bothY_3(2)=weak_strong_bothY_3(2)+1;
            end
        end

        if (frY.zscore_group_EandMmerge(i)==2 | frY.zscore_group_EandMmerge(i)==20 ) & ( frY.zscore_group_late(i)==2 | frY.zscore_group_late(i)==20)
            CategoryY(4)=CategoryY(4)+1;
            a=1;
        end

        if (frY.zscore_group_EandMmerge(i)==2 | frY.zscore_group_EandMmerge(i)==20 ) & ( frY.zscore_group_late(i)==0)
            CategoryY(5)=CategoryY(5)+1;
            a=1;
        end

        if (frY.zscore_group_EandMmerge(i)==0 ) & ( frY.zscore_group_late(i)==2 | frY.zscore_group_late(i)==20)
            CategoryY(6)=CategoryY(6)+1;
            a=1;
        end
        if (frY.zscore_group_EandMmerge(i)==0 ) & ( frY.zscore_group_late(i)==0)
            CategoryY(7)=CategoryY(7)+1;
            a=1;
        end
    else
        CategoryY(7)=CategoryY(7)+1;
    end
end


CategoryY(8)= (numel(IdxoptotaggedPtfaY)- sum(CategoryY));
StatPtf1a=[];


% cat 1: EL, 2 = E, 3= L, 4= ELd, 5= Ed, 6=Ld, 7= NP, 8= Other
StatPtf1a(1) = chisquare(sum(CategoryY(2)), sum(CategoryY(1:7)), sum(CategoryM(2)), sum(CategoryM(1:7)));
StatPtf1a(2) = chisquare(sum(CategoryY(1)), sum(CategoryY(1:7)), sum(CategoryM(1)), sum(CategoryM(1:7)));
StatPtf1a(3) = chisquare(sum(CategoryY(3)), sum(CategoryY(1:7)), sum(CategoryM(3)), sum(CategoryM(1:7)));
StatPtf1a(4) = chisquare(sum(CategoryY(5)), sum(CategoryY(1:7)), sum(CategoryM(5)), sum(CategoryM(1:7)));
StatPtf1a(5) = chisquare(sum(CategoryY(4)), sum(CategoryY(1:7)), sum(CategoryM(4)), sum(CategoryM(1:7)));
StatPtf1a(6) = chisquare(sum(CategoryY(6)), sum(CategoryY(1:7)), sum(CategoryM(6)), sum(CategoryM(1:7)));
StatPtf1a(7) = chisquare(sum(CategoryY(7)), sum(CategoryY(1:7)), sum(CategoryM(7)), sum(CategoryM(1:7)));
StatPtf1a(2,3) = chisquare(numel(find(CategoryY <=3)), sum(CategoryY(1:7)), numel(find(CategoryM <=3)), sum(CategoryM(1:7)));
StatPtf1a(2,6) = chisquare(numel(find(CategoryY >3 & CategoryY <7)), sum(CategoryY(1:7)), numel(find(CategoryM >3 & CategoryM <7)), sum(CategoryM(1:7)));


%% Spider chart


figure
PhasesM = [sum(CategoryM(1:3))/sum(CategoryM(1:7)) sum(CategoryM(4:6))/sum(CategoryM(1:7)) sum(CategoryM(7))/sum(CategoryM(1:7))];   
PhasesY = [sum(CategoryY(1:3))/sum(CategoryY(1:7)) sum(CategoryY(4:6))/sum(CategoryY(1:7)) sum(CategoryY(7))/sum(CategoryY(1:7))];   

Pphase = [PhasesM; PhasesY];

% Spider plot
spider_plot(Pphase, 'AxesLabels', {'Increase', 'Decrease', 'Non Changing'}, 'AxesLimits', [0, 0, 0; 0.8, 0.8, 0.8],'AxesOffset',0);
legend('Learner', 'Control',  'Location', 'southoutside');
title('Ptf1aON')

