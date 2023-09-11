%% Characterisation of the movement
%%

clc
clear all
close all
load ThrshDay.mat

ExpNb=[1];
DayAnalysed=1;
Ndpair{ExpNb(1)}=[1 2];

y=[];
for DayAnalysed=DayAnalysed
    
    iter=0;
    for i=ExpNb
        
        for j=Ndpair{i}
            
            iter=iter+1;
            
            filetoimporteXM=['/Users/simonlavaud/Desktop/Github_Kinematics/ConditioningExample/DataExp' num2str(i) '/AnalyseDay' num2str(DayAnalysed) '/pair' num2str(j) '/XmasterD' num2str(DayAnalysed) 'P' num2str(j) '.mat'];
            arrayXM=struct2array(load('-mat', filetoimporteXM));
            XM(1:length(arrayXM),1:5,iter)=arrayXM;
            
            filetoimporteYM=['/Users/simonlavaud/Desktop/Github_Kinematics/ConditioningExample/DataExp' num2str(i) '/AnalyseDay' num2str(DayAnalysed) '/pair' num2str(j) '/YmasterD' num2str(DayAnalysed) 'P' num2str(j) '.mat'];
            arrayYM=struct2array(load('-mat', filetoimporteYM));
            YM(1:length(arrayYM),1:5,iter)=arrayYM;
            
            filetoimporteZM=['/Users/simonlavaud/Desktop/Github_Kinematics/ConditioningExample/DataExp' num2str(i) '/AnalyseDay' num2str(DayAnalysed) '/pair' num2str(j) '/ZmasterD' num2str(DayAnalysed) 'P' num2str(j) '.mat'];
            arrayZM=struct2array(load('-mat', filetoimporteZM));
            ZM(1:length(arrayZM),1:5,iter)=arrayZM-ThrshDay(DayAnalysed,j*2-1,i)+3;
            
            
            filetoimporteXY=['/Users/simonlavaud/Desktop/Github_Kinematics/ConditioningExample/DataExp' num2str(i) '/AnalyseDay' num2str(DayAnalysed) '/pair' num2str(j) '/XyokedD' num2str(DayAnalysed) 'P' num2str(j) '.mat'];
            arrayXY=struct2array(load('-mat', filetoimporteXY));
            XY(1:length(arrayXY),1:5,iter)=arrayXY;
            
            filetoimporteYY=['/Users/simonlavaud/Desktop/Github_Kinematics/ConditioningExample/DataExp' num2str(i) '/AnalyseDay' num2str(DayAnalysed) '/pair' num2str(j) '/YyokedD' num2str(DayAnalysed) 'P' num2str(j) '.mat'];
            arrayYY=struct2array(load('-mat', filetoimporteYY));
            YY(1:length(arrayYY),1:5,iter)=arrayYY;
            
            filetoimporteZY=['/Users/simonlavaud/Desktop/Github_Kinematics/ConditioningExample/DataExp' num2str(i) '/AnalyseDay' num2str(DayAnalysed) '/pair' num2str(j) '/ZyokedD' num2str(DayAnalysed) 'P' num2str(j) '.mat'];
            arrayZY=struct2array(load('-mat', filetoimporteZY));
            ZY(1:length(arrayZY),1:5,iter)=arrayZY-ThrshDay(DayAnalysed,j*2,i)+3;
            
        end
        
        T(iter)=ThrshDay(DayAnalysed,j*2-1,i);
        TY(iter)=ThrshDay(DayAnalysed,j*2,i);
        
        length_exp(iter)=length(arrayXM);
        
    end
    
end

for i=1:iter
    shockinter=0;
    for k=1:18000
        if ZM(k,1,i)<3
            shockinter=shockinter+1;
        end
    end
    shock_received(i,DayAnalysed)=shockinter;
end

%Down-sample for vizualisation
clearvars y yY
for i=1:iter
    y(1:600,i) = downsample(ZM(1:18000,1,i),30);
    yY(1:600,i) = downsample(ZY(1:18000,1,i),30);
end

clearvars median_ZM median_ZY up_med down_med up_medY down_medY
mean_ZM=mean(y,2); 
mean_ZY=mean(yY,2);

for i=1:600
x = y(i,:); 
xY = yY(i,:);

up_med(i)=mean_ZM(i)+1*std(x);
down_med(i)=mean_ZM(i)-1*std(x);

up_medY(i)=mean_ZY(i)+1*std(xY);
down_medY(i)=mean_ZY(i)-1*std(xY);
end


% percentage time spend below threshold
clearvars perc_below_thres_master perc_below_thres_yoked
perc_below_thres_master(:)=sum(ZM(1:18000,1,:)<3)/18000;
perc_below_thres_yoked(:)=sum(ZY(1:18000,1,:)<3)/18000;



%%  INDIVIDUAL PLOT FOR LEARNER - RAW 30Hz
close all
figure(1)
for i= 1:iter
    subplot(iter,1,i)
    plot(ZM(1:18000,1,i))
    line([0 18000],[3 3],'Color','r')
    if i==1
    title('Individual foot position - Learner')
    end
end

close all
figure(2)
for i= 1:iter
    subplot(iter,1,i)
    plot(ZY(1:18000,1,i))
    line([0 18000],[3 3],'Color','r')
    if i==1
    title('Individual foot position - Control')
    end
end



%% MEAN OF INDIVIDUAL PLOT - DOWNSAMPLED TO 1HZ

figure
subplot(2,1,1)
title('Learner')
hold on
plot(mean_ZM,'r','LineWidth',2)
plot(up_med,'-r')
plot(down_med,'-r')
line([0 600],[3 3],'LineWidth',2)
xticks([0 150 300 450 600])
xticklabels({'0','150','300','450','600'})
xlabel('Time (s)')
axis([0 600 -10 16])
ylabel('Foot position (mm)')
grid on
grid minor
hold off

subplot(2,1,2)
hold on
title('Control')
plot(mean_ZY,'b','LineWidth',2)
plot(up_medY,'-b')
plot(down_medY,'-b')
line([0 600],[3 3],'LineWidth',2)
xticks([0 150 300 450 600])
xticklabels({'0','150','300','450','600'})
xlabel('Time (s)')
ylabel('Foot position (mm)')
axis([0 600 -10 16])
grid on
grid minor
hold off

%% FIRST X SECOND MASTER VS YOKED

mouse=2; %Mouse to selec


start=1; %Start time
time=10; %Time in second
figure
subplot(2,1,1)
plot(ZM(start:time*30+start,1,mouse))
for i=start:numel(ZM(1:time*30,1))+start
    if ZM(i,1,mouse)<3
        line([i-start+1 i-start+1],[-2 2],'LineWidth',1,'Color',[0.5 0.5 0.5])
    end
end
line([0 time*30],[3 3],'LineWidth',2)
title('first selected seconds - Learner')
axis([0 time*30 -5 15 ])

subplot(2,1,2)
plot(ZY(start:time*30+start,1,mouse))
for i=start:numel(ZM(1:time*30,1))+start
    if ZM(i,1,mouse)<3
        line([i-start+1 i-start+1],[-2 2],'LineWidth',1,'Color',[0.5 0.5 0.5])
    end
end
line([0 time*30],[3 3],'LineWidth',2)
title('first selected seconds - Control')
axis([0 time*30 -5 15 ])
