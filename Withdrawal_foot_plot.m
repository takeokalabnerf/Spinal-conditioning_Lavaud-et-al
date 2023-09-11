% % % %
%
% Foot position during the foot withdrawal
%

clear all
close all
w = warning ('off','all');

load ThrshDay.mat

% - - - - - - - 
% Choose a dataset
% - - - - - - -  

ExpNb=[1 0];
Ndays = [1];
Npairs=[1 2 0];

% - - - - - - - 
% Parameters
% - - - - - - -  
NbofPoint=1;
DivideBy=1;
nbframes=15;
a=0;
framed=20;
alpha=0.4;

initNd=1;
initExp=1;



while ExpNb(initExp)>0
    
    initPair=1;
    while Npairs(initExp,initPair)>0
        
        a=a+1
        
        XM4Amplitude=load(['/Users/simonlavaud/Desktop/Github_Kinematics/ConditioningExample/DataExp' num2str(ExpNb(initExp)) '/AnalyseDay' num2str(Ndays) '/pair' num2str(Npairs(initExp,initPair)) '/XmasterD' num2str(Ndays) 'P' num2str(Npairs(initExp,initPair)) '.mat']);
        YM4Amplitude=load(['/Users/simonlavaud/Desktop/Github_Kinematics/ConditioningExample/DataExp' num2str(ExpNb(initExp)) '/AnalyseDay' num2str(Ndays) '/pair' num2str(Npairs(initExp,initPair)) '/YmasterD' num2str(Ndays) 'P' num2str(Npairs(initExp,initPair)) '.mat']);
        ZM4Amplitude=load(['/Users/simonlavaud/Desktop/Github_Kinematics/ConditioningExample/DataExp' num2str(ExpNb(initExp)) '/AnalyseDay' num2str(Ndays) '/pair' num2str(Npairs(initExp,initPair)) '/ZmasterD' num2str(Ndays) 'P' num2str(Npairs(initExp,initPair)) '.mat']);
        
        XY4Amplitude=load(['/Users/simonlavaud/Desktop/Github_Kinematics/ConditioningExample/DataExp' num2str(ExpNb(initExp)) '/AnalyseDay' num2str(Ndays) '/pair' num2str(Npairs(initExp,initPair)) '/XyokedD' num2str(Ndays) 'P' num2str(Npairs(initExp,initPair)) '.mat']);
        YY4Amplitude=load(['/Users/simonlavaud/Desktop/Github_Kinematics/ConditioningExample/DataExp' num2str(ExpNb(initExp)) '/AnalyseDay' num2str(Ndays) '/pair' num2str(Npairs(initExp,initPair)) '/YyokedD' num2str(Ndays) 'P' num2str(Npairs(initExp,initPair)) '.mat']);
        ZY4Amplitude=load(['/Users/simonlavaud/Desktop/Github_Kinematics/ConditioningExample/DataExp' num2str(ExpNb(initExp)) '/AnalyseDay' num2str(Ndays) '/pair' num2str(Npairs(initExp,initPair)) '/ZyokedD' num2str(Ndays) 'P' num2str(Npairs(initExp,initPair)) '.mat']);
        
        
        XM4Amplitude=cell2mat(struct2cell(XM4Amplitude));
        YM4Amplitude=cell2mat(struct2cell(YM4Amplitude));
        ZM4Amplitude=cell2mat(struct2cell(ZM4Amplitude));
        
        XY4Amplitude=cell2mat(struct2cell(XY4Amplitude));
        YY4Amplitude=cell2mat(struct2cell(YY4Amplitude));
        ZY4Amplitude=cell2mat(struct2cell(ZY4Amplitude));
        
        
        T(a)=ThrshDay(Ndays,1+(2*(Npairs(initExp,initPair)-1)),ExpNb(initExp));
        TY(a)=ThrshDay(Ndays,(2*(Npairs(initExp,initPair))),ExpNb(initExp));
        L(a)=floor(length(ZM4Amplitude));
        
        ZM(1:length(XM4Amplitude),a,5)=ZM4Amplitude(1:length(XM4Amplitude),5)-T(a)+3;
        ZM(1:length(XM4Amplitude),a,4)=ZM4Amplitude(1:length(XM4Amplitude),4)-T(a)+3;
        ZM(1:length(XM4Amplitude),a,3)=ZM4Amplitude(1:length(XM4Amplitude),3)-T(a)+3;
        ZM(1:length(XM4Amplitude),a,2)=ZM4Amplitude(1:length(XM4Amplitude),2)-T(a)+3;
        ZM(1:length(XM4Amplitude),a,1)=ZM4Amplitude(1:length(XM4Amplitude),1)-T(a)+3;
        
        ZY(1:length(XM4Amplitude),a,5)=ZY4Amplitude(1:length(XM4Amplitude),5)-TY(a)+3;
        ZY(1:length(XM4Amplitude),a,4)=ZY4Amplitude(1:length(XM4Amplitude),4)-TY(a)+3;
        ZY(1:length(XM4Amplitude),a,3)=ZY4Amplitude(1:length(XM4Amplitude),3)-TY(a)+3;
        ZY(1:length(XM4Amplitude),a,2)=ZY4Amplitude(1:length(XM4Amplitude),2)-TY(a)+3;
        ZY(1:length(XM4Amplitude),a,1)=ZY4Amplitude(1:length(XM4Amplitude),1)-TY(a)+3;
        
        initPair=initPair+1;
    end
    initExp=initExp+1;
end


Itermice=a;
pulsenumber=1;

for i=1:Itermice
    for j=2:(L(i)-nbframes)
        
        
        if ZM(j-1,i,1)>3 && ZM(j,i,1)<3
            
            
            pulsesZMF(1:16,pulsenumber,i)=ZM(j:j+nbframes,i,1);
            pulsesZMA(1:16,pulsenumber,i)=ZM(j:j+nbframes,i,2);
            pulsesZMK(1:16,pulsenumber,i)=ZM(j:j+nbframes,i,3);
            pulsesZMH(1:16,pulsenumber,i)=ZM(j:j+nbframes,i,4);
            pulsesZMI(1:16,pulsenumber,i)=ZM(j:j+nbframes,i,5);
            
            pulsesZYF(1:16,pulsenumber,i)=ZY(j:j+nbframes,i,1);
            pulsesZYA(1:16,pulsenumber,i)=ZY(j:j+nbframes,i,2);
            pulsesZYK(1:16,pulsenumber,i)=ZY(j:j+nbframes,i,3);
            pulsesZYH(1:16,pulsenumber,i)=ZY(j:j+nbframes,i,4);
            pulsesZYI(1:16,pulsenumber,i)=ZY(j:j+nbframes,i,5);
            
            pulsenumber=pulsenumber+1;
            
        end
        
        
    end
    
    pulsenumberstart(i)=pulsenumber;
    pulsenumber=1;
    
end


% - - - - - - - 
% Upsample 10 times to 300Hz
% - - - - - - -  

x = 1:nbframes;
xq1 = .1:0.1:nbframes;  
sizexq1=size(xq1);

for i=1:Itermice
    for j=1:pulsenumberstart(i)-1
        
        
        yZMF=pulsesZMF(x,j,i);
        sZMF=spline(x,yZMF,xq1);
        pulsessortedZMFS(1:sizexq1(2),j,i)=sZMF(1,1:sizexq1(2));
        
        yZMA=pulsesZMA(x,j,i);
        sZMA=spline(x,yZMA,xq1);
        pulsessortedZMAS(1:sizexq1(2),j,i)=sZMA(1,1:sizexq1(2));
        
        yZMK=pulsesZMK(x,j,i);
        sZMK=spline(x,yZMK,xq1);
        pulsessortedZMKS(1:sizexq1(2),j,i)=sZMK(1,1:sizexq1(2));
        
        yZMH=pulsesZMH(x,j,i);
        sZMH=spline(x,yZMH,xq1);
        pulsessortedZMHS(1:sizexq1(2),j,i)=sZMH(1,1:sizexq1(2));
        
        yZMI=pulsesZMI(x,j,i);
        sZMI=spline(x,yZMI,xq1);
        pulsessortedZMIS(1:sizexq1(2),j,i)=sZMI(1,1:sizexq1(2));
        
        
        yZYF=pulsesZYF(x,j,i);
        sZYF=spline(x,yZYF,xq1);
        pulsessortedZYFS(1:sizexq1(2),j,i)=sZYF(1,1:sizexq1(2));
        
        yZYA=pulsesZYA(x,j,i);
        sZYA=spline(x,yZYA,xq1);
        pulsessortedZYAS(1:sizexq1(2),j,i)=sZYA(1,1:sizexq1(2));
        
        yZYK=pulsesZYK(x,j,i);
        sZYK=spline(x,yZYK,xq1);
        pulsessortedZYKS(1:sizexq1(2),j,i)=sZYK(1,1:sizexq1(2));
        
        yZYH=pulsesZYH(x,j,i);
        sZYH=spline(x,yZYH,xq1);
        pulsessortedZYHS(1:sizexq1(2),j,i)=sZYH(1,1:sizexq1(2));
        
        yZYI=pulsesZYI(x,j,i);
        sZYI=spline(x,yZYI,xq1);
        pulsessortedZYIS(1:sizexq1(2),j,i)=sZYI(1,1:sizexq1(2));
        
    end
    
end

% - - - - - - - 
% Start of the movement
% - - - - - - -  

StartCalc=squeeze(mean(pulsessortedZMFS(1,1:10,:)));
StartMasterF=mean(StartCalc);

StartCalc=squeeze(mean(pulsessortedZMAS(1,1:10,:)));
StartMasterA=mean(StartCalc);

StartCalc=squeeze(mean(pulsessortedZMKS(1,1:10,:)));
StartMasterK=mean(StartCalc);

StartCalc=squeeze(mean(pulsessortedZMHS(1,1:10,:)));
StartMasterH=mean(StartCalc);

StartCalc=squeeze(mean(pulsessortedZMIS(1,1:10,:)));
StartMasterI=mean(StartCalc);


StartCalc=squeeze(mean(pulsessortedZMFS(1,1:10,:)));
StartYokedF=mean(StartCalc);

StartCalc=squeeze(mean(pulsessortedZYAS(1,1:10,:)));
StartYokedA=mean(StartCalc);

StartCalc=squeeze(mean(pulsessortedZYKS(1,1:10,:)));
StartYokedK=mean(StartCalc);

StartCalc=squeeze(mean(pulsessortedZYHS(1,1:10,:)));
StartYokedH=mean(StartCalc);

StartCalc=squeeze(mean(pulsessortedZYIS(1,1:10,:)));
StartYokedI=mean(StartCalc);

% - - - - - - - 
% Alignement of the movement
% - - - - - - -  

[M I]=max(pulsessortedZMFS);
[MY IY]=max(pulsessortedZYFS);
b=1;
bY=1;

window_end=numel(xq1)/2;
% try to have 15 frames
window_start=numel(xq1)/10;
% Now 25 frames

for i=1:Itermice
    
    for j=1:pulsenumberstart(i)-1
        %before 6*nbframes
        if I(1,j,i)<window_end && I(1,j,i)>window_start
            pulsessortedZMF(1:91,b,i)=pulsessortedZMFS((I(1,j,i)-window_start):(I(1,j,i)+window_end),j,i)-pulsessortedZMFS(I(1,j,i)-window_start,j,i);%-ZtotCalib(i,1);%
            pulsessortedZMA(1:91,b,i)=pulsessortedZMAS((I(1,j,i)-window_start):(I(1,j,i)+window_end),j,i)-pulsessortedZMAS(I(1,j,i)-window_start,j,i);%-ZtotCalib(i,2);%
            pulsessortedZMK(1:91,b,i)=pulsessortedZMKS((I(1,j,i)-window_start):(I(1,j,i)+window_end),j,i)-pulsessortedZMKS(I(1,j,i)-window_start,j,i);%-ZtotCalib(i,3);%
            pulsessortedZMH(1:91,b,i)=pulsessortedZMHS((I(1,j,i)-window_start):(I(1,j,i)+window_end),j,i)-pulsessortedZMHS(I(1,j,i)-window_start,j,i);%-ZtotCalib(i,4);%
            pulsessortedZMI(1:91,b,i)=pulsessortedZMIS((I(1,j,i)-window_start):(I(1,j,i)+window_end),j,i)-pulsessortedZMIS(I(1,j,i)-window_start,j,i);%-ZtotCalib(i,5);%
            
            
            b=b+1;
        end
        
    end
    nbpulsessortedM(i)=b;
    b=1;
end


for i=1:Itermice
    
    for j=1:pulsenumberstart(i)-1
        
        if IY(1,j,i)<window_end && IY(1,j,i)>window_start
            pulsessortedZYF(1:91,bY,i)=pulsessortedZYFS((IY(1,j,i)-window_start):(IY(1,j,i)+window_end),j,i)-pulsessortedZYFS(IY(1,j,i)-window_start,j,i);;%-ZtotCalib(i,6);%
            pulsessortedZYA(1:91,bY,i)=pulsessortedZYAS((IY(1,j,i)-window_start):(IY(1,j,i)+window_end),j,i)-pulsessortedZYAS(IY(1,j,i)-window_start,j,i);;%-ZtotCalib(i,7);%
            pulsessortedZYK(1:91,bY,i)=pulsessortedZYKS((IY(1,j,i)-window_start):(IY(1,j,i)+window_end),j,i)-pulsessortedZYKS(IY(1,j,i)-window_start,j,i);;%-ZtotCalib(i,8);%
            pulsessortedZYH(1:91,bY,i)=pulsessortedZYHS((IY(1,j,i)-window_start):(IY(1,j,i)+window_end),j,i)-pulsessortedZYHS(IY(1,j,i)-window_start,j,i);;%-ZtotCalib(i,9);%
            pulsessortedZYI(1:91,bY,i)=pulsessortedZYIS((IY(1,j,i)-window_start):(IY(1,j,i)+window_end),j,i)-pulsessortedZYIS(IY(1,j,i)-window_start,j,i);;%-ZtotCalib(i,10);%
            
            
            bY=bY+1;
        end
        
    end
    nbpulsessortedY(i)=bY;
    bY=1;
end


% - - - - - - - 
% Remove pulse that are too low
% - - - - - - - 

clearvars pulsessortedZMF2 pulsessortedZMA2 pulsessortedZMK2 pulsessortedZMH2 pulsessortedZMI2
clearvars pulsessortedZYF2 pulsessortedZYA2 pulsessortedZYK2 pulsessortedZYH2 pulsessortedZYI2

b=1;
bY=1;

for i=1:Itermice
    for j=1:nbpulsessortedM(i)-1
        if max(pulsessortedZMF(:,j,i))>0.5
            
            pulsessortedZMF2(:,b,i)=pulsessortedZMF(:,j,i);
            pulsessortedZMA2(:,b,i)=pulsessortedZMA(:,j,i);
            pulsessortedZMK2(:,b,i)=pulsessortedZMK(:,j,i);
            pulsessortedZMH2(:,b,i)=pulsessortedZMH(:,j,i);
            pulsessortedZMI2(:,b,i)=pulsessortedZMI(:,j,i);
            
            b=b+1;
        end
        
    end
    
    nbpulsessorted2(i)=b;
    b=1;
    
end

for i=1:Itermice
    for j=1:nbpulsessortedY(i)-1
        if max(pulsessortedZYF(framed,j,i))>0.5
            
            pulsessortedZYF2(:,bY,i)=pulsessortedZYF(:,j,i);
            pulsessortedZYA2(:,bY,i)=pulsessortedZYA(:,j,i);
            pulsessortedZYK2(:,bY,i)=pulsessortedZYK(:,j,i);
            pulsessortedZYH2(:,bY,i)=pulsessortedZYH(:,j,i);
            pulsessortedZYI2(:,bY,i)=pulsessortedZYI(:,j,i);
            
            bY=bY+1;
            
        end
        
    end
    nbpulsessorted2Y(i)=bY;
    bY=1;
    
    
end

%%
close all
figure
hold on
meantestM=[];
meantestY=[];

number=1;
numbermean=1;
clearvars StoreAll

for test=1:Itermice

    clearvars StoreIndiv
    number=1;
    if nbpulsessorted2(test)>1
        for i=1:nbpulsessorted2(test)-1
            %plot(pulsessortedZMF2(1:90,i,test))
            StoreIndiv(1:90,number)=(pulsessortedZMF2(1:90,i,test));
            number=number+1;
        end
    
    
        StoreAll(1:90,numbermean)=mean(StoreIndiv(1:90,:),2);
        numbermean=numbermean+1
        plot(StoreAll)
    
    end
end

Withdraw_mean=mean(StoreAll(1:90,:),2);
Withdraw_std=std(StoreAll(1:90,:)'); %/sqrt(number-1);

    
plot(Withdraw_mean,'k','LineWidth',4)
plot(Withdraw_mean+Withdraw_std','k','LineWidth',2)
plot(Withdraw_mean-Withdraw_std','k','LineWidth',2)
hold off

xlim([0 45])
ylim([-1 5])

xlabel('time')
ylabel('mm')
xticks([0 15 30 45])
xticklabels({'0 ms','50 ms','100 ms','150 ms'})