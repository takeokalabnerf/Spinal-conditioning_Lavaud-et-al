%% - - - - Starting code:
tic
clear all
close all
w = warning ('off','all');


%% PARAMETERS INITIALIZATION
load ThrshDay.mat
ZM=zeros(3);

NbofPoint=1;
DivideBy=1;
nbframes=50;
a=0;
framed=20;
alpha=0.1;


%% - - - - Access to the different data

ExpNb=1;
Ndays = [1 0];
Npairs = [1 2 0];


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%% - - - EXTRACTING RAW DATA

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


initNd=1;
while (Ndays(initNd))>0
    
    initNp=1;
    while (Npairs(initNd,initNp))>0
        
        a=a+1;
        
        XM4Amplitude=load(['/Users/simonlavaud/Desktop/Github_Kinematics/ConditioningExample/DataExp' num2str(ExpNb) '/AnalyseDay' num2str(Ndays(initNd)) '/pair' num2str(Npairs(initNd,initNp)) '/XmasterD' num2str(Ndays(initNd)) 'P' num2str(Npairs(initNd,initNp)) '.mat']);
        YM4Amplitude=load(['/Users/simonlavaud/Desktop/Github_Kinematics/ConditioningExample/DataExp' num2str(ExpNb) '/AnalyseDay' num2str(Ndays(initNd)) '/pair' num2str(Npairs(initNd,initNp)) '/YmasterD' num2str(Ndays(initNd)) 'P' num2str(Npairs(initNd,initNp)) '.mat']);
        ZM4Amplitude=load(['/Users/simonlavaud/Desktop/Github_Kinematics/ConditioningExample/DataExp' num2str(ExpNb) '/AnalyseDay' num2str(Ndays(initNd)) '/pair' num2str(Npairs(initNd,initNp)) '/ZmasterD' num2str(Ndays(initNd)) 'P' num2str(Npairs(initNd,initNp)) '.mat']);
        
        XY4Amplitude=load(['/Users/simonlavaud/Desktop/Github_Kinematics/ConditioningExample/DataExp' num2str(ExpNb) '/AnalyseDay' num2str(Ndays(initNd)) '/pair' num2str(Npairs(initNd,initNp)) '/XyokedD' num2str(Ndays(initNd)) 'P' num2str(Npairs(initNd,initNp)) '.mat']);
        YY4Amplitude=load(['/Users/simonlavaud/Desktop/Github_Kinematics/ConditioningExample/DataExp' num2str(ExpNb) '/AnalyseDay' num2str(Ndays(initNd)) '/pair' num2str(Npairs(initNd,initNp)) '/YyokedD' num2str(Ndays(initNd)) 'P' num2str(Npairs(initNd,initNp)) '.mat']);
        ZY4Amplitude=load(['/Users/simonlavaud/Desktop/Github_Kinematics/ConditioningExample/DataExp' num2str(ExpNb) '/AnalyseDay' num2str(Ndays(initNd)) '/pair' num2str(Npairs(initNd,initNp)) '/ZyokedD' num2str(Ndays(initNd)) 'P' num2str(Npairs(initNd,initNp)) '.mat']);
        
        
        XM4Amplitude=cell2mat(struct2cell(XM4Amplitude));
        YM4Amplitude=cell2mat(struct2cell(YM4Amplitude));
        ZM4Amplitude=cell2mat(struct2cell(ZM4Amplitude));
        
        XY4Amplitude=cell2mat(struct2cell(XY4Amplitude));
        YY4Amplitude=cell2mat(struct2cell(YY4Amplitude));
        ZY4Amplitude=cell2mat(struct2cell(ZY4Amplitude));
        
        % Exctract lenght of the recording if longer than 600 sec
        L(a)=floor(length(ZM4Amplitude));
        
        
        XM4Amplitude=XM4Amplitude(1:(floor(length(XM4Amplitude)/(DivideBy))),:);
        YM4Amplitude=YM4Amplitude(1:(floor(length(YM4Amplitude)/(DivideBy))),:);
        ZM4Amplitude=ZM4Amplitude(1:(floor(length(ZM4Amplitude)/(DivideBy))),:);
        
        XY4Amplitude=XY4Amplitude(1:(floor(length(XY4Amplitude)/(DivideBy))),:);
        YY4Amplitude=YY4Amplitude(1:(floor(length(YY4Amplitude)/(DivideBy))),:);
        ZY4Amplitude=ZY4Amplitude(1:(floor(length(ZY4Amplitude)/(DivideBy))),:);
        
        XM(1:length(XM4Amplitude),a,5)=XM4Amplitude(1:length(XM4Amplitude),5);%-XM4Amplitude(1,5);%-XM4Amplitude(1,1)
        XM(1:length(XM4Amplitude),a,4)=XM4Amplitude(1:length(XM4Amplitude),4);%-XM4Amplitude(1,4)-XM(1:length(XM4Amplitude),a,5);
        XM(1:length(XM4Amplitude),a,3)=XM4Amplitude(1:length(XM4Amplitude),3);%-XM4Amplitude(1,3)-XM(1:length(XM4Amplitude),a,4);
        XM(1:length(XM4Amplitude),a,2)=XM4Amplitude(1:length(XM4Amplitude),2);%-XM4Amplitude(1,2)-XM(1:length(XM4Amplitude),a,3);
        XM(1:length(XM4Amplitude),a,1)=XM4Amplitude(1:length(XM4Amplitude),1);%-XM4Amplitude(1,1)-XM(1:length(XM4Amplitude),a,2);
        
        XY(1:length(XY4Amplitude),a,5)=XY4Amplitude(1:length(XY4Amplitude),5);%-XY4Amplitude(1,5);
        XY(1:length(XY4Amplitude),a,4)=XY4Amplitude(1:length(XY4Amplitude),4);%-XY4Amplitude(1,4)-XY(1:length(XM4Amplitude),a,5);
        XY(1:length(XY4Amplitude),a,3)=XY4Amplitude(1:length(XY4Amplitude),3);%-XY4Amplitude(1,3)-XY(1:length(XM4Amplitude),a,4);
        XY(1:length(XY4Amplitude),a,2)=XY4Amplitude(1:length(XY4Amplitude),2);%-XY4Amplitude(1,2)-XY(1:length(XM4Amplitude),a,3);
        XY(1:length(XY4Amplitude),a,1)=XY4Amplitude(1:length(XY4Amplitude),1);%-XY4Amplitude(1,1)-XY(1:length(XM4Amplitude),a,2);
        
        YM(1:length(XM4Amplitude),a,5)=YM4Amplitude(1:length(XM4Amplitude),5);%-YM4Amplitude(1,5);
        YM(1:length(XM4Amplitude),a,4)=YM4Amplitude(1:length(XM4Amplitude),4);%-YM4Amplitude(1,4)-YM(1:length(XM4Amplitude),a,5);
        YM(1:length(XM4Amplitude),a,3)=YM4Amplitude(1:length(XM4Amplitude),3);%-YM4Amplitude(1,3)-YM(1:length(XM4Amplitude),a,4);
        YM(1:length(XM4Amplitude),a,2)=YM4Amplitude(1:length(XM4Amplitude),2);%-YM4Amplitude(1,2)-YM(1:length(XM4Amplitude),a,3);
        YM(1:length(XM4Amplitude),a,1)=YM4Amplitude(1:length(XM4Amplitude),1);%-YM4Amplitude(1,1)-YM(1:length(XM4Amplitude),a,2);
        
        YY(1:length(XY4Amplitude),a,5)=YY4Amplitude(1:length(XY4Amplitude),5);%-YY4Amplitude(1,5);
        YY(1:length(XY4Amplitude),a,4)=YY4Amplitude(1:length(XY4Amplitude),4);%-YY4Amplitude(1,4)-YY(1:length(XM4Amplitude),a,5);
        YY(1:length(XY4Amplitude),a,3)=YY4Amplitude(1:length(XY4Amplitude),3);%-YY4Amplitude(1,3)-YY(1:length(XM4Amplitude),a,4);
        YY(1:length(XY4Amplitude),a,2)=YY4Amplitude(1:length(XY4Amplitude),2);%-YY4Amplitude(1,2)-YY(1:length(XM4Amplitude),a,3);
        YY(1:length(XY4Amplitude),a,1)=YY4Amplitude(1:length(XY4Amplitude),1);%-YY4Amplitude(1,1)-YY(1:length(XM4Amplitude),a,2);
        
        
        ZM(1:length(XM4Amplitude),a,5)=ZM4Amplitude(1:length(XM4Amplitude),5);%-ZM4Amplitude(1,5);
        ZM(1:length(XM4Amplitude),a,4)=ZM4Amplitude(1:length(XM4Amplitude),4);%-ZM4Amplitude(1,4)-ZM(1:length(XM4Amplitude),a,5);
        ZM(1:length(XM4Amplitude),a,3)=ZM4Amplitude(1:length(XM4Amplitude),3);%-ZM4Amplitude(1,3)-ZM(1:length(XM4Amplitude),a,4);
        ZM(1:length(XM4Amplitude),a,2)=ZM4Amplitude(1:length(XM4Amplitude),2);%-ZM4Amplitude(1,2)-ZM(1:length(XM4Amplitude),a,3);
        ZM(1:length(XM4Amplitude),a,1)=ZM4Amplitude(1:length(XM4Amplitude),1);%-ZM4Amplitude(1,1)-ZM(1:length(XM4Amplitude),a,2);
        
        ZY(1:length(XY4Amplitude),a,5)=ZY4Amplitude(1:length(XY4Amplitude),5);%-ZY4Amplitude(1,5);
        ZY(1:length(XY4Amplitude),a,4)=ZY4Amplitude(1:length(XY4Amplitude),4);%-ZY4Amplitude(1,4)-ZY(1:length(XM4Amplitude),a,5);
        ZY(1:length(XY4Amplitude),a,3)=ZY4Amplitude(1:length(XY4Amplitude),3);%-ZY4Amplitude(1,3)-ZY(1:length(XM4Amplitude),a,4);
        ZY(1:length(XY4Amplitude),a,2)=ZY4Amplitude(1:length(XY4Amplitude),2);%-ZY4Amplitude(1,2)-ZY(1:length(XM4Amplitude),a,3);
        ZY(1:length(XY4Amplitude),a,1)=ZY4Amplitude(1:length(XY4Amplitude),1);%-ZY4Amplitude(1,1)-ZY(1:length(XM4Amplitude),a,2);
        
        
        % - - - Extract threshsold position from database
        T(a)=ThrshDay(Ndays(initNd),1+(2*(Npairs(initNd,initNp)-1)),ExpNb);
        TY(a)=ThrshDay(Ndays(initNd),(2*(Npairs(initNd,initNp))),ExpNb);
        
        % - - - Extract the foot position for the last 60 sec
        XSec(a)=L(a)/60;
        
        LastXSecMZF(a)=mean(ZM4Amplitude((L(a)-XSec(a)):L(a),1))-T(a);%-mean(ZM4Amplitude(1:XSec(a),1));
        LastXSecYZF(a)=mean(ZY4Amplitude((L(a)-XSec(a)):L(a),1))-TY(a);%-mean(ZY4Amplitude(1:XSec(a),1));

        if LastXSecMZF(a)>0
            ControlFootPos(a)=1;
        else
            ControlFootPos(a)=0;
        end
        
        if LastXSecYZF(a)>0
            ControlFootPosY(a)=1;
        else
            ControlFootPosY(a)=0;
        end
                
        % - - - Extract mean position with a time bin of 1 min = 60 sec
        for i=1:10
            dec=floor(L(a)/10);
            inter=(1+dec*(i-1)):(dec*(i));
            
            for joint=1:5
                
                ZMmean10(a,i,joint)=mean(ZM4Amplitude(inter,joint))-T(a)+3;
                ZYmean10(a,i,joint)=mean(ZY4Amplitude(inter,joint))-TY(a)+3;
                XMmean10(a,i,joint)=mean(XM4Amplitude(inter,joint))-mean(XM4Amplitude(1:100,joint));
                XYmean10(a,i,joint)=mean(XY4Amplitude(inter,joint))-mean(XY4Amplitude(1:100,joint));
                YMmean10(a,i,joint)=mean(YM4Amplitude(inter,joint))-mean(YM4Amplitude(1:100,joint));
                YYmean10(a,i,joint)=mean(YY4Amplitude(inter,joint))-mean(YY4Amplitude(1:100,joint));
                
            end
            
        end
        
        
        % - - - Global mean of foot position during 600 sec
        ZMmean(Npairs(initNd,initNp),initNd)=mean(ZM4Amplitude(:,1))-T(a)+3;
        ZYmean(Npairs(initNd,initNp),initNd)=mean(ZY4Amplitude(:,1))-TY(a)+3;
        
        
        
        initNp=initNp+1;
        
    end
    
    initNd=initNd+1;
    
end
Itermice=a;

% - - - Exclude nul values

for i=1:Itermice
    for j=1:L(i)
        for k=1:5
            if ZM(j,i,k)==0
                ZM(j,i,k)=NaN;
            end
            if ZY(j,i,k)==0
                ZY(j,i,k)=NaN;
            end
            
            if XY(j,i,k)==0
                XY(j,i,k)=NaN;
            end
            if XM(j,i,k)==0
                XM(j,i,k)=NaN;
            end
            
            if YM(j,i,k)==0
                YM(j,i,k)=NaN;
            end
            if YY(j,i,k)==0
                YY(j,i,k)=NaN;
            end
        end
    end
end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% - - - COMPUTING KINEMATIC PARAMETERS - ADAPTIVE PHASE

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


% - - -  Mean

for a=1:Itermice
    
    for joint=1:5
        
        if joint==1
            MeanZM(a,joint)=mean(ZM(1:L(a),a,joint),'omitnan')-T(a)+3;
            MeanZY(a,joint)=mean(ZY(1:L(a),a,joint),'omitnan')-TY(a)+3;

        else
            MeanZM(a,joint)=mean(ZM(1:L(a),a,joint),'omitnan')-mean(ZM(1:(L(a)/10),a,joint),'omitnan');
            MeanZY(a,joint)=mean(ZY(1:L(a),a,joint),'omitnan')-mean(ZY(1:(L(a)/10),a,joint),'omitnan');

        end
        
        MeanXM(a,joint)=mean(XM(1:L(a),a,joint),'omitnan')-mean(XM(1:(L(a)/10),a,joint),'omitnan');
        MeanXY(a,joint)=mean(XY(1:L(a),a,joint),'omitnan')-mean(XY(1:(L(a)/10),a,joint),'omitnan');
 
        MeanYM(a,joint)=mean(YM(1:L(a),a,joint),'omitnan')-mean(YM(1:(L(a)/10),a,joint),'omitnan');
        MeanYY(a,joint)=mean(YY(1:L(a),a,joint),'omitnan')-mean(YY(1:(L(a)/10),a,joint),'omitnan');

    end
end


% - - -  Amplitude
% To avoid noise artifact, amplitude is sorted from low to high
%Delta_amplitude is amplitude_90% - amplitude_10%

for a=1:Itermice
    
    clearvars VXM VYM VZM VXY VYY VZY
    
    [VXM IXM]=sort(XM(1:L(a),:,:));
    [VYM IYM]=sort(YM(1:L(a),:,:));
    [VZM IZM]=sort(ZM(1:L(a),:,:));
    
    [VXY IXY]=sort(XY(1:L(a),:,:));
    [VYY IYY]=sort(YY(1:L(a),:,:));
    [VZY IZY]=sort(ZY(1:L(a),:,:));
    
    for joint=1:5
        
        MaxAmpXM(a,joint)=VXM(L(a)-floor(L(a)/10),a,joint);
        MinAmpXM(a,joint)=VXM(floor(L(a)/10),a,joint);
        AmpXM(a,joint)=MaxAmpXM(a,joint)-MinAmpXM(a,joint) ;
        
        MaxAmpYM(a,joint)=VYM(L(a)-floor(L(a)/10),a,joint);
        MinAmpYM(a,joint)=VYM(floor(L(a)/10),a,joint);
        AmpYM(a,joint)=MaxAmpYM(a,joint)-MinAmpYM(a,joint);
        
        MaxAmpZM(a,joint)=VZM(L(a)-floor(L(a)/10),a,joint);
        MinAmpZM(a,joint)=VZM(floor(L(a)/10),a,joint);
        AmpZM=MaxAmpZM-MinAmpZM;
        
        
        MaxAmpXY(a,joint)=VXY(L(a)-floor(L(a)/10),a,joint);
        MinAmpXY(a,joint)=VXY(floor(L(a)/10),a,joint);
        AmpXY=MaxAmpXY-MinAmpXY;
        
        MaxAmpYY(a,joint)=VYY(L(a)-floor(L(a)/10),a,joint);
        MinAmpYY(a,joint)=VYY(floor(L(a)/10),a,joint);
        AmpYY=MaxAmpYY-MinAmpYY;
        
        MaxAmpZY(a,joint)=VZY(L(a)-floor(L(a)/10),a,joint);
        MinAmpZY(a,joint)=VZY(floor(L(a)/10),a,joint);
        AmpZY=MaxAmpZY-MinAmpZY;
        
        
    end
    
    
end


% - - - Calculation of shock events and time below threshold

for a=1:Itermice
    
    for i=2:L(a)
        if ZM(i,a,1)<T(a)
            BeingS(i,a)=1;
        else
            BeingS(i,a)=0;
        end
        
        if ZY(i,a,1)<TY(a)
            BeingSY(i,a)=1;
        else
            BeingSY(i,a)=0;
        end
        
    end
    
    TimeBeingS(a)=sum(BeingS(:,a));
    TimeBeingSY(a)=sum(BeingSY(:,a));
    
    PercentageTimeBeingS(a)=TimeBeingS(a)/L(a);
    PercentageTimeBeingSY(a)=TimeBeingSY(a)/L(a);
    
end


% - - - Mean for 2D and mean + std for 3D angles

index=1;

for a=1:Itermice
    
    for i=2:L(a)
        
        A=[YM(i,a,1) ZM(i,a,1)];
        B=[YM(i,a,2) ZM(i,a,2)];
        C=[YM(i,a,3) ZM(i,a,3)];
        D=[YM(i,a,4) ZM(i,a,4)];
        E=[YM(i,a,5) ZM(i,a,5)];
        
        AY=[YY(i,a,1) ZY(i,a,1)];
        BY=[YY(i,a,2) ZY(i,a,2)];
        CY=[YY(i,a,3) ZY(i,a,3)];
        DY=[YY(i,a,4) ZY(i,a,4)];
        EY=[YY(i,a,5) ZY(i,a,5)];
        
        BA=A-B;
        BC=C-B;
        BAY=AY-BY;
        BCY=CY-BY;
        
        CB=B-C;
        CD=D-C;
        CBY=BY-CY;
        CDY=DY-CY;
        
        DC=C-D;
        DE=E-D;
        DCY=CY-DY;
        DEY=EY-DY;
        
        Angle2DYZ(i,a,1)=acos(dot(BA,BC)/(norm(BA)*norm(BC)));%*180/pi;
        Angle2DYZ(i,a,2)=acos(dot(CB,CD)/(norm(CB)*norm(CD)));%*180/pi;
        Angle2DYZ(i,a,3)=acos(dot(DC,DE)/(norm(DC)*norm(DE)));%*180/pi;
        
        Angle2DYYZ(i,a,1)=acos(dot(BAY,BCY)/(norm(BAY)*norm(BCY)));%*180/pi;
        Angle2DYYZ(i,a,2)=acos(dot(CBY,CDY)/(norm(CBY)*norm(CDY)));%*180/pi;
        Angle2DYYZ(i,a,3)=acos(dot(DCY,DEY)/(norm(DCY)*norm(DEY)));%*180/pi;
        
        
        
        
        A=[XM(i,a,1) YM(i,a,1)];
        B=[XM(i,a,2) YM(i,a,2)];
        C=[XM(i,a,3) YM(i,a,3)];
        D=[XM(i,a,4) YM(i,a,4)];
        E=[XM(i,a,5) YM(i,a,5)];
        
        AY=[XY(i,a,1) YY(i,a,1)];
        BY=[XY(i,a,2) YY(i,a,2)];
        CY=[XY(i,a,3) YY(i,a,3)];
        DY=[XY(i,a,4) YY(i,a,4)];
        EY=[XY(i,a,5) YY(i,a,5)];
        
        BA=A-B;
        BC=C-B;
        BAY=AY-BY;
        BCY=CY-BY;
        
        CB=B-C;
        CD=D-C;
        CBY=BY-CY;
        CDY=DY-CY;
        
        DC=C-D;
        DE=E-D;
        DCY=CY-DY;
        DEY=EY-DY;
        
        Angle2DXY(i,a,1)=acos(dot(BA,BC)/(norm(BA)*norm(BC)));%*180/pi;
        Angle2DXY(i,a,2)=acos(dot(CB,CD)/(norm(CB)*norm(CD)));%*180/pi;
        Angle2DXY(i,a,3)=real(acos(dot(DC,DE)/(norm(DC)*norm(DE))));%*180/pi);
        
        Angle2DYXY(i,a,1)=acos(dot(BAY,BCY)/(norm(BAY)*norm(BCY)));%*180/pi;
        Angle2DYXY(i,a,2)=acos(dot(CBY,CDY)/(norm(CBY)*norm(CDY)));%*180/pi;
        Angle2DYXY(i,a,3)=real(acos(dot(DCY,DEY)/(norm(DCY)*norm(DEY))));%*180/pi);
        
        
        
        A=[ZM(i,a,1) XM(i,a,1)];
        B=[ZM(i,a,2) XM(i,a,2)];
        C=[ZM(i,a,3) XM(i,a,3)];
        D=[ZM(i,a,4) XM(i,a,4)];
        E=[ZM(i,a,5) XM(i,a,5)];
        
        AY=[ZY(i,a,1) XY(i,a,1)];
        BY=[ZY(i,a,2) XY(i,a,2)];
        CY=[ZY(i,a,3) XY(i,a,3)];
        DY=[ZY(i,a,4) XY(i,a,4)];
        EY=[ZY(i,a,5) XY(i,a,5)];
        
        
        
        BA=A-B;
        BC=C-B;
        
        CB=B-C;
        CD=D-C;
        
        DC=C-D;
        DE=E-D;
        
        BAY=AY-BY;
        BCY=CY-BY;
        
        CBY=BY-CY;
        CDY=DY-CY;
        
        DCY=CY-DY;
        DEY=EY-DY;
        
        Angle2DZX(i,a,1)=acos(dot(BA,BC)/(norm(BA)*norm(BC)));%*180/pi;
        Angle2DZX(i,a,2)=acos(dot(CB,CD)/(norm(CB)*norm(CD)));%*180/pi;
        Angle2DZX(i,a,3)=acos(dot(DC,DE)/(norm(DC)*norm(DE)));%*180/pi;
        
        Angle2DYZX(i,a,1)=acos(dot(BAY,BCY)/(norm(BAY)*norm(BCY)));%*180/pi;
        Angle2DYZX(i,a,2)=acos(dot(CBY,CDY)/(norm(CBY)*norm(CDY)));%*180/pi;
        Angle2DYZX(i,a,3)=acos(dot(DCY,DEY)/(norm(DCY)*norm(DEY)));%*180/pi;
        
        
        
        A=[XM(i,a,1) YM(i,a,1) ZM(i,a,1)];
        B=[XM(i,a,2) YM(i,a,2) ZM(i,a,2)];
        C=[XM(i,a,3) YM(i,a,3) ZM(i,a,3)];
        D=[XM(i,a,4) YM(i,a,4) ZM(i,a,4)];
        E=[XM(i,a,5) YM(i,a,5) ZM(i,a,5)];
        
        AY=[XY(i,a,1) YY(i,a,1) ZY(i,a,1)];
        BY=[XY(i,a,2) YY(i,a,2) ZY(i,a,2)];
        CY=[XY(i,a,3) YY(i,a,3) ZY(i,a,3)];
        DY=[XY(i,a,4) YY(i,a,4) ZY(i,a,4)];
        EY=[XY(i,a,5) YY(i,a,5) ZY(i,a,5)];
        
        
        BA=A-B;
        BC=C-B;
        
        CB=B-C;
        CD=D-C;
        
        DC=C-D;
        DE=E-D;
        
        
        BAY=AY-BY;
        BCY=CY-BY;
        
        CBY=BY-CY;
        CDY=DY-CY;
        
        DCY=CY-DY;
        DEY=EY-DY;
        
        Angle3D(i,a,1)= acos(dot(BA,BC)/(norm(BA)*norm(BC)))*180/pi;
        Angle3D(i,a,2)= acos(dot(CB,CD)/(norm(CB)*norm(CD)))*180/pi;
        Angle3D(i,a,3)= acos(dot(DC,DE)/(norm(DC)*norm(DE)))*180/pi;
        
        Angle3DY(i,a,1)= acos(dot(BAY,BCY)/(norm(BAY)*norm(BCY)))*180/pi;
        Angle3DY(i,a,2)= acos(dot(CBY,CDY)/(norm(CBY)*norm(CDY)))*180/pi;
        Angle3DY(i,a,3)= acos(dot(DCY,DEY)/(norm(DCY)*norm(DEY)))*180/pi;
        
        
    end
    
    
    for angle=1:3
        MeanAngle2DYZ(a,angle)=mean(Angle2DYZ(1:L(a),a,angle),'omitnan');%-Angle2DYZ(2,a,angle);
        MeanAngle2DZX(a,angle)=mean(Angle2DZX(1:L(a),a,angle),'omitnan');%-Angle2DZX(2,a,angle);
        MeanAngle2DXY(a,angle)=mean(Angle2DXY(1:L(a),a,angle),'omitnan');%-Angle2DXY(2,a,angle);
        MeanAngle2DYYZ(a,angle)=mean(Angle2DYYZ(1:L(a),a,angle),'omitnan');%-Angle2DYYZ(2,a,angle);
        MeanAngle2DYZX(a,angle)=mean(Angle2DYZX(1:L(a),a,angle),'omitnan');%-Angle2DYZX(2,a,angle);
        MeanAngle2DYXY(a,angle)=mean(Angle2DYXY(1:L(a),a,angle),'omitnan');%-Angle2DYXY(2,a,angle);
        
        MeanAngle3D(a,angle)=mean(Angle3D(1:L(a),a,angle),'omitnan');
        MeanAngle3DY(a,angle)=mean(Angle3DY(1:L(a),a,angle),'omitnan');
        

    end
    
    
end


% - - -  Mean speed of the joints

for a=1:Itermice
    
    for joint=1:5
        
        speedZM(a,joint)=mean(diff(ZM(1:L(a),a,joint)));
        speedZY(a,joint)=mean(diff(ZY(1:L(a),a,joint)));
        
        speedXM(a,joint)=mean(diff(XM(1:L(a),a,joint)));
        speedXY(a,joint)=mean(diff(XY(1:L(a),a,joint)));
        
        speedYM(a,joint)=mean(diff(YM(1:L(a),a,joint)));
        speedYY(a,joint)=mean(diff(YY(1:L(a),a,joint)));
        
    end
    
end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% - - - COMPUTING KINEMATIC PARAMETERS - DYNAMIC PHASE

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


% - - - Extract the movement after cues 5 frame before and 15 frame after
% Working at 30Hz :     5 frames = 5/30s ~= 16ms
%                       15 frames = 15/30 = 500ms

% Add in a tab all the isolate joint traces for each mice

nbframes=15;
pulsenumber=1;

for i=1:Itermice
    for j=6:(L(i)-nbframes)
        
        if ZM(j-1,i,1)>T(i) && ZM(j,i,1)<T(i) %&& (ZM(j+5,i,1)-ZM(j,i,1))>0.1 && (ZY(j+5,i,1)-ZY(j,i,1))>0.1
            
            pulsesZMF(:,pulsenumber,i)=ZM(j-5:j+nbframes,i,1) - ZM(j,i,1);
            pulsesZMA(:,pulsenumber,i)=ZM(j-5:j+nbframes,i,2) - ZM(j,i,1);
            pulsesZMK(:,pulsenumber,i)=ZM(j-5:j+nbframes,i,3) - ZM(j,i,1);
            pulsesZMH(:,pulsenumber,i)=ZM(j-5:j+nbframes,i,4) - ZM(j,i,1);
            pulsesZMI(:,pulsenumber,i)=ZM(j-5:j+nbframes,i,5) - ZM(j,i,1);
            
            pulsesZYF(:,pulsenumber,i)=ZY(j-5:j+nbframes,i,1) - ZY(j,i,1);
            pulsesZYA(:,pulsenumber,i)=ZY(j-5:j+nbframes,i,2) - ZY(j,i,1);
            pulsesZYK(:,pulsenumber,i)=ZY(j-5:j+nbframes,i,3) - ZY(j,i,1);
            pulsesZYH(:,pulsenumber,i)=ZY(j-5:j+nbframes,i,4) - ZY(j,i,1);
            pulsesZYI(:,pulsenumber,i)=ZY(j-5:j+nbframes,i,5) - ZY(j,i,1);
            
            pulsesXMF(:,pulsenumber,i)=XM(j-5:j+nbframes,i,1) - XM(j,i,1);
            pulsesXMA(:,pulsenumber,i)=XM(j-5:j+nbframes,i,2) - XM(j,i,1);
            pulsesXMK(:,pulsenumber,i)=XM(j-5:j+nbframes,i,3) - XM(j,i,1);
            pulsesXMH(:,pulsenumber,i)=XM(j-5:j+nbframes,i,4) - XM(j,i,1);
            pulsesXMI(:,pulsenumber,i)=XM(j-5:j+nbframes,i,5) - XM(j,i,1);
            
            pulsesXYF(:,pulsenumber,i)=XY(j-5:j+nbframes,i,1) - XY(j,i,1);
            pulsesXYA(:,pulsenumber,i)=XY(j-5:j+nbframes,i,2) - XY(j,i,1);
            pulsesXYK(:,pulsenumber,i)=XY(j-5:j+nbframes,i,3) - XY(j,i,1);
            pulsesXYH(:,pulsenumber,i)=XY(j-5:j+nbframes,i,4) - XY(j,i,1);
            pulsesXYI(:,pulsenumber,i)=XY(j-5:j+nbframes,i,5) - XY(j,i,1);
            
            pulsesYMF(:,pulsenumber,i)=YM(j-5:j+nbframes,i,1) - YM(j,i,1);
            pulsesYMA(:,pulsenumber,i)=YM(j-5:j+nbframes,i,2) - YM(j,i,1);
            pulsesYMK(:,pulsenumber,i)=YM(j-5:j+nbframes,i,3) - YM(j,i,1);
            pulsesYMH(:,pulsenumber,i)=YM(j-5:j+nbframes,i,4) - YM(j,i,1);
            pulsesYMI(:,pulsenumber,i)=YM(j-5:j+nbframes,i,5) - YM(j,i,1);
            
            pulsesYYF(:,pulsenumber,i)=YY(j-5:j+nbframes,i,1) - YY(j,i,1);
            pulsesYYA(:,pulsenumber,i)=YY(j-5:j+nbframes,i,2) - YY(j,i,1);
            pulsesYYK(:,pulsenumber,i)=YY(j-5:j+nbframes,i,3) - YY(j,i,1);
            pulsesYYH(:,pulsenumber,i)=YY(j-5:j+nbframes,i,4) - YY(j,i,1);
            pulsesYYI(:,pulsenumber,i)=YY(j-5:j+nbframes,i,5) - YY(j,i,1);
            
            
            AnglePulse3D_Ankle(:,pulsenumber,i)=(Angle3D(j-5:j+nbframes,i,1));
            AnglePulse3D_Knee(:,pulsenumber,i)=(Angle3D(j-5:j+nbframes,i,2));
            AnglePulse3D_Hip(:,pulsenumber,i)=(Angle3D(j-5:j+nbframes,i,3));
            
            AnglePulse3DY_Ankle(:,pulsenumber,i)=(Angle3DY(j-5:j+nbframes,i,1));
            AnglePulse3DY_Knee(:,pulsenumber,i)=(Angle3DY(j-5:j+nbframes,i,2));
            AnglePulse3DY_Hip(:,pulsenumber,i)=(Angle3DY(j-5:j+nbframes,i,3));
            
            pulsenumber=pulsenumber+1;
            
        end
        
        
    end
    
    pulsenumberstart(i)=pulsenumber;
    pulsenumber=1;
    
end


PXMF=squeeze(mean(pulsesXMF(:,:,:),2));
PXMA=squeeze(mean(pulsesXMA(:,:,:),2));
PXMK=squeeze(mean(pulsesXMK(:,:,:),2));
PXMH=squeeze(mean(pulsesXMH(:,:,:),2));
PXMC=squeeze(mean(pulsesXMI(:,:,:),2));

PYMF=squeeze(mean(pulsesYMF(:,:,:),2));
PYMA=squeeze(mean(pulsesYMA(:,:,:),2));
PYMK=squeeze(mean(pulsesYMK(:,:,:),2));
PYMH=squeeze(mean(pulsesYMH(:,:,:),2));
PYMC=squeeze(mean(pulsesYMF(:,:,:),2));

PZMF=squeeze(mean(pulsesZMF(:,:,:),2));
PZMA=squeeze(mean(pulsesZMA(:,:,:),2));
PZMK=squeeze(mean(pulsesZMK(:,:,:),2));
PZMH=squeeze(mean(pulsesZMH(:,:,:),2));
PZMC=squeeze(mean(pulsesZMI(:,:,:),2));

PXMF=PXMF-PXMF(1,:);
PXMA=PXMA-PXMA(1,:);
PXMK=PXMK-PXMK(1,:);
PXMH=PXMH-PXMH(1,:);
PXMC=PXMC-PXMC(1,:);

PYMF=PYMF-PYMF(1,:);
PYMA=PYMA-PYMA(1,:);
PYMK=PYMK-PYMK(1,:);
PYMH=PYMH-PYMH(1,:);
PYMC=PYMC-PYMC(1,:);

PZMF=PZMF-PZMF(1,:);
PZMA=PZMA-PZMA(1,:);
PZMK=PZMK-PZMK(1,:);
PZMH=PZMH-PZMH(1,:);
PZMC=PZMC-PZMC(1,:);

% % % % % % % % % % % % % % % % % %

PXYF=squeeze(mean(pulsesXYF(:,:,:),2));
PXYA=squeeze(mean(pulsesXYA(:,:,:),2));
PXYK=squeeze(mean(pulsesXYK(:,:,:),2));
PXYH=squeeze(mean(pulsesXYH(:,:,:),2));
PXYC=squeeze(mean(pulsesXYI(:,:,:),2));

PYYF=squeeze(mean(pulsesYYF(:,:,:),2));
PYYA=squeeze(mean(pulsesYYA(:,:,:),2));
PYYK=squeeze(mean(pulsesYYK(:,:,:),2));
PYYH=squeeze(mean(pulsesYYH(:,:,:),2));
PYYC=squeeze(mean(pulsesYYI(:,:,:),2));

PZYF=squeeze(mean(pulsesZYF(:,:,:),2));
PZYA=squeeze(mean(pulsesZYA(:,:,:),2));
PZYK=squeeze(mean(pulsesZYH(:,:,:),2));
PZYH=squeeze(mean(pulsesZYH(:,:,:),2));
PZYC=squeeze(mean(pulsesZYI(:,:,:),2));

PXYF=PXYF-PXYF(1,:);
PXYA=PXYA-PXYA(1,:);
PXYK=PXYK-PXYK(1,:);
PXYH=PXYH-PXYH(1,:);
PXYC=PXYC-PXYC(1,:);

PYYF=PYYF-PYYF(1,:);
PYYA=PYYA-PYYA(1,:);
PYYK=PYYK-PYYK(1,:);
PYYH=PYYH-PYYH(1,:);
PYYC=PYYC-PYYC(1,:);

PZYF=PZYF-PZYF(1,:);
PZYA=PZYA-PZYA(1,:);
PZYK=PZYK-PZYK(1,:);
PZYH=PZYH-PZYH(1,:);
PZYC=PZYC-PZYC(1,:);

x = 1:15;
xq1 = 0.1:0.1:15;
xq2 = 1:0.1:14.9;


% - - - - - - - Computde the 3D angle

for i=1:Itermice
    for j=1:nbframes+1
        
        
        meanAngle3D_Ankle(j,i)=mean(AnglePulse3D_Ankle(j,1:pulsenumberstart(i)-1,i));
        meanAngle3D_Knee(j,i)=mean(AnglePulse3D_Knee(j,1:pulsenumberstart(i)-1,i));
        meanAngle3D_Hip(j,i)=mean(AnglePulse3D_Hip(j,1:pulsenumberstart(i)-1,i));
        
        meanAngle3DY_Ankle(j,i)=mean(AnglePulse3DY_Ankle(j,1:pulsenumberstart(i)-1,i));
        meanAngle3DY_Knee(j,i)=mean(AnglePulse3DY_Knee(j,1:pulsenumberstart(i)-1,i));
        meanAngle3DY_Hip(j,i)=mean(AnglePulse3DY_Hip(j,1:pulsenumberstart(i)-1,i));
        
    end
end



% - - - - Compute the different parameter for Dynamic phase
% Here the mean

for i=1:Itermice
    

        PXMF_M(i)=mean(PXMF(1:end,i));
        PXYF_M(i)=mean(PXYF(1:end,i));
        PYMF_M(i)=mean(PYMF(1:end,i));
        PYYF_M(i)=mean(PYYF(1:end,i));
        PZMF_M(i)=mean(PZMF(1:end,i));
        PZYF_M(i)=mean(PZYF(1:end,i));
        
        PXMA_M(i)=mean(PXMA(1:end,i));
        PXYA_M(i)=mean(PXYA(1:end,i));
        PYMA_M(i)=mean(PYMA(1:end,i));
        PYYA_M(i)=mean(PYYA(1:end,i));
        PZMA_M(i)=mean(PZMA(1:end,i));
        PZYA_M(i)=mean(PZYA(1:end,i));
        
        PXMK_M(i)=mean(PXMK(1:end,i));
        PXYK_M(i)=mean(PXYK(1:end,i));
        PYMK_M(i)=mean(PYMK(1:end,i));
        PYYK_M(i)=mean(PYYK(1:end,i));
        PZMK_M(i)=mean(PZMK(1:end,i));
        PZYK_M(i)=mean(PZYK(1:end,i));
        
        PXMH_M(i)=mean(PXMH(1:end,i));
        PXYH_M(i)=mean(PXYH(1:end,i));
        PYMH_M(i)=mean(PYMH(1:end,i));
        PYYH_M(i)=mean(PYYH(1:end,i));
        PZMH_M(i)=mean(PZMH(1:end,i));
        PZYH_M(i)=mean(PZYH(1:end,i));
        
        PXMC_M(i)=mean(PXMC(1:end,i));
        PXYC_M(i)=mean(PXYC(1:end,i));
        PYMC_M(i)=mean(PYMC(1:end,i));
        PYYC_M(i)=mean(PYYC(1:end,i));
        PZMC_M(i)=mean(PZMC(1:end,i));
        PZYC_M(i)=mean(PZYC(1:end,i));
        
    
end


% - - - - Area under the curve : Energy of the movement

for i=1:Itermice
    
    A_PXMF(i)=trapz(PXMF(:,i))-mean(PXMF(1:5,i));
    A_PXYF(i)=trapz(PXYF(:,i))-mean(PXYF(1:5,i));
    A_PYMF(i)=trapz(PYMF(:,i))-mean(PYMF(1:5,i));
    A_PYYF(i)=trapz(PYYF(:,i))-mean(PYYF(1:5,i));
    A_PZMF(i)=trapz(PZMF(:,i))-mean(PZMF(1:5,i));
    A_PZYF(i)=trapz(PZYF(:,i))-mean(PZYF(1:5,i));
    
    A_PXMA(i)=trapz(PXMA(:,i))-mean(PXMA(1:5,i));
    A_PXYA(i)=trapz(PXYA(:,i))-mean(PXYA(1:5,i));
    A_PYMA(i)=trapz(PYMA(:,i))-mean(PYMA(1:5,i));
    A_PYYA(i)=trapz(PYYA(:,i))-mean(PYYA(1:5,i));
    A_PZMA(i)=trapz(PZMA(:,i))-mean(PZMA(1:5,i));
    A_PZYA(i)=trapz(PZYA(:,i))-mean(PZYA(1:5,i));
    
    A_PXMK(i)=trapz(PXMK(:,i))-mean(PXMK(1:5,i));
    A_PXYK(i)=trapz(PXYK(:,i))-mean(PXYK(1:5,i));
    A_PYMK(i)=trapz(PYMK(:,i))-mean(PYMK(1:5,i));
    A_PYYK(i)=trapz(PYYK(:,i))-mean(PYYK(1:5,i));
    A_PZMK(i)=trapz(PZMK(:,i))-mean(PZMK(1:5,i));
    A_PZYK(i)=trapz(PZYK(:,i))-mean(PZYK(1:5,i));
    
    A_PXMH(i)=trapz(PXMH(:,i))-mean(PXMH(1:5,i));
    A_PXYH(i)=trapz(PXYH(:,i))-mean(PXYH(1:5,i));
    A_PYMH(i)=trapz(PYMH(:,i))-mean(PYMH(1:5,i));
    A_PYYH(i)=trapz(PYYH(:,i))-mean(PYYH(1:5,i));
    A_PZMH(i)=trapz(PZMH(:,i))-mean(PZMH(1:5,i));
    A_PZYH(i)=trapz(PZYH(:,i))-mean(PZYH(1:5,i));
    
    A_PXMC(i)=trapz(PXMC(:,i))-mean(PXMC(1:5,i));
    A_PXYC(i)=trapz(PXYC(:,i))-mean(PXYC(1:5,i));
    A_PYMC(i)=trapz(PYMC(:,i))-mean(PYMC(1:5,i));
    A_PYYC(i)=trapz(PYYC(:,i))-mean(PYYC(1:5,i));
    A_PZMC(i)=trapz(PZMC(:,i))-mean(PZMC(1:5,i));
    A_PZYC(i)=trapz(PZYC(:,i))-mean(PZYC(1:5,i));
end


% - - - - Speed and acceleration


for i=1:Itermice
   
        PZMF_speed(1:20,i)=diff(PZMF(:,i));
        PZYF_speed(1:20,i)=diff(PYMF(:,i));
        PXMF_speed(1:20,i)=diff(PXMF(:,i));
        PXYF_speed(1:20,i)=diff(PXMF(:,i));
        PYMF_speed(1:20,i)=diff(PYMF(:,i));
        PYYF_speed(1:20,i)=diff(PYMF(:,i));
        
        PZMA_speed(1:20,i)=diff(PZMA(:,i));
        PZYA_speed(1:20,i)=diff(PYMA(:,i));
        PXMA_speed(1:20,i)=diff(PXMA(:,i));
        PXYA_speed(1:20,i)=diff(PXMA(:,i));
        PYMA_speed(1:20,i)=diff(PYMA(:,i));
        PYYA_speed(1:20,i)=diff(PYMA(:,i));
        
        PZMK_speed(1:20,i)=diff(PZMK(:,i));
        PZYK_speed(1:20,i)=diff(PYMK(:,i));
        PXMK_speed(1:20,i)=diff(PXMK(:,i));
        PXYK_speed(1:20,i)=diff(PXMK(:,i));
        PYMK_speed(1:20,i)=diff(PYMK(:,i));
        PYYK_speed(1:20,i)=diff(PYMK(:,i));
        
        PZMC_speed(1:20,i)=diff(PZMC(:,i));
        PZYC_speed(1:20,i)=diff(PYMC(:,i));
        PXMC_speed(1:20,i)=diff(PXMC(:,i));
        PXYC_speed(1:20,i)=diff(PXMC(:,i));
        PYMC_speed(1:20,i)=diff(PYMC(:,i));
        PYYC_speed(1:20,i)=diff(PYMC(:,i));
        
        PZMH_speed(1:20,i)=diff(PZMH(:,i));
        PZYH_speed(1:20,i)=diff(PYMH(:,i));
        PXMH_speed(1:20,i)=diff(PXMH(:,i));
        PXYH_speed(1:20,i)=diff(PXMH(:,i));
        PYMH_speed(1:20,i)=diff(PYMH(:,i));
        PYYH_speed(1:20,i)=diff(PYMH(:,i));
        
        
        
        PZMF_acc(1:19,i)=diff(PZMF_speed(:,i));
        PZYF_acc(1:19,i)=diff(PZYF_speed(:,i));
        PXMF_acc(1:19,i)=diff(PXMF_speed(:,i));
        PXYF_acc(1:19,i)=diff(PXYF_speed(:,i));
        PYMF_acc(1:19,i)=diff(PYMF_speed(:,i));
        PYYF_acc(1:19,i)=diff(PYYF_speed(:,i));
        
        PZMA_acc(1:19,i)=diff(PZMA_speed(:,i));
        PZYA_acc(1:19,i)=diff(PZYA_speed(:,i));
        PXMA_acc(1:19,i)=diff(PXMA_speed(:,i));
        PXYA_acc(1:19,i)=diff(PXYA_speed(:,i));
        PYMA_acc(1:19,i)=diff(PYMA_speed(:,i));
        PYYA_acc(1:19,i)=diff(PYYA_speed(:,i));
        
        PZMK_acc(1:19,i)=diff(PZMK_speed(:,i));
        PZYK_acc(1:19,i)=diff(PZYK_speed(:,i));
        PXMK_acc(1:19,i)=diff(PXMK_speed(:,i));
        PXYK_acc(1:19,i)=diff(PXYK_speed(:,i));
        PYMK_acc(1:19,i)=diff(PYMK_speed(:,i));
        PYYK_acc(1:19,i)=diff(PYYK_speed(:,i));
        
        PZMC_acc(1:19,i)=diff(PZMC_speed(:,i));
        PZYC_acc(1:19,i)=diff(PZYC_speed(:,i));
        PXMC_acc(1:19,i)=diff(PXMC_speed(:,i));
        PXYC_acc(1:19,i)=diff(PXYC_speed(:,i));
        PYMC_acc(1:19,i)=diff(PYMC_speed(:,i));
        PYYC_acc(1:19,i)=diff(PYYC_speed(:,i));
        
        PZMH_acc(1:19,i)=diff(PZMH_speed(:,i));
        PZYH_acc(1:19,i)=diff(PZYH_speed(:,i));
        PXMH_acc(1:19,i)=diff(PXMH_speed(:,i));
        PXYH_acc(1:19,i)=diff(PXYH_speed(:,i));
        PYMH_acc(1:19,i)=diff(PYMH_speed(:,i));
        PYYH_acc(1:19,i)=diff(PYYH_speed(:,i));
        
    
end

% - - - Amplitude


for i=1:Itermice
    
    meanAmpPulseZMF(i)=max(PZMF(10:20,i))-mean(PZMF(1:10,i));
    meanAmpPulseZYF(i)=max(PZYF(10:20,i))-mean(PZYF(1:10,i));
    meanAmpPulseXMF(i)=max(PXMF(10:20,i))-mean(PXMF(1:10,i));
    meanAmpPulseXYF(i)=max(PXYF(10:20,i))-mean(PXYF(1:10,i));
    meanAmpPulseYMF(i)=max(PYMF(10:20,i))-mean(PYMF(1:10,i));
    meanAmpPulseYYF(i)=max(PYYF(10:20,i))-mean(PYYF(1:10,i));
    
    meanAmpPulseZMA(i)=max(PZMA(10:20,i))-mean(PZMA(1:10,i));
    meanAmpPulseZYA(i)=max(PZYA(10:20,i))-mean(PZYA(1:10,i));
    meanAmpPulseXMA(i)=max(PXMA(10:20,i))-mean(PXMA(1:10,i));
    meanAmpPulseXYA(i)=max(PXYA(10:20,i))-mean(PXYA(1:10,i));
    meanAmpPulseYMA(i)=max(PYMA(10:20,i))-mean(PYMA(1:10,i));
    meanAmpPulseYYA(i)=max(PYYA(10:20,i))-mean(PYYA(1:10,i));
    
    meanAmpPulseZMK(i)=max(PZMK(10:20,i))-mean(PZMK(1:10,i));
    meanAmpPulseZYK(i)=max(PZYK(10:20,i))-mean(PZYK(1:10,i));
    meanAmpPulseXMK(i)=max(PXMK(10:20,i))-mean(PXMK(1:10,i));
    meanAmpPulseXYK(i)=max(PXYK(10:20,i))-mean(PXYK(1:10,i));
    meanAmpPulseYMK(i)=max(PYMK(10:20,i))-mean(PYMK(1:10,i));
    meanAmpPulseYYK(i)=max(PYYK(10:20,i))-mean(PYYK(1:10,i));
    
    meanAmpPulseZMH(i)=max(PZMH(10:20,i))-mean(PZMH(1:10,i));
    meanAmpPulseZYH(i)=max(PZYH(10:20,i))-mean(PZYH(1:10,i));
    meanAmpPulseXMH(i)=max(PXMH(10:20,i))-mean(PXMH(1:10,i));
    meanAmpPulseXYH(i)=max(PXYH(10:20,i))-mean(PXYH(1:10,i));
    meanAmpPulseYMH(i)=max(PYMH(10:20,i))-mean(PYMH(1:10,i));
    meanAmpPulseYYH(i)=max(PYYH(10:20,i))-mean(PYYH(1:10,i));
    
    meanAmpPulseZMI(i)=max(PZMC(10:20,i))-mean(PZMC(1:10,i));
    meanAmpPulseZYI(i)=max(PZYC(10:20,i))-mean(PZYC(1:10,i));
    meanAmpPulseXMI(i)=max(PXMC(10:20,i))-mean(PXMC(1:10,i));
    meanAmpPulseXYI(i)=max(PXYC(10:20,i))-mean(PXYC(1:10,i));
    meanAmpPulseYMI(i)=max(PYMC(10:20,i))-mean(PYMC(1:10,i));
    meanAmpPulseYYI(i)=max(PYYC(10:20,i))-mean(PYYC(1:10,i));
    
    Angle3D_Ankle(i)=max(meanAngle3D_Ankle(5:end,i))-min(meanAngle3D_Ankle(5:end,i));
    Angle3D_Knee(i)=max(meanAngle3D_Knee(5:end,i))-min(meanAngle3D_Knee(5:end,i));
    Angle3D_Hip(i)=max(meanAngle3D_Hip(5:end,i))-min(meanAngle3D_Hip(5:end,i));
    
    Angle3DY_Ankle(i)=max(meanAngle3DY_Ankle(5:end,i))-min(meanAngle3DY_Ankle(5:end,i));
    Angle3DY_Knee(i)=max(meanAngle3DY_Knee(5:end,i))-min(meanAngle3DY_Knee(5:end,i));
    Angle3DY_Hip(i)=max(meanAngle3DY_Hip(5:end,i))-min(meanAngle3DY_Hip(5:end,i));
    
    
    meanDiffAmpPulseZMF(i)=mean(PZMF_speed(5:end,i));
    meanDiffAmpPulseZYF(i)=mean(PZYF_speed(5:end,i));
    meanDiffAmpPulseXMF(i)=mean(PXMF_speed(5:end,i));
    meanDiffAmpPulseXYF(i)=mean(PXYF_speed(5:end,i));
    meanDiffAmpPulseYMF(i)=mean(PYMF_speed(5:end,i));
    meanDiffAmpPulseYYF(i)=mean(PYYF_speed(5:end,i));
    
    meanDiffAmpPulseZMA(i)=mean(PZMA_speed(5:end,i));
    meanDiffAmpPulseZYA(i)=mean(PZYA_speed(5:end,i));
    meanDiffAmpPulseXMA(i)=mean(PXMA_speed(5:end,i));
    meanDiffAmpPulseXYA(i)=mean(PXYA_speed(5:end,i));
    meanDiffAmpPulseYMA(i)=mean(PYMA_speed(5:end,i));
    meanDiffAmpPulseYYA(i)=mean(PYYA_speed(5:end,i));
    
    meanDiffAmpPulseZMK(i)=mean(PZMK_speed(5:end,i));
    meanDiffAmpPulseZYK(i)=mean(PZYK_speed(5:end,i));
    meanDiffAmpPulseXMK(i)=mean(PXMK_speed(5:end,i));
    meanDiffAmpPulseXYK(i)=mean(PXYK_speed(5:end,i));
    meanDiffAmpPulseYMK(i)=mean(PYMK_speed(5:end,i));
    meanDiffAmpPulseYYK(i)=mean(PYYK_speed(5:end,i));
    
    meanDiffAmpPulseZMH(i)=mean(PZMH_speed(5:end,i));
    meanDiffAmpPulseZYH(i)=mean(PZYH_speed(5:end,i));
    meanDiffAmpPulseXMH(i)=mean(PXMH_speed(5:end,i));
    meanDiffAmpPulseXYH(i)=mean(PXYH_speed(5:end,i));
    meanDiffAmpPulseYMH(i)=mean(PYMH_speed(5:end,i));
    meanDiffAmpPulseYYH(i)=mean(PYYH_speed(5:end,i));
    
    meanDiffAmpPulseZMI(i)=mean(PZMC_speed(5:end,i));
    meanDiffAmpPulseZYI(i)=mean(PZYC_speed(5:end,i));
    meanDiffAmpPulseXMI(i)=mean(PXMC_speed(5:end,i));
    meanDiffAmpPulseXYI(i)=mean(PXYC_speed(5:end,i));
    meanDiffAmpPulseYMI(i)=mean(PYMC_speed(5:end,i));
    meanDiffAmpPulseYYI(i)=mean(PYYC_speed(5:end,i));
    
    
end



%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% - - - FILL THE TABLE FOR THE PCA

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


PCAtab=cell(Itermice*2,54);

PCAtab{1,1}='STATUS';
PCAtab{1,2}='DAY';
PCAtab{1,3}='PAIR';

Ndays(nnz(Ndays))=1;

% STATUS
for i=1:Itermice
    PCAtab{i+1,1}=1;
    PCAtab{i+1+Itermice,1}=0;
end

% DAY
a=2;
for i=1:Ndays(nnz(Ndays))
    for j=1:nnz(Npairs(i,:))
        
        PCAtab{a,2}=i;
        PCAtab{a+Itermice,2}=i;
        a=a+1;
        
    end
end

% PAIR
a=2;
for i=1:Ndays(nnz(Ndays))
    for j=1:nnz(Npairs(i,:))
        
        PCAtab{a,3}=Npairs(i,j);
        PCAtab{a+Itermice,3}=Npairs(i,j);
        a=a+1;
        
    end
end

Ndays(nnz(Ndays))=1;


% CONTROL
PCAtab{1,4}='Control';
a=2;
for i=1:Ndays(nnz(Ndays))
    for j=1:nnz(Npairs(i,:))
        
        PCAtab{a,4}=ControlFootPos(a-1);
        PCAtab{a+Itermice,4}=ControlFootPosY(a-1);
        a=a+1;
        
    end
end
Col=4;

% TIME SPEND UNDER THE THRESHOLD
for fold = 1
    Col=Col+1;
    Joint=1;
    PCAtab{1,Col}='Percentage_time_spend_under_the_threshold';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(PercentageTimeBeingS(a-1));
            PCAtab{a+Itermice,Col}=+(PercentageTimeBeingSY(a-1));
            a=a+1;
            
        end
    end
    
end

% MEAN GLOBAL MOVEMENT
for fold = 1
    Col=Col+1;
    Joint=1;
    PCAtab{1,Col}='MeanGlobal_XF';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(MeanXM(a-1,Joint));
            PCAtab{a+Itermice,Col}=+(MeanXY(a-1,Joint));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=2;
    PCAtab{1,Col}='MeanGlobal_XA';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(MeanXM(a-1,Joint));
            PCAtab{a+Itermice,Col}=+(MeanXY(a-1,Joint));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=3;
    PCAtab{1,Col}='MeanGlobal_XK';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(MeanXM(a-1,Joint));
            PCAtab{a+Itermice,Col}=+(MeanXY(a-1,Joint));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=4;
    PCAtab{1,Col}='MeanGlobal_XH';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(MeanXM(a-1,Joint));
            PCAtab{a+Itermice,Col}=+(MeanXY(a-1,Joint));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=5;
    PCAtab{1,Col}='MeanGlobal_XI';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(MeanXM(a-1,Joint));
            PCAtab{a+Itermice,Col}=+(MeanXY(a-1,Joint));
            a=a+1;
            
        end
    end
    
    % MEAN Y
    Col=Col+1;
    Joint=1;
    PCAtab{1,Col}='MeanGlobal_YF';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(MeanYM(a-1,Joint));
            PCAtab{a+Itermice,Col}=+(MeanYY(a-1,Joint));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=2;
    PCAtab{1,Col}='MeanGlobal_YA';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(MeanYM(a-1,Joint));
            PCAtab{a+Itermice,Col}=+(MeanYY(a-1,Joint));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=3;
    PCAtab{1,Col}='MeanGlobal_YK';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(MeanYM(a-1,Joint));
            PCAtab{a+Itermice,Col}=+(MeanYY(a-1,Joint));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=4;
    PCAtab{1,Col}='MeanGlobal_YH';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(MeanYM(a-1,Joint));
            PCAtab{a+Itermice,Col}=+(MeanYY(a-1,Joint));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=5;
    PCAtab{1,Col}='MeanGlobal_YI';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(MeanYM(a-1,Joint));
            PCAtab{a+Itermice,Col}=+(MeanYY(a-1,Joint));
            a=a+1;
            
        end
    end
    
    % MEAN Z
    Col=Col+1;
    Joint=1;
    PCAtab{1,Col}='MeanGlobal_ZF';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(MeanZM(a-1,Joint));
            PCAtab{a+Itermice,Col}=+(MeanZY(a-1,Joint));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=2;
    PCAtab{1,Col}='MeanGlobal_ZA';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(MeanZM(a-1,Joint));
            PCAtab{a+Itermice,Col}=+(MeanZY(a-1,Joint));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=3;
    PCAtab{1,Col}='MeanGlobal_ZK';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(MeanZM(a-1,Joint));
            PCAtab{a+Itermice,Col}=+(MeanZY(a-1,Joint));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=4;
    PCAtab{1,Col}='MeanGlobal_ZH';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(MeanZM(a-1,Joint));
            PCAtab{a+Itermice,Col}=+(MeanZY(a-1,Joint));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=5;
    PCAtab{1,Col}='MeanGlobal_ZI';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(MeanZM(a-1,Joint));
            PCAtab{a+Itermice,Col}=+(MeanZY(a-1,Joint));
            a=a+1;
            
        end
    end
end

% ANGLE 3D MEAN
for fold = 1
    Col=Col+1;
    Joint=1;
    PCAtab{1,Col}='Angle_3D_Ankle';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(MeanAngle3D(a-1,Joint));
            PCAtab{a+Itermice,Col}=+(MeanAngle3DY(a-1,Joint));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=2;
    PCAtab{1,Col}='Angle_3D_Knee';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(MeanAngle3D(a-1,Joint));
            PCAtab{a+Itermice,Col}=+(MeanAngle3DY(a-1,Joint));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=3;
    PCAtab{1,Col}='Angle_3D_Hip';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(MeanAngle3D(a-1,Joint));
            PCAtab{a+Itermice,Col}=+(MeanAngle3DY(a-1,Joint));
            a=a+1;
            
        end
    end
    
end

% AMP GLOBAL MOVEMENT
for fold = 1
    Col=Col+1;
    Joint=1;
    PCAtab{1,Col}='AmpGlobal_XF';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(AmpXM(a-1,Joint));
            PCAtab{a+Itermice,Col}=+(AmpXY(a-1,Joint));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=2;
    PCAtab{1,Col}='AmpGlobal_XA';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(AmpXM(a-1,Joint));
            PCAtab{a+Itermice,Col}=+(AmpXY(a-1,Joint));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=3;
    PCAtab{1,Col}='AmpGlobal_XK';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(AmpXM(a-1,Joint));
            PCAtab{a+Itermice,Col}=+(AmpXY(a-1,Joint));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=4;
    PCAtab{1,Col}='AmpGlobal_XH';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(AmpXM(a-1,Joint));
            PCAtab{a+Itermice,Col}=+(AmpXY(a-1,Joint));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=5;
    PCAtab{1,Col}='AmpGlobal_XI';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(AmpXM(a-1,Joint));
            PCAtab{a+Itermice,Col}=+(AmpXY(a-1,Joint));
            a=a+1;
            
        end
    end
    
    % AMP Y
    Col=Col+1;
    Joint=1;
    PCAtab{1,Col}='AmpGlobal_YF';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(AmpYM(a-1,Joint));
            PCAtab{a+Itermice,Col}=+(AmpYY(a-1,Joint));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=2;
    PCAtab{1,Col}='AmpGlobal_YA';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(AmpYM(a-1,Joint));
            PCAtab{a+Itermice,Col}=+(AmpYY(a-1,Joint));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=3;
    PCAtab{1,Col}='AmpGlobal_YK';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(AmpYM(a-1,Joint));
            PCAtab{a+Itermice,Col}=+(AmpYY(a-1,Joint));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=4;
    PCAtab{1,Col}='AmpGlobal_YH';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(AmpYM(a-1,Joint));
            PCAtab{a+Itermice,Col}=+(AmpYY(a-1,Joint));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=5;
    PCAtab{1,Col}='AmpGlobal_YI';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(AmpYM(a-1,Joint));
            PCAtab{a+Itermice,Col}=+(AmpYY(a-1,Joint));
            a=a+1;
            
        end
    end
    
    % AMP Z
    Col=Col+1;
    Joint=1;
    PCAtab{1,Col}='AmpGlobal_ZF';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(AmpZM(a-1,Joint));
            PCAtab{a+Itermice,Col}=+(AmpZY(a-1,Joint));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=2;
    PCAtab{1,Col}='AmpGlobal_ZA';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(AmpZM(a-1,Joint));
            PCAtab{a+Itermice,Col}=+(AmpZY(a-1,Joint));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=3;
    PCAtab{1,Col}='AmpGlobal_ZK';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(AmpZM(a-1,Joint));
            PCAtab{a+Itermice,Col}=+(AmpZY(a-1,Joint));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=4;
    PCAtab{1,Col}='AmpGlobal_ZH';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(AmpZM(a-1,Joint));
            PCAtab{a+Itermice,Col}=+(AmpZY(a-1,Joint));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=5;
    PCAtab{1,Col}='AmpGlobal_ZI';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(AmpZM(a-1,Joint));
            PCAtab{a+Itermice,Col}=+(AmpZY(a-1,Joint));
            a=a+1;
            
        end
    end
end

%Speed
for fold = 1
    Col=Col+1;
    Joint=1;
    PCAtab{1,Col}='Speed_XF';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=speedXM(a-1,Joint);
            %(ParticipationLowJointXM_thigh(a-1));
            PCAtab{a+Itermice,Col}=speedXY(a-1,Joint);
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=2;
    PCAtab{1,Col}='Speed_XA';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=speedXM(a-1,Joint);
            %(ParticipationLowJointXM_thigh(a-1));
            PCAtab{a+Itermice,Col}=speedXY(a-1,Joint);
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=3;
    PCAtab{1,Col}='Speed_XK';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=speedXM(a-1,Joint);
            %(ParticipationLowJointXM_thigh(a-1));
            PCAtab{a+Itermice,Col}=speedXY(a-1,Joint);
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=4;
    PCAtab{1,Col}='Speed_XH';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=speedXM(a-1,Joint);
            %(ParticipationLowJointXM_thigh(a-1));
            PCAtab{a+Itermice,Col}=speedXY(a-1,Joint);
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=5;
    PCAtab{1,Col}='Speed_XC';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=speedXM(a-1,Joint);
            %(ParticipationLowJointXM_thigh(a-1));
            PCAtab{a+Itermice,Col}=speedXY(a-1,Joint);
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=1;
    PCAtab{1,Col}='Speed_YF';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=speedYM(a-1,Joint);
            %(ParticipationLowJointXM_thigh(a-1));
            PCAtab{a+Itermice,Col}=speedYY(a-1,Joint);
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=2;
    PCAtab{1,Col}='Speed_YA';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=speedYM(a-1,Joint);
            %(ParticipationLowJointXM_thigh(a-1));
            PCAtab{a+Itermice,Col}=speedYY(a-1,Joint);
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=3;
    PCAtab{1,Col}='Speed_YK';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=speedYM(a-1,Joint);
            %(ParticipationLowJointXM_thigh(a-1));
            PCAtab{a+Itermice,Col}=speedYY(a-1,Joint);
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=4;
    PCAtab{1,Col}='Speed_YH';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=speedYM(a-1,Joint);
            %(ParticipationLowJointXM_thigh(a-1));
            PCAtab{a+Itermice,Col}=speedYY(a-1,Joint);
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=5;
    PCAtab{1,Col}='Speed_YC';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=speedYM(a-1,Joint);
            %(ParticipationLowJointXM_thigh(a-1));
            PCAtab{a+Itermice,Col}=speedYY(a-1,Joint);
            a=a+1;
            
        end
    end
    
    
    
    Col=Col+1;
    Joint=1;
    PCAtab{1,Col}='Speed_ZF';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=speedZM(a-1,Joint);
            %(ParticipationLowJointXM_thigh(a-1));
            PCAtab{a+Itermice,Col}=speedZY(a-1,Joint);
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=2;
    PCAtab{1,Col}='Speed_ZA';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=speedZM(a-1,Joint);
            %(ParticipationLowJointXM_thigh(a-1));
            PCAtab{a+Itermice,Col}=speedZY(a-1,Joint);
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=3;
    PCAtab{1,Col}='Speed_ZK';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=speedZM(a-1,Joint);
            %(ParticipationLowJointXM_thigh(a-1));
            PCAtab{a+Itermice,Col}=speedZY(a-1,Joint);
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=4;
    PCAtab{1,Col}='Speed_ZH';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=speedZM(a-1,Joint);
            %(ParticipationLowJointXM_thigh(a-1));
            PCAtab{a+Itermice,Col}=speedZY(a-1,Joint);
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=5;
    PCAtab{1,Col}='Speed_ZC';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=speedZM(a-1,Joint);
            %(ParticipationLowJointXM_thigh(a-1));
            PCAtab{a+Itermice,Col}=speedZY(a-1,Joint);
            a=a+1;
            
        end
    end
    
    
end

%Area under the curve steady
for fold = 1
    Col=Col+1;
    Joint=1;
    PCAtab{1,Col}='AreaUnder_XF';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            %PCAtab{a,Col}=trapz(XMmean10(j,:,Joint)-XMmean10(j,1,Joint));
            PCAtab{a,Col}=trapz(XM(1:18000,j,Joint)-mean(XM(1:300,j,Joint)));
            PCAtab{a+Itermice,Col}=trapz(XY(1:18000,j,Joint)-mean(XY(1:300,j,Joint)));
            %(ParticipationLowJointXM_thigh(a-1));
            %PCAtab{a+Itermice,Col}=trapz(XYmean10(j,:,Joint)-XYmean10(j,1,Joint));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=2;
    PCAtab{1,Col}='AreaUnder_XA';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            %PCAtab{a,Col}=trapz(XMmean10(j,:,Joint)-XMmean10(j,1,Joint));
            %(ParticipationLowJointXM_thigh(a-1));
            %PCAtab{a+Itermice,Col}=trapz(XYmean10(j,:,Joint)-XYmean10(j,1,Joint));
            PCAtab{a,Col}=trapz(XM(1:18000,j,Joint)-mean(XM(1:300,j,Joint)));
            PCAtab{a+Itermice,Col}=trapz(XY(1:18000,j,Joint)-mean(XY(1:300,j,Joint)));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=3;
    PCAtab{1,Col}='AreaUnder_XK';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            %PCAtab{a,Col}=trapz(XMmean10(j,:,Joint)-XMmean10(j,1,Joint));
            %(ParticipationLowJointXM_thigh(a-1));
            %PCAtab{a+Itermice,Col}=trapz(XYmean10(j,:,Joint)-XYmean10(j,1,Joint));
            PCAtab{a,Col}=trapz(XM(1:18000,j,Joint)-mean(XM(1:300,j,Joint)));
            PCAtab{a+Itermice,Col}=trapz(XY(1:18000,j,Joint)-mean(XY(1:300,j,Joint)));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=4;
    PCAtab{1,Col}='AreaUnder_XH';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            %PCAtab{a,Col}=trapz(XMmean10(j,:,Joint)-XMmean10(j,1,Joint));
            %(ParticipationLowJointXM_thigh(a-1));
            %PCAtab{a+Itermice,Col}=trapz(XMmean10(j,:,Joint)-XYmean10(j,1,Joint));
            PCAtab{a,Col}=trapz(XM(1:18000,j,Joint)-mean(XM(1:300,j,Joint)));
            PCAtab{a+Itermice,Col}=trapz(XY(1:18000,j,Joint)-mean(XY(1:300,j,Joint)));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=5;
    PCAtab{1,Col}='AreaUnder_XC';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            %PCAtab{a,Col}=trapz(XMmean10(j,:,Joint)-XMmean10(j,1,Joint));
            %(ParticipationLowJointXM_thigh(a-1));
            %PCAtab{a+Itermice,Col}=trapz(XYmean10(j,:,Joint)-XYmean10(j,1,Joint));
            PCAtab{a,Col}=trapz(XM(1:18000,j,Joint)-mean(XM(1:300,j,Joint)));
            PCAtab{a+Itermice,Col}=trapz(XY(1:18000,j,Joint)-mean(XY(1:300,j,Joint)));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=1;
    PCAtab{1,Col}='AreaUnder_YF';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            %PCAtab{a,Col}=trapz(YMmean10(j,:,Joint)-YMmean10(j,1,Joint));
            %(ParticipationLowJointXM_thigh(a-1));
            %PCAtab{a+Itermice,Col}=trapz(YYmean10(j,:,Joint)-YYmean10(j,1,Joint));
            PCAtab{a,Col}=trapz(YM(1:18000,j,Joint)-mean(YM(1:300,j,Joint)));
            PCAtab{a+Itermice,Col}=trapz(YY(1:18000,j,Joint)-mean(YY(1:300,j,Joint)));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=2;
    PCAtab{1,Col}='AreaUnder_YA';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            %PCAtab{a,Col}=trapz(YMmean10(j,:,Joint)-YMmean10(j,1,Joint));
            %(ParticipationLowJointXM_thigh(a-1));
            %PCAtab{a+Itermice,Col}=trapz(YYmean10(j,:,Joint)-YYmean10(j,1,Joint));
            PCAtab{a,Col}=trapz(YM(1:18000,j,Joint)-mean(YM(1:300,j,Joint)));
            PCAtab{a+Itermice,Col}=trapz(YY(1:18000,j,Joint)-mean(YY(1:300,j,Joint)));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=3;
    PCAtab{1,Col}='AreaUnder_YK';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            %PCAtab{a,Col}=trapz(YMmean10(j,:,Joint)-YMmean10(j,1,Joint));
            %(ParticipationLowJointXM_thigh(a-1));
            %PCAtab{a+Itermice,Col}=trapz(YYmean10(j,:,Joint)-YYmean10(j,1,Joint));
            PCAtab{a,Col}=trapz(YM(1:18000,j,Joint)-mean(YM(1:300,j,Joint)));
            PCAtab{a+Itermice,Col}=trapz(YY(1:18000,j,Joint)-mean(YY(1:300,j,Joint)));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=4;
    PCAtab{1,Col}='AreaUnder_YH';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            %PCAtab{a,Col}=trapz(YMmean10(j,:,Joint)-YMmean10(j,1,Joint));
            %(ParticipationLowJointXM_thigh(a-1));
            %PCAtab{a+Itermice,Col}=trapz(YYmean10(j,:,Joint)-YYmean10(j,1,Joint));
            PCAtab{a,Col}=trapz(YM(1:18000,j,Joint)-mean(YM(1:300,j,Joint)));
            PCAtab{a+Itermice,Col}=trapz(YY(1:18000,j,Joint)-mean(YY(1:300,j,Joint)));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=5;
    PCAtab{1,Col}='AreaUnder_YC';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            %PCAtab{a,Col}=trapz(YMmean10(j,:,Joint)-YMmean10(j,1,Joint));
            %(ParticipationLowJointXM_thigh(a-1));
            %PCAtab{a+Itermice,Col}=trapz(YYmean10(j,:,Joint)-YYmean10(j,1,Joint));
            PCAtab{a,Col}=trapz(YM(1:18000,j,Joint)-mean(YM(1:300,j,Joint)));
            PCAtab{a+Itermice,Col}=trapz(YY(1:18000,j,Joint)-mean(YY(1:300,j,Joint)));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=1;
    PCAtab{1,Col}='AreaUnder_ZF';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            %PCAtab{a,Col}=trapz(ZMmean10(j,:,Joint)-ZMmean10(j,1,Joint));
            %(ParticipationLowJointXM_thigh(a-1));
            %PCAtab{a+Itermice,Col}=trapz(ZYmean10(j,:,Joint)-ZYmean10(j,1,Joint));
            PCAtab{a,Col}=trapz(ZM(1:18000,j,Joint)-T(j));
            PCAtab{a+Itermice,Col}=trapz(ZY(1:18000,j,Joint)-TY(j));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=2;
    PCAtab{1,Col}='AreaUnder_ZA';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            %PCAtab{a,Col}=trapz(ZMmean10(j,:,Joint)-ZMmean10(j,1,Joint));
            %(ParticipationLowJointXM_thigh(a-1));
            %PCAtab{a+Itermice,Col}=trapz(ZYmean10(j,:,Joint)-ZYmean10(j,1,Joint));
            PCAtab{a,Col}=trapz(ZM(1:18000,j,Joint)-T(j));
            PCAtab{a+Itermice,Col}=trapz(ZY(1:18000,j,Joint)-TY(j));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=3;
    PCAtab{1,Col}='AreaUnder_ZK';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            %PCAtab{a,Col}=trapz(ZMmean10(j,:,Joint)-ZMmean10(j,1,Joint));
            %(ParticipationLowJointXM_thigh(a-1));
            %PCAtab{a+Itermice,Col}=trapz(ZYmean10(j,:,Joint)-ZYmean10(j,1,Joint));
            PCAtab{a,Col}=trapz(ZM(1:18000,j,Joint)-T(j));
            PCAtab{a+Itermice,Col}=trapz(ZY(1:18000,j,Joint)-TY(j));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=4;
    PCAtab{1,Col}='AreaUnder_ZH';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            %PCAtab{a,Col}=trapz(ZMmean10(j,:,Joint)-ZMmean10(j,1,Joint));
            %(ParticipationLowJointXM_thigh(a-1));
            %PCAtab{a+Itermice,Col}=trapz(ZYmean10(j,:,Joint)-ZYmean10(j,1,Joint));
            PCAtab{a,Col}=trapz(ZM(1:18000,j,Joint)-T(j));
            PCAtab{a+Itermice,Col}=trapz(ZY(1:18000,j,Joint)-TY(j));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=5;
    PCAtab{1,Col}='AreaUnder_ZC';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            %PCAtab{a,Col}=trapz(ZMmean10(j,:,Joint)-ZMmean10(j,1,Joint));
            %(ParticipationLowJointXM_thigh(a-1));
            %PCAtab{a+Itermice,Col}=trapz(ZYmean10(j,:,Joint)-ZYmean10(j,1,Joint));
            PCAtab{a,Col}=trapz(ZM(1:18000,j,Joint)-T(j));
            PCAtab{a+Itermice,Col}=trapz(ZY(1:18000,j,Joint)-TY(j));
            a=a+1;
            
        end
    end
    
end


%% DYNAMIC PARAMETERS


% AMP SPIKES MOVEMENT
for fold = 1
    Col=Col+1;
    Joint=1;
    PCAtab{1,Col}='AmpSpikes_XF';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(meanAmpPulseXMF(a-1));
            PCAtab{a+Itermice,Col}=+(meanAmpPulseXYF(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=2;
    PCAtab{1,Col}='AmpSpikes_XA';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(meanAmpPulseXMA(a-1));
            PCAtab{a+Itermice,Col}=+(meanAmpPulseXYA(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=3;
    PCAtab{1,Col}='AmpSpikes_XK';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(meanAmpPulseXMK(a-1));
            PCAtab{a+Itermice,Col}=+(meanAmpPulseXYK(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=4;
    PCAtab{1,Col}='AmpSpikes_XH';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(meanAmpPulseXMH(a-1));
            PCAtab{a+Itermice,Col}=+(meanAmpPulseXYH(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=5;
    PCAtab{1,Col}='AmpSpikes_XI';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(meanAmpPulseXMI(a-1));
            PCAtab{a+Itermice,Col}=+(meanAmpPulseXYI(a-1));
            a=a+1;
            
        end
    end
    
    % MEAN Y
    Col=Col+1;
    Joint=1;
    PCAtab{1,Col}='AmpSpikes_YF';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(meanAmpPulseYMF(a-1));
            PCAtab{a+Itermice,Col}=+(meanAmpPulseYYF(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=2;
    PCAtab{1,Col}='AmpSpikes_YA';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(meanAmpPulseYMA(a-1));
            PCAtab{a+Itermice,Col}=+(meanAmpPulseYYA(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=3;
    PCAtab{1,Col}='AmpSpikes_YK';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(meanAmpPulseYMK(a-1));
            PCAtab{a+Itermice,Col}=+(meanAmpPulseYYK(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=4;
    PCAtab{1,Col}='AmpSpikes_YH';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(meanAmpPulseYMH(a-1));
            PCAtab{a+Itermice,Col}=+(meanAmpPulseYYH(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=5;
    PCAtab{1,Col}='AmpSpikes_YI';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(meanAmpPulseYMI(a-1));
            PCAtab{a+Itermice,Col}=+(meanAmpPulseYYI(a-1));
            a=a+1;
            
        end
    end
    
    % MEAN Z
    Col=Col+1;
    Joint=1;
    PCAtab{1,Col}='AmpSpikes_ZF';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(meanAmpPulseZMF(a-1));
            PCAtab{a+Itermice,Col}=+(meanAmpPulseZYF(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=2;
    PCAtab{1,Col}='AmpSpikes_ZA';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(meanAmpPulseZMA(a-1));
            PCAtab{a+Itermice,Col}=+(meanAmpPulseZYA(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=3;
    PCAtab{1,Col}='AmpSpikes_ZK';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(meanAmpPulseZMK(a-1));
            PCAtab{a+Itermice,Col}=+(meanAmpPulseZYK(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=4;
    PCAtab{1,Col}='AmpSpikes_ZH';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(meanAmpPulseZMH(a-1));
            PCAtab{a+Itermice,Col}=+(meanAmpPulseZYH(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=5;
    PCAtab{1,Col}='AmpSpikes_ZI';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(meanAmpPulseZMI(a-1));
            PCAtab{a+Itermice,Col}=+(meanAmpPulseZYI(a-1));
            a=a+1;
            
        end
    end
end

% SPEED SPIKES MOVEMENT
for fold = 1
    Col=Col+1;
    Joint=1;
    PCAtab{1,Col}='SpeedSpikes_XF';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(meanDiffAmpPulseXMF(a-1));
            PCAtab{a+Itermice,Col}=+(meanDiffAmpPulseXYF(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=2;
    PCAtab{1,Col}='SpeedSpikes_XA';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(meanDiffAmpPulseXMA(a-1));
            PCAtab{a+Itermice,Col}=+(meanDiffAmpPulseXYA(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=3;
    PCAtab{1,Col}='SpeedSpikes_XK';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(meanDiffAmpPulseXMK(a-1));
            PCAtab{a+Itermice,Col}=+(meanDiffAmpPulseXYK(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=4;
    PCAtab{1,Col}='SpeedSpikes_XH';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(meanDiffAmpPulseXMH(a-1));
            PCAtab{a+Itermice,Col}=+(meanDiffAmpPulseXYH(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=5;
    PCAtab{1,Col}='SpeedSpikes_XI';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(meanDiffAmpPulseXMI(a-1));
            PCAtab{a+Itermice,Col}=+(meanDiffAmpPulseXYI(a-1));
            a=a+1;
            
        end
    end
    
    % MEAN Y
    Col=Col+1;
    Joint=1;
    PCAtab{1,Col}='SpeedSpikes_YF';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(meanDiffAmpPulseYMF(a-1));
            PCAtab{a+Itermice,Col}=+(meanDiffAmpPulseYYF(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=2;
    PCAtab{1,Col}='SpeedSpikes_YA';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(meanDiffAmpPulseYMA(a-1));
            PCAtab{a+Itermice,Col}=+(meanDiffAmpPulseYYA(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=3;
    PCAtab{1,Col}='SpeedSpikes_YK';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(meanDiffAmpPulseYMK(a-1));
            PCAtab{a+Itermice,Col}=+(meanDiffAmpPulseYYK(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=4;
    PCAtab{1,Col}='SpeedGlobal_YH';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(meanDiffAmpPulseYMH(a-1));
            PCAtab{a+Itermice,Col}=+(meanDiffAmpPulseYYH(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=5;
    PCAtab{1,Col}='SpeedSpikes_YI';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(meanDiffAmpPulseYMI(a-1));
            PCAtab{a+Itermice,Col}=+(meanDiffAmpPulseYYI(a-1));
            a=a+1;
            
        end
    end
    
    % MEAN Z
    Col=Col+1;
    Joint=1;
    PCAtab{1,Col}='SpeedSpikes_ZF';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(meanDiffAmpPulseZMF(a-1));
            PCAtab{a+Itermice,Col}=+(meanDiffAmpPulseZYF(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=2;
    PCAtab{1,Col}='SpeedSpikes_ZA';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(meanDiffAmpPulseZMA(a-1));
            PCAtab{a+Itermice,Col}=+(meanDiffAmpPulseZYA(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=3;
    PCAtab{1,Col}='SpeedSpikes_ZK';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(meanDiffAmpPulseZMK(a-1));
            PCAtab{a+Itermice,Col}=+(meanDiffAmpPulseZYK(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=4;
    PCAtab{1,Col}='SpeedSpikes_ZH';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(meanDiffAmpPulseZMH(a-1));
            PCAtab{a+Itermice,Col}=+(meanDiffAmpPulseZYH(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=5;
    PCAtab{1,Col}='SpeedSpikes_ZI';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(meanDiffAmpPulseZMI(a-1));
            PCAtab{a+Itermice,Col}=+(meanDiffAmpPulseZYI(a-1));
            a=a+1;
            
        end
    end
end

% AREA UNDER THE CURVE
for fold = 1
    Col=Col+1;
    Joint=1;
    PCAtab{1,Col}='A_PXMF';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(A_PXMF(a-1));
            PCAtab{a+Itermice,Col}=+(A_PXYF(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=2;
    PCAtab{1,Col}='A_PXMA';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(A_PXMA(a-1));
            PCAtab{a+Itermice,Col}=+(A_PXYA(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=3;
    PCAtab{1,Col}='A_PXMK';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(A_PXMK(a-1));
            PCAtab{a+Itermice,Col}=+(A_PXYK(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=4;
    PCAtab{1,Col}='A_PXMH';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(A_PXMH(a-1));
            PCAtab{a+Itermice,Col}=+(A_PXYH(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=5;
    PCAtab{1,Col}='A_PXMC';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(A_PXMC(a-1));
            PCAtab{a+Itermice,Col}=+(A_PXYC(a-1));
            a=a+1;
            
        end
    end
    
    % MEAN Y
    Col=Col+1;
    Joint=1;
    PCAtab{1,Col}='A_PYMF';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(A_PYMF(a-1));
            PCAtab{a+Itermice,Col}=+(A_PYYF(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=2;
    PCAtab{1,Col}='A_PYMA';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(A_PYMA(a-1));
            PCAtab{a+Itermice,Col}=+(A_PYYA(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=3;
    PCAtab{1,Col}='A_PYMK';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(A_PYMK(a-1));
            PCAtab{a+Itermice,Col}=+(A_PYYK(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=4;
    PCAtab{1,Col}='A_PYMH';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(A_PYMH(a-1));
            PCAtab{a+Itermice,Col}=+(A_PYYH(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=5;
    PCAtab{1,Col}='A_PYMC';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(A_PYMC(a-1));
            PCAtab{a+Itermice,Col}=+(A_PYYC(a-1));
            a=a+1;
            
        end
    end
    
    % MEAN Z
    Col=Col+1;
    Joint=1;
    PCAtab{1,Col}='A_PZMF';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(A_PZMF(a-1));
            PCAtab{a+Itermice,Col}=+(A_PZYF(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=2;
    PCAtab{1,Col}='A_PZMA';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(A_PZMA(a-1));
            PCAtab{a+Itermice,Col}=+(A_PZYA(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=3;
    PCAtab{1,Col}='A_PZMK';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(A_PZMK(a-1));
            PCAtab{a+Itermice,Col}=+(A_PZYK(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=4;
    PCAtab{1,Col}='A_PZMH';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(A_PZMH(a-1));
            PCAtab{a+Itermice,Col}=+(A_PZYH(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=5;
    PCAtab{1,Col}='A_PZMC';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(A_PZMC(a-1));
            PCAtab{a+Itermice,Col}=+(A_PZYC(a-1));
            a=a+1;
            
        end
    end
end

% ANGLE 3D DYN
for fold = 1
    Col=Col+1;
    Joint=1;
    PCAtab{1,Col}='Angle_3D_Ankle_Dyn';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(Angle3D_Ankle(a-1));
            PCAtab{a+Itermice,Col}=+(Angle3DY_Ankle(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=2;
    PCAtab{1,Col}='Angle_3D_Knee_Dyn';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(Angle3D_Knee(a-1));
            PCAtab{a+Itermice,Col}=+(Angle3DY_Knee(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=3;
    PCAtab{1,Col}='Angle3D_Hip_Dyn';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(Angle3D_Hip(a-1));
            PCAtab{a+Itermice,Col}=+(Angle3DY_Hip(a-1));
            a=a+1;
            
        end
    end
    
end

% MEAN 3D DYN
for fold = 1
    Col=Col+1;
    Joint=1;
    PCAtab{1,Col}='MeanDyn_XF';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(PXMF_M(a-1));
            PCAtab{a+Itermice,Col}=+(PXYF_M(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=1;
    PCAtab{1,Col}='MeanDyn_XA';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(PXMA_M(a-1));
            PCAtab{a+Itermice,Col}=+(PXYA_M(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=1;
    PCAtab{1,Col}='MeanDyn_XK';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(PXMK_M(a-1));
            PCAtab{a+Itermice,Col}=+(PXYK_M(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=1;
    PCAtab{1,Col}='MeanDyn_XH';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(PXMH_M(a-1));
            PCAtab{a+Itermice,Col}=+(PXYH_M(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=1;
    PCAtab{1,Col}='MeanDyn_XC';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(PXMC_M(a-1));
            PCAtab{a+Itermice,Col}=+(PXYC_M(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=1;
    PCAtab{1,Col}='MeanDyn_YF';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(PYMF_M(a-1));
            PCAtab{a+Itermice,Col}=+(PYYF_M(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=1;
    PCAtab{1,Col}='MeanDyn_YA';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(PYMA_M(a-1));
            PCAtab{a+Itermice,Col}=+(PYYA_M(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=1;
    PCAtab{1,Col}='MeanDyn_YK';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(PYMK_M(a-1));
            PCAtab{a+Itermice,Col}=+(PYYK_M(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=1;
    PCAtab{1,Col}='MeanDyn_YH';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(PYMH_M(a-1));
            PCAtab{a+Itermice,Col}=+(PYYH_M(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=1;
    PCAtab{1,Col}='MeanDyn_YC';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(PYMC_M(a-1));
            PCAtab{a+Itermice,Col}=+(PYYC_M(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=1;
    PCAtab{1,Col}='MeanDyn_ZF';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(PZMF_M(a-1));
            PCAtab{a+Itermice,Col}=+(PZYF_M(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=1;
    PCAtab{1,Col}='MeanDyn_ZA';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(PZMA_M(a-1));
            PCAtab{a+Itermice,Col}=+(PZYA_M(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=1;
    PCAtab{1,Col}='MeanDyn_ZK';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(PZMK_M(a-1));
            PCAtab{a+Itermice,Col}=+(PZYK_M(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=1;
    PCAtab{1,Col}='MeanDyn_ZH';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(PZMH_M(a-1));
            PCAtab{a+Itermice,Col}=+(PZYH_M(a-1));
            a=a+1;
            
        end
    end
    
    Col=Col+1;
    Joint=1;
    PCAtab{1,Col}='MeanDyn_ZC';
    a=2;
    for i=1:Ndays(nnz(Ndays))
        for j=1:nnz(Npairs(i,:))
            
            PCAtab{a,Col}=+(PZMC_M(a-1));
            PCAtab{a+Itermice,Col}=+(PZYC_M(a-1));
            a=a+1;
            
        end
    end
    
end


