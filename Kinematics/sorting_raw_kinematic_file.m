%% Load raw kimematic file
% You need the function importFileHorridge to make this section works
clc
close all
clear all

%Choose the file you want to import
ExpNb=1;
Day=1;
Npairs=1;

for Pair=Npairs

nameX=strcat(['ConditioningExample/DataExp' num2str(ExpNb) '/Day' num2str(Day) '/pair' num2str(Pair) '/Xpair' num2str(Pair) 'day' num2str(Day) '.txt']);
nameY=strcat(['ConditioningExample/DataExp',num2str(ExpNb),'/Day',num2str(Day),'/pair',num2str(Pair),'/Ypair',num2str(Pair),'day',num2str(Day),'.txt']);
nameZ=strcat(['ConditioningExample/DataExp' num2str(ExpNb) '/Day' num2str(Day) '/pair' num2str(Pair) '/Zpair' num2str(Pair) 'day' num2str(Day) '.txt']);
        
XRaw = importFileHorridge(nameX,1);
YRaw = importFileHorridge(nameY,1);
ZRaw = importFileHorridge(nameZ,1);

end

% Separating master and yoked
matrixParadigm = []; % if = 1: 5 first are master, if 0: 5 first are yoked
nonZeros = find(YRaw(:,1)); 
firstNonZero = nonZeros(1);
for i =  firstNonZero:length(YRaw) 
    if mean(YRaw(i,1:5)) < 200
        matrixParadigm = [matrixParadigm; 1];
    else
        matrixParadigm = [matrixParadigm; 0];
    end
end

for i = 1: length(matrixParadigm)
    if matrixParadigm(i) == 1
        XMaster_uncorrected(i,:) = XRaw(i+(firstNonZero - 1),1:5);
        XYoked_uncorrected(i,:) = XRaw(i+(firstNonZero - 1),6:10);
        YMaster_uncorrected(i,:) = YRaw(i+(firstNonZero - 1),1:5);
        YYoked_uncorrected(i,:) = YRaw(i+(firstNonZero - 1),6:10);
        ZMaster_uncorrected(i,:) = ZRaw(i+(firstNonZero - 1),1:5);
        ZYoked_uncorrected(i,:) = ZRaw(i+(firstNonZero - 1),6:10);
    else
        XMaster_uncorrected(i,:) = XRaw(i+(firstNonZero - 1),6:10);
        XYoked_uncorrected(i,:) = XRaw(i+(firstNonZero - 1),1:5);
        YMaster_uncorrected(i,:) = YRaw(i+(firstNonZero - 1),6:10);
        YYoked_uncorrected(i,:) = YRaw(i+(firstNonZero - 1),1:5);
        ZMaster_uncorrected(i,:) = ZRaw(i+(firstNonZero - 1),6:10);
        ZYoked_uncorrected(i,:) = ZRaw(i+(firstNonZero - 1),1:5);
    end
end



%% SEPARATE MASTER YOKED

threshold=450;

for i=1:18000
     
    %Master
    if YMaster_uncorrected(i,1)<threshold 
        
        buff=YMaster_uncorrected(i,:);
        YMaster_uncorrected(i,:)=YYoked_uncorrected(i,:);
        YYoked_uncorrected(i,:)=buff;
    
        buff=ZMaster_uncorrected(i,:);
        ZMaster_uncorrected(i,:)=ZYoked_uncorrected(i,:);
        ZYoked_uncorrected(i,:)=buff;
        
        buff=XMaster_uncorrected(i,:);
        XMaster_uncorrected(i,:)=XYoked_uncorrected(i,:);
        XYoked_uncorrected(i,:)=buff;
        
    end
    
    
    
end


%% BASIC SORTING BASED ON THE Z-AXIS
% Lower Z value label as the foot marker, etc,...

tic
matrixIndexSortingMaster = [];
ZMaster = [];
matrixIndexSortingYoked = [];
ZYoked = [];

for i = 1:length(ZMaster_uncorrected)
    lignToSort = ZMaster_uncorrected(i,1:5);
    [sortedDist_master,idx_master] = sort(lignToSort,'ascend');
    matrixIndexSortingMaster(i,:) = idx_master;
    ZMaster(i,:) = sortedDist_master;
end

for i = 1:length(ZYoked_uncorrected)
    lignToSort = ZYoked_uncorrected(i,1:5);
    [sortedDist_yoked,idx_yoked] = sort(lignToSort,'ascend');
    matrixIndexSortingYoked(i,:) = idx_yoked;
    ZYoked(i,:) = sortedDist_yoked;
end
toc

%% REMOVE 0 OUTLIERS
% 0 oultiers can happen if a marker was missed by the camera detection

matrixResetMaster = [];
matrixResetYoked = [];
for i = 1:length(ZMaster)
    for j = 1:5
        if ZMaster(i,j) < 5
            if i == 1
                idx = find(ZMaster(:,j)~=0, 1, 'first'); 
                ZMaster(i,j) = ZMaster(idx,j); 
            else
                ZMaster(i,:) = ZMaster(i-1,:);
                matrixResetMaster = [matrixResetMaster; i];
            end
        end
    end
end

for i = 1:length(ZYoked)
    for j = 1:5
        if ZYoked(i,j) < 5
            if i == 1
                idx = find(ZYoked(:,j)~=0, 1, 'first');
                if isempty(idx)
                    ZYoked(i,j)=0;
                else
                ZYoked(i,j) = ZYoked(idx,j);
                end
            else
                ZYoked(i,:) = ZYoked(i-1,:);
                matrixResetYoked = [matrixResetYoked; i];
            end
        end
    end
end

%% RESTRUCTURE THE X AND Y MARKERS BASED ON THE Z LABELING

for i = 1:length(XMaster_uncorrected)
    XMaster(i,:) = XMaster_uncorrected(i,matrixIndexSortingMaster(i,:));
    XYoked(i,:) = XYoked_uncorrected(i,matrixIndexSortingYoked(i,:));
    YMaster(i,:) = YMaster_uncorrected(i,matrixIndexSortingMaster(i,:));
    YYoked(i,:) = YYoked_uncorrected(i,matrixIndexSortingYoked(i,:));
end

for i = matrixResetMaster
    XMaster(i,:) = XMaster(i-1,:);
    YMaster(i,:) = YMaster(i-1,:);
end

for i = matrixResetYoked
    XYoked(i,:) = XYoked(i-1,:);
    YYoked(i,:) = YYoked(i-1,:);
end

for i = 1:length(ZMaster)
    for j = 1:5
        if XMaster(i,j) < 5
            if i == 1
                idx = find(XMaster(:,j)~=0, 1, 'first');
                XMaster(i,j) = XMaster(idx,j);
            else
                XMaster(i,:) = XMaster(i-1,:);
            end
        end
        if YMaster(i,j) < 5
            if i == 1
                idx = find(YMaster(:,j)~=0, 1, 'first');
                YMaster(i,j) = YMaster(idx,j);
            else
                YMaster(i,:) = YMaster(i-1,:);
            end
        end
    end
end

for i = 1:length(ZMaster)
    for j = 1:5
        if XYoked(i,j) < 5
            if i == 1
                idx = find(XYoked(:,j)~=0, 1, 'first');
               
             if isempty(idx)
                    XYoked(i,j)=0;
                else
                XYoked(i,j) = ZYoked(idx,j);
                end   
                
            else
                XYoked(i,:) = XYoked(i-1,:);
            end
        end
        if YYoked(i,j) < 5
            if i == 1
                idx = find(YYoked(:,j)~=0, 1, 'first');
                if isempty(idx)
                    YYoked(i,j)=0;
                else
                YYoked(i,j) = ZYoked(idx,j);
                end 
                
            else
                YYoked(i,:) = YYoked(i-1,:);
            end
        end
    end
end


%% MAKE SURE THE X-Y-Z IS CONSERVED FOR MARKER

for joint = 1:5
    potential_joint = 1:5;
    potential_joint = potential_joint(find(potential_joint~=joint));
    for i = 2:length(XMaster)
        for j = potential_joint
            pointNow = [XMaster(i,joint) YMaster(i,joint) ZMaster(i,joint)];
            pointBefore = [XMaster(i-1,joint) YMaster(i-1,joint) ZMaster(i-1,joint)];
            pointPotential = [XMaster(i,j) YMaster(i,j) ZMaster(i,j)];
            jointPotentialPointBefore = [XMaster(i-1,j) YMaster(i-1,j) ZMaster(i-1,j)];
            if norm(pointPotential-pointBefore)<norm(pointNow-pointBefore) && norm(pointNow-jointPotentialPointBefore)<norm(pointPotential-jointPotentialPointBefore)
                tempX = XMaster(i,joint);
                tempY = YMaster(i,joint);
                tempZ = ZMaster(i,joint);
                
                XMaster(i,joint) = XMaster(i,j);
                XMaster(i,j) = tempX;
                YMaster(i,joint) = YMaster(i,j);
                YMaster(i,j) = tempY;
                ZMaster(i,joint) = ZMaster(i,j);
                ZMaster(i,j) = tempZ;
            end
        end
    end
end

for joint = 1:5
    potential_joint =1:5;
    potential_joint = potential_joint(find(potential_joint~=joint));
    for i = 2:length(XMaster)
        for j = potential_joint
            pointNow = [XYoked(i,joint) YYoked(i,joint) ZYoked(i,joint)];
            pointBefore = [XYoked(i-1,joint) YYoked(i-1,joint) ZYoked(i-1,joint)];
            pointPotential = [XYoked(i,j) YYoked(i,j) ZYoked(i,j)];
            jointPotentialPointBefore = [XYoked(i-1,j) YYoked(i-1,j) ZYoked(i-1,j)];
            if norm(pointPotential-pointBefore)<norm(pointNow-pointBefore) && norm(pointNow-jointPotentialPointBefore)<norm(pointPotential-jointPotentialPointBefore)
                tempX = XYoked(i,joint);
                tempY = YYoked(i,joint);
                tempZ = ZYoked(i,joint);
                
                XYoked(i,joint) = XYoked(i,j);
                XYoked(i,j) = tempX;
                YYoked(i,joint) = YYoked(i,j);
                YYoked(i,j) = tempY;
                ZYoked(i,joint) = ZYoked(i,j);
                ZYoked(i,j) = tempZ;
            end
        end
    end
end



%% REMOVE OUTLIER PEAKS IN CASE MARKER ARE MISSLABELED

ThresholdHigh=10;

XMaster2=XMaster;
ZMaster2=ZMaster;
YMaster2=YMaster;

XYoked2=XYoked;
ZYoked2=ZYoked;
YYoked2=YYoked;


for joint=1:5
    for i=100:length(XMaster)
    
        if XMaster(i,joint)>(ThresholdHigh+mean(XMaster(:,joint))) || XMaster(i,joint)<(mean(XMaster(:,joint)-ThresholdHigh))
            Buff=mean(XMaster(i-50:i-25,joint));
            XMaster(i,joint)=Buff;
        end
        
        if YMaster(i,joint)>(ThresholdHigh+mean(YMaster(:,joint))) || YMaster(i,joint)<(mean(YMaster(:,joint)-ThresholdHigh))
            Buff=mean(YMaster(i-50:i-25,joint));
            YMaster(i,joint)=Buff;
        end
        
        if ZMaster(i,joint)>(ThresholdHigh+mean(ZMaster(:,joint))) || ZMaster(i,joint)<(mean(ZMaster(:,joint)-ThresholdHigh))
            Buff=mean(ZMaster(i-50:i-25,joint));
            ZMaster(i,joint)=Buff;
        end
        
    end
end


ThresholdHigh=5;

for joint=1:5
    for i=100:length(YYoked)
    
        if XYoked(i,joint)>(ThresholdHigh+mean(XYoked(:,joint))) || XYoked(i,joint)<(mean(XYoked(:,joint)-ThresholdHigh))
            Buff=mean(XYoked(i-50:i-25,joint));
            XYoked(i,joint)=Buff;
        end
        
        if YYoked(i,joint)>(ThresholdHigh+mean(YYoked(:,joint))) || YYoked(i,joint)<(mean(YYoked(:,joint)-ThresholdHigh))
            Buff=mean(YYoked(i-50:i-25,joint));
            YYoked(i,joint)=Buff;
        end
        
        if ZYoked(i,joint)>(ThresholdHigh+mean(ZYoked(:,joint))) || ZYoked(i,joint)<(mean(ZYoked(:,joint)-ThresholdHigh))
            Buff=mean(ZYoked(i-50:i-25,joint));
            ZYoked(i,joint)=Buff;
        end
        
    end
end


%% Z plot
for i = 1:5
    figure(1)
    plot(ZYoked(:,i))
    hold on
    if i == 5
        legend('F','A','K','H','C')
        title('Z - Yoked');
        zoom on
    end
    figure(2)
    plot(ZMaster(:,i))
    hold on
    if i == 5
        legend('F','A','K','H','C')
        title('Z - Master');
        zoom on
    end
end



%% Y plot
for i = 1:5
    figure(3)
    plot(YYoked(:,i))
    hold on
    if i == 5
        legend('F','A','K','H','C')
        title('Y - Yoked');
        zoom on
    end
    figure(4)
    plot(YMaster(:,i))
    hold on
    if i == 5
        legend('F','A','K','H','C')
        title('Y - Master');
        zoom on
    end
end

%% X plot
for i = 1:5
    figure(5)
    plot(XYoked(:,i))
    hold on
    if i == 5
        legend('F','A','K','H','C')
        title('X - Yoked');
        zoom on
    end
    figure(6)
    plot(XMaster(:,i))
    hold on
    if i == 5
        legend('F','A','K','H','C')
        title('X - Master');
        zoom on
    end
end

%% saving data

pathname = fileparts(['/Users/simonlavaud/Desktop/Github_Kinematics/ConditioningExample/DataExp' num2str(ExpNb) '/AnalyseDay' num2str(Day) '/pair' num2str(Pair) '/Xpair' num2str(Pair) 'day' num2str(Day)])

mkdir(pathname)

save([pathname '/ZmasterD' num2str(Day) 'P' num2str(Pair)],'ZMaster');
save([pathname '/ZyokedD' num2str(Day) 'P' num2str(Pair)],'ZYoked');

save([pathname '/YmasterD' num2str(Day) 'P' num2str(Pair)],'YMaster');
save([pathname '/YyokedD' num2str(Day) 'P' num2str(Pair)],'YYoked');

save([pathname '/XmasterD' num2str(Day) 'P' num2str(Pair)],'XMaster');
save([pathname '/XyokedD' num2str(Day) 'P' num2str(Pair)],'XYoked');
