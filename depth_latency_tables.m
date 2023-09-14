% Takeoka Lab - NERF empowered by imec, KU Leuven and VIB
% Author: Charlotte Bichara 
% 2023

%% Tables required, in : GitHub\takeokalab\Updated fr table
% load('frM_20220802_11corrected.mat')
% load('frM_UOI_nnmf_Cha.mat')
% load('frY_20220802_orrected.mat')
% load('frY_UOI_nnmf_Cha.mat')


% Determine the velocity of Primary afferent fiber and size of leg. 
%-----------Has to be fully run for any analysis of second order neurons


% Parameters to set :

% Choose to separate Ad and C (1) --Like figS7C-- or fuse them (0)
C_Ad_separated = 0;

% Velocities and distance of conduction parameters: -----------------------
% -------- > Used for paper
sABmin = 12; sABmax = 20;
sAdmin = 0.1; sAdmax = 12; 
sCmin = 0.01; sCmax = 0.09;
sPromin = 20; sPromax = 120; 
if C_Ad_separated == 1
    sAdmin = 5; sAdmax = 12; 
    sCmin = 0.1; sCmax = 2.7;
end

Lenght_leg=5; % average length from recording site to the furthest stimulation electrode (data not shown)
% -------------------------------------------------------------------------

% Compute latencies : Includes 1ms of synaptic delay
ABstart = (Lenght_leg/(sABmax*100))+0.001; ABstop = (Lenght_leg/(sABmin*100))+0.001;
Adstart = (Lenght_leg/(sAdmax*100))+0.001; Adstop = (Lenght_leg/(sAdmin*100))+0.001;
Cstart = (Lenght_leg/(sCmax*100))+0.001; Cstop = (Lenght_leg/(sCmin*100))+0.001;
Prostart = (Lenght_leg/(sPromax*100))+0.001; Prostop = (Lenght_leg/(sPromin*100))+0.001;

%% Create tables with all responses (not grouped by unit nb., info in the table)
% Learner

Resp_M =find(frM_UOI.Recording); % UOI on frM_UOI (tables of second order units)
[RespPtf1aON,pos]=intersect(find(frM.optotagged == 1),frM_UOI.Line); % find Ptf1aON neurons in second order
All_latencies = cat(2, frM_UOI.latency_all{Resp_M,:});  
tableau=table; % Create table for learners

Depths=[];
Units=[];
Rec=[];
Ptf1aON=[];
Category=[];
d=0;
WF = [];

for i=1:length(Resp_M)
    j=Resp_M(i);
    k=frM_UOI.Line(j);
    
    
    depths=frM_UOI.depth(j)*ones(1,length(frM_UOI.latency_all{j,:}));  
    Depths(d+1:d+length(depths))=depths;
     
    category=frM_UOI.Category(j)*ones(1,length(frM_UOI.latency_all{j,:}));  
    Category(d+1:d+length(category))=category;

    units=k*ones(1,length(frM_UOI.latency_all{j,:})); 
    Units(d+1:d+length(units))=units;
    
    rec=frM_UOI.Recording(j)*ones(1,length(frM_UOI.latency_all{j,:})); 
    Rec(d+1:d+length(rec))=rec;
    
    if frM_UOI.Recording(j) > 10 && frM_UOI.Recording(j) <100
        ptf1aON= frM.optotagged(frM_UOI.Line(j))*ones(1,length(frM_UOI.latency_all{j,:})); 
        Ptf1aON(d+1:d+length(ptf1aON))=ptf1aON;

        ptf1aInh= frM.inhibition(frM_UOI.Line(j))*ones(1,length(frM_UOI.latency_all{j,:})); 
        Ptf1aInh(d+1:d+length(ptf1aInh))=ptf1aInh;
    else
        ptf1aON= 0*ones(1,length(frM_UOI.latency_all{j,:})); 
        Ptf1aON(d+1:d+length(ptf1aON))=ptf1aON;
        ptf1aInh= 0*ones(1,length(frM_UOI.latency_all{j,:})); 
        Ptf1aInh(d+1:d+length(ptf1aInh))=ptf1aInh;
    end
    
    d=d+length(depths);
end

tableau.Latence=All_latencies';
tableau.Unit=Units';
tableau.Depth=-Depths';
tableau.Recording=Rec';
tableau.Category=Category';
tableau.Ptf1aON=Ptf1aON';
tableau.Ptf1aInh=Ptf1aInh';



%% Create tables with all responses (not grouped by unit nb., info in the table)
% Control

Resp_Y =find(frY_UOI.Recording);
All_latenciesY = cat(2, frY_UOI.latency_all{Resp_Y,:});
[RespPtf1aON_Y,pos]=intersect(find(frY.optotagged == 1),frY_UOI.Line); 

tableauY=table;

DepthsY=[];
UnitsY=[];
RecY=[];
Ptf1aON_Y=[];
WF_Y=[];
CategoryY=[];
d=0;
for i=1:length(Resp_Y)
    j=Resp_Y(i);
    k=frY_UOI.Line(j);
    depthsY=frY_UOI.depth(j)*ones(1,length(frY_UOI.latency_all{j,:}));  
    categoryY=frY_UOI.Category(j)*ones(1,length(frY_UOI.latency_all{j,:}));  
    unitsY=k*ones(1,length(frY_UOI.latency_all{j,:})); 
    DepthsY(d+1:d+length(depthsY))=depthsY;
    CategoryY(d+1:d+length(categoryY))=categoryY;
    UnitsY(d+1:d+length(unitsY))=unitsY;
    recY=frY_UOI.Recording(j)*ones(1,length(frY_UOI.latency_all{j,:})); 
    RecY(d+1:d+length(recY))=recY;

    if frY_UOI.Recording(j) > 9 && frY_UOI.Recording(j) <100
        ptf1aONY= frY.optotagged(frY_UOI.Line(j))*ones(1,length(frY_UOI.latency_all{j,:})); 
        Ptf1aONY(d+1:d+length(ptf1aONY))=ptf1aONY;
        ptf1aInhY= frY.inhibition(frY_UOI.Line(j))*ones(1,length(frY_UOI.latency_all{j,:})); 
        Ptf1aInhY(d+1:d+length(ptf1aInhY))=ptf1aInhY;
    else
        ptf1aONY= 0*ones(1,length(frY_UOI.latency_all{j,:})); 
        Ptf1aONY(d+1:d+length(ptf1aONY))=ptf1aONY;
        ptf1aInhY= 0*ones(1,length(frY_UOI.latency_all{j,:})); 
        Ptf1aInhY(d+1:d+length(ptf1aInhY))=ptf1aInhY;
    end

    d=d+length(depthsY);
end

tableauY.Unit=UnitsY';
tableauY.Latence=All_latenciesY';
tableauY.Depth=-DepthsY';
tableauY.Recording=RecY';
tableauY.Category=CategoryY';
tableauY.Ptf1aON=Ptf1aONY';
tableauY.Ptf1aInh=Ptf1aInhY';


%% Pie tuning fiber type
% Identify units in each second order groups as lines of the created
% "tableau(Y)"
% WARNING : In the current version of the paper, Ad and C are fused
% together by shifting the windows of velocities, So C groups are empty.

C = find(tableau.Latence >= Cstart);% & tableau.Category ~= 0);% & tableau.Depth>-700);
Cy = find(tableauY.Latence>= Cstart);% & tableauY.Category ~= 0);% & tableauY.Depth>-700);
Pro = find(tableau.Latence < Prostop);% & tableau.Category  ~= 0);% & tableau.Depth>-700);
Proy = find(tableauY.Latence < Prostop);% & tableauY.Category ~= 0);% & tableauY.Depth>-700);
Ab = find(tableau.Latence >= ABstart & tableau.Latence < ABstop);% & tableau.Category ~= 0); %& tableau.Depth>-700);
Aby = find(tableauY.Latence >= ABstart & tableauY.Latence < ABstop);% & tableauY.Category ~= 0);%& tableauY.Depth>-700);
Ad = find(tableau.Latence >= Adstart & tableau.Latence < Adstop);%  & tableau.Category ~= 0);%& tableau.Depth>-700);
Ady = find(tableauY.Latence >= Adstart & tableauY.Latence < Adstop);% & tableauY.Category ~= 0); %& tableauY.Depth>-700);


%% Find phase tuning category (Up/down... compared to the baseline)

FibersM = {Ab, Ad, C, Pro};
FibersY = {Aby, Ady, Cy, Proy};
Fibres = {'AB', 'Ad', 'C', 'P'};

DepthsNP_M = [];
DepthsSI_M = [];
DepthsEI_M = [];
DepthsLI_M = [];

DepthsNP_Y = [];
DepthsSI_Y = [];
DepthsEI_Y = [];
DepthsLI_Y = [];

% Run with i = 1 to 4 to plot your desired Second order type
for i = 1:length(FibersM)
    R=[];Ry=[];
    R = FibersM{i}; 
    Ry=FibersY{i};

    
    if i == 1
        ABunits = tableau.Unit(R);
        ABunitsY = tableauY.Unit(Ry);
    elseif i == 2
        Adunits = tableau.Unit(R);
        AdunitsY = tableauY.Unit(Ry);
        
    elseif i == 3
        Cunits = tableau.Unit(R);
        CunitsY = tableauY.Unit(Ry);
        
    elseif i == 4
        Punits = tableau.Unit(R);
        PunitsY = tableauY.Unit(Ry);
    end

    DepthsNP_M = tableau.Depth(find(tableau.Category(R) == 10 | isnan(tableau.Category(R))));
    DepthsSI_M = tableau.Depth(find(tableau.Category(R) == 1 | tableau.Category(R) == 2));
    DepthsEI_M = tableau.Depth(find(tableau.Category(R) == 3| tableau.Category(R) == 4));
    DepthsLI_M = tableau.Depth(find(tableau.Category(R) == 6 | tableau.Category(R) == 5));

    DepthsNP_Y = tableau.Depth(find(tableauY.Category(Ry) == 10 | isnan(tableauY.Category(Ry))));
    DepthsSI_Y = tableau.Depth(find(tableauY.Category(Ry) == 1 | tableauY.Category(Ry) == 2));
    DepthsEI_Y = tableau.Depth(find(tableauY.Category(Ry) == 3| tableauY.Category(Ry) == 4));
    DepthsLI_Y = tableau.Depth(find(tableauY.Category(Ry) == 6 | tableauY.Category(Ry) == 5));

end

%% Histogram latencies

Lat_hist=[];
figure;
for i=1:numel(tableau.Recording)
   Lat_iterM=tableau.Latence(i);
    Lat_hist=[Lat_hist Lat_iterM];
end
hist(Lat_hist,60, 'm')
title('Learner')
xlabel('latency [s]')

rectangle('Position',[Prostart,40,Prostop-Prostart,2],'FaceColor','b')
rectangle('Position',[ABstart,40,ABstop-ABstart,2],'FaceColor','r')
rectangle('Position',[Adstart,40,Adstop-Adstart,2],'FaceColor','g')
if C_Ad_separated == 1
    rectangle('Position',[Cstart,40,Cstop-Cstart,2],'FaceColor','cyan')
end
xlim([0 0.04])

figure;
Lat_histY=[];
for i=1:numel(tableauY.Recording)
    Lat_iter=tableauY.Latence(i);
    Lat_histY=[Lat_histY Lat_iter];
end
hist(Lat_histY,60, 'b')
title('Control')
xlabel('latency [s]')

rectangle('Position',[Prostart,40,Prostop-Prostart,2],'FaceColor','b')
rectangle('Position',[ABstart,40,ABstop-ABstart,2],'FaceColor','r')
rectangle('Position',[Adstart,40,Adstop-Adstart,2],'FaceColor','g')
if C_Ad_separated == 1
    rectangle('Position',[Cstart,40,Cstop-Cstart,2],'FaceColor','cyan')
end
xlim([0 0.04])

%%  Histogram Depths

dep_hist=[];
figure;
for i=1:numel(tableau.Recording)
   dep_iterM=-tableau.Depth(i);
    dep_hist=[dep_hist dep_iterM];
end

h1 = histogram(dep_hist);
h1.Normalization = 'count';
h1.BinWidth = 100;
h1.FaceColor ='m';
title('Learner')
xlabel('deths [μm]')
xlim([0 2000])

figure;
dep_histY=[];
for i=1:numel(tableauY.Recording)
    dep_iter=-tableauY.Depth(i);
    dep_histY=[dep_histY dep_iter];
end
h2 = histogram(dep_histY);
h2.Normalization = 'count';
h2.BinWidth = 100;
h2.FaceColor ='b';
title('Control')
xlabel('deths [μm]')
xlim([0 2000])

% find % of each Ftype in each depth bin
a=0;b=-1;
for i = 1:20
    a=b; b=-(i*100);
    %P
    binM(1,i) = (numel(find(tableau.Depth(Pro)<=a & tableau.Depth(Pro)>=b))*100)/numel(find(tableau.Depth<=a & tableau.Depth>=b));
    %Ab
    binM(2,i) = (numel(find(tableau.Depth(Ab)<=a & tableau.Depth(Ab)>=b))*100)/numel(find(tableau.Depth<=a & tableau.Depth>=b));
    %Ad
    binM(3,i) =(numel(find(tableau.Depth(Ad)<=a & tableau.Depth(Ad)>=b))*100)/numel(find(tableau.Depth<=a & tableau.Depth>=b));
    binY(1,i) = (numel(find(tableauY.Depth(Proy)<=a & tableauY.Depth(Proy)>=b))*100)/numel(find(tableauY.Depth<=a & tableauY.Depth>=b));
    binY(2,i) = (numel(find(tableauY.Depth(Aby)<=a & tableauY.Depth(Aby)>=b))*100)/numel(find(tableauY.Depth<=a & tableauY.Depth>=b));
    binY(3,i) = (numel(find(tableauY.Depth(Ady)<=a & tableauY.Depth(Ady)>=b))*100)/numel(find(tableauY.Depth<=a & tableauY.Depth>=b));
end
