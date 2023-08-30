%% Run PCA 
% - - - - - - - - 
% - - - - - - - - 

%%
close all
clc
w = warning ('off','all');

% % % % % % % % % % % %  
% % - - Load ALL the mice
% % % % % % % % % % % %  

clearvars -except PCAKinematicsWTandAblatedCellType 
ImportTable= PCAKinematicsWTandAblatedCellType;
Data=table2array(ImportTable);
[n p]=size(Data);

% % % % % % % % % % % %  
% % - - Selection of the mice to display. Day 1 M/Y no modification
% % % % % % % % % % % %  

iter=1;
iterY=1;

for i=1:n
    
    if Data(i,1)==1 || Data(i,1)==2 || Data(i,1)==5 || Data(i,1)==12 || Data(i,1)==13 || Data(i,1)==14 || Data(i,1)==15 || Data(i,1)==19 || Data(i,1)==27 %Experience 1 5 12 13 14 15 19 27
    	% Extract Day 1
        if (Data(i,4)==1) 
            %Sort Learner and Control
            if Data(i,3)==1
            
            Data_Day1_Learner(iter,:)=Data(i,:);
            iter=iter+1;
            else
                
           	Data_Day1_Control(iterY,:)=Data(i,:);
            iterY=iterY+1;

            end
            
    end
    
    end

end


% % % % % % % % % % % %       
% % - - Sort By Experiment type
% % % % % % % % % % % %  

iterMSwitch=1;
iterYSwitch=1;
iterMSwitch2=1;
iterYSwitch2=1;

iterMnaive=1;
iterYnaive=1;

for i=1:n
    
    %Cond 1 - Mice with switch on day 2
    if Data(i,2)==1 && Data(i,4)==1
        if Data(i,3)==1
        Data_Day1_Switch_Learner(iterMSwitch,:)=Data(i,:);
        iterMSwitch=iterMSwitch+1;
        else 
      	Data_Day1_Switch_Control(iterYSwitch,:)=Data(i,:);
        iterYSwitch=iterYSwitch+1;
        end  
    end
    
    if Data(i,2)==1 && Data(i,4)==2
        if Data(i,3)==1
        Data_Day2_Switch_Learner(iterMSwitch2,:)=Data(i,:);
        iterMSwitch2=iterMSwitch2+1;
        else 
      	Data_Day2_Switch_Control(iterYSwitch2,:)=Data(i,:);
        iterYSwitch2=iterYSwitch2+1;
        end  
    end
    
    
    %Cond 10 - Mice just with day 1 training
    if Data(i,2)==10 && Data(i,4)==1
        if Data(i,3)==1
        Data_Day1_training_Learner(iterMnaive,:)=Data(i,:);
        iterMnaive=iterMnaive+1;
        else 
      	Data_Day1_training_Control(iterYnaive,:)=Data(i,:);
        iterYnaive=iterYnaive+1;
        end  
    end
    
end


%%
% % % % % % % % % % % %     
% % - - 1/ PCA TO USE WITH KINEMATIC STANDARDISATION 
% % - - - Based on full 10min movement
% % % % % % % % % % % % 

clearvars StandKine DataPCA SCORE COEFF

StandKine=[];
DataPCA_origin=[Data_Day1_Learner(:,:);Data_Day1_Control(:,:)];

%Export size
[n_norm p_norm] = size(DataPCA_origin);
e = ones(n_norm,1);

% - Binary number related to the position of the foot at the end of the
%   experiment
% - Percentage of time spend below the treshold
inter=[6;7];
[a b]=size(inter);
StandKine=[];
TabOfIn=[];
StandKine_notnormalized=[];

for i=1:a
    TabOfIn=[];
    TabOfIn=DataPCA_origin(:,inter(i,:));
    
    [n_norm p_norm] = size(TabOfIn);
    e = ones(n_norm,1);

    StandKine=[StandKine (TabOfIn-min(min(TabOfIn)))./(max(max(TabOfIn))-min(min(TabOfIn)))];
    StandKine_notnormalized=[StandKine_notnormalized TabOfIn];
end

% Kinematic 3D Angle
inter=[23:25];
[a b]=size(inter);
TabOfIn=[];
for i=1:a
    TabOfIn=[];
    TabOfIn=DataPCA_origin(:,inter(i,:));
    
    [n_norm p_norm] = size(TabOfIn);
    e = ones(n_norm,1);

    StandKine=[StandKine (TabOfIn-min(min(TabOfIn)))./(max(max(TabOfIn))-min(min(TabOfIn)))];
    StandKine_notnormalized=[StandKine_notnormalized TabOfIn];
end


% Mean position of the joint kinematic 
% Mean amplitude of the joint kinematic
% Mean speed of the joint kinematic
inter=[8:12;13:17;18:22; ...
    26:30;31:35;36:40; ...
    41:45;46:50;51:55];

[a b]=size(inter);
TabOfIn=[];

for i=1:a
    TabOfIn=[];
    TabOfIn=DataPCA_origin(:,inter(i,:));
    
    [n_norm p_norm] = size(TabOfIn);
    e = ones(n_norm,1);

    StandKine=[StandKine (TabOfIn-min(min(TabOfIn)))./(max(max(TabOfIn))-min(min(TabOfIn)))];
    StandKine_notnormalized=[StandKine_notnormalized TabOfIn];
end


% Mean area under the curve (absement) of the joint kinematic
inter=[56:60;61:65;66:70];

[a b]=size(inter);
TabOfIn=[];

for i=1:a
    TabOfIn=[];
    TabOfIn=DataPCA_origin(:,inter(i,:));
    
    [n_norm p_norm] = size(TabOfIn);
    e = ones(n_norm,1);

    StandKine=[StandKine (TabOfIn-min(min(TabOfIn)))./(max(max(TabOfIn))-min(min(TabOfIn)))];
    StandKine_notnormalized=[StandKine_notnormalized TabOfIn];
end

% Run the PCA
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED]=pca(StandKine);
RecoDataCentered=SCORE*COEFF';


% Plot the results 
day1_noretrain=1:25;
numel(day1_noretrain);
scatter3(SCORE(1:numel(day1_noretrain),1),SCORE(1:numel(day1_noretrain),2),SCORE(1:numel(day1_noretrain),3),'r')
hold on
scatter3(SCORE((1:numel(day1_noretrain))+numel(day1_noretrain),1),SCORE((1:numel(day1_noretrain))+numel(day1_noretrain),2),SCORE((1:numel(day1_noretrain))+numel(day1_noretrain),3),'b')

Display_1C_Master=SCORE(1:25,1:3);
Display_1C_Yoked=SCORE(26:50,1:3);

% Input for the 95% ellipsioid
y1 = SCORE(1:25,1); 
y2 = SCORE(1:25,2); 
y1_prime = SCORE(26:50,1); 
y2_prime = SCORE(26:50,2); 
data = [y1 y2];
data_prime = [y1_prime y2_prime];

% Calculate the eigenvectors and eigenvalues
covariance = cov(data);
[eigenvec, eigenval ] = eig(covariance);
covariance_prime = cov(data_prime);
[eigenvec_p, eigenval_p ] = eig(covariance_prime);

% Get the index of the largest eigenvector
[largest_eigenvec_ind_c, r] = find(eigenval == max(max(eigenval)));
largest_eigenvec = eigenvec(:, largest_eigenvec_ind_c);
[largest_eigenvec_ind_c_p, r_p] = find(eigenval_p == max(max(eigenval_p)));
largest_eigenvec_p = eigenvec_p(:, largest_eigenvec_ind_c_p);

% Get the largest eigenvalue
largest_eigenval = max(max(eigenval));
largest_eigenval_p = max(max(eigenval_p));

% Get the smallest eigenvector and eigenvalue
if(largest_eigenvec_ind_c == 1)
    smallest_eigenval = max(eigenval(:,2))
    smallest_eigenvec = eigenvec(:,2);
else
    smallest_eigenval = max(eigenval(:,1))
    smallest_eigenvec = eigenvec(1,:);
end


if(largest_eigenvec_ind_c_p == 1)
    smallest_eigenval_p = max(eigenval_p(:,2))
    smallest_eigenvec_p = eigenvec_p(:,2);
else
    smallest_eigenval_p = max(eigenval_p(:,1))
    smallest_eigenvec_p = eigenvec_p(1,:);
end

% Calculate the angle between the x-axis and the largest eigenvector
angle = atan2(largest_eigenvec(2), largest_eigenvec(1));
angle_p = atan2(largest_eigenvec_p(2), largest_eigenvec_p(1));
if(angle_p < 0)
    angle_p = angle_p + 2*pi;
end

if(angle_p < 0)
    angle_p = angle_p + 2*pi;
end

% Get the coordinates of the data mean
avg = mean(data);
avg_p = mean(data_prime);

% Get the 95% confidence interval error ellipse
chisquare_val = 2.4477;
theta_grid = linspace(0,2*pi);

phi = angle;
phi_p = angle_p;
X0=avg(1);
Y0=avg(2);
X0_p=avg_p(1);
Y0_p=avg_p(2);

a=chisquare_val*sqrt(largest_eigenval);
b=chisquare_val*sqrt(smallest_eigenval);
a_p=chisquare_val*sqrt(largest_eigenval_p);
b_p=chisquare_val*sqrt(smallest_eigenval_p);

% the ellipse in x and y coordinates 
ellipse_x_r  = a*cos( theta_grid );
ellipse_y_r  = b*sin( theta_grid );

ellipse_x_r_p  = a_p*cos( theta_grid );
ellipse_y_r_p  = b_p*sin( theta_grid );

%Define a rotation matrix
R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];
R_p = [ cos(phi_p) sin(phi_p); -sin(phi_p) cos(phi_p) ];

%let's rotate the ellipse to some angle phi
r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;
r_ellipse_p = [ellipse_x_r_p;ellipse_y_r_p]' * R_p;

% Draw the error ellipse
plot(r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0,'-b')
hold on;
plot(r_ellipse_p(:,1) + X0_p,r_ellipse_p(:,2) + Y0_p,'-r')


hold on;
hXLabel = xlabel('x');
hYLabel = ylabel('y');


%%
% % % % % % % % % % % %     
% % - - 2/ PCA TO USE WITH KINEMATIC STANDARDISATION 
% % - - - Based on withdrawing movement - 500ms
% % % % % % % % % % % % 

clearvars StandKine DataPCA 

StandKine=[];
DataPCA_origin=[Data_Day1_Learner(:,:);Data_Day1_Control(:,:)];


% Mean position of the joint kinematic 
% Mean amplitude of the joint kinematic
% Mean speed of the joint kinematic
%
inter=[119:123;124:128;129:133; ...
    71:75;76:80;81:85; ...
    86:90;91:95;96:100;...
    101:105;106:110;111:115];

[a b]=size(inter);
TabOfIn=[];
StandKine_notnormalized=[];

for i=1:a
    TabOfIn=[];
    TabOfIn=DataPCA_origin(:,inter(i,:));
    
    [n_norm p_norm] = size(TabOfIn);
    e = ones(n_norm,1);

    StandKine=[StandKine (TabOfIn-min(min(TabOfIn)))./(max(max(TabOfIn))-min(min(TabOfIn)))];
    StandKine_notnormalized=[StandKine_notnormalized TabOfIn];

end


inter=[116:118];
[a , ~]=size(inter);
TabOfIn=[];
for i=1:a
    TabOfIn=[];
    TabOfIn=DataPCA_origin(:,inter(i,:));
    
    [n_norm p_norm] = size(TabOfIn);
    e = ones(n_norm,1);

    StandKine=[StandKine (TabOfIn-min(min(TabOfIn)))./(max(max(TabOfIn))-min(min(TabOfIn)))];
    StandKine_notnormalized=[StandKine_notnormalized TabOfIn];

end


clearvars COEFF SCORE LATENT TSQUARED EXPLAINED
[COEFF_Dyn, SCORE_Dyn, LATENT, TSQUARED, EXPLAINED_Dyn]=pca(StandKine);%,'Centered',true);%,'VariableWeights','variance');
RecoDataCentered=SCORE_Dyn*COEFF_Dyn';

%Normalize COEFF
COEFF_Dyn;
[n_coeff p_coeff] = size(COEFF_Dyn);
e = ones(n_coeff,1);
COEFF_Norm=(2.*COEFF_Dyn-e*(min(COEFF_Dyn)+max(COEFF_Dyn)))./(e*max(COEFF_Dyn)-e*min(COEFF_Dyn));
COEFF_Norm=(COEFF_Norm).*(COEFF_Norm);

day1_noretrain=1:25;
numel(day1_noretrain);
scatter3(SCORE_Dyn(1:numel(day1_noretrain),1),SCORE_Dyn(1:numel(day1_noretrain),2),SCORE_Dyn(1:numel(day1_noretrain),3),'r')
hold on
scatter3(SCORE_Dyn((1:numel(day1_noretrain))+numel(day1_noretrain),1),SCORE_Dyn((1:numel(day1_noretrain))+numel(day1_noretrain),2),SCORE_Dyn((1:numel(day1_noretrain))+numel(day1_noretrain),3),'b')

Display_1E_Master=SCORE_Dyn(1:25,1:3);
Display_1E_Yoked=SCORE_Dyn(26:50,1:3);


% Create some random data
y1 = SCORE_Dyn(1:numel(day1_noretrain),1); %normrnd(s(1).*x,1);
y2 = SCORE_Dyn(1:numel(day1_noretrain),2); %normrnd(s(2).*x,1);

[a b]=ttest2(y1,y2);

y1_prime = SCORE_Dyn((1:numel(day1_noretrain))+numel(day1_noretrain),1); %normrnd(s(1).*x,1);
y2_prime =SCORE_Dyn((1:numel(day1_noretrain))+numel(day1_noretrain),2); %normrnd(s(2).*x,1);
data = [y1 y2];
data_prime = [y1_prime y2_prime];
% Calculate the eigenvectors and eigenvalues
covariance = cov(data);
[eigenvec, eigenval ] = eig(covariance);

covariance_prime = cov(data_prime);
[eigenvec_p, eigenval_p ] = eig(covariance_prime);

% Get the index of the largest eigenvector
[largest_eigenvec_ind_c, r] = find(eigenval == max(max(eigenval)));
largest_eigenvec = eigenvec(:, largest_eigenvec_ind_c);

[largest_eigenvec_ind_c_p, r_p] = find(eigenval_p == max(max(eigenval_p)));
largest_eigenvec_p = eigenvec_p(:, largest_eigenvec_ind_c_p);


% Get the largest eigenvalue
largest_eigenval = max(max(eigenval));
largest_eigenval_p = max(max(eigenval_p));

% Get the smallest eigenvector and eigenvalue
if(largest_eigenvec_ind_c == 1)
    smallest_eigenval = max(eigenval(:,2))
    smallest_eigenvec = eigenvec(:,2);
else
    smallest_eigenval = max(eigenval(:,1))
    smallest_eigenvec = eigenvec(1,:);
end


if(largest_eigenvec_ind_c_p == 1)
    smallest_eigenval_p = max(eigenval_p(:,2))
    smallest_eigenvec_p = eigenvec_p(:,2);
else
    smallest_eigenval_p = max(eigenval_p(:,1))
    smallest_eigenvec_p = eigenvec_p(1,:);
end

% Calculate the angle between the x-axis and the largest eigenvector
angle = atan2(largest_eigenvec(2), largest_eigenvec(1));

angle_p = atan2(largest_eigenvec_p(2), largest_eigenvec_p(1));

% This angle is between -pi and pi.
% Let's shift it such that the angle is between 0 and 2pi
if(angle_p < 0)
    angle_p = angle_p + 2*pi;
end

if(angle_p < 0)
    angle_p = angle_p + 2*pi;
end

% Get the coordinates of the data mean
avg = mean(data);

avg_p = mean(data_prime);

% Get the 95% confidence interval error ellipse
chisquare_val = 2.4477;
theta_grid = linspace(0,2*pi);

phi = angle;
phi_p = angle_p;
X0=avg(1);
Y0=avg(2);
X0_p=avg_p(1);
Y0_p=avg_p(2);

a=chisquare_val*sqrt(largest_eigenval);
b=chisquare_val*sqrt(smallest_eigenval);

a_p=chisquare_val*sqrt(largest_eigenval_p);
b_p=chisquare_val*sqrt(smallest_eigenval_p);

% the ellipse in x and y coordinates 
ellipse_x_r  = a*cos( theta_grid );
ellipse_y_r  = b*sin( theta_grid );

ellipse_x_r_p  = a_p*cos( theta_grid );
ellipse_y_r_p  = b_p*sin( theta_grid );

%Define a rotation matrix
R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];
R_p = [ cos(phi_p) sin(phi_p); -sin(phi_p) cos(phi_p) ];

%let's rotate the ellipse to some angle phi
r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;
r_ellipse_p = [ellipse_x_r_p;ellipse_y_r_p]' * R_p;


% Draw the error ellipse
plot(r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0,'-b')
hold on;
plot(r_ellipse_p(:,1) + X0_p,r_ellipse_p(:,2) + Y0_p,'-r')

hold on;

% Set the axis labels
hXLabel = xlabel('x');
hYLabel = ylabel('y');



%%
% % % % % % % % % % % %     
% % - - 3/ PCA FOR THE SWITCH EXPERIMENT 
% % - - - Based on 10 minutes experiment
% % % % % % % % % % % % 


clearvars StandKine DataPCA SCORE COEFF DataPCA_origin 
StandKine=[];
StandKine_notnormalized=[];

DataPCA_origin=[...
Data_Day1_Switch_Learner; ...
Data_Day2_Switch_Control;...
Data_Day1_Switch_Control;...
Data_Day2_Switch_Learner;...
];

%Export size
[n_norm p_norm] = size(DataPCA_origin);
e = ones(n_norm,1);

% - Binary number related to the position of the foot at the end of the
%   experiment
% - Percentage of time spend below the treshold
inter=[6;7];
[a b]=size(inter);
StandKine=[];
TabOfIn=[];

for i=1:a
    TabOfIn=[];
    TabOfIn=DataPCA_origin(:,inter(i,:));
    
    [n_norm p_norm] = size(TabOfIn);
    e = ones(n_norm,1);

    StandKine=[StandKine (TabOfIn-min(min(TabOfIn)))./(max(max(TabOfIn))-min(min(TabOfIn)))];

end

% Angle 3D 
inter=[23:25];
[a b]=size(inter);
TabOfIn=[];
for i=1:a
    TabOfIn=[];
    TabOfIn=DataPCA_origin(:,inter(i,:));
    
    [n_norm p_norm] = size(TabOfIn);
    e = ones(n_norm,1);

    StandKine=[StandKine (TabOfIn-min(min(TabOfIn)))./(max(max(TabOfIn))-min(min(TabOfIn)))];

end

% Mean position of the joint kinematic 
% Mean amplitude of the joint kinematic
% Mean speed of the joint kinematic
inter=[8:12;13:17;18:22; ...
    26:30;31:35;36:40; ...
    41:45;46:50;51:55];

[a b]=size(inter);
TabOfIn=[];

for i=1:a
    TabOfIn=[];
    TabOfIn=DataPCA_origin(:,inter(i,:));
    
    [n_norm p_norm] = size(TabOfIn);
    e = ones(n_norm,1);

    StandKine=[StandKine (TabOfIn-min(min(TabOfIn)))./(max(max(TabOfIn))-min(min(TabOfIn)))];
    StandKine_notnormalized=[StandKine_notnormalized TabOfIn];
end


% Mean area under the curve (absement) of the joint kinematic
inter=[56:70];

[a b]=size(inter);
TabOfIn=[];

for i=1:a
    TabOfIn=[];
    TabOfIn=DataPCA_origin(:,inter(i,:));
    
    [n_norm p_norm] = size(TabOfIn);
    e = ones(n_norm,1);

    StandKine=[StandKine (TabOfIn-min(min(TabOfIn)))./(max(max(TabOfIn))-min(min(TabOfIn)))];

end


[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED]=pca(StandKine);%,'Centered',true);%,'VariableWeights','variance');
RecoDataCentered=SCORE*COEFF';

%Normalize COEFF
COEFF;
[n_coeff p_coeff] = size(COEFF);
e = ones(n_coeff,1);
COEFF_Norm=(2.*COEFF-e*(min(COEFF)+max(COEFF)))./(e*max(COEFF)-e*min(COEFF));
COEFF_Norm=(COEFF_Norm).*(COEFF_Norm);



figure
% Master Yoked Switch day 1 and 2
plot(SCORE(1:7,1),SCORE(1:7,2),'or')
hold on 
plot(SCORE(8:14,1),SCORE(8:14,2),'ob')

plot(SCORE(15:21,1),SCORE(15:21,2),'xb')
hold on 
plot(SCORE(22:28,1),SCORE(22:28,2),'xr')
axis([-2 2 -2 2])

for i=1:7
    
    plot([SCORE(i,1) SCORE(i+7,1)],[SCORE(i,2) SCORE(i+7,2)],'k')
    hold on
    plot([SCORE(i+14,1) SCORE(i+21,1)],[SCORE(i+14,2) SCORE(i+21,2)],'k')  
end




%%
% % % % % % % % % % % %     
% % - - 4/ PCA ABLATED POPULATION V.S WT
% % - - - Based on 10 minutes experiment
% % % % % % % % % % % % 

clearvars Data_Day1_DMRT3 Data_Day1_En1 Data_Day1_Ptf1a Data_Day1_Shox2 Data_Day1_Sim1 Data_Day1_Tlx3

% - - Extract ablation label
% 3 -- Ptf1a
% 4 -- Tlx3
% 5 -- DMRT3
% 6 -- En1
% 7 -- Shox2
% 8 -- Sim1

Ablat=2;

iterLbx1=1;
iterPtf1a=1;
iterTlx3=1;
iterDMRT3=1;
iterEn1=1;
iterShox2=1;
iterSim1=1;

[n p]=size(Data);

for i=1:n
            
    if Data(i,2)==3 && Data(i,3)==1 && Data(i,4)==1
        Data_Day1_Ptf1a(iterPtf1a,:)=Data(i,:);
        iterPtf1a=iterPtf1a+1; 
    end
    
    if Data(i,2)==4 && Data(i,3)==1 && Data(i,4)==1
        Data_Day1_Tlx3(iterTlx3,:)=Data(i,:);
        iterTlx3=iterTlx3+1;  
    end
    
    if Data(i,2)==5 && Data(i,3)==1 && Data(i,4)==1
        Data_Day1_DMRT3(iterDMRT3,:)=Data(i,:);
        iterDMRT3=iterDMRT3+1;  
    end
    
    if Data(i,2)==6 && Data(i,3)==1 && Data(i,4)==1
        Data_Day1_Sim1(iterSim1,:)=Data(i,:);
        iterSim1=iterSim1+1;  
    end
    
    if Data(i,2)==7 && Data(i,3)==1 && Data(i,4)==1
        Data_Day1_En1(iterEn1,:)=Data(i,:);
        iterEn1=iterEn1+1;  
    end
    
    if Data(i,2)==8 && Data(i,3)==1 && Data(i,4)==1
        Data_Day1_Shox2(iterShox2,:)=Data(i,:);
        iterShox2=iterShox2+1;  
    end
    
    if Data(i,2)==9 && Data(i,3)==1 && Data(i,4)==1
        Data_Day1_Gad2(iterGad2,:)=Data(i,:);
        iterGad2=iterGad2+1;    
    end
    
end


clearvars StandKine DataPCA SCORE_Ablat COEFF_Ablat
StandKine=[];
StandKine_notnormalized=[];

DataPCA_origin=[Data_Day1_Learner; ...
    Data_Day1_Control; ...
    Data_Day1_Ptf1a; ...
    Data_Day1_Tlx3; ...
    Data_Day1_DMRT3; ...
    Data_Day1_Sim1; ...
    Data_Day1_En1; ...
    Data_Day1_Shox2];

%Export size
[n_norm p_norm] = size(DataPCA_origin);
e = ones(n_norm,1);


% - Binary number related to the position of the foot at the end of the
%   experiment
% - Percentage of time spend below the treshold
inter=[6;7];
[a b]=size(inter);
StandKine=[];
TabOfIn=[];

for i=1:a
    TabOfIn=[];
    TabOfIn=DataPCA_origin(:,inter(i,:));
    
    [n_norm p_norm] = size(TabOfIn);
    e = ones(n_norm,1);

    StandKine=[StandKine (TabOfIn-min(min(TabOfIn)))./(max(max(TabOfIn))-min(min(TabOfIn)))];

end

% Angle 3D 
inter=[23:25];
[a b]=size(inter);
TabOfIn=[];
for i=1:a
    TabOfIn=[];
    TabOfIn=DataPCA_origin(:,inter(i,:));
    
    [n_norm p_norm] = size(TabOfIn);
    e = ones(n_norm,1);

    StandKine=[StandKine (TabOfIn-min(min(TabOfIn)))./(max(max(TabOfIn))-min(min(TabOfIn)))];

end

% Mean position of the joint kinematic 
% Mean amplitude of the joint kinematic
% Mean speed of the joint kinematic
inter=[8:12;13:17;18:22; ...
    26:30;31:35;36:40; ...
    41:45;46:50;51:55];

[a b]=size(inter);
TabOfIn=[];

for i=1:a
    TabOfIn=[];
    TabOfIn=DataPCA_origin(:,inter(i,:));
    
    [n_norm p_norm] = size(TabOfIn);
    e = ones(n_norm,1);

    StandKine=[StandKine (TabOfIn-min(min(TabOfIn)))./(max(max(TabOfIn))-min(min(TabOfIn)))];
    StandKine_notnormalized=[StandKine_notnormalized TabOfIn];
end


% Mean area under the curve (absement) of the joint kinematic
inter=[56:60;61:65;66:70];

[a b]=size(inter);
TabOfIn=[];

for i=1:a
    TabOfIn=[];
    TabOfIn=DataPCA_origin(:,inter(i,:));
    
    [n_norm p_norm] = size(TabOfIn);
    e = ones(n_norm,1);

    StandKine=[StandKine (TabOfIn-min(min(TabOfIn)))./(max(max(TabOfIn))-min(min(TabOfIn)))];

end

clearvars COEFF_Ablat SCORE_Ablat

[COEFF_Ablat, SCORE_Ablat, LATENT_Ablat, TSQUARED_Ablat, EXPLAINED_Ablat]=pca(StandKine);%,'Centered',true);%,'VariableWeights','variance');


interMaster=[1:25];
interYoked=[26:50];
interPtf1a=[51:55];
interTlx3=[56:59];
interDMRT3=[60:64];
interSim1=[65:69];
interEn1=[70:74];
interShox2=[75:79];


%Master
plot(SCORE_Ablat(interMaster,1),SCORE_Ablat(interMaster,2),'.k','MarkerSize',20)
hold on
center = [mean(SCORE_Ablat(interMaster,1)), mean(SCORE_Ablat(interMaster,2))];
center_Master=center;
stdev = [std(SCORE_Ablat(interMaster,1)), std(SCORE_Ablat(interMaster,2))];
stdev_Master=stdev;
axh = gca(); 
llc = center(:)-stdev(:);
wh = stdev(:)*2; 
plot(center(1),center(2),'xk')


% Yoked
plot(SCORE_Ablat(interYoked,1),SCORE_Ablat(interYoked,2),'.k','MarkerSize',10)
center = [mean(SCORE_Ablat(interYoked,1)), mean(SCORE_Ablat(interYoked,2))];
center_Yoked=center;
stdev = [std(SCORE_Ablat(interYoked,1)), std(SCORE_Ablat(interYoked,2))];
stdev_Yoked=stdev;
axh = gca(); 
llc = center(:)-stdev(:);
wh = stdev(:)*2; 
plot(center(1),center(2),'xk')


% %Ptf1a
plot(SCORE_Ablat(interPtf1a,1),SCORE_Ablat(interPtf1a,2),'.b','MarkerSize',20)
center = [mean(SCORE_Ablat(interPtf1a,1)), mean(SCORE_Ablat(interPtf1a,2))];
center_Ptf1a=center;
stdev = [std(SCORE_Ablat(interPtf1a,1)), std(SCORE_Ablat(interPtf1a,2))];
stdev_Ptf1a=stdev;
axh = gca(); 
llc = center(:)-stdev(:);
wh = stdev(:)*2; 
plot(center(1),center(2),'xb')


% %Tlx3
plot(SCORE_Ablat(interTlx3,1),SCORE_Ablat(interTlx3,2),'.r','MarkerSize',20)
center = [mean(SCORE_Ablat(interTlx3,1)), mean(SCORE_Ablat(interTlx3,2))];
center_Tlx3=center;
stdev = [std(SCORE_Ablat(interTlx3,1)), std(SCORE_Ablat(interTlx3,2))];
stdev_Tlx3=stdev;
axh = gca(); 
llc = center(:)-stdev(:);
wh = stdev(:)*2; 
plot(center(1),center(2),'xr')


% %DMRT3
plot(SCORE_Ablat(interDMRT3,1),SCORE_Ablat(interDMRT3,2),'.g','MarkerSize',20)
center = [mean(SCORE_Ablat(interDMRT3,1)), mean(SCORE_Ablat(interDMRT3,2))];
center_DMRT3=center;
stdev = [std(SCORE_Ablat(interDMRT3,1)), std(SCORE_Ablat(interDMRT3,2))];
stdev_DMRT3=stdev;
axh = gca(); 
llc = center(:)-stdev(:);
wh = stdev(:)*2; 
plot(center(1),center(2),'xg')
 


% %SIM1
plot(SCORE_Ablat(interSim1,1),SCORE_Ablat(interSim1,2),'.y','MarkerSize',20)
center = [mean(SCORE_Ablat(interSim1,1)), mean(SCORE_Ablat(interSim1,2))];
center_Sim1=center;
stdev = [std(SCORE_Ablat(interSim1,1)), std(SCORE_Ablat(interSim1,2))];
stdev_Sim1=stdev;
axh = gca(); 
llc = center(:)-stdev(:);
wh = stdev(:)*2; 
plot(center(1),center(2),'xy')


% %EN1
plot(SCORE_Ablat(interEn1,1),SCORE_Ablat(interEn1,2),'.c','MarkerSize',20)
center = [mean(SCORE_Ablat(interEn1,1)), mean(SCORE_Ablat(interEn1,2))];
center_En1=center;
stdev = [std(SCORE_Ablat(interEn1,1)), std(SCORE_Ablat(interEn1,2))];
stdev_En1=stdev;
axh = gca(); 
llc = center(:)-stdev(:);
wh = stdev(:)*2; 
plot(center(1),center(2),'xc')


% %SHOX2
plot(SCORE_Ablat(interShox2,1),SCORE_Ablat(interShox2,2),'.m','MarkerSize',20)
center = [mean(SCORE_Ablat(interShox2,1)), mean(SCORE_Ablat(interShox2,2))];
center_Shox2=center;
stdev = [std(SCORE_Ablat(interShox2,1)), std(SCORE_Ablat(interShox2,2))];
stdev_Shox2=stdev;
axh = gca(); 
llc = center(:)-stdev(:);
wh = stdev(:)*2; 
plot(center(1),center(2),'xm')