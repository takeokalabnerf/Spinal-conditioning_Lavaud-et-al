%% PCA version beg May 2021

% CODE TO USE 2021 09
%%
close all
clc
w = warning ('off','all');

% % % % % % % % % % % %  
% % - - Load ALL the mice
% % % % % % % % % % % %  

clearvars -except UpdatePCAAugust2023 %UpdatePCAallJuly2023v2 %PCAjunegroundtothehip2 %PCAEnCSV %PCAjunegroundtothehip2%PCAaddedParamNewAnalJune %PCAAblat8min %PCAjunegroundtothehip2 %PCAaddedParamNewAnalJune %NewPCAfile092021 %PCAaddedParamNewAnalJune %PCAjunegroundtothehip2 %PCAaddedParamNewAnalJune %PCAjunegroundtothehip2 %PCAaddedParamNewAnalJune %PCAAblat8min %PCAaddedParamApril %PCAEndApril %BookPCAv10
ImportTable= UpdatePCAAugust2023; %UpdatePCAallJuly2023v2; %PCAjunegroundtothehip2; %PCAEnCSV; %PCAjunegroundtothehip2; %NewPCAfile092021; % PCAaddedParamNewAnalJune; %PCAaddedParamNewAnalJune; %PCAaddedParamApril; %PCAEndApril; %BookPCAv10;
Data=table2array(ImportTable);
[n p]=size(Data);

% %%
% 
% TAB=95;
% Data1=table2array(PCAjunegroundtothehip2);
% compare(1:100,1)=Data1(1:100,TAB);
% Data2=table2array(UpdatePCAallJuly2023);
% compare(1:100,2)=Data2(1:100,TAB);


%%
% % % % % % % % % % % %  
% % - - Exclusion of mice not performing well
% % - - Exclusion rules have been selected to be:  % <4*SEM %     
% % - - Selection of the mice to display. Day 1 M/Y no modification
% % % % % % % % % % % %  


iter=1;
iterY=1;

for i=1:n
    
    if Data(i,1)==1 || Data(i,1)==2 || Data(i,1)==5 || Data(i,1)==12 || Data(i,1)==13 || Data(i,1)==14 || Data(i,1)==15 || Data(i,1)==19 || Data(i,1)==27 %Experience 1 5 12 13 14 15 19 27
    	if (Data(i,4)==1) %|| (Data(i,4)<5 && Data(i,1)==27) || (Data(i,4)<5 && Data(i,1)==19) || (Data(i,4)<4 && Data(i,1)==15))
            if Data(i,3)==1
            
            Data_Day1_Master(iter,:)=Data(i,:);
            iter=iter+1;
            else
                
           	Data_Day1_Yoked(iterY,:)=Data(i,:);
            iterY=iterY+1;

            end
            
    end
    
    end

end


[n_Data_Day1_Master p_Data_Day1_Master]=size(Data_Day1_Master);
FailMatrice=zeros(n_Data_Day1_Master,p_Data_Day1_Master);

for i=1:p_Data_Day1_Master
        
    M=mean(Data_Day1_Master(:,i));
    S=std(Data_Day1_Master(:,i));
    
    for j=1:n_Data_Day1_Master
          
        if Data_Day1_Master(j,i)<(M-2*S) || Data_Day1_Master(j,i)>(M+2*S)         
            FailMatrice(j,i)=1;
        end
        
        if Data_Day1_Master(j,i)<(M-4*S) || Data_Day1_Master(j,i)>(M+4*S)         
            FailMatrice(j,i)=2;
        end
        
     end
    
end

%
% % % % % % % % % % % %  
% % - - Exclusion if sum > 20 % -> 13 14 28
% % % % % % % % % % % %  

FailMatriceSum=sum(FailMatrice');
a=1;
% for i=n_Data_Day1_Master:-1:1
%     
%     if FailMatriceSum(i)>=20
%         Data_Day1_Master(i,:)=[];
%         Data_Day1_Yoked(i,:)=[];
%         
%         Master_That_Fail(a)=i;
%         a=a+1;
%     end
%     
% end
%Data_Day1_Master(30,:)=[]
%Data_Day1_Yoked(30,:)=[]


Data_Day1_Master(29,:)=[]
Data_Day1_Yoked(29,:)=[]

Data_Day1_Master(20,:)=[]
Data_Day1_Yoked(20,:)=[]

Data_Day1_Master(15,:)=[]
Data_Day1_Yoked(15,:)=[]

Data_Day1_Master(6,:)=[]
Data_Day1_Yoked(6,:)=[]

Data_Day1_Master(5,:)=[]
Data_Day1_Yoked(5,:)=[]



% Master distribution
% 1 - 1 2 3 4 
% 2 - 5 6 
% 5 - 7 8 9
% 12 - 10 11
% 13 - 12 13 14 15
% 14 - 16 17 18 19 20
% 15 - 21 22 23
% 19 - 24 25 26
% 27 - 27 28 29 30

% Remove :
% 5 6 15 20 29
% 2-1 / 2-2 / 13-4 / 14-5 / 27-3
% Fail = [30,16,15,11,6,5]
% Fail_new = [6,11,15,16,30];

% 12 - 1 3 
% 13 - 1 2

% Master distribution
% 1 - 1 2 3 4               4   4 
% 2 -                       0   4
% 5 - 5 6 7                 3   7
% 12 - 8 9                  2   9
% 13 - 10 11                2   11
% 14 - 12 13 14 15 16       5   16
% 15 - 17 18 19             3   19
% 19 - 20 21 22             3   22
% 27 - 23 24 25             3   25

% Total of : 25  - 16 good / 9 bad normally
% % % % % % % % % % % %       
% % - - 1/ Sort DATA_ALL By Included/Excluded
% % % % % % % % % % % %  

clearvars Data_All
iter=1;

% for i=1:n 
%     
%     %Excluse data from exp 2 / Exclude 12.2 13.3 13.4 27.3 
%     %+ if not above 15 : 14.3 12.1
%     if Data(i,1)==1 || Data(i,1)==5 || Data(i,1)==12 || Data(i,1)==13 || Data(i,1)==14 || Data(i,1)==15 || Data(i,1)==19 || Data(i,1)==27 %Experience 1 5 12 13 14 15 19 27
%         if (Data(i,4)>0)
%         %Data to exclude BECAUSE STD EXCLUSION
%             if (Data(i,1)==13 && (Data(i,5)==3 || Data(i,5)==4)) || (Data(i,1)==12 && Data(i,5)==2) %|| (Data(i,1)==14 && Data(i,5)==3) || (Data(i,1)==12 && Data(i,5)==1) || (Data(i,1)==27 && Data(i,5)==4)
%                 continue
%             else
%                 
%             Data_All(iter,:)=Data(i,:);
%             iter=iter+1;
%             
%             end
%             
%         end
%     end
%     
% end
    



for i=1:n 
    
    %Excluse data from exp 2 / Exclude 12.2 13.3 13.4 27.3 
    %+ if not above 15 : 14.3 12.1
    if Data(i,1)==1 || Data(i,1)==5 || Data(i,1)==12 || Data(i,1)==13 || Data(i,1)==14 || Data(i,1)==15 || Data(i,1)==19 || Data(i,1)==27 %Experience 1 5 12 13 14 15 19 27
        if (Data(i,4)>0)
        %Data to exclude BECAUSE STD EXCLUSION
            if (Data(i,1)==2 && (Data(i,5)==2)) || (Data(i,1)==2 && Data(i,5)==1) || (Data(i,1)==12 && Data(i,5)==2) || (Data(i,1)==13 && Data(i,5)==3) || (Data(i,1)==13 && Data(i,5)==4) || (Data(i,1)==27 && Data(i,5)==4)
                continue
            else
                
            Data_All(iter,:)=Data(i,:);
            iter=iter+1;
            
            end
            
        end
    end
    
end




% % % % % % % % % % % %       
% % - - 2/ Sort By Experiment type
% % % % % % % % % % % %  

iterMSwitch=1;
iterYSwitch=1;
iterMSwitch2=1;
iterYSwitch2=1;
iterMSwitch3=1;
iterYSwitch3=1;
iterMSwitch4=1;
iterYSwitch4=1;
iterMSwitch5=1;
iterYSwitch5=1;
iterMSwitch6=1;
iterYSwitch6=1;

iterMSwitch48=1;
iterYSwitch48=1;
iterMSwitch48_1=1;
iterYSwitch48_1=1;

iterMSwitch72=1;
iterYSwitch72=1;
iterMSwitch72_1=1;
iterYSwitch72_1=1;

iterMSwitch48_310=1;
iterYSwitch48_310=1;
iterMSwitch48_310_1=1;
iterYSwitch48_310_1=1;
iterMSwitch48_310_2=1;
iterYSwitch48_310_2=1;
iterMSwitch48_310_3=1;
iterYSwitch48_310_3=1;

iterMSwitch48_410=1;
iterYSwitch48_410=1;
iterMSwitch48_410_1=1;
iterYSwitch48_410_1=1;
iterMSwitch48_410_2=1;
iterYSwitch48_410_2=1;
iterMSwitch48_410_3=1;
iterYSwitch48_410_3=1;
iterMSwitch48_410_4=1;
iterYSwitch48_410_4=1;
iterMSwitch48_410_5=1;
iterYSwitch48_410_5=1;

[p n]=size(Data_All)

for i=1:p
    
    %Cond 1
    if Data_All(i,2)==1 && Data_All(i,4)==1
        if Data_All(i,3)==1
        Data_Day1_Switch_Master(iterMSwitch,:)=Data_All(i,:);
        iterMSwitch=iterMSwitch+1;
        else 
      	Data_Day1_Switch_Yoked(iterYSwitch,:)=Data_All(i,:);
        iterYSwitch=iterYSwitch+1;
        end  
    end
    
    if Data_All(i,2)==1 && Data_All(i,4)==2
        if Data_All(i,3)==1
        Data_Day2_Switch_Master(iterMSwitch2,:)=Data_All(i,:);
        iterMSwitch2=iterMSwitch2+1;
        else 
      	Data_Day2_Switch_Yoked(iterYSwitch2,:)=Data_All(i,:);
        iterYSwitch2=iterYSwitch2+1;
        end  
    end
    
    if Data_All(i,2)==1 && Data_All(i,4)==3
        if Data_All(i,3)==1
        Data_Day3_Switch_Master(iterMSwitch3,:)=Data_All(i,:);
        iterMSwitch3=iterMSwitch3+1;
        else 
      	Data_Day3_Switch_Yoked(iterYSwitch3,:)=Data_All(i,:);
        iterYSwitch3=iterYSwitch3+1;
        end  
    end
    
    if Data_All(i,2)==1 && Data_All(i,4)==4
        if Data_All(i,3)==1
        Data_Day4_Switch_Master(iterMSwitch4,:)=Data_All(i,:);
        iterMSwitch4=iterMSwitch4+1;
        else 
      	Data_Day4_Switch_Yoked(iterYSwitch4,:)=Data_All(i,:);
        iterYSwitch4=iterYSwitch4+1;
        end  
    end
    
    if Data_All(i,2)==1 && Data_All(i,4)==5
        if Data_All(i,3)==1
        Data_Day5_Switch_Master(iterMSwitch5,:)=Data_All(i,:);
        iterMSwitch5=iterMSwitch5+1;
        else 
      	Data_Day5_Switch_Yoked(iterYSwitch5,:)=Data_All(i,:);
        iterYSwitch5=iterYSwitch5+1;
        end  
    end
    
   if Data_All(i,2)==1 && Data_All(i,4)==6
        if Data_All(i,3)==1
        Data_Day6_Switch_Master(iterMSwitch6,:)=Data_All(i,:);
        iterMSwitch6=iterMSwitch6+1;
        else 
      	Data_Day6_Switch_Yoked(iterYSwitch6,:)=Data_All(i,:);
        iterYSwitch6=iterYSwitch6+1;
        end  
    end
    
    %Cond 10 - Switch 48h
    if Data_All(i,2)==10 && Data_All(i,4)==1
        if Data_All(i,3)==1
        Data_Day1_Switch_48_Master(iterMSwitch48,:)=Data_All(i,:);
        iterMSwitch48=iterMSwitch48+1;
        else 
      	Data_Day1_Switch_48_Yoked(iterYSwitch48,:)=Data_All(i,:);
        iterYSwitch48=iterYSwitch48+1;
        end  
    end
    
    if Data_All(i,2)==10 && Data_All(i,4)==2
        if Data_All(i,3)==1
        Data_Day2_Switch_48_Master(iterMSwitch48_1,:)=Data_All(i,:);
        iterMSwitch48_1=iterMSwitch48_1+1;
        else 
      	Data_Day2_Switch_48_Yoked(iterYSwitch48_1,:)=Data_All(i,:);
        iterYSwitch48_1=iterYSwitch48_1+1;
        end  
    end
    

    %Cond 100 - Switch 72h
    if Data_All(i,2)==100 && Data_All(i,4)==1
        if Data_All(i,3)==1
        Data_Day1_Switch_72_Master(iterMSwitch72,:)=Data_All(i,:);
        iterMSwitch72=iterMSwitch72+1;
        else 
      	Data_Day1_Switch_72_Yoked(iterYSwitch72,:)=Data_All(i,:);
        iterYSwitch72=iterYSwitch72+1;
        end  
    end
    
    if Data_All(i,2)==100 && Data_All(i,4)==2
        if Data_All(i,3)==1
        Data_Day2_Switch_72_Master(iterMSwitch72_1,:)=Data_All(i,:);
        iterMSwitch72_1=iterMSwitch72_1+1;
        else 
      	Data_Day2_Switch_72_Yoked(iterYSwitch72_1,:)=Data_All(i,:);
        iterYSwitch72_1=iterYSwitch72_1+1;
        end  
    end
    
    %Cond 310 - Switch 48h after 3days training
    if Data_All(i,2)==310 && Data_All(i,4)==1
        if Data_All(i,3)==1
        Data_Day1_TrainSwitch_Master(iterMSwitch48_310,:)=Data_All(i,:);
        iterMSwitch48_310=iterMSwitch48_310+1;
        else 
      	Data_Day1_TrainSwitch_Yoked(iterYSwitch48_310,:)=Data_All(i,:);
        iterYSwitch48_310=iterYSwitch48_310+1;
        end  
    end
    
    if Data_All(i,2)==310 && Data_All(i,4)==2
        if Data_All(i,3)==1
        Data_Day2_TrainSwitch_Master(iterMSwitch48_310_1,:)=Data_All(i,:);
        iterMSwitch48_310_1=iterMSwitch48_310_1+1;
        else 
      	Data_Day2_TrainSwitch_Yoked(iterYSwitch48_310_1,:)=Data_All(i,:);
        iterYSwitch48_310_1=iterYSwitch48_310_1+1;
        end  
    end
    
    if Data_All(i,2)==310 && Data_All(i,4)==3
        if Data_All(i,3)==1
        Data_Day3_TrainSwitch_Master(iterMSwitch48_310_2,:)=Data_All(i,:);
        iterMSwitch48_310_2=iterMSwitch48_310_2+1;
        else 
      	Data_Day3_TrainSwitch_Yoked(iterYSwitch48_310_2,:)=Data_All(i,:);
        iterYSwitch48_310_2=iterYSwitch48_310_2+1;
        end  
    end
    
    if Data_All(i,2)==310 && Data_All(i,4)==4
        if Data_All(i,3)==1
        Data_Day4_TrainSwitch_Master(iterMSwitch48_310_3,:)=Data_All(i,:);
        iterMSwitch48_310_3=iterMSwitch48_310_3+1;
        else 
      	Data_Day4_TrainSwitch_Yoked(iterYSwitch48_310_3,:)=Data_All(i,:);
        iterYSwitch48_310_3=iterYSwitch48_310_3+1;
        end  
    end

    
    %Cond 410 - Switch 48h after 4days training
    if Data_All(i,2)==410 && Data_All(i,4)==1
        if Data_All(i,3)==1
        Data_Day1_Train4Switch_Master(iterMSwitch48_410,:)=Data_All(i,:);
        iterMSwitch48_410=iterMSwitch48_410+1;
        else 
      	Data_Day1_Train4Switch_Yoked(iterYSwitch48_410,:)=Data_All(i,:);
        iterYSwitch48_410=iterYSwitch48_410+1;
        end  
    end
    
    if Data_All(i,2)==410 && Data_All(i,4)==2
        if Data_All(i,3)==1
        Data_Day2_Train4Switch_Master(iterMSwitch48_410_1,:)=Data_All(i,:);
        iterMSwitch48_410_1=iterMSwitch48_410_1+1;
        else 
      	Data_Day2_Train4Switch_Yoked(iterYSwitch48_410_1,:)=Data_All(i,:);
        iterYSwitch48_410_1=iterYSwitch48_410_1+1;
        end  
    end
    
    if Data_All(i,2)==410 && Data_All(i,4)==3
        if Data_All(i,3)==1
        Data_Day3_Train4Switch_Master(iterMSwitch48_410_2,:)=Data_All(i,:);
        iterMSwitch48_410_2=iterMSwitch48_410_2+1;
        else 
      	Data_Day3_Train4Switch_Yoked(iterYSwitch48_410_2,:)=Data_All(i,:);
        iterYSwitch48_410_2=iterYSwitch48_410_2+1;
        end  
    end
    
    if Data_All(i,2)==410 && Data_All(i,4)==4
        if Data_All(i,3)==1
        Data_Day4_Train4Switch_Master(iterMSwitch48_410_3,:)=Data_All(i,:);
        iterMSwitch48_410_3=iterMSwitch48_410_3+1;
        else 
      	Data_Day4_Train4Switch_Yoked(iterYSwitch48_410_3,:)=Data_All(i,:);
        iterYSwitch48_410_3=iterYSwitch48_410_3+1;
        end  
    end
    
    if Data_All(i,2)==410 && Data_All(i,4)==5
        if Data_All(i,3)==1
        Data_Day5_Train4Switch_Master(iterMSwitch48_410_4,:)=Data_All(i,:);
        iterMSwitch48_410_4=iterMSwitch48_410_4+1;
        else 
      	Data_Day5_Train4Switch_Yoked(iterYSwitch48_410_4,:)=Data_All(i,:);
        iterYSwitch48_410_4=iterYSwitch48_410_4+1;
        end  
    end
    
    if Data_All(i,2)==410 && Data_All(i,4)==6
        if Data_All(i,3)==1
        Data_Day6_Train4Switch_Master(iterMSwitch48_410_5,:)=Data_All(i,:);
        iterMSwitch48_410_5=iterMSwitch48_410_5+1;
        else 
      	Data_Day6_Train4Switch_Yoked(iterYSwitch48_410_5,:)=Data_All(i,:);
        iterYSwitch48_410_5=iterYSwitch48_410_5+1;
        end  
    end
    
end

%
% % % % % % % % 
%% % - -PCA TO USE WITH KINEMATIC STANDARDISATION - STEADY
% % % % % % % %


clearvars StandKine DataPCA SCORE COEFF

StandKine=[];

% Steady - 74 param
Param=[6 7   8:22 41:55 80:94 95:109   23:25]; %56:64
% Dyn - 72 param
%Param=[110:124 125:139 155:163 179:193 194:196 197:211];

%8:22 MEAN
%23:25 ANGLE
%26:40 STD
%41:55 AMPLITUDE
%56:64 PARTICIPATION
%65:79 LINEAR EVOLUTION
%80:94 SPEED
%95:109 AREA

%1-2 Experiment param
%3-17 Mean
%18-32 Amplitude

DataPCA_origin=[Data_Day1_Master(:,:);Data_Day1_Yoked(:,:)];

DataPCA=DataPCA_origin(:,Param);

%

%Export size
[n_norm p_norm] = size(DataPCA);
e = ones(n_norm,1);


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



inter=[8:12;13:17;18:22; ...
    41:45;46:50;51:55; ...
    65:69;70:74;75:79];


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



inter=[80:84;85:89;90:94];


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


% inter=[56:58;59:61;62:64];
% [a b]=size(inter);
% TabOfIn=[];
% 
% for i=1:a
%     TabOfIn=[];
%     TabOfIn=DataPCA_origin(:,inter(i,:));
%     
%     [n_norm p_norm] = size(TabOfIn);
%     e = ones(n_norm,1);
% 
%     StandKine=[StandKine (TabOfIn-min(min(TabOfIn)))./(max(max(TabOfIn))-min(min(TabOfIn)))];
% 
% end





[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED]=pca(StandKine);%,'Centered',true);%,'VariableWeights','variance');
RecoDataCentered=SCORE*COEFF';

%Normalize COEFF
COEFF;
[n_coeff p_coeff] = size(COEFF);
e = ones(n_coeff,1);

factor_pos= 1/max(COEFF(:,1));
factor_neg= -1/min(COEFF(:,1));
    
for i=1:length(COEFF)
   
    if COEFF(i,1)>0 
        
        COEFF_Norm(i,1)=COEFF(i,1)*factor_pos;
        
    else
        
        COEFF_Norm(i,1)=COEFF(i,1)*factor_neg;
        
    end
    
end
        

COEFF_Norm=(COEFF-e*min(COEFF))./(e*max(COEFF)-e*min(COEFF));


COEFF_Norm=(2.*COEFF-e*(min(COEFF)+max(COEFF)))./(e*max(COEFF)-e*min(COEFF));
%COEFF_Norm=(COEFF_Norm).*(COEFF_Norm);

day1_noretrain=1:25;
numel(day1_noretrain);
scatter3(SCORE(1:numel(day1_noretrain),1),SCORE(1:numel(day1_noretrain),2),SCORE(1:numel(day1_noretrain),3),'r')
hold on
scatter3(SCORE((1:numel(day1_noretrain))+numel(day1_noretrain),1),SCORE((1:numel(day1_noretrain))+numel(day1_noretrain),2),SCORE((1:numel(day1_noretrain))+numel(day1_noretrain),3),'b')

Display_1C_Master=SCORE(1:25,1:3);
Display_1C_Yoked=SCORE(26:50,1:3);

% Create some random data
y1 = SCORE(1:25,1); %normrnd(s(1).*x,1);
y2 = SCORE(1:25,2); %normrnd(s(2).*x,1);
y1_prime = SCORE(26:50,1); %normrnd(s(1).*x,1);
y2_prime = SCORE(26:50,2); %normrnd(s(2).*x,1);
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

% Plot the original data
%plot(data(:,1), data(:,2), '.b');
%plot(data_prime(:,1), data_prime(:,2), '.r');

% mindata = min(min(data));
% maxdata = max(max(data));
% Xlim([mindata-3, maxdata+3]);
% Ylim([mindata-3, maxdata+3]);
hold on;

% Plot the eigenvectors
%quiver(X0, Y0, largest_eigenvec(1)*sqrt(largest_eigenval), largest_eigenvec(2)*sqrt(largest_eigenval), '-m', 'LineWidth',2);
%quiver(X0, Y0, smallest_eigenvec(1)*sqrt(smallest_eigenval), smallest_eigenvec(2)*sqrt(smallest_eigenval), '-g', 'LineWidth',2);
hold on;

% Set the axis labels
hXLabel = xlabel('x');
hYLabel = ylabel('y');


%
% % % % % % % % 
%% % - -PCA TO USE WITH KINEMATIC STANDARDISATION - DYN
% % % % % % % %

clearvars StandKine DataPCA 

StandKine=[];
DataPCA_origin=[Data_Day1_Master(:,:);Data_Day1_Yoked(:,:)];

% Steady - 74 param
%Param=[6 7 8:22 41:55 80:94 95:109  56:64 23:25];
% Dyn - 72 param
% Param=[110:124 125:139 155:163 179:193 194:196 197:211];

%110:124 AMPLITUDE SPIKE
%125:139 SPEED
%140:154 ACC
%155:163 PARTICIPATION
%164:178 MAX MIN ENDPOINT
%179:193 AREA
%194:196 ANGLE
%197:211 MEAN 



% inter=[110:114;115:119;120:124; ...
%     125:129;130:134;135:139; ...
%     179:183;184:188;189:193;...
%     197:201;202:206;207:211];

inter=[182:186;187:191;192:196; ...
    95:99;100:104;105:109; ...
    110:114;115:119;120:124;...
    164:168;169:173;174:178];

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

% inter=[155:157;158:160;161:163];
% [a b]=size(inter);
% TabOfIn=[];
% 
% for i=1:a
%     TabOfIn=[];
%     TabOfIn=DataPCA_origin(:,inter(i,:));
%     
%     [n_norm p_norm] = size(TabOfIn);
%     e = ones(n_norm,1);
% 
%     StandKine=[StandKine (TabOfIn-min(min(TabOfIn)))./(max(max(TabOfIn))-min(min(TabOfIn)))];
% 
% end


%inter=[194:196];
inter=[179:181];
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

% Plot the original data
%plot(data(:,1), data(:,2), '.b');
%plot(data_prime(:,1), data_prime(:,2), '.r');

% mindata = min(min(data));
% maxdata = max(max(data));
% Xlim([mindata-3, maxdata+3]);
% Ylim([mindata-3, maxdata+3]);
hold on;

% Plot the eigenvectors
%quiver(X0, Y0, largest_eigenvec(1)*sqrt(largest_eigenval), largest_eigenvec(2)*sqrt(largest_eigenval), '-m', 'LineWidth',2);
%quiver(X0, Y0, smallest_eigenvec(1)*sqrt(smallest_eigenval), smallest_eigenvec(2)*sqrt(smallest_eigenval), '-g', 'LineWidth',2);
hold on;

% Set the axis labels
hXLabel = xlabel('x');
hYLabel = ylabel('y');





% % %
% % 
% % % inter=[110:114;115:119;120:124; ... %Amp
% % %     125:129;130:134;135:139; ... %Speed
% % %     179:183;184:188;189:193;... %Area
% % %     197:201;202:206;207:211]; %Mean
% % 
% % 
% % % Dyn - 72 param
% % %Param=[110:124 125:139 155:163 179:193 194:196 197:211];
% % 
% % %110:124 AMPLITUDE SPIKE
% % %125:139 SPEED
% % %140:154 ACC
% % %155:163 PARTICIPATION
% % %179:193 AREA
% % %194:196 ANGLE
% % %197:211 MEAN 
% % 
% % 
% % CoefToStudy=COEFF(:,1);
% % Joint_2=[];
% % Joint_3=[];
% % Joint_5=[];
% % 
% % % Mean / Amplitude / Speed / Area / Partici / Angle 
% % 
% % % % Dynamic
% % % Joint_5=[ ...
% % %     46:50; 51:55; 56:60;
% % %     1:5; 6:10; 11:15; ...
% % %     16:20; 21:25; 26:30; ...
% % %     31:35; 36:40; 41:45];
% % % 
% % % Joint_3=[61:63; 64:66; 67:69; 70:72];
% % % 
% % % Joint_2=0;
% % 
% % % Steady
% % % inter=[8:12;13:17;18:22;41:45;46:50; ...
% % %     51:55;80:84;85:89;90:94; ...
% % %     95:99;100:104;105:109];
% % 
% % Joint_2=[1:2];
% % 
% % Joint_5=[ ...
% %     1:5; 6:10; 11:15; ...
% %     16:20; 21:25; 26:30; ...
% %     31:35; 36:40; 41:45; ...
% %     46:50; 51:55; 56:60]+2;
% % 
% % Joint_3=[61:63; 64:66; 67:69; 70:72]+2;
% % 
% % 
% % %Joint5
% % [a b]=size(Joint_5);
% % CoeffVisu=[];
% % Col=1;
% % 
% % for i=1:a
% % 
% %     CoeffVisu_iter=sum(abs(COEFF(Joint_5(i,:),Col)));
% %     
% %     CoeffVisu=[CoeffVisu CoeffVisu_iter];
% %     
% % end
% %     
% % %Joint3
% % [a b]=size(Joint_3);
% % for i=1:a
% % 
% %     CoeffVisu_iter=sum(abs(COEFF(Joint_3(i,:),Col)));
% %     
% %     CoeffVisu=[CoeffVisu CoeffVisu_iter];
% %     
% % end
% %     
% % %Joint2
% % clearvars CoeffVisu_iter
% % [a b]=size(Joint_2);
% % if numel(Joint_2)==2
% %     for i=1:a
% %         
% %         CoeffVisu_iter=sum(abs(COEFF(Joint_2(i,:),Col)));
% %         
% %         CoeffVisu=[CoeffVisu CoeffVisu_iter];
% %         
% %     end
% % end


%%
% % % % % % % % % % % %     
% % - - COMPUTE DATA FOR GOOD/BAD MASTER + STIMULATION 
% % % % % % % % % % % % 
% Steady for 1E 
% Dynamic for 1C

iterM=1;
iterY=1;
[n_Data_All p_Data_All]=size(Data_All)


clearvars StimulationMaster
for i=1:25
    StimulationMaster(iterM,1:20)=Data_Day1_Master(i,212:231);
    iterM=iterM+1;
end
% 
% for i=1:n_Data_All
%    
%     if Data_All(i,1)==1 || Data_All(i,1)==5 || Data_All(i,1)==12 || Data_All(i,1)==13 || Data_All(i,1)==14 || Data_All(i,1)==15 || Data_All(i,1)==19 || Data_All(i,1)==27 %Experience 1 5 12 13 14 15 19 27
%             if Data_All(i,4)==1 %Day 1
%                 if Data_All(i,3)==1 % status 1
%                     
%                 %Display_1C_Master(iterM,1:3)=SCORE(i,1:3);
%                 StimulationMaster(iterM,1:20)=Data_All(i,212:231);
%                 %TESTPCASORT(iterM,1:p_Data_All)=Data_All(i,1:p_Data_All);
%                 %MasterDataPCA2Compare(iterM,:)=Data_All(i,Param);
%                 iterM=iterM+1;
%                 
%                else
%                     %Display_1C_Yoked(iterY,1:3)=SCORE(i,1:3);
%                     %WhichIsYoked(iterY)=i;
%                     %YokedDataPCA2Compare(iterY,:)=Data_All(i,Param);
% 
%                     iterY=iterY+1;
%                 end
%             
%             end
%             
%     end
%     
%     
% end

%
[n_Data_Master p_Data_Master]=size(Display_1C_Yoked);

for i=1:10
    for j=1:n_Data_Master
        
    TotalStim(j,i)=StimulationMaster(j,i)*60*60/100;
    TotalStimNorm(j,i)=(StimulationMaster(j,i)/StimulationMaster(j,1))*60*60/100*100/36;
    
    end
    
    TotalStim_Sum_ByTB(i)=sum(TotalStim(:,i));
end

for i=1:n_Data_Master
    TotalStim_Sum_ByMice(i,1)=sum(TotalStim(i,:));
end
 
%Total amount of stimulation
Sca_TSSBM=[ones(n_Data_Master,1) TotalStim_Sum_ByMice];   

%Time touching the threshold 
for i=1:10
    for j=1:n_Data_Master    
    TotalTouch(j,i)=StimulationMaster(j,i+10);
    end
    
    TotalTouch_Sum_ByTB(i)=sum(TotalTouch(:,i));
end


for i=1:n_Data_Master
    TotalTouch_Sum_ByMice(i,1)=sum(TotalTouch(i,:));
end



%%
% % % % % % % % % % % %     
% % - - CORRELATION + GOOD vs. BAD MASTER 
% % % % % % % % % % % % 
%sqeuclidean distance

%X=Display_1C_Master(1:n_Data_Master,1:2);
X=[Display_1C_Master(1:25,1) TotalStim_Sum_ByMice(1:25)];
%X=[yfit(1:25) ones(25,1)];

clearvars Silh

for i=2 %NumberOfClusters

    opts = statset('Display','final');
	[idx,C,sumd,D] = kmeans(X,i,'Distance','sqeuclidean',...
    'Replicates',100,'Options',opts);

    [silh,h] = silhouette(X,idx,'sqeuclidean');

    Silh(i)=mean(silh)

end

%mean(silh)

figure; 
plot(Display_1C_Master(idx==1,1),TotalStim_Sum_ByMice(idx==1,1),'r.','MarkerSize',12)
hold on
plot(Display_1C_Master(idx==2,1),TotalStim_Sum_ByMice(idx==2,1),'b.','MarkerSize',12)
%plot(X(idx==3,1),X(idx==3,1),'g.','MarkerSize',12)

[t i]=ttest2(TotalStim_Sum_ByMice(idx==1,1),TotalStim_Sum_ByMice(idx==2,1))
[t i]=ttest2(Display_1C_Master(idx==1,1),Display_1C_Master(idx==2,1))



% Good 
goodmaster=find(idx==2);
% Bad 
badmaster=find(idx==1);

%%
% % % % % % % % % % % %     
% % - - ANALYSIS DIFFERENCE GOOD vs. BAD MASTER 
% % % % % % % % % % % % 


%Param=[6 7 8:22 41:55 80:94 95:109 56:64 23:25];
%Param=[110:124 125:139 155:163 179:193 194:196 197:211];

%Param=[6 7 8:22 41:55 80:94 95:109 56:64 23:25 110:124 125:139 155:163 179:193 194:196 197:211];

% All param
% Param=[6 7 8:10 13:15 18:20 ...
%     23:25 ...
%     41:43 46:48 51:53 ...
%     56:58 59:61 62:64 ...
%     80:82 85:87 90:92 ...
%     95:97 100:102 105:107 ...
%     110:112 115:117 120:122 ...
%     125:127 130:132 135:137 ...
%     179:181 184:186 189:191 ...
%     194:196 ...
%     197:199 202:204 207:209]; %155:163 ...

% Only dyna
Param=[ ...
    197:199 202:204 207:209 ...
    95:97 100:102 105:107 ...
    110:112 115:117 120:122 ...
    125:127 130:132 135:137 ...
    179:181 184:186 189:191 ...
    194:196]; %155:163 ...

a=1;
b=1;

clearvars Bad Good
for i=1:25
    
    if idx(i)==2
    Bad(a)=i;
    a=a+1;
    else
    Good(b)=i;
    b=b+1;
    end
    
end

%  4  7  8  11  14  18  20  21  24


%Bad=[4 7 8 11 14 18 20 21 24];
%Good=[1 2 3 5 6 9 10 12 13 15 16 17 19 22 23 25];

iter=1;
iter_signi=1;

figure(1)
clearvars TTEST_Abl TabMean1 TabMean5 TabSTD1 TabSTD5 TabDivide 
clearvars TabMean1_Signi TabMean5_Signi TabSTD1_Signi TabSTD5_Signi Tab_Signi_Impact TabSigniWho
clearvars TabData1 TabData5
clc

clearvars Tabreal1 Tabreal5
hold on
for i=Param
    
    A=(Data_Day1_Master(Bad,i));
    B=(Data_Day1_Master(Good,i));
    
    TabToScale2=([A;B]); % [9*Bad 16*Good]
    %TabToScale=(TabToScale2);
    TabToScale=(TabToScale2-min(TabToScale2)) / (max(TabToScale2)-min(TabToScale2));
    
    TabMean1(iter)=mean(TabToScale(1:length(A),1)); %Bad
    TabMean5(iter)=mean(TabToScale((length(A)+1):(length(B)+length(A)),1)); %Good
    

    TabSTD1(iter)=std(TabToScale(1:length(A),1))/sqrt(length(A)-1);
    TabSTD5(iter)=std(TabToScale((length(A)+1):(length(B)+length(A)),1))/sqrt(length(B)-1);

    TabDivide(iter)=(TabMean1(iter)/TabMean5(iter)); %Bad/Good
    
    [a b]=ttest2(TabToScale(1:length(A),1),TabToScale((length(A)+1):(length(B)+length(A)),1));
    TTEST_Abl(iter)=(b);
    
    if b<0.05
        
    Tabreal1(iter_signi)=mean(abs(Data_Day1_Master(Bad,i)));
    Tabreal5(iter_signi)=mean(abs(Data_Day1_Master(Good,i)));
    
    
    TabMean1_Signi(iter_signi)=TabMean1(iter);
    TabMean5_Signi(iter_signi)=TabMean5(iter);
    TabSTD1_Signi(iter_signi)=TabSTD1(iter);
    TabSTD5_Signi(iter_signi)=TabSTD5(iter);
    Tab_Signi_Impact(iter_signi)=b;
    TabSigniWho(iter_signi)=NewPCAfile092021.Properties.VariableNames(([i]));
    
    TabData1(iter_signi,1:length(A))=TabToScale(1:length(A),1);
    TabData5(iter_signi,1:numel((length(A)+1):(length(B)+length(A))))=TabToScale((length(A)+1):(length(B)+length(A)),1);

    TabData1_notnormal(iter_signi,1:length(A))=A;
    TabData5_notnormal(iter_signi,1:numel((length(A)+1):(length(B)+length(A))))=B;
    
    iter_signi=iter_signi+1;
    
    end
    
    iter=iter+1;
    
%     0.05
%     0.01
%     0.001
end


errorbar(TabMean1_Signi,TabSTD1_Signi)
errorbar(TabMean5_Signi,TabSTD5_Signi)

[a b]=size(TabData1)
for i=1:b
    
    plot(i*ones(1,length(A))-0.2,TabData1(i,:),'or')
    plot(0.2+i*ones(1,numel((length(A)+1):(length(B)+length(A)))),TabData5(i,:),'ob')

end

%axis([-1 length(TabMean1_Signi)+1 0 1])

figure(2)
Impact = [TabMean1_Signi 0];
theta = linspace(0,2*pi,(length(Impact)));
polarplot(theta,Impact,'-r');
polarplot(theta,Impact,'or');
hold on

Impact = [TabMean1_Signi-TabSTD1_Signi 0];
theta = linspace(0,2*pi,(length(Impact)));
polarplot(theta,Impact,'--r');
polarplot(theta,Impact,'xr');

Impact = [TabMean1_Signi+TabSTD1_Signi 0];
theta = linspace(0,2*pi,(length(Impact)));
polarplot(theta,Impact,'--r');
polarplot(theta,Impact,'xr');

Impact = [TabMean5_Signi 0];
theta = linspace(0,2*pi,(length(Impact)));
polarplot(theta,Impact,'ob');
polarplot(theta,Impact,'-b');
hold on

Impact = [TabMean5_Signi-TabSTD5_Signi 0];
theta = linspace(0,2*pi,(length(Impact)));
polarplot(theta,Impact,'--b');
polarplot(theta,Impact,'xb');

Impact = [TabMean5_Signi+TabSTD5_Signi 0];
theta = linspace(0,2*pi,(length(Impact)));
polarplot(theta,Impact,'--b');
polarplot(theta,Impact,'xb');

ax=gca
ax.ThetaTick = [0:(360/length(TabMean5_Signi)):(360)]; 
ax.ThetaTickLabel=TabSigniWho;


%%

%%

hold on
plot(ones(1,9),TabData1_notnormal(1,:),'or')
plot(2*ones(1,16),TabData5_notnormal(1,:),'ob')





%% Isolating amplitude / speed
%close all
clearvars TabMean_Bad TabMean_Good TabMean_Bad_mdl TabMean_Good_mdl

hold on
Param_AmpSteady=41:55;
Param_AmpDyn=110:124;

Param_SpeedSteady=80:94;
Param_SpeedDyn=125:139;

Param_AmpDyn_X=110:114;
Param_AmpDyn_Y=115:119;
Param_AmpDyn_Z=120:124;

Param_SpeedDyn_X=125:129;
Param_SpeedDyn_Y=130:134;
Param_SpeedDyn_Z=135:139;

Param_AmpSteady_X=41:45;
Param_AmpSteady_Y=46:50;
Param_AmpSteady_Z=51:55;

Param_SpeedSteady_X=80:84;
Param_SpeedSteady_Y=85:89;
Param_SpeedSteady_Z=90:94;


% for i=Param_SpeedSteady_X
%     A=Data_Day1_Master(Bad,i);
%     B=Data_Day1_Master(Good,i); 
%     TabToScale2=([A;B]);
%     TabToScale=(TabToScale2-min(TabToScale2)) / (max(TabToScale2)-min(TabToScale2));   
%     TabMean_Bad=mean(TabToScale(1:length(A),1)); %Bad
%     TabMean_Good=mean(TabToScale((length(A)+1):(length(B)+length(A)),1)); %Good    
%     figure(1)
%     plot([1 2],[TabMean_Bad TabMean_Good],'-r')
%     plot([1 2],[TabMean_Bad TabMean_Good],'xk')  
%     axis([0 3 0 1])
% end
% 
% for i=Param_SpeedSteady_Y
%     A=Data_Day1_Master(Bad,i);
%     B=Data_Day1_Master(Good,i); 
%     TabToScale2=([A;B]);
%     TabToScale=(TabToScale2-min(TabToScale2)) / (max(TabToScale2)-min(TabToScale2));   
%     TabMean_Bad=mean(TabToScale(1:length(A),1)); %Bad
%     TabMean_Good=mean(TabToScale((length(A)+1):(length(B)+length(A)),1)); %Good    
%     figure(1)
%     plot([3 4],[TabMean_Bad TabMean_Good],'-b')
%     plot([3 4],[TabMean_Bad TabMean_Good],'xk')  
%     axis([0 7 0 1])
% end
% 
% for i=Param_SpeedSteady_Z
%     A=Data_Day1_Master(Bad,i);
%     B=Data_Day1_Master(Good,i); 
%     TabToScale2=([A;B]);
%     TabToScale=(TabToScale2-min(TabToScale2)) / (max(TabToScale2)-min(TabToScale2));   
%     TabMean_Bad=mean(TabToScale(1:length(A),1)); %Bad
%     TabMean_Good=mean(TabToScale((length(A)+1):(length(B)+length(A)),1)); %Good    
%     figure(1)
%     plot([5 6],[TabMean_Bad TabMean_Good],'-g')
%     plot([5 6],[TabMean_Bad TabMean_Good],'xk')  
%     axis([0 7 0 1])
% end

iter=1;
for i=Param_SpeedSteady

    A=Data_Day1_Master(Bad,i);
    B=Data_Day1_Master(Good,i);
    
    TabToScale2=([A;B]);
    TabToScale=(TabToScale2-min(TabToScale2)) / (max(TabToScale2)-min(TabToScale2));
    
    TabMean_Bad=mean(TabToScale(1:length(A),1)); %Bad
    TabMean_Good=mean(TabToScale((length(A)+1):(length(B)+length(A)),1)); %Good
    
    TabMean_Bad_mdl(iter)=mean(TabToScale(1:length(A),1));
    TabMean_Good_mdl(iter)=mean(TabToScale((length(A)+1):(length(B)+length(A)),1)); %Good
    iter=iter+1;
    
    figure(1)
    plot([1 2],[TabMean_Bad TabMean_Good],'-k')
    plot([1 2],[TabMean_Bad TabMean_Good],'.k')
    axis([0 3 0 1])
    
end

plot([1 2],[mean(TabMean_Bad_mdl) mean(TabMean_Good_mdl)],'-k','LineWidth',4)
[a]=ranksum(TabMean_Bad_mdl,TabMean_Good_mdl)




%%
% % % % % % % 
% % -- Look at individual param
% % % % % % %

a=1;
clearvars ParamDiffBad ParamDiffGood

for i=Param
    
    A=Data_Day1_Master(Bad,i);
    B=Data_Day1_Master(Good,i);
    
    TabToScale2=([A;B]); % [10*Bad 16*Good]
    %TabToScale=(TabToScale2);
    TabToScale=(TabToScale2-min(TabToScale2)) / (max(TabToScale2)-min(TabToScale2));
    
    A=TabToScale(1:numel(Bad));
    B=TabToScale(1+numel(Bad):(numel(Bad)+numel(Good)));
    
    [c b]=ttest2(A,B);

    if b<0.05
        
        ParamDiffBad(1:numel(Bad),a)=A;
        ParamDiffGood(1:numel(Good),a)=B;
        PCAaddedParamApril.Properties.VariableNames(([i]));
        a=a+1;
        
    end
    
end



%%
% % % % % % % % % % % %     
% % - - 2/ Sort Master/Yoked Switch day2 - then training
% % % % % % % % % % % % 



clearvars StandKine DataPCA SCORE COEFF DataPCA_origin

StandKine=[];

% Steady - 74 param
Param=[6 7 8:22 41:55 80:94 95:109  56:64 23:25];
% Dyn - 72 param
%Param=[110:124 125:139 155:163 179:193 194:196 197:211];

%8:22 MEAN
%23:25 ANGLE
%26:40 STD
%41:55 AMPLITUDE
%56:64 PARTICIPATION
%65:79 LINEAR EVOLUTION
%80:94 SPEED
%95:109 AREA

% DataPCA_origin=[ ...
% Data_Day1_Master(1:4,:); ...
% Data_Day2_Switch_Yoked(1:4,:); ...
% Data_Day3_Switch_Yoked(:,:); ...
% Data_Day4_Switch_Yoked(:,:); ...
% Data_Day5_Switch_Yoked(:,:); ...
% Data_Day6_Switch_Yoked(:,:); ...
% Data_Day1_Yoked(1:4,:); ...
% Data_Day2_Switch_Master(1:4,:); ...
% Data_Day3_Switch_Master(:,:); ...
% Data_Day4_Switch_Master(:,:); ...
% Data_Day5_Switch_Master(:,:); ...
% Data_Day6_Switch_Master(:,:) ...
% ];


DataPCA_origin=[...
Data_Day1_Switch_Master; ...
Data_Day2_Switch_Yoked;...
Data_Day1_Switch_Yoked;...
Data_Day2_Switch_Master;...
];



% DataPCA_origin=[ ...
% Data_Day1_Master(1:7,:); ...
% Data_Day2_Switch_Yoked(:,:); ...
% Data_Day3_Switch_Yoked(:,:); ...
% Data_Day4_Switch_Yoked(:,:); ...
% Data_Day5_Switch_Yoked(:,:); ...
% Data_Day6_Switch_Yoked(:,:); ...
% Data_Day1_Yoked(1:7,:); ...
% Data_Day2_Switch_Master(:,:); ...
% Data_Day3_Switch_Master(:,:); ...
% Data_Day4_Switch_Master(:,:); ...
% Data_Day5_Switch_Master(:,:); ...
% Data_Day6_Switch_Master(:,:) ...
% ];


clearvars StandKine

DataPCA=DataPCA_origin(:,Param);

%Export size
[n_norm p_norm] = size(DataPCA);
e = ones(n_norm,1);


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

inter=[8:12;13:17;18:22; ...
    41:45;46:50;51:55; ...
    65:69;70:74;75:79];

[a b]=size(inter);
TabOfIn=[];

for i=1:a
    TabOfIn=[];
    TabOfIn=DataPCA_origin(:,inter(i,:));
    
    [n_norm p_norm] = size(TabOfIn);
    e = ones(n_norm,1);

    StandKine=[StandKine (TabOfIn-min(min(TabOfIn)))./(max(max(TabOfIn))-min(min(TabOfIn)))];

end

inter=[80:94];

[a b]=size(inter);
TabOfIn=[];

for i=1:a
    TabOfIn=[];
    TabOfIn=DataPCA_origin(:,inter(i,:));
    
    [n_norm p_norm] = size(TabOfIn);
    e = ones(n_norm,1);

    StandKine=[StandKine (TabOfIn-min(min(TabOfIn)))./(max(max(TabOfIn))-min(min(TabOfIn)))];

end

% inter=[56:58;59:61;62:64];
% [a b]=size(inter);
% TabOfIn=[];
% 
% for i=1:a
%     TabOfIn=[];
%     TabOfIn=DataPCA_origin(:,inter(i,:));
%     
%     [n_norm p_norm] = size(TabOfIn);
%     e = ones(n_norm,1);
% 
%     StandKine=[StandKine (TabOfIn-min(min(TabOfIn)))./(max(max(TabOfIn))-min(min(TabOfIn)))];
% 
% end

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
%7 Master/Yoked for day 1 and day 2
%4 Master/Yoked for day 1 / 2 / 3 / 4 / 5 / 6

%Master and switch
[a1 b1]=size(Data_Day1_Master(1:4,:));
MasterDay1=[1:a1];
MasterDay1_Switch=[1:7];
MasterDay1_ConsDays=[1:4];

[a2 b2]=size(Data_Day2_Switch_Yoked(1:4,:));
YokedDay2=[a1+1:a1+a2];
YokedDay2_ConsDays=[a1+1:a1+4];

[a3 b3]=size(Data_Day3_Switch_Yoked);
YokedDay3=[a1+a2+1:a1+a2+a3];

[a4 b4]=size(Data_Day4_Switch_Yoked);
YokedDay4=[a1+a2+a3+1:a1+a2+a3+a4];

[a5 b5]=size(Data_Day5_Switch_Yoked);
YokedDay5=[a1+a2+a3+a4+1:a1+a2+a3+a4+a5];

[a6 b6]=size(Data_Day6_Switch_Yoked);
YokedDay6=[a1+a2+a3+a4+a5+1:a1+a2+a3+a4+a5+a6];

new_start=a1+a2+a3+a4+a5+a6;

%Yoked and switch
[a1 b1]=size(Data_Day1_Yoked(1:4,:));
YokedDay1=[new_start+1:new_start+a1];
YokedDay1_Switch=[new_start+1:new_start+7];
YokedDay1_ConsDays=[new_start+1:new_start+4];

[a2 b2]=size(Data_Day2_Switch_Master(1:4,:));
MasterDay2=[new_start+a1+1:new_start+a1+a2];
MasterDay2_ConsDays=[new_start+a1+1:new_start+a1+4];

[a3 b3]=size(Data_Day3_Switch_Master);
MasterDay3=[new_start+a1+a2+1:new_start+a1+a2+a3];

[a4 b4]=size(Data_Day4_Switch_Master);
MasterDay4=[new_start+a1+a2+a3+1:new_start+a1+a2+a3+a4];

[a5 b5]=size(Data_Day5_Switch_Master);
MasterDay5=[new_start+a1+a2+a3+a4+1:new_start+a1+a2+a3+a4+a5];

[a6 b6]=size(Data_Day6_Switch_Master);
MasterDay6=[new_start+a1+a2+a3+a4+a5+1:new_start+a1+a2+a3+a4+a5+a6];


figure
% Master Yoked Switch day 1 and 2
plot(SCORE(MasterDay1_Switch,1),SCORE(MasterDay1_Switch,2),'or')
hold on 
plot(SCORE(YokedDay1_Switch,1),SCORE(YokedDay1_Switch,2),'ob')

plot(SCORE(MasterDay2,1),SCORE(MasterDay2,2),'xb')
hold on 
plot(SCORE(YokedDay2,1),SCORE(YokedDay2,2),'xr')
axis([-2 2 -2 2])

Tab_switch_master(1:4,1)=SCORE(MasterDay1_ConsDays(:),1)
Tab_switch_master(1:4,2)=SCORE(YokedDay2_ConsDays(:),1)
Tab_switch_master(1:4,3)=SCORE(YokedDay3(:),1)
Tab_switch_master(1:4,4)=SCORE(YokedDay4(:),1)
Tab_switch_master(1:4,5)=SCORE(YokedDay5(:),1)
Tab_switch_master(1:4,6)=SCORE(YokedDay6(:),1)

Tab_switch_yoked(1:4,1)=SCORE(YokedDay1_ConsDays(:),1)
Tab_switch_yoked(1:4,2)=SCORE(MasterDay2_ConsDays(:),1)
Tab_switch_yoked(1:4,3)=SCORE(MasterDay3(:),1)
Tab_switch_yoked(1:4,4)=SCORE(MasterDay4(:),1)
Tab_switch_yoked(1:4,5)=SCORE(MasterDay5(:),1)
Tab_switch_yoked(1:4,6)=SCORE(MasterDay6(:),1)




figure
test=4;
% Master Yoked Consecutive days
plot(SCORE(MasterDay1_ConsDays(test),1),SCORE(MasterDay1_ConsDays(test),2),'or')
hold on 
plot(SCORE(YokedDay1_ConsDays(test),1),SCORE(YokedDay1_ConsDays(test),2),'ob')
% Day 2
plot(SCORE(MasterDay2_ConsDays(test),1),SCORE(MasterDay2_ConsDays(test),2),'xb')
plot(SCORE(YokedDay2_ConsDays(test),1),SCORE(YokedDay2_ConsDays(test),2),'xr')
% Day 3
plot(SCORE(MasterDay3(test),1),SCORE(MasterDay3(test),2),'>b')
plot(SCORE(YokedDay3(test),1),SCORE(YokedDay3(test),2),'>r')
% Day 4
plot(SCORE(MasterDay4(test),1),SCORE(MasterDay4(test),2),'<b')
plot(SCORE(YokedDay4(test),1),SCORE(YokedDay4(test),2),'<r')
% Day 5
plot(SCORE(MasterDay5(test),1),SCORE(MasterDay5(test),2),'.b')
plot(SCORE(YokedDay5(test),1),SCORE(YokedDay5(test),2),'.r')
% Day 6
plot(SCORE(MasterDay6(test),1),SCORE(MasterDay6(test),2),'sb')
plot(SCORE(YokedDay6(test),1),SCORE(YokedDay6(test),2),'sr')
axis([-1 1 -1 1])



%%
% % % % % % % % % % % %     
% % - - 3/ Sort Master/Yoked Switch day3 - after 48h rest
% % % % % % % % % % % % 
iterM1=1;
iterY1=1;

iterM2=1;
iterY2=1;

clearvars StandKine DataPCA SCORE COEFF DataPCA_origin

StandKine=[];

DataPCA_origin=[ ...
Data_Day1_Switch_48_Master; ...
Data_Day1_Switch_48_Yoked; ...
Data_Day2_Switch_48_Master; ...
Data_Day2_Switch_48_Yoked; ...
Data_Day1_Master; ...
Data_Day1_Yoked; ...
];


clearvars StandKine

DataPCA=DataPCA_origin(:,Param);

%Export size
[n_norm p_norm] = size(DataPCA);
e = ones(n_norm,1);


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

inter=[8:12;13:17;18:22; ...
    41:45; 46:50;51:55; ...
    80:84; 85:89;90:94; ...
    95:99;100:104;105:109];

[a b]=size(inter);
TabOfIn=[];

for i=1:a
    TabOfIn=[];
    TabOfIn=DataPCA_origin(:,inter(i,:));
    
    [n_norm p_norm] = size(TabOfIn);
    e = ones(n_norm,1);

    StandKine=[StandKine (TabOfIn-min(min(TabOfIn)))./(max(max(TabOfIn))-min(min(TabOfIn)))];

end

% inter=[56:58;59:61;62:64];
% [a b]=size(inter);
% TabOfIn=[];
% 
% for i=1:a
%     TabOfIn=[];
%     TabOfIn=DataPCA_origin(:,inter(i,:));
%     
%     [n_norm p_norm] = size(TabOfIn);
%     e = ones(n_norm,1);
% 
%     StandKine=[StandKine (TabOfIn-min(min(TabOfIn)))./(max(max(TabOfIn))-min(min(TabOfIn)))];
% 
% end

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


[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED]=pca(StandKine);%,'Centered',true);%,'VariableWeights','variance');
RecoDataCentered=SCORE*COEFF';

%Normalize COEFF
COEFF;
[n_coeff p_coeff] = size(COEFF);
e = ones(n_coeff,1);
COEFF_Norm=(2.*COEFF-e*(min(COEFF)+max(COEFF)))./(e*max(COEFF)-e*min(COEFF));
COEFF_Norm=(COEFF_Norm).*(COEFF_Norm);




figure(4)
plot(SCORE(1:4,1),SCORE(1:4,2),'or')
hold on
plot(SCORE(5:8,1),SCORE(5:8,2),'ob')

plot(SCORE(9:12,1),SCORE(9:12,2),'xr')
plot(SCORE(13:16,1),SCORE(13:16,2),'xb')
title('2d - Master & Yoked')

%axis([-2 1 -1.5 1])


%%
% % % % % % % % % % % %     
% % - - 3bis/ Sort Master/Yoked Switch day4 - after 72h rest
% % % % % % % % % % % % 

clearvars StandKine DataPCA SCORE COEFF DataPCA_origin

StandKine=[];

DataPCA_origin=[ ...
Data_Day1_Switch_72_Master; ...
Data_Day1_Switch_72_Yoked; ...
Data_Day2_Switch_72_Master; ...
Data_Day2_Switch_72_Yoked; ...
Data_Day1_Master; ...
Data_Day1_Yoked; ...
];


clearvars StandKine

DataPCA=DataPCA_origin(:,Param);

%Export size
[n_norm p_norm] = size(DataPCA);
e = ones(n_norm,1);


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

inter=[8:12;13:17;18:22; ...
    41:45; 46:50;51:55; ...
    80:84; 85:89;90:94; ...
    95:99;100:104;105:109];

[a b]=size(inter);
TabOfIn=[];

for i=1:a
    TabOfIn=[];
    TabOfIn=DataPCA_origin(:,inter(i,:));
    
    [n_norm p_norm] = size(TabOfIn);
    e = ones(n_norm,1);

    StandKine=[StandKine (TabOfIn-min(min(TabOfIn)))./(max(max(TabOfIn))-min(min(TabOfIn)))];

end

% inter=[56:58;59:61;62:64];
% [a b]=size(inter);
% TabOfIn=[];
% 
% for i=1:a
%     TabOfIn=[];
%     TabOfIn=DataPCA_origin(:,inter(i,:));
%     
%     [n_norm p_norm] = size(TabOfIn);
%     e = ones(n_norm,1);
% 
%     StandKine=[StandKine (TabOfIn-min(min(TabOfIn)))./(max(max(TabOfIn))-min(min(TabOfIn)))];
% 
% end

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


[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED]=pca(StandKine);%,'Centered',true);%,'VariableWeights','variance');
RecoDataCentered=SCORE*COEFF';

%Normalize COEFF
COEFF;
[n_coeff p_coeff] = size(COEFF);
e = ones(n_coeff,1);
COEFF_Norm=(2.*COEFF-e*(min(COEFF)+max(COEFF)))./(e*max(COEFF)-e*min(COEFF));
COEFF_Norm=(COEFF_Norm).*(COEFF_Norm);


plot(SCORE(1:5,1),SCORE(1:5,2),'or')
hold on
plot(SCORE(6:10,1),SCORE(6:10,2),'ob')

plot(SCORE(11:15,1),SCORE(11:15,2),'xr')
plot(SCORE(16:20,1),SCORE(16:20,2),'xb')

axis([-1 1 -1 1])
title('2d - 72h supp fig - Master & Yoked')


%%
% % % % % % % % % % % %     
% % - - 4/ Sort Master/Yoked Switch after 4 days training
% % % % % % % % % % % %

clearvars StandKine DataPCA SCORE COEFF DataPCA_origin

StandKine=[];

DataPCA_origin=[ ...
Data_Day1_Train4Switch_Master; ...
Data_Day1_Train4Switch_Yoked; ...
Data_Day2_Train4Switch_Master; ...
Data_Day2_Train4Switch_Yoked; ...
Data_Day3_Train4Switch_Master; ...
Data_Day3_Train4Switch_Yoked; ...
Data_Day4_Train4Switch_Master; ...
Data_Day4_Train4Switch_Yoked; ...
Data_Day5_Train4Switch_Master; ...
Data_Day5_Train4Switch_Yoked; ...
%Data_Day1_Master; ...
%Data_Day1_Yoked; ...
];

clearvars StandKine

% DataPCA=DataPCA_origin(:,Param);
% 
% %Export size
% [n_norm p_norm] = size(DataPCA);
% e = ones(n_norm,1);


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

inter=[8:12;13:17;18:22; ...
    41:45;46:50;51:55; ...
    65:69;70:74;75:79;...
    80:84;85:89;90:94];

[a b]=size(inter);
TabOfIn=[];

for i=1:a
    TabOfIn=[];
    TabOfIn=DataPCA_origin(:,inter(i,:));
    
    [n_norm p_norm] = size(TabOfIn);
    e = ones(n_norm,1);

    StandKine=[StandKine (TabOfIn-min(min(TabOfIn)))./(max(max(TabOfIn))-min(min(TabOfIn)))];

end

% inter=[56:58;59:61;62:64];
% [a b]=size(inter);
% TabOfIn=[];
% 
% for i=1:a
%     TabOfIn=[];
%     TabOfIn=DataPCA_origin(:,inter(i,:));
%     
%     [n_norm p_norm] = size(TabOfIn);
%     e = ones(n_norm,1);
% 
%     StandKine=[StandKine (TabOfIn-min(min(TabOfIn)))./(max(max(TabOfIn))-min(min(TabOfIn)))];
% 
% end

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


[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED]=pca(StandKine);%,'Centered',true);%,'VariableWeights','variance');
RecoDataCentered=SCORE*COEFF';

%Normalize COEFF
COEFF;
[n_coeff p_coeff] = size(COEFF);
e = ones(n_coeff,1);
COEFF_Norm=(2.*COEFF-e*(min(COEFF)+max(COEFF)))./(e*max(COEFF)-e*min(COEFF));
COEFF_Norm=(COEFF_Norm).*(COEFF_Norm);




% % iterM1=1;
% % iterY1=1;
% % iterM2=1;
% % iterY2=1;
% % iterM3=1;
% % iterY3=1;
% % iterM4=1;
% % iterY4=1;
% % iterM5=1;
% % iterY5=1;
% % iterM6=1;
% % iterY6=1;
% % 
% % [n_Data_All p_Data_All]=size(Data_All)
% % for i=1:n_Data_All
% %     
% %     %if Data_All(i,1)==19 || Data_All(i,1)==27 %All the experiment
% %   	if Data_All(i,1)==27 %Only experiment having a day after switch
% %         
% %         if Data_All(i,4)==1
% %             
% %             if Data_All(i,3)==1
% %                 Display_2E_1_Master(iterM1,1:3)=SCORE(i,1:3);
% %                 iterM1=iterM1+1;
% %                 i
% %             else
% %                 Display_2E_1_Yoked(iterY1,1:3)=SCORE(i,1:3);
% %                 iterY1=iterY1+1;
% %             end
% %             
% %         elseif Data_All(i,4)==2
% %             
% %             if Data_All(i,3)==1
% %                 Display_2E_2_Master(iterM2,1:3)=SCORE(i,1:3);
% %                 iterM2=iterM2+1;
% %             else
% %                 Display_2E_2_Yoked(iterY2,1:3)=SCORE(i,1:3);
% %                 iterY2=iterY2+1;
% %             end
% %             
% %         elseif Data_All(i,4)==3
% %             
% %             if Data_All(i,3)==1
% %                 Display_2E_3_Master(iterM3,1:3)=SCORE(i,1:3);
% %                 iterM3=iterM3+1;
% %             else
% %                 Display_2E_3_Yoked(iterY3,1:3)=SCORE(i,1:3);
% %                 iterY3=iterY3+1;
% %             end
% %             
% %         elseif Data_All(i,4)==4
% %             
% %             if Data_All(i,3)==1
% %                 Display_2E_4_Master(iterM4,1:3)=SCORE(i,1:3);
% %                 iterM4=iterM4+1;
% %             else
% %                 Display_2E_4_Yoked(iterY4,1:3)=SCORE(i,1:3);
% %                 iterY4=iterY4+1;
% %             end
% %             
% %         elseif Data_All(i,4)==5
% %             
% %             if Data_All(i,3)==1
% %                 Display_2E_5_Master(iterM5,1:3)=SCORE(i,1:3);
% %                 iterM5=iterM5+1;
% %             else
% %                 Display_2E_5_Yoked(iterY5,1:3)=SCORE(i,1:3);
% %                 iterY5=iterY5+1;
% %             end
% %             
% %         elseif Data_All(i,4)==6
% %             
% %             if Data_All(i,3)==1
% %                 Display_2E_6_Master(iterM6,1:3)=SCORE(i,1:3);
% %                 iterM6=iterM6+1;
% %             else
% %                 Display_2E_6_Yoked(iterY6,1:3)=SCORE(i,1:3);
% %                 iterY6=iterY6+1;
% %             end
% %             
% %         end
% %         
% %     end
% %     
% % end
% % 
% % iter=3;
% % 
% % figure(5)
% % plot(Display_2E_1_Master(iter,1),Display_2E_1_Master(iter,2),'or')
% % hold on
% % plot(Display_2E_1_Yoked(iter,1),Display_2E_1_Yoked(iter,2),'ob')
% % 
% % plot(Display_2E_2_Yoked(iter,1),Display_2E_2_Yoked(iter,2),'oc')
% % plot(Display_2E_2_Master(iter,1),Display_2E_2_Master(iter,2),'om')
% % 
% % plot(Display_2E_3_Yoked(iter,1),Display_2E_3_Yoked(iter,2),'og')
% % plot(Display_2E_3_Master(iter,1),Display_2E_3_Master(iter,2),'og')
% % plot(Display_2E_3_Master(iter,1),Display_2E_3_Master(iter,2),'xg')
% % 
% % plot(Display_2E_4_Yoked(iter,1),Display_2E_4_Yoked(iter,2),'ok')
% % plot(Display_2E_4_Master(iter,1),Display_2E_4_Master(iter,2),'ok')
% % plot(Display_2E_4_Master(iter,1),Display_2E_4_Master(iter,2),'xk')
% % % SWITCH
% % plot(Display_2E_5_Master(iter,1),Display_2E_5_Master(iter,2),'.c','MarkerSize',20)
% % plot(Display_2E_5_Yoked(iter,1),Display_2E_5_Yoked(iter,2),'.m','MarkerSize',20)
% % 
% % plot(Display_2E_6_Master(iter,1),Display_2E_6_Master(iter,2),'xc','MarkerSize',20)
% % plot(Display_2E_6_Yoked(iter,1),Display_2E_6_Yoked(iter,2),'xm','MarkerSize',20)
% % title('2E - Master & Yoked')
% % 
% % axis([-1 1 -1.5 1.5])

%%
% % % % % % % % % % % %     
% % - - 4bis/ Sort Master/Yoked Switch after 3 ONLY days training
% % % % % % % % % % % %


clearvars StandKine DataPCA SCORE COEFF DataPCA_origin

StandKine=[];

DataPCA_origin=[ ...
Data_Day1_TrainSwitch_Master; ...
Data_Day1_TrainSwitch_Yoked; ...
Data_Day2_TrainSwitch_Master; ...
Data_Day2_TrainSwitch_Yoked; ...
Data_Day3_TrainSwitch_Master; ...
Data_Day3_TrainSwitch_Yoked; ...
Data_Day4_TrainSwitch_Master; ...
Data_Day4_TrainSwitch_Yoked; ...
Data_Day1_Master; ...
Data_Day1_Yoked; ...
];

clearvars StandKine

DataPCA=DataPCA_origin(:,Param);

%Export size
[n_norm p_norm] = size(DataPCA);
e = ones(n_norm,1);


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

inter=[8:12;13:17;18:22; ...
    41:45; 46:50;51:55; ...
    80:84; 85:89;90:94; ...
    95:99;100:104;105:109];

[a b]=size(inter);
TabOfIn=[];

for i=1:a
    TabOfIn=[];
    TabOfIn=DataPCA_origin(:,inter(i,:));
    
    [n_norm p_norm] = size(TabOfIn);
    e = ones(n_norm,1);

    StandKine=[StandKine (TabOfIn-min(min(TabOfIn)))./(max(max(TabOfIn))-min(min(TabOfIn)))];

end

% inter=[56:58;59:61;62:64];
% [a b]=size(inter);
% TabOfIn=[];
% 
% for i=1:a
%     TabOfIn=[];
%     TabOfIn=DataPCA_origin(:,inter(i,:));
%     
%     [n_norm p_norm] = size(TabOfIn);
%     e = ones(n_norm,1);
% 
%     StandKine=[StandKine (TabOfIn-min(min(TabOfIn)))./(max(max(TabOfIn))-min(min(TabOfIn)))];
% 
% end

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


[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED]=pca(StandKine);%,'Centered',true);%,'VariableWeights','variance');
RecoDataCentered=SCORE*COEFF';

%Normalize COEFF
COEFF;
[n_coeff p_coeff] = size(COEFF);
e = ones(n_coeff,1);
COEFF_Norm=(2.*COEFF-e*(min(COEFF)+max(COEFF)))./(e*max(COEFF)-e*min(COEFF));
COEFF_Norm=(COEFF_Norm).*(COEFF_Norm);



plot(SCORE(1:3,1),SCORE(1:3,2),'or')
hold on
plot(SCORE(4:6,1),SCORE(4:6,2),'ob')

plot(SCORE(7:9,1),SCORE(7:9,2),'xr')
plot(SCORE(10:12,1),SCORE(10:12,2),'xb')

plot(SCORE(13:15,1),SCORE(13:15,2),'<r')
plot(SCORE(16:18,1),SCORE(16:18,2),'<b')

plot(SCORE(19:21,1),SCORE(19:21,2),'>r')
plot(SCORE(22:24,1),SCORE(22:24,2),'>b')



% % %%
% % iterM1=1;
% % iterY1=1;
% % iterM2=1;
% % iterY2=1;
% % iterM3=1;
% % iterY3=1;
% % iterM4=1;
% % iterY4=1;
% % iterM5=1;
% % iterY5=1;
% % iterM6=1;
% % iterY6=1;
% % 
% % [n_Data_All p_Data_All]=size(Data_All)
% % for i=1:n_Data_All
% %     
% %     if Data_All(i,1)==15
% %         
% %         if Data_All(i,4)==1
% %             
% %             if Data_All(i,3)==1
% %                 Display_2E_3DAYS_1_Master(iterM1,1:3)=SCORE(i,1:3);
% %                 iterM1=iterM1+1;
% %             else
% %                 Display_2E_3DAYS_1_Yoked(iterY1,1:3)=SCORE(i,1:3);
% %                 iterY1=iterY1+1;
% %             end
% %             
% %         elseif Data_All(i,4)==2
% %             
% %             if Data_All(i,3)==1
% %                 Display_2E_3DAYS_2_Master(iterM2,1:3)=SCORE(i,1:3);
% %                 iterM2=iterM2+1;
% %             else
% %                 Display_2E_3DAYS_2_Yoked(iterY2,1:3)=SCORE(i,1:3);
% %                 iterY2=iterY2+1;
% %             end
% %             
% %         elseif Data_All(i,4)==3
% %             
% %             if Data_All(i,3)==1
% %                 Display_2E_3DAYS_3_Master(iterM3,1:3)=SCORE(i,1:3);
% %                 iterM3=iterM3+1;
% %             else
% %                 Display_2E_3DAYS_3_Yoked(iterY3,1:3)=SCORE(i,1:3);
% %                 iterY3=iterY3+1;
% %             end
% %             
% %             
% %         elseif Data_All(i,4)==4
% %             
% %             if Data_All(i,3)==1
% %                 Display_2E_3DAYS_4_Master(iterM4,1:3)=SCORE(i,1:3);
% %                 iterM4=iterM4+1;
% %             else
% %                 Display_2E_3DAYS_4_Yoked(iterY4,1:3)=SCORE(i,1:3);
% %                 iterY4=iterY4+1;
% %             end
% %             
% %         elseif Data_All(i,4)==5
% %             
% %             if Data_All(i,3)==1
% %                 Display_2E_3DAYS_5_Master(iterM5,1:3)=SCORE(i,1:3);
% %                 iterM5=iterM5+1;
% %             else
% %                 Display_2E_3DAYS_5_Yoked(iterY5,1:3)=SCORE(i,1:3);
% %                 iterY5=iterY5+1;
% %             end
% %             
% %         elseif Data_All(i,4)==6
% %             
% %             if Data_All(i,3)==1
% %                 Display_2E_3DAYS_6_Master(iterM6,1:3)=SCORE(i,1:3);
% %                 iterM6=iterM6+1;
% %             else
% %                 Display_2E_3DAYS_6_Yoked(iterY6,1:3)=SCORE(i,1:3);
% %                 iterY6=iterY6+1;
% %             end
% %             
% %         end
% %         
% %     end
% %     
% % end
% % 
% % figure(5)
% % plot(Display_2E_3DAYS_1_Master(:,1),Display_2E_3DAYS_1_Master(:,2),'or')
% % hold on
% % plot(Display_2E_3DAYS_1_Yoked(:,1),Display_2E_3DAYS_1_Yoked(:,2),'ob')
% % 
% % plot(Display_2E_3DAYS_2_Yoked(:,1),Display_2E_3DAYS_2_Yoked(:,2),'oc')
% % plot(Display_2E_3DAYS_2_Master(:,1),Display_2E_3DAYS_2_Master(:,2),'om')
% % 
% % plot(Display_2E_3DAYS_3_Yoked(:,1),Display_2E_3DAYS_3_Yoked(:,2),'og')
% % plot(Display_2E_3DAYS_3_Master(:,1),Display_2E_3DAYS_3_Master(:,2),'og')
% % plot(Display_2E_3DAYS_3_Master(:,1),Display_2E_3DAYS_3_Master(:,2),'xg')
% % 
% % % plot(Display_2E_3DAYS_4_Yoked(:,1),Display_2E_3DAYS_4_Yoked(:,2),'ok')
% % % plot(Display_2E_3DAYS_4_Master(:,1),Display_2E_3DAYS_4_Master(:,2),'ok')
% % % plot(Display_2E_3DAYS_4_Master(:,1),Display_2E_3DAYS_4_Master(:,2),'xk')
% % % SWITCH
% % plot(Display_2E_3DAYS_4_Master(:,1),Display_2E_3DAYS_4_Master(:,2),'.c','MarkerSize',20)
% % plot(Display_2E_3DAYS_4_Yoked(:,1),Display_2E_3DAYS_4_Yoked(:,2),'.m','MarkerSize',20)
% % 
% % %plot(Display_2E_6_Master(:,1),Display_2E_6_Master(:,2),'xc','MarkerSize',20)
% % %plot(Display_2E_6_Yoked(:,1),Display_2E_6_Yoked(:,2),'xm','MarkerSize',20)
% % title('2E - 3DAYS - Master & Yoked')


%%
% % % % % % % % % % % %     
% % - -  6/ Ablation versus Master
% % % % % % % % % % % %

clearvars Data_Day1_Lbx1 Data_Day1_DMRT3 Data_Day1_En1 Data_Day1_Ptf1a Data_Day1_Shox2 Data_Day1_Sim1 Data_Day1_Tlx3
clearvars Data_Day2_Ptf1a Data_Day3_Ptf1a Data_Day4_Ptf1a Data_Day5_Ptf1a
clearvars Data_Day2_Tlx3 Data_Day3_Tlx3 Data_Day4_Tlx3 Data_Day5_Tlx3


% Master Data_Day1_Master

%Extract ablation
% 2 -- Lbx1
% 3 -- Ptf1a
% 4 -- Tlx3
% 5 -- DMRT3
% 6 -- En1
% 7 -- Shox2
% 8 -- Sim1
% 9 -- Gad2


Ablat=2;

iterLbx1=1;
iterPtf1a=1;
iterTlx3=1;
iterDMRT3=1;
iterEn1=1;
iterShox2=1;
iterSim1=1;
iterGad2=1;

iterLbx1_2=1;
iterLbx1_3=1;
iterLbx1_4=1;
iterLbx1_5=1;

iterPtf1a_2=1;
iterPtf1a_3=1;
iterPtf1a_4=1;
iterPtf1a_5=1;


iterTlx3_2=1;
iterTlx3_3=1;
iterTlx3_4=1;
iterTlx3_5=1;

iterGad2_2=1;
iterGad2_3=1;
iterGad2_4=1;
iterGad2_5=1;

[n p]=size(Data);

for i=1:n

    if Data(i,2)==2 && Data(i,3)==1 && Data(i,4)==1
        Data_Day1_Lbx1(iterLbx1,:)=Data(i,:);
        iterLbx1=iterLbx1+1;    
    elseif Data(i,2)==2 && Data(i,3)==1 && Data(i,4)==2
        Data_Day2_Lbx1(iterLbx1_2,:)=Data(i,:);
        iterLbx1_2=iterLbx1_2+1;
    elseif Data(i,2)==2 && Data(i,3)==1 && Data(i,4)==3
        Data_Day3_Lbx1(iterLbx1_3,:)=Data(i,:);
        iterLbx1_3=iterLbx1_3+1;
     elseif Data(i,2)==2 && Data(i,3)==1 && Data(i,4)==4
        Data_Day4_Lbx1(iterLbx1_4,:)=Data(i,:);
        iterLbx1_4=iterLbx1_4+1;   
     elseif Data(i,2)==2 && Data(i,3)==1 && Data(i,4)==5
        Data_Day5_Lbx1(iterLbx1_5,:)=Data(i,:);
        iterLbx1_5=iterLbx1_5+1;    
    end
        
            
    if Data(i,2)==3 && Data(i,3)==1 && Data(i,4)==1
        Data_Day1_Ptf1a(iterPtf1a,:)=Data(i,:);
        iterPtf1a=iterPtf1a+1; 
    elseif Data(i,2)==3 && Data(i,3)==1 && Data(i,4)==2
        Data_Day2_Ptf1a(iterPtf1a_2,:)=Data(i,:);
        iterPtf1a_2=iterPtf1a_2+1; 
    elseif Data(i,2)==3 && Data(i,3)==1 && Data(i,4)==3
        Data_Day3_Ptf1a(iterPtf1a_3,:)=Data(i,:);
        iterPtf1a_3=iterPtf1a_3+1; 
    elseif Data(i,2)==3 && Data(i,3)==1 && Data(i,4)==4
        Data_Day4_Ptf1a(iterPtf1a_4,:)=Data(i,:);
        iterPtf1a_4=iterPtf1a_4+1;     
    elseif Data(i,2)==3 && Data(i,3)==1 && Data(i,4)==5
        Data_Day5_Ptf1a(iterPtf1a_5,:)=Data(i,:);
        iterPtf1a_5=iterPtf1a_5+1;     
    end
    
    
    if Data(i,2)==4 && Data(i,3)==1 && Data(i,4)==1
        Data_Day1_Tlx3(iterTlx3,:)=Data(i,:);
        iterTlx3=iterTlx3+1;  
    elseif Data(i,2)==4 && Data(i,3)==1 && Data(i,4)==2
        Data_Day2_Tlx3(iterTlx3_2,:)=Data(i,:);
        iterTlx3_2=iterTlx3_2+1;  
    elseif Data(i,2)==4 && Data(i,3)==1 && Data(i,4)==3
        Data_Day3_Tlx3(iterTlx3_3,:)=Data(i,:);
        iterTlx3_3=iterTlx3_3+1;  
    elseif Data(i,2)==4 && Data(i,3)==1 && Data(i,4)==4
        Data_Day4_Tlx3(iterTlx3_4,:)=Data(i,:);
        iterTlx3_4=iterTlx3_4+1;  
    elseif Data(i,2)==4 && Data(i,3)==1 && Data(i,4)==5
        Data_Day5_Tlx3(iterTlx3_5,:)=Data(i,:);
        iterTlx3_5=iterTlx3_5+1;    
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
    elseif Data(i,2)==9 && Data(i,3)==1 && Data(i,4)==2
        Data_Day2_Gad2(iterGad2_2,:)=Data(i,:);
        iterGad2_2=iterGad2_2+1;  
    elseif Data(i,2)==9 && Data(i,3)==1 && Data(i,4)==3
        Data_Day3_Gad2(iterGad2_3,:)=Data(i,:);
        iterGad2_3=iterGad2_3+1;  
    elseif Data(i,2)==9 && Data(i,3)==1 && Data(i,4)==4
        Data_Day4_Gad2(iterGad2_4,:)=Data(i,:);
        iterGad2_4=iterGad2_4+1;  
    elseif Data(i,2)==9 && Data(i,3)==1 && Data(i,4)==5
        Data_Day5_Gad2(iterGad2_5,:)=Data(i,:);
        iterGad2_5=iterGad2_5+1;    
    end
    
end

%% Warning ! One tlx3 below 50%
Data_Day1_Tlx3(2,:)=[];
% Data_Day2_Tlx3(2,:)=[];
% Data_Day3_Tlx3(2,:)=[];
% Data_Day4_Tlx3(2,:)=[];
% Data_Day5_Tlx3(2,:)=[];

test=4;
Data_Day1_Ptf1a(test,:)=[];
% Data_Day2_Ptf1a(test,:)=[];
% Data_Day3_Ptf1a(test,:)=[];
% Data_Day4_Ptf1a(test,:)=[];
% Data_Day5_Ptf1a(test,:)=[];

%Remove point not good
%Data_Day1_En1(4,:)=[];
%Data_Day1_En1(3,:)=[];

%Data_Day1_Shox2(4,:)=[];

% % 
% interMaster=[1:25];
%  interYoked=[26:50];
% % %interPtf1a=[51:55];
% % %interTlx3=[56:59];
% % %interDMRT3=[60:64];
%  interSim1=[51:55];
% % %interEn1=[70:73];
% % %interShox2=[51:55];
% % 


interMaster=[1:25];
interYoked=[26:50];
interPtf1a=[51:55];
interTlx3=[56:59];
interDMRT3=[60:64];
interSim1=[65:69];
interEn1=[70:74];
interShox2=[75:79];
%interGad2=[79:79];

% % If ablation only >50% 
% Data_Day1_En1(4,:)=[];
% Data_Day1_En1(2,:)=[];
% 
% Data_Day1_Shox2(5,:)=[];
% Data_Day1_Shox2(1,:)=[];
% 
% Data_Day1_Sim1(2,:)=[];

close all

% Ablation steady


clearvars StandKine


clearvars StandKine DataPCA SCORE_Ablat COEFF_Ablat

StandKine=[];

% Steady - 74 param
Param=[6 7   8:22 41:55 80:94 95:109  23:25];
% Dyn - 72 param
%Param=[110:124 125:139 155:163 179:193 194:196 197:211];

%8:22 MEAN
%23:25 ANGLE
%26:40 STD
%41:55 AMPLITUDE
%56:64 PARTICIPATION
%65:79 LINEAR EVOLUTION
%80:94 SPEED
%95:109 AREA

%DataPCA_origin=[Data_Day1_Master;Data_Day1_Yoked;Data_Day1_Ptf1a;Data_Day1_Tlx3;Data_Day1_DMRT3;Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2; Data_Day2_Ptf1a ; Data_Day3_Ptf1a; Data_Day4_Ptf1a; Data_Day5_Ptf1a; Data_Day2_Tlx3; Data_Day3_Tlx3; Data_Day4_Tlx3; Data_Day5_Tlx3];%;Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2;Data_Day1_Yoked;Data_Day1_Master];%;Data_Day1_Yoked];Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2];
DataPCA_origin=[Data_Day1_Master;Data_Day1_Yoked;Data_Day1_Ptf1a;Data_Day1_Tlx3;Data_Day1_DMRT3;Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2];%;Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2;Data_Day1_Yoked;Data_Day1_Master];%;Data_Day1_Yoked];Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2];


%DataPCA_origin=[Data_Day1_Master;Data_Day1_Ptf1a];%;Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2;Data_Day1_Yoked;Data_Day1_Master];%;Data_Day1_Yoked];Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2];

DataPCA=DataPCA_origin(:,Param);

%Export size
[n_norm p_norm] = size(DataPCA);
e = ones(n_norm,1);


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

inter=[8:12;13:17;18:22; ...
    41:45;46:50;51:55; ...
    65:69;70:74;75:79;...
    80:84;85:89;90:94];



[a b]=size(inter);
TabOfIn=[];

for i=1:a
    TabOfIn=[];
    TabOfIn=DataPCA_origin(:,inter(i,:));
    
    [n_norm p_norm] = size(TabOfIn);
    e = ones(n_norm,1);

    StandKine=[StandKine (TabOfIn-min(min(TabOfIn)))./(max(max(TabOfIn))-min(min(TabOfIn)))];

end

% inter=[56:58;59:61;62:64];
% [a b]=size(inter);
% TabOfIn=[];
% 
% for i=1:a
%     TabOfIn=[];
%     TabOfIn=DataPCA_origin(:,inter(i,:));
%     
%     [n_norm p_norm] = size(TabOfIn);
%     e = ones(n_norm,1);
% 
%     StandKine=[StandKine (TabOfIn-min(min(TabOfIn)))./(max(max(TabOfIn))-min(min(TabOfIn)))];
% 
% end


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



clearvars COEFF_Ablat SCORE_Ablat

[COEFF_Ablat, SCORE_Ablat, LATENT_Ablat, TSQUARED_Ablat, EXPLAINED_Ablat]=pca(StandKine);%,'Centered',true);%,'VariableWeights','variance');
%RecoDataCentered=SCORE*COEFF';

%Normalize COEFF
%COEFF;
%[n_coeff p_coeff] = size(COEFF);
%e = ones(n_coeff,1);
%COEFF_Norm=(2.*COEFF-e*(min(COEFF)+max(COEFF)))./(e*max(COEFF)-e*min(COEFF));
%COEFF_Norm=(COEFF_Norm).*(COEFF_Norm);

% If ablation only >30%
interMaster=[1:25];
interYoked=[26:50];
interPtf1a=[51:55];
interTlx3=[56:59];
interDMRT3=[60:64];
interSim1=[65:69];
interEn1=[70:74];
interShox2=[75:79];

% If ablation only >50% 
% interMaster=[1:25];
% interYoked=[26:50];
% interPtf1a=[51:55];
% interTlx3=[56:59];
% interDMRT3=[60:64];
% interSim1=[65:68];
% interEn1=[69:71];
% interShox2=[72:74];


% % If ablation only >30% and no Yoked 
% clearvars DataPCA
% DataPCA=[Data_Day1_Master;Data_Day1_Ptf1a;Data_Day1_Tlx3;Data_Day1_DMRT3;Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2]%;Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2;Data_Day1_Yoked;Data_Day1_Master];%;Data_Day1_Yoked];Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2];
%interMaster=[1:25];
%interPtf1a=[26:30];
% interTlx3=[31:34];
% interDMRT3=[35:39];
% interSim1=[40:44];
% interEn1=[45:48];
% interShox2=[49:52];


% % If ablation only >50% and no Yoked 
% clearvars DataPCA
% DataPCA=[Data_Day1_Master;Data_Day1_Ptf1a;Data_Day1_Tlx3;Data_Day1_DMRT3;Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2]%;Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2;Data_Day1_Yoked;Data_Day1_Master];%;Data_Day1_Yoked];Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2];
% interMaster=[1:25];
% interPtf1a=[26:30];
% interTlx3=[31:34];
% interDMRT3=[35:39];
% interSim1=[40:43];
% interEn1=[44:46];
% interShox2=[47:49];

% Master Yoked + One pop
%clearvars DataPCA
% %DataPCA=[Data_Day1_Master;Data_Day1_Ptf1a;Data_Day1_Tlx3;Data_Day1_DMRT3;Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2]%;Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2;Data_Day1_Yoked;Data_Day1_Master];%;Data_Day1_Yoked];Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2];
% DataPCA=[Data_Day1_Master;...
%     Data_Day1_Yoked;...
%     %Data_Day1_Ptf1a;...
%     %Data_Day1_Tlx3;...
%     %Data_Day1_DMRT3;...
%     %Data_Day1_Sim1;...
%     %Data_Day1_En1;...
%     Data_Day1_Shox2...
%     ]
% 
% interMaster=[1:25];
% interYoked=[26:50];
% interPtf1a=[51:55];
% interTlx3=[51:54];
% interDMRT3=[51:55];
% interSim1=[51:55];
% interEn1=[51:55];
% interShox2=[51:55];

% interTlx3=[56:59];
% interDMRT3=[60:64];
% interSim1=[65:69];
% interEn1=[70:74];
% interShox2=[75:79];
% %%
% 


%Master
plot(SCORE_Ablat(interMaster,1),SCORE_Ablat(interMaster,2),'.k','MarkerSize',20)
hold on
%plotell([SCORE_Ablat(1:24,1) SCORE_Ablat(1:24,2)],'-r','ok');
center = [mean(SCORE_Ablat(interMaster,1)), mean(SCORE_Ablat(interMaster,2))];
center_Master=center;
stdev = [std(SCORE_Ablat(interMaster,1)), std(SCORE_Ablat(interMaster,2))];
stdev_Master=stdev;
axh = gca(); 
llc = center(:)-stdev(:);
wh = stdev(:)*2; 
%rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1]); 
plot(center(1),center(2),'xk')

% Yoked
plot(SCORE_Ablat(interYoked,1),SCORE_Ablat(interYoked,2),'.k','MarkerSize',10)
%plotell([SCORE_Ablat(interYoked,1) SCORE_Ablat(interYoked,2)],'--k','ok');
center = [mean(SCORE_Ablat(interYoked,1)), mean(SCORE_Ablat(interYoked,2))];
center_Yoked=center;
stdev = [std(SCORE_Ablat(interYoked,1)), std(SCORE_Ablat(interYoked,2))];
stdev_Yoked=stdev;
axh = gca(); 
llc = center(:)-stdev(:);
wh = stdev(:)*2; 
%rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1]); 
plot(center(1),center(2),'xk')

% %%
% DataPCA=[Data_Day1_Master;Data_Day1_Yoked;Data_Day1_Ptf1a;Data_Day1_Tlx3;Data_Day1_DMRT3;Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2]%;Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2;Data_Day1_Yoked;Data_Day1_Master];%;Data_Day1_Yoked];Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2];
% 
% 
% yfitMaster = trainedModel_Para_SVM.predictFcn(Data_Day1_Master(:,Param)) 
% yfitYoked = trainedModel_Para_SVM.predictFcn(Data_Day1_Yoked(:,Param)) 
% yfitPtf1a = trainedModel_Para_SVM.predictFcn(Data_Day1_Ptf1a(:,Param)) 
% yfitTlx3 = trainedModel_Para_SVM.predictFcn(Data_Day1_Tlx3(:,Param)) 
% yfitDMRT3 = trainedModel_Para_SVM.predictFcn(Data_Day1_DMRT3(:,Param)) 
% yfitSim1 = trainedModel_Para_SVM.predictFcn(Data_Day1_Sim1(:,Param)) 
% yfitEn1 = trainedModel_Para_SVM.predictFcn(Data_Day1_En1(:,Param)) 
% yfitShox2 = trainedModel_Para_SVM.predictFcn(Data_Day1_Shox2(:,Param)) 
% 
% plot(1,yfitMaster,'xk')
% hold on
% plot(2,yfitYoked,'xk')
% plot(3,yfitPtf1a,'xk')
% plot(4,yfitTlx3,'xk')
% plot(5,yfitDMRT3,'xk')
% plot(6,yfitSim1,'xk')
% plot(7,yfitEn1,'xk')
% plot(8,yfitShox2,'xk')

%

% yfitMaster = trainedModel.predictFcn(Data_Day1_Master(:,Param)) 
% yfitYoked = trainedModel.predictFcn(Data_Day1_Yoked(:,Param)) 
% yfitPtf1a = trainedModel.predictFcn(Data_Day1_Ptf1a(:,Param)) 
% yfitTlx3 = trainedModel.predictFcn(Data_Day1_Tlx3(:,Param)) 
% yfitDMRT3 = trainedModel.predictFcn(Data_Day1_DMRT3(:,Param)) 
% yfitSim1 = trainedModel.predictFcn(Data_Day1_Sim1(:,Param)) 
% yfitEn1 = trainedModel.predictFcn(Data_Day1_En1(:,Param)) 
% yfitShox2 = trainedModel.predictFcn(Data_Day1_Shox2(:,Param)) 
% 
% plot(1,yfitMaster,'xk')
% hold on
% plot(2,yfitYoked,'xk')
% plot(3,yfitPtf1a,'xk')
% plot(4,yfitTlx3,'xk')
% plot(5,yfitDMRT3,'xk')
% plot(6,yfitSim1,'xk')
% plot(7,yfitEn1,'xk')
% plot(8,yfitShox2,'xk')
% 


% %Ptf1a
plot(SCORE_Ablat(interPtf1a,1),SCORE_Ablat(interPtf1a,2),'.b','MarkerSize',20)
%plotell([SCORE_Ablat(49:53,1) SCORE_Ablat(49:53,2)],'-b','ob');
center = [mean(SCORE_Ablat(interPtf1a,1)), mean(SCORE_Ablat(interPtf1a,2))];
center_Ptf1a=center;
stdev = [std(SCORE_Ablat(interPtf1a,1)), std(SCORE_Ablat(interPtf1a,2))];
stdev_Ptf1a=stdev;
axh = gca(); 
llc = center(:)-stdev(:);
wh = stdev(:)*2; 
%rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1],'EdgeColor','b'); 
plot(center(1),center(2),'xb')
% 
% 
% %Tlx3
plot(SCORE_Ablat(interTlx3,1),SCORE_Ablat(interTlx3,2),'.r','MarkerSize',20)
%plotell([SCORE_Ablat(54:77,1) SCORE_Ablat(54:77,2)],'-r','or');
center = [mean(SCORE_Ablat(interTlx3,1)), mean(SCORE_Ablat(interTlx3,2))];
center_Tlx3=center;
stdev = [std(SCORE_Ablat(interTlx3,1)), std(SCORE_Ablat(interTlx3,2))];
stdev_Tlx3=stdev;
axh = gca(); 
llc = center(:)-stdev(:);
wh = stdev(:)*2; 
%rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1],'EdgeColor','r'); 
plot(center(1),center(2),'xr')
% 


% %DMRT3
plot(SCORE_Ablat(interDMRT3,1),SCORE_Ablat(interDMRT3,2),'.g','MarkerSize',20)
%plotell([SCORE_Ablat(58:62,1) SCORE_Ablat(58:62,2)],'-g','og');
center = [mean(SCORE_Ablat(interDMRT3,1)), mean(SCORE_Ablat(interDMRT3,2))];
center_DMRT3=center;
stdev = [std(SCORE_Ablat(interDMRT3,1)), std(SCORE_Ablat(interDMRT3,2))];
stdev_DMRT3=stdev;
axh = gca(); 
llc = center(:)-stdev(:);
wh = stdev(:)*2; 
%rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1],'EdgeColor','g'); 
plot(center(1),center(2),'xg')
% 


%SIM1
plot(SCORE_Ablat(interSim1,1),SCORE_Ablat(interSim1,2),'.y','MarkerSize',20)
%plotell([SCORE_Ablat(63:67,1) SCORE_Ablat(63:67,2)],'-y','oy');
center = [mean(SCORE_Ablat(interSim1,1)), mean(SCORE_Ablat(interSim1,2))];
center_Sim1=center;
stdev = [std(SCORE_Ablat(interSim1,1)), std(SCORE_Ablat(interSim1,2))];
stdev_Sim1=stdev;
axh = gca(); 
llc = center(:)-stdev(:);
wh = stdev(:)*2; 
%rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1],'EdgeColor','y'); 
plot(center(1),center(2),'xy')


% %EN1
plot(SCORE_Ablat(interEn1,1),SCORE_Ablat(interEn1,2),'.c','MarkerSize',20)
%plotell([SCORE_Ablat(68:72,1) SCORE_Ablat(68:72,2)],'-c','oc');
center = [mean(SCORE_Ablat(interEn1,1)), mean(SCORE_Ablat(interEn1,2))];
center_En1=center;
stdev = [std(SCORE_Ablat(interEn1,1)), std(SCORE_Ablat(interEn1,2))];
stdev_En1=stdev;
axh = gca(); 
llc = center(:)-stdev(:);
wh = stdev(:)*2; 
%rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1],'EdgeColor','c'); 
plot(center(1),center(2),'xc')
% 
% 
% %SHOX2
plot(SCORE_Ablat(interShox2,1),SCORE_Ablat(interShox2,2),'.m','MarkerSize',20)
%plotell([SCORE_Ablat(73:77,1) SCORE_Ablat(73:77,2)],'-m','om');
center = [mean(SCORE_Ablat(interShox2,1)), mean(SCORE_Ablat(interShox2,2))];
center_Shox2=center;
stdev = [std(SCORE_Ablat(interShox2,1)), std(SCORE_Ablat(interShox2,2))];
stdev_Shox2=stdev;
axh = gca(); 
llc = center(:)-stdev(:);
wh = stdev(:)*2; 
%rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1],'EdgeColor','m'); 
plot(center(1),center(2),'xm')

%axis([-1.5 1 -1 1.5])



%%
% 
% clearvars p_valuemaster p_valueyoked
% 
% inter_m=1;
% inter_y=1;
% 
% for i=1:63
%     
%     [a b]=ttest2(StandKine(interMaster,i),StandKine(interPtf1a,i));
%     [c d]=ttest2(StandKine(interYoked,i),StandKine(interPtf1a,i));
% 
%     if b<0.05
%         
%         p_valuemaster(i)=b;
%         inter_m=inter_m+1;
%         
%     else 
%         
%         p_valuemaster(i)=1;
%         
%     end
%     
%     if d<0.05
%         
%         p_valueyoked(i)=d;
%         inter_y=inter_y+1;
%         
%     else 
%         p_valueyoked(i)=1;
%     end
% end
% 
% 
% %%
% 
% clearvars p_valuemaster p_valueyoked
% close all
% inter_m=1;
% inter_y=1;
% hold on
% plo=1;
% 
% for i=1:63
%         
% 	[a b]=ttest2(StandKine(interMaster,i),StandKine(interPtf1a,i));
%     
%     if b<0.05
%     plot(plo,mean(StandKine(interMaster,i)),'or')
%     plot(plo,mean(StandKine(interYoked,i)),'ob')
%     
%     plot(plo,mean(StandKine(interPtf1a,i)),'og')
%     
%     p_valuemaster(i)=b;
%     inter_m=inter_m+1;
%     
%     else 
%         
%         p_valuemaster(i)=1;
%         
%     end
%     
%     plo=plo+1;
%     
% end
% 
% 
% % Higher general amplitude
% % Tlx3 High Speed Y - Low Z - Low Y area
% % DMRT3 High Y Area - High Y Mean
% % En1 / Ptf1a look the same
% % Shox2 lower amplitude - Lower mean X and Y
% 
% %Ptf1a
% 
% axis([0 65 0 1])
% grid on
% grid minor
% 
% 
% 
% %% 
% 
% 
% clearvars p_valuemaster p_valueyoked
% close all
% inter_m=1;
% inter_y=1;
% hold on
% plo=1;
% 
% toshuffle=[interMaster,interYoked];
% 
% interMaster_shuffle=interMaster(randperm(length(interMaster)));
% interMaster_shuffle_1=interMaster_shuffle(1:floor(length(interMaster)/2))
% interMaster_shuffle_2=interMaster_shuffle(floor(length(interMaster)/2):length(interMaster))
% 
% inter_shuffle=toshuffle(randperm(length(toshuffle)));
% 
% inter_shuffle_1=inter_shuffle(1:floor(length(toshuffle)/2));
% inter_shuffle_2=inter_shuffle(floor(length(toshuffle)/2):length(toshuffle));
% 
% 
% for k=1:1000
% 
% inter_shuffle=toshuffle(randperm(length(toshuffle)));
% 
% inter_shuffle_1=inter_shuffle(1:floor(length(toshuffle)/2));
% inter_shuffle_2=inter_shuffle(floor(length(toshuffle)/2):length(toshuffle));
% 
% 
% for i=1:63
%         
% 	[a b]=ttest2(StandKine(inter_shuffle_1,i),StandKine(inter_shuffle_2,i));
%     
%     if b<0.05
%     %plot(plo,mean(StandKine(toshuffle,i)),'or')
%     %plot(plo,mean(StandKine(interYoked,i)),'ob')
%     
%     %plot(plo,mean(StandKine(interPtf1a,i)),'og')
%     
%     p_valuemaster(i)=b;
%     inter_m=inter_m+1;
%     
%     else 
%         
%         p_valuemaster(i,k)=1;
%         
%     end
%     
%     plo=plo+1;
%     
% end
% 
% end
% 
% p_valuemaster_alldata=sum(p_valuemaster(:,:),2)/1000;
% 
% 
% %%
% 
% 
% toshuffle=[interMaster,interYoked];
% 
% 
% for k=1:100
% 
% 
% inter_shuffle=toshuffle(randperm(length(toshuffle)));
% 
% inter_shuffle_1=inter_shuffle;%inter_shuffle(1:floor(length(toshuffle)/2));
% inter_shuffle_2=interPtf1a;%inter_shuffle(floor(length(toshuffle)/2):length(toshuffle));
% 
% 
% 
% p_valuemaster_group=zeros(1,13);
% 
% for i=1:13
%     
%     if i==13 %Groupe Angle
%         
%         for j=1:3
%             
%             val=5*(i-1)+j;
%             [a b]=ttest2(StandKine(inter_shuffle_1,val),StandKine(inter_shuffle_2,val));
%             if b<0.05
%                 p_valuemaster_group(i)=1+p_valuemaster_group(i);
%             else
%                 p_valuemaster_group(i)=0+p_valuemaster_group(i);
%             end
%             
%         end
%         
%     else
%         
%         for j=1:5
%             
%             val=5*(i-1)+j;
%             
%             [a, b]=ttest2(StandKine(inter_shuffle_1,val),StandKine(inter_shuffle_2,val));
%             
%             if b<0.05
%                 
%                 p_valuemaster_group(i)=1+p_valuemaster_group(i);
%                 
%             else
%                 
%                 p_valuemaster_group(i)=0+p_valuemaster_group(i);
%                 
%             end
%             
%         end
%     end
%     
% end
% 
% p_valuemaster_group(1:12)=p_valuemaster_group(1:12)*100/5;
% p_valuemaster_group(13)=p_valuemaster_group(13)*100/3;
% 
% p_valuemaster_group_shuffle(:,k)=p_valuemaster_group(:);
% 
% end 
% 
% 
% p_valuemaster_group_perc(:)=sum(p_valuemaster_group_shuffle(:,:),2)/100


%%


% %%
% 
% clearvars p_valuemaster p_valueyoked
% close all
% inter_m=1;
% inter_y=1;
% hold on
% plo=1;
% 
% toshuffle=[interMaster];
% 
% 
% for ts=1:1000
%     
%     interMaster_shuffle=toshuffle(randperm(length(toshuffle)));
%     
%     interMaster_shuffle_1=interMaster_shuffle(1:5);%interMaster_shuffle(1:floor(length(interMaster)/2));
%     
%     interMaster_shuffle_2=interPtf1a;%interMaster_shuffle(floor(length(interMaster)/2):length(interMaster));
%     
%     
%     for i=1:63
%         
%         [a b]=ttest2(StandKine(interMaster_shuffle_1,i),StandKine(interMaster_shuffle_2,i));
%         
%         if b<0.05
%             %plot(plo,mean(StandKine(interMaster,i)),'or')
%             %plot(plo,mean(StandKine(interYoked,i)),'ob')
%             
%             %plot(plo,mean(StandKine(interShox2,i)),'og')
%             
%             p_valuemaster(i,ts)=0;
%             inter_m=inter_m+1;
%             
%         else
%             
%             p_valuemaster(i,ts)=1;
%             
%         end
%         
%         plo=plo+1;
%         
%     end
%     
%     
%     
%     
% end
% 
% for i=1:63
% 
%     Sh(i)=sum(p_valuemaster(i,:))/1000;
% end


%% 
%% --ABLATION DYNAMIC
%
% % % % % % % %

clearvars StandKine DataPCA_origin SCORE_Ablat COEFF_Ablat
StandKine=[];

DataPCA_origin=[Data_Day1_Master;Data_Day1_Yoked;Data_Day1_Ptf1a;Data_Day1_Tlx3;Data_Day1_DMRT3;Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2]%;Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2;Data_Day1_Yoked;Data_Day1_Master];%;Data_Day1_Yoked];Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2];
%DataPCA_origin=[Data_Day1_Master;Data_Day1_Yoked;Data_Day1_Ptf1a;Data_Day1_DMRT3;Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2]%;Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2;Data_Day1_Yoked;Data_Day1_Master];%;Data_Day1_Yoked];Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2];


interMaster=[1:25];
interYoked=[26:50];
interPtf1a=[51:55];
interTlx3=[56:59];
interDMRT3=[60:64];
interSim1=[65:69];
interEn1=[70:74];
interShox2=[75:79];

% 
% DataPCA_origin=[ ...
%     Data_Day1_Master;...
%     Data_Day1_Yoked;...
%     %Data_Day1_Ptf1a;...
%      %Data_Day1_Tlx3;...
%      Data_Day1_DMRT3;...
%      %Data_Day1_Sim1;...
%    %Data_Day1_En1;...
%     %Data_Day1_Shox2...
%     ];
% 
%  Col=2;
%  
% interMaster=[1:25];
% interYoked=[26:50];
% interPtf1a=[51:55];
%  interTlx3=[51:54];
%  interDMRT3=[51:55];
%  interSim1=[51:55];
%   interEn1=[51:55];
%  interShox2=[51:55];

%interMaster=[1:25];
% interPtf1a=[26:30]-25;
% interTlx3=[31:34]-25;
% interDMRT3=[35:39]-25;
% interSim1=[40:44]-25;
% interEn1=[45:48]-25;
% interShox2=[49:52]-25;


% Steady - 74 param
%Param=[6 7 8:22 41:55 80:94 95:109  56:64 23:25];
% Dyn - 72 param
%Param=[110:124 179:193 194:196]; %179:193 194:196

%110:124 AMPLITUDE SPIKE
%125:139 SPEED
%140:154 ACC
%155:163 PARTICIPATION
%164:178 MAX MIN ENDPOINT
%179:193 AREA
%194:196 ANGLE
%197:211 MEAN 

% %%
% 
% for i=130
%     
%     [a b]=ttest2(Data_Day1_Master(:,i),Data_Day1_Ptf1a(:,i))
% 
% end
% 


% %

inter=[110:114;115:119;120:124; ...
    125:129;130:134;135:139; ...
    179:183;184:188;189:193;...
    197:201;202:206;207:211; ...
    ];

[a b]=size(inter);
TabOfIn=[];

for i=1:a
    TabOfIn=[];
    TabOfIn=DataPCA_origin(:,inter(i,:));
    
    [n_norm p_norm] = size(TabOfIn);
    e = ones(n_norm,1);

    StandKine=[StandKine (TabOfIn-min(min(TabOfIn)))./(max(max(TabOfIn))-min(min(TabOfIn)))];

end
% 
% inter=[155:157;158:160;161:163];
% [a b]=size(inter);
% TabOfIn=[];
% 
% for i=1:a
%     TabOfIn=[];
%     TabOfIn=DataPCA_origin(:,inter(i,:));
%     
%     [n_norm p_norm] = size(TabOfIn);
%     e = ones(n_norm,1);
% 
%     StandKine=[StandKine (TabOfIn-min(min(TabOfIn)))./(max(max(TabOfIn))-min(min(TabOfIn)))];
% 
% end
% % 
% 
inter=[194:196];
[a b]=size(inter);
TabOfIn=[];
for i=1:a
    TabOfIn=[];
    TabOfIn=DataPCA_origin(:,inter(i,:));
    
    [n_norm p_norm] = size(TabOfIn);
    e = ones(n_norm,1);

    StandKine=[StandKine (TabOfIn-min(min(TabOfIn)))./(max(max(TabOfIn))-min(min(TabOfIn)))];

end

% inter=[164:168;169:173;174:178];
% [a b]=size(inter);
% TabOfIn=[];
% 
% for i=1:a
%     TabOfIn=[];
%     TabOfIn=DataPCA_origin(:,inter(i,:));
%     
%     [n_norm p_norm] = size(TabOfIn);
%     e = ones(n_norm,1);
% 
%     StandKine=[StandKine (TabOfIn-min(min(TabOfIn)))./(max(max(TabOfIn))-min(min(TabOfIn)))];
% 
% end

clearvars COEFF_Ablat SCORE_Ablat CoeffVisu

[COEFF_Ablat, SCORE_Ablat, LATENT_Ablat, TSQUARED_Ablat, EXPLAINED_Ablat]=pca(StandKine);%,'Centered',true);%,'VariableWeights','variance');



%CoefToStudy=COEFF_Ablat;

Joint_3=[];
Joint_5=[];

% Mean / Amplitude / Speed / Area / Partici / Angle 

% % Dynamic
% Joint_5=[ ...
%     46:50; 51:55; 56:60;
%     1:5; 6:10; 11:15; ...
%     16:20; 21:25; 26:30; ...
%     31:35; 36:40; 41:45];
% 
% Joint_3=[61:63; 64:66; 67:69; 70:72];
% 
% Joint_2=0;

% Steady
% inter=[8:12;13:17;18:22;41:45;46:50; ...
%     51:55;80:84;85:89;90:94; ...
%     95:99;100:104;105:109];

% Joint_5=[ ...
%     1:5; 6:10; 11:15; ...
%     16:20; 21:25; 26:30; ...
%     31:35; 36:40; 41:45; ...
%     46:50; 51:55; 56:60];
% 
% Joint_3=[61:63; 64:66; 67:69; 70:72];
% 
% 
% %Joint5
% [a b]=size(Joint_5);
% CoeffVisu=[];
% CoeffVisu_iter=[];
% 
% for i=1:a
% 
%     CoeffVisu_iter=sum(abs(COEFF_Ablat(Joint_5(i,:),Col)));
%     
%     CoeffVisu=[CoeffVisu CoeffVisu_iter];
%     
% end
%     
% %Joint3
% [a b]=size(Joint_3);
% for i=1:a
% 
%     CoeffVisu_iter=sum(abs(COEFF_Ablat(Joint_3(i,:),Col)));
%     
%     CoeffVisu=[CoeffVisu CoeffVisu_iter];
%     
% end
% 





% % % %  If ablation only >30%
% interMaster=[1:25];
% interYoked=[26:50];
% interPtf1a=[51:55];
% interTlx3=[56:59];
% interDMRT3=[60:64];
% interSim1=[65:69];
% interEn1=[70:73];
% interShox2=[74:78];

% If ablation only >50% 
% interMaster=[1:25];
% interYoked=[26:50];
% interPtf1a=[51:55];
% interTlx3=[56:59];
% interDMRT3=[60:64];
% interSim1=[65:68];
% interEn1=[69:71];
% interShox2=[72:74];


% % If ablation only >30% and no Yoked 
% clearvars DataPCA
% DataPCA=[Data_Day1_Master;Data_Day1_Ptf1a;Data_Day1_Tlx3;Data_Day1_DMRT3;Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2]%;Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2;Data_Day1_Yoked;Data_Day1_Master];%;Data_Day1_Yoked];Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2];
% interMaster=[1:25];
% interPtf1a=[26:30];
% interTlx3=[31:34];
% interDMRT3=[35:39];
% interSim1=[40:44];
% interEn1=[45:49];
% interShox2=[50:54];


% % If ablation only >50% and no Yoked 
% clearvars DataPCA
% DataPCA=[Data_Day1_Master;Data_Day1_Ptf1a;Data_Day1_Tlx3;Data_Day1_DMRT3;Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2]%;Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2;Data_Day1_Yoked;Data_Day1_Master];%;Data_Day1_Yoked];Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2];
% interMaster=[1:25];
% interPtf1a=[26:30];
% interTlx3=[31:34];
% interDMRT3=[35:39];
% interSim1=[40:43];
% interEn1=[44:46];
% interShox2=[47:49];

% Master Yoked + One pop
% clearvars DataPCA
% DataPCA=[Data_Day1_Master;Data_Day1_Ptf1a;Data_Day1_Tlx3;Data_Day1_DMRT3;Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2]%;Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2;Data_Day1_Yoked;Data_Day1_Master];%;Data_Day1_Yoked];Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2];
% DataPCA=[Data_Day1_Master;...
%     Data_Day1_Yoked;...
%     Data_Day1_Ptf1a;...
%     Data_Day1_Tlx3;...
%     %Data_Day1_DMRT3;...
%     %Data_Day1_Sim1;...
%     %Data_Day1_En1;...
%     %Data_Day1_Shox2...
%     ];
% 
%  interMaster=[1:5];
%  interYoked=[6:10];
%  interPtf1a=[11:15];

 %interTlx3=[31:34];
 
% interMaster=[1:25];
% interYoked=[26:50];
% % interPtf1a=[51:55];
% % % interPtf1a=[51:55];
%  % interTlx3=[51:54];
%  % interDMRT3=[51:55];
% interSim1=[51:55];
% % % interEn1=[51:55];
% % % interShox2=[51:55];
% 
% errorbar(-2,mean(SCORE_Ablat(interMaster,1)),std(SCORE_Ablat(interMaster,1))/sqrt(numel(interMaster)))
% hold on
% errorbar(-2,mean(SCORE_Ablat(interYoked,1)),std(SCORE_Ablat(interYoked,1))/sqrt(numel(interMaster)))
% errorbar(-1,mean(SCORE_Ablat(interPtf1a,1)),std(SCORE_Ablat(interPtf1a,1))/sqrt(numel(interPtf1a)))
% 
% 
% errorbar(1,mean(SCORE_Ablat(interMaster,2)),std(SCORE_Ablat(interMaster,2))/sqrt(numel(interMaster)))
% hold on
% errorbar(1,mean(SCORE_Ablat(interYoked,2)),std(SCORE_Ablat(interYoked,2))/sqrt(numel(interMaster)))
% errorbar(2,mean(SCORE_Ablat(interPtf1a,2)),std(SCORE_Ablat(interPtf1a,2))/sqrt(numel(interPtf1a)))
% 
% 
% axis([-4 4 -1 1])

% %

% Param=[110:124 125:139 155:163 179:193 194:196 197:211];
% c=1;
% clearvars AA
% 
%  interMaster=[14:19];
%  interYoked=[26:50];
%  interPtf1a=[51:55];
%  
% for i=1:72
%     
%     [a b]=ttest2(StandKine(interMaster,i),StandKine(interPtf1a,i));
%     
%     if b<0.05
%         
%      AA([c])=PCAaddedParamApril.Properties.VariableNames(([Param(i)]))
%      c=c+1;
%      
%     end
%     
% end

% %

% interTlx3=[56:59];
% interDMRT3=[60:64];
% interSim1=[65:69];
% interEn1=[70:74];
% interShox2=[75:79];

% Only abalted popu
% DataPCA=[ ...
%     %Data_Day1_Master;...
%     %Data_Day1_Yoked;...
%     Data_Day1_Ptf1a;...
%     Data_Day1_Tlx3;...
%     Data_Day1_DMRT3;...
%     Data_Day1_Sim1;...
%     Data_Day1_En1;...
%     Data_Day1_Shox2...
%     ];

% interPtf1a=[1:5];
% interTlx3=[6:9];
% interDMRT3=[10:14];
% interSim1=[15:19];
% interEn1=[20:24];
% interShox2=[25:29];

% 
% interMaster=[1:25];
% interYoked=[26:50];
% interDMRT3=[51:55];
% interSim1=[56:60];
% interEn1=[61:65];
% interShox2=[66:69];

%
% % Master
plot(SCORE_Ablat(interMaster,1),SCORE_Ablat(interMaster,2),'.k','MarkerSize',20)
hold on
%plotell([SCORE_Ablat(1:24,1) SCORE_Ablat(1:24,2)],'-r','ok');
center = [mean(SCORE_Ablat(interMaster,1)), mean(SCORE_Ablat(interMaster,2))];
center_Master=center;
stdev = [std(SCORE_Ablat(interMaster,1)), std(SCORE_Ablat(interMaster,2))];
stdev_Master=stdev;
axh = gca(); 
llc = center(:)-stdev(:);
wh = stdev(:)*2; 
rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1]); 
plot(center(1),center(2),'xk')

% %Yoked
plot(SCORE_Ablat(interYoked,1),SCORE_Ablat(interYoked,2),'ok')
%plotell([SCORE_Ablat(interYoked,1) SCORE_Ablat(interYoked,2)],'--k','ok');
center = [mean(SCORE_Ablat(interYoked,1)), mean(SCORE_Ablat(interYoked,2))];
center_Yoked=center;
stdev = [std(SCORE_Ablat(interYoked,1)), std(SCORE_Ablat(interYoked,2))];
stdev_Yoked=stdev;
axh = gca(); 
llc = center(:)-stdev(:);
wh = stdev(:)*2; 
rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1]); 
plot(center(1),center(2),'xk')

hold on
% %Ptf1a
plot(SCORE_Ablat(interPtf1a,1),SCORE_Ablat(interPtf1a,2),'.b','MarkerSize',20)
%plotell([SCORE_Ablat(49:53,1) SCORE_Ablat(49:53,2)],'-b','ob');
center = [mean(SCORE_Ablat(interPtf1a,1)), mean(SCORE_Ablat(interPtf1a,2))];
center_Ptf1a=center;
stdev = [std(SCORE_Ablat(interPtf1a,1)), std(SCORE_Ablat(interPtf1a,2))];
stdev_Ptf1a=stdev;
axh = gca(); 
llc = center(:)-stdev(:);
wh = stdev(:)*2; 
rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1],'EdgeColor','b'); 
plot(center(1),center(2),'xb')
% 
% 
% Tlx3
plot(SCORE_Ablat(interTlx3,1),SCORE_Ablat(interTlx3,2),'.r','MarkerSize',20)
%%plotell([SCORE_Ablat(interTlx3,1) SCORE_Ablat(interTlx3,2)],'-r','or');
center = [mean(SCORE_Ablat(interTlx3,1)), mean(SCORE_Ablat(interTlx3,2))];
center_Tlx3=center;
stdev = [std(SCORE_Ablat(interTlx3,1)), std(SCORE_Ablat(interTlx3,2))];
stdev_Tlx3=stdev;
axh = gca(); 
llc = center(:)-stdev(:);
wh = stdev(:)*2; 
rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1],'EdgeColor','r'); 
plot(center(1),center(2),'xr')

% 
% 
% % DMRT3
plot(SCORE_Ablat(interDMRT3,1),SCORE_Ablat(interDMRT3,2),'.g','MarkerSize',20)
%%plotell([SCORE_Ablat(interDMRT3,1) SCORE_Ablat(interDMRT3,2)],'-g','og');
center = [mean(SCORE_Ablat(interDMRT3,1)), mean(SCORE_Ablat(interDMRT3,2))];
center_DMRT3=center;
stdev = [std(SCORE_Ablat(interDMRT3,1)), std(SCORE_Ablat(interDMRT3,2))];
stdev_DMRT3=stdev;
axh = gca(); 
llc = center(:)-stdev(:);
wh = stdev(:)*2; 
rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1],'EdgeColor','g'); 
plot(center(1),center(2),'xg')



% % % % SIM1
plot(SCORE_Ablat(interSim1,1),SCORE_Ablat(interSim1,2),'.y','MarkerSize',20)
plotell([SCORE_Ablat(interSim1,1) SCORE_Ablat(interSim1,2)],'-y','oy');
center = [mean(SCORE_Ablat(interSim1,1)), mean(SCORE_Ablat(interSim1,2))];
center_Sim1=center;
stdev = [std(SCORE_Ablat(interSim1,1)), std(SCORE_Ablat(interSim1,2))];
stdev_Sim1=stdev;
axh = gca(); 
llc = center(:)-stdev(:);
wh = stdev(:)*2; 
rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1],'EdgeColor','y'); 
plot(center(1),center(2),'xy')


% % % % % EN1
plot(SCORE_Ablat(interEn1,1),SCORE_Ablat(interEn1,2),'.c','MarkerSize',20)
plotell([SCORE_Ablat(interEn1,1) SCORE_Ablat(interEn1,2)],'-c','oc');
center = [mean(SCORE_Ablat(interEn1,1)), mean(SCORE_Ablat(interEn1,2))];
center_En1=center;
stdev = [std(SCORE_Ablat(interEn1,1)), std(SCORE_Ablat(interEn1,2))];
stdev_En1=stdev;
axh = gca(); 
llc = center(:)-stdev(:);
wh = stdev(:)*2; 
rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1],'EdgeColor','c'); 
plot(center(1),center(2),'xc')


% % % % % %SHOX2
plot(SCORE_Ablat(interShox2,1),SCORE_Ablat(interShox2,2),'.m','MarkerSize',20)
plotell([SCORE_Ablat(interShox2,1) SCORE_Ablat(interShox2,2)],'-m','om');
center = [mean(SCORE_Ablat(interShox2,1)), mean(SCORE_Ablat(interShox2,2))];
center_Shox2=center;
stdev = [std(SCORE_Ablat(interShox2,1)), std(SCORE_Ablat(interShox2,2))];
stdev_Shox2=stdev;
axh = gca(); 
llc = center(:)-stdev(:);
wh = stdev(:)*2; 
rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1],'EdgeColor','m'); 
plot(center(1),center(2),'xm')
% 


D_yoked = pdist([center_Master;center_Yoked],'euclidean');
D_Ptf1a = pdist([center_Master;center_Ptf1a],'euclidean');
D_Tlx3 = pdist([center_Master;center_Tlx3],'euclidean');
D_DMRT3 = pdist([center_Master;center_DMRT3],'euclidean');
D_Sim1 = pdist([center_Master;center_Sim1],'euclidean');
D_Shox2 = pdist([center_Master;center_Shox2],'euclidean');
D_En1 = pdist([center_Master;center_En1],'euclidean');



%Centroid distance
%bar([D_yoked D_Ptf1a D_Tlx3 D_DMRT3 D_Sim1 D_Shox2 D_En1])

%%

clearvars ZZZZ Norm

DataPCA_origin=[Data_Day1_Master;Data_Day1_Yoked;Data_Day1_Ptf1a;Data_Day1_Tlx3;Data_Day1_DMRT3;Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2]%;Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2;Data_Day1_Yoked;Data_Day1_Master];%;Data_Day1_Yoked];Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2];


interMaster=[1:25];
interYoked=[26:50];
interPtf1a=[51:55];
interTlx3=[56:59];
interDMRT3=[60:64];
interSim1=[65:69];
interEn1=[70:74];
interShox2=[75:79];

inter=[110:114 115:119 120:124 ...
    125:129 130:134 135:139 ...
    179:183 184:188 189:193 ...
    197:201 202:206 207:211 ...
    194:196];

toshuffle=[interMaster];
tostudy=interSim1;

for k=1:1000
    
    inter_shuffle=toshuffle(randperm(length(toshuffle)));
    tostudy=toshuffle(randperm(length(toshuffle)));

    for i=1:63
    
        %h1=adtest(DataPCA_origin(inter_shuffle(1:5),inter(i)));
        %h2=adtest(DataPCA_origin(tostudy,inter(i)));
        %h=kstest(DataPCA_origin(inter_shuffle(1:5),inter(i)),DataPCA_origin(tostudy,inter(i)));
        x=DataPCA_origin(inter_shuffle(1:5),inter(i));
        y=DataPCA_origin(tostudy,inter(i));
        
        h1=adtest(x);
        h2=adtest(y);
        
        if h1==0 && h2==0
            [a b]=ttest2(DataPCA_origin(inter_shuffle(1:5),inter(i)),DataPCA_origin(tostudy,inter(i)));
            Norm(i,k)=1;
        else   
            [a b]=kstest2(DataPCA_origin(inter_shuffle(1:5),inter(i)),DataPCA_origin(tostudy,inter(i)));      
            Norm(i,k)=0;
        end
        
        if b<0.05
            ZZZZ(i,k)=1;
        else
            ZZZZ(i,k)=0;
        end
             
    end

end

ZZZZ_sum=sum(ZZZZ,2)/10;
Norm_sum=sum(Norm,2)/10;

figure(1)
bar(ZZZZ_sum)
figure(2)
bar(Norm_sum)



for i=1:13
    
    if i==13
        
    ZZZZ_inter(i)=mean(ZZZZ_sum((12*5):(12*5)+3));
    
    else
        
    %interval= 1+((i-1)*3) : (i*3);
    interval= 1+((i-1)*5) : (i*5)+3;
    
    ZZZZ_inter(i)=mean(ZZZZ_sum(interval));
    
    end
    
end

% close all
% 
% hold on
% plot(DataPCA_origin(inter_shuffle(1:5),200),1,'or')
% plot(DataPCA_origin(interPtf1a,200),2,'ob')
% hold off

%% Consecutive day of training PTF1A

DataPCA_origin=[   Data_Day1_Master;Data_Day1_Yoked; ...
            Data_Day1_Ptf1a;Data_Day2_Ptf1a; ...
            Data_Day3_Ptf1a;Data_Day4_Ptf1a; ...
            Data_Day5_Ptf1a];


clearvars StandKine SCORE_Ablat

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

inter=[8:12;13:17;18:22;41:45;46:50; ...
    51:55;80:84;85:89;90:94; ...
    95:99;100:104;105:109];

[a b]=size(inter);
TabOfIn=[];

for i=1:a
    TabOfIn=[];
    TabOfIn=DataPCA_origin(:,inter(i,:));
    
    [n_norm p_norm] = size(TabOfIn);
    e = ones(n_norm,1);

    StandKine=[StandKine (TabOfIn-min(min(TabOfIn)))./(max(max(TabOfIn))-min(min(TabOfIn)))];

end


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


clearvars COEFF_Ablat SCORE_Ablat
[COEFF_Ablat, SCORE_Ablat, LATENT_Ablat, TSQUARED_Ablat, EXPLAINED_Ablat]=pca(StandKine);%,'Centered',true);%,'VariableWeights','variance');


interMaster=[1:25];
interYoked=[26:50];
interPtf1a_Day1=[51:55];
interPtf1a_Day2=[56:60];
interPtf1a_Day3=[61:65];
interPtf1a_Day4=[66:70];
interPtf1a_Day5=[71:75];

%%

x=SCORE_Ablat(interMaster,1);
y=SCORE_Ablat(interTlx3_Day5,1);

h1=adtest(x);
h2=adtest(y);

if h1==0 && h2==0
    [a b]=ttest2(x,y);
    'ttest'
    %Norm(i,k)=1;
else
    [a b]=kstest2(x,y);
    'wilkinson'
    %Norm(i,k)=0;
end

b

%%

t = table(species,meas(:,1),meas(:,2),meas(:,3),meas(:,4),...
'VariableNames',{'species','meas1','meas2','meas3','meas4'});
Meas = table([1 2 3 4]','VariableNames',{'Measurements'});


%% Consecutive day of training TLX3

DataPCA_origin=[   Data_Day1_Master;Data_Day1_Yoked; ...
            Data_Day1_Tlx3;Data_Day2_Tlx3; ...
            Data_Day3_Tlx3;Data_Day4_Tlx3; ...
            Data_Day5_Tlx3];


clearvars StandKine SCORE_Ablat

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

inter=[8:12;13:17;18:22;41:45;46:50; ...
    51:55;80:84;85:89;90:94; ...
    95:99;100:104;105:109];

[a b]=size(inter);
TabOfIn=[];

for i=1:a
    TabOfIn=[];
    TabOfIn=DataPCA_origin(:,inter(i,:));
    
    [n_norm p_norm] = size(TabOfIn);
    e = ones(n_norm,1);

    StandKine=[StandKine (TabOfIn-min(min(TabOfIn)))./(max(max(TabOfIn))-min(min(TabOfIn)))];

end


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


clearvars COEFF_Ablat SCORE_Ablat
[COEFF_Ablat, SCORE_Ablat, LATENT_Ablat, TSQUARED_Ablat, EXPLAINED_Ablat]=pca(StandKine);%,'Centered',true);%,'VariableWeights','variance');


interMaster=[1:25];
interYoked=[26:50];
interTlx3_Day1=[51:54];
interTlx3_Day2=[55:58];
interTlx3_Day3=[59:62];
interTlx3_Day4=[63:66];
interTlx3_Day5=[67:70];































































%% ANALYSIS WITHIN ABLATED GROUP

clc

TableComp=zeros(8);

interMaster=[1:25];
interYoked=[26:50];
interPtf1a=[51:55];
interTlx3=[56:59];
interDMRT3=[60:64];
interSim1=[65:69];
interEn1=[70:74];
interShox2=[75:79];


interMice=[1 25;26 50;51 55;56 59;60 64;65 69;70 74;75 79];

DataPCA_origin=[ ...
    Data_Day1_Master;...
    Data_Day1_Yoked;...
    Data_Day1_Ptf1a;...
    Data_Day1_Tlx3;...
    Data_Day1_DMRT3;...
    Data_Day1_Sim1;...
    Data_Day1_En1;...
    Data_Day1_Shox2...
    ];

iter1etoile=0;
iter2etoile=0;
iter1=[];
iter2=[];
    
for j=1:8
    
    for k=j:8
        
        for i=1:63
            
            iter1=interMice(j,1):interMice(j,2);
            iter2=interMice(k,1):interMice(k,2);
            
            [a b]=ttest2(StandKine(iter1,i),StandKine(iter2,i));
            
            if b<0.05
                iter1etoile=iter1etoile+1;
            end
            if b<0.01
                iter2etoile=iter2etoile+1;
            end
            
        end
        
        TableComp(j,k)=iter1etoile;
     	TableComp(k,j)=iter2etoile;
        iter1etoile=0;
        iter2etoile=0;
   
    end
    
end

colormap cool
yvalues = {'Master','Yoked','Ptf1a','Tlx3','DMRT3','Sim1','En1','Shox2'};


h=heatmap(yvalues,yvalues,TableComp)
h.Colormap = parula

h.Title = 'Dynamic State parameters=72';
h.XLabel = 'p<0.05';
h.YLabel = 'p<0.01'; 





%%
clearvars AAA AAA2
% 
% AAA=[ ...
%     SCORE_Ablat(interMaster,1) ....
%     SCORE_Ablat(interYoked,1) ...
%     [SCORE_Ablat(interPtf1a,1);zeros(25-numel(interPtf1a),1)]...
%     [SCORE_Ablat(interTlx3,1);zeros(25-numel(interTlx3),1)]...
%     [SCORE_Ablat(interDMRT3,1);zeros(25-numel(interDMRT3),1)]...
%     [SCORE_Ablat(interSim1,1);zeros(25-numel(interSim1),1)]...
%     [SCORE_Ablat(interEn1,1);zeros(25-numel(interEn1),1)]...
%     [SCORE_Ablat(interShox2,1);zeros(25-numel(interShox2),1)]]
% 
% AAA2=[ ...
%     SCORE_Ablat(interMaster,2) ....
%     SCORE_Ablat(interYoked,2) ...
%     [SCORE_Ablat(interPtf1a,2);zeros(25-numel(interPtf1a),1)]...
%     [SCORE_Ablat(interTlx3,2);zeros(25-numel(interTlx3),1)]...
%     [SCORE_Ablat(interDMRT3,2);zeros(25-numel(interDMRT3),1)]...
%     [SCORE_Ablat(interSim1,2);zeros(25-numel(interSim1),1)]...
%     [SCORE_Ablat(interEn1,2);zeros(25-numel(interEn1),1)]...
%     [SCORE_Ablat(interShox2,2);zeros(25-numel(interShox2),1)]]


AAA_noyoked=[ ...
    SCORE_Ablat(interMaster,1) ....
    [SCORE_Ablat(interPtf1a,1);zeros(25-numel(interPtf1a),1)]...
    [SCORE_Ablat(interTlx3,1);zeros(25-numel(interTlx3),1)]...
    [SCORE_Ablat(interDMRT3,1);zeros(25-numel(interDMRT3),1)]...
    [SCORE_Ablat(interSim1,1);zeros(25-numel(interSim1),1)]...
    [SCORE_Ablat(interEn1,1);zeros(25-numel(interEn1),1)]...
    [SCORE_Ablat(interShox2,1);zeros(25-numel(interShox2),1)]]

AAA2_noyoked=[ ...
    SCORE_Ablat(interMaster,2) ....
    [SCORE_Ablat(interPtf1a,2);zeros(25-numel(interPtf1a),1)]...
    [SCORE_Ablat(interTlx3,2);zeros(25-numel(interTlx3),1)]...
    [SCORE_Ablat(interDMRT3,2);zeros(25-numel(interDMRT3),1)]...
    [SCORE_Ablat(interSim1,2);zeros(25-numel(interSim1),1)]...
    [SCORE_Ablat(interEn1,2);zeros(25-numel(interEn1),1)]...
    [SCORE_Ablat(interShox2,2);zeros(25-numel(interShox2),1)]]




%%
% clc
% close all
% cl = fitcsvm(X_training(:,2:3),X_training(:,1),'KernelFunction','gaussian',...
%     'BoxConstraint',Inf,'ClassNames',[0,1]);
% 
% data3=X_training(:,2:3);
% theclass=X_training(:,1);
% d = 0.02;
% [x1Grid,x2Grid] = meshgrid(min(data3(:,1)):d:max(data3(:,1)),...
%     min(data3(:,2)):d:max(data3(:,2)));
% xGrid = [x1Grid(:),x2Grid(:)];
% [~,scores] = predict(cl,xGrid);
% 
% % Plot the data and the decision boundary
% figure;
% h(1:2) = gscatter(data3(:,1),data3(:,2),theclass,'rb','.');
% hold on
% %ezpolar(@(x)1);
% h(3) = plot(data3(cl.IsSupportVector,1),data3(cl.IsSupportVector,2),'ko');
% contour(x1Grid,x2Grid,reshape(scores(:,2),size(x1Grid)),[0 0],'k');
% legend(h,{'-1','+1','Support Vectors'});
% axis equal
% hold off
% 
% %%
% X=X_training_Param(:,2:end);
% y=X_training_Param(:,1);
% 
% SVMModel = fitcsvm(X,y)
% 
% sv = SVMModel.SupportVectors;
% figure
% gscatter(X_training(:,2),X_training(:,3),y)
% hold on
% plot(sv(:,1),sv(:,2),'ko','MarkerSize',10)
% legend('versicolor','virginica','Support Vector')
% hold off

%% Individual distance to Master centroid


plot(SCORE_Ablat(25:48,1),SCORE_Ablat(25:48,2),'ok')
clearvars D_yoked_ind

for i=1:24
    
    X=[center_Master; SCORE_Ablat(i,1) SCORE_Ablat(i,2)]
    D_master_ind(i)=pdist(X,'euclidean');
    
    D_master_1(i)=abs(center_Master(1)-SCORE_Ablat(i,1));
    D_master_2(i)=abs(center_Master(2)-SCORE_Ablat(i,2));
    
end

for i=1:24
    
    X=[center_Master; SCORE_Ablat(24+i,1) SCORE_Ablat(24+i,2)]
    D_yoked_ind(i)=pdist(X,'euclidean');
    
    D_yoked_1(i)=abs(center_Master(1)-SCORE_Ablat(24+i,1));
    D_yoked_2(i)=abs(center_Master(2)-SCORE_Ablat(24+i,2));
    
end

for i=1:5
    
    X=[center_Master; SCORE_Ablat(48+i,1) SCORE_Ablat(48+i,2)]
    D_ptf1a_ind(i)=pdist(X,'euclidean');
    
    D_ptf1a_1(i)=abs(center_Master(1)-SCORE_Ablat(48+i,1));
    D_ptf1a_2(i)=abs(center_Master(2)-SCORE_Ablat(48+i,2));
    
end

for i=1:4
    
    X=[center_Master; SCORE_Ablat(53+i,1) SCORE_Ablat(53+i,2)]
    D_tlx3_ind(i)=pdist(X,'euclidean');
    
    D_tlx3_1(i)=abs(center_Master(1)-SCORE_Ablat(53+i,1));
    D_tlx3_2(i)=abs(center_Master(2)-SCORE_Ablat(53+i,2));
    
end

for i=1:5
    
    X=[center_Master; SCORE_Ablat(57+i,1) SCORE_Ablat(57+i,2)]
    D_dmrt3_ind(i)=pdist(X,'euclidean');
    
    D_dmrt3_1(i)=abs(center_Master(1)-SCORE_Ablat(57+i,1));
    D_dmrt3_2(i)=abs(center_Master(2)-SCORE_Ablat(57+i,2));
    
    
end

for i=1:5
    
    X=[center_Master; SCORE_Ablat(62+i,1) SCORE_Ablat(62+i,2)]
    D_sim1_ind(i)=pdist(X,'euclidean');
    
    D_sim1_1(i)=abs(center_Master(1)-SCORE_Ablat(62+i,1));
    D_sim1_2(i)=abs(center_Master(2)-SCORE_Ablat(62+i,2));
    
end

for i=1:5
    
    X=[center_Master; SCORE_Ablat(67+i,1) SCORE_Ablat(67+i,2)]
    D_en1_ind(i)=pdist(X,'euclidean');
    
    D_en1_1(i)=abs(center_Master(1)-SCORE_Ablat(67+i,1));
    D_en1_2(i)=abs(center_Master(2)-SCORE_Ablat(67+i,2));
    
end

for i=1:5
    
    X=[center_Master; SCORE_Ablat(72+i,1) SCORE_Ablat(72+i,2)]
    D_shox2_ind(i)=pdist(X,'euclidean');
    
    D_shox2_1(i)=abs(center_Master(1)-SCORE_Ablat(72+i,1));
    D_shox2_2(i)=abs(center_Master(2)-SCORE_Ablat(72+i,2));
    
end


%%

Cat=[zeros(24,1);ones(24,1);2*ones(5,1);3*ones(4,1);4*ones(5,1);5*ones(5,1);6*ones(5,1);7*ones(5,1)];

Cat_cell = num2cell(Cat)

scatterhist(SCORE_Ablat(1:77,1),SCORE_Ablat(1:77,2),'Group',Cat,'Kernel','on');
figure(gcf)
hold on

%%

h = scatterhist(SCORE_Ablat(1:77,1),SCORE_Ablat(1:77,2),'Group',Cat,'Kernel','on');
hold on;
clr = get(h(1),'colororder');
boxplot(h(2),SCORE_Ablat(1:77,1),Cat,'orientation','horizontal','color',clr);
boxplot(h(3),SCORE_Ablat(1:77,2),Cat,'orientation','horizontal','color',clr);
set(h(2:3),'XTickLabel','');
view(h(3),[270,90]);  % Rotate the Y plot
axis(h(1),'auto');  % Sync axes

center = [mean(SCORE_Ablat(1:24,1)), mean(SCORE_Ablat(1:24,2))];
stdev = [std(SCORE_Ablat(1:24,1)), std(SCORE_Ablat(1:24,2))];
llc = center(:)-stdev(:);
axh = gca(); 
wh = stdev(:)*2; 
rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1],'EdgeColor',clr(1,:)); 

center = [mean(SCORE_Ablat(25:48,1)), mean(SCORE_Ablat(25:48,2))];
stdev = [std(SCORE_Ablat(25:48,1)), std(SCORE_Ablat(25:48,2))];
llc = center(:)-stdev(:);
axh = gca(); 
wh = stdev(:)*2; 
rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1],'EdgeColor',clr(2,:)); 

center = [mean(SCORE_Ablat(49:53,1)), mean(SCORE_Ablat(49:53,2))];
stdev = [std(SCORE_Ablat(49:53,1)), std(SCORE_Ablat(49:53,2))];
llc = center(:)-stdev(:);
axh = gca(); 
wh = stdev(:)*2; 
rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1],'EdgeColor',clr(3,:)); 

center = [mean(SCORE_Ablat(54:57,1)), mean(SCORE_Ablat(54:57,2))];
stdev = [std(SCORE_Ablat(54:57,1)), std(SCORE_Ablat(54:57,2))];
llc = center(:)-stdev(:);
axh = gca(); 
wh = stdev(:)*2; 
rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1],'EdgeColor',clr(4,:)); 

center = [mean(SCORE_Ablat(58:62,1)), mean(SCORE_Ablat(58:62,2))];
stdev = [std(SCORE_Ablat(58:62,1)), std(SCORE_Ablat(58:62,2))];
llc = center(:)-stdev(:);
axh = gca(); 
wh = stdev(:)*2; 
rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1],'EdgeColor',clr(5,:)); 

center = [mean(SCORE_Ablat(63:67,1)), mean(SCORE_Ablat(63:67,2))];
stdev = [std(SCORE_Ablat(63:67,1)), std(SCORE_Ablat(63:67,2))];
llc = center(:)-stdev(:);
axh = gca(); 
wh = stdev(:)*2; 
rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1],'EdgeColor',clr(6,:)); 

center = [mean(SCORE_Ablat(68:72,1)), mean(SCORE_Ablat(68:72,2))];
stdev = [std(SCORE_Ablat(68:72,1)), std(SCORE_Ablat(68:72,2))];
llc = center(:)-stdev(:);
axh = gca(); 
wh = stdev(:)*2; 
rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1],'EdgeColor',clr(7,:)); 

center = [mean(SCORE_Ablat(72:77,1)), mean(SCORE_Ablat(72:77,2))];
stdev = [std(SCORE_Ablat(72:77,1)), std(SCORE_Ablat(72:77,2))];
llc = center(:)-stdev(:);
axh = gca(); 
wh = stdev(:)*2; 
rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1],'EdgeColor',clr(1,:)); 






hold off;



%%
































%% ABLATION SORT >50% 

% 
% clearvars Data_Day1_Lbx1 Data_Day1_DMRT3 Data_Day1_En1 Data_Day1_Ptf1a Data_Day1_Shox2 Data_Day1_Sim1 Data_Day1_Tlx3
%  
% % Master Data_Day1_Master
% 
% %Extract ablation
% % 2 -- Lbx1
% % 3 -- Ptf1a
% % 4 -- Tlx3
% % 5 -- DMRT3
% % 6 -- En1
% % 7 -- Shox2
% % 8 -- Sim1
% 
% Ablat=2;
% 
% iterLbx1=1;
% iterPtf1a=1;
% iterTlx3=1;
% iterDMRT3=1;
% iterEn1=1;
% iterShox2=1;
% iterSim1=1;
% 
% iterLbx1_2=1;
% iterLbx1_3=1;
% iterLbx1_4=1;
% iterLbx1_5=1;
% 
% iterPtf1a_2=1;
% iterPtf1a_3=1;
% iterPtf1a_4=1;
% iterPtf1a_5=1;
% 
% iterTlx3_2=1;
% iterTlx3_3=1;
% iterTlx3_4=1;
% iterTlx3_5=1;
% 
% [n p]=size(Data);
% 
% for i=1:n
% 
% 	if ((Data(i,1)==4&&Data(i,5)==5) || (Data(i,1)==22&&Data(i,5)==1) || (Data(i,1)==22&&Data(i,5)==2) || (Data(i,1)==25&&Data(i,5)==1) || (Data(i,1)==25&&Data(i,5)==3)) && Data(i,3)==1 && Data(i,4)==1
%         Data_Day1_Ptf1a(iterPtf1a,:)=Data(i,:);
%     end
%     
% 	if ((Data(i,1)==9&&Data(i,5)==1) || (Data(i,11)==22&&Data(i,9)==1) || (Data(i,1)==22&&Data(i,5)==2) || (Data(i,1)==25&&Data(i,5)==1) || (Data(i,1)==25&&Data(i,5)==3)) && Data(i,3)==1 && Data(i,4)==1
%         Data_Day1_Ptf1a(iterPtf1a,:)=Data(i,:);
%     end
%     
%     
% 
%   
% end
% 
% 
% 



%% FIGURE 3B -  dyn and stat


%% Parameter analysis

Param=[80:94 95:109 110:124 125:129 130:134 135:139 140:154 155:157]; %
%Param=[6 7 8:22 23:25 26:40 41:55 56:58 59:61 62:64 65:79]; 

% 1:15 Amp Spike
% 16:30 Speed
% 31:45  Acc
% 46:50 Max spike
% 51:55 Min spike
% 56:60 Endpoint 
% 61:75 Absement
% 76:78 Angle


% iter=1;
% 
% figure(1)
% clearvars TTEST_Abl
% 
% for i=Param
%     
%     TabToScale=abs([Data_Day1_Lbx1(:,i);Data_Day1_Master(:,i)]);
%     TabToScale=(TabToScale);
%     
%     TabMean1(iter)=mean(TabToScale(1:5,1));
%     TabMean5(iter)=mean(TabToScale(6:31,1));
% 
%     %TabD15(iter)=mean( TabToScale(1:5,1) ./ TabToScale(6:31,1) )-1;
%     
%     TabSTD1(iter)=std(TabToScale(1:5,1))/sqrt(5);
%     TabSTD5(iter)=std(TabToScale(6:31,1))/sqrt(26);
% 
%     TabDivide(iter)=TabMean1(iter)/TabMean5(iter)-1;
%     
% 
%     [a b]=ttest2(TabToScale(1:5,1),TabToScale(6:31,1));
%     TTEST_Abl(iter)=(b);
%     
%     iter=iter+1;
%     
%     
% end
% 
%  %hold on 
%     %bar(TabDivide)
%     %plot(TabSTD5,'bx')
%     
%     %axis([0 length(Param) -1 1])
% % errorbar(TabMean1,TabSTD1)
% % hold on
% % errorbar(TabMean5,TabSTD5)
% 
% % 1:15 Amp Spike
% % 16:30 Speed
% % 31:45  Acc
% % 46:50 Max spike
% % 51:55 Min spike
% % 56:60 Endpoint 
% % 61:75 Absement
% % 76:78 Angle
% 
% %To put in a hexagone
% CatDivi(1)=mean(abs(TabDivide(1:15)))
% CatDivi(2)=mean(abs(TabDivide(16:30)))
% CatDivi(3)=mean(abs(TabDivide(31:45)))
% CatDivi(4)=mean(abs(TabDivide(46:60)))
% CatDivi(5)=mean(abs(TabDivide(61:75)))
% CatDivi(6)=mean(abs(TabDivide(76:78)))
% 
% CatDivi_std(1)=std(abs(TabDivide(1:15)))/sqrt(15)
% CatDivi_std(2)=std(abs(TabDivide(16:30)))/sqrt(15)
% CatDivi_std(3)=std(abs(TabDivide(31:45)))/sqrt(15)
% CatDivi_std(4)=std(abs(TabDivide(46:60)))/sqrt(15)
% CatDivi_std(5)=std(abs(TabDivide(61:75)))/sqrt(15)
% CatDivi_std(6)=std(abs(TabDivide(76:78)))/sqrt(3)
% 
% errorbar(CatDivi,CatDivi_std)


















%%
%% Parameter analysis GOOD versus BAD
%% COMPARE INDIV POPULATION VS. MASTER

%Param=[80:94 95:109 110:124 125:129 130:134 135:139 140:154 155:157]; %
%Param=[6 7 8:22 23:25 26:40 41:55 56:58 59:61 62:64 65:79]; 

%Param_kin
%Param=[8:22 26:40 41:55 56:79 80:94 95:109 110:124 125:129 130:134 135:139 140:154 23:25 155:157];

% Param_inter_stat=[5 5 5  5 5 5  5 5 5  3 3 3  5 5 5  3]; %16
% Param_inter_dyn=[5 5 5  5 5 5  5 5 5  5 5 5  5 5 5  3]; %16
% Param_inter=[Param_inter_stat Param_inter_dyn];

% 
% Param=[8:10 13:15 18:20  26:28 31:33 36:38  41:43 46:48 51:53 ... 
%         56:58 59:61 62:64  65:67 70:72 75:77  23:25 ...  %16 * FAK
%         80:82 85:87 90:92   95:97 100:102 105:107 ... 
%         110:112 115:117 120:122  125:127 130:132 135:137 ...
%         140:142 145:147 150:152  155:157]; %16 * FAK
% 
% Param_inter_stat=[3 3 3  3 3 3  3 3 3  3 3 3  3 3 3]; %15
% Param_inter_dyn=[3 3 3  3 3 3  3 3 3  3 3 3]; %12

Param=[8:10 13:15 18:20  26:28 31:33 36:38  41:43 46:48 51:53 ... 
         56:58 59:61 62:64  65:67 70:72 75:77  23:25];   %16 * FAK
        
Param=[ Param 80:82 85:87 90:92   95:97 100:102 105:107 ... 
        110:112 115:117 120:122  125:127 130:132 135:137 ...
        140:142 145:147 150:152  155:157]; %16 * FAK


% Dyn
% 1:15 Amp Spike
% 16:30 Speed
% 31:45  Acc
% 46:50 Max spike
% 51:55 Min spike
% 56:60 Endpoint 
% 61:75 Absement
% 76:78 Angle

% Bad=[4 6 8 12 15 19 20 21 22 26];
% Good=[1 2 3 5 7 9 10 11 13 14 16 17 18 23 24 25];

Bad=[4 7 8 11 14 18 20 21 24];
Good=[1 2 3 5 6 9 10 12 13 15 16 17 19 22 23 25];


iter=1;
iter_signi=1;

figure(1)
clearvars TTEST_Abl TabMean1 TabMean5 TabSTD1 TabSTD5 TabDivide 
clearvars TabMean1_Signi TabMean5_Signi TabSTD1_Signi TabSTD5_Signi Tab_Signi_Impact TabSigniWho
clc

hold on
for i=Param
    
    A=Data_Day1_Master(Bad,i);
    B=Data_Day1_Master(Good,i);
    
    TabToScale2=([A;B]); % [10*Bad 16*Good]
    %TabToScale=(TabToScale2);
    TabToScale=(TabToScale2-min(TabToScale2)) / (max(TabToScale2)-min(TabToScale2));
    
    TabMean1(iter)=mean(TabToScale(1:length(A),1)); %Bad
    TabMean5(iter)=mean(TabToScale((length(A)+1):(length(B)+length(A)),1)); %Good
    
    TabSTD1(iter)=std(TabToScale(1:length(A),1))/sqrt(length(A)-1);
    TabSTD5(iter)=std(TabToScale((length(A)+1):(length(B)+length(A)),1))/sqrt(length(B)-1);

    TabDivide(iter)=(TabMean1(iter)/TabMean5(iter)); %Bad/Good
    
    [a b]=ttest2(TabToScale(1:length(A),1),TabToScale((length(A)+1):(length(B)+length(A)),1));
    TTEST_Abl(iter)=(b);
    
    if b<0.05
        
    TabMean1_Signi(iter_signi)=TabMean1(iter);
    TabMean5_Signi(iter_signi)=TabMean5(iter);
    TabSTD1_Signi(iter_signi)=TabSTD1(iter);
    TabSTD5_Signi(iter_signi)=TabSTD5(iter);
    Tab_Signi_Impact(iter_signi)=b;
    TabSigniWho(iter_signi)=NewPCAfile092021.Properties.VariableNames(([i]));
    
    iter_signi=iter_signi+1;
    
    end
    
    iter=iter+1;
    
%     0.05
%     0.01
%     0.001
end


errorbar(TabMean1_Signi,TabSTD1_Signi)
errorbar(TabMean5_Signi,TabSTD5_Signi)

%axis([-1 length(TabMean1_Signi)+1 0 1])

figure(2)
Impact = [TabMean1_Signi 0];
theta = linspace(0,2*pi,(length(Impact)));
polarplot(theta,Impact,'-r');
polarplot(theta,Impact,'or');
hold on

Impact = [TabMean1_Signi-TabSTD1_Signi 0];
theta = linspace(0,2*pi,(length(Impact)));
polarplot(theta,Impact,'--r');
polarplot(theta,Impact,'xr');

Impact = [TabMean1_Signi+TabSTD1_Signi 0];
theta = linspace(0,2*pi,(length(Impact)));
polarplot(theta,Impact,'--r');
polarplot(theta,Impact,'xr');

Impact = [TabMean5_Signi 0];
theta = linspace(0,2*pi,(length(Impact)));
polarplot(theta,Impact,'ob');
polarplot(theta,Impact,'-b');
hold on

Impact = [TabMean5_Signi-TabSTD5_Signi 0];
theta = linspace(0,2*pi,(length(Impact)));
polarplot(theta,Impact,'--b');
polarplot(theta,Impact,'xb');

Impact = [TabMean5_Signi+TabSTD5_Signi 0];
theta = linspace(0,2*pi,(length(Impact)));
polarplot(theta,Impact,'--b');
polarplot(theta,Impact,'xb');

ax=gca
ax.ThetaTick = [0:(360/length(TabMean5_Signi)):(360)]; 
ax.ThetaTickLabel=TabSigniWho;



%%
%% COMPARE RAW DATA OF GROUP OF INTEREST.
%%

Param=[80:94 95:109 110:124 125:129 130:134 135:139 140:154 155:157]; %
%Param=[6 7 8:22 23:25 26:40 41:55 56:58 59:61 62:64 65:79]; 

%Param_kin
%Param=[8:22 26:40 41:55 56:79 80:94 95:109 110:124 125:129 130:134 135:139 140:154 23:25 155:157];

% Param_inter_stat=[5 5 5  5 5 5  5 5 5  3 3 3  5 5 5  3]; %16
% Param_inter_dyn=[5 5 5  5 5 5  5 5 5  5 5 5  5 5 5  3]; %16
% Param_inter=[Param_inter_stat Param_inter_dyn];

% 
%  Param=[8:10 13:15 18:20  26:28 31:33 36:38  41:43 46:48 51:53 ... 
%         56:58 59:61 62:64  65:67 70:72 75:77  23:25 ...  %16 * FAK
%         80:82 85:87 90:92   95:97 100:102 105:107 ... 
%         110:112 115:117 120:122  125:127 130:132 135:137 ...
%         140:142 145:147 150:152  155:157]; %16 * FAK

% %Dyn Param
% Param=[ 80:82 85:87 90:92   95:97 100:102 105:107 ... 
%         110:112 115:117 120:122  125:127 130:132 135:137 ...
%         140:142 145:147 150:152  155:157]; %16*3 FAK  % 5*9 +3


% 
% Param_inter_stat=[3 3 3  3 3 3  3 3 3  3 3 3  3 3 3]; %15
% Param_inter_dyn=[3 3 3  3 3 3  3 3 3  3 3 3]; %12


% Dyn
% 1:15 Amp Spike
% 16:30 Speed
% 31:45  Acc
% 46:50 Max spike
% 51:55 Min spike
% 56:60 Endpoint 
% 61:75 Absement
% 76:78 Angle

%
%vbls(i)=BookPCAv6.Properties.VariableNames(Param(I(i)));
clearvars A B C D E F TabToScale TabToScalebe
clearvars TabMean1_SigniEn1 TabMean1_SigniLbx1 TabMean1_SigniSim1 TabMean1_SigniShox2
clearvars TabMean5_SigniEn1 TabMean5_SigniLbx1 TabMean5_SigniSim1 TabMean5_SigniShox2
clearvars TabMean5 TabSTD1 TabSTD1_SigniEn1 TabSTD1_SigniLbx1 TabSTD1_SigniShox2 TabSTD1_SigniSim1 TabSTD5 TabSTD5_SigniEn1 TabSTD5_SigniLbx1 TabSTD5_SigniShox2 TabSTD5_SigniSim1
clearvars Tab_Signi_ImpactLbx1 Tab_Signi_ImpactEn1 Tab_Signi_ImpactSim1 Tab_Signi_ImpactShox2
clc

iter=1;
for i=Param
    
    A=Data_Day1_Lbx1(:,i);
    B=Data_Day1_En1(:,i);
    C=Data_Day1_Sim1(:,i);
    D=Data_Day1_Shox2(:,i);
    
    E=Data_Day1_Master(:,i);
    
    F=Data_Day1_Yoked(:,i);
    
    TabToScalebe=([A;B;C;D;E]);
    TabToScale(:,i)=zscore(TabToScalebe);
    %TabToScale(:,i)=(TabToScalebe-min(TabToScalebe)) / (max(TabToScalebe)-min(TabToScalebe));
    %TabToScale_Cluster(1:length(TabToScalebe),iter)=(TabToScalebe-min(TabToScalebe)) / (max(TabToScalebe)-min(TabToScalebe));
    TabToScale_Cluster(1:length(TabToScalebe),iter)=(TabToScalebe);
    
    iter=iter+1;
end


clearvars MasterCluster Lbx1Cluster Lbx1Cluster2 Lbx1Cluster_STD2 Tab_Lbx1VSMaster
clearvars Tab_Lbx1VSMaster Tab_Master

Tab_Lbx1VSMaster(1,:)=abs(mean(TabToScale_Cluster(1:5,:)));
Tab_En1VSMaster(1,:)=abs(mean(TabToScale_Cluster(6:10,:)));
Tab_Sim1VSMaster(1,:)=abs(mean(TabToScale_Cluster(11:15,:)));
Tab_Shox2VSMaster(1,:)=abs(mean(TabToScale_Cluster(16:20,:)));

Tab_Master(1,:)=abs(mean(TabToScale_Cluster(21:46,:)));

    
iter=1;
for i=1:6
    
    if i==6
    inter = 1+(i-1)*3 : (i)*3 ;
    inter=inter+4*15;
    else
    inter = 1+(i-1)*15 : (i)*15 ;    
    end
    
    Tab_Lbx1VSMaster_cluster(iter)=mean(Tab_Lbx1VSMaster(1,inter));
    Tab_En1VSMaster_cluster(iter)=mean( Tab_En1VSMaster(1,inter));
    Tab_Sim1VSMaster_cluster(iter)=mean(Tab_Sim1VSMaster(1,inter));
    Tab_Shox2VSMaster_cluster(iter)=mean(Tab_Shox2VSMaster(1,inter));
    
	Tab_Master_cluster(iter)=mean(Tab_Master(1,inter));
    
    iter=iter+1;
 
end

    Lbx1VSMaster=Tab_Lbx1VSMaster_cluster./Tab_Master_cluster
    En1VSMaster=Tab_En1VSMaster_cluster./Tab_Master_cluster
    Sim1VSMaster=Tab_Sim1VSMaster_cluster./Tab_Master_cluster
    Shox2VSMaster=Tab_Shox2VSMaster_cluster./Tab_Master_cluster


Impact = [Lbx1VSMaster 0];
theta = linspace(0,2*pi,(length(Impact)));
polarplot(theta,Impact,'or');
hold on   
Impact = [En1VSMaster 0];
theta = linspace(0,2*pi,(length(Impact)));
polarplot(theta,Impact,'ob');

Impact = [Sim1VSMaster 0];
theta = linspace(0,2*pi,(length(Impact)));
polarplot(theta,Impact,'og');

Impact = [Shox2VSMaster 0];
theta = linspace(0,2*pi,(length(Impact)));
polarplot(theta,Impact,'ok');


ax=gca
ax.ThetaTick = [0:(360/length(Lbx1VSMaster)):(360)]; 
    

%% COMPARE ALL DATA INDIV VS. MASTER AND 
%%  DISPLAY THE PVALUE < 0.05 FOR EACH ON THE SAME GRAPH
% % %%
% % 
% % iter=1;
% % for i=Param
% %     
% %     A=Data_Day1_Ptf1a(:,i);
% %     C=Data_Day1_DMRT3(:,i);
% %     B=Data_Day1_Tlx3(:,i);
% %     
% %     D=Data_Day1_Master(:,i);
% %     
% %     F=Data_Day1_Yoked(:,i);
% %     
% %     TabToScalebe=([A;B;C;D]);
% %     %TabToScale(:,i)=zscore(TabToScalebe);
% %     %TabToScale(:,i)=(TabToScalebe-min(TabToScalebe)) / (max(TabToScalebe)-min(TabToScalebe));
% %     %TabToScale_Cluster(1:length(TabToScalebe),iter)=(TabToScalebe-min(TabToScalebe)) / (max(TabToScalebe)-min(TabToScalebe));
% %     TabToScale_Cluster(1:length(TabToScalebe),iter)=(TabToScalebe);
% %     
% %     iter=iter+1;
% % end
% %
% % 
% % clearvars MasterCluster Lbx1Cluster Lbx1Cluster2 Lbx1Cluster_STD2 Tab_Lbx1VSMaster
% % clearvars Tab_Lbx1VSMaster Tab_Master
% % 
% % Tab_Ptf1aVSMaster(1,:)=abs(mean(TabToScale_Cluster(1:5,:)));
% % Tab_Tlx3VSMaster(1,:)=abs(mean(TabToScale_Cluster(6:9,:)));
% % Tab_DMRT3VSMaster(1,:)=abs(mean(TabToScale_Cluster(10:14,:)));
% % %Tab_Shox2VSMaster(1,:)=abs(mean(TabToScale_Cluster(16:20,:)));
% % 
% % Tab_Master(1,:)=abs(mean(TabToScale_Cluster(15:40,:)));
% % 
% %     
% % iter=1;
% % for i=1:6
% %     
% %     if i==6
% %     inter = 1+(i-1)*3 : (i)*3 ;
% %     inter=inter+4*15;
% %     else
% %     inter = 1+(i-1)*15 : (i)*15 ;    
% %     end
% %     
% %     Tab_Ptf1aVSMaster_cluster(iter)=mean(Tab_Ptf1aVSMaster(1,inter));
% %     Tab_Tlx3VSMaster_cluster(iter)=mean( Tab_Tlx3VSMaster(1,inter));
% %     Tab_DMRT3VSMaster_cluster(iter)=mean(Tab_DMRT3VSMaster(1,inter));
% %     %Tab_Shox2VSMaster_cluster(iter)=mean(Tab_Shox2VSMaster(1,inter));
% %     
% % 	Tab_Master_cluster(iter)=mean(Tab_Master(1,inter));
% %     
% %     iter=iter+1;
% %  
% % end
% % 
% %     Ptf1aVSMaster=Tab_Ptf1aVSMaster_cluster./Tab_Master_cluster
% %     Tlx3VSMaster=Tab_Tlx3VSMaster_cluster./Tab_Master_cluster
% %     DMRT3VSMaster=Tab_DMRT3VSMaster_cluster./Tab_Master_cluster
% %     %Shox2VSMaster=Tab_Shox2VSMaster_cluster./Tab_Master_cluster
% % 
% % 
% % Impact = [Ptf1aVSMaster 0];
% % theta = linspace(0,2*pi,(length(Impact)));
% % polarplot(theta,Impact,'or');
% % hold on   
% % Impact = [Tlx3VSMaster 0];
% % theta = linspace(0,2*pi,(length(Impact)));
% % polarplot(theta,Impact,'ob');
% % 
% % Impact = [DMRT3VSMaster 0];
% % theta = linspace(0,2*pi,(length(Impact)));
% % polarplot(theta,Impact,'og');
% % 
% % % Impact = [Shox2VSMaster 0];
% % % theta = linspace(0,2*pi,(length(Impact)));
% % % polarplot(theta,Impact,'ok');
% % 
% % 
% % ax=gca
% % ax.ThetaTick = [0:(360/length(Lbx1VSMaster)):(360)]; 
% %     
% % %rlim([0 2.5])
% %     
% % %%
% % 
% % 
% % Tab_Lbx1VSMaster=mean(TabToScale_Cluster(1:5,:))./mean(TabToScale_Cluster(6:31,:));
% % Tab_En1VSMaster=mean(TabToScale_Cluster(6:10,:))./mean(TabToScale_Cluster(21:46,:));
% % Tab_Sim1VSMaster=mean(TabToScale_Cluster(11:15,:))./mean(TabToScale_Cluster(21:46,:));
% % Tab_Shox2VSMaster=mean(TabToScale_Cluster(16:20,:))./mean(TabToScale_Cluster(21:46,:));
% % 
% % iter=1;
% % for i=1:6
% %     
% %     if i==6
% %     inter = 1+(i-1)*3 : (i)*3 ;
% %     inter=inter+4*15;
% %     else
% %     inter = 1+(i-1)*15 : (i)*15 ;    
% %     end
% %     
% %     
% %     Lbx1VSMaster(iter)=mean(Tab_Lbx1VSMaster(inter));
% %     Lbx1VSMaster_std(iter)=std(Tab_Lbx1VSMaster(inter));
% %     
% %     En1VSMaster(iter)=mean(Tab_En1VSMaster(inter));
% %     En1VSMaster_std(iter)=std(Tab_En1VSMaster(inter));
% %     
% %     Sim1VSMaster(iter)=mean(Tab_Sim1VSMaster(inter));
% %     Sim1VSMaster_std(iter)=std(Tab_Sim1VSMaster(inter));
% %     
% %     Shox2VSMaster(iter)=mean(Tab_Shox2VSMaster(inter));
% %     Shox2VSMaster_std(iter)=std(Tab_Shox2VSMaster(inter));
% %     
% %     iter=iter+1;
% %     
% % end
% % 
% % %%
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % iter=1;
% % iter_signi=1;
% % 
% % clearvars TTEST_Abl TabMean1 TabMean5 TabSTD1 TabSTD5 TabDivide 
% % clearvars TabMean1_Signi TabMean5_Signi TabSTD1_Signi TabSTD5_Signi Tab_Signi_Impact TabSigniWho TabParamInMemory
% % clearvars TabParamInMemoryEn1 TabParamInMemoryLbx1 TabParamInMemoryShox2 TabParamInMemorySim1
% % clearvars Shox2_control Shox2_control_std
% % clearvars En1_control En1_control_std
% % clearvars Lbx1_control Lbx1_control_std
% % clearvars Sim1_control Sim1_control_std
% % clearvars Master_control Master_control_std
% % clc
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % for i=Param
% %     
% % 
% %     TabMean1(iter)=mean(TabToScale(1:length(A),i));
% %     TabMean5(iter)=mean(TabToScale(21:46,i));
% %     
% %     TabSTD1(iter)=std(TabToScale(1:length(A),i))/sqrt(length(A));
% %     TabSTD5(iter)=std(TabToScale(21:46,i))/sqrt(length(21:46));
% % 
% %     TabDivide(iter)=(TabMean1(iter)/TabMean5(iter));
% %     
% % 
% %     [a b]=ttest2(TabToScale(1:length(A),i),TabToScale(21:46,i));
% %     TTEST_Abl(iter)=(b);
% %     
% %     if b<0.05
% %         
% %     TabMean1_SigniLbx1(iter_signi)=TabMean1(iter);
% %     TabMean5_SigniLbx1(iter_signi)=TabMean5(iter);
% %     TabSTD1_SigniLbx1(iter_signi)=TabSTD1(iter);
% %     TabSTD5_SigniLbx1(iter_signi)=TabSTD5(iter);
% %     Tab_Signi_ImpactLbx1(iter_signi)=b;
% %     TabSigniWho(iter_signi)=BookPCAv10.Properties.VariableNames(([i]));
% %     TabParamInMemoryLbx1(iter_signi)=i;
% %     
% %     iter_signi=iter_signi+1;
% %     
% %     end
% %     
% %     iter=iter+1;
% %     
% %     
% % end
% % 
% % iter=1;
% % iter_signi=1;
% % 
% % for i=Param
% %       
% %     interEn1=length(A)+1:length(A)+length(B);
% %     
% %     TabMean1(iter)=mean(TabToScale(interEn1,i));
% %     
% %     TabSTD1(iter)=std(TabToScale(interEn1,i))/sqrt(length(interEn1));    
% % 
% %     [a b]=ttest2(TabToScale(interEn1,i),TabToScale(21:46,i));
% %     TTEST_Abl(iter)=(b);
% %     
% %     if b<0.05
% %         
% %     TabMean1_SigniEn1(iter_signi)=TabMean1(iter);
% %     TabMean5_SigniEn1(iter_signi)=TabMean5(iter);
% %     TabSTD1_SigniEn1(iter_signi)=TabSTD1(iter);
% %     TabSTD5_SigniEn1(iter_signi)=TabSTD5(iter);
% %     Tab_Signi_ImpactEn1(iter_signi)=b;
% %     TabSigniWho(iter_signi)=BookPCAv10.Properties.VariableNames(([i]));
% %     TabParamInMemoryEn1(iter_signi)=i;
% % 
% %     
% %     iter_signi=iter_signi+1;
% %     
% %     end
% %     
% %     iter=iter+1;
% %     
% %     
% % end
% % 
% % iter=1;
% % iter_signi=1;
% % 
% % for i=Param
% %     
% %     interSim1=1+length(A)+length(B) : length(A)+length(B)+length(C);
% %     
% %     TabMean1(iter)=mean(TabToScale(interSim1,i));
% %     
% %     TabSTD1(iter)=std(TabToScale(interSim1,i))/sqrt(length(interSim1));    
% % 
% %     [a b]=ttest2(TabToScale(interSim1,i),TabToScale(21:46,i));
% %     TTEST_Abl(iter)=(b);
% %     
% %     if b<0.05
% %         
% %     TabMean1_SigniSim1(iter_signi)=TabMean1(iter);
% %     TabMean5_SigniSim1(iter_signi)=TabMean5(iter);
% %     TabSTD1_SigniSim1(iter_signi)=TabSTD1(iter);
% %     TabSTD5_SigniSim1(iter_signi)=TabSTD5(iter);
% %     Tab_Signi_ImpactSim1(iter_signi)=b;
% %     TabSigniWho(iter_signi)=BookPCAv10.Properties.VariableNames(([i]));
% %     TabParamInMemorySim1(iter_signi)=i;
% %     
% %     iter_signi=iter_signi+1;
% %     
% %     end
% %     
% %     iter=iter+1;
% %     
% %     
% % end
% % 
% % 
% % iter=1;
% % iter_signi=1;
% % for i=Param
% %     
% %     interShox2= 1+length(A)+length(B)+length(C) : length(A)+length(B)+length(C)+length(D);
% %     
% %     TabMean1(iter)=mean(TabToScale(interShox2,i));
% % 
% %     TabSTD1(iter)=std(TabToScale(interShox2,i))/sqrt(length(interShox2));
% % 
% %     
% % 
% %     [a b]=ttest2(TabToScale(interShox2,i),TabToScale(21:46,i));
% %     TTEST_Abl(iter)=(b);
% %     
% %     if b<0.05
% %         
% %     TabMean1_SigniShox2(iter_signi)=TabMean1(iter);
% %     TabMean5_SigniShox2(iter_signi)=TabMean5(iter);
% %     TabSTD1_SigniShox2(iter_signi)=TabSTD1(iter);
% %     TabSTD5_SigniShox2(iter_signi)=TabSTD5(iter);
% %     Tab_Signi_ImpactShox2(iter_signi)=b;
% %     TabSigniWho(iter_signi)=BookPCAv10.Properties.VariableNames(([i]));
% %     TabParamInMemoryShox2(iter_signi)=i;
% % 
% %     iter_signi=iter_signi+1;
% %     
% %     end
% %     
% %     iter=iter+1;
% %     
% %     
% % end
% % 
% % 
% % TabParamInMemory=[TabParamInMemoryLbx1 TabParamInMemoryEn1 TabParamInMemorySim1 TabParamInMemoryShox2]
% % C=unique(TabParamInMemory)
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % % errorbar(TabMean1_Signi,TabSTD1_Signi)
% % % errorbar(TabMean5_Signi,TabSTD5_Signi)
% % % axis([-1 length(TabMean1_Signi)+1 0 1])
% % 
% % 
% % figure(2)
% % 
% % iter=1;
% % for i=C
% %     
% % TabSigniWhoAblationControl(iter)=BookPCAv10.Properties.VariableNames(([i]));
% % iter=iter+1;
% %     
% % end
% % 
% % 
% % iter=1;
% % iterL=1;
% % iterE=1;
% % iterSi=1;
% % iterSh=1;
% % 
% % for i=C
% %     
% %     if ismember(i,TabParamInMemoryLbx1)
% %         Lbx1_control(iter)=TabMean1_SigniLbx1(iterL);
% %         Lbx1_control_std(iter)=TabSTD1_SigniLbx1(iterL);
% %         iterL=iterL+1;   
% %     else
% %         Lbx1_control(iter)=NaN;
% %         Lbx1_control_std(iter)=NaN;   
% %     end
% %     
% %     if ismember(i,TabParamInMemoryEn1)
% %         En1_control(iter)=TabMean1_SigniEn1(iterE);
% %         En1_control_std(iter)=TabSTD1_SigniEn1(iterE);
% %         iterE=iterE+1;   
% %     else
% %         En1_control(iter)=NaN;
% %         En1_control_std(iter)=NaN;   
% %     end
% %     
% %     if ismember(i,TabParamInMemorySim1)
% %         Sim1_control(iter)=TabMean1_SigniSim1(iterSi);
% %         Sim1_control_std(iter)=TabSTD1_SigniSim1(iterSi);
% %         iterSi=iterSi+1;   
% %     else
% %         Sim1_control(iter)=NaN;
% %         Sim1_control_std(iter)=NaN;   
% %     end
% %     
% %     if ismember(i,TabParamInMemoryShox2)
% %         Shox2_control(iter)=TabMean1_SigniShox2(iterSh);
% %         Shox2_control_std(iter)=TabSTD1_SigniShox2(iterSh);
% %         iterSh=iterSh+1;   
% %     else
% %         Shox2_control(iter)=NaN;
% %         Shox2_control_std(iter)=NaN;   
% %     end
% %     
% %         Master_control(iter)=mean(TabToScale(21:46,i));
% %         Master_control_std(iter)=std(TabToScale(21:46,i))/sqrt(25);
% %     
% %     iter=iter+1;
% %     
% %     
% % end
% % 
% %         
% % Impact = [Lbx1_control 0];
% % theta = linspace(0,2*pi,(length(Impact)));
% % polarplot(theta,Impact,'or');
% % hold on   
% % Impact = [Lbx1_control-Lbx1_control_std 0];
% % theta = linspace(0,2*pi,(length(Impact)));
% % polarplot(theta,Impact,'xr');
% % Impact = [Lbx1_control+Lbx1_control_std 0];
% % theta = linspace(0,2*pi,(length(Impact)));
% % polarplot(theta,Impact,'xr');
% % 
% % 
% % Impact = [En1_control 0];
% % theta = linspace(0,2*pi,(length(Impact)));
% % polarplot(theta,Impact,'ob');
% % Impact = [En1_control-En1_control_std 0];
% % theta = linspace(0,2*pi,(length(Impact)));
% % polarplot(theta,Impact,'xb');
% % Impact = [En1_control+En1_control_std 0];
% % theta = linspace(0,2*pi,(length(Impact)));
% % polarplot(theta,Impact,'xb');
% %  
% % Impact = [Sim1_control 0];
% % theta = linspace(0,2*pi,(length(Impact)));
% % polarplot(theta,Impact,'og');
% % Impact = [Sim1_control-Sim1_control_std 0];
% % theta = linspace(0,2*pi,(length(Impact)));
% % polarplot(theta,Impact,'xg');
% % Impact = [Sim1_control+Sim1_control_std 0];
% % theta = linspace(0,2*pi,(length(Impact)));
% % polarplot(theta,Impact,'xg');
% % 
% % Impact = [Shox2_control 0];
% % theta = linspace(0,2*pi,(length(Impact)));
% % polarplot(theta,Impact,'ok');
% % Impact = [Shox2_control-Shox2_control_std 0];
% % theta = linspace(0,2*pi,(length(Impact)));
% % polarplot(theta,Impact,'xk');
% % Impact = [Shox2_control+Shox2_control_std 0];
% % theta = linspace(0,2*pi,(length(Impact)));
% % polarplot(theta,Impact,'xk');
% % 
% % 
% % Impact = [Master_control 0];
% % theta = linspace(0,2*pi,(length(Impact)));
% % polarplot(theta,Impact,'-c');
% % Impact = [Master_control-Master_control_std 0];
% % theta = linspace(0,2*pi,(length(Impact)));
% % polarplot(theta,Impact,'--c');
% % Impact = [Master_control+Master_control_std 0];
% % theta = linspace(0,2*pi,(length(Impact)));
% % polarplot(theta,Impact,'--c');
% % 
% % 
% % 
% % ax=gca
% % ax.ThetaTick = [0:(360/length(C)):(360)]; 
% %     ax.ThetaTickLabel=TabSigniWhoAblationControl;
% %     
% % 
% %     
% % 
% %     
%%    
    %% Comp Subpopu of LBX1
%%    
    
    
    %Dyn Param
Param=[ 80:82 85:87 90:92   95:97 100:102 105:107 ... 
        110:112 115:117 120:122  125:127 130:132 135:137 ...
        140:142 145:147 150:152  155:157]; %16 * FAK


%
%vbls(i)=BookPCAv6.Properties.VariableNames(Param(I(i)));
clearvars A B C D E TabToScale TabToScalebe
clearvars TabMean1_SigniEn1 TabMean1_SigniLbx1 TabMean1_SigniSim1 TabMean1_SigniShox2
clearvars TabMean5_SigniEn1 TabMean5_SigniLbx1 TabMean5_SigniSim1 TabMean5_SigniShox2
clearvars TabMean5 TabSTD1 TabSTD1_SigniEn1 TabSTD1_SigniLbx1 TabSTD1_SigniShox2 TabSTD1_SigniSim1 TabSTD5 TabSTD5_SigniEn1 TabSTD5_SigniLbx1 TabSTD5_SigniShox2 TabSTD5_SigniSim1

for i=Param
    
    A=Data_Day1_Ptf1a(:,i);
    B=Data_Day1_Tlx3(:,i);
    C=Data_Day1_DMRT3(:,i);
    
    D=Data_Day1_Master(:,i);
    
    
    
    TabToScalebe=([A;B;C;D]);
    %TabToScale(:,i)=zscore(TabToScalebe);
    TabToScale(:,i)=(TabToScalebe-min(TabToScalebe)) / (max(TabToScalebe)-min(TabToScalebe));

    
end

%%
iter=1;
iter_signi=1;

clearvars TTEST_Abl TabMean1 TabMean5 TabSTD1 TabSTD5 TabDivide 
clearvars TabMean1_Signi TabMean5_Signi TabSTD1_Signi TabSTD5_Signi Tab_Signi_Impact TabSigniWho TabParamInMemory
clearvars TabParamInMemoryEn1 TabParamInMemoryLbx1 TabParamInMemoryShox2 TabParamInMemorySim1
clearvars Ptf1a_control Ptf1a_control_std
clearvars Tlx3_control Tlx3_control_std
clearvars DMRT3_control DMRT3_control_std
clearvars Master_control Master_control_std


clc

for i=Param
    

    TabMean1(iter)=mean(TabToScale(1:length(A),i));
    TabMean5(iter)=mean(TabToScale(14:39,i));
    
    TabSTD1(iter)=std(TabToScale(1:length(A),i))/sqrt(length(A));
    TabSTD5(iter)=std(TabToScale(14:39,i))/sqrt(length(14:39));

    TabDivide(iter)=(TabMean1(iter)/TabMean5(iter));
    

    [a b]=ttest2(TabToScale(1:length(A),i),TabToScale(14:39,i));
    TTEST_Abl(iter)=(b);
    
    if b<0.05
        
    TabMean1_SigniPtf1a(iter_signi)=TabMean1(iter);
    TabMean5_SigniPtf1a(iter_signi)=TabMean5(iter);
    TabSTD1_SigniPtf1a(iter_signi)=TabSTD1(iter);
    TabSTD5_SigniPtf1a(iter_signi)=TabSTD5(iter);
    Tab_Signi_ImpactPtf1a(iter_signi)=b;
    TabSigniWho(iter_signi)=BookPCAv7.Properties.VariableNames(([i]));
    TabParamInMemoryPtf1a(iter_signi)=i;
    
    iter_signi=iter_signi+1;
    
    end
    
    iter=iter+1;
    
    
end

iter=1;
iter_signi=1;

for i=Param
      
    interTlx3=length(A)+1:length(A)+length(B);
    
    TabMean1(iter)=mean(TabToScale(interTlx3,i));
    
    TabSTD1(iter)=std(TabToScale(interTlx3,i))/sqrt(length(interTlx3));    

    [a b]=ttest2(TabToScale(interTlx3,i),TabToScale(14:39,i));
    TTEST_Abl(iter)=(b);
    
    if b<0.05
        
    TabMean1_SigniTlx3(iter_signi)=TabMean1(iter);
    TabMean5_SigniTlx3(iter_signi)=TabMean5(iter);
    TabSTD1_SigniTlx3(iter_signi)=TabSTD1(iter);
    TabSTD5_SigniTlx3(iter_signi)=TabSTD5(iter);
    Tab_Signi_ImpactTlx3(iter_signi)=b;
    TabSigniWho(iter_signi)=BookPCAv7.Properties.VariableNames(([i]));
    TabParamInMemoryTlx3(iter_signi)=i;

    
    iter_signi=iter_signi+1;
    
    end
    
    iter=iter+1;
    
    
end

iter=1;
iter_signi=1;

for i=Param
    
    interDMRT3=1+length(A)+length(B) : length(A)+length(B)+length(C);
    
    TabMean1(iter)=mean(TabToScale(interDMRT3,i));
    
    TabSTD1(iter)=std(TabToScale(interDMRT3,i))/sqrt(length(interDMRT3));    

    [a b]=ttest2(TabToScale(interDMRT3,i),TabToScale(14:39,i));
    TTEST_Abl(iter)=(b);
    
    if b<0.05
        
    TabMean1_SigniDMRT3(iter_signi)=TabMean1(iter);
    TabMean5_SigniDMRT3(iter_signi)=TabMean5(iter);
    TabSTD1_SigniDMRT3(iter_signi)=TabSTD1(iter);
    TabSTD5_SigniDMRT3(iter_signi)=TabSTD5(iter);
    Tab_Signi_ImpactDMRT3(iter_signi)=b;
    TabSigniWho(iter_signi)=BookPCAv7.Properties.VariableNames(([i]));
    TabParamInMemoryDMRT3(iter_signi)=i;
    
    iter_signi=iter_signi+1;
    
    end
    
    iter=iter+1;
    
    
end

clearvars C TabParamInMemory
TabParamInMemory=[TabParamInMemoryPtf1a TabParamInMemoryTlx3 TabParamInMemoryDMRT3]
C=unique(TabParamInMemory)







% errorbar(TabMean1_Signi,TabSTD1_Signi)
% errorbar(TabMean5_Signi,TabSTD5_Signi)
% axis([-1 length(TabMean1_Signi)+1 0 1])


figure(2)

iter=1;
for i=C
    
TabSigniWhoAblationControl(iter)=BookPCAv7.Properties.VariableNames(([i]));
iter=iter+1;
    
end


iter=1;
iterP=1;
iterT=1;
iterD=1;

for i=C
    
    if ismember(i,TabParamInMemoryPtf1a)
        Ptf1a_control(iter)=TabMean1_SigniPtf1a(iterP);
        Ptf1a_control_std(iter)=TabSTD1_SigniPtf1a(iterP);
        iterP=iterP+1;   
    else
        Ptf1a_control(iter)=NaN;
        Ptf1a_control_std(iter)=NaN;   
    end
    
    if ismember(i,TabParamInMemoryTlx3)
        Tlx3_control(iter)=TabMean1_SigniTlx3(iterT);
        Tlx3_control_std(iter)=TabSTD1_SigniTlx3(iterT);
        iterT=iterT+1;   
    else
        Tlx3_control(iter)=NaN;
        Tlx3_control_std(iter)=NaN;   
    end
    
    if ismember(i,TabParamInMemoryDMRT3)
        DMRT3_control(iter)=TabMean1_SigniDMRT3(iterD);
        DMRT3_control_std(iter)=TabSTD1_SigniDMRT3(iterD);
        iterD=iterD+1;   
    else
        DMRT3_control(iter)=NaN;
        DMRT3_control_std(iter)=NaN;   
    end
    

    
        Master_control(iter)=mean(TabToScale(14:39,i));
        Master_control_std(iter)=std(TabToScale(14:39,i))/sqrt(25);
    
    iter=iter+1;
    
    
end

        
Impact = [Ptf1a_control 0];
theta = linspace(0,2*pi,(length(Impact)));
polarplot(theta,Impact,'or');
hold on   
Impact = [Ptf1a_control-Ptf1a_control_std 0];
theta = linspace(0,2*pi,(length(Impact)));
polarplot(theta,Impact,'xr');
Impact = [Ptf1a_control+Ptf1a_control_std 0];
theta = linspace(0,2*pi,(length(Impact)));
polarplot(theta,Impact,'xr');


Impact = [Tlx3_control 0];
theta = linspace(0,2*pi,(length(Impact)));
polarplot(theta,Impact,'ob');
Impact = [Tlx3_control-Tlx3_control_std 0];
theta = linspace(0,2*pi,(length(Impact)));
polarplot(theta,Impact,'xb');
Impact = [Tlx3_control+Tlx3_control_std 0];
theta = linspace(0,2*pi,(length(Impact)));
polarplot(theta,Impact,'xb');
 
Impact = [DMRT3_control 0];
theta = linspace(0,2*pi,(length(Impact)));
polarplot(theta,Impact,'og');
Impact = [DMRT3_control-DMRT3_control_std 0];
theta = linspace(0,2*pi,(length(Impact)));
polarplot(theta,Impact,'xg');
Impact = [DMRT3_control+DMRT3_control_std 0];
theta = linspace(0,2*pi,(length(Impact)));
polarplot(theta,Impact,'xg');


Impact = [Master_control 0];
theta = linspace(0,2*pi,(length(Impact)));
polarplot(theta,Impact,'-c');
Impact = [Master_control-Master_control_std 0];
theta = linspace(0,2*pi,(length(Impact)));
polarplot(theta,Impact,'--c');
Impact = [Master_control+Master_control_std 0];
theta = linspace(0,2*pi,(length(Impact)));
polarplot(theta,Impact,'--c');



ax=gca
ax.ThetaTick = [0:(360/length(C)):(360)]; 
    ax.ThetaTickLabel=TabSigniWhoAblationControl;
    
    
    
    
    
    
    %%
    
    
    
    
    
    
    
    
    
    









































clearvars CatDivi CatDivi_std
clearvars TTEST_Abl_1
clc

interval=0;
interval_iter=0;

for i=1:32
    
    interval= 3*(i-1)+1 : 3*(i);
    
% interval=(interval_iter+1):(interval_iter+Param_inter(i))
%     
%CatDivi(i)=(mean((TabDivide(interval))));
%CatDivi_std(i)=std((TabDivide(interval)))/sqrt(Param_inter(i)-1);

CatDivi_1(i)=(mean((TabMean1(interval))));
CatDivi_std_1(i)=std((TabMean1(interval)))/sqrt(3-1);
CatDivi_5(i)=(mean((TabMean5(interval))));
CatDivi_std_5(i)=std((TabMean5(interval)))/sqrt(3-1);

[a b]=ttest2(TabMean1(interval),TabMean5(interval));
TTEST_Abl_1(i)=(b);
    
    
% 
% interval_iter=interval(Param_inter(i));

end


bar(TTEST_Abl_1)
hold on

errorbar(CatDivi_1,CatDivi_std_1)
errorbar(CatDivi_5,CatDivi_std_5)



axis([0 33 0 1])


figure(2)
Impact = [TabMean1 0];
theta = linspace(0,2*pi,(length(Impact)));
polarplot(theta,Impact,'-b');
hold on

Impact = [TabMean1+TabSTD1 0];
theta = linspace(0,2*pi,(length(Impact)));
polarplot(theta,Impact,'--b');

Impact = [TabMean1-TabSTD1 0];
theta = linspace(0,2*pi,(length(Impact)));
polarplot(theta,Impact,'--b');


Impact = [TabMean5 0];
theta = linspace(0,2*pi,(length(Impact)));
polarplot(theta,Impact,'-r');

Impact = [TabMean5+TabSTD5 0];
theta = linspace(0,2*pi,(length(Impact)));
polarplot(theta,Impact,'--r');

Impact = [TabMean5-TabSTD5 0];
theta = linspace(0,2*pi,(length(Impact)));
polarplot(theta,Impact,'--r');


Impact = [0.5-TTEST_Abl/2 0];
theta = linspace(0,2*pi,(length(Impact)));
polarplot(theta,Impact,'-k');


ax=gca
ax.ThetaTick = [0:(360/length(TTEST_Abl)):(360)]; 
 %ax.ThetaTickLabel = {'F','A','K','F','A','K','F','A','K','F','A','K','F','A','K','F','A','K','F','A','K',' '}
 %ax.ThetaTickLabel = {'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27',' '}
% ax.ThetaTickLabel = {'X','Y','Z' , 'X','Y','Z' , 'X','Y','Z' , ... 
%                      'X','Y','Z' , 'X','Y','Z' , '3D Angles', ...
%                      'X','Y','Z' , 'X','Y','Z' , 'X','Y','Z' , ... 
%                      'X','Y','Z' , 'X','Y','Z' , '3D Angles', ...
%                      ' '}

ax.ThetaTickLabel = ...
   {'F','A','K','F','A','K','F','A','K', ...
    'F','A','K','F','A','K','F','A','K', ...
    'F','A','K','F','A','K','F','A','K', ...
    'F','A','K','F','A','K','F','A','K', ...
    'F','A','K','F','A','K','F','A','K', ...
    'F','A','K', ...
    'F','A','K','F','A','K','F','A','K', ...
    'F','A','K','F','A','K','F','A','K', ...
    'F','A','K','F','A','K','F','A','K', ...
    'F','A','K','F','A','K','F','A','K', ...
    'F','A','K','F','A','K','F','A','K', ...
    'F','A','K'}
    


rlim([0 0.8])


% %%
% figure(1)
% ImpactBad = [TTEST_Abl_1 0];
% theta = linspace(0,2*pi,(length(TTEST_Abl_1)));
% polarplot(theta,TTEST_Abl_1,'-b');
% 
% %rlim([0 2])
% % 
% ax=gca
% ax.ThetaTick = [0:(360/33):360]; 
% %ax.ThetaTickLabel = {'F','A','K','F','A','K','F','A','K','F','A','K','F','A','K','F','A','K','F','A','K',' '}
%  %ax.ThetaTickLabel = {'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27',' '}
% 
% %  
 
 
% for i=1:29
%     
%     interval=
% 
%     CatDivi(i)=(mean((TabDivide(interval))));
%     CatDivi_std(i)=std((TabDivide(interval)))/sqrt(5);
%     
%     
% 
% end



%
 %hold on 
    %bar(TabDivide)
    %plot(TabSTD5,'bx')
    
    %axis([0 length(Param) -1 1])
% errorbar(TabMean1,TabSTD1)
% hold on
% errorbar(TabMean5,TabSTD5)

% 1:15 Amp Spike
% 16:30 Speed
% 31:45  Acc
% 46:50 Max spike
% 51:55 Min spike
% 56:60 Endpoint 
% 61:75 Absement
% 76:78 Angle

%%

Param_CoordX_FAK=[8:10 26:28 41:43 56:58 65:67 80:82 95:97 110:112 140:142];
Param_CoordY_FAK=[13:15 31:33 46:48 59:61 70:72 85:87 100:102 115:117 145:147];
Param_CoordZ_FAK=[18:20 36:38 51:53 62:64 75:77 90:92 105:107 120:122 150:152];

clc
%Param_kin
Param=[8:22 26:40 41:55 56:79 80:94 95:109 110:124 125:129 130:134 135:139 140:154 23:25 155:157];

Param_inter_stat=[5 5 5  5 5 5  5 5 5  3 3 3  5 5 5  3]; %16
Param_inter_dyn=[5 5 5  5 5 5  5 5 5  5 5 5  5 5 5  3]; %16

Param_inter=[Param_inter_stat Param_inter_dyn];

%%
































%% Ablation only dynamic
%Param=[6 7 8:22 23:25 26:40 41:55 56:58 59:61 62:64 65:79];
Param=[80:94 95:109 110:124 125:129 130:134 135:139 140:154 155:157]; %
%Param=[6 7 8:22 23:25 26:40 41:55 56:58 59:61 62:64 65:79 80:94 95:109 110:124 125:129 130:134 135:139 140:154 155:157];

% 1:15 Amp Spike
% 16:30 Speed
% 31:45  Acc
% 46:50 Max spike
% 51:55 Min spike
% 56:60 Endpoint 
% 61:75 Absement
% 76:78 Angle

%DataPCA=[Data_Day1_Lbx1;Data_Day1_Ptf1a;Data_Day1_Tlx3;Data_Day1_DMRT3]%;Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2;Data_Day1_Yoked;Data_Day1_Master];%;Data_Day1_Yoked];Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2];
DataPCA=[Data_Day1_Master;Data_Day1_Ptf1a;Data_Day1_Tlx3;Data_Day1_DMRT3;Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2]%;Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2;Data_Day1_Yoked;Data_Day1_Master];%;Data_Day1_Yoked];Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2];


DataPCA=DataPCA(:,Param);

[n_norm p_norm] = size(DataPCA);
e = ones(n_norm,1);
DataStand=(DataPCA-e*min(DataPCA))./(e*max(DataPCA)-e*min(DataPCA));

% Perform PCA
[COEFF_Ablat, SCORE_Ablat, LATENT_Ablat, TSQUARED_Ablat, EXPLAINED_Ablat]=pca(DataStand);%,'Centered',true);%,'VariableWeights','variance');
RecoDataCentered_Ablat=SCORE_Ablat*COEFF_Ablat';

%Master
plot(SCORE_Ablat(1:24,1),SCORE_Ablat(1:24,2),'.k','MarkerSize',20)
hold on
%plotell([SCORE_Ablat(1:24,1) SCORE_Ablat(1:24,2)],'-r','ok');
center = [mean(SCORE_Ablat(1:24,1)), mean(SCORE_Ablat(1:24,2))];
center_Master=center;
stdev = [std(SCORE_Ablat(1:24,1)), std(SCORE_Ablat(1:24,2))];
stdev_Master=stdev;
axh = gca(); 
llc = center(:)-stdev(:);
wh = stdev(:)*2; 
rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1]); 
plot(center(1),center(2),'xk')

% % Yoked
% plot(SCORE_Ablat(25:48,1),SCORE_Ablat(25:48,2),'ok')
% plotell([SCORE_Ablat(25:48,1) SCORE_Ablat(25:48,2)],'--k','ok');
% center = [mean(SCORE_Ablat(25:48,1)), mean(SCORE_Ablat(25:48,2))];
% center_Yoked=center;
% stdev = [std(SCORE_Ablat(25:48,1)), std(SCORE_Ablat(25:48,2))];
% stdev_Yoked=stdev;
% axh = gca(); 
% llc = center(:)-stdev(:);
% wh = stdev(:)*2; 
% rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1]); 
% plot(center(1),center(2),'xk')


%Ptf1a
plot(SCORE_Ablat(25:29,1),SCORE_Ablat(25:29,2),'.b','MarkerSize',20)
%plotell([SCORE_Ablat(49:53,1) SCORE_Ablat(49:53,2)],'-b','ob');
center = [mean(SCORE_Ablat(25:29,1)), mean(SCORE_Ablat(25:29,2))];
center_Ptf1a=center;
stdev = [std(SCORE_Ablat(25:29,1)), std(SCORE_Ablat(25:29,2))];
stdev_Ptf1a=stdev;
axh = gca(); 
llc = center(:)-stdev(:);
wh = stdev(:)*2; 
rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1],'EdgeColor','b'); 
plot(center(1),center(2),'xb')


%Tlx3
plot(SCORE_Ablat(30:33,1),SCORE_Ablat(30:33,2),'.r','MarkerSize',20)
%plotell([SCORE_Ablat(54:77,1) SCORE_Ablat(54:77,2)],'-r','or');
center = [mean(SCORE_Ablat(30:33,1)), mean(SCORE_Ablat(30:33,2))];
center_Tlx3=center;
stdev = [std(SCORE_Ablat(30:33,1)), std(SCORE_Ablat(30:33,2))];
stdev_Tlx3=stdev;
axh = gca(); 
llc = center(:)-stdev(:);
wh = stdev(:)*2; 
rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1],'EdgeColor','r'); 
plot(center(1),center(2),'xr')


%DMRT3
plot(SCORE_Ablat(34:38,1),SCORE_Ablat(34:38,2),'.g','MarkerSize',20)
%plotell([SCORE_Ablat(58:62,1) SCORE_Ablat(58:62,2)],'-g','og');
center = [mean(SCORE_Ablat(34:38,1)), mean(SCORE_Ablat(34:38,2))];
center_DMRT3=center;
stdev = [std(SCORE_Ablat(34:38,1)), std(SCORE_Ablat(34:38,2))];
stdev_DMRT3=stdev;
axh = gca(); 
llc = center(:)-stdev(:);
wh = stdev(:)*2; 
rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1],'EdgeColor','g'); 
plot(center(1),center(2),'xg')

%SIM1
plot(SCORE_Ablat(39:42,1),SCORE_Ablat(39:42,2),'.y','MarkerSize',20)
%plotell([SCORE_Ablat(63:67,1) SCORE_Ablat(63:67,2)],'-y','oy');
center = [mean(SCORE_Ablat(39:42,1)), mean(SCORE_Ablat(39:42,2))];
center_Sim1=center;
stdev = [std(SCORE_Ablat(39:42,1)), std(SCORE_Ablat(39:42,2))];
stdev_Sim1=stdev;
axh = gca(); 
llc = center(:)-stdev(:);
wh = stdev(:)*2; 
rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1],'EdgeColor','y'); 
plot(center(1),center(2),'xy')

%EN1
plot(SCORE_Ablat(43:45,1),SCORE_Ablat(43:45,2),'.c','MarkerSize',20)
%plotell([SCORE_Ablat(68:72,1) SCORE_Ablat(68:72,2)],'-c','oc');
center = [mean(SCORE_Ablat(43:45,1)), mean(SCORE_Ablat(43:45,2))];
center_En1=center;
stdev = [std(SCORE_Ablat(43:45,1)), std(SCORE_Ablat(43:45,2))];
stdev_En1=stdev;
axh = gca(); 
llc = center(:)-stdev(:);
wh = stdev(:)*2; 
rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1],'EdgeColor','c'); 
plot(center(1),center(2),'xc')

%SHOX2
plot(SCORE_Ablat(46:48,1),SCORE_Ablat(46:48,2),'.m','MarkerSize',20)
%plotell([SCORE_Ablat(73:77,1) SCORE_Ablat(73:77,2)],'-m','om');
center = [mean(SCORE_Ablat(46:48,1)), mean(SCORE_Ablat(46:48,2))];
center_Shox2=center;
stdev = [std(SCORE_Ablat(46:48,1)), std(SCORE_Ablat(46:48,2))];
stdev_Shox2=stdev;
axh = gca(); 
llc = center(:)-stdev(:);
wh = stdev(:)*2; 
rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1],'EdgeColor','m'); 
plot(center(1),center(2),'xm')

 
D_yoked = pdist([center_Master;center_Yoked],'euclidean');
D_Ptf1a = pdist([center_Master;center_Ptf1a],'euclidean');
D_Tlx3 = pdist([center_Master;center_Tlx3],'euclidean');
D_DMRT3 = pdist([center_Master;center_DMRT3],'euclidean');
D_Sim1 = pdist([center_Master;center_Sim1],'euclidean');
D_Shox2 = pdist([center_Master;center_Shox2],'euclidean');
D_En1 = pdist([center_Master;center_En1],'euclidean');

%Centroid distance
%bar([D_yoked D_Ptf1a D_Tlx3 D_DMRT3 D_Sim1 D_Shox2 D_En1])



%% Individual distance to Master centroid


plot(SCORE_Ablat(25:48,1),SCORE_Ablat(25:48,2),'ok')
clearvars D_yoked_ind

for i=1:24
    
    X=[center_Master; SCORE_Ablat(i,1) SCORE_Ablat(i,2)]
    D_master_ind(i)=pdist(X,'euclidean');
    
    D_master_1(i)=abs(center_Master(1)-SCORE_Ablat(i,1));
    D_master_2(i)=abs(center_Master(2)-SCORE_Ablat(i,2));
    
end

for i=1:24
    
    X=[center_Master; SCORE_Ablat(24+i,1) SCORE_Ablat(24+i,2)]
    D_yoked_ind(i)=pdist(X,'euclidean');
    
    D_yoked_1(i)=abs(center_Master(1)-SCORE_Ablat(24+i,1));
    D_yoked_2(i)=abs(center_Master(2)-SCORE_Ablat(24+i,2));
    
end

for i=1:5
    
    X=[center_Master; SCORE_Ablat(48+i,1) SCORE_Ablat(48+i,2)]
    D_ptf1a_ind(i)=pdist(X,'euclidean');
    
    D_ptf1a_1(i)=abs(center_Master(1)-SCORE_Ablat(48+i,1));
    D_ptf1a_2(i)=abs(center_Master(2)-SCORE_Ablat(48+i,2));
    
end

for i=1:4
    
    X=[center_Master; SCORE_Ablat(53+i,1) SCORE_Ablat(53+i,2)]
    D_tlx3_ind(i)=pdist(X,'euclidean');
    
    D_tlx3_1(i)=abs(center_Master(1)-SCORE_Ablat(53+i,1));
    D_tlx3_2(i)=abs(center_Master(2)-SCORE_Ablat(53+i,2));
    
end

for i=1:5
    
    X=[center_Master; SCORE_Ablat(57+i,1) SCORE_Ablat(57+i,2)]
    D_dmrt3_ind(i)=pdist(X,'euclidean');
    
    D_dmrt3_1(i)=abs(center_Master(1)-SCORE_Ablat(57+i,1));
    D_dmrt3_2(i)=abs(center_Master(2)-SCORE_Ablat(57+i,2));
    
    
end

for i=1:5
    
    X=[center_Master; SCORE_Ablat(62+i,1) SCORE_Ablat(62+i,2)]
    D_sim1_ind(i)=pdist(X,'euclidean');
    
    D_sim1_1(i)=abs(center_Master(1)-SCORE_Ablat(62+i,1));
    D_sim1_2(i)=abs(center_Master(2)-SCORE_Ablat(62+i,2));
    
end

for i=1:5
    
    X=[center_Master; SCORE_Ablat(67+i,1) SCORE_Ablat(67+i,2)]
    D_en1_ind(i)=pdist(X,'euclidean');
    
    D_en1_1(i)=abs(center_Master(1)-SCORE_Ablat(67+i,1));
    D_en1_2(i)=abs(center_Master(2)-SCORE_Ablat(67+i,2));
    
end

for i=1:5
    
    X=[center_Master; SCORE_Ablat(72+i,1) SCORE_Ablat(72+i,2)]
    D_shox2_ind(i)=pdist(X,'euclidean');
    
    D_shox2_1(i)=abs(center_Master(1)-SCORE_Ablat(72+i,1));
    D_shox2_2(i)=abs(center_Master(2)-SCORE_Ablat(72+i,2));
    
end


%%

Cat=[zeros(24,1);ones(24,1);2*ones(5,1);3*ones(4,1);4*ones(5,1);5*ones(5,1);6*ones(5,1);7*ones(5,1)];

Cat_cell = num2cell(Cat)

scatterhist(SCORE_Ablat(1:77,1),SCORE_Ablat(1:77,2),'Group',Cat,'Kernel','on');
figure(gcf)
hold on

%%

h = scatterhist(SCORE_Ablat(1:77,1),SCORE_Ablat(1:77,2),'Group',Cat,'Kernel','on');
hold on;
clr = get(h(1),'colororder');
boxplot(h(2),SCORE_Ablat(1:77,1),Cat,'orientation','horizontal','color',clr);
boxplot(h(3),SCORE_Ablat(1:77,2),Cat,'orientation','horizontal','color',clr);
set(h(2:3),'XTickLabel','');
view(h(3),[270,90]);  % Rotate the Y plot
axis(h(1),'auto');  % Sync axes

center = [mean(SCORE_Ablat(1:24,1)), mean(SCORE_Ablat(1:24,2))];
stdev = [std(SCORE_Ablat(1:24,1)), std(SCORE_Ablat(1:24,2))];
llc = center(:)-stdev(:);
axh = gca(); 
wh = stdev(:)*2; 
rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1],'EdgeColor',clr(1,:)); 

center = [mean(SCORE_Ablat(25:48,1)), mean(SCORE_Ablat(25:48,2))];
stdev = [std(SCORE_Ablat(25:48,1)), std(SCORE_Ablat(25:48,2))];
llc = center(:)-stdev(:);
axh = gca(); 
wh = stdev(:)*2; 
rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1],'EdgeColor',clr(2,:)); 

center = [mean(SCORE_Ablat(49:53,1)), mean(SCORE_Ablat(49:53,2))];
stdev = [std(SCORE_Ablat(49:53,1)), std(SCORE_Ablat(49:53,2))];
llc = center(:)-stdev(:);
axh = gca(); 
wh = stdev(:)*2; 
rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1],'EdgeColor',clr(3,:)); 

center = [mean(SCORE_Ablat(54:57,1)), mean(SCORE_Ablat(54:57,2))];
stdev = [std(SCORE_Ablat(54:57,1)), std(SCORE_Ablat(54:57,2))];
llc = center(:)-stdev(:);
axh = gca(); 
wh = stdev(:)*2; 
rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1],'EdgeColor',clr(4,:)); 

center = [mean(SCORE_Ablat(58:62,1)), mean(SCORE_Ablat(58:62,2))];
stdev = [std(SCORE_Ablat(58:62,1)), std(SCORE_Ablat(58:62,2))];
llc = center(:)-stdev(:);
axh = gca(); 
wh = stdev(:)*2; 
rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1],'EdgeColor',clr(5,:)); 

center = [mean(SCORE_Ablat(63:67,1)), mean(SCORE_Ablat(63:67,2))];
stdev = [std(SCORE_Ablat(63:67,1)), std(SCORE_Ablat(63:67,2))];
llc = center(:)-stdev(:);
axh = gca(); 
wh = stdev(:)*2; 
rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1],'EdgeColor',clr(6,:)); 

center = [mean(SCORE_Ablat(68:72,1)), mean(SCORE_Ablat(68:72,2))];
stdev = [std(SCORE_Ablat(68:72,1)), std(SCORE_Ablat(68:72,2))];
llc = center(:)-stdev(:);
axh = gca(); 
wh = stdev(:)*2; 
rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1],'EdgeColor',clr(7,:)); 

center = [mean(SCORE_Ablat(72:77,1)), mean(SCORE_Ablat(72:77,2))];
stdev = [std(SCORE_Ablat(72:77,1)), std(SCORE_Ablat(72:77,2))];
llc = center(:)-stdev(:);
axh = gca(); 
wh = stdev(:)*2; 
rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1],'EdgeColor',clr(1,:)); 






hold off;


%% BOOTSTRAP

D_dmrt3_1_boot = bootstrp(24,@(x)[mean(x)],D_dmrt3_1);










%%

%Other
%plot3(SCORE_Ablat(20:34,1),SCORE_Ablat(20:34,2),SCORE_Ablat(20:34,3),'om')



% 
    %d2 = mahal(SCORE_Ablat(15:19,1:1),SCORE_Ablat(1:5,1:1))
% 
% 
% 
% 
% d2 = mahal(SCORE_Ablat(6:10,1:3),SCORE_Ablat(35:60,1:3))
% d2 = mahal(SCORE_Ablat(6:10,1:3),SCORE_Ablat(71:end,1:3))





xlabel('X')
ylabel('Y')
zlabel('Z')


% 
% clearvars AblatedScore
% 
% AblatedScore(1:5,1)=(SCORE_Ablat(1:5,1))
% AblatedScore(1:5,2)=(SCORE_Ablat(6:10,1))
% AblatedScore(1:4,3)=(SCORE_Ablat(11:14,1))
% AblatedScore(1:5,4)=(SCORE_Ablat(15:19,1))
% AblatedScore(1:15,5)=(SCORE_Ablat(20:34,1))


%% Indiv Abla + Master Yoked

Param=[80:94 95:109 110:124 125:129 130:134 135:139 140:154 155:157]; 
%Param=[95:109 110:124 125:129]; 

DataPCA=[Data_Day1_Tlx3;Data_Day1_Yoked;Data_Day1_Master];%;Data_Day1_Yoked];Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2];
DataPCA=DataPCA(:,Param);

[n_norm p_norm] = size(DataPCA);
e = ones(n_norm,1);
DataStand=(DataPCA-e*min(DataPCA))./(e*max(DataPCA)-e*min(DataPCA));

% Perform PCA
[COEFF_Ablat, SCORE_Ablat, LATENT_Ablat, TSQUARED_Ablat, EXPLAINED_Ablat]=pca(DataStand);%,'Centered',true);%,'VariableWeights','variance');
RecoDataCentered_Ablat=SCORE_Ablat*COEFF_Ablat';

%Lbx1
plot3(SCORE_Ablat(1:4,1),SCORE_Ablat(1:4,2),SCORE_Ablat(1:4,3),'om')
hold on

plot3(SCORE_Ablat(5:5+26,1),SCORE_Ablat(5:5+26,2),SCORE_Ablat(5:5+26,3),'ob')
plot3(SCORE_Ablat(5+26:5+51,1),SCORE_Ablat(5+26:5+51,2),SCORE_Ablat(5+26:5+51,3),'oc')







%% PCA ablation + all
DataPCA=[Data_Day1_Lbx1;Data_Day2_Lbx1;Data_Day3_Lbx1;Data_Day4_Lbx1;Data_Day5_Lbx1;Data_Day1_Master;Data_Day1_Yoked];%;Data_Day1_Sim1;Data_Day1_Shox2;Data_Day1_En1;Data_Day1_DMRT3;Data_Day1_Tlx3;Data_Day2_Tlx3;Data_Day3_Tlx3;Data_Day4_Tlx3;Data_Day5_Tlx3;Data_Day1_Ptf1a;Data_Day2_Ptf1a;Data_Day3_Ptf1a;Data_Day4_Ptf1a;Data_Day5_Ptf1a];
DataPCA=DataPCA(:,Param);

[n_norm p_norm] = size(DataPCA);
e = ones(n_norm,1);
DataStand=(DataPCA-e*min(DataPCA))./(e*max(DataPCA)-e*min(DataPCA));

% Perform PCA
[COEFF_Ablat, SCORE_Ablat, LATENT_Ablat, TSQUARED_Ablat, EXPLAINED_Ablat]=pca(DataStand);%,'Centered',true);%,'VariableWeights','variance');
RecoDataCentered_Ablat=SCORE_Ablat*COEFF_Ablat';

plot(SCORE_Ablat(1:5,1),SCORE_Ablat(1:5,2),'og')
hold on
plot(SCORE_Ablat(6:10,1),SCORE_Ablat(6:10,2),'ok')
plot(SCORE_Ablat(11:15,1),SCORE_Ablat(11:15,2),'oc')
plot(SCORE_Ablat(16:20,1),SCORE_Ablat(16:20,2),'om')
plot(SCORE_Ablat(21:25,1),SCORE_Ablat(21:25,2),'oy')


MasterFromData_All=[1 2 3 4 49 50 51 61 62 63 73 74 81 82 83 84 85 101 102 103 125 126 127 155 156 157]+25;
YokedFromData_All=[5 6 7 8 52 53 54 64 65 66 75 76 86 87 88 89 90 104 105 106 128 129 130 158 159 160]+25;

plot(SCORE_Ablat(MasterFromData_All,1),SCORE_Ablat(MasterFromData_All,2),'or')
plot(SCORE_Ablat(YokedFromData_All,1),SCORE_Ablat(YokedFromData_All,2),'ob')



%% Perfom PCA Ptf1a
DataPCA=[Data_Day1_Master;Data_Day1_Yoked;Data_Day1_Lbx1];
DataPCA=DataPCA(:,Param);

[n_norm p_norm] = size(DataPCA);
e = ones(n_norm,1);
DataStand=(DataPCA-e*min(DataPCA))./(e*max(DataPCA)-e*min(DataPCA));

% Perform PCA
[COEFF_Ablat, SCORE_Ablat, LATENT_Ablat, TSQUARED_Ablat, EXPLAINED_Ablat]=pca(DataStand);%,'Centered',true);%,'VariableWeights','variance');
RecoDataCentered_Ablat=SCORE_Ablat*COEFF_Ablat';

plot3(SCORE_Ablat(1:26,1),SCORE_Ablat(1:26,2),SCORE_Ablat(1:26,3),'or')
hold on
plot3(SCORE_Ablat(27:52,1),SCORE_Ablat(27:52,2),SCORE_Ablat(27:52,3),'ob')
plot3(SCORE_Ablat(53:end,1),SCORE_Ablat(53:end,2),SCORE_Ablat(53:end,3),'og')


%% Dynamic over 7 days
%Param=[80:94 95:109 110:124 125:129 130:134 135:139 140:154 155:157]; %
%Param=[80:94 110:124 125:129 130:134 135:139 140:154 155:157]; %



% Param=[80:82 85:87 90:92 95:97 100:102 105:107 110:112 115:117 120:122  ...
%       125:129 130:134 135:139 140:142 145:147 150:152 155:157]; 
%     


% DataPCA=[Data_Day1_Lbx1;Data_Day2_Lbx1;Data_Day3_Lbx1;Data_Day4_Lbx1;Data_Day5_Lbx1;Data_Day1_Master;Data_Day1_Yoked];
% DataPCA=[DataPCA;Data_Day1_Ptf1a;Data_Day2_Ptf1a;Data_Day3_Ptf1a;Data_Day4_Ptf1a;Data_Day5_Ptf1a];
% DataPCA=[DataPCA;Data_Day1_Tlx3;Data_Day2_Tlx3;Data_Day3_Tlx3;Data_Day4_Tlx3;Data_Day5_Tlx3];
% DataPCA=[DataPCA;Data_Day1_DMRT3;Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2]

% DataPCA=[Data_Day1_Ptf1a;Data_Day2_Ptf1a;Data_Day3_Ptf1a;Data_Day4_Ptf1a;Data_Day5_Ptf1a;Data_Day1_Master;Data_Day1_Yoked];
% DataPCA=[DataPCA;Data_Day1_Lbx1;Data_Day2_Lbx1;Data_Day3_Lbx1;Data_Day4_Lbx1;Data_Day5_Lbx1];
% DataPCA=[DataPCA;Data_Day1_Tlx3;Data_Day2_Tlx3;Data_Day3_Tlx3;Data_Day4_Tlx3;Data_Day5_Tlx3];
% DataPCA=[DataPCA;Data_Day1_DMRT3;Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2]

% DataPCA=[Data_Day1_Tlx3;Data_Day2_Tlx3;Data_Day3_Tlx3;Data_Day4_Tlx3;Data_Day5_Tlx3;Data_Day1_Master;Data_Day1_Yoked];
% DataPCA=[DataPCA;Data_Day1_Ptf1a;Data_Day2_Ptf1a;Data_Day3_Ptf1a;Data_Day4_Ptf1a;Data_Day5_Ptf1a];
% DataPCA=[DataPCA;Data_Day1_Lbx1;Data_Day2_Lbx1;Data_Day3_Lbx1;Data_Day4_Lbx1;Data_Day5_Lbx1];
% DataPCA=[DataPCA;Data_Day1_DMRT3;Data_Day1_Sim1;Data_Day1_En1;Data_Day1_Shox2]
% 

Param=[6 7 8:22 23:25 26:40 41:55 56:58 59:61 62:64 65:79]; %65:79 


DataPCA=[Data_Day1_Master;Data_Day1_Yoked;Data_Day1_Tlx3;];
DataPCA=[DataPCA;Data_Day1_Ptf1a];
DataPCA=[DataPCA;Data_Day1_DMRT3;Data_Day1_En1;Data_Day1_Shox2;Data_Day1_Sim1];

DataPCA=[DataPCA;Data_All];



DataPCA=DataPCA(:,Param);

[n_norm p_norm] = size(DataPCA);
e = ones(n_norm,1);
DataStand=(DataPCA-e*min(DataPCA))./(e*max(DataPCA)-e*min(DataPCA));

% Perform PCA
[COEFF_Ablat, SCORE_Ablat, LATENT_Ablat, TSQUARED_Ablat, EXPLAINED_Ablat]=pca(DataStand);%,'Centered',true);%,'VariableWeights','variance');
RecoDataCentered_Ablat=SCORE_Ablat*COEFF_Ablat';


% 
% plot3(SCORE_Ablat(1:5,1),SCORE_Ablat(1:5,2),SCORE_Ablat(1:5,3),'om')
% hold on
% plot3(SCORE_Ablat(6:10,1),SCORE_Ablat(6:10,2),SCORE_Ablat(6:10,3),'oy')
% plot3(SCORE_Ablat(11:15,1),SCORE_Ablat(11:15,2),SCORE_Ablat(11:15,3),'og')
% plot3(SCORE_Ablat(16:20,1),SCORE_Ablat(16:20,2),SCORE_Ablat(16:20,3),'ok')
% plot3(SCORE_Ablat(21:25,1),SCORE_Ablat(21:25,2),SCORE_Ablat(21:25,3),'oc')
% plot3(SCORE_Ablat(26:51,1),SCORE_Ablat(26:51,2),SCORE_Ablat(26:51,3),'or')
% plot3(SCORE_Ablat(52:77,1),SCORE_Ablat(52:77,2),SCORE_Ablat(52:77,3),'ob')


%plot3(SCORE_Ablat(1:26,1),SCORE_Ablat(1:26,2),SCORE_Ablat(1:26,3),'ok')
hold on
%plot3(SCORE_Ablat(27:52,1),SCORE_Ablat(27:52,2),SCORE_Ablat(27:52,3),'ok')

plotell([SCORE_Ablat(1:24,1) SCORE_Ablat(1:24,2)],'-k','ok');
plotell([SCORE_Ablat(25:48,1) SCORE_Ablat(25:48,2)],'-k','ok');


plot3(SCORE_Ablat(49:52,1),SCORE_Ablat(49:52,2),SCORE_Ablat(49:52,3),'ob')
plot3(SCORE_Ablat(53:57,1),SCORE_Ablat(53:57,2),SCORE_Ablat(53:57,3),'or')
plot3(SCORE_Ablat(58:62,1),SCORE_Ablat(58:62,2),SCORE_Ablat(58:62,3),'og')
plot3(SCORE_Ablat(63:67,1),SCORE_Ablat(63:67,2),SCORE_Ablat(63:67,3),'oc')
plot3(SCORE_Ablat(68:72,1),SCORE_Ablat(68:72,2),SCORE_Ablat(68:72,3),'om')
plot3(SCORE_Ablat(73:77,1),SCORE_Ablat(73:77,2),SCORE_Ablat(73:77,3),'oy')
clearvars AblatedScore

DMRT3TEST=[1.02162879709906
0.817177074918871
0.635612903521269
1.39543251197253
0.89349441838626];

rng(0,'twister');
%Create a vector of 1000 random values drawn from a normal distribution with a mean of 500 and a standard deviation of 5.

a = std(SCORE_Ablat(49:52,1));
b = mean(SCORE_Ablat(49:52,1));
y = a.*randn(24,1) + b;

a = std(DMRT3TEST);
b = mean(DMRT3TEST);
y_dmrt3_abla = a.*randn(24,1) + b;

SCORE_Ablat(53:56,1)
Tlx3_boot_abla = bootstrp(24,@(x)[mean(x)],SCORE_Ablat(49:52,1));
Ptf1a_boot_abla = bootstrp(24,@(x)[mean(x)+std(x)],SCORE_Ablat(53:57,1));
DMRT3_boot_abla = bootstrp(24,@(x)[mean(x)+std(x)],SCORE_Ablat(58:62,1));
En1_boot_abla = bootstrp(24,@(x)[mean(x)+std(x)],SCORE_Ablat(63:67,1));
Shox2_boot_abla = bootstrp(24,@(x)[mean(x)+std(x)],SCORE_Ablat(68:72,1));
Sim1_boot_abla = bootstrp(24,@(x)[mean(x)+std(x)],SCORE_Ablat(73:77,1));


% %FOR TLX3
% plot3(SCORE_Ablat(1:4,1),SCORE_Ablat(1:4,2),SCORE_Ablat(1:4,3),'or')
% hold on
% plot3(SCORE_Ablat(5:8,1),SCORE_Ablat(5:8,2),SCORE_Ablat(5:8,3),'ob')
% plot3(SCORE_Ablat(9:12,1),SCORE_Ablat(9:12,2),SCORE_Ablat(9:12,3),'og')
% plot3(SCORE_Ablat(13:16,1),SCORE_Ablat(13:16,2),SCORE_Ablat(13:16,3),'ok')
% plot3(SCORE_Ablat(17:20,1),SCORE_Ablat(17:20,2),SCORE_Ablat(17:20,3),'oc')
% 
% 
% clearvars AblatedScore
% 
% AblatedScore(1:4,1)=(SCORE_Ablat(1:4,1))
% AblatedScore(1:4,2)=(SCORE_Ablat(5:8,1))
% AblatedScore(1:4,3)=(SCORE_Ablat(9:12,1))
% AblatedScore(1:4,4)=(SCORE_Ablat(13:16,1))
% AblatedScore(1:4,5)=(SCORE_Ablat(17:20,1))




xlabel('X')
ylabel('Y')
zlabel('Z')

%%

    iterpolar=1;
for i=Param
    [a b]=ttest2([Data_Day1_Lbx1(:,i)],[Data_Day5_Lbx1(:,i)]);
    TTEST_Abl(iterpolar)=(b);
    iterpolar=iterpolar+1;   
end

errorbar(mean(Data_Day1_Lbx1(:,Param)),std(Data_Day1_Lbx1(:,Param))./sqrt(5))
hold on
errorbar(mean(Data_Day5_Lbx1(:,Param)),std(Data_Day5_Lbx1(:,Param))./sqrt(5))
%AAaaa=[mean(Data_Day1_Lbx1(:,i));s]

axis([0 55 -3 3])
%%
[p,t,stats] = anova1(MPG,Origin,'off');
[c,m,h,nms] = multcompare(stats);

%% Dynamic over 7 days - tlx3 and PTF1A
%Param=[80:94 95:109 110:124 125:129 130:134 135:139 140:154 155:157]; %
Param=[6 7 8:22 23:25 26:40 41:55 56:58 59:61 62:64 65:79 80:94 95:109 110:124 125:129 130:134 135:139 140:154 155:157]; %65:79 

% 1:15 Amp Spike
% 16:30 Speed
% 31:45  Acc
% 46:50 Max spike
% 51:55 Min spike
% 56:60 Endpoint 
% 61:75 Absement
% 76:78 Angle



%Param=[80:94 110:124 125:129 130:134 135:139 140:154 155:157]; %


%DataPCA=[Data_Day1_Master;Data_Day1_Yoked;Data_Day1_Ptf1a;Data_Day2_Ptf1a;Data_Day3_Ptf1a;Data_Day4_Ptf1a;Data_Day5_Ptf1a];
%DataPCA=[DataPCA;Data_Day1_Lbx1;Data_Day2_Lbx1;Data_Day3_Lbx1;Data_Day4_Lbx1;Data_Day5_Lbx1];
DataPCA=[Data_Day1_Master;Data_Day1_Yoked;Data_Day1_Ptf1a;Data_Day2_Ptf1a;Data_Day3_Ptf1a;Data_Day4_Ptf1a;Data_Day5_Ptf1a];
%DataPCA=[Data_Day1_Master;Data_Day1_Yoked;Data_Day1_Ptf1a;Data_Day2_Ptf1a;Data_Day3_Ptf1a;Data_Day4_Ptf1a;Data_Day5_Ptf1a];
DataPCA=[DataPCA;Data_Day1_DMRT3;Data_Day1_En1;Data_Day1_Shox2;Data_Day1_Sim1];
DataPCA=[DataPCA;Data_All];

DataPCA=DataPCA(:,Param);

[n_norm p_norm] = size(DataPCA);
e = ones(n_norm,1);
DataStand=(DataPCA-e*min(DataPCA))./(e*max(DataPCA)-e*min(DataPCA));

% Perform PCA
[COEFF_Ablat, SCORE_Ablat, LATENT_Ablat, TSQUARED_Ablat, EXPLAINED_Ablat]=pca(DataStand);%,'Centered',true);%,'VariableWeights','variance');
%RecoDataCentered_Ablat=SCORE_Ablat*COEFF_Ablat';


hold on

% plot3(SCORE_Ablat(1:25,1),SCORE_Ablat(1:25,2),SCORE_Ablat(1:25,3),'ob')
% plot3(SCORE_Ablat(1:5,1),SCORE_Ablat(1:5,2),SCORE_Ablat(1:5,3),'xb')
% plot3(SCORE_Ablat(21:25,1),SCORE_Ablat(21:25,2),SCORE_Ablat(21:25,3),'xm')

plot(SCORE_Ablat(1:26,1),SCORE_Ablat(1:26,2),'ok')
plot(SCORE_Ablat(27:52,1),SCORE_Ablat(27:52,2),'oy')

% plot(SCORE_Ablat(53:57,1),SCORE_Ablat(53:57,2),'ob')
% plot(SCORE_Ablat(58:62,1),SCORE_Ablat(58:62,2),'or')
% plot(SCORE_Ablat(63:67,1),SCORE_Ablat(63:67,2),'oc')
% plot(SCORE_Ablat(68:72,1),SCORE_Ablat(68:72,2),'om')
% plot(SCORE_Ablat(73:77,1),SCORE_Ablat(73:77,2),'og')

plot(SCORE_Ablat(53:56,1),SCORE_Ablat(53:56,2),'ob')
plot(SCORE_Ablat(57:60,1),SCORE_Ablat(57:60,2),'or')
plot(SCORE_Ablat(61:64,1),SCORE_Ablat(61:64,2),'oc')
plot(SCORE_Ablat(65:68,1),SCORE_Ablat(65:68,2),'om')
plot(SCORE_Ablat(69:72,1),SCORE_Ablat(69:72,2),'og')

xlabel('X')
ylabel('Y')
zlabel('Z')



%%



%% Dynamic/Static comparaison

%Param=[80:94 95:109 110:124 125:129 130:134 135:139 140:154 155:157]; 
Param=[6 7 8:22 23:25 26:40 41:55 56:58 59:61 62:64 65:79]; 


DataPCA=[Data_All;Data_Day1_Tlx3;Data_Day1_Ptf1a;Data_Day1_DMRT3;Data_Day1_En1;Data_Day1_Shox2;Data_Day1_Sim1];

%DataPCA=[DataPCA;Data_Day2_Lbx1;Data_Day3_Lbx1;Data_Day4_Lbx1;Data_Day5_Lbx1];
%DataPCA=[DataPCA;Data_Day2_Ptf1a;Data_Day3_Ptf1a;Data_Day4_Ptf1a;Data_Day5_Ptf1a];
%DataPCA=[DataPCA;Data_Day2_Tlx3;Data_Day3_Tlx3;Data_Day4_Tlx3;Data_Day5_Tlx3];
%DataPCA=[DataPCA;Data_All];


DataPCA=DataPCA(:,Param);

[n_norm p_norm] = size(DataPCA);
e = ones(n_norm,1);
DataStand=(DataPCA-e*min(DataPCA))./(e*max(DataPCA)-e*min(DataPCA));

% Perform PCA
[COEFF_Ablat, SCORE_Ablat, LATENT_Ablat, TSQUARED_Ablat, EXPLAINED_Ablat]=pca(DataStand);%,'Centered',true);%,'VariableWeights','variance');
RecoDataCentered_Ablat=SCORE_Ablat*COEFF_Ablat';

plot(SCORE_Ablat(WhichIsMaster,1),SCORE_Ablat(WhichIsMaster,2),'ok')
hold on
plot(SCORE_Ablat(WhichIsYoked,1),SCORE_Ablat(WhichIsYoked,2),'ok')

plot(SCORE_Ablat(53+138:56+138,1),SCORE_Ablat(53+138:56+138,2),'ob')
plot(SCORE_Ablat(57+138:61+138,1),SCORE_Ablat(57+138:61+138,2),'or')
plot(SCORE_Ablat(62+138:66+138,1),SCORE_Ablat(62+138:66+138,2),'og')
plot(SCORE_Ablat(67+138:71+138,1),SCORE_Ablat(67+138:71+138,2),'om')
plot(SCORE_Ablat(72+138:76+138,1),SCORE_Ablat(72+138:76+138,2),'oc')
plot(SCORE_Ablat(77+138:81+138,1),SCORE_Ablat(77+138:81+138,2),'oy')


axis([-1 1 -1 2])
grid on 


