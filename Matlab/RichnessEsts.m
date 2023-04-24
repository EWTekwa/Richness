function [Richness_raw,Chao1,GP,Chao2,ACE,JK_a,JK_i,Richness_omega,Richness_taylor,Richness_taylor_0,Omega_T_terms,meanStates] = RichnessEsts(TransectAbundance)

%RichnessEsts.m
%Eden Tekwa Feb 8, 2022 - Apr 11, 2022
%function returns richness estimates based on
%the spatial TransectAbundance data: rows=transects, columns=species,
%values=individual counts

TransectAbundance=TransectAbundance(:,sum(TransectAbundance,1)>0); %take out empty species columns

%write observational process model
syms mn k P Dixo
Dixo=(1-exp(-mn/P)); %obersvation probability in occupied patch x for species i
Dix=Dixo*P; %observation probability in patch x
Di=1-(1-Dix)^k; %detection probability of species i across all sampled sites x in community

%second-order partial derivatives for 2nd order Taylor expansion of overall detection probability of any species in the community
%     d2Di_dmn2=diff(Di,'mn',2);
%     d2Di_dP2=diff(Di,'P',2);
%     d2Di_dmnP=diff(diff(Di,'mn'),'P');

%hard coded results from above for speed:
d2Di_dmn2 = - (k*exp(-mn/P)*(P*(exp(-mn/P) - 1) + 1)^(k - 1))/P - k*exp(-(2*mn)/P)*(P*(exp(-mn/P) - 1) + 1)^(k - 2)*(k - 1);
d2Di_dP2 = - k*(P*(exp(-mn/P) - 1) + 1)^(k - 2)*(k - 1)*(exp(-mn/P) + (mn*exp(-mn/P))/P - 1)^2 - (k*mn^2*exp(-mn/P)*(P*(exp(-mn/P) - 1) + 1)^(k - 1))/P^3;
d2Di_dmnP = (k*mn*exp(-mn/P)*(P*(exp(-mn/P) - 1) + 1)^(k - 1))/P^2 + k*exp(-mn/P)*(P*(exp(-mn/P) - 1) + 1)^(k - 2)*(k - 1)*(exp(-mn/P) + (mn*exp(-mn/P))/P - 1);

numTrans=size(TransectAbundance,1); %get number of transects
Richness_raw=sum(sum(TransectAbundance,1)>0); %get raw richness
D=zeros(Richness_raw,1); %store observed species' observation probabilities
P_detected=zeros(1,Richness_raw); %array to record occupancy for each species
for species=1:Richness_raw

    P_detected(species)=sum(TransectAbundance(:,species)>0)/numTrans; %occupancy as number of transects occupied divided by number of transects
end
P_detected(P_detected==0)=NaN;
%compute on full dataset:
mean_P_detected=nanmean(P_detected);
var_P_detected=nanvar(P_detected);
mn_detected=mean(TransectAbundance,1);
if isnan(mn_detected)
    mn_detected=[];
end
mean_mn_detected=mean(mn_detected);
var_mn_detected=var(mn_detected);
cov_mn_P_detected=nancov(mn_detected,P_detected);
if size(cov_mn_P_detected,2)>1
    cov_mn_P_detected=cov_mn_P_detected(1,2);
else
    cov_mn_P_detected=0;
end

%compute Chao1 estimate
TransectAbundance=ceil(TransectAbundance); %convert data to discrete count
f1=sum(sum(TransectAbundance,1)==1); %number of singleton species
f2=sum(sum(TransectAbundance,1)==2); %number of doubleton species
Chao1=Richness_raw+f1*(f1-1)/(2*(f2+1)); %Chao1 richness estimator

%compute GP (Gamma-Poisson mixture) estimate (Chiu 2023, peerJ)
f3=sum(sum(TransectAbundance,1)==3); %number of tripleton species
if f3==0
    f3c=1;
else
    f3c=f3;
end
if f1==0
    f1c=1;
else
    f1c=f1;
end
A=2-(2*f2^2)/(3*f1c*f3c);
if f2>0
    f0Chao1=f1c^2/(2*f2);
else
    f0Chao1=f1c*(f1c-1)/2;
end
if A<1
    GP=Richness_raw+f0Chao1*max(.5,A);
else
    GP=Richness_raw+f0Chao1;
end


%compute Chao2 estimate
q1=sum(sum(TransectAbundance>0)==1); %number of species occurring in one transect only
q2=sum(sum(TransectAbundance>0)==2); %number of species occurring in two transect only
m=sum(sum(TransectAbundance>0));
Chao2=Richness_raw+((m-1)/m)*q1*(q1-1)/(2*(q2+1)); %Chao2 richness estimator

%compute abundance-based coverage estimator (ACE)
S_rare=sum(sum(TransectAbundance,1)<=10); %number of rare species (<=10 individuals)
S_abund=sum(sum(TransectAbundance,1)>10); %number of rare species (<=10 individuals)
n_rare=sum(TransectAbundance(:,sum(TransectAbundance,1)<=10),'all'); %total number of individuals in rare species
C_ACE=1-f1/n_rare; %sample coverage estimate
wsumfa=0;
for a=1:10
    wsumfa=wsumfa+a*(a-1)*sum(sum(TransectAbundance,1)==a);
end
V2=max(((S_rare/C_ACE)*wsumfa/(n_rare*(n_rare-1))-1),0); %coefficient of variation
if C_ACE>0
    ACE=S_abund+S_rare/C_ACE+(f1/C_ACE)*V2;
else
    ACE=Chao1;
end

%compute jackknife abundance estimator
JK_a=Richness_raw+2*f1-f2;

%compute jackknife incidence estimator
if m==1
    JK_i=Richness_raw;
else
    JK_i=Richness_raw+(q1*(2*m-3)/m-q2*((m-2)^2)/(m*(m-1)));
end

%compute correction terms for omega_T
k=numTrans;
mn=mean_mn_detected;
P=mean_P_detected;
Omega_T_terms=([eval(Di)
    eval(d2Di_dmn2)*var_mn_detected/2
    eval(d2Di_dP2)*var_P_detected/2
    eval(d2Di_dmnP)*cov_mn_P_detected]);

if sum(Omega_T_terms)>0.1 %if sum of correction terms is greater than a threshold
    D_omega=sum(Omega_T_terms); %use full correction
else
    D_omega=Omega_T_terms(1); %else, use 0th order correction
end
if D_omega>1
    D_omega=1;
end
if D_omega<0.1
    D_omega=0.1;
end
Richness_taylor=Richness_raw/D_omega;
Richness_taylor_0=Richness_raw/Omega_T_terms(1);
meanStates=[mean_mn_detected,mean_P_detected,var_mn_detected,var_P_detected,cov_mn_P_detected]';

%compute omega
mn=mn_detected;
P=P_detected;
D_means=eval(Di);
D_mean=mean(D_means);
Richness_omega=Richness_raw/D_mean;

if Richness_raw==0
    Richness_raw=NaN;
    Richness_taylor=NaN;
    Richness_taylor_0=NaN;
    Richness_omega=NaN;
    Chao1=NaN;
    GP=NaN;
    Chao2=NaN;
    ACE=NaN;
    JK_a=NaN;
    JK_i=NaN;
end