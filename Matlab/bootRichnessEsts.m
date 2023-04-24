%function [expectedRichness_raw,expectedChao1,expectedGP,expectedChao2,expectedACE,expectedJK_a,expectedJK_i,expectedRichness_omega,expectedRichness_taylor,expectedRichness_taylor_0] = bootRichnessEsts_all(TransectAbundance,numBoot)

function [Richness_raw,Chao1,GP,Chao2,ACE,JK_a,JK_i,Richness_omega,Richness_taylor,Richness_taylor_0,expectedRichness_raw,expectedChao1,expectedGP,expectedChao2,expectedACE,expectedJK_a,expectedJK_i,expectedRichness_omega,expectedRichness_taylor,expectedRichness_taylor_0,meanStates] = bootRichnessEsts(TransectAbundance,numBoot)
%bootRichnessEsts.m
%Eden Tekwa Feb 8, 2022, modified Mar 16, 2023
%function returns bootstrapped richness estimates based on
%the spatial TransectAbundance data: rows=transects, columns=species,
%values=individual counts

TransectAbundance=TransectAbundance(:,sum(TransectAbundance,1)>0); %take out empty species columns
[Richness_raw,Chao1,GP,Chao2,ACE,JK_a,JK_i,Richness_omega,Richness_taylor,Richness_taylor_0,~,meanStates] = RichnessEsts(TransectAbundance);

%store bootstrapped estimates for the 7 estimators:
expectedRichness_raw=zeros(numBoot,1); %raw
expectedChao1=zeros(numBoot,1); %Chao1
expectedGP=zeros(numBoot,1); %Gamma-Poisson
expectedChao2=zeros(numBoot,1); %Chao2
expectedACE=zeros(numBoot,1); %ACE
expectedJK_a=zeros(numBoot,1); %Jackknife (abundance)
expectedJK_i=zeros(numBoot,1); %Jackknife (incidence)
expectedRichness_omega=zeros(numBoot,1); %Omega
expectedRichness_taylor=zeros(numBoot,1); %Omega_taylor
expectedRichness_taylor_0=zeros(numBoot,1); %Omega_0

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


%hard coded results from above for speed:d2Di_dmn2 = - (k*exp(-mn/P)*(P*(exp(-mn/P) - 1) + 1)^(k - 1))/P - k*exp(-(2*mn)/P)*(P*(exp(-mn/P) - 1) + 1)^(k - 2)*(k - 1);
d2Di_dP2 = - k*(P*(exp(-mn/P) - 1) + 1)^(k - 2)*(k - 1)*(exp(-mn/P) + (mn*exp(-mn/P))/P - 1)^2 - (k*mn^2*exp(-mn/P)*(P*(exp(-mn/P) - 1) + 1)^(k - 1))/P^3;
d2Di_dmnP = (k*mn*exp(-mn/P)*(P*(exp(-mn/P) - 1) + 1)^(k - 1))/P^2 + k*exp(-mn/P)*(P*(exp(-mn/P) - 1) + 1)^(k - 2)*(k - 1)*(exp(-mn/P) + (mn*exp(-mn/P))/P - 1);

numTrans=size(TransectAbundance,1); %get number of transects

for resample=1:numBoot
    %resample transects with replacement first, for original transect number of times:
    sampleSet_init=TransectAbundance(randi(numTrans,numTrans,1),:);
    %then resample species with replacement, for "point-estimated richness
    %number of times" (pseudospecies unique for each estimator):
    sampleSetRaw=sampleSet_init(:,randi(Richness_raw,1,Richness_raw));
    sampleSetChao1=sampleSet_init(:,randi(Richness_raw,1,round(Chao1)));
    sampleSetGP=sampleSet_init(:,randi(Richness_raw,1,round(GP)));
    sampleSetChao2=sampleSet_init(:,randi(Richness_raw,1,round(Chao2)));
    sampleSetACE=sampleSet_init(:,randi(Richness_raw,1,round(ACE)));
    sampleSetJK_a=sampleSet_init(:,randi(Richness_raw,1,round(JK_a)));
    sampleSetJK_i=sampleSet_init(:,randi(Richness_raw,1,round(JK_i)));
    sampleSetOmega=sampleSet_init(:,randi(Richness_raw,1,round(Richness_omega)));
    sampleSetTaylor=sampleSet_init(:,randi(Richness_raw,1,round(Richness_taylor)));
    sampleSetTaylor0=sampleSet_init(:,randi(Richness_raw,1,round(Richness_taylor_0)));

    sampleSetRaw=sampleSetRaw(:,sum(sampleSetRaw,1)>0); %take out empty species columns
    sampleSetChao1=sampleSetChao1(:,sum(sampleSetChao1,1)>0); %take out empty species columns
    sampleSetGP=sampleSetGP(:,sum(sampleSetGP,1)>0); %take out empty species columns
    sampleSetChao2=sampleSetChao2(:,sum(sampleSetChao2,1)>0); %take out empty species columns
    sampleSetACE=sampleSetACE(:,sum(sampleSetACE,1)>0); %take out empty species columns
    sampleSetJK_a=sampleSetJK_a(:,sum(sampleSetJK_a,1)>0); %take out empty species columns
    sampleSetJK_i=sampleSetJK_i(:,sum(sampleSetJK_i,1)>0); %take out empty species columns
    sampleSetOmega=sampleSetOmega(:,sum(sampleSetOmega,1)>0); %take out empty species columns
    sampleSetTaylor=sampleSetTaylor(:,sum(sampleSetTaylor,1)>0); %take out empty species columns
    sampleSetTaylor0=sampleSetTaylor0(:,sum(sampleSetTaylor0,1)>0); %take out empty species columns
    
    %compute raw estimate
    Richness_raw_boot=sum(sum(sampleSetRaw,1)>0); %observed richness in raw boot
    
    %compute Chao1 estimate
    raw=sum(sum(sampleSetChao1,1)>0); %observed richness in Chao1 boot
    f1=sum(sum(sampleSetChao1,1)==1); %number of singleton species
    f2=sum(sum(sampleSetChao1,1)==2); %number of doubleton species
    Chao1_boot=raw+f1*(f1-1)/(2*(f2+1)); %Chao1 richness estimator
    
    %compute GP (Gamma-Poisson mixture) estimate (Chiu 2023, peerJ)
    raw=sum(sum(sampleSetGP,1)>0); %observed richness in GP boot
    f1=sum(sum(sampleSetGP,1)==1); %number of singleton species
    f2=sum(sum(sampleSetGP,1)==2); %number of doubleton species
    f3=sum(sum(sampleSetGP,1)==3); %number of tripleton species
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
        GP_boot=raw+f0Chao1*max(.5,A);
    else
        GP_boot=raw+f0Chao1;
    end
    
    %compute Chao2 estimate
    raw=sum(sum(sampleSetChao2,1)>0); %observed richness in Chao2 boot
    q1=sum(sum(sampleSetChao2>0)==1); %number of species occurring in one transect only
    q2=sum(sum(sampleSetChao2>0)==2); %number of species occurring in two transect only
    m=sum(sum(sampleSetChao2>0));
    Chao2_boot=raw+((m-1)/m)*q1*(q1-1)/(2*(q2+1)); %Chao2 richness estimator
    
    %compute abundance-based coverage estimator (ACE)
    raw=sum(sum(sampleSetACE,1)>0); %observed richness in ACE boot
    f1=sum(sum(sampleSetACE,1)==1); %number of singleton species
    f2=sum(sum(sampleSetACE,1)==2); %number of doubleton species
    S_rare=sum(sum(sampleSetACE,1)<=10); %number of rare species (<=10 individuals)
    S_abund=sum(sum(sampleSetACE,1)>10); %number of rare species (<=10 individuals)
    n_rare=sum(sampleSetACE(:,sum(sampleSetACE,1)<=10),'all'); %total number of individuals in rare species
    C_ACE=1-f1/n_rare; %sample coverage estimate
    wsumfa=0;
    for a=1:10
        wsumfa=wsumfa+a*(a-1)*sum(sum(sampleSetACE,1)==a);
    end
    V2=max(((S_rare/C_ACE)*wsumfa/(n_rare*(n_rare-1))-1),0); %coefficient of variation
    if C_ACE>0
        ACE_boot=S_abund+S_rare/C_ACE+(f1/C_ACE)*V2;
    else
        ACE_boot=raw+f1*(f1-1)/(2*(f2+1));
    end
    
    %compute jackknife abundance estimator
    raw=sum(sum(sampleSetJK_a,1)>0); %observed richness in jackknife abundance boot
    f1=sum(sum(sampleSetJK_a,1)==1); %number of singleton species
    f2=sum(sum(sampleSetJK_a,1)==2); %number of doubleton species
    JK_a_boot=raw+2*f1-f2;
    
    %compute jackknife incidence estimator
    raw=sum(sum(sampleSetJK_i,1)>0); %observed richness in jackknifte incidence boot
    q1=sum(sum(sampleSetJK_i>=0)==1); %number of species occurring in one transect only
    q2=sum(sum(sampleSetJK_i>=0)==2); %number of species occurring in two transect only
    m=sum(sum(sampleSetJK_i>=0));
    
    if m==1
        JK_i_boot=raw;
    else
        JK_i_boot=raw+(q1*(2*m-3)/m-q2*((m-2)^2)/(m*(m-1)));
    end
    
    %compute correction terms for omega_T
    raw=sum(sum(sampleSetTaylor,1)>0); %observed richness in omega_T boot
    P_detected=zeros(1,raw); %array to record occupancy for each species
    for species=1:raw
        P_detected(species)=sum(sampleSetTaylor(:,species)>0)/numTrans; %occupancy as number of transects occupied divided by number of transects
    end
    P_detected(P_detected==0)=NaN;
    %compute on full dataset:
    mean_P_detected=nanmean(P_detected);
    var_P_detected=nanvar(P_detected);
    mn_detected=mean(sampleSetTaylor(:,sum(sampleSetTaylor)>0));
    mean_mn_detected=mean(mn_detected);
    var_mn_detected=var(mean(sampleSetTaylor(:,sum(sampleSetTaylor)>0)));
    cov_mn_P_detected=nancov(mean(sampleSetTaylor(:,sum(sampleSetTaylor)>0)),P_detected);
    if size(cov_mn_P_detected,2)>1
        cov_mn_P_detected=cov_mn_P_detected(1,2);
    else
        cov_mn_P_detected=0;
    end
    
    k=numTrans;
    mn=mean_mn_detected;
    P=mean_P_detected;
    %add "eval" before "(" below if using ke in place of k
    Taylor_detectP_terms=([eval(Di)
        eval(d2Di_dmn2)*var_mn_detected/2
        eval(d2Di_dP2)*var_P_detected/2
        eval(d2Di_dmnP)*cov_mn_P_detected]);
    
    if sum(Taylor_detectP_terms)>0.1 %if sum of correction terms is positive and greater than a threshold
        Ds_omega=sum(Taylor_detectP_terms); %use full correction
    else
        Ds_omega=Taylor_detectP_terms(1); %else, use 0th order correction
    end
    if sum(Taylor_detectP_terms)>1
        Ds_omega=1;
    end
    Richness_taylor_boot=raw/Ds_omega;
    
    %omega_0:
    raw=sum(sum(sampleSetTaylor0,1)>0); %observed richness in omega_0 boot
    P_detected=zeros(1,raw); %array to record occupancy for each species
    for species=1:raw
        P_detected(species)=sum(sampleSetTaylor0(:,species)>0)/numTrans; %occupancy as number of transects occupied divided by number of transects
    end
    P_detected(P_detected==0)=NaN;
    %compute on full dataset:
    mean_P_detected=nanmean(P_detected);
    mn_detected=mean(sampleSetTaylor0(:,sum(sampleSetTaylor0)>0));
    mean_mn_detected=mean(mn_detected);

    k=numTrans;
    mn=mean_mn_detected;
    P=mean_P_detected;
    Richness_taylor_0_boot=raw/eval(Di);
    
    %omega:
    raw=sum(sum(sampleSetOmega,1)>0); %observed richness in omega boot
    P_detected=zeros(1,raw); %array to record occupancy for each species
    for species=1:raw
        P_detected(species)=sum(sampleSetOmega(:,species)>0)/numTrans; %occupancy as number of transects occupied divided by number of transects
    end
    P_detected(P_detected==0)=NaN;
    mn_detected=mean(sampleSetOmega(:,sum(sampleSetOmega)>0));
    mn=mn_detected;
    P=P_detected;
    Ds_means=eval(Di);
    Ds_mean=mean(Ds_means);
    Richness_omega_boot=raw/Ds_mean;
    
    expectedRichness_raw(resample)=Richness_raw_boot;
    expectedChao1(resample)=Chao1_boot;
    expectedGP(resample)=GP_boot;
    expectedChao2(resample)=Chao2_boot;
    expectedACE(resample)=ACE_boot;
    expectedJK_a(resample)=JK_a_boot;
    expectedJK_i(resample)=JK_i_boot;
    expectedRichness_omega(resample)=Richness_omega_boot;
    expectedRichness_taylor(resample)=Richness_taylor_boot;
    expectedRichness_taylor_0(resample)=Richness_taylor_0_boot;
end