function [Richness_omegaC,Richness_taylorC,Richness_taylorC_0,Omega_taylorC_terms] = RichnessEsts_ideal(TransectAbundance,k,Richness_raw)

%RichnessEsts_ideal.m
%Eden Tekwa Feb 8, 2022 - Apr 11, 2022
%function returns omegaC and omega_TC versions of omega and omega_T with survivorship bias eliminated for richness means based on
%the spatial TransectAbundance data: rows=transects, columns=species,
%values=individual counts

TransectAbundance=TransectAbundance(:,sum(TransectAbundance,1)>0); %take out empty species columns

%write observational process model
syms mn P
Dix=(1-exp(-mn/P))*P; %local detection probability of species i at sampled site x
Di=1-(1-Dix)^k; %detection probability of species i across all sampled sites x in community s

%second-order partial derivatives for 2nd order Taylor expansion of overall detection probability of any species in the community
%     d2Di_dmn2=diff(Di,'mn',2);
%     d2Di_dP2=diff(Di,'P',2);
%     d2Di_dmnP=diff(diff(Di,'mn'),'P');

%hard coded results from above for speed:
d2Di_dmn2 = - (k*exp(-mn/P)*(P*(exp(-mn/P) - 1) + 1)^(k - 1))/P - k*exp(-(2*mn)/P)*(P*(exp(-mn/P) - 1) + 1)^(k - 2)*(k - 1);
d2Di_dP2 = - k*(P*(exp(-mn/P) - 1) + 1)^(k - 2)*(k - 1)*(exp(-mn/P) + (mn*exp(-mn/P))/P - 1)^2 - (k*mn^2*exp(-mn/P)*(P*(exp(-mn/P) - 1) + 1)^(k - 1))/P^3;
d2Di_dmnP = (k*mn*exp(-mn/P)*(P*(exp(-mn/P) - 1) + 1)^(k - 1))/P^2 + k*exp(-mn/P)*(P*(exp(-mn/P) - 1) + 1)^(k - 2)*(k - 1)*(exp(-mn/P) + (mn*exp(-mn/P))/P - 1);

numTrans=size(TransectAbundance,1); %get number of transects
Richness_true=sum(sum(TransectAbundance,1)>0); %get raw richness
mn_true=mean(TransectAbundance,1);
P_true=zeros(1,Richness_true); %array to record occupancy for each species
for species=1:Richness_true
    P_true(species)=sum(TransectAbundance(:,species)>0)/numTrans; %occupancy as number of transects occupied divided by number of transects
end
P_true(P_true==0)=NaN;
mean_P_true=nanmean(P_true);
var_P_true=nanvar(P_true);
mean_mn_true=mean(mean(TransectAbundance(:,sum(TransectAbundance)>0)));
var_mn_true=var(mean(TransectAbundance(:,sum(TransectAbundance)>0)));
cov_mn_P_true=nancov(mean(TransectAbundance(:,sum(TransectAbundance)>0)),P_true);
if size(cov_mn_P_true,2)>1
    cov_mn_P_true=cov_mn_P_true(1,2);
else
    cov_mn_P_true=0;
end


%compute correction terms for omegaC_T
mn=mean_mn_true;
P=mean_P_true;
Omega_taylorC_terms=[eval(Di)
    eval(d2Di_dmn2)*var_mn_true/2
    eval(d2Di_dP2)*var_P_true/2
    eval(d2Di_dmnP)*cov_mn_P_true];

if sum(Omega_taylorC_terms)>0.1 %if sum of correction terms is greater than a threshold
    Ds_taylorC=sum(Omega_taylorC_terms); %use full correction
else
    Ds_taylorC=Omega_taylorC_terms(1); %else, use 0th order correction
end
if Ds_taylorC>1
    Ds_taylorC=1;
end
if Ds_taylorC<0.1
    Ds_taylorC=0.1;
end
Richness_taylorC=Richness_raw/Ds_taylorC;

%compute omegaC_0
Richness_taylorC_0=Richness_raw/Omega_taylorC_terms(1);

%compute omegaC
mn=mn_true;
P=P_true;
Ds_means=eval(Di);
Ds_mean=mean(Ds_means);
Richness_omegaC=Richness_raw/Ds_mean;
if Richness_raw==0
    Richness_omegaC=NaN;
    Richness_taylorC=NaN;
    Richness_taylorC_0=NaN;
end