%biodiversitySamplingSim.m
%Eden Tekwa Mar 22, 2023
%run simulations to evaluate all estimators

clear
set(0,'DefaultAxesFontSize',22)
scrsz = get(0,'ScreenSize');
figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/2.2 scrsz(4)/1.5]);

%define dataset
rng(50); %set random number generator seed
simReps=2000; %number of communities to simulation (2000)
numTests=2000; %number of random pairwise community comparisons (2000)
minSpecies=1; %minimum number of species
maxSpecies=100; %maximum number of species
k_min=2; %minimum number of sampled sites within region
k_max=10; %maximum number of sampled sites within region
k_mean=(k_max-k_min)/2+k_min;
num_sites=100; %true number of sites within region
maxRichDiff_all=[2 10 20]; %maximum richness difference between community pairs to be compared
numBoot=50; %number of bootstraps for each community estimate

%comment out one case at a time:
%case 1: main text 2 scenarios, a. low portion of inviduals observable, b. high spatial heterogeneitis
m_scenarios=[.4 1]; %m is the portion of individuals observable, or equivalently 1/m is the abundance needed for a species to be detected 63% of the time within a site
P_scenarios=[1 0.4]; %mean occupancy (variance=mean)
meanAbundance_scenarios=[2 2]; %mean abundances

% %case 2: 2 scenarios, a. mixed conditions, b. worse mixed conditions
% m_scenarios=[0.4 0.2]; %m is the portion of individuals observable, or equivalently 1/m is the abundance needed for a species to be detected 63% of the time within a site
% P_scenarios=[0.4 0.2]; %mean occupancy (variance=mean)
% meanAbundance_scenarios=[2 2]; %mean abundances

sdLogAbundance=1; %standard deviation of lognormal distribution for species abundance
%write observational process model
% syms mn k P Dixo
% Dixo=(1-exp(-mn/P)); %obersvation probability in occupied patch x for species i
% Dix=Dixo*P; %observation probability in patch x
% Di=1-(1-Dix)^k; %detection probability of species i across all sampled sites x in community
% 
% d2Dis_dmn2 = - (k*exp(-mn/P)*(P*(exp(-mn/P) - 1) + 1)^(k - 1))/P - k*exp(-(2*mn)/P)*(P*(exp(-mn/P) - 1) + 1)^(k - 2)*(k - 1);
% d2Dis_dP2 = - k*(P*(exp(-mn/P) - 1) + 1)^(k - 2)*(k - 1)*(exp(-mn/P) + (mn*exp(-mn/P))/P - 1)^2 - (k*mn^2*exp(-mn/P)*(P*(exp(-mn/P) - 1) + 1)^(k - 1))/P^3;
% d2Dis_dmnP = (k*mn*exp(-mn/P)*(P*(exp(-mn/P) - 1) + 1)^(k - 1))/P^2 + k*exp(-mn/P)*(P*(exp(-mn/P) - 1) + 1)^(k - 2)*(k - 1)*(exp(-mn/P) + (mn*exp(-mn/P))/P - 1);

%run simulations
ScenarioEvals=zeros(13,6,length(meanAbundance_scenarios)); %store 6 evaluation outcomes (columns) for all estimators (rows) and scenarios (3rd dimension)
ScenarioEvalsRel=zeros(13,6,length(meanAbundance_scenarios));
ScenarioEvalRankIndex=zeros(13,6,length(meanAbundance_scenarios));
ScenarioEvalRanks=zeros(13,7,length(meanAbundance_scenarios));
ScenarioEvalsRelDiff=zeros(13,7,length(meanAbundance_scenarios));
splot=1; %track subplot label
for scenario=1:length(meanAbundance_scenarios)
    meanAbundance=meanAbundance_scenarios(scenario);
    m_mean=m_scenarios(scenario);
    P_mean=P_scenarios(scenario);
    all_k=zeros(1,simReps);
    all_Species=zeros(1,simReps);
    all_True=zeros(1,simReps);
    all_omega=zeros(1,simReps);
    all_T=zeros(1,simReps);
    all_0=zeros(1,simReps);
    all_omegaC=zeros(1,simReps);
    all_TC=zeros(1,simReps);
    all_0C=zeros(1,simReps);
    all_Chao1=zeros(1,simReps);
    all_GP=zeros(1,simReps);
    all_Chao2=zeros(1,simReps);
    all_ACE=zeros(1,simReps);
    all_JK_a=zeros(1,simReps);
    all_JK_i=zeros(1,simReps);
    all_mean_n=zeros(1,simReps);
    all_var_n=zeros(1,simReps);
    all_mean_mn=zeros(1,simReps);
    all_var_mn=zeros(1,simReps);
    all_mean_P=zeros(1,simReps);
    all_var_P=zeros(1,simReps);
    all_cov_mn_P=zeros(1,simReps);
    all_mean_mn_ratio=zeros(1,simReps);
    all_var_mn_ratio=zeros(1,simReps);
    all_mean_P_ratio=zeros(1,simReps);
    all_var_P_ratio=zeros(1,simReps);
    all_true_n_estRelDiff=zeros(1,simReps);
    all_detected_n_estRelDiff=zeros(1,simReps);
    all_Omega_T_detectP_terms=zeros(4,simReps);
    all_Omega_TC_detectP_terms=zeros(4,simReps);
    all_bootStd_raw=zeros(1,simReps);
    all_bootStd_Chao1=zeros(1,simReps);
    all_bootStd_Omega=zeros(1,simReps);
    
    for i=1:simReps %(community i)
        m_community=m_mean+(rand-0.5)*(m_mean*(1-m_mean));
        k=randi([k_min k_max]); %mean k in community
        k(k==0)=1;
        P_community=P_mean+(rand-0.5)*(P_mean*(1-P_mean));
        Abundance_community=rand*2*meanAbundance;
    
        ni=random('lognormal',log(Abundance_community)-(sdLogAbundance^2)/2,sdLogAbundance,1,randi([minSpecies,maxSpecies]))/P_community; %generate random array of species densities
        %ni<1 will have inflated average abundance because at least one
        %individual has to be added to occupied patches (all species must
        %have average abundance of at least P). Thus the average abundance
        %is adjusted downward for all species:
        ni_totAdj=sum(ni(ni<1))-sum(ni<1);
        ni=ni*(sum(ni)+ni_totAdj)/sum(ni);
        
        nk_true=zeros(k,length(ni)); %store true abundance at sites k
        nk=zeros(k,length(ni)); %store observed abundance at sites k
        pi=P_community+(rand(1,length(ni))-0.5)*(P_community*(1-P_community));
        for site=1:num_sites %create true species local abundances at sites
             nk_true(site,:)=(max(random('poisson',ni,1,length(ni)),1)).*(rand(1,length(ni))<pi);
        end
        m_species=m_community+(rand(1,length(ni))-0.5)*(m_community*(1-m_community));
        rand_sites=randperm(num_sites,k); %pick k random site without replacement
        for site=1:k %sample from true sites
            rand_site=rand_sites(site);
            nk(site,:)=random('poisson',nk_true(rand_site,:).*m_species); %generate independent random observed abundance per sample site k
        end
        
        %record true abundance data
        n_s=mean(nk_true(:,sum(nk_true,1)>0));
        mn_true=(nk_true(:,sum(nk_true,1)>0)).*m_species(:,sum(nk_true,1)>0); %true observable abundance
        mn_s=mean(mn_true);
        num_species_true=sum(sum(nk_true)>0);
        sites_occupied_true=sum(nk_true(:,sum(nk_true,1)>0)>0);
        mean_n_true=mean(n_s);
        mean_mn_true=mean(mn_s);
        var_n_true=var(n_s);
        var_mn_true=var(mn_s);
        P_true=sites_occupied_true/num_sites;
        mean_P_true=mean(P_true);
        mean_1P_true=mean(1./P_true);
        var_P_true=var(P_true);
        var_1P_true=var(1./P_true);
        cov_mn_P_true=nancov(mn_s,P_true);
        if size(cov_mn_P_true,2)>1
            cov_mn_P_true=cov_mn_P_true(1,2);
        else
            cov_mn_P_true=0;
        end
        
        %record observation data
        mn_s_detected=mean(nk(:,sum(nk,1)>0),1);
        sites_occupied_detected=sum(nk(:,sum(nk,1)>0)>0);
        num_detected=sum(sum(nk,1)>0);
        mean_detected=mean(mn_s_detected);
        var_detected=var(mn_s_detected);
        P_detected=sites_occupied_detected/k;
        mean_P_detected=mean(P_detected);
        mean_1P_detected=mean(1./P_detected);
        var_P_detected=var(P_detected);
        var_1P_detected=var(1./P_detected);
        cov_mn_P_detected=nancov(mn_s_detected,P_detected);
        if size(cov_mn_P_detected,2)>1
            cov_mn_P_detected=cov_mn_P_detected(1,2);
        else
            cov_mn_P_detected=0;
        end
        
        [Richness_raw,Chao1,GP,Chao2,ACE,JK_a,JK_i,Richness_omega,Richness_T,Richness_0,Omega_T_detectP_terms,~] = RichnessEsts(nk); %compute estimates
        if Richness_raw==0 || isnan(Richness_raw) %if no species observed, rerun simulation replicate
            i=i-1;
            continue
        end
        [~,~,~,~,~,~,~,~,~,~,expectedRichness_raw,expectedChao1,~,~,~,~,~,expectedRichness_omega,~,~,~] = bootRichnessEsts(nk,numBoot);
        [Richness_omegaC,Richness_TC,Richness_0C,Omega_TC_detectP_terms] = RichnessEsts_ideal(mn_true,k,Richness_raw); %compute idealized approximation
        
        all_bootStd_raw(i)=nanstd(expectedRichness_raw);
        all_bootStd_Chao1(i)=nanstd(expectedChao1);
        all_bootStd_Omega(i)=nanstd(expectedRichness_omega);
        
        all_k(i)=k;
        all_mean_n(i)=mean_n_true;
        all_var_n(i)=var_n_true;
        all_mean_mn(i)=mean_mn_true;
        all_var_mn(i)=var_mn_true;
        all_mean_P(i)=mean_P_true;
        all_var_P(i)=var_P_true;
        all_cov_mn_P(i)=cov_mn_P_true;
        
        all_Omega_T_detectP_terms(:,i)=Omega_T_detectP_terms;
        all_Omega_TC_detectP_terms(:,i)=Omega_TC_detectP_terms;
        
        all_mean_mn_ratio(i)=mean_detected/mean_mn_true;
        all_var_mn_ratio(i)=var_detected/var_mn_true;
        all_mean_P_ratio(i)=mean_P_detected/mean_P_true;
        all_var_P_ratio(i)=var_P_detected/var_P_true;
        
        all_Species(i)=num_species_true;
        all_True(i)=num_detected;
        all_omega(i)=Richness_omega;
        all_T(i)=Richness_T;
        all_0(i)=Richness_0;
        all_omegaC(i)=Richness_omegaC;
        all_TC(i)=Richness_TC;
        all_0C(i)=Richness_0C;
        all_Chao1(i)=Chao1;
        all_GP(i)=GP;
        all_Chao2(i)=Chao2;
        all_ACE(i)=ACE;
        all_JK_a(i)=JK_a;
        all_JK_i(i)=JK_i;
        
        %%%%%%%%
        %diagnostic plots for current community:
        %
        %     diagfig=figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/3.5 scrsz(4)/4]);
        %     %set(diagfig,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
        %     subplot(1,2,1)
        %     hold on
        %     [mn_true_sorted,n_index]=sort(mean(mn_true(:,sum(nk_true,1)>0)),'descend');
        %     scatter([1:num_species_true],mn_true_sorted,'m','filled') %plot rank-true m/n curve
        %     scatter([1:num_detected],sort(mean(nk(:,sum(nk,1)>0),1),'descend'),'k','filled') %plot rank-observed m/n curve
        %     ax=gca;
        %     text(ax.XLim(2)*.1,ax.YLim(2)*.9,{['mean true=' num2str(mean_mn_true,2) ', observed=' num2str(mean_detected,2)],['var true=' num2str(var_mn_true,2) ', observed=' num2str(var_detected,2)]},'fontsize',16)
        %     xlabel 'species abundance rank'
        %     ylabel 'observable abundance (mn)'
        %     legend('true','observed','location','east')
        %
        %     subplot(1,2,2)
        %     hold on
        %     [P_sorted,P_index]=sort(P_true,'descend');
        %     [P_detected_sorted,P_sorted_index]=sort(P_detected,'descend');
        %     scatter([1:num_species_true],sites_occupied_true(P_index)./num_sites,'m','filled')
        %     scatter([1:num_detected],sites_occupied_detected(P_sorted_index)./k,'k','filled')
        %     ax=gca;
        %     text(ax.XLim(2)*.1,ax.YLim(2)*.90,{['mean true=' num2str(mean_P_true,2) ', observed=' num2str(mean_P_detected,2)],['var true=' num2str(var_P_true,2) ', observed=' num2str(var_P_detected,2)]},'fontsize',16)
        %     xlabel 'species occupancy rank'
        %     ylabel 'occupancy (P)'
        
    end
    
    %pairwise richness comparison (between simulation repicates)
    numRichDiff=length(maxRichDiff_all);
    numCorrect_detected=zeros(1,numRichDiff);
    numCorrect_Chao1=zeros(1,numRichDiff);
    numCorrect_GP=zeros(1,numRichDiff);
    numCorrect_Chao2=zeros(1,numRichDiff);
    numCorrect_ACE=zeros(1,numRichDiff);
    numCorrect_JK_a=zeros(1,numRichDiff);
    numCorrect_JK_i=zeros(1,numRichDiff);
    numCorrect_unbiased_omega=zeros(1,numRichDiff);
    numCorrect_unbiased_T=zeros(1,numRichDiff);
    numCorrect_unbiased_0=zeros(1,numRichDiff);
    numCorrect_unbiased_omegaC=zeros(1,numRichDiff);
    numCorrect_unbiased_TC=zeros(1,numRichDiff);
    numCorrect_unbiased_0C=zeros(1,numRichDiff);
    for test=1:numTests
        [~,temppos]=sort(all_Species);
        NonZeroPos=temppos(all_True(temppos)>0); %find indices of communities with detected species, ordered from low to high abundance
        CommunityOrder=randi(sum(all_True>0),1,2); %draw only from communities with detected species
        CommunityPair=NonZeroPos(CommunityOrder); %match picked community orders to indices
        rich1=all_Species(CommunityPair(1)); %richness of first random community
        for testvar=1:numRichDiff
            pos2=find(all_Species(NonZeroPos)<rich1+maxRichDiff_all(testvar) & all_Species(NonZeroPos)>rich1-maxRichDiff_all(testvar) & all_Species(NonZeroPos)~=rich1); %pick a second community with richness within the defined range when compared to the first community, but does not have identical richness
            if testvar==1 && isempty(pos2)
                test=test-1; %redo test set
                testvar=numRichDiff; %trigger exit from current test loop
            else
                CommunityPair(2)=NonZeroPos(pos2(randperm(length(pos2),1)));
                CommDiff_true=diff(all_Species(CommunityPair)); %difference (richness2-richness1)
                CommDiff_detected=diff(all_True(CommunityPair));
                CommDiff_Chao1=diff(all_Chao1(CommunityPair));
                CommDiff_GP=diff(all_GP(CommunityPair));
                CommDiff_Chao2=diff(all_Chao2(CommunityPair));
                CommDiff_ACE=diff(all_ACE(CommunityPair));
                CommDiff_JK_a=diff(all_JK_a(CommunityPair));
                CommDiff_JK_i=diff(all_JK_i(CommunityPair));
                CommDiff_unbiased_omega=diff(all_omega(CommunityPair));
                CommDiff_unbiased_T=diff(all_T(CommunityPair));
                CommDiff_unbiased_0=diff(all_0(CommunityPair));
                CommDiff_unbiased_omegaC=diff(all_omegaC(CommunityPair));
                CommDiff_unbiased_TC=diff(all_TC(CommunityPair));
                CommDiff_unbiased_0C=diff(all_0C(CommunityPair));
                numCorrect_detected(testvar)=numCorrect_detected(testvar)+((0^CommDiff_true)==(0^CommDiff_detected)); %add 1 if correct
                numCorrect_Chao1(testvar)=numCorrect_Chao1(testvar)+((0^CommDiff_true)==(0^CommDiff_Chao1)); %add 1 if correct
                numCorrect_GP(testvar)=numCorrect_GP(testvar)+((0^CommDiff_true)==(0^CommDiff_GP)); %add 1 if correct
                numCorrect_Chao2(testvar)=numCorrect_Chao2(testvar)+((0^CommDiff_true)==(0^CommDiff_Chao2)); %add 1 if correct
                numCorrect_ACE(testvar)=numCorrect_ACE(testvar)+((0^CommDiff_true)==(0^CommDiff_ACE)); %add 1 if correct
                numCorrect_JK_a(testvar)=numCorrect_JK_a(testvar)+((0^CommDiff_true)==(0^CommDiff_JK_a)); %add 1 if correct
                numCorrect_JK_i(testvar)=numCorrect_JK_i(testvar)+((0^CommDiff_true)==(0^CommDiff_JK_i)); %add 1 if correct
                numCorrect_unbiased_omega(testvar)=numCorrect_unbiased_omega(testvar)+((0^CommDiff_true)==(0^CommDiff_unbiased_omega)); %add 1 if correct
                numCorrect_unbiased_T(testvar)=numCorrect_unbiased_T(testvar)+((0^CommDiff_true)==(0^CommDiff_unbiased_T)); %add 1 if correct
                numCorrect_unbiased_0(testvar)=numCorrect_unbiased_0(testvar)+((0^CommDiff_true)==(0^CommDiff_unbiased_0)); %add 1 if correct
                numCorrect_unbiased_omegaC(testvar)=numCorrect_unbiased_omegaC(testvar)+((0^CommDiff_true)==(0^CommDiff_unbiased_omegaC)); %add 1 if correct
                numCorrect_unbiased_TC(testvar)=numCorrect_unbiased_TC(testvar)+((0^CommDiff_true)==(0^CommDiff_unbiased_TC)); %add 1 if correct
                numCorrect_unbiased_0C(testvar)=numCorrect_unbiased_0C(testvar)+((0^CommDiff_true)==(0^CommDiff_unbiased_0C)); %add 1 if correct
                
                
            end
        end
    end
    
    scenario
    mean_n=nanmean(all_mean_n)
    meamn=nanmean(all_mean_mn./all_mean_n)
    mean_P=nanmean(all_mean_P)
    
    %regress estimated richness to true richness:
    all_True(all_True==0)=NaN;
    mdl_detected=fitlm(all_Species,all_True);
    mdl_Chao1=fitlm(all_Species,all_Chao1);
    mdl_GP=fitlm(all_Species,all_GP);
    mdl_Chao2=fitlm(all_Species,all_Chao2);
    mdl_ACE=fitlm(all_Species,all_ACE);
    mdl_JK_a=fitlm(all_Species,all_JK_a);
    mdl_JK_i=fitlm(all_Species,all_JK_i);
    mdl_correction_Omega=fitlm(all_Species,all_omega);
    mdl_correction_T=fitlm(all_Species,all_T);
    mdl_correction_0=fitlm(all_Species,all_0);
    mdl_correction_OmegaC=fitlm(all_Species,all_omegaC);
    mdl_correction_TC=fitlm(all_Species,all_TC);
    mdl_correction_0C=fitlm(all_Species,all_0C);
    
    RSS_detected=nansum((all_Species-all_True).^2);
    RSS_Chao1=nansum((all_Species-all_Chao1).^2);
    RSS_GP=nansum((all_Species-all_GP).^2);
    RSS_Chao2=nansum((all_Species-all_Chao2).^2);
    RSS_ACE=nansum((all_Species-all_ACE).^2);
    RSS_JK_a=nansum((all_Species-all_JK_a).^2);
    RSS_JK_i=nansum((all_Species-all_JK_i).^2);
    RSS_correction_Omega=nansum((all_Species-all_omega).^2);
    RSS_correction_T=nansum((all_Species-all_T).^2);
    RSS_correction_0=nansum((all_Species-all_0).^2);
    RSS_correction_OmegaC=nansum((all_Species-all_omegaC).^2);
    RSS_correction_TC=nansum((all_Species-all_TC).^2);
    RSS_correction_0C=nansum((all_Species-all_0C).^2);
    TSS=nansum((all_Species-nanmean(all_Species)).^2);
    
    [~,coeffSEs_detected,coeffs_detected]=hac(mdl_detected);
    [~,coeffSEs_Chao1,coeffs_Chao1]=hac(mdl_Chao1);
    [~,coeffSEs_GP,coeffs_GP]=hac(mdl_GP);
    [~,coeffSEs_Chao2,coeffs_Chao2]=hac(mdl_Chao2);
    [~,coeffSEs_ACE,coeffs_ACE]=hac(mdl_ACE);
    [~,coeffSEs_JK_a,coeffs_JK_a]=hac(mdl_JK_a);
    [~,coeffSEs_JK_i,coeffs_JK_i]=hac(mdl_JK_i);
    [~,coeffSEs_correction_Omega,coeffs_correction_Omega]=hac(mdl_correction_Omega);
    [~,coeffSEs_correction_T,coeffs_correction_T]=hac(mdl_correction_T);
    [~,coeffSEs_correction_0,coeffs_correction_0]=hac(mdl_correction_0);
    [~,coeffSEs_correction_OmegaC,coeffs_correction_OmegaC]=hac(mdl_correction_OmegaC);
    [~,coeffSEs_correction_TC,coeffs_correction_TC]=hac(mdl_correction_TC);
    [~,coeffSEs_correction_0C,coeffs_correction_0C]=hac(mdl_correction_0C);
    
    slope_detected=coeffs_detected(2);
    slope_Chao1=coeffs_Chao1(2);
    slope_GP=coeffs_GP(2);
    slope_Chao2=coeffs_Chao2(2);
    slope_ACE=coeffs_ACE(2);
    slope_JK_a=coeffs_JK_a(2);
    slope_JK_i=coeffs_JK_i(2);
    slope_correction_Omega=coeffs_correction_Omega(2);
    slope_correction_T=coeffs_correction_T(2);
    slope_correction_0=coeffs_correction_0(2);
    slope_correction_OmegaC=coeffs_correction_OmegaC(2);
    slope_correction_TC=coeffs_correction_TC(2);
    slope_correction_0C=coeffs_correction_0C(2);
    slopeSD_detected=coeffSEs_detected(2)*sqrt(simReps);
    slopeSD_Chao1=coeffSEs_Chao1(2)*sqrt(simReps);
    slopeSD_GP=coeffSEs_GP(2)*sqrt(simReps);
    slopeSD_Chao2=coeffSEs_Chao2(2)*sqrt(simReps);
    slopeSD_ACE=coeffSEs_ACE(2)*sqrt(simReps);
    slopeSD_JK_a=coeffSEs_JK_a(2)*sqrt(simReps);
    slopeSD_JK_i=coeffSEs_JK_i(2)*sqrt(simReps);
    slopeSD_correction_Omega=coeffSEs_correction_Omega(2)*sqrt(simReps);
    slopeSD_correction_T=coeffSEs_correction_T(2)*sqrt(simReps);
    slopeSD_correction_0=coeffSEs_correction_0(2)*sqrt(simReps);
    slopeSD_correction_OmegaC=coeffSEs_correction_OmegaC(2)*sqrt(simReps);
    slopeSD_correction_TC=coeffSEs_correction_TC(2)*sqrt(simReps);
    slopeSD_correction_0C=coeffSEs_correction_0C(2)*sqrt(simReps);
    R2_detected=1-RSS_detected/TSS;
    R2_Chao1=1-RSS_Chao1/TSS;
    R2_GP=1-RSS_GP/TSS;
    R2_Chao2=1-RSS_Chao2/TSS;
    R2_ACE=1-RSS_ACE/TSS;
    R2_JK_a=1-RSS_JK_a/TSS;
    R2_JK_i=1-RSS_JK_i/TSS;
    R2_correction_Omega=1-RSS_correction_Omega/TSS;
    R2_correction_T=1-RSS_correction_T/TSS;
    R2_correction_0=1-RSS_correction_0/TSS;
    R2_correction_OmegaC=1-RSS_correction_OmegaC/TSS;
    R2_correction_TC=1-RSS_correction_TC/TSS;
    R2_correction_0C=1-RSS_correction_0C/TSS;
    
    %store all evaluation outcomes
    ScenarioEvals(1,:,scenario)=[slope_detected,slopeSD_detected,R2_detected,numCorrect_detected(1)/numTests,numCorrect_detected(2)/numTests,numCorrect_detected(3)/numTests];
    ScenarioEvals(2,:,scenario)=[slope_Chao1,slopeSD_Chao1,R2_Chao1,numCorrect_Chao1(1)/numTests,numCorrect_Chao1(2)/numTests,numCorrect_Chao1(3)/numTests];
    ScenarioEvals(3,:,scenario)=[slope_GP,slopeSD_GP,R2_GP,numCorrect_GP(1)/numTests,numCorrect_GP(2)/numTests,numCorrect_GP(3)/numTests];
    ScenarioEvals(4,:,scenario)=[slope_Chao2,slopeSD_Chao2,R2_Chao2,numCorrect_Chao2(1)/numTests,numCorrect_Chao2(2)/numTests,numCorrect_Chao2(3)/numTests];
    ScenarioEvals(5,:,scenario)=[slope_ACE,slopeSD_ACE,R2_ACE,numCorrect_ACE(1)/numTests,numCorrect_ACE(2)/numTests,numCorrect_ACE(3)/numTests];
    ScenarioEvals(6,:,scenario)=[slope_JK_a,slopeSD_JK_a,R2_JK_a,numCorrect_JK_a(1)/numTests,numCorrect_JK_a(2)/numTests,numCorrect_JK_a(3)/numTests];
    ScenarioEvals(7,:,scenario)=[slope_JK_i,slopeSD_JK_i,R2_JK_i,numCorrect_JK_i(1)/numTests,numCorrect_JK_i(2)/numTests,numCorrect_JK_i(3)/numTests];
    ScenarioEvals(8,:,scenario)=[slope_correction_Omega,slopeSD_correction_Omega,R2_correction_Omega,numCorrect_unbiased_omega(1)/numTests,numCorrect_unbiased_omega(2)/numTests,numCorrect_unbiased_omega(3)/numTests];
    ScenarioEvals(9,:,scenario)=[slope_correction_T,slopeSD_correction_T,R2_correction_T,numCorrect_unbiased_T(1)/numTests,numCorrect_unbiased_T(2)/numTests,numCorrect_unbiased_T(3)/numTests];
    ScenarioEvals(10,:,scenario)=[slope_correction_0,slopeSD_correction_0,R2_correction_0,numCorrect_unbiased_0(1)/numTests,numCorrect_unbiased_0(2)/numTests,numCorrect_unbiased_0(3)/numTests];
    ScenarioEvals(11,:,scenario)=[slope_correction_OmegaC,slopeSD_correction_OmegaC,R2_correction_OmegaC,numCorrect_unbiased_omegaC(1)/numTests,numCorrect_unbiased_omegaC(2)/numTests,numCorrect_unbiased_omegaC(3)/numTests];
    ScenarioEvals(12,:,scenario)=[slope_correction_TC,slopeSD_correction_TC,R2_correction_TC,numCorrect_unbiased_TC(1)/numTests,numCorrect_unbiased_TC(2)/numTests,numCorrect_unbiased_TC(3)/numTests];
    ScenarioEvals(13,:,scenario)=[slope_correction_0C,slopeSD_correction_0C,R2_correction_0C,numCorrect_unbiased_0C(1)/numTests,numCorrect_unbiased_0C(2)/numTests,numCorrect_unbiased_0C(3)/numTests];
    
    ScenarioEvalsRel(:,:,scenario)=ScenarioEvals(:,:,scenario);
    ScenarioEvalsRel(:,1,scenario)=abs(ScenarioEvals(:,1,scenario)-1);
    ScenarioEvalsRel(:,1:2,scenario)=1-ScenarioEvalsRel(:,1:2,scenario);
    [~,ScenarioEvalRankIndex(:,:,scenario)]=sort(ScenarioEvalsRel(:,:,scenario),1,'descend');
    [~,ScenarioEvalRanks(:,1:end-1,scenario)]=sort(ScenarioEvalRankIndex(:,:,scenario));
    ScenarioEvalRanks(:,end,scenario)=mean(ScenarioEvalRanks(:,1:end-1,scenario),2);
    
    ScenarioEvalsRelSD=std(ScenarioEvalsRel(:,:,scenario),0,1);
    ScenarioEvalsRelDiff(:,1:end-1,scenario)=(ScenarioEvalsRel(:,:,scenario)-mean(ScenarioEvalsRel(:,:,scenario),1))./ScenarioEvalsRelSD;
    ScenarioEvalsRelDiff(:,end,scenario)=mean(ScenarioEvalsRelDiff(:,1:end-1,scenario),2);
    
    subplot(2,2,scenario)
    title(['E[m]=' num2str(meamn,2) ', E[n]=' num2str(mean_n,2) ', E[P]=' num2str(mean_P,2) ', E[K]=' num2str(k_mean,1)],'fontweight','normal')
    hold on
    xvals=[1:100];
    window=10;
    sdBounds_obs=NaN(2,100);
    sdBounds_Chao1=NaN(2,100);
    sdBounds_Omega=NaN(2,100);
    sdBoot_obs=NaN(1,100);
    sdBoot_Chao1=NaN(1,100);
    sdBoot_Omega=NaN(1,100);
    for richCtr=1:100
        %get standard deviation bounds from estimates across communities
        indices=find(all_Species>=richCtr-(window-1)/2 & all_Species<=richCtr+(window-1)/2);
        meanEst=nanmean(all_True(indices));
        stdEst=nanstd(all_True(indices));
        sdBounds_obs(:,richCtr)=[meanEst,stdEst];
        meanEst=nanmean(all_Chao1(indices));
        stdEst=nanstd(all_Chao1(indices));
        sdBounds_Chao1(:,richCtr)=[meanEst,stdEst];
        meanEst=nanmean(all_omega(indices));
        stdEst=nanstd(all_omega(indices));
        sdBounds_Omega(:,richCtr)=[meanEst,stdEst];
        
        %get mean standard deviations from bootstrapped estimates within
        %communities
        sdBoot_obs(richCtr)=nanmean(all_bootStd_raw(indices));
        sdBoot_Chao1(richCtr)=nanmean(all_bootStd_Chao1(indices));
        sdBoot_Omega(richCtr)=nanmean(all_bootStd_Omega(indices));
    end
    plot(xvals,sdBounds_Chao1(1,:)+sdBounds_Chao1(2,:),'b','LineWidth',2);
    plot(xvals,sdBounds_Chao1(1,:)-sdBounds_Chao1(2,:),'b','LineWidth',2);
    plot(xvals,sdBounds_Omega(1,:)+sdBounds_Omega(2,:),'r','LineWidth',2);
    plot(xvals,sdBounds_Omega(1,:)-sdBounds_Omega(2,:),'r','LineWidth',2);
    plot(xvals,sdBounds_Chao1(1,:)+sdBoot_Chao1,':b','LineWidth',2);
    plot(xvals,sdBounds_Chao1(1,:)-sdBoot_Chao1,':b','LineWidth',2);
    plot(xvals,sdBounds_Omega(1,:)+sdBoot_Omega,':r','LineWidth',2);
    plot(xvals,sdBounds_Omega(1,:)-sdBoot_Omega,':r','LineWidth',2);
    
    xlabel 'true species richness'
    ylabel 'estimated species richness'
    ylims=ylim;
    %     text(minSpecies+5*(maxSpecies-minSpecies)/100,ylims(2)*0.95,['Raw slope=' num2str(slope_detected,2) '\pm' num2str(slopeSD_detected,2) ', R^2*=' num2str(R2_detected,2) ', pairwise_{' num2str(maxRichDiff_all(1)) ',' num2str(maxRichDiff_all(2)) ',' num2str(maxRichDiff_all(3)) '}=' num2str(numCorrect_detected(1)/numTests,2) ',' num2str(numCorrect_detected(2)/numTests,2) ',' num2str(numCorrect_detected(3)/numTests,2)])
    %     text(minSpecies+5*(maxSpecies-minSpecies)/100,ylims(2)*0.9,['Chao1 slope=' num2str(slope_Chao1,2) '\pm' num2str(slopeSD_Chao1,2) ', R^2*=' num2str(R2_Chao1,2) ', pairwise_{' num2str(maxRichDiff_all(1)) ',' num2str(maxRichDiff_all(2)) ',' num2str(maxRichDiff_all(3)) '}=' num2str(numCorrect_Chao1(1)/numTests,2) ',' num2str(numCorrect_Chao1(2)/numTests,2) ',' num2str(numCorrect_Chao1(3)/numTests,2)],'Color','b')
    %     text(minSpecies+5*(maxSpecies-minSpecies)/100,ylims(2)*0.85,['GP slope=' num2str(slope_GP,2) '\pm' num2str(slopeSD_GP,2) ', R^2*=' num2str(R2_GP,2) ', pairwise_{' num2str(maxRichDiff_all(1)) ',' num2str(maxRichDiff_all(2)) ',' num2str(maxRichDiff_all(3)) '}=' num2str(numCorrect_GP(1)/numTests,2) ',' num2str(numCorrect_GP(2)/numTests,2) ',' num2str(numCorrect_GP(3)/numTests,2)],'Color',[0.5 0.5 0.5])
    %     text(minSpecies+5*(maxSpecies-minSpecies)/100,ylims(2)*0.8,['Chao2 slope=' num2str(slope_Chao2,2) '\pm' num2str(slopeSD_Chao2,2) ', R^2*=' num2str(R2_Chao2,2) ', pairwise_{' num2str(maxRichDiff_all(1)) ',' num2str(maxRichDiff_all(2)) ',' num2str(maxRichDiff_all(3)) '}=' num2str(numCorrect_Chao2(1)/numTests,2) ',' num2str(numCorrect_Chao2(2)/numTests,2) ',' num2str(numCorrect_Chao2(3)/numTests,2)],'Color',[0.5 0.5 0.5])
    %     text(minSpecies+5*(maxSpecies-minSpecies)/100,ylims(2)*0.75,['ACE slope=' num2str(slope_ACE,2) '\pm' num2str(slopeSD_ACE,2) ', R^2*=' num2str(R2_ACE,2) ', pairwise_{' num2str(maxRichDiff_all(1)) ',' num2str(maxRichDiff_all(2)) ',' num2str(maxRichDiff_all(3)) '}=' num2str(numCorrect_ACE(1)/numTests,2) ',' num2str(numCorrect_ACE(2)/numTests,2)  ',' num2str(numCorrect_ACE(3)/numTests,2)],'Color',[0.5 0.5 0.5])
    %     text(minSpecies+5*(maxSpecies-minSpecies)/100,ylims(2)*0.7,['JK_a slope=' num2str(slope_JK_a,2) '\pm' num2str(slopeSD_JK_a,2) ', R^2*=' num2str(R2_JK_a,2) ', pairwise_{' num2str(maxRichDiff_all(1)) ',' num2str(maxRichDiff_all(2)) ',' num2str(maxRichDiff_all(3)) '}=' num2str(numCorrect_JK_a(1)/numTests,2) ',' num2str(numCorrect_JK_a(2)/numTests,2) ',' num2str(numCorrect_JK_a(3)/numTests,2)],'Color',[0.5 0.5 0.5])
    %     text(minSpecies+5*(maxSpecies-minSpecies)/100,ylims(2)*0.65,['JK_i slope=' num2str(slope_JK_i,2) '\pm' num2str(slopeSD_JK_i,2) ', R^2*=' num2str(R2_JK_i,2) ', pairwise_{' num2str(maxRichDiff_all(1)) ',' num2str(maxRichDiff_all(2)) ',' num2str(maxRichDiff_all(3)) '}=' num2str(numCorrect_JK_i(1)/numTests,2) ',' num2str(numCorrect_JK_i(2)/numTests,2) ',' num2str(numCorrect_JK_i(3)/numTests,2)],'Color',[0.5 0.5 0.5])
    %     text(minSpecies+5*(maxSpecies-minSpecies)/100,ylims(2)*0.6,['\Omega slope=' num2str(slope_correction_Omega,2) '\pm' num2str(slopeSD_correction_Omega,2) ', R^2*=' num2str(R2_correction_Omega,2) ', pairwise_{' num2str(maxRichDiff_all(1)) ',' num2str(maxRichDiff_all(2)) ',' num2str(maxRichDiff_all(3)) '}=' num2str(numCorrect_unbiased_omega(1)/numTests,2) ',' num2str(numCorrect_unbiased_omega(2)/numTests,2) ',' num2str(numCorrect_unbiased_omega(3)/numTests,2)],'Color','r')
    %     text(minSpecies+5*(maxSpecies-minSpecies)/100,ylims(2)*0.55,['\Omega_T slope=' num2str(slope_correction_T,2) '\pm' num2str(slopeSD_correction_T,2) ', R^2*=' num2str(R2_correction_T,2) ', pairwise_{' num2str(maxRichDiff_all(1)) ',' num2str(maxRichDiff_all(2)) ',' num2str(maxRichDiff_all(3)) '}=' num2str(numCorrect_unbiased_T(1)/numTests,2) ',' num2str(numCorrect_unbiased_T(2)/numTests,2) ',' num2str(numCorrect_unbiased_T(3)/numTests,2)],'Color',[0.5 0.5 0.5])
    %     text(minSpecies+5*(maxSpecies-minSpecies)/100,ylims(2)*0.5,['\Omega_0 slope=' num2str(slope_correction_0,2) '\pm' num2str(slopeSD_correction_0,2) ', R^2*=' num2str(R2_correction_0,2) ', pairwise_{' num2str(maxRichDiff_all(1)) ',' num2str(maxRichDiff_all(2)) ',' num2str(maxRichDiff_all(3)) '}=' num2str(numCorrect_unbiased_0(1)/numTests,2) ',' num2str(numCorrect_unbiased_0(2)/numTests,2) ',' num2str(numCorrect_unbiased_0(3)/numTests,2)],'Color',[0.5 0.5 0.5])
    %     text(minSpecies+5*(maxSpecies-minSpecies)/100,ylims(2)*0.45,['\Omega_C slope=' num2str(slope_correction_OmegaC,2) '\pm' num2str(slopeSD_correction_OmegaC,2) ', R^2*=' num2str(R2_correction_OmegaC,2) ', pairwise_{' num2str(maxRichDiff_all(1)) ',' num2str(maxRichDiff_all(2)) ',' num2str(maxRichDiff_all(3)) '}=' num2str(numCorrect_unbiased_omegaC(1)/numTests,2) ',' num2str(numCorrect_unbiased_omegaC(2)/numTests,2) ',' num2str(numCorrect_unbiased_omegaC(3)/numTests,2)],'Color',[0.5 0.5 0.5])
    %     text(minSpecies+5*(maxSpecies-minSpecies)/100,ylims(2)*0.4,['\Omega_{TC} slope=' num2str(slope_correction_TC,2) '\pm' num2str(slopeSD_correction_TC,2) ', R^2*=' num2str(R2_correction_TC,2) ', pairwise_{' num2str(maxRichDiff_all(1)) ',' num2str(maxRichDiff_all(2)) ',' num2str(maxRichDiff_all(3)) '}=' num2str(numCorrect_unbiased_TC(1)/numTests,2) ',' num2str(numCorrect_unbiased_TC(2)/numTests,2) ',' num2str(numCorrect_unbiased_TC(3)/numTests,2)],'Color',[0.5 0.5 0.5])
    %     text(minSpecies+5*(maxSpecies-minSpecies)/100,ylims(2)*0.35,['\Omega_{0C} slope=' num2str(slope_correction_0C,2) '\pm' num2str(slopeSD_correction_0C,2) ', R^2*=' num2str(R2_correction_0C,2) ', pairwise_{' num2str(maxRichDiff_all(1)) ',' num2str(maxRichDiff_all(2)) ',' num2str(maxRichDiff_all(3)) '}=' num2str(numCorrect_unbiased_0C(1)/numTests,2) ',' num2str(numCorrect_unbiased_0C(2)/numTests,2) ',' num2str(numCorrect_unbiased_0C(3)/numTests,2)],'Color',[0.5 0.5 0.5])
    %
    lgd=legend('observed','Chao1','\Omega');
    lgd.AutoUpdate='off';
    xlims=xlim;
    text(xlims(1)-diff(xlims)*0.18,ylims(2)+diff(ylims)*0.12,char(64+splot),'Fontsize',20)
    splot=splot+2;
    
%     %diagnostic plot for species abundance distributions in scenario:
%     figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/3.5 scrsz(4)/6]);
%     subplot(1,2,1)
%     hold on
%     histogram(all_mean_n,80,'facecolor','w');
%     xlim([0,8])
%     xlabel 'mn across communities'
%     ylabel 'count'
%     
%     subplot(1,2,2)
%     hold on
%     histogram(mean(nk(:,sum(nk,1)>0),1),'facecolor','w');
%     xlabel 'mn within a community'
%     ylabel 'count'
end

set(0,'DefaultAxesFontSize',14)
scrsz = get(0,'ScreenSize');
figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/2.2 scrsz(4)/4]);

for sc=1:scenario
    subplot(1,2,sc)
    estOrder=[8 9 10 11 12 13 2 3 4 5 6 7 1];
    %heatmap(ScenarioEvals(:,:,1));
    colormap 'gray'
    % Grid
    [X,Y]=meshgrid(1:7,1:13);
    % Two data sets to compare
    % A=rand(10,10);
    % B=-10.*rand(10,10);
    % Plot grid
    h=imagesc(X(:),Y(:),ScenarioEvalsRelDiff(estOrder,:,sc),[-5.5 1.5])
    % Text
    txt=sprintfc('%.2f',[ScenarioEvals(estOrder,:,sc) ScenarioEvalsRelDiff(estOrder,end,sc)])
    text(X(:),Y(:),txt,'horizontalalignment','center','verticalalignment','middle','fontsize',14)
    colorbar
    xlabels={'slope','\pmS.D','1:1 R^2*','<2','<10','<20','score'};
    ylabels={'\Omega','\Omega_T','\Omega_0','\Omega_C','\Omega_{TC}','\Omega_{0C}','Chao1','GP','Chao2','ACE','JK_a','JK_i','S_{obs}'};
    set(gca,'xaxisLocation','top','XTick',[1:7],'XTickLabel',xlabels,'YTick',[1:13],'YTickLabel',ylabels);
end