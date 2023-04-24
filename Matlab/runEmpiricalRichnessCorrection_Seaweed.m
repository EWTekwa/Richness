%runEmpirricalRichnessCorrection_Seaweed.m
%Eden Tekwa Mar 23, 2023
%analyze BC seaweed dataset

scrsz = get(0,'ScreenSize');
rng(6); %set random number generator seed

load('BC.seaweed_cover.mat'); %load .mat file containing the cell Data_all with species % cover data from Martone seaweed survey

% Data2012=importdata('community_2012.csv');
% Data2013=importdata('community_2013.csv');
% Data2014=importdata('community_2014.csv');
% Data2015=importdata('community_2015.csv');
% Data2016=importdata('community_2016.csv');
% Data2017=importdata('community_2017.csv');
% Data2018=importdata('community_2018.csv');
% Data2019=importdata('community_2019.csv');


%collect data from individual years (note: columns are species, rows are
%quadrats sequenced according to quadrats (10 consecutive rows are
%quadrats belonging to one quadrat)
Data_all={Data2012,Data2013,Data2014,Data2015,Data2016,Data2017,Data2018,Data2019};

ksub_perms=[9 2]; %pick out of 9 transects
ksub_sub_perms=[2 9]; %pick out of 10 quadrats in each transect

m_perms=[1 0.1]; %downsampling experiments within each quadrats (fraction of individuals observed)
%m_perms=[1]; %downsampling experiments within each quadrats (fraction of individuals observed)

NumYears=length(Data_all);
numTran=9; %number of transects
numQuad=10; %number of quadrats in each transect
YearLabels=[12:19];
numResample=40; %number of resamples for subsampled estimates (40)
numBoot=2000; %total number of bootstramps for CI (2000): numBoot/numResample per subsample
maxRich=300; %for display only
minRich=0; %for display only
diffRich=maxRich-minRich;

allSub_Chao1=zeros(length(m_perms)*(length(ksub_perms)+1),NumYears);
allSub_Chao2=zeros(length(m_perms)*(length(ksub_perms)+1),NumYears);
allSub_Richness_omega=zeros(length(m_perms)*(length(ksub_perms)+1),NumYears);
allSub_Richness_raw=zeros(length(m_perms)*(length(ksub_perms)+1),NumYears);

allSub_Chao1SD=zeros(length(m_perms)*(length(ksub_perms)+1),NumYears);
allSub_Chao2SD=zeros(length(m_perms)*(length(ksub_perms)+1),NumYears);
allSub_Richness_omegaSD=zeros(length(m_perms)*(length(ksub_perms)+1),NumYears);
allSub_Richness_rawSD=zeros(length(m_perms)*(length(ksub_perms)+1),NumYears);

allSubSlope_Chao1=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopeSD_Chao1=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopeP_Chao1=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopePPt_Chao1=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopeSDP_Chao1=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopeBeta_Chao1=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopeBetaPt_Chao1=zeros(length(m_perms)*(length(ksub_perms)+1),1);

allSubSlope_Chao2=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopeSD_Chao2=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopeP_Chao2=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopePPt_Chao2=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopeSDP_Chao2=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopeBeta_Chao2=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopeBetaPt_Chao2=zeros(length(m_perms)*(length(ksub_perms)+1),1);

allSubSlope_omega=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopeSD_omega=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopeP_omega=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopePPt_omega=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopeSDP_omega=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopeBeta_omega=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopeBetaPt_omega=zeros(length(m_perms)*(length(ksub_perms)+1),1);

allSubSlope_raw=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopeSD_raw=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopeP_raw=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopePPt_raw=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopeSDP_raw=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopeBeta_raw=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopeBetaPt_raw=zeros(length(m_perms)*(length(ksub_perms)+1),1);

figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/2 scrsz(4)/1.2],'DefaultAxesFontSize',16,'defaultAxesColorOrder',[[0 0 0]; [1 0 1]]);

i=1;
%prepare subsampling before simulations (so same set of resampled quadrats are used for paired downsampling experiments) 
transectsubIDs=cell(length(ksub_perms),numResample);
quadratsubIDs=cell(numQuad,length(ksub_sub_perms),numResample);
for testperm=1:length(ksub_perms)
    ksub=ksub_perms(testperm);
    ksub_sub=ksub_sub_perms(testperm);
    for resample=1:numResample
        transectsubIDs{testperm,resample}=randperm(numTran,ksub); %pick random transects for subsampling, same across years        
        for trans=1:ksub_perms(testperm)
            quadratsubIDs{trans,testperm,resample}=randperm(numQuad,ksub_sub); %pick different quadrats, same across year and transects      
        end
    end
end

Richness_raw_orig=zeros(1,NumYears);
min_cover=zeros(1,NumYears);
Samples=zeros(1,NumYears);
for yr=1:NumYears
    Richness_raw_orig(yr)=sum(sum(Data_all{yr}.data)>0);
    min_cover(yr)=min(Data_all{yr}.data(Data_all{yr}.data>0));
    Samples(yr)=size(Data_all{yr}.data,1);
end

%run simulations for all subsampling and paired downsampling scenarios
for mperm=1:length(m_perms)
    Richness_raw=zeros(1,NumYears);
    Chao1_all=zeros(1,NumYears);
    Chao2_all=zeros(1,NumYears);
    Richness_omega=zeros(1,NumYears);
    Richness_taylor=zeros(1,NumYears);
    min_cover=zeros(1,NumYears);
    meanStates=zeros(5,NumYears); %mean mn,P, var mn,P cov mn,P
    expectedRichness_raw=zeros(numBoot,NumYears);
    expectedChao1=zeros(numBoot,NumYears);
    expectedChao2=zeros(numBoot,NumYears);
    expectedRichness_omega=zeros(numBoot,NumYears);
    expectedRichness_taylor=zeros(numBoot,NumYears);
    ke_yr=zeros(1,NumYears);
    for yr=1:NumYears
        Data_all{yr}.data(Data_all{yr}.data<0.5 & Data_all{yr}.data>0)=0.5; %set nonzerocovers to minimum value of 0.5
        SpeciesIDs=find(sum(Data_all{yr}.data)>0);
        TransectAbundance=zeros(Samples(yr)/numQuad,Richness_raw_orig(yr));
        for species=1:Richness_raw_orig(yr)
            for transect=1:numTran %counting transect as a sample of the community
                TransectAbundance(transect,species)=round(sum(Data_all{yr}.data((transect-1)*10+1:(transect-1)*10+10,SpeciesIDs(species)))/0.5); %convert cover to count
            end
        end
        QuadratAbundance_m=zeros(size(TransectAbundance));
        if m_perms(mperm)==1
            QuadratAbundance_m(:,:,yr)=TransectAbundance; %use original dataset
        else
            QuadratAbundance_m(:,:,yr)=poissrnd(TransectAbundance*m_perms(mperm)); %downsample within each quadrat to a fraction of true abundances per species
        end
       %correct bootstrapping accounting for non-independence of individuals belonging to same species and quadrat:
        [Richness_raw(yr),Chao1_all(yr),~,Chao2_all(yr),~,~,~,Richness_omega(yr),Richness_taylor(yr),~,expectedRichness_raw(:,yr),expectedChao1(:,yr),~,expectedChao2(:,yr),~,~,~,expectedRichness_omega(:,yr),expectedRichness_taylor(:,yr),~,meanStates(:,yr)] = bootRichnessEsts(QuadratAbundance_m(:,:,yr),numBoot);
        %[Richness_raw(yr),Chao1(yr),GP(yr),Chao2(yr),ACE(yr),JK_a(yr),JK_i(yr),Richness_omega(yr),Richness_taylor(yr),Richness_taylor_0(yr),Omega_T_terms,meanStates] = RichnessEsts(QuadratAbundance_m(:,:,yr));

     end
    
    %record estimates for full dataset:
    allSub_Chao1((mperm-1)*(length(ksub_perms)+1)+1,:)=Chao1_all;
    allSub_Chao2((mperm-1)*(length(ksub_perms)+1)+1,:)=Chao2_all;
    allSub_Richness_omega((mperm-1)*(length(ksub_perms)+1)+1,:)=Richness_omega;
    allSub_Richness_raw((mperm-1)*(length(ksub_perms)+1)+1,:)=Richness_raw;
    
    %get temporal richness slope for centred Chao1 and omega boots
    slopeBoot=numBoot;
    RawDiff=(Richness_raw-mean(expectedRichness_raw));
    Chao1Diff=(Chao1_all-mean(expectedChao1));
    Chao2Diff=(Chao2_all-mean(expectedChao2));
    OmegaDiff=(Richness_omega-mean(expectedRichness_omega));  
    TaylorDiff=(Richness_taylor-mean(expectedRichness_taylor));
    RichnessSlope_raw_boot=zeros(slopeBoot,1);
    RichnessSlopeSD_raw_boot=zeros(slopeBoot,1);
    RichnessSlopeP_raw_boot=zeros(slopeBoot,1);
    RichnessSlope_Chao1_boot=zeros(slopeBoot,1);
    RichnessSlopeSD_Chao1_boot=zeros(slopeBoot,1);
    RichnessSlopeP_Chao1_boot=zeros(slopeBoot,1);
    RichnessSlope_Chao2_boot=zeros(slopeBoot,1);
    RichnessSlopeSD_Chao2_boot=zeros(slopeBoot,1);
    RichnessSlopeP_Chao2_boot=zeros(slopeBoot,1);
    RichnessSlope_Omega_boot=zeros(slopeBoot,1);
    RichnessSlopeSD_Omega_boot=zeros(slopeBoot,1);
    RichnessSlopeP_Omega_boot=zeros(slopeBoot,1);
    RichnessSlope_Taylor_boot=zeros(slopeBoot,1);
    RichnessSlopeSD_Taylor_boot=zeros(slopeBoot,1);
    RichnessSlopeP_Taylor_boot=zeros(slopeBoot,1);
    %get temporal slopes:
    %regress based on point richness estimates:
    mdl_Richness_raw=fitlm([1:NumYears],Richness_raw);
    mdl_Chao1=fitlm([1:NumYears],Chao1_all);
    mdl_Chao2=fitlm([1:NumYears],Chao2_all);
    mdl_Omega=fitlm([1:NumYears],Richness_omega);
    mdl_Taylor=fitlm([1:NumYears],Richness_taylor);
    RichnessSlope_raw=mdl_Richness_raw.Coefficients.Estimate(2);
    RichnessSlopeSD_raw=mdl_Richness_raw.Coefficients.SE(2)*sqrt(NumYears);
    RichnessSlopeP_raw=mdl_Richness_raw.Coefficients.pValue(2);
    RichnessSlope_Chao1=mdl_Chao1.Coefficients.Estimate(2);
    RichnessSlopeSD_Chao1=mdl_Chao1.Coefficients.SE(2)*sqrt(NumYears);
    RichnessSlopeP_Chao1=mdl_Chao1.Coefficients.pValue(2);
    RichnessSlope_Chao2=mdl_Chao2.Coefficients.Estimate(2);
    RichnessSlopeSD_Chao2=mdl_Chao2.Coefficients.SE(2)*sqrt(NumYears);
    RichnessSlopeP_Chao2=mdl_Chao2.Coefficients.pValue(2);
    RichnessSlope_Omega=mdl_Omega.Coefficients.Estimate(2);
    RichnessSlopeSD_Omega=mdl_Omega.Coefficients.SE(2)*sqrt(NumYears);
    RichnessSlopeP_Omega=mdl_Omega.Coefficients.pValue(2);
    RichnessSlope_Taylor=mdl_Taylor.Coefficients.Estimate(2);
    RichnessSlopeSD_Taylor=mdl_Taylor.Coefficients.SE(2)*sqrt(NumYears);
    RichnessSlopeP_Taylor=mdl_Taylor.Coefficients.pValue(2);
    %regress based on bootstrapped richness estimates:
    for boot=1:slopeBoot
        mdl_Richness_raw_boot=fitlm([1:NumYears],expectedRichness_raw(boot,:)+RawDiff);
        mdl_Chao1_boot=fitlm([1:NumYears],expectedChao1(boot,:)+Chao1Diff);
        mdl_Chao2_boot=fitlm([1:NumYears],expectedChao2(boot,:)+Chao2Diff);
        mdl_Omega_boot=fitlm([1:NumYears],expectedRichness_omega(boot,:)+OmegaDiff);
        mdl_Taylor_boot=fitlm([1:NumYears],expectedRichness_taylor(boot,:)+TaylorDiff);
        %end
        RichnessSlope_raw_boot(boot)=mdl_Richness_raw_boot.Coefficients.Estimate(2);
        RichnessSlopeSD_raw_boot(boot)=mdl_Richness_raw_boot.Coefficients.SE(2)*sqrt(NumYears);
        RichnessSlopeP_raw_boot(boot)=mdl_Richness_raw_boot.Coefficients.pValue(2);
        RichnessSlope_Chao1_boot(boot)=mdl_Chao1_boot.Coefficients.Estimate(2);
        RichnessSlopeSD_Chao1_boot(boot)=mdl_Chao1_boot.Coefficients.SE(2)*sqrt(NumYears);
        RichnessSlopeP_Chao1_boot(boot)=mdl_Chao1_boot.Coefficients.pValue(2);
        RichnessSlope_Chao2_boot(boot)=mdl_Chao2_boot.Coefficients.Estimate(2);
        RichnessSlopeSD_Chao2_boot(boot)=mdl_Chao2_boot.Coefficients.SE(2)*sqrt(NumYears);
        RichnessSlopeP_Chao2_boot(boot)=mdl_Chao2_boot.Coefficients.pValue(2);
        RichnessSlope_Omega_boot(boot)=mdl_Omega_boot.Coefficients.Estimate(2);
        RichnessSlopeSD_Omega_boot(boot)=mdl_Omega_boot.Coefficients.SE(2)*sqrt(NumYears);
        RichnessSlopeP_Omega_boot(boot)=mdl_Omega_boot.Coefficients.pValue(2);
        RichnessSlope_Taylor_boot(boot)=mdl_Taylor_boot.Coefficients.Estimate(2);
        RichnessSlopeSD_Taylor_boot(boot)=mdl_Taylor_boot.Coefficients.SE(2)*sqrt(NumYears);
        RichnessSlopeP_Taylor_boot(boot)=mdl_Taylor_boot.Coefficients.pValue(2);
    end
    
    
    %record estimated richness SD from full dataset:
    allSub_Chao1SD((mperm-1)*(length(ksub_perms)+1)+1,:)=std(expectedChao1);
    allSub_Chao2SD((mperm-1)*(length(ksub_perms)+1)+1,:)=std(expectedChao2);
    allSub_Richness_omegaSD((mperm-1)*(length(ksub_perms)+1)+1,:)=std(expectedRichness_omega);
    allSub_Richness_rawSD((mperm-1)*(length(ksub_perms)+1)+1,:)=std(expectedRichness_raw);
    
    %record trend estimates from full dataset:
    allSubSlope_Chao1((mperm-1)*(length(ksub_perms)+1)+1,:)=mean(RichnessSlope_Chao1_boot(:));
    allSubSlopeSD_Chao1((mperm-1)*(length(ksub_perms)+1)+1,:)=mean(RichnessSlopeSD_Chao1_boot(:));
    allSubSlopeP_Chao1((mperm-1)*(length(ksub_perms)+1)+1,:)=mean(RichnessSlopeP_Chao1_boot(:));
    allSubSlopePPt_Chao1((mperm-1)*(length(ksub_perms)+1)+1,:)=mean(RichnessSlopeP_Chao1);
    allSubSlopeSDP_Chao1((mperm-1)*(length(ksub_perms)+1)+1,:)=std(RichnessSlopeP_Chao1_boot(:));
    allSubSlopeBeta_Chao1((mperm-1)*(length(ksub_perms)+1)+1,:)=sum(RichnessSlopeP_Chao1_boot(:)>0.05)/numBoot;
    allSubSlopeBetaPt_Chao1((mperm-1)*(length(ksub_perms)+1)+1,:)=sum(RichnessSlopeP_Chao1(:)>0.05)/numResample;
    
    allSubSlope_Chao2((mperm-1)*(length(ksub_perms)+1)+1,:)=mean(RichnessSlope_Chao2_boot(:));
    allSubSlopeSD_Chao2((mperm-1)*(length(ksub_perms)+1)+1,:)=mean(RichnessSlopeSD_Chao2_boot(:));
    allSubSlopeP_Chao2((mperm-1)*(length(ksub_perms)+1)+1,:)=mean(RichnessSlopeP_Chao2_boot(:));
    allSubSlopePPt_Chao2((mperm-1)*(length(ksub_perms)+1)+1,:)=mean(RichnessSlopeP_Chao2);
    allSubSlopeSDP_Chao2((mperm-1)*(length(ksub_perms)+1)+1,:)=std(RichnessSlopeP_Chao2_boot(:));
    allSubSlopeBeta_Chao2((mperm-1)*(length(ksub_perms)+1)+1,:)=sum(RichnessSlopeP_Chao2_boot(:)>0.05)/numBoot;
    allSubSlopeBetaPt_Chao2((mperm-1)*(length(ksub_perms)+1)+1,:)=sum(RichnessSlopeP_Chao2(:)>0.05)/numBoot;
    
    allSubSlope_omega((mperm-1)*(length(ksub_perms)+1)+1,:)=mean(RichnessSlope_Omega_boot(:));
    allSubSlopeSD_omega((mperm-1)*(length(ksub_perms)+1)+1,:)=mean(RichnessSlopeSD_Omega_boot(:));
    allSubSlopeP_omega((mperm-1)*(length(ksub_perms)+1)+1,:)=mean(RichnessSlopeP_Omega_boot(:));
    allSubSlopePPt_omega((mperm-1)*(length(ksub_perms)+1)+1,:)=mean(RichnessSlopeP_Omega);
    allSubSlopeSDP_omega((mperm-1)*(length(ksub_perms)+1)+1,:)=std(RichnessSlopeP_Omega_boot(:));
    allSubSlopeBeta_omega((mperm-1)*(length(ksub_perms)+1)+1,:)=sum(RichnessSlopeP_Omega_boot(:)>0.05)/numBoot;
    allSubSlopeBetaPt_omega((mperm-1)*(length(ksub_perms)+1)+1,:)=sum(RichnessSlopeP_Omega(:)>0.05)/numBoot;
    
    allSubSlope_raw((mperm-1)*(length(ksub_perms)+1)+1,:)=mean(RichnessSlope_raw_boot(:));
    allSubSlopeSD_raw((mperm-1)*(length(ksub_perms)+1)+1,:)=mean(RichnessSlopeSD_raw_boot(:));
    allSubSlopeP_raw((mperm-1)*(length(ksub_perms)+1)+1,:)=mean(RichnessSlopeP_raw_boot(:));
    allSubSlopePPt_raw((mperm-1)*(length(ksub_perms)+1)+1,:)=mean(RichnessSlopeP_raw);
    allSubSlopeSDP_raw((mperm-1)*(length(ksub_perms)+1)+1,:)=std(RichnessSlopeP_raw_boot(:));
    allSubSlopeBeta_raw((mperm-1)*(length(ksub_perms)+1)+1,:)=sum(RichnessSlopeP_raw_boot(:)>0.05)/numBoot;
    allSubSlopeBetaPt_raw((mperm-1)*(length(ksub_perms)+1)+1,:)=sum(RichnessSlopeP_raw(:)>0.05)/numBoot;

    %plot selected richness estimates across years on full spatial samples:
    %plotIDs=[1:]
    subplot(length(m_perms)*2,length(ksub_perms)+1,(mperm-1)*(length(ksub_perms)+1)+1)
    %subplot(length(m_perms)*2,length(ksub_perms)+2,(mperm-1)*(length(ksub_perms)+2)+1)
    %yyaxis left
    hold on
    %centre confidence bounds to match bootstrapped mean to spot estimates
    %on the original dataset:
%     boundedline([1:NumYears]', Chao2_all',[max(mean(expectedChao2)'-prctile(expectedChao2,2.5)',0),prctile(expectedChao1,97.5)'-mean(expectedChao1)'],'--c','alpha','transparency', 0.1);
%     boundedline([1:NumYears]', Chao1_all',[max(mean(expectedChao1)'-prctile(expectedChao1,2.5)',0),prctile(expectedChao1,97.5)'-mean(expectedChao1)'],'-b','alpha','transparency', 0.1);
%     boundedline([1:NumYears]', Richness_raw',[max(mean(expectedRichness_raw)'-prctile(expectedRichness_raw,2.5)',0),prctile(expectedRichness_raw,97.5)'-mean(expectedRichness_raw)'],'--k','alpha','transparency', 0.1);
%     boundedline([1:NumYears]', Richness_omega',[max(mean(expectedRichness_omega)'-prctile(expectedRichness_omega,2.5)',0),prctile(expectedRichness_omega,97.5)'-mean(expectedRichness_omega)'],'-r','alpha','transparency', 0.1);
    boundedline([1:NumYears]', Chao2_all',std(expectedChao2)','--c','alpha','transparency', 0.1);
    boundedline([1:NumYears]', Chao1_all',std(expectedChao1)','-b','alpha','transparency', 0.1);
    boundedline([1:NumYears]', Richness_omega',std(Richness_omega)','-r','alpha','transparency', 0.1);
    boundedline([1:NumYears]', Richness_raw',std(Richness_raw)','--k','alpha','transparency', 0.1);
    
    plot(Chao1_all,'-b','LineWidth',2);
    plot(Chao2_all,'--','LineWidth',2,'Color',[0 0.7 0.7]);
    plot(Richness_omega,'-r','LineWidth',2);
    plot(Richness_raw,'--k','LineWidth',2);

    xlabel 'year'
    ylabel 'richness'
    xticks([1:NumYears])
    xticklabels(YearLabels)
    xlim([1,NumYears])
    title({['9 transects x 10 quadrats'];[num2str(m_perms(mperm)*100) '% individuals observed']},'fontweight','normal')
    %plotLabels{(mperm-1)*(length(ksub_perms)+1)+1}=['9x10x' num2str(m_perms(mperm))];
    ylim([minRich maxRich])
%     text(1.2,minRich+diffRich*.95,['\DeltaRaw=' num2str(mean(RichnessSlope_raw_boot(:)),1) '\pm' num2str(mean(RichnessSlopeSD_raw_boot(:)),1) ', p=' num2str(mean(RichnessSlopeP_raw_boot(:)),1) '(' num2str(mean(RichnessSlopeP_raw),1) ')\pm' num2str(std(RichnessSlopeP_raw_boot(:)),1) ', \alpha=' num2str(sum(RichnessSlopeP_raw_boot(:)<=0.05)/numBoot,1)],'fontsize',14);
%     text(1.2,minRich+diffRich*.85,['\DeltaChao1=' num2str(mean(RichnessSlope_Chao1_boot(:)),1) '\pm' num2str(mean(RichnessSlopeSD_Chao1_boot(:)),1) ', p=' num2str(mean(RichnessSlopeP_Chao1_boot(:)),1) '(' num2str(mean(RichnessSlopeP_Chao1),1) ')\pm' num2str(std(RichnessSlopeP_Chao1_boot(:)),1) ', \alpha=' num2str(sum(RichnessSlopeP_Chao1_boot(:)<=0.05)/numBoot,1)],'fontsize',14,'color','b');
%     text(1.2,minRich+diffRich*.75,['\DeltaChao2=' num2str(mean(RichnessSlope_Chao2_boot(:)),1) '\pm' num2str(mean(RichnessSlopeSD_Chao2_boot(:)),1) ', p=' num2str(mean(RichnessSlopeP_Chao2_boot(:)),1) '(' num2str(mean(RichnessSlopeP_Chao2),1) ')\pm' num2str(std(RichnessSlopeP_Chao2_boot(:)),1) ', \alpha=' num2str(sum(RichnessSlopeP_Chao2_boot(:)<=0.05)/numBoot,1)],'fontsize',14,'color',[0 0.7 0.7]);
%     text(1.2,minRich+diffRich*.63,['\Delta\Omega_T=' num2str(mean(RichnessSlope_Taylor_boot(:)),1) '\pm' num2str(mean(RichnessSlopeSD_Taylor_boot(:)),1) ', p=' num2str(mean(RichnessSlopeP_Taylor_boot(:)),1) '(' num2str(mean(RichnessSlopeP_Taylor),1) ')\pm' num2str(std(RichnessSlopeP_Taylor_boot(:)),1) ', \alpha=' num2str(sum(RichnessSlopeP_Taylor_boot(:)<=0.05)/numBoot,1)],'fontsize',14,'Color','r');
%     text(1,2,minRich+diffRich*.95,['bootstrapped estimates:'],'fontsize',14)
%     text(1.2,minRich+diffRich*.85,['raw slope=' num2str(mean(RichnessSlope_raw_boot),1) '\pm' num2str(mean(RichnessSlopeSD_raw_boot),1) ', p=' num2str(mean(RichnessSlopeP_raw_boot),1) '\pm' num2str(std(RichnessSlopeP_raw_boot),1) ', \beta=' num2str(sum(RichnessSlopeP_raw_boot(:)>0.05)/numBoot,1)],'fontsize',14);
%     text(1.2,minRich+diffRich*.75,['Chao1 slope=' num2str(mean(RichnessSlope_Chao1_boot),1) '\pm' num2str(mean(RichnessSlopeSD_Chao1_boot),1) ', p=' num2str(mean(RichnessSlopeP_Chao1_boot),1) '\pm' num2str(std(RichnessSlopeP_Chao1_boot),1) ', \beta=' num2str(sum(RichnessSlopeP_Chao1_boot(:)>0.05)/numBoot,1)],'fontsize',14,'color','b');
%     text(1.2,minRich+diffRich*.63,['\Omega_o slope=' num2str(mean(RichnessSlope_Omega_boot),1) '\pm' num2str(mean(RichnessSlopeSD_Omega_boot),1) ', p=' num2str(mean(RichnessSlopeP_Omega_boot),1) '\pm' num2str(std(RichnessSlopeP_Omega_boot),1) ', \beta=' num2str(sum(RichnessSlopeP_Omega_boot(:)>0.05)/numBoot,1)],'fontsize',14,'Color','r');
%     text(1,2,minRich+diffRich*.45,['point estimates:'],'fontsize',14)
%     text(1.2,minRich+diffRich*.35,['raw slope=' num2str(mean(RichnessSlope_raw),1) '\pm' num2str(mean(RichnessSlopeSD_raw),1) ', p=' num2str(mean(RichnessSlopeP_raw),1)],'fontsize',14);
%     text(1.2,minRich+diffRich*.25,['Chao1 slope=' num2str(mean(RichnessSlope_Chao1),1) '\pm' num2str(mean(RichnessSlopeSD_Chao1),1) ', p=' num2str(mean(RichnessSlopeP_Chao1),1)],'fontsize',14,'color','b');
%     text(1.2,minRich+diffRich*.13,['\Omega_o slope=' num2str(mean(RichnessSlope_Omega),1) '\pm' num2str(mean(RichnessSlopeSD_Omega),1) ', p=' num2str(mean(RichnessSlopeP_Omega),1)],'fontsize',14,'Color','r');
    text(1-(NumYears-1)*0.15,maxRich*1.15,char(64+i),'Fontsize',16)
    
    %plot fractions of species recovered in subsampled and downsampled
    %experiments
    
    %plot observed mean states
    subplot(4,length(ksub_perms)+1,(mperm-1)*(length(ksub_perms)+1)+1+6)
    hold on
    plot(meanStates(1,:),'-','LineWidth',2,'color','b'); %E[mn]
    plot(meanStates(2,:),'-','LineWidth',2,'color',[255 139 37]/255); %E[P]
    plot(meanStates(3,:),'--','LineWidth',2,'color','b'); %var[mn]
    plot(meanStates(4,:),'--','LineWidth',2,'color',[255 139 37]/255); %var[P]
    plot(meanStates(5,:),'--','LineWidth',2,'color',[69 181 80]/255); %cov[mn,P]
    refline(0,1)
    ax=gca;
    ax.YAxis.Scale= 'log';
    ylim([0.01 2*10^5])
    yticks([0.01 1 10^2 10^4])
    ylabel('state')
    ylims=ylim;
    xlabel 'year'
    xticks([1:NumYears])
    xticklabels(YearLabels)
    xlim([1,NumYears])
    title({['9 transects x 10 quadrats'];[num2str(m_perms(mperm)*100) '% individuals observed']},'fontweight','normal')
    text(1-(NumYears-1)*0.15,ylims(2)*3,char(64+i),'Fontsize',16)
    i=i+1;
    
    for testperm=1:length(ksub_perms)
        subplot(length(m_perms)*2,length(ksub_perms)+1,(mperm-1)*(length(ksub_perms)+1)+testperm+1) %plot estimated richness by year
        hold on
        ksub=ksub_perms(testperm);
        ksub_sub=ksub_sub_perms(testperm);
        
        Richness_raw_sub=zeros(numResample,NumYears); %point estimates
        Chao1_sub=zeros(numResample,NumYears); %point estimates
        Chao2_sub=zeros(numResample,NumYears); %point estimates
        Richness_omega_sub=zeros(numResample,NumYears); %point estimates
        Richness_taylor_sub=zeros(numResample,NumYears); %point estimates
        Richness_raw_sub_mean=zeros(numResample,NumYears); %bootstrap means across resamples
        Chao1_sub_mean=zeros(numResample,NumYears);
        Chao2_sub_mean=zeros(numResample,NumYears);
        Richness_omega_sub_mean=zeros(numResample,NumYears);
        Richness_taylor_sub_mean=zeros(numResample,NumYears);
        Richness_raw_sub_lo=zeros(numResample,NumYears); %bootstrap lower bound across resamples (95%)
        Chao1_sub_lo=zeros(numResample,NumYears);
        Chao2_sub_lo=zeros(numResample,NumYears);
        Richness_omega_sub_lo=zeros(numResample,NumYears);
        Richness_taylor_sub_lo=zeros(numResample,NumYears);
        Richness_raw_sub_up=zeros(numResample,NumYears); %bootstrap upper bounds across resamples (95%)
        Chao1_sub_up=zeros(numResample,NumYears);
        Chao2_sub_up=zeros(numResample,NumYears);
        Richness_omega_sub_up=zeros(numResample,NumYears);
        Richness_taylor_sub_up=zeros(numResample,NumYears);
        meanStates=zeros(5,numResample,NumYears);
        ke_yr=zeros(1,NumYears);
        
 %get temporal richness slope for raw and centred Chao1 and omega
        %boots in subsampling experiments:
        subBoots=round(numBoot/numResample);
        RichnessSlope_raw=zeros(numResample,1);
        RichnessSlopeSD_raw=zeros(numResample,1);
        RichnessSlopeP_raw=zeros(numResample,1);
        RichnessSlope_Chao1=zeros(numResample,1);
        RichnessSlopeSD_Chao1=zeros(numResample,1);
        RichnessSlopeP_Chao1=zeros(numResample,1);
        RichnessSlope_Chao2=zeros(numResample,1);
        RichnessSlopeSD_Chao2=zeros(numResample,1);
        RichnessSlopeP_Chao2=zeros(numResample,1);
        RichnessSlope_Omega=zeros(numResample,1);
        RichnessSlopeSD_Omega=zeros(numResample,1);
        RichnessSlopeP_Omega=zeros(numResample,1);
        RichnessSlope_Taylor=zeros(numResample,1);
        RichnessSlopeSD_Taylor=zeros(numResample,1);
        RichnessSlopeP_Taylor=zeros(numResample,1);
        
        RichnessSlope_raw_boot=zeros(numResample,subBoots);
        RichnessSlopeSD_raw_boot=zeros(numResample,subBoots);
        RichnessSlopeP_raw_boot=zeros(numResample,subBoots);
        RichnessSlope_Chao1_boot=zeros(numResample,subBoots);
        RichnessSlopeSD_Chao1_boot=zeros(numResample,subBoots);
        RichnessSlopeP_Chao1_boot=zeros(numResample,subBoots);
        RichnessSlope_Chao2_boot=zeros(numResample,subBoots);
        RichnessSlopeSD_Chao2_boot=zeros(numResample,subBoots);
        RichnessSlopeP_Chao2_boot=zeros(numResample,subBoots);
        RichnessSlope_Omega_boot=zeros(numResample,subBoots);
        RichnessSlopeSD_Omega_boot=zeros(numResample,subBoots);
        RichnessSlopeP_Omega_boot=zeros(numResample,subBoots);
        RichnessSlope_Taylor_boot=zeros(numResample,subBoots);
        RichnessSlopeSD_Taylor_boot=zeros(numResample,subBoots);
        RichnessSlopeP_Taylor_boot=zeros(numResample,subBoots);
        
        expectedRichness_raw_yrs=zeros(numBoot,NumYears);
        expectedChao1_yrs=zeros(numBoot,NumYears);
        expectedChao2_yrs=zeros(numBoot,NumYears);
        expectedRichness_omega_yrs=zeros(numBoot,NumYears);
        expectedRichness_taylor_yrs=zeros(numBoot,NumYears);
        
        for resample=1:numResample
            expectedRichness_raw=zeros(subBoots,NumYears);
            expectedChao1=zeros(subBoots,NumYears);
            expectedChao2=zeros(subBoots,NumYears);
            expectedRichness_omega=zeros(subBoots,NumYears);
            expectedRichness_taylor=zeros(subBoots,NumYears);
            transectsubID=transectsubIDs{testperm,resample}; %load transet set

            for yr=1:NumYears
                %TransectAbundance_sub=QuadratAbundance_m(:,:,yr);
                SpeciesIDs=find(sum(Data_all{yr}.data)>0);
                quadratAbundance_sub=zeros(ksub,Richness_raw_orig(yr));
                for species=1:Richness_raw_orig(yr)
                    for transect_sub=1:ksub %counting transect as a sample of the community
                        quadratsubID=quadratsubIDs{transect_sub,testperm,resample}; %load quadrat set
                        quadratAbundance_sub(transect_sub,species)=round(sum(Data_all{yr}.data((transectsubID(transect_sub)-1)*10+quadratsubID,SpeciesIDs(species)))/0.5);
                    end
                end
                if m_perms(mperm)~=1
                    quadratAbundance_sub=poissrnd(quadratAbundance_sub*m_perms(mperm)); %downsample within each quadrat to a fraction of true abundances per species
                end
                [Richness_raw_sub(resample,yr),Chao1_sub(resample,yr),~,Chao2_sub(resample,yr),~,~,~,Richness_omega_sub(resample,yr),Richness_taylor_sub(resample,yr),~,expectedRichness_raw_boot,expectedChao1_boot,~,expectedChao2_boot,~,~,~,expectedRichness_omega_boot,expectedRichness_taylor_boot,~,meanStates(:,resample,yr)] = bootRichnessEsts(quadratAbundance_sub,subBoots);            
                expectedRichness_raw(:,yr)=expectedRichness_raw_boot;
                expectedChao1(:,yr)=expectedChao1_boot;
                expectedChao2(:,yr)=expectedChao2_boot;
                expectedRichness_omega(:,yr)=expectedRichness_omega_boot;
                expectedRichness_taylor(:,yr)=expectedRichness_taylor_boot;
                expectedRichness_raw_yrs((resample-1)*subBoots+1:resample*subBoots,yr)=expectedRichness_raw_boot;
                expectedChao1_yrs((resample-1)*subBoots+1:resample*subBoots,yr)=expectedChao1_boot;
                expectedChao2_yrs((resample-1)*subBoots+1:resample*subBoots,yr)=expectedChao2_boot;
                expectedRichness_omega_yrs((resample-1)*subBoots+1:resample*subBoots,yr)=expectedRichness_omega_boot;
                expectedRichness_taylor_yrs((resample-1)*subBoots+1:resample*subBoots,yr)=expectedRichness_taylor_boot;
                Richness_raw_sub_mean(resample,yr)=mean(expectedRichness_raw_boot);
                Chao1_sub_mean(resample,yr)=mean(expectedChao1_boot);
                Chao2_sub_mean(resample,yr)=mean(expectedChao2_boot);
                Richness_omega_sub_mean(resample,yr)=mean(expectedRichness_omega_boot);
                Richness_taylor_sub_mean(resample,yr)=mean(expectedRichness_taylor_boot);
                Richness_raw_sub_lo(resample,yr)=prctile(expectedRichness_raw_boot,2.5);
                Chao1_sub_lo(resample,yr)=prctile(expectedChao1_boot,2.5);
                Chao2_sub_lo(resample,yr)=prctile(expectedChao2_boot,2.5);
                Richness_omega_sub_lo(resample,yr)=prctile(expectedRichness_omega_boot,2.5);
                Richness_taylor_sub_lo(resample,yr)=prctile(expectedRichness_taylor_boot,2.5);
                Richness_raw_sub_up(resample,yr)=prctile(expectedRichness_raw_boot,97.5);
                Chao1_sub_up(resample,yr)=prctile(expectedChao1_boot,97.5);
                Chao2_sub_up(resample,yr)=prctile(expectedChao2_boot,97.5);
                Richness_omega_sub_up(resample,yr)=prctile(expectedRichness_omega_boot,97.5);
                Richness_taylor_sub_up(resample,yr)=prctile(expectedRichness_taylor_boot,97.5);
    
            end
            RawDiff=(Richness_raw_sub(resample,:)-mean(expectedRichness_raw));
            Chao1Diff=(Chao1_sub(resample,:)-mean(expectedChao1));
            Chao2Diff=(Chao2_sub(resample,:)-mean(expectedChao2));
            OmegaDiff=(Richness_omega_sub(resample,:)-mean(expectedRichness_omega));
            TaylorDiff=(Richness_taylor_sub(resample,:)-mean(expectedRichness_taylor));
            mdl_Richness_raw=fitlm([1:NumYears],Richness_raw_sub(resample,:));
            mdl_Chao1=fitlm([1:NumYears],Chao1_sub(resample,:));
            mdl_Chao2=fitlm([1:NumYears],Chao2_sub(resample,:));
            mdl_Omega=fitlm([1:NumYears],Richness_omega_sub(resample,:));
            mdl_Taylor=fitlm([1:NumYears],Richness_taylor_sub(resample,:));
            RichnessSlope_raw(resample)=mdl_Richness_raw.Coefficients.Estimate(2);
            RichnessSlopeSD_raw(resample)=mdl_Richness_raw.Coefficients.SE(2)*sqrt(NumYears);
            RichnessSlopeP_raw(resample)=mdl_Richness_raw.Coefficients.pValue(2);
            RichnessSlope_Chao1(resample)=mdl_Chao1.Coefficients.Estimate(2);
            RichnessSlopeSD_Chao1(resample)=mdl_Chao1.Coefficients.SE(2)*sqrt(NumYears);
            RichnessSlopeP_Chao1(resample)=mdl_Chao1.Coefficients.pValue(2);
            RichnessSlope_Chao2(resample)=mdl_Chao2.Coefficients.Estimate(2);
            RichnessSlopeSD_Chao2(resample)=mdl_Chao2.Coefficients.SE(2)*sqrt(NumYears);
            RichnessSlopeP_Chao2(resample)=mdl_Chao2.Coefficients.pValue(2);
            RichnessSlope_Omega(resample)=mdl_Omega.Coefficients.Estimate(2);
            RichnessSlopeSD_Omega(resample)=mdl_Omega.Coefficients.SE(2)*sqrt(NumYears);
            RichnessSlopeP_Omega(resample)=mdl_Omega.Coefficients.pValue(2);
            RichnessSlope_Taylor(resample)=mdl_Taylor.Coefficients.Estimate(2);
            RichnessSlopeSD_Taylor(resample)=mdl_Taylor.Coefficients.SE(2)*sqrt(NumYears);
            RichnessSlopeP_Taylor(resample)=mdl_Taylor.Coefficients.pValue(2);
            for boot=1:subBoots
                mdl_Richness_raw_boot=fitlm([1:NumYears],expectedRichness_raw(boot,:)+RawDiff);
                mdl_Chao1_boot=fitlm([1:NumYears],expectedChao1(boot,:)+Chao1Diff);
                mdl_Chao2_boot=fitlm([1:NumYears],expectedChao2(boot,:)+Chao2Diff);
                mdl_Omega_boot=fitlm([1:NumYears],expectedRichness_omega(boot,:)+OmegaDiff);
                mdl_Taylor_boot=fitlm([1:NumYears],expectedRichness_taylor(boot,:)+TaylorDiff);
                %end
                RichnessSlope_raw_boot(resample,boot)=mdl_Richness_raw_boot.Coefficients.Estimate(2);
                RichnessSlopeSD_raw_boot(resample,boot)=mdl_Richness_raw_boot.Coefficients.SE(2)*sqrt(NumYears);
                RichnessSlopeP_raw_boot(resample,boot)=mdl_Richness_raw_boot.Coefficients.pValue(2);
                RichnessSlope_Chao1_boot(resample,boot)=mdl_Chao1_boot.Coefficients.Estimate(2);
                RichnessSlopeSD_Chao1_boot(resample,boot)=mdl_Chao1_boot.Coefficients.SE(2)*sqrt(NumYears);
                RichnessSlopeP_Chao1_boot(resample,boot)=mdl_Chao1_boot.Coefficients.pValue(2);
                RichnessSlope_Chao2_boot(resample,boot)=mdl_Chao2_boot.Coefficients.Estimate(2);
                RichnessSlopeSD_Chao2_boot(resample,boot)=mdl_Chao2_boot.Coefficients.SE(2)*sqrt(NumYears);
                RichnessSlopeP_Chao2_boot(resample,boot)=mdl_Chao2_boot.Coefficients.pValue(2);
                RichnessSlope_Omega_boot(resample,boot)=mdl_Omega_boot.Coefficients.Estimate(2);
                RichnessSlopeSD_Omega_boot(resample,boot)=mdl_Omega_boot.Coefficients.SE(2)*sqrt(NumYears);
                RichnessSlopeP_Omega_boot(resample,boot)=mdl_Omega_boot.Coefficients.pValue(2);
                RichnessSlope_Taylor_boot(resample,boot)=mdl_Taylor_boot.Coefficients.Estimate(2);
                RichnessSlopeSD_Taylor_boot(resample,boot)=mdl_Taylor_boot.Coefficients.SE(2)*sqrt(NumYears);
                RichnessSlopeP_Taylor_boot(resample,boot)=mdl_Taylor_boot.Coefficients.pValue(2);
            end
        end
        %record estimates for subsampled+downsampled dataset:
        allSub_Chao1((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=mean(Chao1_sub);
        allSub_Chao2((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=mean(Chao2_sub);
        allSub_Richness_omega((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=mean(Richness_omega_sub);
        allSub_Richness_raw((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=mean(Richness_raw_sub);
        
        
        allSub_Chao1SD((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=std(expectedChao1_yrs);
        allSub_Chao2SD((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=std(expectedChao2_yrs);
        allSub_Richness_omegaSD((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=std(expectedRichness_omega_yrs);
        allSub_Richness_rawSD((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=std(expectedRichness_raw_yrs);
        
        %record trend estimates for subsampled+downsampled dataset:
        allSubSlope_Chao1((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=mean(RichnessSlope_Chao1_boot(:));
        allSubSlopeSD_Chao1((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=mean(RichnessSlopeSD_Chao1_boot(:));
        allSubSlopeP_Chao1((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=mean(RichnessSlopeP_Chao1_boot(:));
        allSubSlopePPt_Chao1((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=mean(RichnessSlopeP_Chao1);
        allSubSlopeSDP_Chao1((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=std(RichnessSlopeP_Chao1_boot(:));
        allSubSlopeBeta_Chao1((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=sum(RichnessSlopeP_Chao1_boot(:)>0.05)/numBoot;
        allSubSlopeBetaPt_Chao1((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=sum(RichnessSlopeP_Chao1(:)>0.05)/numResample;
        
        allSubSlope_Chao2((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=mean(RichnessSlope_Chao2_boot(:));
        allSubSlopeSD_Chao2((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=mean(RichnessSlopeSD_Chao2_boot(:));
        allSubSlopeP_Chao2((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=mean(RichnessSlopeP_Chao2_boot(:));
        allSubSlopePPt_Chao2((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=mean(RichnessSlopeP_Chao2);
        allSubSlopeSDP_Chao2((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=std(RichnessSlopeP_Chao2_boot(:));
        allSubSlopeBeta_Chao2((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=sum(RichnessSlopeP_Chao2_boot(:)>0.05)/numBoot;
        allSubSlopeBetaPt_Chao2((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=sum(RichnessSlopeP_Chao2(:)>0.05)/numResample;
        
        allSubSlope_omega((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=mean(RichnessSlope_Omega_boot(:));
        allSubSlopeSD_omega((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=mean(RichnessSlopeSD_Omega_boot(:));
        allSubSlopeP_omega((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=mean(RichnessSlopeP_Omega_boot(:));
        allSubSlopePPt_omega((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=mean(RichnessSlopeP_Omega);
        allSubSlopeSDP_omega((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=std(RichnessSlopeP_Omega_boot(:));
        allSubSlopeBeta_omega((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=sum(RichnessSlopeP_Omega_boot(:)>0.05)/numBoot;
        allSubSlopeBetaPt_omega((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=sum(RichnessSlopeP_Omega(:)>0.05)/numResample;
        
        allSubSlope_raw((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=mean(RichnessSlope_raw_boot(:));
        allSubSlopeSD_raw((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=mean(RichnessSlopeSD_raw_boot(:));
        allSubSlopeP_raw((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=mean(RichnessSlopeP_raw_boot(:));
        allSubSlopePPt_raw((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=mean(RichnessSlopeP_raw);
        allSubSlopeSDP_raw((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=std(RichnessSlopeP_raw_boot(:));
        allSubSlopeBeta_raw((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=sum(RichnessSlopeP_raw_boot(:)>0.05)/numBoot;
        allSubSlopeBetaPt_raw((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=sum(RichnessSlopeP_raw(:)>0.05)/numResample;
        
        %centre confidence bounds to match bootstrapped mean to spot estimates
        %on the original dataset:
%         boundedline([1:NumYears]', mean(Chao2_sub)',[max(mean(Chao2_sub_mean-Chao2_sub_lo)',0),mean(Chao2_sub_up-Chao2_sub_mean)'],'--c','alpha','transparency', 0.1);
%         boundedline([1:NumYears]', mean(Chao1_sub)',[max(mean(Chao1_sub_mean-Chao1_sub_lo)',0),mean(Chao1_sub_up-Chao1_sub_mean)'],'-b','alpha','transparency', 0.1);
%         boundedline([1:NumYears]', mean(Richness_omega_sub)',[max(mean(Richness_omega_sub_mean-Richness_omega_sub_lo)',0),mean(Richness_omega_sub_up-Richness_omega_sub_mean)'],'-r','alpha','transparency', 0.1);
%         boundedline([1:NumYears]', mean(Richness_raw_sub)',[max(mean(Richness_raw_sub_mean-Richness_raw_sub_lo)',0),mean(Richness_raw_sub_up-Richness_raw_sub_mean)'],'--k','alpha','transparency', 0.1);
        boundedline([1:NumYears]', mean(Chao2_sub)',std(Chao2_sub)','--c','alpha','transparency', 0.1);
        boundedline([1:NumYears]', mean(Chao1_sub)',std(Chao1_sub)','-b','alpha','transparency', 0.1);
        boundedline([1:NumYears]', mean(Richness_omega_sub)',std(Richness_omega_sub)','-r','alpha','transparency', 0.1);
        boundedline([1:NumYears]', mean(Richness_raw_sub)',std(Richness_raw_sub)','--k','alpha','transparency', 0.1);
   
        plot(mean(Chao1_sub),'-b','LineWidth',2);
        plot(mean(Chao2_sub),'--','LineWidth',2,'Color',[0 0.7 0.7]);
        plot(mean(Richness_omega_sub),'-r','LineWidth',2);
        plot(mean(Richness_raw_sub),'--k','LineWidth',2);
        
        xlabel 'year'
        ylim([minRich maxRich])
        ylabel 'richness'
        title({[num2str(ksub) ' transects x ' num2str(ksub_sub) ' quadrats'];[num2str(m_perms(mperm)*100) '% individuals observed']},'fontweight','normal')
        xticks([1:NumYears])
        xticklabels(YearLabels)
        xlim([1,NumYears])
%         text(1.2,minRich+diffRich*.95,['\DeltaRaw=' num2str(mean(RichnessSlope_raw_boot(:)),1) '\pm' num2str(mean(RichnessSlopeSD_raw_boot(:)),1) ', p=' num2str(mean(RichnessSlopeP_raw_boot(:)),1) '(' num2str(mean(RichnessSlopeP_raw),1) ')\pm' num2str(std(RichnessSlopeP_raw_boot(:)),1) ', \alpha=' num2str(sum(RichnessSlopeP_raw_boot(:)<=0.05)/numBoot,1) '(' num2str(sum(RichnessSlopeP_raw(:)<=0.05)/numResample,1) ')' ],'fontsize',14);
%         text(1.2,minRich+diffRich*.85,['\DeltaChao1=' num2str(mean(RichnessSlope_Chao1_boot(:)),1) '\pm' num2str(mean(RichnessSlopeSD_Chao1_boot(:)),1) ', p=' num2str(mean(RichnessSlopeP_Chao1_boot(:)),1) '(' num2str(mean(RichnessSlopeP_Chao1),1) ')\pm' num2str(std(RichnessSlopeP_Chao1_boot(:)),1) ', \alpha=' num2str(sum(RichnessSlopeP_Chao1_boot(:)<=0.05)/numBoot,1) '(' num2str(sum(RichnessSlopeP_Chao1(:)<=0.05)/numResample,1) ')'],'fontsize',14,'color','b');
%         text(1.2,minRich+diffRich*.75,['\DeltaChao2=' num2str(mean(RichnessSlope_Chao2_boot(:)),1) '\pm' num2str(mean(RichnessSlopeSD_Chao2_boot(:)),1) ', p=' num2str(mean(RichnessSlopeP_Chao2_boot(:)),1) '(' num2str(mean(RichnessSlopeP_Chao2),1) ')\pm' num2str(std(RichnessSlopeP_Chao2_boot(:)),1) ', \alpha=' num2str(sum(RichnessSlopeP_Chao2_boot(:)<=0.05)/numBoot,1) '(' num2str(sum(RichnessSlopeP_Chao2(:)<=0.05)/numResample,1) ')'],'fontsize',14,'color',[0 0.7 0.7]);
%         text(1.2,minRich+diffRich*.63,['\Delta\Omega_T=' num2str(mean(RichnessSlope_Taylor_boot(:)),1) '\pm' num2str(mean(RichnessSlopeSD_Taylor_boot(:)),1) ', p=' num2str(mean(RichnessSlopeP_Taylor_boot(:)),1) '(' num2str(mean(RichnessSlopeP_Taylor),1) ')\pm' num2str(std(RichnessSlopeP_Taylor_boot(:)),1) ', \alpha=' num2str(sum(RichnessSlopeP_Taylor_boot(:)<=0.05)/numBoot,1) '(' num2str(sum(RichnessSlopeP_Taylor(:)<=0.05)/numResample,1) ')'],'fontsize',14,'Color','r');
%         text(1,2,minRich+diffRich*.95,['bootstrapped estimates:'],'fontsize',14)
%         text(1.2,minRich+diffRich*.85,['raw slope=' num2str(mean(RichnessSlope_raw_boot(:)),1) '\pm' num2str(mean(RichnessSlopeSD_raw_boot(:)),1) ', p=' num2str(mean(RichnessSlopeP_raw_boot(:)),1) '\pm' num2str(std(RichnessSlopeP_raw_boot(:)),1) ', \beta=' num2str(sum(RichnessSlopeP_raw_boot(:)>0.05)/numBoot,1)],'fontsize',14);
%         text(1.2,minRich+diffRich*.75,['Chao1 slope=' num2str(mean(RichnessSlope_Chao1_boot(:)),1) '\pm' num2str(mean(RichnessSlopeSD_Chao1_boot(:)),1) ', p=' num2str(mean(RichnessSlopeP_Chao1_boot(:)),1) '\pm' num2str(std(RichnessSlopeP_Chao1_boot(:)),1) ', \beta=' num2str(sum(RichnessSlopeP_Chao1_boot(:)>0.05)/numBoot,1)],'fontsize',14,'color','b');
%         text(1.2,minRich+diffRich*.63,['\Omega_o slope=' num2str(mean(RichnessSlope_Omega_boot(:)),1) '\pm' num2str(mean(RichnessSlopeSD_Omega_boot(:)),1) ', p=' num2str(mean(RichnessSlopeP_Omega_boot(:)),1) '\pm' num2str(std(RichnessSlopeP_Omega_boot(:)),1) ', \beta=' num2str(sum(RichnessSlopeP_Omega_boot(:)>0.05)/numBoot,1)],'fontsize',14,'Color','r');
%         text(1,2,minRich+diffRich*.45,['point estimates:'],'fontsize',14)
%         text(1.2,minRich+diffRich*.35,['raw slope=' num2str(mean(RichnessSlope_raw),1) '\pm' num2str(mean(RichnessSlopeSD_raw),1) ', p=' num2str(mean(RichnessSlopeP_raw),1) '\pm' num2str(std(RichnessSlopeP_raw),1) ', \beta=' num2str(sum(RichnessSlopeP_raw(:)>0.05)/numBoot,1)],'fontsize',14);
%         text(1.2,minRich+diffRich*.25,['Chao1 slope=' num2str(mean(RichnessSlope_Chao1),1) '\pm' num2str(mean(RichnessSlopeSD_Chao1),1) ', p=' num2str(mean(RichnessSlopeP_Chao1),1) '\pm' num2str(std(RichnessSlopeP_Chao1),1) ', \beta=' num2str(sum(RichnessSlopeP_Chao1(:)>0.05)/numBoot,1)],'fontsize',14,'color','b');
%         text(1.2,minRich+diffRich*.13,['\Omega_o slope=' num2str(mean(RichnessSlope_Omega),1) '\pm' num2str(mean(RichnessSlopeSD_Omega),1) ', p=' num2str(mean(RichnessSlopeP_Omega),1) '\pm' num2str(std(RichnessSlopeP_Omega),1) ', \beta=' num2str(sum(RichnessSlopeP_Omega(:)>0.05)/numBoot,1)],'fontsize',14,'Color','r');
        text(1-(NumYears-1)*0.15,maxRich*1.15,char(64+i),'Fontsize',16)
        
        %plot state variable means and variances over years (averaged over
        %replicates)
        subplot(length(m_perms)*2,length(ksub_perms)+1,(mperm-1)*(length(ksub_perms)+1)+testperm+1+6) %plot estimated richness by year
        hold on
        plot(meanStates(1,:),'-','LineWidth',2,'color','b'); %E[mn]
        plot(meanStates(2,:),'-','LineWidth',2,'color',[255 139 37]/255); %E[P]
        plot(meanStates(3,:),'--','LineWidth',2,'color','b'); %var[mn]
        plot(meanStates(4,:),'--','LineWidth',2,'color',[255 139 37]/255); %var[P]
        plot(meanStates(5,:),'--','LineWidth',2,'color',[69 181 80]/255); %cov[mn,P]
        refline(0,1)
        ax=gca;
        ax.YAxis.Scale= 'log';
        ylim([0.01 2*10^5])
        yticks([0.01 1 10^2 10^4])
        ylabel('state')
        ylims=ylim;
        xlabel 'year'
        xticks([1:NumYears])
        xticklabels(YearLabels)
        xlim([1,NumYears])
        title({[num2str(ksub) ' transects x ' num2str(ksub_sub) ' quadrats'];[num2str(m_perms(mperm)*100) '% individuals observed']},'fontweight','normal')
        text(1-(NumYears-1)*0.15,ylims(2)*3,char(64+i),'Fontsize',16)
        if i==6
            lgd=legend('E[$$\hat{mn}$$]','E[$$\hat{P}$$]','var[$$\hat{mn}$$]','var[$$\hat{P}$$]','cov[$$\hat{mn}$$,$$\hat{P}$$]','Interpreter','Latex');
            lgd.AutoUpdate='off';
        end
        i=i+1;
    end
end
figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/2 scrsz(4)/6],'DefaultAxesFontSize',16,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);

plotOrder=[1 2 3 4 5 6];
reordered_Chao1=mean(allSub_Chao1(plotOrder,:),2);
reordered_Chao2=mean(allSub_Chao2(plotOrder,:),2);
reordered_omega=mean(allSub_Richness_omega(plotOrder,:),2);
reordered_raw=mean(allSub_Richness_raw(plotOrder,:),2);
reordered_Chao1SD=mean(allSub_Chao1SD(plotOrder,:),2);
reordered_Chao2SD=mean(allSub_Chao2SD(plotOrder,:),2);
reordered_omegaSD=mean(allSub_Richness_omegaSD(plotOrder,:),2);
reordered_rawSD=mean(allSub_Richness_rawSD(plotOrder,:),2);
reordered_Chao1CV=mean(allSub_Chao1SD(plotOrder,:)./allSub_Chao1(plotOrder,:),2);
reordered_Chao2CV=mean(allSub_Chao2SD(plotOrder,:)./allSub_Chao2(plotOrder,:),2);
reordered_omegaCV=mean(allSub_Richness_omegaSD(plotOrder,:)./allSub_Richness_omega(plotOrder,:),2);
reordered_rawCV=mean(allSub_Richness_rawSD(plotOrder,:)./allSub_Richness_raw(plotOrder,:),2);

subplot(1,5,1)
hold on
plot(reordered_Chao1(1:3),'-b','LineWidth',5);
plot(reordered_Chao2(1:3),'--','LineWidth',5,'Color',[0 0.7 0.7]);
plot(reordered_omega(1:3),'-r','LineWidth',5);
plot(reordered_raw(1:3),'--k','LineWidth',5);

plot(reordered_Chao1(4:6),'-b','LineWidth',2);
plot(reordered_Chao2(4:6),'--','LineWidth',2,'Color',[0 0.7 0.7]);
plot(reordered_omega(4:6),'-r','LineWidth',2);
plot(reordered_raw(4:6),'--k','LineWidth',2);

ylabel('mean richness')
xlim([1 3]);
xticks([1:3]);
plotLabels=["9x10" "9x2" "2x9"];
xticklabels(plotLabels)

subplot(1,5,2)
hold on

plot(reordered_Chao1CV(1:3),'-b','LineWidth',5);
plot(reordered_Chao2CV(1:3),'--','LineWidth',5,'Color',[0 0.7 0.7]);
plot(reordered_omegaCV(1:3),'-r','LineWidth',5);
plot(reordered_rawCV(1:3),'--k','LineWidth',5);
plot(reordered_Chao1CV(4:6),'-b','LineWidth',2);
plot(reordered_Chao2CV(4:6),'--','LineWidth',2,'Color',[0 0.7 0.7]);
plot(reordered_omegaCV(4:6),'-r','LineWidth',2);
plot(reordered_rawCV(4:6),'--k','LineWidth',2);
ylabel('mean richness C.V.')

xlim([1 3]);
xticks([1:3]);
plotLabels=["9x10" "9x2" "2x9"];
xticklabels(plotLabels)

subplot(1,5,3)
hold on
plotOrder=[1 2 3 4 5 6];
reordered_Chao1=mean(allSub_Chao1(plotOrder,:),2);
reordered_Chao2=mean(allSub_Chao2(plotOrder,:),2);
reordered_Richness_omega=mean(allSub_Richness_omega(plotOrder,:),2);
reordered_Richness_raw=mean(allSub_Richness_raw(plotOrder,:),2);
plot(100*reordered_Chao1(1:3)./reordered_Chao1(1),'-b','LineWidth',5);
plot(100*reordered_Chao2(1:3)./reordered_Chao2(1),'--','LineWidth',5,'Color',[0 0.7 0.7]);
plot(100*reordered_Richness_omega(1:3)./reordered_Richness_omega(1),'-r','LineWidth',5);
plot(100*reordered_Richness_raw(1:3)./reordered_Richness_raw(1),'--k','LineWidth',5);

plot(100*reordered_Chao1(4:6)./reordered_Chao1(1),'-b','LineWidth',2);
plot(100*reordered_Chao2(4:6)./reordered_Chao2(1),'--','LineWidth',2,'Color',[0 0.7 0.7]);
plot(100*reordered_Richness_omega(4:6)./reordered_Richness_omega(1),'-r','LineWidth',2);
plot(100*reordered_Richness_raw(4:6)./reordered_Richness_raw(1),'--k','LineWidth',2);
ylabel('% of full estimates')
ylim([0 100])
xlim([1 3]);
xticks([1:3]);
plotLabels=["9x10" "9x2" "2x9"];
xticklabels(plotLabels)
xlabel('transects x quadrats sampled')

subplot(1,5,4)
hold on
reordered_Chao1_slope=mean(allSubSlope_Chao1(plotOrder,:),2);
reordered_Chao2_slope=mean(allSubSlope_Chao2(plotOrder,:),2);
reordered_omega_slope=mean(allSubSlope_omega(plotOrder,:),2);
reordered_raw_slope=mean(allSubSlope_raw(plotOrder,:),2);
reordered_Chao1_slopeSD=mean(allSubSlopeSD_Chao1(plotOrder,:),2);
reordered_Chao2_slopeSD=mean(allSubSlopeSD_Chao2(plotOrder,:),2);
reordered_omega_slopeSD=mean(allSubSlopeSD_omega(plotOrder,:),2);
reordered_raw_slopeSD=mean(allSubSlopeSD_raw(plotOrder,:),2);

plot(reordered_Chao1_slope(1:3),'-b','LineWidth',5);
plot(reordered_Chao2_slope(1:3),'--','LineWidth',5,'Color',[0 0.7 0.7]);
plot(reordered_omega_slope(1:3),'-r','LineWidth',5);
plot(reordered_raw_slope(1:3),'--k','LineWidth',5);

plot(reordered_Chao1_slope(4:6),'-b','LineWidth',2);
plot(reordered_Chao2_slope(4:6),'--','LineWidth',2,'Color',[0 0.7 0.7]);
plot(reordered_omega_slope(4:6),'-r','LineWidth',2);
plot(reordered_raw_slope(4:6),'--k','LineWidth',2);

ylabel('\Deltarichness/\DeltaT')
xticks([1:3]);
plotLabels=["9x10" "9x2" "2x9"];
xticklabels(plotLabels)

subplot(1,5,5)
hold on
reordered_Chao1_slopeP=mean(allSubSlopePPt_Chao1(plotOrder,:),2);
reordered_Chao2_slopeP=mean(allSubSlopePPt_Chao2(plotOrder,:),2);
reordered_omega_slopeP=mean(allSubSlopePPt_omega(plotOrder,:),2);
reordered_raw_slopeP=mean(allSubSlopePPt_raw(plotOrder,:),2);

plot(reordered_Chao1_slopeP(1:3),'-b','LineWidth',5);
plot(reordered_Chao2_slopeP(1:3),'--','LineWidth',5,'Color',[0 0.7 0.7]);
plot(reordered_omega_slopeP(1:3),'-r','LineWidth',5);
plot(reordered_raw_slopeP(1:3),'--k','LineWidth',5);

plot(reordered_Chao1_slopeP(4:6),'-b','LineWidth',2);
plot(reordered_Chao2_slopeP(4:6),'--','LineWidth',2,'Color',[0 0.7 0.7]);
plot(reordered_omega_slopeP(4:6),'-r','LineWidth',2);
plot(reordered_raw_slopeP(4:6),'--k','LineWidth',2);

ylabel('p_{\Deltarichness}')
xticks([1:3]);
plotLabels=["9x10" "9x2" "2x9"];
xticklabels(plotLabels)

