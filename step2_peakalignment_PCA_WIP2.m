%%% Step2 PeakAlignment&PCA %%%
%%% By: Mohammed Sarikahya %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dependencies: (Addon within Matlab)
%
% Signal Processing Toolbox, Signal Integrity Toolbox, PADCAT
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% THIS IS A WIP VERSION 09.20.2023 %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% DO NOT USE AS IS, EXPERIMENTAL %%%%%%%%%
%%%%%%% DO NOT USE AS IS, EXPERIMENTAL %%%%%%%%%
%%%%%%% DO NOT USE AS IS, EXPERIMENTAL %%%%%%%%%
%%%%%%% DO NOT USE AS IS, EXPERIMENTAL %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear ;


% Set these Directories for the Script to work
ExcelDirectory = 'I:\Omega-3 + PCE\Lipids\PD45\Region Selections\FEMALE\New\Excels\'; %EDIT' here you send the processed excels
AnalysisDirectory = "I:\Omega-3 + PCE\Lipids\PD45\Region Selections\FEMALE\New\Excels\Analysis2\"; %EDIT' dir for analysis

cd(fullfile(ExcelDirectory))
mkdir Analysis2

% Set Regions of interest, sex, and whether FA, Lipids, or etc.
Sex = 'Female'; % MUST add "_" before Males, else, will be also get Females
ROI = 'DS';

% PRL IL VS DS CA1 core shell
%% Turn on list as needed
% MetaboliteListDirectory = "C:\Users\msari\OneDrive\Desktop\!!!! MALDI SCRIPTS !!!!\codeforYeungLab\Metabolite List PUFA\"; %" Metabolite List - PUFA
%  Type = 'PUFA';
%  tol = 0.001; % PUFA

MetaboliteListDirectory = "C:\Users\msari\OneDrive\Desktop\!!!! MALDI SCRIPTS !!!!\codeforYeungLab\Metabolite List Lipids\"; %" Metabolite List - LIPIDS
Type = 'Lipids';
tol = 0.0008; % Lipid

% MetaboliteListDirectory = "C:\Users\msari\OneDrive\Desktop\!!!! MALDI SCRIPTS !!!!\codeforYeungLab\Metabolite List ZnO\"; %" Metabolite List - Metabolites (ZnO specific)
% Type = 'ZnO';
% tol = 0.0018; % ZnO

%%%How to calculate tolerance value; can only run iff data has been loaded in
%% x = peakmz_Treated;
%% y = peakmz_Veh;
%% tol = 0.00005;
%%%within tolderance 
%% b = tol*max(abs([x(:);y(:)]))


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% DO NOT ALTER ANYTHING BELOW THIS LINE %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NameofFILE = strcat(Sex,'_', ROI,'peakAli+Matched');
NameofFILEx = strcat(Sex,'_', ROI,'_');

cd(fullfile(ExcelDirectory)); % DO NOT EDIT
xx = (strcat('*', Sex,'*', ROI,'*'));
BrainRegion = dir(xx);

%%% can run tolerance value check for Matrix after running up to line 60 %%%

sheet = 'Sheet1'; 

ExcelDIR = ExcelDirectory;

ps = dir(strcat(ExcelDIR));
ps = {ps.name};
ps = ps(~ismember(ps, {'.', '..', '.DS_Store'}));

for i = 1:size(ps,2)
        cd(fullfile(ExcelDIR));
        tmpfile = BrainRegion; 
        tmpfile = {tmpfile.name};
        tmpfile = tmpfile(~ismember(tmpfile,{'.','..','.DS_Store'}));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMMENT EVERYTHING BELOW THIS LINE IF RUNNING .mat FILE ONLY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scanIdx=[];
AllPeaks=[];
h=cell.empty;

for tt = 1:size(tmpfile,2)
    Firstfile = char(tmpfile(1,tt));  
    match = '.xls';
    file1st = erase(Firstfile,match);
    filename = strcat(ExcelDIR, file1st);
    [numeric, text, raw] = xlsread(filename); 
    X =  numeric(:,1:7);
%     X(:,7) = []; 
    
    if X(1,1) == 0
          X(1,:) = [];
    end
    
    X(~any(~isnan(X), 2),:)=[];
        
       peakmz_First = X(:,1);
       
       peak_auc_First = X(:,7);
       
       First_combined = [{peakmz_First}, {peak_auc_First}];


%     mz =  (First_combined(:,1));
%     AUC = First_combined(:,2);
% 
%     [Peaklist, PFWHH] = mspeaks(cell2mat(mz), cell2mat(AUC))
%     
% h{tt,1} = Peaklist(:,:)

 mz =  cell2mat(First_combined(:,1));
    AUC = cell2mat(First_combined(:,2));

h{tt,1} = [mz AUC]

end

    peaks=cellfun(@(x) x(x(:,1)>1 ,:), h, 'uniformoutput',false);

    MetaCompare = dir(strcat(MetaboliteListDirectory));
            MetaCompare = {MetaCompare.name};
            MetaCompare = MetaCompare(~ismember(MetaCompare, {'.', '..', '.DS_Store'}));
            cd(MetaboliteListDirectory)
            [numeric, text, raw] = xlsread(cell2mat(MetaCompare));
            cd(ExcelDIR);
            Metabolites =  numeric(:,:);
            Metabolites(:,1) = [];
            Metabolites = Metabolites';
            MetaNames = raw(:,:);
            MetaNames(1,:) = [];
            MetaNames(:,1) = [];
f=0;
h=0;
fs=0;

for x = 1
for i=1:length(peaks)         
     [max_size, max_index] = max(cellfun('size', peaks, 1));
                maxsizepeak = peaks{max_index};
     [xa,newPeakMz] = alignsignals(peaks{max_index}(:,1)',peaks{i}(:,1)')
     
%      for m = 1:length(newPeakMz)
%          for k = 1:length(peaks{i}(:,1))
%               x = newPeakMz(1,m)';
%               y = peaks{i}(k,1);
%               if ismembertol(x,y, tol, 'ByRows', true)
%                  f = newPeakMz(m);
%                  h = peaks{i}(k);
%               end
%                newidentities = padcat(f,h);
%          end
%       end
%     MatchedPeaks = uniquetol(newidentities,'ByRows',true); 
%     Metabolites = Metabolites';
    MatchedIden = 0;

    %Matching the peak to the intensity for treated
     for x = 1
        for k = 1:length(newPeakMz)
          for bb = 1:length(peaks{i}(:,1))
                    t = peaks{i}(bb,2);
               if newPeakMz(1,k) ==  peaks{i}(bb,1) 
                   MatchedIden(k,:) = t;
               end
          end
        end
     end

    MatchedIden = MatchedIden.'; %'
    unAlignedPeakswAUC{i} = [newPeakMz', MatchedIden']; 
     
 save(strcat(ExcelDIR, NameofFILE,Type,'.mat'),'unAlignedPeakswAUC', '-v7.3');
    save(strcat(ExcelDIR, NameofFILE,Type,'.mat'),'tmpfile', '-append');

    
    movefile((strcat(ExcelDIR, NameofFILE,Type,'.mat')),AnalysisDirectory)

end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMMENT EVERYTHING ABOVE THIS LINE IF RUNNING .mat FILE ONLY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ExcelDIR = ExcelDirectory;

cd(fullfile(ExcelDirectory)); % DO NOT EDIT
NameofFILExx = dir(strcat('*', NameofFILE ,'*'));

%%%%% uncomment next line if running .mat only, otherwise, nothing will load in, adjust file name as needed%%%%%%%%
% load(strcat(AnalysisDirectory,'Male_NASHpeakAli+MatchedLipids.mat')) % change ' ' to name of .mat you want to load


for i = 1:size(ps,2)
    cd(fullfile(ExcelDIR));
    tmpfile = NameofFILExx; 
    tmpfile = {tmpfile.name};
    tmpfile = tmpfile(~ismember(tmpfile,{'.','..','.DS_Store'}));
end
      

for tt = 1:size(unAlignedPeakswAUC,2)      
       peakmz_First = unAlignedPeakswAUC{tt}(:,1);
       peak_auc_First= unAlignedPeakswAUC{tt}(:,2);       
       First_combined = [peakmz_First, peak_auc_First];

     cd(fullfile(ExcelDIR));


    peakmz_First = peakmz_First.';        
        
    f = 0;
    h = 0;
    
    MetaCompare = dir(strcat(MetaboliteListDirectory));
    MetaCompare = {MetaCompare.name};
    MetaCompare = MetaCompare(~ismember(MetaCompare, {'.', '..', '.DS_Store'}));
    cd(MetaboliteListDirectory)
    [numeric, text, raw] = xlsread(cell2mat(MetaCompare));
    cd(ExcelDIR);
    Metabolites =  numeric(:,:);
    Metabolites(:,1) = [];
    Metabolites = Metabolites';
    MetaNames = raw(:,:);
    MetaNames(1,:) = [];
    MetaNames(:,1) = [];

     for i = 1:length(peakmz_First)
         for k = 1:length(Metabolites)
              x = peakmz_First(1,i);
              y = Metabolites(1,k);
              if ismembertol(x,y, tol, 'ByRows', true)
                 f = peakmz_First(i);
                 h = Metabolites(k);
              end
            newidentities(1,i) = f;
            newidentities(2,i) = h;
         end
      end

    newidentities = newidentities.'; %'
    if newidentities(1,1) == 0
        newidentities(1,:) = [];
    end
    if newidentities(1,1) == 0
        newidentities(1,:) = [];
    end
    if newidentities(1,1) == 0
        newidentities(1,:) = [];
    end
    if newidentities(1,1) == 0
        newidentities(1,:) = [];
    end
    if newidentities(1,1) == 0
        newidentities(1,:) = [];
    end
    
    MatchedPeaks = uniquetol(newidentities,'ByRows',true); 

    Metabolites = Metabolites';
    MatchedIden = cell.empty;

    %Matching the peak to the intensity for treated
     for x = 1
        for k = 1:length(MatchedPeaks)
          for bb = 1:length(Metabolites)
                    t = MetaNames(bb,:);
               if MatchedPeaks(k,2) ==  Metabolites(bb,1) 
                   MatchedIden(k,:) = t;
               end
          end
        end
     end
     
     if cell2mat(MatchedIden(1,1)) == []
           MatchedIden(1,:) = [];
     end

        MatchedIden = MatchedIden.'; %'
 
   MatchedIden = MatchedIden.'; %'


    %Matching the peak to the intensity for Vehicle
    for x = 1
        for k = 1:length(MatchedPeaks)
          for i = 1:length(peakmz_First)
                    t = First_combined(i,2);
               if MatchedPeaks(k,1) ==  First_combined(i,1) 
                   MatchedInten_V(k) = t;
               end
          end
        end
    end

    MatchedInten_V = MatchedInten_V.'; %'
    Matched_Peak_Inten_Veh = [MatchedPeaks(:,2), MatchedPeaks(:,1), MatchedInten_V]; 
    
    Metabolites = Metabolites.';
    
    
    
    Matched =  [];
    
    
    Brain_File.Matched_Peak_Inten_Veh = Matched_Peak_Inten_Veh(:,:);
    Brain_File.MetaboliteList = MatchedIden;
    MatchedIden = MatchedIden.';

    Brain_Files = struct2cell(Brain_File);
    Brain_Files = (Brain_Files');
    BrainFiles(tt,:) = (Brain_Files(1,:));
   

    mz_comp = [];
    mz_Treated = [];
    Matched_ratio = [];
    MatchedPeaks = [];
    MatchedInten_V  = [];
    Matched_Peak_Inten_Treated = [];
    Matched_Peak_Inten_Veh = [];
    peakInten_Treated = [];
    mz_Treated = [];
    peakInten_Veh = [];
    mz_comp = [];
    newidentities = [];
    Brain_File = [];
    MatchedAUC = [];
    MatchedAUC_V = [];
       
        peakmz_First_T = [];
        peakInten_First_T = [];
        peak_auc_First_T = [];
        IntenRatio_F = [];
        AUCratio_F = [];
        peakmz_Compared_T = [];
        peakInten_Compared_T = [];
        peak_auc_Compared_T = [];
        IntenRatio_C = [];
        AUCratio_C = [];
        
    First_combined = [];
    Treated_combined = [];
    
%     MatchedIden = [];
    mz_comp = [];
    MetaNames = [];
    Metabolites = [];
    MatchedIden = "";
    MatchedIden = [];


end

%% Save the values

     MatchedFiles = BrainFiles(:,:)';

for i = 1:size(ps,2)
        cd(fullfile(ExcelDIR));
        tmpfile = BrainRegion; 
        tmpfile = {tmpfile.name};
        tmpfile = tmpfile(~ismember(tmpfile,{'.','..','.DS_Store'}));
end



for x = 1
        for mmmmmmmm = 1:size(BrainRegion',2)
             theor = (MatchedFiles{2,mmmmmmmm});
       theor2= (cell2mat(theor(:,1)));
        theorr = uniquetol(theor2,'ByRows',true); 

                       exp = (MatchedFiles{1,mmmmmmmm});
                              expp =  uniquetol(exp,'ByRows',true); 
         T = table(expp);
            writetable(T,strcat(NameofFILE,'_',Type,'.xls'),'WriteMode','Append') 
        end
end


  theorMZ = 'theorMZ';
    theorName = 'theorName';
     expMZ = 'expMZ';
      expAUC = 'expAUC';

 theorMZ = repmat( {theorMZ}, 1, 1)
 theorName = repmat( {theorName}, 1, 1)
 expMZ = repmat( {expMZ}, 1, 1)
 expAUC = repmat( {expAUC}, 1, 1)

sheet = "sheet1";

xlswrite(strcat(ExcelDIR,NameofFILE,'_',Type,'.xls'), tmpfile','Filenames');
xlswrite(strcat(ExcelDIR, NameofFILE,'_',Type,'.xls'), theorMZ,sheet,'A1');     
xlswrite(strcat(ExcelDIR, NameofFILE,'_',Type,'.xls'), expMZ,sheet,'B1');     
xlswrite(strcat(ExcelDIR, NameofFILE,'_',Type,'.xls'), expAUC,sheet,'C1');     
      


movefile((strcat(ExcelDIR, NameofFILE,'_',Type,'.xls')),AnalysisDirectory)



% %%%%%%%%%%%%%%%
% %% PCA CODE  %%
% %%%%%%%%%%%%%%%
% 
folder = ExcelDirectory;

    cd(ExcelDirectory)
    fileXXname = dir("*.xls*");
    fullFileName = fileXXname.name;
    

    msi_table=readtable(fullFileName);
    msi_table=msi_table{:,:};
    
msi_mean = []
    %     msi_mean = msi_table - mean(msi_table);
     
     for ii = 1:length(unAlignedPeakswAUC)
         SelectAUCs(ii,:) =  unAlignedPeakswAUC{1, ii}(:,2) 
         NewAUCs =  vertcat(SelectAUCs, SelectAUCs);
     end
        msi_mean = NewAUCs - mean(NewAUCs);

    [coeff,score,latent,~,explained] = pca(msi_mean);
    [PCALoadings,PCAScores,PCAVar] = pca(msi_mean);
    

    %Explained Variance by Data 
    bar(explained);
    title('Explained Variance: More than X% explained by first two principal components');
    ylabel('PC');
    
    %Loading Plot; load PCALoadings .mat and then run these two lines to
    %revisualize the PCA
    PCA_plot=plot(PCALoadings); %
    mapcaplot(PCALoadings) %PCA visualization tool


save(strcat(NameofFILEx,'PCA_loadings'),'PCALoadings')
movefile((strcat(ExcelDIR, NameofFILEx,'PCA_loadings','.mat')),AnalysisDirectory)
  cd(AnalysisDirectory)
