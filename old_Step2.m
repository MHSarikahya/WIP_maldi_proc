%%%% Match Veh/Treated across all brain regions %%%%%
%%%% By: Mohammed Sarikahya %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to match the matched m/z
%%% across all brain regions

%% How to use: Read the lines 16 to 36, and fill in the places with 'EDIT'.
% Then click hit the F5 button to Run the script, or click Run
% in the Editor tab above

function [] = old_Step2()
clc;
clear all;

% Set these Directories for the Script to work
ExcelDirectory = 'C:\Users\msari\OneDrive\Desktop\pd120 lipid selections\Male\Excels\'; %EDIT' here you send the processed excels
MetaboliteListDirectory = "C:\Users\msari\OneDrive\Desktop\!!!! MALDI SCRIPTS !!!!\codefromYeungLab\Metabolite List PUFA\"; %" Metabolite List - LIPIDS

% Set the brain region of interest, can signify Right/Left as well, if needed
cd(fullfile(ExcelDirectory)); % DO NOT EDIT
BrainRegion = dir('*prl*'); % EDIT; Make sure to maintain the '*x*', otherwise it will not wo

% Name of Output excel file:
name = 'PUFA__dd__';  %EDIT
% Set the name of the sheet where you will save the data to:
sheet = 'Sheet1'; 

%%How to calculate tolerance value
% x = peakmz_Treated;
% y = peakmz_Veh;
% tol = 0.00005;
%%within tolderance 
% b = tol*max(abs([x(:);y(:)]))
tol = 0.0009; % Lipid

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% DO NOT ALTER ANYTHING BELOW THIS LINE %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%

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
       peakInten_First= X(:,2);
       peak_auc_First = X(:,7);
       
       First_combined = [peakmz_First, peakInten_First, peak_auc_First];

     cd(fullfile(ExcelDIR));


    peakmz_First = peakmz_First.';
%     peakmz_Compared = peakmz_Compared.';
        
        
    f = 0;
    h = 0;
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MetaCompare = dir(strcat(MetaboliteListDirectory));
    MetaCompare = {MetaCompare.name};
    MetaCompare = MetaCompare(~ismember(MetaCompare, {'.', '..', '.DS_Store'}));
    cd(MetaboliteListDirectory)
    [numeric, ~, raw] = xlsread(cell2mat(MetaCompare));
    cd(ExcelDIR);
    Metabolites =  numeric(:,:);
%     Metabolites(1,:) = [];
    Metabolites(:,1) = [];
    Metabolites = Metabolites';
    MetaNames = raw(:,:);
    MetaNames(1,:) = [];
    MetaNames(:,1) = [];

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   %Matching treated values to vehicle values
%       for i = 1:length(Metabolites)
%          for k = 1:length(peakmz_First)
%               x = Metabolites(1,i);
%               y = peakmz_First(1,k);
%               if ismembertol(y, x, tol, 'ByRows', true)
%                  f = Metabolites(i);
%                  h = peakmz_First(k);
%               end
%             newidentities(1,i) = h;
%             newidentities(2,i) = f;
%          end
%       end
      
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
%         MatchedPeaks = uniquetol(MatchedPeaks,'ByRows',true); 

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

            
%      for x = 1
%         for k = 1:length(MatchedIden)
%           for bb = 1:length(Metabolites)
%                     t = MetaNames(bb,1);
%                if char(MatchedIden(k,1)) ==  char(MetaNames(bb,2))
%                    MatchedIdenNumb(k) = t;
%                end
%           end
%         end
%     end
    
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
    Matched_Peak_Inten_Veh = [MatchedPeaks(:,1), MatchedInten_V]; 

    
    %Matching the peak to the intensity for Vehicle
    for x = 1
        for k = 1:length(MatchedPeaks)
          for i = 1:length(peakmz_First)
                    t = First_combined(i,3);
               if MatchedPeaks(k,1) ==  First_combined(i,1) 
                   MatchedAUC_V(k) = t;
               end
          end
        end
    end

    MatchedAUC_V = MatchedAUC_V.'; %'
    Matched_Peak_Inten_Veh = [Matched_Peak_Inten_Veh, MatchedAUC_V]; 

    
    Metabolites = Metabolites.';
    
    
    
    Matched =  [];
    
    
    Brain_File.Matched_Peak_Inten_Veh = Matched_Peak_Inten_Veh(:,:);
%     MatchedIden = MatchedIden.';
    Brain_File.MetaboliteList = MatchedIden;
%     Brain_File.Matched_Peak_Inten_Treated = Matched_Peak_Inten_Treated(:,:);
    MatchedIden = MatchedIden.';

    Brain_Files = struct2cell(Brain_File);
    Brain_Files = (Brain_Files');
    BrainFiles(tt,:) = (Brain_Files(1,:));
   
%     MatchedFiles(1:1,:) = []
    
% %     Brain_Files = [Brain_Files{1,1};Brain_Files{1,2}];
% %     BrainFiles(tt) = mat2cell(Brain_Files)
% %     
% % % for i = 1:size(subjects,2)
% % %     a.curSubj = subjects{i};  %structure 'a', filed with all the files from the first line of the loop
% % %     a.dataDir = strcat(DATADIR,a.curSubj,'/');
% % %     a.files = files(:,1);               %where the 4D files are
% % %     a.SPM = SPM_files{1,i};
% %  
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

  %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Save the values %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     MatchedFiles = BrainFiles(:,:)';
%      MatchedFiles(:,2) = BrainFiles(:,6)';
     %MatchedFiles(1:1,:) = [];
     
 xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles(2,1),sheet,'C1');     %Write column header
     
        col_header1 = (tmpfile(1,1));
          col_header2 = (tmpfile(1,2));
            col_header3 = (tmpfile(1,3));
             col_header4 = (tmpfile(1,4));
                col_header5 = (tmpfile(1,5));
                 col_header6 = (tmpfile(1,6));
                 col_header7 = (tmpfile(1,7));
                 col_header8 = (tmpfile(1,8));
                 col_header9 = (tmpfile(1,9));
                 col_header10 = (tmpfile(1,10));
%                              col_header11 = (tmpfile(1,11));

%% 1   
 xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{2,1}(:,1), sheet, 'A2');
   xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{2,1}(:,2), sheet, 'B2');   
      xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), col_header1,sheet,'C1');     
        xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{1,1}(:,1), sheet, 'C2');
           xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{1,1}(:,2), sheet, 'D2');
               xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{1,1}(:,3), sheet, 'E2');

 %  2
 xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{2,2}(:,1), sheet, 'G2');
   xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{2,2}(:,2), sheet, 'H2');   
      xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), col_header2,sheet,'I1');     
        xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{1,2}(:,1), sheet, 'I2');
           xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{1,2}(:,2), sheet, 'J2');
               xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{1,2}(:,3), sheet, 'K2');
                                    
% 
% 
%  3
 xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{2,3}(:,1), sheet, 'M2');
   xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{2,3}(:,2), sheet, 'N2');   
      xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), col_header3,sheet,'O1');     
        xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{1,3}(:,1), sheet, 'O2');
           xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{1,3}(:,2), sheet, 'P2');
               xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{1,3}(:,3), sheet, 'Q2');

% 4
xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{2,4}(:,1), sheet, 'S2');
  xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{2,4}(:,2), sheet, 'T2');   
     xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), col_header4,sheet,'U1');     
       xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{1,4}(:,1), sheet, 'U2');
          xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{1,4}(:,2), sheet, 'V2');
              xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{1,4}(:,3), sheet, 'W2') 
 
%  5
 xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{2,5}(:,1), sheet, 'Y2');
   xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{2,5}(:,2), sheet, 'Z2');   
      xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), col_header5,sheet,'AA1');     
        xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{1,5}(:,1), sheet, 'AA2');
           xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{1,5}(:,2), sheet, 'AB2');
               xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{1,5}(:,3), sheet, 'AC2');
%  6
 xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{2,6}(:,1), sheet, 'AE2');
   xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{2,6}(:,2), sheet, 'AF2');   
      xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), col_header6,sheet,'AG1');     
        xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{1,6}(:,1), sheet, 'AG2');
           xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{1,6}(:,2), sheet, 'AH2');
               xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{1,6}(:,3), sheet, 'AI2');
%  7              
 xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{2,7}(:,1), sheet, 'AK2');
   xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{2,7}(:,2), sheet, 'AL2');   
      xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), col_header7,sheet,'AM1');     
        xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{1,7}(:,1), sheet, 'AN2');
           xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{1,7}(:,2), sheet, 'AO2');
               xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{1,7}(:,3), sheet, 'AP2');
%  8             
  xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{2,8}(:,1), sheet, 'AR2');
   xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{2,8}(:,2), sheet, 'AS2');   
      xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), col_header8,sheet,'AT1');     
        xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{1,8}(:,1), sheet, 'AU2');
           xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{1,8}(:,2), sheet, 'AV2');
               xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{1,8}(:,3), sheet, 'AW2');
% 9              
 xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{2,9}(:,1), sheet, 'AY2');
   xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{2,9}(:,2), sheet, 'AZ2');   
      xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), col_header9,sheet,'BA1');     
        xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{1,9}(:,1), sheet, 'BB2');
           xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{1,9}(:,2), sheet, 'BC2');
               xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{1,9}(:,3), sheet, 'BD2'); 
% 10                
 xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{2,10}(:,1), sheet, 'BF2');
   xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{2,10}(:,2), sheet, 'BG2');   
      xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), col_header10,sheet,'BH1');     
        xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{1,10}(:,1), sheet, 'BI2');
           xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{1,10}(:,2), sheet, 'BJ2');
               xlswrite(strcat(ExcelDIR, name, 'AllExpMatch'), MatchedFiles{1,10}(:,3), sheet, 'BK2'); 


% %     save(strcat(ExcelDIR, name),'BrainFiles');
    save(strcat(ExcelDIR,name,'_MetaMatch'),'MatchedFiles');

      cd(fullfile(ExcelDirectory))

%     xlswrite(strcat(ExcelDIR, name, '_ANALYZED'), Matched_ratio, sheet, 'M2');

end