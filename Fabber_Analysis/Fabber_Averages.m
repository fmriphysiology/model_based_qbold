function Fabber_Averages(data_dir)
% Fabber_Averages.m usage:
%
%       Fabber_Averages(data_dir)
%
% Loops through a set of FABBER analysis results (contained in DATA_DIR) and
% calculates, and plots, the average grey-matter parameter values for R2p and
% DBV. Requires data in the BIDS format. Requires the following scripts:
%       sigstar.m 
%       read_avw.m
%
% This script is designed to be used with the following dataset:
%       Cherukara MT, Stone AJ, Chappell MA, Blockley NP. Data acquired to
%       demonstrate model-based Bayesian inference of brain oxygenation using
%       quantitative BOLD, Oxford University Research Archive 2018. doi: <Please
%       see ORA entry for DOI> 
%
% 
%       Copyright (C) University of Oxford, 2018
%
% 
% Created by MT Cherukara, 12 March 2018
%
% CHANGELOG:
%
% 2018-08-27. Various changes.


%% Initial
close all;

% Hardcoded parameters
nsubs = 7;          % number of subjects
slicenum = 4:9;     % slices that we want
thr_R2p = 20;       % R2' threshold value
thr_DBV = 1;        % DBV threshold value

% Task names
tasks = {'qbold-linear','qbold-1C','qbold-2C'};
ntasks = length(tasks);


%% Load Data and Calculate Averages

% Pre-allocate data arrays
av_R2p = zeros(nsubs,ntasks);
sd_R2p = zeros(nsubs,ntasks);
av_DBV = zeros(nsubs,ntasks);
sd_DBV = zeros(nsubs,ntasks);


% loop through subjects
for ss = 1:nsubs
    
     % load grey matter mask
    mask_gm = LoadSlice(strcat(data_dir,'sub-0',num2str(ss),'/pve/sub-0',...
                        num2str(ss),'_mask_greymatter.nii.gz'), slicenum);
                    
	% loop through tasks
    for tt = 1:ntasks
        
        task_name = tasks{tt};
    
        % directory and beginning of the subject name
        sub_dir = strcat(data_dir,'sub-0',num2str(ss),'/param/task-',task_name, ...
                         '/sub-0',num2str(ss),'_',task_name,'_');

        % load data
        dat_R2p = LoadSlice(strcat(sub_dir,'param-R2p_abs.nii.gz'), slicenum);
        dat_DBV = LoadSlice(strcat(sub_dir,'param-DBV_abs.nii.gz'), slicenum);
        std_R2p = LoadSlice(strcat(sub_dir,'std-R2p.nii.gz'), slicenum);
        std_DBV = LoadSlice(strcat(sub_dir,'std-DBV.nii.gz'), slicenum);


        % mask and vectorize
        dat_R2p = dat_R2p(:).*mask_gm(:);
        dat_DBV = dat_DBV(:).*mask_gm(:);
        std_R2p = std_R2p(:).*mask_gm(:);
        std_DBV = std_DBV(:).*mask_gm(:);

        % create a mask of values to remove
        %    we remove, voxels that are zero (masked), voxels which have non-finite
        %    values, voxels whose values are greater than the threshold, voxels,
        %    voxels whose standard deviation is greater than the threshold, or less
        %    than a very small fraction (0.1%) of the threshold

        bad_R2p =  (dat_R2p == 0) + ~isfinite(dat_R2p) + ~isfinite(std_R2p) + ...
                   (dat_R2p > thr_R2p) + (std_R2p > thr_R2p) + ...
                   (std_R2p < (thr_R2p.*1e-4));

        bad_DBV =  (dat_DBV == 0) + ~isfinite(dat_DBV) + ~isfinite(std_DBV) + ...
                   (dat_DBV > thr_DBV) + (std_DBV > thr_R2p) + ...
                   (std_DBV < (thr_DBV.*1e-4));      

        % Remove the bad values
        dat_R2p(bad_R2p ~= 0) = [];
        std_R2p(bad_R2p ~= 0) = [];
        dat_DBV(bad_DBV ~= 0) = [];
        std_DBV(bad_DBV ~= 0) = [];

        % Calculate and store means and standard deviations
        av_R2p(ss,tt) = mean(dat_R2p);
        sd_R2p(ss,tt) = mean(std_R2p);
        av_DBV(ss,tt) = mean(dat_DBV).*100; % convert DBV values to percentage
        sd_DBV(ss,tt) = mean(std_DBV).*100; 
        
    end % for tt = 1:tasks
    
end % for ss = 1:nsubs

%% Save data
save('fabber_averages_data.mat','av_R2p','av_DBV','sd_R2p','sd_DBV');


%% Plot Histograms

% Plot R2p
figure(1); hold on; box on;
bar(1:ntasks,mean(av_R2p),0.6);
errorbar(1:ntasks,mean(av_R2p),mean(sd_R2p),'k.','LineWidth',2,'MarkerSize',1);

axis([0.5,ntasks+0.5,0,5.6]);
ylabel('R_2'' (s^-^1)');
xticks(1:ntasks);
xticklabels(tasks); 

% Plot DBV
figure(2); hold on; box on;
bar(1:ntasks,mean(av_DBV),0.6);
errorbar(1:ntasks,mean(av_DBV),mean(sd_DBV),'k.','LineWidth',2,'MarkerSize',1);

axis([0.5,ntasks+0.5,0,12.8]);
ylabel('DBV (%)');
xticks(1:ntasks);
xticklabels(tasks); 


%% Perform Statistical Analyses

% Perform ANOVA
[~,~,stat_R2p] = anova2(av_R2p,1,'off');
[~,~,stat_DBV] = anova2(av_DBV,1,'off');

% Multiple comparison
c_R2p = multcompare(stat_R2p,'display','off');
c_DBV = multcompare(stat_DBV,'display','off');


% plot significance stars
figure(1);
HR = sigstar({[1,2],[1,3],[2,3]},c_R2p(:,6),1);
set(HR,'Color','k');

figure(2);
HD = sigstar({[1,2],[1,3],[2,3]},c_DBV(:,6),1);
set(HD,'Color','k');

end % function


%% LoadSlice function
function slicedata = LoadSlice(filename,slicenum)
    % Loads the data from a specified NIFTI file and cuts it down to particular
    % slices

    % Load the selected NIFTY into dataset
    dataset = read_avw(filename);

    % for now
    slicedata = squeeze(dataset(:,:,slicenum,:));
    
end
