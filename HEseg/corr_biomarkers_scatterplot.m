%% summary of biomarker
%% beta catenin vs. COX2
data_dir = 'D:\Documents\multiplex';
coordinate_dir = 'D:\Documents\multiplex\coordinates_all_bm';
image_dir = 'D:\Documents\colon_cancer_data\H&E_virtual';
seg_dir = 'D:\Documents\multiplex\seg_output_50';

%filelist = dir(fullfile(coordinate_dir,'*.mat'));
%filelist = {filelist.name}';

% read in the clincial data file
clinical_data = readtable(fullfile('D:\Documents\multiplex','clinical_data_all_spots.csv'),...
    'Delimiter',',');
filelist = clinical_data.spot_name;

% find the index of SMA
tmp = load(fullfile(coordinate_dir,filelist{1}));
bm_names = cellstr(tmp.bm_names);

% corr_bm_names = {'Median.Cell.BetaCatenin','Median.Cell.COX2'};
% indx = find(ismember(bm_names,corr_bm_names));
% corr_bm_names = {'BetaCatenin','COX2'};
% output_dir = 'D:\Documents\multiplex\pairwise_corr\BetaCat_COX2_2';

corr_bm_names = {'Median.Cell.BetaCatenin','Median.Cell.pMAPKAPK2'};
indx = find(ismember(bm_names,corr_bm_names));
corr_bm_names = {'BetaCat','MAPK'};
output_dir = 'D:\Documents\multiplex\pairwise_corr\BetaCat_MAPK_2';

if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

%dist_vec = 0:15:120;
dist_vec = 0:7.5:60;
%{
pair_bm_vec = cell(length(filelist),length(dist_vec)-1);

tic;
for i = 1:length(filelist)
    imname = filelist{i};
    
    if ~exist(fullfile(coordinate_dir, [imname '.mat']),'file') || ~exist(fullfile(seg_dir, [imname '.jpg']),'file')
        fprintf('Do not have this file: %s\n',imname);
        continue;
    end
    
    im = imread(fullfile(seg_dir,[imname '.jpg']));
  
    tmp = load(fullfile(coordinate_dir, [imname '.mat']));
    x = double(max(1,tmp.x));
    y = double(max(1,tmp.y));
    
    bm_data = tmp.bm_data(:,indx);
    pairwise_dist = squareform(pdist([x;y]'));
    
    for j = 1:(length(dist_vec)-1)
       [r,c] = find(pairwise_dist > dist_vec(j) & ...
           pairwise_dist <= dist_vec(j+1));
       ind0 = r<c;
       r = r(ind0); c = c(ind0);
       pair_bm_vec{i,j} = cat(2,bm_data(r,:),bm_data(c,:));
    end
end
toc;

save(fullfile(output_dir,'pair_bm_vec.mat'),'pair_bm_vec');
%}
% read in the clincial data file
% group data by stage
%load(fullfile(output_dir,'pair_bm_vec.mat'));

%{
for stage = 1:3
   indx_stage = clinical_data.stages == stage;
   %indx_stage = clinical_data.grades == stage;
   %indx_stage = clinical_data.recurrent_5yr == (stage - 1);
   for j  = 1:(length(dist_vec) - 1)
       bm_data = pair_bm_vec(indx_stage,j);
       bm_data = cat(1, bm_data{:});
       %sample_bm_data = datasample(bm_data,1000);
       sample_bm_data = bm_data;
       [~,indx1] = sort(sample_bm_data(:,1));
       [~,indx2] = sort(sample_bm_data(:,2));
       
       mld1 = fitlm(bm_data(:,1),bm_data(:,3)); % 1 v 1
       [ypred1,yci1] = predict(mld1,sample_bm_data(:,1), 'Prediction','observation');
       
       mld2 = fitlm(bm_data(:,1),bm_data(:,4)); % 1 v 2
       [ypred2,yci2] = predict(mld2,sample_bm_data(:,1), 'Prediction','observation');
       
       mld3 = fitlm(bm_data(:,2),bm_data(:,3)); % 2 v 1
       [ypred3,yci3] = predict(mld3,sample_bm_data(:,2), 'Prediction','observation');
       
       mld4 = fitlm(bm_data(:,2),bm_data(:,4)); % 2 v 2
       [ypred4,yci4] = predict(mld4,sample_bm_data(:,2), 'Prediction','observation');
       
       figure('position',[100 100 800 600]); 
       ha = tight_subplot(2,2,[.03 .03],[.08 .01],[.07 .01]);
     
       %subplot(2,2,1)
       axes(ha(1));
       plot(sample_bm_data(:,1), sample_bm_data(:,3),'.');
       hold on;
       plot(sample_bm_data(:,1), ypred1,'-','LineWidth',3);
       plot(sample_bm_data(indx1,1), yci1(indx1,1),'r--');
       plot(sample_bm_data(indx1,1), yci1(indx1,2),'r--');
       hold off;
       xlim([-1 15]); ylim([-1 15]);
       ylabel(corr_bm_names{1});
       formula = sprintf('y = %.2f + %.2f x\nAdj. Rsq = %.2f', ...
           mld1.Coefficients.Estimate(1),mld1.Coefficients.Estimate(2),mld1.Rsquared.Adjusted);
       text(3,12,formula)
       set(gca, 'LooseInset', get(gca,'TightInset'))

       %subplot(2,2,3)
       axes(ha(3));
       plot(sample_bm_data(:,1), sample_bm_data(:,4),'.');
       hold on;
       plot(sample_bm_data(:,1), ypred2,'-','LineWidth',3);
       plot(sample_bm_data(indx1,1), yci2(indx1,1),'r--');
       plot(sample_bm_data(indx1,1), yci2(indx1,2),'r--');
       hold off;
       xlim([-1 15]); ylim([-1 15]);
       ylabel(corr_bm_names{2}); xlabel(corr_bm_names{1});
       formula = sprintf('y = %.2f + %.2f x\nAdj. Rsq = %.2f', ...
           mld2.Coefficients.Estimate(1),mld2.Coefficients.Estimate(2),mld2.Rsquared.Adjusted);
       text(3,12,formula)
       
       %subplot(2,2,2)
       axes(ha(2));
       plot(sample_bm_data(:,2), sample_bm_data(:,3),'.');
       hold on;
       plot(sample_bm_data(:,2), ypred3,'-','LineWidth',3);
       plot(sample_bm_data(indx2,2), yci3(indx2,1),'r--');
       plot(sample_bm_data(indx2,2), yci3(indx2,2),'r--');
       hold off;
       xlim([-1 15]); ylim([-1 15]);
       formula = sprintf('y = %.2f + %.2f x\nAdj. Rsq = %.2f', ...
           mld3.Coefficients.Estimate(1),mld3.Coefficients.Estimate(2),mld3.Rsquared.Adjusted);
       text(3,12,formula)
       
       %subplot(2,2,4)
       axes(ha(4));
       plot(sample_bm_data(:,2), sample_bm_data(:,4),'.');
       hold on;
       plot(sample_bm_data(:,2), ypred4,'-','LineWidth',3);
       plot(sample_bm_data(indx2,2), yci4(indx2,1),'r--');
       plot(sample_bm_data(indx2,2), yci4(indx2,2),'r--');
       hold off;
       xlim([-1 15]); ylim([-1 15]);
       xlabel(corr_bm_names{2});
       formula = sprintf('y = %.2f + %.2f x\nAdj. Rsq = %.2f', ...
           mld4.Coefficients.Estimate(1),mld4.Coefficients.Estimate(2),mld4.Rsquared.Adjusted);
       text(3,12,formula)
       %set(gca,'FontSize',16);
       
       suptitle(['Stage ', num2str(stage), ', D = ', num2str(dist_vec(j+1))]);
       print(fullfile(output_dir,['Stage ', num2str(stage), ', D = ', num2str(dist_vec(j+1)*2)]) ,'-dpng');
       %suptitle(['Grade ', num2str(stage), ', D = ', num2str(dist_vec(j+1))]);
       %print(fullfile(output_dir,['Grade ', num2str(stage), ', D = ', num2str(dist_vec(j+1)*2)]) ,'-dpng'); 
       %suptitle(['Recur = ', num2str(stage-1), ', D = ', num2str(dist_vec(j+1))]);
       %print(fullfile(output_dir,['Recur = ', num2str(stage-1), ', D = ', num2str(dist_vec(j+1)*2)]) ,'-dpng'); 
       close all;
   end
end
%}

for stage = 1:3
   indx_stage = clinical_data.stages == stage;
   for j  = 1:length(dist_vec) 
       bm_data = pair_bm_vec(indx_stage,j);
       bm_data = cat(1, bm_data{:});
       %bm_data = datasample(bm_data,300);
       sample_indx = randperm(size(bm_data,1),1000);
       %sample_bm_data = bm_data(sample_indx,:);
       sample_bm_data = bm_data;
       
       figure; 
       [S,AX,BigAx,H,HAx] = plotmatrix(sample_bm_data(:,1:2),sample_bm_data(:,3:4));
       ylabel(AX(1,1),corr_bm_names{1});
       ylabel(AX(2,1),corr_bm_names{2});
       xlabel(AX(2,1),corr_bm_names{1});
       xlabel(AX(2,2),corr_bm_names{2});
       xlim(AX(1,1),[-1,15]);xlim(AX(2,2),[-1,15]);xlim(AX(1,2),[-1,15]);xlim(AX(2,1),[-1,15]);
       ylim(AX(1,1),[-1,15]);ylim(AX(2,2),[-1,15]);ylim(AX(1,2),[-1,15]);ylim(AX(2,1),[-1,15]);
       title(['Stage ', num2str(stage), ', D = ', num2str(dist_vec(j))]) 
       print(fullfile(output_dir,['Stage ', num2str(stage), ', D = ', num2str(dist_vec(j)*2)]) ,'-dpng');    
       close all;
   end
end

for grade = 1:3
   indx_grade = clinical_data.grades == grade;
   for j  = 1:length(dist_vec)
       bm_data = pair_bm_vec(indx_stage,j);
       bm_data = cat(1, bm_data{:});
       sample_indx = randperm(size(bm_data,1),500);
       %sample_bm_data = bm_data(sample_indx,:);
       sample_bm_data = bm_data;
       
       figure; 
       [S,AX,BigAx,H,HAx] = plotmatrix(sample_bm_data(:,1:2),sample_bm_data(:,3:4));
       ylabel(AX(1,1),corr_bm_names{1});
       ylabel(AX(2,1),corr_bm_names{2});
       xlabel(AX(2,1),corr_bm_names{1});
       xlabel(AX(2,2),corr_bm_names{2});
       xlim(AX(1,1),[-1,15]);xlim(AX(2,2),[-1,15]);xlim(AX(1,2),[-1,15]);xlim(AX(2,1),[-1,15]);
       ylim(AX(1,1),[-1,15]);ylim(AX(2,2),[-1,15]);ylim(AX(1,2),[-1,15]);ylim(AX(2,1),[-1,15]);

       title(['Grade ', num2str(grade), ', D = ', num2str(dist_vec(j))]) 
       print(fullfile(output_dir,['Grade ', num2str(grade), ', D = ', num2str(dist_vec(j)*2)]) ,'-dpng');    
       close all;
   end
end

for rec = 0:1
   indx_grade = clinical_data.recurrent_5yr == rec;
   for j  = 1:length(dist_vec) 
       bm_data = pair_bm_vec(indx_stage,j);
       bm_data = cat(1, bm_data{:});
       sample_indx = randperm(size(bm_data,1),500);
       %sample_bm_data = bm_data(sample_indx,:);
       sample_bm_data = bm_data;
       
       figure; 
       [S,AX,BigAx,H,HAx] = plotmatrix(sample_bm_data(:,1:2),sample_bm_data(:,3:4));
       ylabel(AX(1,1),corr_bm_names{1});
       ylabel(AX(2,1),corr_bm_names{2});
       xlabel(AX(2,1),corr_bm_names{1});
       xlabel(AX(2,2),corr_bm_names{2});
       xlim(AX(1,1),[-1,15]);xlim(AX(2,2),[-1,15]);xlim(AX(1,2),[-1,15]);xlim(AX(2,1),[-1,15]);
       ylim(AX(1,1),[-1,15]);ylim(AX(2,2),[-1,15]);ylim(AX(1,2),[-1,15]);ylim(AX(2,1),[-1,15]);

       title(['Recur ', num2str(rec), ', D = ', num2str(dist_vec(j))]) 
       print(fullfile(output_dir,['Recur =', num2str(rec), ', D = ', num2str(dist_vec(j)*2)]) ,'-dpng');    
       close all;
   end
end


%{
num_pairs = cellfun(@(x) size(x,1), pair_bm_vec,'UniformOutput',false);
num_pairs = cell2mat(num_pairs);
figure;
for stage = 1:3
   indx_stage = clinical_data.stages == stage;
   num_pairs_stage = num_pairs(indx_stage,:);
   num_pair_mean = mean(num_pairs_stage, 1);
   num_pair_se = std(num_pairs_stage,1)./sqrt(length(indx_stage));
   
   errorbar(dist_vec(2:end),num_pair_mean,num_pair_se,'-s',...
       'LineWidth',3)
   xlim([0, max(dist_vec)+1]);
   hold on;
end
legend('Stage 1','Stage 2', 'Stage 3');
xlabel('Distances'); ylabel('Number of pairs');
set(gca,'FontSize',16);

figure;
for grade = 1:3
   indx_stage = clinical_data.grades == grade;
   num_pairs_stage = num_pairs(indx_stage,:);
   num_pair_mean = mean(num_pairs_stage, 1);
   num_pair_se = std(num_pairs_stage,1)./sqrt(length(indx_stage));
   
   errorbar(dist_vec(2:end),num_pair_mean,num_pair_se,'-s','MarkerSize',5,...
       'LineWidth',3)
   xlim([0, max(dist_vec)+1]);
   hold on;
end
legend('Grade 1','Grade 2', 'Grade 3');
xlabel('Distances'); ylabel('Number of pairs');
set(gca,'FontSize',16);

figure;
for rec = 0:1
   indx_stage = clinical_data.recurrent_5yr == rec;
   num_pairs_stage = num_pairs(indx_stage,:);
   num_pair_mean = mean(num_pairs_stage, 1);
   num_pair_se = std(num_pairs_stage,1)./sqrt(length(indx_stage));
   
   errorbar(dist_vec(2:end),num_pair_mean,num_pair_se,'-s','MarkerSize',5,...
       'LineWidth',3)
   xlim([0, max(dist_vec)+1]);
   hold on;
end
legend('Not rec','Rec');
xlabel('Distances'); ylabel('Number of pairs');
set(gca,'FontSize',16);
%}

%{
self_bm_vec = cell(length(filelist),1);

for i = 1:length(filelist)
    imname = filelist{i};
    
    if ~exist(fullfile(coordinate_dir, [imname '.mat']),'file') || ~exist(fullfile(seg_dir, [imname '.jpg']),'file')
        fprintf('Do not have this file: %s i = %d\n',imname,i);
        continue;
    end
    
    im = imread(fullfile(seg_dir,[imname '.jpg']));
  
    tmp = load(fullfile(coordinate_dir, [imname '.mat']));
    x = double(max(1,tmp.x));
    y = double(max(1,tmp.y));
    
    bm_data = tmp.bm_data(:,indx);
    self_bm_vec{i} = cat(2,bm_data,bm_data);
end
%}