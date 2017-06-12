% t-sne presentation
data_dir = 'D:\Documents\multiplex';
coordinate_dir = 'D:\Documents\multiplex\coordinates_all_bm';
%coordinate_dir = 'C:\Users\luong_nguyen\Box Sync\Colon Cancer Projects\data\raw_data_filtered';
image_dir = 'D:\Documents\colon_cancer_data\H&E_virtual';
seg_dir = 'D:\Documents\multiplex\seg_output_50';

%filelist = dir(fullfile(coordinate_dir,'*.mat'));
%filelist = {filelist.name}';

% read in the clincial data file
clinical_data = readtable(fullfile('D:\Documents\multiplex','clinical_data_all_spots.csv'),...
    'Delimiter',',');
filelist = clinical_data.spot_name;

output_dir = 'D:\Documents\multiplex\spots_corrcoef_raw_filtered';

if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

tmp = load(fullfile(coordinate_dir,filelist{1}));
bm_names = cellstr(tmp.bm_names);
bm_names = cellfun(@(x) strsplit(x,'.'),bm_names,'UniformOutput',false);
bm_names = unique(cat(1,bm_names{:}));
bm_names = cellfun(@(x) ['Median.Cell.' x], bm_names,'UniformOutput',false);
%indx = find(ismember(tmp.bm_names,bm_names));
indx = ~(ismember(tmp.bm_names,bm_names));

%{
spot_names = {};
corr_matrics = {};
stages = {}; grades = {}; rec = {}; chemo = {};
count = 0;
for i = 1:length(filelist)
    imname = filelist{i};
    
    if ~exist(fullfile(coordinate_dir, [imname '.mat']),'file') || ~exist(fullfile(seg_dir, [imname '.jpg']),'file')
        fprintf('Do not have this file: %s\n',imname);
        continue;
    end
    count = count + 1;    
    spot_names{count} = imname;
    tmp = load(fullfile(coordinate_dir, [imname '.mat']));
    bm_data = tmp.bm_data(:,indx);
    %corr_matrics{count} = corrcoef(bm_data);
    corr_matrices{count} = corr(bm_data);
    stages{count} = clinical_data.stages(i);
    grades{count} = clinical_data.grades(i);
    rec{count} = clinical_data.recurrent_5yr(i);
    chemo{count} = clinical_data.chemo(i);
end

save(fullfile(output_dir, 'correlation_spots_2.mat'),'spot_names','stages',...
    'grades','rec','chemo','bm_names');    
%}


%load(fullfile(output_dir, 'correlation_spots_2.mat'));
% calculate the Frobenious norm

%corr_matrix_dir = 'C:\Users\luong_nguyen\Box Sync\Colon Cancer Projects\data\corr_tensor_168_bm_addnoise';
corr_matrix_dir = 'D:\Documents\multiplex\corr_tensor_168_bm_addnoise';
spot_names = dir(fullfile(corr_matrix_dir,'*.mat'));
spot_names = {spot_names.name}';

stages = zeros(length(spot_names),1);
grades = zeros(length(spot_names),1);
recurs = zeros(length(spot_names),1);

for i =1:length(spot_names)
    stages(i) = clinical_data(strcmpi(spot_names{i}(1:end-4),clinical_data.spot_name),:).stages;
    grades(i) = clinical_data(strcmpi(spot_names{i}(1:end-4),clinical_data.spot_name),:).grades;
    recurs(i) = clinical_data(strcmpi(spot_names{i}(1:end-4),clinical_data.spot_name),:).recurrent_5yr;
end

distance_matrix = zeros(length(spot_names),length(spot_names));
tic;
for i = 1:(length(spot_names) - 1)
    imname = spot_names{i}(1:end-4);
    tmp = load(fullfile(coordinate_dir, [imname '.mat']));
    mean_i = mean(tmp.bm_data(:,indx),1);
    for j = (i+1):length(spot_names)
        imname = spot_names{j}(1:end-4);
        tmp = load(fullfile(coordinate_dir, [imname '.mat']));
        mean_j = mean(tmp.bm_data(:,indx),1);
        distance_matrix(i,j) = norm(mean_i-mean_j);
        distance_matrix(j,i) = distance_matrix(i,j);
    end
end
toc

%{
for i = 1:length(spot_names) - 1
    for j = (i+1):length(spot_names)
       diff_mat =  corr_matrices{i} - corr_matrices{j};
       distance_matrix(i,j) = norm(diff_mat,'fro');
       distance_matrix(j,i) = distance_matrix(i,j);
    end 
end
%}

%{
corr_matrices = cell(length(spot_names),1);
for i = 1:length(spot_names) 
    tmp = load(fullfile(corr_matrix_dir,spot_names{i}));
    corr_matrices{i} = tmp.corr_tensor;
end

layer_indx = 1:9;
distance_matrix = zeros(length(spot_names),length(spot_names),length(layer_indx));
tic;
for i = 1:(length(spot_names) - 1)
    %t = tic;
    for j = (i+1):length(spot_names)
        for k = 1:length(layer_indx)
           distance_matrix(i,j,k) = distance_matrix(i,j,k) + ...
               norm(corr_matrices{i}(:,:,layer_indx(k))-corr_matrices{j}(:,:,layer_indx(k)),'fro');
           distance_matrix(j,i,k) = distance_matrix(i,j,k);
        end
    end
    %toc(t)
end
toc
%}
dist_vec = 0:15:120;
%tmp = load('C:\Users\luong_nguyen\Box Sync\Colon Cancer Projects\data\pmiDistance_frobenious_30D_15atom.mat');
tmp = load('C:\Users\luong_nguyen\Box Sync\Colon Cancer Projects\data\corr_tensor_168_bm_addnoise\fro_distances.mat');
curr_distance_matrix = sum(tmp.distance_matrix(:,:,[1,5,9]),3);
labels = tmp.stages;

%curr_distance_matrix = distance_matrix;
%labels = stages;

%curr_distance_matrix = tmp.N;
%
%indx_stage13 = (labels == 1) | (labels == 3); 
%curr_distance_matrix = curr_distance_matrix(indx_stage13,indx_stage13);
%labels = labels(indx_stage13);

%labels = labels;
% for some reason entry 718 has value 495. why?
%labels = labels + 1;
cmap = [1 0 0; 0 1 0;0 0 1];
P = d2p(curr_distance_matrix);
figure;
Y = tsne_p(P,labels);

%figure; gscatter(Y(:,1),Y(:,2),labels);
%close all;

sz = 35;

figure; scatter(Y(labels==1,1),Y(labels==1,2),sz,'rs','filled'); hold on
scatter(Y(labels==2,1),Y(labels==2,2),sz,'go','filled');
scatter(Y(labels==3,1),Y(labels==3,2),sz,'b>','filled');
hold off;
legend('Stage 1','Stage 2', 'Stage 3');
set(gca,'FontSize',14)

% 2d projection
%{
figure; 
for i = 1:length(labels), 
	plot(Y(i,1),Y(i,2),'o','color',cmap(labels(i),:)); 
	hold on; 
end
%}
% 3d projection
Y = tsne_p(P,labels,3);
figure; 
for i = 1:length(labels), 
	plot3(Y(i,1),Y(i,2),Y(i,3),'o','color',cmap(labels(i),:)); 
	hold on; 
end
grid on;

