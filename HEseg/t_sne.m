% t-sne presentation
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

output_dir = 'D:\Documents\multiplex\spots_corrcoef';

if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

% find the index of SMA
tmp = load(fullfile(coordinate_dir,filelist{1}));
bm_names = cellstr(tmp.bm_names);
bm_names = cellfun(@(x) strsplit(x,'.'),bm_names,'UniformOutput',false);
bm_names = unique(cat(1,bm_names{:}));
bm_names = cellfun(@(x) ['Median.Nuc.' x], bm_names,'UniformOutput',false);
indx = find(ismember(tmp.bm_names,bm_names));

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
    corr_matrics{count} = corrcoef(bm_data);
    stages{count} = clinical_data.stages(i);
    grades{count} = clinical_data.grades(i);
    rec{count} = clinical_data.recurrent_5yr(i);
    chemo{count} = clinical_data.chemo(i);
end

save(fullfile(output_dir, 'correlation_spots.mat'),'spot_names','stages',...
    'grades','rec','chemo','bm_names');    
%}

% calculate the Frobenious norm
distance_matrix = zeros(length(spot_names),length(spot_names));

for i = 1:length(spot_names) - 1
    for j = (i+1):length(spot_names)
       diff_mat =  abs(corr_matrics{i} - corr_matrics{j});
       distance_matrix(i,j) = norm(diff_mat,'fro');
    end 
end

labels = cell2mat(recur5);
labels = labels';
% for some reason entry 718 has value 495. why?
labels(labels>0) = 1;
labels = labels + 1;

cmap = [0 0 1; 1 0 0];


P = d2p(N);
Y = tsne_p(P,labels);

% 2d projection
figure; 
for i = 1:length(labels), 
	plot(Y(i,1),Y(i,2),'o','color',cmap(labels(i),:)); 
	hold on; 
end

% 3d projection
Y = tsne_p(P,labels,3);
figure; 
for i = 1:length(labels), 
	plot3(Y(i,1),Y(i,2),Y(i,3),'o','color',cmap(labels(i),:)); 
	hold on; 
end
grid on;

