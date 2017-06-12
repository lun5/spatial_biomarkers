%% classification of fibroblast vs. immune cells
% training data
coordinate_dir = 'D:\Documents\multiplex\NN_data\test_biomarker_intensity\nuc_morph';
image_dir = 'D:\Documents\multiplex\H&E_renamed';
train_dir = fullfile(coordinate_dir,'train');
if ~exist(train_dir,'dir')
    mkdir(train_dir);
end

filelist = dir(fullfile(train_dir,'*.mat'));
filelist = {filelist.name}';

morph_feature_names = {'x','y','Area', 'Perimeter','boundingbox_width','boundingbox_height',...
    'MajorAxisLength','MinorAxisLength','Extent','EulerNumber', 'Eccentricity', 'Orientation',...
     'Solidity','EquivDiameter'};
channel_names = {'r','g','b','hsv','lab','luv','he','br'};
first_order_stats_names = {'mean_intensity', 'median_intensity', 'var_intensity',...
    'kurtosis_intensity','skewness_intensity'};
all_feat_names = {};
for i =1:length(channel_names)
    all_feat_names{i} = cellfun(@(x) [x '_' channel_names{i}], first_order_stats_names,'UniformOutput',false);
end
all_feat_names = cat(2, morph_feature_names,all_feat_names{:});

feature_names = {'x','y','Area','MajorAxisLength','MinorAxisLength','Eccentricity'};
indx_features = find(ismember(all_feat_names,feature_names));

load(fullfile(train_dir,filelist{1}));
bm_names = cellstr(tmp.bm_names);
lineage_bm_names = {'Median.Nuc.CD20','Median.Nuc.CD8','Median.Nuc.CD79',...
    'Median.Nuc.CD3','Median.Nuc.EPCAM', 'Median.Nuc.SMA'};
indx_bm = find(ismember(bm_names,lineage_bm_names));

features = cell(length(filelist),1);
labels = cell(length(filelist),1);
cell_segs = cell(length(filelist),1);
for i = 1:length(filelist)
    fname = filelist{i}(1:end-4);
    im = imread(fullfile(image_dir,[fname,'.jpg'])); % H&E image
    
    load(fullfile(train_dir,[fname '.mat']));
    chosen_cells = find((tmp.fibroblast | tmp.immune ) | tmp.til);
    
    curr_cell_seg = cell(length(chosen_cells),1);
    for j = 1:length(chosen_cells)
        bb = regionprops(tmp.filtered_mask_im == chosen_cells(j),'BoundingBox');
        curr_cell_seg{j} = imcrop(im.*repmat(uint8(...
            tmp.filtered_mask_im==chosen_cells(j)),1,1,3),bb.BoundingBox);
    end
    cell_segs{i} = curr_cell_seg;
    labels{i} = cat(1,tmp.fibroblast(chosen_cells), tmp.immune(chosen_cells),...
        tmp.til(chosen_cells))';
    features{i} = cat(2,tmp.features(chosen_cells,3:end),...
        tmp.bm_data(chosen_cells,indx_bm));
end

cell_segs = cat(1,cell_segs{:});
labels = cat(1,labels{:});
features = cat(1,features{:});
curr_labels = labels(:,1) + 2*labels(:,2) + 2*labels(:,3);

% calculate the distance matrix
rand_indx = randperm(size(features,1),1000);
distance_matrix = squareform(pdist(features(rand_indx,13:52),'correlation'));
curr_labels = curr_labels(rand_indx);
curr_segs = cell_segs(rand_indx);
%labels = labels;
% for some reason entry 718 has value 495. why?
%labels = labels + 1;
cmap = [0 0 1; 1 0 0;0 1 1];

P = d2p(distance_matrix);
figure;
Y = tsne_p(P,curr_labels);

figure; gscatter(Y(:,1),Y(:,2),curr_labels);
legend('fibroblast','immune');set(gca,'FontSize',16);

figure; gscatter(Y(:,1),Y(:,2),curr_labels);
hold on;
for i = 1:length(rand_indx)
   imagesc([Y(i,1),Y(i,1)+3],[Y(i,2),Y(i,2)+3],curr_segs{i});     
end
legend('fibroblast','immune');set(gca,'FontSize',16);


