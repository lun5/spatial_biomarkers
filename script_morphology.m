%% script to calculate the morphology of stroma cells
% 5/10/2017
% We want to extract the morphology data to fibroblast vs. immune cells
clear all
filtered_raw_dir = 'C:\Users\luong_nguyen\Box Sync\Colon Cancer Projects\data\raw_data_filtered';
raw_img_dir = 'D:\Documents\multiplex\NN_data\test_biomarker_intensity\from_raw_images';
image_dir = 'D:\Documents\multiplex\H&E_renamed';

data_dir = 'C:\Users\luong_nguyen\Box Sync\Colon Cancer Projects\data';
csv_dir = fullfile(data_dir,'coordinates_all_bm');
T = readtable(fullfile(data_dir,'clinical_data_all_spots.csv'),'Delimiter',',');

filelist = T.spot_name;

bm_compare_dir = 'D:\Documents\multiplex\NN_data\test_biomarker_intensity\compare_image_extracted_bm';
if ~exist(bm_compare_dir,'dir')
    mkdir(bm_compare_dir);
end

input_csv_dir = 'C:\Users\luong_nguyen\Box Sync\Colon Cancer Projects\data\csv_data_filtered';

morph_code = 'C:\Users\luong_nguyen\Documents\GitHub\tissue-component-classification\';
addpath(genpath(morph_code));

nuc_moroph_dir = 'D:\Documents\multiplex\NN_data\test_biomarker_intensity\nuc_morph';
if ~exist(nuc_moroph_dir,'dir')
    mkdir(nuc_moroph_dir);
end

for i = 1:length(filelist)
    fname = filelist{i};
    %spot_name = 'AGA_260_3_115';
    if ~exist(fullfile(filtered_raw_dir,[fname,'.mat']),'file') 
        fprintf('No %s in raw data\n',fname);
        continue
    elseif ~exist(fullfile(raw_img_dir,[fname,'.mat']),'file')
        fprintf('No %s in image extracted data\n',fname);
        continue
    elseif ~exist(fullfile(image_dir,[fname,'.jpg']),'file')
        fprintf('No %s H&E image\n',fname);
        continue;
    end

    csv_file = load(fullfile(filtered_raw_dir,[fname,'.mat']));
    matlab_file = load(fullfile(raw_img_dir,[fname,'.mat']));
   
    clear xy areas biomarkers
    xy{1} = [csv_file.x;csv_file.y]';
    xy{2} = [matlab_file.x; matlab_file.y]';
    
    im = imread(fullfile(image_dir,[fname,'.jpg'])); % H&E image
    seg_im = matlab_file.seg_im;
     
    %figure; subplot(1,2,1); imshow(im); subplot(1,2,2); imshow(mask_im);
    
    area{1} = csv_file.area;
    area{2} = matlab_file.area;
    bm_data{1} = csv_file.bm_data;
    bm_data{2} = matlab_file.bm_data;
    
    bm_names = csv_file.bm_names;
    cell_radius = sqrt(area{1}/pi);
    threshold = prctile(cell_radius*2,20);
    distances = pdist2(xy{1},xy{2});
    [min_val, min_indx] = min(distances,[],2);
    min_indx(min_val > threshold) = 0;
    indx{1} = find(min_indx > 0);
    indx{2} = min_indx(min_indx > 0);
    [c,ia,ic] = unique(indx{2});
    indx{2} = c; indx{1} = indx{1}(ia);
    
    clear csv_file matlab_file threshold cell_radius min_indx min_val c ia ic
%     figure; histogram(cell_radius*2,30,'Normalization','probability','FaceColor',[.8,.8,.8]);
%     xlabel('Cell diameter');ylabel('probability');
%     ax = gca;
%     set(gca,'FontSize',16); hold on; plot([threshold threshold], ax.YLim,'LineWidth',3); hold off
%     
%     figure; histogram(min_val,30,'Normalization','probability','FaceColor',[.8,.8,.8]);
%     xlabel('Minimum distances to matched cells');ylabel('probability');
%     ax = gca;
%     set(gca,'FontSize',16);hold on; plot([threshold threshold], ax.YLim,'LineWidth',3); hold off
    
    bm_data{1} = bm_data{1}(indx{1},:);
    bm_data{2} = bm_data{2}(indx{2},:);
    
    % scatter plot for all 56
    figure('position',[100 100 1200 1000]);
    ha = tight_subplot(7,8,[.03 .03],[.08 .01],[.07 .01]);
   
    for j = 1:length(bm_names)-1
        axes(ha(j));
        plot(bm_data{1}(:,j), bm_data{2}(:,j),'.');
        xlim([-1 15]); ylim([-1 15]);
        text(3,10,bm_names{j});grid on;
        set(gca, 'LooseInset', get(gca,'TightInset'))
    end
    print(fullfile(bm_compare_dir,fname ) ,'-dpng');
    close all
    clear bm_data area
    tic; filtered_mask_im = arrayfun(@(x) seg_im == x, indx{2}, 'UniformOutput',false);toc
    tic; filtered_mask_im = cat(3,filtered_mask_im{:}); toc
    filtered_mask_im = uint16(filtered_mask_im);
    %tic;filtered_seg_im = arrayfun(@(x) (seg_im == indx{2}(x)).*x, 1:length(indx{2}), 'UniformOutput',false); toc
    %tic;filtered_seg_im = cat(3,filtered_seg_im{:});toc
    tic;
    for m =1:length(indx{2})
       filtered_mask_im(:,:,m) = filtered_mask_im(:,:,m)*m; 
    end
    toc
    filtered_mask_im = sum(filtered_mask_im,3);
    
    figure; imshow(im.*repmat(uint8(filtered_mask_im > 0),1,1,3)); hold on;
    plot(xy{1}(indx{1},1), xy{1}(indx{1},2),'rx','MarkerSize',10);
    plot(xy{2}(indx{2},1), xy{2}(indx{2},2),'b.','MarkerSize',10);
    hold off;
    print(fullfile(bm_compare_dir,[fname '_overlay'] ) ,'-dpng');
    close all;
    
    % update csv file
    tmp = load(fullfile(input_csv_dir,[fname,'.mat']));
    x = tmp.x(indx{1});
    y = tmp.y(indx{1});
    area = tmp.area(indx{1});
    epithelial = tmp.epithelial(indx{1});
    bm_data = tmp.bm_data(indx{1},:);
    bm_names = tmp.bm_names;
    features = nucleiFeatures1(filtered_mask_im, im);
    
    save(fullfile(nuc_moroph_dir,[fname '.mat']),'x','y','area','bm_data','bm_names','epithelial','filtered_mask_im','features');
    clear tmp x y area epithelial bm_data bm_names features filtered_mask_im
end


morph_feature_names = {'x','y','Area', 'Perimeter','boundingbox_width','boundingbox_height',...
    'MajorAxisLength','MinorAxisLength','Extent','EulerNumber', 'Eccentricity', 'Orientation',...
     'Solidity','EquivDiameter'};

channel_names = {'r','g','b','hsv','lab','luv','he','br'};

first_order_stats_names = {'mean_intensity', 'median_intensity', 'var_intensity',...
    'kurtosis_intensity','skewness_intensity'};
harr_features_names = {'correllation', 'clusterShape', 'clusterProminence', 'energy', ...
    'entropy', 'contrast','homogeneity'};
runlength_features_names = {'SRE', 'LRE', 'GLN', 'RLN',  'RP', 'LGRE', 'HGRE',...
    'SGLGE', 'SRHGE', 'LRLGE',  'LRHGE'};

texture_features = cat(2,first_order_stats_names, harr_features_names, runlength_features_names);

all_feat_names = {};

for i =1:length(channel_names)
    all_feat_names{i} = cellfun(@(x) [x '_' channel_names{i}], first_order_stats_names,'UniformOutput',false);
end

all_feat_names = cat(2, morph_feature_names,all_feat_names{:});
% filtered_seg_im = seg_im;
% for k = 1:length(indx{2})
%     filtered_seg_im(seg_im == indx{2}(k)) = -10;
% end
% 
% filtered_seg_im(filtered_seg_im > 0) = 0;
% filtered_seg_im = filtered_seg_im./(-10);
% filtered_mask_im = im.*repmat(uint8(filtered_seg_im),1,1,3);
