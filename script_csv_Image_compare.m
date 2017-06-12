%% identify the correct segmentation mask
%{
imdir = 'X:\ColonCancerStudy\AGA_260_3\Results\20140430_015736_305\SingleCellAnalysis\001';
imlist = dir(fullfile(imdir, '*.tif'));
imlist = {imlist.name}';
imname = 'NUCLEI_SEG.TIF';

seg_im = imread(fullfile(imdir, imname));
figure; imshow(imadjust(seg_im));

CC = bwconncomp(seg_im);
%label_im = bwlabel(bin_im);

figure; imshow(label2rgb(label_im));
bdry_im = seg2bdry(label_im);


for i = 1:length(imlist)
    im = imread(fullfile(imdir, imlist{i}));
    if size(im,3) == 1
        figure; imshow(imadjust(im));
    else
        figure; imshow(im);
    end
end
%}
%% registered images
%{
reg_im_dir = 'X:\ColonCancerStudy\AGA_260_3\RegisteredImages\S002_BetaCatenin_S6';
cy3_imname = 'BetaCatenin_S6_AGA_260_3_S002_P001_cy3.tif';
cy3_im = imread(fullfile(reg_im_dir,cy3_imname));
whos cy3_im
figure; imshow(imadjust(cy3_im));

cy5_imname = 'BetaCatenin_S6_AGA_260_3_S002_P001_cy5.tif';
cy5_im = imread(fullfile(reg_im_dir,cy5_imname));
whos cy5_im
figure; imshow(imadjust(cy5_im));

corr_imname = 'BetaCatenin_S6_AGA_260_3_S002_P001_correl.tif';
corr_im = imread(fullfile(reg_im_dir, corr_imname));
whos corr_im
figure; imshow(imadjust(corr_im));

dapi_imname = 'BetaCatenin_S6_AGA_260_3_S002_P001_dapi.tif';
dapi_im = imread(fullfile(reg_im_dir,dapi_imname));
figure; imshow(imadjust(dapi_im));

% now take the mean or median of 
cy3_im_log = log2(double(cy3_im));
cy3_im_log(isinf(cy3_im_log)) = 0;
stats = regionprops(CC, cy3_im_log,'PixelValues','centroid');

for k = 1:length(stats)
    stats(k).MedianIntensity = median(stats(k).PixelValues);
end

centroids = cat(1,stats.Centroid);
medianIntensity = cat(1, stats.MedianIntensity);

figure; imshow(imadjust(dapi_im));
hold on; plot(centroids(:,1), centroids(:,2), 'bx');
hold off;
%}
%% auto removed images
%{
auto_dir = 'X:\ColonCancerStudy\AGA_260_3\AFRemoved';
imname = 'BetaCatenin_AFRemoved_001.tif';
afr_im = imread(fullfile(auto_dir,imname));

figure; imshow(imadjust(afr_im));
hold on; plot(centroids(:,1), centroids(:,2), 'bo');
hold off;

figure; imshow(imadjust(cy3_im));
hold on; plot(centroids(:,1), centroids(:,2), 'bo');
hold off;

prod_im = uint16((double(afr_im.*cy3_im)).^(0.5));
figure; imshow(imadjust(prod_im));
figure; imshowpair(imadjust(afr_im),imadjust(cy3_im),'diff')
figure; imshowpair(imadjust(afr_im),imadjust(cy3_im))

%
imname = 'S6_AFRemoved_001.tif';
afr_im = imread(fullfile(auto_dir,imname));

figure; imshow(imadjust(afr_im));
hold on; plot(centroids(:,1), centroids(:,2), 'bo');
hold off;

figure; imshow(imadjust(cy5_im));
hold on; plot(centroids(:,1), centroids(:,2), 'bo');
hold off;

prod_im = uint16((double(afr_im.*cy5_im)).^(0.5));
figure; imshow(imadjust(prod_im));
%}
%% print out details for these biomarkers
% bm_names = {'pERK', 'CD31', 'BetaCatenin', 'S6', 'pS6235', 'beta_actin', 'pck26',...
%     'Glut1', 'NaKATPase', 'SMA', 'Albumin', 'EGFR', 'p4EBP1', 'pNDRG1', 'MLH1',...
%     'LaminA_C', 'pGSK3beta', 'EZH2', 'Claudin1', 'NDRG1', 'CD68', 'CD8', 'PTEN',...
%     'pMAPKAPK2', 'Akt', 'CA9', 'CleavedCaspase3', 'ERK', 'EPCAM', 'CD3', 'MSH2',...
%     '4EBP1', 'COX2', 'p53', 'ColIV'};
% 
data_dir = 'X:\ColonCancerStudy\';
table_dir = 'C:\Users\luong_nguyen\Box Sync\Colon Cancer Projects\data';
output_dir = 'D:\Documents\multiplex\NN_data\test_biomarker_intensity\from_raw_images';

if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

T = readtable(fullfile(table_dir,'clinical_data_all_spots.csv'),'Delimiter',',');
filelist = T.spot_name;

output_raw_dir = fullfile(table_dir,'raw_data_filtered');
tmp1 = load(fullfile(output_raw_dir,[filelist{1} '.mat']));
bm_names = tmp1.bm_names(1:end-1); % There is NO DAPI

single_cell_dir = {'20140430_015736_305','20140501_083300_738','20140501_083645_537'};
for i = 1:length(filelist)
   spot_name = filelist{i};
   split_name = strsplit(spot_name,'_');
   slide_name = strjoin(split_name(1:3),'_');
   pos = str2num(split_name{4});
%    if exist(fullfile(output_dir,[spot_name '.mat']),'file')
%        fprintf('Already calculated for spot %s\n',spot_name);
%        continue;
%    end
   tic;
   seg_dir = fullfile(data_dir,slide_name,'Results');
   switch slide_name
       case 'AGA_260_3'
           seg_dir = fullfile(seg_dir,single_cell_dir{1},'SingleCellAnalysis');
       case 'AGTA_264_3'
           seg_dir = fullfile(seg_dir,single_cell_dir{2},'SingleCellAnalysis');
       otherwise
           seg_dir = fullfile(seg_dir,single_cell_dir{3},'SingleCellAnalysis');
   end
   imname = 'NUCLEI_SEG.TIF';
   seg_im = imread(fullfile(seg_dir,sprintf('%03d',pos), imname));    
   auto_dir = fullfile(data_dir,slide_name,'AFRemoved');
   stats = regionprops(seg_im, 'Centroid','Area');
   xy = cat(1,stats.Centroid); x = xy(:,1)'; y = xy(:,2)';
   area = cat(1,stats.Area);
   bm_data = zeros(max(seg_im(:)),length(bm_names));
   for j = 1:(length(bm_names)-1)
       afr_im_name = sprintf('%s_AFRemoved_%03d.tif',bm_names{j},pos);
       if ~exist(fullfile(auto_dir, afr_im_name),'file')
           if strcmp(bm_names{j},'COLIV')
              afr_im_name = strrep(afr_im_name,'COLIV','collagenIV');
           else
              continue; 
           end
       end
       
       afr_im = imread(fullfile(auto_dir, afr_im_name));
       afr_im_log = log2(double(afr_im));
       afr_im_log(isinf(afr_im_log)) = 0;
       stats = regionprops(seg_im, afr_im_log,'PixelValues');
       for k = 1:length(stats)
           stats(k).MedianIntensity = median(stats(k).PixelValues);
       end
       bm_data(:,j) = cat(1, stats.MedianIntensity);
   end
   fprintf('Done with %s in %.2f seconds\n',spot_name, toc);
   save(fullfile(output_dir,[spot_name '.mat']),'x','y','area',...
       'bm_data','bm_names','seg_im');   
end


%% matching between the csv files and the matlab output
csv_output_dir = fullfile(data_dir,'NN_data\test_biomarker_intensity\from_csv');
matlab_output_dir = fullfile(data_dir,'NN_data\test_biomarker_intensity\from_matlab');

for i = 1:length(test_spot_names)
    spot_name = test_spot_names{i};
    csv_file = load(fullfile(csv_output_dir,[spot_name,'.mat']));
    matlab_file = load(fullfile(matlab_output_dir,[spot_name,'.mat']));
    
    clear xy areas biomarkers
    xy{1} = csv_file.xy;
    xy{2} = matlab_file.xy;    
    area{1} = csv_file.areas;
    area{2} = matlab_file.areas;
    bm_data{1} = csv_file.biomarkers;
    bm_data{2} = matlab_file.biomarkers;
    
    cell_radius = sqrt(area{1}/pi);
    threshold = prctile(cell_radius*2,95);
    match_id = zeros(size(area{1}));
    distances = pdist2(xy{1},xy{2});
    [min_val, min_indx] = min(distances,[],2);
    min_indx(min_val > threshold) = 0;
    indx{1} = find(min_indx > 0);
    indx{2} = min_indx(min_indx > 0);
    
    figure; histogram(cell_radius*2,30,'Normalization','probability','FaceColor',[.8,.8,.8]);
    xlabel('Cell diameter');ylabel('probability');
    ax = gca;
    set(gca,'FontSize',16); hold on; plot([threshold threshold], ax.YLim,'LineWidth',3); hold off
    
    figure; histogram(min_val,30,'Normalization','probability','FaceColor',[.8,.8,.8]);
    xlabel('Minimum distances to matched cells');ylabel('probability');
    ax = gca;
    set(gca,'FontSize',16);hold on; plot([threshold threshold], ax.YLim,'LineWidth',3); hold off
    
    bm_data{1} = bm_data{1}(indx{1},:);
    bm_data{2} = bm_data{2}(indx{2},:);
    
    figure;
    for j = 1:9
        x = bm_data{1}(:,j);
        y = bm_data{2}(:,j);
        mdl = fitlm(x,y);
        subplot(3,3,j); 
        plot(x, y,'bx'); lsline
        title(sprintf('%s (R^2 = %.2f)',bm_names{j},mdl.Rsquared.Adjusted));
        xlim([0 14]); ylim([0,14]);
    end
    
    figure;
    for j = 10:18
        x = bm_data{1}(:,j);
        y = bm_data{2}(:,j);
        mdl = fitlm(x,y);
        subplot(3,3,j-9); 
        plot(x, y,'bx'); lsline
        title(sprintf('%s (R^2 = %.2f)',bm_names{j},mdl.Rsquared.Adjusted));
        xlim([0 14]); ylim([0,14]);
    end
    
    figure;
    for j = 19:27
        x = bm_data{1}(:,j);
        y = bm_data{2}(:,j);
        mdl = fitlm(x,y);
        subplot(3,3,j-18); 
        plot(x, y,'bx'); lsline
        title(sprintf('%s (R^2 = %.2f)',bm_names{j},mdl.Rsquared.Adjusted));
        xlim([0 14]); ylim([0,14]);
    end
    
    figure;
    for j = 28:35
        x = bm_data{1}(:,j);
        y = bm_data{2}(:,j);
        mdl = fitlm(x,y);
        subplot(3,3,j-27); 
        plot(x, y,'bx'); lsline
        title(sprintf('%s (R^2 = %.2f)',bm_names{j},mdl.Rsquared.Adjusted));
        xlim([0 14]); ylim([0,14]);
    end
    
    split_name = strsplit(spot_name,'_');
    slide_name = strjoin(split_name(1:3),'_');
    pos = str2num(split_name{4});
    imname = 'NUCLEI_SEG.TIF';
    
    seg_dir = fullfile(data_dir,slide_name,'Results');
    switch slide_name
        case 'AGA_260_3'
            seg_dir = fullfile(seg_dir,single_cell_dir{1},'SingleCellAnalysis');
        case 'AGTA_264_3'
            seg_dir = fullfile(seg_dir,single_cell_dir{2},'SingleCellAnalysis');
        otherwise
            seg_dir = fullfile(seg_dir,single_cell_dir{3},'SingleCellAnalysis');
    end
    seg_im = imread(fullfile(seg_dir,sprintf('%03d',pos), imname));
    bin_im = seg_im > 0;
    
    figure; imshow(bin_im); hold on;
    plot(xy{1}(indx{1},1), xy{1}(indx{1},2),'rx','MarkerSize',10);
    %plot(xy{1}(indx{1},2), xy{1}(indx{1},1),'rx');
    plot(xy{2}(indx{2},1), xy{2}(indx{2},2),'b.','MarkerSize',10);
    hold off;
    
    corrplot(bm_data{1}(:,1:5),'varNames',bm_names(1:5))
    %set(gca,'FontSize',16)
    corrplot(bm_data{2}(:,1:5),'varNames',bm_names(1:5))
end


