data_dir = 'D:\Documents\multiplex';
coordinate_dir = 'C:\Users\luong_nguyen\Box Sync\Colon Cancer Projects\data\csv_data_filtered';
image_dir = 'D:\Documents\colon_cancer_data\H&E_virtual';
seg_dir = 'D:\Documents\multiplex\seg_output_50';

clinical_data = readtable(fullfile('D:\Documents\multiplex','clinical_data_all_spots.csv'),...
    'Delimiter',',');
filelist = clinical_data.spot_name;

% find the index of SMA
output_dir = fullfile(data_dir,'corr_tensor_168_bm_addnoise');

if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

dist_vec = 0:15:120;
%dist_vec = 0:7.5:60;

output_raw_dir = 'C:\Users\luong_nguyen\Box Sync\Colon Cancer Projects\data\raw_data_filtered';
tmp1 = load(fullfile(output_raw_dir,[filelist{1} '.mat']));   
tmp2 = load(fullfile(coordinate_dir,[filelist{1} '.mat']));

bm_names = tmp1.bm_names;
bm_names{12} = 'ColIV';
%long_bm_names = cellfun(@(x) ['Median.Cell.' x], bm_names,'UniformOutput',false);
long_bm_names = cellfun(@(x) {['Median.Nuc.' x],['Median.Memb.' x],['Median.Cyt.' x]}, ...
    bm_names,'UniformOutput',false);
long_bm_names = cat(2, long_bm_names{:})';

indx_bm = cellfun(@(x) find(ismember(cellstr(tmp2.bm_names)',x)),...
    long_bm_names,'UniformOutput',false);
indx_bm = cat(1,indx_bm{:});

tic;
for i = 1:length(filelist)
    imname = filelist{i};
    
    if ~exist(fullfile(coordinate_dir, [imname '.mat']),'file') || ~exist(fullfile(seg_dir, [imname '.jpg']),'file')
        fprintf('Do not have this file: %s\n',imname);
        continue;
    end
     
    tmp = load(fullfile(coordinate_dir, [imname '.mat']));
    x = double(max(1,tmp.x));
    y = double(max(1,tmp.y));
    
    bm_data = tmp.bm_data(:,indx_bm);
    pairwise_dist = squareform(pdist([x;y]'));
    
    corr_tensor = zeros(length(long_bm_names),length(long_bm_names),length(dist_vec));
    
    for j = 1:length(dist_vec)
        if j == 1
            corr_tensor(:,:,j) = corr(bm_data + randn(size(bm_data))*1e-9);
        else
            [r,c] = find(pairwise_dist > dist_vec(j-1) & ...
                pairwise_dist <= dist_vec(j));
            ind0 = r<c;
            r = r(ind0); c = c(ind0);
            corr_tensor(:,:,j) = corr(bm_data(r,:)+randn(size(bm_data(r,:)))*1e-9,...
                bm_data(c,:)+randn(size(bm_data(r,:)))*1e-9);
        end
        if sum(sum(isnan(corr_tensor(:,:,j)))) > 0
            a = corr_tensor(:,:,j);
            a(isnan(a)) = 0;
            corr_tensor(:,:,j) = a;
            fprintf('There is NAN in image %s\n',imname);
        end
    end
    save(fullfile(output_dir,[imname '.mat']),'corr_tensor');
end
toc
for j = 1:length(dist_vec)
    figure; imagesc(corr_tensor(:,:,j)); axis square; colorbar; caxis([-1 1]);
end
