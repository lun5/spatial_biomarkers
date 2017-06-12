%% script to segment the data
data_dir = 'D:\Documents\multiplex\coordinates';
image_dir = 'D:\Documents\multiplex\H&E_virtual';

filelist = dir(fullfile(data_dir,'*.mat'));
filelist = {filelist.name}';

threshold_vec = [50,75];
renamed_dir = 'D:\Documents\multiplex\H&E_renamed';
if ~exist(renamed_dir,'dir')
    mkdir(renamed_dir);
end

%{
for i = 1:length(filelist)
    t1 = tic;
    fname = filelist{i};
    fname_split = strsplit(fname,{'_','.'});
    slide_name = strjoin(fname_split(1:3),'_');    
    pos = str2double(fname_split(4));
    if strcmp(slide_name,'AGTA_269_3')
        imname = fullfile(image_dir,slide_name,'VHE',['VHE' sprintf('%04d',pos) '.jpg']);
    else
        imname = fullfile(image_dir,slide_name,'VHE',['VHE' sprintf('%03d',pos) '.jpg']);
    end
    
    if ~exist(imname,'file')
        fprintf('image %s does not exist\n', fname)
        continue;
    end
    copyfile(imname, fullfile(renamed_dir,[fname(1:end-4) '.jpg']));
end
%}
se = strel('disk',7,4);
output_dir = 'D:\Documents\multiplex\seg_output_50';
for i = 1:length(filelist)
    fname = filelist{i};
    if ~exist(fullfile(output_dir,fname),'file')
        continue;
    end
    tmp = load(fullfile(output_dir,fname));
    seg_result = tmp.data;
    save(fullfile(output_dir,fname),'seg_result');
%     tmp = load(fullfile(output_dir,fname));
%     seg = tmp.data;
%     bdry = seg2bdry(double(seg),'imageSize');
%     bdry = imdilate(bdry,se);
%     bdry = logical(bdry);
%     seg_result = struct('seg',seg,'bdry',bdry');
%     parsave(fullfile(output_dir,fname),seg_result);
end
%}

%{
for m = 1:length(threshold_vec)
threshold = threshold_vec(m);
output_dir = ['D:\Documents\multiplex\seg_output' '_' num2str(threshold)];
if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

for i = 1:length(filelist)
    t1 = tic;
    fname = filelist{i};
    fname_split = strsplit(fname,{'_','.'});
    slide_name = strjoin(fname_split(1:3),'_');    
    pos = str2double(fname_split(4));
    if strcmp(slide_name,'AGTA_269_3')
        imname = fullfile(image_dir,slide_name,'VHE',['VHE' sprintf('%04d',pos) '.jpg']);
    else
        imname = fullfile(image_dir,slide_name,'VHE',['VHE' sprintf('%03d',pos) '.jpg']);
    end
    
    if ~exist(imname,'file')
        fprintf('image %s does not exist\n', fname)
        continue;
    end
    
    %if exist(fullfile(output_dir,fname),'file')
    %    continue;
    %    fprintf('Already done with image %s in %.2f seconds\n',fname, toc(t1));
    %end
    im = imread(imname);
    %figure; imshow(im);
    
    tmp = load(fullfile(data_dir,fname));
    x = tmp.x;
    y = tmp.y;
    epi_stroma = tmp.epithelial == 1;
    area = tmp.area;
    
    %figure; imshow(im);
    %hold on;
    %plot(x,y,'b.','MarkerSize',10);
    %hold off;
    
    % epithelial only
    x = double(x(epi_stroma));
    y = double(y(epi_stroma));
    area = double(area(epi_stroma));
    
    % find pairwise distances
    pairwise_dist = pdist(cat(1,x,y)');
    % filter below median
    indx_filter = pairwise_dist > threshold;
    pairwise_dist(indx_filter) = 0;
    
    pairwise_dist_sparse = sparse(squareform(pairwise_dist));
    %view(biograph(pairwise_dist_sparse));
    [~,C] = graphconncomp(pairwise_dist_sparse);
    
    components_areas = zeros(max(C),1);
    for k = 1:length(components_areas)
        indx_cl = C == k;
        components_areas(k) = sum(area(indx_cl));
    end
    
    [~,indx] = sort(components_areas,'descend');
    coordinates = cat(1,x,y)';
    % display output
    redc = im(:,:,1); greenc = im(:,:,2); bluec = im(:,:,3);
    nrow = size(im,1); ncol = size(im,2);
    mask_obj = zeros(nrow, ncol);
    indx_obj = sub2ind(size(mask_obj),y,x);
    %indx_obj = sub2ind(size(mask_obj),x,y);
    mask_obj(indx_obj) = 1;
    se = strel('disk',7,4);
    mask_obj = imdilate(mask_obj,se)>0;
    redc(mask_obj) = 255;
    greenc(mask_obj) = 140;
    bluec(mask_obj) = 0;
    
    % make sure the area is > 1000
    min_area = 500;
    num_comps = sum(components_areas > min_area);
    top_centers = cell(num_comps,1);
    top_radii = cell(num_comps,1);
    indx_sp_cl = cell(num_comps,1); % index of superpixel
    for k = 1:num_comps
        indx_cl = C == indx(k);
        top_centers{k} = coordinates(indx_cl,:);
        top_radii{k} = sqrt(area(indx_cl)/pi)';
    end
    
    %if plot_flag
    
    colors = distinguishable_colors(num_comps)*255;
    seg = zeros(nrow,ncol);
    for j = 1:num_comps
        x = top_centers{j}(:,1); y = top_centers{j}(:,2);
        if length(x) ==1
            ang=0:0.01:2*pi;
            xp=floor(top_radii{j}*cos(ang)'+x);
            yp=floor(top_radii{j}*sin(ang)'+y);
            k = boundary(xp,yp);
            mask = poly2mask(xp(k),yp(k),nrow, ncol);
        elseif length(x) <= 10
            extended_x = [x; min(x+top_radii{j},nrow); max(x - top_radii{j},0)];
            extended_y = [y; min(y+top_radii{j},nrow); max(y - top_radii{j},0)];
            k = boundary(extended_x,extended_y,0.5);
            mask = poly2mask(extended_x(k),extended_y(k),nrow, ncol);
        else
            k = boundary(x,y,0.5);
            mask =  poly2mask(x(k),y(k),nrow, ncol);
        end
        % avoid overlapping
        area_mask = sum(mask(:));
        if area_mask > 8000
            bdry = seg2bdry(double(mask),'imageSize');
            bdry = imdilate(bdry,se);
            bdry = logical(bdry);
            redc(bdry) = colors(j,1);
            greenc(bdry) = colors(j,2);
            bluec(bdry) = colors(j,3);
        end
        if sum(seg(mask)>0)/sum(mask(:)) < 0.3
            seg(mask) = j;
        else
            mask_diff = ((seg - mask) == -1);
            seg(mask_diff) = j;
        end
        %figure; imshow(label2rgb(seg));
    end
    I2 = cat(3,redc, greenc, bluec);
    imwrite(I2, fullfile(output_dir,[fname(1:end-4) '.jpg']));
    parsave(fullfile(output_dir,fname),seg);    
    fprintf('Done with image %s in %.2f seconds\n',fname, toc(t1));
end

end
%}
% calculate the semivariance

seg_dir = 'D:\Documents\multiplex\seg_output_50';
coordinate_dir = 'D:\Documents\multiplex\coordinates';
distance_shell = 5:10:105;
output_dir = 'D:\Documents\multiplex\semovariograms';
if ~exist(output_dir,'dir')
    mkdir(output_dir);
end
bm_names = {'Memb.CD8','Cyt.COX2','Cell.NaKATPase', 'Cell.MSH2',...
             'Cell.pERK','Cyt.BetaCatenin','Cyt.pMAPKAPK2','Nuc.p21',...
             'Memb.BetaCatenin','Nuc.CD8','Nuc.CD3','Cell.PI3Kp110a'};
         
colors = distinguishable_colors(length(bm_names));


%{
count = 0;
semivar = {};
output_names = {};
for i = 1:length(filelist)
    imname = filelist{i}(1:end-4);
    
    if ~exist(fullfile(seg_dir, [imname '.mat']),'file');
        fprintf('File %s does not exist\n',imname);
        continue;
    end
    
    tmp = load(fullfile(seg_dir, [imname '.mat']));
    seg = tmp.data;
    numSegs = max(seg(:));
    im = imread(fullfile(seg_dir,[imname '.jpg']));
    
    for j = 1:numSegs
       area_seg = sum(seg(:) == j);
       if area_seg < 35000
           continue;
       end
       mask = seg == j;
       imseg = im.*repmat(uint8(mask),[1 1 3]);
       if ~exist(fullfile(output_dir,[imname '_seg_' num2str(j) '.mat']),'file');
           fprintf('File %s does not exist\n',[imname '_seg_' num2str(j)]);
           continue;
       end
       tmp = load(fullfile(output_dir,[imname '_seg_' num2str(j) '.mat']));
       count = count + 1;
       output_names{count} = [imname '_seg_' num2str(j)];
       semivar{count} = mean(tmp.bm_diffs,1);
    end
end

semivar = cat(1,semivar{:});

for k = 2:1:6

[idx,C, sumd] = kmeans(semivar,k, 'Options',statset('UseParallel',1));
figure;
[s,h] = silhouette(semivar,idx);
hold on;
plot([mean(s) mean(s)], [0,k],'r-');
hold off;
fprintf('Clusters %d has mean s = %.2f\n',k, mean(s));
end

k = 4;
[idx,C, sumd] = kmeans(semivar,k, 'Options',statset('UseParallel',1));

for j = 1:max(idx)
   idx_cls = idx == j;
   image_cls = output_names(idx_cls);
   if ~exist([output_dir '_' num2str(j)],'dir')
       mkdir([output_dir '_' num2str(j)]);
   end
   for m = 1:length(image_cls)
       copyfile(fullfile(output_dir,[image_cls{m} '.png']),[output_dir '_' num2str(j)])
   end
end
%}
se = strel('disk',8,4);
nclusters = 10;
colors = uint8(distinguishable_colors(nclusters).*255);

for i = 1%:length(filelist)
    imname = filelist{i}(1:end-4);
    
    fname_split = strsplit(imname,{'_','.'});
    slide_name = strjoin(fname_split(1:3),'_');    
    pos = str2double(fname_split(4));
    if strcmp(slide_name,'AGTA_269_3')
        he_imname = fullfile(image_dir,slide_name,'VHE',['VHE' sprintf('%04d',pos) '.jpg']);
    else
        he_imname = fullfile(image_dir,slide_name,'VHE',['VHE' sprintf('%03d',pos) '.jpg']);
    end
    
    if ~exist(he_imname,'file')
        fprintf('image %s does not exist\n', imname)
        continue;
    end
    
    %if exist(fullfile(output_dir,fname),'file')
    %    continue;
    %    fprintf('Already done with image %s in %.2f seconds\n',fname, toc(t1));
    %end
    im = imread(he_imname);
    
    tmp = load(fullfile(seg_dir, [imname '.mat']));
    seg = tmp.data;
    numSegs = max(seg(:));
    
    %im = imread(fullfile(seg_dir,[imname '.jpg']));
    tmp = load(fullfile(coordinate_dir, [imname '.mat']));
    x = tmp.x;
    y = tmp.y;
    
    indx_cells = sub2ind(size(seg),y,x);
    epi_stroma = tmp.epithelial == 1;
    area = tmp.area;
    bm_data = tmp.bm_data;
    
    for j = 1%:numSegs
       area_seg = sum(seg(:) == j);
       if area_seg < 15000
           continue;
       end
       mask = seg == j;
       imseg = im.*repmat(uint8(mask),[1 1 3]);
       cells_in_mask = mask(indx_cells);
       coord_mask = [x(cells_in_mask); y(cells_in_mask)]';
       bm_mask = bm_data(cells_in_mask,:);
       area_mask = area(cells_in_mask);
       
       pairwise_distances = pdist(double(coord_mask));
       max_dist = max(pairwise_distances);
       pairwise_distances = squareform(pairwise_distances);
       threshold = 100;
       cell_variances = zeros(size(pairwise_distances,1),length(bm_names));
       for k = 1:size(pairwise_distances, 1)
          cell_in_shell = find(pairwise_distances(k,:) <= threshold);
          bm_in_shell = bm_mask(cell_in_shell,:);
          cell_variances(k,1:length(bm_names)) = std(bm_in_shell, 1);  
          %cell_variances(k,1:length(bm_names)) = std(bm_in_shell, 1); 
       end
       % do hierachical clustering
       %T = clusterdata(cell_variances, 1);
       indx_keep = v(:,1) < -.02;
       data_keep = cell_variances(indx_keep,:);
       [u1,s1,v1] = svd(data_keep','econ');
       figure; plot(v1(:,1), v1(:,2),'.');
       [indx,C] = kmeans(data_keep,nclusters);
       coord_mask = coord_mask(indx_keep);

       % display the cells
       
       redc = imseg(:,:,1); greenc = imseg(:,:,2); bluec = imseg(:,:,3);
       nrow = size(im,1); ncol = size(im,2);
       for k = 1:max(indx)
          indx_cls = find(indx == k);
          % display output
          cls_x = coord_mask(indx_cls,1);
          cls_y = coord_mask(indx_cls,2);    
          mask_obj = zeros(nrow, ncol);
          indx_obj = sub2ind(size(mask_obj),cls_y,cls_x);
          
          mask_obj(indx_obj) = 1;
          mask_obj = imdilate(mask_obj,se)>0;
          redc(mask_obj) = colors(k,1);
          greenc(mask_obj) = colors(k,2);
          bluec(mask_obj) = colors(k,3);
       end       
       output_im = uint8(cat(3,redc, greenc, bluec));
       figure; imshow(output_im);
    end
end

distance_to_centers = pdist2(cell_variances, C);
[~,indx] = min(distance_to_centers,[],2);

indx_keep = v(:,1) < -.02;
data_keep = cell_variances(indx_keep,:);
[u1,s1,v1] = svd(data_keep','econ');
figure; plot(v1(:,1), v1(:,2),'.');

data_left = cell_variances(~indx_keep,:);
[u2,s2,v2] = svd(data_left','econ');
figure; plot(v2(:,1), v2(:,2),'.');


%{
for i = 1:length(filelist)
    imname = filelist{i}(1:end-4);
    tmp = load(fullfile(seg_dir, [imname '.mat']));
    seg = tmp.data;
    numSegs = max(seg(:));
    
    im = imread(fullfile(seg_dir,[imname '.jpg']));
    tmp = load(fullfile(coordinate_dir, [imname '.mat']));
    x = tmp.x;
    y = tmp.y;
    
    indx_cells = sub2ind(size(seg),y,x);
    epi_stroma = tmp.epithelial == 1;
    area = tmp.area;
    bm_data = tmp.bm_data;
    
    for j = 1:numSegs
       area_seg = sum(seg(:) == j);
       if area_seg < 15000
           continue;
       end
       mask = seg == j;
       imseg = im.*repmat(uint8(mask),[1 1 3]);
       cells_in_mask = mask(indx_cells);
       coord_mask = [x(cells_in_mask); y(cells_in_mask)]';
       bm_mask = bm_data(cells_in_mask,:);
       area_mask = area(cells_in_mask);
       
       pairwise_distances = pdist(double(coord_mask));
       max_dist = max(pairwise_distances);
       pairwise_distances = squareform(pairwise_distances);
       
       curr_dist_shell = distance_shell(distance_shell < max_dist);
       
       bm_diffs = cell(length(curr_dist_shell) - 1,1);
       bm_se = cell(length(curr_dist_shell) - 1,1);
       omit_dist = [];
       for k = 1:length(curr_dist_shell) -1
           [r,c] = find(pairwise_distances >= curr_dist_shell(k) & ...
               pairwise_distances < curr_dist_shell(k+1));
           if length(r) < 10
               omit_dist = cat(1,omit_dist,[k+1]);
               continue;
           end
           cells_in_shell = [r,c];
           cells_in_shell = cells_in_shell(r < c,:);
           curr_diff = bm_mask(cells_in_shell(:,1),:) - bm_mask(cells_in_shell(:,2),:);
           bm_diffs{k} = mean((curr_diff).^2,1);
           bm_se{k} = std(curr_diff,1)/sqrt(size(curr_diff,1));
       end
       curr_dist_shell(omit_dist) = [];
       bm_diffs = cat(1,bm_diffs{:});
       bm_se = cat(1,bm_se{:});
       if isempty(bm_diffs); continue;end
       figure('Position', [100, 100, 800, 800]); 
       subplot(3,3,[1 2 4 5])
       imshow(imseg);% hold on; plot(coord_mask(:,1), coord_mask(:,2),'.'); hold off;
       subplot(3,3,[7 8 9]);
       for m = 1:length(bm_names)
          errorbar(curr_dist_shell(2:end),bm_diffs(:,m), bm_se(:,m),...
              'Color',colors(m,:), 'LineWidth',3); hold on
       end
       hold off
       legend(bm_names{:},'Location','bestoutside');%legend('boxoff')
       xlim([0, max(distance_shell)]); set(gca,'FontSize',12);
       xlabel('Distances'); ylabel('semivariance');
       
       print(fullfile(output_dir,[imname '_seg_' num2str(j)]) ,'-dpng');
       close all
       save(fullfile(output_dir,[imname '_seg_' num2str(j) '.mat']),'curr_dist_shell','bm_diffs','bm_se');
    end
end
%}

