% imdir = 'X:\ColonCancerStudy\AGA_260_3\AFRemoved';
%bm = 'CK19';
%bm = 'CD68';
% bm = 'ALDH1';
imlist = dir(fullfile(imdir,[bm '*.tif']));
imlist = {imlist.name}';

fprintf('There are %d image with %s biomarker.\n', length(imlist), bm);

% outputdir = fullfile('X:\ColonCancerStudy\AGA_260_3\',bm);
% if ~exist(outputdir,'dir'); mkdir(outputdir); end

for i = 1:20%length(imlist)
    im = imread(fullfile(imdir,imlist{i}));
    %figure(1); imagesc(log2(double(im)));axis off; axis equal; colorbar; caxis([0,16])
    imtool(im); 
    figure(2); plot(1:size(im,1), sum(im,1));
    figure(3); plot(1:size(im,1), sum(im,2));
    pause;
    
end 

slide_name = 'AGTA_264_3';
pos = 129;

imdir = 'X:\ColonCancerStudy';
im = imread(fullfile(imdir,slide_name, 'AFRemoved',['_AFRemoved_' , num2str(pos),'.tif']));
imad = imadjust(im);
figure; imshow(imad);

tmp = load(fullfile(imdir,[slide_name '_' num2str(pos) '.mat']));
% parse the biomarker names
bm_name = tmp.biomarkers(:,12:end);

[max_resid, indx_max] = max(tmp.residual(:,1));
[min_resid, indx_min] = min(tmp.residual(:,1));

figure(1); imshow(imad)
hold on; plot(tmp.xy(:,1), tmp.xy(:,2),'bo');
hold off

figure(2); imshow(imad);
hold on; plot(tmp.xy(tmp.neighbors(indx_max,:)+1,1),...
    tmp.xy(tmp.neighbors(indx_max,:)+1,2),'bx');
hold off

figure(3); imshow(imad);
hold on; plot(tmp.xy(tmp.neighbors(indx_min,:)+1,1),...
    tmp.xy(tmp.neighbors(indx_min,:)+1,2),'gx');
hold off

%tmp.biomarkers = tmp.biomarkers(:,12:end);
% show the cells in biomarkers with highest/lowest residuals
% for the max
bm_max = strtrim(bm_name(tmp.residual(indx_max,3),:));
bm_im_max = imread(fullfile(imdir,slide_name, 'AFRemoved',[bm_max, '_AFRemoved_' , num2str(pos),'.tif']));
bm_im_max = imadjust(bm_im_max);

figure(4); imshow(bm_im_max);
hold on; plot(tmp.xy(tmp.neighbors(indx_max,:)+1,1),...
    tmp.xy(tmp.neighbors(indx_max,:)+1,2),'bx');
hold off

bm_min = strtrim(bm_name(tmp.residual(indx_max,5),:));
bm_im_min = imread(fullfile(imdir,slide_name, 'AFRemoved',[bm_min, '_AFRemoved_' , num2str(pos),'.tif']));
bm_im_min = imadjust(bm_im_min);

figure(5); imshow(bm_im_min);
hold on; plot(tmp.xy(tmp.neighbors(indx_max,:)+1,1),...
    tmp.xy(tmp.neighbors(indx_max,:)+1,2),'gx');
hold off

% for the min
bm_max = strtrim(bm_name(tmp.residual(indx_min,3),:));
bm_im_max = imread(fullfile(imdir,slide_name, 'AFRemoved',[bm_max, '_AFRemoved_' , num2str(pos),'.tif']));
bm_im_max = imadjust(bm_im_max);

figure(6); imshow(bm_im_max);
hold on; plot(tmp.xy(tmp.neighbors(indx_min,:)+1,1),...
    tmp.xy(tmp.neighbors(indx_min,:)+1,2),'bx');
hold off

bm_min = strtrim(bm_name(tmp.residual(indx_min,5),:));
bm_im_min = imread(fullfile(imdir,slide_name, 'AFRemoved',[bm_min, '_AFRemoved_' , num2str(pos),'.tif']));
bm_im_min = imadjust(bm_im_min);

figure(7); imshow(bm_im_min);
hold on; plot(tmp.xy(tmp.neighbors(indx_min,:)+1,1),...
    tmp.xy(tmp.neighbors(indx_min,:)+1,2),'gx');
hold off

%% draw circles based on numbers of 
slide_name = 'AGTA_264_3';
pos = 231;
NN_dir = fullfile(imdir, 'NN_radius_100_withRes');
NN_dir = fullfile(imdir, 'NN_data', 'Elastic_NN_radius_50_alpha_1_L1ratio_0.5');

filelist = dir(fullfile(NN_dir,'AG*.mat'));
filelist = {filelist.name}';
%colors = {'yellow', 'magenta','cyan','red', 'green'};
colors = {'blue','cyan','yellow','magenta','red'};
r = 30;

% viz_dir = fullfile(imdir, 'NN_data','NN_radius_50_alpha_1_viz');
% if ~exist(viz_dir,'dir'); mkdir(viz_dir);end
% 
stageT = readtable(fullfile(imdir,'spots_stages.csv'),'Delimiter',',','ReadVariableNames',false);
stageT.Properties.VariableNames = {'spot_name','stage'};

thres = 0.15;
num_nb = 12;
count_type = zeros(length(filelist),num_nb);

for ff = f1:length(stageT.spot_name)
    spot_name = stageT.spot_name{ff};
    split_name = strsplit(spot_name,'_');
    slide_name = strjoin(split_name(1:3),'_');
    pos = split_name{4};
    tmp = load(fullfile(NN_dir,[spot_name '.mat']));%,'residuals','num_nn');
    %indx = tmp.residuals < thres;
    %xy = tmp.xy(tmp.residuals < thres,:);
    % eliminate cells based on areas, number of nuclei
    cell_id = tmp.residuals < thres & tmp.num_nuclei == 1 & ...
        tmp.areas > 50 & tmp.areas < 2000; % only epithelial
    num_nn = tmp.num_nn(cell_id);
    xy = tmp.xy(cell_id,:);
    %num_nn = ceil( tmp.num_nn(cell_id)/2);
    
    %BetaCatenin
    im = imread(fullfile(imdir,slide_name, 'AFRemoved',['BetaCatenin_AFRemoved_' ,...
        sprintf('%3.3d',str2num(pos)),'.tif']));
    imad = imadjust(im);
%     
    xy_stroma = tmp.xy(tmp.epi_stroma == 2,:);
    xy_epi = tmp.xy(tmp.epi_stroma == 1,:);
    figure; imshow(imad)
    hold on
    viscircles(xy_stroma,ones(size(xy_stroma,1),1)*10,'Color','b');
    viscircles(xy_epi,10*ones(size(xy_epi,1),1),'Color','y');
    hold off; 
        
    
%         
%     figure; imshow(imad);
%     pause
%     
%     hold on
%     for i = 1%:num_nb
%         
% %         if i < num_nb
% %             centers = xy(num_nn == i,:);
% %             count_type(ff,i) = sum(num_nn == i);
% %         else
% %             count_type(ff,i) = sum(num_nn >= i);
% %             centers = xy(num_nn >= i,:);
% %         end
%         
%         neighbor_id = tmp.neighbor_id(num_nn == i);
%         neighbor_id = cat(2,neighbor_id{:});
%         neighbor_centers = tmp.xy(neighbor_id+1,:);
%         pos_rect_neigh = cat(2,neighbor_centers(:,1)-r, neighbor_centers(:,2)-r,...
%             r*ones(size(neighbor_centers(:,1))),r*ones(size(neighbor_centers(:,1))));
%         pos_rect = cat(2,centers(:,1) - r,centers(:,2) - r,...
%             r*ones(size(centers(:,1))),r*ones(size(centers(:,1))));
%         for j = 1:size(pos_rect,1)            
%             rectangle('Position',pos_rect(j,:),'Curvature',[1 1], 'FaceColor',colors{i});
%         end                
%         for j = 1:size(neighbor_centers,1)
%             rectangle('Position',pos_rect_neigh(j,:),'Curvature',[1 1], 'FaceColor','y');
%         end 
%     end
    
%     hold off
    
%     f=getframe(gca);
%     [X, map] = frame2im(f);
%     imwrite(X,fullfile(viz_dir,[fname '.tif']));
%     close all;
end


count_type_prop = count_type./repmat(sum(count_type,2),[1 size(count_type,2)]);
[grpMin,grpMed,grpMean, grpStd, grpMax, grpNumel] = grpstats(count_type_prop,stageT.stage, ...
{'min','median','mean','std','max','numel'});
figure; bar(grpMean');
xlabel('Communication index');
set(gca,'FontSize',16);
legend('Stage 1','Stage 2','Stage 3');
ylabel('Probability');
ylim([0 0.2]);xlim([0.5 12.5]);

stage1 = sum(count_type(stageT.stage == 1,:),1);
stage2 = sum(count_type(stageT.stage == 2,:),1);
stage3 = sum(count_type(stageT.stage == 3,:),1);
count_numel = [stage1;stage2;stage3];
prop_numel = count_numel./repmat(sum(count_numel,2),[1 5]);
figure; bar(prop_numel');
xlabel('Communication index');
set(gca,'FontSize',16);
legend('Stage 1','Stage 2','Stage 3');
ylabel('Probability');
ylim([0 0.2]);xlim([0.5 12.5]);



[grpMin,grpMed,grpMean, grpStd, grpMax, grpNumel] = grpstats(count_type_prop,grps, ...
{'min','median','mean','std','max','numel'});
figure; bar(grpMean');
xlabel('Number of dependent neighbors');
set(gca,'FontSize',16);
legend('Stage 1','Stage 2+3','Stage 3');
ylabel('Average proportion of cells/spot');

% There isn't any significant different between 3 group
[grpMin_count,grpMed_count,grpMean_count, grpStd_count, grpMax_count, grpNumel_count] =...
    grpstats(count_type,stageT.stage, {'min','median','mean','std','max','numel'});

entr = sum(count_type_prop.*log2(count_type_prop),2);
grpMean_entr = grpstats(entr, stageT.stage);

for i =1 :num_nb
    [p,tbl,stats] = anova1(count_type_prop(:,i), stageT.stage);
    [c,m] = multcompare(stats);
end

for i =1:num_nb
    [p,tbl,stats] = anova1(count_type_prop(:,i), grps);
    [c,m] = multcompare(stats);
end

for i =1:num_nb
    [p,tbl,stats] = anova1(prop_numel(:,i), grps);
    [c,m] = multcompare(stats);
end

%% Loop through all possible combination (when I can get into idli) 
imdir = 'X:\ColonCancerStudy\';
stageT = readtable(fullfile(imdir,'spots_stages.csv'),'Delimiter',',','ReadVariableNames',false);
stageT.Properties.VariableNames = {'spot_name','stage'};

nb_radii = [25,50,75,100];
alphas = [0.1,1];
L1_ratios = [0, 0.3, 0.5, 1];
thresholds = [0.15, 0.2];

all_combs = combvec(nb_radii, alphas, L1_ratios, thresholds);
output_dir = fullfile(imdir,'NN_data','comm_index');
if ~exist(output_dir,'dir'); mkdir(output_dir,'dir'); end
for i = 1:length(all_combs)
   T1 = tic;
   radius = all_combs(1,i);
   alpha = all_combs(2,i);
   L1_ratio = all_combs(3,i);
   thres = all_combs(4,i);
   
   NN_dir = fullfile(imdir, 'NN_data',[ 'Elastic_NN_radius_' num2str(radius),...
       '_alpha_', num2str(alpha), '_L1ratio_', num2str(L1_ratio)]);
   if ~exist(NN_dir,'dir');
       fprintf('No data for r = %d, alpha = %.1f, L1_rat = %.1f, thres = %.2f combo\n',...
           radius, alpha, L1_ratio, thres);
       continue;
   end
   name_split = strsplit(strrep(NN_dir,'.',''),filesep);
   matfile_name = [name_split{end},'_thres', strrep(num2str(thres),'.',''), '.mat'];
   if exist(fullfile(output_dir, matfile_name), 'file')
       fprintf('already done with r = %d, alpha = %.1f, L1_rat = %.1f, thres = %.2f combo\n',...
       radius, alpha, L1_ratio, thres);
       continue;
   end
   [count_type,count_type_prop, grpMean,count_numel] = ...
       calculateCommIndex(NN_dir, output_dir, stageT, thres);
   fprintf('Done with r = %d, alpha = %.1f, L1_rat = %.1f, thres = %.2f combo in %.2f seconds\n',...
       radius, alpha, L1_ratio, thres, toc(T1));
end

%NN_dir = 'X:\ColonCancerStudy\NN_data\Elastic_NN_radius_50_alpha_1_l1_ratio_0.5_scramble';
%NN_dir = 'X:\ColonCancerStudy\NN_data\Elastic_NN_radius_50_alpha_0.1_l1_ratio_1_scramble';
NN_dir = 'X:\ColonCancerStudy\NN_data\Elastic_NN_radius_50_alpha_1_l1_ratio_0.5_scramble_onlypos';
%NN_dir = 'X:\ColonCancerStudy\NN_data\Elastic_NN_radius_50_alpha_0.1_L1ratio_1';
name_split = strsplit(strrep(NN_dir,'.',''),filesep);
matfile_name = [name_split{end},'_thres', strrep(num2str(thres),'.',''), '.mat'];
thres = 0.1;
[count_type,count_type_prop, grpMean,count_numel] = ...
       calculateCommIndex(NN_dir, NN_dir, stageT, thres);

