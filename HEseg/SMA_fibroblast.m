%% script to examine SMA near and faraway from tumor cells
% Assumption is that epithelial cells are from the tumor areas
% April 17, 2017

data_dir = 'D:\Documents\multiplex';
coordinate_dir = 'D:\Documents\multiplex\coordinates_all_bm';
image_dir = 'D:\Documents\colon_cancer_data\H&E_virtual';

filelist = dir(fullfile(coordinate_dir,'*.mat'));
filelist = {filelist.name}';
seg_dir = 'D:\Documents\multiplex\seg_output_50';
%output_dir = 'D:\Documents\multiplex\SMA_stroma_50';

output_dir = 'D:\Documents\multiplex\Immuno_bm';
if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

% find the index of SMA
tmp = load(fullfile(coordinate_dir,filelist{1}));
bm_names = cellstr(tmp.bm_names);


immune_bm_names = {'Median.Nuc.CD20','Median.Nuc.CD8',...
    'Median.Nuc.CD79','Median.Nuc.CD3','Median.Nuc.EPCAM'};
indx = find(ismember(bm_names,immune_bm_names));

%filelist = {'AGA_260_3_35.mat','AGA_260_3_46.mat','AGA_260_3_71.mat',...
%    'AGA_260_3_106.mat','AGA_260_3_110.mat','AGA_260_3_166.mat',...
%    'AGA_260_3_35.mat','AGA_260_3_100.mat','AGA_260_3_149.mat',...
%    'AGA_260_3_162.mat','AGA_260_3_191.mat','AGA_260_3_199.mat','AGA_260_3_217.mat'};


for i = 1:length(filelist)
    imname = filelist{i}(1:end-4);
    %tmp = load(fullfile(seg_dir, [imname '.mat']));
    %seg_result = tmp.seg_result;
    %seg = seg_result.seg;
    %numSegs = max(seg(:));
    
    if ~exist(fullfile(seg_dir, [imname '.mat']),'file')
        fprintf('Done have this file: %s\n',imname);
        continue;
    end
    
    im = imread(fullfile(seg_dir,[imname '.jpg']));
    
    tmp = load(fullfile(coordinate_dir, [imname '.mat']));
    x = max(1,tmp.x);
    y = max(1,tmp.y);
    
    %indx_cells = sub2ind(size(seg),y,x);
    epi_stroma = tmp.epithelial == 1;
    area = tmp.area;
    % stroma, val = 2
    bm_data = tmp.bm_data(~epi_stroma,indx);
    num_cells = size(tmp.bm_data,1);
    
    [counts1, binCenters1] = histcounts(bm_data(:,1),16);
    [counts2, ~] = histcounts(bm_data(:,2), binCenters1);
    [counts3, ~] = histcounts(bm_data(:,3), binCenters1);
    [counts4, ~] = histcounts(bm_data(:,4), binCenters1);
    [counts5, ~] = histcounts(bm_data(:,5), binCenters1);
    
    figure('Position', [100, 100, 1500, 500]);
    subplot(1,3,1); imshow(im);
    subplot(1,3,3);
    plot(binCenters1(2:end), counts1/num_cells, 'r-','LineWidth',3);
    hold on;
    plot(binCenters1(2:end), counts2/num_cells, 'b-','LineWidth',3);
    plot(binCenters1(2:end), counts3/num_cells, 'r-v','LineWidth',3);
    plot(binCenters1(2:end), counts4/num_cells, 'b-o','LineWidth',3);
    plot(binCenters1(2:end), counts5/num_cells, 'g-x','LineWidth',3);
    legend({'CD20','CD8','CD79','CD3','EpCAM'})
    hold off;
    ylim([0 1]); ylabel('Probability')
    xlim([0 10]); xlabel('Biomarker Exp. (log scale)')
    title('Stroma');set(gca,'FontSize',16);axis square;
    
    %print(fullfile(output_dir,[imname '_stroma']) ,'-dpng');
    
    % epithelial, val = 1
    bm_data = tmp.bm_data(epi_stroma,indx);
    num_cells = size(tmp.bm_data,1);
    
    [counts1, binCenters1] = histcounts(bm_data(:,1),16);
    [counts2, ~] = histcounts(bm_data(:,2), binCenters1);
    [counts3, ~] = histcounts(bm_data(:,3), binCenters1);
    [counts4, ~] = histcounts(bm_data(:,4), binCenters1);
    [counts5, ~] = histcounts(bm_data(:,5), binCenters1);
    
    %figure;%('Position', [100, 100, 800, 800]);
    subplot(1,3,2);
    plot(binCenters1(2:end), counts1/num_cells, 'r-','LineWidth',3);
    hold on;
    plot(binCenters1(2:end), counts2/num_cells, 'b-','LineWidth',3);
    plot(binCenters1(2:end), counts3/num_cells, 'r-v','LineWidth',3);
    plot(binCenters1(2:end), counts4/num_cells, 'b-o','LineWidth',3);
    plot(binCenters1(2:end), counts5/num_cells, 'g-x','LineWidth',3);
    
    legend({'CD20','CD8','CD79','CD3','EpCAM'})
    hold off;
    ylim([0 1]); ylabel('Probability')
    xlim([0 10]); xlabel('Biomarker Exp. (log scale)')
    title('Epithelial'); set(gca,'FontSize',16); axis square;
    
    print(fullfile(output_dir,imname) ,'-dpng');
    close all;
end



indx = find(ismember(bm_names,'Median.Cell.SMA'));

r = 75;
se = strel('disk',r,8);
summary_stroma_cells = zeros(length(filelist),4);

%output_dir = ['D:\Documents\multiplex\SMA_stroma_', num2str(r)];

%output_dir = 'D:\Documents\multiplex\Immuno_bm';
if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

clinical_data = readtable(fullfile(data_dir,'clinical_data_all_spots.csv'),...
    'Delimiter',',');

%{
parfor i = 1:length(clinical_data.spot_name)
   imname = clinical_data.spot_name{i};
   if ~exist(fullfile(output_dir,[imname '.png']),'file')
       continue;
   end
   
   curr_stage = clinical_data.recurrent_5yr(i);
   if ~exist(fullfile(output_dir,['recurr' num2str(curr_stage)]),'dir')
       mkdir(fullfile(output_dir,['recurr' num2str(curr_stage)]))
   end
   
   copyfile(fullfile(output_dir,[imname '.png']),fullfile(output_dir,['recurr' num2str(curr_stage)]));
end
%}
%{
for i = 1:length(filelist)
    imname = filelist{i}(1:end-4);
    if ~exist(fullfile(seg_dir, [imname '.mat']),'file')
        fprintf('Done have this file: %s\n',imname);
        continue;
    end
    
    %if exist(fullfile(output_dir,[imname '.png']),'file')
    %    continue;
    %end
    
    tmp = load(fullfile(seg_dir, [imname '.mat']));
    seg_result = tmp.seg_result;
    seg = seg_result.seg;
    numSegs = max(seg(:));
    
    im = imread(fullfile(seg_dir,[imname '.jpg']));
    tmp = load(fullfile(coordinate_dir, [imname '.mat']));
    x = max(1,tmp.x);
    y = max(1,tmp.y);
    
    indx_cells = sub2ind(size(seg),y,x);
    epi_stroma = tmp.epithelial == 1;
    area = tmp.area;
    bm_data = tmp.bm_data;
    
    mask = seg > 0;
    dilated_mask = imdilate(mask,se);
    
    imseg = im.*repmat(uint8(mask),[1 1 3]);
    dilated_imseg = im.*repmat(uint8(dilated_mask),[1 1 3]);
    
    %figure(1);
    figure('Position', [100, 100, 1200, 600]); 
    subplot(1,2,1);
    imshow(dilated_imseg);
    
    cells_in_mask = dilated_mask(indx_cells);
    coord_mask = [x(cells_in_mask); y(cells_in_mask)]';
    
    stroma_in_mask = cells_in_mask & (~epi_stroma);
    stroma_out_mask = (~cells_in_mask) & (~epi_stroma);
    %fprintf('There are %d stroma cells near tumor cells with mean = %.3f\n',...
    %    sum(stroma_in_mask),mean(bm_data(stroma_in_mask,indx)));
    
    %fprintf('There are %d stroma cells outside tumors with mean = %.3f\n',...
    %    sum(stroma_out_mask),mean(bm_data(stroma_out_mask,indx)));
    
    summary_stroma_cells(i,1) = sum(stroma_in_mask);
    summary_stroma_cells(i,2) = mean(bm_data(stroma_in_mask,indx));
    summary_stroma_cells(i,3) = sum(stroma_out_mask);
    summary_stroma_cells(i,4) = mean(bm_data(stroma_out_mask,indx));
    
    [counts1, binCenters1] = histcounts(bm_data(stroma_in_mask, indx),0:1:16);
    [counts2, ~] = histcounts(bm_data(stroma_out_mask,indx), binCenters1);
    [counts3, ~] = histcounts(bm_data(epi_stroma,indx), binCenters1);
    num_stroma_cells = length(~epi_stroma);
    
    %figure(2);
    subplot(1,2,2);
    semilogy(binCenters1(2:end), counts1/num_stroma_cells + eps, 'r-','LineWidth',3);
    hold on;
    semilogy(binCenters1(2:end), counts2/num_stroma_cells + eps, 'b-','LineWidth',3);
    semilogy(binCenters1(2:end), counts3/num_stroma_cells + eps, 'g-','LineWidth',3);
    legend('near tumors','far outside','epithelial','Location','Best')
    axis square
    hold off; set(gca,'FontSize',16);
    ylim([0 1]); ylabel('Probability')
    xlim([0 16]); xlabel('SMA (log scale)')
    
    print(fullfile(output_dir,imname) ,'-dpng');
    close all;
end

save(fullfile(output_dir,'summary.mat'),'summary_stroma_cells');
%}