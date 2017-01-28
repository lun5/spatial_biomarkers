%% function to output the communication indices for different setting
% input: NN_dir, nearest neighbor directory, threshold for residual
% output: histogram plot, count neighbor type for each spot
function [count_type,count_type_prop, grpMean,count_numel] = calculateCommIndex(NN_dir, output_dir, stageT, thres)

tmp = load(fullfile(NN_dir,'all_spots.mat'));
num_nb = prctile(tmp.nn,98); % take 98 percentile 
count_type = zeros(length(stageT.spot_name),num_nb);

for ff = 1:length(stageT.spot_name)
    spot_name = stageT.spot_name{ff};
    tmp = load(fullfile(NN_dir,[spot_name '.mat']));%,'residuals','num_nn');
    cell_id = tmp.residuals < thres & tmp.num_nuclei == 1 & ...
        tmp.areas > 50 & tmp.areas < 2000; % only epithelial
    num_nn = tmp.num_nn(cell_id);
    for i = 1:num_nb
        if i < num_nb
            count_type(ff,i) = sum(num_nn == i);
        else
            count_type(ff,i) = sum(num_nn >= i);
        end
    end   
end

%% proportion of count type/spot
count_type_prop = count_type./repmat(sum(count_type,2),[1 size(count_type,2)]);
grpMean = grpstats(count_type_prop,stageT.stage);
figure; bar(grpMean');
xlabel('Communication index');
set(gca,'FontSize',16);
legend('Stage 1','Stage 2','Stage 3');
ylabel('Probability');
ylim([0 ceil(max(grpMean(:))*10)/10]);
xlim([0 num_nb + 0.5]);
set(gcf,'color','white');
set(gcf,'PaperPositionMode','auto')
name_split = strsplit(strrep(NN_dir,'.',''),filesep);
print(fullfile(output_dir,[name_split{end}, '_thres', strrep(num2str(thres),'.',''), '_spot']),'-dtiff','-r300');

%% group all cells together, regardless of spot
num_stages = max(stageT.stage);
count_numel = zeros(num_stages,num_nb);
for i =1:num_stages
    count_numel(i,:) = sum(count_type(stageT.stage == i,:),1);
end
prop_numel = count_numel./repmat(sum(count_numel,2),[1 num_nb]);
figure; bar(prop_numel');
xlabel('Communication index');
set(gca,'FontSize',16);
legend('Stage 1','Stage 2','Stage 3');
ylabel('Probability');
ylim([0 ceil(max(grpMean(:))*10)/10]);
xlim([0 num_nb + 0.5]);
set(gcf,'color','white');
set(gcf,'PaperPositionMode','auto')
print(fullfile(output_dir,[name_split{end}, '_thres', strrep(num2str(thres),'.','')]),'-dtiff','-r300');
close all
save(fullfile(output_dir,[name_split{end},'_thres', strrep(num2str(thres),'.',''), '.mat']),...
    'count_type','count_numel','grpMean','count_type_prop');
end