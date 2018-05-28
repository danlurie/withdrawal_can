function res = network_analysis_nothresh(mat_path, out_prefix, gamma)

%% Load BCT, NCT, and GenLouvain
addpath('/home/despoB/dlurie/Software/matlab/BCT/2017_01_15_BCT')
addpath('/home/despoB/dlurie/Software/matlab/NCT/2017_10_22_NCT')
addpath('/home/despoB/dlurie/Software/matlab/GenLouvain')
addpath('/home/despoB/dlurie/Projects/withdrawal_CAN/code')

%% Load the adjacency matrix
corr_tab = readtable(mat_path, 'ReadVariableNames', false, 'Delimiter', ' ');
corr_mat = table2array(corr_tab);
corr_mat = weight_conversion(corr_mat, 'autofix');

%% Run the Louvain algorithm 250 times
n_nodes = size(corr_mat,1);
n_iters = 250;
part_mat = zeros(n_nodes,n_iters);
for idx = 1:n_iters;
        %iter_report = sprintf('Iteration %d of 250.', idx);
        %disp(iter_report);
        n  = n_nodes;    % number of nodes
        M  = 1:n;                   % initial community affiliations
        Q0 = -1; Q1 = 0;            % initialize modularity values
        while Q1-Q0>1e-5;           % while modularity increases
            Q0 = Q1;                % perform community detection
            [M, Q1] = community_louvain(corr_mat, gamma, M, 'negative_asym');
        end   
        part_mat(:,idx) = M;
end

%% Generate consensus partition
[S2, Q2, X_new3, qpc] = consensus_iterative(part_mat.');
consensus_partition = S2(1,:);
consensus_q = Q2(1)

%% Save consensus partition
cpart_path = strcat(out_prefix, '_GraphPartition.txt');
fileID = fopen(cpart_path, 'w');
fprintf(fileID, '%i \r', consensus_partition);
fclose(fileID);

%% Save consensus Q value
modq_path = strcat(out_prefix, '_ModularityQ.txt');
fileID = fopen(modq_path, 'w');
fprintf(fileID, '%d', consensus_q);
fclose(fileID);

%% Save correlation matrix re-ordered by module
[corr_mat_order, corr_mat_reordered] = reorder_mod(corr_mat, consensus_partition);
reorder_path = strcat(out_prefix, '_ReorderedMatrix.txt');
fileID = fopen(reorder_path, 'w');
fprintf(fileID, '%.15e \r', corr_mat_reordered);
fclose(fileID);

%% Save the labels of re-ordered nodes.
reorder_order_path = strcat(out_prefix, '_ReorderedMatrix_order.txt');
fileID = fopen(reorder_order_path, 'w');
fprintf(fileID, '%.15e \r', corr_mat_order);
fclose(fileID);

%% Calculate nodal metrics
[z_all, z_pos, z_neg] = module_degree_zscore_sign(corr_mat, consensus_partition);
[pc_pos, pc_neg] = participation_coef_sign(corr_mat, consensus_partition);

%% Save Within Module Degree Z-scores
wmdz_all_path = strcat(out_prefix, '_WMDz_all.txt');
fileID = fopen(wmdz_all_path, 'w');
fprintf(fileID, '%.15e \r', z_all);
fclose(fileID);
wmdz_pos_path = strcat(out_prefix, '_WMDz_pos.txt');
fileID = fopen(wmdz_pos_path, 'w');
fprintf(fileID, '%.15e \r', z_pos);
fclose(fileID);
wmdz_neg_path = strcat(out_prefix, '_WMDz_neg.txt');
fileID = fopen(wmdz_neg_path, 'w');
fprintf(fileID, '%.15e \r', z_neg);
fclose(fileID);

%% Save Participation Coefficients
pc_pos_path = strcat(out_prefix, '_PC_pos.txt');
fileID = fopen(pc_pos_path, 'w');
fprintf(fileID, '%.15e \r', pc_pos);
fclose(fileID);
pc_neg_path = strcat(out_prefix, '_PC_neg.txt');
fileID = fopen(pc_neg_path, 'w');
fprintf(fileID, '%.15e \r', pc_neg);
fclose(fileID);

end