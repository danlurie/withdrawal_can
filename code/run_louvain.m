function res = run_louvain(mat_path, out_prefix)

%% Load BCT, NCT, and GenLouvain
addpath('/home/despoB/dlurie/Software/matlab/BCT/2017_01_15_BCT')
addpath('/home/despoB/dlurie/Software/matlab/NCT/2017_10_22_NCT')
addpath('/home/despoB/dlurie/Software/matlab/GenLouvain')
addpath('/home/despoB/dlurie/Software/matlab/versatility/')

%% Load the adjacency matrix
corr_tab = readtable(mat_path, 'ReadVariableNames', false, 'Delimiter', ' ');
corr_mat = table2array(corr_tab);
corr_mat = weight_conversion(corr_mat, 'autofix');

%% Evaluate versatility across resolution parameters (and save the results).
versatility = find_optimal_gamma_curve(corr_mat, @find_louvain_communities_neg_asym, 0.5:0.05:4.0, 250);
versatility_path = strcat(out_prefix, '_Versatility.txt');
fileID = fopen(versatility_path, 'w');
fprintf(fileID, '%.15e \r', versatility);
fclose(fileID);
% Set gamma yourself based on the versatility graph.
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
reorder_order_path = strcat(out_prefix, '_ReorderedMatrix_order.txt');
fileID = fopen(reorder_order_path, 'w');
fprintf(fileID, '%.15e \r', corr_mat_order);
fclose(fileID);
%% Calculate nodal metrics
wmd_z = module_degree_zscore(corr_mat, consensus_partition, 0);
[pc_pos, pc_neg] = participation_coef_sign(corr_mat, consensus_partition);
%% Save Within Module Degree Z-scores
wmdz_path = strcat(out_prefix, '_WMDz.txt');
fileID = fopen(wmdz_path, 'w');
fprintf(fileID, '%.15e \r', wmd_z);
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