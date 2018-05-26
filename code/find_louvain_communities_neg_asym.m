function [M,Q]=find_louvain_communities_neg_asym(W, gamma)

[M, Q] = community_louvain(W, gamma, [], 'negative_asym');