function [Z_all, Z_pos, Z_neg]=module_degree_zscore_sign(W,Ci)
%MODULE_DEGREE_ZSCORE_SIGN       Signed within-module degree z-score
%
%   Z=module_degree_zscore_sign(W,Ci);
%
%   The within-module degree z-score is a within-module version of degree
%   centrality. This script is an extension of the traditional signed 
%   within-module degree z-score. Instead of calculating WMDz over all
%   edges to produce a single metric, positive and negative WMDz are 
%   calculated separately
%
%   Inputs:     W,      undirected, signed connection matrix
%               Ci,     community affiliation vector
%               
%
%   Output:     Z_all,      within-module degree z-score.
%               Z_pos,      within-module degree z-score, positive edges.
%               Z_neg,      within-module degree z-score, negative edges.
%
%   Reference: Guimera R, Amaral L. Nature (2005) 433:895-900.
%
%
%   Original script by Mika Rubinov, UNSW, 2008-2010
%   Extension by Dan Lurie, UC Berkeley, 2018


% Traditional WMDz

n=length(W);                        %number of vertices
Z_all=zeros(n,1);
for i=1:max(Ci)
    Koi_all=sum(W(Ci==i,Ci==i),2);
    Z_all(Ci==i)=(Koi_all-mean(Koi_all))./std(Koi_all);
end

Z_all(isnan(Z_all))=0;

% Positive WMDz

W_pos = W;
W_pos(W_pos < 0)=0;                   % zero negative edges

Z_pos=zeros(n,1);
for i=1:max(Ci)
    Koi_pos=sum(W_pos(Ci==i,Ci==i),2);
    Z_pos(Ci==i)=(Koi_pos-mean(Koi_pos))./std(Koi_pos);
end

Z_pos(isnan(Z_pos))=0;

% Negative WMDz

W_neg = W;
W_neg(W_neg > 0)=0;                   % zero positive edges

Z_neg=zeros(n,1);
for i=1:max(Ci)
    Koi_neg=sum(W_neg(Ci==i,Ci==i),2);
    Z_neg(Ci==i)=(Koi_neg-mean(Koi_neg))./std(Koi_neg);
end

Z_neg(isnan(Z_neg))=0;