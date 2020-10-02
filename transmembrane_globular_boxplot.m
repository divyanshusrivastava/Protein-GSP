clear;
close all;
clc;

%%
gsp_path = '/home/divyanshu/MTech/MTech-Thesis/dependencies/gspbox/gsp_start.m';
run(gsp_path);
load('/home/divyanshu/MTech/MTech-Thesis/data/NB-Project/membrane_globular.mat');
%globular_graphs = globular;
%membrane_graphs = membrane;
G_lin = [globular_graphs membrane_graphs];

%%
low_freq_comps = zeros(1, length(G_lin));
clc;
%disp('Name - Assortativity - Clust_Coeff - Low_Frew_comp');
for i = 1:length(G_lin)
    disp(i);
    G = G_lin(i);
    G.gft = gsp_gft(G, G.signal);
    temp = abs(G.gft);
    low_freq_component = sum(temp(G.e < (G.e(end)/2))) / sum(temp);
    low_freq_comps(i) = low_freq_component;    
end

%% Mean 
m1 = mean(low_freq_comps(1:length(globular_graphs)));
m2 = mean(low_freq_comps(length(globular_graphs)+1:length(G_lin)));

%% figure
figure; hold on;
title('hydrophobicity signal on transmembrane and globular proteins')
ylabel('Low Frequency Component');
xlabel('Protein #')
scatter(1:length(globular_graphs), ...
    low_freq_comps(1:length(globular_graphs)), '.r');
scatter(length(globular_graphs)+1:length(G_lin), ...
    low_freq_comps(length(globular_graphs)+1:length(G_lin)), 'xb');
hold off;

%% figure
figure; hold on;
title('Low Frequency Component Box Plot');
ylabel('Low Frequency Component')
groups = cell(1, length(G_lin));
groups(1:length(globular_graphs)) = {'Globular Proteins'};
groups(length(globular_graphs)+1:length(G_lin)) = {'Transmembrane Proteins'};
boxplot(low_freq_comps, groups);
plot([mean(low_freq_comps(1:length(globular_graphs))) ...
    mean(low_freq_comps(length(globular_graphs)+1:length(G_lin)))], ...
    'dg');
hold off;

