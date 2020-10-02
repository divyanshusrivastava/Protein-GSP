%% Plotting transmembrane alpha helix and beta sheet proteins spectrum

%%
clear;
close all;
clc;

%% Adding paths
data_dir = 'data/';
gsp_path = 'dependencies/gspbox/gsp_start.m';
run(gsp_path)
clc;
 
%% Parameters
alpha_path = strcat(data_dir, 'transmembrane/');
beta_path = strcat(data_dir, 'globular/');
alphafilestruct = dir(alpha_path);
alphagft = cell(1, length(alphafilestruct)-2);
betafilestruct = dir(beta_path);
betagft = cell(1, length(betafilestruct)-2);
threshold = 7;

%% Generating FFT
clc;
for i = 3: length(alphafilestruct)  % ignoring . and .. directories
    disp(strcat('Processing--', alphafilestruct(i).name));
    filepath = strcat(alpha_path, alphafilestruct(i).name);
    G = get_rig_from_pdb(filepath, threshold, true);
    alphagft{i-2} = gsp_gft(G, G.signal);
    if i == 3
        membrane_graphs = G;
    else
        membrane_graphs(i-2) = G;
    end
end
for i = 3: length(betafilestruct)  % ignoring . and .. directories
    disp(strcat('Processing--', betafilestruct(i).name));
    filepath = strcat(beta_path, betafilestruct(i).name);
    G = get_rig_from_pdb(filepath, threshold, true);
    betagft{i-2} = gsp_gft(G, G.signal);
    if i == 3
       globular_graphs = G;
    else
       globular_graphs(i-2) = G;
    end
end

%% Finding Correlations and HEAT MAP
size_alphagtf = zeros(1, length(alphagft));
size_betagtf = zeros(1, length(betagft));
n_alpha = length(size_alphagtf);
n_beta = length(size_betagtf);
for i = 1:n_alpha
    size_alphagtf(i) = size(alphagft{i}, 2);
end
for i = 1:n_beta
    size_betagtf(i) = size(betagft{i}, 2);
end

%% Take only 200 to 500 bp
l1 = [];
l2 = [];
for i = 1:n_alpha
   if size_alphagtf(i) < 500 && size_alphagtf(i) > 200
       l1 = [l1 i];
   end
end
for i = 1:n_beta
   if size_betagtf(i) < 500 && size_betagtf(i) > 200
       l2 = [l2 i];
   end
end

n_alpha = length(l1);
n_beta = length(l2);


%%
n_proteins = n_alpha + n_beta;
n_spectrum = min(min(size_alphagtf(l1)), min(size_betagtf(l2)));
all_gft = zeros(n_proteins, n_spectrum);

for i = 1:n_alpha
    all_gft(i, :) = alphagft{l1(i)}(1:n_spectrum);
end
for i = 1:n_beta
    all_gft(n_alpha + i, :) = betagft{l2(i)}(1:n_spectrum);
end


%% Display
figure;
colormap();
all_gft = log(all_gft);
imagesc(all_gft);
colorbar;

%%
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


