##############################################################
### GNU-PLOTTING THE GRAPH ANALYSIS DATA
### PLOT01: Degree, CC, and Avg Shortest Path of each node
### PLOT02: Degree Distribution (/Normalized Degree Distn.) 
### PLOT03: Degree Correlations: All nodes and Avg. for a degree 
### PLOT04: CC x k  and CC(k) x k  
### PLOT05: Clustering Coefficient Distribution (with `cc_dist_bin' bins) 
### PLOT06: Shortest Path Distribution (Non-normalized or Normalized) 
##############################################################

## DATA INPUT ##
n=104; n_breath=5 #Breating space for the n-axis (aa#)
max_degree=16; k_breath=5 #Breating space for the k-axis (degree)
diameter=4; dia_breath=5 #Breating space for the diameter-axis 

set terminal postscript eps enhanced color
##############################################################
### PLOT01: Degree, CC, and Avg Shortest Path of each node ###
set output "node_k_cc_avgsp.ps"

set pointsize 0.3 ## Scalable Point Size - For Clarity of Plots
set nologscale 

## ENTER: MULTI-PLOT MODE ##
set size 1,1
set origin 0,0

set multiplot
set size 1,0.3333
set origin 0,0.6666
set xlabel "Amino Acid Number (i)" font "Helvetica,15"
set ylabel "Degree (k_{i})" 1,0 font "Helvetica,15"
set xtics 50; set mxtics 5
set ytics 5; set mytics 5
plot [:n+n_breath][] 'degree.txt' u 1:2 notitle w linesp 

set size 1,0.3333
set origin 0,0.3333
set xlabel "Amino Acid Number (i)" font "Helvetica,15"
set ylabel "Clustering Coeff. (c_{i})" 1,0 font "Helvetica,15"
set xtics 50; set mxtics 5
set ytics 0.2; set mytics 2
plot [:n+n_breath][] 'cc.txt' u 1:2 notitle w linesp 

set size 1,0.3333
set origin 0,0
set xlabel "Amino Acid Number (i)" font "Helvetica,15"
set ylabel "Avg. Shortest Path (L_{i})" 1,0 font "Helvetica,15"
set xtics 50; set mxtics 5
set ytics 5; set mytics 5
plot [:n+n_breath][] 'avg_shortest_paths.txt' u 1:2 notitle w linesp 

replot
set nomultiplot
## EXIT: MULTI-PLOT MODE ##

set nogrid
##############################################################

##############################################################
### PLOT02: Degree Distribution (/Normalized Degree Distn.) ###

set output "degree_dist.ps"

set pointsize 0.8 ## Scalable Point Size - For Clarity of Plots
set nologscale 
set size 1,1
set origin 0,0

set xlabel "Degree (k)" font "Helvetica,25"
set ylabel "N(k)" 1,0 font "Helvetica,25"
set xtics 1; set mxtics 1
set ytics 0.1; set mytics 0.1
plot [:max_degree+k_breath][:1.1] 'degree_dist_norm.txt' u 1:2 title 'Degree Distribution' w linesp 7
replot

set nogrid
##############################################################

##############################################################
### PLOT03: Degree Correlations: All nodes and Avg. for a degree ###

set output "degree_corr.ps"

set pointsize 0.8 ## Scalable Point Size - For Clarity of Plots
#set logscale xy 10
#set format "10^{%T}"
set size 1,1
set origin 0,0

set xlabel "Degree (k)" font "Helvetica,25"
set ylabel "<K_{nn}>" 1,0 font "Helvetica,25"
set xtics 5; set mxtics 5
set ytics 5; set mytics 5
plot [:max_degree+k_breath][:] 'degree_corr.txt' u 1:2 notitle w points 6,\
	 'degree_corr_avg.txt' u 1:2 notitle w points 7
replot
	
set nogrid
##############################################################

##############################################################
### PLOT04: CC x k  and CC(k) x k  ###

set output "ck.ps"

set pointsize 0.8 ## Scalable Point Size - For Clarity of Plots
set nologscale 
set size 1,1
set origin 0,0

set xlabel "k (Degree)" font "Helvetica,25"
set ylabel "c(k)" 1,0 font "Helvetica,25"
set xtics 5; set mxtics 5
set ytics 0.2; set mytics 2
plot [:max_degree+k_breath][0:1] 'cc_degree.txt' u 2:1 notitle w points 6,\
	 'cc_degree_corr.txt' u 1:2 notitle w points 7
replot

set nogrid
##############################################################

##############################################################
### PLOT05: Clustering Coefficient Distribution (with `cc_dist_bin' bins) ###
###		Non-normalized or Normalized			###
set output "cc_dist.ps"

set pointsize 0.8 ## Scalable Point Size - For Clarity of Plots
set nologscale 
set size 1,1
set origin 0,0

set xlabel "Clustering Coefficient (cc)" font "Helvetica,25"
set ylabel "N (cc)" 1,0 font "Helvetica,25"
set xtics 1; set mxtics 0
set ytics 10; set mytics 2
plot [:][:] 'cc_dist.txt' u 1:2 notitle w linesp 7

replot

set nogrid
##############################################################

##############################################################
### PLOT06: Shortest Path Distribution (Non-normalized or Normalized) ###

set output "shortest_path_dist.ps"

set pointsize 0.8 ## Scalable Point Size - For Clarity of Plots
set nologscale 
set size 1,1
set origin 0,0

set xlabel "Shortest Path (L)" font "Helvetica,25"
set ylabel "N (L)" 1,0 font "Helvetica,25"
set xtics 5; set mxtics 5
set ytics 500; set mytics 5
plot [:diameter+dia_breath][:] 'shortest_path_dist.txt' u 1:2 notitle w linesp 7

replot

set nogrid
##############################################################
set output
set terminal X11
