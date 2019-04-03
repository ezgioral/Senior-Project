clc;
close all;
clear all;
f1=figure;
montemax = 1e3;
P_reach = zeros(montemax,1);
min_node = 2;
max_node = 6;
number = (min_node:max_node)';  
names = int2str(number); 

%%SNR
SNR = 0:20;
snrR = (10.^(SNR/10));

%% p-q 
probs = exp(-1./snrR);  
probf = 1- exp(-1./snrR); 
q = probf;
p = 1-q; 

%% simulations
for nodes = min_node:1:max_node 
    no_of_edge = nchoosek(nodes,2);  
    
for SNR_counter = 1:length(SNR)     
tic    
for montecarlo = 1:montemax 
   
P = rand(no_of_edge,1);
%%
Ptoplama_giden=zeros(no_of_edge,1);
for k = 1:no_of_edge
   
if P(k) <= p(SNR_counter)
   Ptoplama_giden(k,1) = 1;
else 
   Ptoplama_giden(k,1) = 0; 
end
end

%% Prob container
[ii,jj] = ndgrid(1:nodes);             
A = zeros(nodes);                       
A(jj>ii) =  Ptoplama_giden; 
A = A + A';
D = diag(sum(A));  
L = D - A ;
rank(L);
if rank(L) == nodes - 1 ;
P_reach(montecarlo,1) =1;   %% ****
else
P_reach(montecarlo,1) =0;    
end

end
Pavg = sum(P_reach)/montemax;
Probability_Container(1,SNR_counter) =  Pavg;
P_reach = zeros(montemax,1);
toc
end
figure(f1);
hold on
semilogy(SNR,1-Probability_Container,'x -');
hold off
legend(names);
title('Probability of NOT Reching Consensus vs. Different SNR (random)');
savefig('SNR_RANDOM.fig')
end
