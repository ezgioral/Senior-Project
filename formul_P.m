%%%% elle yazma
clc;
clearvars;
close all;

SNR = [0:1:20];
f1 = figure;
min_node = 2;
max_node = 6;
% q = 0:0.001:1
SNR = 0:20 ; 
snrR = (10.^(SNR/10));

probs = exp(-1./snrR);  %% s for success?
probf = 1- exp(-1./snrR); %% f for fail?  veeee bu formül do?ru mu?
q = probf; %%% q = probf olsa ve 1-P ye göre plot edilse daha iyi gibi sorr!!
number = [min_node:1:max_node]';  
names = int2str(number);                                   

for nodes = min_node:1:max_node    
switch nodes
      case 2
          P = 1-q;
      case 3
          P = 1-3*q.^(2)+2*q.^(3);  
      case 4
          P = 1-4*q.^(3)-3*q.^(4)+12*q.^(5)-6*q.^(6);
      case 5
          P = 1-5*q.^(4)-10*q.^(6)+20*q.^(7)+30*q.^(8)-60*q.^(9)+24*q.^(10);
      case 6
          P = 1-6*q.^(5)-15*q.^(8)+20*q.^(9)+120*q.^(11)-90*q.^(12)-270*q.^(13)+360*q.^(14)-120*q.^(15);  
end 
figure(f1);
hold on
semilogy(SNR,1-P,'-');
hold off
legend(names);
title('Probability of NOT Reching Consensus vs. Different SNR');
end
