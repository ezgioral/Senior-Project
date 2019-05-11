%%%% Rician Switching

clc;
clear all;
close all;
p = 1; 
MaxIt = 1000;
SNR = [0:1:20];
epsilon = 1/50;
f1 = figure;
f2 = figure;

min_node = 3;
max_node = 7;

montemax=1e6; 

counter = 0;  %%% iterasyon sayisini hesaplamak icin
countermtx = zeros(montemax,1);    %%% countermtx = zeros(montemax);

Avg = zeros(size(SNR));
number = [min_node:1:max_node]'; 
names = int2str(number);                              
Allowed_Error = 1/100;

Ranks = zeros(montemax,1);   %% L nin ranklerini kaydedecek

for nodes = min_node:1:max_node
N = nodes^2 ;  
Edge_num = nchoosek(nodes,2);

switch nodes
      case 3
          Measured = [-6 ; -3 ; 7];   
      case 4
          Measured = [-6 ; -3 ; 7 ; 14 ];
      case 5
          Measured = [-6 ; -3 ; 7 ; 14; 21];
      case 6
          Measured = [-6; -3; 7; 14; 21; -12.5];
      case 7
          Measured = [-6; -3; 7; 14; 17; -12.5; 8.5];
      case 8
          Measured = [-6; -3; 7; 14; 17; -12.5; 8.5; -9.5];
    end   

StateContainer = zeros(size(Measured,1), MaxIt + 1); 

    for SNR_counter = 1:length(SNR)
        tic
        for monte = 1:montemax
            %%monte kontrol icin
            IterMeasured = Measured;
            StateContainer(:,1) = IterMeasured(:);
            
            for k = 1:MaxIt               
            
      
                 R = 1 ;                          
                 K=3;					
                 mu = sqrt( K/(2*(K+1)) );
                 s = sqrt( 1/(2*(K+1)) );
      
                gama = 10^(SNR(SNR_counter)/10);  %gama değişmiyor (shannon, kapasite)
                threshold = (2^R - 1) / gama ;           % R=1 ?
                h_Rician=abs( s*randn(1,Edge_num) + mu ) + 1i*( s*randn(1,Edge_num) + mu ); 
                Edge_con = abs(h_Rician).^2 > threshold ;   
      
                [ii,jj] = ndgrid(1:nodes);              
                A = zeros(nodes);   
                
                A(jj>ii) =  Edge_con;                   
                A = A + A';                            
                D = diag(sum(A));   %% sayi dizisinin diagonal olarak bastiriyor                
                L = D - A  ;  
      
                P_epsilon = eye(nodes) - epsilon * L;
                IterMeasured = P_epsilon * IterMeasured ;
                StateContainer(1:end,k+1) = IterMeasured(1:end); %% k+1 olmas?n?n sebebi ilk de?eri de görmek istemem         
            %%% Iterasyon sayisi hesaplarken -1 yapmak gerekebilir. Buraya dikkat et
            %%% StateContainer 157 ye kadar doluysa 156 iterasyon yapilmis demektir.
            
                if rank(L) >= nodes-1
                Ranks(monte,1)=1;
                else
                Ranks(monte,1)=0;
                end
                
                if abs(max(IterMeasured) - min(IterMeasured)) <= Allowed_Error 
                    counter = counter + 1; %% en son yaptigi iterasyonu da dahil ediyor
                    counter;   %%% Statecontainer ile dogrulugunu kontrol etmek icin
                    break
                else
                    counter = counter + 1 ;
                end
                
            end
            countermtx(monte) = counter;
            SNR_counter;
            counter = 0;
        end
        %countermtx
        %%% ilerde belki daha az degisken kullanilarak sum elde edilebilir        
        sum = sum(countermtx); %%% matrix halinde veriyor montemax*1 boyutunda
        %%summ = sum(1);  %% mtrxin ilk elemanini seciyor
        Avg(SNR_counter) = sum/montemax; 
        clear sum
        countermtx = zeros(montemax,1);  %%% önemli
        
        Succ_Probability = sum(Ranks(:,1)) / montemax;
        Probability_Container(1,SNR_counter) =  Succ_Probability;
        Ranks = zeros(montemax,1);
        toc
    end    
figure(f1);
hold on;
semilogy(SNR,1-Probability_Container,'x -');
legend(names);
title('Propability of Reaching Consensus vs. Different SNR');

figure(f2);
hold on;
semilogy(SNR,Avg,'x -');
legend(names);
title('Average Number of Iterations vs. Different SNR');

end