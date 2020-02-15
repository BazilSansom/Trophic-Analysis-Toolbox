% EXAMPLES

% This script will add examples/tests for different functions as the are added to the toolbox.

% (1) Show some small examples to illustrate levels and incoherence

% Define some small networks with different structures, then
% obatains trophic levels and incoherence stats (F0, eta) and plots results.

% (1.1) Define some simple networks to use as examples/illustrations

NETS=cell(9,1);

%dyad=[0,1;0,0];                          NETS{1,1}=dyad;
strong_edges=[1,2;2,1;1,3;3,1;2,3;3,2]; 
strong=edgelist2adj(strong_edges);       NETS{1,1}=strong;
star=[0,1,1;0,0,0;0,0,0];                NETS{2,1}=star;
chain=[0,1,0;0,0,1;0,0,0];               NETS{3,1}=chain;
ffl=[0,1,1;0,0,1;0,0,0];                 NETS{4,1}=ffl;
fbl=[0,1,0; 0,0,1;1,0,0];                NETS{5,1}=fbl;
fbl2=[0,1,0,0; 0,0,1,0;1,0,0,0;1,0,0,0]; NETS{6,1}=fbl2;
coherent_edges=[1,4; 1,5; 2,5; 3,5; 3,6; 4,7; 5,7; 5,8; 5,9; 6,9];
coherent=edgelist2adj(coherent_edges);   NETS{7,1}=coherent;
semi_coherent_edges=[1,4; 1,5; 2,5; 3,5; 3,6; 4,7; 5,7; 5,8; 5,9; 6,9; 1,7; 2,8; 3,9; 1,8;2,7;2,9];
semi_coherent=edgelist2adj(semi_coherent_edges); NETS{8,1}=semi_coherent;
g1=[1,2;2,1;1,3;3,1;2,3;3,2]; g2=g1+max(max(g1)); link=[3,5]; gfull=[g1;g2;link];
strong_comp=edgelist2adj(gfull);         NETS{9,1}=strong_comp;


% (1.2) Analyse and plot

figure(1)

for i=1:size(NETS,1)
    
    A=NETS{i,1};
    [h,Fmin,eta] = incoherence(A);
    %figure(100)
    %p=plot(digraph(A),'layout','layered');
    %xdata=get(p,'XData');
    %close
    figure(1)
    subplot(3,3,i)
    %plot(digraph(A),'XData',xdata,'YData',h)
    G=digraph(A);
    plot(G,'NodeLabel',round(h,1))
    set(gca,'xtick',[])
    title({['F0=',num2str(round(Fmin,2))],['eta=',num2str(round(eta,1))]})

end

