% This is Monte Carlo to run funbarnew several times
clear;
clc;
mcit = 10; % Monte Carlo simulation rounds
carlinbc = fileread('CARLIN_raw.txt'); % get CARLIN ref barcode


it = 11; % cell division rounds
ss = 1; % sample size proportion, a number between 0 and 1.

propm = 0.8; % proportion of nonzero counts of a barcode to be considered as a matched pair
mupb = 0.4; % probability of a cut at each target's hotspot, suggested 0.2---0.4
% editP - edit/mutation/cut probability of the first active site in each target
% of DNA tape % suggested value around 0.13
editP = 0.13; 
%editP = 0.9; 
%editP = 0.3; % for drawing illustration
usecomb = 0; % an indicator of whether use combined methods (CARLIN and DNA Tape) to rebuild lineage

ins_sub =  [0.8 0.9 0.95 0.975 0.99]; % probabilities of perfect repair, deleting 1,2,3,4,5 bp respectively, on each end of cut locus
lgdelprob = 0.4; % this is the probability to have a large deletion, given more than 2 cut sites.

pulse = 0; % if pulse = 1, pulse induction; if pulse = 0, constant dox level

divp = 1; % probability of division
        %if division occurs, get two children, then check probability of
        %survival, then mutate.
        %if not dividing, get one child, check probability of surviving,
        %then simulate mutation.
clive = 1; % probability of cell survival/live. NOTE: program may generate error message if all leaf barcodes are empty.

% if DNA type writer is used===============
% twtarn - type writer target number
twtarn = 1;
% tanlen - tandem length, 1-5
tanlen = 2;
% insBC - number of different insertBCs possible
insBC = 10;

% A,C,G,T = 1,2,3,4
tref = zeros(length(carlinbc),1); % convert to 1,2,3,4
for k=1:length(carlinbc)
    switch carlinbc(k)
        case 'A'
            tref(k)=1;
        case 'C'
            tref(k)=2;
        case 'G'
            tref(k)=3;
        case 'T'
            tref(k)=4;
        otherwise
    end
end
carlinref = trim(tref);
n=length(carlinref);% carlin barcode length
tarn = 10; % target number
%trbk = it-1; % this is the number of generations to track back, trbk<it, trbk is at most it-1.
trbk = 4;


 %%for detailed report of CARLIN alleles such as remaining CARLIN
 %%potential and target mutation statistics

 % sumallele = zeros(3,tarn); 
 % 
 %        for i = 1:mcit
 %            i
 %                         [~,~,avelivetar(i,:),tarallele]=funbarNBJNF_det(n,it,propm,ss,mupb,ins_sub,lgdelprob,divp,clive, pulse, trbk, carlinref,tarn); % NBJNF method,detailed report
 %             for k=0:2^it-1
 %                sumallele = sumallele+tarallele(:,:,2^it+k);
 %             end
 %         end
 %    figure
 %    plot(0:it,[10 mean(avelivetar)*10],'black'); % plot average CARLIN potential (# of unedited targets for each generation)
 %    ylim([0,10]);
 %    title('CARLIN Potential');
 %    ylabel('# of unedited targets');
 %    xlabel('Division rounds');
 %    sumallele
 %    %reverse order for bar plot
 %    barsum = sumallele;
 %    barsum(1,:) = sumallele(3,:);
 %    barsum(3,:) = sumallele(1,:);
 %    figure
 %    bar(barsum','stacked');
 %    xlabel('Targets')
 %    ylabel('# of cells');
 %    disp('row 1 - unedited, row 2- intersite deletion with possible indel at both ends, row 3 - intrasite indel');

%===================iterations================================
for twtarn = 2:5:17

    for tanlen=2:5



       %for insBC = 2:4:18
    
        if usecomb == 1
            accuracyNBJNF=zeros(mcit,6);
        else
            accuracyNBJNF=zeros(mcit,4);
        end
        triplet = zeros(mcit,3);
    
        for i = 1:mcit
            i
                       
            if usecomb == 1
                [accuracyNBJNF(i,1), accuracyNBJNF(i,2),  accuracyNBJNF(i,3),  accuracyNBJNF(i,4),accuracyNBJNF(i,5),accuracyNBJNF(i,6),Ctree,TWtree,Botree]=funbarNBJNF_TW(n,it,propm,ss,mupb,ins_sub,lgdelprob,divp,clive, pulse, trbk, carlinref,tarn,twtarn,tanlen,insBC,editP,usecomb); % NBJNF method with ,DNA Type writer
                %triplet(i,3) = tripaccu(Botree,trbk,it);
            else
                [accuracyNBJNF(i,1), accuracyNBJNF(i,2),  accuracyNBJNF(i,3),  accuracyNBJNF(i,4),~,~,Ctree,TWtree,Botree]=funbarNBJNF_TW(n,it,propm,ss,mupb,ins_sub,lgdelprob,divp,clive, pulse, trbk, carlinref,tarn,twtarn,tanlen,insBC,editP,usecomb); % NBJNF method with ,DNA Type writer
            end
            %triplet(i,1) = tripaccu(Ctree,trbk,it);
            %triplet(i,2) = tripaccu(TWtree,trbk,it);
    
         end    
     
        if usecomb ==0
            writematrix([accuracyNBJNF triplet],strcat('MC',num2str(mcit),'div',num2str(it),'trbk',num2str(trbk),'ss',num2str(ss),'mupb',num2str(mupb),'lgdelprob',num2str(lgdelprob), 'twtarn',num2str(twtarn),'tanlen',num2str(tanlen),'insBC',num2str(insBC),'editP',num2str(editP),'divp',num2str(divp),'clive',num2str(clive),'.csv')); %record result
        end
        % 
        [accuracyNBJNF triplet]
        disp('average');
        disp(mean([accuracyNBJNF triplet]))
        disp('max')
        disp(max([accuracyNBJNF triplet]))
        disp('Column 1- CARLIN accuracy of all nodes');
        disp('Column 2- CARLIN accuracy of dividing nodes');
        disp('Column 3- DNA tape accuracy of all nodes');
        disp('Column 4-  DNA tape accuracy of dividing nodes');
    
        if usecomb == 1
            disp('Column 5-  Accuracy of all nodes using both algorithms');
            disp('Column 6-   Accuracy of dividing nodes using both algorithms');
            writematrix([accuracyNBJNF triplet],strcat('MC',num2str(mcit),'div',num2str(it),'trbk',num2str(trbk),'ss',num2str(ss),'mupb',num2str(mupb),'lgdelprob',num2str(lgdelprob), 'twtarn',num2str(twtarn),'tanlen',num2str(tanlen),'insBC',num2str(insBC),'editP',num2str(editP),'divp',num2str(divp),'clive',num2str(clive),'.csv')); %record result
        end
    
        disp('Triplet accuracy (last three columns)-  CARLIN, Type Writer, Dual Algorithm')
    
        %end
    end

end



