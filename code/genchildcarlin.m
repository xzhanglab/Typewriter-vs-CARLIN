% This program generates a child barcode from a parent barcode
function [y, nstp, nendp,nlivetar,nhotspots] = genchildcarlin(n,ldtl,stp,endp,parbarc,tarn,livetar, hotspot, mupb,lgdelprob,ins_sub)
% k is division round
% j is cell index in this generation
% stp =sted(1,2^(k-1)+j-1) (for first child) is the start position in the barcode
% endp = sted(2,2^(k-1)+j-1) (for first child)is the end position
% parbarc = barM(:,2^(k-1)+j-1) is the parent barcode
% livetar records any unedited targets
% hotspot records cutting hotspot of each target
% tarn = target number

% lgdelprob is the base inter target deletion probability, given more than
% 1 cuts. The probability will be higher is more cuts occur.
cutsites = []; % record all cut sites
editedtar = []; % indicate which target is cut
for i=1:tarn
    if livetar(i) % editable target
         mut = rand; % indicator of mutation/cut
         if (mut<=mupb) % mutation/cut occurs
             editedtar = [editedtar i];
             % determine the cut locus
             mut = rand;
             if mut<=0.6
                 cutloc = hotspot(i);
             elseif mut<=0.85
                 cutloc = hotspot(i)-1;
             elseif mut<=0.95
                 cutloc = hotspot(i)-2;
             else
                 cutloc = hotspot(i)-3;
             % elseif mut<=0.95
             %     cutloc = hotspot(i)-4;
             % else
             %     cutloc = hotspot(i)-5;
             end
             if cutloc>0
                cutsites = [cutsites cutloc]; %dynamically extend cutsites
             end
         end
    end
end

y=parbarc; % full parent barcode


                    if (length(editedtar)>=2)% more than 2 cuts
                        lgdelprob = lgdelprob+(length(editedtar)-2)*0.03;
                        if (rand<lgdelprob) % perform one large deletion
                            randi = randsample(length(editedtar),2); % randomly select two sites to have a large deletion
                            y(cutsites((min(randi))):cutsites(max(randi))-1) = 0; % segment deletion
                            y = proend(y,cutsites((min(randi))),1,ins_sub); % left end, process toward right=1
                            y = proend(y,cutsites(max(randi)),0,ins_sub); % right end, process toward left=0
                            cutsites((min(randi)):(max(randi))) = 0; % clear all the involved cutsites
                            livetar(editedtar(min(randi)):editedtar(max(randi))) = 0; % indicate that those targets are no longer editable
                            editedtar(min(randi):max(randi)) = 0; % these edited targets are removed
                        end
                    end
                    
                    newcut = cutsites(cutsites~=0);% remaining cut sites for targets that have been cut
                    newedtar = editedtar(editedtar~=0);
% process all remaining cut sites=============================
% delp = delete probability
% inpseq - input sequence
% cutp - cut position/locus
% indlive - indicate if this target is still editable(1) or not (0)
% Insert probability is half of its deletion counterparts
% insp - insert probability
% nstp - new start point
% nendp - new end point
insp = ins_sub/2;
for i=1:length(newcut) 
    cutp = newcut(i);
totaldel = 0; % total deletion length
totalins = 0; % total insertion length
delseq1 = 0; % indicator of deletion length
ldel = [];
rdel = [];
lins = [];
rins = [];
     mut = rand; 
     if (mut>ins_sub(1) && mut<=ins_sub(2))
            delseq1 =1;% delseq
        elseif (mut>ins_sub(2)&&mut<=ins_sub(3))
            delseq1 =2;
        elseif (mut>ins_sub(3)&&mut<=ins_sub(4))
            delseq1 =3; 
        elseif (mut>ins_sub(4)&&mut<=ins_sub(5))
            delseq1 =4;
        elseif (mut>ins_sub(5))
            delseq1 =5;
     else
     end
    
if delseq1 > 0
      ldel = y(cutp-delseq1:cutp-1); % left delete piece
      totaldel = totaldel+delseq1;
end

delseq2 = 0; % indicator of deletion length
     mut = rand; 
     if (mut>ins_sub(1) && mut<=ins_sub(2))
            delseq2 =1;% delseq
        elseif (mut>ins_sub(2)&&mut<=ins_sub(3))
            delseq2 =2;
        elseif (mut>ins_sub(3)&&mut<=ins_sub(4))
            delseq2 =3; 
        elseif (mut>ins_sub(4)&&mut<=ins_sub(5))
            delseq2 =4;
        elseif (mut>ins_sub(5))
            delseq2 =5;
     else
     end

if delseq2 > 0
      rdel = y(cutp:cutp+delseq2-1); % right delete piece
      totaldel = totaldel+delseq2;
end
%  =============insertion===================
     insseq1 = 0; % indicator of insertion length
      mut = rand; 
     if (mut>insp(1)&&mut<=insp(2))
            insseq1 =1;% delseq
        elseif (mut>insp(2)&&mut<=insp(3))
            insseq1 =2;
        elseif (mut>insp(3)&&mut<=insp(4))
            insseq1 =3; 
        elseif (mut>insp(4)&&mut<=insp(5))
            insseq1 =4;
        elseif (mut>insp(5) && mut<=0.5)
            insseq1 =5;
     else
     end
 if insseq1 >0
    lins = ceil(4*rand(insseq1,1)); % left insertion piece
    totalins = totalins+insseq1;
 end

      insseq2 = 0; % indicator of insertion length
      mut = rand; 
     if (mut>insp(1)&&mut<=insp(2))
            insseq2 =1;% delseq
        elseif (mut>insp(2)&&mut<=insp(3))
            insseq2 =2;
        elseif (mut>insp(3)&&mut<=insp(4))
            insseq2 =3; 
        elseif (mut>insp(4)&&mut<=insp(5))
            insseq2 =4;
        elseif (mut>insp(5) && mut<=0.5)
            insseq2 =5;
     else
     end
 if insseq2 >0
    rins = ceil(4*rand(insseq2,1)); % right insertion piece
    totalins = totalins+insseq2;
 end

    if (totaldel>0 || totalins > 0) % deletion or insertion occurs
    % showedit = zeros(length(y),2);
    % showedit(:,1)=y;
      
     if totaldel>= totalins % no need to shift entries in y
         % delete entries
         if delseq1>0
            y(cutp-delseq1:cutp-1) = 0;
         end
         if delseq2>0
            y(cutp:cutp+delseq2-1) = 0;
         end
        % insert entries
        if insseq1>0
            y(cutp-delseq1:cutp-delseq1+insseq1-1) = lins;
        end
        if insseq2>0
            y(cutp+delseq2-insseq2:cutp+delseq2-1) = rins;
        end
        % check if mutated
        if totaldel== totalins
            if [ldel; rdel] ~= [lins; rins]% mutated target, set uneditable any more
                livetar(newedtar(i)) = 0;
            end
        else
            livetar(newedtar(i)) = 0;
        end
     else % need to shift entries in y
         livetar(newedtar(i)) = 0; % for sure this target has been mutated
        if (stp<n+2*ldtl-endp) % insert downstream
            y(cutp+totalins-totaldel:endp+totalins-totaldel)=y(cutp:endp);
            endp=endp+totalins-totaldel;
            y(cutp-delseq1:cutp+totalins-delseq1-1)=[lins; rins]; % insert this piece
            newcut(i:end)=newcut(i:end)+totalins-totaldel; % update rest cut sites, only for downstream insertion because upstream cut sites have been processed
            hotspot(newedtar(i):end)=hotspot(newedtar(i):end)+totalins-totaldel; % update hotspots for downstream targets
        else % insert upstream
            y(stp-totalins+totaldel:cutp-totalins+totaldel)=y(stp:cutp);
            if (newedtar(i)>1) % shift upstream sites
                hotspot(1:newedtar(i)-1)=hotspot(1:newedtar(i)-1)-totalins+totaldel;
            end
            stp=stp-totalins+totaldel;
            % insert the segment uptream
            y(cutp-totalins+delseq2:cutp+delseq2-1)=[lins; rins];
        end
    % showedit(:,2)=y;
    % ldel
    % rdel
    % lins
    % rins
    % showedit(cutp-10:cutp+10,:)
     end
    end
end

nstp = stp;
nendp = endp;
nlivetar = livetar;
nhotspots = hotspot;