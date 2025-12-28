function y = proend(inpseq,cutp,lr,delp)
% this function processes the cut of a DNA, i.e., indel
% delp = delete probability
% inpseq - input sequence
% cutp - cut position/locus
% lr - process toward left (0) or right (1)

% Insert probability is half of its deletion counterparts
% insp - insert probability
insp = delp/2;
delseq = 0; % indicator of deletion length
insseq = 0; % indicator of insertion length
% showedit = zeros(length(inpseq),2);
% showedit(:,1)=inpseq;

     mut = rand; 
     if (mut>delp(1) && mut<=delp(2))
            delseq =1;% delseq
        elseif (mut>delp(2)&&mut<=delp(3))
            delseq =2;
        elseif (mut>delp(3)&&mut<=delp(4))
            delseq =3; 
        elseif (mut>delp(4)&&mut<=delp(5))
            delseq =4;
        elseif (mut>delp(5))
            delseq =5;
     else
     end
    
      mut = rand; 
     if (mut>insp(1)&&mut<=insp(2))
            insseq =1;% delseq
        elseif (mut>insp(2)&&mut<=insp(3))
            insseq =2;
        elseif (mut>insp(3)&&mut<=insp(4))
            insseq =3; 
        elseif (mut>insp(4)&&mut<=insp(5))
            insseq =4;
        elseif (mut>insp(5) && mut<=0.5)
            insseq =5;
     else
     end

if (delseq ==0) && (insseq == 0) % no more editing on either end
    y = inpseq;
    return;
end

% if delseq > 0
%  if (lr) % toward right but delete bp on left end
%       delpiece = inpseq(cutp-delseq:cutp-1); %delete piece
%  else % toward left but delete bp on right end
%      delpiece = inpseq(cutp:cutp+delseq-1); %delete piece
%  end
% end

if delseq > 0
    if (lr)
        inpseq(cutp-delseq:cutp-1) = 0;
        cutp = cutp-delseq;
    else
        inpseq(cutp:cutp+delseq-1) = 0;
        cutp = cutp+delseq;
    end

end

if insseq >0
    inspiece = ceil(4*rand(insseq,1));
    if (lr)            
        inpseq(cutp:cutp+insseq-1) = inspiece;
    else
        inpseq(cutp-insseq:cutp-1) = inspiece;
    end

end

y = inpseq;
%showedit(:,2)=inpseq;
%showedit(cutp-10:cutp+10,:)