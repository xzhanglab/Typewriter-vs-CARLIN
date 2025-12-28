% This program generates type writer alleles
function [y,nactsite] = twalleles(partwallele, twtarn,paractsite,tanlen,insBC,editP,maxsitelen,tarvital)

% insert BC potential? - that is, if one certain insertBC has been inserted
% many times, its chance to be inserted in subsequent mutations will be
% smaller and smaller
% editP - edit/mutation/cut probability of the first active site in each
% target of DNA tape

% partwallele - parent type writer allele
% paractsite - parent active sites

% impf - impact factor of further edits in one target due to different
% site positions, suggest impf>=1 because if the preceding site is edited,
% there is a better chance to have more edits because the target is
% relatively more active than other targets
impf = 1;

% cell vitality - indicator of cell activity level
cellvital = 2.5*rand;

% with a small probability, there is a second round to insert BC at
% active sites
secP = 0.1; % second round probability

for i=1:twtarn
    if paractsite(i)>=1 && paractsite(i)<=maxsitelen(i) % not a corrupted target, so it is still editable
        % if rand<editP*tarvital(i) % insert a BC at this site
        %     partwallele((i-1)*tanlen+paractsite(i)) = ceil(insBC*rand);
        %     paractsite(i) = paractsite(i)+1;
        % end


        if paractsite(i)==1 % use threshold editP
            if rand<editP*tarvital(i)*cellvital
             partwallele((i-1)*tanlen+paractsite(i)) = ceil(insBC*rand);
             paractsite(i) = paractsite(i)+1;
            end
        elseif paractsite(i)>1 % use threshold editP*discf
             if rand<editP*impf*tarvital(i)*cellvital
             partwallele((i-1)*tanlen+paractsite(i)) = ceil(insBC*rand);
             paractsite(i) = paractsite(i)+1;
            end
        end
    end
end


for i=1:twtarn
    if paractsite(i)>=1 && paractsite(i)<=maxsitelen(i) % not a corrupted target, so it is still editable
        % if rand<editP*tarvital(i)*secP % insert a BC at this site
        %     partwallele((i-1)*tanlen+paractsite(i)) = ceil(insBC*rand);
        %     paractsite(i) = paractsite(i)+1;
        % end

        if paractsite(i)==1 % use threshold editP*secP
            if rand<editP*secP*tarvital(i)*cellvital
             partwallele((i-1)*tanlen+paractsite(i)) = ceil(insBC*rand);
             paractsite(i) = paractsite(i)+1;
            end
        elseif paractsite(i)>1 % use threshld editP*discf*secP
             if rand<editP*impf*secP*tarvital(i)*cellvital
             partwallele((i-1)*tanlen+paractsite(i)) = ceil(insBC*rand);
             paractsite(i) = paractsite(i)+1;
            end
        end
    end
end

y = partwallele;
nactsite = paractsite;