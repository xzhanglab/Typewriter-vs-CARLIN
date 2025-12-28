function y = twscore(u,v,twtarn,tanlen)
% This program calculates similarity score between two dna tape alleles
y = 0;
% matched entry award
matawd = 1;
for i = 1:twtarn
    j=1;
    uentry = u((i-1)*tanlen+j);
    ventry = v((i-1)*tanlen+j);
    while uentry==ventry && uentry~=0
        y = y+matawd;
        j = j+1;
        if j<=tanlen
            uentry = u((i-1)*tanlen+j);
            ventry = v((i-1)*tanlen+j);
        else
            break;
        end
    end
end