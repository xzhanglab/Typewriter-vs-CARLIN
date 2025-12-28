function y = twparent(u,v,twtarn,tanlen)
% This program rebuilds parent dna tape from two children dna tape alleles
y = zeros(twtarn*tanlen,1);
for i = 1:twtarn
    j=1;
    uentry = u((i-1)*tanlen+j);
    ventry = v((i-1)*tanlen+j);
    while uentry==ventry && uentry~=0
        y((i-1)*tanlen+j)=uentry;
        j = j+1;
        if j<=tanlen
            uentry = u((i-1)*tanlen+j);
            ventry = v((i-1)*tanlen+j);
        else
            break;
        end
    end
end