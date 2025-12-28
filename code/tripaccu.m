function y = tripaccu(intree,trbk,it,ss)
% this function calculates tree accuracy using triple nodes
A=[];% record all triple-node combinations, remove duplication later
for i=it-1:-1:it-trbk
    L1 = size(intree{i},1);
    L2 = size(intree{i},2);
    for j=1:L1
        for k=1:L1
            if (j~=k)
                for j1=1:L2
                    if intree{i}(j,j1)>0
                        for k1=1:L2-1
                            if intree{i}(k,k1)>0
                                for k2=k1+1:L2
                                    if intree{i}(k,k2)>0 % find one triplet
                                        temptrip = sort([intree{i}(j,j1) intree{i}(k,k1) intree{i}(k,k2)]);
                                        A=[A; temptrip intree{i}(j,j1)];
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
B = unique(A,'rows')+2^it-1;
allnodes = ceil(ss*2^it);
allcomb = allnodes*(allnodes-1)*(allnodes-2)/6;  % n choose 3

corrtrip = 0; % correct triplets count
for i = 1:length(B)
        tempcomp1=dec2bin(B(i,1));
        tempcomp2=dec2bin(B(i,2));
        tempcomp3=dec2bin(B(i,3));
        for j=2:it+1
            if(tempcomp1(j)==tempcomp2(j) && tempcomp2(j)~=tempcomp3(j))
                if (B(i,3)==B(i,4))
                    corrtrip = corrtrip+1;
                end
                break;
            elseif (tempcomp2(j)==tempcomp3(j) && tempcomp1(j)~=tempcomp2(j))
                if (B(i,1)==B(i,4))
                    corrtrip = corrtrip+1;
                end
                break;
            end
        end
end

y=100*corrtrip/allcomb;