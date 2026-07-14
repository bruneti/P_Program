Inferno_JiLM=zeros(3,4,4,8,7,7,15); %il,l1,l2,l3,m1,m2,m3
CL=Full_CL; %This is needed for JiLM.
%JiLM(CL,il,l1,l2,l3,m1,m2,m3) as we remember
tic
for il=-1:1
    for l1=0:3
        for l2=0:3
            for l3=0:7
                 for m1=-l1:l1
                    for m2=-l2:l2
                        for m3=-l3:l3
                        
                        Inferno_JiLM(il+2,l1+1,l2+1,l3+1,m1+l1+1,m2+l2+1,m3+l3+1)=JiLM(CL,il,l1,l2,l3,m1,m2,m3);         
                              
                        end
                    end
                end
            end
        end
    end
end

save('Inferno','Inferno_JiLM')
toc

% We compile every possible result of JiLM on a tensor 7D, for faster calculations. 