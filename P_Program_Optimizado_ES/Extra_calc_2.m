function [Extra2]=Extra_calc_2(il,l1,l2,m1,m2)
CL=Full_CL;

if l1==0
    CL1=1;
else 
    CL1=CL(1:2*l1+1,1:2*l1+1,l1);
end

if l2==0
    CL2=1;
else 
    CL2=CL(1:2*l2+1,1:2*l2+1,l2);
end    

l3=1;
m3=il;

CL3=CL(1:2*l3+1,1:2*l3+1,1);


m1c=m1+l1+1;
m2c=m2+l2+1;
m3c=m3+l3+1;
%m4c=m4+l4+1;      




P1S=0;
for mu1=-l1:l1
    for mu2=-l2:l2
        for mu3=-l3:l3
            
                if mu1+mu2+mu3==0 
              
                 % if (abs(l1-l2)<=l3) && l3<=l1+l2 
                P1S=P1S+sqrt((2*l1+1)*(2*l2+1)*(2*l3+1)/(4*pi))*CL1(m1c,mu1+l1+1)*CL2(m2c,mu2+l2+1)*CL3(m3c,mu3+l3+1)*Wigner3j([l1,l2,l3],[0,0,0])*Wigner3j([l1,l2,l3],[mu1,mu2,mu3]);%;
                end
              
                 
                    
         end
            
     end
end
Extra2=P1S;%P1
end