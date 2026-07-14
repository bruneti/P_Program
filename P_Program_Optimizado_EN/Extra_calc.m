% Delta_R=(0,0,0); ILL' implementation


function [Extra]=Extra_calc(v,k,SLZ,l1,l2,s1,s2,z1,z2)


slz1=[s1,l1,z1];
slz2=[s2,l2,z2];



VL1=find(ismember(SLZ',slz1,'rows'));

VL2=find(ismember(SLZ',slz2,'rows'));
if size(VL1,1)==0 || size(VL1,2)==0 || size(VL2,1)==0 || size(VL2,2)==0  %Identify the correct V_l to use.

   Extra_1=0;

else


          Extra_1=trapz(k,k.^3.*conj(v(VL1,:)).*v(VL2,:));
        
end



Extra=Extra_1;%P2

end