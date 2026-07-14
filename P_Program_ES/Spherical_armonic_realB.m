% Tabulación del armónico esférico real.

    function [Yrml]=Spherical_armonic_realB(x,y,z,l,m)
%[phi,~,r]=cart2sph(x,y,z);
phi=atan2(y,x);
phi=mod(phi,2*pi);
r=sqrt(x.^2+y.^2+z.^2);
theta=acos(z/r);

if abs(m)>abs(l)
    fprintf('Mal escogido l o m')
end
if l<0
    fprintf('Mal escogido l')
end
if l>7
    fprintf('No tenemos este valor tabulado')
end
if l==0 
   if m==0
       Yrml=1/2*sqrt(1/pi);
   end

end
if l==1
    if m==-1
       Yrml= sqrt(3/(4*pi))*sin(theta)*sin(phi);
    end
    if m==0
        Yrml=sqrt(3/(4*pi))*cos(theta);

    end
    if m==1
        Yrml=sqrt(3/(4*pi))*sin(theta)*cos(phi);
    end
end

if l==2
   if m==-2
       Yrml=1/4*sqrt(15/pi)*sin(theta)^2*sin(2*phi);
   end
   if m==-1
       Yrml=1/4*sqrt(15/pi)*sin(2*theta)*sin(phi);
   end
   if m==0
       Yrml=1/4*sqrt(5/pi)*(3*cos(theta)^2-1);
   end
   if m==1
Yrml=1/4*sqrt(15/pi)*sin(2*theta)*cos(phi);
   end
if m==2
    Yrml=1/4*sqrt(15/pi)*sin(theta)^2*cos(2*phi);
end
end
if l==3
if m==-3 
    Yrml=1/4*sqrt(35/(2*pi))*y*(3*x^2-y^2)/r^3;
end
if m==-2
    Yrml=1/2*sqrt(105/pi)*x*y*z/r^3;
end
if m==-1
    Yrml=1/4*sqrt(21/(2*pi))*y*(5*z^2-r^2)/r^3;
end
if m==0 
    Yrml= 1/4*sqrt(7/pi)*(5*z^3-3*z*r^2)/r^3;
end
if m==1
    Yrml= 1/4*sqrt(21/(2*pi))*x*(5*z^2-r^2)/r^3;
end
if m==2
    Yrml=1/4*sqrt(105/pi)*(x^2-y^2)*z/r^3;
end
if m==3
    Yrml=1/4*sqrt(35/(2*pi))*x*(x^2-3*y^2)/r^3;
end
end
%%Falta l=4;
if l==4
    if m==-4
        Yrml=3/4*sqrt(35/pi)*x*y*(x^2-y^2)/r^4;
    end
    if m==-3
        Yrml=3/4*sqrt(35/(2*pi))*y*(3*x^2-y^2)*z/r^4;
    end
    if m==-2
        Yrml=3/4*sqrt(5/pi)*x*y*(7*z^2-r^2)/r^4;

    end
    if m==-1
        Yrml=3/4*sqrt(5/(2*pi))*y*(7*z^3-3*z*r^2)/r^4;

    end
    if m==0
        Yrml=3/16*sqrt(1/pi)*(35*z^4-30*z^2*r^2+3*r^4)/r^4;
    end
    if m==1
        Yrml=3/4*sqrt(5/(2*pi))*x*(7*z^3-3*z*r^2)/r^4;

    end
    if m==2
        Yrml=3/8*sqrt(5/pi)*(x^2-y^2)*(7*z^2-r^2)/r^4;
    end

    if m==3
        Yrml=3/4*sqrt(35/(2*pi))*x*(x^2-3*y^2)*z/r^4;
    end
    
    if m==4 
        Yrml=3/16*sqrt(35/pi)*(x^2*(x^2-3*y^2)-y^2*(3*x^2-y^2))/r^4;

    end

end
% Caso l=5
 if l == 5
 if m == -5
                Yrml = (3/16) * sqrt(77/pi) * sin(theta)^5 * sin(5*phi);
            elseif m == -4
                Yrml = (3/16) * sqrt(385/(2*pi)) * sin(theta)^4 * cos(theta) * sin(4*phi);
            elseif m == -3
                Yrml = (1/32) * sqrt(385/pi) * sin(theta)^3 * (9*cos(theta)^2 - 1) * sin(3*phi);
            elseif m == -2
                Yrml = (1/8) * sqrt(1155/pi) * sin(theta)^2 * (3*cos(theta)^3 - cos(theta)) * sin(2*phi);
            elseif m == -1
                Yrml = (1/16) * sqrt(165/(2*pi)) * sin(theta) * (21*cos(theta)^4 - 14*cos(theta)^2 + 1) * sin(phi);
            elseif m == 0
                Yrml = (1/16) * sqrt(11/pi) * (63*cos(theta)^5 - 70*cos(theta)^3 + 15*cos(theta));
            elseif m == 1
                Yrml = (1/16) * sqrt(165/(2*pi)) * sin(theta) * (21*cos(theta)^4 - 14*cos(theta)^2 + 1) * cos(phi);
            elseif m == 2
                Yrml = (1/8) * sqrt(1155/pi) * sin(theta)^2 * (3*cos(theta)^3 - cos(theta)) * cos(2*phi);
            elseif m == 3
                Yrml = (1/32) * sqrt(385/pi) * sin(theta)^3 * (9*cos(theta)^2 - 1) * cos(3*phi);
            elseif m == 4
                Yrml = (3/16) * sqrt(385/(2*pi)) * sin(theta)^4 * cos(theta) * cos(4*phi);
            elseif m == 5
                Yrml = (3/16) * sqrt(77/pi) * sin(theta)^5 * cos(5*phi);
            else
                Yrml = 0;
 end
 end   
 if l==6
    if m == -6
                Yrml = (1/32) * sqrt(3003/pi) * sin(theta)^6 * sin(6*phi);
            elseif m == -5
                Yrml = (3/32) * sqrt(1001/pi) * sin(theta)^5 * cos(theta) * sin(5*phi);
            elseif m == -4
                Yrml = (3/32) * sqrt(91/(2*pi)) * sin(theta)^4 * (11*cos(theta)^2 - 1) * sin(4*phi);
            elseif m == -3
                Yrml = (1/32) * sqrt(1365/pi) * sin(theta)^3 * (11*cos(theta)^3 - 3*cos(theta)) * sin(3*phi);
            elseif m == -2
                Yrml = (1/64) * sqrt(1365/pi) * sin(theta)^2 * (33*cos(theta)^4 - 18*cos(theta)^2 + 1) * sin(2*phi);
            elseif m == -1
                Yrml = (1/16) * sqrt(273/(2*pi)) * sin(theta) * (33*cos(theta)^5 - 30*cos(theta)^3 + 5*cos(theta)) * sin(phi);
            elseif m == 0
                Yrml = (1/32) * sqrt(13/pi) * (231*cos(theta)^6 - 315*cos(theta)^4 + 105*cos(theta)^2 - 5);
            elseif m == 1
                Yrml = (1/16) * sqrt(273/(2*pi)) * sin(theta) * (33*cos(theta)^5 - 30*cos(theta)^3 + 5*cos(theta)) * cos(phi);
            elseif m == 2
                Yrml = (1/64) * sqrt(1365/pi) * sin(theta)^2 * (33*cos(theta)^4 - 18*cos(theta)^2 + 1) * cos(2*phi);
            elseif m == 3
                Yrml = (1/32) * sqrt(1365/pi) * sin(theta)^3 * (11*cos(theta)^3 - 3*cos(theta)) * cos(3*phi);
            elseif m == 4
                Yrml = (3/32) * sqrt(91/(2*pi)) * sin(theta)^4 * (11*cos(theta)^2 - 1) * cos(4*phi);
            elseif m == 5
                Yrml = (3/32) * sqrt(1001/pi) * sin(theta)^5 * cos(theta) * cos(5*phi);
            elseif m == 6
                Yrml = (1/32) * sqrt(3003/pi) * sin(theta)^6 * cos(6*phi);
            else
                Yrml = 0;
   end
end
    
    % l = 7
 if l==7
    if m == -7
                Yrml = (1/64) * sqrt(715/(2*pi)) * sin(theta)^7 * sin(7*phi);
            elseif m == -6
                Yrml = (1/64) * sqrt(9009/(2*pi)) * sin(theta)^6 * cos(theta) * sin(6*phi);
            elseif m == -5
                Yrml = (1/64) * sqrt(1155/(2*pi)) * sin(theta)^5 * (13*cos(theta)^2 - 1) * sin(5*phi);
            elseif m == -4
                Yrml = (3/64) * sqrt(385/(2*pi)) * sin(theta)^4 * (13*cos(theta)^3 - 3*cos(theta)) * sin(4*phi);
            elseif m == -3
                Yrml = (3/64) * sqrt(385/(2*pi)) * sin(theta)^3 * (143*cos(theta)^4 - 66*cos(theta)^2 + 3) * sin(3*phi);
            elseif m == -2
                Yrml = (3/64) * sqrt(35/(2*pi)) * sin(theta)^2 * (143*cos(theta)^5 - 110*cos(theta)^3 + 15*cos(theta)) * sin(2*phi);
            elseif m == -1
                Yrml = (1/64) * sqrt(105/(2*pi)) * sin(theta) * (429*cos(theta)^6 - 495*cos(theta)^4 + 135*cos(theta)^2 - 5) * sin(phi);
            elseif m == 0
                Yrml = (1/32) * sqrt(15/pi) * (429*cos(theta)^7 - 693*cos(theta)^5 + 315*cos(theta)^3 - 35*cos(theta));
            elseif m == 1
                Yrml = (1/64) * sqrt(105/(2*pi)) * sin(theta) * (429*cos(theta)^6 - 495*cos(theta)^4 + 135*cos(theta)^2 - 5) * cos(phi);
            elseif m == 2
                Yrml = (3/64) * sqrt(35/(2*pi)) * sin(theta)^2 * (143*cos(theta)^5 - 110*cos(theta)^3 + 15*cos(theta)) * cos(2*phi);
            elseif m == 3
                Yrml = (3/64) * sqrt(385/(2*pi)) * sin(theta)^3 * (143*cos(theta)^4 - 66*cos(theta)^2 + 3) * cos(3*phi);
            elseif m == 4
                Yrml = (3/64) * sqrt(385/(2*pi)) * sin(theta)^4 * (13*cos(theta)^3 - 3*cos(theta)) * cos(4*phi);
            elseif m == 5
                Yrml = (1/64) * sqrt(1155/(2*pi)) * sin(theta)^5 * (13*cos(theta)^2 - 1) * cos(5*phi);
            elseif m == 6
                Yrml = (1/64) * sqrt(9009/(2*pi)) * sin(theta)^6 * cos(theta) * cos(6*phi);
            elseif m == 7
                Yrml = (1/64) * sqrt(715/(2*pi)) * sin(theta)^7 * cos(7*phi);
            else
                Yrml = 0;
    end
 end 
    % Si no se encontró ningún caso válido
    %Yrml = 0;



if r==0
    Yrml=0;
end
end


