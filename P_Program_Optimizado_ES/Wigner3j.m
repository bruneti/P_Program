function w = Wigner3j( j123, m123 )

% Calcula el símbolo 3j de Wigner mediante la fórmula de Racah.
%
% W = Wigner3j( J123, M123 )
%
% J123 = [J1, J2, J3].
% M123 = [M1, M2, M3].
% Todos los Ji y Mi deben ser enteros o semienteros (según corresponda).
%
% Por las reglas de selección, W = 0 salvo que:
%   |Ji - Jj| <= Jk <= (Ji + Jj)    (i,j,k permutaciones de 1,2,3)
%   |Mi| <= Ji    (i = 1,2,3)
%    M1 + M2 + M3 = 0
%
% Referencia:
% Entrada "Wigner 3j-Symbol" de Mathworld (Eric Weinstein):
% http://mathworld.wolfram.com/Wigner3j-Symbol.html
%
% Inspirado en Wigner3j.m de David Terr, Raytheon, 6-17-04
%  (disponible en www.mathworks.com/matlabcentral/fileexchange).
%
% Por Kobi Kraus, Technion, 25-6-08.
% Actualizado 1-8-13.

j1 = j123(1); j2 = j123(2); j3 = j123(3);
m1 = m123(1); m2 = m123(2); m3 = m123(3);

% Comprobación de errores de entrada
if any( j123 < 0 ),
    error( 'Los j deben ser no negativos' )
elseif any( rem( [j123, m123], 0.5 ) ),
    error( 'Todos los argumentos deben ser enteros o semienteros' )
elseif any( rem( (j123 - m123), 1 ) )
    error( 'j123 y m123 no concuerdan' );
end

% Reglas de selección
if ( j3 > (j1 + j2) ) || ( j3 < abs(j1 - j2) ) ... % j3 fuera de intervalo
   || ( m1 + m2 + m3 ~= 0 ) ... % no conserva el momento angular
   || any( abs( m123 ) > j123 ), % m mayor que j
    w = 0;
    return
end
    
% Caso común sencillo
if ~any( m123 ) && rem( sum( j123 ), 2 ), % m1 = m2 = m3 = 0 y j1 + j2 + j3 impar
    w = 0;
    return
end

% Evaluación
t1 = j2 - m1 - j3;
t2 = j1 + m2 - j3;
t3 = j1 + j2 - j3;
t4 = j1 - m1;
t5 = j2 + m2;

tmin = max( 0,  max( t1, t2 ) );
tmax = min( t3, min( t4, t5 ) );

t = tmin : tmax;
w = sum( (-1).^t .* exp( -ones(1,6) * gammaln( [t; t-t1; t-t2; t3-t; t4-t; t5-t] +1 ) + ...
                         gammaln( [j1+j2+j3+1, j1+j2-j3, j1-j2+j3, -j1+j2+j3, j1+m1, j1-m1, j2+m2, j2-m2, j3+m3, j3-m3] +1 ) ...
                         * [-1; ones(9,1)] * 0.5 ) ) * (-1)^( j1-j2-m3 );
         
% Avisos
if isnan( w )
    warning( 'MATLAB:Wigner3j:NaN', '¡Wigner3J es NaN!' )
elseif isinf( w )
    warning( 'MATLAB:Wigner3j:Inf', '¡Wigner3J es Inf!' )
end