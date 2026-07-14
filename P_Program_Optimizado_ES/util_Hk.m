%%
%% Calcula H_k y S_k para un k-point dado.
%% La componente de espín tiene sus dimensiones adecuadas:
%% para Nspin = 1,  Hk = Hk(1:nou,1:nou);
%%      Nspin = 2,  Hk = Hk(1:nou,1:nou, 1:2);
%%      Nspin = 4,8 Hk = Hk(1:2*nou,1:2*nou);
%%
%% Escrito por J. Ferrer, octubre 2020
%

function [ Hk, Ok ] = util_Hk( k, siesta )

%% Abreviaturas
O   = siesta.Osc;
H   = siesta.Hsc;
R   = siesta.lattice.supercell.R;
nsc = siesta.lattice.supercell.Ncells;

%% Fase de Bloch para todas las celdas unidad
eikr = reshape( exp(-1i*k*R), [ 1 1 nsc ] );

%% Calcula Ok
Ok = sum( bsxfun( @times, O, eikr), 3);
Ok = (Ok+Ok')/2;

%% Calcula Hk
if siesta.spin.Nspin == 2
  Hk(:,:,1) = sum(bsxfun(@times,H(1).sp,eikr), 3);
  Hk(:,:,1) = (Hk(:,:,1)+Hk(:,:,1)')/2;
  Hk(:,:,2) = sum(bsxfun(@times,H(2).sp,eikr), 3);
  Hk(:,:,2) = (Hk(:,:,2)+Hk(:,:,2)')/2;
else
  Hk = sum(bsxfun(@times,H,eikr),3);
  Hk = (Hk+Hk')/2;
end

end
