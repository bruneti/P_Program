%%
%% This function computes H_k and S_k for a given k-point, 
%% The spin component has its adequate dimensions:
%% for Nspin = 1,  Hk = Hk(1:nou,1:nou);
%%     Nspin = 2,  Hk = Hk(1:nou,1:nou, 1:2);
%%     Nspin = 4,8 Hk = Hk(1:2*nou,1:2*nou);
%%
%% Written by J. Ferrer, October 2020
%

function [ Hk, Ok ] = util_Hk( k, siesta )

%% Define a few short-hands
O   = siesta.Osc;
H   = siesta.Hsc;
R   = siesta.lattice.supercell.R;
nsc = siesta.lattice.supercell.Ncells;

%% Define Bloch phase for all unit cells
eikr = reshape( exp(-1i*k*R), [ 1 1 nsc ] );

%% Compute Ok
Ok = sum( bsxfun( @times, O, eikr), 3);
Ok = (Ok+Ok')/2;

%% Compute Hk
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
