%%
%% This function computes the Kohn-Sham bands along all lines 
%% each two consecutive k-points.
%%
%% Reference energy is EF
%%
%% Written by J. Ferrer, October 2020
%

function [ post ] = post_bands( post, siesta, flags )

kline.vectors = flags.bands.kline.vectors;
Nk            = flags.bands.kline.Nk;
EF            = post.EF;

k = [];
for i = 1:size(kline.vectors,1)-1
  k1     = kline.vectors(i,:);
  k2     = kline.vectors(i+1,:);
  lambda = (0:1/(Nk(i)-1):1)';
  k      = [k; k1+lambda*(k2-k1)];
end
k = 2*pi*k;
post.k=k;
%% Get bands along the asked line
for ik = 1:size(k,1)
  if strcmp( flags.full_OHDM, 'yes')
    [ Hk, Ok ] = util_Hk( k(ik,:), siesta );
  else
    [ Hk, Ok ] = util_Hk_sparse( k(ik,:), siesta );  
  end    
  if siesta.spin.Nspin == 2
  	bands.up  (:,ik) = eig(Hk(:,:,1),Ok)-EF;
    bands.down(:,ik) = eig(Hk(:,:,2),Ok)-EF;
  else
    bands(:,ik) = eig(Hk,Ok)-EF;
  end
end

if siesta.spin.Nspin == 2
  bands.up         = permute( bands.up, [ 2 1 ] );
  bands.down       = permute( bands.down, [ 2 1 ] );
else
  bands            = permute( bands, [ 2 1 ] );
end

post.bands = bands;

end
