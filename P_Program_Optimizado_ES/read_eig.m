%%
%% Lee del fichero 'label.EIG' la energía de Fermi y las autoenergías, y las
%% guarda en sub-estructuras.
%%
%% EF                   = energía de Fermi
%% eigen_energies(ik,1:no_s*Nspin) = autovalores para cada k
%%
%% Escrito por J. Ferrer, octubre 2020
%%

function [ siesta ] = read_eig( flags, siesta )

feig   = [ siesta.system_label '.EIG' ];

%% Lee el fichero
opt = {'CollectOutput',true};
out = {};
[ fid ] = fopen(feig,'rt');

siesta.EF = str2double(fgetl(fid));
nn     = strsplit(strtrim(fgetl(fid)));
no_s   = str2double(nn{1});
spin   = str2double(nn{2});
neigen = no_s*spin;
if spin > 2, neigen = no_s; end

while ~feof(fid)
    out(end+1) = textscan(fid,'%f',opt{:});
end
fclose(fid);

for ik = 1:siesta.kpoints.Nk
  siesta.eigen_energies(ik,1:neigen)=out{1}((neigen+1)*(ik-1)+(2:neigen+1));
end

end