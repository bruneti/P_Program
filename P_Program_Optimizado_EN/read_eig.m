%%
%% This function reads from file 'label.EIG' the Fermi energy and 
%% eigen-energies and retrieves them in sub-structures
%%
%% EF                   = Fermi energy
%% eigen_energies(ik,1:no_s*Nspin) = eigenvalues for each k
%%
%% Written by J. Ferrer, October 2020
%%

function [ siesta ] = read_eig( flags, siesta )

feig   = [ siesta.system_label '.EIG' ];

%% Read file
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