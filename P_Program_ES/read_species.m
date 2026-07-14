%%
%% Busca en el fichero de salida principal de SIESTA ('logfile') la
%% información relevante de especies y la guarda en las sub-estructuras
%% correspondientes (símbolo, Z, masa, PAOs, ...).
%%
%% Escrito por J. Ferrer, octubre 2020, enero 2022
%%

function [ siesta ] = read_species( logfile, siesta, flags );

%% Lista de masas atómicas
  amass =[  1.01,   4.00,   6.94,   9.01,  10.81,  12.01, ...
           14.01,  16.00,  19.00,  20.18,  22.99,  24.31,  ...
           26.98,  28.09,  30.97,  32.07,  35.45,  39.95, ...
           39.10,  40.08,  44.96,  47.88,  50.94,  52.00, ...
           54.94,  55.85,  58.93,  58.69,  63.55,  65.39, ...
           69.72,  72.61,  74.92,  78.96,  79.90,  83.80, ...
           85.47,  87.62,  88.91,  91.22,  92.91,  95.94, ...
           98.91, 101.07, 102.91, 106.42, 107.87, 112.41, ...
          114.82, 118.71, 121.75, 127.60, 126.90, 131.29, ...
          132.91, 137.33, 138.91, 140.12, 140.91, 144.24, ...
          146.92, 150.36, 151.97, 157.25, 158.93, 162.50, ...
          164.93, 167.26, 168.93, 173.04, 174.97, 178.49, ...
          180.95, 183.85, 186.21, 190.20, 192.22, 195.08, ...
          196.97, 200.59, 204.38, 207.20, 208.98, 208.98, ...
          209.99, 222.02, 223.02, 226.03, 227.03, 232.04, ...
          231.04, 238.03, 237.05, 244.06];

%% Número de especies y código interno de especie de SIESTA
  foundat   = strfind(logfile,'NumberOfSpecies');
  dummy     = strsplit(logfile(foundat:foundat+40));
  Nspecies  = str2double(dummy{2});
  foundat1  = strfind(logfile,'%block ChemicalSpeciesLabel');
  foundat2  = strfind(logfile,'%endblock ChemicalSpeciesLabel');
  dummy     = strsplit(strtrim(logfile(foundat1+29:foundat2)),'\n');

%% Lee la información de especies
  siesta.Nspecies = Nspecies;
  for i = 1:Nspecies
    dummy2            = strsplit(strtrim(dummy{i}));
    species(i).symbol = string(dummy2{3});
    species(i).number = str2double(dummy2{1});
    species(i).Z      = str2double(dummy2{2});
    species(i).mass   = amass( species(i).Z );
  end

%% Lee los datos de PAO de cada especie
  if strcmp( flags.netcdf, 'no' )
    for i = 1:siesta.Nspecies
      ionfile = fileread( [ char(species(i).symbol) '.ion' ] );
      foundat = strfind( ionfile, 'Self energy' );
      dummy   = splitlines( ionfile( foundat:foundat+35));
      dummy   = strsplit( dummy{2});
      species(i).PAO.nPAO = str2num(dummy{3});
      foundat = strfind( ionfile, '# npts, delta, cutoff' );
      data = [];
      dummy2 = splitlines(ionfile(foundat(i)+1:foundat(i+1)-100));
      tmp2 = [];
      for j=1:500
        dummy = strsplit(dummy2{j+1});
        tmp2 = [ tmp2; str2num(dummy{2})];
      end
      for j=1:species(i).PAO.nPAO
        dummy1 = splitlines( ionfile(foundat(j)-150:foundat(j)));
        dummy  = strsplit(dummy1{2});
        n      = str2num(dummy{3});
        l(j)   = str2num(dummy{2});
        z      = str2num(dummy{4});
        pol    = str2num(dummy{5});
        pop    = str2num(dummy{6});
        data   = [ data [ n; l(j); z; pol; pop ] ];
        dummy  = strsplit(dummy1{3});
        Rc(j)  = str2num(dummy{4});
        dummy2 = splitlines(ionfile(foundat(j)+1:foundat(j+1)-100));
        tmp1 = [];
        for k=1:500
          dummy = strsplit(dummy2{k+1});
          tmp1 = [ tmp1; str2num(dummy{3})];
        end
        tmp2 = [ tmp2 tmp1 ];
      end
      species(i).PAO.PAOs   = [ unique(l); histc(l,unique(l)) ];
      species(i).PAO.ordering = l;
      species(i).PAO.data = data;
      species(i).PAO.Rc   = Rc;
      species(i).PAO.Rlz  = tmp2;
    end
  else
    for i = 1:siesta.Nspecies
      fid  = netcdf.open( [ char(species(i).symbol) '.ion.nc' ] ,'NC_NOWRITE');
      l    = netcdf.getVar( fid, 0 );
      n    = netcdf.getVar( fid, 1 );
      z    = netcdf.getVar( fid, 2 );
      pol  = netcdf.getVar( fid, 3 );
      pop  = netcdf.getVar( fid, 4, 'double' );
      Rc   = netcdf.getVar( fid, 5, 'double' );
      dr   = netcdf.getVar( fid, 6 );
      Rlz  = netcdf.getVar( fid, 7 );
      nPAO = size(n,1);
      species(i).PAO.nPAO     = size(l,1);
      species(i).PAO.PAOs     = [ unique(l) histc(l,unique(l)) ];
      species(i).PAO.ordering = l';
      species(i).PAO.data     = [ n l z pol pop ];
      species(i).PAO.grid     = [ Rc dr];
      for iPAO = 1:nPAO
        r(:,iPAO)    = (0:dr(iPAO):Rc(iPAO));
        species(i).PAO.Rlz(iPAO).Rlz = [ r(:,iPAO) r(:,iPAO).^double(l(iPAO)).*Rlz(:,iPAO) ];
      end
    end
  end

%% Vuelca todo en la estructura siesta.species
  siesta.species = species; 

end
