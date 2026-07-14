% 
% Reads H, O and DM from Siesta.
%
%  Written by D. Visontai 2017  & J. Ferrer, October 2020:
%

function [ siesta ] = read_hodm( siesta, constants, flags )

%% AUXILIARY VARIABLES
  tol  = flags.tol;

%% READ H, O & DM
  if strcmp(flags.netcdf,'no')
    fho     = [ siesta.system_label '.HSX' ];
    fdm     = [ siesta.system_label '.DM' ];
  
    %% OPEN 'HSX' FILE & READ SOME STUFF. ONLY nou, nos AND Nspin ARE USED
    fid   = fopen(fho, 'rb');

    hdr   = fread(fid, 1, 'int32');
    nou   = fread(fid, 1, 'int32');
    nos   = fread(fid, 1, 'int32');
    Nspin = fread(fid, 1, 'int32');  
    nh    = fread(fid, 1, 'int32');
    trl   = fread(fid, 1, 'int32');

    nsc = nos/nou;

    hdr   = fread(fid, 1, 'int32');
    gamma = fread(fid, 1, 'int32');
    trl   = fread(fid, 1, 'int32');

    if ~ gamma
      hdr    = fread(fid, 1,   'int32');
      indxuo = fread(fid, nos, 'int32');
      trl    = fread(fid, 1,   'int32');
    end

    hdr  = fread(fid, 1,   'int32');
    numh = fread(fid, nou, 'int32');
    trl  = fread(fid, 1,   'int32');

    %% BUILD MATRIX SUBINDICES IUO & JO
    iuo=[];
    for io=1:nou
	    iuo = [iuo; io.*ones(numh(io),1)];
    end

    jo = [];
    for  io = 1:nou
      hdr   = fread(fid, 1, 'int32') ;
      ibuff = fread(fid, numh(io), 'int32');
      trl   = fread(fid, 1, 'int32') ;
      jo    = [ jo; ibuff ];
    end

    %% READ HAMILTONIAN
    for is = 1:Nspin
      haux = [];
      for io=1:nou
	      hdr   = fread(fid, 1, 'int32');
	      hbuff = fread(fid, numh(io), flags.precision.HO );
	      trl   = fread(fid, 1, 'int32') ;
        haux  = [ haux; hbuff ];
      end
      h(:,is) = haux;
    end
    h(abs(h)<tol) = 0;
    h = constants.RytoeV*h;

    %% READ OVERLAP
    o = [];
    for io = 1:nou
      hdr   = fread(fid, 1, 'int32');
      obuff = fread(fid, numh(io), flags.precision.HO );
      trl   = fread(fid, 1, 'int32') ;
      o     = [ o; obuff ];
    end
    o(abs(o)<tol) = 0;

    hdr  = fread(fid, 1, 'int32');
    Qtot = fread(fid, 1, 'float64');
    T    = constants.RytoeV*constants.eVtoK*fread(fid, 1, 'float64');
    trl  = fread(fid, 1, 'int32');

    fclose(fid); 

    %% OPEN 'DM' FILE & READ SOME STUFF. ONLY nou AND Nspin ARE USED
    fid = fopen(fdm, 'rb');

    hdr   = fread(fid, 1, 'int32');
    nou   = fread(fid, 1, 'int32');
    Nspin = fread(fid, 1, 'int32');
    trl   = fread(fid, 1, 'int32');

    %% WARNING: 3 OF THE 4 NEXT HDR SHOULD BE COMMENTED OUT FOR SIESTA 
    %% VERSIONS EARLIER THAN SIESTA-4.1.5
    hdr   = fread(fid, 1,   'int32');
    hdr   = fread(fid, 1,   'int32');
    hdr   = fread(fid, 1,   'int32');
    hdr   = fread(fid, 1,   'int32');
    numd  = fread(fid, nou, 'int32');
    trl   = fread(fid, 1,   'int32');

    %% BUILD MATRIX SUBINDICES IUO & JO
    iuo=[];
    for io=1:nou
      iuo = [iuo; io.*ones(numd(io),1)];
    end

    jo = [];
    for  io = 1:nou
      hdr   = fread(fid, 1,        'int32');
      ibuff = fread(fid, numd(io), 'int32');
      trl   = fread(fid, 1,        'int32');
      jo    = [jo; ibuff ];
    end

    %% READ DENSITY MATRIX
    for is=1:Nspin
      dmaux = [];
      for io=1:nou
        hdr = fread(fid, 1,          'int32') ;
        hbuff = fread(fid, numd(io), flags.precision.DM) ;
        trl = fread(fid, 1,          'int32') ;
        dmaux = [ dmaux; hbuff ];
      end
      dm(:,is) = dmaux;
    end
    dm(abs(dm)<tol) = 0;

    fclose(fid); 

  else

    fdmho  = 'DMHS.nc';

    %% OPEN 'DMHS.nc' FILE & READ nou, nos, Nspin, AND nh
    fileid = netcdf.open(fdmho,'NC_NOWRITE');

    [ dummy, nou ]   = netcdf.inqDim( fileid, 0 );
    [ dummy, nos ]   = netcdf.inqDim( fileid, 1 );
    [ dummy, Nspin ] = netcdf.inqDim( fileid, 2 );
    [ dummy, nh ]    = netcdf.inqDim( fileid, 3 );

    nsc = nos/nou;

    numh = netcdf.getVar(fileid,0);
    rowp = netcdf.getVar(fileid,1);
    jo   = netcdf.getVar(fileid,2,'double');

    %% BUILD MATRIX SUBINDICES IUO
    iuo = [];
    for io = 1:nou
      iuo = [iuo; io.*ones(numh(io),1)];
    end

    %% READ O, H and DM
    o               = netcdf.getVar(fileid,3,'double');
    dm              = netcdf.getVar(fileid,5,'double');
    h               = netcdf.getVar(fileid,6,'double');
    o(abs(o)<tol)   = 0;
    dm(abs(dm)<tol) = 0;
    h(abs(h)<tol)   = 0;
    h               = constants.RytoeV*h;

    netcdf.close(fileid);
  end

%% DEFINE SPIN STUFF
if     Nspin == 1
    spinflavor = 'This simulation is Non magnetic';
elseif Nspin == 2
    spinflavor = 'This simulation is Spin Polarized';
elseif Nspin == 4
    spinflavor = 'This simulation is Spin Non-Collinear';
elseif Nspin == 8
    spinflavor = 'This simulation includes Spin-Orbit interactions';
end

spin.spinflavor = spinflavor;
spin.Nspin      = Nspin;

%% CONVERT OVERLAP MATRIX FROM SIESTA'S TO REGULAR RECTANGULAR FORMAT &
%% ADD SPIN STRUCTURE
O = sparse(jo, iuo, o, nos, nou);
if Nspin > 2, O = kron( O, eye(2) ); end

%% PAULI-LIKE MATRICES
tau1 =     [1 0; 0 0];
tau2 =     [0 0; 0 1];
tau3 =     [0 1; 0 0]; 
tau7 =     [0 0; 1 0];
tau4 = -1i * tau3;
tau5 =  1i * tau1;
tau6 =  1i * tau2;
tau8 =  1i * tau7;

%% CONVERT HAMILTONIAN MATRIX FROM SIESTA'S TO REGULAR RECTANGULAR FORMAT
%% ADD SPIN STRUCTURE
if Nspin == 1
  H = sparse(jo, iuo, h, nos, nou);
elseif Nspin == 2
  H(1).sp = sparse(jo, iuo, h(:,1), nos, nou);
  H(2).sp = sparse(jo, iuo, h(:,2), nos, nou);
elseif Nspin == 4
  H1 = sparse(jo, iuo, h(:,1), nos, nou);
  H2 = sparse(jo, iuo, h(:,2), nos, nou);
  H3 = sparse(jo, iuo, h(:,3), nos, nou);
  H4 = sparse(jo, iuo, h(:,4), nos, nou);
  H  = kron(H1,tau1) + kron(H2,tau2) + kron(H3,tau3+tau7) + kron(H4,tau4+tau8); 
elseif Nspin == 8
  H1 = sparse(jo, iuo, h(:,1), nos, nou);
  H2 = sparse(jo, iuo, h(:,2), nos, nou);
  H3 = sparse(jo, iuo, h(:,3), nos, nou);
  H4 = sparse(jo, iuo, h(:,4), nos, nou);
  H5 = sparse(jo, iuo, h(:,5), nos, nou);
  H6 = sparse(jo, iuo, h(:,6), nos, nou);
  H7 = sparse(jo, iuo, h(:,7), nos, nou);
  H8 = sparse(jo, iuo, h(:,8), nos, nou);
  H  = kron(H1,tau1) + kron(H2,tau2) + kron(H3,tau3) + kron(H4,tau4) + ...
       kron(H5,tau5) + kron(H6,tau6) + kron(H7,tau7) + kron(H8,tau8);
end

%% CONVERT DENSITY MATRIX FROM SIESTA'S TO REGULAR RECTANGULAR FORMAT 
%% THE AWKWARD SPIN CONVERSION IS DUE TO SIESTA'S AWKWARD DEFINITIONS
if Nspin == 1
  DM = sparse(jo, iuo, dm, nos, nou);
elseif Nspin == 2
  DM(1).sp = sparse(jo, iuo, dm(:,1), nos, nou);
  DM(2).sp = sparse(jo, iuo, dm(:,2), nos, nou);
elseif Nspin == 4
  DM1 = sparse(jo, iuo, dm(:,1), nos, nou);
  DM2 = sparse(jo, iuo, dm(:,2), nos, nou);
  DM3 = sparse(jo, iuo, dm(:,3), nos, nou);
  DM4 = sparse(jo, iuo, dm(:,4), nos, nou);
  DM  = kron(DM1,tau1) + kron(DM2,tau2) + kron(DM3,tau3) + kron(DM4,tau4);
elseif Nspin == 8
  DM1 =   sparse(jo, iuo, dm(:,1), nos, nou);
  DM2 =   sparse(jo, iuo, dm(:,2), nos, nou);
  DM3 =   sparse(jo, iuo, dm(:,7), nos, nou);
  DM4 =   sparse(jo, iuo, dm(:,8), nos, nou);
  DM5 = - sparse(jo, iuo, dm(:,5), nos, nou);
  DM6 = - sparse(jo, iuo, dm(:,6), nos, nou);
  DM7 =   sparse(jo, iuo, dm(:,3), nos, nou);
  DM8 =   sparse(jo, iuo, dm(:,4), nos, nou);
  DM  =   kron(DM1,tau1) + kron(DM2,tau2) + kron(DM3,tau3) + kron(DM4,tau4) + ...
          kron(DM5,tau5) + kron(DM6,tau6) + kron(DM7,tau7) + kron(DM8,tau8);
end

%% CONVERT TALL TO WIDE MATRIX
if Nspin > 2, nou = 2*nou; end;
if Nspin == 2
  O        = sparse( reshape( permute( reshape( full(O).'       , [ nou nou nsc ] ), [ 2 1 3 ] ), [ nou nou*nsc ] ) );
  H(1).sp  = sparse( reshape( permute( reshape( full(H(1).sp).' , [ nou nou nsc ] ), [ 2 1 3 ] ), [ nou nou*nsc ] ) );
  H(2).sp  = sparse( reshape( permute( reshape( full(H(2).sp).' , [ nou nou nsc ] ), [ 2 1 3 ] ), [ nou nou*nsc ] ) );
  DM(1).sp = sparse( reshape( permute( reshape( full(DM(1).sp).', [ nou nou nsc ] ), [ 2 1 3 ] ), [ nou nou*nsc ] ) );
  DM(2).sp = sparse( reshape( permute( reshape( full(DM(2).sp).', [ nou nou nsc ] ), [ 2 1 3 ] ), [ nou nou*nsc ] ) );
else
  O =  sparse( reshape( permute( reshape( full(O).',  [ nou nou nsc ] ), [ 2 1 3 ] ), [ nou nou*nsc ] ) );
  H  = sparse( reshape( permute( reshape( full(H).',  [ nou nou nsc ] ), [ 2 1 3 ] ), [ nou nou*nsc ] ) );
  DM = sparse( reshape( permute( reshape( full(DM).', [ nou nou nsc ] ), [ 2 1 3 ] ), [ nou nou*nsc ] ) );
end

%% MAKE THEM FULL WITH DIMENSIONS (nou,nou,nsc) IF DESIRED
if strcmp(flags.full_OHDM,'yes')
  if Nspin == 2
    siesta.Osc = reshape( full(O), [ nou nou nsc ] );
    for ispin = 1:2
      siesta.Hsc(ispin).sp  = reshape( full(H(ispin).sp), [ nou nou nsc ] );
      siesta.DMsc(ispin).sp = reshape( full(DM(ispin).sp), [ nou nou nsc ] );
    end
  else
    siesta.Osc  = reshape( full(O),  [ nou nou nsc ] );
    siesta.Hsc  = reshape( full(H),  [ nou nou nsc ] );
    siesta.DMsc = reshape( full(DM), [ nou nou nsc ] );
  end

end

siesta.O    = O;
siesta.H    = H;
siesta.DM   = DM;
siesta.spin = spin;

end
