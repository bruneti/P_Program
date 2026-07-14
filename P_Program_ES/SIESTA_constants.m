%% SIESTA_constants.m
%
% Constantes físicas (SI salvo indicación) utilizadas a lo largo del
% pipeline de tasa de transición dipolar.
%
% Extraído de GROGU (J. Ferrer, 2020).

function constants = SIESTA_constants

constants.q      = 1.6022e-19;                      % Coulomb
constants.h      = 6.6261e-34;                      % Joule*sec
constants.hbar   = constants.h/(2*pi);              % Joule*sec
constants.kB     = 1.3807e-23;                      % Joule/Kelvin
constants.eV     = 1.6022e-19;                      % Joule
constants.G0     = 2*constants.q^2/constants.h;     % Siemen
constants.RytoeV = 13.60580;                        % Rydberg to eV
constants.eVtoK  = 1e5/8.617;                       % eV to K
constants.Bohr   = 0.529177;                        % Bohr to Angstrom

end
