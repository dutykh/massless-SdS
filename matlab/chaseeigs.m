% Polynomial eigenvalue solver for massless SdS QNMs.
% Quadratic problem  M0 + i*Omega*M1 + Omega^2*M2 = 0  solved with Advanpix
% multiprecision  polyeig.  Reads parameters and resolutions from shared
% config files (params.txt, resolutions.txt), so the same script runs both
% in single/ and inside any runs/<case>/ directory.
%
% Author: Dr. Denys Dutykh (Khalifa University of Science and Technology,
%         Abu Dhabi, UAE)

close
clear
format longE

addpath('/home/user/Soft/Advanpix/');
maxNumCompThreads(30);

list = load('resolutions.txt')';
N    = length(list);

% Read parameters (s, L, Lambda) from shared config:
fid = fopen('params.txt', 'r');
s_val      = mp(strtrim(fgetl(fid)));
L_val      = mp(strtrim(fgetl(fid)));
Lambda_val = mp(strtrim(fgetl(fid)));
fclose(fid);

for idx = 1:N
  n    = list(idx);
  nstr = num2str(n);

  fprintf('Computing eigs for n = %3d ... ', n);

  mp.Digits(n);

  M0 = mp.read(strcat('assemble/M0_', nstr, '.mat'));
  M1 = mp.read(strcat('assemble/M1_', nstr, '.mat'));
  M2 = mp.read(strcat('assemble/M2_', nstr, '.mat'));

  % Quadratic polynomial eigenvalue problem (massless case):
  e = polyeig(M0, mp('1i')*M1, M2);

  % Save QNMs to a platform-agnostic text file at full multiprecision:
  fname = strcat('results/eigs_', nstr, '.dat');
  fid = fopen(fname, 'w');
  for k = 1:length(e)
    fprintf(fid, '%s %s\n', num2str(real(e(k)), n), num2str(imag(e(k)), n));
  end
  fclose(fid);

  fprintf('done.\n');
end % for
