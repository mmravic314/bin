%  -----------------------------------------------------------------------
%  Copyright (C) 2011 Gevorg Grigoryan
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%  -----------------------------------------------------------------------
%
%FCRICK Fit Crick parameters to an input structure.
%   [err xyz p] = fcrick(cFile, Nc, pType, cType, oType, IP, LB, UB, mask)
%   Fits Crick parameters given the structure in cFile. If you use this
%   program in your research, please cite G. Grigoryan, W. F. DeGrado,
%   "Probing Designability via a Generalized Model of Helical Bundle
%   Geometry", J. Mol. Biol., 405(4): 1079-1100 (2011). See also
%   http://grigoryanlab.org/cccp/.
%
%   === Outputs:
%   err -- RMSD between ideal and input structure
%   xyz -- an N-by-3 matrix of coordinates of the ideal structure
%   p   -- a vector of structs representing best-fit Crick parameters.
%   Field 'name' contains the name of each parameter (some are derivative
%   parameters and some are primary ones). Field 'val' contains parameter
%   values. Fields 'UB' and 'LB' show upper and lower bounds (usually not
%   used) and field 'pri' indicates which parameters are primary.
%
%   === Inputs:
%   If cType is 1, assumes that cFile is a PDB file and uses the unix tools
%   grep and gawk to extract coordinates of CA atoms. NOTE: this will only
%   work on unix/linux, and this assumes that the PDB file is formatted
%   correctly and has only one CA atom per residue (e.g. no alternative
%   locations for CA atom). If cType is set to 2, cFile is assumed to be a
%   flat text file with coordinates of CA atoms (one row per atom, space
%   separated). This is a universal option and will work on any OS.
%   Nc is the number of chains in the file. The number of coordinates has to
%   be divisible by Nc.
%   pType is the parameterization type. All fits include superhelical
%   radius (R0), helical radius (R1), superhelical frequency (w0), helical
%   frequency (w1), and pitch angle (alpha) as global parameters. pType
%   then controls which of the remaining parameters are local to each
%   chain and which are global. Choices are:
%   'GENERAL' -- this is the most general option and means that each helix
%   gets its own helical phase (ph1), superhelical phase ffset (dph0), and
%   Z offset (Zoff), regarless of orientation.
%   'SYMMETRIC' -- the most symmetric option. ph1 is still individual to
%   each chain, but there is only one Zoff denoting the offset between all
%   parallel and all antiparallel chains (i.e. all chains sharing direction
%   are assumed to have zero offset with respect to each other). Also, for''
%   parallel topologies, dph0 is fixed at i*2*pi/Nc for the i-th chain (so
%   dph0 is not a parameter). For strictly alternating topologies
%   (up, down, up, down, ...), one dph0 parameter is allowed to designate
%   the relative rotation of all parallel chains with respect to all
%   anti-parallel chains. NOTE: this reqiures that chains alternate
%   direction in the order listed in the input file.
%   'ZOFF-SYMM' -- same as GENERAL, but Zoff is treated as in SYMMETRIC.
%   'DPH0-SYMM' -- same as GENERAL, but dph0 is treated as in SYMMETRIC.
%   'GENERAL-HLXPH' -- same as GENERAL, but all chains share the same ph1.
%   'SYMMETRIC-HLXPH' -- same as SYMMETRIC, but all chains share ph1.
%   'ZOFF-SYMM-HLXPH' -- same as ZOFF-SYMM, but all chains share ph1.
%   'DPH0-SYMM-HLXPH' -- same as DPH0-SYMM, but all chains share ph1.
%
%   oType determines the type of output. If a string, final parameters will
%   be written out to a file with that name and extension ".par". If 1,
%   will display parameters on the screen and also plot a figure to show
%   the difference between ideal and real coordinates. If 0.5, will only
%   display parameters without the plot.
%   IP -- vector with initial parameters for the search. Must be either
%   empty (i.e. []), in which case initialization is done automatically
%   (recommended) or have six entries, signifying starting values for
%   parameters is R0, R1, w0, w1, alpha, ph1 in that order (one value for
%   each, even if various non-symmetric options are specified).
%   LB, UB -- lower and upper bounds for parameters. Can be left empty
%   (recommended) or must be of the same length as IP.
%   mask -- if not empty, removes the contribution of certain atoms to the
%   total error (and hence the overall fit). Must have the same number of
%   elements as the total number of residues in the structure. Elements
%   should be either 0 (indicating this residues does not contribute) or 1
%   (indicating that the residue does contribute). The order of residues
%   starts with the first residue of the first chain and goes on to the
%   last residue of the last chain, as listed in cFIle. So, for example, if
%   we have a trimer with 28 residues in each, and we'd like to remove the
%   contribution of the 6 C-terminal residues of the second chain, imask
%   would be set to [ones(1, 50) zeros(1, 6) ones(1, 28)].
%

function [err XYZ pret] = fcrick(file, chains, parType, coorType, outType, iIP, iLB, iUB, imask)
global M x0 p0 sym co cr cparType shw extr mask;

% check input params
parType = upper(parType);
if (strcmp(parType, 'GENERAL') + strcmp(parType, 'SYMMETRIC') + strcmp(parType, 'ZOFF-SYMM') + strcmp(parType, 'DPH0-SYMM') + ...
    strcmp(parType, 'GENERAL-HLXPH') + strcmp(parType, 'SYMMETRIC-HLXPH') + strcmp(parType, 'ZOFF-SYMM-HLXPH') + strcmp(parType, 'DPH0-SYMM-HLXPH')== 0)
    error('Unknown parameterization type "%s"', parType);
end
if (isempty(iIP))
    IP = [5.0 2.26 -2*pi()/100 4*pi()/7 -12*pi()/180 0];
elseif (size(iIP, 1)*size(iIP, 2) ~= 6)
    error('Unexpected number of initial parameters!');
else
    IP = iIP;
end

if (coorType == 1)
    tmpxyzf = '.tmp.xyz';
    [status, result] = system(sprintf('grep " CA " %s | gawk ''{print $7"\t"$8"\t"$9}'' > %s', file, tmpxyzf));
    if (~strcmp(result, ''))
        error(sprintf('Error extracting CA atoms (%s)', result));
    end
    M = load(tmpxyzf);
else
    M = load(file);
end

n = size(M, 1);
if (mod(n, chains) ~= 0) error(sprintf('number of coordinates (%d) is not divisible by the number of chains (%d)', n, chains)); end
sym = struct('olig', chains, 'type', 'GENERAL');
shw = 0;
mask = imask;

% determine chain order (clock-wise)
if (sym.olig > 2)
    nc = n/sym.olig;
    ind = setdiff(1:nc, find(mask==0));
    c = M(ind(round(length(ind)/2)), :); % coordinates of the middle CA in each layer
    pz = M(ind(end), :) - M(ind(1), :); % the first chain defines the positive axis direction
    for i = 2:sym.olig
        ind = setdiff((i-1)*nc+1:i*nc, find(mask==0));
        [m mi] = min(sum((M(ind, :) - repmat(c(1, :), length(ind), 1)).^2, 2));
        c(end+1, :) = M(ind(1)+mi-1, :);
    end
    cc = mean(c); % center of the bundle
    an = [0]; % rotation angles of all chains
    for i = 2:sym.olig
        an(i) = acos(dot(cc - c(1, :), cc - c(i, :))/norm(c(1, :) - cc)/norm(cc - c(i, :)));
        if (dot(cross(cc - c(1, :), cc - c(i, :)), pz) < 0)
            an(i) = 2*pi - an(i);
        end
    end
    [an co] = sort(an);
else
    % for dimers (and of course monomers) chain order does not matter
    co = 1:sym.olig;
end

% determine chain orientation and handedness (sign of crossing angle)
if (sym.olig >= 2)
    nc = n/sym.olig;
    % the first chain defines the positive axis direction
    ind = setdiff(1:nc, find(mask==0));
    pz = M(ind(end), :) - M(ind(1), :); cr = 1;
    crsg = zeros(sym.olig-1, 1);
    for i = 2:sym.olig
        ind = setdiff((i-1)*nc+1:i*nc, find(mask==0));
        cr(i) = (dot(M(ind(end), :) - M(ind(1), :), pz) > 0);
        % determine the sign of crossing angle as the dihedral angle
        % between vectors representing the first heptad on each chain
        crsg(i-1) = crossingAngle(M(setdiff(1:nc, find(mask==0)), :), M(setdiff((i-1)*nc+1:i*nc, find(mask==0)), :), cr(i));
    end
    % dial in the right handedness, unless initial params were set
    if (isempty(iIP))
        hand = sign(mean(crsg));
        if (hand == 0) hand = 1; end % resolve in case the average comes out to be zero for some reason
        IP(3) = abs(IP(3))*sign(hand);
%        IP(5) = mean(crsg)/1.3; % guess pitch angle as function of crossing angle
        IP(5) = abs(IP(5))*sign(hand);
%        IP(3) = 1.51*sin(IP(5))/IP(1);
    end
else
    % for monomers chain orientation does not matter
    cr = 1;
end
ap = find(cr == 0);
nap = length(ap); % number of chains anti-parallel to the first

%%%%%%%% Minimization

% initial parameters
p0 = struct('val', [], 'name', {}, 'LB', [], 'UB', []); 
p0 = addparam(p0, IP(1), 'R0 (A)', 0, 30, 'pri', 1); % super-helical radius
p0 = addparam(p0, IP(2), 'R1 (A)', 2, 3, 'pri', 1); % alpha-helical radius
p0 = addparam(p0, IP(3), 'w0 (rad/res)', -10*pi()/100, 10*pi()/100, 'pri', 1); % uper-helical frequency
p0 = addparam(p0, IP(4), 'w1 (rad/res)', 1*pi()/100, 8*pi()/100, 'pri', 1); % alpha-helical frequency
p0 = addparam(p0, IP(5), 'alpha (rad)', -pi()/6, pi()/6, 'pri', 1); % crossing angle
% alpha-helical phase
if (isempty(strfind(parType, 'HLXPH')))
    p0 = addparam(p0, IP(6), 'ph1 (rad)', -pi(), pi(), 'pri', 1);
else
    for i = 1:sym.olig
        p0 = addparam(p0, IP(6), sprintf('ph1 for chain %d (rad)', i), -pi(), pi(), 'pri', 1);
    end
end
% z-shifts
if (~isempty(strfind(parType, 'SYMMETRIC')) || ~isempty(strfind(parType, 'ZOFF-SYMM')))
    % one z-offset (all parallel vs. all anti-parallel)
    if (~isempty(find(cr == 0)))
        p0 = addparam(p0, max(M(:, 3))-min(M(:, 3)), 'absolute ap zoff (A)', -max(M(:, 3))+min(M(:, 3)), max(M(:, 3))-min(M(:, 3)), 'pri', 4);
    end
else
    % all chains have a z-offset (except the first, of course)
    for i = 2:sym.olig
        if (cr(i))
            p0 = addparam(p0, 0, sprintf('absolute zoff_%d (A)', i), -max(M(:, 3))+min(M(:, 3)), max(M(:, 3))-min(M(:, 3)), 'pri', 4);
        else
            p0 = addparam(p0, max(M(:, 3))-min(M(:, 3)), sprintf('absolute zoff_%d (A)', i), -max(M(:, 3))+min(M(:, 3)), max(M(:, 3))-min(M(:, 3)), 'pri', 4);
        end
    end
end
% super-helical phase offsets
if (isempty(strfind(parType, 'SYMMETRIC')) && isempty(strfind(parType, 'DPH0-SYMM')))
    for i = 2:sym.olig
        p0 = addparam(p0, -(co(i)-1)*2*pi()/sym.olig, sprintf('dph0_%d (rad)', i), -2*pi(), 0, 'pri', 1);
    end
elseif (~isempty(strfind(parType, 'SYMMETRIC')) && (all(cr(co(1:2:end)) == 1) && all(cr(co(2:2:end)) == 0)))
    % for symmetric cases with alternating up/down orientations, assume Dn symmetry, so
    % there is only one superhelical phase offset - all parallels relative to all antiparallels
	p0 = addparam(p0, -2*pi()/sym.olig, 'dph0_p_ap (rad)', -2*pi(), 0, 'pri', 1);
end

progV = ver;
% optimization settings
%opts = optimset('MaxFunEvals', 30000, 'MaxIter', 30000, 'Jacobian','on', 'DerivativeCheck', 'on');
%opts = optimset('MaxFunEvals', 30000, 'MaxIter', 30000, 'Jacobian','on');
%x = lsqnonlin(@crickssd, x0, [4 1.5 -2*pi()/50 2*pi()/7 0 0 0 -pi -pi -pi -100 -100 -100], [10 3 -2*pi()/400 6*pi()/7 pi()/6 pi() pi() pi pi pi 100 100 100], opts);
opts = optimset('MaxFunEvals', 30000, 'MaxIter', 30000);

% alternate between minimizing over orientation only or over everything
err = inf; x0 = [p0.val];
for it =1:100 % do at most 100 iterations
    % minimize over each parameter separately
    cparType = sprintf('single%s', parType);
%    ordi = setdiff(randperm(length(x0)), [mystrfind('R1', {p0.name}) mystrfind('w1', {p0.name})]);
%    ordi = randperm(length(x0));
%    ordi = [mystrfind('zoff', {p0.name})' setdiff(randperm(length(x0)), [mystrfind('R1', {p0.name})' mystrfind('w1', {p0.name})'])];
%    ordi = [mystrfind('zoff', {p0.name})' setdiff(1:length(x0), [mystrfind('R1', {p0.name})' mystrfind('w1', {p0.name})'])];
    ordi = mymat2cell([mystrfind('zoff', {p0.name}) mystrfind('dph0', {p0.name}) mystrfind('ph1', {p0.name})]);
    ordi = [ordi mymat2cell(setdiff(1:length(x0), [[ordi{:}] mystrfind('R1', {p0.name})' mystrfind('w1', {p0.name})']))];
    for i = 1:length(ordi)
        extr.vary = ordi{i};
        x1 = fminunc(@crickssd, x0(ordi{i}), opts);
%        x1 = fminbnd(@crickssd, p0(ordi{i}).LB, p0(ordi{i}).UB, opts);
%        x1 = fmincon(@crickssd, x0(i), [], [], [], [], p0(i).LB, p0(i).UB, [], opts);
        x0(ordi{i}) = x1;
    end

    % minimize over everything together
    cparType = parType;
    if (isempty(iLB))
        if (~isempty(mystrfind('MATLAB', {progV.Name})))
            x = lsqnonlin(@crickssd, x0);
        else
            x = fminunc(@crickssd, x0);
        end
        x = fminunc(@crickssd, x, opts);
%        x = fminunc(@crickssd, x0, opts);
    else
        x = fmincon(@crickssd, x0, [], [], [], [], iLB, iUB, [], opts);
    end
    x0 = x;

    % has the error dicreased significantly?
    errn = crickssd(x);
    disp(errn);
    if (err - errn < 0.00001) break; end
    err = errn;
end
for i = 1:length(x)
    p0(i).val = x(i);
end

% remove the effect of 2*pi additives
p0(mystrfind('w0', {p0.name})).val = angle_pmp(p0(mystrfind('w0', {p0.name})).val);
p0(mystrfind('w1', {p0.name})).val = angle_pmp(p0(mystrfind('w1', {p0.name})).val);
p0(mystrfind('alpha', {p0.name})).val = angle_pmp(p0(mystrfind('alpha', {p0.name})).val);
for i = mystrfind('ph1', {p0.name})
    p0(i).val = angle_pmp(p0(i).val);
end
for i = mystrfind('dph0', {p0.name})
    p0(i).val = angle_pmp(p0(i).val);
end
% alpha has to be the same sign as w0 (the objective function adjusts inside, but the answer can come out having different signs)
p0(mystrfind('alpha', {p0.name})).val = abs(p0(mystrfind('alpha', {p0.name})).val)*sign(p0(mystrfind('w0', {p0.name})).val);

% -- if R1 comes out negative, adjust the helical phase to make it positive
r1ind = mystrfind('R1', {p0.name});
if (p0(r1ind).val < 0)
    ph1ind = mystrfind('ph1', {p0.name});
    for i = 1:length(ph1ind)
        p0(ph1ind(i)).val = p0(ph1ind(i)).val + pi();
    end
    p0(r1ind).val = -p0(r1ind).val;
end

% -- display results/compute final error and coordinates
shw = 0; [err XYZi] = crickssd(x); % coordinates without global transformations
if ((~ischar(outType)) && (outType == 1)) shw = 1;
else shw = 0.5; end
[err XYZ] = crickssd(x);

% -- for anti-parallel chains, compute Z-offset between a and the closest d
% first, find the starting heptad position
ph1ind = mystrfind('ph1', {p0.name});
if (length(ph1ind) > 1)
    for i = 1:length(ph1ind)
        ph1(i) = p0(ph1ind(i)).val;
        fhp(i) = getHeptadPosition(ph1(i));
        p0 = addparam(p0, fhp(i), sprintf('starting heptad position for chain %d', i), fhp(i), fhp(i), 'pri', 1);
    end
else
    ph1 = p0(ph1ind).val*ones(1, sym.olig);
    fhp = char(getHeptadPosition(ph1(1))*ones(1, sym.olig));
    p0 = addparam(p0, fhp(1), 'starting heptad position', fhp(1), fhp(1), 'pri', 1);
end
% do heptad positions persist for at least one heptad (depends on w1)?
if (abs(p0(mystrfind('w1', {p0.name})).val - 4*pi()/7) > pi()/7)
    p0 = addparam(p0, 'The fit structure does not appear to be a canonical 7-residue repeat coiled coil. Be cautious when interpreting heptad based parameters!', 'message', '', '');
end

% extract useful values
R0 = p0(mystrfind('R0', {p0.name})).val;
R1 = p0(mystrfind('R1', {p0.name})).val;
w0 = p0(mystrfind('w0', {p0.name})).val;
w1 = p0(mystrfind('w1', {p0.name})).val;
a = p0(mystrfind('alpha', {p0.name})).val;

% compute some additional commonly-reported values
p0 = addparam(p0, abs(2*pi()*R0/tan(a)), 'pitch (A)', abs(2*pi()*R0/tan(a)), abs(2*pi()*R0/tan(a)), 'pri', 1); % pitch (defined always to be positive)
p0 = addparam(p0, R0*w0/sin(a), 'rise per residue (A)', R0*w0/sin(a), R0*w0/sin(a), 'pri', 1); % pitch

% fetch absolute z offsets
zoff = zeros(sym.olig, 1);
if (~isempty(strfind(parType, 'SYMMETRIC')) || ~isempty(strfind(parType, 'ZOFF-SYMM')))
    if (~isempty(find(cr == 0)))
        zo = [p0(mystrfind('absolute ap zoff', {p0.name})).val];
        zoff(find(cr == 0)) = zo;
    end
else
    zoff = [0 [p0(mystrfind('absolute zoff_', {p0.name})).val]];
end

% Compute corrected Z-offsets
% Find the Z-offset between points on the helices that have a helical phase
% of pi (e.g. point directly into the superhelix). Values of parameter t
% that correspond to such points are t = a*n + b, with n an integer, a =
% 2*pi/w1, and b = (pi - ph1). The slope, a, is the same for all chains,
% while intercept b can be different, if helical phase offsets are varied
% differently.
aa1 = 2*pi/w1; b1 = (pi - ph1(1))/w1; z1 = (R0*w0/tan(a))*(aa1*1 + b1);
if (~((~isempty(strfind(parType, 'SYMMETRIC')) || ~isempty(strfind(parType, 'ZOFF-SYMM'))) && isempty(find(cr == 0))))
    for ci = 2:sym.olig
%        dz = absoluteToRegisterZoff(zoff(ci), R0, w0, a, w1, ph1(1), ph1(ci), cr(ci));
%        p0 = addparam(p0, dz, sprintf('register Z-offset for chain %d (A)', ci), dz, dz, 'pri', 1);
        dz = absoluteToZoff_aa(zoff(ci), R0, w0, a, R1, w1, ph1(1), ph1(ci), cr(ci));
        p0 = addparam(p0, dz, sprintf('Z_aa for chain %d (A)', ci), dz, dz, 'pri', 1);
    end
end

% correct absolute Z-offsets for AP chains (N to C distance, not N to N distance)
for ci = 2:sym.olig
    if (cr(ci) == 1) continue; end
    if (~isempty(strfind(parType, 'SYMMETRIC')) || ~isempty(strfind(parType, 'ZOFF-SYMM')))
        zoffi = mystrfind('absolute ap zoff', {p0.name});
%        p0(zoffi).val = XYZi((ci-1)*nc + 1, 3) - XYZi(1,3);
        p0(zoffi).val = XYZi(ci*nc, 3) - XYZi(1,3);
        break;
    else
        zoffi = mystrfind('absolute zoff_', {p0.name});
%        p0(zoffi(ci-1)).val = XYZi((ci-1)*nc + 1, 3) - XYZi(1,3);
        p0(zoffi(ci-1)).val = XYZi(ci*nc, 3) - XYZi(1,3);
    end
end

% create a structure summary
p0 = addparam(p0, sprintf('%d (fit %d residues in each)', sym.olig, nc), 'structure summary: chains', [], [], 'pri', 0);
for ci = 1:sym.olig
    ss = sprintf('%.3f Angstrom in Z-direction', abs(XYZi((ci-1)*nc + 1, 3) - XYZi(ci*nc, 3)));
    if (ci > 1)
        if (cr(ci)) ss = sprintf('%s; parallel to chain 1', ss);
        else ss = sprintf('%s; anti-parallel to chain 1', ss); end
    end
    p0 = addparam(p0, ss, sprintf('structure summary: chain %d', ci), [], [], 'pri', 0);
end
pret = p0;

% -- print results
for i = 1:length(p0)
    if (ischar(p0(i).val)) disp(sprintf('%s = %s', p0(i).name, p0(i).val));
    else disp(sprintf('%s = %f', p0(i).name, p0(i).val)); end
end
if (ischar(outType))
    p0(end+1).name = 'error'; p0(end).val = err;
    saveParms(p0, sprintf('%s.par', outType));
end

%%%%%%%%
function p = addparam(p, val, name, lb, ub, varargin)

n = length(p);
p(n+1).val = val;
p(n+1).name = name;
if (isempty(lb))
    p(n+1).LB = val;
else
    p(n+1).LB = lb;
end
if (isempty(ub))
    p(n+1).UB = val;
else
    p(n+1).UB = ub;
end
if (nargin > 5)
    if (mod(length(varargin), 2) ~= 0) error('Expected an even number of ''property'', ''value'' pairs!'); end
    eval(sprintf('p(n+1).%s = %f;', varargin{1}, varargin{2}));
end


%%%%%%%%
function p = saveParms(p, file)

fp = fopen(file, 'w');
if (fp < 0)
    error('could not open file %s', file);
end
% order parameters by priority
pri = sort(unique([p.pri]));
for ii = 1:length(pri)+1
    for i = 1:length(p)
        % first priorities by order, and then those without a priority
        % (some elements can have no priority so be careful in test)
        if (((ii <= length(pri)) && (sum([(p(i).pri == pri(ii)) 1]) ~= 2)) || ((ii == length(pri)+1) && (~isempty(p(i).pri))))
            continue;
        end
        if (ischar(p(i).val))
            fprintf(fp, '%s = %s\n', p(i).name, p(i).val);
        else
            fprintf(fp, '%s = %f\n', p(i).name, p(i).val);
        end
    end
end
fclose(fp);

%%%%%%%%
function hp = getHeptadPosition(ph1)
ph1 = ph1 + pi()/7;
in = find(ph1 < 0);
if (~isempty(in)) ph1(in) = ph1(in) + 2*pi()*ceil(abs(ph1(in)/2/pi())); end

in = find(ph1 > 2*pi());
if (~isempty(in)) ph1(in) = ph1(in) - 2*pi()*floor(abs(ph1(in)/2/pi())); end
    
fhpi = floor(7*ph1/(2*pi())) + 1;
if (any(any(fhpi <= 0) + any(fhpi > 7)))
    disp('ph1 is:'); disp(ph1);
end
hps = 'fcgdaeb'; hp = hps(fhpi);

%%%%%%%%
% computes all pairwise distances between coordinates in the first argument
% and the second argument
function p = myvdist(X1, X2)

p = [];
for i = 1:size(X1, 1)
    for j = 1:size(X2, 1)
        p(end+1) = sqrt(sum((X1(i, :) - X2(j, :)).^2));
    end
end

%%%%%%%%
% looks for matches of string in first argument in all of elements of cell
% array in the second argument
function m = mystrmatch(query, array)

m = find(~cellfun(@isempty, strfind(array, query)));

%%%%%%%%
% puts an angle between +pi and -pi
function a = angle_pmp(a)
a = mod(a, 2*pi);
i = find(a > pi);
a(i) = a(i) - 2*pi;

%%%%%%%%
% Computes the crossing angle between two chains of any size. First looks
% for a pair of segment of 8 residues, one on each chain, such that they
% are closest to each other. Then computes the crossing angle for these two
% segments
function a = crossingAngleold(A, B, pap)
if ((size(A, 2) ~= 3) || (size(B, 2) ~= 3))
    error('Unexpected matrix size: A is [%d x %d], B is [%d x %d]!', size(A), size(B));
end

if (pap == 0)
    B = B(end:-1:1, :);
end

dnc = min([7 size(A, 1)-1 size(B, 1)-1]);
b1 = 1; b2 = 1;
bd = sqrt(mean(sum((A(b1:b1+dnc, :) - B(b2:b2+dnc, :)).^2, 2)));
for i = 1:size(A, 1)-dnc
    for j = 1:size(B, 1)-dnc
        if (bd > sqrt(mean(sum((A(i:i+dnc, :) - B(j:j+dnc, :)).^2, 2))))
            b1 = i; b2 = j;
            bd = sqrt(mean(sum((A(i:i+dnc, :) - B(j:j+dnc, :)).^2, 2)));
        end
    end
end

a = -dihe(A(b1, :), A(b1+dnc, :), B(b2+dnc, :), B(b2, :));

%%%%%%%%
% Computes the crossing angle between two helical chains by considering the
% helical axis
function a = crossingAngle(A, B, pap)
if ((size(A, 2) ~= 3) || (size(B, 2) ~= 3))
    error('Unexpected matrix size: A is [%d x %d], B is [%d x %d]!', size(A), size(B));
end

if (pap == 0)
    B = B(end:-1:1, :);
end

if (min(size(A, 1), size(B, 1)) < 3)
    a = 0; return;
end

% find helical axes
axsA = helicalAxisPoints(A);
axsB = helicalAxisPoints(B);

a = -dihe(axsA(1, :), axsA(end, :), axsB(end, :), axsB(1, :));

%%%%%%%%
% Computes the crossing angle between two helical chains by considering the
% helical axis
function axs = helicalAxisPoints(H)

axs = [];
if (size(H, 1) < 3)
    return;
end

for i = 2:size(H, 1)-1
    r = (H(i-1, :) - H(i, :)) + (H(i+1, :) - H(i, :));
    r = 2.26*r/norm(r);
    axs(i-1, :) = 2.26*r/norm(r) + H(i, :);
end

%%%%%%%%
% Converts a given 2D matrix into a cell matrix
function c = mymat2cell(m)
c = mat2cell(m, ones(size(m, 1), 1), ones(size(m, 2), 1));


%%%%%%%%
% absoluteToZoff_aa Converts from an absolute Z offset to a Zaa' Z offset,
% defined as the Z distance between closest 'a' positions on opposite chains
function zaa = absoluteToZoff_aa(zoff, R0, w0, a, R1, w1, ph1_1, ph1_2, p_ap)

assert((p_ap == 1) || (p_ap == 0), 'p_ap flag expected to be either 1 or -1');

% find first a-positions on chain 1
rng = 0:6;
[m mi] = min(abs(angleDiff(ph1_1 + w1*rng, canonicalPhases(1)*pi/180))); % the position closest to the canonical 'a'-position phase in the first heptad plus a little
aph1_1 = mod(ph1_1 + w1*rng(mi), 2*pi);
az1 = w0*rng(mi)*R0/tan(a) - R1*sin(a)*sin(w1*rng(mi) + ph1_1);

% * keep going through 'a' positions on chain 2, looking for smallest distance with
% the first a-position on chain 1, until the sign of the distance switches
% * start with the residue on chain 2 that is "close-ish" to the first 'a'
% on chain 2
if (p_ap == 0)
    n = round((zoff - az1)*tan(a)/w0/R0);
else
    n = round((az1 - zoff)*tan(a)/w0/R0);
end
sgn = nan; zaa = inf; f = 0;
for count = 0:99 % if the loop does not happen within a couple of iterations, then we have a very bad fit and it does not matter anyway
    % try up and down
    for ni = [n-count, n+count]
        aph1_2 = ph1_2 + w1*ni;
        if (getHeptadPos(aph1_2, 1) ~= 1) % though phase changes sign, to determine whether something is an 'a' or not, we still need the original phase
            continue;
        end
        if (p_ap == 0)
            az2 = zoff - (w0*ni*R0/tan(a) - R1*sin(a)*sin(w1*ni + ph1_2));
        else
            az2 = zoff + w0*ni*R0/tan(a) - R1*sin(a)*sin(w1*ni + ph1_2);
        end
        if (abs(zaa) > abs(az2 - az1))
            zaa = az2 - az1;
        end
        if (isnan(sgn) && (sign(az2 - az1) ~= 0))
            sgn = sign(az2 - az1);
        end
        if (sgn*sign(az2 - az1) <= 0) % make sure both positive and negative Zaa' are tried, or if zero is found, that's obviously the lowest possible value
            f = 1; break;
        end
    end
    if (f) break; end
end


%%%%%%%%
% Returns the error and its gradient between the coordinate set and the
% ideal backbone given Crick parameters
function [ssd, J] = crickssd(p)
global M x0 sym co cr cparType shw extr mask;

if (strfind(upper(cparType), 'SINGLE'))
    ptemp = p; p = x0; p(extr.vary) = ptemp;
end

if (strfind(upper(cparType), 'GENERAL'))
    [r0 p] = shift(p);
    [r1 p] = shift(p);
    [w0 p] = shift(p);
    [w1 p] = shift(p);
    [a p] = shift(p);
    if (strfind(upper(cparType), 'HLXPH'))
        [ph1 p] = shift(p, sym.olig);
    else
        [ph1 p] = shift(p);
    end
    [zoff p] = shift(p, sym.olig - 1);
    [dph0 p] = shift(p, sym.olig - 1);
end

if (strfind(upper(cparType), 'SYMMETRIC'))
    [r0 p] = shift(p);
    [r1 p] = shift(p);
    [w0 p] = shift(p);
    [w1 p] = shift(p);
    [a p] = shift(p);
    if (strfind(upper(cparType), 'HLXPH'))
        [ph1 p] = shift(p, sym.olig);
    else
        [ph1 p] = shift(p);
    end
    zoff = zeros(1, sym.olig);
    ap = find(cr == 0);
    if (~isempty(ap)) [zoff(ap) p] = shift(p); end
    if (all(cr(co(1:2:end)) == 1) && all(cr(co(2:2:end)) == 0))
        [dph0_pap p] = shift(p);
        dph0 = zeros(sym.olig, 1);
        par_ind = find(cr); [tmp si] = sort(co(par_ind)); [tmp si] = sort(si); % clockwise order of parallel chains
        dph0(par_ind) = -(si-1)*2*pi()/length(par_ind); % the first chain (the one listed first in chain order) is always "parallel"
        apar_ind = find(cr == 0); [tmp si] = sort(co(apar_ind)); [tmp si] = sort(si); % clockwise order of anti-parallel chains
        dph0(apar_ind) = dph0_pap - (si-1)*2*pi()/length(apar_ind);
    else
        dph0 = -(co-1)*2*pi()/sym.olig;
    end
end

if (strfind(upper(cparType), 'ZOFF-SYMM'))
    [r0 p] = shift(p);
    [r1 p] = shift(p);
    [w0 p] = shift(p);
    [w1 p] = shift(p);
    [a p] = shift(p);
    if (strfind(upper(cparType), 'HLXPH'))
        [ph1 p] = shift(p, sym.olig);
    else
        [ph1 p] = shift(p);
    end
    zoff = zeros(1, sym.olig);
    ap = find(cr == 0);
    if (~isempty(ap)) [zoff(ap) p] = shift(p); end
    [dph0 p] = shift(p, sym.olig-1);
end

if (strfind(upper(cparType), 'DPH0-SYMM'))
    [r0 p] = shift(p);
    [r1 p] = shift(p);
    [w0 p] = shift(p);
    [w1 p] = shift(p);
    [a p] = shift(p);
    if (strfind(upper(cparType), 'HLXPH'))
        [ph1 p] = shift(p, sym.olig);
    else
        [ph1 p] = shift(p);
    end
    [zoff p] = shift(p, sym.olig-1);
    dph0 = -(co-1)*2*pi()/sym.olig;
end
a = abs(a).*sign(w0); % make sure pitch angle and frequency have the same sign

switch upper(sym.type)
    case {'NONE'}
        n = size(M, 1);
        t = 0:n-1;
        [x y z] = crickEQ(r0, r1, w0, w1, a, dph0, ph1, t);
    case {'C2'}
        n = size(M, 1)/2;
        t = 0:n-1;
        [x1 y1 z1] = crickEQ(r0, r1, w0, w1, a, 0, ph1, t);
        [x2 y2 z2] = crickEQ(r0, r1, w0, w1, a, 0 + pi, ph1, t);
        x = [x1 x2]; y = [y1 y2]; z = [z1 z2];
    case {'C3'}
        n = size(M, 1)/3;
        t = 0:n-1;
        [x1 y1 z1] = crickEQ(r0, r1, w0, w1, a, 0, ph1, t);
        [x2 y2 z2] = crickEQ(r0, r1, w0, w1, a, 0 + 2*pi/3, ph1, t);
        [x3 y3 z3] = crickEQ(r0, r1, w0, w1, a, 0 + 4*pi/3, ph1, t);
        x = [x1 x2 x3]; y = [y1 y2 y3]; z = [z1 z2 z3];
    case {'C'}
        n = size(M, 1); nc = n/sym.olig;
        x = zeros(1, n); y = x; z = x;
        t = 0:nc-1;
        for i = 1:sym.olig
            [x1 y1 z1] = crickEQ(r0, r1, w0, w1, a, 0 + 2*pi()*(i-1)/sym.olig, ph1, t);
            x((co(i)-1)*nc+1:co(i)*nc) = x1;
            y((co(i)-1)*nc+1:co(i)*nc) = y1;
            z((co(i)-1)*nc+1:co(i)*nc) = z1;
        end
    case {'GENERAL.old'}
        n = size(M, 1); nc = n/sym.olig;
        x = zeros(1, n); y = x; z = x;
        t = -floor(nc/2):floor(nc/2) + mod(nc, 2) - 1;
        for i = 1:sym.olig
            if (cr(i) == 0) lt = t(end:-1:1);
            else lt = t; end
            [x1 y1 z1] = crickEQ(r0, r1, w0, w1, a, 0 + dph0(i), ph1, lt);
            x((co(i)-1)*nc+1:co(i)*nc) = x1;
            y((co(i)-1)*nc+1:co(i)*nc) = y1;
            z((co(i)-1)*nc+1:co(i)*nc) = z1 + zoff(i);
        end
    case {'GENERAL'}
        n = size(M, 1); nc = n/sym.olig;
        XYZ = generateCrickBB(sym.olig, nc, r0, r1, w0, w1, a, ph1, cr, dph0, zoff);
end

% sometimes, some very weird parameters are passed (like inf), causing the
% coordinates not to make sense. In this case, set RMSD to infinity.
if (any(isnan(XYZ(:))))
    ssd = 10^5; mat = eye(3,3);
    return;
end

%ssd = sqrt(sum((M - [x' y' z']).^2, 2)); % for least squares
if (isempty(mask))
%    ssd = sum(sum((M - [x' y' z']).^2)); % for others
    [ssd mat] = superimpose(M', XYZ');
else
%    ssd = sum(sum((M - [x' y' z']).^2, 2).*mask); % for others
%    [ssd mat] = superimpose(M'.*repmat(mask.^0.5, 1, 3)', XYZ'.*repmat(mask.^0.5, 1, 3)');
    [ssd mat] = superimpose(M(mask == 1, :)', XYZ(mask == 1, :)');
end
% constrain super-helical phase offset for dimers to be realistic
if ((sym.olig == 2) && (shw == 0)) % only during fitting
    dd = abs(-(co(2:sym.olig)-1)*2*pi()/sym.olig - dph0(end));
    % difference between ideal values for this olig state vs the next
    % highest (for dimers/trimers this is 60 degrees)
%    ind = find(dd > (2*pi/sym.olig() - 2*pi/(sym.olig()+1)));
%    ssd = ssd + sum(dd(ind));
    ssd = ssd + (dd/(pi/3)/50)*n*nc/42; % penalize by 0.02 RMSD for every pi/3, for every three heptads of a dimer
end
% constrain Ca-Ca distances to be realistic
%d = sqrt((x(1) - x(2)).^2 + (y(1) - y(2)).^2 + (z(1) - z(2)).^2);
%if ((d > 3.9) || (d < 3.6))
%    ssd = ssd*(1+100*(min(abs([d-3.9, d-3.6]))));
%end
if (isnan(ssd) || isinf(ssd)) ssd = 10^10; end


if ((shw == 0.5) || (shw == 1))
    if (isempty(mask))
        XYZ = (XYZ - repmat(mean(XYZ), size(XYZ, 1), 1))*mat + repmat(mean(M), size(M, 1), 1); % super-impose fit onto target
    else
        XYZ = (XYZ - repmat(mean(XYZ(mask == 1, :)), size(XYZ, 1), 1))*mat + repmat(mean(M(mask == 1, :)), size(M, 1), 1); % super-impose fit onto target
    end
    if (shw == 1)
        plot3(XYZ(:, 1), XYZ(:, 2), XYZ(:, 3), '-o'); hold on; plot3(M(:, 1), M(:, 2), M(:, 3), '-ro');
    end
end
if (nargout > 1)
    J = XYZ;
    return;
end

if nargout > 1
    dxd = [
        cos(w0*t+dph0);
        (cos(w0*t+dph0).*cos(w1*t+ph1) - cos(a)*sin(w0*t+dph0).*sin(w1*t+ph1));
        (-r0*sin(w0*t+dph0).*t - r1*sin(w0*t+dph0).*cos(w1*t+ph1).*t - r1*cos(a)*cos(w0*t+dph0).*sin(w1*t+ph1).*t);
        (-r1*cos(w0*t+dph0).*sin(w1*t+ph1).*t - r1*cos(a)*sin(w0*t+dph0).*cos(w1*t+ph1).*t);
        0*ones(1, nc);
        (-r1*cos(w0*t+dph0).*sin(w1*t+ph1) - r1*cos(a)*sin(w0*t+dph0).*cos(w1*t+ph1));
        (-r0*sin(w0*t+dph0) - r1*sin(w0*t+dph0).*cos(w1*t+ph1) - r1*cos(a)*cos(w0*t+dph0).*sin(w1*t+ph1))
    ]';
    dyd = [
        sin(w0*t+dph0);
        (sin(w0*t+dph0).*cos(w1*t+ph1) + cos(a)*cos(w0*t+dph0).*sin(w1*t+ph1));
        (r0*cos(w0*t+dph0).*t + r1*cos(w0*t+dph0).*cos(w1*t+ph1).*t - r1*cos(a)*sin(w0*t+dph0).*sin(w1*t+ph1).*t);
        (-r1*sin(w0*t+dph0).*sin(w1*t+ph1).*t + r1*cos(a)*cos(w0*t+dph0).*cos(w1*t+ph1).*t)
        0*ones(1, nc);
        (-r1*sin(w0*t+dph0).*sin(w1*t+ph1) + r1*cos(a)*cos(w0*t+dph0).*cos(w1*t+ph1))
         (r0*cos(w0*t+dph0) + r1*cos(w0*t+dph0).*cos(w1*t+ph1) - r1*cos(a)*sin(w0*t+dph0).*sin(w1*t+ph1))
    ]';
    dzd = [
        -w0*t/tan(a);
        (-sin(a)*sin(w1*t+ph1));
        (-t*r0/tan(a));
        (-r1*sin(a)*cos(w1*t+ph1).*t);
        w0*t*r0/((sin(a))^2);
        (-r1*sin(a)*cos(w1*t+ph1));
        0*ones(1, nc)
    ]';
    
    for i = 1:length(p)
%        J(:, i) = ((x' - M(:, 1)).*dxd(:, i) + (y' - M(:, 2)).*dyd(:, i) + (z' - M(:, 3)).*dzd(:, i))./ssd; % for least squares
        J(i) = sum(2*((x' - M(:, 1)).*dxd(:, i) + (y' - M(:, 2)).*dyd(:, i) + (z' - M(:, 3)).*dzd(:, i))); % for others
    end
end


%%%%%%%%
% DIHE calculates the dihedral angle given four points
%   DIHE(p1, p2, p3, p4) returns the dihedral angle p1p2p3p4.
%   Author: Gevorg Grigoryan
function d = dihe(p1, p2, p3, p4)

dim = size(p1, 1);
v12 = p1 - p2;
v23 = p2 - p3;
v43 = p4 - p3;

px1 = cross(v12, v23);
px1 = px1./(ones(dim, 1)*sqrt(sum(px1.^2)));

px2 = cross(v43, v23);
px2 = px2./(ones(dim, 1)*sqrt(sum(px2.^2)));

dp12 = dot(px1, px2);
sin2 = 1 - dp12.^2;

d = pi()/2.0 - atan(dp12./sqrt(sin2));

px3 = cross(px1, px2);
ind = find(dot(px3, v23) > 0);
if (~isempty(ind)) d(ind) = -d(ind); end


%%%%%%%%
function d = angleDiff (a, b)

d = mod((mod(a, 2*pi) - mod(b, 2*pi)), 2*pi);
in = find(d > pi);
d(in) = d(in) - 2*pi;


%%%%%%%%
function ph = canonicalPhases(ind)

% canonical phases are:
% [41.0, 95.0, 146.0, 197.0, 249.0, 300.0, 351.0]
% corresponding to {'c', 'g', 'd', 'a', 'e', 'b', 'f'}, respectively

% in the order a-g:
median_phases = [197.0, 300.0, 41.0, 146.0, 249.0, 351.0 95.0];
ph = median_phases(ind);

function hp = getHeptadPos(ph1, varargin)
%getHeptadPos(ph1) Returns the heptad position corresponding to phase ph1.
%getHeptadPos(ph1, 1) Same, but returns an int from 1 to 7, corresponding to
% 'a' through 'g', instead of characters.

meds = canonicalPhases(1:7)*pi/180;
hps = {'a', 'b', 'c', 'd', 'e', 'f', 'g'};

% sort phases in the order they appear on the helical wheel
[meds si] = sort(meds);
hps = hps(si);

ph1 = mod(ph1, 2*pi);
for i = 1:length(meds)
    pin = i - 1;
    if (i == 1) pin = length(meds); end
    nin = i + 1;
    if (i == length(meds)) nin = 1; end

    lb = mod(angleDiff(meds(pin), meds(i))/2 + meds(i), 2*pi);
    ub = mod(angleDiff(meds(nin), meds(i))/2 + meds(i), 2*pi);
    if ((angleDiff(ph1, lb) > 0) && (angleDiff(ub, ph1) > 0))
        if ((nargin > 1) && varargin{1})
            hp = hps{i} - 'a' + 1;
        else
            hp = hps{i};
        end
        break;
    end
end


function [x, y, z] = crickEQ(r0, r1, w0, w1, a, ph0, ph1, t)

%ph1 = -ph1;
x = r0.*cos(w0.*t+ph0) + r1.*cos(w0.*t+ph0).*cos(w1.*t+ph1) - r1.*cos(a).*sin(w0.*t+ph0).*sin(w1.*t+ph1);
y = (r0.*sin(w0.*t+ph0) + r1.*sin(w0.*t+ph0).*cos(w1.*t+ph1) + r1.*cos(a).*cos(w0.*t+ph0).*sin(w1.*t+ph1));
z = w0*t*r0/tan(a) - r1*sin(a)*sin(w1*t+ph1); % Crick exactly, as published in Acta Cryst. (1953). 6, 685
%z = -w0.*t.*r0./tan(a) + r1.*sin(a).*sin(w1.*t+ph1); % this is exactly negative relative to what is published in Acta Cryst. (1953). 6, 685
%z = abs(w0)*t*r0/tan(a) + r1*sin(a)*sin(w1*t+ph1);


function XYZ = generateCrickBB(chains, chL, r0, r1, w0, w1, a, ph1, cr, dph0, zoff, varargin)

if (~isempty(varargin))
    opts = varargin{1};
else
    opts = struct();
end

if (numel(cr) == chains-1) cr = [1 cr]; end
if (numel(dph0) == chains-1) dph0 = [0 dph0]; end
if (numel(zoff) == chains-1) zoff = [0 zoff]; end
if (cr(1) == 0) error('The first entry of the orientation vector can not be 0.'); end
if (numel(ph1) == 1)
    ph1 = ph1*ones(chains, 1);
elseif (numel(ph1) ~= chains)
    error('%d helical phases specified for %d chains', numel(ph1), chains);
end

XYZ = zeros(chL*chains, 3);
t = 0:chL-1;
for i = 1:chains
    if (cr(i) == 0)
        [x y z] = crickEQ(r0, r1, -w0, -w1, a, 0, -ph1(i), t);
        if (isfield(opts, 'apNNzoff'))
            zoff(i) = zoff(i) + XYZ((i-1)*chL+1, 3) - z(end);
        elseif (isfield(opts, 'registerzoff'))
            zo = XYZ(1, 3) - z(end); % start with z-offset that brings termini together. What register z-offset does this give?
            dz = absoluteToRegisterZoff(zo, r0, w0, a, w1, ph1(1), ph1(i), cr(i));
            % correct initial z-offset to give one with a 0 register
            % z-offset - this will be the basal z-offset, to which the
            % specified offset will be added.
            zo = zo - dz;
            zoff(i) = zoff(i) + zo;
        elseif (isfield(opts, 'zoffaa'))
            zo = XYZ(1, 3) - z(end); % start with z-offset that brings termini together. What Zaa does this give?
            dz = absoluteToZoff_aa(zo, r0, w0, a, r1, w1, ph1(1), ph1(i), cr(i));
            % correct initial z-offset to give one with a 0 Z_aa
            % offset - this will be the basal z-offset, to which the
            % specified offset will be added.
            zo = zo - dz;
            zoff(i) = zoff(i) + zo;
        end
        T = [cos(dph0(i) - zoff(i)*tan(a)/r0)  sin(dph0(i) - zoff(i)*tan(a)/r0)  0;
             -sin(dph0(i) - zoff(i)*tan(a)/r0)  cos(dph0(i) - zoff(i)*tan(a)/r0)  0;
              0        0        1];
    else
        [x y z] = crickEQ(r0, r1, w0, w1, a, 0, ph1(i), t);
        if (isfield(opts, 'registerzoff'))
            zo = 0; % start with z-offset that brings termini together. What register z-offset does this give?
            dz = absoluteToRegisterZoff(zo, r0, w0, a, w1, ph1(1), ph1(i), cr(i));
            % correct initial z-offset to give one with a 0 register
            % z-offset - this will be the basal z-offset, to which the
            % specified offset will be added.
            zo = zo - dz;
            zoff(i) = zoff(i) + zo;
        elseif (isfield(opts, 'zoffaa'))
            zo = 0; % start with z-offset that brings termini together. What Zaa does this give?
            dz = absoluteToZoff_aa(zo, r0, w0, a, r1, w1, ph1(1), ph1(i), cr(i));
            % correct initial z-offset to give one with a 0 Z_aa
            % offset - this will be the basal z-offset, to which the
            % specified offset will be added.
            zo = zo - dz;
            zoff(i) = zoff(i) + zo;
        end
        T  = [cos(dph0(i) - zoff(i)*tan(a)./r0)  sin(dph0(i) - zoff(i)*tan(a)./r0)  repmat(0, 1, length(r0));
             -sin(dph0(i) - zoff(i)*tan(a)./r0)  cos(dph0(i) - zoff(i)*tan(a)./r0)  repmat(0, 1, length(r0));
              repmat(0, 1, length(r0)) repmat(0, 1, length(r0)) repmat(1, 1, length(r0))];
    end
    if (length(r0) == chL)
        xyz = [x; y; z];
        for k = 1:length(r0)
            xyz(:, k) = T(:, k:chL:end)*xyz(:, k);
        end
    else
        xyz = T*([x; y; z]);
    end
    XYZ((i-1)*chL+1:i*chL, 1) = xyz(1, :)';
    XYZ((i-1)*chL+1:i*chL, 2) = xyz(2, :)';
    XYZ((i-1)*chL+1:i*chL, 3) = xyz(3, :)' + zoff(i);
end

function [v x] = shift(x, varargin)

n = 1;
if (nargin > 1)
    n = varargin{1};
end
if (n > length(x))
    error('vector does not have %d values...', n);
end
v = x(1:n); x = x(n+1:end);


% [rmsd, M] = superimpose(X, Y)
% For two 3-by-N matrices, X and Y, returns the best fit root mean square
% deviation between the two, resulting from an optimal alignment
% (a translation and a rotation by some angle around some axis) of X onto
% Y. The matrix performing this optimal linear transformation is returned
% as well. Finally, optionally, the third output argument is the array of
% residual distances from Y onto the optimal superposition of X.
%
% Implements the SVD method by Kabsh et al. (Acta Crystallogr 1976, A32,
% 922) and taken from Coutsias et al. J Comp Chem 2004, 25(15) 1849.
%
function [rmsd, M, varargout] = superimpose(X, Y)

if ((size(X, 1) ~= 3) || (size(Y, 1) ~= 3)) error('Matrices X and Y must be 3-by-N'); end
if (size(X, 2) ~= size(Y, 2)) error('Matrices X and Y must be of the same size'); end
N = size(X, 2);

X = X - repmat(mean(X, 2), 1, N);
Y = Y - repmat(mean(Y, 2), 1, N);
R = X*(Y');
[V S W] = svd(R);
I = eye(3,3);
sgn = sign(det(R));
I(3,3) = sgn;
M = W*I*V';
rmsd = sqrt((sum(sum(X.^2) + sum(Y.^2)) - 2*(S(1,1) + S(2,2) + sgn*S(3,3)))/N);
if (nargout > 2)
    varargout{1} = sqrt(sum((Y - M*X).^2));
end


function rzoff = absoluteToRegisterZoff(zoff, R0, w0, a, w1, ph1_1, ph1_2, p_ap)
%absoluteToRegisterZoff Converts from an absolute Z offset to a register Z
%offset (defined as the Z distance between points most innert to the
%bundle - i.e. points with phase pi)

assert((p_ap == 1) || (p_ap == 0), 'p_ap flag expected to be either 1 or -1');

aa1 = 2*pi/w1; b1 = (pi - ph1_1)/w1; z1 = (R0*w0/tan(a))*(aa1*1 + b1);
if (p_ap == 0)
    % for a chain running in the opposite orientation, the meaning
    % of clockwise-ness changes, so w0, w1 and phase flip sign
    aa2 = 2*pi/(-w1);
    b2 = (pi + ph1_2)/(-w1);
    n = ((z1 - zoff)/(-R0*w0/tan(a)) - b2)/aa2;
    dz = -(R0*w0/tan(a))*(aa2*floor(n) + b2) + zoff - z1;
    dz1 = -(R0*w0/tan(a))*(aa2*ceil(n) + b2) + zoff - z1;
else
    b2 = (pi - ph1_2)/w1;
    n = ((z1 - zoff)/(R0*w0/tan(a)) - b2)/aa1;
    dz = (R0*w0/tan(a))*(aa1*floor(n) + b2) + zoff - z1;
    dz1 = (R0*w0/tan(a))*(aa1*ceil(n) + b2) + zoff - z1;
end

if (abs(dz1) < abs(dz))
    rzoff = dz1;
else
    rzoff = dz;
end


function ind = mystrfind(str, list)
ind = find(strncmp(str, list, length(str)));
