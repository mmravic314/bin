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
%generateCrickBB Generate a coiled-coil backbone from Crick parameters.
%   xyz = generateCrickBB(cN, chL, R0, R1, w0, w1, A, ph1, cr, dph0, Zoff)
%   Generates an ideal coiled-coil backbone from given Crick parameters.
%   The output variable is a N-by-3 array of coordinates of backbone atoms
%   (one per residue). If you use this program in your research, please
%   cite G. Grigoryan, W. F. DeGrado, "Probing Designability via a
%   Generalized Model of Helical Bundle Geometry", J. Mol. Biol., 405(4):
%   1079-1100 (2011). See also http://grigoryanlab.org/cccp/.
%
%   cN is the number of chains. chL is the length of each chain. R0 and R1
%   are the superhelical and helical radii, respectively (in Angstroms). w0
%   and w1 are superhelical and helical frequencies, respectively, in
%   degrees per residue. A is the pith angle in degrees.
%   ph1 -- the helical phase angle in degrees. Can be either a scalar, in
%   which case all chains have the same helical phase angle, or vector of
%   length cN, in which case chain k will have helical phase angle ph1(k).
%   Parameters cr, dph0, and Zoff have to be vectors of length cN-1, and
%   the value in the k-th element of the vector corresponds to chain k+1.
%   dph0 -- superhelical chain offset, in degrees. For chain k (k > 1),
%   dph0(k-1) will indicate its superhelical phase offset with respect to
%   chain 1.
%   cr   -- chain orientation. For each chain k (k > 1), its orientation is
%   parallel to the first chain if cr(k-1) is 1, and anti-parallel to the
%   first chain if cr(k-1) is 0.
%   Zoff -- individual chain Z-offsets. For each chain k (k > 1), Zoff(k-1)
%   denotes the Z offset of chain k with respect to the first chain. By
%   default, the Z offset refers to an absolute offset between chains
%   following the application of Crick equations. This is almost always not
%   useful in generation/fitting, and the specific meaning of Z offset can
%   be controlled by the optional last parameter.
%
%   generateCrickBB(cN, chL, R0, R1, w0, w1, A, ph1, cr, dph0, Zoff, zType)
%   zType indicates the type of Z offset to use. It should be a struct with
%   one field (i.e. struct('field', 1)). Possible field values and their
%   meanings are (see Grigoryan, DeGrado, JMB 405:1079, 2011 for exact
%   definitions and details):
%   'zoffaa' -- the Z offset is interpreted as the offset between 'a'
%   positions on opposing chains.
%   'registerzoff' -- the Z offset refers to the offset between points on
%   the Crick curves of opposing chains that point directly into the
%   interface.
%   'apNNzoff' -- the Z offset is interpreted as between the N termini of
%   chains 1 and k if k runs parallel to 1, and the N terminus of chain 1
%   and the C terminus of chain k if k is anti-parallel to 1.
% 
function XYZ = generateCrickBB(chains, chL, r0, r1, w0, w1, a, ph1, cr, dph0, zoff, varargin)

if (~isempty(varargin))
    opts = varargin{1};
else
    opts = struct();
end

% convert degrees to radians
w0 = w0*pi()/180;
w1 = w1*pi()/180;
a = a*pi()/180;
ph1 = ph1*pi()/180;
dph0 = dph0*pi()/180;

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



function zaa = absoluteToZoff_aa(zoff, R0, w0, a, R1, w1, ph1_1, ph1_2, p_ap)
%absoluteToZoff_aa Converts from an absolute Z offset to a Zaa' Z offset,
%defined as the Z distance between closest 'a' positions on opposite chains

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


function d = angleDiff (a, b)

d = mod((mod(a, 2*pi) - mod(b, 2*pi)), 2*pi);
in = find(d > pi);
d(in) = d(in) - 2*pi;


function ph = canonicalPhases(ind)

% canonical phases are:
% [41.0, 95.0, 146.0, 197.0, 249.0, 300.0, 351.0]
% corresponding to {'c', 'g', 'd', 'a', 'e', 'b', 'f'}, respectively

% in the order a-g:
median_phases = [197.0, 300.0, 41.0, 146.0, 249.0, 351.0 95.0];
ph = median_phases(ind);


function [x, y, z] = crickEQ(r0, r1, w0, w1, a, ph0, ph1, t)

%ph1 = -ph1;
x = r0.*cos(w0.*t+ph0) + r1.*cos(w0.*t+ph0).*cos(w1.*t+ph1) - r1.*cos(a).*sin(w0.*t+ph0).*sin(w1.*t+ph1);
y = (r0.*sin(w0.*t+ph0) + r1.*sin(w0.*t+ph0).*cos(w1.*t+ph1) + r1.*cos(a).*cos(w0.*t+ph0).*sin(w1.*t+ph1));
z = w0*t*r0/tan(a) - r1*sin(a)*sin(w1*t+ph1); % Crick exactly, as published in Acta Cryst. (1953). 6, 685
%z = -w0.*t.*r0./tan(a) + r1.*sin(a).*sin(w1.*t+ph1); % this is exactly negative relative to what is published in Acta Cryst. (1953). 6, 685
%z = abs(w0)*t*r0/tan(a) + r1*sin(a)*sin(w1*t+ph1);


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
