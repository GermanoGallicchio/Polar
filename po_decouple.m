function [theta_decoupled, r_decoupled, offset_pnt ] = po_decouple(theta, r, po_cfg)
% SYNTAX:
%   [theta_decoupled, r_decoupled, offset_pnt] = po_decouple(theta, r, po_cfg)
%
% DESCRIPTION:
%   Decouple theta (angles) and r (magnitudes) by circularly shifting r relative to theta
%   along the first dimension by a random offset within a specified range.
%   Useful for breaking within-sample phase-amplitude coupling while preserving
%   marginal distributions and temporal/cycle structure within each vector.
%
% INPUT:
%   theta   - numeric array of angles (radians). Samples are along the first dimension
%   r       - numeric array of magnitudes, same size as theta.
%   po_cfg  - struct with field 'decoupling':
%               .offsetRange   [1x2 numeric] percentage range of possible offsets
%                               relative to the length of first dimension. For example, [20 80]
%                               means choose a random offset (by sampling from a uniform distribution) between
%                               ceil(0.20*N) to floor(0.80*N), where N=size(theta,1).
%               .randomSeed    [numeric] (optional) seed for reproducible offset
%               .perColumn     [logical] (optional, default=true) if true, each column
%                               gets an independent random shift; if false, all columns
%                               shifted by the same offset
%
% OUTPUT:
%   theta_decoupled - same as theta (no change)
%   r_decoupled     - r circularly shifted by offset_pnt along dim=1
%   offset_pnt      - integer shift(s) (in samples) applied to r.
%                     Scalar if perColumn=false; vector [1 x nCols] if perColumn=true
%
% AUTHOR:
%   Germano Gallicchio (germano.gallicchio@gmail.com)

%% input checks and defaults


% theta and r are numeric arrays of the same size
if ~isnumeric(theta) || ~isnumeric(r)
    error('theta and r must be numeric arrays. See po_documentation() for details');
end
if ~isequal(size(theta), size(r))
    error('theta and r must have the same size. See po_documentation() for details');
end

% theta and r have at least 2 values along the first dimension
if size(theta,1) < 2
    error('theta and r must have at least 2 values (i.e., samples) along the first dimension. See po_documentation() for details');
end

% po_cfg.decoupling field exists
if ~isfield(po_cfg,'decoupling')
    error('po_cfg must contain a decoupling field. See po_documentation() for details');
end

% po_cfg.decoupling.offsetRange field exists, otherwise set default
if ~isfield(po_cfg.decoupling, 'offsetRange')
    po_cfg.decoupling.offsetRange = [20 80]; % percentage
    warning(['po_cfg.decoupling.offsetRange was not defined. Setting to ' num2str(po_cfg.decoupling.offsetRange) ' (percentage of the sample length). Specify explicitly to avoid this warning.'])
end

% po_cfg.decoupling.offsetRange must be a 1x2 numeric vector
if ~isnumeric(po_cfg.decoupling.offsetRange) || numel(po_cfg.decoupling.offsetRange)~=2
    error('po_cfg.decoupling.offsetRange must be a 1x2 numeric vector of percentages [min max]. See po_documentation() for details');
end

% offset range must be sorted [min max]
if ~issorted(po_cfg.decoupling.offsetRange)
    error('po_cfg.decoupling.offsetRange must be sorted in ascending order [min max]. See po_documentation() for details');
end

% offset range must be within [0,100]
if any(po_cfg.decoupling.offsetRange < 0) || any(po_cfg.decoupling.offsetRange > 100)
    error('po_cfg.decoupling.offsetRange values must be between 0 and 100 (percentages). See po_documentation() for details');
end

% if randomSeed not set by user, set it to 42
% (randomSeed used to make stochastic simulations replicable)
if ~isfield(po_cfg.decoupling,'randomSeed')
    po_cfg.decoupling.randomSeed = 42;
    warning(['po_cfg.decoupling.randomSeed not defined by user. i am setting it to ' num2str(po_cfg.decoupling.randomSeed) ' , but you can specify it explicitly to avoid this warning. This is not a bug. See po_documentation() for details'])
end

% if perColumn not set by user, default to true (column-wise shift--each column gets an independent shift)
if ~isfield(po_cfg.decoupling,'perColumn')
    po_cfg.decoupling.perColumn = true;
    warning(['po_cfg.decoupling.perColumn not defined by user. i am setting it to true , but you can specify it explicitly to avoid this warning. This is not a bug. See po_documentation() for details'])
end

%% shortcuts

nSamples = size(theta,1); % number of samples (from first dimension)
nCols = size(theta,2);     % number of columns

offsetRange_prc = po_cfg.decoupling.offsetRange;

randomSeed = po_cfg.decoupling.randomSeed;

perColumn = po_cfg.decoupling.perColumn;

%% implementation

dimShift = 1; % dimension where operating the circular shift

% convert to integer sample bounds
lo = ceil(offsetRange_prc(1) * nSamples / 100);
hi = floor(offsetRange_prc(2) * nSamples / 100);

% random seed reproducibility
prevRng = rng; % store whichever was set before
rng(randomSeed); % set new seed for this function

% determine offset point(s)
if perColumn
    % Independent offset per column
    offset_pnt = randi([lo, hi], 1, nCols);
    
    % apply column-wise circular shifts
    r_decoupled = zeros(size(r)); % initialize
    for colIdx = 1:nCols
        r_decoupled(:, colIdx) = circshift(r(:, colIdx), offset_pnt(colIdx), dimShift);
    end
else
    % Single global offset for all columns
    offset_pnt = randi([lo, hi], 1, 1);
    
    r_decoupled = circshift(r, offset_pnt, dimShift);
end

theta_decoupled = theta; % no changes

% restore random seed
rng(prevRng);


