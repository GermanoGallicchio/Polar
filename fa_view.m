function fa_view(theta, r, fa_cfg)
% SYNTAX
%       fa_view(theta, r, fa_Cfg)
%
% DESCRIPTION: 
%       Visualize polar data in angular and linear coordinates
%
% INPUT:        
%       theta       - angles (radians)
%       r           - amplitudes
%       fa_cfg structure with field "viewParams"
%       viewParams  - structure with optional fields
%                       .thetaLim   (default: [-pi pi])
%                       .thetaStep  (default: pi/4)
%                       .thetaTickValues (numeric vector, optional)
%                       .thetaTickLabels (string/cellstr, optional)
%                       .thetaTickLabelInterpreter ('tex'|'latex'|'none', optional)
%                       .xLabel     (char/string scalar, default: 'Argument')
%                       .yLabel     (char/string scalar, default: 'Modulus')
%                       .units (structure, optional)
%                           .name  (default: 'radians')
%                           .lim   (default: thetaLim; e.g., [0 1000] for ms labels)
%                       .type       ('line' or 'point', default: 'line')
%                       .drawMeanResultant (logical, default: false)
%
% AUTHOR:
%   Germano Gallicchio (germano.gallicchio@gmail.com)

%% input check and defaults

% sanity check: fa_cfg.viewParams field exists
if ~isfield(fa_cfg,'viewParams')
    error('fa_cfg structure needs to have a "viewParams" field. see fa_documentation() for details')
end


if ~isfield(fa_cfg.viewParams,'thetaLim')
    fa_cfg.viewParams.thetaLim = [-pi pi];
end
if ~isfield(fa_cfg.viewParams,'thetaStep')
    fa_cfg.viewParams.thetaStep = pi/4;
end
if ~isfield(fa_cfg.viewParams,'type')
    fa_cfg.viewParams.type = 'line';
end
if ~isfield(fa_cfg.viewParams,'drawMeanResultant')
    fa_cfg.viewParams.drawMeanResultant = false;
end
if ~isfield(fa_cfg.viewParams,'xLabel')
    fa_cfg.viewParams.xLabel = 'Argument';
end
if ~isfield(fa_cfg.viewParams,'yLabel')
    fa_cfg.viewParams.yLabel = 'Modulus';
end
if ~islogical(fa_cfg.viewParams.drawMeanResultant) || numel(fa_cfg.viewParams.drawMeanResultant)~=1
    error('fa_cfg.viewParams.drawMeanResultant must be a logical scalar. see fa_documentation() for details')
end
if ~(ischar(fa_cfg.viewParams.xLabel) || (isstring(fa_cfg.viewParams.xLabel) && numel(fa_cfg.viewParams.xLabel)==1))
    error('fa_cfg.viewParams.xLabel must be a char or string scalar. see fa_documentation() for details')
end
if ~(ischar(fa_cfg.viewParams.yLabel) || (isstring(fa_cfg.viewParams.yLabel) && numel(fa_cfg.viewParams.yLabel)==1))
    error('fa_cfg.viewParams.yLabel must be a char or string scalar. see fa_documentation() for details')
end
if ~isfield(fa_cfg.viewParams,'units')
    fa_cfg.viewParams.units = struct();
end
if ~isfield(fa_cfg.viewParams.units,'name')
    fa_cfg.viewParams.units.name = 'radians';
end
if ~isfield(fa_cfg.viewParams.units,'lim')
    fa_cfg.viewParams.units.lim = fa_cfg.viewParams.thetaLim;
end
if ~(ischar(fa_cfg.viewParams.units.name) || (isstring(fa_cfg.viewParams.units.name) && numel(fa_cfg.viewParams.units.name)==1))
    error('fa_cfg.viewParams.units.name must be a char or string scalar. see fa_documentation() for details')
end
if ~(isnumeric(fa_cfg.viewParams.units.lim) && numel(fa_cfg.viewParams.units.lim)==2)
    error('fa_cfg.viewParams.units.lim must be a numeric 1x2 vector. see fa_documentation() for details')
end

%% shortcuts

nSamples  = numel(theta);
thetaLim  = fa_cfg.viewParams.thetaLim;
thetaStep = fa_cfg.viewParams.thetaStep;
unitsName = char(fa_cfg.viewParams.units.name);
unitsLim = fa_cfg.viewParams.units.lim;
xLabel = char(fa_cfg.viewParams.xLabel);
yLabel = char(fa_cfg.viewParams.yLabel);

if isfield(fa_cfg.viewParams,'thetaTickValues')
    thetaTicks = fa_cfg.viewParams.thetaTickValues;
else
    thetaTicks = thetaLim(1):thetaStep:thetaLim(2);
end
if ~isnumeric(thetaTicks) || isempty(thetaTicks)
    error('fa_cfg.viewParams.thetaTickValues must be a non-empty numeric vector. see fa_documentation() for details')
end
thetaTicks = thetaTicks(:)';

if isfield(fa_cfg.viewParams,'thetaTickLabels')
    thetaTickLabels = fa_cfg.viewParams.thetaTickLabels;

    if isstring(thetaTickLabels)
        thetaTickLabels = cellstr(thetaTickLabels(:));
    elseif ischar(thetaTickLabels)
        thetaTickLabels = {thetaTickLabels};
    elseif iscell(thetaTickLabels)
        thetaTickLabels = cellfun(@char, thetaTickLabels, 'UniformOutput', false);
    else
        error('fa_cfg.viewParams.thetaTickLabels must be string, char, or cellstr. see fa_documentation() for details')
    end

    if numel(thetaTickLabels) ~= numel(thetaTicks)
        error('fa_cfg.viewParams.thetaTickLabels must have the same length as thetaTickValues. see fa_documentation() for details')
    end

    if isfield(fa_cfg.viewParams,'thetaTickLabelInterpreter')
        thetaTickLabelInterpreter = fa_cfg.viewParams.thetaTickLabelInterpreter;
    else
        thetaTickLabelInterpreter = 'none';
    end
else
    % Default labels: preserve radians style unless user requests a different unit/range.
    if strcmpi(unitsName,'radians') && isequal(unitsLim, thetaLim)
        thetaTickLabels = arrayfun(@(x) sprintf('%.2g\\pi', x/pi), thetaTicks, 'UniformOutput', false);
        thetaTickLabelInterpreter = 'tex';
    else
        if diff(thetaLim)==0
            error('fa_cfg.viewParams.thetaLim must span a non-zero range when auto-mapping unit labels. see fa_documentation() for details')
        end
        tickValsUnits = unitsLim(1) + (thetaTicks - thetaLim(1)) * (diff(unitsLim) / diff(thetaLim));
        thetaTickLabels = arrayfun(@(x) sprintf('%.4g', x), tickValsUnits, 'UniformOutput', false);
        thetaTickLabelInterpreter = 'none';
    end
end

% wrap the thetaTickLabels around to show coincidence of first and last point, only if theta ticks span a full circle (i.e., 2*pi). 
% This is to highlight that the first and last points are the same angle.
span = thetaTicks(end) - thetaTicks(1);
tol = 1e-10;
isFullCircle = numel(thetaTicks) > 1 && abs(mod(span, 2*pi)) < tol;
if isFullCircle
    thetaTickLabels_first = thetaTickLabels{1};
    thetaTickLabels_last  = thetaTickLabels{end};
    thetaTickLabels{1}  = [thetaTickLabels_first '=' thetaTickLabels_last];
    thetaTickLabels{end} = [thetaTickLabels_last '=' thetaTickLabels_first];
end



if fa_cfg.viewParams.drawMeanResultant
    fa_cfg_metrics = struct();
    fa_cfg_metrics.metrics.requests = ["meanResultantLength" "meanResultantAngle"];
    metrics = fa_meanResultant(theta(:), r(:), fa_cfg_metrics);
    meanResultantAngle = metrics.meanResultantAngle;
    meanResultantLength = metrics.meanResultantLength;

    % map angle to the displayed linear-space range (assumes a 2*pi span)
    meanResultantAngleLinear = mod(meanResultantAngle - thetaLim(1), 2*pi) + thetaLim(1);
end


%% general figure settings

tld = tiledlayout(2,1); 
tld.TileSpacing = 'compact';
%tld.Padding = 'compact';

%% panel 1: angular space

nexttile(tld,1,[1 1])
c = lines(1); % base color


switch fa_cfg.viewParams.type
    case 'point'
        polarplot(theta,r,'.','Color',c);

    case 'line'
        if nSamples <= 200
            polarplot([zeros(size(theta)) theta]', ...
                [zeros(size(theta)) r]', ...
                '-', 'Color', c);
        else
            alphaVal = min(1,200/nSamples);
            polarplot([zeros(size(theta)) theta]', ...
                [zeros(size(theta)) r]', ...
                '-', 'Color', [c alphaVal]);
        end

end

        axPolar = gca;
        axPolar.ThetaAxisUnits = "radians"; % plotting is always radians under the hood
        axPolar.ThetaLim = thetaLim;
        axPolar.ThetaTick = thetaTicks;
        axPolar.ThetaTickLabel = thetaTickLabels;
        axPolar.ThetaZeroLocation = 'right';
        axPolar.ThetaDir = 'counterclockwise';

% optional overlay: mean resultant vector
if fa_cfg.viewParams.drawMeanResultant
    hold on
    polarplot([meanResultantAngle meanResultantAngle], [0 meanResultantLength], ...
        '-', 'LineWidth', 2, 'Color', [1 0.2 0]);
end


%% panel 2: linear space

nexttile(tld,2,[1 1])

switch fa_cfg.viewParams.type
    case 'point'
        
        plot(theta, r, '.', 'Color', c);
    case 'line'

        if nSamples <= 200
            stem(theta, r, 'Marker','.', 'Color', c);
        else
            alphaVal = min(1,200/nSamples);
            plot([theta theta]', [zeros(size(r)) r]', ...
                '-', 'Color', [c alphaVal]);
        end
end

        ax = gca;
        ax.XAxis.Label.String = xLabel;
        ax.XAxis.Limits = thetaLim;
        ax.XAxis.TickValues = thetaTicks;
        ax.XAxis.TickLabels = thetaTickLabels;
        ax.XAxis.TickLabelInterpreter = thetaTickLabelInterpreter;
        ax.YAxis.Label.String = yLabel;
        ax.YAxis.Limits = unique([0 max(r(:))]) + [0 1]*0.1;
        ax.YAxis.TickValues = unique([0 max(r(:))]);

if fa_cfg.viewParams.drawMeanResultant
    hold(ax, 'on')
    mr = plot(ax, [meanResultantAngleLinear meanResultantAngleLinear], [0 meanResultantLength], ...
        '-', 'LineWidth', 2, 'Color', [1 0.2 0]);

    legend(mr, 'mean resultant', 'Location', 'northoutside');

    
end
