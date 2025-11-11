function po_view(theta, r, viewParams)
% SYNTAX
%       po_view(theta, r, viewParams)
%
% DESCRIPTION: 
%       Visualize polar data in angular and linear coordinates
%
% INPUT:        
%       theta       - angles (radians)
%       r           - amplitudes
%       viewParams  - structure with optional fields
%                       .thetaLim   (default: [-pi pi])
%                       .thetaStep  (default: pi/4)
%                       .type       ('line' or 'point', default: 'line')
%
% AUTHOR:
%   Germano Gallicchio (germano.gallicchio@gmail.com)

%% set defaults

if ~isfield(viewParams,'thetaLim')
    viewParams.thetaLim = [-pi pi];
end
if ~isfield(viewParams,'thetaStep')
    viewParams.thetaStep = pi/4;
end
if ~isfield(viewParams,'type')
    viewParams.type = 'line';
end

%% shortcuts

nSamples  = numel(theta);
thetaLim  = viewParams.thetaLim;
thetaStep = viewParams.thetaStep;


%% general figure settings

tld = tiledlayout(2,1);
thetaTicks = thetaLim(1):thetaStep:thetaLim(2);

%% panel 1: angular space

%% panel 1: angular space

nexttile(tld,1,[1 1])
c = lines(1); % base color


switch viewParams.type
    case 'point'
        pp = polarplot(theta,r,'.','Color',c);
        pp.Parent.ThetaAxisUnits = "radians";
        pp.Parent.ThetaLim = thetaLim;
        pp.Parent.ThetaTick = thetaTicks;
        pp.Parent.ThetaZeroLocation = 'right';
        pp.Parent.ThetaDir = 'counterclockwise';

    case 'line'
        if nSamples <= 200
            pp = polarplot([zeros(size(theta)) theta]', ...
                [zeros(size(theta)) r]', ...
                '-', 'Color', c);
        else
            alphaVal = min(1,200/nSamples);
            pp = polarplot([zeros(size(theta)) theta]', ...
                [zeros(size(theta)) r]', ...
                '-', 'Color', [c alphaVal]);
        end
        pp(1).Parent.ThetaAxisUnits = "radians";
        pp(1).Parent.ThetaLim = thetaLim;
        pp(1).Parent.ThetaTick = thetaTicks;
        pp(1).Parent.ThetaZeroLocation = 'right';
        pp(1).Parent.ThetaDir = 'counterclockwise';
end


%% panel 2: linear space

nexttile(tld,2,[1 1])

switch viewParams.type
    case 'point'
        
        lp = plot(theta, r, '.', 'Color', c);
        ax = lp.Parent;
        ax.XAxis.Label.String = 'Phase [radians]';
        ax.XAxis.Limits = thetaLim;
        ax.XAxis.TickValues = thetaTicks;
        ax.XAxis.TickLabels = arrayfun(@(x) sprintf('%.2g\\pi', x/pi), ...
            thetaTicks, 'UniformOutput', false);
        ax.XAxis.TickLabelInterpreter = 'tex';
        ax.YAxis.Label.String = 'Density';
        ax.YAxis.Limits = unique([0 prctile(r(:),[0 100])]);
        ax.YAxis.TickValues = unique([0 prctile(r(:),[0 100])]);

        ax.YAxis.Limits = unique([0 max(r(:))]) + [0 1]*0.1;
        ax.YAxis.TickValues = unique([0 max(r(:))]);

        

        

    case 'line'

        if nSamples <= 200
            lp = stem(theta, r, 'Marker','.', 'Color', c);
        else
            alphaVal = min(1,200/nSamples);
            lp = plot([theta theta]', [zeros(size(r)) r]', ...
                '-', 'Color', [c alphaVal]);
        end
        ax = lp(1).Parent;
        ax.XAxis.Label.String = 'Phase [radians]';
        ax.XAxis.Limits = thetaLim;
        ax.XAxis.TickValues = thetaTicks;
        ax.XAxis.TickLabels = arrayfun(@(x) sprintf('%.2g\\pi', x/pi), ...
            thetaTicks, 'UniformOutput', false);
        ax.XAxis.TickLabelInterpreter = 'tex';
        ax.YAxis.Label.String = 'Density';

        ax.YAxis.Limits = unique([0 max(r(:))]) + [0 1]*0.1;
        ax.YAxis.TickValues = unique([0 max(r(:))]);
end
