function po_view(theta, rho, po_cfg, viewParams)

% DESCRIPTION:
%
% INPUT:        
%
%
%
% AUTHOR:
%   Germano Gallicchio (germano.gallicchio@gmail.com)

%% set defaults

fieldLbl = fieldnames(viewParams);


varLbl = 'thetaLim';
if ~any(strcmp(fieldLbl,varLbl))
    viewParams.(varLbl) = [-pi  pi];
end

varLbl = 'thetaStep';
if ~any(strcmp(fieldLbl,varLbl))
    viewParams.(varLbl) = pi/4;
end


%% shortcuts

nSamples = size(theta,1);


thetaLim  = viewParams.thetaLim;
thetaStep = viewParams.thetaStep;

%% general figure settings

tld = tiledlayout(2,1);
thetaTicks = thetaLim(1) : thetaStep : thetaLim(2);



%% panel 1: angular space

%
nexttile(tld,1,[1 1])

switch viewParams.type
    case 'point'
        pp = polarplot(theta,rho,'.');
        pp.Parent.ThetaAxisUnits = "radians";
        pp.Parent.ThetaLim = thetaLim;
        pp.Parent.ThetaTick = thetaTicks;

    case 'line'
         if nSamples<=200
            pp = polarplot([zeros(size(theta)) theta]',[zeros(size(theta)) rho]','-','Color',lines(1));
            pp(1).Parent.ThetaAxisUnits = "radians";
            pp(1).Parent.ThetaLim = thetaLim;
            pp(1).Parent.ThetaTick = thetaTicks;
         else
             colAlpha = length(theta)/(length(theta)/200*length(theta));
            pp = polarplot([zeros(size(theta)) theta]',[zeros(size(theta)) rho]','-','Color',[lines(1) colAlpha]);
            pp(1).Parent.ThetaAxisUnits = "radians";
            pp(1).Parent.ThetaLim = thetaLim;
            pp(1).Parent.ThetaTick = thetaTicks;
        end
end


%% panel 2: linear space

nexttile(tld,2,[1 1])

switch viewParams.type
    case 'point'
        
        lp = plot(theta, rho, '.');
        lp.Parent.XAxis.Label.String = 'Phase [radians]';
        lp.Parent.XAxis.Limits = thetaLim;
        lp.Parent.XAxis.TickValues = thetaTicks;
        lp.Parent.XAxis.TickLabels = lp.Parent.XAxis.TickValues/pi + "\pi";
        lp.Parent.YAxis.Label.String = 'Density';

    case 'line'

        if nSamples<=200
            lp = stem(theta, rho,'Marker','.');
            lp.Parent.XAxis.Label.String = 'Phase [radians]';
            lp.Parent.XAxis.Limits = thetaLim;
            lp.Parent.XAxis.TickValues = thetaTicks;
            lp.Parent.XAxis.TickLabels = lp.Parent.XAxis.TickValues/pi + "\pi";
            lp.Parent.YAxis.Label.String = '';
            lp.Parent.YAxis.TickValues = [0 1];
        else

            colAlpha = length(theta)/(length(theta)/200*length(theta));
            lp = plot([theta theta]', [zeros(size(rho)) rho]','-','Color',[lines(1) colAlpha]);
            lp(1).Parent.XAxis.Label.String = 'Phase [radians]';
            lp(1).Parent.XAxis.Limits = thetaLim;
            lp(1).Parent.XAxis.TickValues = thetaTicks;
            lp(1).Parent.XAxis.TickLabels = lp(1).Parent.XAxis.TickValues/pi + "\pi";
            lp(1).Parent.YAxis.Label.String = '';
            lp(1).Parent.YAxis.TickValues = [0 1];
        end
end
