function [ControlPoints, InducedPoints, Normals, InfiniteVortices, Vortices, internalMesh, WingExtremes2Export] = createStructure(config, iCorpo)

YAngle = config.RotationAngle_Y(iCorpo);
YRot = [ cosd(YAngle)     0        sind(YAngle);
              0           1            0
        -sind(YAngle)     0        cosd(YAngle)];
    
XAngle = config.DihedralAngle(iCorpo);
XRot = [ 1         0              0      ;
         0    cosd(XAngle)  -sind(XAngle);
         0    sind(XAngle)   cosd(XAngle) ];

WingExtremes.RootLE = [0 0 0];
WingExtremes.RootTE = WingExtremes.RootLE + (YRot*[config.RootChord(iCorpo) 0 0]')';
% WingExtremes.TipLE = [0 config.SemiSpan(iCorpo)*cosd(config.DihedralAngle(iCorpo)) 0];
WingExtremes.TipLE = [0 config.SemiSpan(iCorpo) 0];
% Adding Sweep angle
WingExtremes.TipLE(1) = WingExtremes.RootLE(1) + config.RootChord(iCorpo)/4 + config.SemiSpan(iCorpo) * tand(config.SweepAngle(iCorpo)) - config.TipChord(iCorpo)/4;

% Adding Dihedral angle
WingExtremes.TipLE = (XRot*WingExtremes.TipLE')';


% Add Y-axis rotation
WingExtremes.TipLE = (YRot*WingExtremes.TipLE')';


WingExtremes.TipTE = WingExtremes.TipLE +  (YRot*[config.TipChord(iCorpo) 0 0]')';

WingExtremes.RootLE = WingExtremes.RootLE + [config.LEPosition_X(iCorpo) config.LEPosition_Y(iCorpo) config.LEPosition_Z(iCorpo)];
WingExtremes.RootTE = WingExtremes.RootTE + [config.LEPosition_X(iCorpo) config.LEPosition_Y(iCorpo) config.LEPosition_Z(iCorpo)];
WingExtremes.TipLE = WingExtremes.TipLE + [config.LEPosition_X(iCorpo) config.LEPosition_Y(iCorpo) config.LEPosition_Z(iCorpo)];
WingExtremes.TipTE = WingExtremes.TipTE + [config.LEPosition_X(iCorpo) config.LEPosition_Y(iCorpo) config.LEPosition_Z(iCorpo)];

Extremes = [WingExtremes.RootLE; WingExtremes.RootTE; WingExtremes.TipTE; WingExtremes.TipLE; WingExtremes.RootLE];


% LELength = norm(WingExtremes.RootLE - WingExtremes.TipLE);
% TELength = norm(WingExtremes.RootTE - WingExtremes.TipTE);


% Spanwise discretization
LEDiscr = linspace(0, 1, config.SemiSpanwiseDiscr(iCorpo)+1);
TEDiscr = linspace(0, 1, config.SemiSpanwiseDiscr(iCorpo)+1);

LELine = @(x) (WingExtremes.TipLE - WingExtremes.RootLE ).*x + WingExtremes.RootLE;
TELine = @(x) (WingExtremes.TipTE - WingExtremes.RootTE ).*x + WingExtremes.RootTE;

LEPoints = zeros(length(LEDiscr), 3);
TEPoints = zeros(length(TEDiscr), 3);

for i = 1:length(LEDiscr)
    LEPoints(i, :) = LELine(LEDiscr(i));
    TEPoints(i, :) = TELine(TEDiscr(i));
end


% Chordwise discretization

RootDiscr = linspace(0, 1, config.ChordwiseDiscr(iCorpo)+1);
TipDiscr = linspace(0, 1, config.ChordwiseDiscr(iCorpo)+1);

RootLine = @(x) (WingExtremes.RootTE - WingExtremes.RootLE ).*x + WingExtremes.RootLE;
TipLine = @(x) (WingExtremes.TipTE - WingExtremes.TipLE ).*x + WingExtremes.TipLE;

RootPoints = zeros(length(RootDiscr), 3);
TipPoints = zeros(length(RootDiscr), 3);

for i = 1:length(RootDiscr)
    RootPoints(i, :) = RootLine(RootDiscr(i));
    TipPoints(i, :) = TipLine(TipDiscr(i));
end


%% Create internal structured mesh

internalMesh = cell( config.ChordwiseDiscr(iCorpo), config.SemiSpanwiseDiscr(iCorpo)*2 );
Vortices = cell( config.ChordwiseDiscr(iCorpo), config.SemiSpanwiseDiscr(iCorpo)*2 );
ControlPoints = cell( config.ChordwiseDiscr(iCorpo), config.SemiSpanwiseDiscr(iCorpo)*2 );
InducedPoints = cell( config.ChordwiseDiscr(iCorpo), config.SemiSpanwiseDiscr(iCorpo)*2 );
Normals = cell( config.ChordwiseDiscr(iCorpo), config.SemiSpanwiseDiscr(iCorpo)*2 );

InfiniteVortices = cell( config.ChordwiseDiscr(iCorpo), config.SemiSpanwiseDiscr(iCorpo)*2 );

InfiniteVorticesLength = 50 * config.RootChord(iCorpo);

colormap = jet(length(RootDiscr)-1);

% Going in the chordwise direction
for i = 2:length(RootDiscr)
    LocalLE = @(x) ( TipPoints(i-1, :) - RootPoints(i-1, :) ).*x + RootPoints(i-1, :);
    LocalTE = @(x) ( TipPoints(i, :) - RootPoints(i, :) ).*x + RootPoints(i, :);
    
    LocalQuarterTip = (TipPoints(i, :) + 3.*TipPoints(i-1, :)) ./ 4;
    LocalQuarterRoot = (RootPoints(i, :) + 3.*RootPoints(i-1, :)) ./ 4;
    VortexFun = @(x) ( LocalQuarterTip  - LocalQuarterRoot ) .* x + LocalQuarterRoot;
    
    LocalThreeQuarterTip = 3.*(TipPoints(i, :) - TipPoints(i-1, :)) ./ 4 + TipPoints(i-1, :);
    LocalThreeQuarterRoot = 3.*(RootPoints(i, :) - RootPoints(i-1, :)) ./ 4 + RootPoints(i-1, :);
    ControlPointsFun = @(x) ( LocalThreeQuarterTip  - LocalThreeQuarterRoot ) .* x + LocalThreeQuarterRoot;
    
    % Going in the spanwise direction
    for j = 2:length(LEDiscr)
        
%         LocalRoot = @(x) ( TEPoints(j-1, :) - LEPoints(j-1, :) ) .*x + LEPoints(j-1, :);
%         LocalTip = @(x) ( TEPoints(j, :) - LEPoints(j, :) ) .*x + LEPoints(j, :);
        
        
        % Construct finite vortex extremes
        Vortices{i-1, j-1}.Root = VortexFun(LEDiscr(j-1));
        Vortices{i-1, j-1}.Tip = VortexFun(LEDiscr(j));
        
        % Construct infinite vortex extremes
        InfiniteVortices{i-1, j-1}.Root.onWing = Vortices{i-1, j-1}.Root;
        InfiniteVortices{i-1, j-1}.Root.toInfty = Vortices{i-1, j-1}.Root + [InfiniteVorticesLength, 0, 0];
        
        InfiniteVortices{i-1, j-1}.Tip.onWing = Vortices{i-1, j-1}.Tip;
        InfiniteVortices{i-1, j-1}.Tip.toInfty = Vortices{i-1, j-1}.Tip + [InfiniteVorticesLength, 0, 0];
        
        % Construct control points
        ControlPoints{i-1, j-1}.Coords = ControlPointsFun((LEDiscr(j) + LEDiscr(j-1)) / 2);
        
        % Construct inducing points
        InducedPoints{i-1, j-1}.Coords = VortexFun((LEDiscr(j) + LEDiscr(j-1)) / 2);
        

    end
end




%% Now create the symmetric portion

% Normals

yShift = config.LEPosition_Y(iCorpo);

% Going in the chordwise direction
for i = 2:length(RootDiscr)
    
    % Going in the spanwise direction
    for j = length(LEDiscr)+1:(2*length(LEDiscr)-1)
        
        %% First shift everything to make Root points to the y = 0 axis
        
        % Construct finite vortex extremes
        Vortices{i-1, j-1}.Root = Vortices{i-1, j- length(LEDiscr)}.Tip - [0 yShift 0];
        Vortices{i-1, j-1}.Tip = Vortices{i-1, j- length(LEDiscr)}.Root - [0 yShift 0];
        
        % Construct infinite vortex extremes
        InfiniteVortices{i-1, j-1}.Tip.onWing = InfiniteVortices{i-1, j- length(LEDiscr)}.Root.onWing - [0 yShift 0];
        InfiniteVortices{i-1, j-1}.Tip.toInfty = InfiniteVortices{i-1, j- length(LEDiscr)}.Root.toInfty - [0 yShift 0];
        
        InfiniteVortices{i-1, j-1}.Root.onWing = InfiniteVortices{i-1, j- length(LEDiscr)}.Tip.onWing - [0 yShift 0];
        InfiniteVortices{i-1, j-1}.Root.toInfty = InfiniteVortices{i-1, j- length(LEDiscr)}.Tip.toInfty - [0 yShift 0];
        
        % Construct control points
        ControlPoints{i-1, j-1}.Coords = ControlPoints{i-1, j- length(LEDiscr)}.Coords - [0 yShift 0];
        
        % Construct inducing points
        InducedPoints{i-1, j-1}.Coords = InducedPoints{i-1, j- length(LEDiscr)}.Coords - [0 yShift 0];
        
        
        %% Then create symmetric counterparts
        
        % Vortex extremes
        Vortices{i-1, j-1}.Root(2) = -Vortices{i-1, j-1}.Root(2);
        Vortices{i-1, j-1}.Tip(2) = -Vortices{i-1, j-1}.Tip(2);
        
        % Infinite vortex extremes
        InfiniteVortices{i-1, j-1}.Root.onWing(2) = -InfiniteVortices{i-1, j-1}.Root.onWing(2);
        InfiniteVortices{i-1, j-1}.Root.toInfty(2) = -InfiniteVortices{i-1, j-1}.Root.toInfty(2);
        
        InfiniteVortices{i-1, j-1}.Tip.onWing(2) = -InfiniteVortices{i-1, j-1}.Tip.onWing(2);
        InfiniteVortices{i-1, j-1}.Tip.toInfty(2) = -InfiniteVortices{i-1, j-1}.Tip.toInfty(2);
        
        % Control points
        ControlPoints{i-1, j-1}.Coords(2) = -ControlPoints{i-1, j-1}.Coords(2);
        
        % Inducing points
        InducedPoints{i-1, j-1}.Coords(2) = -InducedPoints{i-1, j-1}.Coords(2);
        
        
        
         %% Then shift everything to return back to original LE_Root_y values
        
        % Construct finite vortex extremes
        Vortices{i-1, j-1}.Root = Vortices{i-1, j-1}.Root + [0 yShift 0];
        Vortices{i-1, j-1}.Tip = Vortices{i-1, j-1}.Tip + [0 yShift 0];
        
        % Construct infinite vortex extremes
        InfiniteVortices{i-1, j-1}.Root.onWing = InfiniteVortices{i-1, j-1}.Root.onWing + [0 yShift 0];
        InfiniteVortices{i-1, j-1}.Root.toInfty = InfiniteVortices{i-1, j-1}.Root.toInfty + [0 yShift 0];
        
        InfiniteVortices{i-1, j-1}.Tip.onWing = InfiniteVortices{i-1, j-1}.Tip.onWing + [0 yShift 0];
        InfiniteVortices{i-1, j-1}.Tip.toInfty = InfiniteVortices{i-1, j-1}.Tip.toInfty + [0 yShift 0];
        
        % Construct control points
        ControlPoints{i-1, j-1}.Coords = ControlPoints{i-1, j-1}.Coords + [0 yShift 0];
        
        % Construct inducing points
        InducedPoints{i-1, j-1}.Coords = InducedPoints{i-1, j-1}.Coords + [0 yShift 0];
        
   

    end
end


% %% Perform eventual rotations
% 
% % For now only Y rotations
% for i = 2:length(RootDiscr)
%     
%     % Going in the spanwise direction
%     for j = 2:length(LEDiscr)
% 
%         % First translate everything such that LE is at (0, 0, 0)
% 
%         % Vortex extremes
%         Vortices{i-1, j-1}.Root = Vortices{i-1, j-1}.Root - [config.LEPosition_X(iCorpo) config.LEPosition_Y(iCorpo) config.LEPosition_Z(iCorpo)];
%         Vortices{i-1, j-1}.Tip = Vortices{i-1, j-1}.Tip;
%         
%         % Infinite vortex extremes
%         InfiniteVortices{i-1, j-1}.Root.onWing(2) = -InfiniteVortices{i-1, j-1}.Root.onWing(2);
%         InfiniteVortices{i-1, j-1}.Root.toInfty(2) = -InfiniteVortices{i-1, j-1}.Root.toInfty(2);
%         
%         InfiniteVortices{i-1, j-1}.Tip.onWing(2) = -InfiniteVortices{i-1, j-1}.Tip.onWing(2);
%         InfiniteVortices{i-1, j-1}.Tip.toInfty(2) = -InfiniteVortices{i-1, j-1}.Tip.toInfty(2);
%         
%         % Control points
%         ControlPoints{i-1, j-1}.Coords(2) = -ControlPoints{i-1, j-1}.Coords(2);
%         
%         % Inducing points
%         InducedPoints{i-1, j-1}.Coords(2) = -InducedPoints{i-1, j-1}.Coords(2);
% 
%     end
% end


%% Compute the normals

for i = 2:length(RootDiscr)
    LocalLE = @(x) ( TipPoints(i-1, :) - RootPoints(i-1, :) ).*x + RootPoints(i-1, :);
    LocalTE = @(x) ( TipPoints(i, :) - RootPoints(i, :) ).*x + RootPoints(i, :);
    
    % Going in the spanwise direction
    for j = 2:length(LEDiscr)
        
%         LocalRoot = @(x) ( TEPoints(j-1, :) - LEPoints(j-1, :) ) .*x + LEPoints(j-1, :);
%         LocalTip = @(x) ( TEPoints(j, :) - LEPoints(j, :) ) .*x + LEPoints(j, :);
        
        % Construct the extremes of the panel
        internalMesh{i-1, j-1}.LERoot = LocalLE(LEDiscr(j-1));
        internalMesh{i-1, j-1}.LEtip = LocalLE(LEDiscr(j));
        internalMesh{i-1, j-1}.TEtip = LocalTE(TEDiscr(j));
        internalMesh{i-1, j-1}.TERoot = LocalTE(TEDiscr(j-1));
        
        % Construct normal vector
        Normals{i-1, j-1}.Coords = cross(internalMesh{i-1, j-1}.TEtip - internalMesh{i-1, j-1}.TERoot, internalMesh{i-1, j-1}.LERoot - internalMesh{i-1, j-1}.TERoot);
        Normals{i-1, j-1}.Coords = Normals{i-1, j-1}.Coords/norm(Normals{i-1, j-1}.Coords);
        

    end
end
% Going in the chordwise direction
for i = 2:length(RootDiscr)
    
    % Going in the spanwise direction
    for j = length(LEDiscr)+1:(2*length(LEDiscr)-1)
        
        %% First shift everything to make Root points to the y = 0 axis
        
        % Construct the extremes of the panel
        internalMesh{i-1, j-1}.LERoot = internalMesh{i-1, j- length(LEDiscr)}.LERoot - [0 yShift 0];
        internalMesh{i-1, j-1}.LEtip = internalMesh{i-1, j- length(LEDiscr)}.LEtip - [0 yShift 0];
        internalMesh{i-1, j-1}.TEtip = internalMesh{i-1, j- length(LEDiscr)}.TEtip - [0 yShift 0];
        internalMesh{i-1, j-1}.TERoot = internalMesh{i-1, j- length(LEDiscr)}.TERoot - [0 yShift 0];

        
        %% Then create symmetric counterparts
        
        % Extremes of the panel
        internalMesh{i-1, j-1}.LERoot(2) = -internalMesh{i-1, j-1}.LERoot(2);
        internalMesh{i-1, j-1}.LEtip(2) = -internalMesh{i-1, j-1}.LEtip(2);
        internalMesh{i-1, j-1}.TEtip(2) = -internalMesh{i-1, j-1}.TEtip(2);
        internalMesh{i-1, j-1}.TERoot(2) = -internalMesh{i-1, j-1}.TERoot(2);
        
        
        
         %% Then shift everything to return back to original LE_Root_y values
        
        % Construct the extremes of the panel
        internalMesh{i-1, j-1}.LERoot = internalMesh{i-1, j-1}.LERoot + [0 yShift 0];
        internalMesh{i-1, j-1}.LEtip = internalMesh{i-1, j-1}.LEtip + [0 yShift 0];
        internalMesh{i-1, j-1}.TEtip = internalMesh{i-1, j-1}.TEtip + [0 yShift 0];
        internalMesh{i-1, j-1}.TERoot = internalMesh{i-1, j-1}.TERoot + [0 yShift 0];
        
        
        
        % Construct normal vector
        Normals{i-1, j-1}.Coords = cross(internalMesh{i-1, j-1}.LEtip - internalMesh{i-1, j-1}.LERoot, internalMesh{i-1, j-1}.TERoot - internalMesh{i-1, j-1}.LERoot);
        Normals{i-1, j-1}.Coords = Normals{i-1, j-1}.Coords/norm(Normals{i-1, j-1}.Coords);
        

    end
end



%% Now flip the original vectors to have tip -> root -> root -> tip

copia_internalMesh = internalMesh;
copia_Vortices = Vortices;
copia_InfiniteVortices = InfiniteVortices;
copia_ControlPoints = ControlPoints;
copia_InducedPoints = InducedPoints;
copia_Normals = Normals;
% Going in the chordwise direction
for i = 1:length(RootDiscr)-1
    for j = 1:length(LEDiscr)-1
        
        internalMesh{i, j} = copia_internalMesh{i, length(LEDiscr)-1 - (j-1)};
        Vortices{i, j} = copia_Vortices{i, length(LEDiscr)-1 - (j-1)};
        InfiniteVortices{i, j} = copia_InfiniteVortices{i, length(LEDiscr)-1 - (j-1)};
        ControlPoints{i, j} = copia_ControlPoints{i, length(LEDiscr)-1 - (j-1)};
        InducedPoints{i, j} = copia_InducedPoints{i, length(LEDiscr)-1 - (j-1)};
        Normals{i, j} = copia_Normals{i, length(LEDiscr)-1 - (j-1)};
    
    end
end

WingExtremes2Export = cell(3, 1);
WingExtremes2Export{1}.LE = internalMesh{1, end}.LEtip;
WingExtremes2Export{1}.TE = internalMesh{end, end}.TEtip;
WingExtremes2Export{2}.LE = WingExtremes.RootLE;
WingExtremes2Export{2}.TE = WingExtremes.RootTE;
WingExtremes2Export{3}.LE = WingExtremes.TipLE;
WingExtremes2Export{3}.TE = WingExtremes.TipTE;
