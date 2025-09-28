%[text] # Analyze Expression Data
% Standard MATLAB initialization - ensure reproducibility
clear; clc; clf;
rng default;
warning('off', 'MATLAB:print:ContentTypeImageSuggested');
ver % display version information %[output:01e7e2ea]
% MATLAB Code Conventions
%  - local variable and function names begin with an uppercase letter to
%    distinquish them from those built into MATLAB
%  - functions use named parameters instead of positional parameters
%%
%[text] ## Structure the Raw Data
Data = ReadData(FileName="example.csv");
MeanData = CalculateMeanExpression(Data=Data);
%%
%[text] ## Verify Peak Time
PeakTime = FindTimeOfPeakExpression(MeanData=MeanData);
figure; PlotTimeCourse(MeanData=MeanData, PeakTime=PeakTime, TileParam="Signal1");
figure; Figure1 = PlotTimeCourse(MeanData=MeanData, PeakTime=PeakTime, TileParam="Signal2");
SaveFigure(Figure=Figure1, FigureName="Figure1");
%%
%[text] ## Fit Models
% Use only the peak expression time for the model
PeakData = Data(Data.Time == PeakTime, :);
PeakData.Time = [];
%%
%[text] ### Single Signal Models
%[text] $E(S) = \\alpha\_0 + \\alpha\_i \\frac{\[S\_i\]}{\[S\_i\] + K\_i}\n$
% For simple estimates, we only need the mean value of expression across
% all replicates.
MeanPeakData = renamevars( ...
    groupsummary(PeakData, ["Signal1", "Signal2"], "mean", "Expression"), ...
    "mean_Expression", "Expression" ...
);

% Before fitting the model for real, calculate simple estimates of the
% model parameters. These estimates can serve as initial guesses for the
% actual model fitting algorithm, and they can provide a quick check on
% its results.
EstS1 = EstimateSingleSignal(MeanPeakData=MeanPeakData, Signal="Signal1");
EstS2 = EstimateSingleSignal(MeanPeakData=MeanPeakData, Signal="Signal2");

% With the estimates as initial guesses, fit models for each signal.
S1Params = FitSingleSignal(PeakData=PeakData, Signal="Signal1", InitialGuesses=EstS1);
S2Params = FitSingleSignal(PeakData=PeakData, Signal="Signal2", InitialGuesses=EstS2);

figure; PlotSingleSignal(PeakData=PeakData, Signal="Signal1", Params=S1Params);
Figure2 = figure; PlotSingleSignal(PeakData=PeakData, Signal="Signal2", Params=S2Params);
SaveFigure(Figure=Figure2, FigureName="Figure2");
%%
%[text] ### Multi-Signal Model
%[text] $E(S) = \\alpha\_0 + \\alpha\_1 \\frac{\[S\_1\]}{\[S\_1\] + K\_1} + \\alpha\_2 \\frac{\[S\_2\]}{\[S\_2\] + K\_2} + \\alpha\_{1,2} \\frac{\[S\_1\]\[S\_2\]}{(\[S\_1\] \[S\_2\] + {K\_{1,2}^2})}$
% As with single signal models, first calculate simple estimates for the
% parameter. Note that we use the computed model parameters for the single
% signal models, not the initial estimates, since the computed model is
% assumed to be more accurate.
EstModel = EstimateMultiSignal(MeanPeakData=MeanPeakData, EstS1=S1Params, EstS2=S2Params);

% Check estimates
Figure3 = figure; PlotMultiSignal( ...
    Signal1Values=unique(MeanPeakData.Signal1), ...
    Signal2Values=unique(MeanPeakData.Signal2), ...
    Params=EstModel, ...
    Observations=PeakData ...
);
SaveFigure(Figure=Figure3, FigureName="Figure3", Complex=true);

Params = FitMultiSignal(PeakData=PeakData, InitialGuesses=EstModel);

figure; PlotMultiSignal( ...
    Signal1Values=unique(MeanPeakData.Signal1), ...
    Signal2Values=unique(MeanPeakData.Signal2), ...
    Params=Params, ...
    Observations=PeakData ...
);
%%
%[text] ### Dynamical Models
%[text] $\\frac{\\mathrm{d}S\_i}{\\mathrm{dt}} = \nc\_i E\_i(\\mathbf{S})\\cdot N - \n \\delta\_i \\cdot S\_i  -  \n m \\cdot S\_i$
Synthase1 = EstimateFullModel(FileName="example.csv");
Synthase2 = EstimateFullModel(FileName="example2.csv");

D1_D2 = 1.7;
C_D2 = [1.2652e-05, 2.5366e-05];
NRange = 0 : 0.005 : 1.0;
MRange = 0 : 0.025 : 5.0;

Columns = [
    ["Density",      "double"]; ... 
    ["MassTransfer", "double"]; ...
    ["Signal1Star",  "double"]; ...
    ["Signal2Star",  "double"]  ...
];
NMTable = table( ...
    'Size', [length(MRange) * length(NRange), height(Columns)], ...
    'VariableNames', Columns(:,1), ...
    'VariableTypes', Columns(:,2) ...
);

idx = 1;
for N = NRange
    S0 = [1, 1];
    for M = MRange
        SStar = Equilibrium(N=N, d1_d2=D1_D2, m_d2=M, c_d2=C_D2, Synthase1=Synthase1, Synthase2=Synthase2, S0=S0);
        NMTable.Density(idx) = N;
        NMTable.MassTransfer(idx) = M;
        NMTable.Signal1Star(idx) = SStar(1);
        NMTable.Signal2Star(idx) = SStar(2);
        idx = idx + 1;
        S0 = SStar;
    end
end

Figure4 = figure; PlotSignalHeatmap(Data=NMTable, Signal="Signal1");
SaveFigure(Figure=Figure4, FigureName="Figure4");

figure; PlotSignalHeatmap(Data=NMTable, Signal="Signal2");

Effector = EstimateFullModel(FileName="effector.csv");
NMTable.Effector = Effector.ExpFn(NMTable.Signal1Star, NMTable.Signal2Star);

Figure5 = figure; PlotEffectorHeatmap(Data=NMTable);
SaveFigure(Figure=Figure5, FigureName="Figure5");
%%
%[text] ## Additional Examples
%[text] ### Extended Peak Expression Time
LongPeak = EstimateFullModel(FileName="longpeak.csv");
Figure6 = figure; PlotMaxTimeCourse( ...
    Data=LongPeak.Data, ...
    MeanData=LongPeak.MeanData, ...
    PeakTimes=[LongPeak.PeakTime - 1, LongPeak.PeakTime, LongPeak.PeakTime + 1]);
SaveFigure(Figure=Figure6, FigureName="Figure6");
%%
%[text] ### Multi-Signal Response
Additive = EstimateFullModel(FileName="additive.csv");
Figure7 = tiledlayout(2, 2);

nexttile; CompareMultiSignal( ...
    Signal1Values=unique(Synthase1.MeanPeakData.Signal1), ...
    Signal2Values=unique(Synthase1.MeanPeakData.Signal2), ...
    Params=Synthase1.Params, ...
    Observations=Synthase1.PeakData, ...
    Tiled=2, ...
    Title="Gene A Response" ...
);
text(-0.15, 1.05, "A", FontWeight="bold", Units="normalized");
nexttile; CompareMultiSignal( ...
    Signal1Values=unique(Additive.MeanPeakData.Signal1), ...
    Signal2Values=unique(Additive.MeanPeakData.Signal2), ...
    Params=Additive.Params, ...
    Observations=Additive.PeakData, ...
    Tiled=2, ...
    Title="Gene B Response" ...
);
text(-0.15, 1.05, "B", FontWeight="bold", Units="normalized");
nexttile; CompareMultiSignal( ...
    Signal1Values=unique(Effector.MeanPeakData.Signal1), ...
    Signal2Values=unique(Effector.MeanPeakData.Signal2), ...
    Params=Effector.Params, ...
    Observations=Effector.PeakData, ...
    Tiled=2, ...
    Title="Gene C Response" ...
);
text(-0.15, 1.05, "C", FontWeight="bold", Units="normalized");
nexttile; PlotMultiSignal( ...
    Signal1Values=unique(Effector.MeanPeakData.Signal1), ...
    Signal2Values=unique(Effector.MeanPeakData.Signal2), ...
    Params=Effector.Params, ...
    Observations=Effector.PeakData, ...
    Tiled=2, ...
    Title="Gene C Model" ...
);
text(-0.15, 1.05, "D", FontWeight="bold", Units="normalized");
% Cleanup some of MATLAB's errors
Figure7.Children(1).Position = [0.8, 0.42, 0.1077, 0.0608];
Figure7.Children(3).Position = [0.3, 0.42, 0.1491, 0.0608];
Figure7.Children(5).Position = [0.75, 0.9, 0.1491, 0.0608];
Figure7.Children(7).Position = [0.3, 0.9, 0.1491, 0.0608];
Figure7.Children(2).FontSize = 5;
Figure7.Children(4).FontSize = 5;
Figure7.Children(6).FontSize = 5;
Figure7.Children(8).FontSize = 5;

SaveFigure(Figure=Figure7, FigureName="Figure7", Complex=true);
%%
%[text] ### Diverse Signals
NMTable2 = table( ...
    'Size', [length(MRange) * length(NRange), height(Columns)], ...
    'VariableNames', Columns(:,1), ...
    'VariableTypes', Columns(:,2) ...
);

idx = 1;
for N = NRange
    S0 = [1, 1];
    for M = MRange
        SStar = Equilibrium(N=N, d1_d2=4, m_d2=M, c_d2=C_D2, Synthase1=Synthase1, Synthase2=Synthase2, S0=S0);
        NMTable2.Density(idx) = N;
        NMTable2.MassTransfer(idx) = M;
        NMTable2.Signal1Star(idx) = SStar(1);
        NMTable2.Signal2Star(idx) = SStar(2);
        idx = idx + 1;
        S0 = SStar;
    end
end

NMTable2.Effector = Effector.ExpFn(NMTable2.Signal1Star, NMTable2.Signal2Star);

Figure8 = figure;
subplot(2,1,1);
PlotEffectorHeatmap(Data=NMTable, Title="Example A", SignalContours=true);
text(-0.15, 1.05, "A", FontWeight="bold", Units="normalized");

subplot(2,1,2);
PlotEffectorHeatmap(Data=NMTable2, Title="Example B", SignalContours=true);
text(-0.15, 1.05, "B", FontWeight="bold", Units="normalized");

Figure8.Position(4) = Figure8.Position(4) * 2;
SaveFigure(Figure=Figure8, FigureName="Figure8");
%%
%[text] ## Local Functions
%[text] ### Data Manipulation
function Data = ReadData(Args)
    arguments
        Args.FileName string
    end

    Data = readtable("../Example Data/" + Args.FileName);

    % Ensure all required columns are present
    MissingColumns = setdiff( ...
        {'Time', 'Signal1', 'Signal2', 'Replicate'}, ...
        Data.Properties.VariableNames ...
    );
    if (~isempty(MissingColumns))
        error("File is missing the following columns: %s", strjoin(MissingColumns, ", "))
    end

    % Sort rows into specific order for ease of processing
    Data = sortrows(Data, {'Time', 'Signal1', 'Signal2'});
end

function MeanData = CalculateMeanExpression(Args)
    arguments
        Args.Data table
    end

    % Preserve the Time, Signal1, and Signal2 values, but replace
    % Expression with its mean (across replicates).
    MeanData = renamevars( ...
        groupsummary(Args.Data, ["Time", "Signal1", "Signal2"], "mean", "Expression"), ...
        "mean_Expression", "Expression" ...
    );
    % Remove the additional column added by group summary
    MeanData.GroupCount = [];
end

function [PeakTime, PeakTimes] = FindTimeOfPeakExpression(Args)
    arguments
        Args.MeanData table
    end

    Data = Args.MeanData;

    % Split the data into groups for each combination of signal values.
    % SignalGroups will indicate the group number of each row in the data
    % table, and UniqueS1 and UniqueS2 will hold the unique values of each
    % signal.
    [SignalGroups, UniqueS1, UniqueS2] = findgroups(Data.Signal1, Data.Signal2);

    % Find the row index that corresponds to the maximum expression within
    % each group.
    MaxExpressionIdx = splitapply( ...
        @(idx) idx(find(Data.Expression(idx) == max(Data.Expression(idx)), 1)), ...
        (1:height(Data))', ...
        SignalGroups ...
    );

    % Create a new table with a row for each combination of signal values
    % and the time at which peak expression occurs for that combination
    PeakTimes = table(UniqueS1, UniqueS2, Data.Time(MaxExpressionIdx), ...
                    'VariableNames', {'Signal1', 'Signal2', 'Time'});

    % The most common peak time value in the table is the single value for
    % time of peak expression.
    PeakTime = mode(PeakTimes.Time);
end
%[text] ### Model Computations
function EstSingleSignalParams = EstimateSingleSignal(Args)
    arguments
        Args.MeanPeakData table
        Args.Signal string {mustBeMember(Args.Signal, ["Signal1", "Signal2"])}
    end

    Data = Args.MeanPeakData;
    Signal = Args.Signal;
    if Signal == "Signal1"
        NullSignal = "Signal2";
    else
        NullSignal = "Signal1";
    end

    SignalValues = unique(Data.(Signal));

    % Extract the rows we need for a single signal
    SignalData = Data(Data.(NullSignal) == 0, :);

    % A0: mean value of expression with no signal
    A0 = SignalData{SignalData.(Signal) == 0, "Expression"};

    % A: difference between extreme expression and A0; it could be
    % positive or negative
    ExtremeExpression = SignalData{SignalData.(Signal) == max(SignalValues), "Expression"};

    if (ExtremeExpression / A0 > 1.05)
        A = ExtremeExpression - A0; % A is positive
        HalfExpression = A0 + (ExtremeExpression - A0) / 2;
        SignalHalfIdx = find(SignalData.Expression > HalfExpression, 1, 'first');
    elseif (ExtremeExpression / A0 < 0.95)
        A = ExtremeExpression - A0; % A is negative
        HalfExpression = A0 - (A0 - ExtremeExpression) / 2;
        SignalHalfIdx = find(SignalData.Expression < HalfExpression, 1, 'first');
    else
        A = 0;
        HalfExpression = A0;
        SignalHalfIdx = 2;
    end

    % K: concentration that yields expression half way between A0 and
    % max/min expression; use linear interpolation between the two available
    % signal values that bracket that half-way point
    K = interp1( ...
        [SignalData{SignalHalfIdx - 1, "Expression"}, SignalData{SignalHalfIdx, "Expression"}], ...
        [SignalData{SignalHalfIdx - 1, Signal}, SignalData{SignalHalfIdx, Signal}], ...
        HalfExpression ...
    );

    EstSingleSignalParams = struct;
    EstSingleSignalParams.A0 = A0;
    EstSingleSignalParams.A = A;
    EstSingleSignalParams.K = K;
end

function Params = FitSingleSignal(Args)
    arguments
        Args.PeakData table
        Args.Signal string {mustBeMember(Args.Signal, ["Signal1", "Signal2"])}
        Args.InitialGuesses struct
    end

    Data = Args.PeakData;
    Signal = Args.Signal;
    if Signal == "Signal1"
        NullSignal = "Signal2";
    else
        NullSignal = "Signal1";
    end

    X = Data{Data.(NullSignal) == 0, Signal};
    Y = Data{Data.(NullSignal) == 0, "Expression"};
    Beta0 = [Args.InitialGuesses.A0, Args.InitialGuesses.A, Args.InitialGuesses.K];
    ModelFn = @(b, x) b(1) + b(2) * (x ./ (b(3) + x));

    % Maintain the sign of the initial guess A values, but keep A0 and
    % no less than the minimum observed value.
    LowerBounds = [min(Data.Expression), 0, 0];
    UpperBounds = [inf, inf, inf];
    if Beta0(2) < 0
       LowerBounds(2) = -inf;
       UpperBounds(2) = 0;
    end

    % Display off to supress warnings about local optimum. If the initial
    % estimates are good, these warnings are expected. If uncertain of the
    % initial estimates, change this option.
    Beta = lsqcurvefit(ModelFn, Beta0, X, Y, LowerBounds, UpperBounds, optimset('Display','off'));

    Params = struct;
    Params.A0 = Beta(1);
    Params.A = Beta(2);
    Params.K = Beta(3);
end

function EstParams = EstimateMultiSignal(Args)
    arguments
        Args.MeanPeakData table
        Args.EstS1 struct
        Args.EstS2 struct
    end

    Data = Args.MeanPeakData;
    EstS1 = Args.EstS1;
    EstS2 = Args.EstS2;

    A0 = mean([EstS1.A0, EstS2.A0]);
    A1 = EstS1.A;
    K1 = EstS1.K;
    A2 = EstS2.A;
    K2 = EstS2.K;

    % Calculate how much of the expression value is not explained by the
    % single signal models at each combination of signal values.
    Data.Multisignal = Data.Expression ...
        - A0 ...
        - A1 .* Data.Signal1 ./ (Data.Signal1 + K1) ...
        - A2 .* Data.Signal2 ./ (Data.Signal2 + K2);

    Signal1Values = unique(Data.Signal1);
    Signal2Values = unique(Data.Signal2);
    ExtremeExpression = Data{Data.Signal1 == max(Signal1Values) & Data.Signal2 == max(Signal2Values), "Expression"};

    SingleSignalExpression = A0 + A1 + A2;

    % Guess for A12 is how much additional expression is needed to reach
    % the maximum or minimum value
    if (ExtremeExpression / SingleSignalExpression > 1.05)
        A12 = ExtremeExpression - SingleSignalExpression; % A12 is positive
        Data.Multisignal = max(Data.Multisignal, 0);
    elseif (ExtremeExpression / SingleSignalExpression < 0.95)
        A12 = ExtremeExpression - SingleSignalExpression; % A12 is negative
        Data.Multisignal = min(Data.Multisignal, 0);
    else
        A12 = 0;
        K12 = sqrt(realmin);
    end

    if A12 ~= 0
        % If there is a nonzero A12, use a simple nonlinear least squares fit
        % to estimate K12. By assuming the estimates for all other parameters
        % are accurate, the model only needs to fit the single parameter K12.
        X = [Data.Signal1, Data.Signal2];
        Y = Data.Multisignal;
        ModelFn = @(k12,x) A12 * x(:,1) .* x(:,2) ./ (x(:,1) .* x(:,2) + k12^2);
        K12 = nlinfit(X, Y, ModelFn, 0.001);
    end
 
    EstParams = struct;
    EstParams.A0 = A0;
    EstParams.A1 = A1;
    EstParams.K1 = K1;
    EstParams.A2 = A2;
    EstParams.K2 = K2;
    EstParams.A12 = A12;
    EstParams.K12 = K12;
end

function Params = FitMultiSignal(Args)
    arguments
        Args.PeakData table
        Args.InitialGuesses struct
    end

    Data = Args.PeakData;

    X = [Data.Signal1, Data.Signal2];
    Y = Data.Expression;
    Beta0 = [ ...
        Args.InitialGuesses.A0, ...
        Args.InitialGuesses.A1, ...
        Args.InitialGuesses.K1, ...
        Args.InitialGuesses.A2, ...
        Args.InitialGuesses.K2, ...
        Args.InitialGuesses.A12, ...
        Args.InitialGuesses.K12 ...
    ];
    ModelFn = @(b,x) ...
        b(1) + ...
        b(2) * x(:,1) ./ (x(:,1) + b(3)) + ...
        b(4) * x(:,2) ./ (x(:,2) + b(5)) + ...
        b(6) * x(:,1) .* x(:,2) ./ (x(:,1) .* x(:,2) + b(7) * b(7));

    % Maintain the sign of the initial guess A values, but keep A0 and
    % no less than the minimum observed value.
    LowerBounds = [min(Data.Expression), 0, 0, 0, 0, 0, 0];
    UpperBounds = [inf, inf, inf, inf, inf, inf, inf];
    for idx = [2, 4, 6]
        if Beta0(idx) < 0
           LowerBounds(idx) = -inf;
           UpperBounds(idx) = 0;
        end
    end

    % Display off to supress warnings about local optimum. If the initial
    % estimates are good, these warnings are expected. If uncertain of the
    % initial estimates, change this option.
    Beta = lsqcurvefit(ModelFn, Beta0, X, Y, LowerBounds, UpperBounds, optimset('Display','off'));

    Params = struct;
    Params.A0 = Beta(1);
    Params.A1 = Beta(2);
    Params.K1 = Beta(3);
    Params.A2 = Beta(4);
    Params.K2 = Beta(5);
    Params.A12 = Beta(6);
    Params.K12 = Beta(7);
end

function Model = EstimateFullModel(Args)
    arguments
        Args.FileName string
    end
    
    Model = struct;
    Model.Data = ReadData(FileName=Args.FileName);
    Model.MeanData = CalculateMeanExpression(Data=Model.Data);
    Model.PeakTime = FindTimeOfPeakExpression(MeanData=Model.MeanData);
    Model.PeakData = Model.Data(Model.Data.Time == Model.PeakTime, :);
    Model.PeakData.Time = [];
    Model.MeanPeakData = renamevars( ...
        groupsummary(Model.PeakData, ["Signal1", "Signal2"], "mean", "Expression"), ...
        "mean_Expression", "Expression" ...
    );
    Model.EstS1 = EstimateSingleSignal(MeanPeakData=Model.MeanPeakData, Signal="Signal1");
    Model.EstS2 = EstimateSingleSignal(MeanPeakData=Model.MeanPeakData, Signal="Signal2");
    Model.S1Params = FitSingleSignal(PeakData=Model.PeakData, Signal="Signal1", InitialGuesses=Model.EstS1);
    Model.S2Params = FitSingleSignal(PeakData=Model.PeakData, Signal="Signal2", InitialGuesses=Model.EstS2);

    Model.EstModel = EstimateMultiSignal(MeanPeakData=Model.MeanPeakData, EstS1=Model.S1Params, EstS2=Model.S2Params);
    Model.Params = FitMultiSignal(PeakData=Model.PeakData, InitialGuesses=Model.EstModel);

    Model.ExpFn = @(s1, s2) Model.Params.A0 + ...
        Model.Params.A1 .* s1 ./ (s1 + Model.Params.K1) + ...
        Model.Params.A2 .* s2 ./ (s2 + Model.Params.K2) + ...
        Model.Params.A12 .* s1 .* s2 ./ (s1 .* s2 + Model.Params.K12 .* Model.Params.K12);

end
%%
%[text] ### Dynamical Models
function DeltaS_d2 = dS_dt(Args)
% dS_dt() calculates the dynamics of the system of two
% signals as the change in concentration for each signal
% over time interval of C4 decay rate.
    arguments
        Args.S (2, 1) double;      % [Signal1, Signal2] concentration
        Args.N double;             % population size
        Args.d1_d2 double;         % Signal1 / Signal2 decay rate
        Args.m_d2 double = 0;      % mass transfer / Signal2 decay rate
        Args.c_d2 (2, 1) double;   % [c1_d2, c2_d2] constants
        Args.Synthase1 struct;     % Synthase1 parameters
        Args.Synthase2 struct;     % Synthase2 parameters
    end

    S = Args.S;
    N = Args.N;
    d1_d2 = Args.d1_d2;
    m_d2 = Args.m_d2;
    c_d2 = Args.c_d2;
    Synthase1 = Args.Synthase1;
    Synthase2 = Args.Synthase2;

    Synthase1Exp = Synthase1.ExpFn(S(1), S(2));
    Synthase2Exp = Synthase2.ExpFn(S(1), S(2));

    DeltaS_d2 = [ ...
        c_d2(1) * Synthase1Exp * N - S(1) * (d1_d2 + m_d2);
        c_d2(2) * Synthase2Exp * N - S(2) * (1.0 + m_d2)
    ];

end

function Sstar = Equilibrium(Args)
% Equilibrium() returns equilibrium concentrations of signals [Signal1, Signal2]
    arguments (Input)
        Args.N double;                   % population size
        Args.d1_d2 double;               % Signal1 / Signal2 decay rate
        Args.m_d2 double = 0;            % mass transfer / Signal2 decay rate
        Args.c_d2 (2, 1) double;         % [c1_d2, c2_d2] constants 
        Args.S0 (2, 1) double = [1, 1]   % Optional starting guess; see note below
        Args.Synthase1 struct;           % Synthase1 parameters
        Args.Synthase2 struct;           % Synthase2 parameters
    end

    % Note: Setting S0 to different values can help prevent the
    % optimization algorithm from converging on alternate solutions in
    % which the concentrations are negative. A possible enhancement would
    % be to use non-linear optimization instead of a fsolve so that
    % constraints can be defined, but this is simpler for now. As a
    % safety check to make sure that negative values are not returned,
    % we use output validation.
    arguments (Output)
        Sstar (2, 1) double {mustBeNonnegative}
    end

    N = Args.N;
    d1_d2 = Args.d1_d2;
    m_d2 = Args.m_d2;
    c_d2 = Args.c_d2;
    S0 = Args.S0;
    Synthase1 = Args.Synthase1;
    Synthase2 = Args.Synthase2;

    SolveOptions = optimoptions('fsolve', 'Display', 'none');
    dS = @(S) dS_dt(S=S, N=N, d1_d2=d1_d2, m_d2=m_d2, c_d2=c_d2, Synthase1=Synthase1, Synthase2=Synthase2);
    Sstar = fsolve(dS, S0, SolveOptions);
end
%%
%[text] ### Plotting
function TimeCoursePlot = PlotTimeCourse(Args)
    arguments
        Args.MeanData table
        Args.PeakTime double
        Args.TileParam string {mustBeMember(Args.TileParam, ["Signal1", "Signal2"])}
    end

    Data = Args.MeanData;
    PeakTime = Args.PeakTime;
    MaxTime = max(Data.Time);

    Signal1Values = unique(Data.Signal1)';
    Signal2Values = unique(Data.Signal2)';

    if Args.TileParam == "Signal1"
        TileValues = Signal2Values;
        LineParam = "Signal2";
        LineValues = Signal1Values;
    else
        TileValues = Signal1Values;
        LineParam = "Signal1";
        LineValues = Signal2Values;
    end

    TileColumns = 3;
    TimeCoursePlot = tiledlayout(ceil(length(TileValues)/TileColumns), TileColumns);
    XLimits = [0, MaxTime];
    YLimits = [0, max(Data.Expression)];
    X = unique(Data.Time);

    for TileIdx = 1 : length(TileValues)
        TileValue = TileValues(TileIdx);
        nexttile;
        hold on

        % Use a single-value bar plot to indicate the peak time value. It's
        % plotted first so that the line plots will appear "on top of" the
        % bar. Note that the bar plot should be excluded from any legend
        % since it's not really actual data.
        bar(PeakTime, YLimits(2) * 1.1, 1, FaceColor="black", EdgeA=0, FaceA=0.10);

        % Reset the color orders since the bar plot did not use a standard
        % value.
        ax = gca;
        ax.ColorOrderIndex = 1;

        for LineValue = LineValues
            if Args.TileParam == "Signal1"
                y = Data{Data.Signal1 == TileValue & Data.Signal2 == LineValue, "Expression"};
            else
                y = Data{Data.Signal1 == LineValue & Data.Signal2 == TileValue, "Expression"};
            end
            plot(X, y, LineWidth=1.5);
            yticks({});
        end

        xlim(XLimits);
        ylim(YLimits);

        if (mod(TileIdx, TileColumns) == 1)
            % ylabel("Expression");
        end
        if ((length(TileValues) - TileIdx) < TileColumns)
            xlabel("Time");
        end

        text(1, 1, Args.TileParam + ": " + num2str(TileValue), ...
            FontSize=7, Units="normalized", HorizontalAlignment="right", VerticalAlignment="bottom");

        hold off
    end

    % Add a legend that shows the signal values corresponding to each line.
    % Only one legend is needed for the entire layout since those signal
    % values are the same for each plot. Note that whe insert an empty
    % string first to exclude the bar chart used to highlight the time of
    % peak expression.
    LegendValues = [{''}, arrayfun(@num2str, LineValues, UniformOutput=false)];
    LTimeCourse = legend(LegendValues, Location="eastoutside", Direction="reverse");
    LTimeCourse.Layout.Tile = 'east';
    LTimeCourse.Title.String = LineParam;

    % Use a two-line title to force some separation from the plot itself.
    title(TimeCoursePlot, {'Mean Expression', ''});
end

function PlotMaxTimeCourse(Args)
    arguments
        Args.Data table
        Args.MeanData table
        Args.PeakTimes double
    end

    Data = Args.Data;
    MeanData = Args.MeanData;
    PeakTimes = Args.PeakTimes;

    MaxSignal1 = max(MeanData.Signal1);
    MaxSignal2 = max(MeanData.Signal2);
    PeakValues = ones(size(PeakTimes)) * max(Data.Expression) * 1.2;

    % Use a bar chart to show the peak time range in a muted color. After
    % plotting, reset the color order since the bar plot isn't using a
    % color from that range
    bar(PeakTimes, PeakValues, 1, FaceColor="black", EdgeA=0, FaceA=0.10);
    hold on;
    ax = gca;
    ax.ColorOrderIndex = 1;

    % Plot the mean expression as a line
    X = MeanData{MeanData.Signal1 == MaxSignal1 & MeanData.Signal2 == MaxSignal2, "Time"};
    Y = MeanData{MeanData.Signal1 == MaxSignal1 & MeanData.Signal2 == MaxSignal2, "Expression"};

    plot(X, Y, LineWidth=2);
    yticks({});
    xticks(X);
    xlabel("Time");
    ylabel("Expression");

    % Plot the observations
    Observations = Data(Data.Signal1 == MaxSignal1 & Data.Signal2 == MaxSignal2, :);
    plot(Observations.Time, Observations.Expression, 'o', MarkerSize=8, LineWidth=1.5);
    if mean(PeakTimes) > max(Data.Time)/2
        Location = "northwest";
    else
        Location = "northeast";
    end
    legend({'Peak Range', 'Mean Value', 'Observations'}, Location=Location)
end

function PlotSingleSignal(Args)
    arguments
        Args.PeakData table
        Args.Signal string {mustBeMember(Args.Signal, ["Signal1", "Signal2"])}
        Args.Params struct
    end

    Data = Args.PeakData;
    Signal = Args.Signal;
    if Signal == "Signal1"
        NullSignal = "Signal2";
    else
        NullSignal = "Signal1";
    end
    Params = Args.Params;

    ModelFn = @(x) Params.A0 + Params.A * (x / (x + Params.K));

    Normalize = mean(Data{Data.Signal1 == 0 & Data.Signal2 == 0, "Expression"});

    X = Data{Data.(NullSignal) == 0, Signal};
    Y = Data{Data.(NullSignal) == 0, "Expression"};
    Y = Y / Normalize;
    YPred = arrayfun(ModelFn, X) ./ Normalize;
    semilogx(X, YPred, LineWidth=2);
    colororder("glow");
    hold on;
    semilogx(X, Y, "o", MarkerSize=8, LineWidth=1.5);
    hold off;
    ylabel("N-Fold Change in Expression");
    xlabel(Signal);
    title("Single Signal Model");
    legend("Predicted", "Observed", Location="northwest");
end

function PlotMultiSignal(Args)
    arguments
        Args.Signal1Values
        Args.Signal2Values
        Args.Params struct
        Args.Observations table
        Args.Tiled double = 1
        Args.Title string = "Multi-Signal Model"
    end

    A0 = Args.Params.A0;
    A1 = Args.Params.A1;
    K1 = Args.Params.K1;
    A2 = Args.Params.A2;
    K2 = Args.Params.K2;
    A12 = Args.Params.A12;
    K12 = Args.Params.K12;

    Scale = 1/Args.Tiled;

    Observations = Args.Observations;

    Normalize = mean(Observations{Observations.Signal1 == 0 & Observations.Signal2 == 0, "Expression"});
    Observations.FoldChange = Observations.Expression / Normalize;

    Predict = @(s1,s2) ...
        A0 ...
        + A1 .* s1 ./ (s1 + K1) ...
        + A2 .* s2 ./ (s2 + K2) ...
        + A12 .* s1 .* s2 ./ (s1 .* s2 + K12^2);

    S1 = Args.Signal1Values;
    S2 = Args.Signal2Values;
    [X, Y] = meshgrid(S1, S2);
    Z = reshape(Predict(X,Y) ./ Normalize, length(S1), length(S2));

    surf(X, Y, Z, FaceColor="none", EdgeColor="black");
    hold on

    colororder("glow");

    plot3(Observations, "Signal1", "Signal2", "FoldChange", ...
        LineStyle="none", Marker=".", MarkerSize=15*Scale);
    set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
    hold off;

    xlabel("Signal 1");
    ylabel("Signal 2");
    zlabel("N-Fold Change in Expression");
    title(Args.Title);
    Legend = legend("Predicted", "Observed", Location="northeast", FontSize=10*Scale);
    if (Scale ~= 1)
        ax = gca;
        ax.TitleHorizontalAlignment = "left";
        ax.FontSize = ax.FontSize * Scale;
        Legend.ItemTokenSize = Legend.ItemTokenSize .* [Scale; 1];
    end
end

function CompareMultiSignal(Args)
    arguments
        Args.Signal1Values
        Args.Signal2Values
        Args.Params struct
        Args.Observations table
        Args.Tiled double = 1
        Args.Title string = "Multi-Signal Response"
    end

    A0 = Args.Params.A0;
    A1 = Args.Params.A1;
    K1 = Args.Params.K1;
    A2 = Args.Params.A2;
    K2 = Args.Params.K2;

    Scale = 1/Args.Tiled;

    Observations = Args.Observations;

    Normalize = mean(Observations{Observations.Signal1 == 0 & Observations.Signal2 == 0, "Expression"});
    MaxExpression = max(Observations.Expression);
    
    PredictS1 = @(s1) A0 + A1 .* s1 ./ (s1 + K1);
    PredictS2 = @(s2) A0 + A2 .* s2 ./ (s2 + K2);

    S1 = Args.Signal1Values;
    S2 = Args.Signal2Values;

    [X, Y] = meshgrid(S1, S2);
    Z = zeros(size(X));

    % Surface plot (actually, just a horizontal plane, i.e. constant value
    % for Z) to show the sum of single signals
    colororder("glow");
    Colors = colororder;
    Color = Colors(1, :);
    surf(X, Y, (Z + A0 + A1 + A2)/Normalize, FaceColor=Color, FaceAlpha=0.5, EdgeColor="none");
    hold on

    % Line plots in color to emphasize the single signal predictions that
    % make up the plane above
    plot3(S1, S2 * 0, PredictS1(S1)/Normalize, LineWidth=3.5*Scale, Color=Color);
    plot3(S1 * 0, S2, PredictS2(S2)/Normalize, LineWidth=3.5*Scale, Color=Color);

    % Second plane showing maximum expression
    Color2 = Colors(2, :);
    surf(X, Y, (Z + MaxExpression)/Normalize, FaceColor=Color2, FaceAlpha=0.75, EdgeColor="none");

    hold off;

    xlabel("Signal 1");
    ylabel("Signal 2");
    zlabel("N-Fold Change in Expression");
    title(Args.Title);
    Legend = legend( ...
        '\Sigma Single Signals', '', '', 'Max Expression', ...
        Location="northeast", FontSize=10*Scale);
    Legend.Direction = "reverse";
    Legend.Position(2) = 0.875;
    ax = gca;
    ax.TitleHorizontalAlignment = "left";
    if (Scale ~= 1)
        ax.FontSize = ax.FontSize * Scale;
        Legend.ItemTokenSize = Legend.ItemTokenSize .* [Scale; 1];
    end
end

function PlotSignalHeatmap(Args)
    arguments
        Args.Data table
        Args.Signal string {mustBeMember(Args.Signal, ["Signal1", "Signal2"])}
    end

    Data = Args.Data;
    
    MRange = sort(unique(Data.MassTransfer));
    XLabels = string(MRange);
    XLabels(mod(MRange, 1) ~= 0) = "";
    NRange = sort(unique(Data.Density));
    YLabels = string(NRange);
    YLabels(mod(NRange, 0.2) ~= 0) = "";

    ColorVariable = Args.Signal + "Star";
    Title = Args.Signal + " Equilibrium";

    HMap = heatmap(Data, "MassTransfer", "Density", ...
        ColorVariable=ColorVariable, ...
        XDisplayLabels=XLabels, YDisplayLabels=YLabels);
    HMap.YDisplayData=flip(HMap.YDisplayData);
    grid off;

    title(Title);
    xlabel("Mass Transfer Rate");
    ylabel("Population Density");

    originalUnits = HMap.Units;  % save original units (probaly normalized)
    HMap.Units = 'centimeters';  % any unit that will result in squares
    % Get number of rows & columns
    sz = size(HMap.ColorData); 
    % Change axis size & position;
    originalPos = HMap.Position; 
    % make axes square (not the table cells, just the axes)
    HMap.Position(3:4) = min(HMap.Position(3:4))*[1,1]; 
    if sz(1)>sz(2)
        % make the axis size more narrow and re-center
        HMap.Position(3) = HMap.Position(3)*(sz(2)/sz(1)); 
        HMap.Position(1) = (originalPos(1)+originalPos(3)/2)-(HMap.Position(3)/2);
    else
        % make the axis size shorter and re-center
        HMap.Position(4) = HMap.Position(4)*(sz(1)/sz(2));
        HMap.Position(2) = (originalPos(2)+originalPos(4)/2)-(HMap.Position(4)/2);
    end
    % Return axis to original units
    HMap.Units = originalUnits;

    % Code below is not working with current release.
    HMapAxes=HMap.NodeChildren(3);
    HMapAxes.XAxis.TickLabelRotationMode = 'manual';
    HMapAxes.XAxis.TickLabelRotation = 0;
end

function PlotEffectorHeatmap(Args)
    arguments
        Args.Data table
        Args.SignalContours logical = false
        Args.Title string = ""
    end

    Data = Args.Data;

    if Args.Title == ""
        Title = "N-Fold Change in Expression";
    else
        Title = Args.Title;
    end
    
    MRange = sort(unique(Data.MassTransfer));
    NRange = sort(unique(Data.Density));
    MinExpression = min(Data.Effector);

    [X, Y] = meshgrid(MRange, NRange);

    Z = griddata(Data.MassTransfer, Data.Density, Data.Effector, X, Y);

    pcolor(X, Y, Z / MinExpression);
    colormap("sky");
    axis square;
    shading interp;
    Colorbar = colorbar;

    if Args.Title ~= ""
        Colorbar.Label.String = "N-Fold Change in Expression";
    end

    hold on;

    contour(X, Y, Z / MinExpression, 1, LineWidth=2, LineColor="white");

    if Args.SignalContours
        Colors = orderedcolors("gem");
     
        Z1 = griddata(Data.MassTransfer, Data.Density, Data.Signal1Star, X, Y);
        contour(X, Y, Z1, 1, LineWidth=2, LineColor=Colors(2,:));
    
        Z2 = griddata(Data.MassTransfer, Data.Density, Data.Signal2Star, X, Y);
        contour(X, Y, Z2, 1, LineWidth=2, LineColor=Colors(3,:));
    end

    hold off;

    title(Title);
    xlabel("Mass Transfer Rate");
    ylabel("Population Density");

    if Args.SignalContours
        Legend = legend(["", "", "Signal 1", "Signal 2"], Location="southeast");
        title(Legend, ["50% Maximum", "Concentration"], FontWeight="normal");
    end

end

function SaveFigure(Args)
    arguments
        Args.Figure
        Args.FigureName string
        Args.Complex logical = false;
    end

    Figure = Args.Figure;
    FigureName = Args.FigureName;

    if (Args.Complex)
        % With complex figures, the current version of MATLAB (2025a update
        % 1) goes nuts, generating SVG files as large as a GigaByte with 
        % hundreds of thousands of duplicate elements. For those figures, 
        % fall back to a PNG format, but wrap it inside an SVG file so that
        % we can use high resolution PNG and specify a manageable intrinsic
        % size for markdown.
        exportgraphics( ...
            Figure, ...
            "../Figures/" + FigureName + ".png", ...
            Resolution=600 ...
        );

        % Unfortunately, we can't simple reference the resulting PNG file
        % inside of an SVG file. For security reasons, browsers won't
        % display an externally referenced <image> element. To work around
        % this, convert the external PNG file into an internal data URI
        % using base64 encoding. To do that, we have to read it as a binary
        % file and extract the raw data.
        PNGFile = fopen("../Figures/" + FigureName + ".png", 'rb');
        PNGData = fread(PNGFile, '*uint8');
        fclose(PNGFile);
        Base64String = matlab.net.base64encode(PNGData);

        % Retrieve the figure's intrinsic dimensions in pixels
        Figure.Units = "pixels";
        Width = Figure.OuterPosition(3);
        Height = Figure.OuterPosition(4);

        SVG = [
            "<?xml version='1.0' encoding='UTF-8'?>", ...
            "<svg width='" + Width + "px' height='" + Height + "px' " ...
            +      "xmlns='http://www.w3.org/2000/svg' " ...
            +      "xmlns:xlink='http://www.w3.org/1999/xlink'  version='1.2'" ...
            + ">" , ...
            "<image href='data:image/png;base64," + Base64String + "' " ...
            +      "width='" + Width + "' height='" + Height + "' />", ...
            "</svg>" ...
        ];
        writelines(SVG, "../Figures/" + FigureName + ".svg");
    else
        exportgraphics( ...
            Figure, ...
            "../Figures/" + FigureName + ".svg", ...
            ContentType="vector" ...
        );
    end
    exportgraphics( ...
        Figure, ...
        "../Submission/Figures/" + FigureName + ".eps", ...
        ContentType="vector" ...
    );
end


%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
%[output:01e7e2ea]
%   data: {"dataType":"text","outputData":{"text":"-----------------------------------------------------------------------------------------\nMATLAB Version: 25.1.0.2973910 (R2025a) Update 1\nMATLAB License Number: 41249549\nOperating System: macOS  Version: 15.7 Build: 24G222 \nJava Version: Java 11.0.25+9-LTS with Amazon.com Inc. OpenJDK 64-Bit Server VM mixed mode\n-----------------------------------------------------------------------------------------\nMATLAB                                                Version 25.1        (R2025a)\nCurve Fitting Toolbox                                 Version 25.1        (R2025a)\nOptimization Toolbox                                  Version 25.1        (R2025a)\nStatistics and Machine Learning Toolbox               Version 25.1        (R2025a)\nSymbolic Math Toolbox                                 Version 25.1        (R2025a)\n","truncated":false}}
%---
