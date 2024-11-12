function CalPlotting(myBayesianAnalysis)
% ========================================================================= %
% DESCRIPTION:
%   Plot histograms and posterior probability density functions (PDFs) for 
%   the Bayesian calibration analysis results and visualize maximum 
%   likelihood estimate (MLE) for each model parameter.
%
% INPUT:
%   myBayesianAnalysis:     Struct containing the results of the Bayesian
%                           calibration analysis.
%
% OUTPUT:
%   None
%
% USAGE:
%   CalPlotting(myBayesianAnalysis);
%
% NOTES:
%   - This function assumes that the myBayesianAnalysis struct contains
%     necessary fields such as PriorDist, Results, and Internal for
%     extracting data related to the Bayesian calibration analysis.
%   - The function generates histograms for prior samples and plots the
%     posterior PDFs obtained from the analysis results.
% ========================================================================= %

    % Extract the number of non-constant model parameters
    nVals = myBayesianAnalysis.Internal.nNonConstModelParams;
    
    % varNames of model parameters (Need to automise!)
    varNames = {'\alpha','\beta','k_a','k_b','k_{f,a}',...
            'k_{f,b}','k_{f,c}','k_{f,d}','d','\eta',...
            'p_p', 'TI_{exp}'};

    % Create empty storage arrays
    means = zeros(nVals,1);
    MLE = zeros(nVals,1);
    
    % Create a new figure for plotting
    figure('Position', [50,50,1350,850]);
    
    % Loop through each model parameter
    for i = 1:nVals

        % Extract information for the current parameter
	    range = myBayesianAnalysis.PriorDist.Marginals(1,i);
	    postRuns = myBayesianAnalysis.Results.PostProc.PostSample(:,i,:);
	    priorRuns = myBayesianAnalysis.Results.PostProc.PriorSample(:,i);

        % Compute the posterior PDF using KDE
	    xi = linspace(range.Parameters(1), range.Parameters(2), 100);
	    f = ksdensity(reshape(postRuns,[],1), xi);

% 	    test = myBayesianAnalysis.Results.PostProc.PostSample(end,i,:);
% 	    f1 = ksdensity(reshape(priorRuns,[],1), xi);
% 	    f2 = ksdensity(reshape(test,[],1), xi);

        % Find the peak of the posterior PDF (MLE)
	    [Pmax, Imax] = max(f);
	    
        % Calculate mean and MLE values (Compare with uqlab outputs)
	    means(i) = mean(postRuns, 'all');
	    MLE(i) = xi(Imax);
    
	    % Plotting
        subplot(3,4,i)
	    histogram(priorRuns, 50, 'Normalization', 'pdf', ...
            'Edgecolor', 'none', 'Displayname', 'Prior Samples'); hold on;
	    area(xi, f, 'FaceAlpha', 0.3, 'Displayname', 'Posterior PDF'); hold on;
 	    % histogram(postRuns, 50, 'Normalization', 'pdf', 'FaceAlpha', 0.3, 'Edgecolor', 'none', 'Displayname', 'Posterior PDF', 'HandleVisibility','off'); hold on;
	    % area(xi, f1, 'FaceAlpha', 0.3); hold on;
	    % plot(xi, f, 'k-', 'Displayname', ''); hold on; 
	    % plot(xi, f2, 'k-'); hold on; 
	    scatter(xi(Imax), Pmax, 30, 'o', 'filled', 'Displayname', 'MLE'); 
        plot([xi(Imax) xi(Imax)], [0 Pmax], 'k--', 'Displayname', '', ...
            'HandleVisibility','off'); 
	    xline(means(i), 'Displayname', '', 'HandleVisibility','off'); 
        grid on; xlim([xi(1) xi(end)]);
        xlabel(sprintf('%s range', varNames{i}))
        ylabel(sprintf('PDF (%s)', varNames{i}))

        % Adjust legend position
	    if i == nVals
		    leg = legend(); 
		    
		    ax = subplot(3,4,12,'Visible','off');
		    axPos = ax.Position;
		    axPos(2) = axPos(2) + axPos(4)/4;
		    axPos(1) = axPos(1) + axPos(3)/4;
		    leg.Position(1:2) = axPos(1:2);
	    end
    end

end

