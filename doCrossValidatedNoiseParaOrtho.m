function [vecNoiseParallelCV,vecNoiseOrthogonalCV,vecNoiseTotalCV] = doCrossValidatedNoiseParaOrtho(matData,vecTrialTypes,intTypeCV)
	%getCrossValidatedNoiseParaOrtho Cross-validated decomposition of
	%			neuronal activity noise in coding and non-coding directions
	%
	%[vecNoiseParallelCV,vecNoiseOrthogonalCV,vecNoiseTotalCV] = ...
	%	getCrossValidatedNoiseParaOrtho(matData,vecTrialTypes,intTypeCV)
	%
	%Inputs:
	% - matData; [n x p]  Matrix of n observations/trials of p predictors/neurons
	% - vecTrialTypes; [n x 1] Trial indexing vector of c classes in n observations/trials
	% - intTypeCV; [int or vec] Integer switch 0-2 or trial repetition vector. 
	%				Val=0, no CV; val=1, leave-one-out CV, val=2 (or
	%				vector), leave-repetition-out. 
	%
	%Outputs:
	% - vecNoiseParallelCV; 
	% - vecNoiseOrthogonalCV;
	% - vecNoiseTotalCV;
	%
	%Version History:
	%2020-12-18 Created function [by Jorrit Montijn]
	
	%% check which kind of cross-validation
	if nargin < 3 || isempty(intTypeCV)
		intTypeCV = 2;
	end
	
	%% prepare
	intVerbose = 0;
	
	%get number of trials
	vecTrialTypes = vecTrialTypes(:);
	intTrials = numel(vecTrialTypes);
	
	%check if matData is [trial x neuron] or [neuron x trial]
	if size(matData,1) == intTrials && size(matData,2) == intTrials
		%number of neurons and trials is the same
		warning([mfilename ':SameNeuronsTrials'],'Number of neurons and trials is identical; please double check the proper orientation of [intNeurons x intTrials]');
	elseif size(matData,1) == intTrials
		%rotate
		matData = matData';
	elseif size(matData,2) == intTrials
		%size is correct
	else
		error([mfilename ':SameNeuronsTrials'],'Size of matData and vecTrialTypes do not match');
	end
	intNeurons = size(matData,1);
	vecUniqueTrialTypes = unique(vecTrialTypes);
	intStimTypes = length(vecUniqueTrialTypes);
	vecTrialTypeIdx = label2idx(vecTrialTypes);
	
	%pre-allocate output
	vecNoiseParallelCV = zeros(1,intTrials);
	vecNoiseOrthogonalCV = zeros(1,intTrials);
	vecNoiseTotalCV = zeros(1,intTrials);
	ptrTic = tic;
	
	%% cross-validate
	if numel(intTypeCV) == intTrials
		%third input is trial repetition index
		vecTrialRepetition = intTypeCV;
		
		%remove repetition
		intRepNum = max(vecTrialRepetition);
		for intRep=1:intRepNum
			%msg
			if toc(ptrTic) > 5
				ptrTic = tic;
				pause(eps);
				if intVerbose > 0,fprintf('Decomposing; now at rep %d/%d [%s]\n',intRep,intRepNum,getTime);end
			end
			
			%remove trials
			indThisRep = vecTrialRepetition==intRep;
			matTrainData = matData(:,~indThisRep);
			vecTrainTrialType = vecTrialTypeIdx(~indThisRep);
			matTestData = matData(:,indThisRep);
			vecTestTrialType = vecTrialTypeIdx(indThisRep);
			
			%get noise decomposition
			[vecNoiseParallelTest,vecNoiseOrthogonalTest,vecNoiseTotalTest] = ...
				getNoiseParaOrthoCV(matTrainData,vecTrainTrialType,matTestData,vecTestTrialType);
	
			%assign
			vecNoiseParallelCV(indThisRep) = vecNoiseParallelTest;
			vecNoiseOrthogonalCV(indThisRep) = vecNoiseOrthogonalTest;
			vecNoiseTotalCV(indThisRep) = vecNoiseTotalTest;
		end
		
		
	elseif intTypeCV == 0
		%no CV
		[vecNoiseParallelCV,vecNoiseOrthogonalCV,vecNoiseTotalCV] = getNoiseParaOrtho(matData,vecTrialTypes);

	elseif intTypeCV == 1
		
		%leave one out
		for intLeaveOut=1:intTrials
			%get info on to-be-left-out trial
			indSelect = ~isnan(vecTrialTypeIdx);
			indSelect(intLeaveOut) = false;
			
			%remove trials
			matTrainData = matData(:,indSelect);
			vecTrainTrialType = vecTrialTypeIdx(indSelect);
			matTestData = matData(:,~indSelect);
			vecTestTrialType = vecTrialTypeIdx(~indSelect);
			
			%get noise decomposition
			[vecNoiseParallelTest,vecNoiseOrthogonalTest,vecNoiseTotalTest] = ...
				getNoiseParaOrthoCV(matTrainData,vecTrainTrialType,matTestData,vecTestTrialType);
	
			%assign
			vecNoiseParallelCV(~indSelect) = vecNoiseParallelTest;
			vecNoiseOrthogonalCV(~indSelect) = vecNoiseOrthogonalTest;
			vecNoiseTotalCV(~indSelect) = vecNoiseTotalTest;
			
			
			%msg
			if toc(ptrTic) > 5
				ptrTic = tic;
				pause(eps);
				if intVerbose > 0,fprintf('Decomposing; now at trial %d/%d [%s]\n',intLeaveOut,intTrials,getTime);end
			end
		end
	elseif intTypeCV == 2
		%remove repetition
		intRepNum = intTrials/intStimTypes;
		if round(intRepNum) ~= intRepNum,error([mfilename ':IncompleteRepetitions'],'Number of repetitions is not an integer');end
		for intRep=1:intRepNum
			%msg
			if toc(ptrTic) > 5
				ptrTic = tic;
				pause(eps);
				if intVerbose > 0,fprintf('Decomposing; now at rep %d/%d [%s]\n',intRep,intRepNum,getTime);end
			end
			
			%remove trials
			intTrialStopRep = intRep*intStimTypes;
			intTrialStartRep = intTrialStopRep - intStimTypes + 1;
			indSelect = true(1,intTrials);
			indSelect(intTrialStartRep:intTrialStopRep) = false;
			
			%remove trials
			matTrainData = matData(:,indSelect);
			vecTrainTrialType = vecTrialTypeIdx(indSelect);
			matTestData = matData(:,~indSelect);
			vecTestTrialType = vecTrialTypeIdx(~indSelect);
			
			%get noise decomposition
			[vecNoiseParallelTest,vecNoiseOrthogonalTest,vecNoiseTotalTest] = ...
				getNoiseParaOrthoCV(matTrainData,vecTrainTrialType,matTestData,vecTestTrialType);
	
			%assign
			vecNoiseParallelCV(~indSelect) = vecNoiseParallelTest;
			vecNoiseOrthogonalCV(~indSelect) = vecNoiseOrthogonalTest;
			vecNoiseTotalCV(~indSelect) = vecNoiseTotalTest;
			
		end
	else
		error([mfilename ':SyntaxError'],'CV type not recognized');
	end
end