function [vecNoiseParallelTest,vecNoiseOrthogonalTest,vecNoiseTotalTest] = getNoiseParaOrthoCV(matSpikeCountsTrain,vecOrientationTrain,matSpikeCountsTest,vecOrientationTest)
	%UNTITLED2 Summary of this function goes here
	%   Detailed explanation goes here
		
	[vecUniqueOris,dummy,vecOriIdx] = unique(vecOrientationTrain);
	intTrainTrials = numel(vecOrientationTrain);
	intTestTrials = numel(vecOrientationTest);
	intOriNum = numel(vecUniqueOris);

	vecNoiseParallelTest = nan(1,intTestTrials);
	vecNoiseOrthogonalTest = nan(1,intTestTrials);
	vecNoiseTotalTest = nan(1,intTestTrials);
	for intOriIdx=1:intOriNum
		%stim1
		matTrainAct1 = matSpikeCountsTrain(:,vecOriIdx==intOriIdx);
		vecTrainMu1 = mean(matTrainAct1,2);
		
		%stim2
		intOriIdx2 = modx(intOriIdx+1,intOriNum);
		matTrainAct2 = matSpikeCountsTrain(:,vecOriIdx==intOriIdx2);
		vecTrainMu2 = mean(matTrainAct2,2);
		
		%get ref & project
		vecTrainRef = vecTrainMu2 - vecTrainMu1;
		vecTrainRef = vecTrainRef / norm(vecTrainRef);
		matRecenteredTestPoints = matSpikeCountsTest  - vecTrainMu1;
		
		%pre-alloc
		vecTestTrialsOri1 = find(vecOrientationTest==intOriIdx);
		intTestRepNum = numel(vecTestTrialsOri1);
		vecTempNoiseParallel = nan(1,intTestRepNum);
		vecTempNoiseOrthogonal = nan(1,intTestRepNum);
		vecTempNoiseTotal = nan(1,intTestRepNum);
		for intTestRep = 1:intTestRepNum
			intTrial = vecTestTrialsOri1(intTestRep);
			vecPoint = matRecenteredTestPoints(:,intTrial);
			
			%project vecPoint onto vecRef and decompose in parallel & orthogonal
			[vecParallel,vecOrtho] = getProjectPointOnVector(vecTrainRef,vecPoint);
			vecTempNoiseParallel(intTestRep) = norm(vecParallel);
			vecTempNoiseOrthogonal(intTestRep) = norm(vecOrtho);
			vecTempNoiseTotal(intTestRep) = norm(vecPoint - vecTrainRef);
		end
		
		%assign
		vecNoiseParallelTest(vecTestTrialsOri1) = vecTempNoiseParallel;
		vecNoiseOrthogonalTest(vecTestTrialsOri1) = vecTempNoiseOrthogonal;
		vecNoiseTotalTest(vecTestTrialsOri1) = vecTempNoiseTotal; %sqrt(p^2 + o^2)
	end
end

