function [vecParallel,vecOrtho] = getProjectPointOnVector(vecRefVector,vecPoint)
	
	
	vecParallel = ((vecRefVector'*vecPoint)/(norm(vecRefVector).^2)).*vecRefVector;
	vecOrtho = vecPoint - vecParallel;
	
end