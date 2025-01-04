function [FRF, FBB] = Novel_Algorithm(Fopt, NRF, Steering_Vec1)
    [~, Ns] = size(Fopt);
    FRF = [];
    Fres = Fopt;

    %  top NRF directions
numDirectionsToSelect = NRF;



% Loop to find the top directions
for i = 1:numDirectionsToSelect
    % Calculate Epsi
    Epsi = Steering_Vec1' * Fres;
    
    % Find the maximum index
    [~, Ind_Direction] = max(diag(Epsi * Epsi'));
    
    % Find FRF
        FRF =[FRF Steering_Vec1(:, Ind_Direction)];

 % Update Fres  (Gram-Shmidt)
  Fres=Fres-((Steering_Vec1(:, Ind_Direction) * (Steering_Vec1(:, Ind_Direction)' * Fres))/norm((Steering_Vec1(:, Ind_Direction) * Steering_Vec1(:, Ind_Direction)'),'fro'));

end

%% FBB Calculation 
    [U,S,V] = svd(Fopt'*FRF);
    FBB = V(:,[1:Ns])*U';

    
    end