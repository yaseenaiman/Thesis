function FRF = DiscreteRFGeneration(Fopt,Steering_Vec,NtRf)
         [~, Ns] = size(Fopt);
    FRF = [];
    Fres = Fopt;

    % Assuming you want the top NtRf directions
numDirectionsToSelect = NtRf;

% Initialize arrays to store the top directions and their indices
topDirections = zeros(size(Steering_Vec, 1), numDirectionsToSelect);
topIndices = zeros(1, numDirectionsToSelect);

% Loop to find the top directions
for i = 1:numDirectionsToSelect
    % Calculate Epsi
    Epsi_ZF = Steering_Vec' * Fres;
    
    % Find the maximum index
    [~, Ind_Direction2] = max(diag(Epsi_ZF * Epsi_ZF'));
    
    % Storing the result
    topDirections(:, i) = Steering_Vec(:, Ind_Direction2);
    topIndices(i) = Ind_Direction2;
        FRF =[FRF Steering_Vec(:, Ind_Direction2) ];


 % Update Fres based on the selected direction (GramShmidt)
  Fres=Fres-(FRF * (FRF' * Fres)/norm(FRF * (FRF'),'fro'));

end

end