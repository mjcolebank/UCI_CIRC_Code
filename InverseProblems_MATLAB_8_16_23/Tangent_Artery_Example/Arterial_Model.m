%% Arterial pressure-diameter model - algebraic equation
% the general equation is 
%       P(D) = stiffness*tan(pi/gamma * (D/Dref - 1))
% where stiffness (mmHg) describes the material properties of the artery,
% gamma (dimensionless) describes the inflection point in the tangent
% curver, Dref (cm) is the reference diameter of the artery and D (cm) is
% the current diameter

function P = Arterial_Model(D,param,Dref)

stiff     = param(1);
gamma     = param(2);
P0        = param(3);

P = stiff.*tan(pi.*(D./Dref- 1.0)./gamma) + P0;
end