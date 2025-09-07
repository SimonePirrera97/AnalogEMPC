q = 1;
Q = q*eye(size(Ad));
P = dlyap(Ad', Q);  % P solves A' * P * A - P = -Q
% Check positive definiteness
eigenvalues = eig(P);
is_pos_def = all(eigenvalues > 0);
disp(['P is positive definite: ', string(is_pos_def)])
Ad'*P*Ad-P  % should give -Q

normP(Ad,P) + normP(T_sw*0.5*[Ad^2-eye(2)]*inv(A)*V_in,P)*normP(B2*F_c{2}(1:2),P)

function nPvec = normPvec(v,P)
    nPvec = sqrt(v'*P*v);
end

function nP = normP(A,P)
    nP = sqrt(max(eig(P\(A'*P*A))));
end