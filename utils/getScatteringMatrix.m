function [S] = getScatteringMatrix(B,Q,Z)


%Get scattering matrix S
%from matrices B and Q

n = size(Q, 2);
q = size(Q, 1);
p = size(B, 1);


    if (q <= p)
       Z_inv = inv(Z);
       S = 2*Q'*inv(Q*Z_inv*Q')*Q*Z_inv - eye(n);

       %According to MATLAb it's faster and more accurate like this
       %But also much less readable
       %S = 2*Q'*(Q*(Z\Q')\Q)/Z - eye(n);
     else
       S = eye(n) - 2*Z*B'*inv(B*Z*B')*B;

       %Same as above
       %S = eye(n) - 2*Z*B'*((B*Z*B')\B);

    end


end

