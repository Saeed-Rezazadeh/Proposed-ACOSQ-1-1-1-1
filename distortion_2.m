function [D_2] = distortion_2(f , Pr , T_2 , codebook , numLevel , delta)
%% Overall adaptive distortion at step 2
summation = 0 ;
parfor x = 1 : numLevel
    for y_1 = 1 : 2
        for y_prime = 1 : 8
            y = (y_1 - 1) * 8 + y_prime ;
            u_index = find(T_2(: , 1 + y_1) == x) ;
            
            u = T_2(u_index , 1) ;
            
            summation = summation + delta * Pr(x , y) * sum(f(u_index) .* (u - codebook(y)) .^2 ) ;
        end
    end
end
D_2 = summation ;
end