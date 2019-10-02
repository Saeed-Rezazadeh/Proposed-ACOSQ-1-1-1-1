function [D_3] = distortion_3(f , Pr , T_3 , codebook , numLevel , delta)
%% Overall adaptive distortion at step 3
summation = 0 ;
parfor y_1_y_2 = 1 : 4
    for y_prime = 1 : 4
        y = (y_1_y_2 - 1) * 4 + y_prime ;
        for x = 1 : numLevel
            u_index = find(T_3(: , 1 + y_1_y_2) == x ) ;
            u = T_3(u_index , 1) ;
            

            summation = summation + Pr(x , y) * delta * sum(f(u_index) .* (u - codebook(y)) .^ 2) ;
        end
    end
end
D_3 = summation ;
end