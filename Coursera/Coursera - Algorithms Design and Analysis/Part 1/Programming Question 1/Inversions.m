function [inv, array] = Inversions( array )
    array = array(:);
    N = numel(array);
    
    % deal with cases of elementary arrays
    if N <= 1
        inv = 0;
    elseif N == 2
        inv = array(1) > array(2);
        if inv
            temp = array(2);
            array(2) = array(1);
            array(1) = temp;
        end
    else
        % divide and conquer so that both subarrays are counted and sorted
        half = round(N/2);
        array1 = array(1:half);
        array2 = array(half+1:N);
        [inv1, array1] = Inversions(array1);
        [inv2, array2] = Inversions(array2);
        
        % find the inversions between two subarrays
        % and put them together into a whole array
        i = 1; j = 1; inv3 = 0;
        while i <= half && j <= N-half
            if array1(i) > array2(j)
                inv3 = inv3 + half - i + 1;
                array(i+j-1) = array2(j);
                j = j + 1;
            else
                array(i+j-1) = array1(i);
                i = i + 1;
            end
        end
        while i<= half
            array(i+j-1) = array1(i);
            i = i + 1;
        end
        while j<=N-half
            array(i+j-1) = array2(j);
            j = j + 1;
        end
        
        inv = inv1+inv2+inv3;
    end
end