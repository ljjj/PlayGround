function [comp, array] = QuickSort( array, pivotChoice)
    array = array(:);
    N = numel(array);
    % deal with elementary cases
    if N <= 1
        comp = 0;
    elseif N == 2
        comp = 1;
        if array(1) > array(2)
            swap(1,2);
        end
    else
        % choose the pivot and swap to the 1st index
        if pivotChoice == 1
        elseif pivotChoice == 2
            swap(1, N);
        elseif pivotChoice == 3
            mot = [array(1), array(N), array(ceil(N/2))];
            medInd = find(mot == median(mot));
            if medInd == 3
                swap(1, ceil(N/2));
            elseif medInd == 2
                swap(1,N);
            end
        else
            error('Wrong parameter input for pivotChoice')
        end
        
        % implement partition
        i = 2;
        for j = 2:N
            if array(j) < array(1)
                swap(i, j);
                i = i + 1;
            end
        end
        swap(1, i-1);
        
        % divide and conquer so subarrays are sorted
        array1 = array(1:(i-2));
        array2 = array(i:N);
        [comp1, array1] = QuickSort(array1, pivotChoice);
        [comp2, array2] = QuickSort(array2, pivotChoice);
        comp = comp1 + comp2 + N-1;
        array = [array1; array(i-1); array2];
    end
    
    % swap the array with specified indices
    function swap(p, q)
        if p == q
            return;
        end
        temp = array(q);
        array(q) = array(p);
        array(p) = temp;
    end
end