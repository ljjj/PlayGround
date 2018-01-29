function nTarget = TwoSum(array, left, right)
% compute the number of target values t in the interval [left,right]
% (inclusive) such that there are distinct numbers x,y in the input file that satisfy x+y=t
% assuming all entries integers

    % initialization
    if left > right % make sure left <= right
        temp = left;
        left = right;
        right = temp;
    end
    array = unique(sort(array)); % processed array
    targetFound = false(1,right-left+1); % whether a target value is found
    N = numel(array);
    
    % for each x find y between [left, right]
    for i = 1:N
        if mod(i,10000) == 0
            disp(i);
        end
        x = array(i);
        L = left - x; % y search left bound
        R = right - x; % y search right bound
        indL = ceilInd(L, i+1, N);
        if indL > N || array(indL) > R % no y can be found
            continue;
        end
        indR = ceilInd(R, indL+1, N) - 1; % if ceilInd returns key position exceeding R
        if indR < N && array(indR+1) == R % if ceilInd returns key position exactly R
            indR = indR + 1;
        end
        % sanity check
        if indL < 1 || (indL > i+1 && x+array(indL-1) > left)
            error(['indL wrong at i = ',num2str(i)])
        end
        if indR > N || (indR < N && x + array(indR+1) < right)
            error(['indR wrong at i = ',num2str(i)])
        end
        % loop over y values to find targets
        for j = indL:indR
            y = array(j);
            t = x + y;
            targetFound(t - left + 1) = true;
        end
    end
    
    nTarget = sum(targetFound);
    
    % find the index of the smallest value that is not below key
    % in the unique array between index ll and rr
    % if ind = N+1 no entry exceeds key
    function ind = ceilInd(key, ll, rr)
        if ll > N || array(ll) > key
            ind = ll;
            return
        end
        if array(rr) < key
            ind = rr+1;
            return
        end
        
        m = floor((ll+rr)/2);
        if array(m) < key
            ind = ceilInd(key, m+1, rr);
        elseif array(m) > key
            ind = ceilInd(key, ll, m-1);
        else
            ind = m;
        end
    end
end