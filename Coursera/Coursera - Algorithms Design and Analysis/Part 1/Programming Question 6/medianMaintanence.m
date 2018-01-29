function medians = medianMaintenance(array)
% Given an array with no repetitive elements, find the medians of the first
% (unsorted) k elements, for k from 1 to N.
% The median is defined as the ceil(k/2)th smallest element.
% The algorithm is called median maintenance, where two heaps are
% dynamically maintained as more elements are added from the array.

    % initialize global variables
    N = numel(array);
    heapLow = zeros(1,ceil(N/2)); % store ceil(k/2) to k smallest numbers, the first element smallest
    NHeapLow = 0; % size of heapLow
    heapHigh = zeros(1,ceil(N/2)); % store 1 to ceil(k/2) smallest numbers, the first element largest
    NHeapHigh = 0; % size of heapHigh
    medians = zeros(1,N); % medians(k) = median of the first k elements
    
    % dynamic heap maintenance, keeping the median in heapHigh
    hHInsert(array(1));
    medians(1) = array(1);
    for i = 2:N
        key = array(i);
        if key < heapHigh(1)
            hHInsert(key);
            if NHeapHigh > NHeapLow + 1 % rebalance the heaps, heapHigh can have 1 element more than heapLow
                reb = hHExtract();
                hLInsert(reb);
            end
        else
            hLInsert(key);
            if NHeapLow > NHeapHigh % rebalance the heaps
                reb = hLExtract();
                hHInsert(reb);
            end
        end
        medians(i) = heapHigh(1); % median is at the root of heapHigh
    end
    
    % insert key in heapLow
    function hLInsert(k)
        NHeapLow = NHeapLow+1;
        heapLow(NHeapLow) = k;
        swim = NHeapLow;
        parent = floor(swim/2);
        % swim the inserted element towards the root if smaller than parent
        while swim > 1 && heapLow(swim) < heapLow(parent)
            temp = heapLow(swim);
            heapLow(swim) = heapLow(parent);
            heapLow(parent) = temp;            
            swim = parent;
            parent = floor(swim/2);
        end
    end

    % extract root in heapLow
    function r = hLExtract()
        r = heapLow(1);
        % put the last element at root and sink towards leaf if larger than the smaller child
        heapLow(1) = heapLow(NHeapLow);
        NHeapLow = NHeapLow - 1;
        sink = 1;
        if sink*2 > NHeapLow % no child
            return
        elseif sink*2+1 > NHeapLow % single child
            child = sink*2;
        elseif heapLow(sink*2) < heapLow(sink*2+1) % find the smaller child to ensure heap structure after swap
            child = sink*2;
        else
            child = sink*2+1;
        end
        while heapLow(sink) > heapLow(child) % keep sinking
            temp = heapLow(sink);
            heapLow(sink) = heapLow(child);
            heapLow(child) = temp;
            sink = child;
            if sink*2 > NHeapLow
                return
            elseif sink*2+1 > NHeapLow
                child = sink*2;
            elseif heapLow(sink*2) < heapLow(sink*2+1)
                child = sink*2;
            else
                child = sink*2+1;
            end
        end
    end

    % insert key in heapHigh
    function hHInsert(k)
        NHeapHigh = NHeapHigh+1;
        heapHigh(NHeapHigh) = k;
        swim = NHeapHigh;
        parent = floor(swim/2);
        % swim the inserted element towards the root if larger than parent
        while swim > 1 && heapHigh(swim) > heapHigh(parent)
            temp = heapHigh(swim);
            heapHigh(swim) = heapHigh(parent);
            heapHigh(parent) = temp;            
            swim = parent;
            parent = floor(swim/2);
        end
    end

    % extract root in heapHigh
    function r = hHExtract()
        r = heapHigh(1);
        % put the last element at root and sink towards leaf if smaller than the larger child
        heapHigh(1) = heapHigh(NHeapHigh);
        NHeapHigh = NHeapHigh - 1;
        sink = 1;
        if sink*2 > NHeapHigh % no child
            return
        elseif sink*2+1 > NHeapHigh % single child
            child = sink*2;
        elseif heapHigh(sink*2) > heapHigh(sink*2+1) % find the larger child to ensure heap structure after swap
            child = sink*2;
        else
            child = sink*2+1;
        end
        while heapHigh(sink) < heapHigh(child) % keep sinking
            temp = heapHigh(sink);
            heapHigh(sink) = heapHigh(child);
            heapHigh(child) = temp;
            sink = child;
            if sink*2 > NHeapHigh
                return
            elseif sink*2+1 > NHeapHigh
                child = sink*2;
            elseif heapHigh(sink*2) > heapHigh(sink*2+1)
                child = sink*2;
            else
                child = sink*2+1;
            end
        end
    end
end