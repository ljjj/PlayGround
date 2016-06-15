function [N,NS1,NS2] = Ponder_This_2016_Jan(S,m1,m2)

if nargin < 3 
    m2 = 3*S^2;
end
if nargin < 2
    m1 = S;
end

N = 0;
NS1 = [];
NS2 = [];

% there are a total of S^2 products and each covers at most 3 gears, so N
% will not exceed 3*S^2. If the maximum of S1 and S2 exceeds that it is useless
for m = m1:m2
    disp(['Testing m = ',num2str(m)])
    tic
    testMaxWithSetSize(m,S,S);
    toc
end

%     function testMax(M) % find maximal N with S1 S2 constrained to be <= M
%         if M < S % the numbers in S1 should all be distinct, same for S2
%             error('M < S in testMax')
%         end
%         S1 = flip((M-S+1):M); % linspace(M,M-S+1,S); % front chainrings
%         % loop over all possible sets of S1 keeping S1(1) = M
%         while S1(1) == M
%             % loop over all possible sets of S2 starting with S1
%             S2 = S1;
%             while 1 == 1
%                 Ncurr = AllGears(S1,S2);
%                 if Ncurr > N % replace by the new maximal N and the sets
%                     N = Ncurr;
%                     NS1 = S1;
%                     NS2 = S2;
%                 elseif Ncurr == N % register equivalent solutions
%                     NS1 = [NS1; S1];
%                     NS2 = [NS2; S2];
%                 end
%                 if S2(1) == S
%                     break
%                 end
%                 S2 = advanceSet(S2);
%             end
%             if S1(1) == S
%                 break
%             end
%             S1 = advanceSet(S1);
%         end
%     end

    function testMaxWithSetSize(M, size1, size2) % find maximal N with size1 size2 sets constrained to be <= M
        if M < size1 || M < size2 % the numbers in the sets should all be distinct, same for S2
            error('M < size1 or M < size2 in testMaxSetWithSetSize')
        end
        si1 = flip((M-size1+1):M);
        % loop over all possible sets of S1 keeping S1(1) = M
        while si1(1) == M
            % loop over all possible sets of S2 starting with S1
            if size2 == size1 % exploit symmetry
                si2 = si1;
            else
                si2 = flip((M-size2+1):M);
            end
            while 1 == 1
                Ncurr = AllGears(si1,si2);
                if Ncurr > N % replace by the new maximal N and the sets
                    N = Ncurr;
                    NS1 = si1;
                    NS2 = si2;
                    disp(N)
                    disp(NS1)
                    disp(NS2)
                elseif Ncurr == N % register equivalent solutions
                    NS1 = [NS1; si1];
                    NS2 = [NS2; si2];
                end
                if si2(1) == size2
                    break
                end
                si2 = advanceSet(si2);
            end
            if si1(1) == size1
                break
            end
            si1 = advanceSet(si1);
        end
    end

    function Si = advanceSet(Si) % find the next set
        Ns = numel(Si);
        if Si(1) == Ns % assuming Si is sorted in reverse order
            error('cannot advance')
        end
        k = Ns;
        l = 1;
        while Si(k) == l % find the number to change
            k = k - 1; % k > 0 or else the above error catches it
            l = l + 1;
        end
        Si(k) = Si(k) - 1; % the k'th number is changed
        Si((k+1):Ns) = flip((Si(k)-Ns+k):Si(k)-1); % linspace(Si(k)-1,Si(k)-Ns+k,Ns-k); % all numbers after it are consecutive decreasing integers
    end

    function maxN = AllGears(s1,s2) % find the maximal N such that every gear value between 1 and maxN is produced by S1 and S2
        gears = zeros(1,max(s1(:))*max(s2(:))+2);
        for i = 1:numel(s1)
            for j = 1:numel(s2)
                p = s1(i) * s2(j);
                gears(p) = 1;
                gears(p+1) = 1;
                if p > 1
                    gears(p-1) = 1;
                end
            end
        end
        br = find(gears == 0);
        maxN = br(1) - 1;
%         gears = sort(unique(s1(:)*s2(:)')); % available gears
%         gap = diff([-1;gears]) > 3;
%         if ~any(gap)
%             maxN = gears(end) + 1;
%         elseif gap(1)
%             maxN = 0;
%         else
%             gear = gears(gap);
%             gear = gears(gears < gear(1));
%             maxN = gear(end) + 1;
%         end
    end
end