function States = Puzzle_2016_Jan(N) % Professor Random picks two numbers from 1 to N independently

% student constants
Daphne = 0.1;
Max = 0.2;
Mindy = 0.3;
Sam = 0.4;
Tim = 0.5;

% initialize
States = zeros(1,N*(N+1)/2); % Day + Student;
% Difs = zeros(1,N*(N+1)/2);
% Maxs = zeros(1,N*(N+1)/2);
% Mins = zeros(1,N*(N+1)/2);
% Sums = zeros(1,N*(N+1)/2);
% Pros = zeros(1,N*(N+1)/2);
% for i = 1:N
%     for j = 1:i
%         Difs(SUB2IND(i,j)) = i-j;
%         Maxs(SUB2IND(i,j)) = i;
%         Mins(SUB2IND(i,j)) = j;
%         Sums(SUB2IND(i,j)) = i+j;
%         Pros(SUB2IND(i,j)) = i*j;
%     end
% end
[i,j] = IND2SUB(1:N*(N+1)/2);
Difs = i-j;
Maxs = i;
Mins = j;
Sums = i+j;
Pros = i.*j;

day = 0;
noChange = false; % game finishes if there's no more change (whether all combinations completed or not)
while ~noChange
    day = day + 1;
    noChange = true;
    % Daphne
    unsolved = find(States == 0); % unsolved cases
    id = Only(Difs(unsolved));
    if ~isempty(id)
        noChange = false;
    end
    States(unsolved(id)) = day + Daphne; % solved by Daphne on this day
    % Max
    unsolved = find(States == 0); % unsolved cases
    id = Only(Maxs(unsolved));
    if ~isempty(id)
        noChange = false;
    end
    States(unsolved(id)) = day + Max; % solved by Max on this day
    % Mindy
    unsolved = find(States == 0); % unsolved cases
    id = Only(Mins(unsolved));
    if ~isempty(id)
        noChange = false;
    end
    States(unsolved(id)) = day + Mindy; % solved by Mindy on this day
    % Sam
    unsolved = find(States == 0); % unsolved cases
    id = Only(Sums(unsolved));
    if ~isempty(id)
        noChange = false;
    end
    States(unsolved(id)) = day + Sam; % solved by Sam on this day
    % Tim
    unsolved = find(States == 0); % unsolved cases
    id = Only(Pros(unsolved));
    if ~isempty(id)
        noChange = false;
    end
    States(unsolved(id)) = day + Tim; % solved by Tim on this day
end

    function [I,J] = IND2SUB(IND)
        I = N+1-round(N+1-sqrt(1+8*IND)/2);
        J = IND - I.*(I-1)/2;
    end

    function IND = SUB2IND(I,J)
        if I < J % swap to make sure I >= J
            J = I + J;
            I = J - I;
            J = J - I;
        end
        IND = I.*(I-1)/2 + J;
    end

    function ID = Only(Vector) % find the index where the value only appears once in the Vector
        ID = [];
        uni = unique(Vector);
        for k = 1:numel(uni)
            ids = Vector == uni(k);
            if sum(ids) == 1
                ID = [ID find(ids == 1)];
            end
        end
    end
end