    function result = double_factorial(n)
        % Check if n is a non-negative integer
        if n < 0 || floor(n) ~= n
            error('Input must be a non-negative integer.');
        end

        % Calculate double factorial
        if n == 0 || n == 1
            result = 1; % Base case: 0!! = 1 and 1!! = 1
        else
            result = 1; % Initialize result
            for k = n:-2:1
                result = result * k; % Multiply by every second integer
            end
        end
    end
