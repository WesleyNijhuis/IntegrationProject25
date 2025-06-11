function S_function_mpc(block)
    setup(block);
end

function setup(block)
    % Define the number of input ports and output ports
    block.NumInputPorts  = 2;  % One for states, one for reference
    block.NumOutputPorts = 1;  % Output for control input

    % Set up port widths (dimensions of input and output)
    block.InputPort(1).Dimensions = 4;  % Assuming state is a 4-dimensional vector
    block.InputPort(2).Dimensions = 4;  % Assuming reference is a 4-dimensional vector
    block.OutputPort(1).Dimensions = 1;  % Output control input

    % Configure the sample time
    block.SampleTimes = [-1, 0];  % Set to appropriate sample time

    % Register the methods for this block
    block.RegBlockMethod('Outputs', @Output);

    % Register the dialog parameter for `fast_opt`
    block.NumDialogPrms = 1;
    %block.DialogPrm = {'fast_opt'};  % Pass the fast_opt as a dialog parameter
    
    % Set the tunability for the dialog parameter (make it non-tunable)
    block.DialogPrmsTunable = {'Nontunable'};  % Specify that it's non-tunable
end

function Output(block)
    %tic
    % Get the input data (states and reference)
    states = block.InputPort(1).Data;     % Current states
    reference = block.InputPort(2).Data;  % Reference signal

    % Get the precompiled optimizer from the block parameter
    fast_opt = block.DialogPrm(1).Data;  % Access the optimizer passed from Simulink

    % Call the MPC function with the states, reference, and the precompiled optimizer (fast_opt)
    [mu,flag] = fast_opt{[states; reference]};
    if flag ~= 0
        u0 = -[-1.8453,4.7438,-0.3246,0.4136]*states; % LQR feedback as fail-safe
        disp('mpc returned NaN')
    else
        u0=mu(1);
    end

    % Set the control output
    block.OutputPort(1).Data = u0;
    %toc
end
