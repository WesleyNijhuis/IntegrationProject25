%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Crane: RT blocks paramaters, I/O conversions
% Adapted on 20-01-2011 by Mernout Burger
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Values for the angle sensors, and the joystick input
%                       [angle  -angle  joy_y   joy_x]
adinchannels =          [ 1      2       3       4    ]; % input channel
adinoffs =              [ 0.011  0.075  -0.013   0.005]; % input offset for angle and joystick 
adingain =              [-0.851 -0.851   1       1    ]; % input gain (machine units -> degree) 

% Values for the hoist and trolley sensors
%                       [trolley hoist   ]
adinchannelsencoders =  [ 1       2      ];              % encoder output[hoist trolley ]
adinoffsencoders =      [ 0       0      ];              % input offset for encoder inputs
adingainencoders =      [-4.5e-6  2.35e-5];              % input gain (machine units -> meters) [trolley hoist]


% no change
dinchannels = [17:19];                                   % input channel [left right top]

daoutchannels = [1 2];                                   % output channel [trolley hoist]
daoutoffs = [-1e-4 -2e-3];                               % output offset  [trolley hoist]
daoutgain = [0.8 2.5];                                   % output gain    [trolley hoist]

ticklost = 100;                                          % max number of ticks lost

