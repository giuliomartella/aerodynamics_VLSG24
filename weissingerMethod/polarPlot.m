function polarPlot(alpha_range, wing, tail)

cl_values_wing = zeros(size(alpha_range));
cl2d_wing = zeros(length(alpha_range), wing.discretize(2));
span_wing = linspace(-0.5 * wing.span, 0.5 * wing.span, wing.discretize(2));

if nargin == 3
    % Initialize variables for tail if provided
    cl_values_tail = zeros(size(alpha_range));
    cl2d_tail = zeros(length(alpha_range), tail.discretize(2));
    span_tail = linspace(-0.5 * tail.span, 0.5 * tail.span, tail.discretize(2));
end

for i = 1:length(alpha_range)
    alpha = alpha_range(i);
    uInf = [cos(alpha), -sin(alpha), 0; sin(alpha), cos(alpha), 0; 0, 0, 1] * [1.0; 0.0; 0.0];

    if nargin == 3
        % Solve the system with both wing and tail
        [wing, tail] = buildLinearSystem(uInf, wing, tail);
    else
        % Solve the system only for the wing
        wing = buildLinearSystem(uInf, wing);
    end

    % Post-process wing data
    wing = postprocessing(wing);
    cl_values_wing(i) = wing.cL;
    cl2d_wing(i, :) = wing.cL2D;

    if nargin == 3
        % Post-process tail data if provided
        tail = postprocessing(tail);
        cl_values_tail(i) = tail.cL;
        cl2d_tail(i, :) = tail.cL2D;
    end
end

% Plot results
figure;

% Wing polar plot
subplot(2, 2, 1);
plot(alpha_range * 180/pi, cl_values_wing, 'o-', 'LineWidth', 2);
xlabel('Angle of Attack (degrees)');
ylabel('C_L');
title('Wing Polar Plot');
grid on;

if nargin == 3
    % Tail polar plot
    subplot(2, 2, 2);
    plot(alpha_range * 180/pi, cl_values_tail, 'o-', 'LineWidth', 2);
    xlabel('Angle of Attack (degrees)');
    ylabel('C_L');
    title('Tail Polar Plot');
    grid on;
end

% Wing span-wise lift distribution
subplot(2, 2, 3);
for i = 1:length(alpha_range)
    plot(span_wing, cl2d_wing(i, :), 'LineWidth', 2);
    hold on;
end
xlabel('Span (m)');
ylabel('C_L2D');
title('Wing Lift Distribution');
grid on;

if nargin == 3
    % Tail span-wise lift distribution
    subplot(2, 2, 4);
    for i = 1:length(alpha_range)
        plot(span_tail, cl2d_tail(i, :), 'LineWidth', 2);
        hold on;
    end
    xlabel('Span (m)');
    ylabel('C_L2D');
    title('Tail Lift Distribution');
    grid on;
end

end
