% Mobility-MLEACH Implementation in MATLAB with Implemntation of Energy and Mobility Model 
% Mobility Induced Multi-Hop LEACH Protocol in Heterogeneous Mobile Network
% Implemented by LARAIB AZMAT 
% Github ->  Laraib-Azmat

clear all;

% Parameters
numNodes = 20; % Total number of nodes in the network
numSuperNodes = 1; % Number of super nodes
numAdvancedNodes = 5; % Number of advanced nodes
numNormalNodes = numNodes - numSuperNodes - numAdvancedNodes; % Remaining nodes are normal nodes
numClusters = min(3, numAdvancedNodes); % Number of clusters, minimum of 3 or the number of advanced nodes
areaWidth = 100; % Width of the simulation area
areaHeight = 100; % Height of the simulation area
initialEnergy = 100; % Initial energy level of each node
maxIterations = 100; % Maximum number of iterations for the simulation
RSSI = 10; % RSSI (Received Signal Strength Indication) threshold (signal strength between nodes)

% Initialize Nodes
nodes = struct('id', 1:numNodes, ... % Assign IDs to nodes
    'x', randi(areaWidth, 1, numNodes), ... % Randomly assign x coordinates within the area
    'y', randi(areaHeight, 1, numNodes), ... % Randomly assign y coordinates within the area
    'energy', initialEnergy, ... % Set initial energy level
    'clusterHead', false, ... % Initialize as non-cluster heads
    'clusterID', 0, ... % Initialize cluster ID
    'nodeType', 'normal'); % Set node type as normal initially

% Assign node types
nodeTypes = {'super', 'advanced', 'normal'}; % Define node types
nodeTypeCounts = [numSuperNodes, numAdvancedNodes, numNormalNodes]; % Count of each node type
startIndex = 1; % Start index for node assignment
for i = 1:length(nodeTypes)
    endIndex = startIndex + nodeTypeCounts(i) - 1; % Calculate end index for node assignment
    [nodes(startIndex:endIndex).nodeType] = deal(nodeTypes{i}); % Assign node types to nodes
    startIndex = endIndex + 1; % Update start index for next node type
end

% Cluster Initialization - Select cluster heads only from advanced nodes
advancedNodeIndices = find(strcmp({nodes.nodeType}, 'advanced')); % Find indices of advanced nodes
if ~isempty(advancedNodeIndices) % Check if there are advanced nodes available
    clusterHeadIndices = advancedNodeIndices(randperm(length(advancedNodeIndices), numClusters)); % Randomly select cluster heads from advanced nodes
    for i = 1:length(clusterHeadIndices)
        nodes(clusterHeadIndices(i)).clusterHead = true; % Mark selected nodes as cluster heads
        nodes(clusterHeadIndices(i)).clusterID = i; % Assign cluster ID to selected cluster heads
    end
else
    disp('Error: No advanced nodes available to select as cluster heads.'); % Display error message if no advanced nodes are available
end

% Simulation Loop
for iter = 1:maxIterations % Iterate through each time step
    disp(['Iteration: ' num2str(iter)]); % Display current iteration number
    
    % Mobility Model (Random Waypoint)
    nodes = updateNodePositions(nodes, areaWidth, areaHeight); % Update node positions based on mobility model
    
    % Energy Model (Simple) - Update energy calculations
    nodes = updateEnergy(nodes, RSSI); % Update energy levels of nodes
    
    % Cluster Formation
    nodes = leachClustering(nodes, numClusters); % Perform clustering based on LEACH protocol
    
    % Data Transmission (Simple) - Display data transmission information
    disp('Data Transmission:');
    for i = 1:length(nodes)
        if ~nodes(i).clusterHead % If node is not a cluster head
            clusterHeadID = find([nodes.clusterID] == nodes(i).clusterID & [nodes.clusterHead], 1); % Find the cluster head for the node
            if ~isempty(clusterHeadID) % If cluster head is found
                distance = sqrt((nodes(i).x - nodes(clusterHeadID).x)^2 + (nodes(i).y - nodes(clusterHeadID).y)^2); % Calculate distance between node and cluster head
                energyConsumption = distance * 0.1; % Calculate energy consumption for data transmission
                disp(['Node ' num2str(nodes(i).id) ' transmitted data to Cluster Head ' num2str(clusterHeadID) ...
                      '. Distance: ' num2str(distance) ', Energy Consumed: ' num2str(energyConsumption)]); % Display data transmission information
            end
        end
    end
    nodes = dataTransmission(nodes); % Perform data transmission
    
    % M-LEACH Routing Protocol
    [nodes, clusterHeads, networkEnergy] = MLEACHRouting(nodes, numNodes, numSuperNodes, numAdvancedNodes, initialEnergy); % Perform M-LEACH routing
    
    % Energy Calculation - Display energy calculations after updating
    disp('Energy Calculation:');
    for i = 1:length(nodes)
        disp(['Node ' num2str(nodes(i).id) ' Energy updated! ' num2str(nodes(i).energy)]); % Display energy level of each node
    end
    
    % Visualization
    visualizeSimulation(nodes, areaWidth, areaHeight); % Visualize the current state of the network
    pause(0.5); % Pause for visualization
end

% Function Definitions

function nodes = updateNodePositions(nodes, areaWidth, areaHeight)
    % Implement mobility model here
    for i = 1:length(nodes)
        % Example: Random waypoint model
        if rand < 0.05 % Probability of node movement
            nodes(i).x = randi(areaWidth); % Randomly update x coordinate
            nodes(i).y = randi(areaHeight); % Randomly update y coordinate
        end
    end
end

function nodes = updateEnergy(nodes, RSSI)
    % Implement energy model here
    for i = 1:length(nodes)
        % Example: Simple energy depletion
        nodes(i).energy = nodes(i).energy - randi(5); % Randomly deplete energy level
        
        % Calculate distance between advanced nodes and normal nodes
        advancedNodes = nodes(strcmp({nodes.nodeType}, 'advanced')); % Extract advanced nodes
        normalNodes = nodes(strcmp({nodes.nodeType}, 'normal')); % Extract normal nodes
        for j = 1:length(normalNodes)
            for k = 1:length(advancedNodes)
                distance = sqrt((normalNodes(j).x - advancedNodes(k).x)^2 + (normalNodes(j).y - advancedNodes(k).y)^2); % Calculate distance between node i and cluster head j
                if distance < RSSI % If distance is less than RSSI threshold
                    % Calculate Emob
                    q = 0.5; % Constant for Emob calculation
                    d = distance; % Distance between nodes
                    m = 0.5; % Constant for Emob calculation
                    v = 1; % Velocity of the node
                    Emob = q * d + 0.5 * m * v; % Calculate Emob
                    % Calculate Ec using Ech, Ench
                    Ech = 0; % Placeholder for Ech calculation
                    Ench = 0; % Placeholder for Ench calculation
                    Ec = Ech + Ench; % Calculate Ec
                    % Calculate Egrp
                    Egrp = 0; % Placeholder for Egrp calculation
                    % Calculate Enet
                    Enet = 0; % Placeholder for Enet calculation
                end
            end
        end
    end
end

function nodes = leachClustering(nodes, numClusters)
    % Implement clustering logic here (e.g., LEACH protocol)
    % Calculate distances between nodes and cluster heads
    for i = 1:length(nodes)
        if ~nodes(i).clusterHead % If node is not a cluster head
            minDistance = inf; % Initialize minimum distance
            closestClusterID = 0; % Initialize closest cluster ID
            for j = 1:length(nodes)
                if nodes(j).clusterHead % If node j is a cluster head
                    % Check if cluster head has valid coordinates
                    if ~isempty(nodes(j).x) && ~isempty(nodes(j).y)
                        % Calculate distance between node i and cluster head j
                        distance = sqrt((nodes(i).x - nodes(j).x)^2 + (nodes(i).y - nodes(j).y)^2);
                        if distance < minDistance % If distance is smaller than current minimum distance
                            minDistance = distance; % Update minimum distance
                            closestClusterID = nodes(j).clusterID; % Update closest cluster ID
                        end
                    end
                end
            end
            % Assign cluster ID to the node
            nodes(i).clusterID = closestClusterID; % Assign the closest cluster ID to the node
        end
    end
    
end

function nodes = dataTransmission(nodes)
    % Implement data transmission logic here
    % Example: Nodes communicate with their cluster heads
    % Update energy levels accordingly
    for i = 1:length(nodes)
        if ~nodes(i).clusterHead % If node is not a cluster head
            % Find the cluster head for the node
            clusterHeadID = find([nodes.clusterID] == nodes(i).clusterID & [nodes.clusterHead], 1);
            if ~isempty(clusterHeadID) % If cluster head is found
                % Calculate distance between node and cluster head
                distance = sqrt((nodes(i).x - nodes(clusterHeadID).x)^2 + (nodes(i).y - nodes(clusterHeadID).y)^2);
                % Example: Energy consumption for data transmission
                energyConsumption = distance * 0.1; % Calculate energy consumption
                nodes(i).energy = nodes(i).energy - energyConsumption; % Update energy level of the node
            end
        end
    end
end

function [nodes, clusterHeads, networkEnergy] = MLEACHRouting(nodes, numNodes, numSuperNodes, numAdvancedNodes, initialEnergy)
    % Initialize cluster heads list
    clusterHeads = [];
    
    % Step 1: Select Cluster Heads
    % Example: Randomly select cluster heads based on residual energy
    advancedNodeIndices = find(strcmp({nodes.nodeType}, 'advanced')); % Find indices of advanced nodes
    advancedNodes = nodes(advancedNodeIndices); % Extract advanced nodes
    
    % Sort advanced nodes by residual energy
    [~, sortedIndices] = sort([advancedNodes.energy], 'descend'); % Sort advanced nodes by energy level
    sortedAdvancedNodes = advancedNodes(sortedIndices); % Get sorted list of advanced nodes
    
    % Select cluster heads from advanced nodes
    numClusters = min(3, length(sortedAdvancedNodes)); % Number of clusters is the minimum of 3 or the number of advanced nodes
    numClusterHeads = min(numClusters, length(sortedAdvancedNodes)); % Number of cluster heads is the minimum of numClusters or the number of advanced nodes
    clusterHeadCandidates = sortedAdvancedNodes(1:numClusterHeads); % Select candidate cluster heads
    clusterHeadIndices = [clusterHeadCandidates.id]; % Extract indices of cluster head candidates
    
    % Update nodes with selected cluster heads
    for i = 1:length(clusterHeadCandidates)
        nodeIndex = find([nodes.id] == clusterHeadCandidates(i).id); % Find index of cluster head candidate in nodes array
        nodes(nodeIndex).clusterHead = true; % Mark node as cluster head
        nodes(nodeIndex).clusterID = i; % Assign cluster ID sequentially
        clusterHeads = [clusterHeads, nodes(nodeIndex)]; % Add node to cluster heads list
    end
    
    % Step 2: Energy Dissipation Model
    % Implement energy dissipation due to mobility and data transmission
    for i = 1:length(nodes)
        % Example: Energy consumption due to node movements
        energyConsumptionMovement = randi(5); % Random energy consumption due to node movement
        nodes(i).energy = nodes(i).energy - energyConsumptionMovement; % Update energy level
        
        %  Energy consumption due to data transmission
        if ~nodes(i).clusterHead % If node is not a cluster head
            % Find the cluster head for the node
            clusterHeadID = find([nodes.clusterID] == nodes(i).clusterID & [nodes.clusterHead], 1);
            if ~isempty(clusterHeadID) % If cluster head is found
                % Calculate distance between node and cluster head
                distance = sqrt((nodes(i).x - nodes(clusterHeadID).x)^2 + (nodes(i).y - nodes(clusterHeadID).y)^2);
                energyConsumptionTransmission = distance * 0.1; % Calculate energy consumption for data transmission
                nodes(i).energy = nodes(i).energy - energyConsumptionTransmission; % Update energy level
            end
        end
    end
    
    % Step 3: Network Energy Calculation
    % Update network energy considering mobility-induced energy consumption
    networkEnergy = sum([nodes.energy]); % Calculate total network energy
    
    % Return updated nodes, cluster heads, and network energy
end

function visualizeSimulation(nodes, areaWidth, areaHeight)
    clf; % Clear figure
    hold on; % Hold current plot
    axis([0 areaWidth 0 areaHeight]); % Set axis limits
    
    % Define colors and markers for different node types
    nodeTypes = {'super', 'advanced', 'normal'}; % Define node types
    nodeMarkers = {'o', 's', '^'}; % Define markers for node types circle square and triangle
    nodeColors = {'g', 'r', 'b'}; % Define colors for node types green red and blue
    markerSize = 6; % Marker size for normal nodes
    
    % Plot nodes
    for typeIndex = 1:length(nodeTypes)
        typeNodes = nodes(strcmp({nodes.nodeType}, nodeTypes{typeIndex})); % Get nodes of current type
        if ~isempty(typeNodes) % If there are nodes of this type
            plot([typeNodes.x], [typeNodes.y], [nodeColors{typeIndex} nodeMarkers{typeIndex}], 'MarkerSize', markerSize); % Plot nodes with appropriate color and marker
        end
    end
    
    % Plot cluster heads with black color
    clusterHeadIndices = find([nodes.clusterHead]); % Find indices of cluster heads
    clusterHeadNodes = nodes(clusterHeadIndices); % Extract cluster head nodes
    if ~isempty(clusterHeadNodes) % If there are cluster heads
        plot([clusterHeadNodes.x], [clusterHeadNodes.y], 'ko', 'MarkerSize', 12); % Plot cluster heads in black color
    end
    
    % Add legend
    legend('Super Node', 'Advanced Node', 'Normal Node', 'Cluster Head'); % Add legend with node types
    
    title('Mobility M-LEACH Simulation'); % Set title
    xlabel('X'); % Set x-axis label
    ylabel('Y'); % Set y-axis label
    hold off; % Release current plot
    drawnow; % Update figure
end