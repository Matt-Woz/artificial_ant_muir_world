clear;

population_size = 200;


population = generateInitialPopulation(population_size);
newPop = runProgram(population, population_size, 300);

%This function generates the intiial population represented as a
%(population_sizex30) matrix with randomly generated genes
function original_population = generateInitialPopulation(population_size)
    original_population = zeros(population_size, 30);
    for i = 1:population_size
        chromosome = zeros(1, 30);
        for k = 1:30/3
            chromosome((k*3)-2) = randi([1 4]);
            chromosome((k*3)-1) = randi([0 9]);
            chromosome(k*3) = randi([0 9]);
        end
        original_population(i, :) = chromosome(:);
    end
    newColumn = zeros(size(original_population,1),1);
    original_population = [original_population, newColumn];
end

function [population, trail] = calculateFitness(population, worldGrid, population_size)
    for i = 1:population_size
        [fitness, trail] = simulate_ant(worldGrid, population(i, 1:31));
        population(i, 31) = fitness;
    end
    population = sortrows(population, 31);
end

function finalPop = runProgram(currentPopulation, population_size, iterations)
    fitness_data = zeros(iterations, 1);
    for i = 1:iterations
        [currentPopulation(:, 1:30), maxFitness] = createNewPopulation(currentPopulation, population_size);
        fitness_data(i) = maxFitness;
    end
    [finalPop, trail] = calculateFitness(currentPopulation, dlmread('muir_world.txt',' '), population_size);
    maxFitness = max(finalPop(:,end));
    fitness_data(iterations) = maxFitness; %Update final fitness score
    plotFitnessScore(fitness_data, iterations);
    plotBestTrail(maxFitness, iterations, trail);
end

function [newPopulation, maxFitness] = createNewPopulation(currentPopulation, population_size)
    %Find fitness score
    [currentPopulation, ~] = calculateFitness(currentPopulation, dlmread('muir_world.txt',' '), population_size);
    maxFitness = max(currentPopulation(:, end));
    %Elitism - Finds best 2 chromosomes to keep in new population
    newPopulation = zeros(population_size, 30);
    newPopulation(1:2, :) = currentPopulation(population_size-1:population_size, 1:30);
    populationNewNum = 2;
    %Selection process
    while (populationNewNum < population_size)
        weights = currentPopulation(:,31)/sum(currentPopulation(:,31));
        choice1 = rouletteWheelSelection(weights);
        choice2 = rouletteWheelSelection(weights);
        temp_chromosome1 = currentPopulation(choice1,1:30);
        temp_chromosome2 = currentPopulation(choice2,1:30);
        %Crossover process
        if (rand < 0.8)
            [temp_chromosome1, temp_chromosome2] = onePointCrossover(temp_chromosome1, temp_chromosome2);
        end
        %Mutation process
        if (rand < 0.3)
            temp_chromosome1 = swapMutation(temp_chromosome1);
        end
        if (rand < 0.3)
            temp_chromosome2 = swapMutation(temp_chromosome2);
        end
        populationNewNum = populationNewNum + 1;
        newPopulation(populationNewNum,:) = temp_chromosome1;
        if (populationNewNum < population_size)
            populationNewNum = populationNewNum + 1;
            newPopulation(populationNewNum,:) = temp_chromosome2;
        end
    end
end


function swappedChromosome = swapMutation(chromosome)
    point1 = randi([1, length(chromosome)]);
    if mod(point1,3) == 1
        point2 = randi([1, length(chromosome)]);
        while(point2 == point1 || mod(point2,3) ~= 1)
            point2 = randi([1, length(chromosome)]);
        end
    else
        point2 = randi([2, length(chromosome)]);
        while(point2 == point1 || mod(point2,3) == 1)
            point2 = randi([2, length(chromosome)]);
        end
    end
    [chromosome(point1), chromosome(point2)] = swap(chromosome(point1), chromosome(point2));
    swappedChromosome = chromosome;
end


%Mutation option 2 -> Scramble???


%Mutation option 3 -> Inversion??

function [b, a] = swap(a, b)
end


function [temp_chromosome1, temp_chromosome2] = onePointCrossover(temp_chromosome1, temp_chromosome2)
    temp_i = randi([1,29]);
    temp_chromosome1 = [temp_chromosome1(1:temp_i) temp_chromosome2(temp_i+1:end)];
    temp_chromosome2 = [temp_chromosome2(1:temp_i) temp_chromosome1(temp_i+1:end)];
end

%crossover option 2 -> 2 point crossover?


%crossover option 3 -> uniform crossover? 

function choice = rouletteWheelSelection(weights)
  accumulation = cumsum(weights);
  p = rand();
  chosen_index = -1;
  for index = 1 : length(accumulation)
    if (accumulation(index) > p)
      chosen_index = index;
      break;
    end
  end
  choice = chosen_index;
end

%Selection option 2 -> Tournament selection
function choice = TournamentSelection(fitness)
    [fitness, index] = sortrows(fitness);
    chr1 = randi([1, length(fitness)]);
    chr2 = randi([1, length(fitness)]);
    if(fitness(chr1)> fitness(chr2))
        choice = index(chr1);
    else
        choice = index(chr2);
    end
 end


%Selection Option 3 -> ???



function plotFitnessScore(fitness_data, Ngen)
    hf = figure(1); set(hf,'Color',[1 1 1]);
    hp = plot(1:Ngen,100*fitness_data/89,'r');
    set(hp,'LineWidth',2);
    axis([0 Ngen 0 100]); grid on;
    xlabel('Generation number');
    ylabel('Ant fitness [%]');
    title('Ant fitness as a function of generation');
end

function plotBestTrail(best_fitness, Ngen, trail)
    % read the John Moir Trail (world)
    filename_world = 'muir_world.txt'; 
    world_grid = dlmread(filename_world,' ');
    % display the John Moir Trail (world)
    world_grid = rot90(rot90(rot90(world_grid)));
    xmax = size(world_grid,2);
    ymax = size(world_grid,1);
    hf = figure(2); set(hf,'Color',[1 1 1]);
    for y=1:ymax
        for x=1:xmax
            if(world_grid(x,y) == 1)
                h1 = plot(x,y,'sk');
                hold on
            end
        end
    end
    grid on
    % display the fittest Individual trail
    for k=1:size(trail,1)
        h2 = plot(trail(k,2),33-trail(k,1),'*m');
        hold on
    end
    axis([1 32 1 32])
    title_str = sprintf('John Muri Trail - Hero Ant fitness %d%% in %d generation ',uint8(100*best_fitness/89), Ngen);
    title(title_str)
    lh = legend([h1 h2],'Food cell','Ant movement');
    set(lh,'Location','SouthEast');
end

    
    
    
    
    



