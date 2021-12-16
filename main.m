clear;
UI();

%This is the function for the UI which sets the default paramters, and asks
%the user if they want to run the default parameters, or try out different
%ones.
function UI()
    %Set default values
    selectionType = 1; %1= roulette, 2=tournament, 3 = truncation
    crossoverType = 2; %1 = 1k-point, 2= 2k-point, 3 = uniform
    mutationType = 2; %1 = swap mutation, 2 = random mutation, 3 = scramble mutation
    iterations = 2000;
    populationSize = 100;
  
    option = input('Select 1 for default options, select 2 for custom ');
    if option == 1
        fprintf('Running with default values: \n Selection type: roulette \n Crossover type: 2k-point \n Mutation type: Random \n Generations: 2000 \n Population size: 100\n');
        runProgram(populationSize, iterations, selectionType, crossoverType, mutationType);
    elseif option == 2
        selectionType = input('Selection type: 1 for roulette wheel, 2 for tournament, 3 for truncation ');
        crossoverType = input('Crossover type: 1 for 1k-point, 2 for 2k-point, 3 for uniform ');
        mutationType = input('Mutation type: 1 for swap mutation, 2 for random mutation, 3 for scramble mutation ');
        iterations = input('Please enter the number of generations you want to simulate ');
        populationSize = input('Please enter the population size of each generation ');
        
        runProgram(populationSize, iterations, selectionType, crossoverType, mutationType);
    end
end

%This is the main loop function which runs the simulation 'iterations'
%number of times, it also calls the functions for graphing the results.
function runProgram(population_size, iterations, selectionType, crossoverType, mutationType)
    currentPopulation = generateInitialPopulation(population_size); %Generates first population
    fitness_data = zeros(iterations, 1);
    tic
    for i = 1:iterations
        [currentPopulation(:, 1:30), maxFitness] = createNewPopulation(currentPopulation, population_size, selectionType, crossoverType, mutationType);
        fitness_data(i) = maxFitness;
    end
    toc %Average of 10 runs at default values: elapsed time = 6.24 seconds, fitness = 93%
    [finalPop, trail] = calculateFitness(currentPopulation, dlmread('muir_world.txt',' '), population_size);
    maxFitness = max(finalPop(:,end));
    fitness_data(iterations) = maxFitness; %Update final fitness score
    %Call graphing functions
    plotFitnessScore(fitness_data, iterations);
    plotBestTrail(maxFitness, iterations, trail);
end

%This function generates the intiial population represented as a
%(population_sizex30) matrix with randomly generated genes, in line with
%the FSM controller specification
function original_population = generateInitialPopulation(population_size)
    original_population = zeros(population_size, 30);
    for i = 1:population_size
        chromosome = zeros(1, 30);
        for k = 1:30/3
            chromosome((k*3)-2) = randi([1, 4]);
            chromosome((k*3)-1) = randi([0, 9]);
            chromosome(k*3) = randi([0, 9]);
        end
        original_population(i, :) = chromosome(:);
    end
    newColumn = zeros(size(original_population,1),1);
    original_population = [original_population, newColumn];
end

%This function runs the simulate_ant function for each chromosome in the
%population. It appends the population matrix to include the fitness of
%each chromosome. It returns the appended population and the trail (used for
%graphing)
function [population, trail] = calculateFitness(population, worldGrid, population_size)
    for i = 1:population_size
        [fitness, trail] = simulate_ant(worldGrid, population(i, 1:31));
        population(i, 31) = fitness;
    end
    population = sortrows(population, 31);
end



%This is the function which generates a new population for each generation. 
%It calls the methods required to select, crossover, and mutate the
%chromosomes depending on the user choice/default options.
%It returns the new population as a matrix and fitness value of the fittest
%ant of the population.
function [newPopulation, maxFitness] = createNewPopulation(currentPopulation, population_size, selectionChoice, crossoverChoice, mutationChoice)
    %Find fitness score of population
    [currentPopulation, ~] = calculateFitness(currentPopulation, dlmread('muir_world.txt',' '), population_size);
    maxFitness = max(currentPopulation(:, end));
    %Elitism - Finds best 2 chromosomes to keep in new population
    newPopulation = zeros(population_size, 30);
    newPopulation(1:2, :) = currentPopulation(population_size-1:population_size, 1:30);
    populationNewNum = 2;
    while (populationNewNum < population_size) %Do until population is filled out
        %Selection process
        if selectionChoice == 1 %Roulette wheel selection
            weights = currentPopulation(:,31)/sum(currentPopulation(:,31)); %Weights for roulette wheel selection
            choice1 = rouletteWheelSelection(weights);
            choice2 = rouletteWheelSelection(weights);
        elseif selectionChoice == 2 %Tournament selection
            choice1 = TournamentSelection(currentPopulation(:,31));
            choice2 = TournamentSelection(currentPopulation(:,31));
        elseif selectionChoice == 3 %Truncation selection
            multiplier = 0.1; %Takes top 10% of ants
            choice1 = TruncationSelection(currentPopulation(:,31), multiplier);
            choice2 = TruncationSelection(currentPopulation(:,31), multiplier);
        end
        %Gets parent chromsomes from population with index from selection
        %algorithm
        temp_chromosome1 = currentPopulation(choice1,1:30);
        temp_chromosome2 = currentPopulation(choice2,1:30);
        %Crossover process
        if (rand < 0.8)
            if crossoverChoice == 1
                [temp_chromosome1, temp_chromosome2] = onePointCrossover(temp_chromosome1, temp_chromosome2);
            elseif crossoverChoice == 2
                [temp_chromosome1, temp_chromosome2] = twoPointCrossover(temp_chromosome1, temp_chromosome2);
            elseif crossoverChoice == 3
                [temp_chromosome1, temp_chromosome2] = uniformCrossover(temp_chromosome1, temp_chromosome2);
            end
        end
        %Mutation process
        if mutationChoice == 1
            if (rand < 0.4)
                temp_chromosome1 = swapMutation(temp_chromosome1);
            end
            if (rand < 0.4)
                temp_chromosome2 = swapMutation(temp_chromosome2);
            end
        elseif mutationChoice == 2
            if (rand < 0.4)
                temp_chromosome1 = randomMutation(temp_chromosome1);
            end
            if (rand < 0.4)
                temp_chromosome2 = randomMutation(temp_chromosome2);
            end
        elseif mutationChoice == 3
            if (rand < 0.4)
                temp_chromosome1 = scrambleMutation(temp_chromosome1);
            end
            if (rand < 0.4)
                temp_chromosome2 = scrambleMutation(temp_chromosome2);
            end
        end
        %Adds child chromosomes to population
        populationNewNum = populationNewNum + 1;
        newPopulation(populationNewNum,:) = temp_chromosome1;
        if (populationNewNum < population_size)
            populationNewNum = populationNewNum + 1;
            newPopulation(populationNewNum,:) = temp_chromosome2;
        end
    end
end

%Selection option 1 -> Roulette wheel selection, this is a fitness
%proportionate technique which weighs each ant on their fitness, and is
%more likely to choose one with a higher weight. It returns the index of
%the chosen ant
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

%Selection option 2 -> Tournament selection, this function selects two ants
%at random, and compares their fitness to see which one is higher, it then
%chooses the ant with the higher fitness
function choice = TournamentSelection(fitness)
    [fitness, index] = sortrows(fitness);
    chr1 = randi([1, length(fitness)]);
    chr2 = randi([1, length(fitness)]);
    while chr1 == chr2 %Ensures they are not the same
        chr2 = randi([1, length(fitness)]);
    end
    if(fitness(chr1)> fitness(chr2)) %Chooses higher fitness index
        choice = index(chr1);
    else
        choice = index(chr2);
    end
 end


%Selection Option 3 -> Truncation selection, this function truncates the
%overall pool of ants by applying a multiplier (E.g. 0.1, truncating to
%10%), meaning that only the top (e.g. 10%) of ants are available, it then chooses
%at random an ant from this truncated list, and returns its index.
function choice=TruncationSelection(fitness, multiplier)
    number=round(length(fitness)- (length(fitness)* multiplier));
    choice=randi([number, length(fitness)]);
end

%Crossover option 1 -> 1k point crossover - This function chooses a random
%point of the genetic material of the ants, and then swaps the genes of the
%parents after this point, it returns the child ants.
function [temp_chromosome1, temp_chromosome2] = onePointCrossover(temp_chromosome1, temp_chromosome2)
    point = randi([1, ((30/3))])*3; %choosing between 1 and 10 and multiplying by 3 ensures the point is at a start of a state
    temp_chromosome1 = [temp_chromosome1(1:point), temp_chromosome2(point+1:end)];
    temp_chromosome2 = [temp_chromosome2(1:point), temp_chromosome1(point+1:end)];
end

%crossover option 2 -> 2k point crossover - This function chooses two
%points of the genetic material of the ants at random, and swaps the genes
%that are between these two points. it returns the child ants.
function [temp_chromosome1, temp_chromosome2] = twoPointCrossover(temp_chromosome1, temp_chromosome2)
  point1 = randi([1, ((30/3)-2)])*3; % -2 ensures that the first point is not at the end of the genetic material, with space for second point
  point2 = randi([(point1/3), ((30/3)-1)])*3;
  temp_chromosome1 = [temp_chromosome1(1:point1), temp_chromosome2((point1+1):point2), temp_chromosome1(point2+1:30)];
  temp_chromosome2 = [temp_chromosome2(1:point1), temp_chromosome1((point1+1):point2), temp_chromosome2(point2+1:30)];
end

%crossover option 3 -> Uniform crossover - This function chooses each
%state of the parent chromosomes at random to give to the children (With a
%50% chance of state passing on from each parent). It returns the child ants.
function [temp_chromosome1, temp_chromosome2] = uniformCrossover(temp_chromosome1, temp_chromosome2)
    for index = 1 : 3: 30
        if rand < 0.5
            temp = (temp_chromosome1(index:index+2));
            [temp_chromosome1(index:index+2)] = (temp_chromosome2(index:index+2));
            [temp_chromosome2(index:index+2)] = temp;
        end
    end
end

%Mutation option 1 - Swap mutation - This function takes the chromosome and
%randomly chooses two points of its genetic material. It then swaps these
%two points (Ensuring valid results i.e. first value of each state must be between 1-4,
%rest is 0-4) It returns the mutated chromosome.
function mutatedChromosome = swapMutation(chromosome)
    point1 = randi([1, length(chromosome)]);
    if mod(point1,3) == 1 %Ensures if a 1-4 gene is selected, it is only swapped with other 1-4 gene
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
    mutatedChromosome = chromosome;
end


%Mutation option 2 -> random reset - This function takes a random gene
%within the chromosome and replaces it with another random (valid) number.
%It returns the mutated chromosome
function mutatedChromosome = randomMutation(chromosome)
    digit1 = randi([1, 30]);
    if mod(digit1, 3) == 1
        chromosome(digit1) = randi([1,4]);
    else
        chromosome(digit1) = randi([0,9]);
    end
    mutatedChromosome = chromosome;
end


%Mutation option 3 -> scramble - This function groups the genes into a cell
%array of each state each containing the 3 gene digits, and then chooses two random points in the array, it
%then scrambles each element between those two points, and returns the
%scrambled chromosome as a matrix.
function mutatedChromosome = scrambleMutation(chromosome)
    c = cell(1,10);
    for i = 1:3:30
        newArray = chromosome(i:i+2);
        c{(i/3) + (2/3)} = newArray; %index is function mapping i to a 1-10 index
    end
    for k = 1:5 %Looped 5 times to scramble between different points in chromosomes
        r1 = randi([1,10]);
        r2 = randi([1,10]);
        while(r1>=r2)
            r1 = randi([1,10]);
            r2 = randi([1,10]);
        end
        for i = 1:10 %Repetead 10 times to ensure it is well scrambled
            i1 = randi([r1, r2]);
            i2 = randi([r1, r2]);
            a = c(i1);
            c(i1) = c(i2);
            c(i2) = a;
        end
    end
    mutatedChromosome = cell2mat(c); %Converts grouped cell array back to matrix
end


function [b, a] = swap(a, b)
end
    
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