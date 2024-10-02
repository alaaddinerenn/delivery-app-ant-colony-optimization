import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;
/**
 * Program to find the quickest route for a delivery car with 2 different methods. Brute force method and ant colony optimization approach
 * @author Alaaddin Eren Namlı
 * @since Date: 11.05.2024
 */
public class Main {
    public static double minDistance1 = Double.MAX_VALUE; // global min distance variable
    public static void main(String[] args) throws FileNotFoundException {

        int chosenMethod = 2; // Chosen method variable to use different methods. Set 1 for brute force method, set 2 for ant colony optimization approach
        int chosenMapType = 2; // Type of the output map. Set 1 for shortest path, set 2 for pheromone intensities graph
        String fileName = "input04.txt"; // Name of the file

        // Hyperparameters of the ant colony optimization approach
        double iteration; // Number of iterations
        double ants; // Number of ants per iteration
        double degradationFactor; // Pheromone degradation factor
        double alpha; // Power of the pheromone level in the edge value formula
        double beta; // Power of the distance in the edge value formula
        double initialPheromoneIntensity;
        double q; // Constant value for calculating delta

        if (fileName.equals("input01.txt")) {
            iteration = 100; // Number of iterations
            ants = 100; // Number of ants per iteration
            degradationFactor = 0.9; // Pheromone degradation factor
            alpha = 0.7; // Power of the pheromone level in the edge value formula
            beta = 1.6; // Power of the distance in the edge value formula
            initialPheromoneIntensity = 0.1;
            q = 0.0001; // Constant value for calculating delta
        } else if (fileName.equals("input02.txt")) {
            iteration = 100; // Number of iterations
            ants = 100; // Number of ants per iteration
            degradationFactor = 0.9; // Pheromone degradation factor
            alpha = 0.7; // Power of the pheromone level in the edge value formula
            beta = 1.6; // Power of the distance in the edge value formula
            initialPheromoneIntensity = 0.1;
            q = 0.0001; // Constant value for calculating delta
        } else if (fileName.equals("input03.txt")) {
            iteration = 100; // Number of iterations
            ants = 100; // Number of ants per iteration
            degradationFactor = 0.9; // Pheromone degradation factor
            alpha = 0.7; // Power of the pheromone level in the edge value formula
            beta = 1.5; // Power of the distance in the edge value formula
            initialPheromoneIntensity = 0.1;
            q = 0.0001; // Constant value for calculating delta
        } else if (fileName.equals("input04.txt")) {
            iteration = 100; // Number of iterations
            ants = 100; // Number of ants per iteration
            degradationFactor = 0.9; // Pheromone degradation factor
            alpha = 1.05; // Power of the pheromone level in the edge value formula
            beta = 1.05; // Power of the distance in the edge value formula
            initialPheromoneIntensity = 0.1;
            q = 0.0001; // Constant value for calculating delta
        } else {
            iteration = 100; // Number of iterations
            ants = 100; // Number of ants per iteration
            degradationFactor = 0.9; // Pheromone degradation factor
            alpha = 1.0; // Power of the pheromone level in the edge value formula
            beta = 3.5; // Power of the distance in the edge value formula
            initialPheromoneIntensity = 0.1;
            q = 0.0001; // Constant value for calculating delta
        }

        File file = new File(fileName); // Create file

        if (!file.exists()) { // If file does not exist, exit the program
            System.out.printf("%s can not be found.", fileName); // Prints error message
            System.exit(1); // Exit the program
        }

        Scanner inputFile = new Scanner(file); // Create scanner object
        ArrayList<Double> xCoordinates = new ArrayList<>(); // Arraylist of x coordinates
        ArrayList<Double> yCoordinates = new ArrayList<>(); // Arraylist of y coordinates

        while (inputFile.hasNextLine()) { // Reads the file
            String line = inputFile.nextLine(); // Reads the next line
            String[] lineParts = line.split(","); // Split the line
            double xCoordinate = Double.parseDouble(lineParts[0]); // First element is the x coordinate
            double yCoordinate = Double.parseDouble(lineParts[1]); // Second element is the x coordinate
            xCoordinates.add(xCoordinate);
            yCoordinates.add(yCoordinate);
        }
        inputFile.close(); // Close the file

        if (chosenMethod == 1) { // If chosen method is 1 use the brute force method
            bruteForceMethod(xCoordinates,yCoordinates); // Brute force method
        } else { // else use the ant colony optimization approach
            antColonyOptimization(iteration, ants, degradationFactor,alpha, beta, initialPheromoneIntensity, q,xCoordinates, yCoordinates,chosenMapType); // ant colony optimization approach
        }
    }

    /**
     * Implements ant colony optimization approach
     * @param iteration Number of iterations
     * @param ants Number of ants per iteration
     * @param degradationFactor Pheromone degradation factor
     * @param alpha Power of the pheromone level in the edge value formula
     * @param beta Power of the distance in the edge value formula
     * @param initialPheromoneIntensity Initial pheromone intensities
     * @param q  Constant value for calculating delta
     * @param xCoordinates Arraylist of x coordinates
     * @param yCoordinates Arraylist of y coordinates
     * @param chosenMapType Type of the output map
     */
    public static void antColonyOptimization(double iteration, double ants, double degradationFactor, double alpha, double beta, double initialPheromoneIntensity, double q, ArrayList<Double> xCoordinates, ArrayList<Double> yCoordinates,int chosenMapType) {

        double[][] pheromoneIntensitiesMatrix = new double[xCoordinates.size()][xCoordinates.size()]; // Matrix to store pheromone levels
        // This for loop fills matrix with initial pheromone intensity
        // There is no pheromone from the node n to the node n, we use 0.0
        for (int i = 0; i < pheromoneIntensitiesMatrix.length; i++) {
            for (int j = 0; j < pheromoneIntensitiesMatrix.length; j++) {
                if (i != j) {
                    pheromoneIntensitiesMatrix[i][j] = initialPheromoneIntensity;
                } else {
                    pheromoneIntensitiesMatrix[i][j] = 0.0;
                }
            }
        }

        ArrayList<Integer> nodes1 = new ArrayList<>(); // Arraylist of nodes
        for (int i = 0; i < xCoordinates.size(); i++) { // Creates the node arraylist
            nodes1.add(i);
        }

        ArrayList<ArrayList<Integer>> bestRoutes = new ArrayList<>(); // Arraylist to store routes with minimum distances
        ArrayList<Double> shortestDistances = new ArrayList<>(); // Arraylist to store minimum distances

        double initialTime = System.currentTimeMillis(); // Get the time before the process

        for (int iterations = 0; iterations < iteration ; iterations++) { // For loop for iterations

            ArrayList<Integer> bestRoute = new ArrayList<>(); // Arraylist to store the route with minimum distance
            double minDistance = Double.MAX_VALUE;
            Random random = new Random(); // Create random object to get random numbers

            for (int ant = 0; ant < ants; ant++) { // For loop for ants

                int randomSourceIndex = random.nextInt(0, xCoordinates.size()); // Get a random source index to start the route
                ArrayList<Integer> visitedNodes = new ArrayList<>(); // Arraylist of visited nodes, to keep track of the nodes on the route
                visitedNodes.add(randomSourceIndex); // Add the first node, random source node, to the visited nodes
                ArrayList<Integer> route = new ArrayList<>(); // Arraylist of nodes to store the route
                route.add(randomSourceIndex);
                double totalCycleDistance = 0;

                for (int vns = visitedNodes.size(); vns < nodes1.size(); vns++) { // For loop for 1 ant, 1 cycle

                    ArrayList<Double> distances = new ArrayList<>(); // Arraylist to store distances from 1 node to the others
                    for (int i = 0; i < nodes1.size(); i++) {
                        if (!visitedNodes.contains(i)) { // If the node is visited, add 0.0 else calculate the distance between the nodes
                            double distance = Math.sqrt(Math.pow(xCoordinates.get(randomSourceIndex) - xCoordinates.get(i), 2) + Math.pow(yCoordinates.get(randomSourceIndex) - yCoordinates.get(i), 2)); // Calculate the distance
                            distances.add(distance);
                        } else {
                            distances.add(0.0);
                        }
                    }

                    ArrayList<Double> edgeValues = new ArrayList<>(); // Arraylist to store edge values of the nodes
                    double sumOfEdgeValues = 0;
                    for (int i = 0; i < nodes1.size(); i++) { // If the node is visited add 0.0 else calculate the edge value between the nodes
                        if (!visitedNodes.contains(i)) {
                            double edgeValue = (Math.pow(pheromoneIntensitiesMatrix[randomSourceIndex][i], alpha)) / (Math.pow(distances.get(i), beta)); // Calculate the edge value
                            sumOfEdgeValues = sumOfEdgeValues + edgeValue; // Sum the edge values
                            edgeValues.add(edgeValue);
                        } else {
                            edgeValues.add(0.0);
                        }
                    }

                    ArrayList<Double> probabilities = new ArrayList<>(); // Arraylist to store probabilities of the nodes
                    // Calculate the probability of the node by dividing the edge value of the node by the sum of edge values
                    for (int i = 0; i < nodes1.size(); i++) { // If the node is visited add 0.0 else calculate the probabilities of the nodes
                        if (!visitedNodes.contains(i)) {
                            double probability = edgeValues.get(i) / sumOfEdgeValues; // Calculate the probability
                            probabilities.add(probability);
                        } else {
                            probabilities.add(0.0);
                        }
                    }

                    // To choose which node to go to, we get a random number between 0 and 1
                    // Then we divide the space between 0 and 1 into regions
                    // The region in which the random number is located is the next node
                    double randomNumber = random.nextDouble(); // Create the random number
                    double sumOfProbabilities = 0;
                    int nextNodeIndex = 0;

                    for (int i = 0; i < nodes1.size(); i++) {
                        if (!visitedNodes.contains(i)) {
                            sumOfProbabilities = sumOfProbabilities + probabilities.get(i); // Add probabilities
                            if (randomNumber <= sumOfProbabilities) {
                                nextNodeIndex = i;
                                break; // When the next node is found, exit the loop
                            }
                        }
                    }

                    visitedNodes.add(nextNodeIndex); // Add the next node to the visited nodes
                    randomSourceIndex = nextNodeIndex; // Update the node
                    route.add(nextNodeIndex); // Add the next node to the route
                    totalCycleDistance = totalCycleDistance + distances.get(nextNodeIndex); // Add the distances to get total cycle distance
                 }

                // At the end of the cycles

                // If the total distance of the cycle is less than min distance, keep the route and the total cycle distance
                if (totalCycleDistance < minDistance) {
                    minDistance = totalCycleDistance;
                    bestRoute = (ArrayList<Integer>) route.clone();
                }

                // This for loop updates the pheromone intensities
                double delta = q / totalCycleDistance; // Delta is the amount that we update the pheromone intensities
                for (int i = 0; i < route.size() - 1; i++) {
                    int j = route.get(i); // Get the first index
                    int k = route.get(i + 1); // Get the second index
                    pheromoneIntensitiesMatrix[j][k] = pheromoneIntensitiesMatrix[j][k] + delta;
                    pheromoneIntensitiesMatrix[k][j] = pheromoneIntensitiesMatrix[k][j] + delta;
                }
            }

            // At the end of each iteration

            // Add the distance between the first node and second node of the route to get the distance of the cycle
            int index1 = bestRoute.getFirst(); // Get first node
            int index2 = bestRoute.getLast(); // Get last node
            double d = Math.sqrt(Math.pow(xCoordinates.get(index1) - xCoordinates.get(index2),2) + Math.pow(yCoordinates.get(index1) - yCoordinates.get(index2),2)); // Calculate the distance
            minDistance = minDistance + d; // Update the distance

            // Re-order the route so that starting node is 0th index
            ArrayList<Integer> tempArraylist = new ArrayList<>(); // Create temporary arraylist
            int sourceIndex = bestRoute.indexOf(0);
            for (int i = sourceIndex; i < bestRoute.size() ; i++) { // This for loop adds indexes after the 0 index
                int value = bestRoute.get(i);
                tempArraylist.add(value);
            }
            for (int i : bestRoute) { // This for loop adds indexes before the 0 index
                if (!tempArraylist.contains(i)) {
                    tempArraylist.add(i);
                }
            }

            // Increment the indexes so that we get the node numbers
            ArrayList<Integer> tempArraylist2 = new ArrayList<>(); // Create temporary arraylist
            for (int i = 0; i < bestRoute.size(); i++) {
                tempArraylist2.add(tempArraylist.get(i) + 1); // Increment the value by 1
            }
            tempArraylist2.addLast(1); // Complete the cycle
            bestRoute = (ArrayList<Integer>) tempArraylist2.clone();

            bestRoutes.add(bestRoute); // Add the best route of the all ants
            shortestDistances.add(minDistance); // Add the minimum distance off all ants' cycles

            // Update the pheromone intensities with degradation factor
            for (int i = 0; i < pheromoneIntensitiesMatrix.length; i++) {
                for (int j = 0; j < pheromoneIntensitiesMatrix.length; j++) {
                    if (i != j) {
                        pheromoneIntensitiesMatrix[i][j] = pheromoneIntensitiesMatrix[i][j] * degradationFactor; // Update the pheromone level
                        pheromoneIntensitiesMatrix[j][i] = pheromoneIntensitiesMatrix[j][i] * degradationFactor; // Update the pheromone level
                    } else {
                        pheromoneIntensitiesMatrix[i][j] = 0.0;
                    }
                }
            }
        }

        double minDistance2 = Double.MAX_VALUE;
        int minIndex = -1;
        // Iterate through all relatively short distances and find the shortest one
        for (int i = 0; i < shortestDistances.size(); i++) {
            if (shortestDistances.get(i) < minDistance2) { // If the distance is shorter than current shortest distance
                minDistance2 = shortestDistances.get(i); // Update the min distance
                minIndex = i; // Update the index of minimum distance
            }
        }
        double lastTime = System.currentTimeMillis(); // Get the time after the process finished
        double time = (lastTime - initialTime) / 1000; // Calculate the time it takes to find the shortest path in seconds

        ArrayList<Integer> bestPath = bestRoutes.get(minIndex); // Get the shortest path
        // Print the messages on the console
        System.out.println("Method: Ant Colony Optimization");
        System.out.printf(Locale.US,"Shortest Distance: %.5f\n",shortestDistances.get(minIndex));
        System.out.println("Shortest Path: "+bestPath);
        System.out.printf(Locale.US,"Time it takes to find the shortest path: %.2f seconds.",time);

        // STD DRAW
        StdDraw.setXscale(0,1);
        StdDraw.setYscale(0,1);
        StdDraw.enableDoubleBuffering(); // Enables double buffering for a faster process

        if (chosenMapType == 1) { // If the chosen map type is the shortest path type
            for (int i = 0; i < bestPath.size() - 1; i++) { // This for loop draws lines
                int j = bestPath.get(i) - 1; // Get first index
                int k = bestPath.get(i + 1) - 1; // Get second index
                StdDraw.setPenRadius(0.007);
                StdDraw.line(xCoordinates.get(j), yCoordinates.get(j), xCoordinates.get(k), yCoordinates.get(k)); // Draw the line between nodes
            }

            for (int i = 1; i < xCoordinates.size(); i++) { // This for loop draws circles
                StdDraw.setPenColor(StdDraw.GRAY);
                StdDraw.filledCircle(xCoordinates.get(i), yCoordinates.get(i), 0.025); // Draws the circles of the nodes
                StdDraw.setPenColor(StdDraw.BLACK);
                StdDraw.text(xCoordinates.get(i), yCoordinates.get(i), String.valueOf(i + 1)); // Writes the number of the circles
            }
            StdDraw.setPenColor(StdDraw.PRINCETON_ORANGE); // Set color to princeton orange for migros
            StdDraw.filledCircle(xCoordinates.getFirst(), yCoordinates.getFirst(), 0.025); // Draws the circle of the migros
            StdDraw.setPenColor(StdDraw.BLACK);
            StdDraw.text(xCoordinates.getFirst(), yCoordinates.getFirst(), "1"); // Writes "1" on the migros circle

        } else { // If the chosen map type is the pheromone intensities graph type

            // This nested for loop creates the pheromone intensities graph
            for (int i = 0; i < pheromoneIntensitiesMatrix.length; i++) {
                for (int j = 0; j < pheromoneIntensitiesMatrix.length; j++) {
                    double penRadius = pheromoneIntensitiesMatrix[i][j] * 0.7; // Update pen radius according to pheromone intensities
                    StdDraw.setPenRadius(penRadius); // Set pen radius according to pheromone intensities
                    double xCoordinate1 = xCoordinates.get(i);
                    double xCoordinate2 = xCoordinates.get(j);
                    double yCoordinate1 = yCoordinates.get(i);
                    double yCoordinate2 = yCoordinates.get(j);
                    StdDraw.line(xCoordinate1,yCoordinate1,xCoordinate2,yCoordinate2); // Draws the lines between the nodes
                }
            }

            for (int i = 0; i < xCoordinates.size(); i++) { // This for loop draws circles
                StdDraw.setPenColor(StdDraw.GRAY);
                StdDraw.filledCircle(xCoordinates.get(i), yCoordinates.get(i), 0.025); // Draws the circles of the nodes
                StdDraw.setPenColor(StdDraw.BLACK);
                StdDraw.text(xCoordinates.get(i), yCoordinates.get(i), String.valueOf(i + 1)); // Writes the number of the circles
            }
        }
        StdDraw.show(); // Show the drawing
    }

    /**
     * Implements the brute force method
     * @param xCoordinates Arraylist of x coordinates
     * @param yCoordinates Arraylist of y coordinates
     */
    public static void bruteForceMethod(ArrayList<Double> xCoordinates,ArrayList<Double> yCoordinates) {

        int[] initialRoute = new int[xCoordinates.size()-1]; // Create the initial route array
        for (int i = 0; i < xCoordinates.size()-1; i++) { // Fill the array with nodes
            initialRoute[i] = i+1;
        }

        double initialTime = System.currentTimeMillis(); // Get the time before the process
        ArrayList<int[]> result = new ArrayList<>(); // Initialize result arraylist
        int[] tempArray = null; // Initialize temporary array
        permutationsAndDistances(result,initialRoute,0,xCoordinates,yCoordinates); // Get the result of the permutations
        // Shortest path is in the result arraylist
        tempArray = result.getFirst(); // Get the shortest path
        double shortestDistance = calculateDistance(tempArray,xCoordinates,yCoordinates); // Calculate the shortest distance
        double lastTime = System.currentTimeMillis(); // Get the time after the process finished
        for (int i = 0; i < tempArray.length; i++) { // Add 1 to indexes in shortest path array to get number of the nodes
            tempArray[i] = tempArray[i] + 1;
        }
        ArrayList<Integer> bestRoute = ArrayToArraylist(tempArray); // Convert temporary array to arraylist and get the shortest route arraylist
        bestRoute.addFirst(1); // Add 1 to the first index of the shortest route arraylist
        bestRoute.addLast(1); // Add 1 to the last index of the shortest route arraylist
        double time = (lastTime - initialTime) / 1000; // Calculate the time it takes to find the shortest path in seconds
        // Print the messages on the console
        System.out.println("Method: Brute Force Method");
        System.out.printf(Locale.US,"Shortest Distance: %.5f\n",shortestDistance);
        System.out.println("Shortest Path: "+bestRoute);
        System.out.printf(Locale.US,"Time it takes to find the shortest path: %.2f seconds.",time);

        //STD DRAW
        StdDraw.setXscale(0,1);
        StdDraw.setYscale(0,1);

        for (int i = 0; i < bestRoute.size()-1; i++) { // This for loop draws lines
            int j = bestRoute.get(i) - 1; // Get first index
            int k = bestRoute.get(i+1) - 1; // Get second index
            StdDraw.setPenRadius(0.007);
            StdDraw.line(xCoordinates.get(j),yCoordinates.get(j),xCoordinates.get(k),yCoordinates.get(k));
        }

        for (int i = 1; i < xCoordinates.size(); i++) { // This for loop draws circles
            StdDraw.setPenColor(StdDraw.GRAY);
            StdDraw.filledCircle(xCoordinates.get(i),yCoordinates.get(i),0.025); // Draws the circles of the nodes
            StdDraw.setPenColor(StdDraw.BLACK);
            StdDraw.text(xCoordinates.get(i),yCoordinates.get(i),String.valueOf(i+1)); // Writes the number of the circles
        }
        StdDraw.setPenColor(StdDraw.PRINCETON_ORANGE); // Set color to princeton orange for migros
        StdDraw.filledCircle(xCoordinates.getFirst(),yCoordinates.getFirst(),0.025); // Draws the circle of the migros
        StdDraw.setPenColor(StdDraw.BLACK);
        StdDraw.text(xCoordinates.getFirst(),yCoordinates.getFirst(),"1"); // Writes "1" on the migros circle
        StdDraw.show(); // Show the drawing
    }

    /**
     * A recursive method that finds all the permutations of the nodes and calculate the distances of all path to get the shortest route
     * @param result Arraylist that contains shortest route array
     * @param route Array that contains the nodes on the route
     * @param c Counter
     * @param xCoordinates Arraylist of x coordinates
     * @param yCoordinates Arraylist of y coordinates
     */
    public static void permutationsAndDistances(ArrayList<int[]> result, int[] route, int c, ArrayList<Double> xCoordinates,ArrayList<Double> yCoordinates) {
        if (c == route.length) {
            double distance = calculateDistance(route,xCoordinates,yCoordinates); // Calculate the distance
            if (distance < minDistance1) { // Compare the distance with the current minimum distance
                result.clear();
                minDistance1 = distance; // Update the minimum distance
                result.add(route.clone()); // Add the shortest route
            }
        } else {
            for (int i = c; i < route.length; i++) {
                // Swap the value at index i with the one at index c
                int temp = route[i];
                route[i] = route[c];
                route[c] = temp;
                permutationsAndDistances(result, route, c + 1, xCoordinates, yCoordinates); // Recursive call
                // Again swap the value at index i with the one at index c
                temp = route[c];
                route[c] = route[i];
                route[i] = temp;
            }
        }
    }

    /**
     * Calculates the total distance of the cycle
     * @param route Array that contains the nodes on the route
     * @param xCoordinates Arraylist of x coordinates
     * @param yCoordinates Arraylist of y coordinates
     * @return The total distance of the cycle
     */
    public static double calculateDistance(int[] route,ArrayList<Double> xCoordinates,ArrayList<Double> yCoordinates) {
        int migrosIndex = 0;
        int previousIndex = migrosIndex;
        double totalDistance = 0;

        for (int j = 0; j < route.length; j++) {
            double xCoordinate1 = xCoordinates.get(route[j]); // Get first x coordinate
            double yCoordinate1 = yCoordinates.get(route[j]); // Get first y coordinate
            double xCoordinate2 = xCoordinates.get(previousIndex); // Get second x coordinate
            double yCoordinate2 = yCoordinates.get(previousIndex); // Get second y coordinate
            double distance = Math.sqrt(Math.pow(xCoordinate1-xCoordinate2,2) + Math.pow(yCoordinate1-yCoordinate2,2)); // Calculate the distance
            totalDistance = totalDistance + distance; // Add distances
            previousIndex = route[j]; // Update the previous node
        }
        double migrosX = xCoordinates.get(0);
        double migrosY = yCoordinates.get(0);
        double lastX = xCoordinates.get(route[route.length-1]); // Get the last x coordinate
        double lastY = yCoordinates.get(route[route.length-1]); // Get the last y coordinate
        double lastDistance = Math.sqrt(Math.pow(migrosX-lastX,2) + Math.pow(migrosY-lastY,2)); // Calculate the distance between the migros and the last node
        totalDistance = totalDistance + lastDistance; // Add the distance
        return totalDistance;
    }

    /**
     * Converts arrays to arraylists
     * @param array An array
     * @return Arraylist version of the array
     */
    public static ArrayList<Integer> ArrayToArraylist(int[] array) {
        ArrayList<Integer> arrayList = new ArrayList<>(); // Initialize arraylist
        for (int i : array) { // Add all elements of the array to the arraylist
            arrayList.add(i);
        }
        return arrayList;
    }

}
