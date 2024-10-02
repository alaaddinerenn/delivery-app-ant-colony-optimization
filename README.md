# DELIVERY APP

## What It Does
This program finds the quickest route for a delivery car by using 2 different approaches according to coordinates of the input delivery locations.
Taking the advantage of the **Ant Colony Optimization** approach, this application finds the best route more than **941 times** faster as the input size increases, compared to the **brute force** approach.

## Report
You can check the report file for the hyper parameters used with the sample inputs for the Ant Colony Optimization approach. Also, you can find the comparison of the approaches in the report.  

## Features
- **Fast:** Thanks to the **Ant Colony Optimization** approach, this program finds nearly the best route faster than the brute force approach.

- **Comparability:** You can switch between two approach and compare their running times with chosenMethod variable. Set 1 for brute force method, set 2 for ant colony optimization approach.

- **Different Input Files** You can switch between the sample input files with the fileName variable to observe the running times of the two approach for different input sizes.

- **Timer:** You can see the running time of the program on the console.

- **Visualization** By using the StdDraw library, delivery locations, the warehouse and the best route are drawn.  

- **Map Types** You can switch between the types of maps for the **Ant Colony Optimization** approach with the chosenMapType variable. Set 1 for the shortest path, set 2 for the pheromone intensities graph.

## Getting Started

### Prerequisites
Before starting, ensure that you have installed the StdDraw library. Also, you need to have the input files in the input-files folder.

### Usage
When the library is installed, you can start to use the program. Give the coordinates of the warehouse and the delivery locations in an input file and choose the approach and the map type. Then, the program finds the best route and draws your graph. 
