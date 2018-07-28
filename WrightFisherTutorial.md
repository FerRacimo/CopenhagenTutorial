Exercises using the Wright-Fisher application
===============

Download and open the Wright-Fisher application from www.coalescent.dk. It is developed by Anders M. Mikkelsen, Jotun Hein and Mikkel Schierup. It is possible to use the Wright-Fisher animator to follow the reproduction process forward in time and the coalescent process backwards in time one generation at a time. After you have followed the reproduction process a number of generations forwards in time it is possible to “untangle” the genealogy, and then to follow both how many descendants each of the original genes leave over generations (click on upper row), and to follow the ancestors to the sequence in the bottom row (by pressing the circles in the bottom row).

A new simulation is done by setting the parameters and pressing the new bottom. The simulation can then be controlled by the buttons in the bottom (right) part of the window, e.g. one generation at a time. One button enables you to untangle the resulting genealogy (i.e. rearranging individuals so that lines do not cross).

Set N = number of gene copies = 10, and G = number of generations = 15.

## 1 - Thinking forwards in time

- a) For each of the initial gene copies, record how long they persist in the population.
- b) Calculate the mean and variance in the number of descendants of the initial gene copies in the next generation.
- c) Calculate the mean and variance in the number of descendants of the initial gene copies in the present (bottom) generation.

If you need the formula for the sample mean and sample variance, you can find them here: http://www.math.uah.edu/stat/sample/Variance.html 

## 2- Thinking backwards in time (coalescence)

- a) Start a new simulation. Before you untangle the genealogy, choose 3 gene copies at random in the present generation. Do all these copies in the present generation coalesce within this time span? If so, record the number of generations to coalescence (the first time they have a common ancestor) of all 3 gene copies. If not, record that the copies did not coalesce and write down how many ancestors there are in the first generation.
- b) Repeat exercise a 5 times, writing down the number of ancestors and the time to coalescence. This should give you an idea about the variance of the process.
- c) Repeat exercises a-b above with N=5 and N=30. How does the time to coalescence for the 3 gene copies in the present scale with N? In other words, do they tend to coalesce faster when N is large or when N is small?
