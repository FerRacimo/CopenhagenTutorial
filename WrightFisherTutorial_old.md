Exercises using the Wright-Fisher application
===============

Download and open the Wright-Fisher application from www.coalescent.dk. It is possible to use the Wright-Fisher animator to follow the reproduction process forward in time and the coalescent process backwards in time one generation at a time. After you have followed the reproduction process a number of generations forwards in time, it is possible to “untangle” the genealogy, and then to follow both how many descendants each of the original genes leave over generations (click on upper row), and to follow the ancestors to the sequence in the bottom row (by pressing the circles in the bottom row).

A new simulation is done by setting the parameters and pressing the new bottom. The simulation can then be controlled by the buttons in the bottom-right part of the window. The left-most button in the bottom-right corner enables you to untangle the resulting genealogy (i.e. rearranging individuals so that lines do not cross).

Set N = number of gene copies = 10, and G = number of generations = 15.

## 1 - Thinking forwards in time

If you click on an individual in the upper-most row, you can see the entire genealogy of descendants of that particular individual. The number of descendants of an individual "X" at a particular generation is just the number of individuals in that generation that can trace their ancestry back to "X". 

- a) Calculate the mean and variance in the number of descendants of the first-generation individuals in the second generation (second top-most row).
- b) Calculate the mean and variance in the number of descendants of the first-generation individuas in the present-day (bottom-most) generation.

If you need the formula for the sample mean and sample variance, you can find them here: http://www.math.uah.edu/stat/sample/Variance.html 

Does the sample mean differ by a lot in exercises a) and b)? What about the variance?


## 2- Thinking backwards in time (coalescence)

- a) Start a new simulation. Before you untangle the genealogy, choose 3 gene copies at random in the present generation. Do all these copies in the present generation coalesce within this time span? If so, record the number of generations to coalescence (the first time they have a common ancestor) of all 3 gene copies. If not, record that the copies did not coalesce and write down how many ancestors there are in the first generation.
- b) Repeat exercise a 4 times, writing down the number of ancestors and the time to coalescence. This should give you an idea about the variance of the process.
- c) Repeat exercises a & b above, but with N=5 and with N=30. How does the time to coalescence for the 3 gene copies in the present scale with N? In other words, do lineages tend to coalesce faster when N is large or when N is small?
