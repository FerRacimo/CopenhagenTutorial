Exercises using the Hudson animator - Copenhagen 2017 course
===============

Download and open the Hudson animator program from www.coalescent.dk. It is developed by Anders M. Mikkelsen, Jotun Hein and Mikkel Schierup as a tool for the visualization of the following continuous time processes.
 
- The basic coalescent and the Coalescent with recombination
- Coalescent with exponential growth
- Coalescent with migration
 
Please consult the manual before doing the exercises, it can be found under help at the start page. It briefly describes how to control the applet.
 
## The basic Coalescent
 
Choose coalescent with recombination (we will set recombination to zero).

a)	The basic coalescent. Choose n=5 sequences and rho=0 (no recombination). Press recalc and the animation starts. The speed can be controlled with speed. Details regarding each node can be seen in the small window at the right when moving the mouse over the node.
- What is the time to the first coalescence event? Write it down.
- What is the time to the most recent common ancestor? Write it down.
- Repeat a and b 5 times (by pressing recalc). How does the time to first coalescence and time to most recent common ancestor vary?

b)	Try with 10 and 20 sequences. What are the times to the first coalescence and the most recent common ancestor in these cases? Write it down for 5 different runs.


## Coalescent with exponential growth
 
Now choose coalescent with exponential growth. This is controlled with the parameter exp, which is equal to Nb. This parameter measure how many times the present population is larger than the population 2N (N=size of present day population) ago. In studies of human mitochondria (there is no recombination in mitochondria) all estimates suggest that exp>100.
 
c)	Try to simulate n=10 sequences and different runs with exp=1, 10, 100, 1000
- How does the shape of the genealogical tree depend on the value of exp?
- How can that be?
- How would this altered shape be visible in a set of sequences evolving on the tree? Would there be fewer or more “singletons”? How would Tajima’s D be affected?
Hint: If you push trees the same tree will appear without any crossing branches
 
