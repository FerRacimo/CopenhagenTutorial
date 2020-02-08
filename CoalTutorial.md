Exercises in coalescent theory
===============

Based on notes by Ida Moltke.

## Exercise	A:	Simulating	a	coalescence	tree	assuming	a	constant	population	size

The	purpose	of	this first	exercise	is	to	make	sure	it	is	clear	how	a	coalescence tree	is	simulated. We will use R so a little familiarity with this language will help. First, let	us try to	simulate	a	coalescence tree	for	five gene	copies by	hand:

1. Start	by	drawing	on a piece of paper a	small circle	for	each	of	the	five	gene	copies. They should be lined up on	an	invisible horizontal line	and you should leave enough	space	above the circles for	
drawing	a	tree	above them (which we will do shortly).	We will henceforth call these five circles "nodes" and label them 1,2,3,4,5

2. Also,	make	a	list	of	the	node	names.	You	can	either	do	this	by	hand	or	you	can	do	it	in	R	by	
simply	writing:

```
R
nodes = c(1,2,3,4,5) # make the list and call it nodes
nodes # print the list
```

3. Sample	which	two	nodes	will	coalesce	first	(going	back	in	time)	by	randomly	picking	two	of	the	
nodes.	You	can	either	do	this	by	hand	or	you	can	do	it	in	R	by	typing:

```
nodecount = length(nodes) # save the number of nodes in the variable nodecount
tocoalesce = sample(1:nodecount, size=2) # sample 2 different nodes in node list
nodes[tocoalesce[1]] # print the first node sampled
nodes[tocoalesce[2]] # print the second node sampled
```

If	you	used	R	then	make	sure	you	understand	what	the	R	code	does	before	moving	on.

4. Sample	the	time	it	takes	before	these	two	nodes	coalesce	(measured	from	previous	
coalescence	event in	units	of	2N)	by	sampling	from	an	exponential	distribution	with	rate	equal	
to	nodecount*(nodecount-1)/2	where	nodecount	is	the	number	of	nodes	in	your node	list.	Do	
this	in	R	by	typing:

```
coalescencerate = nodecount*(nodecount-1)/2 # calculate the coalescent rate
coalescencetime = rexp(1, rate=coalescencerate) # sample from exponential w. that rate
coalescencetime
```

Make sure	you	understand	what	the	R	code	does	before	moving	on.

5. Now	draw	a	node	that	is	the	sampled amount	of	time	further	up	in	the	tree	than	the currently	
highest node	(so	if	the	currently	highest	node	is	drawn	at	height	T	then	draw	the	new	one	at	
height	T plus the	sampled	coalescence	time)	and	draw	a	branch	from	each	of	the	nodes	you	
sampled	in	step	3	to	this	new	node	indicating	that	these	two	nodes	coalesce at	this	time.	

6. Next,	make	an	updated	list	of	the	nodes	that	are	left	by	removing	the two	nodes	that	
coalesced	and	instead	adding	the	newly	drawn	node	that represents	their	common	ancestor.	
You	can	call	the	new	node	the	next	number	not	used	as	a	name	yet	(e.g. if	this	is	the	first	
coalescence event you	can	call	it	6, if	it	is	the	second	coalescence	event you	can	call	it	7	etc.).	
You	can	either	do	this	by	hand	or	in	R.	If	you	want	to	do	it	R	you	can	do	it	as	follows:

```
nodes <- nodes[-tocoalesce] # remove the two nodes that coalesced
nodes <- c(nodes,2*5-length(nodes)-1) # add the new node
nodes # print the new list
```

If	you	used	R	then	make	sure	you	understand	what	the	R	code	does	before	moving	on.

7. If	you	only	have	one	node	left	in	your	list	of	remaining	nodes	you	are	done.	If	not,	go	back	to	
step	3.	

In	the	end you	should	have	a	tree,	which	is	a	simulation	of	a	coalescence	tree. Try	to	do	this	a	
couple times	until	you	feel	like	you	know	how	it	is	done	and	understand	how the coalescence process works	(if after	drawing	a	few	trees still	don’t	understand,	then	feel	free	to	ask for	help!).

## Exercise	B:	Exploring	the	basic	properties	of	a	standard	coalescence tree	

Doing	this	by	hand	is	obviously a	bit	tedious.	So	based	on	the	R	code	snippets	you	already	got, we have built a function	that	allows	you	to	do	this	automatically	(it	even	makes	a	drawing	of	the	tree).	You	
can	use	it from the course server	by	typing the	following	in	R:

```
R
source("/home/fernando/simulatecoalescencetrees.R")
```

Once	you	have	done	this	you can	simulate	and	draw	trees just like	you	just	did	by hand by	typing the code below, which will print out ten trees on the screen:

```
par (mfrow=c(2,5))
for (i in c(1:10)){
         print("New Tree")
         yourtree <-simtree(5) # simulate tree with 5 nodes
         ct<-read.tree(text=yourtree);plot(ct,cex=1.5);add.scale.bar(y=1.2,x=0.2,cex = 2,col = "red",lcol="red",lwd=3)
         print(" ")
}

```

You should see several trees printed out in the screen. If this doesn't happen, try downloading the R script from this github website, and then running it locally in your machine (after you cd to the folder in which you downloaded the script).

```
R
install.packages("ape")
source("simulatecoalescencetrees.R")
```

Note that the code	also	prints	the	simulated	coalescence	times.	

Based	on	the	results you	get answer	the	following	questions:

1) Which	coalescence event takes	the	longest on	average (the	first coalescence event,	the	
second,	…,	or	the	last)?	And	which	event	takes	the	shortest on	average?

2) Is	that	what	you	would	expect? Recall that the	mean	of	an	exponential	distribution	with rate	lambda	is	1/lambda	
and	the	coalescence rate	when	there	are	n nodes	left	is	n(n-1)/2.	So	the	mean	is	2/(n(n-1)),	so
for	instance	for	when	there	are	5	nodes	left	the	mean	coalescent	time	is	2/(5(5-1))=0.1

3) Which	coalescence event	time	seems to	vary	the	most?

4) Is	that	what	you	would	expect? Recall that, if we have a random variable that follows an exponential distribution with rate lambda, then its variance is equal 1/(lambda^2).

5) Finally,	imagine	the	following	case:	a	researcher	has	estimated	the	structure	of	a	tree	for	
mtDNA	from	a	species	sampled	in	a	single	location.	She	obtains	a	tree	looking	as	follows:

![alt text](https://github.com/FerRacimo/CopenhagenTutorial/blob/master/Tree0.png)

Based	on	the	structure	of	the	tree,	i.e.	two	groups	of	related	individuals	separated	by	long	
branches	down	to the	root	of	the	tree,	she	concludes	that	there	must	be population	
subdivision	with	two	clearly	differentiated	groups.	 Based	on	what	you	have	learned	from	
the	simulations,	do	you	agree	with	this	conclusion?
