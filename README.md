
# Chouca: a probabilistic cellular automaton engine 

Probabilistic cellular automata are a class of models that describe a system in space the 
dynamics of cells on a 2D regular grid. Each cell can be in one of several states. At each 
time step, they can transition to other states with some probability. This probability 
typically depends on the neighbors of the cell and the global state of the landscape. 

You probably already know Conway's game of life -- a probabilistic cellular automaton 
is identical, except that cell transitions do not always occur when a rule is satisfied, 
but with a given probability.

Right now, this project is highly experimental, and interface may change at any moment. Do 
not use this for serious business (but watch when this notice disappears!). 



## What this package implements 

This package is an *engine* for cellular automata. The goal is that you declare your 
transition rules and starting conditions, then the package handles the rest of the 
simulation for you. 

For example, Kubo's forest model (Kubo, 1996), which describes how gaps create and expand in forests, can be implemented using the following few lines of code: 

```r
kubo <- camodel( 
  transition(from = "TREE", 
             to   = "EMPTY", 
             prob = ~ d + delta * q["EMPTY"] ), 
  transition(from = "EMPTY", 
             to   = "TREE", 
             prob = ~ alpha), 
  parms = list(d = 0.125, 
               delta = 0.5, 
               alpha = 0.1), 
  all_states = c("EMPTY", "TREE")
)
```

Running the model for 200 iterations on a 100x100 grid is another few lines of code: 

```r
initmat <- generate_initmat(kubo, c(0.5, 0.5), 100, 100)
run_camodel(kubo, initmat, niter = 200)
```

At the moment `chouca` only runs cellular automata for which the transition probabilities linear combinations of the global covers of each state `p` and the local covers of each state `q`. In other words a transition rule must follow the following pattern: 

```
P(A -> B) = a + b1 p1 + ... + bn pn + c1 q1 + ... + cn qn 
```

where n is the total number of states. 



## Motivation and objectives

Probabilistic cellular automata are widely used in ecology to describe the dynamics of 
organisms in the landscape, and investigate how local interactions between organisms 
may affect the dynamic as a whole. 

`chouca` wants to be user-friendly, but stay reasonably fast. The core of the engine is in C++ (using RcppArmadillo), which helps getting results fast. The goal is to spend
more time on thinking about your results, rather than implementing your model ;)



## Authors and acknowledgements 

`chouca` is mainly developed by Alexandre Génin, but contributions and discussion are 
welcome. 

*A.G. has received funding from the European Union’s Horizon 2020 research and innovation programme under the Marie Sklodowska-Curie grant agreement N°896159.*



## References 

Kubo, Takuya, Yoh Iwasa, and Naoki Furumoto. 1996. “Forest Spatial Dynamics with Gap Expansion: Total Gap Area and Gap Size Distribution.” Journal of Theoretical Biology 180 (3): 229–46.
