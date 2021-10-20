# Mesh-generation <br>
## A simple algorithm for maximal Poisson-disk sampling 
We provide a simple algorithm and data structures implementation  for 2D unbiased maximal Poisson-disk sampling.This algorithm use an order of magnitude less memory and time 
than the alternatives. 
Our results become more favorable as the dimension increases. This allows us to produce bigger samplings. Domains may be non-convex with holes.
The generated point cloud is maximal up to round-off error. The serial algorithm is provably bias-free. For an output sampling of size n in fixed dimension d,
we use a linear memory budget and empirical Θ(n) runtime.


The algorithm proceeds through a finite sequence of uniform grids. The grids guide the dart throwing and track
the remaining disk-free area. The top-level grid provides an efficient way to test if a candidate dart is disk-free.
Our uniform grids are like quadtrees, except we delay splits and refine all leaves at once. Since the quadtree is
flat it can be represented using very little memory: we just need the indices of the active leaves and a global level.
Also it is very simple to sample from leaves with uniform probability.

![image](https://user-images.githubusercontent.com/67281513/127719944-d3bf3a3f-8258-4ee4-b14e-de5f29e9b5f7.png)

                                                  fig.1
In case of figure 1 domain width is 1cm and length is 1cm  and the disk radius is .01cm .

## Voronoi diagram
providing a voronoi diagram based on neighbor cells " we get that cells from previous algorithm".
That algorithm uses hyperplan idea for creating the voronoi cell that has a corner in equal disance from all neighbours cells.
![image](https://user-images.githubusercontent.com/67281513/138172654-ef4ecb27-df28-4a0e-95ed-0de293df90ee.png)              
                                               
                                               fig.2
For high quality voronoi cells we use an uniform distribution of seeds "between each seed and another same distance" 
![image](https://user-images.githubusercontent.com/67281513/138172952-c4605243-0c9a-492c-94e3-b426e9cae3f7.png)
                                                
                                                fig.3
in figure 3 , the relation between number of corners and its angle .

## Referances

EBEIDA M., MiITCHELL S., PATNEY A., DAVIDSON A., OWENS J. : A simple algorithm for maximal Poisson-disk sampling in high dimensions. Computer Graphics Forum 31, 2pt4 (2012), 785-794. 
https://onlinelibrary.wiley.com/doi/full/10.1111/j.1467-8659.2012.03059.x
