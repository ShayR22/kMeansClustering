run()

systemInit():
	- master will read file
	- master will broadcast to all slaves:
		* number of points
		* points
		* number of clusters
		* delta t of between iterations
		* LIMIT (the limit to the kMeans iterations number)
	- every machine allocate memory clusters, clusterHelpers and points

distributePoints():
	- master send to each slave 2 numbers:
		* the number from which to start calculating points within the points array
		* the number of points to work on from the starting point
	- slave recv those numbers

loop until: (time > T or q < QM) (meaining no more iterations required or a quality meassre good enough has been achevied)
	kMeans() // stabelize the system or at least try to within the LIMIT of iterations that has been given from the file
	for full description look down in the this README file.

	distributeHelpers():
	- master assign clusterHelpers to each machine
	- slave recv data from master
	
	calcDiameter():
	-  each machine calculating the diameter of the helpers recv from distributeHelpers() method, the method
	go over each helper and find the maximum distance between any 2 points that belongs to the helper, O(n^2)

	foldHelpers():
	- slaves send diameter results to master
	- master put results in the right place

	calcQ():
	- master calculate the q.



kMeans(): (full explanation)
	- master broadcast all the clusters
	- reset clusterHelpers, go over its helper and make its diameter, number of points and sum of location all to 0.
	- each machine assign for each point in its range, its nearest cluster.
	- each machine go over its helpers and for each one go over all the points and assign each point that belong to it, in an 
	index (IDs) array. 
		NOTE: this calculation might be O(k*n) but with this split finding diameters later would be less computioan heavy cause we would be able
		splitted between both CPU and GPU. (which is O(n^2)) 
		and also avoid a critical section (this happen due to the fact when going over points if 2 points belong
	    to the same cluster who will get in which spot in the array).

	- each machine go over its helpers and sum the all of its points (from its IDs array) into the helper
	- all slaves send thier helpers to master
	- master recv every k helpers into a seperate space inside its helpers (master helpers space = numOfMachines*k)
	- master fold the results of all its k sets of helpers into one1 set of k helpers.
	- master find the centers of helpers by: dividing the sumLocation with its number of points.
	- master comprare each helper center its correlative cluster center, if they differ
	the method will update the cluster center.
	- if any changes in any cluster occurred and the number of iteration is less then LIMIT repeat the entire process.   





/// ---------------------------------- PROBABLY END HERE FOR NOW -------------------------------------------























OVERVIEW:
	having T/dt intervals in each interval all the clusters will be calculated using the LOGIC below.
	after all the clusters have been calculated, the QM (quality measure) will be calculated, if the QM that
	was calculated is less the const of QM given the program is done, if not go to next iteration 


LOGIC:
	- Each iteration the cluster need to be assigned all the nearest points
	to it, this can be paralleled if we let the point to decide where to go (can this be paralleled?).
	HOW-TO
		- In order to parallel the cluster, we can have each cluster have a static array of pointers to points
		in the size of all the points.
		this way memory consumption is low (array of pointers) and parallel can be achivied by having each point an
		id 1..N, each point will assign itself to its id location in the cluster's array, NULL would be indication 
		of an empty slot. 

	- After each iteration if there were any changes, meaning one of the centers of the clusters has been changed,
	(or points reassigned (maybe its the same)). 

	- In case there were changes in one of the centers => clean each cluster pointer array and repeat clusters
	

PARALLEL:
	MPI:
	- each process needs to know all the clusters so he can identify points correctly, therefore points should be
	parallelized:

		--- kMeansIteration Start ---
		* resetAllClusters()
			+ basically every process can do this, instead of just master do and re-send, this way no need to send the cluster
			  all over.

		* addPointsToNearestCluster()
			+ each process will calculate numPoints/numProcs, and assign it to the right cluster, therfore
			  master will go over the points and re-assign points based on id to the cluster (cant copy the pointer
			  of recv cluster as pointer dont have to be same in the other machine).
		
		* calculateNewCenter()
			+ after completing each point assignment master can scatter the work of each cluster to another machine
			  by sending the cluster and let it send  

		--- kMeansIteration End -----

		* increment point
			+ each process knows all its points and can just increment those points
		
		* QM (quality measure) can it be done using different machines?


	



	----------------------------------------------------------------------------------------------

	TEST

	cluster is just a
	 - id (maybe not needed) 
	 - center (location vector)
	 - max distance
	 point is: 

	 - id (maybe not needed)
	 - loacation vector
	 - speed vector
	 - cluster number belonging
	 - distance from cluster

	 maybe: clusterHelper
	 - sumOfPoints (location vector)
	 - numOfPoints
	 
	  
	 stages for kMeansAlgorithm
	  1) master broadcast its clusters (DONE)
	  2) master spread points across processes (DONE)
	  3) each process (including master) go over its own points and calculate for each point to which cluster
	    its belongs too. and whats its distance from it (DONE)
	  4) each process will determine its own set of clusterHelper 
	  5) send results to master where master will calculate the new means (clusters check QM(mabye with mpi))
	  
	 master broadcast clusters.
	 master scatter all points across machines.
	 each machine go overs it point and compare distances with every cluster and determine the cluster which
	 it is most close too and its distance from it.
	 clusterHelpers could have an array of pointers to points and an int indicating the number of points in the
	 clusterHelper, this could be easily done with openMp but not so much with CUDA if the distances

	 steps:
	 CUDA + openMP determine each point distance and cluster belonging
	 openMP go over and assign each point to a cluster
	
	