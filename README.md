# lattice_CS
group project for UCLA Math 285J spring 2018 on compressive sensing on lattices

cosamp.m and cosampround.m are matlab function files which contain the cosamp algorithm and the cosamp-but-round-after-every-iteration algorithm. workspace.m is a matlab script file in which I've been playing around with the cosamp algorithm

(added by Adam, 4/24): cosamproundend.m is CoSaMP where it only rounds at the end.  testing.m is a script which ranges the error norm from 0 to 5 in steps of 0.1, and at each error norm it runs the three CoSaMPs 50 times and takes the average of norm(reconstruction - signal), and plots those as a function of the error norm.  The picture looks similar to the professor's but for some reason intermediate rounding and end rounding are giving the exact same recontsruction every single time, even though I used all the same parameters as the professor.  Not sure why.
