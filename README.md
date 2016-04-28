# README #

This is an embedding algorithm that takes a graph as input and outputs hyperbolic coordinates.

### Setup ###

* To compile, you need to install the libraries `gflags`, `glog`, and `gsl`. 
* Once that is done, compile with `make`. If it fails, double check that the includes of above libraries can be found in `/usr/local/include`, otherwise change the path in the makefile.
* The input file must be an edge list, i.e. it has two (labels of) nodes per line, delimited by a whitespace. The file may contain comments starting with `#`. Graph files from http://snap.stanford.edu/data/ should typically work.
* Type `./embedder --help` to see an overview over the existing flags.
* To (re-)embed an artificially generated hyperbolic random graph, run e.g. `./embedder --logtostderr --generate="mygraph" --n=5000 --C=-1 --alpha=0.75 --seed=32472351 --embed=myembedding`
* To embed a given graph, run e.g. `./embedder --logtostderr --input="edgelist.txt" --seed=32472351 --embed=myembedding`
* Results are stored in `myembedding-links.txt` and `myembedding-coordinates.txt`, use `--help` for more information

### Who do I talk to? ###

* If you want to use this algorithm, please cite our corresponding paper (currently under review).