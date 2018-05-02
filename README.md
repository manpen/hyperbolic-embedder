# README #

Hyperbolic geometry appears to be intrinsic in many large real networks. We developed and implemented a new maximum likelihood estimation algorithm that embeds scale-free graphs in the hyperbolic space. Our algorithm achieves quasilinear runtime, which makes it the first algorithm that can embed networks with hundreds of thousands of nodes in less than one hour. In typical metrics like Log-likelihood and greedy routing, our algorithm discovers embeddings that are very close to the ground truth.

### Setup ###

* To compile, you need to install the libraries `gflags` (>=2.1.2), `glog` (>=0.3.4), and `gsl` (>=1.16). 
* Once that is done, compile with `make`. If it fails, double check that the includes of above libraries can be found in `/usr/local/include`, otherwise change the path in the makefile. The makefile also looks for the libraries in `/usr/local/lib`.
* The input file must be an edge list, i.e. it has two (labels of) nodes per line, delimited by a whitespace. The file may contain comments starting with `#`. Graph files from http://snap.stanford.edu/data/ should typically work.
* Type `./embedder --help` to see an overview over the existing flags.
* To (re-)embed an artificially generated hyperbolic random graph, run e.g. `./embedder --logtostderr --generate="mygraph" --n=5000 --C=-1 --alpha=0.75 --seed=32472351 --embed=myembedding`. This will create and store a hyperbolic random graph, so it needs no further inputs.
* To embed a given graph, run e.g. `./embedder --logtostderr --input="edgelist.txt" --seed=32472351 --embed=myembedding`
* Results are stored in `myembedding-links.txt` and `myembedding-coordinates.txt`, use `--help` for more information

### Misc. ###

* If you want to use this algorithm, please cite our [corresponding paper](http://dx.doi.org/10.1109/TNET.2018.2810186).	 
* Please contact anton.krohmer@hpi.de for further questions.