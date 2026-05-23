# Secure computation of social prestige and influence from multiple networks

This repository contains the research code for the secure multiparty centrality protocols presented in:

Gilad Asharov, Francesco Bonchi, David Garcia-Soriano and Tamir Tassa:
"Secure Computation of Social Prestige and Influence from Multiple Networks",
in Proceedings of the 26th International World Wide Web Conference (WWW), pages 957-966, 2017.
       [PDF](https://www.cs.upc.edu/~dgarcia/pubs/maxmin-fair-ranking-individual-fairness-under-group-fairness-constraints.pdf)

## License

Everything in this repository is distributed under the Apache License, version
2.0. See the [LICENSE](LICENSE) file and the [NOTICE](NOTICE) file.  

## Building

The code is written in C++11 and requires a C++ compiler, Boost headers, and
Victor Shoup's Number Theory Library (NTL).

Build with:

    make

This produces the `secure_inf` executable.

## Usage

Run:

    ./secure_inf [-u] [-a] [-c] D GRAPH_1 GRAPH_2 ...

Arguments and options:

    D          maximum path length.
    GRAPH_i    edge list for one host/layer.
    -u         treat input graphs as undirected by duplicating every edge.
    -a         also run the slower secure protocols for pi_1 and pi_2.
    -c         compute and print the three reference centrality vectors.

Each graph file is a text edge list. Lines beginning with `#` are ignored. Each
non-comment line has either `source target` or `source target weight`; omitted
weights default to 1.

By default the program computes exact reference scores for pi_2 and pi_3, runs
the efficient secure protocol for pi_3, and prints a comparison. With `-c`, the
centrality vectors are written in Python syntax to standard error for downstream
analysis scripts.

## Contact

For questions or bug reports, please contact David Garcia-Soriano
<david.garcia.soriano@upc.edu>.
