#compile and run single target problem
make SG
./SingleTarget --filepath=Datasets/graph_liquor_updated --method=IP



#compile and run multi targets problem
make MT
./MultiTarget --filepath=Datasets/graph_bankClean_updated --method=IP


