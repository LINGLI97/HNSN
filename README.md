# Source code for the paper: Heavy Nodes in A Small Neighborhood: Exact and Peeling Algorithms with Applications''

Our approach was compiled with g++ version 11.4.0 using the following flags: -std=c++11 and -O3. 
We ran experiments on Ubuntu 22.04.3 LTS.

## Install
1. Before compiling, please install the [Gurobi](https://www.gurobi.com/downloads/) solver.
2. After installing Gurobi, please change the path `/home/ling/opt/gurobi912/linux64` in the Makefile into `\path\to\your\installation\dir`.
3. Install [boost](https://www.boost.org/users/download/) library.
``` bash
sudo apt install libboost-all-dev
```

## Compile and run

```bash
./compile.sh
```

## Reproducibility & Example to Run the Codes

We utilize the parameters --filepath and --method to specify the dataset and algorithm, respectively. Use Macros VERBOSE and PRINT to check the log information.
```bash
./SingleTarget --filepath=Datasets/graph_liquor_updated --method=IP
```
**--------------------------Ouput of Iterative peeling on Liquor dataset-------------------------** \
Now run IP on graph_liquor_updated\
numSource = 4026 , numMiddle  = 52131\
numEdges = 410619\
Score is 66195.3 in iteration 1\
Score is 66195.3 in iteration 2\
Score is 73073 in iteration 3\
Score is 73073 in iteration 4\
Score is 73073 in iteration 5\
Score is 73073 in iteration 6\
Score is 73073 in iteration 7\
Score is 73073 in iteration 8\
Score is 73073 in iteration 9\
Score is 73073 in iteration 10\
Score is 73073 in iteration 11\
Score is 73073 in iteration 12\
Score is 73073 in iteration 13\
Rutime for IP = 11327 ms\
Score for IP = 73073\
Max score was achieved in iteration 3

```bash
./MultiTarget --filepath=Datasets/graph_bankClean_updated --method=IP
```
**--------------------------Ouput of Iterative peeling on Czech Financial dataset-------------------------** \
Score is 0.00794237 in iteration 1\
Score is 0.00794237 in iteration 2\
Score is 0.00794237 in iteration 3\
Score is 0.00794237 in iteration 4\
Score is 0.00794237 in iteration 5\
Score is 0.00794237 in iteration 6\
Score is 0.00794237 in iteration 7\
Score is 0.00794237 in iteration 8\
Score is 0.00794237 in iteration 9\
Score is 0.00794237 in iteration 10\
Score is 0.00794237 in iteration 11\
Score is 0.00794237 in iteration 12\
Score is 0.00794237 in iteration 13\
Score is 0.00794237 in iteration 14\
Score is 0.00794237 in iteration 15\
Score is 0.00794237 in iteration 16\
Score is 0.00794237 in iteration 17\
Score is 0.00794237 in iteration 18\
Score is 0.00794237 in iteration 19\
Score is 0.00901631 in iteration 20\
Score is 0.0220886 in iteration 21\
Score is 0.0220886 in iteration 22\
Score is 0.0220886 in iteration 23\
Score is 0.0220886 in iteration 24\
Score is 0.0220886 in iteration 25\
Score is 0.0220886 in iteration 26\
Score is 0.0220886 in iteration 27\
Score is 0.0220886 in iteration 28\
Score is 0.0220886 in iteration 29\
Score is 0.0220886 in iteration 30\
Score is 0.0220886 in iteration 31\
Score is 0.0220886 in iteration 32\
Score is 0.0220886 in iteration 33\
Score is 0.0220886 in iteration 34\
Score is 0.0220886 in iteration 35\
Score is 0.0220886 in iteration 36\
Score is 0.0220886 in iteration 37\
Score is 0.0220886 in iteration 38\
Score is 0.0220886 in iteration 39\
Score is 0.0220886 in iteration 40\
Score is 0.0220886 in iteration 41\
Rutime for IP = 7992 ms\
Score for IP = 0.0220886\
Max score was achieved in iteration 21


## Information about FlowScope & AA_Smurf

Please refer to [here](https://github.com/csqjxiao/FlowScope) and [here](https://github.com/mengchillee/AutoAudit/tree/master) for the implementation of FL and AA_Smurf, respectively.


## Datasets
We provide 12 bipartite graph datasets and 3 tripartite datasets, see the details in the table .

| Datasets                      | filename                  |
|-------------------------------|---------------------------|
| Foodmart (FM)                 | graph_foodmart_updated    |
| E-commerce (EC)               | graph_ecommerce_updated   |
| Liquor (LI)                   | graph_liquor_updated      |
| Fruithut (FR)                 | graph_fruithut_updated    |
| YooChoose (YC)                | graph_yoochoose_updated   |
| Kosarak (KS)                   | graph_kosarak_updated     |
| Connectious (CN)               | graph_connectious_updated |
| Digg (DI)                      | graph_digg_updated        |
| NotreDame (ND)                 | graph_notredame_updated   |
| IMDB (IM)                      | graph_imdb_updated        |
|NBA Shot (NBA) |           graph_NBA_updated               |
|ACM Citation (ACM)|          graph_ACM_updated|
|Czech Financial Dataset (CFD) | graph_bankClean_updated   |
|$AML_H$| AML_H/wo_labels_updated.txt|
|$AML_L$|AML_L/wo_labels_updated.txt|

The node indices are consecutive, starting from 0. We first index V, then add indices for U, and W (if applicable).

The format of single target datasets: 
- `graph_***_updated`: 
	- The first line is the number of source nodes |U|
	- The second line is the number of middle nodes |V|
	- Then, each of the next |V| lines contains the weight of each v node. 
	- The last |E| lines are the edges in the graph.

The format of multi target dataset `CFD`, `AML_H` and `AML_L`: 
	- The first line is the number of source nodes |U|
	- The second line is the number of middle nodes |V|
	- The third line is the number of target nodes |W|
	- The fourth line is the number of edges between U and V.
	- Then, each of the next |V| lines contains the weight of each v  node. 
	- The last |E| lines are the edges in the graph.

Since `AML_H` and `AML_L` have labels, we prepared two files (`wo_labels_updated.txt` and `w_labels_updated.txt`) under each folder. `w_labels_updated.txt` includes labels, while `wo_labels_updated.txt` does not. Please use `wo_labels_updated.txt` as the  input of our methods because we are unsupervised methods. 



## License

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see [this](http://www.gnu.org/licenses/).

## Citing & Authors

...