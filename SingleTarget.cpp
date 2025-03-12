//
// Created by ling on 02/01/24.
//

#include <iostream>
#include <map>
#include <unordered_map>
#include <set>
#include <float.h>
//#include "matplotlibcpp.h"
#include <dirent.h>

#include <unordered_set>
#include <fstream>

#include <chrono>
#include <vector>
#include<cmath>
#include <algorithm>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include<boost/heap/fibonacci_heap.hpp>
//#include<boost/heap/fibonacci_heap_bp.hpp>

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <queue>
#include "gurobi_c++.h"

using namespace std;
using namespace chrono;
using namespace boost::heap;
#include <filesystem>
#include "cmdline.h"
#define VERBOSE
#define PRINT
//namespace fs = std::filesystem;
#define log(x) std::cout<<x<<endl

struct Tolerance{
    int max = 10;
    int cnt = 1;
    int iteration =1;
    int total=1;
};

void Stringsplit(string str, const char split,vector<string>& res)
{
    istringstream iss(str);
    string token;
    while (getline(iss, token, split))
    {
        res.push_back(token);
    }
}


tuple<vector <unordered_map<int, unordered_set<int>>>,vector <int>,vector <int>,vector <int>, vector <vector<double>>, unordered_map<int, int>> readFile(string filename, double &pre_decomposition_time,int merge_size ,int merge_bound){

    string line;
    ifstream myfile (filename);
//    vector <unordered_map <int, int>> mapping_relation;
    int numSource =0, numMiddle = 0, numEdges = 0;
//    vector <vector<double>> Multi_score;
    unordered_map<int, vector<int>> trans_vec;
    unordered_map< string, unordered_set<int>> trans_string;


    if (myfile.is_open())
    {
        getline (myfile,line);
        numSource =  stoi(line);
        // cout<<"numSource = "<<numSource<<endl;
        getline (myfile,line);
        numMiddle =  stoi(line);

        vector<double> score(numMiddle);


        int counter = 0;

        while(counter<numMiddle){
            getline (myfile,line) ;
            score[counter] = stof(line);
            counter++;
        }

        using namespace boost;
        {

            using std::unordered_set;
            using std::unordered_map;

            typedef adjacency_list<vecS, vecS, undirectedS> Graph;
            Graph G;
            int inNode, outNode;
            double small_time_sum=0.0;

            while (myfile >> inNode >> outNode) {
                auto small_cnt_start = std::chrono::high_resolution_clock::now();

                numEdges++;
                trans_vec[outNode].push_back(inNode);
                add_edge(outNode, inNode, G);
                auto small_cnt_end = std::chrono::high_resolution_clock::now();
                auto small_time=  std::chrono::duration_cast < std::chrono::milliseconds > (small_cnt_end - small_cnt_start).count();
                small_time_sum += small_time;
            }


            auto pre_decomposition_time_start = std::chrono::high_resolution_clock::now();

            for (auto it = trans_vec.begin();it != trans_vec.end() ; it++){
                // radix sort the innodes first then make them become strings
//                (*it).second = radixsort((*it).second, (*it).second.size());

                sort((*it).second.begin(),(*it).second.end());


                stringstream ss;

                copy((*it).second.begin(),(*it).second.end(),ostream_iterator<int>(ss," "));
                string s=ss.str();
                //get rid of the useless " "
                s = s.substr(0, s.length()-1);

                // make the <innodes string, outnode set> map
                auto iter = trans_string.find(s);
                if (iter != trans_string.end()) {
                    trans_string[s].insert((*it).first);
                } else {
                    unordered_set<int> temp;
                    temp.insert((*it).first);
                    trans_string[s] = temp;
                }

            }


            vector<int> component(num_vertices(G));
            int num = connected_components(G, &component[0]);
#ifdef VERBOSE

            cout<<"Before merging, the number of components: "<<num<<endl;
#endif
            vector <int> Multi_numMiddle_old(num);

            for (int i = 0; i<numMiddle; i++){
                Multi_numMiddle_old[component[i]]++;
            }
            vector<int> component_mapping(num);

            bool Merge_flag = false;
            int component_id;
            int new_num =0;
            int merge_component_size= 0;

            for (int i = 0; i<num; i++){
                if (1<Multi_numMiddle_old[i] && Multi_numMiddle_old[i]<merge_size){
                    if (Merge_flag && merge_component_size< merge_bound){
                        component_mapping[i] = component_id;
                        merge_component_size += Multi_numMiddle_old[i];

                    } else{
                        component_mapping[i] = new_num;
                        component_id = new_num;
                        Merge_flag = true;
                        new_num++;
                        merge_component_size = Multi_numMiddle_old[i];

                    }
                } else{

                    component_mapping[i] = new_num;
                    new_num++;

                }
            }
#ifdef VERBOSE
            cout<< "After merging, the number of component: "<<new_num<<endl;

#endif

            vector <unordered_map<int, unordered_set<int>>> Multi_trans(new_num);

            vector <int> Multi_numSource(new_num);
            vector <int> Multi_numMiddle(new_num);
            vector <int> Multi_numEdges(new_num);
            vector <vector<double>> Multi_score(new_num);

            unordered_map<int, int> mapping_relation;

            vector<int> current_counts_per_component(new_num);

            for (auto iter = trans_string.begin(); iter!= trans_string.end(); iter++){
                int id;
                for (auto i :iter->second){
                    id = component_mapping[component[i]];
                    int cnt= current_counts_per_component[id];
                    //map the outnodes in the same set to one node
                    mapping_relation[i]=cnt;
//                    cout << i <<":"<<cnt<<endl;
                }
                current_counts_per_component[id]++;
                Multi_numMiddle[id]++;
                vector<string> strList;
                Stringsplit(iter->first, ' ', strList);
                Multi_numEdges[id]+= strList.size();
            }



            for (int i = numMiddle; i != numMiddle + numSource;++i)
            {
                //build the mapping relation for source nodes
                int id = component_mapping[component[i]];
                int cnt=current_counts_per_component[id];
                mapping_relation[i]=cnt;

//                cout << i <<":"<<cnt<<endl;
                current_counts_per_component[id]++;
                Multi_numSource[id]++;
            }

            vector<bool> component_flags(num, false);
            for (int i = 0; i < numMiddle; ++i)
            {
                //initialization Multi_score
                int id = component_mapping[component[i]];

                if (!component_flags[id]){
                    vector<double> multi_score_for_component(Multi_numMiddle[id]) ;
                    Multi_score[id]=multi_score_for_component;
                    component_flags[id] = true;
                }
                //accumulating the scores for all nodes after combination

                Multi_score[id][mapping_relation[i]] += score[i];


                unordered_set<int> tmp;

                for(auto &it2 : trans_vec[i])  //create a set with the new ids
                {
                    tmp.insert(mapping_relation[it2]);
                }
                Multi_trans[id][mapping_relation[i]]=tmp; //add it into the correct unordered_map
//                }
            }



            auto pre_decomposition_time_end = std::chrono::high_resolution_clock::now();

            myfile.close();

            pre_decomposition_time= std::chrono::duration_cast < std::chrono::milliseconds > (pre_decomposition_time_end - pre_decomposition_time_start).count() + small_time_sum;
//            cout<<"Runtime for Decomposition= "<<decomposition_time<<" ms"<<endl;

            return std::make_tuple(Multi_trans, Multi_numSource, Multi_numMiddle, Multi_numEdges, Multi_score, mapping_relation);



        }



    }//end if open file
    else cout << "Unable to open file";

}

vector<double> readFile(const string &filename, unordered_map <int, unordered_set<int>> &trans,
                        int &numSource, int &numMiddle, int &numEdges){
    string line;
    ifstream myfile (filename);
    if (myfile.is_open()){
        getline (myfile,line);
        numSource =  stoi(line);
        getline (myfile,line);
        numMiddle =  stoi(line);
#ifdef VERBOSE

        cout<<"numSource = "<<numSource<<", numMiddle  = "<<numMiddle<<endl;
#endif
        vector<double> score(numMiddle);
        score.reserve(numMiddle);

        int counter = 0;

        while(counter<numMiddle){
            getline (myfile,line) ;
            score[counter] = stof(line);
            counter++;
        }

        int inNode, outNode;
        while(myfile >>inNode>>outNode){
            numEdges++;
            auto it =  trans.find(outNode);
            if(it!= trans.end()){
                trans[outNode].insert(inNode);
            }
            else{
                unordered_set<int> temp;
                temp.insert(inNode);
                trans.insert({outNode, temp});

            }

        }
#ifdef VERBOSE
        cout<<"numEdges = "<<numEdges<<endl;
#endif
        myfile.close();
        return score;
    }

    else cout << "Unable to open file";

    vector<double> score;

    return score;

}

vector<double> readFile(const string &filename, unordered_map <int, unordered_set<int>> &midSource,
                        int &numSource, int &numMiddle, int &numEdges,unordered_map <int, unordered_set<int>> &sorceMid){
    string line;

//    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;
//    Graph G;
    ifstream myfile (filename);
    if (myfile.is_open()){
        getline (myfile,line);
        numSource =  stoi(line);
        // cout<<"numSource = "<<numSource<<endl;
        getline (myfile,line);
        numMiddle =  stoi(line);
        cout<<"numSource = "<<numSource<<" , numMiddle  = "<<numMiddle<<endl;

        vector<double> score(numMiddle);

        int counter = 0;

        while(counter<numMiddle){
            getline (myfile,line) ;
            score[counter] = stof(line);
            counter++;
        }




        int inNode, outNode;
        while(myfile >>inNode>>outNode){
            numEdges++;
            // cout << "in and out : "<<inNode << " "<<outNode<<endl;

            auto it =  midSource.find(outNode);
            if(it!= midSource.end()){
                midSource[outNode].insert(inNode);
//                add_edge(outNode, inNode, G);
            }
            else{
                unordered_set<int> temp;
                temp.insert(inNode);
                midSource.insert({outNode, temp});
//                add_edge(outNode, inNode, G);
            }

            auto iter_mid =  sorceMid.find(inNode);
            if(iter_mid!= sorceMid.end()){
                sorceMid[inNode].insert(outNode);
                //cout<<"this case insert "<<outNode<<" to "<<inNode<<endl;
            }
            else{
                unordered_set<int> temp;
                temp.insert(outNode);
                sorceMid.insert({inNode, temp});
                //cout<<"insert "<<outNode<<" to "<<inNode<<endl;

            }


            //cout<<"inNode = "<<inNode <<", outNode = "<<outNode<<endl;
        }





//
//        std::vector<int> component_(num_vertices(G));
//        int num = connected_components(G, &component_[0]);
//        component = component_;

//        cout<<"num components= "<<num<<endl;

        cout<<"numEdges = "<<numEdges<<endl;



        myfile.close();
        return score;
    }

    else cout << "Unable to open file";

}




vector<double> readFile(const string &filename, unordered_map <int, unordered_set<int>> &midSource,
                        int &numSource, int &numMiddle, int &numEdges,unordered_map <int, unordered_set<int>> &sorceMid,double &MaxScoreSeparateU,unordered_set<int> &MaxSet){
    string line;

//    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;
//    Graph G;
    ifstream myfile (filename);
    if (myfile.is_open()){
        getline (myfile,line);
        numSource =  stoi(line);
        getline (myfile,line);
        numMiddle =  stoi(line);
#ifdef VERBOSE
        cout<<"numSource = "<<numSource<<" , numMiddle  = "<<numMiddle<<endl;
#endif

        vector<double> score(numMiddle);

        int counter = 0;

        while(counter<numMiddle){
            getline (myfile,line) ;
            score[counter] = stof(line);
            counter++;
        }


//        vector<vector<int>> MidVec;
        vector<vector<int>> SourceVec(numSource);

        vector<int> MidDegreeOne(numMiddle, 1); //all initialized as 1
//        vector<int> SourceVec;


        int inNode, outNode;
        while(myfile >>inNode>>outNode){
            numEdges++;
            // cout << "in and out : "<<inNode << " "<<outNode<<endl;

            auto it =  midSource.find(outNode);
            if(it!= midSource.end()){
                midSource[outNode].insert(inNode);
//                MidVec[outNode].push_back(inNode);
//                add_edge(outNode, inNode, G);
                MidDegreeOne[outNode]=0;  // set value as o if degree >1
            }
            else{
                unordered_set<int> temp;
                temp.insert(inNode);
                midSource.insert({outNode, temp});
//                MidVec[outNode].push_back(inNode);

//                add_edge(outNode, inNode, G);
            }

            auto iter_mid =  sorceMid.find(inNode);
            if(iter_mid!= sorceMid.end()){
                sorceMid[inNode].insert(outNode);
                SourceVec[inNode-numMiddle].push_back(outNode);

                //cout<<"this case insert "<<outNode<<" to "<<inNode<<endl;
            }
            else{
                unordered_set<int> temp;
                temp.insert(outNode);
                sorceMid.insert({inNode, temp});
                SourceVec[inNode-numMiddle].push_back(outNode);

                //cout<<"insert "<<outNode<<" to "<<inNode<<endl;

            }


            //cout<<"inNode = "<<inNode <<", outNode = "<<outNode<<endl;
        }




        for (auto &i: SourceVec){
            unordered_set<int> MaxV;
            double scoretmp = 0.0;
            for (auto &j : i){
                if(MidDegreeOne[j] == 1){
                    MaxV.insert(j);
                    scoretmp += score[j];
                }
            }
            if(scoretmp> MaxScoreSeparateU){
                MaxScoreSeparateU = scoretmp;
                MaxSet = MaxV;

            }
        }
//
//        std::vector<int> component_(num_vertices(G));
//        int num = connected_components(G, &component_[0]);
//        component = component_;

//        cout<<"num components= "<<num<<endl;

        cout<<"numEdges = "<<numEdges<<endl;



        myfile.close();
        return score;
    }

    else cout << "Unable to open file";

}



void ILP(const int &numMiddle, vector<double> &score, const int &numSource,
         unordered_map <int, unordered_set<int>> &trans, vector<int> &finalSetILP,
         double &objRE, const string &filename, const int &numEdges, int &sourceSelect, vector<int> &sourceSelectVec){
    try {

        // Create an environment
        GRBEnv env = GRBEnv(true);
//        env.set("LogFile", "ilp.log");
        env.set("LogFile", "ilp.log");
//        env.set(GRB_IntParam_OutputFlag,0);

        env.start();

        // Create an empty model
        GRBModel model = GRBModel(env);

        // Create variables
        double *w{ new double[numMiddle]{} };

        double *z_ub{ new double[numMiddle]{} };
        double *v_ub{ new double[numSource]{} };

        char* z_type = new char[numMiddle];
        char* v_type = new char[numSource];

        for(int i = 0; i< numMiddle; i++){
            z_type[i] = GRB_CONTINUOUS;
            z_ub[i] = 1.0;
            w[i]=score[i];
        }

        for(int i = 0; i< numSource; i++){
            v_type[i] = GRB_CONTINUOUS;
            v_ub[i] = 1.0;
        }

        GRBVar *v = model.addVars(NULL,v_ub, NULL, v_type, NULL, numSource);
        GRBVar *z = model.addVars(NULL,z_ub, NULL, z_type, NULL, numMiddle);

        GRBLinExpr  myexpr = 0;
        GRBLinExpr  sumV = 0;
        // GRBLinExpr  ysum = 0;

        for(int i=0;i<numMiddle;++i){
            myexpr += z[i]*w[i];

        }
        for(int i=0;i<numSource;++i){
            sumV += v[i];
            // ysum+=y[i];
        }

        // Add constraint: x + 2 y + 3 z <= 4
        model.addConstr(sumV-1 ==0);

        for(auto &it: trans){
            for(auto &socrceNode : it.second){
                model.addConstr(z[it.first]-v[socrceNode-numMiddle]<= 0);
//                model.addConstr(z[it.first]-v[socrceNode]<= 0);

            }
        }


        // Set objective: maximize x + y + 2 z
        model.setObjective(myexpr, GRB_MAXIMIZE);
        // Optimize model
        model.optimize();


        for(int i = 0; i< numMiddle ;i++){
            if(z[i].get(GRB_DoubleAttr_X)>0){


                finalSetILP.push_back(i);
            }

        }


        objRE = model.get(GRB_DoubleAttr_ObjVal);
//        cout << "MaxScore for LP: " << objRE << endl;



        for(int i = 0; i< numSource ;i++){
            if(v[i].get(GRB_DoubleAttr_X)!=0){
                sourceSelectVec.push_back(i);
                sourceSelect++;
            }

        }




//        cout<<"Number of item selected = "<<sourceSelect<<endl;


        delete []z;
        delete []v;
        delete []w;
        delete []z_ub;
        delete []v_ub;

        delete []v_type;
        delete []z_type;


    } catch(GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch(...) {
        cout << "Exception during optimization" << endl;
    }
}



void update_score_mp(multimap<double, int> &sourceScore_id, unordered_map<int, double> &id_sourceScore,const unordered_map<int, double> &node_affected,const unordered_set<int> &node_erased) {
//    multimap<double, int>::iterator itlow, itup;



    for (auto &i: node_erased) {
        double old_score = id_sourceScore[i];
        pair<multimap<double, int>::iterator, multimap<double, int>::iterator> all_score = sourceScore_id.equal_range(old_score);

        for (auto &iter = all_score.first; iter != all_score.second; ++iter) {
            if (iter->second == i) {
//                    double new_score = get_mid_score(i, score, sorceMid);

                sourceScore_id.erase(iter);
                break;

            }
        }

    }



    for (auto &i: node_affected) {
        int node_id = i.first;
        double old_score = id_sourceScore[node_id];
        pair<multimap<double, int>::iterator, multimap<double, int>::iterator> all_score = sourceScore_id.equal_range(old_score);

        if (all_score.first != all_score.second) { // the elements exsits
            for (auto& iter = all_score.first; iter != all_score.second; ++iter) {
                if (iter->second == node_id) {
                    double new_score = old_score - i.second;
                    sourceScore_id.erase(iter);
                    sourceScore_id.insert({new_score, node_id});
                    id_sourceScore[node_id] = new_score;
                    break;

                }
            }

        }

    }

}




void update_mp(unordered_map <int, unordered_set<int>> &midSource,
               unordered_map <int, unordered_set<int>> &sorceMid, vector<double> &score, const int& sourceErase, double &sum_mid, int &denominator, unordered_map<int, double> &node_affected, unordered_set <int> &node_erased){
    denominator -= 1; // for the sourceErase
    unordered_set <int> midNodes = sorceMid[sourceErase]; // get the midnodes set of the source that need to be erased
    for(auto &mid : midNodes){
        int numlink = midSource[mid].size();
        if( numlink!= 1){   // erase it directly if the midnode is only connected to this sourceErase
            for (auto &i :midSource[mid]){ // get every sourcenode that is linked to every midnode
//                auto posVal= find_if(sourceScore.begin(),sourceScore.end(), cmp); // find the source node (weight, id)
                if (i==sourceErase){
                    continue;
                }

                if (sorceMid[i].size() != 1){
                    auto iter_affected = node_affected.find(i);
                    if(iter_affected!= node_affected.end()){
                        iter_affected->second += score[mid];
                    }else{
                        node_affected[i] = score[mid];
                    }
                    sorceMid[i].erase(mid);
                }else{

//                    node_erased.insert({score[*sorceMid[i].begin()],i});// size =1 means: after remove the middle node, the neighbor source node is isoloated
                    auto iter = node_affected.find(i);
                    if (iter!= node_affected.end()){
                        node_affected.erase(iter);
                    }
                    node_erased.insert(i);
                    sorceMid.erase(i);//remove it from sorcemid map
                    denominator -= 1;

                }

            }
            midSource.erase(mid);

        } else{

            midSource.erase(mid);
        }
        sum_mid -= score[mid];
    }
    sorceMid.erase(sourceErase);

}


double GAR(const int& numMiddle, vector<double> &score, const int& numSource,
                           unordered_map <int, unordered_set<int>> &midSource,
                           unordered_map <int, unordered_set<int>> &sorceMid, unordered_set<int> &Pickedset){

//    unordered_set<int> nodeSet; //nodes are selected, at the beginning it select all nodes
    //unordered_set<int> nodeNoneZero; //the denominator (|N(M)|-|N(M\u)|)) none 0, use this to find the node in the multimap above
//    unordered_map<int, int> denominatorAll; //key node, value: (|N(M)|-|N(M\u)|))
    //unordered_map<int, double> suDivedsource; // key: node, value: s_u/|N(M)|
    multimap<double, int> sourceScore_id;
    multimap<double, int> MaxScore_id;

    unordered_map<int, double> id_sourceScore;

//    unordered_set<int> sourceSet;
//    int nonZeo = 0;
//    multimap<double, int> suScore;
//    vector<pair<int, double>> suScore;

    //first loop
    //Count the denominator |N(M)|-|N(M\u)|
    unordered_set<int> maxSet;
    unordered_set<int> nodeSet;

    for(int i = 0; i< numSource; i++){
        int source_id = i + numMiddle;
        double sum_source_i = 0;
        for (auto &it: sorceMid[source_id]){
            sum_source_i += score[it];
        }
        sourceScore_id.insert({sum_source_i,source_id});
        id_sourceScore.insert({source_id,sum_source_i});

//        sourceSet.insert(source_id);
    }

    double sum_mid = 0;
    for(int i = 0; i<numMiddle; i++){
        nodeSet.insert(i);

        sum_mid+=score[i];
    }
    vector<double> new_score_vec;
    int denominator = sorceMid.size();
    double max_numerator = sum_mid;
    double maxScore = (double) max_numerator/denominator;

//    int num_mid = 0;
//    int num_sor = 0;
//
    maxSet = nodeSet;


    double newScore = maxScore;
//    new_score_vec.push_back(log10(newScore));
    double max_denominator = denominator;
//    int cnt_iter = 0;
//    int best_iter = 0;

//    unordered_map <int, unordered_set<int>> best_midSource;
//    unordered_map <int, unordered_set<int>> best_sorceMid;


    while(!sourceScore_id.empty()){
//        cnt_iter++;
//        cout<<cnt_iter<<endl;

        unordered_map<int, double> node_affected; // node-> the score need to be minus
        unordered_set <int> node_erased;

//        cout<<sourceScore_id.size()<<endl;
        double min_source_score = sourceScore_id.begin()->first;
        int min_source_node = sourceScore_id.begin()->second;
        sourceScore_id.erase(sourceScore_id.begin());

        for (auto xx:sorceMid[min_source_node]){
            nodeSet.erase(xx);
        }

        update_mp(midSource, sorceMid,score, min_source_node, sum_mid, denominator, node_affected,node_erased);
        update_score_mp(sourceScore_id,id_sourceScore, node_affected,node_erased);

        if(denominator!=0){
            newScore = (double) sum_mid/ denominator;
//            new_score_vec.push_back(log10(newScore));
        }
        if(newScore >maxScore){
            maxSet = nodeSet;

            max_numerator= sum_mid;
            max_denominator = denominator;
            maxScore = newScore;

            MaxScore_id = sourceScore_id;

//            num_mid = midSource.size();
//            num_sor = sorceMid.size();
//            if(num_mid==2){
//                best_midSource  = midSource;
//            best_sorceMid= sorceMid;
//

        }
//            best_iter = cnt_iter;
    }

//    for (auto i: MaxScore_id){
//        Sourc.insert(i.second);
//    }
    Pickedset = maxSet;


    //    ploter.data.insert({datasetName,new_score_vec});
    return maxScore;

}





bool compareBySecEle(const pair<int, double>& p1, const pair<int, double>& p2){
    return p1.second > p2.second;
}


int updatePA(unordered_map <int, unordered_set<int>> &midSource, int node, vector<int> &PA, int numMiddle){

    int numIsolateSource = 0;
    for(auto &source : midSource[node]){
        if(PA[source-numMiddle] ==1){
            numIsolateSource++;
        }
        PA[source-numMiddle]--;
    }
    midSource.erase(node);
    return numIsolateSource;

}


double FGR(const int& numMiddle, const vector<double>& score, const int &numSource,
           unordered_map <int, unordered_set<int>> &midSource,
           vector<int> &PA, double &scoreFast, int &sourceFast,  int &middleFast, unordered_set<int> &FGR_set){

    vector <int> nodeSort;
    unordered_set<int> nodeSet;
//    unordered_set<int> maxSet;
    //count the score S_u/|N(u)| and sort
    //unordered_map<int, double> suDivedsource;
    //vector<pair<int, double>>
    double sum = 0;
    int nonZeo = 0;
//    unordered_set<int> sourceSet;
    //double suScore[numMiddle] = {0.0};

    vector<pair<int, double>> suScore;

//    suScore.reserve(numMiddle);

    for(int i = 0; i<numMiddle; i++){
        sum+=score[i];

        //int numIsolate = countDenominatorFast(midSource_vec,PA, i);
        //denominatorAll[i] = numIsolate;
        //suScore[i] = (double)score[i]/midSource[i].size();
        if(score[i]!=0){
            nonZeo++;
        }
        nodeSort.push_back(i);
        suScore.push_back({i, (double)score[i]/midSource[i].size()});
        nodeSet.insert(i);

        //}
        //suDivedsource.insert({i,(double)score[i]/midSource[i].size()});
//        for(auto &it :  midSource[i]){
//            sourceSet.insert(it);
//        }
    }
    FGR_set = nodeSet;
    //cout<<"first node suScore "<<suScore[0]<<endl;
    //quickSort(suScore, 0 , numMiddle-1, nodeSort);
    sort(suScore.begin(), suScore.end(), compareBySecEle);

    int i = 0;
    for (const auto& p : suScore) {
        // cout << p.first << ' ' << p.second << endl;
        nodeSort[i] = p.first;
        i++;
    }
    //cout<<"first node "<<nodeSort[0]<<endl;
    int denominator = numSource;

    //nodeSort = sort(suDivedsource);

    //count the value of \sum of M S_u/|N(M)|
    double maxScore = (double) sum/denominator;
    double max_numerator = sum;
    double max_denominator = denominator;

    int curretSize = nodeSort.size();
    int maxNum = nonZeo;
    //cout<< "maxNum = "<<maxNum<<endl;
    double newScore = maxScore;
//    vector < double> new_score_vec;
//    new_score_vec.push_back(newScore);

    //peeling algorithm and remember the set with max value
    while(!nodeSort.empty()){
        int node = nodeSort.back();
        nodeSort.pop_back();
        //	cout<<"remove node "<<node<<endl;

        curretSize--;
        //int numIsolateSource = updateGraph(midSource, sorceMid, node);
        sum = sum - score[node];

        //int numIsolateSource = denominatorAll[node];
        int numIsolateSource = updatePA(midSource, node,PA,numMiddle);
        nodeSet.erase(node);

        if(numIsolateSource != 0){
            denominator = denominator - numIsolateSource;
            if(denominator!=0){
                newScore = (double) sum/ denominator;
//                new_score_vec.push_back(newScore);
            }
            if(newScore >maxScore){
                max_numerator= sum;
                max_denominator = denominator;
                maxScore = newScore;
                maxNum = curretSize;
                FGR_set = nodeSet;


            }

        }

    }

//    cout<<"FGR max set { ";
//    for(auto &it : FGR_set){
//        cout<<it<<", ";
//    }
//    cout<<"}"<<endl;

//    plt::plot(FGR_maxscore_vec,{{"label","FGR"}});
//    plt::title("Score of the FGR objective function - "+ name +" "+ flag);

//    plt::save("./FGR_Score_"+name+"_"+flag+".png");
//    plt::show();

    scoreFast = maxScore;
    sourceFast = max_denominator;
    middleFast = maxNum;




//    cout<<"  MaxScore FGR = "<<maxScore<<endl;
//    ploter.data.insert({datasetName,new_score_vec});

    return maxScore;
}








void IP(const int& numMiddle, vector<double> &score, const int& numSource,
                         unordered_map <int, unordered_set<int>> midSource,
                         unordered_map <int, unordered_set<int>> sorceMid, double & maxScore,vector<double> &L_u, int& t, unordered_set<int> &maxSet){


//    unordered_set<int> nodeSet; //nodes are selected, at the beginning it select all nodes
    //unordered_set<int> nodeNoneZero; //the denominator (|N(M)|-|N(M\u)|)) none 0, use this to find the node in the multimap above
//    unordered_map<int, int> denominatorAll; //key node, value: (|N(M)|-|N(M\u)|))
    //unordered_map<int, double> suDivedsource; // key: node, value: s_u/|N(M)|
    multimap<double, int> sourceScore_id;
    unordered_map<int, double> id_sourceScore;

//    unordered_set<int> sourceSet;
//    int nonZeo = 0;
//    multimap<double, int> suScore;
//    vector<pair<int, double>> suScore;

    //first loop
    //Count the denominator |N(M)|-|N(M\u)|

//    unordered_set<int> maxSet;
    unordered_set<int> nodeSet;

    for(int i = 0; i< numSource; i++){
        int source_id = i + numMiddle;
        double sum_source_i = 0;
        for (auto &it: sorceMid[source_id]){
            sum_source_i += score[it];
        }
        sourceScore_id.insert({sum_source_i+L_u[i],source_id});
        id_sourceScore.insert({source_id,sum_source_i+L_u[i]});

//        sourceSet.insert(source_id);
    }

    double sum_mid = 0;
    for(int i = 0; i<numMiddle; i++){
        nodeSet.insert(i);

        sum_mid+=score[i];
    }
    vector<double> new_score_vec;
    int denominator = sorceMid.size();
    double max_numerator = sum_mid;
    if (t == 1){
        maxSet = nodeSet;

        maxScore = (double) max_numerator/denominator;
    }

//    int num_mid = 0;
//    int num_sor = 0;
//

    double newScore = maxScore;
//    new_score_vec.push_back(log10(newScore));
    double max_denominator = denominator;
//    int cnt_iter = 0;
//    int best_iter = 0;

//    unordered_map <int, unordered_set<int>> best_midSource;
//    unordered_map <int, unordered_set<int>> best_sorceMid;


    while(!sourceScore_id.empty()){
//        cnt_iter++;
//        cout<<cnt_iter<<endl;

        unordered_map<int, double> node_affected; // node-> the score need to be minus
        unordered_set <int> node_erased;

//        cout<<sourceScore_id.size()<<endl;
        double min_source_score = sourceScore_id.begin()->first;
        int min_source_node = sourceScore_id.begin()->second;

        sourceScore_id.erase(sourceScore_id.begin());

        for (auto xx:sorceMid[min_source_node]){
            nodeSet.erase(xx);
        }



//        double sum_mid_tmp = sum_mid;
        update_mp(midSource, sorceMid,score, min_source_node, sum_mid, denominator, node_affected,node_erased);
        L_u[min_source_node -numMiddle]  += min_source_score; // update the l_u according to the marginal value of min_source->second

//        log(min_source.)
//        log(sum_mid_tmp - sum_mid);
        update_score_mp(sourceScore_id,id_sourceScore, node_affected,node_erased);

        if(denominator!=0){
            newScore = (double) sum_mid/ denominator;
//            new_score_vec.push_back(log10(newScore));
        }
        if(newScore >maxScore){
            max_numerator= sum_mid;
            max_denominator = denominator;
            maxScore = newScore;
            maxSet = nodeSet;


        }
//            best_iter = cnt_iter;
    }


//        cout<<"when t = "<<t<<"; IP max set { ";
//        for(auto &it : maxSet){
//            cout<<it<<", ";
//        }
//        cout<<endl;

    //    ploter.data.insert({datasetName,new_score_vec});
//    return maxScore;

}



int countDenominator(unordered_map <int, unordered_set<int>> &midSource,
                     unordered_map <int, unordered_set<int>> &sorceMid, int node){

    int numIsolateSource = 0;

    for(auto &source : midSource[node]){
        if(sorceMid[source].size() == 1){
            //	cout<<"source = "<< source<<endl;
            numIsolateSource++;
        }
    }
    return numIsolateSource;
}




void updateDenominator(unordered_map <int, unordered_set<int>> &midSource,
                       unordered_map <int, unordered_set<int>> &sorceMid, int node,
                       unordered_map<int, int> &denominatorAll,
                       vector<double> &score, multimap<double, int> &noneZeroDenominator){

    int numIsolateSource = 0;
    //unordered_set<int, int> updateNodeDenominator; //key: node, value: |N(M)|-|N(M\u)| update
    vector<int> nodeUpdate;
//    nodeUpdate.reserve
    unordered_set<int> sourceSet;
    for(auto &source : midSource[node]){
        if(sorceMid[source].size() == 1){
            sorceMid.erase(source);
        }
        else{
            sorceMid[source].erase(node);
            if(sorceMid[source].size() == 1){
                //cout<<source<<" sorceMid[source].size()  = "<<sorceMid[source].size() <<endl;
                int node = *sorceMid[source].begin();
                //noneZeroDenominator.erase((double)score[node]/denominatorAll[node]);
                auto iter = denominatorAll.find(node);
                if(iter != denominatorAll.end()){
                    denominatorAll[node]++;
                }
                else{
                    denominatorAll.insert({node,1});
                }
                nodeUpdate.push_back(node);
            }
        }

    }
    midSource.erase(node);
    //cout<<"!!!!!!!!!nodeUpdate size = "<<nodeUpdate.size()<<endl;

    for(auto &it: nodeUpdate){
        int de = denominatorAll[it];
        noneZeroDenominator.insert({(double)score[it]/de, it});
    }

}




double GR(int numMiddle, vector<double> &score, int &numSource,
          unordered_map <int, unordered_set<int>> &midSource,
          unordered_map <int, unordered_set<int>> &sorceMid,
          double &scoreGreedy, int &sourceGreedy, int &middleGreedy, unordered_set<int> &GR_set){

    double sum = 0;
    unordered_set<int> nodeSet; //nodes are selected, at the beginning it select all nodes
    multimap<double, int> noneZeroDenominator; //score for s_u/(|N(M)|-|N(M\u)|)
    //unordered_set<int> nodeNoneZero; //the denominator (|N(M)|-|N(M\u)|)) none 0, use this to find the node in the multimap above
    unordered_map<int, int> denominatorAll; //key node, value: (|N(M)|-|N(M\u)|))
    //unordered_map<int, double> suDivedsource; // key: node, value: s_u/|N(M)|
    unordered_set<int> sourceSet; //source nodes indicate to the selected middle nodes
    int nonZeo = 0;
    vector<pair<int, double>> suScore;
    //first loop
    //Count the denominator |N(M)|-|N(M\u)|
    for(int i = 0; i< numMiddle; i++){
        //if(score[i]!=0){
        sum+=score[i];
        if(score[i]>0){
            nonZeo++;
        }
        int numIsolate = countDenominator(midSource,sorceMid, i);
        if(numIsolate >0){
            //If there exist number that not 0 then put in multimap
            noneZeroDenominator.insert({(double)score[i]/numIsolate, i});
            //nodeNoneZero.insert(i);
            denominatorAll.insert({i,numIsolate});

        }
        else{
            suScore.push_back({i,(double)score[i]/midSource[i].size()});
            //suDivedsource.insert({i,(double)score[i]/midSource[i].size()});
        }
        nodeSet.insert(i);
        //}

        for(auto &it :  midSource[i]){
            sourceSet.insert(it);
        }

    }

    int denominator = sourceSet.size();

    //vector <int> nodeSort = sortMap(suDivedsource);
    sort(suScore.begin(), suScore.end(), compareBySecEle);
    vector<int> nodeSort(suScore.size());
    int i = 0;
    for (const auto& p : suScore) {
        // cout << p.first << ' ' << p.second << endl;
        nodeSort[i] = p.first;
        i++;
    }
    long double maxScore = (double) sum/denominator;

    //unordered_set<int> maxSet;
    //maxSet = nodeSet;
    //int maxSelection = nodeSet.size();
    int maxSelection = nonZeo;
    double max_denominator = denominator;
    double max_numerator = sum;
    //In the first loop count everything, the size of none 0 will increase after deletion
    //it won't decrease
    vector <double> new_score_vec;

    while(!nodeSet.empty()){
        double newScore = 0.0;

        int nodeErase =-1;
        for(auto iter = noneZeroDenominator.begin(); iter!= noneZeroDenominator.end();){
            auto iter_node = nodeSet.find(iter->second);
            if(iter_node != nodeSet.end()){
                nodeErase = iter->second;
                iter = noneZeroDenominator.erase(iter);
                break;
            }
            iter = noneZeroDenominator.erase(iter);
        }

        //noneZeroDenominator.erase(iter);

        if(nodeErase!=-1){
            //cout<<"nodeErase = "<< nodeErase<<endl;
            updateDenominator(midSource, sorceMid, nodeErase,denominatorAll,score,
                              noneZeroDenominator);
            sum -= score[nodeErase];
            denominator -= denominatorAll[nodeErase];
            if(denominator!=0){
                newScore = (long double)sum/ denominator;
                new_score_vec.push_back(newScore);
            }
            //cout<<"denominator = "<<denominator <<"  newScore = "<<newScore<<endl;
            nodeSet.erase(nodeErase);
            //cout<<"node erase = "<<nodeErase<<" , newScore = "<<newScore<<endl;
            if(newScore >maxScore){
                max_numerator = sum;
                max_denominator = denominator;
                maxScore = newScore;
                GR_set = nodeSet;
                maxSelection = nodeSet.size();
            }

        }
            //nodeNoneZero.erase(nodeErase);
            //}
        else{
            while(!nodeSort.empty()){
                int node = nodeSort.back();
                nodeSort.pop_back();

                auto iter = nodeSet.find(node);
                if(iter!= nodeSet.end()){
                    //cout<<"node erase = "<<node<<endl;
                    nodeSet.erase(node);
                    updateDenominator(midSource, sorceMid, node,denominatorAll,score,
                                      noneZeroDenominator);
                    sum = sum - score[node];
                    if(denominator!=0){
                        newScore = (double)sum/ denominator;
                        new_score_vec.push_back(newScore);

                    }
                    if(newScore >maxScore){
                        maxScore = newScore;
                        maxSelection = nodeSet.size();
                        GR_set = nodeSet;
                    }
                    break;
                }
            }
        }


    }

    scoreGreedy = maxScore;
    sourceGreedy = max_denominator;
    middleGreedy = maxSelection;

//    plt::plot(GR_maxscore_vec,{{"label","GR"}});
//    plt::title("Score of the FGR / GR objective function - "+ name+" "+ flag);
//    plt::legend();
//    plt::save("./FGR_GR_Score_"+name+"_"+ flag+".png");
//    plt::show();


//    cout<<"  MaxScore GR = "<<maxScore<<endl;

//    ploter.data.insert({datasetName,new_score_vec});
    return maxScore;


}

int countExtraSourceV2(unordered_map <int, unordered_set<int>> &midSource,
                       unordered_map <int, unordered_set<int>> &sorceMid, int node,
                       unordered_set <int> &sourceSelected){

    int numIsolateSource = 0;

    for(auto &source : midSource[node]){
        auto iter = sourceSelected.find(source);
        if(iter == sourceSelected.end()){
            numIsolateSource++;
        }

    }

    return numIsolateSource;
}



double GRR(int numMiddle, vector<double> &score, int &numSource,
           unordered_map <int, unordered_set<int>> &midSource,
           unordered_map <int, unordered_set<int>> &sorceMid,
           double &scoreBott, int &sourceBott, int &middleBott, unordered_set<int> &GRR_set){

    unordered_set<int> nodeSet;
    int maxDenominator = 0;

    int numSelectedMid = 0;

    double scoreCurt = 0.0;
    int denominator = 0;
    unordered_set <int> sourceSelected;

    for(int i = 0; i<numMiddle; i++){
        if(score[i]!=0){
            nodeSet.insert(i);
        }
    }
    double maxScore = 0.0;
    int counter = 0;
    int maxSelection = 0;
    GRR_set = nodeSet;

    while(!nodeSet.empty()){

        int nodeSelect = -1;
        double maxScoreCurt = DBL_MAX;
        //select one node for current iteration
        for(auto &it: nodeSet){
            //Greg 22/08: what if denominator 0?
            int difference = countExtraSourceV2(midSource,sorceMid, it, sourceSelected);
            double newScore = (double)difference/score[it];
            if(newScore ==0 ){
                nodeSelect = it;
                break;
            }
            //cout<<score[it]<<" "<<it <<" newScore = "<<newScore<<endl;
            if(newScore < maxScoreCurt){
                maxScoreCurt = newScore;
                nodeSelect = it;
            }
        }
        counter++;
        // cout<<counter<<" node select = "<<nodeSelect<<endl;
        //cout<<"maxScoreCurt = "<<maxScoreCurt<<endl;
        for(auto &source : midSource[nodeSelect]){
            sourceSelected.insert(source);
        }
        scoreCurt += score[nodeSelect];
        double currentScore = (double)scoreCurt/sourceSelected.size();
        if(currentScore > maxScore){
            maxScore = currentScore;
            maxSelection = counter;
            maxDenominator = sourceSelected.size();
            GRR_set = nodeSet;

        }

        nodeSet.erase(nodeSelect);
    }

    scoreBott = maxScore;
    sourceBott = maxDenominator;
    middleBott = maxSelection;
    return maxScore;

}



void runLP(string &filename){

#ifdef PRINT
    cout<<"Now run LP on "<<filename<<endl;
#endif

    int sourceSelect = 0;
//    vector <vector<double>> Multi_score;
    double readfile_time = 0.0;
    int numSource = 0;
    int numMiddle = 0;
    int numEdges = 0;
    vector<int> sourceSelectVec;
    vector<int> finalSetILP;

    unordered_map <int, unordered_set<int>> trans;

    vector<double> score = readFile(filename, trans,numSource, numMiddle, numEdges);


    double objRE = 0.0;


    auto ILP_start = std::chrono::high_resolution_clock::now();

    ILP(numMiddle,score, numSource, trans,finalSetILP, objRE,filename, numEdges, sourceSelect,sourceSelectVec);

    auto ILP_end = std::chrono::high_resolution_clock::now();

    double time_ILP = std::chrono::duration_cast < std::chrono::milliseconds > (ILP_end - ILP_start).count();
#ifdef PRINT
        cout<<"Max score = "<<objRE<<endl;
        cout<<"Rutime for LP  = "<<time_ILP<<" ms"<<endl;

#endif





}



void runCD(string &filename){
#ifdef PRINT
    cout<<"Now run CD on "<<filename<<endl;
#endif
    vector<int> finalSetILP;
    double objRE = 0.0;
    vector<int> bestSetILP;


    int sourceSelect = 0;

    int merge_size = 100;
    int merge_bound = 500;

    double readfile_time = 0.0;

    vector<unordered_map<int, unordered_set<int>>> Multi_trans;
    vector<int> Multi_numSource;
    vector<int> Multi_numMiddle;
    vector<int> Multi_numEdges;
    vector<vector<double>> Multi_score;
    unordered_map<int, int> mapping;

    std::tuple<vector<unordered_map<int, unordered_set<int>>>, vector<int>, vector<int>, vector<int>, vector<vector<double>>,unordered_map<int, int>> result = readFile(
            filename, readfile_time, merge_size, merge_bound);


    std::tie(Multi_trans, Multi_numSource, Multi_numMiddle, Multi_numEdges, Multi_score, mapping) = result;

    double best_score = -1.0;
    int best_i = -1;


    auto ILP_start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < Multi_numMiddle.size(); ++i) {
//        cout<<"the size of component "<< i << " is " <<Multi_numMiddle[i]<<endl;

        if (Multi_numMiddle[i] == 1) {
            objRE = Multi_score[i][0]/Multi_numSource[i];
            if (objRE > best_score) {
                best_score = objRE;
                best_i = i;
            }

        } else {
            vector<int> finalSetILP;
            vector<int> sourceSelectVec;

            ILP(Multi_numMiddle[i], Multi_score[i], Multi_numSource[i], Multi_trans[i], finalSetILP, objRE, filename, Multi_numEdges[i], sourceSelect, sourceSelectVec);


            if (objRE > best_score) {
                best_score = objRE;
                best_i = i;
                bestSetILP = finalSetILP;
            }

        }

    }
    auto ILP_end = std::chrono::high_resolution_clock::now();
    double time_ILP = std::chrono::duration_cast<std::chrono::milliseconds>(ILP_end - ILP_start).count();

#ifdef PRINT
    cout << "Max score is:" << best_score << " found in " << best_i +1 << "th component"<<endl;
    cout << "Total time for CD = " << time_ILP + readfile_time << " ms" << endl;

#endif

#ifdef VERBOSE

    cout << "Runtime for ILP: " << time_ILP << " ms" << endl;
    cout << "Runtime for preprocessing and decomposition = " << readfile_time << " ms" << endl;
#endif

}



void runIP(string &filename){
#ifdef PRINT
    cout<<"Now run IP on "<<filename<<endl;
#endif
    int numSource = 0;
    int numMiddle = 0;

    unordered_map <int, unordered_set<int>> midSource;
    unordered_map <int, unordered_set<int>> sorceMid;
    int numEdges = 0;

    vector<double> score = readFile(filename, midSource, numSource, numMiddle, numEdges, sorceMid);



    int T = 100;

    unordered_set<int> maxSet;
    unordered_set<int> maxUSet;

//        double time_greedyAproxV2=0.0;


    Tolerance tolerance;
    double maxScore = -1;
    double newScore;
    vector<double> L_u (numSource,0); // initialize L_u at t = 0

    auto IP_start = std::chrono::high_resolution_clock::now();

    for (int t = 1; t < T; ++t) {

        IP(numMiddle, score, numSource, midSource, sorceMid, newScore, L_u, t, maxSet);


#ifdef VERBOSE
        cout<<"Score is "<<newScore<<" in iteration "<<t<<endl;
#endif
        tolerance.total++;

        if (newScore > maxScore) {
            maxScore = newScore;
            tolerance.cnt = 1;
            tolerance.iteration = t;
        } else {

            tolerance.cnt++;

            if (tolerance.cnt > tolerance.max) {
                break;
            }
        }
    }


    auto IP_end = std::chrono::high_resolution_clock::now();
    double time_IP = std::chrono::duration_cast<std::chrono::milliseconds>(IP_end - IP_start).count();



#ifdef PRINT
    cout << "Rutime for IP = " << time_IP << " ms" << endl;
    cout << "Score for IP = " << maxScore  << endl;
#endif

#ifdef VERBOSE
    cout<< "Max score was achieved in iteration " <<tolerance.iteration<<endl;

#endif

}




void runGAR(string &filename){
#ifdef PRINT
    cout<<"Now run GAR on "<<filename<<endl;
#endif
    int numSource = 0;
    int numMiddle = 0;

    unordered_map<int, unordered_set<int>> midSource;
    unordered_map<int, unordered_set<int>> sorceMid;
//        vector<int> component;
    int numEdges = 0;
    vector<double> score = readFile(filename, midSource, numSource, numMiddle, numEdges, sorceMid);






    int max_degree = -1;
    for (auto const &i : midSource) {
        int degree_i = i.second.size();
        if (degree_i > max_degree){
            max_degree = degree_i;
        }

    }
#ifdef VERBOSE

    cout<<"The 1/max degree(v) in V: "<<(double) 1/max_degree<<endl;
#endif




    unordered_set<int> Pickedset;



    auto GAR_start = std::chrono::high_resolution_clock::now();

    double maxScoreGR = GAR(numMiddle, score, numSource, midSource, sorceMid, Pickedset);

    auto GAR_end = std::chrono::high_resolution_clock::now();

    double time_GAR = std::chrono::duration_cast<std::chrono::milliseconds>(GAR_end - GAR_start).count();

#ifdef PRINT
    cout << "Score for GAR = " << maxScoreGR << endl;
    cout << "Rutime for GAR = " << time_GAR << " ms" << endl;
#endif


}




void runGR(string &filename) {

#ifdef PRINT
    cout<<"Now run GR on "<<filename<<endl;
#endif

    int numSource = 0;
    int numMiddle = 0;

    unordered_map<int, unordered_set<int>> midSource;
    unordered_map<int, unordered_set<int>> sorceMid;
    vector<int> component;

    int numEdges = 0;

    double MaxScoreSeparateU = 0.0;

    auto SeparateU_start = std::chrono::high_resolution_clock::now();

    unordered_set<int> MaxSet;
    vector<double> score = readFile(filename, midSource, numSource, numMiddle, numEdges, sorceMid, MaxScoreSeparateU, MaxSet);
    auto SeparateU_end = std::chrono::high_resolution_clock::now();
    double time_SeparateU = std::chrono::duration_cast < std::chrono::milliseconds > (SeparateU_end - SeparateU_start).count();

#ifdef VERBOSE

    cout<<"Rutime for readfile  = "<<time_SeparateU<<" ms"<<endl;

#endif


    vector<int> PA(numSource);
    PA.reserve(numSource);
    for (int i = 0; i < numSource; i++) {
        PA[i] = sorceMid[i + numMiddle].size();
    }


    double maxscore_GR = 0;
    int denomiGreedy = 0;
    int middleGreedy = 0;




    unordered_set<int> GR_middleset;
//    unordered_set<int> GR_sourceSet;



    auto GR_start = std::chrono::high_resolution_clock::now();


    GR(numMiddle, score, numSource, midSource, sorceMid, maxscore_GR, denomiGreedy, middleGreedy, GR_middleset);
    auto GR_end = std::chrono::high_resolution_clock::now();

    double time_GR = std::chrono::duration_cast<std::chrono::milliseconds>( GR_end - GR_start).count();


    if (MaxScoreSeparateU>maxscore_GR) {
        maxscore_GR = MaxScoreSeparateU;
    }

#ifdef PRINT
    cout<<"MaxScore GR = "<<maxscore_GR<<endl;

    cout<<"Rutime for GR = "<<time_GR<<" ms"<<endl;

#endif





}




void runFGR(string &filename) {

#ifdef PRINT
    cout<<"Now run FGR on "<<filename<<endl;
#endif

    int numSource = 0;
    int numMiddle = 0;

    unordered_map<int, unordered_set<int>> midSourceFast;
    unordered_map<int, unordered_set<int>> sorceMid;

    int numEdges = 0;



    double MaxScoreSeparateU = 0.0;

    auto SeparateU_start = std::chrono::high_resolution_clock::now();

    unordered_set<int> MaxSet;
    vector<double> score = readFile(filename, midSourceFast, numSource, numMiddle, numEdges, sorceMid, MaxScoreSeparateU, MaxSet);
    auto SeparateU_end = std::chrono::high_resolution_clock::now();
    double time_SeparateU = std::chrono::duration_cast < std::chrono::milliseconds > (SeparateU_end - SeparateU_start).count();

#ifdef VERBOSE

    cout<<"Rutime for readfile  = "<<time_SeparateU<<" ms"<<endl;

#endif






    vector<int> PA(numSource);
    PA.reserve(numSource);
    for (int i = 0; i < numSource; i++) {
        PA[i] = sorceMid[i + numMiddle].size();
    }



    unordered_set<int> FGR_set;
    double maxscore_FGR = 0;
    int denomiFast = 0;
    int middleFast = 0;

    auto FGR_start = std::chrono::high_resolution_clock::now();

    FGR(numMiddle, score, numSource, midSourceFast, PA, maxscore_FGR, denomiFast, middleFast, FGR_set);
    auto FGR_end = std::chrono::high_resolution_clock::now();

    double time_FGR = std::chrono::duration_cast<std::chrono::milliseconds>(FGR_end - FGR_start).count();


    if (MaxScoreSeparateU>maxscore_FGR) {

        maxscore_FGR = MaxScoreSeparateU;

    }

#ifdef PRINT

    cout<<"MaxScore FGR = "<<maxscore_FGR<<endl;

    cout<<"Rutime for FGR  = "<<time_FGR<<" ms"<<endl;


#endif



}




void runGRR(string &filename) {

#ifdef PRINT
    cout<<"Now run GRR on "<<filename<<endl;
#endif

    int numSource = 0;
    int numMiddle = 0;
    int numEdges = 0;

    unordered_map<int, unordered_set<int>> midSourceBottom;
    unordered_map<int, unordered_set<int>> sorceMidBottom;

    vector<double> score = readFile(filename, midSourceBottom, numSource, numMiddle, numEdges, sorceMidBottom);



    double maxScore_GRR = 0;
    int denomiBott = 0;
    int middleBott = 0;



    unordered_set<int> GRR_set;
    auto bottom_up_start = std::chrono::high_resolution_clock::now();


    GRR(numMiddle,score, numSource, midSourceBottom, sorceMidBottom, maxScore_GRR, denomiBott, middleBott, GRR_set);

    auto bottom_up_end = std::chrono::high_resolution_clock::now();

    double time_bottom = std::chrono::duration_cast < std::chrono::milliseconds > (bottom_up_end - bottom_up_start).count();
#ifdef PRINT
    cout << "GRR score = " << maxScore_GRR << endl;

    cout<<"Rutime for GRR  = "<<time_bottom<<" ms"<<endl;

#endif



}





int main(int argc, char *argv[]){

    cmdline::parser parser;
    // 1st argument is long name
    // 2nd argument is short name (no short name if '\0' specified)
    // 3rd argument is description
    // 4th argument is mandatory (optional. default is false)
    // 5th argument is default value  (optional. it used when mandatory is false)

    parser.add<string>("filepath", 'f', "the path to input graph", true, "");
    //available options of methodFlag = {LP, CD, IP GAR, GR, FGR, GRR}, default as LP
    parser.add<string>("method", 'm', "choose a method", false, "LP", cmdline::oneof<string>("LP", "CD", "IP", "GAR", "GR", "FGR", "GRR"));




    parser.parse_check(argc, argv);
    string filename = parser.get<string>("filepath");
    string methodFlag = parser.get<string>("method");

    if (methodFlag == "LP"){
        runLP(filename);

    } else if (methodFlag == "CD"){

        runCD(filename);
    }else if (methodFlag == "IP"){
        runIP(filename);

    }else if (methodFlag == "GAR"){
        runGAR(filename);

    }else if (methodFlag == "GR"){
        runGR(filename);

    }else if (methodFlag == "FGR"){
        runFGR(filename);

    }else if (methodFlag == "GRR"){
        runGRR(filename);

    }




    return 0;
}

