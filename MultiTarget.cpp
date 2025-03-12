//
// Created by ling on 03/01/24.
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
#include "cmdline.h"

using namespace std;
using namespace chrono;
using namespace boost::heap;
#include <filesystem>
//namespace fs = std::filesystem;
#define log(x) std::cout<<x<<endl
#define VERBOSE
#define PRINT

struct Tolerance{
    int max = 5;
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



std::tuple<vector<unordered_map<int, unordered_set<int>>>, vector<int> ,vector<int>, vector<int>,vector<vector<pair<int,int>>>,
        vector<unordered_map <int, unordered_set<int>>> ,vector<vector<pair<int,int>>>, vector<vector<double>>> readFile(string &filename, double &readtime,int &merge_size ,int& merge_bound){
    string line;
    ifstream myfile (filename);
    int numSource =0, numMiddle = 0, numEdges = 0, numTarget =0;


    unordered_map <int, unordered_set<int>> trans;
    vector<pair<int,int>> sourceMiddle;
    vector<pair<int,int>> middleTarget;
    unordered_map <int, unordered_set<int>> middleTargetMap;

    if (myfile.is_open()) {
        getline(myfile, line);
        numSource = stoi(line);
        // cout<<"numSource = "<<numSource<<endl;
        getline(myfile, line);
        numMiddle = stoi(line);

        getline(myfile, line);
        numTarget = stoi(line);

        getline(myfile, line);
        numEdges = stoi(line);
#ifdef VERBOSE
        cout<<"numSource = "<<numSource<<", numMiddle  = "<<numMiddle<<" ,  numTarget = "<<numTarget<<endl;
#endif

        vector<double> score(numMiddle);

        int counter = 0;

        //cout << line << '\n';
        while (counter < numMiddle) {
            getline(myfile, line);
            score[counter] = stof(line);
            counter++;
        }


        counter = 0;


        using namespace boost;
        {

            using std::unordered_set;
            using std::unordered_map;

            typedef adjacency_list<vecS, vecS, undirectedS> Graph;
            Graph G;
            double small_time_sum=0.0;

            while (counter < numEdges) {
                int inNode, outNode;

                //getline (myfile,line) ;
                //cout<<"The line is  : "<<line<<endl;
                myfile >> inNode >> outNode;
                auto small_cnt_start = std::chrono::high_resolution_clock::now();

//                inNode = inNode + numMiddle;
                // cout<<"here "<<inNode<<", "<<outNode<<endl;
                auto it = trans.find(outNode);
                if (it != trans.end()) {
                    trans[outNode].insert(inNode);
                    add_edge(outNode, inNode, G);
                } else {
                    unordered_set<int> temp;
                    temp.insert(inNode);
                    trans.insert({outNode, temp});
                    add_edge(outNode, inNode, G);

                }
                counter++;
                sourceMiddle.push_back({inNode, outNode});
                auto small_cnt_end = std::chrono::high_resolution_clock::now();

                auto small_time=  std::chrono::duration_cast < std::chrono::milliseconds > (small_cnt_end - small_cnt_start).count();
                small_time_sum += small_time;


            }

            int inNode, outNode;
            while (myfile >> inNode >> outNode) {
                auto small_cnt_start = std::chrono::high_resolution_clock::now();


                auto it = middleTargetMap.find(inNode);
                if (it != middleTargetMap.end()) {
                    middleTargetMap[inNode].insert(outNode);
                    add_edge(outNode, inNode, G);
                } else {
                    unordered_set<int> temp;
                    temp.insert(outNode);
                    middleTargetMap.insert({inNode, temp});
                    add_edge(outNode, inNode, G);
                }

                middleTarget.push_back({inNode, outNode});

                auto small_cnt_end = std::chrono::high_resolution_clock::now();

                auto small_time=  std::chrono::duration_cast < std::chrono::milliseconds > (small_cnt_end - small_cnt_start).count();
                small_time_sum += small_time;
                //cout<<"inNode = "<<inNode <<", outNode = "<<outNode<<endl;
            }
            // cout<<"middleTarget = "<<middleTarget.size()<<endl;


            auto decomposition_time_start = std::chrono::high_resolution_clock::now();

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
                if (0<Multi_numMiddle_old[i] && Multi_numMiddle_old[i]<merge_size){
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


            vector<unordered_map<int, unordered_set<int>>> Multi_trans(new_num);
            vector<int> Multi_numSource(new_num);
            vector<int> Multi_numMiddle(new_num);
            vector<int> Multi_numTarget(new_num);
            vector<vector<pair<int,int>>> Multi_middleTarget(new_num);
            vector<unordered_map <int, unordered_set<int>>> Multi_middleTargetMap(new_num);
            vector<vector<pair<int,int>>> Multi_sourceMiddle(new_num);
            std::unordered_map<int, int> mapping_relation_decomposition;
            vector<vector<double>> Multi_score(new_num);


            vector<int> current_counts_per_component(new_num);


            for (int i = 0; i < numMiddle; ++i){
                int componentId_new = component_mapping[component[i]];
                Multi_numMiddle[componentId_new]++;
                int cnt = current_counts_per_component[componentId_new];
                mapping_relation_decomposition[i]=cnt;
                current_counts_per_component[componentId_new]++;
//                Multi_numEdges[componentId_new] += trans[i].size();
            }

            for (int i = numMiddle; i != numMiddle + numSource;++i)
            {
                //build the mapping relation for source nodes
                int componentId_new = component_mapping[component[i]];
                int cnt=current_counts_per_component[componentId_new];
                mapping_relation_decomposition[i]=cnt;
//                cout << i <<":"<<cnt<<endl;
                current_counts_per_component[componentId_new]++;
                Multi_numSource[componentId_new]++;
            }

            for (int i = numMiddle + numSource; i != numMiddle + numSource + numTarget;++i)
            {
                //build the mapping relation for source nodes
                int componentId_new = component_mapping[component[i]];
                int cnt=current_counts_per_component[componentId_new];
                mapping_relation_decomposition[i]=cnt;
//                cout << i <<":"<<cnt<<endl;
                current_counts_per_component[componentId_new]++;
                Multi_numTarget[componentId_new]++;
            }





            vector<bool> component_flags(num, false);
            for (int i = 0; i < numMiddle; ++i)
            {
                //initialization Multi_score
                int componentId_new = component_mapping[component[i]];

                if (!component_flags[componentId_new]){
                    vector<double> multi_score_for_component(Multi_numMiddle[componentId_new]) ;
                    Multi_score[componentId_new]=multi_score_for_component;
                    component_flags[componentId_new] = true;
                }
                //accumulating the scores for all nodes after combination

                Multi_score[componentId_new][mapping_relation_decomposition[i]] += score[i];


                std::unordered_set<int> tmp;

                std::unordered_set<int> tar_tmp;

                for(auto &it2 : trans[i])  //create a set with the new ids
                {
                    tmp.insert(mapping_relation_decomposition[it2]);
                    Multi_sourceMiddle[componentId_new].push_back({mapping_relation_decomposition[it2], mapping_relation_decomposition[i]});

                }

                for(auto &it2 : middleTargetMap[i])  //create a set with the new ids
                {
                    tar_tmp.insert(mapping_relation_decomposition[it2]);
                    Multi_middleTarget[componentId_new].push_back({mapping_relation_decomposition[i], mapping_relation_decomposition[it2]});

                }
                Multi_middleTargetMap[componentId_new][mapping_relation_decomposition[i]] = tar_tmp;
                Multi_trans[componentId_new][mapping_relation_decomposition[i]]=tmp; //add it into the correct unordered_map
//                }
            }



            auto decomposition_time_end = std::chrono::high_resolution_clock::now();


#ifdef VERBOSE

            cout<<"After merging, the number of components: "<<new_num<<endl;
#endif


            readtime= std::chrono::duration_cast < std::chrono::milliseconds > (decomposition_time_end - decomposition_time_start).count() + small_time_sum;





            myfile.close();
            return std::make_tuple(Multi_trans,Multi_numSource, Multi_numMiddle, Multi_numTarget, Multi_middleTarget,Multi_middleTargetMap, Multi_sourceMiddle,Multi_score);
        }
    }
    else cout << "Unable to open file";

}


vector<double> readFile(string filename, int &numSource, int &numMiddle,
                        unordered_map <int, unordered_set<int>> &midSource,
                        unordered_map <int, unordered_set<int>> &sorceMid,double &MaxScoreSeparateV){
    string line;
    ifstream myfile (filename);
    int numTarget;
    int numSource_tmp;
    if (myfile.is_open()){
        getline (myfile,line);
        numSource_tmp =  stoi(line);
        // cout<<"numSource = "<<numSource<<endl;
        getline (myfile,line);
        numMiddle =  stoi(line);

        getline (myfile,line);
        numTarget =  stoi(line);
        numSource = numSource_tmp+ numTarget;
        getline (myfile,line);
        int numEdges = stoi(line);


        vector<double> score(numMiddle);

        int counter = 0;

        //cout << line << '\n';
        while(counter<numMiddle){
            getline (myfile,line) ;
            score[counter] = stof(line);
            counter++;
        }


        counter = 0;
        while(counter<numEdges){
            int inNode, outNode;


            myfile >>inNode>>outNode;
            // cout<<"here "<<inNode<<", "<<outNode<<endl;
            auto it =  midSource.find(outNode);
            if(it!= midSource.end()){
                midSource[outNode].insert(inNode);
//                MidDegreeOne[outNode]=0;  // set value as o if degree >1

            }
            else{
                unordered_set<int> temp;
                temp.insert(inNode);
                midSource.insert({outNode, temp});

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

            counter++;

        }

//        cout<<sorceMid.size()<<endl;
//        cout<<midSource.size()<<endl;

        int inNode, outNode;
        while(myfile >>inNode>>outNode){
            //cout<<"inNode = "<<inNode<<", outNode = "<<outNode<<endl;
//            outNode = outNode + numSource_tmp;

            numEdges++;
            auto it =  midSource.find(inNode);
            if(it!= midSource.end()){
                midSource[inNode].insert(outNode);

            }
            else{
                unordered_set<int> temp;
                temp.insert(outNode);
                midSource.insert({inNode, temp});

            }

            auto iter_mid =  sorceMid.find(outNode);
            if(iter_mid!= sorceMid.end()){
                sorceMid[outNode].insert(inNode);
                //cout<<"this case insert "<<outNode<<" to "<<inNode<<endl;
            }
            else{
                unordered_set<int> temp;
                temp.insert(inNode);
                sorceMid.insert({outNode, temp});
                //cout<<"insert "<<outNode<<" to "<<inNode<<endl;


            }
            //cout<<"inNode = "<<inNode <<", outNode = "<<outNode<<endl;
        }

        for (auto &i: midSource){
            double tmp_score = (double)score[i.first]/i.second.size();
            if (tmp_score> MaxScoreSeparateV){
                MaxScoreSeparateV = tmp_score;
            }
        }

        myfile.close();
        return score;
    }

    else cout << "Unable to open file";

    vector<double> score;

    return score;

}




vector<double> readFile(string filename, unordered_map <int, unordered_set<int>> &trans,
                              int &numSource, int &numMiddle, int &numTarget,
                              vector<pair<int,int>> &middleTarget, unordered_map <int, unordered_set<int>> &middleTargetMap,
                              vector<pair<int,int>> &sourceMiddle){
    string line;
    ifstream myfile (filename);
    if (myfile.is_open()){
        getline (myfile,line);
        numSource =  stoi(line);
        // cout<<"numSource = "<<numSource<<endl;
        getline (myfile,line);
        numMiddle =  stoi(line);

        getline (myfile,line);
        numTarget =  stoi(line);

        getline (myfile,line);
        int numEdges = stoi(line);
#ifdef VERBOSE
        cout<<"numSource = "<<numSource<<", numMiddle  = "<<numMiddle<<" ,  numTarget = "<<numTarget<<endl;

#endif

        vector<double> score(numMiddle);

        int counter = 0;

        //cout << line << '\n';
        while(counter<numMiddle){
            getline (myfile,line) ;
            score[counter] = stof(line);
            counter++;
        }


        counter = 0;
        while(counter<numEdges){
            int inNode, outNode;

            //getline (myfile,line) ;
            //cout<<"The line is  : "<<line<<endl;
            myfile >>inNode>>outNode;
            // cout<<"here "<<inNode<<", "<<outNode<<endl;
            auto it =  trans.find(outNode);
            if(it!= trans.end()){
                trans[outNode].insert(inNode);
            }
            else{
                unordered_set<int> temp;
                temp.insert(inNode);
                trans.insert({outNode, temp});

            }
            counter++;
            sourceMiddle.push_back({inNode, outNode});

        }

        int inNode, outNode;
        while(myfile >>inNode>>outNode){
            auto it =  middleTargetMap.find(inNode);
            if(it!= middleTargetMap.end()){
                middleTargetMap[inNode].insert(outNode);
            }
            else{
                unordered_set<int> temp;
                temp.insert(outNode);
                middleTargetMap.insert({inNode, temp});
            }

            middleTarget.push_back({inNode, outNode});
        }
        myfile.close();
        return score;
    }

    else cout << "Unable to open file";

    vector<double> score;

    return score;

}



vector<double> readFile(string filename, int &numSource, int &numMiddle,
                               unordered_map <int, unordered_set<int>> &midSource,
                               unordered_map <int, unordered_set<int>> &sorceMid){
    string line;
    ifstream myfile (filename);
    int numTarget;
    int numSource_tmp;
    if (myfile.is_open()){
        getline (myfile,line);
        numSource_tmp =  stoi(line);
        // cout<<"numSource = "<<numSource<<endl;
        getline (myfile,line);
        numMiddle =  stoi(line);

        getline (myfile,line);
        numTarget =  stoi(line);
        numSource = numSource_tmp+ numTarget;
        getline (myfile,line);
        int numEdges = stoi(line);


        vector<double> score(numMiddle);

        int counter = 0;

        //cout << line << '\n';
        while(counter<numMiddle){
            getline (myfile,line) ;
            score[counter] = stof(line);
            counter++;
        }

        counter = 0;
        while(counter<numEdges){
            int inNode, outNode;


            myfile >>inNode>>outNode;
            // cout<<"here "<<inNode<<", "<<outNode<<endl;
            auto it =  midSource.find(outNode);
            if(it!= midSource.end()){
                midSource[outNode].insert(inNode);
            }
            else{
                unordered_set<int> temp;
                temp.insert(inNode);
                midSource.insert({outNode, temp});

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

            counter++;

        }


//        cout<<sorceMid.size()<<endl;
//        cout<<midSource.size()<<endl;

        int inNode, outNode;
        while(myfile >>inNode>>outNode){
            //cout<<"inNode = "<<inNode<<", outNode = "<<outNode<<endl;
//            outNode = outNode + numSource_tmp;
            numEdges++;
            auto it =  midSource.find(inNode);
            if(it!= midSource.end()){
                midSource[inNode].insert(outNode);

            }
            else{
                unordered_set<int> temp;
                temp.insert(outNode);
                midSource.insert({inNode, temp});

            }

            auto iter_mid =  sorceMid.find(outNode);
            if(iter_mid!= sorceMid.end()){
                sorceMid[outNode].insert(inNode);
                //cout<<"this case insert "<<outNode<<" to "<<inNode<<endl;
            }
            else{
                unordered_set<int> temp;
                temp.insert(inNode);
                sorceMid.insert({outNode, temp});
                //cout<<"insert "<<outNode<<" to "<<inNode<<endl;


            }
            //cout<<"inNode = "<<inNode <<", outNode = "<<outNode<<endl;
        }



/*
        auto a = *max_element(score.begin(), score.end());
        auto b = *min_element(score.begin(), score.end());

        int max_degree= 0;
        int min_degree= 9999999;
        int sum = 0;
        for (auto &i: midSource){
            int t = i.second.size();
            if (t> max_degree){
                max_degree = t;
            }
            if (t< min_degree){
                min_degree = t;
            }
            sum += t;
        }
        double avg = (double) sum/midSource.size();


*/






//        cout<<"After combining T and U, numSource = "<<numSource<<", numMiddle  = "<<numMiddle<<endl;

//        cout<<"numEdges = "<<numEdges<<endl;
        // cout<<"middleTarget = "<<middleTarget.size()<<endl;
        myfile.close();
        return score;
    }

    else cout << "Unable to open file";

    vector<double> score;

    return score;

}





vector<double> readFile(string filename, int &numSource, int &numMiddle, int &numTarget,
                        unordered_map <int, unordered_set<int>> &midSource,
                        unordered_map <int, unordered_set<int>> &sorceMid,
                        unordered_map <int, unordered_set<int>> &midTarget,
                        unordered_map <int, unordered_set<int>> &targetMid){
    string line;
    ifstream myfile (filename);
    if (myfile.is_open()){
        getline (myfile,line);
        numSource =  stoi(line);
        // cout<<"numSource = "<<numSource<<endl;
        getline (myfile,line);
        numMiddle =  stoi(line);

        getline (myfile,line);
        numTarget =  stoi(line);

        getline (myfile,line);
        int numEdges = stoi(line);
#ifdef VERBOSE

        cout<<"numSource = "<<numSource<<", numMiddle  = "<<numMiddle<<" ,  numTarget = "<<numTarget<<endl;
#endif

        vector<double> score(numMiddle);

        int counter = 0;

        //cout << line << '\n';
        while(counter<numMiddle){
            getline (myfile,line) ;
            score[counter] = stof(line);
            counter++;
        }

        counter = 0;
        while(counter<numEdges){
            int inNode, outNode;


            myfile >>inNode>>outNode;
            // cout<<"here "<<inNode<<", "<<outNode<<endl;
            auto it =  midSource.find(outNode);
            if(it!= midSource.end()){
                midSource[outNode].insert(inNode);
            }
            else{
                unordered_set<int> temp;
                temp.insert(inNode);
                midSource.insert({outNode, temp});

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

            counter++;

        }


        int inNode, outNode;
        while(myfile >>inNode>>outNode){
            //cout<<"inNode = "<<inNode<<", outNode = "<<outNode<<endl;
            auto it =  targetMid.find(outNode);
            if(it!= targetMid.end()){
                targetMid[outNode].insert(inNode);
            }
            else{
                unordered_set<int> temp;
                temp.insert(inNode);
                targetMid.insert({outNode, temp});

            }

            auto iter_mid =  midTarget.find(inNode);
            if(iter_mid!= midTarget.end()){
                midTarget[inNode].insert(outNode);
                //cout<<"this case insert "<<outNode<<" to "<<inNode<<endl;
            }
            else{
                unordered_set<int> temp;
                temp.insert(outNode);
                midTarget.insert({inNode, temp});
                //cout<<"insert "<<outNode<<" to "<<inNode<<endl;
            }
            //cout<<"inNode = "<<inNode <<", outNode = "<<outNode<<endl;
        }
        // cout<<"middleTarget = "<<middleTarget.size()<<endl;
        myfile.close();
        return score;
    }

    else cout << "Unable to open file";

    vector<double> score;

    return score;

}





bool compareBySecEle(const pair<int, double>& p1, const pair<int, double>& p2){
    return p1.second > p2.second;
}

void ILP(int numMiddle, vector<double> &score, int numSource, int numTarget,
              unordered_map <int, unordered_set<int>> &trans, vector<int> &finalSetILP,
              vector<pair<int,int>> &middleTarget, vector<pair<int,int>> &sourceMiddle,
              double &objRE){
    try {

        GRBEnv env = GRBEnv(true);
        env.set("LogFile", "ilp.log");
//        env.set(GRB_IntParam_OutputFlag,0);

        env.start();

        // Create an empty model
        GRBModel model = GRBModel(env);

        // Create variables
        double *w{ new double[numMiddle]{} };

        double *z_ub{ new double[numMiddle]{} };
        double *v_ub{ new double[numSource]{} };
        double *s_ub{ new double[numTarget]{} };

        char* z_type = new char[numMiddle];
        char* v_type = new char[numSource];
        char* s_type = new char[numTarget];

        for(int i = 0; i< numMiddle; i++){
            z_type[i] = GRB_CONTINUOUS;
            z_ub[i] = 1.0;
            w[i]=score[i];
        }

        for(int i = 0; i< numSource; i++){
            v_type[i] = GRB_CONTINUOUS;
            v_ub[i] = 1.0;
        }

        for(int i = 0; i< numTarget; i++){
            s_type[i] = GRB_CONTINUOUS;
            s_ub[i] = 1.0;
        }

        GRBVar *v = model.addVars(NULL,v_ub, NULL, v_type, NULL, numSource);
        GRBVar *z = model.addVars(NULL,z_ub, NULL, z_type, NULL, numMiddle);
        GRBVar *s = model.addVars(NULL,s_ub, NULL, s_type, NULL, numTarget);

        GRBLinExpr  myexpr = 0;
        GRBLinExpr  sumV = 0;
        GRBLinExpr  sumS = 0;
        //  GRBLinExpr  temp_expr = 0;
        // GRBLinExpr  ysum = 0;
        for(int i=0;i<numMiddle;++i){
            myexpr += z[i]*w[i];
        }

        for(int i=0;i<numSource;++i){
            sumV += v[i];
            // ysum+=y[i];
        }

        for(int i=0;i<numTarget;++i){
            sumS += s[i];
            // ysum+=y[i];
        }


        model.addConstr(sumV+sumS -1 == 0);


//        for(auto &it: sourceMiddle){
//            model.addConstr(z[it.second]-v[it.first-numMiddle]<= 0);
//        }
//
//        for(auto &it: middleTarget){
//            model.addConstr(z[it.first]-s[it.second-numMiddle-numSource]<= 0);
//        }

        for(auto &it: sourceMiddle){
            model.addConstr(z[it.second]-v[it.first-numMiddle]<= 0);
        }

        for(auto &it: middleTarget){
            model.addConstr(z[it.first]-s[it.second-numMiddle- numSource]<= 0);
        }

        // Set objective: maximize x + y + 2 z
        model.setObjective(myexpr, GRB_MAXIMIZE);
        // Optimize model
        model.optimize();

        //  cout<<"List of Middle nodes"<<endl;
        for(int i = 0; i< numMiddle ;i++){
            if(z[i].get(GRB_DoubleAttr_X)!=0){
                //	cout<<"i = "<<i<<", "<<z[i].get(GRB_DoubleAttr_X)<<endl;
                finalSetILP.push_back(i);
            }

        }


/*
        cout<<"CD max set { ";
        for(auto &it : finalSetILP){
            cout<<it<<", ";
        }
        cout<<"}"<<endl;

*/
        delete []z;
        delete []v;
        delete []s;
        delete []w;
        delete []z_ub;
        delete []v_ub;
        delete []s_ub;

        delete []v_type;
        delete []z_type;
        delete []s_type;


        objRE = model.get(GRB_DoubleAttr_ObjVal);
//        cout << "MaxScore for LP: " << objRE << endl;

    } catch(GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch(...) {
        cout << "Exception during optimization" << endl;
    }
}






int updateGraph(unordered_map <int, unordered_set<int>> &midSource,
                unordered_map <int, unordered_set<int>> &sorceMid, int node,
                unordered_map <int, unordered_set<int>> &midTarget,
                unordered_map <int, unordered_set<int>> &targetMid){

    int numIsolateSourceTarget = 0;

    for(auto &source : midSource[node]	){
        if(sorceMid[source].size() == 1){
            //	cout<<"source = "<< source<<endl;
            numIsolateSourceTarget++;
            sorceMid.erase(source);
        }
        sorceMid[source].erase(node);
    }
    midSource.erase(node);

    for(auto &target:midTarget[node] ){
        if(targetMid[target].size() == 1){
            numIsolateSourceTarget++;
            targetMid.erase(target);
        }
        targetMid[target].erase(node);
    }
    midTarget.erase(node);


    return numIsolateSourceTarget;
}



double FGR_MT(int numMiddle, vector<double> &score, int numSource,
            unordered_map <int, unordered_set<int>> &midSource,
            unordered_map <int, unordered_set<int>> &sorceMid,
            unordered_map <int, unordered_set<int>> &midTarget,
            unordered_map <int, unordered_set<int>> &targetMid){

    vector <int> nodeSort(numMiddle);
    unordered_set<int> nodeSet;
    unordered_set<int> maxSet;


    //count the score S_u/|N(u)| and sort
    //unordered_map<int, double> suDivedsource;
    vector<pair<int, double>> suScore;

    double sum = 0;
//    unordered_set<int> sourceSet; //store the souce incident to selected middle nodes
//    unordered_set<int> targetSet; //store the target incident to selected middle nodes
    //count the score fist and then sort by decreasing order
    for(int i = 0; i<numMiddle ; i++){
        sum+=score[i];
        //if(score[i]!=0){
        suScore.push_back({i, (double)score[i]/(midSource[i].size()+midTarget[i].size())});

        //suDivedsource.insert({i,(double)score[i]/(midSource[i].size()+midTarget[i].size())});
//        for(auto &it :  midSource[i]){
//            sourceSet.insert(it);
//        }
//        for(auto &it :  midTarget[i]){
//            targetSet.insert(it);
//        }

        nodeSet.insert(i);

    }
    maxSet = nodeSet;

    sort(suScore.begin(), suScore.end(), compareBySecEle);
    int i = 0;
    for (const auto& p : suScore) {
        // cout << p.first << ' ' << p.second << endl;
        nodeSort[i] = p.first;
        i++;
    }

//    int denominator = sourceSet.size()+ targetSet.size();

    int denominator = sorceMid.size()+ targetMid.size();

    double maxScore = (double) sum/denominator;

    int curretSize = nodeSort.size();
    int maxNum = nodeSet.size();
    int max_denominator = denominator;
    //maxSet = nodeSet;

    //peeling algorithm and remember the set with max value
    while(!nodeSort.empty()){
        int node = nodeSort.back();
        int numIsolateSourceTarget = updateGraph(midSource, sorceMid, node, midTarget, targetMid);
        curretSize--;
        sum = sum - score[node];
        denominator = denominator - numIsolateSourceTarget;
        double newScore = (double) sum/ denominator;
        //cout<<"remove node "<<node<<" "<<numIsolateSourceTarget<<" sources and target will be removed, new score = "<< newScore<<endl;
        nodeSet.erase(node);
        if(newScore > maxScore){
            maxScore = newScore;
//            maxSet = nodeSet;
            max_denominator = denominator;
            maxNum = curretSize;

        }

        nodeSort.pop_back();
    }

//    cout<<"max set { ";
//    for(auto &it : maxSet){
//        cout<<it<<", ";
//    }
//    cout<<"}"<<endl;
//
//
//    cout<<"maxScore adapted FGR = "<<maxScore<<endl;
    return maxScore;

}




int countDenominatorMT(unordered_map <int, unordered_set<int>> &midSource,
                       unordered_map <int, unordered_set<int>> &sorceMid, int node,
                       unordered_map <int, unordered_set<int>> &midTarget,
                       unordered_map <int, unordered_set<int>> &targetMid){

    int numIsolateSourceTarget = 0;

    for(auto &source : midSource[node]){
        if(sorceMid[source].size() == 1){
            numIsolateSourceTarget++;
        }
    }

    for(auto &target : midTarget[node]){
        if(targetMid[target].size() == 1){
            numIsolateSourceTarget++;
        }
    }

    return numIsolateSourceTarget;
}




void updateDenominatorMT(unordered_map <int, unordered_set<int>> &midSource,
                         unordered_map <int, unordered_set<int>> &sorceMid,
                         unordered_map <int, unordered_set<int>> &midTarget,
                         unordered_map <int, unordered_set<int>> &targetMid, int node,
                         unordered_map<int, int> &denominatorAll,
                         vector<double> &score, multimap<double, int> &noneZeroDenominator){

    int numIsolateSource = 0;
    //unordered_set<int, int> updateNodeDenominator; //key: node, value: |N(M)|-|N(M\u)| update
    vector<int> nodeUpdate;
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

    for(auto &target : midTarget[node]){
        if(targetMid[target].size() == 1){
            targetMid.erase(target);
        }
        else{
            targetMid[target].erase(node);
            if(targetMid[target].size() == 1){
                //cout<<source<<" sorceMid[source].size()  = "<<sorceMid[source].size() <<endl;
                int node = *targetMid[target].begin();
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
    midTarget.erase(node);

    for(auto &it: nodeUpdate){
        int de = denominatorAll[it];
        noneZeroDenominator.insert({(double)score[it]/de, it});
    }

}



double GR_MT(int numMiddle, vector<double> &score, int numSource,
           unordered_map <int, unordered_set<int>> &midSource,
           unordered_map <int, unordered_set<int>> &sorceMid,
           unordered_map <int, unordered_set<int>> &midTarget,
           unordered_map <int, unordered_set<int>> &targetMid){


    double sum = 0;

    unordered_set<int> nodeSet; //nodes are selected, at the beginning it select all nodes
    multimap<double, int> noneZeroDenominator; //score for s_u/(|N(M)|-|N(M\u)|)
    //unordered_set<int> nodeNoneZero; //the denominator (|N(M)|-|N(M\u)|)) none 0, use this to find the node in the multimap above
    unordered_map<int, int> denominatorAll; //key node, value: (|N(M)|-|N(M\u)|))
    unordered_map<int, double> suDivedsource; // key: node, value: s_u/|N(M)|

    unordered_set<int> sourceSet; //source nodes indicate to the selected middle nodes
    unordered_set<int> targetSet; //target nodes indicate to the selected middle nodes

    vector<pair<int, double>> suScore;
    //first loop
    //Count the denominator |N(M)|-|N(M\u)|
    for(int i = 0; i< numMiddle; i++){
        sum+=score[i];
        //if(score[i]!=0){
        int numIsolate = countDenominatorMT(midSource,sorceMid, i, midTarget, targetMid);
        if(numIsolate >0){
            //If there exist number that not 0 then put in multimap
            noneZeroDenominator.insert({(double)score[i]/numIsolate, i});
            //nodeNoneZero.insert(i);
            denominatorAll.insert({i,numIsolate});
        }
        else{
            suScore.push_back({i,(double)score[i]/(midSource[i].size()+midTarget[i].size())});
            //suDivedsource.insert({i,(double)score[i]/(midSource[i].size()+midTarget[i].size())});
        }

        for(auto &it :  midSource[i]){
            sourceSet.insert(it);
        }

        for(auto &it :  midTarget[i]){
            targetSet.insert(it);
        }
        nodeSet.insert(i);
        //}
    }
    unordered_set<int> maxSet;
    maxSet = nodeSet;


    int denominator = sourceSet.size()+ targetSet.size();
    sort(suScore.begin(), suScore.end(), compareBySecEle);
    vector<int> nodeSort(numMiddle);
    int i = 0;
    for (const auto& p : suScore) {
        // cout << p.first << ' ' << p.second << endl;
        nodeSort[i] = p.first;
        i++;
    }

    //vector <int> nodeSort = sort(suDivedsource);
    double maxScore = (double) sum/denominator;

    int maxSelection = nodeSet.size();
    int max_denominator = denominator;

    //In the first loop count everything, the size of none 0 will increase after deletion
    //it won't decrease
    double newScore = 0.0;
    while(!nodeSet.empty()){
        //if the middle node incident to the unique source node
        //we find the middle node with minimum score of S_u/(|N(M)|-|N(M\u)|)
        //if(!noneZeroDenominator.empty()){
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


        if(nodeErase!=-1){
            updateDenominatorMT(midSource, sorceMid, midTarget, targetMid, nodeErase,denominatorAll,score,
                                noneZeroDenominator);
            sum -= score[nodeErase];
            denominator -= denominatorAll[nodeErase];
            newScore = (double)sum/ denominator;
            //cout<<"newScore = "<<newScore<<endl;
            nodeSet.erase(nodeErase);
            if(newScore > maxScore){
                max_denominator = denominator;
                maxScore = newScore;
                maxSelection = nodeSet.size();
                maxSet = nodeSet;

            }

        }

        else{
            while(!nodeSort.empty()){
                int node = nodeSort.back();
                nodeSort.pop_back();
                auto iter = nodeSet.find(node);
                if(iter!= nodeSet.end()){
                    nodeSet.erase(node);
                    updateDenominatorMT(midSource, sorceMid, midTarget, targetMid,node,denominatorAll,score,
                                        noneZeroDenominator);
                    sum = sum - score[node];
                    break;
                }
            }
        }

    }


    return maxScore;
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



double GAR(const int& numMiddle, vector<double> &score, const int& numSource,
                           unordered_map <int, unordered_set<int>> &midSource,
                           unordered_map <int, unordered_set<int>> &sorceMid, unordered_set<int> &Pickedset){


    multimap<double, int> sourceScore_id;
    unordered_map<int, double> id_sourceScore;

//    unordered_set<int> sourceSet;
//    int nonZeo = 0;
//    multimap<double, int> suScore;
//    vector<pair<int, double>> suScore;

    //first loop
    //Count the denominator |N(M)|-|N(M\u)|
    unordered_set<int> nodeSet;
    unordered_set<int> maxSet;

    for(int i = 0; i< numSource; i++){
        int source_id = i + numMiddle;
//        int source_id = i;

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
    maxSet = nodeSet;
//    int num_mid = 0;
//    int num_sor = 0;
//

    double newScore = maxScore;
//    new_score_vec.push_back(log10(newScore));
    double max_denominator = denominator;


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

        }
        if(newScore >maxScore){
            maxSet = nodeSet;
            max_numerator= sum_mid;
            max_denominator = denominator;
            maxScore = newScore;



        }
    }



    Pickedset = maxSet;

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

    unordered_set<int> nodeSet;

//    unordered_set<int> sourceSet;
//    int nonZeo = 0;
//    multimap<double, int> suScore;
//    vector<pair<int, double>> suScore;

    //first loop
    //Count the denominator |N(M)|-|N(M\u)|


    for(int i = 0; i< numSource; i++){
        int source_id = i + numMiddle;
//        int source_id = i ;

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

        maxScore = (double) max_numerator/denominator;
        maxSet = nodeSet;

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

//        L_u[min_source_node]  += min_source_score; // update the l_u according to the marginal value of min_source->second


//        log(min_source.)
//        log(sum_mid_tmp - sum_mid);
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
//            cout<<"Find better! "<<endl;



        }
//            best_iter = cnt_iter;
    }





    //    ploter.data.insert({datasetName,new_score_vec});
//    return maxScore;

}











void runLP (string filename){
//    cout<<"Now run on "<<filename<<endl;

#ifdef PRINT
    cout<<"Now run LP on "<<filename<<endl;
#endif

    int numSource = 0;
    int numMiddle = 0;
    int numTarget = 0;

    unordered_map <int, unordered_set<int>> trans;
    vector<pair<int,int>> sourceMiddle;
    vector<pair<int,int>> middleTarget;
    unordered_map <int, unordered_set<int>> middleTargetMap;

    vector<double> score = readFile(filename, trans,numSource, numMiddle, numTarget,middleTarget,middleTargetMap, sourceMiddle);






    vector<int> finalSetILP;
    double objRE = 0.0;
    auto ILP_start = std::chrono::high_resolution_clock::now();

    ILP(numMiddle,score, numSource, numTarget, trans,finalSetILP, middleTarget, sourceMiddle, objRE);

    auto ILP_end = std::chrono::high_resolution_clock::now();

    double time_ILP = std::chrono::duration_cast < std::chrono::milliseconds > (ILP_end - ILP_start).count();
#ifdef PRINT
    cout<<"Max score =  "<<objRE<<endl;
    cout<<"Rutime for LP = "<<time_ILP<<" ms"<<endl;

#endif

//    cout << "-------------------------------------------------------------------------------" << endl;


}



void runCD(string filename){
#ifdef PRINT
    cout<<"Now run CD on "<<filename<<endl;
#endif
    int merge_size = 100;
    int merge_bound = 500;


    vector<unordered_map<int, unordered_set<int>>> Multi_trans;
    vector<int> Multi_numSource;
    vector<int> Multi_numMiddle;
    vector<int> Multi_numTarget;
    vector<vector<double>> Multi_score;

    vector<vector<pair<int,int>>> Multi_sourceMiddle;
    vector<vector<pair<int,int>>> Multi_middleTarget;

    vector<unordered_map <int, unordered_set<int>>> Multi_middleTargetMap;


    double readfile_time = 0.0;


    std::tuple<vector<unordered_map<int, unordered_set<int>>>, vector<int> ,vector<int>, vector<int>,vector<vector<pair<int,int>>>,
            vector<unordered_map <int, unordered_set<int>>> ,vector<vector<pair<int,int>>>, vector<vector<double>>> result= readFile(filename, readfile_time,merge_size , merge_bound);



    std::tie(Multi_trans, Multi_numSource, Multi_numMiddle, Multi_numTarget,Multi_middleTarget,Multi_middleTargetMap, Multi_sourceMiddle,Multi_score) = result;
    double best_score = -1.0;
    int best_i = -1;

    auto ILP_start = std::chrono::high_resolution_clock::now();


    for (int i = 0; i < Multi_numMiddle.size(); ++i) {
        vector<int> finalSetILP;
        double objRE = 0.0;

        ILP(Multi_numMiddle[i], Multi_score[i], Multi_numSource[i], Multi_numTarget[i], Multi_trans[i], finalSetILP, Multi_middleTarget[i], Multi_sourceMiddle[i], objRE);

        if (objRE > best_score) {
            best_score = objRE;
            best_i = i;
        }
    }
    auto ILP_end = std::chrono::high_resolution_clock::now();

    double time_ILP = std::chrono::duration_cast < std::chrono::milliseconds > (ILP_end - ILP_start).count();


    cout<<time_ILP + readfile_time <<endl;
#ifdef PRINT
    cout << "Max score is:" << best_score << " found in " << best_i +1 << "th component"<<endl;
    cout << "Total time for CD = " << time_ILP + readfile_time << " ms" << endl;

#endif

}



void runIP (string filename){
//    cout<<"Now run on "<<filename<<endl;




    int numSource = 0;
    int numMiddle = 0;
//        int numTarget = 0;

    unordered_map <int, unordered_set<int>> midSource;
    unordered_map <int, unordered_set<int>> sorceMid;

//        unordered_map <int, unordered_set<int>> midTarget;
//        unordered_map <int, unordered_set<int>> targetMid;
    double MaxScoreSeparateV =0.0;
    int numEdges = 0;
//        vector<double> score = readFile(filename, midSource, numSource, numMiddle, numEdges, sorceMid);
//        vector<double> score = readFile_GAR_MT(filename, trans,numSource, numMiddle, numTarget,middleTarget,middleTargetMap, sourceMiddle);
    vector<double> score = readFile(filename, numSource, numMiddle, midSource, sorceMid);



    int T = 100;

    Tolerance tolerance;
    tolerance.max = 20;
    double maxScore = -1;
    double newScore;
    vector<double> L_u (numSource,0); // initialize L_u at t = 0
    unordered_set<int> maxSet;


    auto IP_start = std::chrono::high_resolution_clock::now();

    for (int t = 1; t < T; ++t) {

        IP(numMiddle,score, numSource, midSource, sorceMid, newScore, L_u, t,maxSet);
#ifdef VERBOSE
        cout<<"Score is "<<newScore<<" in iteration "<<t<<endl;
#endif
        tolerance.total++;

        if (newScore>maxScore){
            maxScore = newScore;
            tolerance.cnt= 1;
            tolerance.iteration = t;
        }else{

            tolerance.cnt++;

            if (tolerance.cnt>tolerance.max){
                break;
            }
        }
    }
    auto IP_end = std::chrono::high_resolution_clock::now();
    double time_IP = std::chrono::duration_cast < std::chrono::milliseconds > (IP_end - IP_start).count();

#ifdef PRINT
    cout << "Rutime for IP = " << time_IP << " ms" << endl;
    cout << "Score for IP = " << maxScore  << endl;
#endif

#ifdef VERBOSE
    cout<< "Max score was achieved in iteration " <<tolerance.iteration<<endl;

#endif



}




void runGAR (string filename){
#ifdef PRINT
    cout<<"Now run GAR on "<<filename<<endl;
#endif
    int numSource = 0;
    int numMiddle = 0;


    unordered_map <int, unordered_set<int>> midSource;
    unordered_map <int, unordered_set<int>> sorceMid;


    double MaxScoreSeparateV =0.0;
    int numEdges = 0;

    vector<double> score = readFile(filename, numSource, numMiddle, midSource, sorceMid, MaxScoreSeparateV);

    unordered_set<int> Pickedset;
    auto GAR_start = std::chrono::high_resolution_clock::now();

    double max_GAR = GAR(numMiddle, score, numSource, midSource, sorceMid, Pickedset);

    auto GAR_end = std::chrono::high_resolution_clock::now();

    double time_GAR = std::chrono::duration_cast<std::chrono::milliseconds>(GAR_end - GAR_start).count();

    if (MaxScoreSeparateV>max_GAR){

        max_GAR = MaxScoreSeparateV;

    }

#ifdef PRINT
    cout << "Score for GAR = " << max_GAR << endl;
    cout << "Rutime for GAR = " << time_GAR << " ms" << endl;
#endif
}




void runGR(string filename){


#ifdef PRINT
    cout<<"Now run GR on "<<filename<<endl;
#endif
    string line;
    int numSource = 0;
    int numMiddle = 0;
    int numTarget = 0;

    unordered_map <int, unordered_set<int>> midSource;
    unordered_map <int, unordered_set<int>> sorceMid;

    unordered_map <int, unordered_set<int>> midTarget;
    unordered_map <int, unordered_set<int>> targetMid;

//    int numEdges = 0;
    vector<double> score = readFile(filename, numSource, numMiddle, numTarget,  midSource, sorceMid, midTarget, targetMid);




    auto GR_start = std::chrono::high_resolution_clock::now();

    double maxscore_GR = GR_MT(numMiddle,score, numSource, midSource, sorceMid, midTarget, targetMid);
    auto GR_end = std::chrono::high_resolution_clock::now();

    double time_GR = std::chrono::duration_cast < std::chrono::milliseconds > (GR_end - GR_start).count();

#ifdef PRINT
    cout<<"MaxScore GR = "<<maxscore_GR<<endl;

    cout<<"Rutime for GR = "<<time_GR<<" ms"<<endl;

#endif

}





void runFGR(string filename){


#ifdef PRINT
    cout<<"Now run FGR on "<<filename<<endl;
#endif

    string line;
    int numSource = 0;
    int numMiddle = 0;
    int numTarget = 0;

    unordered_map <int, unordered_set<int>> midSource;
    unordered_map <int, unordered_set<int>> sorceMid;

    unordered_map <int, unordered_set<int>> midTarget;
    unordered_map <int, unordered_set<int>> targetMid;

//    int numEdges = 0;
    vector<double> score = readFile(filename, numSource, numMiddle, numTarget,  midSource, sorceMid, midTarget, targetMid);




    auto FGR_start = std::chrono::high_resolution_clock::now();

    double maxscore_FGR = FGR_MT(numMiddle,score, numSource, midSource, sorceMid, midTarget, targetMid);

    auto FGR_end = std::chrono::high_resolution_clock::now();

    double time_FGR = std::chrono::duration_cast < std::chrono::milliseconds > (FGR_end - FGR_start).count();


#ifdef PRINT
    cout<<"MaxScore FGR = "<<maxscore_FGR<<endl;

    cout<<"Rutime for FGR  = "<<time_FGR<<" ms"<<endl;

#endif


}








int main(int argc, char *argv[]){


    cmdline::parser parser;


    parser.add<string>("filepath", 'f', "the path to input graph", true, "");
    //available options of methodFlag = {LP, CD, IP GAR, GR, FGR, GRR}, default as LP
    parser.add<string>("method", 'm', "choose a method", false, "LP", cmdline::oneof<string>("LP", "CD", "IP", "GAR", "GR", "FGR"));


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

    }




    return 0;
}