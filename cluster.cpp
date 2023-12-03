#include <iostream>
#include <random>
#include <bits/stdc++.h>
#include <vector>
#include "generals.h"
#include "lsh_func.h"
#include "cube_func.h"
#include "cluster_func.h"

int main(int argc, char const *argv[]){
    std::string input_file;
    std::string output_file;
    std::string conf_file;
    std::string method;
    const int NO_IMAGES = 60000;
    const int DIMENSION = 784;
    const long M = 34359738363;  // 2^35
    const float stopfactor = 1.0;
    int L = 3;
    int K = 4;
    int Mcube = 10;
    int dt = 14;
    int probes = 2;
    int number_of_clusters = 10;
    bool complete = false;
    float R;

    // Ελεγχος γραμμης εντολων
    for (int i = 1 ; i < argc ; i++){
        if (strcmp(argv[i], "-i") == 0 && i + 1 < argc){
            input_file = argv[i + 1];
        }
        else if (strcmp(argv[i], "-c") == 0 && i + 1 < argc){
            conf_file = argv[i + 1];
        }
        else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc){
            output_file = argv[i + 1];
        }
        else if (strcmp(argv[i], "-complete") == 0){
            complete = true;
        }
        else if (strcmp(argv[i], "-m") == 0 && i + 1 < argc){
            if (strcmp(argv[i + 1], "Classic") == 0){
                method = "Lloyds";

            }
            else if (strcmp(argv[i + 1], "LSH") == 0){
                method = "Range Search LSH";

            }
            else if (strcmp(argv[i + 1], "Hypercube") == 0){
                method = "Range Search Hypercube";

            }
            else{
                std::cerr << "Wrong method. Enter one of the following methods next time: Classic, LSH, Hypercube.";
                return -1;
            }
        }
    }

    if (method.empty()){
        std::cerr << "Error. No algorithm entered. Next time enter Classic, LSH or Hypercube";
        return -1;
    }

    // Input file (αν δεν δοθηκε στη γραμμη εντολων)
    if (input_file.empty()){
        std::cout << "Enter input file: ";
        std::cin >> input_file;
    }

    // Configuration file (αν δεν δοθηκε στη γραμμη εντολων)
    if (conf_file.empty()){
        std::cout << "Enter configuration file: ";
        std::cin >> conf_file;
    }

    std::cout << "Preprocessing the data..." << std::endl;

    if (readconf(conf_file, number_of_clusters, L, K, Mcube, dt, probes) != 0){
        return -1;
    }

    int** pixels = readfile<int>(input_file, NO_IMAGES, DIMENSION);

    int* clustindexforps = new int[NO_IMAGES];  // Σε ποιο cluster θα ανηκει το καθε διανυσμα

    Cluster<int>** clusters = new Cluster<int>*[number_of_clusters];
    int centers[number_of_clusters];
    for (int i = 0 ; i < number_of_clusters ; i++){
        centers[i] = -1;
    }

    // K-Means++ Initialization
    srand(time(NULL));
    kmeans_plusplus(pixels, clusters, centers, number_of_clusters, NO_IMAGES, DIMENSION, &dist);

    int** w = new int*[L];
    double** t = new double*[L];
    int** rs = new int*[L]; // οι L πινακες που θα εχουν τα r για καθε map
    std::unordered_multimap<int, int>* mm[L]; // empty multimap container
    long* id = new long[NO_IMAGES];
    bool* assigned = new bool[NO_IMAGES];
    std::map<int, int> hypervalues; // empty map container
    std::unordered_multimap<long, int> hypercube; // empty multimap container

    if (method == "Range Search LSH" || method == "Range Search Hypercube"){
        for (int i = 0 ; i < L ; i++){
            w[i] = new int[K];
            t[i] = new double[K];
            
            for (int j = 0 ; j < K ; j++){
                w[i][j] = rand()%5 + 2; // τυχαιο για καθε μαπ απο 2 εως 6
                t[i][j] = ( rand()%(w[i][j]*1000) )/1000.0; // τυχαιο για καθε μαπ στο [0,w)
            }
            if (method == "Range Search LSH"){
                rs[i] = new int[K];
                for (int j = 0 ; j < K ; j++){
                    rs[i][j] = rand();  // τα r ειναι τυχαια
                }
            }
            else{
                delete[] rs;
            }
        }

        if (method == "Range Search LSH"){
            for (int l = 0 ; l < L ; l++){
                mm[l] = new std::unordered_multimap<int, int>();
                for (int j = 0 ; j < NO_IMAGES ; j++){
                    bool already_center = false;
                    for (int k = 0 ; k < number_of_clusters ; k++){
                        if (j == centers[k]){
                            already_center = true;
                            break;
                        }
                    }
                    if (already_center) continue;

                    int key = g(pixels, w[l], t[l], rs[l], id, j, K, M, NO_IMAGES/8, DIMENSION, l);
                    mm[l]->insert({key, j});
                }
            }
        }
        else{
            for (int i = 0 ; i < NO_IMAGES ; i++){
                bool already_center = false;
                for (int k = 0 ; k < number_of_clusters ; k++){
                    if (i == centers[k]){
                        already_center = true;
                        break;
                    }
                }
                if (already_center) continue;

                preprocess_cube(pixels, hypervalues, hypercube, w, t, i, L, K, dt, DIMENSION);
            }
        }

        for (int i = 0 ; i < NO_IMAGES ; i++){
            assigned[i] = false;
        }
    }
    else{
        delete[] w;
        delete[] t;
        delete[] rs;
        delete[] id;
        delete[] assigned;
    }

    if (output_file.empty()){
        std::cout << "Enter output file: ";
        std::cin >> output_file;
    }

    // Δημιουργια αρχειου εξοδου
    std::ofstream Output(output_file);

    if (method == "Range Search LSH" || method == "Range Search Hypercube"){
        // Οριζουμε το R = min(αποσταση δυο cluster centers)
        R = R_init(clusters, number_of_clusters, DIMENSION, &dist);
    }

    std::cout << "Processing the data..." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();

    long query_key[number_of_clusters];
    if (method == "Range Search Hypercube"){
        for (int i = 0 ; i < number_of_clusters ; i++){
            query_key[i] = query_key_init(clusters[i]->get_center(), hypervalues, w, t, dt, K, L, DIMENSION);
        }
    }
    
    if (method == "Lloyds"){
        cluster_lloyds(pixels, clusters, clustindexforps, number_of_clusters, NO_IMAGES, DIMENSION, stopfactor, &dist);
    }
    else if (method == "Range Search LSH"){
        cluster_lsh(pixels, clusters, mm, clustindexforps, centers, assigned, w, t, rs, id, stopfactor, number_of_clusters, L, K, M, NO_IMAGES, DIMENSION, R, &dist);
    }
    else if (method == "Range Search Hypercube"){
        cluster_hypercube(pixels, clusters, hypercube, query_key, clustindexforps, assigned, number_of_clusters, stopfactor, dt, probes, Mcube, R, NO_IMAGES, DIMENSION, &dist);
    }

    if (method == "Range Search LSH" || method == "Range Search Hypercube"){
        // Κανουμε assign τα σημεια που εμειναν unassigned
        assign_unassigned(pixels, clusters, clustindexforps, assigned, number_of_clusters, NO_IMAGES, DIMENSION, &dist);
    }

    auto stop = std::chrono::high_resolution_clock::now();

    Output << "Algorithm: " << method << std::endl;
    for (int i = 0 ; i < number_of_clusters ; i++){
        Output << "CLUSTER-" << i + 1 << " {size : " << clusters[i]->get_size() << ", centroid: [" << clusters[i]->print_center_coordinates() << "] }" << std::endl;
    }

    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);

    Output << "clustering_time: " << duration.count() << " seconds" << std::endl;

    // double ss = s(0,NO_IMAGES,number_of_clusters,clustindexforps,centers,pixels,DIMENSION);

    if (complete){
        std::vector<int>** vects = new std::vector<int>*[number_of_clusters];
        for (int i = 0 ; i < number_of_clusters ; i++){
            vects[i] = new std::vector<int>;
        }
        for (int i = 0 ; i < NO_IMAGES ; i++){
            vects[clustindexforps[i]]->push_back(i);
            
        }
        for (int i = 0 ; i < number_of_clusters ; i++){
            Output << "CLUSTER-" << i + 1 << " {centroid, " << print_vector(*vects[i]) << "}" << std::endl;
        }
    }

    // Κλεισιμο αρχειου εξοδου
    Output.close();

    // Deallocations
    for (int i = 0 ; i < NO_IMAGES ; i++){
        delete[] pixels[i];
    }
    delete[] pixels;
    delete[] clustindexforps;
    delete[] clusters;

    return 0;
}

double s(int v,int no_images,int number_of_clusters,int* clustindexfops,int** centers,int** vects,int DIMENSION){

    int obj = v;
    double a = 0;
    int maincluster;
    int divider = 0;
    for(int i = 0; i < no_images; i++){
        if( clustindexfops[i] == clustindexfops[obj] && i != obj ){
            a += dist(vects[i],vects[obj],2,DIMENSION);
            maincluster = clustindexfops[obj];
            divider++;
        }
    }
    if( divider > 0 ){
        a = a/divider;
    }
    


    double mindcluster;
    int neibcluster;
    int countc = 0;
    double mm;
    for(int i = 0; i < number_of_clusters; i++){
        if( ( mm = dist(centers[i],centers[maincluster],2,DIMENSION) ) != 0 ){
            countc++;
            if( countc == 1 ){
                mindcluster = mm;
                neibcluster = i;
            }
            else{
                if( mm < mindcluster ){
                    mindcluster = mm;
                    neibcluster = i;
                }
            }
        }
    }

    double b = 0;
    int divider2 = 0;
    for(int i = 0; i < no_images; i++){
        if( clustindexfops[i] == neibcluster ){
            b += dist(vects[i],vects[obj],2,DIMENSION);
        }
    }
    if( divider > 0 ){
        b = b/divider;
    }


    if( a < b ){
        1 - a/b;
    }
    if( a > b ){
        b/a - 1;
    }
    return 0;

}