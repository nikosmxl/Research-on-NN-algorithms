#include <iostream>
#include <random>
#include <bits/stdc++.h>
#include <vector>
#include "generals.h"
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
            else if (strcmp(argv[i], "LSH") == 0){
                method = "Range Search LSH";

            }
            else if (strcmp(argv[i], "Hypercube") == 0){
                method = "Range Search Hypercube";

            }
            else{
                std::cerr << "Wrong method. Enter one of the following methods next time: Classic, LSH, Hypercube.";
            }
        }
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

    int** pixels = new int*[NO_IMAGES];

    for (int i = 0 ; i < NO_IMAGES ; i++){
        pixels[i] = new int[DIMENSION];
    }

    readfile(input_file, pixels, NO_IMAGES, DIMENSION);

    Cluster** clusters = new Cluster*[number_of_clusters];
    int centers[number_of_clusters];
    for (int i = 0 ; i < number_of_clusters ; i++){
        centers[i] = -1;
    }

    // K-Means++ Initialization
    srand(time(NULL));
    kmeans_plusplus(pixels, clusters, centers, number_of_clusters, NO_IMAGES, DIMENSION, &dist);

    // Δημιουργια number_of_clusters clusters

    if (output_file.empty()){
        std::cout << "Enter output file: ";
        std::cin >> output_file;
    }

    // Δημιουργια αρχειου εξοδου
    std::ofstream Output(output_file);

    auto start = std::chrono::high_resolution_clock::now();

    // Lloyd's algorithm
    std::cout << "Processing the data..." << std::endl;

    // Assignment (Expectation)
    int* clustindexforps = new int[NO_IMAGES];
    assignment_lloyds(pixels, clusters, clustindexforps, number_of_clusters, NO_IMAGES, DIMENSION, &dist);

    // Update (Maximization)
    double changefactor = 100;
    while( changefactor >= stopfactor ){       // Οσο υπαρχει διαφορα πανω απο 5% σε σχεση με την τελευταια κατασταση συνεχιζουμε να ενημερωνουμε
        update(number_of_clusters, DIMENSION, NO_IMAGES, clustindexforps, pixels, clusters);
        int sumofchanged = 0;           // Αθροισμα διανυσματων που αλλαξαν cluster
        for(int i = 0; i < NO_IMAGES; i++){
            int min = dist(pixels[i], clusters[0]->get_center(), 2, DIMENSION);
            int mincenter = 0;
            for(int j = 1; j < number_of_clusters; j++){
                double ddd = dist(pixels[i], clusters[j]->get_center(), 2, DIMENSION);
                if( ddd < min ){
                    min = ddd;
                    mincenter = j;
                }
            }
            if( clustindexforps[i] != mincenter ){
                clusters[clustindexforps[i]]->size_down();  // Μειωνουμε το size του αφου θα χασει ενα vector
                clustindexforps[i] = mincenter;
                clusters[clustindexforps[i]]->size_up();    // Αυξανουμε το size του cluster που πηρε το vector
                sumofchanged++;
            }
        }
        changefactor = sumofchanged*100/(double)(NO_IMAGES);
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