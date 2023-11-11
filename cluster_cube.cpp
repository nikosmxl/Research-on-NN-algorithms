#include <iostream>
#include <random>
#include <bits/stdc++.h>
#include <vector>
#include "generals.h"
#include "lsh_func.h"
#include "cluster_func.h"

int main(int argc, char const *argv[]){
    std::string input_file;
    std::string output_file;
    std::string conf_file;
    std::string method;
    const int NO_IMAGES = 60000;
    const int DIMENSION = 784;
    const long M = 34359738363;  // 2^35
    const float stopfactor = 5.0;
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

    int** w = new int*[L];
    double** t = new double*[L];

    for (int i = 0 ; i < L ; i++){
        w[i] = new int[K];
        t[i] = new double[K];
        for (int j = 0 ; j < K ; j++){
            w[i][j] = rand()%5 + 2; // τυχαιο για καθε μαπ απο 2 εως 6
            t[i][j] = ( rand()%(w[i][j]*1000) )/1000.0; // τυχαιο για καθε μαπ στο [0,w)
        }
    }

    std::map<int, int> hypervalues; // empty map container
    std::unordered_multimap<long, int> hypercube; // empty multimap container
    
    for (int i = 0 ; i < NO_IMAGES ; i++){
        bool already_center = false;
        for (int k = 0 ; k < number_of_clusters ; k++){
            if (i == centers[k]){
                already_center = true;
                break;
            }
        }
        if (already_center) continue;

        long key = 0;
        int bit = 1;
        for (int j = 0 ; j < dt ; j++){
            int zero_or_one;
            int hi = rand() % K;   // Ποια h θα επιλεξουμε
            int l = rand() % L;
            int num = h(pixels[i], w[l][hi], t[l][hi], hi, DIMENSION, l);
            auto itr = hypervalues.find(num);
            if (itr == hypervalues.end()){      // Αν δεν υπαρχει μες στο map
                zero_or_one = rand() % 2;
                hypervalues.insert({ num, zero_or_one });
            }
            else{
                zero_or_one = itr->second;
            }
            key += zero_or_one * bit;
            bit *= 2;
        }
        
        hypercube.insert({key,i});
    }

    bool* assigned = new bool[NO_IMAGES];
    for (int i = 0 ; i < NO_IMAGES ; i++){
        assigned[i] = false;
    }

    // Δημιουργια number_of_clusters clusters

    if (output_file.empty()){
        std::cout << "Enter output file: ";
        std::cin >> output_file;
    }

    // Δημιουργια αρχειου εξοδου
    std::ofstream Output(output_file);

    
    // Οριζουμε το R = min(αποσταση δυο cluster centers)
    float R = R_init(clusters, number_of_clusters, DIMENSION, &dist);
    
    std::cout << "Processing the data..." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();

    long query_key[number_of_clusters];
    for (int i = 0 ; i < number_of_clusters ; i++){
        query_key[i] = 0;
        int bit = 1;
        for (int j = 0 ; j < dt ; j++){
            int zero_or_one;
            int hi = rand() % K;   // Ποια h θα επιλεξουμε
            int l = rand() % L;
            int num = h(clusters[i]->get_center(), w[l][hi], t[l][hi], hi, DIMENSION, l); 
            auto itr = hypervalues.find(num);
            if (itr != hypervalues.end()){      // Αν δεν υπαρχει μες στο map
                zero_or_one = rand() % 2;
                hypervalues.insert({ num, zero_or_one });
            }
            else{
                zero_or_one = itr->second;
            }
            query_key[i] += zero_or_one * bit;
            bit *= 2;
        }
    }

    // Assignment by Hypercube Range Search
    int* clustindexforps = new int[NO_IMAGES];

    double changefactor = 100;
    while (changefactor <= stopfactor){        // Οσο υπαρχει διαφορα πανω απο 5% σε σχεση με την τελευταια κατασταση συνεχιζουμε να ενημερωνουμε
        int sumofchanged = 0;           // Αθροισμα διανυσματων που αλλαξαν cluster
        for (int i = 0 ; i < number_of_clusters ; i++){
            auto itr = hypercube.equal_range(query_key[i]);

            std::bitset<32> binary_query_key(query_key[i]); // Convert p to binary representation
            int countVerticesRange = 0;
            int countCubeElementsRange = 0;
            int hamming_dist_wantedRange = 0;

            int hd = 0;

            while (countCubeElementsRange < Mcube && countVerticesRange < probes){
                for (auto it = itr.first; it != itr.second; it++) {
                    countCubeElementsRange++;
                    int d = dist(pixels[it->second],clusters[i]->get_center(),2,DIMENSION);
                    if(d >= R) continue;
                    
                    if (assigned[it->second]){
                        if (d < dist(pixels[it->second], clusters[clustindexforps[it->second]]->get_center(), 2, DIMENSION)){
                            clusters[clustindexforps[it->second]]->size_down();
                            clustindexforps[it->second] = i;
                            clusters[clustindexforps[it->second]]->size_up();
                            sumofchanged++;
                        }
                    }
                    else{
                        clustindexforps[it->second] = i;
                        clusters[clustindexforps[it->second]]->size_up();
                        assigned[it->second] = true;
                        sumofchanged++;
                    }
                }
                countVerticesRange++;

                if (hd == query_key[i]){
                    hd++;
                    hd = hd % dt;
                }
                if (hd == 0) hamming_dist_wantedRange++;
                while (1){
                    std::bitset<32> binary_hd(hd); // Convert num to binary representation
                    int differingBits = (binary_query_key ^ binary_hd).count(); // Count differing bits
                    if (differingBits == hamming_dist_wantedRange){
                        itr = hypercube.equal_range(hd);
                        break;
                    }

                    hd++;
                    hd = hd % dt;
                    if (hd == 0) hamming_dist_wantedRange++;
                }
                
            }
        }
        
        R *= 2;
        changefactor = sumofchanged*100/(double)(NO_IMAGES);
    }

    // Κανουμε assign τα σημεια που εμειναν unassigned
    for (int i = 0 ; i < NO_IMAGES ; i++){
        if (assigned[i]) continue;

        int min = dist(pixels[i], clusters[0]->get_center(), 2, DIMENSION);
        int mincenter = 0;
        for(int j = 1; j < number_of_clusters; j++){
            double ddd = dist(pixels[i], clusters[j]->get_center(), 2, DIMENSION);
            if( ddd < min ){
                min = ddd;
                mincenter = j;
            }
        }

        clustindexforps[i] = mincenter;
        clusters[clustindexforps[i]]->size_up();
    }

    auto stop = std::chrono::high_resolution_clock::now();

    // Prints
    Output << "Algorithm: " << method << std::endl;
    for (int i = 0 ; i < number_of_clusters ; i++){
        Output << "CLUSTER-" << i + 1 << " {size : " << clusters[i]->get_size() << ", centroid: [" << clusters[i]->print_center_coordinates() << "] }" << std::endl;
    }

    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);

    Output << "clustering_time: " << duration.count() << " seconds" << std::endl;

    //floyd kmeans

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