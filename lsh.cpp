#include <iostream>
#include <fstream>
#include <random>
#include <unordered_map>
#include <bits/stdc++.h>
#include "generals.h"
#include "lsh_func.h"

int main(int argc, char const *argv[]){
    std::string input_file;
    std::string output_file;
    std::string query_file;
    const int NO_QUERIES = 10;
    const int NO_IMAGES = 60000;
    const int DIMENSION = 784;
    const long M = 34359738363;     // 2^35 - 5 οπως λενε οι διαφανειες
    int L = 5;
    int K = 4;    
    int N = 1;
    float R = 10000.0;

    for (int i = 1 ; i < argc ; i++){
        if (strcmp(argv[i], "-d") == 0 && i + 1 < argc){
            input_file = argv[i + 1];
        }
        else if (strcmp(argv[i], "-q") == 0 && i + 1 < argc){
            query_file = argv[i + 1];
        }
        else if (strcmp(argv[i], "-k") == 0 && i + 1 < argc){
            K = atoi(argv[i + 1]);
        }
        else if (strcmp(argv[i], "-L") == 0 && i + 1 < argc){
            L = atoi(argv[i + 1]);
        }
        else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc){
            output_file = argv[i + 1];
        }
        else if (strcmp(argv[i], "-N") == 0 && i + 1 < argc){
            N = atoi(argv[i + 1]);
        }
        else if (strcmp(argv[i], "-R") == 0 && i + 1 < argc){
            R = atof(argv[i + 1]);
        }
    }

    if (input_file.empty()){
        std::cout << "Enter input file: ";
        std::cin >> input_file;
        std::cout << "Preprocessing the data..." << std::endl;
    }

    // Pixel array
    int** pixels = new int*[NO_IMAGES];
    int** queries = new int*[NO_QUERIES];

    for (int i = 0 ; i < NO_IMAGES ; i++){
        pixels[i] = new int[DIMENSION];
    }

    for (int i = 0 ; i < NO_QUERIES ; i++){
        queries[i] = new int[DIMENSION];
    }

    readfile(input_file, pixels, NO_IMAGES, DIMENSION);

    srand(time(NULL));

    int** w = new int*[L];
    double** t = new double*[L];
    int** rs = new int*[L]; // οι L πινακες που θα εχουν τα r για καθε map

    for (int i = 0 ; i < L ; i++){
        w[i] = new int[K];
        t[i] = new double[K];
        rs[i] = new int[K];
        for (int j = 0 ; j < K ; j++){
            w[i][j] = rand()%5 + 2; // τυχαιο για καθε μαπ απο 2 εως 6
            t[i][j] = ( rand()%(w[i][j]*1000) )/1000.0; // τυχαιο για καθε μαπ στο [0,w)
            rs[i][j] = rand();  // τα r ειναι τυχαια
        }
    }

    std::unordered_multimap<int, int>* mm[L]; // empty multimap container
    long* id = new long[NO_IMAGES];

    for (int l = 0 ; l < L ; l++){
        mm[l] = new std::unordered_multimap<int, int>();
        for (int j = 0 ; j < NO_IMAGES ; j++){
            int key = g(pixels, w[l], t[l], rs[l], id, j, K, M, NO_IMAGES/8, DIMENSION, l);
            mm[l]->insert({key, j});
        }
    }

    if (query_file.empty()){
        std::cout << "Enter query file: ";
        std::cin >> query_file;
    }
    if (output_file.empty()){
        std::cout << "Enter output file: ";
        std::cin >> output_file;
    }

    // Create Output file to write
    std::ofstream Output(output_file);

    while (1){
        // Read from query file
        readfile(query_file, queries, NO_QUERIES, DIMENSION);

        auto startLSH = std::chrono::high_resolution_clock::now();
        
        int query = 0;
        Neibs* lsh = new Neibs(pixels, queries, DIMENSION, N, query, &dist);
        
        int countLSH = 0;
        for(int l = 0; l < L; l++){
            int key = g(pixels, w[l], t[l], rs[l], id, query, K, M, NO_IMAGES/8, DIMENSION, l);
            auto itr = mm[l]->equal_range(key);
            for (auto it = itr.first; it != itr.second; it++) {
                lsh->insertionsort_insert(it->second);
                countLSH++;
                if( countLSH > 10*L ){
                    break;
                }
            }
            
            if( countLSH > 10*L ){
                    break;
            }
        }

        Neibs* rangeSearch = new Neibs(pixels, queries, DIMENSION, NO_IMAGES, query, &dist);
        
        int countRange = 0;
        for(int l = 0; l < L; l++){
            int key = g(pixels, w[l], t[l], rs[l], id, query, K, M, NO_IMAGES/8, DIMENSION, l);
            auto itr = mm[l]->equal_range(key);
            for (auto it = itr.first; it != itr.second; it++) {
                countRange++;
                if ( dist(pixels[it->second],queries[query],2,DIMENSION) >= R ) continue;
                    rangeSearch->insertionsort_insert(it->second);
                if( countRange > 20*L ){
                    break;
                }
            }
            
            if( countRange > 20*L ){
                break;
            }
        }

        auto stopLSH = std::chrono::high_resolution_clock::now();

        auto startReal = std::chrono::high_resolution_clock::now();

        Neibs* real_neighbs = new Neibs(pixels, queries, DIMENSION, N, query, &dist);
        for (int i = 0 ; i < NO_IMAGES ; i++){
            real_neighbs->insertionsort_insert(i);
        }

        auto stopReal = std::chrono::high_resolution_clock::now();

        Output << "Query: " << query << std::endl;
        for (int i = 0 ; i < N ; i++){
            Output << "Nearest neighbor-" << i + 1 << ": " << lsh->givenn(i) << std::endl;
            Output << "distanceLSH: " << lsh->givedist(i) << std::endl;
            Output << "distanceTrue: " << real_neighbs->givedist(i) << std::endl;
        }

        auto durationLSH = std::chrono::duration_cast<std::chrono::milliseconds>(stopLSH - startLSH);
        auto durationReal = std::chrono::duration_cast<std::chrono::milliseconds>(stopReal - startReal);

        Output << "tLSH: " << durationLSH.count() << " milliseconds" << std::endl;
        Output << "tTrue: " << durationReal.count() << " milliseconds" << std::endl;

        int rangecount = rangeSearch->give_size();
        Output << "R-near Neighbors: " << rangecount << std::endl;
        for (int i = 0 ; i < rangecount ; i++){
            Output << rangeSearch->givenn(i) << std::endl;
        }

        delete lsh;
        delete rangeSearch;
        delete real_neighbs;

        std::cout << "Type quit to stop or a different query file name to rerun it with" << std::endl;
        std::cin >> query_file;
        if (query_file == "quit") break;
    }
    
    // Close Output file
    Output.close();

    // Deallocations
    for (int i = 0 ; i < NO_IMAGES ; i++){
        delete[] pixels[i];
    }
    delete[] pixels;

    for (int i = 0 ; i < L ; i++){
        delete[] w[i];
        delete[] t[i];
        delete[] rs[i];
        // delete mm[i];
    }
    delete[] w;
    delete[] t;
    delete[] rs;

    delete[] id;

    return 0;
}