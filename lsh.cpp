#include <iostream>
#include <fstream>
#include <random>
#include <unordered_map>
#include <bits/stdc++.h>
#include "generals.h"
#include "lsh_func.h"

int main(int argc, char const *argv[]){
    std::string input_file = "input.dat";
    std::string output_file;
    std::string query_file = "query.dat";
    std::string encoded_input_file;
    std::string encoded_query_file;
    const int NO_QUERIES = 10;
    const int NO_IMAGES = 60000;
    const int DIMENSION = 784;
    const int LATENT_DIMENSION = 10;
    const long M = 34359738363;     // 2^35 - 5 οπως λενε οι διαφανειες
    int L = 5;
    int K = 4;    
    int N = 1;
    float R = 10000.0;

    for (int i = 1 ; i < argc ; i++){
        if (strcmp(argv[i], "-d") == 0 && i + 1 < argc){
            encoded_input_file = argv[i + 1];
        }
        else if (strcmp(argv[i], "-q") == 0 && i + 1 < argc){
            encoded_query_file = argv[i + 1];
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

    if (encoded_input_file.empty()){
        std::cout << "Enter input file: ";
        std::cin >> encoded_input_file;
    }

    std::cout << "Preprocessing the data..." << std::endl;

    // Pixel array
    int** pixels = readfile<int>(input_file, NO_IMAGES, DIMENSION);

    // Encoded pixel array
    int** encoded_pixels = readfile<int>(encoded_input_file, NO_IMAGES, LATENT_DIMENSION);

    // Pixel array
    int** queries = readfile<int>(query_file, NO_QUERIES, DIMENSION);

    srand(time(NULL));

    int** w = new int*[L];
    double** t = new double*[L];
    int** rs = new int*[L]; // οι L πινακες που θα εχουν τα r για καθε map
    std::unordered_multimap<int, int>* mm[L]; // empty multimap container
    long* id = new long[NO_IMAGES];

    lsh_init(encoded_pixels, w, t, rs, mm, id, L, K, M, NO_IMAGES, LATENT_DIMENSION);

    if (encoded_query_file.empty()){
        std::cout << "Enter query file: ";
        std::cin >> encoded_query_file;
    }
    if (output_file.empty()){
        std::cout << "Enter output file: ";
        std::cin >> output_file;
    }

    // Create Output file to write
    std::ofstream Output(output_file);

    while (1){
        std::cout << "Processing the data..." << std::endl;

        // Read from query file
        int** encoded_queries = readfile<int>(encoded_query_file, NO_QUERIES, LATENT_DIMENSION);
    
        auto startLSH = std::chrono::high_resolution_clock::now();
        
        int query = rand() % NO_QUERIES;
        Neibs<int>* encoded_lsh = new Neibs<int>(encoded_pixels, encoded_queries, LATENT_DIMENSION, N, query, &dist);
        
        lsh_knn(encoded_queries, mm, encoded_lsh, w, t, rs, id, query, L, K, M, NO_IMAGES, LATENT_DIMENSION, false);
        
        Neibs<int>* encoded_rangeSearch = new Neibs<int>(encoded_pixels, encoded_queries, LATENT_DIMENSION, NO_IMAGES, query, &dist);
        
        lsh_rangeSearch(encoded_pixels, encoded_queries, mm, encoded_rangeSearch, w, t, rs, id, query, L, K, M, NO_IMAGES, LATENT_DIMENSION, R);
        
        auto stopLSH = std::chrono::high_resolution_clock::now();

        auto startReal = std::chrono::high_resolution_clock::now();

        Neibs<int>* real_neighbs = new Neibs<int>(pixels, queries, DIMENSION, N, query, &dist);
        for (int i = 0 ; i < NO_IMAGES ; i++){
            real_neighbs->insertionsort_insert(i);
        }

        auto stopReal = std::chrono::high_resolution_clock::now();

        Neibs<int>* lsh = new Neibs<int>(pixels, queries, DIMENSION, N, query, &dist);

        for (int i = 0 ; i < N ; i++){
            lsh->insertionsort_insert(encoded_lsh->givenn(i));
        }

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

        Neibs<int>* rangeSearch = new Neibs<int>(pixels, queries, DIMENSION, NO_IMAGES, query, &dist);
        int rangecount = encoded_rangeSearch->give_size();

        for (int i = 0 ; i < rangecount ; i++){
            rangeSearch->insertionsort_insert(encoded_rangeSearch->givenn(i));
        }
        
        Output << "R-near Neighbors: " << rangecount << std::endl;
        for (int i = 0 ; i < rangecount ; i++){
            Output << rangeSearch->givenn(i) << std::endl;
        }

        delete lsh;
        delete encoded_lsh;
        delete rangeSearch;
        delete encoded_rangeSearch;
        delete real_neighbs;
        for (int i = 0 ; i < NO_QUERIES ; i++){
            delete[] encoded_queries[i];
        }
        delete[] encoded_queries;
        std::cout << "Type quit to stop or a different query file name to rerun it with" << std::endl;
        std::cin >> encoded_query_file;
        if (encoded_query_file == "quit") break;
    }
    
    // Close Output file
    Output.close();

    // Deallocations
    for (int i = 0 ; i < NO_IMAGES ; i++){
        delete[] pixels[i];
        delete[] encoded_pixels[i];
    }

    for (int i = 0 ; i < NO_QUERIES ; i++){
        delete[] queries[i];
    }

    delete[] pixels;
    delete[] encoded_pixels;
    delete[] queries;

    for (int i = 0 ; i < L ; i++){
        delete[] w[i];
        delete[] t[i];
        delete[] rs[i];
        delete mm[i];
    }
    delete[] w;
    delete[] t;
    delete[] rs;
    delete[] id;

    return 0;
}