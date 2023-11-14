#include <iostream>
#include <fstream>
#include <random>
#include <unordered_map>
#include <bits/stdc++.h>
#include "generals.h"
#include "lsh_func.h"
#include "cube_func.h"

int main(int argc, char const *argv[]){
    std::string input_file;
    std::string output_file;
    std::string query_file;
    const int NO_QUERIES = 10;
    const int NO_IMAGES = 60000;
    const int DIMENSION = 784;
    const int K = 4;
    const int L = 5;
    int M = 10;
    int N = 1;
    int dt = 14;    // d'
    int probes = 2;
    float R = 10000.0;

    for (int i = 1 ; i < argc ; i++){
        if (strcmp(argv[i], "-d") == 0 && i + 1 < argc){
            input_file = argv[i + 1];
        }
        else if (strcmp(argv[i], "-q") == 0 && i + 1 < argc){
            query_file = argv[i + 1];
        }
        else if (strcmp(argv[i], "-k") == 0 && i + 1 < argc){
            dt = atoi(argv[i + 1]);
        }
        else if (strcmp(argv[i], "-M") == 0 && i + 1 < argc){
            M = atoi(argv[i + 1]);
        }
        else if (strcmp(argv[i], "-probes") == 0 && i + 1 < argc){
            probes = atoi(argv[i + 1]);
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
    }

    std::cout << "Preprocessing the data..." << std::endl;

    // Pixel array
    int** pixels = readfile<int>(input_file, NO_IMAGES, DIMENSION);

    srand(time(NULL));

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

    if (query_file.empty()){
        std::cout << "Enter query file: ";
        std::cin >> query_file;
    }
    if (output_file.empty()){
        std::cout << "Enter output file: ";
        std::cin >> output_file;
    }

    std::map<int, int> hypervalues; // empty map container
    std::unordered_multimap<long, int> hypercube; // empty multimap container

    for (int i = 0 ; i < NO_IMAGES ; i++){
        preprocess_cube(pixels, hypervalues, hypercube, w, t, i, L, K, dt, DIMENSION);
    }

    // Create Output file to write
    std::ofstream Output(output_file);

    while (1){
        std::cout << "Processing the data..." << std::endl;

        // Read from query file
        int** queries = readfile<int>(query_file, NO_QUERIES, DIMENSION);

        auto startCube = std::chrono::high_resolution_clock::now();
        
        int query = rand() % NO_QUERIES;

        long query_key = query_key_init(queries[query], hypervalues, w, t, dt, K, L, DIMENSION);
        
        // Hypercube knn
        Neibs<int>* cube_neibs = new Neibs<int>(pixels, queries, DIMENSION, N, query, &dist);

        cube_knn(cube_neibs, hypervalues, hypercube, query_key, dt, M, probes);

        // Hypercube Range Search
        Neibs<int>* rangeSearchCube = new Neibs<int>(pixels, queries, DIMENSION, NO_IMAGES, query, &dist);
        
        cube_rangeSearch(rangeSearchCube, hypervalues, hypercube, query_key, dt, M, probes, pixels, queries, query, DIMENSION, R);

        auto stopCube = std::chrono::high_resolution_clock::now();

        auto startReal = std::chrono::high_resolution_clock::now();

        Neibs<int>* real_neighbs = new Neibs<int>(pixels, queries, DIMENSION, N, query, &dist);
        for (int i = 0 ; i < NO_IMAGES ; i++){
            real_neighbs->insertionsort_insert(i);
        }

        auto stopReal = std::chrono::high_resolution_clock::now();

        Output << "Query: " << query << std::endl;

        for (int i = 0 ; i < N ; i++){
            Output << "HyperCube Nearest neighbor-" << i + 1 << ": " << cube_neibs->givenn(i) << std::endl;
            Output << "distance HyperCube: " << cube_neibs->givedist(i) << std::endl;
            Output << "distanceTrue: " << real_neighbs->givedist(i) << std::endl;
        }

        auto durationCube = std::chrono::duration_cast<std::chrono::milliseconds>(stopCube - startCube);
        auto durationReal = std::chrono::duration_cast<std::chrono::milliseconds>(stopReal - startReal);

        Output << "tCube: " << durationCube.count() << " milliseconds" << std::endl;
        Output << "tTrue: " << durationReal.count() << " milliseconds" << std::endl;

        int rangecount = rangeSearchCube->give_size();
        Output << "R-near neighbors: " << rangecount << std::endl;
        for (int i = 0 ; i < rangecount ; i++){
            Output << rangeSearchCube->givenn(i) << std::endl;
        }

        delete cube_neibs;
        delete rangeSearchCube;
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
        // delete mm[i];
    }
    delete[] w;
    delete[] t;

    return 0;
}
