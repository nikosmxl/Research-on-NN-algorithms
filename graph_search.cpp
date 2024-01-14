#include <iostream>
#include <string.h>
#include "graph_search_func.h"

int main(int argc, char const *argv[]) {
    std::string input_file = "input.dat";
    std::string output_file;
    std::string query_file = "query.dat";
    std::string encoded_input_file;
    std::string encoded_query_file;
    std::string method;
    const int NO_QUERIES = 10;
    const int NO_IMAGES = 60000;
    const int DIMENSION = 784;
    const int LATENT_DIMENSION = 10;
    const long M = 34359738363;         // 2^35 - 5 οπως λενε οι διαφανειες
    const int T = 30;                   // Οριο αριθμου βηματων
    unsigned int N = 1;                 // Πλησιεστεροι γειτονες
    unsigned int k = 50;                // Πλησιεστεροι γειτονες στον γραφο k-NN
    unsigned int R = 1;                 // Αριθμος τυχαιων επανεκκινησεων
    unsigned int E = 30;                // Αριθμος επεκτασεων
    std::size_t candidate_size = 20;    // Το πληθος υποψηφιων L που θα δεχτουμε
    unsigned int L = ceil(k *0.125);    // L maps της LSH
    unsigned int K = floor(L *0.625);   // K παραμετρος της LSH
    bool local_optimal = false;

    for (int i = 1 ; i < argc ; i++){
        if (strcmp(argv[i], "-d") == 0 && i + 1 < argc){
            encoded_input_file = argv[i + 1];
        }
        else if (strcmp(argv[i], "-q") == 0 && i + 1 < argc){
            encoded_query_file = argv[i + 1];
        }
        else if (strcmp(argv[i], "-k") == 0 && i + 1 < argc){
            k = atoi(argv[i + 1]);
        }
        else if (strcmp(argv[i], "-E") == 0 && i + 1 < argc){
            E = atoi(argv[i + 1]);
        }
        else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc){
            output_file = argv[i + 1];
        }
        else if (strcmp(argv[i], "-N") == 0 && i + 1 < argc){
            N = atoi(argv[i + 1]);
        }
        else if (strcmp(argv[i], "-R") == 0 && i + 1 < argc){
            R = atoi(argv[i + 1]);
        }
        else if (strcmp(argv[i], "-l") == 0 && i + 1 < argc){
            candidate_size = atoi(argv[i + 1]);
        }
        else if (strcmp(argv[i], "-m") == 0 && i + 1 < argc){
            if (strcmp(argv[i + 1], "1") == 0){
                method = "GNNS Results";

            }
            else if (strcmp(argv[i + 1], "2") == 0){
                method = "MRNG Results";

            }
            else{
                std::cerr << "Wrong method. Enter one of the following methods next time: 1 for GNNS, 2 for MRNG ";
                return -1;
            }
        }
    }

    // If method is empty exit
    if (method.empty()){
        std::cerr << "Error. No method entered. Next time enter 1 for GNNS, 2 for MRNG";
        return -1;
    }

    // N <= l
    if (method == "MRNG Results" && N > candidate_size){
        std::cerr << "Error. N cannot be bigger than L. Try again with valid arguments";
        return -1;
    }

    // If input file is not given
    if (encoded_input_file.empty()){
        std::cout << "Enter input file: ";
        std::getline(std::cin, encoded_input_file);
        if (std::cin.fail() || encoded_input_file.empty()) exit(-1);
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

    std::cout << "Creating the graph..." << std::endl;
    
    Graph* gr = new Graph(NO_IMAGES);
    if (method == "GNNS Results"){
        gnns_construction(encoded_pixels, gr, mm, w, t, rs, id, k, L, K, M, NO_IMAGES, LATENT_DIMENSION, &dist);
    }
    else if (method == "MRNG Results"){
        // threaded_mrng(4, gr, encoded_pixels, mm, w, t, rs, id, k, L, K, M, NO_IMAGES, LATENT_DIMENSION, &dist);
        mrng_optimal(gr, encoded_pixels, mm, w, t, rs, id, k, L, K, M, NO_IMAGES, LATENT_DIMENSION, &dist);
    }
    std::cout << "EDW " << std::endl;
    // Calculating navigation node in case it is MRNG method
    Neibs<int>* navigation_node;
    if (method == "MRNG Results"){
        int** centroid = new int*[1];
        centroid[0] = new int[LATENT_DIMENSION];
        for (int i = 0 ; i < LATENT_DIMENSION ; i++){
            centroid[0][i] = 0;
            for (int j = 0 ; j < NO_IMAGES ; j++){
                centroid[0][i] += encoded_pixels[i][j];
            }
            int avg = round(centroid[0][i] / (double)NO_IMAGES);
            centroid[0][i] = avg;
        }
        
        navigation_node = new Neibs<int>(encoded_pixels, centroid, LATENT_DIMENSION, 1, 0, &dist);
        lsh_knn(centroid, mm, navigation_node, w, t, rs, id, 0, L, K, M, NO_IMAGES, LATENT_DIMENSION, false);
    }

    // If query file is not given
    if (encoded_query_file.empty()){
        std::cout << "Enter query file: ";
        std::getline(std::cin, encoded_query_file);
        if (std::cin.fail() || encoded_query_file.empty()) exit(-1);
    }

    // If output file is not given
    if (output_file.empty()){
        std::cout << "Enter output file: ";
        std::getline(std::cin, output_file);
        if (std::cin.fail() || output_file.empty()) exit(-1);
    }

    // Create Output file to write
    std::ofstream Output(output_file);

    std::vector<double> timesSearch;
    std::vector<double> timesTrue;
    double AAF = 0;

    Output << method << std::endl;

    while (1){
        std::cout << "Processing the data..." << std::endl;
        
        // Read from query file
        int** encoded_queries = readfile<int>(encoded_query_file, NO_QUERIES, LATENT_DIMENSION);

        auto startSearch = std::chrono::high_resolution_clock::now();

        int query = rand() % NO_QUERIES;        // Random query

        // Process
        Neibs<int>* EncodedSearch;
        if (method == "GNNS Results"){
            EncodedSearch = gnns_search(gr, encoded_pixels, encoded_queries, query, N, T, R, E, NO_IMAGES, LATENT_DIMENSION, local_optimal, &dist);
        }
        else if (method == "MRNG Results"){
            EncodedSearch = search_on_graph(gr, encoded_pixels, encoded_queries, query, navigation_node, candidate_size, N, NO_IMAGES, LATENT_DIMENSION, &dist);
        }

        auto stopSearch = std::chrono::high_resolution_clock::now();

        auto startReal = std::chrono::high_resolution_clock::now();

        Neibs<int>* real_neighbs = new Neibs<int>(pixels, queries, DIMENSION, N, query, &dist);
        for (int i = 0 ; i < NO_IMAGES ; i++){
            real_neighbs->insertionsort_insert(i);
        }

        auto stopReal = std::chrono::high_resolution_clock::now();

        Neibs<int>* Search = new Neibs<int>(pixels, queries, DIMENSION, N, query, dist);
        for (int i = 0 ; i < N ; i++){
            Search->insertionsort_insert(EncodedSearch->givenn(i));
        }

        // Prints
        Output << "Query: " << query << std::endl;
        for (int i = 0 ; i < N ; i++){
            Output << "Nearest neighbor-" << i + 1 << ": " << Search->givenn(i) << std::endl;
            Output << "distanceApproximate: " << Search->givedist(i) << std::endl;
            Output << "distanceTrue: " << real_neighbs->givedist(i) << std::endl;
            
            if (i == 0){
                AAF += Search->givedist(i) / real_neighbs->givedist(i);
            }
        }

        // Calculating total time taken by the program.
        double durationSearch = std::chrono::duration_cast<std::chrono::milliseconds>(stopSearch - startSearch).count();
        durationSearch *= 1e-3;
        timesSearch.push_back(durationSearch);

        double durationReal = std::chrono::duration_cast<std::chrono::milliseconds>(stopReal - startReal).count();
        durationReal *= 1e-3;
        timesTrue.push_back(durationReal);

        // Deallocations
        for (int i = 0 ; i < NO_QUERIES ; i++){
            delete[] encoded_queries[i];
        }
        delete[] encoded_queries;
        delete real_neighbs;
        delete EncodedSearch;
        delete Search;

        // Quit or rerun with a different/same file
        std::cout << "Type quit to stop, type a different query file name to rerun it with or press enter to rerun with the same query file" << std::endl;
        std::string input;
        std::getline(std::cin, input);
        if (std::cin.fail() || input.empty()) continue;
        if (input == "quit") break;
        encoded_query_file = input;
    }

    // Calculating average time of the algorithm and average real time
    double averageSearch = 0.0;
    double averageTrue = 0.0;
    for (auto it = timesSearch.begin(); it != timesSearch.end(); it++){
        averageSearch += *it;
    }
    averageSearch /= timesSearch.size();

    for (auto it = timesTrue.begin(); it != timesTrue.end(); it++){
        averageTrue += *it;
    }
    averageTrue /= timesTrue.size();

    AAF /= timesTrue.size();
    
    // Printing the average time of the algorithm and the average real time
    Output << "tAverageApproximate: " << std::fixed << std::setprecision(4) << averageSearch << " seconds" << std::endl;
    Output << "tAverageTrue: " << std::fixed << std::setprecision(4) << averageTrue << " seconds" << std::endl;
    Output << "AAF: " << std::fixed << std::setprecision(4) <<  AAF << std::endl;
    
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
    delete gr;

    return 0;
}
