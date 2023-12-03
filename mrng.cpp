#include <iostream>
#include <string.h>
#include "graph_search_func.h"

int main(int argc, char const *argv[]) {
    std::string input_file;
    std::string output_file;
    std::string query_file;
    std::string method;
    const int NO_QUERIES = 10;
    const int NO_IMAGES = 60000;
    const int DIMENSION = 784;
    const long M = 34359738363;         // 2^35 - 5 οπως λενε οι διαφανειες
    const int T = 30;                   // Οριο αριθμου βηματων
    unsigned int N = 1;                 // Πλησιεστεροι γειτονες
    unsigned int k = 50;                // Πλησιεστεροι γειτονες στον γραφο k-NN
    unsigned int R = 1;                 // Αριθμος τυχαιων επανεκκινησεων
    unsigned int E = 30;                // Αριθμος επεκτασεων
    std::size_t candidate_size = 20;    // Το πληθος υποψηφιων L που θα δεχτουμε
    unsigned int L = ceil(k / 8.0);     // L maps της LSH
    unsigned int K = floor(L *0.625);   // K παραμετρος της LSH

    for (int i = 1 ; i < argc ; i++){
        if (strcmp(argv[i], "-d") == 0 && i + 1 < argc){
            input_file = argv[i + 1];
        }
        else if (strcmp(argv[i], "-q") == 0 && i + 1 < argc){
            query_file = argv[i + 1];
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

    if (method.empty()){
        std::cerr << "Error. No method entered. Next time enter 1 for GNNS, 2 for MRNG";
        return -1;
    }

    if (method == "MRNG Results" && N > candidate_size){
        std::cerr << "Error. N cannot be bigger than L. Try again with valid arguments";
        return -1;
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
    int** rs = new int*[L]; // οι L πινακες που θα εχουν τα r για καθε map
    std::unordered_multimap<int, int>* mm[L]; // empty multimap container
    long* id = new long[NO_IMAGES];

    lsh_init(pixels, w, t, rs, mm, id, L, K, M, NO_IMAGES, DIMENSION);

    std::cout << "Creating the graph..." << std::endl;

    Graph* gr = new Graph(NO_IMAGES);

    std::ofstream Apotelesmata("apotelesmata.txt");

    // auto startdist = std::chrono::high_resolution_clock::now();

    // double* distances = distances_init("apostaseis.txt", NO_IMAGES);

    // auto stopdist = std::chrono::high_resolution_clock::now();

    // double durationdistances = std::chrono::duration_cast<std::chrono::minutes>(stopdist - startdist).count();
    // Apotelesmata << "WHOLE TIME FOR DISTANCES IS : " << durationdistances << " minutes" << std::endl;

    auto startmrng = std::chrono::high_resolution_clock::now();

    threaded_mrng(4, gr, pixels, mm, w, t, rs, id, k, L, K, M, NO_IMAGES, DIMENSION, &dist);
    // graph_init(gr, "graph2.txt");
    
    auto stopmrng = std::chrono::high_resolution_clock::now();

    double durationMRNG = std::chrono::duration_cast<std::chrono::minutes>(stopmrng - startmrng).count();
    Apotelesmata << "WHOLE TIME FOR MRNG WITH 4 THREADS IS : " << durationMRNG << " minutes" << std::endl;

    std::cout << "Made the graph..." << std::endl;
    std::ofstream Graphfile("graph.txt");
    for (int i = 0 ; i < NO_IMAGES ; i++){
        Graphfile << i << " ";
        std::vector<int>* vert = gr->get_node_nn(i);
        for (auto j = vert->begin() ; j != vert->end() ; j++){
            Graphfile << *j << " ";
        }
        Graphfile << std::endl;
    }
    std::cout << "Creating the centroid..." << std::endl;
    int** centroid = new int*[1];
    centroid[0] = new int[DIMENSION];
    for (int i = 0 ; i < DIMENSION ; i++){
        centroid[0][i] = 0;
        for (int j = 0 ; j < NO_IMAGES ; j++){
            centroid[0][i] += pixels[i][j];
        }
    }
    std::cout << "Creating the navigation node..." << std::endl;
    Neibs<int>* navigation_node = new Neibs<int>(pixels, centroid, DIMENSION, 1, 0, &dist);
    lsh_knn(centroid, mm, navigation_node, w, t, rs, id, 0, L, K, M, NO_IMAGES, DIMENSION, false);

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

    std::vector<double> timesSOG;
    std::vector<double> timesTrue;
    double MAF = -1;

    std::string del;
    std::cout << "ta duskola perasan..mallon?" << std::endl;
    std::cin >> del;

    Output << method << std::endl;

    while (1){
        std::cout << "Processing the data..." << std::endl;
        
        // Read from query file
        int** queries = readfile<int>(query_file, NO_QUERIES, DIMENSION);

        auto startSOG = std::chrono::high_resolution_clock::now();

        int query = rand() % NO_QUERIES;

        // Process
        bool checked[NO_IMAGES];
        for (int i = 0 ; i < NO_IMAGES ; i++){
            checked[i] = false;
        }
        std::set<int, DistComparator<int>> Rset(DistComparator<int>(pixels, queries, query, DIMENSION));
        Rset.insert(navigation_node->givenn(0));
        std::cout << "nav node is " << navigation_node->givenn(0) << std::endl;
        int cands_checked = 0;
        while (cands_checked < candidate_size){
            std::cout << "candschecked is " << cands_checked << std::endl;
            for (auto p : Rset){
                std::cout << "p is " << p << " and checked of p is " << checked[p] << std::endl;
                if (!checked[p]){std::cout << "MPAINW GENIKA" << std::endl;
                    checked[p] = true;
                    cands_checked++;
                    std::vector<int>* N = gr->get_node_nn(p);
                    for (auto it = N->begin() ; it < N->end() ; it++){
                        std::cout << "EXEI " << *it << std::endl;
                        Rset.insert(*it);
                    }

                    delete N;
                    break;
                }
            }
            
        }

        Neibs<int>* SOG = new Neibs<int>(pixels, queries, DIMENSION, N, query, &dist);
        int first_k = 0;
        for (auto r : Rset){
            SOG->insertionsort_insert(r);
            first_k++;
            if (first_k == N) break;
        }

        std::cout << "Did the search..." << std::endl;
        auto stopSOG = std::chrono::high_resolution_clock::now();

        auto startReal = std::chrono::high_resolution_clock::now();

        Neibs<int>* real_neighbs = new Neibs<int>(pixels, queries, DIMENSION, N, query, &dist);
        for (int i = 0 ; i < NO_IMAGES ; i++){
            real_neighbs->insertionsort_insert(i);
        }

        auto stopReal = std::chrono::high_resolution_clock::now();

        Output << "Query: " << query << std::endl;
        for (int i = 0 ; i < N ; i++){
            Output << "Nearest neighbor-" << i + 1 << ": " << SOG->givenn(i) << std::endl;
            Output << "distanceApproximate: " << SOG->givedist(i) << std::endl;
            Output << "distanceTrue: " << real_neighbs->givedist(i) << std::endl;
            if (i == 0){
                double AF = SOG->givedist(i) / real_neighbs->givedist(i);
                if (MAF != -1){
                    if (MAF < AF){
                        MAF = AF;
                    }
                }
                else{
                    MAF = AF;
                }
            }
        }

        // Calculating total time taken by the program.
        double durationSOG = std::chrono::duration_cast<std::chrono::milliseconds>(stopSOG - startSOG).count();
        durationSOG *= 1e-3;
        timesSOG.push_back(durationSOG);

        double durationReal = std::chrono::duration_cast<std::chrono::milliseconds>(stopReal - startReal).count();
        durationReal *= 1e-3;
        timesTrue.push_back(durationReal);

        // Deallocations
        for (int i = 0 ; i < NO_QUERIES ; i++){
            delete[] queries[i];
        }
        delete[] queries;
        delete real_neighbs;
        delete SOG;

        // Quit or rerun with a different file
        std::cout << "Type quit to stop, type a different query file name to rerun it with or press enter to rerun with the same query file" << std::endl;
        std::string input;
        std::getline(std::cin, input);
        if (std::cin.fail() || input.empty()) continue;
        if (input == "quit") break;
        query_file = input;
    }

    // Calculating average time of the algorithm and average real time
    double averageSOG = 0.0;
    double averageTrue = 0.0;
    for (auto it = timesSOG.begin(); it != timesSOG.end(); it++){
        averageSOG += *it;
    }
    averageSOG /= timesSOG.size();

    for (auto it = timesTrue.begin(); it != timesTrue.end(); it++){
        averageTrue += *it;
    }
    averageTrue /= timesTrue.size();
    
    // Printing the average time of the algorithm and the average real time
    Output << "tAverageApproximate: " << std::fixed << std::setprecision(4) << averageSOG << " seconds" << std::endl;
    Output << "tAverageTrue: " << std::fixed << std::setprecision(4) << averageTrue << " seconds" << std::endl;
    Output << "MAF: " << std::fixed << std::setprecision(4) <<  MAF << std::endl;
    
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
        delete mm[i];
    }
    delete[] w;
    delete[] t;
    delete[] rs;
    delete[] id;
    delete gr;

    return 0;
}
