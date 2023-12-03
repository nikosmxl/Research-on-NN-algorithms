#include "graph_search_func.h"

Graph::Graph(int vertices) : vertices(vertices) {
    adjacencyList = new std::list<int>[vertices];
}

Graph::~Graph() {
    delete[] adjacencyList;
}

void Graph::addEdge(int source, int destination) {
    adjacencyList[source].push_back(destination);
}

std::vector<int>* Graph::get_node_nn(int node, int E){
    std::vector<int>* vect = new std::vector<int>();
    int i = 0;
    for (const auto& neighbor : adjacencyList[node]) {
        vect->push_back(neighbor);
        if (i == E - 1) break;
    }

    return vect;
}

std::vector<int>* Graph::get_node_nn(int node){
    std::vector<int>* vect = new std::vector<int>();

    for (const auto& neighbor : adjacencyList[node]) {
        vect->push_back(neighbor);
    }

    return vect;
}

void graph_init_file(Graph* gr, std::string filename){
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cerr << "Error opening the file." << std::endl;
        exit(-1); // Return an error code
    }

    // Read each line from the file
    std::string line;

    int index = 0;
    int value;
    while (std::getline(inputFile, line)){
        // Use std::istringstream to extract values
        std::istringstream iss(line);

        iss >> index;

        while (iss >> value) {
            gr->addEdge(index, value);
        }
    }

    // Close the input file
    inputFile.close();
}

unsigned long getIndex(unsigned long i, unsigned long j) {
    if (i == j) {
        return -1; // Diagonal elements
    }
    if (i < j) {
        std::swap(i, j); // Ensure i >= j for upper triangular part
    }
    return i * (i - 1) * 0.5 + j; // Formula to convert 2D indices to 1D index
}

int end_init(long NUM_PS, long start, int numThreads){
    long long equality_approximation;
    unsigned long num = (NUM_PS*(NUM_PS-1)/numThreads);

    int a = start;
    int b = NUM_PS;
    while(1){
        long c = (a + b)/2;
        equality_approximation = c*c + c - num - start*start + start;
        if (abs(equality_approximation) < 0.003 * num){
            unsigned long elements = (NUM_PS - c)*(NUM_PS - 1 + c);
            if (elements < 0.997*num){
                return NUM_PS;
            }
            return c;
        }
        else if (equality_approximation < 0){
            a = c;
        }
        else{
            b = c;
        }
    }

}

double* distances_init_file(std::string filename, long NUM_PS){

    // Calculate the size of the 1D array for the upper triangular part
    unsigned long arraySize = NUM_PS * (NUM_PS - 1) * 0.5;
    double* d = new double[arraySize];

    // Open the input file
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cerr << "Error opening the file." << std::endl;
        return NULL; // Return an error code
    }

    // Read each line from the file
    std::string line;

    int index = 0;
    double value;
    while (std::getline(inputFile, line)){
        // Use std::istringstream to extract values
        std::istringstream iss(line);
        while (iss >> value) {
            // std::cout << "index is " << index << std::endl;
            d[index++] = value;
        }
    }

    // Close the input file
    inputFile.close();

    return d;
}

void set_insert(std::vector<int>& S, std::vector<int>* cands){
    for (auto it = cands->begin(); it != cands->end(); it++){
        S.push_back(*it);
    }
}

int min(int** pixels, int** queries, int query, std::vector<int>* cands, int DIM, double (*dist)(int*, int*, int, int)){
    if (cands->empty()) {
        // Handle the case where the candidates vector is empty
        std::cerr << "Error: Empty candidates vector." << std::endl;
        return -1;  // Return an error code or handle it appropriately
    }

    auto curr_it = cands->begin();
    if (curr_it == cands->end()) {
        // Handle the case where the candidates vector is empty
        std::cerr << "Error: Empty candidates vector." << std::endl;
        return -1;  // Return an error code or handle it appropriately
    }

    int min = *curr_it;

    double min_dist = dist(pixels[min], queries[query], 2, DIM);
    curr_it++;

    for (auto it = curr_it; it != cands->end(); it++) {
        double curr_dist = dist(pixels[*it], queries[query], 2, DIM);

        if (curr_dist < min_dist) {
            min = *it;
            min_dist = curr_dist;
        }
    }

    return min;
}