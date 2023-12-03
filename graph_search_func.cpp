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