#include <set>
#include <list>
#include <vector>
#include <sstream>
#include <thread>
#include "generals.h"
#include "lsh_func.h"

unsigned long getIndex(unsigned long i, unsigned long j);

int end_init(long NUM_PS, long start, int numThreads);

double* distances_init_file(std::string filename, long NUM_PS);

template<typename T>
void distances_init(T** p, double* d, long start, long stop, int DIMENSION, double (*dist)(T*, T*, int, int)){
    for (unsigned long i = start ; i < stop ; i++){
        for (unsigned long j = 0 ; j < i ; j++){
            unsigned long index = i * (i - 1) * 0.5 + j;
            d[index] = dist(p[i], p[j], 2, DIMENSION);
        }
    }
}

template<typename T>
double* threaded_distances_init(int numThreads, T** p, long NUM_PS, int DIMENSION, double (*dist)(T*, T*, int, int)){
    // Calculate the size of the 1D array for the upper triangular part
    unsigned long arraySize = NUM_PS * (NUM_PS - 1) * 0.5;
    double* d = new double[arraySize];

    // Create threads
    std::thread threads[numThreads];

    int start = 1;
    for (int i = 0 ; i < numThreads ; i++){
        int end = end_init(NUM_PS, start, numThreads);
        threads[i] = std::thread(distances_init<T>, p, d, start, end, DIMENSION, dist);
        start = end;
    }

    // Join threads
    for (int thr = 0; thr < numThreads; ++thr) {
        threads[thr].join();
    }

    return d;
}

template<typename T>
class DistComparator {
private:
    T** dataset;
    T** queryset;
    const int point;
    const int DIMENSION;

public:
    DistComparator(T** array, T** search_array, int pt, int dim)
        : dataset(array), queryset(search_array), point(pt), DIMENSION(dim) {}

    bool operator()(T a,T b) const {
        T* y = queryset[point];
        return dist(dataset[a], y, 2, DIMENSION) < dist(dataset[b], y, 2, DIMENSION);
    }
};

class Graph {
private:
    int vertices; // Number of vertices
    std::list<int>* adjacencyList; // Array of linked lists to represent adjacency list

public:
    // Constructor
    Graph(int vertices);

    // Destructor
    ~Graph();

    // Add a directed edge to the graph
    void addEdge(int source, int destination);

    std::vector<int>* get_node_nn(int node, int E);

    std::vector<int>* get_node_nn(int node);
};

void graph_init_file(Graph* gr, std::string filename);

template<typename T>
std::set<int, DistComparator<int>>* lsh_knn_mrng(T** p, std::unordered_multimap<int, int>** mm, int** w, double** t, int** rs, long* id, int query, int k, int L, int K, long M, int NO_IMAGES, int DIMENSION, bool skip){
    int countLSH = 0;
    std::set<int, DistComparator<int>>* m = new std::set<int, DistComparator<int>>(DistComparator<int>(p, p, query, DIMENSION));
    int inserted = 0;
    
    double min = DBL_MAX;
    for(int l = 0; l < L; l++){
        int key = g(p, w[l], t[l], rs[l], id, query, K, M, NO_IMAGES/8, DIMENSION, l);
        auto itr = mm[l]->equal_range(key);
        for (auto it = itr.first; it != itr.second; it++) {
            if (skip && query == it->second) continue;  // Useful if we are trying to find knn of the same file as the lsh init happened
            
            double currdist = dist(p[it->second], p[query], 2, DIMENSION);
            if (currdist < min){
                min = currdist;
                if (!m->empty()){
                    delete m;
                    m = new std::set<int, DistComparator<int>>(DistComparator<int>(p, p, query, DIMENSION));
                }
                
                inserted = 0;
            }
            
            if (currdist == min){
                m->insert(it->second);
                if (inserted == k){
                    auto last = m->rend();
                    m->erase(*last);
                }
                else{
                    inserted++;
                }
            }
            
            countLSH++;
            if( countLSH > 10*L ){
                break;
            }
        }
        
        if( countLSH > 10*L ){
            break;
        }
    }

    return m;
}

template <typename T>
void mrng(int start, int end, Graph* gr, T** pixels, double* distances, std::unordered_multimap<int, int>** mm, int** w, double** t, int** rs, long* id, int k, int L, int K, long M, int NO_IMAGES, int DIMENSION){
    std::unordered_set<int> Rp;
    std::cout << "ALALOUM" << std::endl;
    for (int j = 0 ; j < NO_IMAGES ; j++){
        Rp.insert(j);
    }
    
    for (int i = start ; i < end ; i++){
        std::cout << "i is " << i << std::endl;
        auto startmini = std::chrono::high_resolution_clock::now();
        std::set<int, DistComparator<int>>* Lp = lsh_knn_mrng(pixels, mm, w, t, rs, id, i, k, L, K, M, NO_IMAGES, DIMENSION, true);
        
        Rp.erase(i);    // Rp = S - {p}

        for (auto l : *Lp){
            Rp.erase(l);    // Rp = S - {p} - Lp
        }
        
        for (auto r : Rp){
            bool condition = true;
            double pr = distances[getIndex(i, r)];
            for (auto t = Lp->begin() ; t != Lp->end() ; t++){
                double rt = distances[getIndex(r, *t)];
                double pt = distances[getIndex(i, *t)];
                if ((pr > rt) && (pr > pt)){
                    condition = false;
                    break;
                }
            }
            if (condition == true){
                Lp->insert(r);
                if (Lp->size() > k){
                    auto last = Lp->end();
                    Lp->erase(*last);
                }
            }
        }
        
        for (auto it = Lp->begin(); it != Lp->end(); it++){
            gr->addEdge(i, *it);
        }
        
        Rp.insert(i);
        for (auto l : *Lp){
            Rp.insert(l);
        }
        auto stopmini = std::chrono::high_resolution_clock::now();

        double durationmini = std::chrono::duration_cast<std::chrono::milliseconds>(stopmini - startmini).count();
        durationmini *= 1e-3;
        std::cout << "time is : " << durationmini << " seconds" << std::endl;

        delete Lp;
    }
}

template <typename T>
void threaded_mrng(int numThreads, Graph* gr, T** pixels, std::unordered_multimap<int, int>** mm, int** w, double** t, int** rs, long* id, int k, int L, int K, long M, int NO_IMAGES, int DIMENSION, double (*dist)(T *, T *, int, int)){
    std::cout << "Creating the dist..." << std::endl;
    double* distances = threaded_distances_init(numThreads, pixels, NO_IMAGES, DIMENSION, dist);
    // double* distances = new double[NO_IMAGES];
    
    // Calculate chunk size for each thread
    int chunkSize = NO_IMAGES / numThreads;

    // Create threads
    std::thread threads[numThreads];
    std::cout << "mrng..." << std::endl;
    // Launch threads to fill array in parallel
    for (int thr = 0; thr < numThreads; thr++) {
        int start = thr * chunkSize;
        int end = (thr + 1) * chunkSize;
        threads[thr] = std::thread(mrng<T>, start, end, gr, pixels, distances, mm, w, t, rs, id, k, L, K, M, NO_IMAGES, DIMENSION);
    }

    // Join threads
    for (int thr = 0; thr < numThreads; ++thr) {
        threads[thr].join();
    }

    delete[] distances;
}