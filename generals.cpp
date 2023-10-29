#include "generals.h"

int readfile(std::string filename, int** pixels, int images, int dim){

    std::ifstream rf(filename, std::ios::binary);
    if (rf.is_open()) {
        char buffer[dim];
        rf.read(buffer, 16);    // Τα αρχικα bytes με τις πληροφοριες που ηδη μας ειναι γνωστες και δεν χρειαζομαστε

        // Check if the read operation was successful
        if (rf) {
            for (int i = 0 ; i < images ; i++){
                rf.read(buffer, dim);
                for (int j = 0; j < dim; j++) {
                    pixels[i][j] = static_cast<int>(static_cast<unsigned char>(buffer[j]));
                }
            }
        } else {
            std::cerr << "Error reading file." << std::endl;
        }
        
        // Close the file
        rf.close();
    } else {
        std::cerr << "Error opening file." << std::endl;
    }

    return 0;
}

double dist(int* x,int* y, int k, int DIM){
    unsigned long long sum = 0;
    for(int i = 0; i < DIM; i++){
        sum += pow(abs(x[i] - y[i]), k);
    }
    double d = pow(sum, 1.0/k);
    return d;
}

Neibs::Neibs(int** pixels, int** qs, int dimension, int s,int query, double (*func)(int*, int*, int, int)) : p(pixels), queries(qs), DIM(dimension), size(0), maxsize(s), q(query), distfunc(func) {
    nn = new int[s];
}

Neibs::~Neibs(){
    delete[] nn;
}

int* Neibs::givenn() const{
    return nn;
}

int Neibs::givenn(int i) const{
    return nn[i];
}

int Neibs::give_size() const{
    return size;
}

void Neibs::insertionsort_insert(int v){

    if( size == 0 ){
        nn[size++] = v;
        return;
    }

    if( size == maxsize ){
        if( distfunc(queries[q],p[nn[maxsize - 1]], 2, DIM) <= distfunc(queries[q],p[v], 2, DIM) ){
            return;
        }
    }
    
    for(int i = 0; i < size; i++){
        if( distfunc(queries[q],p[nn[i]], 2, DIM) > distfunc(queries[q],p[v], 2, DIM) ){ // αν ειναι το v κοντινοτερο σε σχεση με το nn[i]?
            int temp = v;
            for(int j = i; j < size; j++){
                int temp2 = nn[j];
                nn[j] = temp;
                temp = temp2;
            }

            if (size < maxsize){
                nn[size] = temp;
                size++;
            }
            
            return;
        }
    }
    if( size < maxsize){
        nn[size] = v;
        size++;
        return;
    }
    return;
}

double Neibs::givedist(int i) const{
        return distfunc(queries[q],p[nn[i]], 2, DIM);
    }