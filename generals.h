#include <iostream>
#include <fstream>
#include <bits/stdc++.h>

#ifndef GENERALS_H
#define GENERALS_H

template<typename T>
T** readfile(std::string filename, std::size_t images, std::size_t dim){
    T** array = new T*[images];

    // Allocate memory for each row
    for (std::size_t i = 0; i < images; ++i) {
        array[i] = new T[dim];
    }

    std::ifstream rf(filename, std::ios::binary);
    if (rf.is_open()) {
        char buffer[dim];
        rf.read(buffer, 16);    // Τα αρχικα bytes με τις πληροφοριες που ηδη μας ειναι γνωστες και δεν χρειαζομαστε

        // Check if the read operation was successful
        if (rf) {
            for (std::size_t i = 0; i < images; i++){
                rf.read(buffer, dim);
                for (std::size_t j = 0; j < dim; j++) {
                    array[i][j] = static_cast<T>(static_cast<unsigned char>(buffer[j]));
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

    return array;
}

template<typename T>
double dist(T* x,T* y, int k, int DIM){
    unsigned long long sum = 0;
    for(int i = 0; i < DIM; i++){
        sum += pow(abs(x[i] - y[i]), k);
    }
    double d = pow(sum, 1.0/k);
    return d;
}

template<typename T>
class Neibs
{
    private:
        T** p;
        T** queries;
        int DIM;
        int* nn;
        std::size_t size;
        std::size_t maxsize;
        int q;
        double (*distfunc)(T*, T*, int, int);
    public:
    Neibs(T** pixels, T** qs, int dimension, std::size_t s,int query, double (*func)(T*, T*, int, int)) : p(pixels), queries(qs), DIM(dimension), size(0), maxsize(s), q(query), distfunc(func) {
        nn = new int[s];
    }

    ~Neibs(){
        delete[] nn;
    }

    int* givenn() const{
        return nn;
    }

    int givenn(int i) const{
        return nn[i];
    }

    std::size_t give_size() const{
        return size;
    }

    void insertionsort_insert(int v){
        if( size == 0 ){
            nn[size++] = v;
            return;
        }

        if( size == maxsize ){
            if( distfunc(queries[q],p[nn[maxsize - 1]], 2, DIM) <= distfunc(queries[q],p[v], 2, DIM) ){
                return;
            }
        }
        
        for(std::size_t i = 0; i < size; i++){
            if( distfunc(queries[q],p[nn[i]], 2, DIM) > distfunc(queries[q],p[v], 2, DIM) ){ // αν ειναι το v κοντινοτερο σε σχεση με το nn[i]?
                int temp = v;
                for(std::size_t j = i; j < size; j++){
                    int temp2 = nn[j];
                    nn[j] = temp;
                    temp = temp2;
                }

                if (size < maxsize){
                    nn[size++] = temp;
                }
                
                return;
            }
        }
        if( size < maxsize){
            nn[size++] = v;
            return;
        }
        return;
    }

    double givedist(int i) const{
        return distfunc(queries[q],p[nn[i]], 2, DIM);
    }
};

#endif