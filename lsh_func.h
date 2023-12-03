#include <iostream>
#include <random>
#include <unordered_map>
#include "generals.h"

#ifndef LSH_FUNC_H
#define LSH_FUNC_H

template<typename T>
int h(T* p, int w, double t, int K, int DIM, int L){
    std::default_random_engine generator;
    generator.seed(K + 10*L + 1);           // Για να σωσουμε μνημη χωρις να χρειαζομαστε τον 3διαστατο πινακα v.
    std::normal_distribution<double> distribution(0.0, 1.0);
    double multiplyvexts = 0;
    for(int i = 0; i < DIM; i++){
        double v = fabs(distribution(generator));  // v vector με τιμες απο κανονικη κατανομη
        multiplyvexts += v*p[i];
    }

    int hp = ((int)(  (multiplyvexts + t)/(w*1.0)  ));  // κατω ακεραιο μερος οπως λεει στον τυπο στις διαφανειες
    
    return hp;
}

template<typename T>
int g(T** pixels, int* w, double* t, int* rs, long* id, int j, int K, long M, int tablesize, int DIM, int L){
    unsigned long long sum = 0;
    for(int i = 0; i < K; i++){
        sum += rs[i]*h(pixels[j], w[i], t[i], i, DIM, L) % M;
    }
    
    id[j] = sum;
    return (id[j] % tablesize);
}

template<typename T>
void lsh_init(T** pixels, int** w, double** t, int** rs, std::unordered_multimap<int, int>** mm, long* id, int L, int K, long M, int NO_IMAGES, int DIMENSION){
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

    for (int l = 0 ; l < L ; l++){
        mm[l] = new std::unordered_multimap<int, int>();
        for (int j = 0 ; j < NO_IMAGES ; j++){
            int key = g(pixels, w[l], t[l], rs[l], id, j, K, M, NO_IMAGES/8, DIMENSION, l);
            mm[l]->insert({key, j});
        }
    }
}

template<typename T>
void lsh_knn(T** p, std::unordered_multimap<int, int>** mm, Neibs<T>* lsh, int** w, double** t, int** rs, long* id, int query, int L, int K, long M, int NO_IMAGES, int DIMENSION, bool skip){
    int countLSH = 0;
    for(int l = 0; l < L; l++){
        int key = g(p, w[l], t[l], rs[l], id, query, K, M, NO_IMAGES/8, DIMENSION, l);
        auto itr = mm[l]->equal_range(key);
        for (auto it = itr.first; it != itr.second; it++) {
            if (skip && query == it->second) continue;  // Useful if we are trying to find knn of the same file as the lsh init happened
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
}

template<typename T>
void lsh_rangeSearch(T** pixels, T** queries, std::unordered_multimap<int, int>** mm, Neibs<T>* rangeSearch, int** w, double** t, int** rs, long* id, int query, int L, int K, long M, int NO_IMAGES, int DIMENSION, float R){
    int countRange = 0;
    for(int l = 0; l < L; l++){
        int key = g(queries, w[l], t[l], rs[l], id, query, K, M, NO_IMAGES/8, DIMENSION, l);
        auto itr = mm[l]->equal_range(key);
        auto it = itr.first;
        for (int i = 0; i != mm[l]->count(key) - 1; i++) {
            countRange++;
            if (it->second == query){
                if (++it == itr.second) break;
            }
            if ( dist(pixels[it->second],queries[query],2,DIMENSION) >= R ) continue;
                rangeSearch->insertionsort_insert(it->second);
            if( countRange > 20*L ){
                break;
            }
            it++;
        }
        
        if( countRange > 20*L ){
            break;
        }
    }
}

#endif