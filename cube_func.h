#include <iostream>
#include <random>
#include <unordered_map>
#include <bits/stdc++.h>
#include "generals.h"
#include "lsh_func.h"

#ifndef CUBE_FUNC_H
#define CUBE_FUNC_H

template <typename T>
void preprocess_cube(T** pixels, std::map<int, int>& hypervalues, std::unordered_multimap<long, int>& hypercube, int** w, double** t, int i, int L, int K, int dt, int DIMENSION){
    long key = 0;
    int bit = 1;
    for (int j = 0 ; j < dt ; j++){
        int zero_or_one;
        int hi = rand() % K;   // Ποια h θα επιλεξουμε
        int l = rand() % L;
        int num = h(pixels[i], w[l][hi], t[l][hi], hi, DIMENSION, l);
        auto itr = hypervalues.find(num);
        if (itr == hypervalues.end()){      // Αν δεν υπαρχει μες στο map
            zero_or_one = rand() % 2;
            hypervalues.insert({ num, zero_or_one });
        }
        else{
            zero_or_one = itr->second;
        }
        key += zero_or_one * bit;
        bit *= 2;
    }
    
    hypercube.insert({key,i});
}

template <typename T>
void cube_knn(Neibs<T>* cube_neibs, std::map<int, int> hypervalues, std::unordered_multimap<long, int> hypercube, int query_key, int dt, int M, int probes){
    int countVertices = 0;
    int countCubeElements = 0;
    int hamming_dist_wanted = 0;
    int i = 0;
    auto itr = hypercube.equal_range(query_key);
    std::bitset<32> binary_query_key(query_key);    // Convert to binary
    while (countCubeElements < M && countVertices < probes){
        for (auto it = itr.first; it != itr.second; it++) {
            cube_neibs->insertionsort_insert(it->second);
            countCubeElements++;
        }
        countVertices++;

        if (i == query_key){
            i++;
            i = i % dt;
        }
        if (i == 0) hamming_dist_wanted++;
        while (1){
            std::bitset<32> binary_i(i); // Convert i to binary representation
            int differingBits = (binary_query_key ^ binary_i).count(); // Count differing bits
            if (differingBits == hamming_dist_wanted){
                itr = hypercube.equal_range(i);
                break;
            }

            i++;
            i = i % dt;
            if (i == 0) hamming_dist_wanted++;
        }
        
    }
}

template <typename T>
void cube_rangeSearch(Neibs<T>* rangeSearchCube, std::map<int, int> hypervalues, std::unordered_multimap<long, int> hypercube, int query_key, int dt, int M, int probes, T** pixels, T** queries, int query, int DIMENSION, float R){
    int countVerticesRange = 0;
    int countCubeElementsRange = 0;
    int hamming_dist_wantedRange = 0;
    int i = 0;
    auto itr = hypercube.equal_range(query_key);
    std::bitset<32> binary_query_key(query_key);    // Convert to binary
    while (countCubeElementsRange < M && countVerticesRange < probes){
        for (auto it = itr.first; it != itr.second; it++) {
            countCubeElementsRange++;
            if( dist(pixels[it->second],queries[query],2,DIMENSION) >= R ) continue;
            
            rangeSearchCube->insertionsort_insert(it->second);
        }
        countVerticesRange++;

        if (i == query_key){
            i++;
            i = i % dt;
        }
        if (i == 0) hamming_dist_wantedRange++;
        while (1){
            std::bitset<32> binary_i(i); // Convert i to binary representation
            int differingBits = (binary_query_key ^ binary_i).count(); // Count differing bits
            if (differingBits == hamming_dist_wantedRange){
                itr = hypercube.equal_range(i);
                break;
            }

            i++;
            i = i % dt;
            if (i == 0) hamming_dist_wantedRange++;
        }
        
    }
}

template <typename T>
long query_key_init(T* vector, std::map<int, int>& hypervalues, int** w, double** t, int dt, int K, int L, int DIMENSION){
    long query_key = 0;
    int bit = 1;
    for (int j = 0 ; j < dt ; j++){
        int zero_or_one;
        int hi = rand() % K;   // Ποια h θα επιλεξουμε
        int l = rand() % L;
        int num = h(vector, w[l][hi], t[l][hi], hi, DIMENSION, l);
        auto itr = hypervalues.find(num);
        if (itr != hypervalues.end()){      // Αν δεν υπαρχει μες στο map
            zero_or_one = rand() % 2;
            hypervalues.insert({ num, zero_or_one });
        }
        else{
            zero_or_one = itr->second;
        }
        query_key += zero_or_one * bit;
        bit *= 2;
    }

    return query_key;
}

#endif