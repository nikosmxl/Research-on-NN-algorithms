#include "lsh_func.h"

int h(int* p, int w, double t, int K, int DIM, int L){
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
int g(int** pixels, int* w, double* t, int* rs, long* id, int j, int K, long M, int tablesize, int DIM, int L){
    unsigned long long sum = 0;
    for(int i = 0; i < K; i++){
        sum += rs[i]*h(pixels[j], w[i], t[i], i, DIM, L);
    }
    
    id[j] = sum % M;
    return (id[j] % tablesize);
}