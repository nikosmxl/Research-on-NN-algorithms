#include "cluster_func.h"

int readconf(std::string conf, int& number_of_clusters, int& L, int& K, int& Mcube, int& dt, int& probes){
    std::ifstream confFile(conf); // Open the file for reading
    if (!confFile) {
        std::cerr << "Error: Could not open the file." << std::endl;
        return 1;
    }

    std::string line;
    int i = 0;
    while (std::getline(confFile, line)) { // Read each line from the file
        char* token = std::strtok(&line[0], " \t\n"); // Tokenize the line using space, tab, and newline as delimiters
        int flag = -1;
        while (token != nullptr) {
            if (strcmp(token, "//") == 0) break;

            if (flag == -1){
                if (strcmp(token, "number_of_clusters:") == 0){
                    flag = 0;
                }
                if (strcmp(token, "number_of_vector_hash_tables:") == 0){
                    flag = 1;
                }
                if (strcmp(token, "number_of_vector_hash_functions:") == 0){
                    flag = 2;
                }
                if (strcmp(token, "max_number_M_hypercube:") == 0){
                    flag = 3;
                }
                if (strcmp(token, "number_of_hypercube_dimensions:") == 0){
                    flag = 4;
                }
                if (strcmp(token, "number_of_probes:") == 0){
                    flag = 5;
                }
            }
            else{
                if (flag == 0){
                    number_of_clusters = atoi(token);
                }
                else if (flag == 1){
                    L = atoi(token);
                }
                else if (flag == 2){
                    K = atoi(token);
                }
                else if (flag == 3){
                    Mcube = atoi(token);
                }
                else if (flag == 4){
                    dt = atoi(token);
                }
                else if (flag == 5){
                    probes = atoi(token);
                }
                flag = -1;
            }
            token = std::strtok(nullptr, " \t\n"); // Get the next token
        }
        i++;
    }

    confFile.close(); // Close the file
    return 0;
}