#include "Data.h"
#include <iostream>
#include <vector>
#include <time.h> 
#include <cstdlib>
#include <algorithm>

using namespace std;

/*Resolver o bestImprovenmentSwap que nao termina a execução.
*/

typedef struct Solucao {
        vector<int> sequence;
        double valorObj= 0.0;
    } Solucao;

void ExibirSolucao(Solucao *s){
    for(int i = 0; i < s->sequence.size() - 1;i++){
        cout << s->sequence[i] << "->";    
    }
    cout << s->sequence.back() <<  endl;
}
double Latencia (int i, Solucao *s, Data& data){
    double latencia = 0;
    for (int j = 0; j < i - 1; j++){
        latencia += data.getDistance(s->sequence[j], s->sequence[j+1]);
        
    }
    return latencia;
}
void CalculaValorObj(Data& data, Solucao *s){
    s->valorObj = 0;
    for(int i = 1; i <= data.getDimension() + 1; i++){
        s->valorObj += Latencia(i, s, data);
    }
}

typedef struct Subsequence {
    double T, C;
    int W;
    int first,last;
    inline static Subsequence Concatenate(Subsequence& sigma_1, Subsequence& sigma_2, Data& data){
        Subsequence sigma;
        double temp =  data.getDistance(sigma_1.last, sigma_2.first);
        sigma.W = sigma_1.W + sigma_2.W;
        sigma.T = sigma_1.T + temp + sigma_2.T;
        sigma.C = sigma_1.C + sigma_2.W * (sigma_1.T + temp) + sigma_2.C;
        sigma.first = sigma_1.first;
        sigma.last = sigma_2.last;
        return sigma;
    }
} Subsequence;
typedef struct InsertionInfo{
    int noInserido;
    int vertceAnterior;
    double custo;
}InsertionInfo;

void UpdateAllSubseq(Solucao *s, vector<vector<Subsequence>>& subseq_matrix, Data& data){
    int n = s->sequence.size();
    for(int i = 0; i < n; i++){
        int v = s->sequence[i];
        subseq_matrix[i][i].W = (i > 0);
        subseq_matrix[i][i].C = 0;
        subseq_matrix[i][i].T = 0;
        subseq_matrix[i][i].first = s->sequence[i];
        subseq_matrix[i][i].last = s->sequence[i];
        }
    for(int i = 0; i < n; i++){
        for(int j = i+1; j < n; j++){
            subseq_matrix[i][j] = Subsequence::Concatenate(subseq_matrix[i][j-1], subseq_matrix[j][j], data);
        }
    }

    for (int i = n - 1; i >= 0; i--){
        for (int j = i - 1; j >= 0; j--){
            subseq_matrix[i][j] = Subsequence::Concatenate(subseq_matrix[i][j+1], subseq_matrix[j][j], data);
        }
    }
    // cout << "4"
}
vector<int> nosRestantes(size_t vertices){
    vector<int> CL;
    for(size_t j = 1; j <= vertices; j++){
        CL.push_back(j);
    }
    return CL;
}

vector<InsertionInfo> CalcularCusto (int r, vector<int>& CL, Data& data){
    vector<InsertionInfo> custoInsercao(CL.size());
    int l = 0;
    for(auto k : CL){
        custoInsercao[l].custo = data.getDistance(k, r);
        custoInsercao[l].noInserido = k;
        l++;
    }
    return custoInsercao;
}
bool myfunction (InsertionInfo i, InsertionInfo j) { return (i.custo<j.custo); }

void inserirNaSolucao (Solucao *s, int c, int r, vector<InsertionInfo>& custoInsercao){
    for(int i = 0; i < s->sequence.size(); i++){
        if (s->sequence[i] == r){
            s->sequence.insert(s->sequence.begin() + i + 1, custoInsercao[c].noInserido);
            break;
        }
    }

}

Solucao Construcao(Data& data){
    Solucao s = {{1,1}, 0.0};
    size_t n = data.getDimension();
    vector<int> CL = nosRestantes(n);
    CL.erase(CL.begin());
    int r = 1;
    while (CL.size() != 0){
        vector <InsertionInfo> custoInsercao = CalcularCusto(r, CL, data);
        sort(custoInsercao.begin(), custoInsercao.end(), myfunction);
        double alpha = (double) rand() / RAND_MAX;
        int c = rand() % ((int) ceil((alpha + 0.000001) * custoInsercao.size()));
        inserirNaSolucao(&s, c, r, custoInsercao);
        r = custoInsercao[c].noInserido;
        for (int i = 0; i < CL.size(); i++){
            if(CL[i] == custoInsercao[c].noInserido){
                CL.erase(CL.begin() + i);
                break;
            }
        }
    }
    return s;
}
bool bestImprovenmentSwap(Solucao *s,vector<vector<Subsequence>>& subseq_matrix, Data& data){
    double bestCost = s->valorObj;
    double cost;
    int best_i,best_j;
    for(int i = 1; i < s->sequence.size() - 1; i++){
   
        for(int j = i+1; j < s->sequence.size() - 1;j++){
            
            if((i + 1)!= j){
                Subsequence sigma_1 = Subsequence :: Concatenate(subseq_matrix[0][i-1], subseq_matrix[j-1][j-1], data);
                
                Subsequence sigma_2 = Subsequence :: Concatenate(sigma_1, subseq_matrix[i+1][j-1], data);
                
                Subsequence sigma_3 = Subsequence :: Concatenate(sigma_2, subseq_matrix[i][i], data);
                
                //s->sequence.size() -1 ou -2???????
                Subsequence sigma_4 = Subsequence :: Concatenate(sigma_3, subseq_matrix[j+1][s->sequence[s->sequence.size()-1]], data);
                
                cost = sigma_4.C;
                
            } else {
                
                Subsequence sigma_1 = Subsequence :: Concatenate(subseq_matrix[0][vi_prev-1],subseq_matrix[vj-1][vi-1], data);
                
                //s->sequence.size() -1 ou -2???????
                Subsequence sigma_2 = Subsequence :: Concatenate(sigma_1, subseq_matrix[vj_next-1][s->sequence.size()-1], data);
                
                cost = sigma_2.C;
            }
            if(cost < bestCost){
                bestCost = cost;
                best_i = i;
                best_j = j;
            }
        }
    }
    
    if(bestCost < s->valorObj){
        swap(s->sequence[best_i], s->sequence[best_j]);
        
        s->valorObj = bestCost;
        UpdateAllSubseq(s, subseq_matrix, data);
        return true;
    }
    return false;
}
bool bestImprovement2Opt(Solucao *s,vector<vector<Subsequence>>& subseq_matrix, Data& data){
    double bestCost = s->valorObj;
    int best_i;
    int best_j;
    for(int i = 1; i < s->sequence.size() - 1;i++){
       //i = primeiro elemento da aresta a ser invertida
        for(int j = i + 3; j < s->sequence.size() - 1; j++){
            cout << "i: " << i << " j: " << j<< endl;
            // j = ultimo elemento da aresta a ser invertida
            Subsequence sigma_1;
            sigma_1 = Subsequence :: Concatenate(subseq_matrix[0][i-1], subseq_matrix[j][i], data);
            Subsequence sigma_2 = Subsequence :: Concatenate(sigma_1, subseq_matrix[j+1][s->sequence.size()-1], data);
            cout << sigma_2.C<< endl;
            if(sigma_2.C < bestCost){
                bestCost = sigma_2.C;
                best_i = i;
                best_j = j;
            }
        }
        
    }
    if(bestCost < s->valorObj){
        int cont = 0;
        while((best_i + cont) < (best_j - cont)){
            swap(s->sequence[best_i + cont], s->sequence[best_j - cont]);
            cont++;
        }
        s->valorObj = bestCost;
        UpdateAllSubseq(s, subseq_matrix, data);
        return true;
    }
    return false;
}
bool bestImprovementOrOpt(Solucao *s, vector<vector<Subsequence>>& subseq_matrix, int tam_bloco){
    double bestDelta = 0;
    int best_i1;
    int best_i2;
    int best_j;
    for(int i = 1; i < s->sequencia.size() - 1;i++){
        int vi1 = s->sequencia[i];
        int vi1_prev = s->sequencia[i-1];
        int vi2 = s->sequencia[i + tam_bloco - 1];
        int vi2_next = s->sequencia[i + tam_bloco];
        for(int j = i + tam_bloco; j < s->sequencia.size() - 1; j++){
            int vj = s->sequencia[j];
            int vj_next = s->sequencia[j+1];
            Subsequence sigma_1 = Subsequence :: Concatenate(subseq_matrix)
            double delta = -data.getDistance(vi1_prev, vi1) - data.getDistance(vi2, vi2_next) - data.getDistance(vj, vj_next) + data.getDistance(vi1_prev, vi2_next) + data.getDistance(vj, vi2) + data.getDistance(vi1, vj_next);
            double delta = -Matriz[vi1_prev-1][vi1-1] - Matriz[vi2-1][vi2_next-1]- Matriz[vj-1][vj_next-1] + Matriz[vi1_prev-1][vi2_next-1] + Matriz[vj-1][vi2-1] + Matriz[vi1-1][vj_next-1];
            if (delta < bestDelta){
                bestDelta = delta;
                best_i1 = i;
                best_i2 = i+tam_bloco -1;
                best_j = j;
            }
        }
    }

    int cont = 0;
    if(bestDelta < 0){
        do{
        s->sequencia.insert(s->sequencia.begin() + best_j + 1, s->sequencia[best_i1]);
        s->sequencia.erase(s->sequencia.begin() + best_i1);
        best_j--;
        cont++;
        }while(cont < tam_bloco);
        s->valorObj = s->valorObj + bestDelta;
        return true;
    }
    return false;
}
void BuscaLocal(Solucao *s, vector<vector<double>>& Matriz){
    vector<int> NL = {1, 2, 3, 4, 5};
    bool improved = false;

    while (NL.empty() == false){
        int n = rand() % NL.size();
        switch (NL[n]){
            case 1:
                improved = bestImprovementSwap(s, Matriz);
                break;
            case 2:
                improved = bestImprovement2Opt(s, Matriz);
                break;
            case 3:
                improved = bestImprovementOrOpt(s, Matriz, 1);
                break;
            case 4: 
                improved = bestImprovementOrOpt(s, Matriz, 2);
                break;
            case 5:
                improved = bestImprovementOrOpt(s, Matriz, 3);
                break;
        }
        if(improved){
            NL = {1,2,3,4,5};
        }else{
            NL.erase(NL.begin() + n);
        }
    }
}

Solucao ILS(int maxIter, int maxIterIls, Data& data){
    Solucao bestOfAll;
    bestOfAll.valorObj = INFINITY;
    for(int i = 0; i < maxIter;i++){
        Solucao s = Construcao(data);
        CalculaValorObj(data, &s);
       /* Solucao best = s;
        int iterIls = 0;
        
        while(iterIls <= maxIterIls){
            
            BuscaLocal(&s, Matriz);
            if(s.valorObj < best.valorObj){
                best = s;
                iterIls = 0;
            }
            
            s = Perturbacao(best, Matriz);
            iterIls++;
        }
    
        if(best.valorObj < bestOfAll.valorObj)
            bestOfAll = best;*/
    }
    return bestOfAll;   
}
int main(int argc, char** argv) {
    auto data = Data(argc, argv[1]);
    data.read(); 
    srand(time(NULL));
    
    Solucao s = Construcao(data);
    ExibirSolucao(&s);
    CalculaValorObj(data, &s);
    cout << s.valorObj << endl;
    int n = s.sequence.size();
    vector<vector<Subsequence>> subseq_matrix(n, vector<Subsequence>(n));
    UpdateAllSubseq(&s, subseq_matrix, data);
    for (int k = 0; k < n; k++){
        for (int l = 0 ; l < n; l++){
            cout <<"i: " << k << " j: "<< l <<  " W: " << subseq_matrix[k][l].W << " T:" << subseq_matrix[k][l].T << " C: " << subseq_matrix[k][l].C << endl;
        }
     }
    cout << endl;
    bool troca = bestImprovement2Opt(&s, subseq_matrix, data);
    //CalculaValorObj(data, &s);
    cout << s.valorObj << endl;
    ExibirSolucao(&s);

    bool troca2 = bestImprovenmentSwap(&s, subseq_matrix, data);
     cout << troca2 << " " << s.valorObj << endl;
    ExibirSolucao(&s);

    /*
    bool troca2 = bestImprovementSwap(&s, subseq_matrix, data);
    CalculaValorObj(data, &s);
    cout << troca2 << " " << s.valorObj << endl;
    ExibirSolucao(&s);*/
    return 0;
}






















/*
vector<int> escolher3NosAleatorios(size_t vertices){
    vector<int> sequencia;
    
    sequencia.push_back(1);
    sequencia.push_back(1);

    sequencia.insert(sequencia.begin() + 1, (rand() % (vertices) + 1));
    while(sequencia[0] == sequencia[1]){
        sequencia.erase(sequencia.begin() + 1);
        sequencia.insert(sequencia.begin() + 1, (rand() % (vertices) + 1));
    }

    sequencia.insert(sequencia.begin() + 2, (rand() % (vertices) + 1));
    while(sequencia[1] == sequencia[2] || sequencia[2] == sequencia[0]){
        sequencia.erase(sequencia.begin() + 2);
        sequencia.insert(sequencia.begin() + 2, (rand() % (vertices) + 1));
    }
    
    sequencia.insert(sequencia.begin() + 3, (rand() % (vertices) + 1));
    while(sequencia[1] == sequencia[3] || sequencia[2] == sequencia[3] || sequencia[0] == sequencia[3]){
        sequencia.erase(sequencia.begin() + 3);
        sequencia.insert(sequencia.begin() + 3, (rand() % (vertices) + 1));
    }
    
    

    return sequencia;
}*/
/*
void ExibirSolucao(Solucao *s){
    for(int i = 0; i < s->sequencia.size() - 1;i++){
        cout << s->sequencia[i] << "->";    
    }
    cout << s->sequencia.back() <<  endl;
}

}*/
/*
vector<InsertionInfo> calcularCustoInsercao(Solucao& s, vector<int>& CL, vector<vector<double>>& Matriz){
    vector<InsertionInfo> custoInsercao((s.sequencia.size() - 1) * CL.size());
    int l = 0;

    for(int a = 0; a < s.sequencia.size() - 1;a++){
        int i = s.sequencia[a];
        int j = s.sequencia[a + 1];
        for (auto k : CL){
            //custoInsercao[l].custo = data.getDistance(i, k) + data.getDistance(k, j) - data.getDistance(i, j);
            custoInsercao[l].custo = Matriz[i-1][k-1] + Matriz[k-1][j-1] - Matriz[i-1][j-1];
            custoInsercao[l].arestaRemovida = a;
            custoInsercao[l].noInserido = k;
            l++;
        }
    }
    return custoInsercao;
}

void inserirNaSolucao(Solucao& s, vector<InsertionInfo>& custoInsercao, vector<int>& CL, int selecionado){
    for(int i = 0; i < s.sequencia.size();i++){
        if(s.sequencia[i] == s.sequencia[custoInsercao[selecionado].arestaRemovida]){
            s.sequencia.insert(s.sequencia.begin() + i + 1, custoInsercao[selecionado].noInserido);
            break;
        }
    }
    for(int i = 0; i < CL.size();i++){
        if(CL[i] == custoInsercao[selecionado].noInserido){
            CL.erase(CL.begin() + i);
        }
    }
}*/
/*
bool myfunction (InsertionInfo i, InsertionInfo j) { return (i.custo<j.custo); }

Solucao Construcao(int vertices, vector<vector<double>>& Matriz){
    double tempo_ordenacao = 0;
    Solucao s = {{}, 0.0};
    s.sequencia = escolher3NosAleatorios(vertices);
    vector<int> CL = nosRestantes(&s, vertices);
    while(!CL.empty()){
        vector<InsertionInfo> custoInsercao = calcularCustoInsercao(s, CL, Matriz);
        sort(custoInsercao.begin(), custoInsercao.end(), myfunction);
        double alpha = (double) rand() / RAND_MAX;
        int selecionado = rand() % ((int) ceil((alpha + 0.000001) * custoInsercao.size()));
        inserirNaSolucao(s, custoInsercao, CL, selecionado); 
    }
    return s;
}*/
/*
bool bestImprovementSwap(Solucao *s, vector<vector<double>>& Matriz){
    double bestDelta = 0;
    int best_i,best_j;
    for(int i = 1; i < s->sequencia.size() - 1; i++){
        int vi = s->sequencia[i];
        int vi_next = s->sequencia[i+1];
        int vi_prev = s->sequencia[i-1];
        for(int j = i+1; j < s->sequencia.size() - 1;j++){
            int vj = s->sequencia[j];
            int vj_next = s->sequencia[j+1];
            int vj_prev = s->sequencia[j-1];
            double delta = 0;
            if(i + 1!= j){
                delta = -data.getDistance(vi_prev,vi) - data.getDistance(vi, vi_next) + data.getDistance(vi_prev, vj) + data.getDistance(vj, vi_next) - data.getDistance(vj_prev, vj) - data.getDistance(vj, vj_next) + data.getDistance(vj_prev, vi) + data.getDistance(vi, vj_next);
                delta = -Matriz[vi_prev-1][vi-1] - Matriz[vi-1][vi_next-1] + Matriz[vi_prev-1][vj-1] + Matriz[vj-1][vi_next-1] - Matriz[vj_prev-1][vj-1] - Matriz[vj-1][vj_next-1] + Matriz[vj_prev-1][vi-1] + Matriz[vi-1][vj_next-1];
            } else {
                delta = -data.getDistance(vi_prev,vi) + data.getDistance(vi_prev, vj) - data.getDistance(vj, vj_next)  + data.getDistance(vi, vj_next);
                delta = -Matriz[vi_prev-1][vi-1] + Matriz[vi_prev-1][vj-1] - Matriz[vj-1][vj_next-1] + Matriz[vi-1][vj_next-1];
            }
            if(delta < bestDelta){
                bestDelta = delta;
                best_i = i;
                best_j = j;
            }
        }
    }
    
    if(bestDelta < 0){
        swap(s->sequencia[best_i], s->sequencia[best_j]);
        s->valorObj = s->valorObj + bestDelta;
        return true;
    }
    return false;
}*/
/*
bool bestImprovement2Opt(Solucao *s, vector<vector<double>>& Matriz){
    double bestDelta = 0;
    int best_i;
    int best_j;
    for(int i = 1; i < s->sequencia.size() - 1;i++){
        i = primeiro elemento da aresta a ser invertida
        int vi = s->sequencia[i];
        int vi_prev = s->sequencia[i-1];
        for(int j = i + 3; j < s->sequencia.size() - 1; j++){
             j = ultimo elemento da aresta a ser invertida
            int vj = s->sequencia[j];
            int vj_next = s->sequencia[j+1];
            double delta = - data.getDistance(vi_prev, vi) - data.getDistance(vj, vj_next) + data.getDistance(vi_prev, vj) + data.getDistance(vi, vj_next);
            double delta = -Matriz[vi_prev-1][vi-1] - Matriz[vj-1][vj_next-1] + Matriz[vi_prev-1][vj-1] + Matriz[vi-1][vj_next-1];
            if(delta < bestDelta){
                bestDelta = delta;
                best_i = i;
                best_j = j;
            }
        }
    }
    if(bestDelta < 0){
        int cont = 0;
        while((best_i + cont) < (best_j - cont)){
            swap(s->sequencia[best_i + cont], s->sequencia[best_j - cont]);
            cont++;
        }
        s->valorObj = s->valorObj + bestDelta;
        return true;
    }
    return false;
}*/
/*
bool bestImprovementOrOpt(Solucao *s, vector<vector<double>>& Matriz, int tam_bloco){
    double bestDelta = 0;
    int best_i1;
    int best_i2;
    int best_j;
    for(int i = 1; i < s->sequencia.size() - 1;i++){
        int vi1 = s->sequencia[i];
        int vi1_prev = s->sequencia[i-1];
        int vi2 = s->sequencia[i + tam_bloco - 1];
        int vi2_next = s->sequencia[i + tam_bloco];
        for(int j = i + tam_bloco; j < s->sequencia.size() - 1; j++){
            int vj = s->sequencia[j];
            int vj_next = s->sequencia[j+1];
            double delta = -data.getDistance(vi1_prev, vi1) - data.getDistance(vi2, vi2_next) - data.getDistance(vj, vj_next) + data.getDistance(vi1_prev, vi2_next) + data.getDistance(vj, vi2) + data.getDistance(vi1, vj_next);
            double delta = -Matriz[vi1_prev-1][vi1-1] - Matriz[vi2-1][vi2_next-1]- Matriz[vj-1][vj_next-1] + Matriz[vi1_prev-1][vi2_next-1] + Matriz[vj-1][vi2-1] + Matriz[vi1-1][vj_next-1];
            if (delta < bestDelta){
                bestDelta = delta;
                best_i1 = i;
                best_i2 = i+tam_bloco -1;
                best_j = j;
            }
        }
    }

    int cont = 0;
    if(bestDelta < 0){
        do{
        s->sequencia.insert(s->sequencia.begin() + best_j + 1, s->sequencia[best_i1]);
        s->sequencia.erase(s->sequencia.begin() + best_i1);
        best_j--;
        cont++;
        }while(cont < tam_bloco);
        s->valorObj = s->valorObj + bestDelta;
        return true;
    }
    return false;
}*/
/*
void BuscaLocal(Solucao *s, vector<vector<double>>& Matriz){
    vector<int> NL = {1, 2, 3, 4, 5};
    bool improved = false;

    while (NL.empty() == false){
        int n = rand() % NL.size();
        switch (NL[n]){
            case 1:
                improved = bestImprovementSwap(s, Matriz);
                break;
            case 2:
                improved = bestImprovement2Opt(s, Matriz);
                break;
            case 3:
                improved = bestImprovementOrOpt(s, Matriz, 1);
                break;
            case 4: 
                improved = bestImprovementOrOpt(s, Matriz, 2);
                break;
            case 5:
                improved = bestImprovementOrOpt(s, Matriz, 3);
                break;
        }
        if(improved){
            NL = {1,2,3,4,5};
        }else{
            NL.erase(NL.begin() + n);
        }
    }
}*/
/*
Solucao Perturbacao (Solucao best, vector<vector<double>>& Matriz){

    Solucao s = best;
    int vertices = best.sequencia.size() - 1;
    int i, j, tam_i, tam_j;
    //escolhe os valores aleatorios para os dois segmentos e seus tamanhos
    do{
        i = (rand() % (vertices - 1)) + 1;
        if(vertices/10 > 2){
            tam_i = (rand() % ((vertices/10) - 2) + 2);
        } else{
            tam_i = 2;
        }

    }while((i + tam_i) > vertices);
    do{
        j = (rand() % (vertices - 1)) + 1;
        if(vertices/10 > 2){
            tam_j = (rand() % (vertices/10 - 2) + 2);
        }
        else{
            tam_j = 2;
        }
    
    }while(((i + tam_i) > j && (j + tam_j) > i) || ((j + tam_j) > vertices));
    
    double delta;
    if((i+ tam_i) == j || (j +tam_j) == i){
        if( i < j){
            delta = - Matriz[s.sequencia[i]-1][s.sequencia[i-1]-1] - Matriz[s.sequencia[i+tam_i -1]-1][s.sequencia[i+tam_i]-1] -Matriz[s.sequencia[j+tam_j-1]-1][s.sequencia[j+tam_j]-1] + Matriz[s.sequencia[i-1]-1][s.sequencia[j]-1] + Matriz[s.sequencia[j+tam_j-1]-1][s.sequencia[i]-1] + Matriz[s.sequencia[i+tam_i-1]-1][s.sequencia[j+tam_j]-1];
        }
        else{
            delta = -Matriz[s.sequencia[j]-1][s.sequencia[j-1]-1] - Matriz[s.sequencia[j+tam_j -1]-1][s.sequencia[j+tam_j]-1] -Matriz[s.sequencia[i+tam_i-1]-1][s.sequencia[i+tam_i]-1] + Matriz[s.sequencia[j-1]-1][s.sequencia[i]-1] + Matriz[s.sequencia[i+tam_i-1]-1][s.sequencia[j]-1] + Matriz[s.sequencia[j+tam_j-1]-1][s.sequencia[i+tam_i]-1];
        }
    }else{
        if(i < j){
            delta = -Matriz[s.sequencia[i]-1][s.sequencia[i-1]-1] - Matriz[s.sequencia[i+tam_i-1]-1][s.sequencia[i+tam_i]-1] - Matriz[s.sequencia[j]-1][s.sequencia[j-1]-1] - Matriz[s.sequencia[j+tam_j]-1][s.sequencia[j+tam_j-1]-1] + Matriz[s.sequencia[i-1]-1][s.sequencia[j]-1] + Matriz[s.sequencia[j+tam_j-1]-1][s.sequencia[i+tam_i]-1] + Matriz[s.sequencia[j-1]-1][s.sequencia[i]-1] + Matriz[s.sequencia[i+tam_i-1]-1][s.sequencia[j+tam_j]-1];
        }
        else{
            delta = -Matriz[s.sequencia[j]-1][s.sequencia[j-1]-1] - Matriz[s.sequencia[j+tam_j-1]-1][s.sequencia[j+tam_j]-1] - Matriz[s.sequencia[i]-1][s.sequencia[i-1]-1] - Matriz[s.sequencia[i+tam_i]-1][s.sequencia[i+tam_i-1]-1] + Matriz[s.sequencia[j-1]-1][s.sequencia[i]-1] + Matriz[s.sequencia[i+tam_i-1]-1][s.sequencia[j+tam_j]-1] + Matriz[s.sequencia[i-1]-1][s.sequencia[j]-1] + Matriz[s.sequencia[j+tam_j-1]-1][s.sequencia[i+tam_i]-1];
        }
    }
    s.valorObj = s.valorObj + delta;
    
    //faz a troca dos segmentos
    int cont = 0;
    if(i < j){
        do{
            s.sequencia.insert(s.sequencia.begin() + j, s.sequencia[i]);
            s.sequencia.erase(s.sequencia.begin() + i);
            cont++;
        }while(cont < tam_i);
        
        cont = 0;
        do{
            s.sequencia.insert(s.sequencia.begin() + i + cont, s.sequencia[j + cont]);
            s.sequencia.erase(s.sequencia.begin() + j + cont + 1);
            cont++;
        }while(cont < tam_j);
        
    }
    else{
        do{
            s.sequencia.insert(s.sequencia.begin() + i, s.sequencia[j]);
            s.sequencia.erase(s.sequencia.begin() + j);
            cont++;
        }while(cont < tam_j);
        
        cont = 0;
        do{
            s.sequencia.insert(s.sequencia.begin() + j + cont, s.sequencia[i + cont]);
            s.sequencia.erase(s.sequencia.begin() + i + cont + 1);
            cont++;
        }while(cont < tam_i);
    }
    
   return s; 
}
*/
/*
Solucao ILS(int maxIter, int maxIterIls, vector<vector<double>>& Matriz){
    Solucao bestOfAll;
    bestOfAll.valorObj = INFINITY;
    int vertices = Matriz.size();
    for(int i = 0; i < maxIter;i++){
        Solucao s = Construcao(vertices, Matriz);
        CalculaValorObj(& s, Matriz);
        Solucao best = s;
        int iterIls = 0;
        
        while(iterIls <= maxIterIls){
            
            BuscaLocal(&s, Matriz);
            if(s.valorObj < best.valorObj){
                best = s;
                iterIls = 0;
            }
            
            s = Perturbacao(best, Matriz);
            iterIls++;
        }
    
        if(best.valorObj < bestOfAll.valorObj)
            bestOfAll = best;
    }
    return bestOfAll;   
}
*/
/*
int main(int argc, char** argv) {
    double custo = 0;
    double tempo = 0;
    srand(time(NULL));
    for(int l = 0; l < 10;l++){
        clock_t start, end;
        start = clock();
        auto data = Data(argc, argv[1]);
        data.read();
        int n = data.getDimension();
        vector<vector<double>> Matriz(n, vector<double>(n));
    
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                Matriz[i][j] = data.getDistance(i+1,j+1);
            }
        }
   
        int maxIter = 50;
        int maxIterIls;
        if(n >= 150)
            maxIterIls = n/2;
        else
            maxIterIls = n;
        Solucao best = ILS(maxIter, maxIterIls, Matriz);
        CalculaValorObj(& best, Matriz);
        custo += best.valorObj;
        end = clock();
        double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
        tempo += time_taken;
    }
    
    double custo_medio = custo/10;
    double tempo_medio = tempo/10;
    cout << "custo: " << custo_medio << endl;
    cout << "media de tempo gasto: " << fixed << tempo_medio;
    cout << " secs" << endl;
    
    
    return 0;

}
 */