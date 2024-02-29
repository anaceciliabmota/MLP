#include "Data.h"
#include <iostream>
#include <vector>
#include <time.h> 
#include <cstdlib>
#include <algorithm>
#include <iomanip>

using namespace std;

/*bestImprovementOrOpt(&s, subseq_matrixs, 2, data); e  bestImprovementOrOpt(&s, subseq_matrixs, 3, data) não calcula o custo corretamente.
Consertar.
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
    s->sequence.insert(s->sequence.end() -1, custoInsercao[c].noInserido);
    /*for(int i = 0; i < s->sequence.size(); i++){
        if (s->sequence[i] == r){
            s->sequence.insert(s->sequence.begin() + i + 1, custoInsercao[c].noInserido);
            break;
        }
    }*/

}

Solucao Construcao(Data& data, double a){
    Solucao s = {{1,1}, 0.0};
    size_t n = data.getDimension();
    vector<int> CL = nosRestantes(n);
    CL.erase(CL.begin());
    int r = 1;
    while (CL.size() != 0){
        vector <InsertionInfo> custoInsercao = CalcularCusto(r, CL, data);
        sort(custoInsercao.begin(), custoInsercao.end(), myfunction);
        int c = rand() % ((int) ceil((a + 0.000001) * custoInsercao.size()));
        inserirNaSolucao(&s, c, r, custoInsercao);
        r = custoInsercao[c].noInserido;
        for (int i = 0; i < CL.size(); i++){
            if(CL[i] == custoInsercao[c].noInserido){
                CL.erase(CL.begin() + i);
                break;
            }
        }
    }
    CalculaValorObj(data, &s);
    return s;
}
bool bestImprovementSwap(Solucao *s,vector<vector<Subsequence>>& subseq_matrix, Data& data){
    double bestCost = s->valorObj;
    
    int best_i,best_j;
    for(int i = 1; i < s->sequence.size() - 1; i++){
   
        for(int j = i+1; j < s->sequence.size() - 1;j++){
            double cost;
            if((i + 1)!= j){
                Subsequence sigma_1 = Subsequence :: Concatenate(subseq_matrix[0][i-1], subseq_matrix[j][j], data);
                Subsequence sigma_2 = Subsequence :: Concatenate(sigma_1, subseq_matrix[i+1][j-1], data);
                Subsequence sigma_3 = Subsequence :: Concatenate(sigma_2, subseq_matrix[i][i], data);
                Subsequence sigma_4 = Subsequence :: Concatenate(sigma_3, subseq_matrix[j+1][s->sequence.size()-1], data);
                
                cost = sigma_4.C;
                
            } else {
                
                Subsequence sigma_1 = Subsequence :: Concatenate(subseq_matrix[0][i-1],subseq_matrix[j][i], data);
                Subsequence sigma_2 = Subsequence :: Concatenate(sigma_1, subseq_matrix[j+1][s->sequence.size()-1], data);
                
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
        //UpdateAllSubseq(s, subseq_matrix, data);
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
            // j = ultimo elemento da aresta a ser invertida
            Subsequence sigma_1;
            sigma_1 = Subsequence :: Concatenate(subseq_matrix[0][i-1], subseq_matrix[j][i], data);
            Subsequence sigma_2 = Subsequence :: Concatenate(sigma_1, subseq_matrix[j+1][s->sequence.size()-1], data);
            if(sigma_2.C < bestCost){
                bestCost = sigma_2.C;
                best_i = i;
                best_j = j;
            }
        }
        
    }
    if(bestCost < s->valorObj){
        int cont = 0;
        //cout << "best i: " << best_i << "best j" << best_j;
        while((best_i + cont) < (best_j - cont)){
            swap(s->sequence[best_i + cont], s->sequence[best_j - cont]);
            cont++;
        }
        s->valorObj = bestCost;
        //UpdateAllSubseq(s, subseq_matrix, data);
        return true;
    }
    return false;
}
bool bestImprovementOrOpt(Solucao *s, vector<vector<Subsequence>>& subseq_matrix, int tam_bloco, Data& data){
    double bestCost = s->valorObj;
    int best_i1;
    int best_i2;
    int best_j;
    for(int i = 1; i < s->sequence.size() - 1;i++){
        /*int vi1 = s->sequencia[i];
        int vi1_prev = s->sequencia[i-1];
        int vi2 = s->sequencia[i + tam_bloco - 1];
        int vi2_next = s->sequencia[i + tam_bloco];*/
        int i2 = i + tam_bloco -1;
        for(int j = i + tam_bloco; j < s->sequence.size() - 1; j++){
            //int vj = s->sequencia[j];
            //int vj_next = s->sequencia[j+1];
            Subsequence sigma_1 = Subsequence :: Concatenate(subseq_matrix[0][i-1], subseq_matrix[i2+1][j], data);
            Subsequence sigma_2 = Subsequence :: Concatenate(sigma_1, subseq_matrix[i][i2], data);
            Subsequence sigma_3 = Subsequence :: Concatenate(sigma_2, subseq_matrix[j+1][s->sequence.size()-1], data);
            //double delta = -data.getDistance(vi1_prev, vi1) - data.getDistance(vi2, vi2_next) - data.getDistance(vj, vj_next) + data.getDistance(vi1_prev, vi2_next) + data.getDistance(vj, vi2) + data.getDistance(vi1, vj_next);
            //double delta = -Matriz[vi1_prev-1][vi1-1] - Matriz[vi2-1][vi2_next-1]- Matriz[vj-1][vj_next-1] + Matriz[vi1_prev-1][vi2_next-1] + Matriz[vj-1][vi2-1] + Matriz[vi1-1][vj_next-1];
            if (sigma_3.C < bestCost){
                bestCost = sigma_3.C;
                best_i1 = i;
                best_i2 = i+tam_bloco -1;
                best_j = j;
            }
        }
    }

    int cont = 0;
    if(bestCost < s->valorObj){
        do{
        s->sequence.insert(s->sequence.begin() + best_j + 1, s->sequence[best_i1]);
        s->sequence.erase(s->sequence.begin() + best_i1);
        best_j--;
        cont++;
        }while(cont < tam_bloco);
        s->valorObj = bestCost;
        //UpdateAllSubseq(s, subseq_matrix, data);
        
        return true;
    }
    return false;
}
void BuscaLocal(Solucao *s, vector<vector<Subsequence>>& subseq_matrix, Data& data){
    vector<int> NL = {1, 2, 3, 4, 5};
    bool improved = false;

    while (NL.empty() == false){
        int n = rand() % NL.size();
        switch (NL[n]){
            case 1:
                improved = bestImprovementSwap(s, subseq_matrix, data);
                if(improved){UpdateAllSubseq(s, subseq_matrix, data);}
                /*cout << "caso swap: " << s->valorObj << " = ";
                CalculaValorObj(data, s);
                cout << s->valorObj << endl;*/
                break;
            case 2:
                improved = bestImprovement2Opt(s, subseq_matrix, data);
                if(improved){UpdateAllSubseq(s, subseq_matrix, data);}
                /*cout << "caso 2Opt: " << s->valorObj << " = ";
                CalculaValorObj(data, s);
                cout << s->valorObj << endl;*/
                break;
            case 3:
                improved = bestImprovementOrOpt(s, subseq_matrix, 1, data);
                if(improved){UpdateAllSubseq(s, subseq_matrix, data);}
                /*cout << "caso OrOpt1: " << s->valorObj << " = ";
                CalculaValorObj(data, s);
                cout << s->valorObj << endl;*/
                break;
            case 4: 
                improved = bestImprovementOrOpt(s, subseq_matrix, 2, data);
                if(improved){UpdateAllSubseq(s, subseq_matrix, data);}
                /*cout << "caso OrOpt2: " << s->valorObj << " = ";
                CalculaValorObj(data, s);
                cout << s->valorObj << endl;*/
                break;
            case 5:
                improved = bestImprovementOrOpt(s, subseq_matrix, 3, data);
                if(improved){UpdateAllSubseq(s, subseq_matrix, data);}
                /* cout << "caso OrOpt3: " << s->valorObj << " = ";
                CalculaValorObj(data, s);
                cout << s->valorObj << endl;*/
                break;
        }
        if(improved){
            NL = {1,2,3,4,5};
        }else{
            NL.erase(NL.begin() + n);
        }
    }
}
Solucao Perturbacao (Solucao best, vector<vector<Subsequence>>& subseq_matrix, Data& data){

    Solucao s = best;
    int vertices = best.sequence.size() - 1;
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
    
    double cost;
    if((i+ tam_i) == j || (j +tam_j) == i){
        if( i < j){
            Subsequence sigma_1 = Subsequence :: Concatenate(subseq_matrix[0][i-1], subseq_matrix[j][j+tam_j-1], data);
            Subsequence sigma_2 = Subsequence :: Concatenate(sigma_1, subseq_matrix[i][i+tam_i-1], data);
            Subsequence sigma_3 = Subsequence :: Concatenate(sigma_2, subseq_matrix[j+tam_j][s.sequence.size()-1], data);
            cost = sigma_3.C;
        }
        else{
            Subsequence sigma_1 = Subsequence :: Concatenate(subseq_matrix[0][j-1], subseq_matrix[i][i+tam_i-1], data);
            Subsequence sigma_2 = Subsequence :: Concatenate(sigma_1, subseq_matrix[j][j+tam_j-1], data);
            Subsequence sigma_3 = Subsequence :: Concatenate(sigma_2, subseq_matrix[i+tam_i][s.sequence.size()-1], data);
            cost = sigma_3.C;
        }
    }else{
        if(i < j){
            Subsequence sigma_1 = Subsequence :: Concatenate(subseq_matrix[0][i-1], subseq_matrix[j][j + tam_j -1], data);
            Subsequence sigma_2 = Subsequence :: Concatenate(sigma_1, subseq_matrix[i+tam_i][j-1], data);
            Subsequence sigma_3 = Subsequence :: Concatenate(sigma_2, subseq_matrix[i][i+tam_i-1], data);
            Subsequence sigma_4 = Subsequence :: Concatenate(sigma_3, subseq_matrix[j+tam_j][s.sequence.size()-1], data);
            cost = sigma_4.C;
        }
        else{
            Subsequence sigma_1 = Subsequence :: Concatenate(subseq_matrix[0][j-1], subseq_matrix[i][i + tam_i -1], data);
            Subsequence sigma_2 = Subsequence :: Concatenate(sigma_1, subseq_matrix[j+tam_j][i-1], data);
            Subsequence sigma_3 = Subsequence :: Concatenate(sigma_2, subseq_matrix[j][j+tam_j-1], data);
            Subsequence sigma_4 = Subsequence :: Concatenate(sigma_3, subseq_matrix[i+tam_i][s.sequence.size()-1], data);
            cost = sigma_4.C;
        }
    }
    s.valorObj = cost;
    
    //faz a troca dos segmentos
    int cont = 0;
    if(i < j){
        do{
            s.sequence.insert(s.sequence.begin() + j, s.sequence[i]);
            s.sequence.erase(s.sequence.begin() + i);
            cont++;
        }while(cont < tam_i);
        
        cont = 0;
        do{
            s.sequence.insert(s.sequence.begin() + i + cont, s.sequence[j + cont]);
            s.sequence.erase(s.sequence.begin() + j + cont + 1);
            cont++;
        }while(cont < tam_j);
        
    }
    else{
        do{
            s.sequence.insert(s.sequence.begin() + i, s.sequence[j]);
            s.sequence.erase(s.sequence.begin() + j);
            cont++;
        }while(cont < tam_j);
        
        cont = 0;
        do{
            s.sequence.insert(s.sequence.begin() + j + cont, s.sequence[i + cont]);
            s.sequence.erase(s.sequence.begin() + i + cont + 1);
            cont++;
        }while(cont < tam_i);
    }
   return s; 
}

Solucao GILSRVND(int maxIter, int maxIterIls, Data& data, vector<double> R){
    Solucao bestOfAll;
    bestOfAll.valorObj = INFINITY;
    for(int i = 0; i < maxIter;i++){
        int indice = rand()%26;
        double a = R[indice];
        Solucao s = Construcao(data, a);
        int n = s.sequence.size();
        vector<vector<Subsequence>> subseq_matrixs(n, vector<Subsequence>(n));
        vector<vector<Subsequence>> subseq_matrixb(n, vector<Subsequence>(n));
        UpdateAllSubseq(&s, subseq_matrixs, data);
        Solucao best = s;
        UpdateAllSubseq(&best, subseq_matrixb, data);
        int iterIls = 0;
        while(iterIls <= maxIterIls){
            BuscaLocal(&s, subseq_matrixs, data);
            if(s.valorObj < best.valorObj){
                best = s;
                UpdateAllSubseq(&best, subseq_matrixb, data);
                iterIls = 0;
            }
            s = Perturbacao(best, subseq_matrixb, data);
            /*cout << "custo: " << s.valorObj;
            CalculaValorObj(data, &s);
            cout << " = " << s.valorObj << endl; */
            UpdateAllSubseq(&s, subseq_matrixs, data);
            iterIls++;
        }
    
        if(best.valorObj < bestOfAll.valorObj)
            bestOfAll = best;

    }
    return bestOfAll;   
}
int main(int argc, char** argv) {
    double custo = 0;
    double tempo = 0;
    srand(time(NULL));
    vector<double> R = {0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.1, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25};
    for(int l = 0; l < 10;l++){
        clock_t start, end;
        start = clock();
        auto data = Data(argc, argv[1]);
        data.read();
        int n = data.getDimension();
        int maxIter = 10;
        int maxIterIls;
        if(n >= 100)
            maxIterIls = 100;
        else
            maxIterIls = n;
        Solucao best = GILSRVND(maxIter, maxIterIls, data, R);
        CalculaValorObj(data, &best);
        //ExibirSolucao(&best);
        custo += best.valorObj;
        end = clock();
        double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
        tempo += time_taken;
        cout << "Execucao " << l+1 << " - custo: " << fixed << setprecision(2) << best.valorObj << endl;
    }
    double custo_medio = custo/10;
    double tempo_medio = tempo/10;
    cout << "custo: " << fixed << setprecision(2) << custo_medio << endl;
    cout << "media de tempo gasto: " << fixed << tempo_medio;
    cout << " secs" << endl;
   //codigo para teste para descobrir onde está o erro do meu código
    /*auto data = Data(argc, argv[1]);
    data.read();
    int indice = rand()%26;
    double a = R[indice];
    Solucao s = Construcao(data, a);
    int n = s.sequence.size();
    vector<vector<Subsequence>> subseq_matrixs(n, vector<Subsequence>(n));
    UpdateAllSubseq(&s, subseq_matrixs, data);
    for(int l = 0; l < 100; l++){
        cout << "antes:" << fixed << setprecision(2) << s.valorObj<< endl;
        bool improved = bestImprovementOrOpt(&s, subseq_matrixs, 3, data);
        if(improved){UpdateAllSubseq(&s, subseq_matrixs, data);}
        cout << "depois: "<< fixed << setprecision(2) << s.valorObj << " ";
        CalculaValorObj(data, &s);
        cout << fixed << setprecision(2) << s.valorObj<< endl;
    }
    */
    return 0;

}





















