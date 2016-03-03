/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: Alfradique
 *
 * Created on 27 de Fevereiro de 2016, 21:51
 */

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
using namespace std;



double f1(double x1, double x2){
	
	return pow(log(x1), 2) + pow(log(x2), 2);
	
}

vector<double> gradf1(double x1, double x2){
	
    vector<double> ret(2);
	ret[0] = (2.0/x1);
	ret[1] = (2.0/x2);
	
	return ret;
		
}


double f2(double x1, double x2){
	
	return pow(pow(log(x1), 2) + pow(log(x2), 2), 0.5);
	
}


vector<double> gradf2(double x1, double x2){
	
	vector<double> ret(2);
	ret[0] = (1)/(x1*(pow(pow(log(x1), 2) + pow(log(x2), 2), 0.5)));
	ret[1] = (1)/(x2*(pow(pow(log(x1), 2) + pow(log(x2), 2), 0.5)));
	
	return ret;
		
}

bool isZeroGradF1(double x1, double x2){
    
    vector<double> ret(2);
    ret = gradf1(x1, x2);
    
    return ret[0] == ret[1] == 0.0;
    
}

bool isZeroGradF2(double x1, double x2){
    
    vector<double> ret(2);
    ret = gradf2(x1, x2);
    
    
    return ret[0] == ret[1] == 0.0;
    
}

vector<double> menosGrad1(double x1, double x2){
    
    vector<double> ret(2);
    ret = gradf1(x1, x2);
    
    ret[0] = -ret[0];
    ret[1] = -ret[1];
    
    return ret;
    
}

vector<double> menosGrad2(double x1, double x2){
    
    vector<double> ret(2);
    ret = gradf2(x1, x2);
    
    ret[0] = -ret[0];
    ret[1] = -ret[1];
    
    return ret;
    
}

vector<double> somaVetor(vector<double> a, vector<double> b){
    
    if(a.size() == b.size()){
        
    vector<double> ret(a.size());
    for(int i = 0; i< a.size() ; i++){
        
        ret[i] = a[i] + b[i];
        
    }
    
    return ret;
    
    } else{
        
        return vector<double>(a.size());
    }
    
}

vector<double> multiplicaConst(vector<double> a, double tk){
    
    vector<double> ret(a.size());
    for(int i = 0; i< a.size(); i++){
        
        ret[i] = tk*a[i];
        
    }
    
    return ret;
    
}


void metodoGradiente(){
    
    
    int k = 0;
    
    vector<double> xk(2);
    
    vector<double> dk(2);
    
    xk[0] = 2.0;
    xk[1] = 3.0;
    
    while(isZeroGradF1(xk[0], xk[1])){
        
        dk = menosGrad1(xk[0], xk[1]);
        vector<double> newvector(2);
        
        double tk = 0.5;
        newvector = multiplicaConst(dk,tk);
        newvector = somaVetor(newvector, xk);
        
        
        while(f1(newvector[0], newvector[1]) - f1(xk[0],xk[1]) <=0.0){
        
            tk = tk/2.0;
        }
        xk = somaVetor(xk, multiplicaConst(dk, tk));
     
    }
}

double multVetor(vector<double> a, vector<double> b){
    
    double ret;
    
    for(int i = 0; i<a.size(); i++ ){
        
        ret += a[i]*b[i];       
        
    }
    
    return ret;
}

vector<double> multiMatrizVetor(vector<vector<double> > matriz, vector<double> vetor){
    
    vector<double> ret(vetor.size());
    
    if(matriz[0].size() == vetor.size()){
        
    
    for(int i = 0; i < vetor.size(); i++){
        
        ret[i] = multVetor(matriz[i], vetor);
    }
    
    } else{
        
        cout << "multiMatrizVetor======Matriz e vetor tamanhos errados " << endl;
    }
    
    return ret;
}

vector< vector<double> > matrizInversa22(vector< vector<double> > matriz){
    
    vector< vector<double> > inversa(matriz[0].size());
    
    for(int i = 0; i< matriz[0].size; i++){
        
        inversa[i].resize(matriz[i].size());
    }
    
    double det = (matriz[0][0]*matriz[1][1] - matriz[0][1]*matriz[1][0]);    
    
    inversa[0][0] = 1/det;
    inversa[0][1] = -1/det;
    inversa[1][0] = -1/det;
    inversa[1][1] = 1/det;
    
    return inversa;
}

vector< <vector<double> > Hessianf1(double x1, double x2){
    
    vector< <vector<double> > ret(2, vector<double>(2));
    
    ret[0][0] = (-2*(log(x1) - 2)*log(x1))/x1*x1
    ret[1][1] = (-2*(log(x2) - 2)*log(x2))/x2*x2
    ret[0][1] = ret[1][0] = 0.0;
    
    return ret;
    
}


void metodoDeNewton(){
    
    int k = 0;
    
    vector<double> xk(2);
    
    vector<double> dk(2);
    
    xk[0] = 2.0;
    xk[1] = 3.0;
    while(isZeroGradF1(xk[0], xk[1])){
        
        dk = multiMatrizVetor(matrizInversa22(Hessianf1(xk[0], xk[1])),gradf1(xk[0], xk[1]));
        
        xk = somaVetor(xk, dk);
        
    }
    
    
}

/*
 * 
 */
int main(int argc, char** argv) {

    
    
    
    
    
    
    return 0;
}

