#include <iostream>
#include <math.h>
#include <random>
#include <time.h>
#include <chrono>
#include <fstream>
#include <unistd.h>
#include "parameters.h"

using namespace std::chrono;
using namespace std;

double Rand() {
    thread_local std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<double> dis(0, 1);
    return dis(gen);
}

int mod(int a,int b){
    return (a%b + b)%b;
}
void shiftx(int b[], int a[], int sign, int v){
    for (int i=0; i<d; i++){
        if (i==v){
            if (v==0){ b[i]=mod(a[i]+sign,Nt);  }
            else { b[i]=mod(a[i]+sign,Ns);  }
        }
        else { b[i]=a[i]; }
    }
}
void random(int a[]){
    for (int i=1; i<d; i++){
        a[i]=rand()%Ns;
    }
    a[0]=rand()%Nt;
}

long double I(int s){
    long double a=0,r=0;
    while(r<inf){
        a+=dr*pow(r,s+1)*exp(-eta*pow(r,2)-lmd*pow(r,4));
        r+=dr;
    }
    return a;
}

bool eql(int x1[],int x2[]){
    bool flag=true;
    for (int i=0;i<d;i++){
        if(x1[i]!=x2[i]){
            flag=false;
            break;
        }
    }
    return flag;
}

int ksum(int k[][Ns][d]){
    int sum=0;
    for (int i=0;i<Nt;i++){
        for (int j=0;j<Ns;j++){
            sum+=k[i][j][0];
        }
    }
    return sum;
}
long int sx(int x[], int k[][Ns][d], int a[][Ns][d]){
    int sum=0;
    int v[d]={0};
    for (int i=0;i<d;i++){
        v[i]=1;

        sum+=abs(k[x[0]][x[1]][i])
        +abs(k[mod(x[0]-v[0],Nt)][mod(x[1]-v[1],Ns)][i])
        +2*(a[x[0]][x[1]][i]
        +a[mod(x[0]-v[0],Nt)][mod(x[1]-v[1],Ns)][i]);

        v[i]=0;
    }
    return sum;
}

void update(){
    int y;
    for (int i=0;i<Nt;i++){
        for (int j=0;j<Ns;j++){
            for (int t=0;t<d;t++){
                y=a_[i][j][t];
                a_[i][j][t]=a[i][j][t]+2*(rand()%2)-1;
                
                if (a_[i][j][t]<0){
                    a_[i][j][t]=y;
                    continue; 
                }
                
                int temp[d], temp1[d];
                temp[0]=i; temp[1]=j;
                shiftx(temp1, temp, 1, t);
                if (a_[i][j][t]>a[i][j][t]){
                    rho=1.0/(abs(k[i][j][t])+a_[i][j][t])
                    /a_[i][j][t]
                    *I_val[sx(temp,k,a_)]*I_val[sx(temp1,k,a_)]
                    /I_val[sx(temp,k,a)]/I_val[sx(temp1,k,a)];
                    
                }
                else {
                    rho=1.0*(abs(k[i][j][t])+a[i][j][t])
                    *a[i][j][t]
                    *I_val[sx(temp,k,a_)]*I_val[sx(temp1,k,a_)]
                    /I_val[sx(temp,k,a)]/I_val[sx(temp1,k,a)];
                }
                if (Rand()<rho){
                    a[i][j][t]=a_[i][j][t];
                }
                else{
                    a_[i][j][t]=y;
                }
            }
        }
    }
    
    bool start=false;
    while (!start){
        del=2*(rand()%2)-1;
        v=rand()%d;
        sign=2*(rand()%2)-1;
        random(x0);
        shiftx(x,x0,sign,v);
        if (sign<0){
            shiftx(xx,x0,-1,v);
        }
        else{
            shiftx(xx,x0,0,v);
        }
        y=k_[xx[0]][xx[1]][v];
        k_[xx[0]][xx[1]][v]=k[xx[0]][xx[1]][v]+del*sign;
        
        if(abs(k_[xx[0]][xx[1]][v])>abs(k[xx[0]][xx[1]][v])){
            rho=exp(sign*mu*del*(v==0))
            /(abs(k_[xx[0]][xx[1]][v])+a[xx[0]][xx[1]][v])
            *A/I_val[sx(x,k,a)]/I_val[sx(x0,k,a)];
            
        }
        else {
            rho=exp(sign*mu*del*(v==0))
            *(abs(k[xx[0]][xx[1]][v])+a[xx[0]][xx[1]][v])
            *A/I_val[sx(x,k,a)]/I_val[sx(x0,k,a)];
        }
        if (Rand()<rho){
            k[xx[0]][xx[1]][v]=k_[xx[0]][xx[1]][v];
            start=true;
        }
        else{
            k_[xx[0]][xx[1]][v]=y;
        }
    }
    
    while (!eql(x,x0)){
        v=rand()%d;
        sign=2*(rand()%2)-1;
        shiftx(x_,x,sign,v);
        
        if (sign<0){
            shiftx(xx,x,-1,v);
        }
        else{
            shiftx(xx,x,0,v);
        }
        y=k_[xx[0]][xx[1]][v];
        k_[xx[0]][xx[1]][v]=k[xx[0]][xx[1]][v]+del*sign;
        
        if(abs(k_[xx[0]][xx[1]][v])>abs(k[xx[0]][xx[1]][v])){
            if (eql(x_,x0)){
                rho=exp(sign*mu*del*(v==0))
                /(abs(k_[xx[0]][xx[1]][v])+a[xx[0]][xx[1]][v])
                *I_val[sx(x,k_,a)]*I_val[sx(x_,k_,a)]/A;
            }
            else {
                rho=exp(sign*mu*del*(v==0))
                /(abs(k_[xx[0]][xx[1]][v])+a[xx[0]][xx[1]][v])
                *I_val[sx(x,k_,a)]/I_val[sx(x_,k,a)];
            }
        }
        else {
            if (eql(x_,x0)){
                rho=exp(sign*mu*del*(v==0))
                *(abs(k[xx[0]][xx[1]][v])+a[xx[0]][xx[1]][v])
                *I_val[sx(x,k_,a)]*I_val[sx(x_,k_,a)]/A;
            }
            else{
                rho=exp(sign*mu*del*(v==0))
                *(abs(k[xx[0]][xx[1]][v])+a[xx[0]][xx[1]][v])
                *I_val[sx(x,k_,a)]/I_val[sx(x_,k,a)];
            }
        }
        
        if (Rand()<rho){
            k[xx[0]][xx[1]][v]=k_[xx[0]][xx[1]][v];
            shiftx(x,x_,0,v);
        }
        else{
            k_[xx[0]][xx[1]][v]=y;
        }
    }
    
}
float phi2(int k[][Ns][d], int a[][Ns][d], long double I[int_val]){
    long double sum=0;
    int x[d];
    for (int i=0;i<Nt;i++){
        for (int j=0;j<Ns;j++){
            x[0]=i; x[1]=j;
            sum+=I[sx(x, k, a)+2]/I[sx(x, k, a)];
        }
    }
    return sum/Nt/Ns;
}
float phi4(int k[][Ns][d], int a[][Ns][d], long double I[int_val]){
    long double sum=0;
    int x[d];
    for (int i=0;i<Nt;i++){
        for (int j=0;j<Ns;j++){
            x[0]=i; x[1]=j;
            sum+=I[sx(x, k, a)+4]/I[sx(x, k, a)];
        }
    }
    return sum/Nt/Ns;
}
float errorjack(float xi[configs]){
    float x_i[configs], x_=0, stddev=0;
    for (int i=0; i<configs; i++){
        for (int j=0; j<configs; j++){
            x_i[i]=(1-(i==j))*xi[j];
        }
        x_i[i]=x_i[i]/(configs-1);
        x_+=x_i[i];
    }
    x_=x_/configs;
    for (int i=0; i<configs; i++){
        stddev+=(x_i[i]-x_)*(x_i[i]-x_);
    }
    stddev=sqrt(stddev*(configs-1)/configs);
    return stddev;
}
int main(int argc, char **argv)
{
    auto begin=high_resolution_clock::now();
    srand(time(NULL));
    for (int i=0; i<int_val; i++){
        I_val[i]=I(i);
    }
    float phi2_avg, phi4_avg, dmu=(mu_max-mu_min)/mu_n;;
    float xi[configs],phi2i[configs];
    mu=mu_min;
    
    ofstream data, data1, data2;
    data.open("mu_vs_n.txt");
    data1.open("mu_vs_phi2.txt");
    data2.open("mu_vs_phi4.txt");
    for (int g=0; g<mu_n; g++){
        
        for (int i=0; i<equil; i++){
            update();
            
        }
        phi2_avg=0;
        n_avg=0;
        phi4_avg=0;
        for (int i=0; i<configs; i++){
            for (int j=0; j<gaps; j++){
                update();
                
            }
            update();
            
            xi[i]=1.0*ksum(k)/Nt/Ns;
            phi2i[i]=phi2(k,a,I_val);
            
            n_avg+=xi[i];
            phi2_avg+=phi2i[i];
            phi4_avg+=phi4(k,a,I_val);
        }
        n_avg=n_avg/configs;
        phi2_avg=phi2_avg/configs;
        phi4_avg=phi4_avg/configs;
        
        data<<mu<<"\t"<<n_avg<<"\t"<<errorjack(xi)<<endl;
        data1<<mu<<"\t"<<phi2_avg<<"\t"<<errorjack(phi2i)<<endl;
        data2<<mu<<"\t"<<phi4_avg<<endl;
        mu+=dmu;
        cout<<g<<endl;
    }
    
    data.close();
    data1.close();
    data2.close();
    auto stop=high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop-begin);
    cout<<duration.count()<<endl;
    
    return 0;
}