#include<iostream>
#include<iomanip>
using namespace std;
float cost[20];
float table[10][20];
int xb[10];
float cb[10];
float b[10];
float min_ratio[10];
float delta[20];
int m,n;
float solution;

void print_table()
{
    cout<<"Table: \n";
    int i,j;
    cout<<"xb  ";
    cout<<"cb    ";
    for(j=0;j<n+m;j++){
        cout<<"x"<<j+1<<"   ";
    }
    cout<<"b"<<"\n";

    for(i=0;i<m;i++){ 
        cout<<"x"<<xb[i]+1<<" ";
        cout<<cb[i]<<"  ";
 	    for(j=0;j<n+m;j++){
    	   cout<<table[i][j]<<" ";
        }
        cout<<" "<<b[i]<<"\n";
    }
    cout<<"         ";    
    for(j=0;j<m+n;j++)
        cout<<delta[j]<<" "; 
    cout<<"\n";     
}



int min(float a[],int n)
{
    int position=0;float min=a[0];
    for(int i=1;i<n;i++){
    	if(a[i]<min){
    		min=a[i];
    		position=i;
    	}
    }
    return position;
}


void gauss_elimination(int p,int q)
{
    int i,j;
    for(j=0;j<m+n;j++){
        if(j!=q){ table[p][j]=table[p][j]/table[p][q]; }     
    }
    b[p]=b[p]/table[p][q];
    table[p][q]=1;

    for(i=0;i<m;i++){
        if(i!=p){
            for(j=0;j<m+n;j++){
            	if(j!=q){
            		table[i][j]=table[i][j]-table[i][q]*table[p][j];
            	}
            }
            b[i]=b[i]-table[i][q]*b[p];
            table[i][q]=0;
        }
    }
}

void simplex()
{
    int i,j,pivot_i,pivot_j,iteration=0;
    solution=0;
    while(1)
    { 
    	//calculating delta
        for(j=0;j<m+n;j++){
        	delta[j]=0;
            for(i=0;i<m;i++){
            	delta[j]+=table[i][j]*cb[i];
            }
            delta[j]=delta[j]-cost[j];
        }
    	cout<<"\n"<<iteration<<" ";
    	print_table();

    	//selecting pivot column
        pivot_j=min(delta,m+n);
        cout<<"minimum delta "<<delta[pivot_j]<<"\n";

        if(delta[pivot_j]>=0){
            cout<<"========== Optimality reached ==========\n";
            print_table();
            cout<<"Solution :\n";
            for(i=0;i<m;i++){
                cout<<"x"<<xb[i]+1<<" = "<<b[i]<<"\n";
                solution+= cost[xb[i]] *b[i];
            }
            cout<<"Answer : "<<solution<<"\n";
            break;
        }
        //calculating minimim ratios
        for(i=0;i<m;i++){ 
            if(table[i][pivot_j]>0)
                min_ratio[i]=b[i]/table[i][pivot_j];
            else
                min_ratio[i]=1000000;
        }
        //selecting pivot row
        pivot_i=min(min_ratio,m);
        cout<<"minimum ration at "<<b[pivot_i]<<"\n";

        if(min_ratio[pivot_i] == 1000000.0){      //when all element in
            cout<<"____infeasibe solution____\n"; //pivot row are -ve
            print_table();
            break;
        }
        cout<<"x"<<pivot_j+1<<" enters "<<"x"<<xb[pivot_i]+1<<" leaves \n";


        //proceeding to next iteration
        xb[pivot_i]=pivot_j;       //updating xb
        cb[pivot_i]=cost[pivot_j]; //updating cb
        gauss_elimination(pivot_i,pivot_j);

        //condition to break loop if iteration > 10
        if(iteration>10){
        	cout<<"iteration limit crossed\n";
        	break;
        }
        iteration++;
    }
}

int main()
{
 std::cout << std::fixed;
 std::cout << std::setprecision(2);  

 cout<<"========== Simplex Method ==========\n";
 cout<<"for Maximization and <= type equations\n"; 
 int i,j;
 cout<<"Enter no. of equations: ";
 cin>>m;
 cout<<"Enter no. of variables: ";
 cin>>n;
 cout<<"Enter space seperated coefficient of equations including 'b's\n";
 for(i=0;i<m;i++)
 {
   	 for(j=0;j<n;j++)
     {
    	cin>>table[i][j];
     }
     cin>>b[i];
     table[i][n+i]=1;
 }
 cout<<"Enter cost:\n";
 for(i=0;i<n;i++)
	 cin>>cost[i];

 
//preparing table for first iteration    
 for(i=0;i<m;i++){//filling xb with slack variable's number  
 	xb[i]=n+i;     
 }


 simplex();
}







