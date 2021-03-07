#include <stdlib.h>
#include <stdio.h>
#include <time.h>

FILE * per_strength;

#define L 150     /* Linear dimension */
#define N (L*L)
#define EMPTY (-N-1) 
#define en 1000

int ptr[N];          /* Array of pointers */
int nn[N][4];        /* Nearest neighbors */
int order[N];        /* Occupation order */
int q = 0,a,b,check_v,check_h;   
double ensemble[en][N];

void boundaries()            // Neighbours - without periodic boundary condition
{
  int i;

  for (i=0; i<N; i++) 
  {
    if((i+1)%L == 0){nn[i][0] = EMPTY;}          //Right
    else {nn[i][0] = (i+1);} 
                   
    if(i%L == 0 || i ==0){nn[i][1] = EMPTY;}      //Left
    else {nn[i][1] = (i-1);}
                   
    
    if(i >= (L-1)*L){nn[i][2] = EMPTY;}
    else{nn[i][2] = (i+L);}                        //Down
   
   
    if(i < L){nn[i][3] = EMPTY;}
    else{nn[i][3] = (i-L);}                        //Up
  }
}


void permutation()
{
  int i,j,seed = 12345;
  int temp;
  //seed = time(NULL);
  //srand(seed);
  double r;
  
  for (i=0; i<N; i++) {order[i] = i;}
  
  for (i=0; i<N; i++) {
    r = ((double) rand() / (RAND_MAX));
    //printf("%lf\n",r);
    
    j = i + (N-i)*r;
    temp = order[i];
    order[i] = order[j];
    order[j] = temp;
  }
  
                               
}




int findroot(int i)
{
  
  if (ptr[i]<0) return i;
  return ptr[i] = findroot(ptr[i]);
}


void percolate()
{
  int i,j,z,s1,s2,r1,r2,k,q,a,b;
  double avg[N],sum[N],sm;
  per_strength = fopen("per_strength_L_400_E_1000_v_0.txt","w");

  for(k = 0; k < en; k++)
  {
    double str = 0;
    int big =0;
    if(k%50==0) printf("Number of iteration : %d\n",k+1);
    
	permutation();
    for (i=0; i<N; i++) ptr[i] = EMPTY;
    for (i=0; i<N; i++) 
      {

        check_v = 0;
        check_h = 0;
       
        r1 = s1 = order[i];
        ptr[s1] = -1;
        for (j=0; j<4; j++) 
          {
            s2 = nn[s1][j];
            if(s2 >= 0){
              if (ptr[s2]!=EMPTY) 
                {
                  r2 = findroot(s2);
                  if (r2!=r1) 
                    {
                    if(ptr[r1]>ptr[r2]) 
                    {
                      ptr[r2] +=  ptr[r1];
                      ptr[r1] = r2;
                      r1 = r2;
                    }    
                  else 
                    {
                      ptr[r1] += ptr[r2];
                      ptr[r2] = r1;
                    }
                if (-ptr[r1]>big) big = -ptr[r1];
                }
            
            }
          }
        }
      
      /*
    for(q = 0; q < N; q++)
    {
      if(ptr[q] >= 0){ptr[q] = findroot(q);}
    }
    */
  /*for(q = 0; q < N; q++)      
    {
      if(ptr[q] == EMPTY)
        {printf("x\t");}
      else{printf("%d\t",ptr[q]);}
      if((q+1)%L == 0 && q>0)
      {printf("\n");}
    }
  */  
   // printf("\n\n");

    for(a = 0; a < L; a++)      
    {
      for(b = N-L; b < N; b++)
      {
        if(ptr[a] > 0 && ptr[a] == ptr[b] )
          {
            check_v++;
            str=(-1)*ptr[findroot(ptr[a])];
          }

      }
    }

    for(a = 0; a < N; a += L)      
        {
          for(b = L-1; b < N; b += L)
            {
              if(ptr[a] > 0 && ptr[a] == ptr[b])
                {
                  check_h++;
                   str=(-1)*ptr[findroot(ptr[a])];
                }

            }
        }
    if(check_v != 0 || check_h != 0)
      {
        //str = double(big);
        ensemble[k][i] = str;
      }
      //ensemble[k][i] = str;
    }
  }

for(a = 0; a < N; a++)
{
  sm = 0.0;
  for(b = 0; b < en; b++)
  {
    sm += ensemble[b][a];
  }
sum[a] = sm;
}
for(a = 0; a < N; a++)
{
  avg[a] = sum[a]/en;
  fprintf(per_strength,"%lf\t  %lf           %lf\n",double(a+1)/N,avg[a]/N,avg[a]/(a+1) );
}

  fclose(per_strength);
}


main()
{
  boundaries();
  percolate();
 
}
