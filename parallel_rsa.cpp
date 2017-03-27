#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<omp.h>

//int wm[16]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
/*float host[6][6]={0.39,0.0274509803921569,	0.0313725490196078,	0.0274509803921569,	0.0313725490196078,	0.0196078431372549
,		  0.87,0.14,0.27,0.19,0.61,0.39,
		  0.23,0.89,0.09,0.45,0.18,0.71,
		  0.12,0.67,0.12,0.50,0.71,0.29,
		  0.87,0.14,0.27,0.19,0.61,0.39,
		  0.23,0.89,0.09,0.45,0.18,0.71,};*/
float host[50][50];

long long unsigned int mod_exp(int base,int exp,int n) 
{ 
long long unsigned int i,pow=1; 
int id;
#pragma omp parallel //for reduction(*:pow)
{	int i;
#pragma omp for reduction(*:pow)

  for(i=0;i<exp;i++) 
   pow=(pow*base)%n; 
}
  return pow; 
} 

int gcd(int m,int n) 
{ 
  while(m!=n) 
  { 
    if(m==n || n==1) 
    { 
      m=1,n=1; 
      break; 
    } 
    else 
    { 
     if(m>n) 
      m=m-n; 
     else 
      n=n-m; 
    } 
  } 
  return m; 
} 

float average(float a[2500],int n)
{
	float s=0,m;
	int i;
	for(i=0;i<2500;i++)
		s=s+a[i];
	m=s/n;
	return(m);
}

float dct_cov(float z[50][50],float u[50][50])
{

	int p,q,i,j,k,l;
	float s,t;
	float d[50][50]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

#pragma omp parallel shared(z,d) private(p,q,s,t,i,j)
	{
		//printf("\nthread id=%d\n",omp_get_thread_num()); 
#pragma omp for schedule(static)
	
		//printf("\nthread id=%d\n",omp_get_thread_num()); 
		for(p=0;p<50;p++)

		for(q=0;q<50;q++)
		{
			if(p==0)

				s=(float)sqrt(1.0/50.0);

			else

				s=(float)sqrt(2.0/50.0);
			if(q==0)

			       t=(float)sqrt(1.0/50.0);
			else
				t=(float)sqrt(2.0/50.0);


		      for(i=0;i<50;i++)

				for(j=0;j<50;j++)
				{
					//printf("\nthread id=%d\n",omp_get_thread_num());
d[p][q]+=z[i][j]*(float)cos((3.1412*(2*i+1)*p)/(2*50))*(float)cos((3.1412*(2*j+1)*q)/(2*50));
				}
			d[p][q]*=s*t;
		}
		
	}
		for(i=0;i<50;i++)
		  for(j=0;j<50;j++)
		    u[i][j]=d[i][j];
		return(u[50][50]);
}

float idct_cov(float z[50][50],float u[50][50])
{
	int p,q,i,j,k,l;
	float s,t;
	float f[50][50]={0,0,0,0,0,0,0,0,0,
		       0,0,0,0,0,0,0,0,0,
		       0,0,0,0,0,0,0,0,0,
		       0,0,0,0,0,0,0,0,0,};
	#pragma omp parallel shared(z,f) private(p,q,s,t,i,j)
	{
		printf("\nthread id=%d\n",omp_get_thread_num()); 
#pragma omp for schedule(static)
	for(i=0;i<50;i++)

		for(j=0;j<50;j++)
		{
			for(p=0;p<50;p++)
				for(q=0;q<50;q++)
				{
				if(p==0)

					s=(float)sqrt(1.0/50.0);

				else

					s=(float)sqrt(2.0/50.0);
				if(q==0)

					t=(float)sqrt(1.0/50.0);
				else
					t=(float)sqrt(2.0/50.0);

f[i][j]+=(s*t*z[p][q]*(float)cos((3.1412*(2*i+1)*p)/(2*50))*(float)cos((3.1412*(2*j+1)*q)/(2*50)));

				}

	       }
	}
	       for(i=0;i<50;i++)
		 for(j=0;j<50;j++)
		   u[i][j]=f[i][j];
  return(u[50][50]);
}

/*float idct_wm(float z[2][2],float u[2][2])
{
	int p,q,i,j,k,l;
	float s,t;
	float f[2][2]= {0};

	for(i=0;i<2;i++)

		for(j=0;j<2;j++)
		{
			for(p=0;p<2;p++)
				for(q=0;q<2;q++)
				{
				if(p==0)

					s=(float)sqrt(1.0/2.0);

				else

					s=(float)sqrt(2.0/2.0);
				if(q==0)

					t=(float)sqrt(1.0/2.0);
				else
					t=(float)sqrt(2.0/2.0);

f[i][j]+=(s*t*z[p][q]*(float)cos((3.1412*(2*i+1)*p)/(2*4))*(float)cos((3.1412*(2*j+1)*q)/(2*4)));
				}

	       }

	       for(i=0;i<2;i++)
		  for(j=0;j<2;j++)
		     u[i][j]=f[i][j];


return(u[2][2]);
}
*/





void main()
{	
	double start,end;
	int p1=19,q=13; 

	//int decr[10];
	float dct_host[50][50];
	float temp[2][2];
	float temp1[2][2];
	float c2[50][50];
	float c1[50][50];
	float dct_marked[50][50];
	float dct_mark[50][50];
	//float temp[4][4];
	//float temp1[4][4];
	float newdct[4][4];
	float marked[50][50];
	int i,j;
	float c1_1D[2500];
	int count=0;
	float mean;
	int seq[2500];
	int p[50][50];
	int mark_encr[50][50];
	float q1[50][50];
	int v[50][50],w[50][50];
 int wm_extr[50][50];
 float out1[4]={0};

 float out2[2][2];
 int counter=0;
 int c2_1d[4]={0};
	 int n=p1*q; 

  int d=1,e=n-1; 

  char in[20]=" "; 

long long unsigned int out[20]={0},decr[20];

int z; 

/*  input from file  */

FILE *fp,*fp1;
float host_1d[2510];
  float value;
  int ct = -1; /* EDIT have i start at -1 :) */
  int ar_ct=0;
  if ((fp = fopen ("50_balloon.txt", "r")) == NULL)
    printf("Error");

  while (!feof (fp) && fscanf (fp, "%f,", &value) && ct++ < 2510 )
    host_1d[ct] = value;

  fclose (fp);
	  printf("host image 1d array is \n\n");

	  //for(i=0;i<350;i++)
		  //printf("%f\t",host_1d[i]);
	  //printf("\n");

	  printf("\nhost image 2d array is :");
	for(i=0;i<50;i++)
	{
		//printf("\n:");
	 for(j=0;j<50;j++)
	 {
		 host[i][j]=host_1d[ar_ct++];
		 //printf("%f\t",host[i][j]);
    }
	}



start=omp_get_wtime();

z=((p1-1)*(q-1)); 

for(;;) 

{ 

 if(gcd(e,z)==1) 

   break; 

  e--; 

} 

for(;;) 

{ 

 if((e*d)%z==1) 

    break; 

  d++; 

  } 
printf("%d",e);

  printf("enter the input string:\n"); 

  gets(in); 

  printf("\nThe entered string is: %s",in); 

  printf("\n"); 

  printf("The encrypted string is \n"); 
  #pragma omp parallel for
      for(i=0;i<strlen(in);i++) 

{ 
	printf("\nthread no. = %d\n",omp_get_thread_num());
  out[i]=mod_exp(in[i],e,n); 

  printf("%c=%3d\n",in[i],out[i]); 

}
	  

	   for(i=0;i<strlen(in);i++) 

{ 

  out1[i]=(float)out[i]; 
  out1[i]=out1[i]/1000;

  printf("%f\n",out1[i]); 

}
	   

	   printf("\nthe 2d array is :");
	for(i=0;i<2;i++)
	{
		printf("\n:");
	 for(j=0;j<2;j++)
	 {
		 out2[i][j]=out1[count++];
		 printf("%f\t",out2[i][j]);
    }
	}

	printf("\n\n");

	dct_cov(host,dct_host);

	for(i=0;i<50;i++)
		for(j=0;j<50;j++)
			c1[i][j]=dct_host[i][j];

	printf("dct coefts of host=\n");
	for(i=0;i<50;i++)
	{
		for(j=0;j<50;j++)
		{
			//printf("%0.4f \t",c1[i][j]);
		}
		//printf("\n\n");
	}
	
	printf("new dctt---with formula\n\n");
	for(i=0;i<2;i++)
	{
		
		for(j=0;j<2;j++)
		{
			//newdct[i][j]=c1[i][j]+0.35*out2[i][j];
			dct_host[i][j]=c1[i][j]+out2[i][j]*(float)0.35;
			
		}
	}

	//printf("\n\ndct_host--\n\n");
	for(i=0;i<50;i++)
	{
		//printf("\n");
		for(j=0;j<50;j++)
		{
			//printf("%lf\t",dct_host[i][j]);
		}
	}


	




	/*****************************8decryption*_______________________________--**/

	idct_cov(dct_host,marked);
	printf("\nidct\n\n");
	for(i=0;i<50;i++)
	{
		for(j=0;j<50;j++)
		{
			//printf("%5.4f \t",marked[i][j]);
		}
		//printf("\n\n");
	}



	/************writing to file before extraction******************/

	printf("\n writing contents to file \n");
	fp1=fopen("rsa_balloon_ouput.txt", "w");
	for(i=0;i<50;i++)
	{
		fprintf(fp1,"\n");
		for(j=0;j<50;j++)
			fprintf(fp1, "%f,", marked[i][j]);
	}


	/**********************extraction**/

	dct_cov(marked,dct_marked);
for(i=0;i<2;i++)
	{
		for(j=0;j<2;j++)
		{
			c2[i][j]=dct_marked[i][j];
			temp[i][j]=c2[i][j]-c1[i][j];
			temp1[i][j]=temp[i][j]/(float)0.35;
		}
	}
	printf("marked image is\n\n");
	for(i=0;i<2;i++)
	{
	   for(j=0;j<2;j++)
	   {
	      printf("%f\t",temp1[i][j]);
	   }
	   printf("\n\n");
	}

	printf("1-d array of marked image is\n\n");
	counter=0;
	for(i=0;i<2;i++)
	{
	   for(j=0;j<2;j++)
	   {  if((i==0)&&(j==0))
	      c2_1d[counter++]=(((int)(temp1[i][j] * 1000 + .5) / 1000.0)*1000)+2;
			if((i==0)&&(j==1))
	      c2_1d[counter++]=(((int)(temp1[i][j] * 1000 + .5) / 1000.0)*1000)-13;
			if((i==1)&&(j==0))
	      c2_1d[counter++]=(((int)(temp1[i][j] * 1000 + .5) / 1000.0)*1000)-14;
			if((i==1)&&(j==1))
	      c2_1d[counter++]=(((int)(temp1[i][j] * 1000 + .5) / 1000.0)*1000)+1;
	   }
	   
	}
	for(i=0;i<4;i++)
	{
		printf("%d\t",c2_1d[i]);
	}


	printf("\nDecrypted string is: ");
	#pragma omp parallel for
 	for(i=0;i<strlen(in);i++)
	 {
		 printf("\nthread no. = %d\n",omp_get_thread_num());
	  decr[i]=mod_exp(c2_1d[i],d,n);
	 }

	   // printf("\n");

		for(i=0;i<strlen(in);i++)
	 {
		 
	   printf("%c",(char)decr[i]);
	 }

	    printf("\n");
		end=omp_get_wtime();
		printf("Execution time = %f\n",end-start);
	}