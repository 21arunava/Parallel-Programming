#include<stdio.h>
#include<math.h>
#include<omp.h>

/*int wm[10][10]=	{1,0,1,1,0,0,0,0,0,0,
		1,0,1,1,0,0,0,0,0,0,
		1,0,1,1,0,0,0,0,0,0,
		1,0,1,1,0,0,0,0,0,0,
		1,0,1,1,0,0,0,0,0,0,
		1,0,1,1,0,0,0,0,0,0,
		1,0,1,1,0,0,0,0,0,0,
		1,0,1,1,0,0,0,0,0,0,
		1,0,1,1,0,0,0,0,0,0,
		1,0,1,1,0,0,0,0,0,0};*/

/*float host[6][6]={0.34,0.12,0.52,0.19,0.74,0.67,
		  0.87,0.14,0.27,0.19,0.61,0.39,
		  0.23,0.89,0.09,0.45,0.18,0.71,
		  0.12,0.67,0.12,0.200,0.71,0.29,
		  0.87,0.14,0.27,0.19,0.61,0.39,
		  0.23,0.89,0.09,0.45,0.18,0.71,};*/
float host[200][200];
int wm[50][50];

void main()
{
	float dct_wm(int z[50][50],float u[50][50]);
	float dct_cov(float z[200][200],float u[200][200]);
	float idct_wm(float z[50][50],float u[50][50]);
	float idct_cov(float z[200][200],float u[200][200]);
	float average(float a[2500],int n);
	float corr(int m1[50][50], int m2[50][50]);
	float mean;
	int p[50][50];
	int mark[50][50];
	float dct_mark[50][50];
	float a=2.75,b=0.15;
	float dct_host[200][200];
	float c1[50][50];
	float newdct[50][50];
	int seq[2500];
	float marked[200][200];
	float dct_marked[200][200];
	float c2[50][50];
	float temp[50][50];
	float temp1[50][50];
	float wm_encr[50][50];
	int wm_extr[50][50];
	float q[50][50];
	int v[50][50],w[50][50];
	int mark_encr[50][50];

	float x0=1.4,y0=0.7;
	int i,j,k;
	float x[2500],y[2500];
	double start,end;


	FILE *fp,*fp1,*fp3,*fp4,*fp5,*fp6;
float host_1d[41000];
  float value;
  int ct = -1; /* EDIT have i start at -1 :) */
  int ar_ct=0;

  int wm_1d[2510];
  int value1;
  int ct1=-1;
  int ar_ct1=0;

  //  reading cover image
  if ((fp = fopen ("200_balloon.txt", "r")) == NULL)
    printf("Error");

  while (!feof (fp) && fscanf (fp, "%f,", &value) && ct++ < 41000 )
    host_1d[ct] = value;

  fclose (fp);
	  //printf("host image 1d array is \n\n");

	  /*for(i=0;i<350;i++)
		  printf("%f\t",host_1d[i]);
	  printf("\n");*/

	  //printf("\nhost image 2d array is :");
	for(i=0;i<200;i++)
	{
		//printf("\n:");
	 for(j=0;j<200;j++)
	 {
		 host[i][j]=host_1d[ar_ct++];
		 //printf("%f\t",host[i][j]);
    }
	}

	//reading watermark image

	if ((fp3 = fopen ("50_tree.txt", "r")) == NULL)
    printf("Error");

  while (!feof (fp3) && fscanf (fp3, "%d,", &value1) && ct1++ < 2510 )
    wm_1d[ct1] = value1;

  fclose (fp3);
	  //printf("watermark image 1d array is \n\n");

	  /*for(i=0;i<350;i++)
		  printf("%f\t",host_1d[i]);
	  printf("\n");*/

	 // printf("\nwatermark image 2d array is :");
	for(i=0;i<50;i++)
	{
		//printf("\n:");
	 for(j=0;j<50;j++)
	 {
		 wm[i][j]=wm_1d[ar_ct1++];
		 //printf("%f\t",host[i][j]);
    }
	}


	start=omp_get_wtime();
  /*----------------------------Chaotic sequence generator---------------------------------*/
#pragma omp parallel for
	for(i=0;i<2500;i++)
	{
		x[i]=y0;
		y[i]=(-b*x0)+(a*y0)-(y0*y0*y0);
		y0=y[i];
		x0=x[i];
		printf("Thread id in chaos parallel is %d\n",omp_get_thread_num());
	}

	/*printf("chaotic sequence=\n");
	for(j=0;j<2500;j++)
		printf("%0.4f \t",x[j]);
		printf("\n\n");*/
/*--------------------------Convert to binary--------------------------------------*/

	mean=average(x,2500);

	for(i=0;i<2500;i++)
	{
		if(x[i]>mean)
		      seq[i]=1;
		else
			seq[i]=0;
	}

	for(i=0;i<50;i++)
	 for(j=0;j<50;j++)
	  p[i][j]=seq[j+(50*i)];

	printf("average=\n");
	printf("%0.4f\n\n",mean);

	printf("binary chaotic sequence=\n");

/*------------------------------Encryption of watermark----------------------------*/

	for(i=0;i<50;i++)

		for(j=0;j<50;j++)
			mark[i][j]=wm[i][j]^p[i][j];


	printf("encrypted watermark=\n\n");
	for(i=0;i<50;i++)
	{
		for(j=0;j<50;j++)
		{
			printf("%d \t",mark[i][j]);
		}
		printf("\n\n");
	}

	printf("\n writing contents to file3 \n");
	fp5=fopen("output_file_encrypted_balloon.txt", "w");
	for(i=0;i<50;i++)
	{
		fprintf(fp5,"\n");
		for(j=0;j<50;j++)
			fprintf(fp5, "%d,", mark[i][j]);
	}

	dct_wm(mark,dct_mark);

	printf("dct coeft sequence=\n");
	for(i=0;i<50;i++)
	{
		for(j=0;j<50;j++)
		{
			printf("%0.4f \t",dct_mark[i][j]);
		}
		printf("\n\n");
	}


/*----------------------------Computing 2-D DCT of host-----------------------------*/

	dct_cov(host,dct_host);
	for(i=0;i<50;i++)
		for(j=0;j<50;j++)
			c1[i][j]=dct_host[i][j];

	printf("dct coefts of host=\n");
	for(i=0;i<50;i++)
	{
		for(j=0;j<50;j++)
		{
			printf("%0.4f \t",c1[i][j]);
		}
		printf("\n\n");
	}

	for(i=0;i<50;i++)
	{
		for(j=0;j<50;j++)
		{
			newdct[i][j]=c1[i][j]+0.15*dct_mark[i][j];
			dct_host[i][j]=newdct[i][j];
		}
	}

	idct_cov(dct_host,marked);


/*-------------------------END OF EMBEDDING PROCESS---------------*/


	/************writing to file before extraction******************/

	printf("\n writing contents to file1 \n");
	fp1=fopen("output_file_chaos_balloon.txt", "w");
	for(i=0;i<200;i++)
	{
		fprintf(fp1,"\n");
		for(j=0;j<200;j++)
			fprintf(fp1, "%f,", marked[i][j]);
	}
	






/*-------------------------EXTRACTION PROCESS---------------------*/

	dct_cov(marked,dct_marked);

	/*printf("\n writing contents to file \n");
	fp1=fopen("output_file_chaos.txt", "w");
	for(i=0;i<200;i++)
	{
		fprintf(fp1,"\n");
		for(j=0;j<200;j++)
			fprintf(fp1, "%f,", dct_marked[i][j]);
	}
	*/
	for(i=0;i<50;i++)
	{
		for(j=0;j<50;j++)
		{
			c2[i][j]=dct_marked[i][j];
			temp[i][j]=c2[i][j]-c1[i][j];
			temp1[i][j]=temp[i][j]/0.15;
		}
	}
	printf("marked image is\n\n");
	for(i=0;i<50;i++)
	{
	   for(j=0;j<50;j++)
	   {
	      printf("%0.4f\t",c2[i][j]);
	   }
	   printf("\n\n");
	}

	idct_wm(temp1,wm_encr);

/*----------------------------------------------------------------------------------------*/

/*---------------------------Decryption of watermark------------------------------*/

	for(i=0;i<50;i++)
	{
		for(j=0;j<50;j++)
		{
		q[i][j]=wm_encr[i][j]+0.5;
		v[i][j]=q[i][j];
		w[i][j]=wm_encr[i][j];

		if(v[i][j]==w[i][j])
			mark_encr[i][j]=(floor(wm_encr[i][j]));
		else
			mark_encr[i][j]=(ceil(wm_encr[i][j]));

		if(wm_encr[i][j]<0)
			mark_encr[i][j]=0;
		}
	}

	printf("\n writing contents to file4 \n");
	fp6=fopen("output_file_encrypted1_balloon.txt", "w");
	for(i=0;i<50;i++)
	{
		fprintf(fp6,"\n");
		for(j=0;j<50;j++)
			fprintf(fp6, "%d,", mark_encr[i][j]);
	}


	for(i=0;i<50;i++)

		for(j=0;j<50;j++)
			wm_extr[i][j]=mark_encr[i][j]^p[i][j];

	/*for(i=0;i<10;i++)
	{
		for(j=0;j<10;j++)
		{
			printf("%d\t",wm_extr[i][j]);
		}
		printf("\n\n");
	}*/

	/************writing to file extracted watermark******************/

	printf("\n writing contents to file2 \n");
	fp4=fopen("output_file_watermark_balloon.txt", "w");
	for(i=0;i<50;i++)
	{
		fprintf(fp4,"\n");
		for(j=0;j<50;j++)
			fprintf(fp4, "%d,", wm_extr[i][j]);
	}

	end=omp_get_wtime();
	printf("Execution time = %f",end-start);
  
}
/*--------------------------------Function for 2-D DCT-----------------------------*/

float average(float a[2500],int n)
{
	float s=0,m;
	int i;
	for(i=0;i<2500;i++)
		s=s+a[i];
	m=s/n;
	return(m);
}


float dct_wm(int z[50][50],float u[50][50])
{

	int p,q,i,j,k,l;
	float s,t;
	float d[50][50]={0};
	#pragma omp parallel shared(z,d) private(p,q,s,t,i,j)
	{
		printf("\nthread id=%d\n",omp_get_thread_num()); 
#pragma omp for schedule(static)
	for(p=0;p<50;p++)

		for(q=0;q<50;q++)
		{
			if(p==0)

				s=sqrt(1.0/50.0);

			else

				s=sqrt(2.0/50.0);
			if(q==0)

			       t=sqrt(1.0/50.0);
			else
				t=sqrt(2.0/50.0);


		      for(i=0;i<50;i++)

				for(j=0;j<50;j++)
				{
d[p][q]+=z[i][j]*cos((3.1412*(2*i+1)*p)/(2*50))*cos((3.1412*(2*j+1)*q)/(2*50));
				}
			d[p][q]*=s*t;

		}
	}
		for(i=0;i<50;i++)
		  for(j=0;j<50;j++)
			u[i][j]=d[i][j];
		return(u[50][50]);

}


float dct_cov(float z[200][200],float u[200][200])
{

	int p,q,i,j,k,l;
	float s,t;
	float d[200][200]={0};
	#pragma omp parallel shared(z,d) private(p,q,s,t,i,j)
	{
		printf("\nthread id=%d\n",omp_get_thread_num()); 
#pragma omp for schedule(static)
	for(p=0;p<200;p++)

		for(q=0;q<200;q++)
		{
			if(p==0)

				s=sqrt(1.0/200.0);

			else

				s=sqrt(2.0/200.0);
			if(q==0)

			       t=sqrt(1.0/200.0);
			else
				t=sqrt(2.0/200.0);


		      for(i=0;i<200;i++)

				for(j=0;j<200;j++)
				{
d[p][q]+=z[i][j]*cos((3.1412*(2*i+1)*p)/(2*200))*cos((3.1412*(2*j+1)*q)/(2*200));
				}
			d[p][q]*=s*t;

		}
	}
		for(i=0;i<200;i++)
		  for(j=0;j<200;j++)
		    u[i][j]=d[i][j];
		return(u[200][200]);
}

/*----------------------IDCT FUNCTIONS--------------------------------*/
float idct_cov(float z[200][200],float u[200][200])
{
	int p,q,i,j,k,l;
	float s,t;
	float f[200][200]={0};
	#pragma omp parallel shared(z,f) private(p,q,s,t,i,j)
	{
		printf("\nthread id=%d\n",omp_get_thread_num()); 
#pragma omp for schedule(static)
	for(i=0;i<200;i++)

		for(j=0;j<200;j++)
		{
			for(p=0;p<200;p++)
				for(q=0;q<200;q++)
				{
				if(p==0)

					s=sqrt(1.0/200.0);

				else

					s=sqrt(2.0/200.0);
				if(q==0)

					t=sqrt(1.0/200.0);
				else
					t=sqrt(2.0/200.0);

f[i][j]+=(s*t*z[p][q]*cos((3.1412*(2*i+1)*p)/(2*200))*cos((3.1412*(2*j+1)*q)/(2*200)));

				}

	       }
	}
	       for(i=0;i<200;i++)
		 for(j=0;j<200;j++)
		   u[i][j]=f[i][j];
  return(u[200][200]);
}


float idct_wm(float z[50][50],float u[50][50])
{
	int p,q,i,j,k,l;
	float s,t;
	float f[50][50]={0};
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

					s=sqrt(1.0/50.0);

				else

					s=sqrt(2.0/50.0);
				if(q==0)

					t=sqrt(1.0/50.0);
				else
					t=sqrt(2.0/50.0);

f[i][j]+=(s*t*z[p][q]*cos((3.1412*(2*i+1)*p)/(2*50))*cos((3.1412*(2*j+1)*q)/(2*50)));
				}

	       }
	}

	       for(i=0;i<50;i++)
		  for(j=0;j<50;j++)
		     u[i][j]=f[i][j];


return(u[50][50]);
}

