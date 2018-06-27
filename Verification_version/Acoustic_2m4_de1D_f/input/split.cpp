#include"stdio.h"
#include"stdlib.h"
#include <string.h>
#define pml 15
#define nz 512
#define ny 512
#define nx 512
#define numpx 16


void num2str(char asc[6],int num);

int main()
{
    int sizex[numpx];
    int xef[numpx];
    int modx,avex;


    int nxp=nx+2*pml;
    modx=nxp%numpx;
    avex=(nxp-modx)/numpx;


    int i,j,k;
    for(i=0;i<numpx;i++)
    {
        sizex[i]=avex;
        if(i<modx){sizex[i]=sizex[i]+1;}

        xef[i]=sizex[i];
        if(i==0){xef[i]=xef[i]-pml;}
        if(i==numpx-1){xef[i]=xef[i]-pml;}

        printf("i=%d size=%d\n",i,xef[i]);

    }

    FILE *fp=fopen("rho512.bin","rb");
    FILE *fp1;
    char name[40],ascx[6],ascy[6];
	float *yoz;

    FILE *fpv=fopen("vp512.bin","rb");
    FILE *fpv1;
    char namev[40];
    float *yozv;
   
    
    for(i=0;i<numpx;i++)
    {
	    yoz=(float*)malloc(sizeof(float)*ny*nz);
		yozv=(float*)malloc(sizeof(float)*ny*nz);

        num2str(ascx,i);
        strcpy(name,"");
        sprintf(name,"x%s-rho.dat",ascx);
        fp1=fopen(name,"wb");

        num2str(ascx,i);
        strcpy(namev,"");
        sprintf(namev,"x%s-vp.dat",ascx);
        fpv1=fopen(namev,"wb");

        k=0;
        do{
           fread(yoz,sizeof(float),ny*nz,fp);
           fwrite(yoz,sizeof(float),ny*nz,fp1);
           fread(yozv,sizeof(float),ny*nz,fpv);
           fwrite(yozv,sizeof(float),ny*nz,fpv1);

           k=k+1;
          }while(k<xef[i]);
         fclose(fp1);
		 fclose(fpv1);
		 free(yoz);
		 free(yozv);
	}
    fclose(fp);fclose(fpv);
    
    
    return 0;
}

void num2str(char asc[6],int num)
{
  char asc1,asc2,asc3,asc4,asc5,asc6;
  asc6='\0';
  if(num<10)
  {
	asc5='0'; // wan
    asc4='0'; // qian
    asc3='0'; // bai
    asc2='0'; // shi
    asc1=(char)(num+48); // ge
  }
  if(num>=10&&num<=99)
  {
    asc5='0'; // wan
    asc4='0'; // qian
    asc3='0'; // bai
    asc2=(char)(num/10+48); // shi
    asc1=(char)(num-num/10*10+48); // ge
  }
  if(num>99&&num<=999)
  {
	asc5='0';                    // wan
    asc4='0';                    // qian
    asc3=(char)(num/100+48);     //bai
    asc2=(char)(num%100/10+48);  //shi
    asc1=(char)(num%10+48);      //ge
  }
  if(num>=1000&&num<=9999)
  {
	asc5='0'; // wan
    asc4=(char)(num/1000+48);      //qian
    asc3=(char)((num/100)%10+48); //bai
    asc2=(char)((num/10)%10+48);  //shi
    asc1=(char)(num%10+48);       //ge

  }
  if(num>=10000&&num<=99999)
  {
	  asc5=(char)(num/10000+48); // wan
      asc4=(char)((num/1000)%10+48); //qian
      asc3=(char)((num/100)%10+48);//bai
      asc2=(char)((num/10)%10+48); //shi
      asc1=(char)(num%10+48);        //ge

  }
   asc[0]=asc5;  // bai
   asc[1]=asc4;  // shi
   asc[2]=asc3;  // ge
   asc[3]=asc2;
   asc[4]=asc1;
   asc[5]=asc6;


   return;
}
