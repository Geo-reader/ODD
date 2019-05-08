/* 
This progranm is designed to split the whole velocity and density models into several sub-domains

Coordinate configuration of seismic data:
 is=ix*ny*nz+iy*nz+iz;
 The fastest dim:        *.z
 The second fastest dim: *.y
 The slowest dim:        *.x
 
 Acknowledgement:

   Copyright (C) 2019 China University of Petroleum, Beijing
   Copyright (C) 2019 Ning Wang

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details: http://www.gnu.org/licenses/

*/

#include"stdio.h"
#include"stdlib.h"
#include <string.h>
#define pml 15            //the number of PML layers
#define nz 512            // z- grid number of model(without PML)
#define ny 512            // y- grid number of model(without PML)
#define nx 512            // x- grid number of model(without PML)
#define numpx 16          //the number of sub-domains at x- direction


void num2str(char asc[6],int num);

int main()
{
    int sizex[numpx]; //arrays for storing the number of grid points (with pml) for sub-domains at x- direction
    int xef[numpx];   //arrays for storing the number of grid points (without pml) for sub-domains at x- direction
    int modx;         //the remainder of the number of whole grid points at x- axis(with pml) divided by the number of sub-domains at x- direction
	int avex;

    int nxp=nx+2*pml;// x- grid number of model(with PML)
    modx=nxp%numpx;
    avex=(nxp-modx)/numpx;


    int i,j,k;
	//define the number of grid-points for every sub-domain at x- direction 
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
   
 	 //spliting the whole velocity & density models into subdomains     
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

// ==========================================================
//  This subroutine is used for transfroming number to string
//  =========================================================
void num2str(char asc[6],int num)
{
  char asc1,asc2,asc3,asc4,asc5,asc6;
  asc6='\0';
  if(num<10)
  {
	asc5='0'; 
    asc4='0'; 
    asc3='0'; 
    asc2='0'; 
    asc1=(char)(num+48); 
  }
  if(num>=10&&num<=99)
  {
    asc5='0'; 
    asc4='0'; 
    asc3='0'; 
    asc2=(char)(num/10+48); 
    asc1=(char)(num-num/10*10+48); 
  }
  if(num>99&&num<=999)
  {
	asc5='0';                    
    asc4='0';                   
    asc3=(char)(num/100+48);     
    asc2=(char)(num%100/10+48);  
    asc1=(char)(num%10+48);      
  }
  if(num>=1000&&num<=9999)
  {
	asc5='0'; 
    asc4=(char)(num/1000+48);      
    asc3=(char)((num/100)%10+48);
    asc2=(char)((num/10)%10+48);  
    asc1=(char)(num%10+48);       

  }
  if(num>=10000&&num<=99999)
  {
	  asc5=(char)(num/10000+48); 
      asc4=(char)((num/1000)%10+48); 
      asc3=(char)((num/100)%10+48);
      asc2=(char)((num/10)%10+48); 
      asc1=(char)(num%10+48);        

  }
   asc[0]=asc5;  
   asc[1]=asc4;  
   asc[2]=asc3;  
   asc[3]=asc2;
   asc[4]=asc1;
   asc[5]=asc6;


   return;
}
