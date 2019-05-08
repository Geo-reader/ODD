
/* 
This progranm is designed to splice the sub-records into a whole seismic record

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

#define pml 15             //the number of PML layers
#define nt 501            // the number of time steps
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

    int ix,iy,it,i,j,is;
    int offset_x;

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

	 //open a file to store the whole seismic record
     FILE *fp=fopen("rec_final.bin","wb");

   float *rec_whole;
   rec_whole=(float*)malloc(sizeof(float)*nx*ny*nt);

   
	for(it=0;it<nt;it++)
	{
		for(iy=0;iy<ny;iy++)
		{
			for(ix=0;ix<nx;ix++)
			{
			is=it*ny*nx+iy*nx+ix;
			rec_whole[is]=0.0;
			}
		}
	}	

 // splicing all sub-records into the whole seismic record
 float *rec_each;

   FILE *fpp[numpx];
   char name[60],ascx[6];
   
			offset_x=0;
			for(i=0;i<numpx;i++)	 
			{ 
				num2str(ascx,i);  				
				rec_each=(float*)malloc(sizeof(float)*xef[i]*ny*nt);
				sprintf(name,"00000-record2m4-x%s.dat",ascx);
                fpp[i]=fopen(name,"rb"); 
				

				for(it=0;it<nt;it++)
                 {
					for(ix=0;ix<xef[i];ix++)
					 {
						 for(iy=0;iy<ny;iy++)
							{
								is=it*xef[i]*ny+ix*ny+iy;
                                fread(&rec_each[is],sizeof(float),1,fpp[i]);
								rec_whole[it*nx*ny+(ix+offset_x)*ny+iy]=rec_each[is];
						 }
					}
				}
				offset_x=offset_x+xef[i];
			}


 fwrite(rec_whole,sizeof(float),nx*ny*nt,fp);

 
    fclose(fp);
    free(rec_whole);           
    
    
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
