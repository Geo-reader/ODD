#include"stdio.h"
#include"stdlib.h"
#include <string.h>
#define pml 15
#define nt 501
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



    int ix,iy,it,i,j,is;
    int offset_x;
    for(i=0;i<numpx;i++)
    {
        sizex[i]=avex;
        if(i<modx){sizex[i]=sizex[i]+1;}

        xef[i]=sizex[i];
        if(i==0){xef[i]=xef[i]-pml;}
        if(i==numpx-1){xef[i]=xef[i]-pml;}

        printf("i=%d size=%d\n",i,xef[i]);

    }



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
