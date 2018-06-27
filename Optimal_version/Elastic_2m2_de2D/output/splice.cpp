#include"stdio.h"
#include"stdlib.h"
#include <string.h>
#define pml 15
#define nt 1001
#define ny 101
#define nx 101
#define numpx 2
#define numpy 2

void num2str(char asc[6],int num);

int main()
{
    int sizex[numpx];
    int xef[numpx];
	int sizey[numpy];
    int yef[numpy];
    int modx,avex;
    int mody,avey;

    int nxp=nx+2*pml;
    modx=nxp%numpx;
    avex=(nxp-modx)/numpx;

    int nyp=ny+2*pml;
    mody=nyp%numpy;
    avey=(nyp-mody)/numpy;

    int ix,iy,it,i,j,is;
    int offset_y,offset_x;
    for(i=0;i<numpx;i++)
    {
        sizex[i]=avex;
        if(i<modx){sizex[i]=sizex[i]+1;}

        xef[i]=sizex[i];
        if(i==0){xef[i]=xef[i]-pml;}
        if(i==numpx-1){xef[i]=xef[i]-pml;}

        printf("i=%d size=%d\n",i,xef[i]);

    }


    for(j=0;j<numpy;j++)
    {
        sizey[j]=avey;
        if(j<mody){sizey[j]=sizey[j]+1;}

        yef[j]=sizey[j];
        if(j==0){yef[j]=yef[j]-pml;}
        if(j==numpy-1){yef[j]=yef[j]-pml;}

      printf("j=%d size=%d\n",j,yef[j]);

    }


     FILE *fp=fopen("rec_whole.bin","wb");

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

   FILE *fpp[numpy][numpx];
   char name[60],ascx[6],ascy[6];
   

        offset_y=0;
        for(j=0;j<numpy;j++)
        {
			num2str(ascy,j);  
			offset_x=0;
			for(i=0;i<numpx;i++)	 
			{ 
				num2str(ascx,i);  				
				rec_each=(float*)malloc(sizeof(float)*xef[i]*yef[j]*nt);
				sprintf(name,"00000-record_elas2m2-x%s-y%s.dat",ascx,ascy);///record_elas2m2
                fpp[j][i]=fopen(name,"rb"); 
				

				for(it=0;it<nt;it++)
                 {
					for(iy=0;iy<yef[j];iy++)
					 {
						 for(ix=0;ix<xef[i];ix++)
							{
								is=it*yef[j]*xef[i]+iy*xef[i]+ix;
                                fread(&rec_each[is],sizeof(float),1,fpp[j][i]);
								rec_whole[it*ny*nx+(iy+offset_y)*nx+(ix+offset_x)]=rec_each[is];
						 }
					}
				}
				offset_x=offset_x+xef[i];
			}
			offset_y=offset_y+yef[j];
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
