#include"stdio.h"
#include"stdlib.h"
#include <string.h>
#define pml 15
#define nz 512
#define ny 512
#define nx 512
#define numpx 4
#define numpy 4

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

    int ix,iy,iz,i,j,is;
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

   float *vp_whole;
   float *vp_each;
   FILE *fp=fopen("vp512.bin","rb");
   FILE *fpp1[numpy][numpx];
   char name[60],ascx[6],ascy[6];
   vp_whole=(float*)malloc(sizeof(float)*nx*ny*nz);

   float *den_whole;
   float *den_each;
   FILE *fp2=fopen("rho512.bin","rb");
   FILE *fpp2[numpy][numpx];
   char name2[60];
   den_whole=(float*)malloc(sizeof(float)*nx*ny*nz);

	for(iz=0;iz<nz;iz++)
	{
		for(iy=0;iy<ny;iy++)
		{
			for(ix=0;ix<nx;ix++)
			{
			is=iz*ny*nx+iy*nx+ix;
			fread(&vp_whole[is],sizeof(float),1,fp);
			fread(&den_whole[is],sizeof(float),1,fp2);
			}
		}
	}	
    
    offset_y=0;
        for(j=0;j<numpy;j++)
        {
         num2str(ascy,j);                  
         
		 offset_x=0;
		for(i=0;i<numpx;i++)
        {
         num2str(ascx,i);          
          sprintf(name,"x%s_vp_y%s.dat",ascx,ascy);
          fpp1[j][i]=fopen(name,"wb"); 
		  vp_each=(float*)malloc(sizeof(float)*xef[i]*yef[j]*nz);

          sprintf(name2,"x%s_rho_y%s.dat",ascx,ascy);
          fpp2[j][i]=fopen(name2,"wb"); 
		  den_each=(float*)malloc(sizeof(float)*xef[i]*yef[j]*nz);
		 
         for(iz=0;iz<nz;iz++)
         
	       {
		    for(iy=0;iy<yef[j];iy++)
		      {
		        for(ix=0;ix<xef[i];ix++)
		        {
			       is=iz*yef[j]*xef[i]+iy*xef[i]+ix;
				   vp_each[is]=vp_whole[iz*ny*nx+(iy+offset_y)*nx+(ix+offset_x)];
			       fwrite(&vp_each[is],sizeof(float),1,fpp1[j][i]);

				   den_each[is]=den_whole[iz*ny*nx+(iy+offset_y)*nx+(ix+offset_x)];
			       fwrite(&den_each[is],sizeof(float),1,fpp2[j][i]);
			    }
		     }  
		   }

			offset_x=offset_x+xef[i];
		 } 
		 offset_y=offset_y+yef[j];  
		}

        for(j=0;j<numpy;j++)
        {

		    for(i=0;i<numpx;i++)
              {
			   fclose(fpp1[j][i]);
			   fclose(fpp2[j][i]);
		      }
		}
  
    fclose(fp);fclose(fp2);
    free(vp_whole);free(vp_each);free(den_whole);free(den_each);    

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
