/* 
3D acoustic modeling using MPI spliting x and y directions(optimized domain decomposition method)
(2M,4) SGFD scheme

Coordinate configuration of seismic data:
 is=iz*ny*nx+iy*nx+ix;
 The fast dim:        *.x
 The second fast dim: *.y
 The slowest dim:     *.z
     

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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include "omp.h"
#include "Myfunctions.h"

#define pi 3.141592653
#define pml 15



int main(int argc,char *argv[])
{

	//=========================================================
	//  MPI Indix
	//  =======================================================
	MPI_Comm comm=MPI_COMM_WORLD;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
    int myid,numprocs,namelen;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(comm,&myid);
    MPI_Comm_size(comm,&numprocs);
    MPI_Get_processor_name(processor_name,&namelen);

	//=========================================================
	//  Parameters of the time of the system...
	//  =======================================================
	double starttime,endtime;

	//=========================================================
	//  input modeling parameters & Bcast
	//  =======================================================

    float rectime;			    // recording time
	float dt;				    // time sample  
	float fm;                   // Riker wavelet domiant frequency
    float dx,dy,dz;             // grid interval of x, y and z
    int shotnum;                // the number of shot 
    int nxo,nyo,nzo;            // x-,y- and z-grid number of model(without PML)
    int nx,ny,nz;               // x-,y- and z-grid number of model(including PML)
    int Lx,Ly,Lz;               // FD order for x, y and z derivatives
	int numx,numy;              // sub-domains in x- and y- direction
    char filepathv[40];         // path of record to be saved 
    char filepathr[40];         // path of velocity & density
	int tmax;                   // the number of time-steps for simulation
	int i,j,k,n,t;
       
	if(myid==0)
	{

	   input_parameters(&rectime,&dt,&fm,&shotnum,&dx,&dy,&dz,
	   &nxo,&nyo,&nzo,&Lx,&Ly,&Lz,&numx,&numy,filepathv,filepathr);

	   tmax=(int)(rectime/dt*1000)+1;
           nx=nxo+2*pml;
           ny=nyo+2*pml;
           nz=nzo+2*(pml+Lz);
           dt=dt/1000.0;

	}
        
	    MPI_Bcast(&tmax,1,MPI_INT,0,comm);
        MPI_Bcast(&nx,1,MPI_INT,0,comm);
        MPI_Bcast(&ny,1,MPI_INT,0,comm);
        MPI_Bcast(&nz,1,MPI_INT,0,comm);
        MPI_Bcast(&nzo,1,MPI_INT,0,comm);
        MPI_Bcast(&nyo,1,MPI_INT,0,comm);
        MPI_Bcast(&nxo,1,MPI_INT,0,comm);   
        MPI_Bcast(&Lx,1,MPI_INT,0,comm);
        MPI_Bcast(&Ly,1,MPI_INT,0,comm);
        MPI_Bcast(&Lz,1,MPI_INT,0,comm);
       
        MPI_Bcast(&shotnum,1,MPI_INT,0,comm);
        MPI_Bcast(&dt,1,MPI_FLOAT,0,comm);
        MPI_Bcast(&dx,1,MPI_FLOAT,0,comm);
        MPI_Bcast(&dy,1,MPI_FLOAT,0,comm);
        MPI_Bcast(&dz,1,MPI_FLOAT,0,comm);
        MPI_Bcast(&fm,1,MPI_FLOAT,0,comm);
	    MPI_Bcast(&numx,1,MPI_INT,0,comm);
	    MPI_Bcast(&numy,1,MPI_INT,0,comm);

        MPI_Bcast(filepathv,40,MPI_CHAR,0,comm);
        MPI_Bcast(filepathr,40,MPI_CHAR,0,comm);


	//=========================================================
	//  input geometry parameters & Bcast
	//  =======================================================	
	int *shotxg,*shotyg,*shotzg;   // location of the sources
	int recg;                      // depth of recievers
	shotxg=(int*)malloc(sizeof(int)*shotnum);
    shotyg=(int*)malloc(sizeof(int)*shotnum);
	shotzg=(int*)malloc(sizeof(int)*shotnum);
	
	if(myid==0)
	{
	   input_shotsposition(shotnum,dx,dy,dz,shotxg,shotyg,shotzg,&recg);
	}
	MPI_Bcast(shotxg,shotnum,MPI_INT,0,comm);
    MPI_Bcast(shotyg,shotnum,MPI_INT,0,comm);
	MPI_Bcast(shotzg,shotnum,MPI_INT,0,comm);
    MPI_Bcast(&recg,1,MPI_INT,0,comm);
 

	//=========================================================
	//  parameters of splitting domain
	//  =======================================================	  
	int avex,avey,insizex,insizey,modx,mody;
	int sizex,sizey;          // samples of x- and y- directions for sub-domains (including ghost areas)
	int xef,yef,zef;          // samples of x- and y- directions for sub-domains (without ghost areas)
	int idx,idy;              // ID number of process for x- and y- directions
	int offsetrowx,offsetrowy,offsetrowxn,offsetrowyn;
	int xl,xr,yl,yr,zu,zd;
	zef=nzo;        
        

	   split_domains(nx,ny,nxo,nyo,Lx,Ly,nz,Lz,myid,&avex,&avey,&insizex,&insizey,
		         &modx,&mody,&sizex,&sizey,&offsetrowx,&offsetrowy,
		         &offsetrowxn,&offsetrowyn,numx,numy,&idx,&idy,
		         &xl,&xr,&yl,&yr,&zu,&zd,&xef,&yef);
        


	//=========================================================
	//  locate the source's rank and grid
	//  =======================================================	
	   int *shotxlc,*shotylc,*shotzlc;   // location of sources
	   int *flagshots;    

	   	flagshots=(int*)malloc(sizeof(int)*shotnum);
	    shotxlc=(int*)malloc(sizeof(int)*shotnum);
	    shotylc=(int*)malloc(sizeof(int)*shotnum);
	    shotzlc=(int*)malloc(sizeof(int)*shotnum);

        for(i=0;i<shotnum;i++)
	    {
	        flagshots[i]=0;
	        shotxlc[i]=0;
	        shotylc[i]=0;
	        shotzlc[i]=0;

         if(shotxg[i]>=offsetrowxn&&shotxg[i]<offsetrowxn+xef&&
	      shotyg[i]>=offsetrowyn&&shotyg[i]<offsetrowyn+yef)
	       {
	    	flagshots[i]=1;
		    shotxlc[i]=(shotxg[i]-offsetrowxn)+xl;
		    shotylc[i]=(shotyg[i]-offsetrowyn)+yl;
		    shotzlc[i]=zu+shotzg[i];
	       }
	    }

	//=========================================================
	//  input velocity and density
	//  =======================================================	
        float *modelv,*modeld;     //velocity & density models for sub-domains (without ghost areas)
		float *vp,*dens;           //velocity & density models for sub-domains (including ghost areas)

        modelv=(float*)malloc(sizeof(float)*xef*yef*nzo);
        modeld=(float*)malloc(sizeof(float)*xef*yef*nzo);
        vp=(float*)malloc(sizeof(float)*sizex*sizey*nz);
        dens=(float*)malloc(sizeof(float)*sizex*sizey*nz);
             
        input_v_dens(filepathv,xef,yef,zef,idx,idy,modelv,modeld);

        input_model_parameters(vp,dens,modelv,modeld,xef,yef,nzo,xl,xr,yl,yr,zu,zd,nz,sizex,sizey);



	//=========================================================
	//  get max and min velocity
	//  =======================================================	
     int Ns=xef*yef*xef;
     float vpmin; float vpmax;               
     get_max_min(modelv,Ns,&vpmax,&vpmin);
     float dv=5.0;
     int Nf=(int)((vpmax-vpmin)/dv+1.5); // number of FD coefficients
     int Lpx=Lx+2,Lpy=Ly+2,Lpz=Lz+2;

     if(dy==dx){Lpy=Lpx;}
     if(dz==dx){Lpz=Lpx;}

     

	//=========================================================
	//  FD coefficients 
	//  =======================================================	
     float *Gp,*Gpy,*Gpz;   // (2m,4) FD coefficients for x- ,y- and z directions

     Gp=(float*)malloc(sizeof(float)*Nf*Lpx);
     Gpy=(float*)malloc(sizeof(float)*Nf*Lpy);
     Gpz=(float*)malloc(sizeof(float)*Nf*Lpz);

     TE_2M4_3dxyz(vpmin,dv,Nf,Lpx,dt,dx,dy,dz,Gp); // for x direction


     free(modelv);free(modeld);

	//=========================================================
	//  PML parameters
	//  =======================================================	
    float *i_a_x,*i_b_x,*i_a_z,*i_b_z,*i_a_y,*i_b_y;
    float *h_a_x,*h_b_x,*h_a_z,*h_b_z,*h_a_y,*h_b_y;
    float *i_ax,*i_bx,*h_ax,*h_bx;
    float *i_ay,*i_by,*h_ay,*h_by;

            i_a_x=(float*)malloc(sizeof(float)*nx);
            i_b_x=(float*)malloc(sizeof(float)*nx);
            h_a_x=(float*)malloc(sizeof(float)*nx);
            h_b_x=(float*)malloc(sizeof(float)*nx);

            i_a_y=(float*)malloc(sizeof(float)*ny);
            i_b_y=(float*)malloc(sizeof(float)*ny);
            h_a_y=(float*)malloc(sizeof(float)*ny);
            h_b_y=(float*)malloc(sizeof(float)*ny);
    
            i_a_z=(float*)malloc(sizeof(float)*nz);
            i_b_z=(float*)malloc(sizeof(float)*nz);
            h_a_z=(float*)malloc(sizeof(float)*nz);
            h_b_z=(float*)malloc(sizeof(float)*nz);

            cal_pml_parameters(i_a_x,i_b_x,h_a_x,h_b_x,i_a_z,i_b_z,h_a_z,h_b_z,
                               i_a_y,i_b_y,h_a_y,h_b_y,nx,ny,nz,Lz,dx,dy,dz,dt,fm);
            
            i_ax=(float*)malloc(sizeof(float)*insizex);
            i_bx=(float*)malloc(sizeof(float)*insizex);
            h_ax=(float*)malloc(sizeof(float)*insizex);
            h_bx=(float*)malloc(sizeof(float)*insizex);

            i_ay=(float*)malloc(sizeof(float)*insizey);
            i_by=(float*)malloc(sizeof(float)*insizey);
            h_ay=(float*)malloc(sizeof(float)*insizey);
            h_by=(float*)malloc(sizeof(float)*insizey);


            for(i=0;i<insizex;i++)
            {
                i_ax[i]=i_a_x[offsetrowx+i];
                i_bx[i]=i_b_x[offsetrowx+i];
                h_ax[i]=h_a_x[offsetrowx+i];
                h_bx[i]=h_b_x[offsetrowx+i];
            }
            for(i=0;i<insizey;i++)
            {
                i_ay[i]=i_a_y[offsetrowy+i];
                i_by[i]=i_b_y[offsetrowy+i];
                h_ay[i]=h_a_y[offsetrowy+i];
                h_by[i]=h_b_y[offsetrowy+i];
            }
            free(i_a_x);free(i_b_x);free(h_a_x);free(h_b_x);
            free(i_a_y);free(i_b_y);free(h_a_y);free(h_b_y);

	//=========================================================
	//  wavefield variables
	//  =======================================================	
        float *p,*vx,*vy,*vz;
        float *fi_px,*fi_py,*fi_pz,*fi_vx,*fi_vy,*fi_vz;

        p=(float*)malloc(sizeof(float)*sizex*sizey*nz);
        vx=(float*)malloc(sizeof(float)*sizex*sizey*nz);
        vy=(float*)malloc(sizeof(float)*sizex*sizey*nz);
        vz=(float*)malloc(sizeof(float)*sizex*sizey*nz);
        fi_px=(float*)malloc(sizeof(float)*sizex*sizey*nz);
        fi_py=(float*)malloc(sizeof(float)*sizex*sizey*nz);
        fi_pz=(float*)malloc(sizeof(float)*sizex*sizey*nz);
        fi_vx=(float*)malloc(sizeof(float)*sizex*sizey*nz);
        fi_vy=(float*)malloc(sizeof(float)*sizex*sizey*nz);
        fi_vz=(float*)malloc(sizeof(float)*sizex*sizey*nz);
            

        
        MPI_Barrier(comm);
        starttime=MPI_Wtime();

        
	//=========================================================
	//  shot number iteration
	//  =======================================================	
        for(n=0;n<shotnum;n++)
        {

               forward3d_c(sizex,sizey,nz,tmax,dx,dy,dz,dt,vpmin,dv,Gp,vp,dens,xl,xr,yl,yr,zu,zd,
	                   insizex,insizey,nxo,nyo,nzo,shotxlc,shotylc,shotzlc,recg,flagshots,
                           i_ax,i_bx,h_ax,h_bx,i_ay,i_by,h_ay,h_by,i_a_z,i_b_z,h_a_z,h_b_z,
                           fm,n,xef,yef,zef,Lpx,Lpy,Lpz,Nf,Lx,Ly,Lz,idx,idy,myid,numx,numy,p,vx,vy,vz,
                           fi_px,fi_py,fi_pz,fi_vx,fi_vy,fi_vz,numprocs,offsetrowxn,offsetrowyn,filepathr);
         
        }
         
        endtime=MPI_Wtime();

        if(myid==0){
        printf("It tooks %f secodes\n", endtime-starttime);}


           free(shotxg);free(shotyg);free(shotzg);
           free(p);free(vx);free(vy);free(vz);
           free(fi_px);free(fi_py);free(fi_pz);
           free(fi_vx);free(fi_vy);free(fi_vz);
           free(i_ax);free(i_bx);free(h_ax);free(h_bx);
           free(i_ay);free(i_by);free(h_ay);free(h_by);
           free(i_a_z);free(i_b_z);free(h_a_z);free(h_b_z);
           free(Gp);free(Gpy);free(Gpz);free(vp);free(dens);
           free(flagshots);free(shotxlc);free(shotylc);free(shotzlc);

        MPI_Finalize();
        return 0;
}

// ==========================================================
//  This subroutine is used for forward modelling ...
//  =========================================================
void forward3d_c(int sizex,int sizey,int nz,int tmax,float dx,float dy,float dz,float dt,
	       float vpmin,float dv,float *Gp,float *vp,float *dens,
	       int xl,int xr,int yl,int yr,int zu,int zd,int insizex,int insizey,
	       int nxo,int nyo,int nzo, int *shotxlc,int *shotylc,int *shotzlc,
           int recg,int *flagshots,float *i_ax,float *i_bx,float *h_ax,float *h_bx,
	       float *i_ay,float *i_by,float *h_ay,float *h_by,
	       float *i_a_z,float *i_b_z,float *h_a_z,float *h_b_z,
	       float fm,int n,int xef,int yef,int zef,int Lpx,int Lpy,int Lpz,int Nf,int Lx,int Ly,int Lz,
	       int idx,int idy,int myid,int numx,int numy,float *p,float *vx,float *vy,
           float *vz,float *fi_px,float *fi_py,float *fi_pz,float *fi_vx,
	       float *fi_vy,float *fi_vz,int numprocs,int offsetrowxn,int offsetrowyn,char filepathr[40])
{
	int ix,iy,iz,t,ii,jj,kk,is,iss,flag1;
	float temp,temp1,temp2,temp3,temp4,t0=1.0/fm;
	float sumx,sumy,sumz;
	int index;

	MPI_Comm comm=MPI_COMM_WORLD;
	MPI_Info info=MPI_INFO_NULL;
	int amodew=MPI_MODE_CREATE+MPI_MODE_WRONLY;
	MPI_Offset offsetrec,offsetsnap;

	MPI_File fs;
	MPI_Request reqx[4],reqy[4];
    MPI_Status statusqx[4],statusqy[4]; 

	char ascx[6],ascy[6],ascn[6],ascp[6],asct[6],tmp1[30],tmp2[30];
	num2str(ascn,n);
	num2str(ascx,idx);num2str(ascy,idy);

	float *rec;
	rec=(float*)malloc(sizeof(float)*xef*yef);
	

	// define nearby process
	int tagx1=1,tagx2=2,tagx3=3,tagx4=4;
    int tagy1=5,tagy2=6,tagy3=7,tagy4=8;
    int up,down,left,right;
       
    if(idy>0){up=myid-numx;}
    else{up=MPI_PROC_NULL;}
       
    if(idx>0){left=myid-1;}
    else{left=MPI_PROC_NULL;}
       
    if(idy<numy-1){down=myid+numx;}
    else{down=MPI_PROC_NULL;}
       
    if(idx<numx-1){right=myid+1;}
    else{right=MPI_PROC_NULL;}


	float *sareaup,*sareadw,*sarealf,*sarearg;
	float *rareaup,*rareadw,*rarealf,*rarearg;
	sareaup=(float*)malloc(sizeof(float)*insizex*Ly*nz);
	sareadw=(float*)malloc(sizeof(float)*insizex*Ly*nz);
	rareaup=(float*)malloc(sizeof(float)*insizex*Ly*nz);
	rareadw=(float*)malloc(sizeof(float)*insizex*Ly*nz);

    sarealf=(float*)malloc(sizeof(float)*Lx*insizey*nz);
	sarearg=(float*)malloc(sizeof(float)*Lx*insizey*nz);
	rarealf=(float*)malloc(sizeof(float)*Lx*insizey*nz);
	rarearg=(float*)malloc(sizeof(float)*Lx*insizey*nz);
    
	
	// initialize the repeated non-blocking communication
	MPI_Send_init(sareaup,insizex*Ly*nz,MPI_FLOAT,up,tagy1,comm,&reqy[0]);
	MPI_Recv_init(rareadw,insizex*Ly*nz,MPI_FLOAT,down,tagy1,comm,&reqy[1]);
	MPI_Send_init(sareadw,insizex*Ly*nz,MPI_FLOAT,down,tagy2,comm,&reqy[2]);
	MPI_Recv_init(rareaup,insizex*Ly*nz,MPI_FLOAT,up,tagy2,comm,&reqy[3]);

	MPI_Send_init(sarealf,Lx*insizey*nz,MPI_FLOAT,left,tagx1,comm,&reqx[0]);
	MPI_Recv_init(rarearg,Lx*insizey*nz,MPI_FLOAT,right,tagx1,comm,&reqx[1]);
	MPI_Send_init(sarearg,Lx*insizey*nz,MPI_FLOAT,right,tagx2,comm,&reqx[2]);
	MPI_Recv_init(rarealf,Lx*insizey*nz,MPI_FLOAT,left,tagx2,comm,&reqx[3]);

    int pmlz=pml+Lz;
 
	  // set OpenMp threads 
    omp_set_num_threads (100); 
        
    #pragma omp parallel for private(ix,iy,iz,is)  
	for(iz=0;iz<nz;iz++)
	{
		for(iy=0;iy<sizey;iy++)
		{
			for(ix=0;ix<sizex;ix++)
			{
                is=iz*sizey*sizex+iy*sizex+ix;
				p[is]=0.0;vx[is]=0.0;vy[is]=0.0;vz[is]=0.0;
                fi_vx[is]=0.0;fi_vy[is]=0.0;fi_vz[is]=0.0;
				fi_px[is]=0.0;fi_py[is]=0.0;fi_pz[is]=0.0;
			}

		}
	}
        
	#pragma omp parallel for private(ix,iy,iz) 
	for(iz=0;iz<nz;iz++)
	{
	    for(iy=0;iy<Ly;iy++)
	    {
	        for(ix=0;ix<insizex;ix++)
	        {
	            rareaup[iz*Ly*insizex+iy*insizex+ix]=0.0;
	            rareadw[iz*Ly*insizex+iy*insizex+ix]=0.0;           
	        }
	    }
	}
	#pragma omp parallel for private(ix,iy,iz) 
	for(iz=0;iz<nz;iz++)
	{
	    for(iy=0;iy<insizey;iy++)
	    {
	        for(ix=0;ix<Lx;ix++)
	        {
	            rarealf[iz*Lx*insizey+iy*Lx+ix]=0.0;
	            rarearg[iz*Lx*insizey+iy*Lx+ix]=0.0;             
	        }
	    }
	}

    char recname[100];
	strcpy(recname,"");
	strcpy(recname,filepathr);
	sprintf(tmp1,"%s-rec_acou2m4-x%s-y%s.dat",ascn,ascx,ascy);
	strcat(recname,tmp1);
	// open files for recording for each area
	FILE *frec;
	frec=fopen(recname,"wb");

	
	int flagout=0,flagin=0,flagoutsnap=0;
	if(shotxlc[n]>=2*Lx&&shotxlc[n]<sizex-2*Lx&&shotylc[n]>=2*Ly&&shotylc[n]<sizey-2*Ly)
        {
	   flagin=1;flagout=0;
	    }
	else
	{
	   flagin=0;flagout=1;
	}
	flagout=flagout*flagshots[n];
	flagin=flagin*flagshots[n];

    int shotpos=(shotzlc[n])*sizey*sizex+shotylc[n]*sizex+shotxlc[n];	
    int nt=0;
	for(t=0;t<tmax;t++)     ///// time iteration starts
	{
        
		if(myid==0&&t%10==0)
		{
		    printf("%d\n",t);
		}

		
		  // (0) output recording (according to areas)
		  for(iy=0;iy<yef;iy++)
		  {
			for(ix=0;ix<xef;ix++)
			{
			     rec[iy*xef+ix]=p[(recg+zu)*sizey*sizex+(iy+yl)*sizex+ix+xl];
			}
		  }
		  fwrite(rec,sizeof(float),xef*yef,frec);

		//calculate the send areas of p
		// (1) up-down-y
		#pragma omp parallel for private(ix,iz,iy,sumx,sumy,sumz,ii,is,index) 
		
		for(iz=Lz;iz<nz-Lz;iz++)
		{
			for(ix=Lx;ix<sizex-Lx;ix++)
			{
			     for(iy=Ly;iy<2*Ly;iy++)   //
			        {
                     is=iz*sizey*sizex+iy*sizex+ix;
                     index=(int)((vp[is]-vpmin)/dv+0.5); 
                     sumx=0.0;sumy=0.0;sumz=0.0;
					 for(ii=0;ii<Lx;ii++)
					 {
				      sumx=sumx+Gp[index*Lpx+ii]*(vx[iz*sizey*sizex+iy*sizex+(ix+ii)]-vx[iz*sizey*sizex+iy*sizex+(ix-ii-1)]);
					 }
					  sumx=sumx+Gp[index*Lpx+Lx]*(vx[iz*sizey*sizex+(iy+1)*sizex+ix]-vx[iz*sizey*sizex+(iy+1)*sizex+(ix-1)]+vx[iz*sizey*sizex+(iy-1)*sizex+ix]-vx[iz*sizey*sizex+(iy-1)*sizex+(ix-1)]);
					  sumx=sumx+Gp[index*Lpx+Lx+1]*(vx[(iz+1)*sizey*sizex+iy*sizex+ix]-vx[(iz+1)*sizey*sizex+iy*sizex+(ix-1)]+vx[(iz-1)*sizey*sizex+iy*sizex+ix]-vx[(iz-1)*sizey*sizex+iy*sizex+(ix-1)]);					 
					 sumx=sumx/dx;
					 for(ii=0;ii<Ly;ii++)
					 {
					  sumy=sumy+Gp[index*Lpy+ii]*(vy[iz*sizey*sizex+(iy+ii)*sizex+ix]-vy[iz*sizey*sizex+(iy-ii-1)*sizex+ix]);
					 }
					  sumy=sumy+Gp[index*Lpy+Ly]*(vy[iz*sizey*sizex+iy*sizex+(ix+1)]-vy[iz*sizey*sizex+(iy-1)*sizex+(ix+1)]+vy[iz*sizey*sizex+iy*sizex+(ix-1)]-vy[iz*sizey*sizex+(iy-1)*sizex+(ix-1)]);
					  sumy=sumy+Gp[index*Lpy+Ly+1]*(vy[(iz+1)*sizey*sizex+iy*sizex+ix]-vy[(iz+1)*sizey*sizex+(iy-1)*sizex+ix]+vy[(iz-1)*sizey*sizex+iy*sizex+ix]-vy[(iz-1)*sizey*sizex+(iy-1)*sizex+ix]);					 
					 sumy=sumy/dy;
					 for(ii=0;ii<Lz;ii++)
					 {
				     sumz=sumz+Gp[index*Lpz+ii]*(vz[(iz+ii)*sizey*sizex+iy*sizex+ix]-vz[(iz-ii-1)*sizey*sizex+iy*sizex+ix]);
					 }
					 sumz=sumz/dz;

					 
					 fi_vx[is]=i_bx[ix-Lx]*fi_vx[is]+i_ax[ix-Lx]*sumx;
					 fi_vy[is]=i_by[iy-Ly]*fi_vy[is]+i_ay[iy-Ly]*sumy;
					 fi_vz[is]=i_b_z[iz]*fi_vz[is]+i_a_z[iz]*sumz;
                     p[is]=p[is]-dens[is]*vp[is]*vp[is]*dt*(fi_vx[is]+fi_vy[is]+fi_vz[is]+sumx+sumy+sumz);
				 }
				 
				 for(iy=sizey-2*Ly;iy<sizey-Ly;iy++)
				 {
					 is=iz*sizey*sizex+iy*sizex+ix;
					 index=(int)((vp[is]-vpmin)/dv+0.5); 
				     sumx=0.0;sumy=0.0;sumz=0.0;
					 for(ii=0;ii<Lx;ii++)
					 {
				      sumx=sumx+Gp[index*Lpx+ii]*(vx[iz*sizey*sizex+iy*sizex+(ix+ii)]-vx[iz*sizey*sizex+iy*sizex+(ix-ii-1)]);
					 }
					  sumx=sumx+Gp[index*Lpx+Lx]*(vx[iz*sizey*sizex+(iy+1)*sizex+ix]-vx[iz*sizey*sizex+(iy+1)*sizex+(ix-1)]+vx[iz*sizey*sizex+(iy-1)*sizex+ix]-vx[iz*sizey*sizex+(iy-1)*sizex+(ix-1)]);
					  sumx=sumx+Gp[index*Lpx+Lx+1]*(vx[(iz+1)*sizey*sizex+iy*sizex+ix]-vx[(iz+1)*sizey*sizex+iy*sizex+(ix-1)]+vx[(iz-1)*sizey*sizex+iy*sizex+ix]-vx[(iz-1)*sizey*sizex+iy*sizex+(ix-1)]);					 
					 sumx=sumx/dx;
					 for(ii=0;ii<Ly;ii++)
					 {
					  sumy=sumy+Gp[index*Lpy+ii]*(vy[iz*sizey*sizex+(iy+ii)*sizex+ix]-vy[iz*sizey*sizex+(iy-ii-1)*sizex+ix]);
					 }
					  sumy=sumy+Gp[index*Lpy+Ly]*(vy[iz*sizey*sizex+iy*sizex+(ix+1)]-vy[iz*sizey*sizex+(iy-1)*sizex+(ix+1)]+vy[iz*sizey*sizex+iy*sizex+(ix-1)]-vy[iz*sizey*sizex+(iy-1)*sizex+(ix-1)]);
					  sumy=sumy+Gp[index*Lpy+Ly+1]*(vy[(iz+1)*sizey*sizex+iy*sizex+ix]-vy[(iz+1)*sizey*sizex+(iy-1)*sizex+ix]+vy[(iz-1)*sizey*sizex+iy*sizex+ix]-vy[(iz-1)*sizey*sizex+(iy-1)*sizex+ix]);					 
					 sumy=sumy/dy;
					 for(ii=0;ii<Lz;ii++)
					 {
				      sumz=sumz+Gp[index*Lpz+ii]*(vz[(iz+ii)*sizey*sizex+iy*sizex+ix]-vz[(iz-ii-1)*sizey*sizex+iy*sizex+ix]);
					 }
					 sumz=sumz/dz;

					 
					 fi_vx[is]=i_bx[ix-Lx]*fi_vx[is]+i_ax[ix-Lx]*sumx;
					 fi_vy[is]=i_by[iy-Ly]*fi_vy[is]+i_ay[iy-Ly]*sumy;
					 fi_vz[is]=i_b_z[iz]*fi_vz[is]+i_a_z[iz]*sumz;
					 p[is]=p[is]-dens[is]*vp[is]*vp[is]*dt*(fi_vx[is]+fi_vy[is]+fi_vz[is]+sumx+sumy+sumz);

				 }
			}
                  		
		}
                
		// (2) left-right-x
		#pragma omp parallel for private(ix,iy,iz,sumx,sumy,sumz,ii,is,index) 
		
		for(iz=Lz;iz<nz-Lz;iz++)
		{
		     for(iy=2*Ly;iy<sizey-2*Ly;iy++) //attention
		     {
			     for(ix=Lx;ix<2*Lx;ix++)
			     {
					 is=iz*sizey*sizex+iy*sizex+ix;
					 index=(int)((vp[is]-vpmin)/dv+0.5); 
					 sumx=0.0;sumy=0.0;sumz=0.0;
					 for(ii=0;ii<Lx;ii++)
					 {
					 sumx=sumx+Gp[index*Lpx+ii]*(vx[iz*sizey*sizex+iy*sizex+(ix+ii)]-vx[iz*sizey*sizex+iy*sizex+(ix-ii-1)]);
					 }
					 sumx=sumx+Gp[index*Lpx+Lx]*(vx[iz*sizey*sizex+(iy+1)*sizex+ix]-vx[iz*sizey*sizex+(iy+1)*sizex+(ix-1)]+vx[iz*sizey*sizex+(iy-1)*sizex+ix]-vx[iz*sizey*sizex+(iy-1)*sizex+(ix-1)]);
					 sumx=sumx+Gp[index*Lpx+Lx+1]*(vx[(iz+1)*sizey*sizex+iy*sizex+ix]-vx[(iz+1)*sizey*sizex+iy*sizex+(ix-1)]+vx[(iz-1)*sizey*sizex+iy*sizex+ix]-vx[(iz-1)*sizey*sizex+iy*sizex+(ix-1)]);					 
					 sumx=sumx/dx;
					 for(ii=0;ii<Ly;ii++)
					 {
					  sumy=sumy+Gp[index*Lpy+ii]*(vy[iz*sizey*sizex+(iy+ii)*sizex+ix]-vy[iz*sizey*sizex+(iy-ii-1)*sizex+ix]);
					 }
					  sumy=sumy+Gp[index*Lpy+Ly]*(vy[iz*sizey*sizex+iy*sizex+(ix+1)]-vy[iz*sizey*sizex+(iy-1)*sizex+(ix+1)]+vy[iz*sizey*sizex+iy*sizex+(ix-1)]-vy[iz*sizey*sizex+(iy-1)*sizex+(ix-1)]);
					  sumy=sumy+Gp[index*Lpy+Ly+1]*(vy[(iz+1)*sizey*sizex+iy*sizex+ix]-vy[(iz+1)*sizey*sizex+(iy-1)*sizex+ix]+vy[(iz-1)*sizey*sizex+iy*sizex+ix]-vy[(iz-1)*sizey*sizex+(iy-1)*sizex+ix]);					 
					  sumy=sumy/dy;
					 for(ii=0;ii<Lz;ii++)
					 {
				      sumz=sumz+Gp[index*Lpz+ii]*(vz[(iz+ii)*sizey*sizex+iy*sizex+ix]-vz[(iz-ii-1)*sizey*sizex+iy*sizex+ix]);
					 }
					 sumz=sumz/dz;

					 
					 fi_vx[is]=i_bx[ix-Lx]*fi_vx[is]+i_ax[ix-Lx]*sumx;
					 fi_vy[is]=i_by[iy-Ly]*fi_vy[is]+i_ay[iy-Ly]*sumy;
					 fi_vz[is]=i_b_z[iz]*fi_vz[is]+i_a_z[iz]*sumz;
					 p[is]=p[is]-dens[is]*vp[is]*vp[is]*dt*(fi_vx[is]+fi_vy[is]+fi_vz[is]+sumx+sumy+sumz);
			     }
				 
			     for(ix=sizex-2*Lx;ix<sizex-Lx;ix++)
			     {
					 is=iz*sizey*sizex+iy*sizex+ix;
					 index=(int)((vp[is]-vpmin)/dv+0.5); 
				     sumx=0.0;sumy=0.0;sumz=0.0;
					 for(ii=0;ii<Lx;ii++)
					 {
					 sumx=sumx+Gp[index*Lpx+ii]*(vx[iz*sizey*sizex+iy*sizex+(ix+ii)]-vx[iz*sizey*sizex+iy*sizex+(ix-ii-1)]);
					 }
					 sumx=sumx+Gp[index*Lpx+Lx]*(vx[iz*sizey*sizex+(iy+1)*sizex+ix]-vx[iz*sizey*sizex+(iy+1)*sizex+(ix-1)]+vx[iz*sizey*sizex+(iy-1)*sizex+ix]-vx[iz*sizey*sizex+(iy-1)*sizex+(ix-1)]);
					 sumx=sumx+Gp[index*Lpx+Lx+1]*(vx[(iz+1)*sizey*sizex+iy*sizex+ix]-vx[(iz+1)*sizey*sizex+iy*sizex+(ix-1)]+vx[(iz-1)*sizey*sizex+iy*sizex+ix]-vx[(iz-1)*sizey*sizex+iy*sizex+(ix-1)]);					 
					 sumx=sumx/dx;
					 for(ii=0;ii<Ly;ii++)
					 {
					  sumy=sumy+Gp[index*Lpy+ii]*(vy[iz*sizey*sizex+(iy+ii)*sizex+ix]-vy[iz*sizey*sizex+(iy-ii-1)*sizex+ix]);
					 }
					  sumy=sumy+Gp[index*Lpy+Ly]*(vy[iz*sizey*sizex+iy*sizex+(ix+1)]-vy[iz*sizey*sizex+(iy-1)*sizex+(ix+1)]+vy[iz*sizey*sizex+iy*sizex+(ix-1)]-vy[iz*sizey*sizex+(iy-1)*sizex+(ix-1)]);
					  sumy=sumy+Gp[index*Lpy+Ly+1]*(vy[(iz+1)*sizey*sizex+iy*sizex+ix]-vy[(iz+1)*sizey*sizex+(iy-1)*sizex+ix]+vy[(iz-1)*sizey*sizex+iy*sizex+ix]-vy[(iz-1)*sizey*sizex+(iy-1)*sizex+ix]);					 
					 sumy=sumy/dy;
					 for(ii=0;ii<Lz;ii++)
					 {
					  sumz=sumz+Gp[index*Lpz+ii]*(vz[(iz+ii)*sizey*sizex+iy*sizex+ix]-vz[(iz-ii-1)*sizey*sizex+iy*sizex+ix]);
					 }
					 sumz=sumz/dz;

					 
					 
					 fi_vx[is]=i_bx[ix-Lx]*fi_vx[is]+i_ax[ix-Lx]*sumx;
					 fi_vy[is]=i_by[iy-Ly]*fi_vy[is]+i_ay[iy-Ly]*sumy;
					 fi_vz[is]=i_b_z[iz]*fi_vz[is]+i_a_z[iz]*sumz;
					 p[is]=p[is]-dens[is]*vp[is]*vp[is]*dt*(fi_vx[is]+fi_vy[is]+fi_vz[is]+sumx+sumy+sumz); 
			     }
				 
		      }	
		}
            
		
		// (3) add excitation or not
		temp=(pi*fm*(t*dt-t0))*(pi*fm*(t*dt-t0));
	    p[shotpos]=p[shotpos]+flagout*(1.0-2.0*temp)*exp(-temp);
		
		#pragma omp parallel for private(iz,iy,ix,iss)
		for(iz=0;iz<nz;iz++)
		{
			for(iy=0;iy<Ly;iy++)
			{
				for(ix=0;ix<insizex;ix++)
				{
					iss=iz*Ly*insizex+iy*insizex+ix;
					sareaup[iss]=p[iz*sizey*sizex+(iy+Ly)*sizex+(ix+Lx)];
					sareadw[iss]=p[iz*sizey*sizex+(iy+sizey-2*Ly)*sizex+(ix+Lx)];
				}
			}
		}
		#pragma omp parallel for private(iz,iy,ix,iss)
		for(iz=0;iz<nz;iz++)
		{
			for(iy=0;iy<insizey;iy++)
			{
				for(ix=0;ix<Lx;ix++)   
				{
					iss=iz*insizey*Lx+iy*Lx+ix;
					sarealf[iss]=p[iz*sizey*sizex+(iy+Ly)*sizex+(ix+Lx)];
					sarearg[iss]=p[iz*sizey*sizex+(iy+Ly)*sizex+(ix+sizex-2*Lx)]; 
				}
			}
		}

		// start non-blocking send and receive
		MPI_Startall(4,reqx);
		MPI_Startall(4,reqy);
		
		// (4) calculate the rest p areas
		#pragma omp parallel for private(ix,iz,iy,sumx,sumy,sumz,ii,is,index) 
		
		for(iz=Lz;iz<nz-Lz;iz++)
		{
			for(iy=2*Ly;iy<sizey-2*Ly;iy++)
			{
				for(ix=2*Lx;ix<sizex-2*Lx;ix++)
				{
					 is=iz*sizey*sizex+iy*sizex+ix;
					 index=(int)((vp[is]-vpmin)/dv+0.5); 
					 sumx=0.0;sumy=0.0;sumz=0.0;
					 for(ii=0;ii<Lx;ii++)
					 {
				      sumx=sumx+Gp[index*Lpx+ii]*(vx[iz*sizey*sizex+iy*sizex+(ix+ii)]-vx[iz*sizey*sizex+iy*sizex+(ix-ii-1)]);
					 }
					  sumx=sumx+Gp[index*Lpx+Lx]*(vx[iz*sizey*sizex+(iy+1)*sizex+ix]-vx[iz*sizey*sizex+(iy+1)*sizex+(ix-1)]+vx[iz*sizey*sizex+(iy-1)*sizex+ix]-vx[iz*sizey*sizex+(iy-1)*sizex+(ix-1)]);
					  sumx=sumx+Gp[index*Lpx+Lx+1]*(vx[(iz+1)*sizey*sizex+iy*sizex+ix]-vx[(iz+1)*sizey*sizex+iy*sizex+(ix-1)]+vx[(iz-1)*sizey*sizex+iy*sizex+ix]-vx[(iz-1)*sizey*sizex+iy*sizex+(ix-1)]);					 
					 sumx=sumx/dx;
					 for(ii=0;ii<Ly;ii++)
					 {
					   sumy=sumy+Gp[index*Lpy+ii]*(vy[iz*sizey*sizex+(iy+ii)*sizex+ix]-vy[iz*sizey*sizex+(iy-ii-1)*sizex+ix]);
					 }
					   sumy=sumy+Gp[index*Lpy+Ly]*(vy[iz*sizey*sizex+iy*sizex+(ix+1)]-vy[iz*sizey*sizex+(iy-1)*sizex+(ix+1)]+vy[iz*sizey*sizex+iy*sizex+(ix-1)]-vy[iz*sizey*sizex+(iy-1)*sizex+(ix-1)]);
					   sumy=sumy+Gp[index*Lpy+Ly+1]*(vy[(iz+1)*sizey*sizex+iy*sizex+ix]-vy[(iz+1)*sizey*sizex+(iy-1)*sizex+ix]+vy[(iz-1)*sizey*sizex+iy*sizex+ix]-vy[(iz-1)*sizey*sizex+(iy-1)*sizex+ix]);					 
					 sumy=sumy/dy;
					 for(ii=0;ii<Lz;ii++)
					 {
					  sumz=sumz+Gp[index*Lpz+ii]*(vz[(iz+ii)*sizey*sizex+iy*sizex+ix]-vz[(iz-ii-1)*sizey*sizex+iy*sizex+ix]);
					 }
					 sumz=sumz/dz;

					 
					 fi_vx[is]=i_bx[ix-Lx]*fi_vx[is]+i_ax[ix-Lx]*sumx;
					 fi_vy[is]=i_by[iy-Ly]*fi_vy[is]+i_ay[iy-Ly]*sumy;
					 fi_vz[is]=i_b_z[iz]*fi_vz[is]+i_a_z[iz]*sumz;
					 p[is]=p[is]-dens[is]*vp[is]*vp[is]*dt*(fi_vx[is]+fi_vy[is]+fi_vz[is]+sumx+sumy+sumz);
				}
                       
			}
		}
		// adding excitation source
	        p[shotpos]=p[shotpos]+flagin*(1.0-2.0*temp)*exp(-temp);

		// completes the non-blocking send and receive
		MPI_Waitall(4,reqx,statusqx);
		MPI_Waitall(4,reqy,statusqy);

		// update the p values of receiving areas
		//(1) up-down areas
		#pragma omp parallel for private(iz,iy,ix,iss) 
		for(iz=0;iz<nz;iz++)
		{
			for(iy=0;iy<Ly;iy++)
			{
				for(ix=0;ix<insizex;ix++)
				{
					iss=iz*Ly*insizex+iy*insizex+ix;
					p[iz*sizey*sizex+iy*sizex+(ix+Lx)]=rareaup[iss];
					p[iz*sizey*sizex+(iy+sizey-Ly)*sizex+(ix+Lx)]=rareadw[iss];
				}
			}
		}
		//(2) left-right areas
		#pragma omp parallel for private(iz,iy,ix,iss) 
		for(iz=0;iz<nz;iz++)  
		{
			for(iy=0;iy<insizey;iy++)
			{
				for(ix=0;ix<Lx;ix++) 
				{
					iss=iz*insizey*Lx+iy*Lx+ix;
					p[iz*sizey*sizex+(iy+Ly)*sizex+ix]=rarealf[iss];
					p[iz*sizey*sizex+(iy+Ly)*sizex+(ix+sizex-Lx)]=rarearg[iss]; 
				}
			}
		}


		// calculate velocity components
		//(1) calculate vx send areas: left and right
		#pragma omp parallel for private(ix,iy,iz,sumx,ii,is,index)
		
		for(iz=Lz;iz<nz-Lz;iz++)
		{
			for(iy=Ly;iy<sizey-Ly;iy++)
			{
				for(ix=Lx;ix<2*Lx;ix++)
				{
					is=iz*sizey*sizex+iy*sizex+ix;
					index=(int)((vp[is]-vpmin)/dv+0.5); 
					sumx=0.0;
					for(ii=0;ii<Lx;ii++)
					{
					  sumx=sumx+Gp[index*Lpx+ii]*(p[iz*sizey*sizex+iy*sizex+(ix+ii+1)]-p[iz*sizey*sizex+iy*sizex+(ix-ii)]);
					}
					 sumx=sumx+Gp[index*Lpx+Lx]*(p[iz*sizey*sizex+(iy+1)*sizex+(ix+1)]-p[iz*sizey*sizex+(iy+1)*sizex+ix]+p[iz*sizey*sizex+(iy-1)*sizex+(ix+1)]-p[iz*sizey*sizex+(iy-1)*sizex+ix]);
					 sumx=sumx+Gp[index*Lpx+Lx+1]*(p[(iz+1)*sizey*sizex+iy*sizex+ix+1]-p[(iz+1)*sizey*sizex+iy*sizex+ix]+p[(iz-1)*sizey*sizex+iy*sizex+ix+1]-p[(iz-1)*sizey*sizex+iy*sizex+ix]);	
					sumx=sumx/dx;
					
					fi_px[is]=h_bx[ix-Lx]*fi_px[is]+h_ax[ix-Lx]*sumx;
					vx[is]=vx[is]-dt/dens[is]*(sumx+fi_px[is]);
				}
				for(ix=sizex-2*Lx;ix<sizex-Lx;ix++)
				{
					is=iz*sizey*sizex+iy*sizex+ix;
					index=(int)((vp[is]-vpmin)/dv+0.5); 
					sumx=0.0;
					for(ii=0;ii<Lx;ii++)
					{
					 sumx=sumx+Gp[index*Lpx+ii]*(p[iz*sizey*sizex+iy*sizex+(ix+ii+1)]-p[iz*sizey*sizex+iy*sizex+(ix-ii)]);
					}
					 sumx=sumx+Gp[index*Lpx+Lx]*(p[iz*sizey*sizex+(iy+1)*sizex+(ix+1)]-p[iz*sizey*sizex+(iy+1)*sizex+ix]+p[iz*sizey*sizex+(iy-1)*sizex+(ix+1)]-p[iz*sizey*sizex+(iy-1)*sizex+ix]);
					 sumx=sumx+Gp[index*Lpx+Lx+1]*(p[(iz+1)*sizey*sizex+iy*sizex+ix+1]-p[(iz+1)*sizey*sizex+iy*sizex+ix]+p[(iz-1)*sizey*sizex+iy*sizex+ix+1]-p[(iz-1)*sizey*sizex+iy*sizex+ix]);						
					sumx=sumx/dx;
					
					fi_px[is]=h_bx[ix-Lx]*fi_px[is]+h_ax[ix-Lx]*sumx;
					vx[is]=vx[is]-dt/dens[is]*(sumx+fi_px[is]);
				}
				
			}
		}
	     	//(2) calculate vy send areas: up and down
	     	#pragma omp parallel for private(iz,iy,ix,sumy,ii,is,index)
	     	
		for(ix=Lx;ix<sizex-Lx;ix++)
		{
			for(iz=Lz;iz<nz-Lz;iz++)
			{
				for(iy=Ly;iy<2*Ly;iy++)
				{
					is=iz*sizey*sizex+iy*sizex+ix;
					index=(int)((vp[is]-vpmin)/dv+0.5); 
					sumy=0.0;
					for(ii=0;ii<Ly;ii++)
					{
					   sumy=sumy+Gp[index*Lpy+ii]*(p[iz*sizey*sizex+(iy+ii+1)*sizex+ix]-p[iz*sizey*sizex+(iy-ii)*sizex+ix]);
					}
					   sumy=sumy+Gp[index*Lpy+Ly]*(p[iz*sizey*sizex+(iy+1)*sizex+(ix+1)]-p[iz*sizey*sizex+iy*sizex+(ix+1)]+p[iz*sizey*sizex+(iy+1)*sizex+(ix-1)]-p[iz*sizey*sizex+iy*sizex+(ix-1)]);
					   sumy=sumy+Gp[index*Lpy+Ly+1]*(p[(iz+1)*sizey*sizex+(iy+1)*sizex+ix]-p[(iz+1)*sizey*sizex+iy*sizex+ix]+p[(iz-1)*sizey*sizex+(iy+1)*sizex+ix]-p[(iz-1)*sizey*sizex+iy*sizex+ix]);					
					sumy=sumy/dy;
					
					fi_py[is]=h_by[iy-Ly]*fi_py[is]+h_ay[iy-Ly]*sumy;
					vy[is]=vy[is]-dt/dens[is]*(sumy+fi_py[is]);
				}
				
				for(iy=sizey-2*Ly;iy<sizey-Ly;iy++)
				{
					is=iz*sizey*sizex+iy*sizex+ix;
					index=(int)((vp[is]-vpmin)/dv+0.5); 
					sumy=0.0;
					for(ii=0;ii<Ly;ii++)
					{
						sumy=sumy+Gp[index*Lpy+ii]*(p[iz*sizey*sizex+(iy+ii+1)*sizex+ix]-p[iz*sizey*sizex+(iy-ii)*sizex+ix]);
					}
					   sumy=sumy+Gp[index*Lpy+Ly]*(p[iz*sizey*sizex+(iy+1)*sizex+(ix+1)]-p[iz*sizey*sizex+iy*sizex+(ix+1)]+p[iz*sizey*sizex+(iy+1)*sizex+(ix-1)]-p[iz*sizey*sizex+iy*sizex+(ix-1)]);
					   sumy=sumy+Gp[index*Lpy+Ly+1]*(p[(iz+1)*sizey*sizex+(iy+1)*sizex+ix]-p[(iz+1)*sizey*sizex+iy*sizex+ix]+p[(iz-1)*sizey*sizex+(iy+1)*sizex+ix]-p[(iz-1)*sizey*sizex+iy*sizex+ix]);						
					sumy=sumy/dy;
					
					fi_py[is]=h_by[iy-Ly]*fi_py[is]+h_ay[iy-Ly]*sumy;
					vy[is]=vy[is]-dt/dens[is]*(sumy+fi_py[is]);
				}
			}
		}
                
		// filling the send areas
		#pragma omp parallel for private(iz,iy,ix,iss)
		for(iz=0;iz<nz;iz++)
		{
			for(iy=0;iy<Ly;iy++)
			{
				for(ix=0;ix<insizex;ix++)
				{
					iss=iz*Ly*insizex+iy*insizex+ix;
					sareaup[iss]=vy[iz*sizey*sizex+(iy+Ly)*sizex+(ix+Lx)];
					sareadw[iss]=vy[iz*sizey*sizex+(iy+sizey-2*Ly)*sizex+(ix+Lx)];
				}
			}
		}
		#pragma omp parallel for private(iz,iy,ix,iss)
		for(iz=0;iz<nz;iz++) 
		{
			for(iy=0;iy<insizey;iy++)
			{
				for(ix=0;ix<Lx;ix++)  
				{
					iss=iz*insizey*Lx+iy*Lx+ix;
					sarealf[iss]=vx[iz*sizey*sizex+(iy+Ly)*sizex+(ix+Lx)];
					sarearg[iss]=vx[iz*sizey*sizex+(iy+Ly)*sizex+(ix+sizex-2*Lx)]; 
				}
			}
		}

		// start non-blocking send and receive
		MPI_Startall(4,reqy);
		MPI_Startall(4,reqx);

		// calculate the rest vx
		#pragma omp parallel for private(ix,iz,iy,sumx,ii,is,index) 
		for(iz=Lz;iz<nz-Lz;iz++)
		{
			for(iy=Ly;iy<sizey-Ly;iy++)
			{
				for(ix=2*Lx;ix<sizex-2*Lx;ix++)
				{
					is=iz*sizey*sizex+iy*sizex+ix;
					index=(int)((vp[is]-vpmin)/dv+0.5); 
					sumx=0.0;
					for(ii=0;ii<Lx;ii++)
					{
					  sumx=sumx+Gp[index*Lpx+ii]*(p[iz*sizey*sizex+iy*sizex+(ix+ii+1)]-p[iz*sizey*sizex+iy*sizex+(ix-ii)]);
					}
					 sumx=sumx+Gp[index*Lpx+Lx]*(p[iz*sizey*sizex+(iy+1)*sizex+(ix+1)]-p[iz*sizey*sizex+(iy+1)*sizex+ix]+p[iz*sizey*sizex+(iy-1)*sizex+(ix+1)]-p[iz*sizey*sizex+(iy-1)*sizex+ix]);
					 sumx=sumx+Gp[index*Lpx+Lx+1]*(p[(iz+1)*sizey*sizex+iy*sizex+ix+1]-p[(iz+1)*sizey*sizex+iy*sizex+ix]+p[(iz-1)*sizey*sizex+iy*sizex+ix+1]-p[(iz-1)*sizey*sizex+iy*sizex+ix]);							
					sumx=sumx/dx;
					
					fi_px[is]=h_bx[ix-Lx]*fi_px[is]+h_ax[ix-Lx]*sumx;
					vx[is]=vx[is]-dt/dens[is]*(sumx+fi_px[is]);
				}
			}
		}

		// calculate the rest vy
		#pragma omp parallel for private(iz,iy,ix,sumy,ii,is,index) 
		for(iz=Lz;iz<nz-Lz;iz++)
		{
			for(iy=2*Ly;iy<sizey-2*Ly;iy++)  
			{
				for(ix=Lx;ix<sizex-Lx;ix++)
				{
					is=iz*sizey*sizex+iy*sizex+ix;
					index=(int)((vp[is]-vpmin)/dv+0.5); 
					sumy=0.0;
					for(ii=0;ii<Ly;ii++)
					{
					 sumy=sumy+Gp[index*Lpx+ii]*(p[iz*sizey*sizex+(iy+ii+1)*sizex+ix]-p[iz*sizey*sizex+(iy-ii)*sizex+ix]);
					}
					 sumy=sumy+Gp[index*Lpy+Ly]*(p[iz*sizey*sizex+(iy+1)*sizex+(ix+1)]-p[iz*sizey*sizex+iy*sizex+(ix+1)]+p[iz*sizey*sizex+(iy+1)*sizex+(ix-1)]-p[iz*sizey*sizex+iy*sizex+(ix-1)]);
					 sumy=sumy+Gp[index*Lpy+Ly+1]*(p[(iz+1)*sizey*sizex+(iy+1)*sizex+ix]-p[(iz+1)*sizey*sizex+iy*sizex+ix]+p[(iz-1)*sizey*sizex+(iy+1)*sizex+ix]-p[(iz-1)*sizey*sizex+iy*sizex+ix]);						
					sumy=sumy/dy;
					
					fi_py[is]=h_by[iy-Ly]*fi_py[is]+h_ay[iy-Ly]*sumy;
					vy[is]=vy[is]-dt/dens[is]*(sumy+fi_py[is]);
				}
			}
		}
		// calculate the whole vz
		#pragma omp parallel for private(iz,iy,ix,sumz,ii,is,index) 
		for(iz=Lz-1;iz<nz-Lz;iz++)
		{
			for(iy=Ly;iy<sizey-Ly;iy++)
			{
				for(ix=Lx;ix<sizex-Lx;ix++)
				{
					is=iz*sizey*sizex+iy*sizex+ix;
					index=(int)((vp[is]-vpmin)/dv+0.5); 
					sumz=0.0;
					for(ii=0;ii<Lz;ii++)
					{
					sumz=sumz+Gp[index*Lpz+ii]*(p[(iz+ii+1)*sizey*sizex+iy*sizex+ix]-p[(iz-ii)*sizey*sizex+iy*sizex+ix]);
					}
					sumz=sumz+Gp[index*Lpz+Lz]*(p[(iz+1)*sizey*sizex+(iy+1)*sizex+ix]-p[iz*sizey*sizex+(iy+1)*sizex+ix]+p[(iz+1)*sizey*sizex+(iy-1)*sizex+ix]-p[iz*sizey*sizex+(iy-1)*sizex+ix]);
					sumz=sumz+Gp[index*Lpz+Lz+1]*(p[(iz+1)*sizey*sizex+iy*sizex+ix+1]-p[iz*sizey*sizex+iy*sizex+(ix+1)]+p[(iz+1)*sizey*sizex+iy*sizex+ix-1]-p[iz*sizey*sizex+iy*sizex+(ix-1)]);						
					sumz=sumz/dz;
					
					fi_pz[is]=h_b_z[iz]*fi_pz[is]+h_a_z[iz]*sumz;
					vz[is]=vz[is]-dt/dens[is]*(sumz+fi_pz[is]);
				}
			}
		}
		// completes the non-blocking send and receive
		MPI_Waitall(4,reqx,statusqx);
		MPI_Waitall(4,reqy,statusqy);

                
		// update the send-receive area vx and vy
		#pragma omp parallel for private(iz,iy,ix,iss) 
		for(iz=0;iz<nz;iz++)
		{
			for(iy=0;iy<Ly;iy++)
			{
				for(ix=0;ix<insizex;ix++)
				{
					iss=iz*Ly*insizex+iy*insizex+ix;
					vy[iz*sizey*sizex+iy*sizex+(ix+Lx)]=rareaup[iss];
					vy[iz*sizey*sizex+(iy+sizey-Ly)*sizex+(ix+Lx)]=rareadw[iss];
				}
			}
		}
		#pragma omp parallel for private(iz,iy,ix,iss) 
		for(iz=0;iz<nz;iz++)
		{
			for(iy=0;iy<insizey;iy++)
			{
				for(ix=0;ix<Lx;ix++)   
				{
					iss=iz*insizey*Lx+iy*Lx+ix;
					vx[iz*sizey*sizex+(iy+Ly)*sizex+ix]=rarealf[iss];
					vx[iz*sizey*sizex+(iy+Ly)*sizex+(ix+sizex-Lx)]=rarearg[iss]; 
				}
			}
		}
                
		
	}// time iteration ends
	fclose(frec);

	free(rec);
	free(sareaup);free(sareadw);free(rareaup);free(rareadw);
	free(sarealf);free(sarearg);free(rarealf);free(rarearg);


	return;
}

// ==========================================================
//  This subroutine is used for input the geometry ...
//  =========================================================
void input_shotsposition(int shotnum,float dx,float dy,float dz,int *shotxg,
						 int *shotyg,int *shotzg,int *recg)
{
    int i;
    float temp1,temp2,temp3,temp4;
    char strtmp[256];
    FILE *fp=fopen("geometry.txt","r");
    if(fp==0)
    {
        printf("Cannot open the geometry1.txt file!\n");
        exit(0);
    }

    fgets(strtmp,256,fp);
    fscanf(fp,"\n");
    for(i=0;i<shotnum;i++)
    {
       fscanf(fp,"%f",&temp1);
       shotxg[i]=(int)(temp1/dx);
    } 
    fscanf(fp,"\n");
    
    fgets(strtmp,256,fp);
    fscanf(fp,"\n");
    for(i=0;i<shotnum;i++)
    {
       fscanf(fp,"%f",&temp2);
       shotyg[i]=(int)(temp2/dy);
    }

    fscanf(fp,"\n");
    fgets(strtmp,256,fp);
    fscanf(fp,"\n");
    for(i=0;i<shotnum;i++)
    {
       fscanf(fp,"%f",&temp3);
       shotzg[i]=(int)(temp3/dz);
    }

	fscanf(fp,"\n");
    fgets(strtmp,256,fp);
    fscanf(fp,"\n");
	fscanf(fp,"%f",&temp4);
	*recg=(int)(temp4/dz);
    fclose(fp);
    
    return;
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
	asc5='0'; // wan
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


// ==========================================================
//  This subroutine is used for input modeling parameters
//  =========================================================
void input_parameters(float *rectime,float *dt,float *fm,int *shotnum,float *dx,float *dy,float *dz,int *nxo,int *nyo,
                      int *nzo,int *Lx,int *Ly,int *Lz,int *numx,int *numy,char filepathv[40],char filepathr[40])
		      
                      
{
    char strtmp[256];
    FILE *fp=fopen("parameters3d.txt","r");
    if(fp==0)
    {
        printf("Cannot open the parameters3d4.txt file!\n");
        exit(0);
    }

    fgets(strtmp,256,fp);
    fscanf(fp,"\n");
    fscanf(fp,"%f",rectime);
    fscanf(fp,"\n");

    fgets(strtmp,256,fp); 
    fscanf(fp,"\n");
    fscanf(fp,"%f",dt);
    fscanf(fp,"\n");

    fgets(strtmp,256,fp); 
    fscanf(fp,"\n");
    fscanf(fp,"%f",fm);
    fscanf(fp,"\n");

    fgets(strtmp,256,fp);
    fscanf(fp,"\n");
    fscanf(fp,"%d",shotnum);
    fscanf(fp,"\n");
  
    fgets(strtmp,256,fp);
    fscanf(fp,"\n");
    fscanf(fp,"%f",dx);
    fscanf(fp,"\n");
  
    fgets(strtmp,256,fp);
    fscanf(fp,"\n");
    fscanf(fp,"%f",dy);
    fscanf(fp,"\n");
  
    fgets(strtmp,256,fp);
    fscanf(fp,"\n");
    fscanf(fp,"%f",dz);
    fscanf(fp,"\n");
  
    fgets(strtmp,256,fp);
    fscanf(fp,"\n");
    fscanf(fp,"%d",nxo);
    fscanf(fp,"\n");
  
    fgets(strtmp,256,fp);
    fscanf(fp,"\n");
    fscanf(fp,"%d",nyo);
    fscanf(fp,"\n");
  
    fgets(strtmp,256,fp);
    fscanf(fp,"\n");
    fscanf(fp,"%d",nzo);
    fscanf(fp,"\n");
  
    fgets(strtmp,256,fp);
    fscanf(fp,"\n");
    fscanf(fp,"%d",Lx);
    fscanf(fp,"\n");
  
    fgets(strtmp,256,fp);
    fscanf(fp,"\n");
    fscanf(fp,"%d",Ly);
    fscanf(fp,"\n");
  
    fgets(strtmp,256,fp);
    fscanf(fp,"\n");
    fscanf(fp,"%d",Lz);
    fscanf(fp,"\n");

    fgets(strtmp,256,fp);
    fscanf(fp,"\n");
    fscanf(fp,"%d",numx);
    fscanf(fp,"\n");

    fgets(strtmp,256,fp);
    fscanf(fp,"\n");
    fscanf(fp,"%d",numy);
    fscanf(fp,"\n");


    fgets(strtmp,256,fp);
    fscanf(fp,"\n");
    fscanf(fp,"%s",filepathr);
    fscanf(fp,"\n");

    
    fgets(strtmp,256,fp);
    fscanf(fp,"\n");
    fscanf(fp,"%s",filepathv);

    fclose(fp);
  
    return;
}

// ==========================================================
//  This subroutine is used for splitting the domain
//  =========================================================
void split_domains(int nx,int ny,int nxo,int nyo,int Lx,int Ly,int nz,int Lz,int myid,
                   int *avex,int *avey,int *insizex,int *insizey,int *modx,int *mody,int *sizex,int *sizey,
	               int *offsetrowx,int *offsetrowy,int *offsetrowxn,int *offsetrowyn,
		           int numx,int numy,int *idx,int *idy,int *xl,int *xr,int *yl,int *yr,
                   int *zu,int *zd,int *xef,int *yef)
{
      *modx=nx%numx;
	  *mody=ny%numy;
	  *avex=(nx-*modx)/numx;
	  *avey=(ny-*mody)/numy;

      *idy=myid/numx;
      *idx=myid-(*idy)*numx;

	  *insizex=*avex;
	  *insizey=*avey;

	  if(*idx<*modx)
          {
               *insizex=*insizex+1;
          }
          if(*idy<*mody)
          {
               *insizey=*insizey+1;
          }

	  *sizex=*insizex+2*Lx;
	  *sizey=*insizey+2*Ly;

	  if(*idx<*modx){
              *offsetrowx=(*idx)*(*avex+1);}
          else{
               *offsetrowx=(*modx)*(*avex+1)+(*idx-*modx)*(*avex);}
	 
          if(*idy<*mody){
               *offsetrowy=(*idy)*(*avey+1);}
          else{
               *offsetrowy=(*mody)*(*avey+1)+(*idy-*mody)*(*avey);}

	  *offsetrowxn=(*offsetrowx)-pml;
	  *offsetrowyn=(*offsetrowy)-pml;
	  if(*idx==0){*offsetrowxn=0;}
	  if(*idy==0){*offsetrowyn=0;}

	  *xl=Lx;*xr=(*sizex)-Lx-1;
	  *yl=Ly;*yr=(*sizey)-Ly-1;
	  *zu=pml+Lz;*zd=nz-pml-Lz-1;

	  if(*idx==0){*xl=*xl+pml;}
	  if(*idx==numx-1){*xr=*xr-pml;}
	  if(*idy==0){*yl=*yl+pml;}
	  if(*idy==numy-1){*yr=*yr-pml;}

	  *xef=*xr-*xl+1;
	  *yef=*yr-*yl+1;

	  return;
}

// ==========================================================
//  This subroutine is used for input velocity & density for sub-domains
//  =========================================================
void input_model_parameters(float *vp,float *dens,float *modelv,float *modeld,
			                int xef,int yef,int nzo,int xl,int xr,int yl,int yr,
		                    int zu,int zd,int nz,int sizex,int sizey)
{
        int i,j,k;
	for(i=xl;i<=xr;i++)
	{
		for(j=yl;j<=yr;j++)
		{
			for(k=zu;k<=zd;k++)
			{
				vp[k*sizey*sizex+j*sizex+i]=modelv[(k-zu)*yef*xef+(j-yl)*xef+i-xl];
				dens[k*sizey*sizex+j*sizex+i]=modeld[(k-zu)*yef*xef+(j-yl)*xef+i-xl];
			}
		}
	}

	for(i=xl;i<=xr;i++)
	{
	        for(j=yl;j<=yr;j++)
		{
			for(k=0;k<zu;k++)
			{
				vp[k*sizey*sizex+j*sizex+i]=vp[zu*sizey*sizex+j*sizex+i];
				dens[k*sizey*sizex+j*sizex+i]=dens[zu*sizey*sizex+j*sizex+i];
			}
			for(k=zd+1;k<nz;k++)
			{
				vp[k*sizey*sizex+j*sizex+i]=vp[zd*sizey*sizex+j*sizex+i];
				dens[k*sizey*sizex+j*sizex+i]=dens[zd*sizey*sizex+j*sizex+i];
			}
		}
		for(k=0;k<nz;k++)
		{
			for(j=0;j<yl;j++)
			{
				vp[k*sizey*sizex+j*sizex+i]=vp[k*sizey*sizex+yl*sizex+i];
				dens[k*sizey*sizex+j*sizex+i]=dens[k*sizey*sizex+yl*sizex+i];
			}

			for(j=yr+1;j<sizey;j++)
			{
				vp[k*sizey*sizex+j*sizex+i]=vp[k*sizey*sizex+yr*sizex+i];
				dens[k*sizey*sizex+j*sizex+i]=dens[k*sizey*sizex+yr*sizex+i];
			}
		}
	}
	for(i=0;i<xl;i++)
	{
		for(j=0;j<sizey;j++)
		{
			for(k=0;k<nz;k++)
			{
				vp[k*sizey*sizex+j*sizex+i]=vp[k*sizey*sizex+j*sizex+xl];
				dens[k*sizey*sizex+j*sizex+i]=dens[k*sizey*sizex+j*sizex+xl];
			}
		}
	}

	for(i=xr+1;i<sizex;i++)
	{
		for(j=0;j<sizey;j++)
		{
			for(k=0;k<nz;k++)
			{
				vp[k*sizey*sizex+j*sizex+i]=vp[k*sizey*sizex+j*sizex+xr];
				dens[k*sizey*sizex+j*sizex+i]=dens[k*sizey*sizex+j*sizex+xr];
			}
		}
	}

	return;
}

// ==========================================================
//  This subroutine is used for calculating PML parameters
//  =========================================================
void cal_pml_parameters(float *i_a_x,float *i_b_x,float *h_a_x,float *h_b_x,
                        float *i_a_z,float *i_b_z,float *h_a_z,float *h_b_z,
                        float *i_a_y,float *i_b_y,float *h_a_y,float *h_b_y,
                        int nx,int ny,int nz,int Lz,float dx,float dy,
			            float dz,float dt,float fm)
{     
        float alpha_max=pi*fm,dmaxx,dmaxz,dmaxy,widthx,widthz,widthy;
        float Re=1e-5,vmax=5000;
        int n1=2,n2=1,n3=1,i,j;
        float temp,temp1,temp2,temp3,temp4;
        
        float *int_dx,*int_dz,*int_alphax,*int_alphaz,*int_dy,*int_alphay;
        float *half_dx,*half_dz,*half_alphax,*half_alphaz,*half_dy,*half_alphay;
        
        int_dz=(float*)malloc(sizeof(float)*nz);
        half_dz=(float*)malloc(sizeof(float)*nz);
        int_alphaz=(float*)malloc(sizeof(float)*nz);
        half_alphaz=(float*)malloc(sizeof(float)*nz);
        
        int_dx=(float*)malloc(sizeof(float)*nx);
        half_dx=(float*)malloc(sizeof(float)*nx);
        int_alphax=(float*)malloc(sizeof(float)*nx);
        half_alphax=(float*)malloc(sizeof(float)*nx);

	    int_dy=(float*)malloc(sizeof(float)*ny);
        half_dy=(float*)malloc(sizeof(float)*ny);
        int_alphay=(float*)malloc(sizeof(float)*ny);
        half_alphay=(float*)malloc(sizeof(float)*ny);
        
        int pmlx=pml,pmly=pml;
        int pmlz=pml+Lz;
          
        widthx=pmlx*dx;widthz=pmlz*dz;widthy=pmly*dy;
        dmaxx=(1+n1+n2)*vmax*log(1.0/Re)/(2.0*widthx);
        dmaxz=(1+n1+n2)*vmax*log(1.0/Re)/(2.0*widthz);
	    dmaxy=(1+n1+n2)*vmax*log(1.0/Re)/(2.0*widthy);
     
        // integer absorbing parameter
        for(j=0;j<pmlz;j++)
        {
            temp1=pow(1.0*(pmlz-1-j)/pmlz,n1+n2);
            int_dz[j]=dmaxz*temp1;
         
            temp3=pow(1.0*j/(pmlz-1),n3);
            int_alphaz[j]=alpha_max*temp3;

            int_dz[nz-1-j]=int_dz[j];
            int_alphaz[nz-1-j]=int_alphaz[j];
        }
        for(j=pmlz;j<nz-pmlz;j++)
        {
            int_dz[j]=0.0;
            int_alphaz[j]=int_alphaz[pmlz-1];
        }
        
        for(i=0;i<pmlx;i++)
        {
            temp1=pow(1.0*(pmlx-1-i)/pmlx,n1+n2);
            int_dx[i]=dmaxx*temp1;
         
            temp3=pow(1.0*i/(pmlx-1),n3);
            int_alphax[i]=alpha_max*temp3;

            int_dx[nx-1-i]=int_dx[i];
            int_alphax[nx-1-i]=int_alphax[i];
        }
        for(i=pmlx;i<nx-pmlx;i++)
        {
            int_dx[i]=0.0;
            int_alphax[i]=int_alphax[pmlx-1];
        }

	for(i=0;i<pmly;i++)
        {
            temp1=pow(1.0*(pmly-1-i)/pmly,n1+n2);
            int_dy[i]=dmaxy*temp1;
         
            temp3=pow(1.0*i/(pmly-1),n3);
            int_alphay[i]=alpha_max*temp3;

            int_dy[ny-1-i]=int_dy[i];
            int_alphay[ny-1-i]=int_alphay[i];
        }
        for(i=pmly;i<ny-pmly;i++)
        {
            int_dy[i]=0.0;
            int_alphay[i]=int_alphay[pmly-1];
        }


        // half absorbing parameter
        for(j=0;j<pmlz-1;j++)
        {
            temp2=pow(1.0*(pmlz-1.5-j)/pmlz,n1+n2);
            half_dz[j]=dmaxz*temp2;
            half_dz[nz-2-j]=half_dz[j];
         
            temp4=pow(1.0*(j+0.5)/(pmlz-1),n3);
            half_alphaz[j]=alpha_max*temp4;
            half_alphaz[nz-2-j]=half_alphaz[j];
        }
        for(j=pmlz-1;j<nz-pmlz;j++)
        {
            half_dz[j]=0.0;
            half_alphaz[j]=half_alphaz[pmlz-2];;
        }
        half_dz[nz-1]=0.0;
        half_alphaz[nz-1]=half_alphaz[nz-2];
        
        
        for(i=0;i<pmlx-1;i++)
        {
            temp2=pow(1.0*(pmlx-1.5-i)/pmlx,n1+n2);
            half_dx[i]=dmaxx*temp2;
            half_dx[nx-2-i]=half_dx[i];
         
            temp4=pow(1.0*(i+0.5)/(pmlx-1),n3);
            half_alphax[i]=alpha_max*temp4;
            half_alphax[nx-2-i]=half_alphax[i];
        }
        for(i=pmlx-1;i<nx-pmlx;i++)
        {
            half_dx[i]=0.0;
            half_alphax[i]=half_alphax[pmlx-2];;
        }
        half_dx[nx-1]=0.0;
        half_alphax[nx-1]=half_alphax[nx-2];


	for(i=0;i<pmly-1;i++)
        {
            temp2=pow(1.0*(pmly-1.5-i)/pmly,n1+n2);
            half_dy[i]=dmaxy*temp2;
            half_dy[ny-2-i]=half_dy[i];
         
            temp4=pow(1.0*(i+0.5)/(pmly-1),n3);
            half_alphay[i]=alpha_max*temp4;
            half_alphay[ny-2-i]=half_alphay[i];
        }
        for(i=pmly-1;i<ny-pmly;i++)
        {
            half_dy[i]=0.0;
            half_alphay[i]=half_alphay[pmly-2];
        }
        half_dy[ny-1]=0.0;
        half_alphay[ny-1]=half_alphay[ny-2];
        
     
        for(j=0;j<nz;j++)
        {
            temp=int_dz[j]+int_alphaz[j];
            i_b_z[j]=exp(-dt*temp);
            i_a_z[j]=int_dz[j]/temp*(i_b_z[j]-1.0);
        
            temp=half_dz[j]+half_alphaz[j];
            h_b_z[j]=exp(-dt*temp);
            h_a_z[j]=half_dz[j]/temp*(h_b_z[j]-1.0);
        }
        
        for(i=0;i<nx;i++)
        {
            temp=int_dx[i]+int_alphax[i];
            i_b_x[i]=exp(-dt*temp);
            i_a_x[i]=int_dx[i]/temp*(i_b_x[i]-1.0);
        
            temp=half_dx[i]+half_alphax[i];
            h_b_x[i]=exp(-dt*temp);
            h_a_x[i]=half_dx[i]/temp*(h_b_x[i]-1.0);
        }
        
        for(i=0;i<ny;i++)
        {
            temp=int_dy[i]+int_alphay[i];
            i_b_y[i]=exp(-dt*temp);
            i_a_y[i]=int_dy[i]/temp*(i_b_y[i]-1.0);
        
            temp=half_dy[i]+half_alphay[i];
            h_b_y[i]=exp(-dt*temp);
            h_a_y[i]=half_dy[i]/temp*(h_b_y[i]-1.0);
        }

        free(int_dz);free(int_alphaz);free(half_dz);free(half_alphaz);
        free(int_dx);free(int_alphax);free(half_dx);free(half_alphax);
	free(int_dy);free(int_alphay);free(half_dy);free(half_alphay);
        
        return;
}

// ==========================================================
//  This subroutine is used for input velocity & density 
//  =========================================================
void input_v_dens(char filepathv[40],int xef,int yef,int zef,int idx,int idy,float *modelv,float *modeld)
{
     int i,j,k;
     char ascx[6],ascy[6];
     char name1[180],name2[180];
     char tmp1[40],tmp2[40];
     
     num2str(ascx,idx);num2str(ascy,idy);

     strcpy(name1,filepathv);
     strcpy(name2,filepathv);
     
    sprintf(tmp1,"x%s_vp_y%s.dat",ascx,ascy);
    sprintf(tmp2,"x%s_rho_y%s.dat",ascx,ascy);   
     
  
     strcat(name1,tmp1);
     strcat(name2,tmp2);
     
     FILE *fp1=fopen(name1,"rb");
     FILE *fp2=fopen(name2,"rb");
     
     fread(modelv,sizeof(float),xef*yef*zef,fp1);
     fread(modeld,sizeof(float),xef*yef*zef,fp2);
     
     fclose(fp1);fclose(fp2);
     
     return;
}

// ==========================================================
//  This subroutine is used for get the max & min velocity
//  =========================================================
void get_max_min(float *vpo,int N,float *vpmax,float *vpmin)
{
     int i;
     *vpmax=-1.0;
     *vpmin=8000;
     for(i=0;i<N;i++)
     {
         if(vpo[i]>*vpmax){*vpmax=vpo[i];}
         if(vpo[i]<*vpmin&&vpo[i]!=0.0){*vpmin=vpo[i];}
     }
     return;
}

// ==========================================================
//  This subroutine is used for calculate the (2m,4) FD coe
//  =========================================================
void TE_2M4_3dxyz(float vpmin,float dv,int N,int Lx,float dt,float dx,float dy,float dz,float *Gp)
{
     int i,ii,Lxo=Lx-2;
     float dtdx=dt/dx;
     float dtdy=dt/dy;
     float dtdz=dt/dz; 

     float *rp;
     float v,gamx,gamy,gamz;
     rp=(float*)malloc(sizeof(float)*Lxo);

     float temp1,temp2,sumx1;
     
     
     for(i=0;i<N;i++)
     {
         v=vpmin+i*dv;
         gamx=v*dtdx;
         gamy=v*dtdy;
         gamz=v*dtdz;
         

         cal_xishur(Lxo,rp,gamx);
         temp1=gamy*gamy/24.0;
         temp2=gamz*gamz/24.0;

         sumx1=0.0;
         for(ii=1;ii<Lxo;ii++)
         {
                sumx1=sumx1+(2*ii+1)*rp[ii];
         }
         Gp[i*Lx+0]=1.0-2.0*temp1-2.0*temp2-sumx1;    

         for(ii=1;ii<Lx-2;ii++)
         {
                Gp[i*Lx+ii]=rp[ii];
         }
         Gp[i*Lx+Lx-2]=temp1;
         Gp[i*Lx+Lx-1]=temp2;
         
     }

     free(rp);

     return;
}

void cal_xishur(int Lx,float *rx,float gam)
{
     int m,i;
     float s1,s2;
     for(m=1;m<=Lx;m++)
     {
        s1=1.0;s2=1.0;
        for(i=1;i<m;i++)
        {
            s1=s1*((2.0*i-1)*(2.0*i-1)-gam*gam);
            s2=s2*((2.0*m-1)*(2.0*m-1)-(2.0*i-1)*(2.0*i-1));
        }
        for(i=m+1;i<=Lx;i++)
        {
            s1=s1*((2.0*i-1)*(2.0*i-1)-gam*gam);
            s2=s2*((2.0*m-1)*(2.0*m-1)-(2.0*i-1)*(2.0*i-1));
        }
        s2=fabs(s2);

        rx[m-1]=pow(-1.0,m+1)*s1/(s2*(2.0*m-1));
     }
     
     return;
} 
