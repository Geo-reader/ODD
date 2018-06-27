void input_parameters(float *rectime,float *dt,float *fm,int *shotnum,float *dx,
		              float *dy,float *dz,int *nxo,int *nyo,int *nzo,int *Lx,
		              int *Ly,int *Lz,
		              int *numx,int *numy,char filepathv[40],
                      char filepathr[40]);

void split_domains(int nx,int ny,int nxo,int nyo,int Lx,int Ly,int nz,int Lz,int myid,
                   int *avex,int *avey,int *insizex,int *insizey,int *modx,int *mody,
		           int *sizex,int *sizey,int *offsetrowx,int *offsetrowy,int *offsetrowxn,
	               int *offsetrowyn,int numx,int numy,int *idx,int *idy,int *xl,
		           int *xr,int *yl,int *yr,
                   int *zu,int *zd,int *xef,int *yef);


void input_v_dens(char filepathv[40],int xef,int yef,int zef,int idx,int idy,
                  float *modelv,float *modeld);

void input_model_parameters(float *vp,float *dens,float *modelv,float *modeld,
			                int xef,int yef,int nzo,int xl,int xr,int yl,int yr,
		                    int zu,int zd,int nz,int sizex,int sizey);

void cal_pml_parameters(float *i_a_x,float *i_b_x,float *h_a_x,float *h_b_x,
                        float *i_a_z,float *i_b_z,float *h_a_z,float *h_b_z,
                        float *i_a_y,float *i_b_y,float *h_a_y,float *h_b_y,
                        int nx,int ny,int nz,int Lz,float dx,float dy,
			            float dz,float dt,float fm);

void num2str(char asc[6],int);

void cal_coef(int Lx,float *rx);

void input_shotsposition(int shotnum,float dx,float dy,float dz,int *shotxg,
			             int *shotyg,int *shotzg,int *recg);


void forward3d_c(int sizex,int sizey,int nz,int tmax,float dx,float dy,float dz,float dt,
	             float *rx,float *ry,float *rz,float *vp,float *dens,
	             int xl,int xr,int yl,int yr,int zu,int zd,int insizex,int insizey,
	             int nxo,int nyo,int nzo, int *shotxlc,int *shotylc,int *shotzlc,
                 int recg,int *flagshots,float *i_ax,float *i_bx,float *h_ax,float *h_bx,
	             float *i_ay,float *i_by,float *h_ay,float *h_by,
	             float *i_a_z,float *i_b_z,float *h_a_z,float *h_b_z,
	             float fm,int n,int xef,int yef,int zef,int Lx,int Ly,int Lz,
	             int idx,int idy,int myid,int numx,int numy,float *p,float *vx,float *vy,
                 float *vz,float *fi_px,float *fi_py,float *fi_pz,float *fi_vx,
	             float *fi_vy,float *fi_vz,int numprocs,
	             int offsetrowxn,int offsetrowyn,char filepathr[40]);

void check_stability(int sizex,int sizey,int nz,float dt,int Lx,int Ly,int Lz,int dx,int dy,int dz,float *rx,float *ry,float *rz,float *vp,float delta);
