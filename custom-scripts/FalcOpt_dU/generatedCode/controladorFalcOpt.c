/* Copyright (c) 2017, ETH Zurich, Automatic Control Laboratory */
/* Main developers: */
/* Giampaolo Torrisi <giampaolo.torrisi@gmail.com> */
/* Damian Frick <falcopt@damianfrick.com> */
/* Tommaso Robbiani <tommasro@student.ethz.ch> */
/* Scientic Contributors: */
/* Sergio Grammatico */
/* Roy S. Smith */
/* Manfred Morari */

#include "math.h"

#include "controladorFalcOpt.h"
/*Static data for Jacobian_x*/
static const double F[13] = {  0.7887,   0.282,   0.6439,   0.1823,   0.8486,   0.2155,   0.0973,   0.4836,   0.6808,   0.9657,   0.7654,   0.7362,   1.0};
/*Static data for Jacobian_u*/
static const double G[7] = {  0.6439,   0.1823,   0.0973,   0.4836,   0.9657,   0.7362,   1.0};
/* Static data for build_amu_1() */
static const double umin_1[1] = {
	-0.5 };
/* Static data for build_umb_1() */
static const double minus_umax_1[1] = {
	-0.5 };
/* Static data for build_vNnc_lb_1() */
static const double alpha_inv_lb_1[1] = {
	2.5 };
/* Static data for build_vNnc_ub_1() */
static const double alpha_inv_ub_1[1] = {
	2.5 };
/* Static data for scale_nu() */
static const double alpha_I[1] = {
	-0.4 };

double my_fmax(double x, double y){
	if (x < y) { return y; }
	else { return x; }
}

double my_fmin(double x, double y){
	if (x > y) { return y; }
	else { return x; }
}

/*Dynamics of the system*/
void model_mpc(const double* x, const double* u, double* xp){

	xp[0] =   u[0]*0.6439+u[1]*0.1823+x[0]*0.7887+x[2]*0.282+x[4]*0.6439+x[5]*0.1823;
	xp[1] =   u[0]*0.0973+u[1]*0.4836+x[1]*0.8486+x[3]*0.2155+x[4]*0.0973+x[5]*0.4836;
	xp[2] =   u[1]*0.9657+x[2]*0.6808+x[5]*0.9657;
	xp[3] =   u[0]*0.7362+x[3]*0.7654+x[4]*0.7362;
	xp[4] =   u[0]+x[4];
	xp[5] =   u[1]+x[5];
}

static void build_amu_1(double* r, const double* u) {
	r[0] = -u[0] + umin_1[0];
	r[1] = -u[1] + umin_1[0];
}


void build_amu(const double* u, const unsigned int k, double* amu){

	if((k >= 0) && (k <= 14)) {
		build_amu_1(&amu[0], &u[0]);
	}
}
static void build_umb_1(double* r, const double* u) {
	r[0] = u[0] + minus_umax_1[0];
	r[1] = u[1] + minus_umax_1[0];
}


void build_umb(const double* u, const unsigned int k, double* umb){

	if((k >= 0) && (k <= 14)) {
		build_umb_1(&umb[0], &u[0]);
	}
}

/* det_x is a forward simulation of model_mpc with initial state x0 
and sequence of predicted inputs u (of dim. N*nu) */
void det_x (const double* x0, const double* u, double* x){

	unsigned int ii = 0;

	model_mpc(x0,u,x);
	for (ii = 1;ii < 15; ++ii)
		model_mpc(x + (ii-1)* 6, u + ii* 2, x + ii* 6);

}

/* It computes Q*x */
static void Qmul(double* r, const double* dx) {
	r[0] = dx[0];
	r[1] = dx[1];
	r[2] = dx[2];
	r[3] = dx[3];
	r[4] = 0.0;
	r[5] = 0.0;
}


/* It computes P*x */
static void Pmul(double* r, const double* dx) {
	r[0] = dx[0];
	r[1] = dx[1];
	r[2] = dx[2];
	r[3] = dx[3];
	r[4] = 0.0;
	r[5] = 0.0;
}


/* It computes R*u */
static void Rmul(double* r, const double* du) {
	r[0] = du[0];
	r[1] = du[1];
}


/* dot product x^top *x */
static void dot_product_nx_nx(double* r, const double* du, const double* R) {
	r[0] = R[0]*du[0] + R[1]*du[1] + R[2]*du[2] + R[3]*du[3] + R[4]*du[4] + R[5]*du[5];
}


/* dot product u^top *u */
static void dot_product_nu_nu(double* r, const double* du, const double* R) {
	r[0] = R[0]*du[0] + R[1]*du[1];
}


/* It computes G^top * x + u (exploiting structure of G) */
static void product_and_sum_nu(double* r, const double* u1, const double* u2, const double* A) {
	r[0] = A[0]*u1[0] + A[2]*u1[1] + A[5]*u1[3] + A[6]*u1[4] + u2[0];
	r[1] = A[1]*u1[0] + A[3]*u1[1] + A[4]*u1[2] + A[6]*u1[5] + u2[1];
}


/* It computes F^top * u1 (exploiting structure of F) + u2 */
static void product_and_sum_nx(double* r, const double* u1, const double* u2, const double* A) {
	r[0] = A[0]*u1[0] + u2[0];
	r[1] = A[4]*u1[1] + u2[1];
	r[2] = A[1]*u1[0] + A[8]*u1[2] + u2[2];
	r[3] = A[5]*u1[1] + A[10]*u1[3] + u2[3];
	r[4] = A[2]*u1[0] + A[6]*u1[1] + A[11]*u1[3] + A[12]*u1[4] + u2[4];
	r[5] = A[3]*u1[0] + A[7]*u1[1] + A[9]*u1[2] + A[12]*u1[5] + u2[5];
}


/* It computes dx = x - xref */
static void diffX(double* dx, const double* x, const double* xref) {
	dx[0] = x[0] - xref[0];
	dx[1] = x[1] - xref[1];
	dx[2] = x[2] - xref[2];
	dx[3] = x[3] - xref[3];
	dx[4] = x[4] - xref[4];
	dx[5] = x[5] - xref[5];
}


/* It computes du = u - uref */
static void diffU(double* du, const double* u, const double* uref) {
	du[0] = u[0] - uref[0];
	du[1] = u[1] - uref[1];
}


/* it copies the content of a variable */
static void copy_nx(double* res, const double* x) {
	res[0] = x[0];
	res[1] = x[1];
	res[2] = x[2];
	res[3] = x[3];
	res[4] = x[4];
	res[5] = x[5];
}


/* it copies the content of a variable */
static void copy_Nnx(double* res, const double* x) {
	res[0] = x[0];
	res[1] = x[1];
	res[2] = x[2];
	res[3] = x[3];
	res[4] = x[4];
	res[5] = x[5];
	res[6] = x[6];
	res[7] = x[7];
	res[8] = x[8];
	res[9] = x[9];
	res[10] = x[10];
	res[11] = x[11];
	res[12] = x[12];
	res[13] = x[13];
	res[14] = x[14];
	res[15] = x[15];
	res[16] = x[16];
	res[17] = x[17];
	res[18] = x[18];
	res[19] = x[19];
	res[20] = x[20];
	res[21] = x[21];
	res[22] = x[22];
	res[23] = x[23];
	res[24] = x[24];
	res[25] = x[25];
	res[26] = x[26];
	res[27] = x[27];
	res[28] = x[28];
	res[29] = x[29];
	res[30] = x[30];
	res[31] = x[31];
	res[32] = x[32];
	res[33] = x[33];
	res[34] = x[34];
	res[35] = x[35];
	res[36] = x[36];
	res[37] = x[37];
	res[38] = x[38];
	res[39] = x[39];
	res[40] = x[40];
	res[41] = x[41];
	res[42] = x[42];
	res[43] = x[43];
	res[44] = x[44];
	res[45] = x[45];
	res[46] = x[46];
	res[47] = x[47];
	res[48] = x[48];
	res[49] = x[49];
	res[50] = x[50];
	res[51] = x[51];
	res[52] = x[52];
	res[53] = x[53];
	res[54] = x[54];
	res[55] = x[55];
	res[56] = x[56];
	res[57] = x[57];
	res[58] = x[58];
	res[59] = x[59];
	res[60] = x[60];
	res[61] = x[61];
	res[62] = x[62];
	res[63] = x[63];
	res[64] = x[64];
	res[65] = x[65];
	res[66] = x[66];
	res[67] = x[67];
	res[68] = x[68];
	res[69] = x[69];
	res[70] = x[70];
	res[71] = x[71];
	res[72] = x[72];
	res[73] = x[73];
	res[74] = x[74];
	res[75] = x[75];
	res[76] = x[76];
	res[77] = x[77];
	res[78] = x[78];
	res[79] = x[79];
	res[80] = x[80];
	res[81] = x[81];
	res[82] = x[82];
	res[83] = x[83];
	res[84] = x[84];
	res[85] = x[85];
	res[86] = x[86];
	res[87] = x[87];
	res[88] = x[88];
	res[89] = x[89];
}


/* it copies the content of a variable */
static void copy_Nnc(double* res, const double* x) {
	res[0] = x[0];
	res[1] = x[1];
	res[2] = x[2];
	res[3] = x[3];
	res[4] = x[4];
	res[5] = x[5];
	res[6] = x[6];
	res[7] = x[7];
	res[8] = x[8];
	res[9] = x[9];
	res[10] = x[10];
	res[11] = x[11];
	res[12] = x[12];
	res[13] = x[13];
	res[14] = x[14];
	res[15] = x[15];
	res[16] = x[16];
	res[17] = x[17];
	res[18] = x[18];
	res[19] = x[19];
	res[20] = x[20];
	res[21] = x[21];
	res[22] = x[22];
	res[23] = x[23];
	res[24] = x[24];
	res[25] = x[25];
	res[26] = x[26];
	res[27] = x[27];
	res[28] = x[28];
	res[29] = x[29];
	res[30] = x[30];
	res[31] = x[31];
	res[32] = x[32];
	res[33] = x[33];
	res[34] = x[34];
	res[35] = x[35];
	res[36] = x[36];
	res[37] = x[37];
	res[38] = x[38];
	res[39] = x[39];
	res[40] = x[40];
	res[41] = x[41];
	res[42] = x[42];
	res[43] = x[43];
	res[44] = x[44];
	res[45] = x[45];
	res[46] = x[46];
	res[47] = x[47];
	res[48] = x[48];
	res[49] = x[49];
	res[50] = x[50];
	res[51] = x[51];
	res[52] = x[52];
	res[53] = x[53];
	res[54] = x[54];
	res[55] = x[55];
	res[56] = x[56];
	res[57] = x[57];
	res[58] = x[58];
	res[59] = x[59];
}


/* it copies the content of a variable */
static void copy_Nnu(double* res, const double* x) {
	res[0] = x[0];
	res[1] = x[1];
	res[2] = x[2];
	res[3] = x[3];
	res[4] = x[4];
	res[5] = x[5];
	res[6] = x[6];
	res[7] = x[7];
	res[8] = x[8];
	res[9] = x[9];
	res[10] = x[10];
	res[11] = x[11];
	res[12] = x[12];
	res[13] = x[13];
	res[14] = x[14];
	res[15] = x[15];
	res[16] = x[16];
	res[17] = x[17];
	res[18] = x[18];
	res[19] = x[19];
	res[20] = x[20];
	res[21] = x[21];
	res[22] = x[22];
	res[23] = x[23];
	res[24] = x[24];
	res[25] = x[25];
	res[26] = x[26];
	res[27] = x[27];
	res[28] = x[28];
	res[29] = x[29];
}

 void det_J_and_dot_J(const double* x0, const double* u, const double* x, const double* xref, const double* uref, double* J, double* dot_J){

	double Px[6], Qx[6], mem_tmp2[6], 
		Ru[2], tmp_x = 0.0, tmp_u = 0.0;
	double dx[6], du[2];
	unsigned int ii = 0;


	/* Compute tmp_x (cost associated to last stage) */
	diffX(dx, x + 84, xref + 84);
	Pmul(Px, dx);
	dot_product_nx_nx(&tmp_x,Px, dx);

	diffU(du, u + 28, uref + 28);
	Rmul(Ru, du);
	dot_product_nu_nu(&tmp_u,Ru,du);
	(*J) = 0.5*(tmp_x + tmp_u);

	/* Start computing dot_J from the bottom */
	product_and_sum_nu(dot_J + 28,Px, Ru, G);

	for (ii=14; ii-->0; ) {
		diffX(dx, x + ii*6, xref + ii*6);
		Qmul(Qx, dx);
		dot_product_nx_nx(&tmp_x,Qx,dx);
		diffU(du, u + ii*2, uref + ii*2);
		Rmul(Ru, du);
		dot_product_nu_nu(&tmp_u,Ru,du);

		/* Increment J */
		(*J) += .50*(tmp_x + tmp_u);


		/* Compute dot_J */
		product_and_sum_nx(mem_tmp2, Px, Qx, F);
		copy_nx(Px,mem_tmp2);
		product_and_sum_nu(dot_J + ii*2, Px, Ru, G);
	}
}

 void det_J(const double* x0, const double* u, const double* x, const double* xref, const double* uref, double* J){

	double Qx[6], Ru[2], tmp_x = 0.0, tmp_u = 0.0;
	double dx[6], du[2];
	unsigned int ii = 0;

	diffX(dx, x + 84, xref + 84);
	Pmul(Qx, dx);
	dot_product_nx_nx(&tmp_x,Qx, dx);
	diffU(du, u + 28, uref + 28);
	Rmul(Ru, du);
	dot_product_nu_nu(&tmp_u,Ru,du);

	(*J) = 0.5*(tmp_x + tmp_u);
	for (ii=14; ii-->0; ) {

		diffX(dx, x + ii*6, xref + ii*6);
		Qmul(Qx, dx);
		dot_product_nx_nx(&tmp_x,Qx,dx);

		diffU(du, u + ii*2, uref + ii*2);
		Rmul(Ru, du);
		dot_product_nu_nu(&tmp_u,Ru,du);

		(*J) += .5*(tmp_x + tmp_u);
	}
}

 void build_sl_slsqr(const double* amu, const double* umb, const unsigned int nc, const unsigned int na, const unsigned int nb, double* sl, double* sl_sqr, double* gps){

	unsigned int jj;

	for (jj=0;jj< na;jj++){
		sl_sqr[jj] = my_fmax(1.0, -2.0* amu[jj]);
		gps[jj]= amu[jj] + 0.5*sl_sqr[jj];
		sl[jj] = sqrt(sl_sqr[jj]);
	}
	for (jj=na;jj< nb;jj++){
		sl_sqr[jj] = my_fmax(1.0, -2.0* umb[jj-na]);
		gps[jj]= umb[jj-na] + 0.5*sl_sqr[jj];
		sl[jj] = sqrt(sl_sqr[jj]);
	}
}


/* It initialize the slack variables sl and its squares sl_sqr 
such that, if possible, gps = g + 0.5*sl_sqr = 0 */
 void initialize_slack( const double* u, double* sl, double* sl_sqr, double* gps){

	double amu[2];
	double umb[2];

	/* Unrolling the for loop: iteration 0 of 14 */
	build_amu(&u[0], 0, &amu[0]);
	build_umb(&u[0],0,&umb[0]);
	build_sl_slsqr( &amu[0],  &umb[0], 4, 2, 4, &sl[0], &sl_sqr[0], &gps[0]);

	/* Unrolling the for loop: iteration 1 of 14 */
	build_amu(&u[2], 1, &amu[0]);
	build_umb(&u[2],1,&umb[0]);
	build_sl_slsqr( &amu[0],  &umb[0], 4, 2, 4, &sl[4], &sl_sqr[4], &gps[4]);

	/* Unrolling the for loop: iteration 2 of 14 */
	build_amu(&u[4], 2, &amu[0]);
	build_umb(&u[4],2,&umb[0]);
	build_sl_slsqr( &amu[0],  &umb[0], 4, 2, 4, &sl[8], &sl_sqr[8], &gps[8]);

	/* Unrolling the for loop: iteration 3 of 14 */
	build_amu(&u[6], 3, &amu[0]);
	build_umb(&u[6],3,&umb[0]);
	build_sl_slsqr( &amu[0],  &umb[0], 4, 2, 4, &sl[12], &sl_sqr[12], &gps[12]);

	/* Unrolling the for loop: iteration 4 of 14 */
	build_amu(&u[8], 4, &amu[0]);
	build_umb(&u[8],4,&umb[0]);
	build_sl_slsqr( &amu[0],  &umb[0], 4, 2, 4, &sl[16], &sl_sqr[16], &gps[16]);

	/* Unrolling the for loop: iteration 5 of 14 */
	build_amu(&u[10], 5, &amu[0]);
	build_umb(&u[10],5,&umb[0]);
	build_sl_slsqr( &amu[0],  &umb[0], 4, 2, 4, &sl[20], &sl_sqr[20], &gps[20]);

	/* Unrolling the for loop: iteration 6 of 14 */
	build_amu(&u[12], 6, &amu[0]);
	build_umb(&u[12],6,&umb[0]);
	build_sl_slsqr( &amu[0],  &umb[0], 4, 2, 4, &sl[24], &sl_sqr[24], &gps[24]);

	/* Unrolling the for loop: iteration 7 of 14 */
	build_amu(&u[14], 7, &amu[0]);
	build_umb(&u[14],7,&umb[0]);
	build_sl_slsqr( &amu[0],  &umb[0], 4, 2, 4, &sl[28], &sl_sqr[28], &gps[28]);

	/* Unrolling the for loop: iteration 8 of 14 */
	build_amu(&u[16], 8, &amu[0]);
	build_umb(&u[16],8,&umb[0]);
	build_sl_slsqr( &amu[0],  &umb[0], 4, 2, 4, &sl[32], &sl_sqr[32], &gps[32]);

	/* Unrolling the for loop: iteration 9 of 14 */
	build_amu(&u[18], 9, &amu[0]);
	build_umb(&u[18],9,&umb[0]);
	build_sl_slsqr( &amu[0],  &umb[0], 4, 2, 4, &sl[36], &sl_sqr[36], &gps[36]);

	/* Unrolling the for loop: iteration 10 of 14 */
	build_amu(&u[20], 10, &amu[0]);
	build_umb(&u[20],10,&umb[0]);
	build_sl_slsqr( &amu[0],  &umb[0], 4, 2, 4, &sl[40], &sl_sqr[40], &gps[40]);

	/* Unrolling the for loop: iteration 11 of 14 */
	build_amu(&u[22], 11, &amu[0]);
	build_umb(&u[22],11,&umb[0]);
	build_sl_slsqr( &amu[0],  &umb[0], 4, 2, 4, &sl[44], &sl_sqr[44], &gps[44]);

	/* Unrolling the for loop: iteration 12 of 14 */
	build_amu(&u[24], 12, &amu[0]);
	build_umb(&u[24],12,&umb[0]);
	build_sl_slsqr( &amu[0],  &umb[0], 4, 2, 4, &sl[48], &sl_sqr[48], &gps[48]);

	/* Unrolling the for loop: iteration 13 of 14 */
	build_amu(&u[26], 13, &amu[0]);
	build_umb(&u[26],13,&umb[0]);
	build_sl_slsqr( &amu[0],  &umb[0], 4, 2, 4, &sl[52], &sl_sqr[52], &gps[52]);

	/* Unrolling the for loop: iteration 14 of 14 */
	build_amu(&u[28], 14, &amu[0]);
	build_umb(&u[28],14,&umb[0]);
	build_sl_slsqr( &amu[0],  &umb[0], 4, 2, 4, &sl[56], &sl_sqr[56], &gps[56]);
}

 void build_gpsl_lowLevel(const double* amu, const double* umb, const unsigned int nc, const unsigned int na, const unsigned int nb, const double* sl_sqr, double* gps){

	unsigned int jj;

	for (jj=0;jj< na;jj++){
		gps[jj]= amu[jj] + 0.5*sl_sqr[jj];
	}
	for (jj=na;jj< nb;jj++){
		gps[jj]= umb[jj-na] + 0.5*sl_sqr[jj];
	}
}


/* It computes gps = g + 0.5 * sl_sqr */
 void build_gpsl(const double* u, const double* sl_sqr, double* gps){

	double amu[2];
	double umb[2];

	/* Unrolling the for loop: iteration 0 of 14 */
	build_amu(&u[0],0,&amu[0]);
	build_umb(&u[0],0,&umb[0]);
	build_gpsl_lowLevel( &amu[0],  &umb[0], 4, 2, 4, &sl_sqr[0], &gps[0]);

	/* Unrolling the for loop: iteration 1 of 14 */
	build_amu(&u[2],1,&amu[0]);
	build_umb(&u[2],1,&umb[0]);
	build_gpsl_lowLevel( &amu[0],  &umb[0], 4, 2, 4, &sl_sqr[4], &gps[4]);

	/* Unrolling the for loop: iteration 2 of 14 */
	build_amu(&u[4],2,&amu[0]);
	build_umb(&u[4],2,&umb[0]);
	build_gpsl_lowLevel( &amu[0],  &umb[0], 4, 2, 4, &sl_sqr[8], &gps[8]);

	/* Unrolling the for loop: iteration 3 of 14 */
	build_amu(&u[6],3,&amu[0]);
	build_umb(&u[6],3,&umb[0]);
	build_gpsl_lowLevel( &amu[0],  &umb[0], 4, 2, 4, &sl_sqr[12], &gps[12]);

	/* Unrolling the for loop: iteration 4 of 14 */
	build_amu(&u[8],4,&amu[0]);
	build_umb(&u[8],4,&umb[0]);
	build_gpsl_lowLevel( &amu[0],  &umb[0], 4, 2, 4, &sl_sqr[16], &gps[16]);

	/* Unrolling the for loop: iteration 5 of 14 */
	build_amu(&u[10],5,&amu[0]);
	build_umb(&u[10],5,&umb[0]);
	build_gpsl_lowLevel( &amu[0],  &umb[0], 4, 2, 4, &sl_sqr[20], &gps[20]);

	/* Unrolling the for loop: iteration 6 of 14 */
	build_amu(&u[12],6,&amu[0]);
	build_umb(&u[12],6,&umb[0]);
	build_gpsl_lowLevel( &amu[0],  &umb[0], 4, 2, 4, &sl_sqr[24], &gps[24]);

	/* Unrolling the for loop: iteration 7 of 14 */
	build_amu(&u[14],7,&amu[0]);
	build_umb(&u[14],7,&umb[0]);
	build_gpsl_lowLevel( &amu[0],  &umb[0], 4, 2, 4, &sl_sqr[28], &gps[28]);

	/* Unrolling the for loop: iteration 8 of 14 */
	build_amu(&u[16],8,&amu[0]);
	build_umb(&u[16],8,&umb[0]);
	build_gpsl_lowLevel( &amu[0],  &umb[0], 4, 2, 4, &sl_sqr[32], &gps[32]);

	/* Unrolling the for loop: iteration 9 of 14 */
	build_amu(&u[18],9,&amu[0]);
	build_umb(&u[18],9,&umb[0]);
	build_gpsl_lowLevel( &amu[0],  &umb[0], 4, 2, 4, &sl_sqr[36], &gps[36]);

	/* Unrolling the for loop: iteration 10 of 14 */
	build_amu(&u[20],10,&amu[0]);
	build_umb(&u[20],10,&umb[0]);
	build_gpsl_lowLevel( &amu[0],  &umb[0], 4, 2, 4, &sl_sqr[40], &gps[40]);

	/* Unrolling the for loop: iteration 11 of 14 */
	build_amu(&u[22],11,&amu[0]);
	build_umb(&u[22],11,&umb[0]);
	build_gpsl_lowLevel( &amu[0],  &umb[0], 4, 2, 4, &sl_sqr[44], &gps[44]);

	/* Unrolling the for loop: iteration 12 of 14 */
	build_amu(&u[24],12,&amu[0]);
	build_umb(&u[24],12,&umb[0]);
	build_gpsl_lowLevel( &amu[0],  &umb[0], 4, 2, 4, &sl_sqr[48], &gps[48]);

	/* Unrolling the for loop: iteration 13 of 14 */
	build_amu(&u[26],13,&amu[0]);
	build_umb(&u[26],13,&umb[0]);
	build_gpsl_lowLevel( &amu[0],  &umb[0], 4, 2, 4, &sl_sqr[52], &gps[52]);

	/* Unrolling the for loop: iteration 14 of 14 */
	build_amu(&u[28],14,&amu[0]);
	build_umb(&u[28],14,&umb[0]);
	build_gpsl_lowLevel( &amu[0],  &umb[0], 4, 2, 4, &sl_sqr[56], &gps[56]);

}

/* It computes x_sqr = x.*x */
 void build_sqr_Nnc( const double *x, double *x_sqr){

	unsigned int ii = 0;

	for (ii=0;ii<60;ii++)
		x_sqr[ii] = x[ii]*x[ii];
}


/* dot_product of dim N*nu */
static void dot_product_Nnu(double* r, const double* du, const double* u) {
	r[0] = u[0]*du[0] + u[1]*du[1] + u[2]*du[2] + u[3]*du[3] + u[4]*du[4] + u[5]*du[5] + u[6]*du[6] + u[7]*du[7] + u[8]*du[8] + u[9]*du[9] + u[10]*du[10] + u[11]*du[11] + u[12]*du[12] + u[13]*du[13] + u[14]*du[14] + u[15]*du[15] + u[16]*du[16] + u[17]*du[17] + u[18]*du[18] + u[19]*du[19] + u[20]*du[20] + u[21]*du[21] + u[22]*du[22] + u[23]*du[23] + u[24]*du[24] + u[25]*du[25] + u[26]*du[26] + u[27]*du[27] + u[28]*du[28] + u[29]*du[29];
}

static void build_vNnc_lb_1(double* z, const double* x1, const double* x2) {
	z[0] = alpha_inv_lb_1[0]*x1[0] + x2[0];
	z[1] = alpha_inv_lb_1[0]*x1[1] + x2[1];
}

 void build_vNnc_lb(const double* gps, const double* dot_J, const unsigned int k, double* res){
	if(k <= 14) {
		build_vNnc_lb_1(&res[0], &gps[0], &dot_J[0]);
	}
}
static void build_vNnc_ub_1(double* z, const double* x1, const double* x2) {
	z[0] = alpha_inv_ub_1[0]*x1[0] - x2[0];
	z[1] = alpha_inv_ub_1[0]*x1[1] - x2[1];
}

 void build_vNnc_ub(const double* gps, const double* dot_J, const unsigned int k, double* res){
	if(k <= 14) {
		build_vNnc_ub_1(&res[0], &gps[0], &dot_J[0]);
	}
}
static void minus_Ina_muG_1(double* z, const double* x1) {
	z[0] = -x1[0];
	z[1] = -x1[1];
}

 void minus_Ina_muG(const double* muG, const unsigned int k, double* res){
	if(k <= 14) {
		minus_Ina_muG_1(&res[0], &muG[0]);
	}
}
static void Inb_muG_1(double* z, const double* x1) {
	z[0] = x1[0];
	z[1] = x1[1];
}

 void Inb_muG(const double* muG, const unsigned int k, double* res){
	if(k <= 14) {
		Inb_muG_1(&res[0], &muG[0]);
	}
}
static void sum_nr_constr(double* r, const double* v1, const double* v2, const double* v3) {
	r[0] = v1[0] + v2[0] + v3[0];
	r[1] = v1[1] + v2[1] + v3[1];
}

static void scale_nu(double* z, const double* x) {
	z[0] = alpha_I[0]*x[0];
	z[1] = alpha_I[0]*x[1];
}

static void scale_nc_1(double* z, const double* x) {
	z[0] = alpha_I[0]*x[0];
	z[1] = alpha_I[0]*x[1];
	z[2] = alpha_I[0]*x[2];
	z[3] = alpha_I[0]*x[3];
}

 void product_matlab_nc_1(const double* x, const double* y, double* z){

	unsigned int ii=0;

	for (ii=0;ii<4;ii++)
		z[ii] = x[ii]*y[ii];

}
 void minus_scale_nc(const double* x, const unsigned int k, double* res){
	if(k <= 14) {
		scale_nc_1( &res[0], &x[0]);
	}
}
 void product_matlab_nc(const double* x, const double* y, const unsigned int k, double* res){
	if(k <= 14) {
		product_matlab_nc_1(&x[0], &y[0], &res[0]);
	}
}
static void solveConstraintSystem_DMult1(double* r, const double* v1, const double* v2, const double* M1, const double* M2) {
	r[0] = M1[0]*v1[0] + M2[0]*v2[0];
	r[1] = M1[1]*v1[1] + M2[1]*v2[1];
	r[2] = M1[2]*v1[0] + M2[2]*v2[0];
	r[3] = M1[3]*v1[1] + M2[3]*v2[1];
}

/** 
 * @brief Processe slacks and generate matrices D1 and D2.
 * @param s The (squared) slacks for stage k. A vector of dimension 4.
 * @param D1 Return value. A matrix with 4 elements.
 * @param D2 Return value. A matrix with 4 elements.
 * @param tmp A vector of dimension 2 used for storing intermediate results.
 */
static  void solveConstraintSystem_processSlacks(const double* s, double* D1, double* D2, double* tmp) {
	/** Compute intermediate results involving divisions **/
	tmp[0] = 1.0/(s[0]*s[2+0] + s[0] + s[2+0]);
	tmp[1] = 1.0/(s[1]*s[2+1] + s[1] + s[2+1]);
	/** Construct M_k := [dp_k]'*diag({S_i}_{i=1}^nu)*[dp_k]+diag({s_{k,i+L+U}}_{i=1}^np)
	     where S_i = 1 if i-th component of u has neither lower nor upper bounds
	               = s_{k,j}/(1+s_{k,j}) if i-th component of u has only lower bounds (j is index of lower bound)
	               = s_{k,L+j}/(1+s_{k,L+j}) if i-th component of u has only upper bounds (j is index of upper bound)
	               = s_{k,j}*s_{k,L+l}/(s_{k,j}*s_{k,L+l} + s_{k,j} + s_{k,L+l})if i-th component of u has lower and upper bounds (j is index of lower bound, l is index of upper bound)	     and nu is the number of inputs, np is the number of constraints, L is the number of lower bounds and U is the number of upper bounds. **/
	/* Nothing to do */
	/* Nothing to do */
	/* Nothing to do */
	/* Nothing to do */
	/** Construct D1 and D2, where **/
	/* Construct D1 */
	D1[0] = tmp[0]*(1.0 + s[2+0]); /* Element #1 of D1 (1,1), corresponds to lower bound #1 and upper bound #1, component (1) */
	D1[1] = tmp[1]*(1.0 + s[2+1]); /* Element #2 of D1 (2,2), corresponds to lower bound #2 and upper bound #2, component (2) */
	D1[2] = tmp[0]; /* Element #3 of D1 (3,1), corresponds to lower bound #1 and upper bound #1 */
	D1[3] = tmp[1]; /* Element #4 of D1 (4,2), corresponds to lower bound #2 and upper bound #2 */
	/* Construct D2 */
	D2[0] = tmp[0]; /* Element #1 of D2 (1,1), corresponds to lower bound #1 and upper bound #1 */
	D2[1] = tmp[1]; /* Element #2 of D2 (2,2), corresponds to lower bound #2 and upper bound #2 */
	D2[2] = tmp[0]*(1.0 + s[0]); /* Element #3 of D2 (3,1), corresponds to upper bound #1 and lower bound #1, component (1) */
	D2[3] = tmp[1]*(1.0 + s[1]); /* Element #4 of D2 (4,2), corresponds to upper bound #2 and lower bound #2, component (2) */
}


/** 
 */
static void solveConstraintSystem(double* r, const double* v, const double* sq) {
	unsigned int i;
	/* Variables to store intermediate results */
	double D1[4]; double D2[4];
	double tmp[2];

	/** k = 0 **/
	/* Build auxiliary matrices */
	solveConstraintSystem_processSlacks(&sq[0], D1, D2, tmp);
	/* Initialize vector r_k^{a,b} */
	solveConstraintSystem_DMult1(&r[0], &v[0], &v[0+2], D1, D2);

	/** k = 1 **/
	/* Build auxiliary matrices */
	solveConstraintSystem_processSlacks(&sq[4], D1, D2, tmp);
	/* Initialize vector r_k^{a,b} */
	solveConstraintSystem_DMult1(&r[4], &v[4], &v[4+2], D1, D2);

	/** k = 2 **/
	/* Build auxiliary matrices */
	solveConstraintSystem_processSlacks(&sq[8], D1, D2, tmp);
	/* Initialize vector r_k^{a,b} */
	solveConstraintSystem_DMult1(&r[8], &v[8], &v[8+2], D1, D2);

	/** k = 3 **/
	/* Build auxiliary matrices */
	solveConstraintSystem_processSlacks(&sq[12], D1, D2, tmp);
	/* Initialize vector r_k^{a,b} */
	solveConstraintSystem_DMult1(&r[12], &v[12], &v[12+2], D1, D2);

	/** k = 4 **/
	/* Build auxiliary matrices */
	solveConstraintSystem_processSlacks(&sq[16], D1, D2, tmp);
	/* Initialize vector r_k^{a,b} */
	solveConstraintSystem_DMult1(&r[16], &v[16], &v[16+2], D1, D2);

	/** k = 5 **/
	/* Build auxiliary matrices */
	solveConstraintSystem_processSlacks(&sq[20], D1, D2, tmp);
	/* Initialize vector r_k^{a,b} */
	solveConstraintSystem_DMult1(&r[20], &v[20], &v[20+2], D1, D2);

	/** k = 6 **/
	/* Build auxiliary matrices */
	solveConstraintSystem_processSlacks(&sq[24], D1, D2, tmp);
	/* Initialize vector r_k^{a,b} */
	solveConstraintSystem_DMult1(&r[24], &v[24], &v[24+2], D1, D2);

	/** k = 7 **/
	/* Build auxiliary matrices */
	solveConstraintSystem_processSlacks(&sq[28], D1, D2, tmp);
	/* Initialize vector r_k^{a,b} */
	solveConstraintSystem_DMult1(&r[28], &v[28], &v[28+2], D1, D2);

	/** k = 8 **/
	/* Build auxiliary matrices */
	solveConstraintSystem_processSlacks(&sq[32], D1, D2, tmp);
	/* Initialize vector r_k^{a,b} */
	solveConstraintSystem_DMult1(&r[32], &v[32], &v[32+2], D1, D2);

	/** k = 9 **/
	/* Build auxiliary matrices */
	solveConstraintSystem_processSlacks(&sq[36], D1, D2, tmp);
	/* Initialize vector r_k^{a,b} */
	solveConstraintSystem_DMult1(&r[36], &v[36], &v[36+2], D1, D2);

	/** k = 10 **/
	/* Build auxiliary matrices */
	solveConstraintSystem_processSlacks(&sq[40], D1, D2, tmp);
	/* Initialize vector r_k^{a,b} */
	solveConstraintSystem_DMult1(&r[40], &v[40], &v[40+2], D1, D2);

	/** k = 11 **/
	/* Build auxiliary matrices */
	solveConstraintSystem_processSlacks(&sq[44], D1, D2, tmp);
	/* Initialize vector r_k^{a,b} */
	solveConstraintSystem_DMult1(&r[44], &v[44], &v[44+2], D1, D2);

	/** k = 12 **/
	/* Build auxiliary matrices */
	solveConstraintSystem_processSlacks(&sq[48], D1, D2, tmp);
	/* Initialize vector r_k^{a,b} */
	solveConstraintSystem_DMult1(&r[48], &v[48], &v[48+2], D1, D2);

	/** k = 13 **/
	/* Build auxiliary matrices */
	solveConstraintSystem_processSlacks(&sq[52], D1, D2, tmp);
	/* Initialize vector r_k^{a,b} */
	solveConstraintSystem_DMult1(&r[52], &v[52], &v[52+2], D1, D2);

	/** k = 14 **/
	/* Build auxiliary matrices */
	solveConstraintSystem_processSlacks(&sq[56], D1, D2, tmp);
	/* Initialize vector r_k^{a,b} */
	solveConstraintSystem_DMult1(&r[56], &v[56], &v[56+2], D1, D2);
}



/* Gradient step function */
 void gradient_step(const double* dot_J, const double* u, const double* sl,
	const double* sl_sqr, const double* gps, double* du, double* dsl, double* muG){

	double v_Nnc[60], tmp_nu[2], tmp_nc_m[4], tmp_contr = 0.0;
	double temp_lb[2];
	double temp_ub[2];

	/* loop unrolling: step 0 of 14 */
	build_vNnc_lb(&gps[0], &dot_J[0], 0, &v_Nnc[0]);
	build_vNnc_ub(&gps[2], &dot_J[0], 0, &v_Nnc[2]);

	/* loop unrolling: step 1 of 14 */
	build_vNnc_lb(&gps[4], &dot_J[2], 1, &v_Nnc[4]);
	build_vNnc_ub(&gps[6], &dot_J[2], 1, &v_Nnc[6]);

	/* loop unrolling: step 2 of 14 */
	build_vNnc_lb(&gps[8], &dot_J[4], 2, &v_Nnc[8]);
	build_vNnc_ub(&gps[10], &dot_J[4], 2, &v_Nnc[10]);

	/* loop unrolling: step 3 of 14 */
	build_vNnc_lb(&gps[12], &dot_J[6], 3, &v_Nnc[12]);
	build_vNnc_ub(&gps[14], &dot_J[6], 3, &v_Nnc[14]);

	/* loop unrolling: step 4 of 14 */
	build_vNnc_lb(&gps[16], &dot_J[8], 4, &v_Nnc[16]);
	build_vNnc_ub(&gps[18], &dot_J[8], 4, &v_Nnc[18]);

	/* loop unrolling: step 5 of 14 */
	build_vNnc_lb(&gps[20], &dot_J[10], 5, &v_Nnc[20]);
	build_vNnc_ub(&gps[22], &dot_J[10], 5, &v_Nnc[22]);

	/* loop unrolling: step 6 of 14 */
	build_vNnc_lb(&gps[24], &dot_J[12], 6, &v_Nnc[24]);
	build_vNnc_ub(&gps[26], &dot_J[12], 6, &v_Nnc[26]);

	/* loop unrolling: step 7 of 14 */
	build_vNnc_lb(&gps[28], &dot_J[14], 7, &v_Nnc[28]);
	build_vNnc_ub(&gps[30], &dot_J[14], 7, &v_Nnc[30]);

	/* loop unrolling: step 8 of 14 */
	build_vNnc_lb(&gps[32], &dot_J[16], 8, &v_Nnc[32]);
	build_vNnc_ub(&gps[34], &dot_J[16], 8, &v_Nnc[34]);

	/* loop unrolling: step 9 of 14 */
	build_vNnc_lb(&gps[36], &dot_J[18], 9, &v_Nnc[36]);
	build_vNnc_ub(&gps[38], &dot_J[18], 9, &v_Nnc[38]);

	/* loop unrolling: step 10 of 14 */
	build_vNnc_lb(&gps[40], &dot_J[20], 10, &v_Nnc[40]);
	build_vNnc_ub(&gps[42], &dot_J[20], 10, &v_Nnc[42]);

	/* loop unrolling: step 11 of 14 */
	build_vNnc_lb(&gps[44], &dot_J[22], 11, &v_Nnc[44]);
	build_vNnc_ub(&gps[46], &dot_J[22], 11, &v_Nnc[46]);

	/* loop unrolling: step 12 of 14 */
	build_vNnc_lb(&gps[48], &dot_J[24], 12, &v_Nnc[48]);
	build_vNnc_ub(&gps[50], &dot_J[24], 12, &v_Nnc[50]);

	/* loop unrolling: step 13 of 14 */
	build_vNnc_lb(&gps[52], &dot_J[26], 13, &v_Nnc[52]);
	build_vNnc_ub(&gps[54], &dot_J[26], 13, &v_Nnc[54]);

	/* loop unrolling: step 14 of 14 */
	build_vNnc_lb(&gps[56], &dot_J[28], 14, &v_Nnc[56]);
	build_vNnc_ub(&gps[58], &dot_J[28], 14, &v_Nnc[58]);

	solveConstraintSystem(&muG[0], &v_Nnc[0], &sl_sqr[0]);

	/* loop unrolling: step 0 of 14 */
	minus_Ina_muG(&muG[0], 0, &temp_lb[0]);
	Inb_muG(&muG[2], 0, &temp_ub[0]);
	sum_nr_constr(&tmp_nu[0], &temp_lb[0], &temp_ub[0], &dot_J[0]);
	scale_nu(&du[0], &tmp_nu[0]);
	product_matlab_nc(&sl[0], &muG[0], 0, &tmp_nc_m[0]);
	minus_scale_nc(&tmp_nc_m[0], 0, &dsl[0]);

	/* loop unrolling: step 1 of 14 */
	minus_Ina_muG(&muG[4], 1, &temp_lb[0]);
	Inb_muG(&muG[6], 1, &temp_ub[0]);
	sum_nr_constr(&tmp_nu[0], &temp_lb[0], &temp_ub[0], &dot_J[2]);
	scale_nu(&du[2], &tmp_nu[0]);
	product_matlab_nc(&sl[4], &muG[4], 1, &tmp_nc_m[0]);
	minus_scale_nc(&tmp_nc_m[0], 1, &dsl[4]);

	/* loop unrolling: step 2 of 14 */
	minus_Ina_muG(&muG[8], 2, &temp_lb[0]);
	Inb_muG(&muG[10], 2, &temp_ub[0]);
	sum_nr_constr(&tmp_nu[0], &temp_lb[0], &temp_ub[0], &dot_J[4]);
	scale_nu(&du[4], &tmp_nu[0]);
	product_matlab_nc(&sl[8], &muG[8], 2, &tmp_nc_m[0]);
	minus_scale_nc(&tmp_nc_m[0], 2, &dsl[8]);

	/* loop unrolling: step 3 of 14 */
	minus_Ina_muG(&muG[12], 3, &temp_lb[0]);
	Inb_muG(&muG[14], 3, &temp_ub[0]);
	sum_nr_constr(&tmp_nu[0], &temp_lb[0], &temp_ub[0], &dot_J[6]);
	scale_nu(&du[6], &tmp_nu[0]);
	product_matlab_nc(&sl[12], &muG[12], 3, &tmp_nc_m[0]);
	minus_scale_nc(&tmp_nc_m[0], 3, &dsl[12]);

	/* loop unrolling: step 4 of 14 */
	minus_Ina_muG(&muG[16], 4, &temp_lb[0]);
	Inb_muG(&muG[18], 4, &temp_ub[0]);
	sum_nr_constr(&tmp_nu[0], &temp_lb[0], &temp_ub[0], &dot_J[8]);
	scale_nu(&du[8], &tmp_nu[0]);
	product_matlab_nc(&sl[16], &muG[16], 4, &tmp_nc_m[0]);
	minus_scale_nc(&tmp_nc_m[0], 4, &dsl[16]);

	/* loop unrolling: step 5 of 14 */
	minus_Ina_muG(&muG[20], 5, &temp_lb[0]);
	Inb_muG(&muG[22], 5, &temp_ub[0]);
	sum_nr_constr(&tmp_nu[0], &temp_lb[0], &temp_ub[0], &dot_J[10]);
	scale_nu(&du[10], &tmp_nu[0]);
	product_matlab_nc(&sl[20], &muG[20], 5, &tmp_nc_m[0]);
	minus_scale_nc(&tmp_nc_m[0], 5, &dsl[20]);

	/* loop unrolling: step 6 of 14 */
	minus_Ina_muG(&muG[24], 6, &temp_lb[0]);
	Inb_muG(&muG[26], 6, &temp_ub[0]);
	sum_nr_constr(&tmp_nu[0], &temp_lb[0], &temp_ub[0], &dot_J[12]);
	scale_nu(&du[12], &tmp_nu[0]);
	product_matlab_nc(&sl[24], &muG[24], 6, &tmp_nc_m[0]);
	minus_scale_nc(&tmp_nc_m[0], 6, &dsl[24]);

	/* loop unrolling: step 7 of 14 */
	minus_Ina_muG(&muG[28], 7, &temp_lb[0]);
	Inb_muG(&muG[30], 7, &temp_ub[0]);
	sum_nr_constr(&tmp_nu[0], &temp_lb[0], &temp_ub[0], &dot_J[14]);
	scale_nu(&du[14], &tmp_nu[0]);
	product_matlab_nc(&sl[28], &muG[28], 7, &tmp_nc_m[0]);
	minus_scale_nc(&tmp_nc_m[0], 7, &dsl[28]);

	/* loop unrolling: step 8 of 14 */
	minus_Ina_muG(&muG[32], 8, &temp_lb[0]);
	Inb_muG(&muG[34], 8, &temp_ub[0]);
	sum_nr_constr(&tmp_nu[0], &temp_lb[0], &temp_ub[0], &dot_J[16]);
	scale_nu(&du[16], &tmp_nu[0]);
	product_matlab_nc(&sl[32], &muG[32], 8, &tmp_nc_m[0]);
	minus_scale_nc(&tmp_nc_m[0], 8, &dsl[32]);

	/* loop unrolling: step 9 of 14 */
	minus_Ina_muG(&muG[36], 9, &temp_lb[0]);
	Inb_muG(&muG[38], 9, &temp_ub[0]);
	sum_nr_constr(&tmp_nu[0], &temp_lb[0], &temp_ub[0], &dot_J[18]);
	scale_nu(&du[18], &tmp_nu[0]);
	product_matlab_nc(&sl[36], &muG[36], 9, &tmp_nc_m[0]);
	minus_scale_nc(&tmp_nc_m[0], 9, &dsl[36]);

	/* loop unrolling: step 10 of 14 */
	minus_Ina_muG(&muG[40], 10, &temp_lb[0]);
	Inb_muG(&muG[42], 10, &temp_ub[0]);
	sum_nr_constr(&tmp_nu[0], &temp_lb[0], &temp_ub[0], &dot_J[20]);
	scale_nu(&du[20], &tmp_nu[0]);
	product_matlab_nc(&sl[40], &muG[40], 10, &tmp_nc_m[0]);
	minus_scale_nc(&tmp_nc_m[0], 10, &dsl[40]);

	/* loop unrolling: step 11 of 14 */
	minus_Ina_muG(&muG[44], 11, &temp_lb[0]);
	Inb_muG(&muG[46], 11, &temp_ub[0]);
	sum_nr_constr(&tmp_nu[0], &temp_lb[0], &temp_ub[0], &dot_J[22]);
	scale_nu(&du[22], &tmp_nu[0]);
	product_matlab_nc(&sl[44], &muG[44], 11, &tmp_nc_m[0]);
	minus_scale_nc(&tmp_nc_m[0], 11, &dsl[44]);

	/* loop unrolling: step 12 of 14 */
	minus_Ina_muG(&muG[48], 12, &temp_lb[0]);
	Inb_muG(&muG[50], 12, &temp_ub[0]);
	sum_nr_constr(&tmp_nu[0], &temp_lb[0], &temp_ub[0], &dot_J[24]);
	scale_nu(&du[24], &tmp_nu[0]);
	product_matlab_nc(&sl[48], &muG[48], 12, &tmp_nc_m[0]);
	minus_scale_nc(&tmp_nc_m[0], 12, &dsl[48]);

	/* loop unrolling: step 13 of 14 */
	minus_Ina_muG(&muG[52], 13, &temp_lb[0]);
	Inb_muG(&muG[54], 13, &temp_ub[0]);
	sum_nr_constr(&tmp_nu[0], &temp_lb[0], &temp_ub[0], &dot_J[26]);
	scale_nu(&du[26], &tmp_nu[0]);
	product_matlab_nc(&sl[52], &muG[52], 13, &tmp_nc_m[0]);
	minus_scale_nc(&tmp_nc_m[0], 13, &dsl[52]);

	/* loop unrolling: step 14 of 14 */
	minus_Ina_muG(&muG[56], 14, &temp_lb[0]);
	Inb_muG(&muG[58], 14, &temp_ub[0]);
	sum_nr_constr(&tmp_nu[0], &temp_lb[0], &temp_ub[0], &dot_J[28]);
	scale_nu(&du[28], &tmp_nu[0]);
	product_matlab_nc(&sl[56], &muG[56], 14, &tmp_nc_m[0]);
	minus_scale_nc(&tmp_nc_m[0], 14, &dsl[56]);
}


/* dot_product of dim N*nc */
static void dot_product_Nnc(double* r, const double* du, const double* u) {
	r[0] = u[0]*du[0] + u[1]*du[1] + u[2]*du[2] + u[3]*du[3] + u[4]*du[4] + u[5]*du[5] + u[6]*du[6] + u[7]*du[7] + u[8]*du[8] + u[9]*du[9] + u[10]*du[10] + u[11]*du[11] + u[12]*du[12] + u[13]*du[13] + u[14]*du[14] + u[15]*du[15] + u[16]*du[16] + u[17]*du[17] + u[18]*du[18] + u[19]*du[19] + u[20]*du[20] + u[21]*du[21] + u[22]*du[22] + u[23]*du[23] + u[24]*du[24] + u[25]*du[25] + u[26]*du[26] + u[27]*du[27] + u[28]*du[28] + u[29]*du[29] + u[30]*du[30] + u[31]*du[31] + u[32]*du[32] + u[33]*du[33] + u[34]*du[34] + u[35]*du[35] + u[36]*du[36] + u[37]*du[37] + u[38]*du[38] + u[39]*du[39] + u[40]*du[40] + u[41]*du[41] + u[42]*du[42] + u[43]*du[43] + u[44]*du[44] + u[45]*du[45] + u[46]*du[46] + u[47]*du[47] + u[48]*du[48] + u[49]*du[49] + u[50]*du[50] + u[51]*du[51] + u[52]*du[52] + u[53]*du[53] + u[54]*du[54] + u[55]*du[55] + u[56]*du[56] + u[57]*du[57] + u[58]*du[58] + u[59]*du[59];
}


/* It computes dx = x - xref */
static void diff_Nnc(double* dx, const double* x, const double* xref) {
	dx[0] = x[0] - xref[0];
	dx[1] = x[1] - xref[1];
	dx[2] = x[2] - xref[2];
	dx[3] = x[3] - xref[3];
	dx[4] = x[4] - xref[4];
	dx[5] = x[5] - xref[5];
	dx[6] = x[6] - xref[6];
	dx[7] = x[7] - xref[7];
	dx[8] = x[8] - xref[8];
	dx[9] = x[9] - xref[9];
	dx[10] = x[10] - xref[10];
	dx[11] = x[11] - xref[11];
	dx[12] = x[12] - xref[12];
	dx[13] = x[13] - xref[13];
	dx[14] = x[14] - xref[14];
	dx[15] = x[15] - xref[15];
	dx[16] = x[16] - xref[16];
	dx[17] = x[17] - xref[17];
	dx[18] = x[18] - xref[18];
	dx[19] = x[19] - xref[19];
	dx[20] = x[20] - xref[20];
	dx[21] = x[21] - xref[21];
	dx[22] = x[22] - xref[22];
	dx[23] = x[23] - xref[23];
	dx[24] = x[24] - xref[24];
	dx[25] = x[25] - xref[25];
	dx[26] = x[26] - xref[26];
	dx[27] = x[27] - xref[27];
	dx[28] = x[28] - xref[28];
	dx[29] = x[29] - xref[29];
	dx[30] = x[30] - xref[30];
	dx[31] = x[31] - xref[31];
	dx[32] = x[32] - xref[32];
	dx[33] = x[33] - xref[33];
	dx[34] = x[34] - xref[34];
	dx[35] = x[35] - xref[35];
	dx[36] = x[36] - xref[36];
	dx[37] = x[37] - xref[37];
	dx[38] = x[38] - xref[38];
	dx[39] = x[39] - xref[39];
	dx[40] = x[40] - xref[40];
	dx[41] = x[41] - xref[41];
	dx[42] = x[42] - xref[42];
	dx[43] = x[43] - xref[43];
	dx[44] = x[44] - xref[44];
	dx[45] = x[45] - xref[45];
	dx[46] = x[46] - xref[46];
	dx[47] = x[47] - xref[47];
	dx[48] = x[48] - xref[48];
	dx[49] = x[49] - xref[49];
	dx[50] = x[50] - xref[50];
	dx[51] = x[51] - xref[51];
	dx[52] = x[52] - xref[52];
	dx[53] = x[53] - xref[53];
	dx[54] = x[54] - xref[54];
	dx[55] = x[55] - xref[55];
	dx[56] = x[56] - xref[56];
	dx[57] = x[57] - xref[57];
	dx[58] = x[58] - xref[58];
	dx[59] = x[59] - xref[59];
}


/* It computes the merit function phi */
 void det_phi (const double J, const double* gps, const double* mu, const double rho, double* phi){

	double tmp[60], pr = 0.0;
	unsigned int ii = 0;

	for (ii=60; ii--; )
		tmp[ii] = mu[ii] + 0.50*rho*gps[ii];
	dot_product_Nnc(&pr, gps,tmp);

	(*phi) = J + pr;
}


/* It computes the derivative of the merit function dot_phi */
 void det_dot_phi (const double* du, const double* DJ, const double rho, const double* gps, 
	 const double* mu, const double* dm, double* dot_phi){

	double tmp_prod[60], prod_1 = 0.0, prod_2 = 0.0;
	unsigned int ii = 0;

	for (ii=60; ii--; )
		tmp_prod[ii] = mu[ii] - dm[ii] + rho*gps[ii];
	dot_product_Nnu(&prod_1, du, DJ);

	dot_product_Nnc(&prod_2, gps,tmp_prod);

	(*dot_phi) = prod_1 - prod_2;
}


/* It checks the decrease condition of the merit function */
 int conditions_rho_PM_simpler (const double dot_phi, const double du_sqr, const double dsl_sqr, const double alpha){

	unsigned int res = 2;

	if (dot_phi <= (-0.5/alpha)*(du_sqr + dsl_sqr))
		res = 1;
	else
		res = 0;

	return res;
}


/*  dx = t*x + xref */
static void weighted_sum_Nnc(double* dx, const double* x, const double* xref, const double* I1) {
	dx[0] = I1[0]*x[0] + xref[0];
	dx[1] = I1[0]*x[1] + xref[1];
	dx[2] = I1[0]*x[2] + xref[2];
	dx[3] = I1[0]*x[3] + xref[3];
	dx[4] = I1[0]*x[4] + xref[4];
	dx[5] = I1[0]*x[5] + xref[5];
	dx[6] = I1[0]*x[6] + xref[6];
	dx[7] = I1[0]*x[7] + xref[7];
	dx[8] = I1[0]*x[8] + xref[8];
	dx[9] = I1[0]*x[9] + xref[9];
	dx[10] = I1[0]*x[10] + xref[10];
	dx[11] = I1[0]*x[11] + xref[11];
	dx[12] = I1[0]*x[12] + xref[12];
	dx[13] = I1[0]*x[13] + xref[13];
	dx[14] = I1[0]*x[14] + xref[14];
	dx[15] = I1[0]*x[15] + xref[15];
	dx[16] = I1[0]*x[16] + xref[16];
	dx[17] = I1[0]*x[17] + xref[17];
	dx[18] = I1[0]*x[18] + xref[18];
	dx[19] = I1[0]*x[19] + xref[19];
	dx[20] = I1[0]*x[20] + xref[20];
	dx[21] = I1[0]*x[21] + xref[21];
	dx[22] = I1[0]*x[22] + xref[22];
	dx[23] = I1[0]*x[23] + xref[23];
	dx[24] = I1[0]*x[24] + xref[24];
	dx[25] = I1[0]*x[25] + xref[25];
	dx[26] = I1[0]*x[26] + xref[26];
	dx[27] = I1[0]*x[27] + xref[27];
	dx[28] = I1[0]*x[28] + xref[28];
	dx[29] = I1[0]*x[29] + xref[29];
	dx[30] = I1[0]*x[30] + xref[30];
	dx[31] = I1[0]*x[31] + xref[31];
	dx[32] = I1[0]*x[32] + xref[32];
	dx[33] = I1[0]*x[33] + xref[33];
	dx[34] = I1[0]*x[34] + xref[34];
	dx[35] = I1[0]*x[35] + xref[35];
	dx[36] = I1[0]*x[36] + xref[36];
	dx[37] = I1[0]*x[37] + xref[37];
	dx[38] = I1[0]*x[38] + xref[38];
	dx[39] = I1[0]*x[39] + xref[39];
	dx[40] = I1[0]*x[40] + xref[40];
	dx[41] = I1[0]*x[41] + xref[41];
	dx[42] = I1[0]*x[42] + xref[42];
	dx[43] = I1[0]*x[43] + xref[43];
	dx[44] = I1[0]*x[44] + xref[44];
	dx[45] = I1[0]*x[45] + xref[45];
	dx[46] = I1[0]*x[46] + xref[46];
	dx[47] = I1[0]*x[47] + xref[47];
	dx[48] = I1[0]*x[48] + xref[48];
	dx[49] = I1[0]*x[49] + xref[49];
	dx[50] = I1[0]*x[50] + xref[50];
	dx[51] = I1[0]*x[51] + xref[51];
	dx[52] = I1[0]*x[52] + xref[52];
	dx[53] = I1[0]*x[53] + xref[53];
	dx[54] = I1[0]*x[54] + xref[54];
	dx[55] = I1[0]*x[55] + xref[55];
	dx[56] = I1[0]*x[56] + xref[56];
	dx[57] = I1[0]*x[57] + xref[57];
	dx[58] = I1[0]*x[58] + xref[58];
	dx[59] = I1[0]*x[59] + xref[59];
}


/*  du = t*u + uref */
static void weighted_sum_Nnu(double* dx, const double* x, const double* xref, const double* I1) {
	dx[0] = I1[0]*x[0] + xref[0];
	dx[1] = I1[0]*x[1] + xref[1];
	dx[2] = I1[0]*x[2] + xref[2];
	dx[3] = I1[0]*x[3] + xref[3];
	dx[4] = I1[0]*x[4] + xref[4];
	dx[5] = I1[0]*x[5] + xref[5];
	dx[6] = I1[0]*x[6] + xref[6];
	dx[7] = I1[0]*x[7] + xref[7];
	dx[8] = I1[0]*x[8] + xref[8];
	dx[9] = I1[0]*x[9] + xref[9];
	dx[10] = I1[0]*x[10] + xref[10];
	dx[11] = I1[0]*x[11] + xref[11];
	dx[12] = I1[0]*x[12] + xref[12];
	dx[13] = I1[0]*x[13] + xref[13];
	dx[14] = I1[0]*x[14] + xref[14];
	dx[15] = I1[0]*x[15] + xref[15];
	dx[16] = I1[0]*x[16] + xref[16];
	dx[17] = I1[0]*x[17] + xref[17];
	dx[18] = I1[0]*x[18] + xref[18];
	dx[19] = I1[0]*x[19] + xref[19];
	dx[20] = I1[0]*x[20] + xref[20];
	dx[21] = I1[0]*x[21] + xref[21];
	dx[22] = I1[0]*x[22] + xref[22];
	dx[23] = I1[0]*x[23] + xref[23];
	dx[24] = I1[0]*x[24] + xref[24];
	dx[25] = I1[0]*x[25] + xref[25];
	dx[26] = I1[0]*x[26] + xref[26];
	dx[27] = I1[0]*x[27] + xref[27];
	dx[28] = I1[0]*x[28] + xref[28];
	dx[29] = I1[0]*x[29] + xref[29];
}


/* it generates a safeguarded quadratic interpolation */
 double quadratic_interp (const double f_l, const double g_l, const double t_u, const double f_u) {

	double t_theo;

	t_theo = -0.5*(g_l*t_u*t_u)/(f_u - f_l - g_l*t_u);
	return my_fmax(my_fmin(t_theo,0.99*t_u), 0.01*t_u);
}


/*  compute maximum of a vector */
 double compute_max_Nnc( const double* x) {
	unsigned int ii = 0;
	double m = -100.0;

	for (ii=60; ii--; ) {
		m = my_fmax(m,x[ii]);
	}
	return m;
}


/* main function of the algorithm */
int proposed_algorithm(const double* x0, double* u, const double* xref, const double* uref, double* x, double* fval, unsigned int* iter, unsigned int* iter_ls, double* optimval, double* feasval, double* meritval){

	unsigned int conditions_f = 1, conditions_x = 1,
		conditions_n = 1, cond_err = 1,
		it = 0, it_ls = 0, ii = 0, jj = 0;
	int cond = -2;
	unsigned int flag_ini = 0;
	double J= 0.0, Jt = 0.0, constr_viol = 1.0;
	double xp[90],
		dot_J[30], du[30], up[30],
		gps[60], gpsp[60],
		sl[60], sl_sqr[60], dsl[60], slp[60], slp_sqr[60], muG[60],
		mu[60], dm[60], mup[60],
		dsl_sqr = 0.0, gps_sqr = 0.0, dm_sqr = 0.0, rho_hat_tmp = 0.0,
		rho = 0.0, rho_hat = 0.0,
		phi0 = 0.0, phit = 0.0, phi0_dot = 0.0,
		t = 0.0, t_u = 0.0,
		du_sqr = 0.0, u_sqr = 0.0;

/* Initialization of the predicted state x (dimension: N*nx) */
	det_x(x0,u,x);

	/* This function initializes the slack variables sl, its squares sl_sqr 
and the function gps = g + 0.5 * diag(sl_sqr) */
	initialize_slack(u, sl, sl_sqr, gps);


	/* Outer loop over the iterates it. The break command interrups the cycle 
	if convergence or error is encountered before termination */
	for (it=0; it < CONTROLADOR_FALC_OPT_MAXIMUM_ITERATIONS; it++) {

		/* Computation of the cost J and its derivative dot_J */
		det_J_and_dot_J(x0, u, x, xref, uref, &J, dot_J);


		/* Compute the gradient step, projected onto the linearization of the constraints 
du: primal variable gradient steps, dsl: slack variables, muG: dual variables */
		gradient_step(dot_J, u, sl, sl_sqr, gps, du, dsl, muG);

		dot_product_Nnu(&du_sqr, du, du);
		dot_product_Nnu(&u_sqr, u, u);
		dot_product_Nnc(&gps_sqr, gps, gps);
		dot_product_Nnc(&dsl_sqr, dsl, dsl);
		if (it==0)
			copy_Nnc(mu, muG);

		/* dm = muG - mu; */
		diff_Nnc(dm,muG,mu);


		/* compute the directional derivative of the merit function, phi0_dot */
		det_dot_phi (du,dot_J, rho, gps, mu, dm, &phi0_dot);

		/* Check the penalty parameter rho via phi0_dot */
		if (conditions_rho_PM_simpler(phi0_dot,du_sqr,dsl_sqr,CONTROLADOR_FALC_OPT_LINE_SEARCH_ALPHA_MAX)==0) {
			dot_product_Nnc(&dm_sqr,dm,dm);
			rho_hat_tmp = dm_sqr / gps_sqr;

			/* Update the penalty parameter rho */
			rho_hat = 2.0 * sqrt(rho_hat_tmp);
			rho = my_fmax(2.0*rho,rho_hat);
			flag_ini = 1;

			/* Recompute phi0_dot */
			det_dot_phi (du,dot_J, rho, gps, mu, dm, &phi0_dot);
		}
		if ((flag_ini == 1)||(it == 0)){
			det_phi (J, gps, mu, rho, &phi0);
			flag_ini = 0;
		}

		if (phi0_dot <= -CONTROLADOR_FALC_OPT_LINE_SEARCH_ALPHA_TOL_SQ) {
			t = 1.0;
			t_u = 1.0;
			for (it_ls = 0; it_ls < CONTROLADOR_FALC_OPT_MAXIMUM_LINE_SEARCH_ITERATIONS; it_ls++) {

				/* up = u + t*du; */
				/* slp = sl + t*dsl; */
				weighted_sum_Nnu(up,du,u,&t);
				weighted_sum_Nnc(slp,dsl,sl,&t);

				/* mup = mu + t*dm; */
				weighted_sum_Nnc(mup,dm,mu,&t);


				/* Compute the predicted state xp (dim. N*nx) using the new predicted input up (dim. N*nu) */
				det_x(x0,up,xp);

				/* Compute the new cost Jt, square the slack variables slp 

				and compute gpsp = gps + 0.5 * slp_sqr */
				det_J(x0, up, xp, xref, uref, &Jt);

				build_sqr_Nnc(slp, slp_sqr);
				build_gpsl(up,slp_sqr, gpsp);

				/* Compute the merit function phit (function of the step size t) */
				det_phi (Jt,gpsp,mup,rho,&phit);

				/* Check the Armijo condition */

				if (phit - phi0 <= 0.3*t*phi0_dot){

					/* step size t accepted. Break the line-search */
					break;
				}
				else {
					t_u = t;

					/* Reduce the step size t via quadratic interpolation */
					t = quadratic_interp (phi0, phi0_dot, t_u, phit);
				}
			}

			/* if t gets too small, output an error */
			if (t_u <= 0.0001)
				cond_err = 0;
		}
		else {

			/* phi0_dot is sufficiently small and we have converged */
			conditions_f = 0;
		}

		iter_ls[it] = it_ls+1;
		if (it_ls == CONTROLADOR_FALC_OPT_MAXIMUM_LINE_SEARCH_ITERATIONS)
			cond_err = 0;


		/* update step */
		copy_Nnc(mu,mup);
		copy_Nnu(u,up);
		copy_Nnc(sl,slp);
		copy_Nnc(sl_sqr,slp_sqr);
		copy_Nnc(gps,gpsp);
		copy_Nnx(x,xp);
		J = Jt;
		phi0 = phit;


		/* check convergence (update step) */

		constr_viol = compute_max_Nnc(gpsp); 
		if ((du_sqr >= CONTROLADOR_FALC_OPT_LINE_SEARCH_ALPHA_SQ_TOL_SQ)||(constr_viol >= CONTROLADOR_FALC_OPT_KKTOPTIMALITY_TOLERANCE))
			conditions_x = 1;
		else
		conditions_x = 0;
		if (it == CONTROLADOR_FALC_OPT_MAXIMUM_ITERATIONS-1)
			conditions_n = 0;
		if(!((conditions_f && cond_err) && (conditions_n && conditions_x)))
			break;

	}


	/* Assign optimality */
	dot_product_Nnc(&dsl_sqr, dsl, dsl);
	(*optimval) = (sqrt(du_sqr + dsl_sqr))/CONTROLADOR_FALC_OPT_LINE_SEARCH_ALPHA_MAX; 
	(*feasval) = constr_viol; 
	(*meritval) = phi0_dot; 

	/* Assign exitflag */
	(*iter) = it+1;
	if (conditions_f == 0)
		cond = 2;
	else {
		if (conditions_x == 0)
			cond = 1;
		else {
			if (conditions_n == 0) /* Maximum number of iterations reached */
				cond = 0;
			else {
				if (phi0_dot < -CONTROLADOR_FALC_OPT_LINE_SEARCH_ALPHA_TOL) 
					cond = -1;
				else
					cond = -10;
			}
		}
	}


		/* Output cost function */
	*fval = J;
	return cond;
} 

