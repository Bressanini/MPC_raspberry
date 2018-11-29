#include "./myBibs/utilidades.h"
#include "./myBibs/modelo.h"

#include "./generatedCode/controladorFalcOpt.h"

int printiter = 0;

void inicializaFalcOpt(int NX, int NU,double** f_x0, double** f_u, double** f_xref, double** f_uref, double** f_x, double** f_fval, unsigned int** f_iter, unsigned int** f_iter_ls, double** f_optimval, double** f_feasval, double** f_meritval){

	int i,j;

	*f_x0 		= (double*		 )malloc(sizeof(double		 )*(NX+NU)	   );
	*f_x 	 	= (double*		 )malloc(sizeof(double		 )*(NX+NU)*(15));
	*f_u  		= (double*		 )malloc(sizeof(double		 )*(NU   )*(15));
	*f_xref  	= (double*		 )malloc(sizeof(double		 )*(NX+NU)*(15));
	*f_uref	 	= (double*		 )malloc(sizeof(double		 )*(NU   )*(15));
	*f_fval  	= (double*		 )malloc(sizeof(double		 )			   );
	*f_iter  	= (unsigned int* )malloc(sizeof(unsigned int )			   );
	*f_iter_ls	= (unsigned int* )malloc(sizeof(unsigned int )*8000		   );
	*f_optimval	= (double*		 )malloc(sizeof(double		 )			   );
	*f_feasval	= (double*		 )malloc(sizeof(double		 )			   );
	*f_meritval	= (double*		 )malloc(sizeof(double		 )			   );

	for (i=0; i<(NX+NU)		; i++) (*f_x0)[i] = 0.0;
	for (i=0; i<(NX+NU)*(15); i++) (*f_x )[i] = 0.0;
	for (i=0; i<(NU   )*(15); i++) (*f_u )[i] = 0.0;

	for (i=0; i<(NX+NU)*(15); i++) (*f_xref )[i] = 0.0;
	for (i=0; i<(NU   )*(15); i++) (*f_uref )[i] = 0.0;

	**f_fval 	 	= 0.0;
	**f_iter 	 	= 0;
	**f_iter_ls 	= 0;
	**f_optimval	= 0.0;
	**f_feasval		= 0.0;
	**f_meritval	= 0.0;

	return;


}

int main(int argc, char **argv){

// Declaração de variáveis ---------------------------------------------------------------------------------

int NX, 	// Número de estados do modelo da planta
	NU, 	// Número de entradas do modelo da planta
	NY; 	// Número de saídas do modelo da planta

double 	*A,	// Matriz A do modelo SS-LIN discreto
		*B,	// Matriz B do modelo SS-LIN discreto
		*U,	// Matriz U com as entradas para a simulação
		*Uref,
		*X,
		*Xref;
		
double TS; 	// Passo de integração do modelo da planta
int	   NIT;	// Número de iterações da simulação

int num_iters; // Número de iterações do solver

int i,j,k;

// Inicialização das variáveis da simulação ----------------------------------------------------------------

inicializa(&A, &B, &NX, &NU, &TS, &NIT, &U); // Alocando matrizes e inicializando variáveis de acordo com arquivo de parâmetros

// Declaração e definicção de entradas e saídas da simulação -----------------------------------------------

X  = mZero(NX,NIT+1); //  Alocando e inicializando com zeros a matriz X que guardará as saídas da simulação + condição inicial

free(U);
U = mZero(NU,NIT);

double *u0, *x0, *x0sim;
csv2mat(&u0	,NU 	,1	,"sim/u0.csv");
csv2mat(&x0 	,NX	,1	,"sim/x0.csv");
csv2mat(&x0sim	,NX	,1	,"sim/x0sim.csv");

double dxdt[4], x[4], u[2]; // Variáveis de apoio usados no modelo não linear

for (i=0;i<NX;i++) X[(NIT+1)*i+0] = x0sim[i]; // Passando condição inicial da simulação

// Controlador ---------------------------------------------------------------------------------------------

double *f_x0, *f_u, *f_xref, *f_uref, *f_x, *f_fval, *f_optimval, *f_feasval, *f_meritval;
unsigned int *f_iter, *f_iter_ls;
int flag;

inicializaFalcOpt(NX, NU, &f_x0, &f_u, &f_xref, &f_uref, &f_x, &f_fval, &f_iter, &f_iter_ls, &f_optimval, &f_feasval, &f_meritval);

csv2mat(&Xref,NX,NIT,"sim/xREF.csv");
csv2mat(&Uref,NU,NIT,"sim/uREF.csv");

// Repetição -----------------------------------------------------------------------------------------------

int num_repet = 10;
int l;

// Timing --------------------------------------------------------------------------------------------------

clock_t begin, end;
double time_spent;
double tproc[NIT*num_repet];

double *X1, *U1, *U2;
X1 = mZero(num_repet,NIT+1);
U1 = mZero(num_repet,NIT);
U2 = mZero(num_repet,NIT);

// Loop da simulação ---------------------------------------------------------------------------------------

for(l=0;l<num_repet;l++){

for (k=0; k<NIT;k++){

	for (i=0; i<NX; i++) f_x0[i] = X[i*(NIT+1)+k]; 						// Passando medição atual
	for (i=0; i<NU; i++){
		if (k<=0) 	f_x0[i+NX] = 0;
		else		f_x0[i+NX] = U[i*(NIT)+k-1];
	}

	for (i=0;i<NX;i++) for (j=0; j<15; j++) f_xref[NX*j+i]  = Xref[NIT*i+k]; // Passando referência de x

	begin = clock();

	flag = proposed_algorithm(f_x0, f_u, f_xref, f_uref, f_x, f_fval, f_iter, f_iter_ls, f_optimval, f_feasval, f_meritval);

	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	tproc[k+l*NIT] = time_spent;

	if (printiter){
		printf("%3d|\tTempo do cálculo do controlador: %lf\n", k, time_spent);
		printf("  t|\t%3.1f\n",k*TS+TS);
		printf("   |\tNúmero de iterações: %d\n", *f_iter);
		printf("\n");
	}

	// Obtendo resultado do MPC
	if (k <= 0)	for (i=0; i<NU; i++) U[i*NIT+k] = f_u[i] + 0;
	else		for (i=0; i<NU; i++) U[i*NIT+k] = f_u[i] + U[i*(NIT)+k-1];

	// Aplicando entradas no modelo linear
	AXmBU(A,X+k,B,U+k,NX,NU,NIT);

	/*for (j=0;j<NX;j++) x[j] = X[(NIT+1)*j+i];
	for (j=0;j<NU;j++) u[j] = U[NIT*j+i];

	quatroTQnl(dxdt,x,u);

	for (j=0;j<NX;j++) X[(NIT+1)*j+i+1] = dxdt[j]*TS + X[(NIT+1)*j+i];*/

}
for (i=0;i<NIT;i++){
	X1[(NIT+1)*l+i] = X[i];
	U1[NIT*l+i] = U[i    ];
	U2[NIT*l+i] = U[NIT+i];
}
X1[(NIT+1)*l+NIT] = X[NIT];
printf("l:\t%d\n",l);
} // Loop da repetição


for (i=0;i<4;i++) printf("X(%d) ini = %lf\n",i,X[(NIT+1)*i+0]); printf("\n");
for (i=0;i<4;i++) printf("X(%d) end = %lf\n",i,X[(NIT+1)*i+NIT]);

mat2csv(tproc,num_repet ,NIT  ,"sim/tproc.csv");
mat2csv(X1,num_repet ,NIT+1  ,"sim/X1.csv");
mat2csv(U1,num_repet ,NIT  ,"sim/U1.csv");
mat2csv(U2,num_repet ,NIT  ,"sim/U2.csv");


}


