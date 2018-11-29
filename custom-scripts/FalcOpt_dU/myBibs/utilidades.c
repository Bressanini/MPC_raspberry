#include "utilidades.h"

void inicializa(double** A, double** B,  int* NX, int* NU, double* TS, int* NIT, double** U){
// Função que carrega os parâmetros da simulação a partir do arquivo dadosSIM.csv
//
//	Parâmetros
//		A	: Ponteiro para ponteiro para double que irá representar a matriz A do modelo SS-LIN
//		B	: Ponteiro para ponteiro para double que irá representar a matriz B do modelo SS-LIN
//		NX	: Ponteiro para inteiro que irá representar o número de estados do modelo SS-LIN
//		NU	: Ponteiro para inteiro que irá representar o número de entradas do modelo SS-LIN
//		TS	: Ponteiro para double que irá representar o tempo de amostragem do modelo SS-LIN discreto
//		NIT : Ponteiro para inteiro que irá representar o número de iterações da simulação
//		U	:  Ponteiro para ponteiro para double que irá representar a matriz U que contém as entradas para a simulação
//
//	Observações:
//		* As matrizes A e B e o vetor U não devem estar alocados
//		* Lembrar de desalocar A, B e U
//

	int  i																	// Contador de linha
		,j;																	// Contador de coluna

	FILE *ARQ = fopen("./sim/dadosSIM.csv","r"); 							// Abrindo o arquivo com os parâmetros da simulação
	int dummy;
	
	if (ARQ==NULL) printf("Erro ao abrir o arquivo ./sim/dadosSIM.csv\n"); 	// Caso tenha acontecido algum erro ao abrir o arquivo, uma mensagem é impressa na tela
	else{ 																	// Caso não tenha ocorrido um erro ao abrir o arquivo
		
		dummy = fscanf(ARQ,"%*[^;];%d;%*[^\n]\n",NX);								// 		Lendo o número de estados
		dummy = fscanf(ARQ,"%*[^;];%d;%*[^\n]\n",NU);								// 		Lendo o número de entradas
		dummy = fscanf(ARQ,"%*[^;];%lf;%*[^\n]\n",TS);								// 		Lendo o tempo de amostragem do modelo
		dummy = fscanf(ARQ,"%*[^;];%d;%*[^\n]\n",NIT);								// 		Lendo o número de iterações da simulação a serem realizadas
		
		*A = (double*)malloc(sizeof(double)*(*NX)*(*NX));					// 		Alocando a matriz A
		*B = (double*)malloc(sizeof(double)*(*NX)*(*NU));					// 		Alocando a matriz B
		*U = (double*)malloc(sizeof(double)*(*NU)*(*NIT));					// 		Alocando a matriz com as entradas a serem aplicadas
		
		if ( (*A==NULL) || (*B==NULL) || (*U==NULL) ){ 						// 		Caso tenha acontecido algum erro na alocação das matrizes
			printf("Erro ao alocar matrizes da simulação\n\n");				//			Uma mensagem é exibida na tela
		}
		else{																// 		Caso não tenha ocorrido erro na alocação das matrizes
			// Leitura da matriz A -------------------------------------------------------------------------------------------------------------------------------
			for (i=0;i<*NX;i++){											//			Para cada linha da matriz A
				dummy = fscanf(ARQ,"%*[^;];");										//				O arquivo é varrido até o próximo caracter ;
				for (j=0;j<*NX;j++) dummy = fscanf(ARQ,"%lf;",*A + i*(*NX) + j);	//				Para cada coluna da matriz A, cada elemento da matriz A é lido
				dummy = fscanf(ARQ,"%*[^\n]\n");									//				O arquivo é varrido até a próxima quebra de linha
			}
			// Leitura da matriz  B ------------------------------------------------------------------------------------------------------------------------------
			for (i=0;i<*NX;i++){											//			Para cada linha da matriz B
				dummy = fscanf(ARQ,"%*[^;];");										//				O arquivo é varrido até o próximo caracter ;
				for (j=0;j<*NU;j++) dummy = fscanf(ARQ,"%lf;",*B + i*(*NU) + j);	//				Para cada coluna da matriz B, cada elemento da matriz A é lido
				dummy = fscanf(ARQ,"%*[^\n]\n");									//				O arquivo é varrido até a próxima quebra de linha
			}
			// Leitura de U --------------------------------------------------------------------------------------------------------------------------------------
			for (i=0;i<*NU;i++){											//			Para cada linha da matriz B
				dummy = fscanf(ARQ,"%*[^;];");										//				O arquivo é varrido até o próximo caracter ;
				for (j=0;j<*NIT;j++) dummy = fscanf(ARQ,"%lf;",*U + i*(*NIT) + j);	//				Para cada coluna da matriz B, cada elemento da matriz A é lido
				dummy = fscanf(ARQ,"%*[^\n]\n");									//				O arquivo é varrido até a próxima quebra de linha
			}
		
		}
	
		fclose(ARQ); 														// Fechando o arquivo
	}
		
	return;
}

void mat2csv(double* M, int NL, int NC, char NOME[20]){
// Função para imprimir em um arquivo CSV o conteúdo de uma matriz
//
//	Parâmetros
//		M 	: Ponteiro para double que representa a matriz a ser salva
//		NL	: Número de linhas da matriz a ser salva
//		NC	: Número de colunas da matriz a ser salva
//		NOME: String contendo o nome do arquivo a ser salvo
//

	FILE *ARQ = fopen(NOME,"w");									//	Abrindo o arquivo
	
	int  i															// Contador de linha
		,j;															// Contador de coluna
	
	if (ARQ==NULL) printf("Erro ao abrir o arquivo %s\n",NOME);		// Caso ocorra erro ao abrir o arquivo, uma mensagem é impressa na tela
	else{															// Caso não tenha ocorrido um erro ao abrir o arquivo

		for (i=0;i<NL;i++){											// 		Para cada linha da matriz 
			for (j=0;j<NC;j++) fprintf(ARQ,"%lf,",*(M+i*NC+j));	//			Para cada coluna da matriz, cada elemento da matriz é gravado no arquivo
			fprintf(ARQ,"\n");										//			É gravado no arquivo uma quebra de linha
		}
		
		fclose(ARQ);												// 		O arquivo é fechado
	}
	
	return;
}

void csv2mat(double** M, int NL, int NC, char NOME[20]){

	FILE *ARQ = fopen(NOME,"r");									//	Abrindo o arquivo
		
		int  i															// Contador de linha
			,j;															// Contador de coluna
		int dummy;

	*M = (double*)malloc(sizeof(double)*NL*NC);	
		
	if (ARQ==NULL) printf("Erro ao abrir o arquivo %s\n",NOME);		// Caso ocorra erro ao abrir o arquivo, uma mensagem é impressa na tela
	else if(*M==NULL) printf("Erro ao alocar memória\n");			//  Caso ocorra erro ao alocar memória, uma mensagem é impressa na tela
	else{															// Caso não tenha ocorrido um erro ao abrir o arquivo ou alocar memória

		for (i=0;i<NL;i++){											// 		Para cada linha da matriz 
			for (j=0;j<NC;j++) dummy =fscanf(ARQ,"%lf,",*M+i*NC+j);	//			Para cada coluna da matriz, cada elemento da matriz é gravado no arquivo
			fprintf(ARQ,"\n");										//			É gravado no arquivo uma quebra de linha
		}

		fclose(ARQ);												// 		O arquivo é fechado
	}
	
	return;
}

double* mZero(int NL, int NC){
// Função que aloca e inicializa uma matriz com zeros, retornando o endereço da matriz
//
//	Parâmetros
//		NL	: Número de linhas da matriz
//		NC	: Número de colunas da matriz
//
//	Observações:
//		* Lembrar de desalocar a memória

	int 	i;														// Contador
	double *M;														// Ponteiro para a matriz a ser alocada
	
	M = (double*)malloc(sizeof(double)*NL*NC); 						// Alocando o espaço da matiz
	
	if (M==NULL) printf("Erro ao alocar memória!\n");				// Caso tenha ocorrido algum erro na alocação da matriz, uma mensagem é impressa na tela
	else{															// Caso não tenha ocorrido algum erro na alocação da matriz
		
		for (i=0; i<NL*NC; i++) *(M+i) = 0.0;							// 		Para cada uma das posições de memória ocupadas pela matriz, cada elemento é inicializado com zero
		
		return(M);													// 		Retornando o endereço da matriz alocada
	}

}

void AXmBU(double* A, double* X, double* B, double* U, int NX, int NU, int NIT){
// Função que avalia a equação de estado X[k+1] = A*X[k] + B*U[k] 
//
// Parâmetros
//		A	: Matriz dinâmica do sistema
//		X	: Vetor com os valores dos estados
//		B	: Matriz entrada-estado
//		U	: Entradas do sistema
//		NX	: Número de estados do sistema
//		NU	: Número de entradas do sistema
//		NIT : Número total de iterações da simulação
//
// Observações:
//		
//		* Assumindo que todas as matrizes e vetores estão na forma Row-major
//		* A função cblas_dgemv() é utilizada em duas etapas:
//			- X[k+1] = alfa*A*X[k] + beta*X[k+1], sendo alfa = 1, beta = 0;
//			- X[k+1] = alfa*B*U[k] + beta*X[k+1], sendo alfa = 1, beta = 1;
//

	cblas_dgemv(CblasRowMajor	, 	// Indicação se as matrizes são Column-major (fortran) ou Row-major (C)
				CblasNoTrans 	,	// Indicação se a matriz A deve ser transposta
				NX           	,	// Número de linhas da mastriz A
				NX           	,	// Número de colunas da mastriz A
				1.0          	,	// Fator alfa
				A				,	// Matriz A
				NX				,	// Tamamnho da primeira dimensão do vetor X. Se Row-major, deve ser nº de col
				X				,	// Vetor dos estados X
				NIT+1			,	// Tamamnho da primeira dimensão do vetor X. Se Row-major, deve ser nº de col
				0.0				,	// Fator Beta
				X+1				,	// Vetor X[k+1] = alfa*A*X[k] + beta*X[k+1]
				NIT+1			);	// Tamamnho da primeira dimensão do vetor X[k+1]. Se Row-major, deve ser nº de col);

	cblas_dgemv(CblasRowMajor	, 	// Indicação se as matrizes são Column-major (fortran) ou Row-major (C)
				CblasNoTrans	,	// Indicação se a matriz A deve ser transposta
				NX           	,	// Número de linhas da matriz B
				NU           	,	// Número de colunas das mastrizes A
				1.0         	,	// Fator alfa
				B				,	// Matriz B
				NU				,	// Tamamnho da primeira dimensão da matriz B. Se Row-major, deve ser nº de col
				U				,	// Vetor U
				NIT				,	// Tamamnho da primeira dimensão do vetor U. Se Row-major, deve ser nº de col
				1.0				,	// Fator Beta
				X+1				,	// Vetor X[k+1] = alfa*B*U[k] + beta*X[k+1]
				NIT+1			);	// Tamamnho da primeira dimensão do vetor X[k+1]. Se Row-major, deve ser nº de col);
	
	return;
	
}
