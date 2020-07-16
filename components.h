float calculateLocalD(int i,mesh m){
    Matrix matrix;
    Vector row1, row2, row3;

    element e = m.getElement(i);
    
    row1.push_back(calcularTenedor(e, EQUIS, 2, 1, m));
    row1.push_back(calcularTenedor(e, YE, 2, 1, m));
    row1.push_back(calcularTenedor(e, ZETA, 2, 1, m));

    row2.push_back(calcularTenedor(e, EQUIS, 3, 1, m));
    row2.push_back(calcularTenedor(e, YE, 3, 1, m));
    row2.push_back(calcularTenedor(e, ZETA, 3, 1, m));

    row3.push_back(calcularTenedor(e, EQUIS, 4, 1, m));
    row3.push_back(calcularTenedor(e, YE, 4, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 4, 1, m));

    matrix.push_back(row1);
    matrix.push_back(row2);
    matrix.push_back(row3);

    return determinant(matrix);
}

void calculateAlpha(int i,Matrix &A,mesh m){
    zeroes(A,3);
    element e = m.getElement(i);
    
    A.at(0).at(0) = OperarRestaTenedor(e, YE, ZETA, 3, 4, m);
    A.at(0).at(1) = OperarRestaTenedor(e, YE, ZETA, 4, 2, m);
    A.at(0).at(2) = OperarRestaTenedor(e, YE, ZETA, 2, 3, m);

    A.at(1).at(0) = OperarRestaTenedor(e, EQUIS, ZETA, 4, 3, m);
    A.at(1).at(1) = OperarRestaTenedor(e, EQUIS, ZETA, 2, 4, m);
    A.at(1).at(2) = OperarRestaTenedor(e, EQUIS, ZETA, 3, 2, m);

    A.at(2).at(0) = OperarRestaTenedor(e, EQUIS, YE, 3, 4, m);
    A.at(2).at(1) = OperarRestaTenedor(e, EQUIS, YE, 4, 2, m);
    A.at(2).at(2) = OperarRestaTenedor(e, EQUIS, YE, 2, 3, m);

}

void calculateBeta(Matrix &B){
    zeroes(B,3,12);

    B.at(0).at(0) = -1; 
    B.at(0).at(1) =  1; 
    B.at(0).at(2) =  0; 
    B.at(0).at(3) =  0; 
    B.at(0).at(4) = -1; 
    B.at(0).at(5) =  1; 
    B.at(0).at(6) =  0; 
    B.at(0).at(7) =  0; 
    B.at(0).at(8) = -1; 
    B.at(0).at(9) =  1; 
    B.at(0).at(10) =  0; 
    B.at(0).at(11) =  0; 
    
    B.at(1).at(0) = -1; 
    B.at(1).at(1) = 0; 
    B.at(1).at(2) = 1;
    B.at(1).at(3) = 0;
    B.at(1).at(4) = -1; 
    B.at(1).at(5) = 0; 
    B.at(1).at(6) = 1;
    B.at(1).at(7) = 0;
    B.at(1).at(8) = -1; 
    B.at(1).at(9) = 0; 
    B.at(1).at(10) = 1;
    B.at(1).at(11) = 0;

    B.at(2).at(0) = -1;
    B.at(2).at(1) = 0;
    B.at(2).at(2) = 0;
    B.at(2).at(3) = 1;
    B.at(2).at(4) = -1;
    B.at(2).at(5) = 0;
    B.at(2).at(6) = 0;
    B.at(2).at(7) = 1;
    B.at(2).at(8) = -1;
    B.at(2).at(9) = 0;
    B.at(2).at(10) = 0;
    B.at(2).at(11) = 1;
}

void calculateOmega(Matrix &C){
    zeroes(C,3,4);
    C.at(0).at(0) = -1; C.at(0).at(1) = 1; C.at(0).at(2) = 0; C.at(0).at(3) = 0;
    C.at(1).at(0) = -1; C.at(1).at(1) = 0; C.at(1).at(2) = 1; C.at(1).at(3) = 0;
    C.at(2).at(0) = -1; C.at(2).at(1) = 0; C.at(2).at(2) = 0; C.at(2).at(3) = 1;
}


void calculateGamma(Matrix &m){
	zeroes(m,12,3);

	m.at(0).at(0) = 1;   m.at(0).at(1) = 0;   m.at(0).at(2) = 0;
	m.at(1).at(0) = 1;   m.at(1).at(1) = 0;   m.at(1).at(2) = 0; 
    m.at(2).at(0) = 1;   m.at(2).at(1) = 0;   m.at(2).at(2) = 0;
	m.at(3).at(0) = 1;   m.at(3).at(1) = 0;   m.at(3).at(2) = 0; 
    m.at(4).at(0) = 0;   m.at(4).at(1) = 1;   m.at(4).at(2) = 0;
	m.at(5).at(0) = 0;   m.at(5).at(1) = 1;   m.at(5).at(2) = 0; 
    m.at(6).at(0) = 0;   m.at(6).at(1) = 1;   m.at(6).at(2) = 0;
	m.at(7).at(0) = 0;   m.at(7).at(1) = 1;   m.at(7).at(2) = 0; 
    m.at(8).at(0) = 0;   m.at(8).at(1) = 0;   m.at(8).at(2) = 1;
	m.at(9).at(0) = 0;   m.at(9).at(1) = 0;   m.at(9).at(2) = 1; 
    m.at(10).at(0) = 0;  m.at(10).at(1) = 0;  m.at(10).at(2) = 1;
	m.at(11).at(0) = 0;  m.at(11).at(1) = 0;  m.at(11).at(2) = 1; 
	
}

float calculateLocalJ(int i,mesh m){
    Matrix matrix;
    Vector row1, row2, row3;

    element e = m.getElement(i);
    
    row1.push_back(calcularTenedor(e, EQUIS, 2, 1, m));
    row1.push_back(calcularTenedor(e, EQUIS, 3, 1, m));
    row1.push_back(calcularTenedor(e, EQUIS, 4, 1, m));

    row2.push_back(calcularTenedor(e, YE, 2, 1, m));
    row2.push_back(calcularTenedor(e, YE, 3, 1, m));
    row2.push_back(calcularTenedor(e, YE, 4, 1, m));

    row3.push_back(calcularTenedor(e, ZETA, 2, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 3, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 4, 1, m));

    matrix.push_back(row1);
    matrix.push_back(row2);
    matrix.push_back(row3);

    return determinant(matrix);
}

//matriz A
void calcularT(int i, mesh m,Matrix t){
	float t1,t2,t3,t4;
	element e = m.getElement(i);
	
	t1=calcularTenedor(e,EQUIS,2,1,m)+(2*calcularTenedor(e,EQUIS,3,1,m))+calcularTenedor(e,EQUIS,4,1,m)+calcularTenedor(e,EQUIS,4,1,m);
	t2=(2*calcularTenedor(e,EQUIS,2,1,m))+calcularTenedor(e,EQUIS,3,1,m)+calcularTenedor(e,EQUIS,4,1,m)+calcularTenedor(e,EQUIS,4,1,m);
	t3=calcularTenedor(e,EQUIS,2,1,m)+calcularTenedor(e,EQUIS,3,1,m)+(2*calcularTenedor(e,EQUIS,4,1,m))+calcularTenedor(e,EQUIS,4,1,m);
	t4=calcularTenedor(e,EQUIS,2,1,m)+calcularTenedor(e,EQUIS,3,1,m)+calcularTenedor(e,EQUIS,4,1,m)+(2*calcularTenedor(e,EQUIS,4,1,m));
	
	t.at(0).at(0)=t1;
	t.at(1).at(0)=t2;
	t.at(2).at(0)=t3;
	t.at(3).at(0)=t4;
	t.at(4).at(0)=0;
	t.at(5).at(0)=0;
	t.at(6).at(0)=0;
	t.at(7).at(0)=0;
	t.at(8).at(0)=0;
	t.at(9).at(0)=0;
	t.at(10).at(0)=0;
	t.at(11).at(0)=0;
	
	t.at(0).at(1)=0;
	t.at(1).at(1)=0;
	t.at(2).at(1)=0;
	t.at(3).at(1)=0;
	t.at(4).at(1)=t1;
	t.at(5).at(1)=t2;
	t.at(6).at(1)=t3;
	t.at(7).at(1)=t4;
	t.at(8).at(1)=0;
	t.at(9).at(1)=0;
	t.at(10).at(1)=0;
	t.at(11).at(1)=0;
	
	t.at(0).at(2)=0;
	t.at(1).at(2)=0;
	t.at(2).at(2)=0;
	t.at(3).at(2)=0;
	t.at(4).at(2)=0;
	t.at(5).at(2)=0;
	t.at(6).at(2)=0;
	t.at(7).at(2)=0;
	t.at(8).at(2)=t1;
	t.at(9).at(2)=t2;
	t.at(10).at(2)=t3;
	t.at(11).at(2)=t4;	
	
}

//matriz G
void calcularT2(int i, mesh m,Matrix t){
	float t1,t2,t3,t4;
	
	element e = m.getElement(i);
	
	node n1= selectNode(1,e,m);
	
		
	t1=calcularTenedor(e,YE,2,1,m)+(2*calcularTenedor(e,YE,3,1,m))+calcularTenedor(e,YE,4,1,m)+calcularTenedor(e,YE,4,1,m);
	t2=(2*calcularTenedor(e,YE,2,1,m))+calcularTenedor(e,YE,3,1,m)+calcularTenedor(e,YE,4,1,m)+calcularTenedor(e,YE,4,1,m);
	t3=calcularTenedor(e,YE,2,1,m)+calcularTenedor(e,YE,3,1,m)+(2*calcularTenedor(e,YE,4,1,m))+calcularTenedor(e,YE,4,1,m);
	t4=calcularTenedor(e,YE,2,1,m)+calcularTenedor(e,YE,3,1,m)+calcularTenedor(e,YE,4,1,m)+(2*calcularTenedor(e,YE,4,1,m));	
	
	
	t.at(0).at(0)=t1;
	t.at(1).at(0)=t2;
	t.at(2).at(0)=t3;
	t.at(3).at(0)=t4;
	t.at(4).at(0)=0;
	t.at(5).at(0)=0;
	t.at(6).at(0)=0;
	t.at(7).at(0)=0;
	t.at(8).at(0)=0;
	t.at(9).at(0)=0;
	t.at(10).at(0)=0;
	t.at(11).at(0)=0;
	
	t.at(0).at(1)=0;
	t.at(1).at(1)=0;
	t.at(2).at(1)=0;
	t.at(3).at(1)=0;
	t.at(4).at(1)=t1;
	t.at(5).at(1)=t2;
	t.at(6).at(1)=t3;
	t.at(7).at(1)=t4;
	t.at(8).at(1)=0;
	t.at(9).at(1)=0;
	t.at(10).at(1)=0;
	t.at(11).at(1)=0;
	
	t.at(0).at(2)=0;
	t.at(1).at(2)=0;
	t.at(2).at(2)=0;
	t.at(3).at(2)=0;
	t.at(4).at(2)=0;
	t.at(5).at(2)=0;
	t.at(6).at(2)=0;
	t.at(7).at(2)=0;
	t.at(8).at(2)=t1;
	t.at(9).at(2)=t2;
	t.at(10).at(2)=t3;
	t.at(11).at(2)=t4;
		
}

//matriz K
float calcularKTermino(int i, mesh m){
	float k;
	node n1;
	element e = m.getElement(i);
	
	n1= selectNode(1,e,m);
	
	k=calcularTenedor(e,YE,4,1,m)+calcularTenedor(e,EQUIS,2,1,m),calcularTenedor(e,EQUIS,3,1,m)+
	calcularTenedor(e,EQUIS,4,1,m)+(4*selectCoord(EQUIS,n1))+calcularTenedor(e,YE,2,1,m)+
	calcularTenedor(e,YE,3,1,m)+(4*selectCoord(YE,n1));
	
	return k;

}

Matrix createLocalM(int e,mesh &m){
    Matrix matrixA,matrixK,matrixG,matrixD,matrixT,matrixT2;
    float u_bar,nu,rho,Ve,J,Determinant;
    
    /* [ A+K  G ]
       [  D   0 ]
    */

    //Matrix A
    Matrix g_matrix, Alpha, Beta;
	zeroes(matrixT,12,3);
	zeroes(matrixT2,12,3);
    Determinant = calculateLocalD(e,m);
    J = calculateLocalJ(e,m);
    calcularT(e,m,matrixT);

    calcularT2(e,m,matrixT2);

    if(Determinant == 0){
        cout << "\n!---CATASTROPHIC FAILURE---!\n";
        exit(EXIT_FAILURE);
    }
    
    float real_a = (float) (J)/(120*Determinant);
    
    calculateGamma(g_matrix);
    calculateAlpha(e,Alpha,m);
    calculateBeta(Beta);
    productRealMatrix(real_a, productMatrixMatrix(matrixT,productMatrixMatrix(Alpha,Beta,3,3,12),12,3,12),matrixA);


    //Matrix K
    Matrix Alpha_t,Beta_t;

    
    float real_k = (float) (calcularKTermino(e,m)*J)/(120*Determinant);

    transpose(Alpha,Alpha_t);
    transpose(Beta,Beta_t);

    productRealMatrix(real_k,productMatrixMatrix(Beta_t,productMatrixMatrix(Alpha_t,productMatrixMatrix(Alpha,Beta,3,3,12),3,3,12),12,3,12),matrixK);

    
    //Matrix G
    Matrix Omega;
    
    float real_g = (float) (J/(120*Determinant));

    calculateOmega(Omega);
    productRealMatrix(real_g,productMatrixMatrix(matrixT2,productMatrixMatrix(Alpha,Omega,3,3,4),12,3,4),matrixG);

    //Matrix D
    Matrix g_matrix_t,Omega_t;
    float real_d = (float)(J/(120*Determinant));

    transpose(Omega, Omega_t);
    transpose(g_matrix,g_matrix_t);
    productRealMatrix(real_d,productMatrixMatrix(Omega_t,productMatrixMatrix(Alpha_t,g_matrix_t,3,3,12),4,3,12),matrixD);

    //Matrix M
    Matrix M;
    zeroes(M,16);
    ubicarSubMatriz(M,0,11,0,11, sumMatrix(matrixA,matrixK,12,12));
    ubicarSubMatriz(M,0,11,12,15,matrixG);
    ubicarSubMatriz(M,12,15,0,11,matrixD);

    return M;
}

void calculateF(Vector &f, mesh &m){
    zeroes(f,3);

    f.at(0) += m.getParameter(EXTERNAL_FORCE_X);
    f.at(1) += m.getParameter(EXTERNAL_FORCE_Y);
    f.at(2) += m.getParameter(EXTERNAL_FORCE_Z);

}

Vector createLocalb(int e,mesh &m){
    float J;
    Vector b,b_aux,f;
    Matrix g_matrix;

    calculateF(f, m);

    calculateGamma(g_matrix);

    J = calculateLocalJ(e,m);

    if(J == 0){
        cout << "\n!---CATASTROPHIC FAILURE---!\n";
        exit(EXIT_FAILURE);
    }
    
    zeroes(b_aux,16);

    productMatrixVector(g_matrix,f,b_aux);

    float val = m.getParameter(HX) + m.getParameter(HY) + m.getParameter(HZ);

    b_aux.at(12) = (float) val;
    b_aux.at(13) = (float) val;
    b_aux.at(14) =(float) val;
    b_aux.at(15) =(float) val;

    productRealVector(J/120,b_aux,b);
    
    cout << b.size();
    cout << endl<<val;

    return b;
}
