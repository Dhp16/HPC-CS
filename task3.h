


void K_effective(double* G_Mm, double* K_effective, const int Nx, const double delta_t){
	// assuming the K_effective that is input is already equal to G_Ke
	double u = 1/(0.250*delta_t);
	Matrix_by_Scalar(G_Mm, u, Nx);
	Matrix_Add_Diagonal(K_effective, G_Mm, Nx, 0);
}




void T3_inputs(){



	return;
}