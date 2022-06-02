
def build_structs_data_file(dir_path, phi_log2, poly_deg, rho_log2, nb_max_add, big_int_name, big_int, small_int, base_beta, ext_pol, ri_poly, neg_inv_ri, phi, gmm, p, F, R, PP, mask_for_redint):
	
	with open(dir_path+"/structs_data.h", "w") as f:
		
		f.write("#ifndef STRUCTS_DATA\n")
		f.write("#define STRUCTS_DATA\n\n\n")
		
		f.write("#define PHI_LOG2 " + str(phi_log2) + "\n")
		f.write("#define POLY_DEG " + str(poly_deg) + "\n")
		f.write("#define NB_COEFF " + str(poly_deg+1) + "\n")
		f.write("#define NB_ADD_MAX " + str(nb_max_add) + "\n\n")
		
		f.write("#define RHO_LOG2 " + str(rho_log2) + "  // rho = 1 << RHO_LOG2.\n\n")
		
		f.write("#define CONV_MASK " + str(base_beta - 1) + "UL  // CONV_MASK = (1 << RHO_LOG2) - 1, for conversion\n\n")
		
		if mask_for_redint != 0 :
			f.write("#define REDINT_MASK " + str(mask_for_redint) + "UL  // REDINT_MASK = (1 << PHI_LOG2) - 1, for internal reduction\n\n")
		
		f.write("typedef " + big_int_name + " " + big_int + ";\n\n")
		
		f.write("//~ representations of polynomials Pi, used for conversion into the AMNS\n")
		f.write("//~ Note: each Pi is a representation of ((1 << RHO_LOG2)^i * phi^2)\n")
		f.write(generate_conv_polys_P(small_int, base_beta, (poly_deg+1), ext_pol, ri_poly, neg_inv_ri, phi, gmm, p, F, R, PP))
		
		f.write("mpz_t modul_p;\n")
		f.write("mpz_t gama_pow[POLY_DEG];\n\n")
			
		f.write("#endif\n\n")
		
		
	


