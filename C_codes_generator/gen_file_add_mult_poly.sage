
def build_add_mult_poly_h_file(dir_path, small_int, big_int):
	with open(dir_path+"/add_mult_poly.h", "w") as f:
		
		f.write("#ifndef POLY_MULT_ADD\n")
		f.write("#define POLY_MULT_ADD\n\n\n")
		
		f.write("void sub_poly(" + small_int + " *rop, " + small_int + " *pa, " + small_int + " *pb);\n")
		f.write("void add_poly(" + small_int + " *rop, " + small_int + " *pa, " + small_int + " *pb);\n")
		f.write("void double_add_poly(" + small_int + " *rop, " + small_int + " *pa, " + small_int + " *pb);\n")
		f.write("void double_sub_poly(" + small_int + " *rop, " + small_int + " *pa, " + small_int + " *pb);\n")
		f.write("void neg_poly(" + small_int + " *rop, " + small_int + " *op);\n")
		f.write("void scalar_mult_poly(" + small_int + " *rop, " + small_int + " *op, " + small_int + " scalar);\n")
		f.write("void double_poly_coeffs(" + small_int + " *rop, " + small_int + " *op);\n")
		f.write("void lshift_poly_coeffs(" + small_int + " *rop, " + small_int + " *op, int nb_pos);\n\n")
		
		f.write("void mult_mod_poly(" + small_int + " *rop, " + small_int + " *pa, " + small_int + " *pb);\n\n")
		f.write("void square_mod_poly(" + small_int + " *rop, " + small_int + " *pa);\n\n")	
		
		f.write("void internal_reduction(" + small_int + " *rop, " + big_int + " *op);\n\n")
		f.write("void exact_coeffs_reduction(" + small_int + " *rop, " + small_int + " *op);\n\n")
		
		f.write("#endif\n\n")


def build_add_mult_poly_c_file(dir_path, n, mont_phi, upPow_distrib, small_int, unsigned_small_int, big_int, red_int_coeff, neg_inv_ri_rep_coeff, mask_for_redint):
	with open(dir_path+"/add_mult_poly.c", "w") as f:
		
		f.write("#include \"add_mult_poly.h\"\n\n\n")
		
		f.write("void add_poly(" + small_int + " *rop, " + small_int + " *pa, " + small_int + " *pb){\n")
		f.write("	int j;\n")
		f.write("	for (j=0; j<NB_COEFF; j++)\n")
		f.write("		rop[j] = pa[j] + pb[j];\n")
		f.write("}\n\n")

		f.write("void sub_poly(" + small_int + " *rop, " + small_int + " *pa, " + small_int + " *pb){\n")
		f.write("	int j;\n")
		f.write("	for (j=0; j<NB_COEFF; j++)\n")
		f.write("		rop[j] = pa[j] - pb[j];\n")
		f.write("}\n\n")
		
		f.write("//~ computes : pa + 2.pb\n")
		f.write("void double_add_poly(" + small_int + " *rop, " + small_int + " *pa, " + small_int + " *pb){\n")
		f.write("	int j;\n")
		f.write("	for (j=0; j<NB_COEFF; j++)\n")
		f.write("		rop[j] = pa[j] + 2*pb[j];\n")
		f.write("}\n\n")
		
		f.write("//~ computes : pa - 2.pb\n")
		f.write("void double_sub_poly(" + small_int + " *rop, " + small_int + " *pa, " + small_int + " *pb){\n")
		f.write("	int j;\n")
		f.write("	for (j=0; j<NB_COEFF; j++)\n")
		f.write("		rop[j] = pa[j] - 2*pb[j];\n")
		f.write("}\n\n")

		f.write("void neg_poly(" + small_int + " *rop, " + small_int + " *op){\n")
		f.write("	int j;\n")
		f.write("	for (j=0; j<NB_COEFF; j++)\n")
		f.write("		rop[j] = -op[j];\n")
		f.write("}\n\n")

		f.write("//~ assumes 'scalar' and/or coeffs of 'op' small enough to avoid an overflow.\n")
		f.write("void scalar_mult_poly(" + small_int + " *rop, " + small_int + " *op, " + small_int + " scalar){\n")
		f.write("	int j;\n")
		f.write("	for (j=0; j<NB_COEFF; j++)\n")
		f.write("		rop[j] = scalar * op[j];\n")
		f.write("}\n\n")
		
		f.write("//~ assumes 'scalar' and/or coeffs of 'op' small enough to avoid an overflow.\n")
		f.write("void double_poly_coeffs(" + small_int + " *rop, " + small_int + " *op){\n")
		f.write("	int j;\n")
		f.write("	for (j=0; j<NB_COEFF; j++)\n")
		f.write("		rop[j] = op[j] << 1;\n")
		f.write("}\n\n")
		
		f.write("//~ assumes 'nb_pos' and/or coeffs of 'op' small enough to avoid an overflow.\n")
		f.write("void lshift_poly_coeffs(" + small_int + " *rop, " + small_int + " *op, int nb_pos){\n")
		f.write("	int j;\n")
		f.write("	for (j=0; j<NB_COEFF; j++)\n")
		f.write("		rop[j] = op[j] << nb_pos;\n")
		f.write("}\n\n")
		
		#~ to check if it is an AMNS or not
		if (type(upPow_distrib[0]) == int) or (type(upPow_distrib[0]) == Integer) : 
			is_AMNS = True
			upPow_distrib = upPow_distrib[0]
		else :
			is_AMNS = False
		
		if is_AMNS :
			prod_funct = build_prod_code
			sqr_funct = build_square_code
		else :
			prod_funct = build_prod_code_v2
			sqr_funct = build_square_code_v2
		
		f.write("//~ Computes pa(X)*pb(X) mod(E)\n")
		f.write("void mult_mod_poly(" + small_int + " *rop, " + small_int + " *pa, " + small_int + " *pb){\n\n")
		f.write("	" + big_int + " tmp_prod_result[NB_COEFF];\n")
		f.write(prod_funct(n, upPow_distrib, big_int))
		f.write("\n	internal_reduction(rop, tmp_prod_result);\n")
		f.write("}\n\n")
		
		f.write("//~ Computes pa(X)^2 mod(E)\n")
		f.write("void square_mod_poly(" + small_int + " *rop, " + small_int + " *pa){\n\n")
		f.write("	" + big_int + " tmp_prod_result[NB_COEFF];\n")
		f.write(sqr_funct(n, upPow_distrib, big_int))
		f.write("\n	internal_reduction(rop, tmp_prod_result);\n")
		f.write("}\n\n")
		
		
		if is_AMNS :
			red_funct = build_red_int_code
			exactRedCoeff_interProd_funct = build_exactRedCoeff_interProd_code
		else:
			red_funct = build_red_int_code_v2
			exactRedCoeff_interProd_funct = build_exactRedCoeff_interProd_code_v2
		
		f.write("//~ performs the internal reduction on 'op' and puts the result in 'rop'\n")
		f.write("//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.\n")
		f.write("void internal_reduction(" + small_int + " *rop, " + big_int + " *op){\n\n")
		f.write(red_funct(unsigned_small_int, big_int, n, upPow_distrib, mont_phi, red_int_coeff, neg_inv_ri_rep_coeff, mask_for_redint))
		f.write("}\n\n")
		
		f.write("void exact_coeffs_reduction(" + small_int + " *rop, " + small_int + " *op){\n\n")
		f.write("	int i;\n")
		f.write("	" + big_int + " tmp[NB_COEFF];\n\n")
		f.write("	for(i=0; i<NB_COEFF; i++)\n")
		f.write("		tmp[i] = (" + big_int + ") op[i];\n")
		f.write("\n")
		f.write("	internal_reduction(rop, tmp);\n")
		f.write(exactRedCoeff_interProd_funct(n, upPow_distrib, big_int))
		f.write("\n")
		f.write("	internal_reduction(rop, tmp);\n")
		f.write("}\n\n")
		


