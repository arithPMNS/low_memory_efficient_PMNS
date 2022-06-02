def build_useful_functs_h_file(dir_path, small_int, unsigned_small_int, big_int):
	with open(dir_path+"/useful_functs.h", "w") as f:
		
		f.write("#ifndef USEFUL_FUNCTS\n")
		f.write("#define USEFUL_FUNCTS\n\n\n")
		
		f.write("void init_data();\n\n")
		f.write("void free_data();\n\n")
		f.write("void from_int_to_amns(" + small_int + " *rop, mpz_t op);\n\n")
		f.write("void from_amns_to_int(mpz_t rop, " + small_int + " *op);\n\n")
		f.write("int cmp_poly_evals(" + small_int + " *pa, " + small_int + " *pb);\n\n")
		f.write("void copy_poly(" + small_int + " *rop, " + small_int + " *op);\n\n")
		f.write("void add_lpoly(" + big_int + " *rop, " + big_int + " *pa, " + big_int + " *pb);\n\n")
		f.write("void scalar_mult_lpoly(" + big_int + " *rop, " + small_int + " *op, " + unsigned_small_int + " scalar);\n\n")
		f.write("void from_mont_domain(" + small_int + " *rop, " + small_int + " *op);\n\n")
		f.write("void print_element(" + small_int + " *poly);\n\n")
		
		f.write("#endif\n\n")


def build_useful_functs_c_file(dir_path, small_int, unsigned_small_int, big_int, p, gmm):
	with open(dir_path+"/useful_functs.c", "w") as f:
		
		f.write("#include \"useful_functs.h\"\n\n\n")
		
		f.write("//~ Assumes allocation already done for 'rop'.\n")
		f.write("//~ IMPORTANT : convertion to montgomery domain will be done here\n")
		f.write("void from_int_to_amns(" + small_int + " *rop, mpz_t op){\n")
		f.write("	int i;\n")
		f.write("	mpz_t tmp;\n")
		f.write("	" + big_int + " tmp_poly[NB_COEFF];\n")
		f.write("	" + big_int + " tmp_sum[NB_COEFF];\n\n")
		f.write("	mpz_init_set(tmp, op);\n\n")
		f.write("	for(i=0; i<NB_COEFF; i++){\n")
		f.write("		rop[i] = 0;\n")
		f.write("		tmp_sum[i] = 0;\n")
		f.write("	}\n\n")
		f.write("	if(tmp->_mp_size == 0)\n")
		f.write("		return;\n\n")
		f.write("	i = 0;\n")
		f.write("	while(tmp->_mp_size && (i < NB_COEFF)){\n")
		f.write("		scalar_mult_lpoly(tmp_poly, polys_P[i++], (tmp->_mp_d[0]) & CONV_MASK);\n")
		f.write("		add_lpoly(tmp_sum, tmp_sum, tmp_poly);\n\n")
		f.write("		mpz_tdiv_q_2exp (tmp, tmp, RHO_LOG2);\n")
		f.write("	}\n\n")
		f.write("	internal_reduction(rop, tmp_sum);\n\n")
		f.write("	mpz_clear(tmp);\n")
		f.write("}\n\n")
		
		f.write("//~ Assumes \"rop\" already initialized.\n")
		f.write("//~ IMPORTANT : convertion from montgomery domain will be done here.\n")
		f.write("void from_amns_to_int(mpz_t rop, " + small_int + " *op){\n")
		f.write("	int i;\n")
		f.write("	mpz_t tmp_sum;\n")
		f.write("	" + small_int + " tmp_conv[NB_COEFF];\n\n")
		f.write("	mpz_init(tmp_sum);\n\n")
		f.write("	//~ convertion out of mont domain\n")
		f.write("	from_mont_domain(tmp_conv, op);\n\n")
		f.write("	mpz_set_si(rop, tmp_conv[0]);\n")
		f.write("	for(i=0; i<POLY_DEG; i++){\n")
		f.write("		mpz_mul_si(tmp_sum, gama_pow[i], tmp_conv[i+1]);\n")
		f.write("		mpz_add(rop, rop, tmp_sum);\n")
		f.write("	}\n")
		f.write("	mpz_mod (rop, rop, modul_p);\n\n")
		f.write("	mpz_clear(tmp_sum);\n")
		f.write("}\n\n")
		
		f.write("//~ computes : op/phi\n")
		f.write("void from_mont_domain(" + small_int + " *rop, " + small_int + " *op){\n\n")
		f.write("	int i;\n")
		f.write("	" + big_int + " tmp[NB_COEFF];\n\n")
		f.write("	for(i=0; i<NB_COEFF; i++)\n")
		f.write("		tmp[i] = (" + big_int + ") op[i];\n")
		f.write("\n")
		f.write("	internal_reduction(rop, tmp);\n")
		f.write("}\n\n")
		
		f.write("void init_data(){\n\n")
		f.write("	int i;\n")
		f.write("	for(i=0; i<POLY_DEG; i++)\n")
		f.write("		mpz_init (gama_pow[i]);\n\n")
		f.write("	mpz_init (modul_p);\n\n\n")
		f.write("	mpz_set_str (modul_p, \"" + str(p) + "\", 10);\n\n")
		f.write("	mpz_set_str (gama_pow[0], \"" + str(gmm) + "\", 10);\n")
		f.write("	for(i=1; i<POLY_DEG; i++){\n")
		f.write("		mpz_mul (gama_pow[i], gama_pow[i-1], gama_pow[0]);\n")
		f.write("		mpz_mod (gama_pow[i], gama_pow[i], modul_p);\n")
		f.write("	}\n")
		f.write("}\n\n")
		
		f.write("void free_data(){\n")
		f.write("	int i;\n")
		f.write("	for(i=0; i<POLY_DEG; i++)\n")
		f.write("		mpz_clear (gama_pow[i]);\n")
		f.write("\n")
		f.write("	mpz_clear (modul_p);\n")
		f.write("}\n\n")
		
		f.write("//~ ----------------------------------------------------------------------------------------\n\n")
		f.write("//~ return a positive value if pa > pb, zero if pa = pb, or a negative value if pa < pb.\n")
		f.write("//~ Important : evaluation is done using the corresponding integers modulo 'p'.\n")
		f.write("int cmp_poly_evals(" + small_int + " *pa, " + small_int + " *pb){\n")
		f.write("	int rep;\n")
		f.write("	mpz_t a, b;\n")
		f.write("	mpz_inits (a, b, NULL);\n")
		f.write("	from_amns_to_int(a, pa);\n")
		f.write("	from_amns_to_int(b, pb);\n")
		f.write("	rep = mpz_cmp (a, b);\n")
		f.write("	mpz_clears (a, b, NULL);\n")
		f.write("	return rep;\n")
		f.write("}\n\n")
		
		f.write("void copy_poly(" + small_int + " *rop, " + small_int + " *op){\n")
		f.write("	int i;\n")
		f.write("	for(i=0; i<NB_COEFF; i++)\n")
		f.write("		rop[i] = op[i];\n")
		f.write("}\n\n")

		f.write("void add_lpoly(" + big_int + " *rop, " + big_int + " *pa, " + big_int + " *pb){\n")
		f.write("	int j;\n")
		f.write("	for (j=0; j<NB_COEFF; j++)\n")
		f.write("		rop[j] = pa[j] + pb[j];\n")
		f.write("}\n\n")

		f.write("//~ assumes 'scalar' and/or coeffs of 'op' small enough to avoid an overflow.\n")
		f.write("void scalar_mult_lpoly(" + big_int + " *rop, " + small_int + " *op, " + unsigned_small_int + " scalar){\n")
		f.write("	int j;\n")
		f.write("	for (j=0; j<NB_COEFF; j++)\n")
		f.write("		rop[j] = (" + big_int + ")op[j] * scalar;\n")
		f.write("}\n\n")
		
		f.write("void print_element(" + small_int + " *poly){\n")
		f.write("	int i;\n")
		f.write("	printf(\"[\");\n")
		f.write("	for (i=0; i<POLY_DEG; i++)\n")
		f.write("		printf(\"%2ld, \", poly[i]);\n")
		f.write("	printf(\"%2ld]\", poly[POLY_DEG]);\n")
		f.write("}\n\n")

















