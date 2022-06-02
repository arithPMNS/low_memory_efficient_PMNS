from os import makedirs as create_dir
from shutil import rmtree as delete_code
from shutil import copy as cpfile
load("c_codes_gen_utils.sage")
load("build_avx_data.sage")
load("gen_file_structs_data.sage")
load("gen_file_add_mult_poly.sage")
load("gen_file_useful_functs.sage")



#~ assumes that polynomial representation form is P(X) = a_0 + ... + a_n.X^n = [a_0, ..., a_n].
#~ assumes 'amns_data' is a list containning (in this order): [nb_max_add, n, redExtPol_coeffs, rhoUp_log2, phi_log2, gmm, M]
#~ assumes 'target_archi_info' is a list containning (in this order): [small_int_name, unsigned_small_int_name, big_int_name, big_int_new_name, target_word_size]
#~ WARNING : 'amns_data' is supposed generated for 'target_archi_info'
def build_amns_c_codes(p, amns_data, target_archi_info, for_avx, p_num=0, amns_num=0):
	
	amnsData = list(amns_data)
	
	redExtPol_coeffs = amnsData[2]
	
	
	if ((len(redExtPol_coeffs) - redExtPol_coeffs.count(0)) == 2) and (redExtPol_coeffs[0] != 0) : # E(X) = X^n - lambda
		
		amnsData.append([-redExtPol_coeffs[0]])
		
		is_amns = True
		
	else :
		
		is_amns = False
		
		R.<x> = ZZ[]
		x = R.gen()
		re_pol = R(redExtPol_coeffs)
		
		n = re_pol.degree()
		
		tmp = x**n
		hcoeffs = []
		for d in range(n - 1):
			ll = list(tmp%re_pol)
			ll += [0] * (n - len(ll))
			hcoeffs.append(ll)
			tmp *= x

		amnsData.append(matrix(hcoeffs))
		
	gen_c_codes(p, amnsData, target_archi_info, p_num, amns_num, for_avx, is_amns)
	
	return;


def compute_neg_inv_ri_poly(n, PP, ri_polyC, ext_polyC):
	
	R.<X> = QQ[]
	
	imy = PP(R(ri_polyC).inverse_mod(R(ext_polyC)))
	
	neg_inv_riP = -imy
	
	neg_inv_riC = []
	for el in list(neg_inv_riP):
		neg_inv_riC.append(Integer(el))

	return (neg_inv_riC, neg_inv_riP)


def gen_c_codes(p, amns_data, target_archi_info, p_num, amns_num, for_avx, is_amns):
	
	small_int = target_archi_info[0]
	unsigned_small_int = target_archi_info[1]
	big_int_name = target_archi_info[2]
	big_int = target_archi_info[3]
	target_word_size = target_archi_info[4]
	
	nb_max_add = amns_data[0]
	n = amns_data[1]
	redExtPol = amns_data[2]
	rho_log2 = amns_data[3]
	phi_log2 = amns_data[4]
	gmm = amns_data[5]
	red_int_pol = amns_data[6]
	upPow_distrib = amns_data[7]
	
	F = GF(p)
	R.<x> = ZZ[]
	mont_phi = 1 << phi_log2
	P = ZZ.quo(mont_phi); PP.<y> = P[]
	
	base_beta = 1 << rho_log2
	ext_pol = R(redExtPol)
	ri_poly = R(red_int_pol)
	
	neg_inv_ri, neg_inv_riP = compute_neg_inv_ri_poly(n, PP, red_int_pol, redExtPol)
	
	mask_for_redint = ((1 << phi_log2) - 1) if target_word_size != phi_log2 else 0
	
	dir_name = 'p' + str(p.nbits()) + '_n'  + str(n) + '_rho'  + str(rho_log2) + '_phi'  + str(phi_log2) + '_d'  + str(nb_max_add)
	if p_num != 0:
		dir_name += '_' + str(p_num)
	if amns_num != 0:
		dir_name += '_' + str(amns_num)
	
	dir_path = "c_codes/" + dir_name
	if is_amns:
		dir_path += "__amns"
	else:
		dir_path += "__pmns"
	
	if for_avx :
		dir_path += "_avx"
		
	try:
		create_dir(dir_path)
	except OSError:  # if this directory already exist
		delete_code(dir_path)
		create_dir(dir_path)
	
	cpfile("constant_files/main__for_nb_cycles.c", dir_path+"/main.c")
	cpfile("constant_files/main__check_ops.c", dir_path+"/main__check_ops.c")
		
	
	cpfile("constant_files/gmp_stuff.c", dir_path)
	cpfile("constant_files/intel_measurement_stuff.c", dir_path)
	
	
	build_structs_data_file(dir_path, phi_log2, (n-1), rho_log2, nb_max_add, big_int_name, big_int, small_int, base_beta, ext_pol, ri_poly, neg_inv_riP, mont_phi, gmm, p, F, R, PP, mask_for_redint)

	build_add_mult_poly_h_file(dir_path, small_int, big_int)
	build_add_mult_poly_c_file(dir_path, n, mont_phi, upPow_distrib, small_int, unsigned_small_int, big_int, red_int_pol, neg_inv_ri, mask_for_redint)

	build_useful_functs_h_file(dir_path, small_int, unsigned_small_int, big_int)
	build_useful_functs_c_file(dir_path, small_int, unsigned_small_int, big_int, p, gmm)
	
	if for_avx :
		ri_mat, iri_mat = compute_redintmatrices(n, mont_phi, red_int_pol, neg_inv_ri, redExtPol)
		build_avx_const_data(dir_path, n, ri_mat, iri_mat)
		#return [ri_mat, iri_mat]
	
	return;


def compute_redintmatrices(n, phi, M_coeffs, invM_coeffs, E_coeffs):
	R.<x> = ZZ[]
	P = ZZ.quo(phi); PP.<y> = P[]
	X = R.gen()
	Y = PP.gen()
	EX = R(E_coeffs)
	EY = PP(E_coeffs)
	tmp = R(M_coeffs)
	tmp2 = PP(invM_coeffs)
	
	lx = M_coeffs + [0]*(n - len(M_coeffs))
	ly = invM_coeffs + [0]*(n - len(invM_coeffs))
	
	ri_mat = [lx]
	iri_mat = [ly]
	for i in range(n-1):
		tmp = (tmp*X)%EX
		tmp2 = (tmp2*Y)%EY
		
		l1 = list(tmp)
		l2 = list(tmp2)
		
		lx = l1 + [0]*(n - len(l1))
		ly = l2 + [0]*(n - len(l2))
		
		ri_mat.append(lx)
		iri_mat.append(ly)
	
	# ~ return [matrix(ri_mat).transpose(), matrix(iri_mat).transpose()]
	return [matrix(ri_mat), matrix(iri_mat)]













