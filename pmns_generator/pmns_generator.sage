from multiprocessing import Process, Queue
load("pmns_gen_utils.sage")

#~ ----------------------------------------------------------------------------------------------

def compute_neg_inv_ri_poly(n, phi, ri_poly, ext_poly):
	
	R.<x> = QQ[]
	P = ZZ.quo(phi)
	PP.<y> = P[]

	e = R(ext_poly) 
	m = R(ri_poly)
	imy = PP(m.inverse_mod(e))

	return (-imy)


#~ returns a representation of 'op/phi'
def pmns_red_int(op, ext_pol, ri_poly, neg_inv_ri, R, PP, phi):
	q = (PP(op)*neg_inv_ri).mod(PP(ext_pol))
	r0 = op + (R(q)*ri_poly).mod(ext_pol)
	return (r0/phi)


#~ returns a representation of '(op1*op2)/phi'
def pmns_mont_mult(op1, op2, ext_pol, ri_poly, neg_inv_ri, R, PP, phi):
	c = (op1*op2).mod(ext_pol)
	q = (PP(c)*neg_inv_ri).mod(PP(ext_pol))
	r0 = c + (R(q)*ri_poly).mod(ext_pol)
	return (r0/phi)

#~ Uses method 2 of conversion to pmns, to compute a representation of 'val', with infinity norm lower than 'rho_min'
def compute_rep_in_pmns(val, n, p, gmm, phi, ext_pol, ri_poly, neg_inv_ri):
	
	F = GF(p)
	R.<x> = ZZ[]
	P = ZZ.quo(phi); PP.<y> = P[]
	
	w = compute_redExtPol_w(n, list(ext_pol))
	rho_min = 2 * w * ri_poly.norm(infinity)
	
	tmp_rep = Integer(F(val * phi.powermod(n-1, p)))
	
	rep = pmns_red_int(R(tmp_rep), ext_pol, ri_poly, neg_inv_ri, R, PP, phi)
	for i in xrange(n-2):
		rep = pmns_red_int(rep, ext_pol, ri_poly, neg_inv_ri, R, PP, phi)
	
	if F(rep(gmm)) != F(val) :
		print("ERROR : Bad conversion !!!")
	
	if rep.norm(infinity) >= rho_min :
		print("ERROR : Element has infinity norm greater than expected !!!")
	
	return rep


#~ ----------------------------------------------------------------------------------------------

def roots_computation(ext_pol, queue):
	res = []
	try:
		res = ext_pol.roots(multiplicities=False)
	except ValueError:  # ValueError, msg:
		#~ print msg
		res = []
	finally:
		queue.put(res)

#~ tries to find a root of 'ext_pol' within 'roots_max_duration_checks' seconds
def find_roots_with_timeout(ext_pol, roots_max_duration_checks):
	
	res = []
	queue = Queue() # used to get the result
	proc = Process(target=roots_computation, args=(ext_pol, queue)) # creation of a process calling longfunction with the specified arguments
	proc.start() # lauching the processus on another thread
	try:
		res = queue.get(timeout=roots_max_duration_checks) # getting the resultat under 'max_duration' seconds or stop
		proc.join() # proper delete if the computation has take less than timeout seconds
	except Exception:  #Exception, msg:
		proc.terminate() # kill the process
		#~ print ("Timed out!")
	return res


#~ ------------------------------------------------------------------------------------------------------------------

#~ checks if 'extPol' produces some PMNS and returns them if so. 
def quick_check_ext_pol(ZX, F, p, n, redExtPol_coeffs, roots_max_duration_checks, resFileName, nb_free_add, phi_log2):
	
	T.<y> = F[]
	redExtPolT = T(redExtPol_coeffs)
	
	wE = compute_redExtPol_w(n, redExtPol_coeffs)
	
	print("Starting: " + str(ZX(redExtPol_coeffs)))
	print("E = " + str(redExtPol_coeffs))
	print("w = " + str(wE))
	
	gmms = find_roots_with_timeout(redExtPolT, roots_max_duration_checks)
	
	if gmms == []:
		print("Done: NO ROOT FOUND!\n")
		return []
	
	pmns_listt = []
	for gmm in gmms :
		
		gmm = Integer(gmm)
		
		if (gmm == 0) or (gmm == 1) or (gmm == (p-1)):
			continue # not useful roots
		
		(redIntPol_coeffs, rhoUp_log2) = find_small_valid_vect__and__compute_rhoUp_log2_and_nb_add_max(p, n, gmm, redExtPol_coeffs, nb_free_add, phi_log2)
		
		if redIntPol_coeffs != -1 :
			
			pmns_found = [nb_free_add, n, redExtPol_coeffs, rhoUp_log2, phi_log2, gmm, redIntPol_coeffs]
			
			pmns_listt.append(pmns_found)
			
			with open(resFileName, "a") as f:
				pmnsres = str(pmns_found)+"\n"
				f.write(pmnsres)
				
	print("Done: " + str(len(pmns_listt)) + " PMNS found. See corresponding results file.\n")
	
	return pmns_listt


#~ checks if 'extPol' produces some PMNS and returns them if so. 
def quick_check_ext_pol__avx(F, p, n, redExtPol_coeffs, roots_max_duration_checks, resFileName, nb_free_add, phi_log2):
	
	T.<y> = F[]
	redExtPolT = T(redExtPol_coeffs)
	
	ZW.<W> = ZZ[]
	EW = ZW(redExtPol_coeffs)
	
	wE = compute_redExtPol_w(n, redExtPol_coeffs)
	
	print("Starting: " + str(ZW(redExtPol_coeffs)))
	print("E = " + str(redExtPol_coeffs))
	print("w = " + str(wE))
	
	gmms = find_roots_with_timeout(redExtPolT, roots_max_duration_checks)
	
	if gmms == []:
		print("Done: NO ROOT FOUND!\n")
		return []
	
	pmns_listt = []
	search_fac_max = 1 << 10
	
	for gmm in gmms :
		
		gmm = Integer(gmm)
		
		if (gmm == 0) or (gmm == 1) or (gmm == (p-1)):
			continue # not useful roots
		
		(redIntPol_coeffs, rhoUp_log2) = find_all_posi_coeff_with_babai_rounding(p, n, gmm, nb_free_add, phi_log2, redExtPol_coeffs, EW, ZW, wE, nb_free_add, search_fac_max)
		
		if redIntPol_coeffs != -1 :
			
			pmns_found = [nb_free_add, n, redExtPol_coeffs, rhoUp_log2, phi_log2, gmm, redIntPol_coeffs]
			
			pmns_listt.append(pmns_found)
			
			with open(resFileName, "a") as f:
				pmnsres = str(pmns_found)+"\n"
				f.write(pmnsres)
				
	print("Done: " + str(len(pmns_listt)) + " PMNS found. See corresponding results file.\n")
	
	return pmns_listt


#~ ----------------------------------------------------------------------------------------------

def gen_extPol_set(n, lambda_max, for_avx):
	
	extPol_set = []
	
	# ~ extPol_set.append([-1]+[0]*(n-1)+[1])
	
	extPol_set.append([-1, -1]+[0]*(n-2)+[1])
	
	if not for_avx:
		
		if Integer(n).is_power_of(2):
			extPol_set.append([1]+[0]*(n-1)+[1])
		
		extPol_set.append([1, -1]+[0]*(n-2)+[1])
		extPol_set.append([-1, 1]+[0]*(n-2)+[1])
		extPol_set.append([1, 1]+[0]*(n-2)+[1])
		extPol_set.append([1]*(n+1))
		
		tpol=[1]
		coeff = -1
		for i in range(n):
			tpol = [coeff] + tpol
			coeff *= -1
		extPol_set.append(tpol)
		
		if n%2 == 0 :
			
			zero = [0]*((n//2)-1)
			extPol_set.append([1] + zero + [1] + zero + [1])
			extPol_set.append([1] + zero + [-1] + zero + [1])
			
			tpol = [1]
			for i in range(0,n,2):
				tpol = [1,0] + tpol
			extPol_set.append(tpol)
			
			tpol = [1]
			coeff = -1
			for i in range(0,n,2):
				tpol = [coeff,0] + tpol
				coeff *= -1
			extPol_set.append(tpol)
		
	for lmbd in range(2, lambda_max+1):
		extPol_set.append([-lmbd]+[0]*(n-1)+[1])
		if not for_avx:
			extPol_set.append([lmbd]+[0]*(n-1)+[1])
	
	
	extPol_set.sort(key=extPolys_sorting_criteria)
	
	return extPol_set


def extPolys_sorting_criteria(ext_pol_coeffs):
		
	return compute_redExtPol_w(len(ext_pol_coeffs) - 1, ext_pol_coeffs)

#~ ----------------------------------------------------------------------------------------------

def build_pmns_candidates_for_n(p, n, phi_log2, nb_free_add, lambda_max, roots_max_duration_checks, for_avx):
	
	F = GF(p)
	
	extPol_set = gen_extPol_set(n, lambda_max, for_avx)
	
	print("STARTING") 
	print("Candidate external polynomials generation done, nbre cdts: " + str(len(extPol_set)) + "\n") 
	
	resFileName = "results/p" + str(p.nbits()) + "_n" + str(n) + "_phi" + str(phi_log2) + "_delta" + str(nb_free_add)
	if for_avx:
		resFileName += "__for_avx"
	
	with open(resFileName, "w") as f:
		pval = "p = " + str(p)+"\n\n"
		phi_size = "phi_log2 = " + str(phi_log2)+"\n\n"
		dt_struc = "data_struct = [nb_free_add, n, redExtPol_coeffs, rhoUp_log2, phi_log2, gmm, redIntPol_coeffs]\n\n"
		f.write(pval)
		f.write(phi_size)
		f.write(dt_struc)
	
	pmns_cands = []	
	
	if not for_avx:
		ZX.<X> = QQ[]
		for extPol in extPol_set:
			rep = quick_check_ext_pol(ZX, F, p, n, extPol, roots_max_duration_checks, resFileName, nb_free_add, phi_log2)
			if rep != [] :
				pmns_cands.append(rep)
	else:
		for extPol in extPol_set:
			rep = quick_check_ext_pol__avx(F, p, n, extPol, roots_max_duration_checks, resFileName, nb_free_add, phi_log2)
			if rep != [] :
				pmns_cands.append(rep)
	
	pmns_cands = flatten(pmns_cands, max_level=1)
	
	print("DONE. Total number of PMNS found: " + str(len(pmns_cands)))
		
	return pmns_cands



#~ ----------------------------------------------------------------------------------------------





















