#~ representation is of the form : a0 + a1.X + ... + a_{n-1}.X^{n-1}
# important: assumes val < base**n
def compute_rep_in_base(val, base, n):
	val = Integer(val)
	rep = []
	while val != 0 :
		rep.append(val % base)
		val //= base
	return rep + [0]*(n-len(rep))

#~ assumes that little-endian representation is used 
def compute_val_from_rep(rep, base):
	l = len(rep)
	if l==0 :
		return 0
	val = rep[0]
	for i in range(1, l):
		val += rep[i]*(base**i)
	return val

#~ ----------------------------------------------------------------------------------------------

#~ returns the infinity norm of 'vect'
def get_inf_norm(vect):
	tmp = [abs(a) for a in vect]
	return max(tmp)

#~ returns an element of 'm' with the smallest infinity norm
#~ important: assumes 'm' is not empty
def get_smallest_vect(m): 
	sml_vect = m[0]
	sv_norm = get_inf_norm(m[0])
	for v in m[1:] :
		tmp = get_inf_norm(v)
		if sv_norm > tmp :
			sv_norm = tmp
			sml_vect = v
	return list(sml_vect)

#~ returns an element of 'm' with the biggest infinity norm
#~ important: assumes 'm' is not empty
def get_biggest_vect(m): 
	bg_vect = m[0]
	bv_norm = get_inf_norm(m[0])
	for v in m[1:] :
		tmp = get_inf_norm(v)
		if bv_norm < tmp :
			bv_norm = tmp
			bg_vect = v
	return list(bg_vect)

#~ ------------------------------------------------------------------------------------------

def compute_red_int_matrix_norm1(redIntPol_coeffs, redExtPol_coeffs, n):
	
	R.<x> = ZZ[]
	X = R.gen()
	
	ext_pol = R(redExtPol_coeffs)
	ri_poly = R(redIntPol_coeffs)
	
	bs = [redIntPol_coeffs]
	
	XiM = ri_poly
	for i in range(1,n):
		XiM = (X * XiM)%ext_pol
		ll = list(XiM)
		bs.append(ll + [0]*(n-len(ll)))
	
	mat = matrix(bs)
	
	return mat.norm(1)


def compute_redExtPol_w(n, redExtPol_coeffs):
	
	R.<x> = ZZ[]
	x = R.gen()
	re_pol = R(redExtPol_coeffs)
	
	c = n-1
	tmp = x**n
	l1 = list(tmp%re_pol)
	l1 = [abs(k) for k in l1]
	V = c * vector(l1 + [0]*(n-len(l1)))
	for d in range(n-2):
		c -= 1
		tmp *= x
		l1 = list(tmp%re_pol)
		l1 = [abs(k) for k in l1]
		V += c * vector(l1 + [0]*(n-len(l1)))
	
	V += vector(range(1, n+1))
	
	return max(V)


#~ ----------------------------------------------------------------------------------------------
#~ Note: we use the equivalence: a0 + a1.X + ... + a_{n-1}.X^{n-1} <==> (a0, a1, ..., a_{n-1})

def build_lattice_base(p, n, gmm):
	b = []
	l = [p] + [0]*(n-1)
	b.append(l)
	m = identity_matrix(n)
	for i in range(1,n):
		t = (-gmm.powermod(i, p))%p
		if t%2 == 1 :
			t += p
		l = [t] + list(m[i][1:])
		b.append(l)
	bb = matrix(b)
	return bb.LLL(delta=1, algorithm='NTL:LLL')


#~ assumes : k lower than 2^n
def build_sumVect(k, n):
	rep = []
	for i in range(n):
		tmp = k & 1
		k >>= 1 
		rep.append(tmp)
	return rep


#~ note: we assume that elements coeffs can be negatives in the PMNS, so one bit will be allocated for that.
def compute_rhoUp_log2(n, redExtPol_coeffs, redIntPol_coeffs, nb_free_add, phi_log2):
	
	phi = 1 << phi_log2
	phi2 = phi**2
	
	w = compute_redExtPol_w(n, redExtPol_coeffs)
	
	redInt_mat_norm1 = compute_red_int_matrix_norm1(redIntPol_coeffs, redExtPol_coeffs, n)
		
	rho_min = 2 * redInt_mat_norm1
	
	rhoUp_log2 = ceil(log(rho_min, 2))
	rho = 1 << rhoUp_log2 
	
	tmp_phi = 2 * w * rho
	tmp_phi *= ((nb_free_add + 1)**2)  # for 'free' adds
	if tmp_phi > phi :
		return -1
		
	tmp_prod1 = w * ((nb_free_add + 1) * (rho - 1))**2
	tmp_prod2 = w * (phi -1) * get_inf_norm(redIntPol_coeffs)
	prod_acc_max = tmp_prod1 + tmp_prod2
	if prod_acc_max >= phi2:
		return -1
	
	return rhoUp_log2

#~ ------------------------------------------------------------------------------------------

#~ returns a small valid vect if any
# ~ def find_small_valid_vect(n, tE, lattice_base, T):
def find_small_valid_vect(rrbase, lb, tE, T):
	
	lattice_base = matrix(rrbase)
	
	vInit = 0*lattice_base[0] # to obtain the good type of vector later
	
	vmax = 2**lb 
	vld_list = []
	for k in range(1, vmax): 
		
		vsum = vInit
		k_vt = build_sumVect(k, lb)
		
		for i in range(lb):
			vsum += k_vt[i]*lattice_base[i]
		
		M = list(vsum)
		tM = T(M)
		
		if (tE.resultant(tM) % 2) == 1:
			vld_list.append(M)
	
	if vld_list != [] :
		return get_smallest_vect(vld_list)
	return []


def exhaustive_search(rrbase, lb, n, redExtPol_coeffs, nb_free_add, phi_log2, tE, T):
	
	M = find_small_valid_vect(rrbase, lb, tE, T)
	
	if M == []:
		return (-1, -1)
	
	rhoUp_log2 = compute_rhoUp_log2(n, redExtPol_coeffs, M, nb_free_add, phi_log2)
			
	if rhoUp_log2 != -1 :
		return (M, rhoUp_log2)
	
	return (-1, -1)


def check_base_rows(rbase, n, redExtPol_coeffs, nb_free_add, phi_log2, tE, T):
	
	vld_zero = (0,0)
	
	for vt in rbase:
		
		M = list(vt)
		tM = T(M)
		
		if (tE.resultant(tM) % 2) == 1:
			rhoUp_log2 = compute_rhoUp_log2(n, redExtPol_coeffs, M, nb_free_add, phi_log2)
			
			if rhoUp_log2 != -1 :
				if (vld_zero == (0,0)) or (rhoUp_log2 < vld_zero[1]): # note: the order is important here
					vld_zero = (M, rhoUp_log2)
	
	return vld_zero


def eliminate_too_big_elemts(rbase, phi_log2):
	res = []
	phi = 1 << phi_log2
	for el in rbase:
		infn = get_inf_norm(el)
		if get_inf_norm(el) < phi:
			res.append(list(el))
	return res


#~ ------------------------------------------------------------------------------------------

#~ Find a small valid vect if any
#~ note: we assume that elements coeffs can be negatives in the PMNS, so one bit will be allocated for that.
#~ NOTE : rho is taken as a power of two
def find_small_valid_vect__and__compute_rhoUp_log2_and_nb_add_max(p, n, gmm, redExtPol_coeffs, nb_free_add, phi_log2):
	
	T.<y> = ZZ[]
	tE = T(redExtPol_coeffs)
	coeff_size_max = phi_log2
	
	rbase = build_lattice_base(p, n, gmm)
	rrbase = eliminate_too_big_elemts(rbase, phi_log2)
	
	lb = len(rrbase)
	if lb == 0:
		return (-1, -1)
		
	exhaustive_search_bound = 10
	if lb <= exhaustive_search_bound: # 'lb' is small enough for an exhaustive search
		return exhaustive_search(rrbase, lb, n, redExtPol_coeffs, nb_free_add, phi_log2, tE, T)
	
	
	(M, rhoUp_log2) = check_base_rows(rrbase, n, redExtPol_coeffs, nb_free_add, phi_log2, tE, T)
	if rhoUp_log2 != 0:
		return (M, rhoUp_log2)
	elif len(rrbase) < 2:
		return (-1, -1)


	print("a random search here ... ")
		
	safe_trans_limit = lb - exhaustive_search_bound
	rand_lshift_pos = randint(0, safe_trans_limit)
	rrbase = rrbase[rand_lshift_pos:(rand_lshift_pos+exhaustive_search_bound)]
	
	vmax = 1 << exhaustive_search_bound
	
	vInit = 0*rbase[0] # to obtain the good type of vector in each iteration
	
	for k in range(3, vmax): 
		
		if (k & (k-1)) == 0: #k.is_power_of(2)
			#already done (with 'check_base_rows' above)
			#print(bin(k)[2:])
			continue
		
		(M, rhoUp_log2) = check_k(k, exhaustive_search_bound, n, vInit, rrbase, redExtPol_coeffs, nb_free_add, phi_log2, tE, T)
		
		if rhoUp_log2 != 0:
			return (M, rhoUp_log2)
		
	return (-1, -1)

#~ ---------------------------------------------------------------------------

def check_k(k, l_rbase, n, vInit, rbase, redExtPol_coeffs, nb_free_add, phi_log2, tE, T):
	
	vsum = vInit
	
	k_vt = build_sumVect(k, l_rbase)
	
	for i in range(l_rbase):
		vsum += k_vt[i]*rbase[i]
	
	M = list(vsum)
	tM = T(M)
	
	if (tE.resultant(tM) % 2) == 1:
		
		rhoUp_log2 = compute_rhoUp_log2(n, redExtPol_coeffs, M, nb_free_add, phi_log2)
		
		if rhoUp_log2 != -1 :
			return (M, rhoUp_log2)
	
	return (0, 0)

#~ ---------------------------------------------------------------------------
#~ ---------------------------------------------------------------------------

def is_base_for_avx(rbase, max_val_log2):
	
	max_val = 1 << max_val_log2
	
	for el in rbase:
		if get_inf_norm(el) > max_val:
			return False
	
	return True


def find_all_posi_coeff_with_babai_rounding(p, n, gmm, dlt, phi_log2, redExtPol_coeffs, ET, ZT, wE, nb_free_add, fac_max):
	
	B = build_lattice_base(p, n, gmm)
	
	if not is_base_for_avx(B, phi_log2):
		return 0
	
	iB = B.inverse()
	
	phi = 1 << phi_log2
	
	k = 2 * wE * (dlt + 1) * (dlt + 1)
	k2 = wE * wE * (dlt + 1) * (dlt + 1)
	
	vect_ubound = floor(phi/k2)
		
	for i in range(1,fac_max):
		(M, rhoUp_log2) = check_vect(1/i, vect_ubound, n, B, iB, redExtPol_coeffs, ET, ZT, wE, k, phi_log2, nb_free_add)
		if rhoUp_log2 != 0:
			return (M, rhoUp_log2)
	
	for i in range(2,fac_max):
		(M, rhoUp_log2) = check_vect(i, vect_ubound, n, B, iB, redExtPol_coeffs, ET, ZT, wE, k, phi_log2, nb_free_add)
		if rhoUp_log2 != 0:
			return (M, rhoUp_log2)
	
	return (-1, -1)


def check_vect(fact, vect_ubound, n, B, iB, redExtPol_coeffs, ET, ZT, wE, k, phi_log2, nb_free_add):
	
	vect = vector([fact*vect_ubound]*n)
	
	# ~ #or
	# ~ vect = []
	# ~ lint = floor(vect_ubound/2)
	# ~ rint = vect_ubound
	# ~ for i in range(n):
		# ~ vect.append(fact*randint(lint, rint))
	# ~ vect = vector(vect)
	
	M = check_vect_step(vect, n, B, iB, ET, ZT)

	if M == 0:
		return (0, 0)
	
	rhoUp_log2 = compute_rhoUp_log2(n, redExtPol_coeffs, M, nb_free_add, phi_log2)
	if rhoUp_log2 != -1 :
		return (M, rhoUp_log2)
	else:
		return (0, 0)


# ~ looking for valid poly M, with positive coeffs
def check_vect_step(vect, n, B, iB, ET, ZT):
	
	rvect = vect*iB
	
	rM = []
	for j in range(n):
		if rvect[j] < 0:
			rM.append(ceil(rvect[j]))
		else:
			rM.append(floor(rvect[j]))
	# ~ for j in range(n):
		# ~ rM.append(round(rvect[j]))
	
	M = list(vector(rM)*B)
	minM, maxM = min(M), max(M)
	
	check = (minM >= 0) and (minM < maxM) and (ET.resultant(ZT(M)) % 2 == 1)
	
	if check:
		return M
	else:
		return 0


#~ ---------------------------------------------------------------------------



