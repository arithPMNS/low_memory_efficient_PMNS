def get_sign(k):
	if k < 0 :
		return '-'
	else :
		return '+'

def is_power_of_two(k) :
	if k <= 0 :
		return False
	return ((k & (k - 1)) == 0)

def same_sign(id_list):
	s = id_list[0].sign()
	for el in id_list[1:]:
		if el.sign() != s :
			return [False, 1]
	return [True, s]

def build_ith_mat(pol_coeffs):
	l = len(pol_coeffs)
	ll = [0]*l
	mat = []
	for i in range(1, l+1):
		ll = [pol_coeffs[l-i]] + ll[:-1]
		mat.append(ll)
	return matrix(mat)

def build_itl_mat(pol_coeffs):
	n = len(pol_coeffs)
	ll = list(pol_coeffs)
	mat = [ll]
	for i in range(1, n):
		ll = [0] + ll[:-1]
		mat.append(ll)
	return matrix(mat)

def build_lambda_mat(lambd, l):
	ll = [lambd] + [0]*l
	mat = [ll]
	for i in range(1, l):
		ll = [0] + ll[:-1]
		mat.append(ll)
	return matrix(mat)

def build_cell_prod(cell_dict):
	
	ch = ""
	for k in cell_dict.keys() : # elements are > 0
		
		id_list = cell_dict[k]
		
		if len(id_list) == 1 :
			
			idd = id_list[0]
			
			if k == 1 :
				tch = ""
			else :
				tch = str(k) + "*"
			
			ch += " " + get_sign(idd) + " " + tch + "c" + str(abs(idd))
		
		elif k == 1 :
			for idd in id_list :
				ch += " " + get_sign(idd) + " c" + str(abs(idd))
		else :
		
			[vr, s] = same_sign(id_list)
			
			if vr :
				ch += " " + get_sign(s) + " " + str(k) + "*("
				
				idd = s*id_list[0] # > 0
				ch += "c" + str(idd)
				for idd in id_list[1:] :
					ch += " + c" + str(s*idd)
				
			else :
				ch += " + " + str(k) + "*("
				
				idd = id_list[0] 
				if idd > 0 :
					ch += "c" + str(idd)
				else :
					ch += "-c" + str(-idd)
				for idd in id_list[1:] :
					ch += " " + get_sign(idd) + " c" + str(abs(idd))
		
			ch += ")"
	
	return ch

def build_cell_tmpZero(cell_dict, big_int):
	
	ch = ""
	for k in cell_dict.keys() : # keys are > 0
		
		id_list = cell_dict[k]
		
		if len(id_list) == 1 :
			
			idd = id_list[0]
			
			if k == 1 :
				ch += " " + get_sign(idd) + " (" + big_int + ")tmpQ[" + str(abs(idd)-1) + "]"
			else :
				ch += " + (" + big_int + ")tmpQ[" + str(abs(idd)-1) + "] * " + str(hex(k*idd.sign())) + "L"
			
		elif k == 1 :
			for idd in id_list :
				ch += " " + get_sign(idd) + " (" + big_int + ")tmpQ[" + str(abs(idd)-1)+ "]"
		else :
		
			[vr, s] = same_sign(id_list)
			
			if vr :
				ch += " + ((" + big_int + ")tmpQ[" + str(s*id_list[0] - 1)+ "]"
				for idd in id_list[1:] :
					ch += " + (" + big_int + ")tmpQ[" + str(s*idd - 1)+ "]"
				k *= s 
			else :
				ch += " + ("
				
				idd = id_list[0] 
				if idd > 0 :
					ch += "(" + big_int + ")tmpQ[" + str(idd - 1)+ "]"
				else :
					ch += "-(" + big_int + ")tmpQ[" + str(-idd - 1)+ "]"
				for idd in id_list[1:] :
					ch += " " + get_sign(idd) + " (" + big_int + ")tmpQ[" + str(abs(idd) - 1)+ "]"
		
			ch += ") * " + str(hex(k)) + "L"
	
	return ch


def build_cell_tmpQ(cell_dict, unsigned_small_int):
	
	ch = ""
	for k in cell_dict.keys() : # keys are > 0
		
		id_list = cell_dict[k] # elements are > 0
		
		if len(id_list) == 1 :
			if k == 1 :
				ch += " + (" + unsigned_small_int + ")op[" + str(id_list[0]-1) + "]"
			else :
				ch += " + (" + unsigned_small_int + ")op[" + str(id_list[0]-1) + "] * " + str(hex(k)) + "UL"
		elif k == 1 :
			for idd in id_list :
				ch += " + (" + unsigned_small_int + ")op[" + str(idd-1) + "]"
		else :
			ch += " + ((" + unsigned_small_int + ")op[" + str(id_list[0]-1) + "]"
			for idd in id_list[1:] :
				ch += " + (" + unsigned_small_int + ")op[" + str(idd-1) + "]"
			ch += ") * " + str(hex(k)) + "UL"
	
	return ch

#~ note : jump > 0
def build_scal_prod(mat, jump, call_for, big_int=None, unsigned_small_int=None):
	
	(n, m) = mat.dimensions()
	
	tmp = []
	for i in range(m):
		
		tmp_dict = {}
		
		for j in range(n):
			k = abs(mat[j][i])
			if k == 0 :
				continue
			s = mat[j][i].sign() 
			try:
				tmp_dict[k].append(s*(j+jump))
			except KeyError :
				tmp_dict[k] = [s*(j+jump)]
		
		if tmp_dict == {} :
			tmp.append("")
		else :
			if call_for == "tmpQ" :
				tmp.append(build_cell_tmpQ(tmp_dict, unsigned_small_int))
			elif call_for == "tmpZero" :
				tmp.append(build_cell_tmpZero(tmp_dict, big_int))
			else :
				tmp.append(build_cell_prod(tmp_dict))
	
	return tmp


#~ --------------------------------------- for E(X) = X^n - lambda ------------------------------------------------------


def build_prod_code(n, lambd, big_int):
	
	c_sign = get_sign(lambd)
	abs_c = abs(lambd)
	abs_c_is_2pow = is_power_of_two(abs_c)
	
	tmp_part1 = ['']*n
	tmp_part2 = ['']*(n-1)
	
	for i in range(n):
		tmp_part1[i] = 'tmp_prod_result['+str(i)+'] ='
	
	for i in range(n):
		for j in range(n):
			pos = i+j
			if pos < n :
				tmp_part1[pos]+= ' (' + big_int + ')pa['+ str(i) + '] * pb['+ str(j) + '] +'
			else:
				tmp_part2[pos%n]+= ' (' + big_int + ')pa['+ str(i) + '] * pb['+ str(j) + '] +'
	
	for i in range(n-1):
		tmp_part1[i] = tmp_part1[i][:-2]
		if abs_c == 1 :
			tmp_part2[i] = ' ' + c_sign + ' (' + tmp_part2[i][1:-2] +')'
		elif abs_c_is_2pow :
			tmp_part2[i] = ' ' + c_sign + ' ((' + tmp_part2[i][1:-2] +') << ' + str(int(log(abs_c, 2))) + ')'
		else :
			tmp_part2[i] = ' ' + c_sign + ' ((' + tmp_part2[i][1:-2] +') * ' + str(abs_c) + ')'
	
	tmp_part1[n-1] = tmp_part1[n-1][:-2]
	
	result = '\n'
	for i in range(n-1):
		result += "	" + tmp_part1[i] + tmp_part2[i] + ';\n'
	result +=  "	" + tmp_part1[n-1] + ';\n'
	
	return result


def build_exactRedCoeff_interProd_code(n, lambd, big_int):
	
	c_sign = get_sign(lambd)
	abs_c = abs(lambd)
	abs_c_is_2pow = is_power_of_two(abs_c)
	
	tmp_part1 = ['']*n
	tmp_part2 = ['']*(n-1)
	
	for i in range(n):
		tmp_part1[i] = 'tmp['+str(i)+'] ='
	
	for i in range(n):
		for j in range(n):
			pos = i+j
			if pos < n :
				tmp_part1[pos]+= ' (' + big_int + ')rop['+ str(i) + '] * polys_P[0]['+ str(j) + '] +'
			else:
				tmp_part2[pos%n]+= ' (' + big_int + ')rop['+ str(i) + '] * polys_P[0]['+ str(j) + '] +'
	
	for i in range(n-1):
		tmp_part1[i] = tmp_part1[i][:-2]
		if abs_c == 1 :
			tmp_part2[i] = ' ' + c_sign + ' (' + tmp_part2[i][1:-2] +')'
		elif abs_c_is_2pow :
			tmp_part2[i] = ' ' + c_sign + ' ((' + tmp_part2[i][1:-2] +') << ' + str(int(log(abs_c, 2))) + ')'
		else :
			tmp_part2[i] = ' ' + c_sign + ' ((' + tmp_part2[i][1:-2] +') * ' + str(abs_c) + ')'
	
	tmp_part1[n-1] = tmp_part1[n-1][:-2]
	
	result = '\n'
	for i in range(n-1):
		result += "	" + tmp_part1[i] + tmp_part2[i] + ';\n'
	result +=  "	" + tmp_part1[n-1] + ';\n'
	
	return result

def build_square_code(n, lambd, big_int):
	
	c_sign = get_sign(lambd)
	abs_c = abs(lambd)
	abs_c_is_2pow = is_power_of_two(abs_c)
	
	tmp_part1 = ['']*n
	tmp_part3 = ['']*n
	tmp_part2 = ['']*(n-1)
	
	for i in range(n):
		tmp_part1[i] = 'tmp_prod_result['+str(i)+'] = '
	
	for i in range(n):
		for j in range(i):
			pos = i+j
			if pos < n :
				tmp_part3[pos] += ' (' + big_int + ')pa['+ str(i) + '] * pa['+ str(j) + '] +'
			else:
				tmp_part2[pos%n] += ' (' + big_int + ')pa['+ str(i) + '] * pa['+ str(j) + '] +'
	
	for i in range(n):
		if tmp_part3[i] != '' :
			tmp_part3[i] = '((' + tmp_part3[i][1:-2] + ') << 1)'
	
	for i in range(n, 2*n-1):
		if tmp_part2[i%n] != '' :
			if i%2 == 0 :
				tmp_part2[i%n] = '((' + tmp_part2[i%n][1:-2] + ') << 1)'
			else :
				tmp_part2[i%n] = tmp_part2[i%n][1:-2]
	
	for i in range(n):
		pos = i+i
		if pos < n :
			if tmp_part3[pos] != '' :
				tmp_part3[pos] += ' + (' + big_int + ')pa['+ str(i) + '] * pa['+ str(i) + ']'
			else :
				tmp_part3[pos] += '(' + big_int + ')pa['+ str(i) + '] * pa['+ str(i) + ']'
		else:
			if tmp_part2[pos%n] != '' :
				tmp_part2[pos%n] += ' + (' + big_int + ')pa['+ str(i) + '] * pa['+ str(i) + ']'
			else :
				tmp_part2[pos%n] += '(' + big_int + ')pa['+ str(i) + '] * pa['+ str(i) + ']'
	
	for i in range(n, 2*n-1):
		if i%2 == 0 :
			if abs_c == 1 :
				tmp_part2[i%n] = ' ' + c_sign + ' (' + tmp_part2[i%n] +')'
			elif abs_c_is_2pow :
				tmp_part2[i%n] = ' ' + c_sign + ' ((' + tmp_part2[i%n] +') << ' + str(int(log(abs_c, 2))) + ')'
			else :
				tmp_part2[i%n] = ' ' + c_sign + ' ((' + tmp_part2[i%n] +') * '  + str(abs_c) + ')'
		else :
			if abs_c_is_2pow :
				tmp_part2[i%n] = ' ' + c_sign + ' ((' + tmp_part2[i%n] +') << ' + str(int(log(abs_c, 2))+1) + ')'
			else :
				tmp_part2[i%n] = ' ' + c_sign + ' ((' + tmp_part2[i%n] +') * '  + str(abs_c*2) + ')'
	
	result = '\n'
	for i in range(n-1):
		result += "	" + tmp_part1[i]+ tmp_part3[i] + tmp_part2[i] + ';\n'
	result += "	" + tmp_part1[n-1]+ tmp_part3[n-1] + ';\n'
	
	return result


#~ --------------------------------------- for any unitary E(X)  ---------------------------------------------------------

def build_prod_code_v2(n, upPow_distrib, big_int):

	tmp_part1 = ['']*n
	tmp_part2 = ['']*(n-1)
	
	for i in range(n):
		tmp_part1[i] = 'tmp_prod_result['+str(i)+'] ='
	
	for i in range(n):
		for j in range(n):
			pos = i+j
			if pos < n :
				tmp_part1[pos]+= ' (' + big_int + ')pa['+ str(i) + '] * pb['+ str(j) + '] +'
			else:
				tmp_part2[pos%n]+= ' (' + big_int + ')pa['+ str(i) + '] * pb['+ str(j) + '] +'
	
	for i in range(n-1):
		tmp_part1[i] = tmp_part1[i][:-2]
		tmp_part2[i] = tmp_part2[i][1:-2]
	tmp_part1[n-1] = tmp_part1[n-1][:-2]
	
	
	result = "	" + big_int
	for i in range(n-2):
		result += " c" + str(n+i) + ","
	result += " c" + str(2*n-2) + ";\n\n"
	
	
	for i in range(n-1):
		result += "	c" + str(n+i) + " = " + tmp_part2[i] + ';\n'
	result += '\n'	
	
	
	upPow_contrib = build_scal_prod(upPow_distrib, n, "prd")
	for i in range(n):
		result += "	" + tmp_part1[i] + upPow_contrib[i] + ';\n'
	
	return result


def build_exactRedCoeff_interProd_code_v2(n, upPow_distrib, big_int):

	tmp_part1 = ['']*n
	tmp_part2 = ['']*(n-1)
	
	for i in range(n):
		tmp_part1[i] = 'tmp['+str(i)+'] ='
	
	for i in range(n):
		for j in range(n):
			pos = i+j
			if pos < n :
				tmp_part1[pos]+= ' (' + big_int + ')rop['+ str(i) + '] * polys_P[0]['+ str(j) + '] +'
			else:
				tmp_part2[pos%n]+= ' (' + big_int + ')rop['+ str(i) + '] * polys_P[0]['+ str(j) + '] +'
	
	for i in range(n-1):
		tmp_part1[i] = tmp_part1[i][:-2]
		tmp_part2[i] = tmp_part2[i][1:-2]
	tmp_part1[n-1] = tmp_part1[n-1][:-2]
	
	
	result = "	" + big_int
	for i in range(n-2):
		result += " c" + str(n+i) + ","
	result += " c" + str(2*n-2) + ";\n\n"
	
	
	for i in range(n-1):
		result += "	c" + str(n+i) + " = " + tmp_part2[i] + ';\n'
	result += '\n'	
	
	
	upPow_contrib = build_scal_prod(upPow_distrib, n, "prd")
	for i in range(n):
		result += "	" + tmp_part1[i] + upPow_contrib[i] + ';\n'
	
	return result


def build_square_code_v2(n, upPow_distrib, big_int):
	
	tmp_part1 = ['']*n
	tmp_part3 = ['']*n
	tmp_part2 = ['']*(n-1)
	
	for i in range(n):
		tmp_part1[i] = 'tmp_prod_result['+str(i)+'] = '
	
	for i in range(n):
		for j in range(i):
			pos = i+j
			if pos < n :
				tmp_part3[pos] += ' (' + big_int + ')pa['+ str(i) + '] * pa['+ str(j) + '] +'
			else:
				tmp_part2[pos%n] += ' (' + big_int + ')pa['+ str(i) + '] * pa['+ str(j) + '] +'
	
	for i in range(n):
		if tmp_part3[i] != '' :
			tmp_part3[i] = '((' + tmp_part3[i][1:-2] + ') << 1)'
	
	for i in range(n, 2*n-1):
		if tmp_part2[i%n] != '' :
			tmp_part2[i%n] = '((' + tmp_part2[i%n][1:-2] + ') << 1)'
	
	for i in range(n):
		pos = i+i
		if pos < n :
			if tmp_part3[pos] != '' :
				tmp_part3[pos] += ' + (' + big_int + ')pa['+ str(i) + '] * pa['+ str(i) + ']'
			else :
				tmp_part3[pos] = '(' + big_int + ')pa['+ str(i) + '] * pa['+ str(i) + ']'
		else:
			if tmp_part2[pos%n] != '' :
				tmp_part2[pos%n] += ' + (' + big_int + ')pa['+ str(i) + '] * pa['+ str(i) + ']'
			else :
				tmp_part2[pos%n] = '(' + big_int + ')pa['+ str(i) + '] * pa['+ str(i) + ']'
	
	result = "	" + big_int
	for i in range(n-2):
		result += " c" + str(n+i) + ","
	result += " c" + str(2*n-2) + ";\n\n"
	
	for i in range(n-1):
		result += "	c" + str(n+i) + " = " + tmp_part2[i] + ';\n'
	result += '\n'
	
	upPow_contrib = build_scal_prod(upPow_distrib, n, "sqr")
	for i in range(n):
		result += "	" + tmp_part1[i] + tmp_part3[i] + upPow_contrib[i] + ';\n'
	
	return result



#~ --------------------------------------- internal reduction methods  ---------------------------------------------------------


def build_red_int_code(unsigned_small_int, big_int, n, lambd, mont_phi, red_int_coeff, neg_inv_ri_rep_coeff, mask_for_redint):
	
	lambda_mat = build_lambda_mat(lambd, n-1)
	
	return build_red_int_code_v2(unsigned_small_int, big_int, n, lambda_mat, mont_phi, red_int_coeff, neg_inv_ri_rep_coeff, mask_for_redint)


def build_red_int_code_v2(unsigned_small_int, big_int, n, upPow_distrib, mont_phi, red_int_coeff, neg_inv_ri_rep_coeff, mask_for_redint):
	
	result = "	" + unsigned_small_int
	result += " tmpQ[" + str(n) + "];\n"
	
	result += "	" + big_int
	result += " tmpZero[" + str(n) + "];\n\n"
	
	result += "	//~ computation of : op*neg_inv_ri_rep_coeff mod(E, mont_phi)\n"
	result += build_redInt_tmpQ_prod_code_v2(n, mont_phi, upPow_distrib, neg_inv_ri_rep_coeff, unsigned_small_int, mask_for_redint)
	
	result += "\n	//~ computation of : tmp_q*red_int_coeff mod(E)\n"
	result += build_redInt_tmpZero_prod_code_v2(n, upPow_distrib, red_int_coeff, big_int)
	
	tmp = ''
	result += "\n	//~ computation of : (op + tmp_zero)/mont_phi\n"
	for i in range(n):
		tmp = "	rop[" + str(i) + "] = (op[" + str(i) + "] + tmpZero[" + str(i) + "]) >> PHI_LOG2;\n"
		result += tmp
	
	return result


#~ note : elements of 'neg_inv_ri_rep_coeff' are >= 0
def build_redInt_tmpQ_prod_code_v2(n, mont_phi, upPow_distrib, neg_inv_ri_rep_coeff, unsigned_small_int, mask_for_redint):
	
	tmp_partStart = ['']*n
	for i in range(n):
		tmp_partStart[i] = 'tmpQ['+str(i)+'] = '
	
	#~ some adjustments, if needed
	neg_inv_ri_rep_coeff += [0] * (n - len(neg_inv_ri_rep_coeff))
	
	low_prod = build_itl_mat(neg_inv_ri_rep_coeff)
	
	#~ accumulation of products
	ivm_prodMat = list(build_ith_mat(neg_inv_ri_rep_coeff[1:])*upPow_distrib)
	ivm_prodMat.insert(0, [0]*n)
	ivm_prodMat = matrix(ivm_prodMat)
	
	sum_mat = list(low_prod + ivm_prodMat)
	for i in range(n):
		for j in range(n):
			sum_mat[i][j] = sum_mat[i][j] % mont_phi
	sum_mat = matrix(sum_mat)
	
	prod_list = build_scal_prod(sum_mat, 1, "tmpQ", unsigned_small_int=unsigned_small_int)
	
	(mask_p1, mask_p2) = ('(', ') & REDINT_MASK') if mask_for_redint != 0 else ('','')
	
	result = ""
	for i in range(n):
		if prod_list[i] == '' :
			result += "	" + tmp_partStart[i] + '0;\n'
		else :
			result += "	" + tmp_partStart[i] + mask_p1 + prod_list[i][3:] + mask_p2 + ';\n'
		
	return result


# IMPORTANT : 'big_int' is assumed big enough so that overflow will not happen on target architecture.
def build_redInt_tmpZero_prod_code_v2(n, upPow_distrib, red_int_coeff, big_int):
	
	tmp_partStart = ['']*n
	for i in range(n):
		tmp_partStart[i] = 'tmpZero['+str(i)+'] = '
	
	#~ some adjustments, if needed
	red_int_coeff += [0] * (n - len(red_int_coeff))
	
	low_prod = build_itl_mat(red_int_coeff)
	
	#~ accumulation of products
	m_prodMat = list(build_ith_mat(red_int_coeff[1:])*upPow_distrib)
	m_prodMat.insert(0, [0]*n)
	m_prodMat = matrix(m_prodMat)
	
	sum_mat = low_prod + m_prodMat
	
	prod_list = build_scal_prod(sum_mat, 1, "tmpZero", big_int=big_int)
	
	result = ""
	for i in range(n):
		if prod_list[i] == '' :
			result += "	" + tmp_partStart[i] + '0;\n'
		else :
			if prod_list[i][1] == '-':
				result += "	" + tmp_partStart[i] + prod_list[i][1:] + ';\n'
			else:
				result += "	" + tmp_partStart[i] + prod_list[i][3:] + ';\n'
	
	return result


#~ --------------------------------------- conversion polynomials P  ---------------------------------------------------------

def compute_redExtPol_w(n, redExtPol, R):
	
	x = R.gen()
	
	c = n-1
	tmp = x**n
	l1 = list(tmp%redExtPol)
	l1 = [abs(k) for k in l1]
	V = c * vector(l1 + [0]*(n-len(l1)))
	for d in range(n-2):
		c -= 1
		tmp *= x
		l1 = list(tmp%redExtPol)
		l1 = [abs(k) for k in l1]
		V += c * vector(l1 + [0]*(n-len(l1)))
	
	V += vector(range(1, n+1))
	
	return max(V)

#~ returns a representation of 'op/phi'
def amns_red_int(op, ext_pol, ri_poly, neg_inv_ri, R, PP, phi):
	q = (PP(op)*neg_inv_ri).mod(PP(ext_pol))
	r0 = op + (R(q)*ri_poly).mod(ext_pol)
	return (r0/phi)


#~ Uses method 2 of conversion to AMNS, to compute a representation of 'val', with infinity norm lower than 'rho_min'
def compute_rep_in_amns(val, n, p, gmm, phi, ext_pol, ri_poly, neg_inv_ri, F, R, PP):
	
	w = compute_redExtPol_w(n, ext_pol, R)
	rho_min = 2 * w * ri_poly.norm(infinity)
	
	tmp_rep = Integer(F(val * phi.powermod(n-1, p)))
	
	rep = amns_red_int(R(tmp_rep), ext_pol, ri_poly, neg_inv_ri, R, PP, phi)
	for i in range(n-2):
		rep = amns_red_int(rep, ext_pol, ri_poly, neg_inv_ri, R, PP, phi)
	
	if F(rep(gmm)) != F(val) :
		print("ERROR : Bad conversion !!!")
	
	if rep.norm(infinity) >= rho_min :
		print("ERROR : Element has infinity norm greater than expected !!!")
	
	return rep


#~ computes the polynomials Pi, for conversion: Pi ~ (beta^i * phi^2)
def compute_conv_polysP(base_beta, n, ext_pol, ri_poly, neg_inv_ri, phi, gmm, p, F, R, PP):
	
	phi2 = phi**2
	
	conv_polys = []
	for i in range(n):
		poly = compute_rep_in_amns(phi2*(base_beta**i), n, p, gmm, phi, ext_pol, ri_poly, neg_inv_ri, F, R, PP)
		conv_polys.append(poly)
	
	return conv_polys


def generate_conv_polys_P(small_int_name, base_beta, n, ext_pol, ri_poly, neg_inv_ri, phi, gmm, p, F, R, PP):
	
	conv_polys = compute_conv_polysP(base_beta, n, ext_pol, ri_poly, neg_inv_ri, phi, gmm, p, F, R, PP)
	
	result = small_int_name + " polys_P[NB_COEFF][NB_COEFF] = {\n"
	
	# ~ for i in range(n-1):
		# ~ result += "	{" + str(list(conv_polys[i]))[1:-1] + "},\n"
	# ~ result += "	{" + str(list(conv_polys[n-1]))[1:-1] + "}};\n\n"
	
	for i in range(n-1):
		result += "	{" + list__from_int_to_str_hex(list(conv_polys[i])) + "},\n"
	result += "	{" + list__from_int_to_str_hex(list(conv_polys[n-1])) + "}};\n"
	
	return result 

def list__from_int_to_str_hex(llist, el_type='L'):
	res = ''
	for el in llist[:-1]:
		res += hex(el) + el_type + ', '
	res += hex(llist[-1]) + el_type
	return res












