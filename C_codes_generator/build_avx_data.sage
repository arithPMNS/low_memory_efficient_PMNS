def build_avx_const_data(dir_path, n, ri_mat, iri_mat, blc_size=8):
	
	poly_deg = n - 1
	
	ref_kq = "#define KQ_"
	ref_kz = "#define KZ_"
	ref_type = " (__m512i){"
	
	with open(dir_path+"/avx_const_data.h", "w") as f:
		
		if n <= blc_size:
			
			for i in range(n):
				l = ref_kq + str(i) + ref_type
				t = ref_kz + str(i) + ref_type
				for j in range(poly_deg):
					l += str(hex(iri_mat[i][j]))+"UL,"
					t += str(hex(ri_mat[i][j]))+"L,"
				l += str(hex(iri_mat[i][poly_deg]))+"UL,}"
				t += str(hex(ri_mat[i][poly_deg]))+"L,}"
				f.write(l+"\n")
				f.write(t+"\n")
		
		else:
			
			nl = floor(n/blc_size)
			rmd = n%blc_size
			
			#print(nl,rmd)
			
			for i in range(n):
				for j in range(nl):
					tmp_kq = ref_kq + str(i) + "_" + str(j) + ref_type
					tmp_kz = ref_kz + str(i) + "_" + str(j) + ref_type
					for k in range(blc_size-1):
						tmp_kq += str(hex(iri_mat[i][blc_size*j+k]))+"UL,"
						tmp_kz += str(hex(ri_mat[i][blc_size*j+k]))+"L,"
					tmp_kq += str(hex(iri_mat[i][blc_size*(j+1)-1]))+"UL}"
					tmp_kz += str(hex(ri_mat[i][blc_size*(j+1)-1]))+"L}"
					f.write(tmp_kq+"\n")
					f.write(tmp_kz+"\n")
				if rmd != 0:
					j = nl
					tmp_kq = ref_kq + str(i) + "_" + str(j) + ref_type
					tmp_kz = ref_kz + str(i) + "_" + str(j) + ref_type
					for k in range(rmd-1):
						tmp_kq += str(hex(iri_mat[i][blc_size*j+k]))+"UL,"
						tmp_kz += str(hex(ri_mat[i][blc_size*j+k]))+"L,"
					tmp_kq += str(hex(iri_mat[i][blc_size*j+rmd-1]))+"UL,}"
					tmp_kz += str(hex(ri_mat[i][blc_size*j+rmd-1]))+"L,}"
					f.write(tmp_kq+"\n")
					f.write(tmp_kz+"\n")
	return;




