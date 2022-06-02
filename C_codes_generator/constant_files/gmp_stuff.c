//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/////// GNU MP elements for Montgomery modular multiplication //////
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#define mpn_redc_1 __MPN(redc_1)
__GMP_DECLSPEC mp_limb_t mpn_redc_1 (mp_ptr, mp_ptr, mp_srcptr, mp_size_t, mp_limb_t);

#define mpn_redc_n __MPN(redc_n)
__GMP_DECLSPEC void mpn_redc_n (mp_ptr, mp_ptr, mp_srcptr, mp_size_t, mp_srcptr);

#define   mpn_binvert __MPN(binvert)
__GMP_DECLSPEC void mpn_binvert (mp_ptr, mp_srcptr, mp_size_t, mp_ptr);

#define binvert_limb_table  __gmp_binvert_limb_table
__GMP_DECLSPEC extern const unsigned char  binvert_limb_table[128];
#define binvert_limb(inv,n)						\
  do {									\
    mp_limb_t  __n = (n);						\
    mp_limb_t  __inv;							\
									\
    __inv = binvert_limb_table[(__n/2) & 0x7F]; /*  8 */		\
    if (GMP_NUMB_BITS > 8)   __inv = 2 * __inv - __inv * __inv * __n;	\
    if (GMP_NUMB_BITS > 16)  __inv = 2 * __inv - __inv * __inv * __n;	\
    if (GMP_NUMB_BITS > 32)  __inv = 2 * __inv - __inv * __inv * __n;	\
									\
    if (GMP_NUMB_BITS > 64)						\
      {									\
	int  __invbits = 64;						\
	do {								\
	  __inv = 2 * __inv - __inv * __inv * __n;			\
	  __invbits *= 2;						\
	} while (__invbits < GMP_NUMB_BITS);				\
      }									\
									\
    (inv) = __inv & GMP_NUMB_MASK;					\
  } while (0)

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void mpn_mod_mult(mp_limb_t *rop, const mp_limb_t *op1, const mp_limb_t *op2, const mp_limb_t *p_limbs, int nb_limbs){
	
	mp_limb_t c_limbs[2*nb_limbs];
	mp_limb_t q_limbs[nb_limbs+1];
	
	mpn_mul_n (c_limbs, op1, op2, nb_limbs); 
	mpn_tdiv_qr (q_limbs, rop, 0, c_limbs, (nb_limbs*2), p_limbs, nb_limbs); 
}

void mpn_mont_mul_red_1(mp_limb_t *rop, mp_limb_t *op1, mp_limb_t *op2, const mp_limb_t *p_limbs, mp_limb_t mip0, int nb_limbs){
	
	mp_limb_t tmp_limbs[2*nb_limbs];
	
	mpn_mul_n (tmp_limbs, op1, op2, nb_limbs);
	
	if (mpn_redc_1(rop, tmp_limbs, p_limbs, nb_limbs, mip0) != 0)
		mpn_sub_n (rop, rop, p_limbs, nb_limbs);
}

void mpn_mont_mul_red_n(mp_limb_t *rop, mp_limb_t *op1, mp_limb_t *op2, const mp_limb_t *p_limbs, mp_limb_t *mip_limbs, int nb_limbs){
	
	mp_limb_t tmp_limbs[2*nb_limbs];
	
	mpn_mul_n (tmp_limbs, op1, op2, nb_limbs);
	mpn_redc_n (rop, tmp_limbs, p_limbs, nb_limbs, mip_limbs);
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
////////////////////////////////////////////////////////////////////
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//~ assumes 'rop' already initialised
void from_limbs_to_mpz_t(mpz_t rop, mp_limb_t *op, int op__nb_limbs){
	
	mpz_set_ui(rop, 0);
	
	if (op__nb_limbs == 0)
		return;
	
	int i, nb_dec;
	
	nb_dec = 8*sizeof(mp_limb_t);
	
	mpz_set_ui(rop, op[op__nb_limbs - 1]);
	
	for(i=(op__nb_limbs - 2); i >=0; i--){
		mpz_mul_2exp (rop, rop, nb_dec);
		mpz_add_ui (rop, rop, op[i]);
	}
}


//~ important: assumes allocation of 'nb_limbs' limbs already done for 'rop'
void copy_limbs(mp_limb_t *rop, mpz_t op, int nb_limbs){
	int i;
	const mp_limb_t *tmp;
	tmp = mpz_limbs_read(op);
	for(i=0; i<nb_limbs; i++)
		rop[i] = tmp[i];
}

void clean_limbs(mp_limb_t *op, int nb_limbs){
	for(int i=0; i<nb_limbs; i++)
		op[i] = 0;
}



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
////////////////////////////////////////////////////////////////////
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
